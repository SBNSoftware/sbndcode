////////////////////////////////////////////////////////////////////////
//
// CalWireSBND class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//  copied over to 1053 - andrzej.szelc@yale.edu
//  copied over and modified to SBND   
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <stdint.h>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "larcore/Geometry/Geometry.h"
//#include "Filters/ChannelFilter.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/Utilities/sparse_vector.h"
#include "lardata/Utilities/AssociationUtil.h"  //--Hec

#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

///creation of calibrated signals on wires
namespace caldata {

  class CalWireSBND : public art::EDProducer {

  public:
    
    typedef lar::sparse_vector<float> RegionsOfInterest_t;

    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireSBND(fhicl::ParameterSet const& pset); 
    virtual ~CalWireSBND();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    
  //  int          fDataSize;          ///< size of raw data on one wire
    int          fPostsample;        ///< number of postsample bins
    int          fBaseSampleBins;        ///< number of postsample bins
    bool         fDoBaselineSub;     ///< baseline subtraction
    float        fBaseVarCut;        ///< baseline variance cut
    std::string  fDigitModuleLabel;  ///< module that made digits
                                                       ///< constants
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

    //void SubtractBaseline(std::vector<float>& holder, int fBaseSampleBins);
    void SubtractBaseline(std::vector<float>& holder);

  protected: 
    
  }; // class CalWireSBND

  DEFINE_ART_MODULE(CalWireSBND)
  
  //-------------------------------------------------
  CalWireSBND::CalWireSBND(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);

    //--Hec if(fSpillName.size()<1) produces< std::vector<recob::Wire> >();
    //--Hec else produces< std::vector<recob::Wire> >(fSpillName);
  
    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
  }
  //-------------------------------------------------
  CalWireSBND::~CalWireSBND()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireSBND::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fPostsample       = p.get< int >        ("PostsampleBins");
    fBaseSampleBins   = p.get< int >        ("BaseSampleBins");
    fBaseVarCut       = p.get< int >        ("BaseVarCut");
    fDoBaselineSub    = p.get< bool >        ("DoBaselineSub");
    
    fSpillName="";
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
  }

  //-------------------------------------------------
  void CalWireSBND::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void CalWireSBND::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWireSBND::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    int transformSize = fFFT->FFTSize();

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceSBND> sss;
    double DeconNorm = sss->GetDeconNorm();

    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    
    // ... and an association set     --Hec
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("CalWireSBND") << "CalWireSBND:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors

    if( (unsigned int)transformSize < dataSize){
      mf::LogWarning("CalWireSBND")<<"FFT size (" << transformSize << ") "
                                    << "is smaller than the data size (" << dataSize << ") "
                                    << "\nResizing the FFT now...";
      fFFT->ReinitializeFFT(dataSize,fFFT->FFTOptions(),fFFT->FFTFitBins());
      transformSize = fFFT->FFTSize();
      mf::LogWarning("CalWireSBND")<<"FFT size is now (" << transformSize << ") "
                                    << "and should be larger than the data size (" << dataSize << ")";
    }

    mf::LogInfo("CalWireSBND") << "Data size is " << dataSize << " and transform size is " << transformSize;

    if(fBaseSampleBins > 0 && dataSize % fBaseSampleBins != 0) {
      mf::LogError("CalWireSBND")<<"Set BaseSampleBins modulo dataSize= "<<dataSize;
    }

    uint32_t     channel(0); // channel number
    unsigned int bin(0);     // time bin loop variable
    
///    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();  

    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    
    // loop over all wires    
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();

      // skip bad channels
      //  if(!chanFilt->BadChannel(channel)) {
      if(true) {

        // resize and pad with zeros
        holder.resize(transformSize, 0.);
        
        // uncompress the data
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        
        // loop over all adc values and subtract the pedestal
        float pdstl = digitVec->GetPedestal();
        
        for(bin = 0; bin < dataSize; ++bin) 
          holder[bin]=(rawadc[bin]-pdstl);

	//fill the remaining bin with data
	for(bin = dataSize; bin < holder.size(); bin++){
	  //  philosophy change - don't repeat data but instead fill extra space with zeros.
          //    not sure that one is better than the other.
	  //	  holder[bin] = (rawadc[bin-dataSize]-pdstl);
	  holder[bin] = 0.0;
	}

        // Do deconvolution.
        sss->Deconvolute(channel, holder);
	  for(bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;
      } // end if not a bad channel 
      
      holder.resize(dataSize,1e-5);

      //This restores the DC component to signal removed by the deconvolution.
      if(fPostsample) {
        float average=0.0;
        for(bin=0; bin < (unsigned short)fPostsample; ++bin) 
          average += holder[holder.size()-1+bin];
        average = average / (float)fPostsample;
        for(bin = 0; bin < holder.size(); ++bin) holder[bin]-=average;
      }  
      // adaptive baseline subtraction
      //if(fBaseSampleBins) SubtractBaseline(holder, fBaseSampleBins);
      if(fDoBaselineSub) SubtractBaseline(holder);

      // Make a single ROI that spans the entire data size
      RegionsOfInterest_t sparse_holder;
      sparse_holder.add_range(0,holder.begin(),holder.end());
      auto view = geom->View(digitVec->Channel());
      wirecol->emplace_back(sparse_holder,digitVec->Channel(),view);
    
      // add an association between the last object in wirecol--Hec
      // (that we just inserted) and digitVec
      if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
        throw cet::exception("CalWireSBND")
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key() << "\n";
      } // if failed to add association
    }


    if(wirecol->size() == 0)
      mf::LogWarning("CalWireSBND") << "No wires made for this event.";

    //--Hec if(fSpillName.size()>0)
    //--Hec   evt.put(std::move(wirecol), fSpillName);
    //--Hec else evt.put(std::move(wirecol));

    evt.put(std::move(wirecol), fSpillName);        //--Hec
    evt.put(std::move(WireDigitAssn), fSpillName);  //--Hec
    
   // delete chanFilt;
    return;
  }


  void CalWireSBND::SubtractBaseline(std::vector<float>& holder)
  {
    float min = 0, max = 0;
    for(unsigned int bin = 0; bin < holder.size(); bin++){
      if (holder[bin] > max) max = holder[bin];
      if (holder[bin] < min) min = holder[bin];
    }
    int nbin = max - min;
    if (nbin!=0) {
      TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
      for(unsigned int bin = 0; bin < holder.size(); bin++){
	h1->Fill(holder[bin]);
      }
      float ped = h1->GetMaximum();
      float ave = 0, ncount = 0;
      for(unsigned int bin = 0; bin < holder.size(); bin++){
	if (fabs(holder[bin]-ped) < 2){
	  ave += holder[bin];
	  ncount++;
	}
      }
      if (ncount==0) ncount = 1;
      ave = ave/ncount;
      for(unsigned int bin = 0; bin < holder.size(); bin++){
	holder[bin] -= ave;
      }
      h1->Delete();
    }
  }
  
  
} // end namespace caldata
