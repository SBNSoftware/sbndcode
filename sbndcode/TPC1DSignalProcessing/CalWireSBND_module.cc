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
#include "art/Utilities/make_tool.h"
#include "lardata/ArtDataHelper/WireCreator.h"

#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "sbndcode/TPC1DSignalProcessing/IROIFinder.h"
#include "larcore/Geometry/Geometry.h"
//#include "Filters/ChannelFilter.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/Utilities/sparse_vector.h"
#include "lardata/Utilities/AssociationUtil.h"  //--Hec
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TComplex.h"
#include "TFile.h"
#include "TH1F.h"

///creation of calibrated signals on wires
namespace caldata {

  class CalWireSBND : public art::EDProducer {

  public:
    
    //typedef lar::sparse_vector<float> RegionsOfInterest_t;

    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireSBND(fhicl::ParameterSet const& pset); 
    virtual ~CalWireSBND();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);

    std::unique_ptr<sbnd_tool::IROIFinder> fROITool;

    // Define the waveform container
    using Waveform        = std::vector<float>;
    
    // Define the ROI and its container
    using CandidateROI    = std::pair<size_t, size_t>;
    using CandidateROIVec = std::vector<CandidateROI>;

  private:
    
    bool          fDoBaselineSub;     ///< subtract baseline to restore DC component post-deconvolution
    bool          fDoAdvBaselineSub;  ///< use interpolation-based baseline subtraction
    int           fBaseSampleBins;    ///< bin grouping size in "interpolate"  method
    float         fBaseVarCut;        ///< baseline variance cut used in "interpolate" method
    //FFT service options
    int fFFTSize;
    std::string fFFTOption;
    int fFFTFitBins;
   
    std::string  fDigitModuleLabel;   ///< module that made digits
                                                       
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data
    
    void          SubtractBaseline(std::vector<float>& holder);
    void          SubtractBaselineAdv(std::vector<float>& holder);
    

  protected: 
    
  }; // class CalWireSBND

  DEFINE_ART_MODULE(CalWireSBND)
  
  //-------------------------------------------------
  CalWireSBND::CalWireSBND(fhicl::ParameterSet const& pset)
  : EDProducer(pset)
  {
    fSpillName="";
    this->reconfigure(pset);
    
    fROITool = art::make_tool<sbnd_tool::IROIFinder>(pset.get<fhicl::ParameterSet>("ROITool"));

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
    fDoBaselineSub    = p.get< bool >       ("DoBaselineSub");
    fDoAdvBaselineSub = p.get< bool >       ("DoAdvBaselineSub");
    fBaseSampleBins   = p.get< int >        ("BaseSampleBins");
    fBaseVarCut       = p.get< int >        ("BaseVarCut");
    fFFTSize          = p.get< int >        ("FFTSize");
    fFFTOption        = p.get< std::string >("FFTOption");
    fFFTFitBins       = p.get< int >        ("FFTFitBins");
    
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

    // reset FFT service for each event
    fFFT->ReinitializeFFT(fFFTSize,fFFTOption,fFFTFitBins);

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
      mf::LogInfo("CalWireSBND")<<"FFT size (" << transformSize << ") "
                                    << "is smaller than the data size (" << dataSize << ") "
                                    << "\nResizing the FFT now...";
      fFFT->ReinitializeFFT(dataSize,fFFT->FFTOptions(),fFFT->FFTFitBins());
      transformSize = fFFT->FFTSize();
      mf::LogInfo("CalWireSBND")<<"FFT size is now (" << transformSize << ") "
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
    
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

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
        sss->Deconvolute(clockData, channel, holder);
	  for(bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;
      } // end if not a bad channel 
      
      holder.resize(dataSize,1e-5);

      // restore DC component through baseline subtraction
      if( fDoBaselineSub ) SubtractBaseline(holder);
      // more advanced, interpolation-based subtraction alg 
      // that uses the BaseSampleBins and BaseVarCut params
      if( fDoAdvBaselineSub ) SubtractBaselineAdv(holder);

      // Make a single ROI that spans the entire data size
      //RegionsOfInterest_t sparse_holder;
      //sparse_holder.add_range(0,holder.begin(),holder.end());
      CandidateROIVec candROIVec;
      fROITool->FindROIs( holder, channel, candROIVec);//calculates ROI and returns it to roiVec.
      //auto view = geom->View(channel);
      recob::Wire::RegionsOfInterest_t roiVec;

      //looping over roiVec to make a RegionOfInterest_t object.
      for(auto const& CandidateROI: candROIVec){
	size_t roiStart = CandidateROI.first;
	size_t roiStop = CandidateROI.second;
	std::vector<float> roiHolder;
	for(size_t i_holder=roiStart; i_holder<=roiStop; i_holder++){
	  roiHolder.push_back(holder[i_holder]);
	}
	roiVec.add_range(roiStart, std::move(roiHolder));
      }
      wirecol->push_back(recob::WireCreator(std::move(roiVec),*digitVec).move());


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
    // Robust baseline calculation that effectively ignores outlier 
    // samples from large pulses:
    //   (1) fill a histogram with every sample's value,
    //   (2) find mode (bin with most entries),
    //   (3) calculate the mean along the entire waveform using
    //       only samples with values close to this mode.
    unsigned int bin(0);  
    float min = 0, max = 0;
    for(bin = 0; bin < holder.size(); bin++){
      if (holder[bin] > max) max = holder[bin];
      if (holder[bin] < min) min = holder[bin];
    }
    int nbin = max - min;
    if (nbin > 0) {
      TH1F h("h","h",nbin,min,max);
      for(bin = 0; bin < holder.size(); bin++) h.Fill(holder[bin]);
      float x_max = h.GetXaxis()->GetBinCenter(h.GetMaximumBin());
      float ped   = x_max;
      float sum   = 0;
      int ncount  = 0;
      for(bin = 0; bin < holder.size(); bin++){
        if( fabs(holder[bin]-x_max) < 2. ) {
          sum += holder[bin];
          ncount++; 
        }
      }
      if (ncount) ped = sum/ncount;
      for(bin = 0; bin < holder.size(); bin++) holder[bin] -= ped;
    }
  }
 
  void CalWireSBND::SubtractBaselineAdv(std::vector<float>& holder)
  {
      // Subtract baseline using linear interpolation between regions defined
      // by the datasize and fBaseSampleBins

      // number of points to characterize the baseline
      unsigned short nBasePts = holder.size() / fBaseSampleBins;

      // the baseline offset vector
      std::vector<float> base;
      for(unsigned short ii = 0; ii < nBasePts; ++ii) base.push_back(0.);
      // find the average value in each region, using values that are
      // similar
      float fbins = fBaseSampleBins;
      unsigned short nfilld = 0;
      for(unsigned short ii = 0; ii < nBasePts; ++ii) {
        unsigned short loBin = ii * fBaseSampleBins;
        unsigned short hiBin = loBin + fBaseSampleBins;
        float ave = 0.;
        float sum = 0.;
        for(unsigned short bin = loBin; bin < hiBin; ++bin) {
          ave += holder[bin];
          sum += holder[bin] * holder[bin];
        } // jj
        ave = ave / fbins;
        float var = (sum - fbins * ave * ave) / (fbins - 1.);
        // Set the baseline for this region if the variance is small
        if(var < fBaseVarCut) {
          base[ii] = ave;
          ++nfilld;
        }
      } // ii
      // fill in any missing points if there aren't too many missing
      if(nfilld < nBasePts && nfilld > nBasePts / 2) {
        bool baseOK = true;
        // check the first region
        if(base[0] == 0) {
          unsigned short ii1 = 0;
          for(unsigned short ii = 1; ii < nBasePts; ++ii) {
            if(base[ii] != 0) {
              ii1 = ii;
              break;
            }
          } // ii
          unsigned short ii2 = 0;
          for(unsigned short ii = ii1 + 1; ii < nBasePts; ++ii) {
            if(base[ii] != 0) {
              ii2 = ii;
              break;
            }
          } // ii
          // failure
          if(ii2 > 0) {
            float slp = (base[ii2] - base[ii1]) / (float)(ii2 - ii1);
            base[0] = base[ii1] - slp * ii1;
          } else {
            baseOK = false;
          }
        } // base[0] == 0
        // check the last region
        if(baseOK && base[nBasePts] == 0) {
          unsigned short ii2 = 0;
          for(unsigned short ii = nBasePts - 1; ii > 0; --ii) {
            if(base[ii] != 0) {
              ii2 = ii;
              break;
            }
          } // ii
          baseOK = false; // assume the worst, hope for better
          unsigned short ii1 = 0;
          if (ii2 >= 1) {
            for(unsigned short ii = ii2 - 1; ii > 0; --ii) {
              if(base[ii] != 0) {
                ii1 = ii;
                baseOK = true;
                break;
              } // if base[ii]
            } // ii
          } // if ii2
          if (baseOK) {
            float slp = (base[ii2] - base[ii1]) / (float)(ii2 - ii1);
            base[nBasePts] = base[ii2] + slp * (nBasePts - ii2);
          }
        } // baseOK && base[nBasePts] == 0
        // now fill in any intermediate points
        for(unsigned short ii = 1; ii < nBasePts - 1; ++ii) {
          if(base[ii] == 0) {
            // find the next non-zero region
            for(unsigned short jj = ii + 1; jj < nBasePts; ++jj) {
              if(base[jj] != 0) {
                float slp = (base[jj] - base[ii - 1]) / (jj - ii + 1);
                base[ii] = base[ii - 1] + slp;
                break;
              }
            } // jj
          } // base[ii] == 0
        } // ii
      } // nfilld < nBasePts

      // interpolate and subtract
      float slp = (base[1] - base[0]) / (float)fBaseSampleBins;
      // bin offset to the origin (the center of the region)
      unsigned short bof = fBaseSampleBins / 2;
      unsigned short lastRegion = 0;
      for(unsigned short bin = 0; bin < holder.size(); ++bin) {
        // in a new region?
        unsigned short region = bin / fBaseSampleBins;
        if(region > lastRegion) {
          // update the slope and offset
          slp = (base[region] - base[lastRegion]) / (float)fBaseSampleBins;
          bof += fBaseSampleBins;
          lastRegion = region;
        }
        holder[bin] -= base[region] + (bin - bof) * slp;
      }
    }
  

} // end namespace caldata
