////////////////////////////////////////////////////////////////////////
//
// CalWireT1053 class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//  copied over to 1053 - andrzej.szelc@yale.edu
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <stdint.h>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

#include "lar1ndcode/Utilities/SignalShapingServiceT1053.h"
#include "larcore/Geometry/Geometry.h"
//#include "Filters/ChannelFilter.h"

#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "lardata/RecoBase/Wire.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/sparse_vector.h"
#include "lardata/Utilities/AssociationUtil.h"  //--Hec

#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

///creation of calibrated signals on wires
namespace caldata {

  class CalWireT1053 : public art::EDProducer {

  public:
    
    typedef lar::sparse_vector<float> RegionsOfInterest_t;

    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireT1053(fhicl::ParameterSet const& pset); 
    virtual ~CalWireT1053();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    
    int          fDataSize;          ///< size of raw data on one wire
    int          fPostsample;        ///< number of postsample bins
    int          fBaseSampleBins;        ///< number of postsample bins
    float        fBaseVarCut;        ///< baseline variance cut
    std::string  fDigitModuleLabel;  ///< module that made digits
                                                       ///< constants
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

    void SubtractBaseline(std::vector<float>& holder, int fBaseSampleBins);

  protected: 
    
  }; // class CalWireT1053

  DEFINE_ART_MODULE(CalWireT1053)
  
  //-------------------------------------------------
  CalWireT1053::CalWireT1053(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);

    //--Hec if(fSpillName.size()<1) produces< std::vector<recob::Wire> >();
    //--Hec else produces< std::vector<recob::Wire> >(fSpillName);
  
    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
  }
  //-------------------------------------------------
  CalWireT1053::~CalWireT1053()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireT1053::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fPostsample       = p.get< int >        ("PostsampleBins");
    fBaseSampleBins   = p.get< int >        ("BaseSampleBins");
    fBaseVarCut       = p.get< int >        ("BaseVarCut");
    
    fSpillName="";
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
  }

  //-------------------------------------------------
  void CalWireT1053::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void CalWireT1053::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWireT1053::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    int transformSize = fFFT->FFTSize();

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceT1053> sss;


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
    mf::LogInfo("CalWireT1053") << "CalWireT1053:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors

    if( (unsigned int)transformSize < dataSize){
      mf::LogWarning("CalWireT1053")<<"FFT size (" << transformSize << ") "
                                    << "is smaller than the data size (" << dataSize << ") "
                                    << "\nResizing the FFT now...";
      fFFT->ReinitializeFFT(dataSize,fFFT->FFTOptions(),fFFT->FFTFitBins());
      transformSize = fFFT->FFTSize();
      mf::LogWarning("CalWireT1053")<<"FFT size is now (" << transformSize << ") "
                                    << "and should be larger than the data size (" << dataSize << ")";
    }

    mf::LogInfo("CalWireT1053") << "Data size is " << dataSize << " and transform size is " << transformSize;

    if(fBaseSampleBins > 0 && dataSize % fBaseSampleBins != 0) {
      mf::LogError("CalWireT1053")<<"Set BaseSampleBins modulo dataSize= "<<dataSize;
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

        // Do deconvolution.
        sss->Deconvolute(channel, holder);

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
      if(fBaseSampleBins) SubtractBaseline(holder, fBaseSampleBins);

      // Make a single ROI that spans the entire data size
      RegionsOfInterest_t sparse_holder;
      sparse_holder.add_range(0,holder.begin(),holder.end());
      auto view = geom->View(digitVec->Channel());
      wirecol->emplace_back(sparse_holder,digitVec->Channel(),view);
    
      // add an association between the last object in wirecol--Hec
      // (that we just inserted) and digitVec
      if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
        throw art::Exception(art::errors::InsertFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
      } // if failed to add association
    }


    if(wirecol->size() == 0)
      mf::LogWarning("CalWireT1053") << "No wires made for this event.";

    //--Hec if(fSpillName.size()>0)
    //--Hec   evt.put(std::move(wirecol), fSpillName);
    //--Hec else evt.put(std::move(wirecol));

    evt.put(std::move(wirecol), fSpillName);        //--Hec
    evt.put(std::move(WireDigitAssn), fSpillName);  //--Hec
    
   // delete chanFilt;
    return;
  }
  
  void CalWireT1053::SubtractBaseline(std::vector<float>& holder,
     int fBaseSampleBins)
  {
    // subtract baseline using linear interpolation between regions defined
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
  } // SubtractBaseline
  
} // end namespace caldata
