////////////////////////////////////////////////////////////////////////
// Class:       CALLOS
// Plugin Type: analyzer (Unknown Unknown)
// File:        CALLOS_module.cc
//
// Generated at Fri Feb  9 19:47:11 2024 by Rodrigo Alvarez Garrote using cetskelgen
// from  version .
// 
// This module is a Calibration Analyzer for Low Light Optical Signals
// CALLOS is made to read raw waveforms from raw and deconvolved stages and 
// compute the charges and average waveform for each optical channel selected.
// The resulting histograms and NTuples are saved as output.
// A SPE will be produced for each channel but ALL events. While charges will be dumped
// for each event to a TTree to prevent memory overflow (hundreds of channels).
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"


namespace callos {
  class CALLOS;
}

class AverageWaveform {
public:
  AverageWaveform(int size) : size_(size){
    array_ = new float[size_];
    wvf_count_ = 1;
  }

  ~AverageWaveform() {
    delete[] array_;
  }

  void addToAverage(float* inputArray) {
    for (int i = 0; i < size_; i++) {
      // divide by wvf_count_ to avoid overflow, maybe for large number of waveforms is better to stack a bunch and divide after
      array_[i] += inputArray[i]/wvf_count_;
    }
    wvf_count_++;
  }

private:
  float* array_;
  int size_;
  int wvf_count_;
};



class callos::CALLOS : public art::EDAnalyzer {
public:
  explicit CALLOS(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CALLOS(CALLOS const&) = delete;
  CALLOS(CALLOS&&) = delete;
  CALLOS& operator=(CALLOS const&) = delete;
  CALLOS& operator=(CALLOS&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // placeholders
  // void beginJob(art::Run& run) override ;
  // void beginRun(art::Run& run) override ;
  // void beginSubRun(art::Run& run) override ;
  
  // void endSubRun(art::SubRun& sr) override ;
  // void endRun(art::Run& run) override ;
  // void endJob(art::SubRun& sr) override ;

  // Tools shall be initialized here 
private:

  // Declare member data here.

  std::string fInputLabel;

  // Size of the Region of Interest(ROI) of each signal, must be smaller than RawWaveforms size
  const int ROI_samples;

  // Number of samples from start of the ROI to peak. 
  const int start_to_peak;

};


callos::CALLOS::CALLOS(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void callos::CALLOS::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // Load the waveforms
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fInputLabel, wfHandle);

  // Map waveforms to channels

  // select only relevant channels (this allows to sepparate diferent XAs, PMTs, etc based on PDS map stored in json format)
  // also prevents checking wvf by wvf

  // Call the ROI(peak/rise...) finder tool for selected channels

  // Compute charges and average waveforms.

  // Store then in job containers.

  // Clean variables needed for next event.

}

DEFINE_ART_MODULE(callos::CALLOS)
