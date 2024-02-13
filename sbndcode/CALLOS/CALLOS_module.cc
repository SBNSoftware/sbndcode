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
// A SPE and a vector of charges is produced for each channel and for ALL events processed.
// ROI tools and PEAK finders are defined as Tools with templates to allow for
// different implementations to be used. 
// 
// The aim of this module is to be as fast as possible, hence the conscious decision of not using 
// if checks when possible. Specially for analyze() function(run on each event).
// 
// Note: the module assumes raw::waveform channels go from 0 to NChannels-1. 
// If needed, change from std::vector to std::map.
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
#include "AverageWaveform.h"

#include <memory>
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"
#include <iostream>
#include <string>
#include <vector>

#include "Tools/ROIFinderAlg.hh"

namespace callos {
  class CALLOS;
}

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
  void endJob() override ;

  // Tools shall be initialized here 
private:

  // Declare member data here.

  // fhicl parameters
  std::string fInputLabel;
  std::string fPDType;
  // std::vector<std::string> fElectronics; //Not needed for now

  //PDS map
  opdet::sbndPDMapAlg pdsmap;

  // Total number of channels
  int fNPDSChannels= pdsmap.size();
  
  // Filter for the loop and mapping
  std::vector<bool> fPDSchannelStatus;
  std::vector<int>  fPDSchannelMap;

  // Selected channels
  std::vector<int> fSelectedChannels;
  // Number of selected channels
  int fNSelectedChannels;

  // Size of the Region of Interest(ROI) of each signal, must be smaller than RawWaveforms size. rodrigoa: maybe should be a tool parameter
  int fROISamples;

  // Number of samples from start of the ROI to peak. rodrigoa: maybe should be a tool parameter
  int fStartToPeak;

  // Average Waveforms container.
  std::vector<AverageWaveform> AverageWaveform_SelectedChannels;
  //Charge container
  std::vector<std::vector<float>> Charge_SelectedChannels;

  // Tool pointers
  std::unique_ptr<callos::ROIFINDERALG> fROIFinderAlgPtr;
};


callos::CALLOS::CALLOS(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fInputLabel = p.get<std::string>("InputLabel","opdaq");
  fPDType = p.get<std::string>("PDType","xarapuca_vuv");
  fROISamples = p.get<int>("ROI_samples", 1000);
  fStartToPeak = p.get<int>("StartToPeak",200);

  // get map info
  std::vector<int> fSelectedChannels = pdsmap.getChannelsOfType(fPDType);
  fNSelectedChannels = fSelectedChannels.size();
  std::cout<<"CALLOS: Selected "<<fNSelectedChannels<<" channels of selected type"<<std::endl;

  // Prepare the filter and map
  fPDSchannelStatus.resize(fNPDSChannels);
  fPDSchannelMap   .resize(fNPDSChannels);
  for (int i=0; i<fNSelectedChannels; i++)
  { 
    fPDSchannelStatus[fSelectedChannels[i]] = true; 
    fPDSchannelMap   [fSelectedChannels[i]] = i;
  }

  // Initialize the average waveforms container
  AverageWaveform_SelectedChannels.reserve(fNSelectedChannels);
  for (int i=0; i<fNSelectedChannels; i++){ AverageWaveform_SelectedChannels.push_back(AverageWaveform(fROISamples));}
  Charge_SelectedChannels.reserve(fNSelectedChannels);
  // Initialize the ROIFinderAlg tool
  fROIFinderAlgPtr = art::make_tool<callos::ROIFINDERALG>(p.get<fhicl::ParameterSet>("ROIFinderAlg"));
}

void callos::CALLOS::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  // Load the waveforms
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fInputLabel, wfHandle);

  if (!wfHandle.isValid()) 
  {
   mf::LogError("SBNDCALLOS")<<"Input waveforms with input label "<<fInputLabel<<" couldn't be loaded..."<<std::endl;
   throw cet::exception("SBNDCALLOS") << "Input waveforms with input label " << fInputLabel << " not found\n";
  }
  
  // Get the Raw waveforms
  auto NWvf = wfHandle->size();
  std::cout<<"CALLOS: Event "<<e.id().event()<<" has      "<<NWvf<<" waveforms"<<std::endl;
  int iWvf=0;
    for(auto const& wf : *wfHandle)
    {
      // Get raw wvf info
      int wfChannel = wf.ChannelNumber();


      // select only relevant channels (this allows to sepparate diferent XAs, PMTs, etc based on PDS map)
      if (fPDSchannelStatus[wfChannel]) // faster than checking if the channel is in the list of selected channels
      {

        // move Raw wvf to float
        size_t wfsize=wf.Waveform().size();
        std::vector<float> wave;
        wave.reserve(wfsize);
        wave.assign(wf.Waveform().begin(), wf.Waveform().end());
        //Prepare ROI container and initialize to 0s
        std::vector<float> ROI;
        ROI.resize(fROISamples,0);

        // Call the tool for selected channels
        fROIFinderAlgPtr->ProcessWaveform(wave, ROI, Charge_SelectedChannels[fPDSchannelMap[wfChannel]]);
        
        // Update the Average wvf of the channel.
        AverageWaveform_SelectedChannels[fPDSchannelMap[wfChannel]].addToAverage(ROI);

        iWvf++;
      }
    }
  std::cout<<"CALLOS: Event "<<e.id().event()<<" analized "<<iWvf<<" waveforms"<<std::endl;
  // Store info in containers.

  // Clean variables needed for next event.

}

void callos::CALLOS::endJob() 
{
  // Pending real implementation
  // Save the histograms and NTuples to the output file.
  // Close the output file

  std::cout<<"CALLOS: End of Job"<<std::endl;
  // free the memory?
}



DEFINE_ART_MODULE(callos::CALLOS)
