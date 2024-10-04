#include "sbndcode/SERCalibration/SERCalibration_module.hh"

opdet::SERCalibration::SERCalibration(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  //
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fInputLabel = p.get< std::string >("InputLabel");
  fPDTypes = p.get< std::vector<std::string> >("PDTypes");
  fElectronics = p.get< std::vector<std::string> >("Electronics");
  fSERPulseFinderPtr = art::make_tool<opdet::SERPulseFinderBase>( p.get< fhicl::ParameterSet >("SERPulseFinder") );
}


// -------- Begi job function --------
void opdet::SERCalibration::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("OpAnaTree", "Optical Analyzer Tree");
  fTree->Branch("eventID", &_eventID, "eventID/i");
  fTree->Branch("runID", &_runID, "runID/i");
  fTree->Branch("subrunID", &_subrunID, "subrunID/i");
  
  // Initialize the SER vector
  for(size_t i=0; i<320; i++)
  {
    TH1D SER(TString::Format("SER for channel %zu", i), "" , fSEREnd-fSERStart, fSERStart, fSEREnd);
    calibratedSER_v.push_back(SER);
  }
}


// -------- Main function --------
void opdet::SERCalibration::analyze(art::Event const& e)
{
  
  // --- Event General Info
  _eventID = e.id().event();
  _runID = e.id().run();
  _subrunID = e.id().subRun();

  if(fVerbosity>0)
    std::cout << " -- Running SERCalibrationModule -- \n Run=" << _runID << " Subrun=" << _subrunID << " Event=" << _eventID << std::endl;

  //Load the waveforms
  art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
  e.getByLabel(fInputLabel, wfHandle);
  if (!wfHandle.isValid()) {
   mf::LogError("SBNDOpDeconvolution")<<"Input waveforms with input label "<<fInputLabel<<" couldn't be loaded..."<<std::endl;
   throw cet::exception("SBNDOpDeconvolution") << "Input waveforms with input label " << fInputLabel << " not found\n";
  }

  std::vector< raw::OpDetWaveform > RawWfVector;
  RawWfVector.reserve(wfHandle->size());
  for(auto const& wf : *wfHandle){
    if((std::find(fPDTypes.begin(), fPDTypes.end(), pdsmap.pdType(wf.ChannelNumber()) ) != fPDTypes.end() ) &&
       (std::find(fElectronics.begin(), fElectronics.end(), pdsmap.electronicsType(wf.ChannelNumber()) ) != fElectronics.end()) ){
      RawWfVector.push_back(wf);
    }
  }
  fSERPulseFinderPtr->RunSERCalibration(RawWfVector, calibratedSER_v);

}


void opdet::SERCalibration::endJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  for(size_t i=0; i<320; i++)
  {
    TH1D *wvfHist = tfs->make< TH1D >(TString::Format("SER for channel %zu", i), "Time [TTicks]", calibratedSER_v[i].GetNbinsX(), fSERStart, fSEREnd);
    for(size_t j=0; j<size_t(fSEREnd-fSERStart); j++)
    {
      std::cout << " Setting bin content of channel " << i << " in time tick " << j << " to " << calibratedSER_v[i].GetBinContent(j+1) << std::endl;
      wvfHist->SetBinContent(j+1, calibratedSER_v[i].GetBinContent(j+1));
    }
  }
}


DEFINE_ART_MODULE(opdet::SERCalibration)



