////////////////////////////////////////////////////////////////////////
// Class:       OnlinePMTAna
// Plugin Type: analyzer (art v3_03_01)
// File:        OnlinePMTAna.cxx
//
// Generated at Sun Dec 15 20:17:27 2019 by Gray Putnam using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "OnlinePMTAna.h"

#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/RawData/OpDetWaveform.h"

sbnom::OnlinePMTAna::OnlinePMTAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  // access fhicl config -- get the name of the "tag" of the input
  // PMT waveform data
  fWaveformTag = p.get<std::string>("waveform_tag", "opdaq"); 

  // Get a handle to the "File Service"
  // This is a way of accessing a regular ROOT file
  // in the art framework
  art::ServiceHandle<art::TFileService> tfs;

  // make a TTree to save stuff -- you can also make TH1's, etc...
  // art will take care of writing the TTree on exit
  fOutput = tfs->make<TTree>("sbnom_pmt", "SBN Online PMT Data");

  // make a branch
  fOutput->Branch("n_waveforms", &fNWaveforms);
}

void sbnom::OnlinePMTAna::analyze(art::Event const& e)
{
  // Get the handle of OpDetWaveforms
  const art::ValidHandle<std::vector<raw::OpDetWaveform>> waveform_handle = e.getValidHandle<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  // make a list
  const std::vector<raw::OpDetWaveform> waveforms = *waveform_handle;

  // print out stuff
  for (const raw::OpDetWaveform &wvf: waveforms) {
    std::cout << "Waveform at channel: " << wvf.ChannelNumber() << " at time: " << wvf.TimeStamp() << " with size: " << wvf.size() << std::endl;;
  }

  // save stuff into a (regular) root file
  //
  // At some point any interesting analysis stuff will have to go to a database
  // but for now TTree's / histograms are fine
  fNWaveforms = waveforms.size();
  fOutput->Fill();
}

