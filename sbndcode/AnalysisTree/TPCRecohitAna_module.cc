////////////////////////////////////////////////////////////////////////
// Class:       TPCRecohitAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        TPCRecohitAna_module.cc
//
// Generated at Wed Jul 17 14:05:22 2024 by Thomas Junk using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
#ifndef TPCRECOHITANA_H
#define TPCRECOHITAN_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "art_root_io/TFileService.h"

#include <TTree.h>

using namespace std;

class TPCRecohitAna;

class TPCRecohitAna : public art::EDAnalyzer {
public:
  explicit TPCRecohitAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCRecohitAna(TPCRecohitAna const&) = delete;
  TPCRecohitAna(TPCRecohitAna&&) = delete;
  TPCRecohitAna& operator=(TPCRecohitAna const&) = delete;
  TPCRecohitAna& operator=(TPCRecohitAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  //std::vector<art::Ptr<recob::Hit>> _recohit_handle;

  std::string         fHitProducer;
  
private:

  TTree *fEventTree;
  int _run;
  int _subrun;
  int _event;
  std::vector<int> _chanid;
  std::vector<double> _peaktime;
  std::vector<double> _rms;
  std::vector<double> _peakamp;
  std::vector<double> _integ;
};


TPCRecohitAna::TPCRecohitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fHitProducer = p.get<std::string> ("recohit_producer", "gaushit");
}

void TPCRecohitAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

  fEventTree = tfs->make<TTree>("EventTree", "event by event info");
  fEventTree->Branch("run", &_run);
  fEventTree->Branch("subrun", &_subrun);
  fEventTree->Branch("event", &_event);
  fEventTree->Branch("chanid", &_chanid);
  fEventTree->Branch("peaktime", &_peaktime);
  fEventTree->Branch("rms", &_rms);
  fEventTree->Branch("peakamp", &_peakamp);
  fEventTree->Branch("integ", &_integ);
}

void TPCRecohitAna::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event = e.id().event();
  _chanid.clear();
  _peaktime.clear();
  _rms.clear();
  _peakamp.clear();
  _integ.clear();

  // Get time stamp
  art::Handle< std::vector<recob::Hit>> _recohit_handle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (e.getByLabel(fHitProducer,_recohit_handle))
    art::fill_ptr_vector(hitlist, _recohit_handle);

  // analyze timestamp distributions
  // collect timestamps
  for(size_t i=0; i<hitlist.size(); i++){
    int this_chanid = hitlist[i]->Channel();
    double this_peadktime = hitlist[i] -> PeakTime();
    double this_rms = hitlist[i] -> RMS();
    double this_peakamp = hitlist[i] -> PeakAmplitude();
    double this_integ = hitlist[i] -> Integral();

    _chanid.push_back(this_chanid);
    _peaktime.push_back(this_peadktime);
    _rms.push_back(this_rms);
    _peakamp.push_back(this_peakamp);
    _integ.push_back(this_integ);
  }

  fEventTree->Fill();
}

void TPCRecohitAna::endJob()
{
  mf::LogVerbatim("TPCRecohitAna") << "TPCRecohitAna finished job";
}

DEFINE_ART_MODULE(TPCRecohitAna)

#endif
