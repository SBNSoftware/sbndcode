////////////////////////////////////////////////////////////////////////
// Class:       SimpleAnalyzer
// Plugin Type: analyzer
// File:        SimpleAnalyzer_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
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
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

namespace sbnd {
  class SimpleAnalyzer;
}

class sbnd::SimpleAnalyzer : public art::EDAnalyzer {
public:
  explicit SimpleAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimpleAnalyzer(SimpleAnalyzer const&) = delete;
  SimpleAnalyzer(SimpleAnalyzer&&) = delete;
  SimpleAnalyzer& operator=(SimpleAnalyzer const&) = delete;
  SimpleAnalyzer& operator=(SimpleAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec);

private:

  std::string fMCTruthModuleLabel, fMCParticleModuleLabel;

  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;
  int _interaction;

  std::vector<int16_t>          _mc_trackid;
  std::vector<int16_t>          _mc_pdg;
  std::vector<int16_t>          _mc_status;
  std::vector<std::string>      _mc_process;
  std::vector<uint16_t>         _mc_nchildren;
  std::vector<std::vector<int>> _mc_children;
  std::vector<double>           _mc_vx;
  std::vector<double>           _mc_vy;
  std::vector<double>           _mc_vz;
  std::vector<double>           _mc_vt;
  std::vector<double>           _mc_endx;
  std::vector<double>           _mc_endy;
  std::vector<double>           _mc_endz;
  std::vector<double>           _mc_endt;
  std::vector<double>           _mc_vpx;
  std::vector<double>           _mc_vpy;
  std::vector<double>           _mc_vpz;
  std::vector<double>           _mc_ve;
  std::vector<double>           _mc_endpx;
  std::vector<double>           _mc_endpy;
  std::vector<double>           _mc_endpz;
  std::vector<double>           _mc_ende;
};

sbnd::SimpleAnalyzer::SimpleAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fMCTruthModuleLabel    = p.get<std::string>("MCTruthModuleLabel");
    fMCParticleModuleLabel = p.get<std::string>("MCParticleModuleLabel");

    art::ServiceHandle<art::TFileService> fs;

    fTree = fs->make<TTree>("tree","");
    fTree->Branch("run", &_run);
    fTree->Branch("subrun", &_subrun);
    fTree->Branch("event", &_event);
    fTree->Branch("interaction", &_interaction);

    fTree->Branch("mc_trackid", "std::vector<int16_t>", &_mc_trackid);
    fTree->Branch("mc_pdg", "std::vector<int16_t>", &_mc_pdg);
    fTree->Branch("mc_status", "std::vector<int16_t>", &_mc_status);
    fTree->Branch("mc_process", "std::vector<std::string>", &_mc_process);
    fTree->Branch("mc_nchildren", "std::vector<uint16_t>", &_mc_nchildren);
    fTree->Branch("mc_children", "std::vector<std::vector<int>>", &_mc_children);
    fTree->Branch("mc_vx", "std::vector<double>", &_mc_vx);
    fTree->Branch("mc_vy", "std::vector<double>", &_mc_vy);
    fTree->Branch("mc_vz", "std::vector<double>", &_mc_vz);
    fTree->Branch("mc_vt", "std::vector<double>", &_mc_vt);
    fTree->Branch("mc_endx", "std::vector<double>", &_mc_endx);
    fTree->Branch("mc_endy", "std::vector<double>", &_mc_endy);
    fTree->Branch("mc_endz", "std::vector<double>", &_mc_endz);
    fTree->Branch("mc_endt", "std::vector<double>", &_mc_endt);
    fTree->Branch("mc_vpx", "std::vector<double>", &_mc_vpx);
    fTree->Branch("mc_vpy", "std::vector<double>", &_mc_vpy);
    fTree->Branch("mc_vpz", "std::vector<double>", &_mc_vpz);
    fTree->Branch("mc_ve", "std::vector<double>", &_mc_ve);
    fTree->Branch("mc_endpx", "std::vector<double>", &_mc_endpx);
    fTree->Branch("mc_endpy", "std::vector<double>", &_mc_endpy);
    fTree->Branch("mc_endpz", "std::vector<double>", &_mc_endpz);
    fTree->Branch("mc_ende", "std::vector<double>", &_mc_ende);
  }

void sbnd::SimpleAnalyzer::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  // Get MCTruths
  art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
  e.getByLabel(fMCTruthModuleLabel, MCTruthHandle);
  if(!MCTruthHandle.isValid()){
    std::cout << "MCTruth product " << fMCTruthModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
  art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

  _interaction = 0;

  for(auto const& mct : MCTruthVec)
    {
      const simb::MCNeutrino mcn = mct->GetNeutrino();

      if(mcn.CCNC() == 0)
	continue;

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles( { mct }, e, fMCParticleModuleLabel);
      std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(0);

      bool signal = false;

      for(auto const& mcp : MCParticleVec)
	{
	  if(mcp->Process() == "primary" && mcp->StatusCode() == 1 && mcp->PdgCode() == 111)
	    signal = true;
	}

      if(signal)
	{
	  // Fill MCParticle variables
	  AnalyseMCParticles(MCParticleVec);

	  // Fill the Tree
	  fTree->Fill();

	  ++_interaction;
	}
    }
}

void sbnd::SimpleAnalyzer::AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec)
{
  const unsigned nMCParticles = MCParticleVec.size();

  _mc_trackid.resize(nMCParticles);
  _mc_pdg.resize(nMCParticles);
  _mc_status.resize(nMCParticles);
  _mc_process.resize(nMCParticles);
  _mc_nchildren.resize(nMCParticles);
  _mc_children.resize(nMCParticles);
  _mc_vx.resize(nMCParticles);
  _mc_vy.resize(nMCParticles);
  _mc_vz.resize(nMCParticles);
  _mc_vt.resize(nMCParticles);
  _mc_endx.resize(nMCParticles);
  _mc_endy.resize(nMCParticles);
  _mc_endz.resize(nMCParticles);
  _mc_endt.resize(nMCParticles);
  _mc_vpx.resize(nMCParticles);
  _mc_vpy.resize(nMCParticles);
  _mc_vpz.resize(nMCParticles);
  _mc_ve.resize(nMCParticles);
  _mc_endpx.resize(nMCParticles);
  _mc_endpy.resize(nMCParticles);
  _mc_endpz.resize(nMCParticles);
  _mc_ende.resize(nMCParticles);
  
  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      const auto mcp = MCParticleVec[i];

      _mc_trackid[i]    = mcp->TrackId();
      _mc_pdg[i]        = mcp->PdgCode();
      _mc_status[i]     = mcp->StatusCode();
      _mc_process[i]    = mcp->Process();
      _mc_nchildren[i] = mcp->NumberDaughters();
      _mc_vx[i]         = mcp->Vx();
      _mc_vy[i]         = mcp->Vy();
      _mc_vz[i]         = mcp->Vz();
      _mc_vt[i]         = mcp->T();
      _mc_endx[i]       = mcp->EndX();
      _mc_endy[i]       = mcp->EndY();
      _mc_endz[i]       = mcp->EndZ();
      _mc_endt[i]       = mcp->EndT();
      _mc_vpx[i]        = mcp->Px();
      _mc_vpy[i]        = mcp->Py();
      _mc_vpz[i]        = mcp->Pz();
      _mc_ve[i]         = mcp->E();
      _mc_endpx[i]      = mcp->EndPx();
      _mc_endpy[i]      = mcp->EndPy();
      _mc_endpz[i]      = mcp->EndPz();
      _mc_ende[i]       = mcp->EndE();

      _mc_children[i].resize(_mc_nchildren[i]);

      for(unsigned ii = 0; ii < _mc_nchildren[i]; ++ii)
        _mc_children[i][ii] = mcp->Daughter(ii);
    }
}

DEFINE_ART_MODULE(sbnd::SimpleAnalyzer)
