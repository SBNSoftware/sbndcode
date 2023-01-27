////////////////////////////////////////////////////////////////////////
// Class:       PiZeroAnalysis
// Plugin Type: analyzer
// File:        PiZeroAnalysis_module.cc
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

#include "TTree.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

namespace sbnd {
  class PiZeroAnalysis;
}

class sbnd::PiZeroAnalysis : public art::EDAnalyzer {
public:
  explicit PiZeroAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PiZeroAnalysis(PiZeroAnalysis const&) = delete;
  PiZeroAnalysis(PiZeroAnalysis&&) = delete;
  PiZeroAnalysis& operator=(PiZeroAnalysis const&) = delete;
  PiZeroAnalysis& operator=(PiZeroAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalyseMCParticles(const std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec, const art::FindManyP<simb::MCTruth> &MCParticlesToMCTruths);

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

  std::string fMCParticleModuleLabel;
  bool fDebug;

  TTree* fPiZeroTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;

  std::vector<int>    _mct_origin;
  std::vector<double> _mc_vx;
  std::vector<double> _mc_vy;
  std::vector<double> _mc_vz;
  std::vector<double> _mc_vt;
  std::vector<bool>   _mc_av;
  std::vector<bool>   _mc_fv;
  std::vector<bool>   _mc_uboone_fv;
  std::vector<double> _mc_endx;
  std::vector<double> _mc_endy;
  std::vector<double> _mc_endz;
  std::vector<double> _mc_endt;
  std::vector<double> _mc_vpx;
  std::vector<double> _mc_vpy;
  std::vector<double> _mc_vpz;
  std::vector<double> _mc_vtheta;
  std::vector<double> _mc_vphi;
  std::vector<double> _mc_ve;
  std::vector<double> _mc_endpx;
  std::vector<double> _mc_endpy;
  std::vector<double> _mc_endpz;
  std::vector<double> _mc_ende;

  std::vector<int>    _mc_gamma1_pdg;
  std::vector<int>    _mc_gamma1_status;
  std::vector<double> _mc_gamma1_vx;
  std::vector<double> _mc_gamma1_vy;
  std::vector<double> _mc_gamma1_vz;
  std::vector<double> _mc_gamma1_vt;
  std::vector<double> _mc_gamma1_endx;
  std::vector<double> _mc_gamma1_endy;
  std::vector<double> _mc_gamma1_endz;
  std::vector<double> _mc_gamma1_endt;
  std::vector<double> _mc_gamma1_vpx;
  std::vector<double> _mc_gamma1_vpy;
  std::vector<double> _mc_gamma1_vpz;
  std::vector<double> _mc_gamma1_vtheta;
  std::vector<double> _mc_gamma1_vphi;
  std::vector<double> _mc_gamma1_ve;
  std::vector<double> _mc_gamma1_endpx;
  std::vector<double> _mc_gamma1_endpy;
  std::vector<double> _mc_gamma1_endpz;
  std::vector<double> _mc_gamma1_ende;

  std::vector<int>    _mc_gamma2_pdg;
  std::vector<int>    _mc_gamma2_status;
  std::vector<double> _mc_gamma2_vx;
  std::vector<double> _mc_gamma2_vy;
  std::vector<double> _mc_gamma2_vz;
  std::vector<double> _mc_gamma2_vt;
  std::vector<double> _mc_gamma2_endx;
  std::vector<double> _mc_gamma2_endy;
  std::vector<double> _mc_gamma2_endz;
  std::vector<double> _mc_gamma2_endt;
  std::vector<double> _mc_gamma2_vpx;
  std::vector<double> _mc_gamma2_vpy;
  std::vector<double> _mc_gamma2_vpz;
  std::vector<double> _mc_gamma2_vtheta;
  std::vector<double> _mc_gamma2_vphi;
  std::vector<double> _mc_gamma2_ve;
  std::vector<double> _mc_gamma2_endpx;
  std::vector<double> _mc_gamma2_endpy;
  std::vector<double> _mc_gamma2_endpz;
  std::vector<double> _mc_gamma2_ende;

  std::vector<double> _mc_open_angle;
};

sbnd::PiZeroAnalysis::PiZeroAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fMCParticleModuleLabel = p.get<std::string>("MCParticleModuleLabel", "largeant");
    fDebug                 = p.get<bool>("Debug", false);

    art::ServiceHandle<art::TFileService> fs;

    fPiZeroTree = fs->make<TTree>("pizeros","");
    fPiZeroTree->Branch("run", &_run);
    fPiZeroTree->Branch("subrun", &_subrun);
    fPiZeroTree->Branch("event", &_event);

    fPiZeroTree->Branch("mct_origin", "std::vector<int>", &_mct_origin);
    fPiZeroTree->Branch("mc_vx", "std::vector<double>", &_mc_vx);
    fPiZeroTree->Branch("mc_vy", "std::vector<double>", &_mc_vy);
    fPiZeroTree->Branch("mc_vz", "std::vector<double>", &_mc_vz);
    fPiZeroTree->Branch("mc_vt", "std::vector<double>", &_mc_vt);
    fPiZeroTree->Branch("mc_av", "std::vector<bool>", &_mc_av);
    fPiZeroTree->Branch("mc_fv", "std::vector<bool>", &_mc_fv);
    fPiZeroTree->Branch("mc_uboone_fv", "std::vector<bool>", &_mc_uboone_fv);
    fPiZeroTree->Branch("mc_endx", "std::vector<double>", &_mc_endx);
    fPiZeroTree->Branch("mc_endy", "std::vector<double>", &_mc_endy);
    fPiZeroTree->Branch("mc_endz", "std::vector<double>", &_mc_endz);
    fPiZeroTree->Branch("mc_endt", "std::vector<double>", &_mc_endt);
    fPiZeroTree->Branch("mc_vpx", "std::vector<double>", &_mc_vpx);
    fPiZeroTree->Branch("mc_vpy", "std::vector<double>", &_mc_vpy);
    fPiZeroTree->Branch("mc_vpz", "std::vector<double>", &_mc_vpz);
    fPiZeroTree->Branch("mc_vtheta", "std::vector<double>", &_mc_vtheta);
    fPiZeroTree->Branch("mc_vphi", "std::vector<double>", &_mc_vphi);
    fPiZeroTree->Branch("mc_ve", "std::vector<double>", &_mc_ve);
    fPiZeroTree->Branch("mc_endpx", "std::vector<double>", &_mc_endpx);
    fPiZeroTree->Branch("mc_endpy", "std::vector<double>", &_mc_endpy);
    fPiZeroTree->Branch("mc_endpz", "std::vector<double>", &_mc_endpz);
    fPiZeroTree->Branch("mc_ende", "std::vector<double>", &_mc_ende);

    fPiZeroTree->Branch("mc_gamma1_pdg", "std::vector<int>", &_mc_gamma1_pdg);
    fPiZeroTree->Branch("mc_gamma1_status", "std::vector<int>", &_mc_gamma1_status);
    fPiZeroTree->Branch("mc_gamma1_vx", "std::vector<double>", &_mc_gamma1_vx);
    fPiZeroTree->Branch("mc_gamma1_vy", "std::vector<double>", &_mc_gamma1_vy);
    fPiZeroTree->Branch("mc_gamma1_vz", "std::vector<double>", &_mc_gamma1_vz);
    fPiZeroTree->Branch("mc_gamma1_vt", "std::vector<double>", &_mc_gamma1_vt);
    fPiZeroTree->Branch("mc_gamma1_endx", "std::vector<double>", &_mc_gamma1_endx);
    fPiZeroTree->Branch("mc_gamma1_endy", "std::vector<double>", &_mc_gamma1_endy);
    fPiZeroTree->Branch("mc_gamma1_endz", "std::vector<double>", &_mc_gamma1_endz);
    fPiZeroTree->Branch("mc_gamma1_endt", "std::vector<double>", &_mc_gamma1_endt);
    fPiZeroTree->Branch("mc_gamma1_vpx", "std::vector<double>", &_mc_gamma1_vpx);
    fPiZeroTree->Branch("mc_gamma1_vpy", "std::vector<double>", &_mc_gamma1_vpy);
    fPiZeroTree->Branch("mc_gamma1_vpz", "std::vector<double>", &_mc_gamma1_vpz);
    fPiZeroTree->Branch("mc_gamma1_vtheta", "std::vector<double>", &_mc_gamma1_vtheta);
    fPiZeroTree->Branch("mc_gamma1_vphi", "std::vector<double>", &_mc_gamma1_vphi);
    fPiZeroTree->Branch("mc_gamma1_ve", "std::vector<double>", &_mc_gamma1_ve);
    fPiZeroTree->Branch("mc_gamma1_endpx", "std::vector<double>", &_mc_gamma1_endpx);
    fPiZeroTree->Branch("mc_gamma1_endpy", "std::vector<double>", &_mc_gamma1_endpy);
    fPiZeroTree->Branch("mc_gamma1_endpz", "std::vector<double>", &_mc_gamma1_endpz);
    fPiZeroTree->Branch("mc_gamma1_ende", "std::vector<double>", &_mc_gamma1_ende);

    fPiZeroTree->Branch("mc_gamma2_pdg", "std::vector<int>", &_mc_gamma2_pdg);
    fPiZeroTree->Branch("mc_gamma2_status", "std::vector<int>", &_mc_gamma2_status);
    fPiZeroTree->Branch("mc_gamma2_vx", "std::vector<double>", &_mc_gamma2_vx);
    fPiZeroTree->Branch("mc_gamma2_vy", "std::vector<double>", &_mc_gamma2_vy);
    fPiZeroTree->Branch("mc_gamma2_vz", "std::vector<double>", &_mc_gamma2_vz);
    fPiZeroTree->Branch("mc_gamma2_vt", "std::vector<double>", &_mc_gamma2_vt);
    fPiZeroTree->Branch("mc_gamma2_endx", "std::vector<double>", &_mc_gamma2_endx);
    fPiZeroTree->Branch("mc_gamma2_endy", "std::vector<double>", &_mc_gamma2_endy);
    fPiZeroTree->Branch("mc_gamma2_endz", "std::vector<double>", &_mc_gamma2_endz);
    fPiZeroTree->Branch("mc_gamma2_endt", "std::vector<double>", &_mc_gamma2_endt);
    fPiZeroTree->Branch("mc_gamma2_vpx", "std::vector<double>", &_mc_gamma2_vpx);
    fPiZeroTree->Branch("mc_gamma2_vpy", "std::vector<double>", &_mc_gamma2_vpy);
    fPiZeroTree->Branch("mc_gamma2_vpz", "std::vector<double>", &_mc_gamma2_vpz);
    fPiZeroTree->Branch("mc_gamma2_vtheta", "std::vector<double>", &_mc_gamma2_vtheta);
    fPiZeroTree->Branch("mc_gamma2_vphi", "std::vector<double>", &_mc_gamma2_vphi);
    fPiZeroTree->Branch("mc_gamma2_ve", "std::vector<double>", &_mc_gamma2_ve);
    fPiZeroTree->Branch("mc_gamma2_endpx", "std::vector<double>", &_mc_gamma2_endpx);
    fPiZeroTree->Branch("mc_gamma2_endpy", "std::vector<double>", &_mc_gamma2_endpy);
    fPiZeroTree->Branch("mc_gamma2_endpz", "std::vector<double>", &_mc_gamma2_endpz);
    fPiZeroTree->Branch("mc_gamma2_ende", "std::vector<double>", &_mc_gamma2_ende);

    fPiZeroTree->Branch("mc_open_angle", "std::vector<double>", &_mc_open_angle);
  }

void sbnd::PiZeroAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  if(fDebug) std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get MCParticles
  art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
  e.getByLabel(fMCParticleModuleLabel, MCParticleHandle);
  if(!MCParticleHandle.isValid()){
    std::cout << "MCParticle product " << fMCParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCParticle>> MCParticleVec;
  art::fill_ptr_vector(MCParticleVec, MCParticleHandle);

  // Get MCTruth to MCParticles Assns
  art::FindManyP<simb::MCTruth> MCParticlesToMCTruths(MCParticleHandle, e, fMCParticleModuleLabel);

  // Fill MCParticle variables
  AnalyseMCParticles(MCParticleVec, MCParticlesToMCTruths);

  // Fill the Tree
  fPiZeroTree->Fill();
}

void sbnd::PiZeroAnalysis::AnalyseMCParticles(const std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec, const art::FindManyP<simb::MCTruth> &MCParticlesToMCTruths)
{
  const unsigned nMCParticles = MCParticleVec.size();
  unsigned nPiZeros           = 0;

  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      const auto mcp = MCParticleVec[i];

      if(mcp->PdgCode() != 111 || mcp->StatusCode() != 1)
	continue;

      ++nPiZeros;
    }

  _mct_origin.resize(nPiZeros);
  _mc_vx.resize(nPiZeros);
  _mc_vy.resize(nPiZeros);
  _mc_vz.resize(nPiZeros);
  _mc_vt.resize(nPiZeros);
  _mc_av.resize(nPiZeros);
  _mc_fv.resize(nPiZeros);
  _mc_uboone_fv.resize(nPiZeros);
  _mc_endx.resize(nPiZeros);
  _mc_endy.resize(nPiZeros);
  _mc_endz.resize(nPiZeros);
  _mc_endt.resize(nPiZeros);
  _mc_vpx.resize(nPiZeros);
  _mc_vpy.resize(nPiZeros);
  _mc_vpz.resize(nPiZeros);
  _mc_vtheta.resize(nPiZeros);
  _mc_vphi.resize(nPiZeros);
  _mc_ve.resize(nPiZeros);
  _mc_endpx.resize(nPiZeros);
  _mc_endpy.resize(nPiZeros);
  _mc_endpz.resize(nPiZeros);
  _mc_ende.resize(nPiZeros);

  _mc_gamma1_pdg.resize(nPiZeros);
  _mc_gamma1_status.resize(nPiZeros);
  _mc_gamma1_vx.resize(nPiZeros);
  _mc_gamma1_vy.resize(nPiZeros);
  _mc_gamma1_vz.resize(nPiZeros);
  _mc_gamma1_vt.resize(nPiZeros);
  _mc_gamma1_endx.resize(nPiZeros);
  _mc_gamma1_endy.resize(nPiZeros);
  _mc_gamma1_endz.resize(nPiZeros);
  _mc_gamma1_endt.resize(nPiZeros);
  _mc_gamma1_vpx.resize(nPiZeros);
  _mc_gamma1_vpy.resize(nPiZeros);
  _mc_gamma1_vpz.resize(nPiZeros);
  _mc_gamma1_vtheta.resize(nPiZeros);
  _mc_gamma1_vphi.resize(nPiZeros);
  _mc_gamma1_ve.resize(nPiZeros);
  _mc_gamma1_endpx.resize(nPiZeros);
  _mc_gamma1_endpy.resize(nPiZeros);
  _mc_gamma1_endpz.resize(nPiZeros);
  _mc_gamma1_ende.resize(nPiZeros);

  _mc_gamma2_pdg.resize(nPiZeros);
  _mc_gamma2_status.resize(nPiZeros);
  _mc_gamma2_vx.resize(nPiZeros);
  _mc_gamma2_vy.resize(nPiZeros);
  _mc_gamma2_vz.resize(nPiZeros);
  _mc_gamma2_vt.resize(nPiZeros);
  _mc_gamma2_endx.resize(nPiZeros);
  _mc_gamma2_endy.resize(nPiZeros);
  _mc_gamma2_endz.resize(nPiZeros);
  _mc_gamma2_endt.resize(nPiZeros);
  _mc_gamma2_vpx.resize(nPiZeros);
  _mc_gamma2_vpy.resize(nPiZeros);
  _mc_gamma2_vpz.resize(nPiZeros);
  _mc_gamma2_vtheta.resize(nPiZeros);
  _mc_gamma2_vphi.resize(nPiZeros);
  _mc_gamma2_ve.resize(nPiZeros);
  _mc_gamma2_endpx.resize(nPiZeros);
  _mc_gamma2_endpy.resize(nPiZeros);
  _mc_gamma2_endpz.resize(nPiZeros);
  _mc_gamma2_ende.resize(nPiZeros);

  _mc_open_angle.resize(nPiZeros);

  unsigned pizero = 0;
  
  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      const auto mcp = MCParticleVec[i];

      if(mcp->PdgCode() != 111 || mcp->StatusCode() != 1)
	continue;

      const auto mcts = MCParticlesToMCTruths.at(mcp.key());
      if(mcts.size() != 1) std::cout << "There are multiple truths" << std::endl;
      const auto mct  = mcts[0];

      _mct_origin[pizero]   = mct->Origin();
      _mc_vx[pizero]        = mcp->Vx();
      _mc_vy[pizero]        = mcp->Vy();
      _mc_vz[pizero]        = mcp->Vz();
      _mc_vt[pizero]        = mcp->T();
      _mc_av[pizero]        = VolumeCheck(mcp->Position().Vect());
      _mc_fv[pizero]        = VolumeCheck(mcp->Position().Vect(), 20., 10., 10., 50.);
      _mc_uboone_fv[pizero] = VolumeCheck(mcp->Position().Vect(), 5., 5., 5., 5.);
      _mc_endx[pizero]      = mcp->EndX();
      _mc_endy[pizero]      = mcp->EndY();
      _mc_endz[pizero]      = mcp->EndZ();
      _mc_endt[pizero]      = mcp->EndT();
      _mc_vpx[pizero]       = mcp->Px();
      _mc_vpy[pizero]       = mcp->Py();
      _mc_vpz[pizero]       = mcp->Pz();
      _mc_vtheta[pizero]    = TMath::RadToDeg() * mcp->Momentum().Theta();
      _mc_vphi[pizero]      = TMath::RadToDeg() * mcp->Momentum().Phi();
      _mc_ve[pizero]        = mcp->E();
      _mc_endpx[pizero]     = mcp->EndPx();
      _mc_endpy[pizero]     = mcp->EndPy();
      _mc_endpz[pizero]     = mcp->EndPz();
      _mc_ende[pizero]      = mcp->EndE();

      if(mcp->NumberDaughters() == 2)
	{
	  const simb::MCParticle* gamma1 = particleInv->TrackIdToParticle_P(mcp->Daughter(0));
	  const simb::MCParticle* gamma2 = particleInv->TrackIdToParticle_P(mcp->Daughter(1));

	  const simb::MCParticle* gammalead    = gamma2->E() > gamma1->E() ? gamma2 : gamma1;
	  const simb::MCParticle* gammasublead = gamma2->E() > gamma1->E() ? gamma1 : gamma2;
	  
	  _mc_gamma1_pdg[pizero]    = gammalead->PdgCode();
	  _mc_gamma1_status[pizero] = gammalead->StatusCode();
	  _mc_gamma1_vx[pizero]     = gammalead->Vx();
	  _mc_gamma1_vy[pizero]     = gammalead->Vy();
	  _mc_gamma1_vz[pizero]     = gammalead->Vz();
	  _mc_gamma1_vt[pizero]     = gammalead->T();
	  _mc_gamma1_endx[pizero]   = gammalead->EndX();
	  _mc_gamma1_endy[pizero]   = gammalead->EndY();
	  _mc_gamma1_endz[pizero]   = gammalead->EndZ();
	  _mc_gamma1_endt[pizero]   = gammalead->EndT();
	  _mc_gamma1_vpx[pizero]    = gammalead->Px();
	  _mc_gamma1_vpy[pizero]    = gammalead->Py();
	  _mc_gamma1_vpz[pizero]    = gammalead->Pz();
	  _mc_gamma1_vtheta[pizero] = TMath::RadToDeg() * gammalead->Momentum().Theta();
	  _mc_gamma1_vphi[pizero]   = TMath::RadToDeg() * gammalead->Momentum().Phi();
	  _mc_gamma1_ve[pizero]     = gammalead->E();
	  _mc_gamma1_endpx[pizero]  = gammalead->EndPx();
	  _mc_gamma1_endpy[pizero]  = gammalead->EndPy();
	  _mc_gamma1_endpz[pizero]  = gammalead->EndPz();
	  _mc_gamma1_ende[pizero]   = gammalead->EndE();
	  
	  _mc_gamma2_pdg[pizero]    = gammasublead->PdgCode();
	  _mc_gamma2_status[pizero] = gammasublead->StatusCode();
	  _mc_gamma2_vx[pizero]     = gammasublead->Vx();
	  _mc_gamma2_vy[pizero]     = gammasublead->Vy();
	  _mc_gamma2_vz[pizero]     = gammasublead->Vz();
	  _mc_gamma2_vt[pizero]     = gammasublead->T();
	  _mc_gamma2_endx[pizero]   = gammasublead->EndX();
	  _mc_gamma2_endy[pizero]   = gammasublead->EndY();
	  _mc_gamma2_endz[pizero]   = gammasublead->EndZ();
	  _mc_gamma2_endt[pizero]   = gammasublead->EndT();
	  _mc_gamma2_vpx[pizero]    = gammasublead->Px();
	  _mc_gamma2_vpy[pizero]    = gammasublead->Py();
	  _mc_gamma2_vpz[pizero]    = gammasublead->Pz();
	  _mc_gamma2_vtheta[pizero] = TMath::RadToDeg() * gammasublead->Momentum().Theta();
	  _mc_gamma2_vphi[pizero]   = TMath::RadToDeg() * gammasublead->Momentum().Phi();
	  _mc_gamma2_ve[pizero]     = gammasublead->E();
	  _mc_gamma2_endpx[pizero]  = gammasublead->EndPx();
	  _mc_gamma2_endpy[pizero]  = gammasublead->EndPy();
	  _mc_gamma2_endpz[pizero]  = gammasublead->EndPz();
	  _mc_gamma2_ende[pizero]   = gammasublead->EndE();

	  const TVector3 gamma1_mom = gammalead->Momentum().Vect();
	  const TVector3 gamma2_mom = gammasublead->Momentum().Vect();
	  _mc_open_angle[pizero]    = TMath::RadToDeg() * gamma1_mom.Angle(gamma2_mom);
	}
      else
	{
	  _mc_gamma1_pdg[pizero]    = -999999;
	  _mc_gamma1_status[pizero] = -999999;
	  _mc_gamma1_vx[pizero]     = -999999.;
	  _mc_gamma1_vy[pizero]     = -999999.;
	  _mc_gamma1_vz[pizero]     = -999999.;
	  _mc_gamma1_vt[pizero]     = -999999.;
	  _mc_gamma1_endx[pizero]   = -999999.;
	  _mc_gamma1_endy[pizero]   = -999999.;
	  _mc_gamma1_endz[pizero]   = -999999.;
	  _mc_gamma1_endt[pizero]   = -999999.;
	  _mc_gamma1_vpx[pizero]    = -999999.;
	  _mc_gamma1_vpy[pizero]    = -999999.;
	  _mc_gamma1_vpz[pizero]    = -999999.;
	  _mc_gamma1_ve[pizero]     = -999999.;
	  _mc_gamma1_endpx[pizero]  = -999999.;
	  _mc_gamma1_endpy[pizero]  = -999999.;
	  _mc_gamma1_endpz[pizero]  = -999999.;
	  _mc_gamma1_ende[pizero]   = -999999.;

	  _mc_gamma2_pdg[pizero]    = -999999;
	  _mc_gamma2_status[pizero] = -999999;
	  _mc_gamma2_vx[pizero]     = -999999.;
	  _mc_gamma2_vy[pizero]     = -999999.;
	  _mc_gamma2_vz[pizero]     = -999999.;
	  _mc_gamma2_vt[pizero]     = -999999.;
	  _mc_gamma2_endx[pizero]   = -999999.;
	  _mc_gamma2_endy[pizero]   = -999999.;
	  _mc_gamma2_endz[pizero]   = -999999.;
	  _mc_gamma2_endt[pizero]   = -999999.;
	  _mc_gamma2_vpx[pizero]    = -999999.;
	  _mc_gamma2_vpy[pizero]    = -999999.;
	  _mc_gamma2_vpz[pizero]    = -999999.;
	  _mc_gamma2_ve[pizero]     = -999999.;
	  _mc_gamma2_endpx[pizero]  = -999999.;
	  _mc_gamma2_endpy[pizero]  = -999999.;
	  _mc_gamma2_endpz[pizero]  = -999999.;
	  _mc_gamma2_ende[pizero]   = -999999.;
	}
      ++pizero;
    }
}

bool sbnd::PiZeroAnalysis::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

DEFINE_ART_MODULE(sbnd::PiZeroAnalysis)
