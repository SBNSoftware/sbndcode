////////////////////////////////////////////////////////////////////////
// Class:       NCPiZeroAnalysis
// Plugin Type: analyzer
// File:        NCPiZeroAnalysis_module.cc
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
#include "canvas/Persistency/Common/FindOneP.h"

#include "TTree.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "NCPiZeroStructs.h"

#include <numeric>

constexpr int def_int       = -std::numeric_limits<int>::max();
constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NCPiZeroAnalysis;
}

class sbnd::NCPiZeroAnalysis : public art::EDAnalyzer {
public:
  explicit NCPiZeroAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NCPiZeroAnalysis(NCPiZeroAnalysis const&) = delete;
  NCPiZeroAnalysis(NCPiZeroAnalysis&&) = delete;
  NCPiZeroAnalysis& operator=(NCPiZeroAnalysis const&) = delete;
  NCPiZeroAnalysis& operator=(NCPiZeroAnalysis&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;

  void ResetVars();
  void ResizeNeutrinoVectors(const int size);

  void ClearMaps(const art::Event &e);

  void SetupMaps(const art::Event &e);

  void AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  double Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  double Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  std::string fMCParticleModuleLabel, fPFParticleModuleLabel, fHitModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fTrackCalorimetryModuleLabel;
  bool fDebug;

  std::map<int,int> fHitsMap;
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap;

  TTree* fEventTree;

  // Tree variables

  int  _run;
  int  _subrun;
  int  _event;

  int _n_nu;

  std::map<std::string, VecVar*> nuVars = {
    { "nu_event_type", new InhVecVar<int>("nu_event_type") },
    { "nu_signal", new InhVecVar<bool>("nu_signal") },
    { "nu_true_en_dep", new InhVecVar<double>("nu_true_en_dep") },
    { "nu_ccnc", new InhVecVar<int>("nu_ccnc") },
    { "nu_mode", new InhVecVar<int>("nu_mode") },
    { "nu_int_type", new InhVecVar<int>("nu_int_type") },
    { "nu_w", new InhVecVar<double>("nu_w") },
    { "nu_x", new InhVecVar<double>("nu_x") },
    { "nu_y", new InhVecVar<double>("nu_y") },
    { "nu_q_sqr", new InhVecVar<double>("nu_q_sqr") },
    { "nu_pt", new InhVecVar<double>("nu_pt") },
    { "nu_theta", new InhVecVar<double>("nu_theta") },
    { "nu_e", new InhVecVar<double>("nu_e") },
    { "nu_vtx_x", new InhVecVar<double>("nu_vtx_x") },
    { "nu_vtx_y", new InhVecVar<double>("nu_vtx_y") },
    { "nu_vtx_z", new InhVecVar<double>("nu_vtx_z") },
  };
  };

sbnd::NCPiZeroAnalysis::NCPiZeroAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fMCParticleModuleLabel       = p.get<std::string>("MCParticleModuleLabel", "largeant");
    fPFParticleModuleLabel       = p.get<std::string>("PFParticleModuleLabel", "pandoraSCE");
    fHitModuleLabel              = p.get<std::string>("HitModuleLabel", "gaushit");
    fTrackModuleLabel            = p.get<std::string>("TrackModuleLabel", "pandoraSCETrack");
    fShowerModuleLabel           = p.get<std::string>("ShowerModuleLabel", "pandoraSCEShower");
    fTrackCalorimetryModuleLabel = p.get<std::string>("TrackCalorimetryModuleLabel", "pandoraSCECalo");
    fDebug                       = p.get<bool>("Debug", false);

    art::ServiceHandle<art::TFileService> fs;

    fEventTree = fs->make<TTree>("events","");
    fEventTree->Branch("run", &_run);
    fEventTree->Branch("subrun", &_subrun);
    fEventTree->Branch("event", &_event);

    fEventTree->Branch("n_nu", &_n_nu);

    for(auto const& [name, nuVar] : nuVars)
      {
	switch(nuVar->Identify())
	  {
	  case kBool:
	    fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<bool>*>(nuVar)->Var()));
	    break;
	  case kInt:
	    fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<int>*>(nuVar)->Var()));
	    break;
	  case kDouble:
	    fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<double>*>(nuVar)->Var()));
	    break;
	  case kUnknownVar:
	    break;
	  }
      }
  }

void sbnd::NCPiZeroAnalysis::analyze(const art::Event &e)
{
  ResetVars();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if(fDebug)
    std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  ClearMaps(e);
  SetupMaps(e);

  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();

  AnalyseNeutrinos(e, MCTruthHandles);

  // Fill the Tree
  fEventTree->Fill();
}

void sbnd::NCPiZeroAnalysis::ClearMaps(const art::Event &e)
{
  fHitsMap.clear();
  fRecoPFPMap.clear();
}

void sbnd::NCPiZeroAnalysis::SetupMaps(const art::Event &e)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  art::Handle<std::vector<recob::Hit> > hitsHandle;
  e.getByLabel(fHitModuleLabel,hitsHandle);

  for(unsigned hit_i = 0; hit_i < hitsHandle->size(); ++hit_i) {
    const art::Ptr<recob::Hit> hit(hitsHandle,hit_i);
    fHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
  }
}

void sbnd::NCPiZeroAnalysis::AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
	{
	  if(mct->Origin() != 1)
	    continue;

	  ++_n_nu;
	}
    }

  ResizeNeutrinoVectors(_n_nu);

  int nuCounter = 0;

  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
	{
	  if(mct->Origin() != 1)
	    continue;

	  const simb::MCNeutrino mcn = mct->GetNeutrino();
	  const simb::MCParticle nu  = mcn.Nu();

	  const bool nc = mcn.CCNC() == 1;
	  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

          unsigned pizeros = 0;

          for(int i = 0; i < mct->NParticles(); ++i)
            {
              const auto mcp = mct->GetParticle(i);

              if(mcp.PdgCode() == 111 && mcp.StatusCode() != 1)
                ++pizeros;
            }

	  const bool pizero = pizeros > 0;

	  if(nc && fv && pizero)
	    {
	      dynamic_cast<InhVecVar<int>*>(nuVars["nu_event_type"])->SetVal(nuCounter, kNCPiZero);
	      dynamic_cast<InhVecVar<bool>*>(nuVars["nu_signal"])->SetVal(nuCounter, true);
	    }
	  else 
	    {
	      dynamic_cast<InhVecVar<bool>*>(nuVars["nu_signal"])->SetVal(nuCounter, false);

	      if(nc && fv)
		  dynamic_cast<InhVecVar<int>*>(nuVars["nu_event_type"])->SetVal(nuCounter, kOtherNC);
	      else if(abs(nu.PdgCode()) == 14 && !nc && fv)
		dynamic_cast<InhVecVar<int>*>(nuVars["nu_event_type"])->SetVal(nuCounter, kCCNuMu);
	      else if(abs(nu.PdgCode()) == 12 && !nc && fv)
		dynamic_cast<InhVecVar<int>*>(nuVars["nu_event_type"])->SetVal(nuCounter, kCCNuE);
	      else if(!fv)
		dynamic_cast<InhVecVar<int>*>(nuVars["nu_event_type"])->SetVal(nuCounter, kNonFV);
	      else
		dynamic_cast<InhVecVar<int>*>(nuVars["nu_event_type"])->SetVal(nuCounter, kUnknownEv);
	    }

	  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
	  double trueEnDep = 0.;

	  for(auto const& mcp : MCParticleVec)
	    {
	      std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

	      for(auto const& ide : ides)
		trueEnDep += ide->energy / 1000.;
	    }

	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_true_en_dep"])->SetVal(nuCounter, trueEnDep);
	  dynamic_cast<InhVecVar<int>*>(nuVars["nu_ccnc"])->SetVal(nuCounter, mcn.CCNC());
	  dynamic_cast<InhVecVar<int>*>(nuVars["nu_mode"])->SetVal(nuCounter, mcn.Mode());
	  dynamic_cast<InhVecVar<int>*>(nuVars["nu_int_type"])->SetVal(nuCounter, mcn.InteractionType());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_w"])->SetVal(nuCounter, mcn.W());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_x"])->SetVal(nuCounter, mcn.X());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_y"])->SetVal(nuCounter, mcn.Y());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_q_sqr"])->SetVal(nuCounter, mcn.QSqr());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_pt"])->SetVal(nuCounter, mcn.Pt());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_theta"])->SetVal(nuCounter, mcn.Theta());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_e"])->SetVal(nuCounter, nu.E());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_vtx_x"])->SetVal(nuCounter, nu.Vx());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_vtx_y"])->SetVal(nuCounter, nu.Vy());
	  dynamic_cast<InhVecVar<double>*>(nuVars["nu_vtx_z"])->SetVal(nuCounter, nu.Vz());

	  ++nuCounter;
	}
    }
}

double sbnd::NCPiZeroAnalysis::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (objectHits.size() == 0) ? def_double : objectHitsMap[trackID]/static_cast<double>(objectHits.size());
}

double sbnd::NCPiZeroAnalysis::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::map<int, int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

  return (fHitsMap[trackID] == 0) ? def_double : objectHitsMap[trackID]/static_cast<double>(fHitsMap[trackID]);
}

bool sbnd::NCPiZeroAnalysis::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

void sbnd::NCPiZeroAnalysis::ResetVars()
{
  _run = -1; _subrun = -1; _event  = -1;

  _n_nu = 0;
}

void sbnd::NCPiZeroAnalysis::ResizeNeutrinoVectors(const int size)
{
  for(auto const& [name, nuVar] : nuVars)
    {
      nuVar->Resize(size);
    }
}

DEFINE_ART_MODULE(sbnd::NCPiZeroAnalysis)
