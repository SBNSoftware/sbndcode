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

#include <numeric>

constexpr int def_int       = -std::numeric_limits<int>::max();
constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NCPiZeroAnalysis;
}

enum EventType
  {
    kNCPiZero,
    kOtherNC,
    kCCNuMu,
    kCCNuE,
    kDirt,
    kNonFV,
    kCosmic,
    kBadRecoSignal,
    kUnknown = -1
  };

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
  std::vector<EventType> _nu_event_type;
  std::vector<bool>      _nu_signal;
  std::vector<double>    _nu_true_en_dep;
  std::vector<int>       _nu_ccnc;
  std::vector<int>       _nu_mode;
  std::vector<int>       _nu_int_type;
  std::vector<double>    _nu_w;
  std::vector<double>    _nu_x;
  std::vector<double>    _nu_y;
  std::vector<double>    _nu_q_sqr;
  std::vector<double>    _nu_pt;
  std::vector<double>    _nu_theta;
  std::vector<double>    _nu_e;
  std::vector<double>    _nu_vtx_x;
  std::vector<double>    _nu_vtx_y;
  std::vector<double>    _nu_vtx_z;
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

    fEventTree = fs->make<TTree>("pizeros","");
    fEventTree->Branch("run", &_run);
    fEventTree->Branch("subrun", &_subrun);
    fEventTree->Branch("event", &_event);

    fEventTree->Branch("_n_nu", &_n_nu);
    fEventTree->Branch("_nu_event_type", &_nu_event_type);
    fEventTree->Branch("_nu_signal", &_nu_signal);
    fEventTree->Branch("_nu_true_en_dep", &_nu_true_en_dep);
    fEventTree->Branch("_nu_ccnc", &_nu_ccnc);
    fEventTree->Branch("_nu_mode", &_nu_mode);
    fEventTree->Branch("_nu_int_type", &_nu_int_type);
    fEventTree->Branch("_nu_w", &_nu_w);
    fEventTree->Branch("_nu_x", &_nu_x);
    fEventTree->Branch("_nu_y", &_nu_y);
    fEventTree->Branch("_nu_q_sqr", &_nu_q_sqr);
    fEventTree->Branch("_nu_pt", &_nu_pt);
    fEventTree->Branch("_nu_theta", &_nu_theta);
    fEventTree->Branch("_nu_e", &_nu_e);
    fEventTree->Branch("_nu_vtx_x", &_nu_vtx_x);
    fEventTree->Branch("_nu_vtx_y", &_nu_vtx_y);
    fEventTree->Branch("_nu_vtx_z", &_nu_vtx_z);
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
	      _nu_event_type[nuCounter] = kNCPiZero;
	      _nu_signal[nuCounter]     = true;
	    }
	  else if(nc && fv)
	    _nu_event_type[nuCounter] = kOtherNC;
	  else if(abs(nu.PdgCode()) == 14 && !nc && fv)
	    _nu_event_type[nuCounter] = kCCNuMu;
	  else if(abs(nu.PdgCode()) == 12 && !nc && fv)
	    _nu_event_type[nuCounter] = kCCNuE;
	  else if(!fv)
	    _nu_event_type[nuCounter] = kNonFV;
	  else
	    _nu_event_type[nuCounter] = kUnknown;

	  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
	  double trueEnDep = 0.;

	  for(auto const& mcp : MCParticleVec)
	    {
	      std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

	      for(auto const& ide : ides)
		trueEnDep += ide->energy / 1000.;
	    }

	  _nu_true_en_dep[nuCounter] = trueEnDep;
	  _nu_ccnc[nuCounter]        = mcn.CCNC();
	  _nu_mode[nuCounter]        = mcn.Mode();
	  _nu_int_type[nuCounter]    = mcn.InteractionType();
	  _nu_w[nuCounter]           = mcn.W();
	  _nu_x[nuCounter]           = mcn.X();
	  _nu_y[nuCounter]           = mcn.Y();
	  _nu_q_sqr[nuCounter]       = mcn.QSqr();
	  _nu_pt[nuCounter]          = mcn.Pt();
	  _nu_theta[nuCounter]       = mcn.Theta();
	  _nu_e[nuCounter]           = nu.E();
	  _nu_vtx_x[nuCounter]       = nu.Vx();
	  _nu_vtx_y[nuCounter]       = nu.Vy();
	  _nu_vtx_z[nuCounter]       = nu.Vz();

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
  _nu_event_type.resize(size);
  _nu_signal.resize(size);
  _nu_true_en_dep.resize(size);
  _nu_ccnc.resize(size);
  _nu_mode.resize(size);
  _nu_int_type.resize(size);
  _nu_w.resize(size);
  _nu_x.resize(size);
  _nu_y.resize(size);
  _nu_q_sqr.resize(size);
  _nu_pt.resize(size);
  _nu_theta.resize(size);
  _nu_e.resize(size);
  _nu_vtx_x.resize(size);
  _nu_vtx_y.resize(size);
  _nu_vtx_z.resize(size);
}

DEFINE_ART_MODULE(sbnd::NCPiZeroAnalysis)
