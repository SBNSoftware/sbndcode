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

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "sbnobj/Common/Reco/CRUMBSResult.h"

#include "NCPiZeroStructs.h"

#include <numeric>

constexpr int def_int       = -std::numeric_limits<int>::max();
constexpr size_t def_size   = std::numeric_limits<size_t>::max();
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
  void ClearMaps();
  void SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                 const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);

  void SetupBranches(VecVarMap &map);

  void ResizeVectors(VecVarMap &map, const int size);
  void ResizeSubVectors(VecVarMap &map, const int pos, const int size);

  void AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                     const art::Handle<std::vector<recob::PFParticle>> &pfpHandle);


  double Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);
  double Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int trackID);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  template<typename T>
  void FillElement(VecVar *vec, const int pos, const T value);
  template<typename T>
  void FillElement(VecVar *vec, const int posA, const int posB, const T value);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  std::string fMCParticleModuleLabel, fSliceModuleLabel, fPFParticleModuleLabel, fVertexModuleLabel,
    fHitModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fTrackCalorimetryModuleLabel,
    fCRUMBSModuleLabel;
  bool fDebug;

  std::map<int,int> fHitsMap;
  std::map<int, art::Ptr<recob::PFParticle>> fPFPMap;
  std::map<int, std::set<art::Ptr<recob::PFParticle>>> fRecoPFPMap;

  TTree* fEventTree;

  // Tree variables

  int  _run;
  int  _subrun;
  int  _event;

  int _n_nu;

  VecVarMap nuVars = {
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

  int _n_slc;

  VecVarMap slcVars = {
    { "slc_n_pfps", new InhVecVar<size_t>("slc_n_pfps") },
    { "slc_primary_pfp_id", new InhVecVar<size_t>("slc_primary_pfp_id") },
    { "slc_primary_pfp_pdg", new InhVecVar<int>("slc_primary_pfp_pdg") },
    { "slc_is_clear_cosmic", new InhVecVar<bool>("slc_is_clear_cosmic") },
    { "slc_n_primary_daughters", new InhVecVar<int>("slc_n_primary_daughters") },
    { "slc_vtx_x", new InhVecVar<double>("slc_vtx_x") },
    { "slc_vtx_y", new InhVecVar<double>("slc_vtx_y") },
    { "slc_vtx_z", new InhVecVar<double>("slc_vtx_z") },
    { "slc_is_fv", new InhVecVar<bool>("slc_is_fv") },
    { "slc_crumbs_score", new InhVecVar<double>("slc_crumbs_score") },
    { "slc_crumbs_nc_score", new InhVecVar<double>("slc_crumbs_nc_score") },
    { "slc_crumbs_ccnue_score", new InhVecVar<double>("slc_crumbs_ccnue_score") },
    { "slc_pfp_id", new InhVecVecVar<size_t>("slc_pfp_id") },
    { "slc_pfp_track_score", new InhVecVecVar<double>("slc_pfp_track_score") },
    { "slc_pfp_good_track", new InhVecVecVar<bool>("slc_pfp_good_track") },
    { "slc_pfp_good_shower", new InhVecVecVar<bool>("slc_pfp_good_shower") },
  };

};

sbnd::NCPiZeroAnalysis::NCPiZeroAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fMCParticleModuleLabel       = p.get<std::string>("MCParticleModuleLabel", "largeant");
    fSliceModuleLabel            = p.get<std::string>("SliceModuleLabel", "pandoraSCE");
    fPFParticleModuleLabel       = p.get<std::string>("PFParticleModuleLabel", "pandoraSCE");
    fVertexModuleLabel           = p.get<std::string>("VertexModuleLabel", "pandoraSCE");
    fHitModuleLabel              = p.get<std::string>("HitModuleLabel", "gaushit");
    fTrackModuleLabel            = p.get<std::string>("TrackModuleLabel", "pandoraSCETrack");
    fShowerModuleLabel           = p.get<std::string>("ShowerModuleLabel", "pandoraSCEShower");
    fTrackCalorimetryModuleLabel = p.get<std::string>("TrackCalorimetryModuleLabel", "pandoraSCECalo");
    fCRUMBSModuleLabel           = p.get<std::string>("CRUMBSModuleLabel", "crumbs");
    fDebug                       = p.get<bool>("Debug", false);

    art::ServiceHandle<art::TFileService> fs;

    fEventTree = fs->make<TTree>("events","");
    fEventTree->Branch("run", &_run);
    fEventTree->Branch("subrun", &_subrun);
    fEventTree->Branch("event", &_event);

    fEventTree->Branch("n_nu", &_n_nu);
    SetupBranches(nuVars);

    fEventTree->Branch("n_slc", &_n_slc);
    SetupBranches(slcVars);
  }

void sbnd::NCPiZeroAnalysis::SetupBranches(VecVarMap &map)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kOneD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<bool>*>(var)->Var()));
              break;
            case kInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<int>*>(var)->Var()));
              break;
            case kUInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<size_t>*>(var)->Var()));
              break;
            case kDouble:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVar<double>*>(var)->Var()));
              break;
            case kUnknownVar:
              break;
            }
        }
      else if(var->IdentifyVec() == kTwoD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<bool>*>(var)->Var()));
              break;
            case kInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<int>*>(var)->Var()));
              break;
            case kUInt:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<size_t>*>(var)->Var()));
              break;
            case kDouble:
              fEventTree->Branch(name.c_str(), &(dynamic_cast<InhVecVecVar<double>*>(var)->Var()));
              break;
            case kUnknownVar:
              break;
            }
        }
    }
}

void sbnd::NCPiZeroAnalysis::analyze(const art::Event &e)
{
  ResetVars();
  ClearMaps();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  if(fDebug)
    std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();


  // Get Hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitModuleLabel, hitHandle);
  if(!hitHandle.isValid()){
    std::cout << "Hit product " << fHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get Slices
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceModuleLabel, sliceHandle);
  if(!sliceHandle.isValid()){
    std::cout << "Slice product " << fSliceModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get PFParticles
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPFParticleModuleLabel, pfpHandle);
  if(!pfpHandle.isValid()){
    std::cout << "PFParticle product " << fPFParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  SetupMaps(e, hitHandle, pfpHandle);
  AnalyseNeutrinos(e, MCTruthHandles);
  AnalyseSlices(e, sliceHandle, pfpHandle);

  // Fill the Tree
  fEventTree->Fill();
}

void sbnd::NCPiZeroAnalysis::ClearMaps()
{
  fHitsMap.clear();
  fPFPMap.clear();
  fRecoPFPMap.clear();
}

void sbnd::NCPiZeroAnalysis::SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle,
                                       const art::Handle<std::vector<recob::PFParticle>> &pfpHandle)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  for(auto const& hit : hitVec)
    fHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;

  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  art::fill_ptr_vector(pfpVec, pfpHandle);

  for(auto const& pfp : pfpVec)
    fPFPMap[pfp->Self()] = pfp;
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

  ResizeVectors(nuVars, _n_nu);

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
          const bool av = VolumeCheck(nu.Position().Vect());
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
              FillElement(nuVars["nu_event_type"], nuCounter, (int) kNCPiZero);
              FillElement(nuVars["nu_signal"], nuCounter, true);
            }
          else
            {
              FillElement(nuVars["nu_signal"], nuCounter, false);

              if(nc && fv)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kOtherNC);
              else if(abs(nu.PdgCode()) == 14 && !nc && fv)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kCCNuMu);
              else if(abs(nu.PdgCode()) == 12 && !nc && fv)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kCCNuE);
              else if(!fv && av)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kNonFV);
              else if(!av)
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kDirt);
              else
                FillElement(nuVars["nu_event_type"], nuCounter, (int) kUnknownEv);
            }

          const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
          double trueEnDep = 0.;

          for(auto const& mcp : MCParticleVec)
            {
              std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

              for(auto const& ide : ides)
                trueEnDep += ide->energy / 1000.;
            }

          FillElement(nuVars["nu_true_en_dep"], nuCounter, trueEnDep);
          FillElement(nuVars["nu_ccnc"], nuCounter, mcn.CCNC());
          FillElement(nuVars["nu_mode"], nuCounter, mcn.Mode());
          FillElement(nuVars["nu_int_type"], nuCounter, mcn.InteractionType());
          FillElement(nuVars["nu_w"], nuCounter, mcn.W());
          FillElement(nuVars["nu_x"], nuCounter, mcn.X());
          FillElement(nuVars["nu_y"], nuCounter, mcn.Y());
          FillElement(nuVars["nu_q_sqr"], nuCounter, mcn.QSqr());
          FillElement(nuVars["nu_pt"], nuCounter, mcn.Pt());
          FillElement(nuVars["nu_theta"], nuCounter, mcn.Theta());
          FillElement(nuVars["nu_e"], nuCounter, nu.E());
          FillElement(nuVars["nu_vtx_x"], nuCounter, nu.Vx());
          FillElement(nuVars["nu_vtx_y"], nuCounter, nu.Vy());
          FillElement(nuVars["nu_vtx_z"], nuCounter, nu.Vz());

          ++nuCounter;
        }
    }
}

void sbnd::NCPiZeroAnalysis::AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                           const art::Handle<std::vector<recob::PFParticle>> &pfpHandle)
{
  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);

  _n_slc = sliceVec.size();
  ResizeVectors(slcVars, _n_slc);

  art::FindManyP<recob::PFParticle>                slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindOneP<recob::Vertex>                     pfpToVertices(pfpHandle, e, fVertexModuleLabel);
  art::FindOneP<sbn::CRUMBSResult>                 sliceToCRUMBS(sliceHandle, e, fCRUMBSModuleLabel);
  art::FindOneP<recob::Track>                      pfpToTrack(pfpHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>                     pfpToShower(pfpHandle, e, fShowerModuleLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpToMeta(pfpHandle, e, fPFParticleModuleLabel);

  for (auto&& [slcCounter, slc] : enumerate(sliceVec))
    {
      const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slc.key());
      FillElement(slcVars["slc_n_pfps"], slcCounter, pfps.size());

      if(pfps.size() == 0)
        {
          FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, true);
          continue;
        }

      const art::Ptr<recob::PFParticle> prim = GetPrimaryPFP(pfps);
      if(prim.isNull())
        continue;

      FillElement(slcVars["slc_primary_pfp_id"], slcCounter, prim->Self());
      FillElement(slcVars["slc_primary_pfp_pdg"], slcCounter, prim->PdgCode());
      FillElement(slcVars["slc_n_primary_daughters"], slcCounter, prim->NumDaughters());

      if(abs(prim->PdgCode()) == 13 || abs(prim->PdgCode()) == 11)
        FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, true);
      else
        FillElement(slcVars["slc_is_clear_cosmic"], slcCounter, false);

      const art::Ptr<recob::Vertex> vtx = pfpToVertices.at(prim.key());
      geo::Point_t vtxPos = vtx.isNonnull() ? vtx->position() : geo::Point_t(def_double, def_double, def_double);
      FillElement(slcVars["slc_vtx_x"], slcCounter, vtxPos.X());
      FillElement(slcVars["slc_vtx_y"], slcCounter, vtxPos.Y());
      FillElement(slcVars["slc_vtx_z"], slcCounter, vtxPos.Z());
      FillElement(slcVars["slc_is_fv"], slcCounter, VolumeCheck(vtxPos, 20., 5., 10., 50.));

      const art::Ptr<sbn::CRUMBSResult> crumbs = sliceToCRUMBS.at(slc.key());
      if(crumbs.isNonnull())
        {
          FillElement<double>(slcVars["slc_crumbs_score"], slcCounter, crumbs->score);
          FillElement<double>(slcVars["slc_crumbs_nc_score"], slcCounter, crumbs->ncscore);
          FillElement<double>(slcVars["slc_crumbs_ccnue_score"], slcCounter, crumbs->ccnuescore);
        }

      ResizeSubVectors(slcVars, slcCounter, prim->NumDaughters());

      for(auto&& [pfpCounter, id] : enumerate(prim->Daughters()))
        {
          const art::Ptr<recob::PFParticle> pfp = fPFPMap[id];
          FillElement(slcVars["slc_pfp_id"], slcCounter, pfpCounter, id);

          const art::Ptr<larpandoraobj::PFParticleMetadata> meta         = pfpToMeta.at(pfp.key());
          const larpandoraobj::PFParticleMetadata::PropertiesMap metaMap = meta->GetPropertiesMap();
          const std::map<std::string, float>::const_iterator scoreIter   = metaMap.find("TrackScore");

          if(scoreIter != metaMap.end())
            FillElement<double>(slcVars["slc_pfp_track_score"], slcCounter, pfpCounter, scoreIter->second);

          const art::Ptr<recob::Track> track   = pfpToTrack.at(pfp.key());
          const art::Ptr<recob::Shower> shower = pfpToShower.at(pfp.key());

          FillElement(slcVars["slc_pfp_good_track"], slcCounter, pfpCounter, track.isNonnull());
          FillElement(slcVars["slc_pfp_good_shower"], slcCounter, pfpCounter, shower.isNonnull());
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

bool sbnd::NCPiZeroAnalysis::VolumeCheck(const geo::Point_t &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const TVector3 posVec(pos.X(), pos.Y(), pos.Z());
  return VolumeCheck(posVec, walls, cath, front, back);
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

  _n_nu = 0; _n_slc  = 0;
}

void sbnd::NCPiZeroAnalysis::ResizeVectors(VecVarMap &map, const int size)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kOneD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              dynamic_cast<InhVecVar<bool>*>(var)->Assign(size, false);
              break;
            case kInt:
              dynamic_cast<InhVecVar<int>*>(var)->Assign(size, def_int);
              break;
            case kUInt:
              dynamic_cast<InhVecVar<size_t>*>(var)->Assign(size, def_size);
              break;
            case kDouble:
              dynamic_cast<InhVecVar<double>*>(var)->Assign(size, def_double);
              break;
            case kUnknownVar:
              break;
            }
        }
      else if(var->IdentifyVec() == kTwoD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              dynamic_cast<InhVecVecVar<bool>*>(var)->Resize(size);
              break;
            case kInt:
              dynamic_cast<InhVecVecVar<int>*>(var)->Resize(size);
              break;
            case kUInt:
              dynamic_cast<InhVecVecVar<size_t>*>(var)->Resize(size);
              break;
            case kDouble:
              dynamic_cast<InhVecVecVar<double>*>(var)->Resize(size);
              break;
            case kUnknownVar:
              break;
            }
        }
    }
}

void sbnd::NCPiZeroAnalysis::ResizeSubVectors(VecVarMap &map, const int pos, const int size)
{
  for(auto const& [name, var] : map)
    {
      if(var->IdentifyVec() == kTwoD)
        {
          switch(var->IdentifyVar())
            {
            case kBool:
              dynamic_cast<InhVecVecVar<bool>*>(var)->Assign(pos, size, false);
              break;
            case kInt:
              dynamic_cast<InhVecVecVar<int>*>(var)->Assign(pos, size, def_int);
              break;
            case kUInt:
              dynamic_cast<InhVecVecVar<size_t>*>(var)->Assign(pos, size, def_size);
              break;
            case kDouble:
              dynamic_cast<InhVecVecVar<double>*>(var)->Assign(pos, size, def_double);
              break;
            case kUnknownVar:
              break;
            }
        }
    }
}

template<typename T>
void sbnd::NCPiZeroAnalysis::FillElement(VecVar *vec, const int pos, const T value)
{
  dynamic_cast<InhVecVar<T>*>(vec)->SetVal(pos, value);
}

template<typename T>
void sbnd::NCPiZeroAnalysis::FillElement(VecVar *vec, const int posA, const int posB, const T value)
{
  dynamic_cast<InhVecVecVar<T>*>(vec)->SetVal(posA, posB, value);
}

art::Ptr<recob::PFParticle> sbnd::NCPiZeroAnalysis::GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps)
{
  for(auto pfp : pfps)
    if(pfp->IsPrimary())
      return pfp;

  return art::Ptr<recob::PFParticle>();
}

DEFINE_ART_MODULE(sbnd::NCPiZeroAnalysis)
