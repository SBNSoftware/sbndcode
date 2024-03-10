////////////////////////////////////////////////////////////////////////
// Class:       NCPiZeroXSecTrees
// Plugin Type: analyzer
// File:        NCPiZeroXSecTrees_module.cc
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
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSpline.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"

#include "NCPiZeroStructs.h"

#include <numeric>

#include "CLHEP/Random/RandEngine.h"
#include "CLHEP/Random/RandGauss.h"

constexpr int def_int       = std::numeric_limits<int>::min();
constexpr float def_float   = -std::numeric_limits<float>::max();
constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NCPiZeroXSecTrees;
}

class sbnd::NCPiZeroXSecTrees : public art::EDAnalyzer {
public:
  explicit NCPiZeroXSecTrees(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NCPiZeroXSecTrees(NCPiZeroXSecTrees const&) = delete;
  NCPiZeroXSecTrees(NCPiZeroXSecTrees&&) = delete;
  NCPiZeroXSecTrees& operator=(NCPiZeroXSecTrees const&) = delete;
  NCPiZeroXSecTrees& operator=(NCPiZeroXSecTrees&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;
  virtual void beginSubRun(const art::SubRun& sr);
  virtual void endSubRun(const art::SubRun& sr);

  void ResetSubRunVars();
  void ResetEventVars();
  void ResetNuVars();
  void ResetSliceVars();

  void ClearMaps();

  void SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle);

  void AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void AnalyseMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct);

  void AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                     const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                     const art::Handle<std::vector<recob::Shower>> &showerHandle);

  void AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc,
                         const art::Handle<std::vector<recob::Slice>> &sliceHandle);

  void AnalysePFPs(const art::Event &e, const art::Ptr<recob::PFParticle> &prim, const std::vector<art::Ptr<recob::PFParticle>> &pfps,
                   const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                   const art::Handle<std::vector<recob::Shower>> &showerHandle);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  int GetTotalGenEvents(const art::Event &e);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps);

  double CorrectEnergy(const double &energy);

  int MomBin(const double &mom);
  int CosThetaBin(const double &costheta);
  int TwoDBins(const double &mom, const double &costheta);

  double ShowerEnergy(const art::Ptr<recob::Shower> &shower, const art::FindManyP<recob::Hit> &showerToHits);
  double TrackEnergy(const art::Ptr<recob::Track> &track, const art::FindManyP<anab::Calorimetry> &trackToCalos);
  TVector3 TrackDir(const art::Ptr<recob::Track> &track);

private:

  CLHEP::HepRandomEngine& fRandomEngine;

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

  art::InputTag fMCParticleModuleLabel, fSliceModuleLabel, fPFParticleModuleLabel, fVertexModuleLabel,
    fHitModuleLabel, fTrackModuleLabel, fCaloModuleLabel, fShowerModuleLabel, fCRUMBSModuleLabel,
    fPOTModuleLabel, fOpT0ModuleLabel, fRazzledModuleLabel;

  std::vector<art::InputTag> fEventWeightModuleLabels;
  bool fBeamOff;

  std::string fShowerEnergyCorrectionFileName;
  TProfile* fShowerEnergyCorrectionHist;

  std::map<const art::Ptr<simb::MCTruth>, int> fNuHitsMap;

  TTree* fSubRunTree;

  double _pot;
  int    _spills, _ngenevts;

  TTree *fNuTree, *fSliceTree;

  int  _run, _subrun, _event, _nu_id, _slc_id;

  int _event_type_incl, _event_type_0p0pi, _event_type_Np0pi,
    _event_type_cc, _event_mode, _true_mom_bin, _true_cos_theta_bin,
    _true_twod_bin, _reco_mom_bin, _reco_cos_theta_bin, _reco_twod_bin,
    _n_primary_razzled_muons, _n_primary_razzled_photons,
    _n_primary_razzled_pions_thresh, _n_primary_razzled_protons_thresh;

  double _pizero_mom, _cos_theta_pizero, _reco_pizero_mom, _reco_cos_theta_pizero,
    _reco_invariant_mass;

  bool _sel_incl, _sel_0p0pi, _sel_Np0pi, _sel_cc, _is_clear_cosmic, _is_fv,
    _best_pzc_good_kinematics, _all_other_trks_contained;

  float _crumbs_nc, _crumbs_ccnumu, _opt0_fracPE, _opt0_score, _comp;

  std::map<std::string, std::map<int, double>> genie_multisigma_universe_weights;

  std::map<std::string, std::vector<float>> _flux_weights, _genie_weights, _geant4_weights;
};

sbnd::NCPiZeroXSecTrees::NCPiZeroXSecTrees(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fRandomEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(createEngine(0, "HepJamesRandom", "SystematicWeightThrows"),
                                                                                     "HepJamesRandom", "SystematicWeightThrows", 101))
  , fMCParticleModuleLabel          (p.get<art::InputTag>("MCParticleModuleLabel"))
  , fSliceModuleLabel               (p.get<art::InputTag>("SliceModuleLabel"))
  , fPFParticleModuleLabel          (p.get<art::InputTag>("PFParticleModuleLabel"))
  , fVertexModuleLabel              (p.get<art::InputTag>("VertexModuleLabel"))
  , fHitModuleLabel                 (p.get<art::InputTag>("HitModuleLabel"))
  , fTrackModuleLabel               (p.get<art::InputTag>("TrackModuleLabel"))
  , fCaloModuleLabel                (p.get<art::InputTag>("CaloModuleLabel"))
  , fShowerModuleLabel              (p.get<art::InputTag>("ShowerModuleLabel"))
  , fCRUMBSModuleLabel              (p.get<art::InputTag>("CRUMBSModuleLabel"))
  , fPOTModuleLabel                 (p.get<art::InputTag>("POTModuleLabel"))
  , fOpT0ModuleLabel                (p.get<art::InputTag>("OpT0ModuleLabel"))
  , fRazzledModuleLabel             (p.get<art::InputTag>("RazzledModuleLabel"))
  , fEventWeightModuleLabels        (p.get<std::vector<art::InputTag>>("EventWeightModuleLabels"))
  , fBeamOff                        (p.get<bool>("BeamOff", false))
  , fShowerEnergyCorrectionFileName (p.get<std::string>("ShowerEnergyCorrectionFileName"))
  {
    cet::search_path sp("FW_SEARCH_PATH");
    std::string showerEnergyCorrectionFileFullPath;
    if (!sp.find_file(fShowerEnergyCorrectionFileName, showerEnergyCorrectionFileFullPath))
      throw cet::exception("NCPiZeroXSecTrees") << "Could not find shower energy correction file: \n"
                                                << fShowerEnergyCorrectionFileName;

    TFile* file = TFile::Open(showerEnergyCorrectionFileFullPath.c_str());
    fShowerEnergyCorrectionHist = (TProfile*) file->Get("hShowerEnergy2DRecoFractionalResolution_pfx");

    art::ServiceHandle<art::TFileService> fs;

    fSubRunTree = fs->make<TTree>("subruns","");

    fSubRunTree->Branch("pot", &_pot);
    fSubRunTree->Branch("spills", &_spills);
    fSubRunTree->Branch("ngenevts", &_ngenevts);

    fNuTree = fs->make<TTree>("neutrinos","");

    fNuTree->Branch("run", &_run);
    fNuTree->Branch("subrun", &_subrun);
    fNuTree->Branch("event", &_event);
    fNuTree->Branch("nu_id", &_nu_id);
    fNuTree->Branch("event_type_incl", &_event_type_incl);
    fNuTree->Branch("event_type_0p0pi", &_event_type_0p0pi);
    fNuTree->Branch("event_type_Np0pi", &_event_type_Np0pi);
    fNuTree->Branch("event_type_cc", &_event_type_cc);
    fNuTree->Branch("event_mode", &_event_mode);
    fNuTree->Branch("true_mom_bin", &_true_mom_bin);
    fNuTree->Branch("true_cos_theta_bin", &_true_cos_theta_bin);
    fNuTree->Branch("true_twod_bin", &_true_twod_bin);
    fNuTree->Branch("pizero_mom", &_pizero_mom);
    fNuTree->Branch("cos_theta_pizero", &_cos_theta_pizero);

    fSliceTree = fs->make<TTree>("slices","");

    fSliceTree->Branch("run", &_run);
    fSliceTree->Branch("subrun", &_subrun);
    fSliceTree->Branch("event", &_event);
    fSliceTree->Branch("slc_id", &_slc_id);
    fSliceTree->Branch("event_type_incl", &_event_type_incl);
    fSliceTree->Branch("event_type_0p0pi", &_event_type_0p0pi);
    fSliceTree->Branch("event_type_Np0pi", &_event_type_Np0pi);
    fSliceTree->Branch("event_type_cc", &_event_type_cc);
    fSliceTree->Branch("event_mode", &_event_mode);
    fSliceTree->Branch("true_mom_bin", &_true_mom_bin);
    fSliceTree->Branch("true_cos_theta_bin", &_true_cos_theta_bin);
    fSliceTree->Branch("true_twod_bin", &_true_twod_bin);
    fSliceTree->Branch("reco_mom_bin", &_reco_mom_bin);
    fSliceTree->Branch("reco_cos_theta_bin", &_reco_cos_theta_bin);
    fSliceTree->Branch("reco_twod_bin", &_reco_twod_bin);
    fSliceTree->Branch("pizero_mom", &_pizero_mom);
    fSliceTree->Branch("cos_theta_pizero", &_cos_theta_pizero);
    fSliceTree->Branch("reco_pizero_mom", &_reco_pizero_mom);
    fSliceTree->Branch("reco_cos_theta_pizero", &_reco_cos_theta_pizero);
    fSliceTree->Branch("reco_invariant_mass", &_reco_invariant_mass);
    fSliceTree->Branch("sel_incl", &_sel_incl);
    fSliceTree->Branch("sel_0p0pi", &_sel_0p0pi);
    fSliceTree->Branch("sel_Np0pi", &_sel_Np0pi);
    fSliceTree->Branch("sel_cc", &_sel_cc);
    fSliceTree->Branch("n_primary_razzled_muons", &_n_primary_razzled_muons);
    fSliceTree->Branch("n_primary_razzled_photons", &_n_primary_razzled_photons);
    fSliceTree->Branch("n_primary_razzled_pions_thresh", &_n_primary_razzled_pions_thresh);
    fSliceTree->Branch("n_primary_razzled_protons_thresh", &_n_primary_razzled_protons_thresh);
    fSliceTree->Branch("is_clear_cosmic", &_is_clear_cosmic);
    fSliceTree->Branch("is_fv", &_is_fv);
    fSliceTree->Branch("best_pzc_good_kinematics", &_best_pzc_good_kinematics);
    fSliceTree->Branch("all_other_trks_contained", &_all_other_trks_contained);
    fSliceTree->Branch("crumbs_nc", &_crumbs_nc);
    fSliceTree->Branch("crumbs_ccnumu", &_crumbs_ccnumu);
    fSliceTree->Branch("opt0_fracPE", &_opt0_fracPE);
    fSliceTree->Branch("opt0_score", &_opt0_score);
    fSliceTree->Branch("comp", &_comp);

    for(auto const& name : flux_weight_names)
      {
        _flux_weights[name] = std::vector<float>();
        fNuTree->Branch(name.c_str(), &_flux_weights[name]);
        fSliceTree->Branch(name.c_str(), &_flux_weights[name]);
      }

    _flux_weights["flux_weights_all"] = std::vector<float>();
    fNuTree->Branch("flux_weights_all", &_flux_weights["flux_weights_all"]);
    fSliceTree->Branch("flux_weights_all", &_flux_weights["flux_weights_all"]);

    for(auto const& name : genie_weight_names)
      {
        _genie_weights[name] = std::vector<float>();
        fNuTree->Branch(name.c_str(), &_genie_weights[name]);
        fSliceTree->Branch(name.c_str(), &_genie_weights[name]);

        if(name.find("multisigma") != std::string::npos)
          {
            _genie_weights[name + "_multisim"] = std::vector<float>();
            fNuTree->Branch((name + "_multisim").c_str(), &_genie_weights[name + "_multisim"]);
            fSliceTree->Branch((name + "_multisim").c_str(), &_genie_weights[name + "_multisim"]);
          }
      }

    _genie_weights["genie_weights_all"] = std::vector<float>();
    fNuTree->Branch("genie_weights_all", &_genie_weights["genie_weights_all"]);
    fSliceTree->Branch("genie_weights_all", &_genie_weights["genie_weights_all"]);

    for(auto const& name : geant4_weight_names)
      {
        _geant4_weights[name] = std::vector<float>();
        fNuTree->Branch(name.c_str(), &_geant4_weights[name]);
        fSliceTree->Branch(name.c_str(), &_geant4_weights[name]);
      }

    _geant4_weights["geant4_weights_all"] = std::vector<float>();
    fNuTree->Branch("geant4_weights_all", &_geant4_weights["geant4_weights_all"]);
    fSliceTree->Branch("geant4_weights_all", &_geant4_weights["geant4_weights_all"]);

    CLHEP::RandGauss randGauss(fRandomEngine, 0., 1.);

    for(auto const& name : genie_weight_names)
      {
        if(name.find("multisigma") != std::string::npos)
          {
            genie_multisigma_universe_weights[name] = std::map<int, double>();

            for(int univ = 0; univ < n_genieweight_univs; ++univ)
              genie_multisigma_universe_weights[name][univ] = randGauss.fire();
          }
      }
  }

void sbnd::NCPiZeroXSecTrees::beginSubRun(const art::SubRun &sr)
{
  ResetSubRunVars();

  if(fBeamOff)
    return;

  // Get POT
  art::Handle<sumdata::POTSummary> potHandle;
  sr.getByLabel(fPOTModuleLabel, potHandle);
  if(!potHandle.isValid()){
    std::cout << "POT product " << fPOTModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  _pot    = potHandle->totpot;
  _spills = potHandle->totspills;
}

void sbnd::NCPiZeroXSecTrees::endSubRun(const art::SubRun &sr)
{
  fSubRunTree->Fill();
}

void sbnd::NCPiZeroXSecTrees::analyze(const art::Event &e)
{
  ResetEventVars();
  ClearMaps();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  // Note this can only be accessed from the event object but is a subrun level quantity.
  // Hence, we override it every event but it is only filled in the subrun tree.
  _ngenevts = GetTotalGenEvents(e);

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

  // Get Tracks
  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackModuleLabel, trackHandle);
  if(!trackHandle.isValid()){
    std::cout << "Track product " << fTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  // Get Showers
  art::Handle<std::vector<recob::Shower>> showerHandle;
  e.getByLabel(fShowerModuleLabel, showerHandle);
  if(!showerHandle.isValid()){
    std::cout << "Shower product " << fShowerModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  SetupMaps(e, hitHandle);
  AnalyseSlices(e, sliceHandle, pfpHandle, trackHandle, showerHandle);
  AnalyseNeutrinos(e, MCTruthHandles);
}

void sbnd::NCPiZeroXSecTrees::ClearMaps()
{
  fNuHitsMap.clear();
}

void sbnd::NCPiZeroXSecTrees::SetupMaps(const art::Event &e, const art::Handle<std::vector<recob::Hit>> &hitHandle)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  for(auto const& hit : hitVec)
    {
      const int trackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      fNuHitsMap[mct]++;
    }
}

void sbnd::NCPiZeroXSecTrees::AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  _nu_id = 0;

  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          ResetNuVars();

          AnalyseMCTruth(e, mct);

          fNuTree->Fill();

          ++_nu_id;
        }
    }
}

void sbnd::NCPiZeroXSecTrees::AnalyseMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct)
{
  if(mct->Origin() == 2)
    {
      _event_type_incl  = NC::kCosmic;
      _event_type_0p0pi = NC::kCosmic;
      _event_type_Np0pi = NC::kCosmic;
      _event_type_cc    = CC::kCosmic;

      return;
    }
  else if(mct->Origin() != 1)
    {
      _event_type_incl  = NC::kUnknownEv;
      _event_type_0p0pi = NC::kUnknownEv;
      _event_type_Np0pi = NC::kUnknownEv;
      _event_type_cc    = CC::kUnknownEv;

      return;
    }

  const simb::MCNeutrino mcn = mct->GetNeutrino();
  const simb::MCParticle nu  = mcn.Nu();

  const bool nc = mcn.CCNC() == 1;
  const bool av = VolumeCheck(nu.Position().Vect());
  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

  art::FindManyP<simb::MCParticle> MCTruthToMCParticles( { mct }, e, fMCParticleModuleLabel);
  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(0);

  int protons = 0, charged_pions = 0, neutral_pions = 0;

  for(auto const& mcp : MCParticleVec)
    {
      if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
        {
          switch(abs(mcp->PdgCode()))
            {
            case 2212:
              if(mcp->P() > .4)
                ++protons;
              break;
            case 211:
              if(mcp->P() > .15)
                ++charged_pions;
              break;
            case 111:
              if(mcp->NumberDaughters() == 2)
                {
                  ++neutral_pions;
                  _pizero_mom = 1e3 * mcp->P();
                  _cos_theta_pizero = mcp->Pz() / mcp->P();

                  _true_mom_bin       = MomBin(_pizero_mom);
                  _true_cos_theta_bin = CosThetaBin(_cos_theta_pizero);
                  _true_twod_bin      = TwoDBins(_pizero_mom, _cos_theta_pizero);
                }
              break;
            default:
              break;
            }
        }
    }

  const bool pizero = neutral_pions == 1;

  if(nc && fv && pizero)
    {
      _event_type_incl = NC::kSignalNCPiZero;
      _event_type_cc   = CC::kNC;

      if(charged_pions == 0)
        {
          if(protons == 0)
            _event_type_0p0pi = NC::kSignalNCPiZero;
          else
            _event_type_0p0pi = NC::kOtherNCPiZero;

          if(protons > 0)
            _event_type_Np0pi = NC::kSignalNCPiZero;
          else
            _event_type_Np0pi = NC::kOtherNCPiZero;
        }
      else
        {
          _event_type_0p0pi = NC::kOtherNCPiZero;
          _event_type_Np0pi = NC::kOtherNCPiZero;
        }
    }
  else
    {
      if(nc && fv)
        {
          _event_type_incl  = NC::kOtherNC;
          _event_type_0p0pi = NC::kOtherNC;
          _event_type_Np0pi = NC::kOtherNC;
          _event_type_cc    = CC::kNC;
        }
      else if(abs(nu.PdgCode()) == 14 && !nc && fv)
        {
          _event_type_incl  = NC::kCCNuMu;
          _event_type_0p0pi = NC::kCCNuMu;
          _event_type_Np0pi = NC::kCCNuMu;

          if(pizero)
            _event_type_cc = CC::kSignalCCPiZero;
          else
            _event_type_cc = CC::kOtherCCNuMu;
        }
      else if(abs(nu.PdgCode()) == 12 && !nc && fv)
        {
          _event_type_incl  = NC::kCCNuE;
          _event_type_0p0pi = NC::kCCNuE;
          _event_type_Np0pi = NC::kCCNuE;
          _event_type_cc    = CC::kCCNuE;
        }
      else if(!fv && av)
        {
          _event_type_incl  = NC::kNonFV;
          _event_type_0p0pi = NC::kNonFV;
          _event_type_Np0pi = NC::kNonFV;
          _event_type_cc    = CC::kNonFV;
        }
      else if(!av)
        {
          _event_type_incl  = NC::kDirt;
          _event_type_0p0pi = NC::kDirt;
          _event_type_Np0pi = NC::kDirt;
          _event_type_cc    = CC::kDirt;
        }
      else
        {
          _event_type_incl  = NC::kUnknownEv;
          _event_type_0p0pi = NC::kUnknownEv;
          _event_type_Np0pi = NC::kUnknownEv;
          _event_type_cc    = CC::kUnknownEv;
        }
    }

  _event_mode = mcn.Mode();

  for(auto const& weightModuleLabel : fEventWeightModuleLabels)
    {
      art::FindManyP<sbn::evwgh::EventWeightMap> MCTruthToWeights( { mct }, e, weightModuleLabel);
      const std::vector<art::Ptr<sbn::evwgh::EventWeightMap>> ewms = MCTruthToWeights.at(0);

      int n_univs = 0;
      if(weightModuleLabel == "fluxweight")
        n_univs = n_fluxweight_univs;
      else if(weightModuleLabel == "systtools")
        n_univs = n_genieweight_univs;
      else if(weightModuleLabel == "geant4weight")
        n_univs = n_geant4weight_univs;

      std::vector<float> all(n_univs, 1.);

      for(auto const& ewm : ewms)
        {
          for(auto const& [ name, weights ] : *ewm)
            {
              if(weightModuleLabel == "fluxweight")
                {
                  _flux_weights[name] = weights;

                  for(int univ = 0; univ < n_univs; ++univ)
                    all[univ] *= weights[univ];
                }
              else if(weightModuleLabel == "geant4weight")
                {
                  _geant4_weights[name] = weights;
                }
              else if(weightModuleLabel == "systtools" && name.find("multisim") != std::string::npos)
                {
                  _genie_weights[name] = weights;

                  for(int univ = 0; univ < n_univs; ++univ)
                    all[univ] *= weights[univ];
                }
              else if(weightModuleLabel == "systtools" && name.find("multisigma") != std::string::npos)
                {
                  _genie_weights[name] = weights;

                  std::vector<float> thrown_weights(n_univs, 1.);

                  if(weights.size() == 6)
                    {
                      double multisigma_sigmas[6] = { -1, 1, -2, 2, -3, 3 };
                      double multisigma_vals[6]   = { weights[0], weights[1], weights[2], weights[3], weights[4], weights[5] };

                      TSpline3 *spline = new TSpline3(Form("%s_spline", name.c_str()), multisigma_sigmas, multisigma_vals, 6);

                      for(int univ = 0; univ < n_univs; ++univ)
                        {
                          thrown_weights[univ] = spline->Eval(genie_multisigma_universe_weights[name][univ]);
                          all[univ] *= thrown_weights[univ];
                        }
                    }
                  else if(weights.size() == 1)
                    {
                      for(int univ = 0; univ < n_univs; ++univ)
                        {
                          thrown_weights[univ] = 1 + (weights[0] - 1) * 2 * genie_multisigma_universe_weights[name][univ];
                          all[univ] *= thrown_weights[univ];
                        }
                    }

                  _genie_weights[name + "_multisim"] = thrown_weights;
                }
            }
        }

      if(weightModuleLabel == "fluxweight")
        _flux_weights["flux_weights_all"] = all;
      else if(weightModuleLabel == "systtools")
        _genie_weights["genie_weights_all"] = all;
    }
}

void sbnd::NCPiZeroXSecTrees::AnalyseSlices(const art::Event &e, const art::Handle<std::vector<recob::Slice>> &sliceHandle,
                                            const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                                            const art::Handle<std::vector<recob::Shower>> &showerHandle)
{
  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);

  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFParticleModuleLabel);
  art::FindOneP<recob::Vertex>      pfpToVertices(pfpHandle, e, fVertexModuleLabel);
  art::FindOneP<sbn::CRUMBSResult>  slicesToCRUMBS(sliceHandle, e, fCRUMBSModuleLabel);
  art::FindManyP<sbn::OpT0Finder>   slicesToOpT0(sliceHandle, e, fOpT0ModuleLabel);
  art::FindManyP<recob::Hit>        slicesToHits(sliceHandle, e, fSliceModuleLabel);

  _slc_id = -1;

  for(auto const &slc : sliceVec)
    {
      ResetSliceVars();
      ++_slc_id;

      const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slc.key());
      const std::vector<art::Ptr<recob::Hit>> hits = slicesToHits.at(slc.key());

      if(pfps.size() == 0)
        {
          _is_clear_cosmic = true;
          continue;
        }

      const art::Ptr<recob::PFParticle> prim = GetPrimaryPFP(pfps);

      if(prim.isNull())
        {
          _is_clear_cosmic = true;
          continue;
        }

      _is_clear_cosmic = abs(prim->PdgCode()) == 13 || abs(prim->PdgCode()) == 11;

      const art::Ptr<recob::Vertex> vtx = pfpToVertices.at(prim.key());
      geo::Point_t vtxPos = vtx.isNonnull() ? vtx->position() : geo::Point_t(def_double, def_double, def_double);
      _is_fv = VolumeCheck(vtxPos, 20., 5., 10., 50.);

      const art::Ptr<sbn::CRUMBSResult> crumbs = slicesToCRUMBS.at(slc.key());
      if(crumbs.isNonnull())
        {
          _crumbs_nc     = crumbs->ncscore;
          _crumbs_ccnumu = crumbs->ccnumuscore;
        }

      std::vector<art::Ptr<sbn::OpT0Finder>> opT0Vec = slicesToOpT0.at(slc.key());
      if(opT0Vec.size() > 0)
        {
          std::sort(opT0Vec.begin(), opT0Vec.end(),
                    [](auto const& a, auto const& b)
                    { return a->score > b->score; });

          _opt0_score  = opT0Vec[0]->score;
          _opt0_fracPE = (opT0Vec[0]->hypoPE - opT0Vec[0]->measPE) / opT0Vec[0]->measPE;
        }

      AnalysePFPs(e, prim, pfps, pfpHandle, trackHandle, showerHandle);

      const bool common_sel = !_is_clear_cosmic && _is_fv && _n_primary_razzled_photons > 1 && _best_pzc_good_kinematics;

      _sel_incl = common_sel && _crumbs_nc > -0.005 && _n_primary_razzled_muons == 0
        && _opt0_fracPE < 0.756 && _opt0_fracPE > -0.7 && _opt0_score > 5 && _all_other_trks_contained;

      _sel_0p0pi = common_sel && _crumbs_nc > -0.005 && _n_primary_razzled_muons == 0
        && _opt0_fracPE < 0.844 && _opt0_fracPE > -0.7 && _opt0_score > 125 && _n_primary_razzled_pions_thresh == 0
        && _n_primary_razzled_protons_thresh == 0;

      _sel_Np0pi = common_sel && _crumbs_nc > -0.005 && _n_primary_razzled_muons == 0
        && _opt0_fracPE < 0.836 && _opt0_fracPE > -0.376 && _opt0_score > 210 && _n_primary_razzled_pions_thresh == 0
        && _n_primary_razzled_protons_thresh > 0;

      _sel_cc = common_sel && _crumbs_ccnumu > 0 && _n_primary_razzled_muons == 1
        && _opt0_fracPE < 0.836 && _opt0_fracPE > -0.376 && _opt0_score > 200;

      AnalyseSliceTruth(e, slc, sliceHandle);

      fSliceTree->Fill();
    }
}

void sbnd::NCPiZeroXSecTrees::AnalysePFPs(const art::Event &e, const art::Ptr<recob::PFParticle> &prim, const std::vector<art::Ptr<recob::PFParticle>> &pfps,
                                          const art::Handle<std::vector<recob::PFParticle>> &pfpHandle, const art::Handle<std::vector<recob::Track>> &trackHandle,
                                          const art::Handle<std::vector<recob::Shower>> &showerHandle)
{
  art::FindOneP<recob::Track>       pfpToTrack(pfpHandle, e, fTrackModuleLabel);
  art::FindOneP<recob::Shower>      pfpToShower(pfpHandle, e, fShowerModuleLabel);
  art::FindManyP<recob::Hit>        showerToHits(showerHandle, e, fShowerModuleLabel);
  art::FindOneP<sbn::MVAPID>        pfpToRazzled(pfpHandle, e, fRazzledModuleLabel);
  art::FindManyP<anab::Calorimetry> trackToCalos(trackHandle, e, fCaloModuleLabel);

  _n_primary_razzled_muons        = 0; _n_primary_razzled_photons        = 0;
  _n_primary_razzled_pions_thresh = 0; _n_primary_razzled_protons_thresh = 0;

  struct RazzPhot
  {
    size_t   index;
    double   shwEnergy;
    TVector3 trkDir;
  };

  std::vector<RazzPhot> razzled_photons;

  for(auto&& [pfp_i, pfp] : enumerate(pfps))
    {
      const bool primary_child = prim->Self() == pfp->Parent();

      const art::Ptr<recob::Track> track   = pfpToTrack.at(pfp.key());
      const art::Ptr<recob::Shower> shower = pfpToShower.at(pfp.key());

      const double shwEnergy = ShowerEnergy(shower, showerToHits);
      const double trkEnergy = TrackEnergy(track, trackToCalos);

      const double pfpenergy = pfp->PdgCode() == 13 ? trkEnergy : pfp->PdgCode() == 11 ? shwEnergy : def_double;

      const art::Ptr<sbn::MVAPID> razzled = pfpToRazzled.at(pfp.key());
      const int razzledpdg = razzled.isNull() ? -1 : razzled->BestPDG();

      if(primary_child)
        {
          if(razzledpdg == 13)
            ++_n_primary_razzled_muons;
          if(razzledpdg == 22)
            {
              ++_n_primary_razzled_photons;

              const TVector3 trkDir = TrackDir(track);
              const RazzPhot phot = { pfp_i, shwEnergy, trkDir };
              razzled_photons.push_back(phot);
            }
          if(razzledpdg == 211 && pfpenergy > 32.1)
            ++_n_primary_razzled_pions_thresh;
          if(razzledpdg == 2212 && pfpenergy > 32.7)
            ++_n_primary_razzled_protons_thresh;
        }
    }

  double bestInvMass = std::numeric_limits<double>::max();
  size_t best_candidate_0 = -1, best_candidate_1 = -1;

  for(auto&& [ i, phot_0 ] : enumerate(razzled_photons))
    {
      for(auto&& [ j, phot_1 ] : enumerate(razzled_photons))
        {
          if(j <= i)
            continue;

          const double energy_0 = CorrectEnergy(phot_0.shwEnergy);
          const double energy_1 = CorrectEnergy(phot_1.shwEnergy);

          const TVector3 dir_0 = phot_0.trkDir;
          const TVector3 dir_1 = phot_1.trkDir;

          const bool goodKinematics       = !(dir_0.X() == -999 || dir_1.X() == -999 || dir_0.X() == def_float || dir_1.X() == def_float
                                              || energy_0 < 0 || energy_1 < 0);
          const double cosThetaGammaGamma = dir_0.Dot(dir_1) / (dir_0.Mag() * dir_1.Mag());
          const TVector3 pizeroDir        = (energy_0 * dir_0) + (energy_1 * dir_1);
          const double invariantMass      = sqrt(2 * energy_0 * energy_1 * (1 - cosThetaGammaGamma));
          const double pizeroMom          = pizeroDir.Mag();
          const double pizeroCosTheta     = pizeroDir.Z() / pizeroMom;

          if(goodKinematics && abs(134.9769 - invariantMass) < bestInvMass)
            {
              _reco_pizero_mom          = pizeroMom;
              _reco_cos_theta_pizero    = pizeroCosTheta;
              _reco_invariant_mass      = invariantMass;
              _best_pzc_good_kinematics = true;
              best_candidate_0          = phot_0.index;
              best_candidate_1          = phot_1.index;


              _reco_mom_bin       = MomBin(_reco_pizero_mom);
              _reco_cos_theta_bin = CosThetaBin(_reco_cos_theta_pizero);
              _reco_twod_bin      = TwoDBins(_reco_pizero_mom, _reco_cos_theta_pizero);
            }
        }
    }

  _all_other_trks_contained = true;

  for(auto&& [pfp_i, pfp] : enumerate(pfps))
    {
      if(pfp_i == best_candidate_0 || pfp_i == best_candidate_1)
        continue;

      if(pfp->PdgCode() != 13)
        continue;

      const art::Ptr<recob::Track> track = pfpToTrack.at(pfp.key());

      if(track.isNull())
        {
          _all_other_trks_contained = false;
          continue;
        }

      const geo::Point_t start = track->Start();
      const geo::Point_t end   = track->End();

      _all_other_trks_contained &= VolumeCheck(start, 10, 0, 10, 10) && VolumeCheck(end, 10, 0, 10, 10);
    }
}

void sbnd::NCPiZeroXSecTrees::AnalyseSliceTruth(const art::Event &e, const art::Ptr<recob::Slice> &slc,
                                                const art::Handle<std::vector<recob::Slice>> &sliceHandle)
{
  art::FindManyP<recob::Hit> slicesToHits(sliceHandle, e, fSliceModuleLabel);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  const std::vector<art::Ptr<recob::Hit>> sliceHits = slicesToHits.at(slc.key());

  std::map<int, int> objectHitMap;
  for(auto const &hit : sliceHits)
    objectHitMap[TruthMatchUtils::TrueParticleID(clockData, hit, true)]++;

  std::map<const art::Ptr<simb::MCTruth>, int> mcTruthHitMap;
  for(auto const& [trackID, nhits] : objectHitMap)
    {
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
      mcTruthHitMap[mct] += nhits;
    }

  int maxHits = def_int;
  art::Ptr<simb::MCTruth> bestMCT = art::Ptr<simb::MCTruth>();

  for(auto const& [mct, nhits] : mcTruthHitMap)
    {
      if(nhits > maxHits)
        {
          maxHits = nhits;
          bestMCT = mct;
        }
    }

  _comp = fNuHitsMap[bestMCT] == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(fNuHitsMap[bestMCT]);

  if(fBeamOff)
    {
      _event_type_incl  = NC::kCosmic;
      _event_type_0p0pi = NC::kCosmic;
      _event_type_Np0pi = NC::kCosmic;
      _event_type_cc    = CC::kCosmic;
    }
  else if(bestMCT.isNonnull())
    AnalyseMCTruth(e, bestMCT);
  else
    {
      _event_type_incl  = NC::kFailedTruthMatch;
      _event_type_0p0pi = NC::kFailedTruthMatch;
      _event_type_Np0pi = NC::kFailedTruthMatch;
      _event_type_cc    = CC::kFailedTruthMatch;
    }
}

bool sbnd::NCPiZeroXSecTrees::VolumeCheck(const geo::Point_t &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const TVector3 posVec(pos.X(), pos.Y(), pos.Z());
  return VolumeCheck(posVec, walls, cath, front, back);
}

bool sbnd::NCPiZeroXSecTrees::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

void sbnd::NCPiZeroXSecTrees::ResetSubRunVars()
{
  _pot = 0.; _spills = 0; _ngenevts = 0;
}

void sbnd::NCPiZeroXSecTrees::ResetEventVars()
{
  _run = -1; _subrun = -1; _event  = -1;
}

void sbnd::NCPiZeroXSecTrees::ResetNuVars()
{
  _event_type_incl    = -1; _event_type_0p0pi = -1;
  _event_type_Np0pi   = -1; _event_type_cc    = -1;
  _event_mode         = -1; _true_mom_bin     = -1;
  _true_cos_theta_bin = -1; _true_twod_bin    = -1;

  _pizero_mom = def_double; _cos_theta_pizero = def_double;

  for(auto&& [ name, vec ] : _flux_weights)
    vec = std::vector<float>();

  for(auto&& [ name, vec ] : _genie_weights)
    vec = std::vector<float>();

  for(auto&& [ name, vec ] : _geant4_weights)
    vec = std::vector<float>();
}

void sbnd::NCPiZeroXSecTrees::ResetSliceVars()
{
  ResetNuVars();

  _reco_mom_bin                     = -1; _reco_cos_theta_bin             = -1;
  _reco_twod_bin                    = -1; _n_primary_razzled_muons        = -1;
  _n_primary_razzled_photons        = -1; _n_primary_razzled_pions_thresh = -1;
  _n_primary_razzled_protons_thresh = -1;

  _reco_pizero_mom     = def_double; _reco_cos_theta_pizero = def_double;
  _reco_invariant_mass = def_double;

  _sel_incl                 = false; _sel_0p0pi                = false;
  _sel_Np0pi                = false; _sel_cc                   = false;
  _is_clear_cosmic          = false; _is_fv                    = false;
  _best_pzc_good_kinematics = false; _all_other_trks_contained = false;

  _crumbs_nc   = def_float; _crumbs_ccnumu = def_float;
  _opt0_fracPE = def_float; _opt0_score    = def_float;
  _comp        = def_float;
}

int sbnd::NCPiZeroXSecTrees::GetTotalGenEvents(const art::Event &e)
{
  int nGenEvt = 0;
  for (const art::ProcessConfiguration &process: e.processHistory()) {
    std::optional<fhicl::ParameterSet> genConfig = e.getProcessParameterSet(process.processName());
    if (genConfig && genConfig->has_key("source") && genConfig->has_key("source.maxEvents") && genConfig->has_key("source.module_type") ) {
      int maxEvents = genConfig->get<int>("source.maxEvents");
      std::string moduleType = genConfig->get<std::string>("source.module_type");
      if (moduleType == "EmptyEvent") {
        nGenEvt += maxEvents;
      }
    }
  }

  return nGenEvt;
}

art::Ptr<recob::PFParticle> sbnd::NCPiZeroXSecTrees::GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfps)
{
  for(auto pfp : pfps)
    if(pfp->IsPrimary())
      return pfp;

  return art::Ptr<recob::PFParticle>();
}

double sbnd::NCPiZeroXSecTrees::CorrectEnergy(const double &energy)
{
  const int bin = fShowerEnergyCorrectionHist->FindBin(energy);

  return energy * (1 - fShowerEnergyCorrectionHist->GetBinContent(bin));
}

int sbnd::NCPiZeroXSecTrees::MomBin(const double &mom)
{
  int bin = -1;

  if(mom < 0.)
    bin = 0;
  else if(mom < 60.)
    bin = 1;
  else if(mom < 120.)
    bin = 2;
  else if(mom < 180.)
    bin = 3;
  else if(mom < 240.)
    bin = 4;
  else if(mom < 300.)
    bin = 5;
  else if(mom < 400.)
    bin = 6;
  else if(mom < 600.)
    bin = 7;
  else if(mom < 1000.)
    bin = 8;
  else
    bin = 9;

  return bin;
}

int sbnd::NCPiZeroXSecTrees::CosThetaBin(const double &costheta)
{
  int bin = -1;

  if(costheta < -1.)
    bin = 0;
  else if(costheta < -0.5)
    bin = 1;
  else if(costheta < 0.)
    bin = 2;
  else if(costheta < 0.2)
    bin = 3;
  else if(costheta < 0.4)
    bin = 4;
  else if(costheta < 0.6)
    bin = 5;
  else if(costheta < 0.8)
    bin = 6;
  else if(costheta < 0.9)
    bin = 7;
  else if(costheta < 0.95)
    bin = 8;
  else if(costheta < 1.)
    bin = 9;
  else
    bin = 10;

  return bin;
}

int sbnd::NCPiZeroXSecTrees::TwoDBins(const double &mom, const double &costheta)
{
  const int momBin      = MomBin(mom);
  const int cosThetaBin = CosThetaBin(costheta);

  return 10 * cosThetaBin + momBin;
}

double sbnd::NCPiZeroXSecTrees::ShowerEnergy(const art::Ptr<recob::Shower> &shower, const art::FindManyP<recob::Hit> &showerToHits)
{
  if(shower.isNull())
    return def_double;

  const std::vector<art::Ptr<recob::Hit>> hits = showerToHits.at(shower.key());

  std::array<int, 3> planeHits = { 0, 0, 0 };

  for(auto const& hit : hits)
    planeHits[hit->WireID().Plane]++;

  int bestPlane = -1;

  if(planeHits[2] >= planeHits[1] && planeHits[2] >= planeHits[0])
    bestPlane = 2;
  else if(planeHits[0] >= planeHits[1])
    bestPlane = 0;
  else
    bestPlane = 1;

  return shower->Energy()[bestPlane];
}

double sbnd::NCPiZeroXSecTrees::TrackEnergy(const art::Ptr<recob::Track> &track, const art::FindManyP<anab::Calorimetry> &trackToCalos)
{
  if(track.isNull())
    return def_double;

  const std::vector<art::Ptr<anab::Calorimetry>> calos = trackToCalos.at(track.key());
  const size_t maxHits = calos.size() != 3 ? -1 : std::max({calos[0]->dEdx().size(), calos[1]->dEdx().size(), calos[2]->dEdx().size()});
  const int bestPlane  = calos.size() != 3 ? -1 : (calos[2]->dEdx().size() == maxHits) ? 2 : (calos[0]->dEdx().size() == maxHits) ? 0 :
    (calos[1]->dEdx().size() == maxHits) ? 1 : -1;

  const std::vector<float> &dEdx  = calos[bestPlane]->dEdx();
  const std::vector<float> &pitch = calos[bestPlane]->TrkPitchVec();

  double ke = 0.;

  for(size_t i = 0; i < dEdx.size(); ++i)
    ke += dEdx[i] * pitch[i];

  return ke;
}

TVector3 sbnd::NCPiZeroXSecTrees::TrackDir(const art::Ptr<recob::Track> &track)
{
  if(track.isNull())
    return TVector3(def_float, def_float, def_float);

  const geo::Vector_t dir = track->StartDirection();

  return TVector3(dir.X(), dir.Y(), dir.Z());
}

DEFINE_ART_MODULE(sbnd::NCPiZeroXSecTrees)
