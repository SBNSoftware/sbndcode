///////////////////////////////////////////////////////////////////////////
// Class:        RecoEff                                                 //
// Module Type:  Analyser                                                //
// File:         RecoEff_module.cc                                       //
// Author:       Henry Lay h.lay@lancaster.ac.uk                         //
//                                                                       //
// Purpose:      The module saves truth information regarding the        //
//               primary products of neutrino interactions. It also      //
//               saves their status in reconstruction. This allows       //
//               for plotting of a series of plots of different          //
//               reconstruction metrics.                                 //
///////////////////////////////////////////////////////////////////////////

//Standard
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//Reco Base
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

//Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

//Root
#include "art_root_io/TFileService.h"
#include "TTree.h"

//SBNDCODE
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

constexpr int def_int     = -999;
constexpr float def_float = -999.0f;

class RecoEff;

class RecoEff : public art::EDAnalyzer {
public:
  explicit RecoEff(fhicl::ParameterSet const& pset);

  RecoEff(RecoEff const&) = delete;
  RecoEff(RecoEff&&) = delete;
  RecoEff& operator=(RecoEff const&) = delete;
  RecoEff& operator=(RecoEff&&) = delete;

  void analyze(art::Event const &e) override;

private:

  void ClearMaps();
  void ResetData();
  void SetupMaps(art::Event const &e);
  void ReconstructionProcessor(art::Event const &e);
  void TruthProcessor(art::Event const &e);
  void FillTrueParticle(const simb::MCParticle *particle, const bool pionChild);

  art::Ptr<recob::PFParticle> GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfpVec);
  art::Ptr<recob::PFParticle> GetPFP(art::Event const &e, long unsigned int const &id);

  float Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID);
  float TrueTrackLength(const simb::MCParticle *particle);


  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  sbnd::TPCGeoAlg fTPCGeo;

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  
  std::string fNuGenModuleLabel, fLArGeantModuleLabel, fSliceModuleLabel,
    fPFParticleModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fHitsModuleLabel;

  float fXEdgeCut, fYEdgeCut, fZFrontCut, fZBackCut,
    fXEdgeCutShowers, fYEdgeCutShowers, fZFrontCutShowers, fZBackCutShowers,
    fCathodeXCut;

  float fXOutEdge, fYEdge, fZFrontEdge, fZBackEdge,
    fXOutEdgeShowers, fYEdgeShowers, fZFrontEdgeShowers, fZBackEdgeShowers,
    fXCathodeEdge;

  TTree *fParticleTree;

  int fMC_trackID, fMC_PDG, fReco_nTracks, fReco_nShowers, fReco_primary_nTracks, fReco_primary_nShowers;
  float fMC_x0, fMC_y0, fMC_z0, fMC_xEnd, fMC_yEnd, fMC_zEnd, fMC_pX0,
    fMC_pY0, fMC_pZ0, fMC_energy0, fMC_momentum, fMC_pXEnd, fMC_pYEnd, fMC_pZEnd,
    fMC_energyEnd, fMC_mass, fMC_theta_xy, fMC_theta_yz, fMC_theta_xz, fMC_length, 
    fReco_showerPurity, fReco_showerCompleteness, fReco_trackPurity, 
    fReco_trackCompleteness, fReco_trackLength, fReco_showerdEdx,
    fReco_primary_showerPurity, fReco_primary_showerCompleteness, fReco_primary_trackPurity,
    fReco_primary_trackCompleteness, fReco_primary_trackLength, fReco_primary_showerdEdx;
  bool fReco_isReconstructed, fReco_primary_isReconstructed;

  std::map<int,int> fNShowersMap, fNTracksMap, fNPrimaryShowersMap, fNPrimaryTracksMap;
  std::map<int,float> fShowerCompMap, fShowerPurMap, fTrackCompMap,
    fTrackPurMap, fTrackLengthMap, fShowerdEdxMap, fPrimaryShowerCompMap,
    fPrimaryShowerPurMap, fPrimaryTrackCompMap, fPrimaryTrackPurMap,
    fPrimaryTrackLengthMap, fPrimaryShowerdEdxMap;
  std::map<int,int> fHitsMap;
};


RecoEff::RecoEff(fhicl::ParameterSet const &pset)
  : EDAnalyzer{pset},
  fNuGenModuleLabel (pset.get<std::string>("NuGenModuleLabel")),
  fLArGeantModuleLabel (pset.get<std::string>("LArGeantModuleLabel")),
  fSliceModuleLabel (pset.get<std::string>("SliceModuleLabel")),
  fPFParticleModuleLabel (pset.get<std::string>("PFParticleModuleLabel")),
  fTrackModuleLabel (pset.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel (pset.get<std::string>("ShowerModuleLabel")),
  fHitsModuleLabel (pset.get<std::string>("HitsModuleLabel")),
  fXEdgeCut (pset.get<float>("XEdgeCut",15.f)),
  fYEdgeCut (pset.get<float>("YEdgeCut",15.f)),
  fZFrontCut (pset.get<float>("ZFrontCut",30.f)),
  fZBackCut (pset.get<float>("ZBackCut",65.f)),
  fXEdgeCutShowers (pset.get<float>("XEdgeCutShowers",25.f)),
  fYEdgeCutShowers (pset.get<float>("YEdgeCutShowers",25.f)),
  fZFrontCutShowers (pset.get<float>("ZFrontCutShowers",30.f)),
  fZBackCutShowers (pset.get<float>("ZBackCutShowers",50.f)),
  fCathodeXCut (pset.get<float>("CathodeXCut",1.5f))

  {
    art::ServiceHandle<art::TFileService> tfs;
    fParticleTree = tfs->make<TTree>("ParticleTree","Particle data TTree");

    fParticleTree->Branch("mc_trackID",&fMC_trackID);
    fParticleTree->Branch("mc_PDG",&fMC_PDG);
    fParticleTree->Branch("mc_x0",&fMC_x0);
    fParticleTree->Branch("mc_y0",&fMC_y0);
    fParticleTree->Branch("mc_z0",&fMC_z0);
    fParticleTree->Branch("mc_xEnd",&fMC_xEnd);
    fParticleTree->Branch("mc_yEnd",&fMC_yEnd);
    fParticleTree->Branch("mc_zEnd",&fMC_zEnd);
    fParticleTree->Branch("mc_pX0",&fMC_pX0);
    fParticleTree->Branch("mc_pY0",&fMC_pY0);
    fParticleTree->Branch("mc_pZ0",&fMC_pZ0);
    fParticleTree->Branch("mc_energy0",&fMC_energy0);
    fParticleTree->Branch("mc_momentum",&fMC_momentum);
    fParticleTree->Branch("mc_pXEnd",&fMC_pXEnd);
    fParticleTree->Branch("mc_pYEnd",&fMC_pYEnd);
    fParticleTree->Branch("mc_pZEnd",&fMC_pZEnd);
    fParticleTree->Branch("mc_energyEnd",&fMC_energyEnd);
    fParticleTree->Branch("mc_mass",&fMC_mass);
    fParticleTree->Branch("mc_theta_xy",&fMC_theta_xy);
    fParticleTree->Branch("mc_theta_yz",&fMC_theta_yz);
    fParticleTree->Branch("mc_theta_xz",&fMC_theta_xz);
    fParticleTree->Branch("mc_length",&fMC_length);
    fParticleTree->Branch("reco_isReconstructed",&fReco_isReconstructed);
    fParticleTree->Branch("reco_nTracks",&fReco_nTracks);
    fParticleTree->Branch("reco_nShowers",&fReco_nShowers);
    fParticleTree->Branch("reco_track_purity",&fReco_trackPurity);
    fParticleTree->Branch("reco_track_completeness",&fReco_trackCompleteness);
    fParticleTree->Branch("reco_shower_purity",&fReco_showerPurity);
    fParticleTree->Branch("reco_shower_completeness",&fReco_showerCompleteness);
    fParticleTree->Branch("reco_track_length",&fReco_trackLength);
    fParticleTree->Branch("reco_shower_dEdx",&fReco_showerdEdx);
    fParticleTree->Branch("reco_primary_isReconstructed",&fReco_primary_isReconstructed);
    fParticleTree->Branch("reco_primary_nTracks",&fReco_primary_nTracks);
    fParticleTree->Branch("reco_primary_nShowers",&fReco_primary_nShowers);
    fParticleTree->Branch("reco_primary_track_purity",&fReco_primary_trackPurity);
    fParticleTree->Branch("reco_primary_track_completeness",&fReco_primary_trackCompleteness);
    fParticleTree->Branch("reco_primary_shower_purity",&fReco_primary_showerPurity);
    fParticleTree->Branch("reco_primary_shower_completeness",&fReco_primary_showerCompleteness);
    fParticleTree->Branch("reco_primary_track_length",&fReco_primary_trackLength);
    fParticleTree->Branch("reco_primary_shower_dEdx",&fReco_primary_showerdEdx);


    fXOutEdge   = fTPCGeo.MaxX() - fXEdgeCut;
    fYEdge      = fTPCGeo.MaxY() - fYEdgeCut;
    fZFrontEdge = fTPCGeo.MinZ() + fZFrontCut;
    fZBackEdge  = fTPCGeo.MaxZ() - fZBackCut;

    fXOutEdgeShowers   = fTPCGeo.MaxX() - fXEdgeCutShowers;
    fYEdgeShowers      = fTPCGeo.MaxY() - fYEdgeCutShowers;
    fZFrontEdgeShowers = fTPCGeo.MinZ() + fZFrontCutShowers;
    fZBackEdgeShowers  = fTPCGeo.MaxZ() - fZBackCutShowers;

    fXCathodeEdge = 0 + fCathodeXCut;

  }

void RecoEff::ClearMaps()
{
  fNTracksMap.clear(); fTrackCompMap.clear(); fTrackPurMap.clear();
  fNShowersMap.clear(); fShowerCompMap.clear(); fShowerPurMap.clear();
  fTrackLengthMap.clear(); fShowerdEdxMap.clear();
  fHitsMap.clear();
}

void RecoEff::ResetData()
{
  fMC_trackID = def_int; fMC_PDG = def_int; 

  fMC_x0 = def_float; fMC_y0 = def_float; fMC_z0 = def_float; fMC_xEnd = def_float; fMC_yEnd = def_float; 
  fMC_zEnd = def_float; fMC_pX0 = def_float; fMC_pY0 = def_float; fMC_pZ0 = def_float; fMC_energy0 = def_float; 
  fMC_momentum = def_float; fMC_pXEnd = def_float; fMC_pYEnd = def_float; fMC_pZEnd = def_float; 
  fMC_energyEnd = def_float; fMC_mass = def_float; fMC_theta_xy = def_float; fMC_theta_yz = def_float;
  fMC_theta_xz = def_float; fMC_length = def_float; 

  fReco_nTracks = def_int; fReco_nShowers = def_int; 
  
  fReco_showerPurity = def_float; fReco_showerCompleteness = def_float; fReco_trackPurity = def_float; 
  fReco_trackCompleteness = def_float; fReco_trackLength = def_float; fReco_showerdEdx = def_float; 

  fReco_isReconstructed = false;

  fReco_primary_nTracks = def_int; fReco_primary_nShowers = def_int;

  fReco_primary_showerPurity = def_float; fReco_primary_showerCompleteness = def_float; fReco_primary_trackPurity = def_float;
  fReco_primary_trackCompleteness = def_float; fReco_primary_trackLength = def_float; fReco_primary_showerdEdx = def_float;

  fReco_primary_isReconstructed = false;
}

void RecoEff::SetupMaps(art::Event const &e)
{
  art::Handle<std::vector<recob::Hit> > handleHits;
  e.getByLabel(fHitsModuleLabel,handleHits);

  for(unsigned hit_i = 0; hit_i < handleHits->size(); ++hit_i)
    {
      const art::Ptr<recob::Hit> hit(handleHits,hit_i);
      fHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
    }
}

void RecoEff::ReconstructionProcessor(art::Event const &e)
{
  art::Handle<std::vector<recob::Slice> > handleSlices;
  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  art::Handle<std::vector<recob::Track> > handleTracks;
  art::Handle<std::vector<recob::Shower> > handleShowers;
  e.getByLabel(fSliceModuleLabel,handleSlices);
  e.getByLabel(fPFParticleModuleLabel,handlePFPs);
  e.getByLabel(fTrackModuleLabel,handleTracks);
  e.getByLabel(fShowerModuleLabel,handleShowers);

  art::FindManyP<recob::PFParticle> slicePFPAssn(handleSlices,e,fPFParticleModuleLabel);
  art::FindManyP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
  art::FindManyP<recob::Shower> pfpShowerAssn(handlePFPs,e,fShowerModuleLabel);
  art::FindManyP<recob::Hit> trackHitAssn(handleTracks,e,fTrackModuleLabel);
  art::FindManyP<recob::Hit> showerHitAssn(handleShowers,e,fShowerModuleLabel);

  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, handleSlices);

  for(auto const &slc : sliceVec)
    {
      const std::vector<art::Ptr<recob::PFParticle> > pfpVec = slicePFPAssn.at(slc.key());

      const art::Ptr<recob::PFParticle> primary = GetPrimaryPFP(pfpVec);

      for(auto const &pfp : pfpVec)
        {
          const bool primaryChild = (primary->PdgCode() == 12 || primary->PdgCode() != 14) && pfp->Parent() == pfp->Self();

          if(pfp->PdgCode() == 13)
            {
              const std::vector<art::Ptr<recob::Track> > tracksVec = pfpTrackAssn.at(pfp.key());
              if(tracksVec.size() != 1)
                continue;
              const art::Ptr<recob::Track> track = tracksVec[0];

              float x = track->Start().X(), y = track->Start().Y(), z = track->Start().Z();

              if(TMath::Abs(x) > fXOutEdge || TMath::Abs(x) < fXCathodeEdge ||
                 TMath::Abs(y) > fYEdge || z < fZFrontEdge || z > fZBackEdge)
                continue;

              std::vector<art::Ptr<recob::Hit> > trackHits = trackHitAssn.at(track.key());
              int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);
              float comp = Completeness(trackHits,trackID);
              float pur = Purity(trackHits,trackID);
              float length = track->Length();

              if(fNTracksMap[trackID] == 0)
                {
                  fTrackCompMap[trackID] = comp;
                  fTrackPurMap[trackID] = pur;
                  fTrackLengthMap[trackID] = length;
                }
              else if(comp > fTrackCompMap[trackID])
                {
                  fTrackCompMap[trackID] = comp;
                  fTrackPurMap[trackID] = pur;
                  fTrackLengthMap[trackID] = length;
                }

              fNTracksMap[trackID]++;

              if(primaryChild)
                {
                  if(fNPrimaryTracksMap[trackID] == 0)
                    {
                      fPrimaryTrackCompMap[trackID] = comp;
                      fPrimaryTrackPurMap[trackID] = pur;
                      fPrimaryTrackLengthMap[trackID] = length;
                    }
                  else if(comp > fPrimaryTrackCompMap[trackID])
                    {
                      fPrimaryTrackCompMap[trackID] = comp;
                      fPrimaryTrackPurMap[trackID] = pur;
                      fPrimaryTrackLengthMap[trackID] = length;
                    }

                  fNPrimaryTracksMap[trackID]++;
                }
            }
          else if(pfp->PdgCode() == 11)
            {
              const std::vector<art::Ptr<recob::Shower> > showersVec = pfpShowerAssn.at(pfp.key());
              if(showersVec.size() != 1)
                continue;
              const art::Ptr<recob::Shower> shower = showersVec[0];

              float x = shower->ShowerStart().X(), y = shower->ShowerStart().Y(), z = shower->ShowerStart().Z();

              if(TMath::Abs(x) > fXOutEdgeShowers || TMath::Abs(x) < fXCathodeEdge ||
                 TMath::Abs(y) > fYEdgeShowers || z < fZFrontEdgeShowers || z > fZBackEdgeShowers)
                continue;

              std::vector<art::Ptr<recob::Hit> > showerHits = showerHitAssn.at(shower.key());
              int trackID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,showerHits,true);
              float comp = Completeness(showerHits,trackID);
              float pur = Purity(showerHits,trackID);

              std::vector<double> dEdxVec = shower->dEdx();
              int best_plane = shower->best_plane();
              float dEdx = def_float;
              if(dEdxVec.size() != 0)
                dEdx = dEdxVec[best_plane];

              if(fNShowersMap[trackID] == 0)
                {
                  fShowerCompMap[trackID] = comp;
                  fShowerPurMap[trackID] = pur;
                  fShowerdEdxMap[trackID] = dEdx;
                }
              else if(comp > fShowerCompMap[trackID])
                {
                  fShowerCompMap[trackID] = comp;
                  fShowerPurMap[trackID] = pur;
                  fShowerdEdxMap[trackID] = dEdx;
                }

              fNShowersMap[trackID]++;

              if(primaryChild)
                {
                  if(fNPrimaryShowersMap[trackID] == 0)
                    {
                      fPrimaryShowerCompMap[trackID] = comp;
                      fPrimaryShowerPurMap[trackID] = pur;
                      fPrimaryShowerdEdxMap[trackID] = dEdx;
                    }
                  else if(comp > fPrimaryShowerCompMap[trackID])
                    {
                      fPrimaryShowerCompMap[trackID] = comp;
                      fPrimaryShowerPurMap[trackID] = pur;
                      fPrimaryShowerdEdxMap[trackID] = dEdx;
                    }

                  fNPrimaryShowersMap[trackID]++;
                }
            }
        }
    }
}

void RecoEff::TruthProcessor(art::Event const &e)
{
  art::Handle<std::vector<simb::MCTruth> > handleNeutrinos;
  e.getByLabel(fNuGenModuleLabel,handleNeutrinos);

  art::FindManyP<simb::MCParticle> nuParticleAssn(handleNeutrinos,e,fLArGeantModuleLabel);

  for(unsigned int nu_i = 0; nu_i < handleNeutrinos->size(); ++nu_i)
    {
      const art::Ptr<simb::MCTruth> truthNeutrino(handleNeutrinos,nu_i);
      if(truthNeutrino.isNull())
        continue;
      std::vector<art::Ptr<simb::MCParticle> > particles = nuParticleAssn.at(truthNeutrino.key());

      for(auto particle : particles)
        {
          if(particle->Process() != "primary" || particle->StatusCode() != 1)
            continue;

          if(particle->PdgCode() == 111)
            {
              for(int child_i = 0; child_i < particle->NumberDaughters(); ++child_i)
                {
                  const int childTrackID = particle->Daughter(child_i);
                  const simb::MCParticle* child = particleInv->TrackIdToParticle_P(childTrackID);
                  FillTrueParticle(child, true);
                }
            }

          FillTrueParticle(particle.get(), false);
        }
    }
}

void RecoEff::FillTrueParticle(const simb::MCParticle *particle, const bool pionChild)
{
  ResetData();

  fMC_trackID = particle->TrackId();

  if(pionChild)
    fMC_PDG = 111e6 + particle->PdgCode();
  else
    fMC_PDG = particle->PdgCode();

  fMC_x0 = particle->Vx();
  fMC_y0 = particle->Vy();
  fMC_z0 = particle->Vz();

  if(std::abs(fMC_PDG) == 11 || std::abs(fMC_PDG) == 22)
    {
      if(TMath::Abs(fMC_x0) > fXOutEdgeShowers || TMath::Abs(fMC_x0) < fXCathodeEdge ||
         TMath::Abs(fMC_y0) > fYEdgeShowers || fMC_z0 < fZFrontEdgeShowers || fMC_z0 > fZBackEdgeShowers)
        return;
    }
  else if(TMath::Abs(fMC_x0) > fXOutEdge || TMath::Abs(fMC_x0) < fXCathodeEdge ||
          TMath::Abs(fMC_y0) > fYEdge || fMC_z0 < fZFrontEdge || fMC_z0 > fZBackEdge)
    return;

  fMC_xEnd = particle->EndX();
  fMC_yEnd = particle->EndY();
  fMC_zEnd = particle->EndZ();
  fMC_pX0 = particle->Px();
  fMC_pY0 = particle->Py();
  fMC_pZ0 = particle->Pz();
  fMC_energy0 = particle->E();
  fMC_momentum = particle->P();
  fMC_pXEnd = particle->EndPx();
  fMC_pYEnd = particle->EndPy();
  fMC_pZEnd = particle->EndPz();
  fMC_energyEnd = particle->EndE();
  fMC_mass = particle->Mass();
  fMC_theta_xy = TMath::RadToDeg() * TMath::ATan(fMC_pX0/fMC_pY0);
  fMC_theta_yz = TMath::RadToDeg() * TMath::ATan(fMC_pY0/fMC_pZ0);
  fMC_theta_xz = TMath::RadToDeg() * TMath::ATan(fMC_pX0/fMC_pZ0);
  fMC_length = TrueTrackLength(particle);

  if(fNTracksMap[fMC_trackID] > 0 || fNShowersMap[fMC_trackID] > 0)
    fReco_isReconstructed = true;

  fReco_nTracks = fNTracksMap[fMC_trackID];
  fReco_nShowers = fNShowersMap[fMC_trackID];
  fReco_showerPurity = fShowerPurMap[fMC_trackID];
  fReco_showerCompleteness = fShowerCompMap[fMC_trackID];
  fReco_trackPurity = fTrackPurMap[fMC_trackID];
  fReco_trackCompleteness = fTrackCompMap[fMC_trackID];
  fReco_trackLength = fTrackLengthMap[fMC_trackID];
  fReco_showerdEdx = fShowerdEdxMap[fMC_trackID];

  if(fNPrimaryTracksMap[fMC_trackID] > 0 || fNPrimaryShowersMap[fMC_trackID] > 0)
    fReco_primary_isReconstructed = true;

  fReco_primary_nTracks = fNPrimaryTracksMap[fMC_trackID];
  fReco_primary_nShowers = fNPrimaryShowersMap[fMC_trackID];
  fReco_primary_showerPurity = fPrimaryShowerPurMap[fMC_trackID];
  fReco_primary_showerCompleteness = fPrimaryShowerCompMap[fMC_trackID];
  fReco_primary_trackPurity = fPrimaryTrackPurMap[fMC_trackID];
  fReco_primary_trackCompleteness = fPrimaryTrackCompMap[fMC_trackID];
  fReco_primary_trackLength = fPrimaryTrackLengthMap[fMC_trackID];
  fReco_primary_showerdEdx = fPrimaryShowerdEdxMap[fMC_trackID];

  fParticleTree->Fill();
}

art::Ptr<recob::PFParticle> RecoEff::GetPrimaryPFP(const std::vector<art::Ptr<recob::PFParticle>> &pfpVec)
{
  for(auto const& pfp : pfpVec)
    {
      if(pfp->IsPrimary())
        return pfp;
    }

  const art::Ptr<recob::PFParticle> nullReturn;
  return nullReturn;
}

art::Ptr<recob::PFParticle> RecoEff::GetPFP(art::Event const &e, long unsigned int const &id)
{
  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticleModuleLabel,handlePFPs);

  const art::Ptr<recob::PFParticle> nullReturn;

  for(unsigned int pfp_i = 0; pfp_i < handlePFPs->size(); ++pfp_i)
    {
      const art::Ptr<recob::PFParticle> pfp(handlePFPs,pfp_i);
      if(pfp->Self() == id)
        return pfp;
    }

  return nullReturn;
}

float RecoEff::Purity(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;

  return (objectHits.size() == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(objectHits.size());
}

float RecoEff::Completeness(std::vector< art::Ptr<recob::Hit> > const &objectHits, int const &trackID)
{
  std::map<int,int> objectHitsMap;

  for(unsigned int i = 0; i < objectHits.size(); ++i)
    objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)]++;

  return (fHitsMap[trackID] == 0) ? def_float : objectHitsMap[trackID]/static_cast<float>(fHitsMap[trackID]);
}

float RecoEff::TrueTrackLength(const simb::MCParticle *particle)
{
  float length = 0;
  unsigned int nTrajPoints = particle->NumberTrajectoryPoints();

  if(nTrajPoints < 2)
    return length;

  for(unsigned int point = 1; point < nTrajPoints; ++point)
    {
      TVector3 l = particle->Position(point).Vect();
      if(l.X() > fTPCGeo.MaxX() || l.X() < fTPCGeo.MinX() || l.Y() > fTPCGeo.MaxY() || l.Y() < fTPCGeo.MinY() ||
         l.Z() > fTPCGeo.MaxZ() || l.Z() < fTPCGeo.MinZ())
        break;

      TVector3 diff = particle->Position(point).Vect() - particle->Position(point-1).Vect();
      length += diff.Mag();
    }

  return length;
}


void RecoEff::analyze(art::Event const &e)
{
  clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  ClearMaps();
  SetupMaps(e);
  ReconstructionProcessor(e);
  TruthProcessor(e);
}

DEFINE_ART_MODULE(RecoEff)