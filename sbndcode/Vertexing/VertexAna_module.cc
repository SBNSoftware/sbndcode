////////////////////////////////////////////////////////////////////////
// Class:       VertexAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        VertexAna_module.cc
//
// Generated at Thu Feb  1 03:03:22 2024 by Henry Lay using cetskelgen
// from  version .
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

#include "TTree.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

constexpr int def_int     = std::numeric_limits<int>::min();
constexpr float def_float = std::numeric_limits<float>::lowest();

namespace sbnd {
  class VertexAna;
}

class sbnd::VertexAna : public art::EDAnalyzer {
public:
  explicit VertexAna(fhicl::ParameterSet const& p);

  VertexAna(VertexAna const&) = delete;
  VertexAna(VertexAna&&) = delete;
  VertexAna& operator=(VertexAna const&) = delete;
  VertexAna& operator=(VertexAna&&) = delete;

  void analyze(art::Event const& e) override;

  void ClearVars();

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

  art::InputTag fSliceLabel, fPFPLabel, fVertexLabel, fHitLabel;

  TTree *fSliceTree;

  int run, subrun, event, slice_id;
  float vtx_x, vtx_y, vtx_z, true_vtx_x, true_vtx_y, true_vtx_z,
    x_correction, comp, pur, dr;
  bool is_fv;
};


sbnd::VertexAna::VertexAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fSliceLabel (p.get<art::InputTag>("SliceLabel"))
  , fPFPLabel   (p.get<art::InputTag>("PFPLabel"))
  , fVertexLabel(p.get<art::InputTag>("VertexLabel"))
  , fHitLabel   (p.get<art::InputTag>("HitLabel"))
  {
    art::ServiceHandle<art::TFileService> fs;
    fSliceTree = fs->make<TTree>("SliceTree", "");

    fSliceTree->Branch("run", &run);
    fSliceTree->Branch("subrun", &subrun);
    fSliceTree->Branch("event", &event);
    fSliceTree->Branch("slice_id", &slice_id);
    fSliceTree->Branch("vtx_x", &vtx_x);
    fSliceTree->Branch("vtx_y", &vtx_y);
    fSliceTree->Branch("vtx_z", &vtx_z);
    fSliceTree->Branch("true_vtx_x", &true_vtx_x);
    fSliceTree->Branch("true_vtx_y", &true_vtx_y);
    fSliceTree->Branch("true_vtx_z", &true_vtx_z);
    fSliceTree->Branch("x_correction", &x_correction);
    fSliceTree->Branch("comp", &comp);
    fSliceTree->Branch("pur", &pur);
    fSliceTree->Branch("dr", &dr);
    fSliceTree->Branch("is_fv", &is_fv);
  }

void sbnd::VertexAna::analyze(art::Event const& e)
{
  run    = e.run();
  subrun = e.subRun();
  event  = e.event();

  slice_id = 0;

  const detinfo::DetectorClocksData clockData    = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
  const detinfo::DetectorPropertiesData propData = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);

  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitLabel, hitHandle);
  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  std::map<const art::Ptr<simb::MCTruth>, int> nuHitsMap;

  for(const art::Ptr<recob::Hit> &hit : hitVec)
    {
      const int trackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
      const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);

      if(nuHitsMap.find(mct) == nuHitsMap.end())
        nuHitsMap[mct] = 0;

      ++nuHitsMap[mct];
    }

  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceLabel, sliceHandle);
  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPFPLabel, pfpHandle);

  art::FindManyP<recob::Hit> slicesToHits(sliceHandle, e, fSliceLabel);
  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFPLabel);
  art::FindOneP<recob::Vertex> pfpsToVertices(pfpHandle, e, fVertexLabel);

  for(const art::Ptr<recob::Slice> &slice : sliceVec)
    {
      ClearVars();

      const std::vector<art::Ptr<recob::Hit>> sliceHits = slicesToHits.at(slice.key());

      std::map<int, int> objectHitMap;

      for(const art::Ptr<recob::Hit> &hit : sliceHits)
        {
          const int trackID = TruthMatchUtils::TrueParticleID(clockData, hit, true);
          if(objectHitMap.find(trackID) == objectHitMap.end())
            objectHitMap[trackID] = 0;

          ++objectHitMap[trackID];
        }

      std::map<const art::Ptr<simb::MCTruth>, int> mcTruthHitMap;

      for(auto const& [ trackID, nhits ] : objectHitMap)
        {
          const art::Ptr<simb::MCTruth> mct = trackID == def_int ? art::Ptr<simb::MCTruth>() : particleInv->TrackIdToMCTruth_P(trackID);
          if(mcTruthHitMap.find(mct) == mcTruthHitMap.end())
            mcTruthHitMap[mct] = 0;

          mcTruthHitMap[mct] += nhits;
        }

      int maxHits = def_int;
      art::Ptr<simb::MCTruth> bestMCT = art::Ptr<simb::MCTruth>();

      for(auto const& [ mct, nhits ] : mcTruthHitMap)
        {
          if(nhits > maxHits)
            {
              maxHits = nhits;
              bestMCT = mct;
            }
        }

      comp = nuHitsMap[bestMCT] == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(nuHitsMap[bestMCT]);
      pur  = sliceHits.size() == 0 ? def_float : mcTruthHitMap[bestMCT] / static_cast<float>(sliceHits.size());

      if(bestMCT.isNull() || bestMCT->Origin() != 1)
        {
          ++slice_id;
          continue;
        }

      const simb::MCParticle nu = bestMCT->GetNeutrino().Nu();
      true_vtx_x = nu.Vx();
      true_vtx_y = nu.Vy();
      true_vtx_z = nu.Vz();

      is_fv = std::abs(true_vtx_x) < 180 && std::abs(true_vtx_x) > 5
        && std::abs(true_vtx_y) < 180 && true_vtx_z > 10 && true_vtx_z < 500;

      x_correction = propData.ConvertTicksToX(clockData.TPCG4Time2Tick(nu.T()) - clockData.Time2Tick(clockData.BeamGateTime()), 0, 0, 0) - propData.ConvertTicksToX(0, 0, 0, 0);

      if(true_vtx_x < 0)
        x_correction *= -1.;

      const std::vector<art::Ptr<recob::PFParticle>> slicePFPs = slicesToPFPs.at(slice.key());

      for(const art::Ptr<recob::PFParticle> pfp : slicePFPs)
        {
          if(pfp->IsPrimary())
            {
              const art::Ptr<recob::Vertex> vtx = pfpsToVertices.at(pfp.key());
              vtx_x = vtx->position().X();
              vtx_y = vtx->position().Y();
              vtx_z = vtx->position().Z();
              break;
            }
        }
      
      dr = TMath::Sqrt(TMath::Power(vtx_x - true_vtx_x, 2)
                       + TMath::Power(vtx_y - true_vtx_y, 2)
                       + TMath::Power(vtx_z - true_vtx_z, 2));

      fSliceTree->Fill();
      ++slice_id;
    }
}

void sbnd::VertexAna::ClearVars()
{
  vtx_x = def_float;
  vtx_y = def_float;
  vtx_z = def_float;

  true_vtx_x = def_float;
  true_vtx_y = def_float;
  true_vtx_z = def_float;

  x_correction = def_float;

  comp = def_float;
  pur  = def_float;
  
  dr = def_float;

  is_fv = false;
}

DEFINE_ART_MODULE(sbnd::VertexAna)
