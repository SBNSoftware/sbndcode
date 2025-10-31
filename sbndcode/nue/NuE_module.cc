////////////////////////////////////////////////////////////////////////
// Class:       NuE
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuE_module.cc
//
// Generated at Fri Oct 31 08:23:53 2025 by Rachel Coackley using cetskelgen
// from cetlib version 3.18.02.
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

// Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// Reco Base
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
// #include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcore/Geometry/WireReadout.h"

#include "sbnobj/Common/Reco/CRUMBSResult.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

// C++
#include <string>
#include <vector>
#include <map>
#include <cmath>

// ROOT
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <TFile.h>


namespace sbnd {
  class NuE;
}


class sbnd::NuE : public art::EDAnalyzer {
public:
  explicit NuE(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuE(NuE const&) = delete;
  NuE(NuE&&) = delete;
  NuE& operator=(NuE const&) = delete;
  NuE& operator=(NuE&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void Slices(art::Event const& e);
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

  const std::string PFParticleLabel;
  const std::string sliceLabel;
  const std::string sliceSCELabel;
  const std::string vertexLabel;
  const std::string nuGenModuleLabel;
  const std::string hitLabel;
  const std::string crumbsLabel;
  const std::string trackLabel;
  const std::string showerLabel;
  const std::string MCTruthLabel;
  double DLCurrent;
  const std::string spacePointLabel;
  const std::string clusterLabel;
  const std::string POTModuleLabel;
  double signal;

};


sbnd::NuE::NuE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  PFParticleLabel(p.get<std::string>("PFParticleLabel")),
  sliceLabel(p.get<std::string>("SliceLabel")),
  sliceSCELabel(p.get<std::string>("SliceSCELabel")),
  vertexLabel(p.get<std::string>("VertexLabel")),
  nuGenModuleLabel(p.get<std::string>("NuGenModuleLabel")),
  hitLabel(p.get<std::string>("HitLabel")),
  crumbsLabel(p.get<std::string>("CRUMBSLabel")),
  trackLabel(p.get<std::string>("TrackLabel")),
  showerLabel(p.get<std::string>("ShowerLabel")),
  MCTruthLabel(p.get<std::string>("MCTruthLabel")),
  DLCurrent(p.get<double>("DLCurrent")),
  spacePointLabel(p.get<std::string>("SpacePointLabel")),
  clusterLabel(p.get<std::string>("ClusterLabel")),
  POTModuleLabel(p.get<std::string>("POTLabel")),
  signal(p.get<double>("Signal"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::NuE::analyze(art::Event const& e)
{
  // Implementation of required member function here.
    Slices(e);
}

void sbnd::NuE::Slices(art::Event const& e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    art::Handle<std::vector<recob::Slice>>  sliceHandle;
    std::vector<art::Ptr<recob::Slice>>     sliceVec;
    if(e.getByLabel(sliceLabel, sliceHandle))
        art::fill_ptr_vector(sliceVec, sliceHandle);

    if(!sliceVec.empty()){
        art::Handle<std::vector<recob::Hit>> hitHandle;
        std::vector<art::Ptr<recob::Hit>> hitVec;

        if(e.getByLabel(hitLabel, hitHandle))
            art::fill_ptr_vector(hitVec, hitHandle);

        if(!hitVec.empty()){
            int sliceID(std::numeric_limits<int>::max());

            for(const art::Ptr<recob::Slice> &slice : sliceVec){
                sliceID = slice->ID();
                std::cout << "_______________________________________________________" << std::endl;
                std::cout << "Slice " << sliceID << std::endl;
                if(sliceID == std::numeric_limits<int>::max()) continue;

                art::FindManyP<recob::Hit> sliceHitAssns(sliceVec, e, sliceLabel);
                const std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssns.at(slice.key()));

                // Gets the true particle ID of the truth particle who owns the most hits in the slice. True is for rollup.
                const int sliceTrueParticleID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sliceHits, true);
                std::cout << "The slice has most hits coming from true particle with ID = " << sliceTrueParticleID << std::endl;
            
                if(sliceTrueParticleID == std::numeric_limits<int>::min()) continue;

                // Get the MCParticle associated with the majority of the hits in the slice
                const simb::MCParticle* sliceMCParticle = particleInv->TrackIdToParticle_P(sliceTrueParticleID);
                
                // Get the MCTruth associated with the majority of the hits in the slice
                const art::Ptr<simb::MCTruth> sliceMCTruth = particleInv->TrackIdToMCTruth_P(sliceTrueParticleID);

                // Get the MCNeutrino associated with the majority of the hits in the slice
                const simb::MCNeutrino sliceMCNeutrino = sliceMCTruth->GetNeutrino();
                const simb::MCParticle sliceMCNeutrinoParticle = sliceMCNeutrino.Nu();

                std::cout << "MCParticle: PDG = " << sliceMCParticle->PdgCode() << std::endl;
                std::cout << "MCTruth: Origin = " << sliceMCTruth->Origin() << ", Num Particles = " << sliceMCTruth->NParticles() << std::endl;
                
                if(sliceMCTruth->Origin() == simb::kBeamNeutrino){
                    // This is a beam neutrino
                    std::cout << "MCNeutrino: Interaction Type = " << sliceMCNeutrino.InteractionType() << ", CCNC = " << sliceMCNeutrino.CCNC() << ", PDG = " << sliceMCNeutrinoParticle.PdgCode() << ", Vertex = (" << sliceMCNeutrinoParticle.Vx() << ", " << sliceMCNeutrinoParticle.Vy() << ", " << sliceMCNeutrinoParticle.Vy() << ")" << std::endl;
                } else if(sliceMCTruth->Origin() == simb::kCosmicRay){
                    // This is a cosmic ray
                }

                int numTruthMatchedHits = 0;
                int numTruthMatchedHitsSlice = 0;
                int numHitsSlice = ;

                for(const art::Ptr<recob::Hit> &hit : hitVec){
                    const int hitTrueParticleID = TruthMatchUtils::TrueParticleID(clockData, hit, true);
                    std::cout << "The hit comes from true particle with ID = " << hitTrueParticleID << std::endl;

                }

                

            }

        }
    }


}

void sbnd::NuE::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::NuE::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::NuE)
