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
  void beginSubRun(const art::SubRun &sr);
  void endSubRun(const art::SubRun &sr); 
  void analyze(art::Event const& e) override;

  void MCParticles(art::Event const& e);
  void Slices(art::Event const& e);
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:


  TTree* SubRunTree;
  double pot;
  int spills, numGenEvents;
  unsigned int subRunNumber;
  unsigned int subRunRun;

  unsigned int eventID; // Event num
  unsigned int runID; // Run num  
  unsigned int subRunID; // Subrun num


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
  const std::string TruthLabel;
  double DLCurrent;
  const std::string spacePointLabel;
  const std::string clusterLabel;
  const std::string POTModuleLabel;
  double signal;

  TFile *outputFile = TFile::Open("NuEAnalyserOutput.root","RECREATE");

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
  TruthLabel(p.get<std::string>("TruthLabel")),
  DLCurrent(p.get<double>("DLCurrent")),
  spacePointLabel(p.get<std::string>("SpacePointLabel")),
  clusterLabel(p.get<std::string>("ClusterLabel")),
  POTModuleLabel(p.get<std::string>("POTLabel")),
  signal(p.get<double>("Signal"))
{

  art::ServiceHandle<art::TFileService> fs;  

  SubRunTree = fs->make<TTree>("SubRun","");   
  SubRunTree->Branch("pot", &pot);
  SubRunTree->Branch("spills", &spills);
  SubRunTree->Branch("numGenEvents", &numGenEvents);
  SubRunTree->Branch("DLCurrent", &DLCurrent);
  SubRunTree->Branch("signal", &signal);
  SubRunTree->Branch("run", &subRunRun);
  SubRunTree->Branch("subRun", &subRunNumber);
}


void sbnd::NuE::beginSubRun(const art::SubRun &sr){
    pot = 0.; spills = 0; numGenEvents = 0;

    art::Handle<sumdata::POTSummary> potHandle;
    sr.getByLabel(POTModuleLabel, potHandle);
    if(!potHandle.isValid()){
        //std::cout << "POT product " << POTModuleLabel << " not found..." << std::endl;
            return;
    } 
  
    pot = potHandle->totpot;
    spills = potHandle->totspills;
    subRunNumber = sr.subRun();
    subRunRun = sr.run();
    printf("Run = %i, SubRun = %i, POT = %f, Spills = %i\n", subRunRun, subRunNumber, pot, spills);

}

void sbnd::NuE::endSubRun(const art::SubRun &sr){
    //printf("POT = %f, Spills = %i, Num Events Generated = %i\n", pot, spills, numGenEvents);
    SubRunTree->Fill();
}

void sbnd::NuE::analyze(art::Event const& e)
{
    eventID = e.id().event();
    runID = e.id().run();
    subRunID = e.id().subRun();

    std::cout << "" << std::endl;
    std::cout << "========================================================================================================" << std::endl; 
    Slices(e);
    MCParticles(e);
}

void sbnd::NuE::MCParticles(art::Event const& e){
    art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
    if(e.getByLabel(TruthLabel, MCTruthHandle))
        art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

    std::cout << "_________ MCTruth _________" << std::endl;
    std::cout << "MCTruthVec.size = " << MCTruthVec.size() << std::endl;
    if(!MCTruthVec.empty()){
        int counter = 0;
        for(auto &MCTruth : MCTruthVec){
            counter++;
            std::cout << "MCTruth " << counter << ": Origin = " << MCTruth->Origin() << ", Number of Particles = " << MCTruth->NParticles() << std::endl;
            if(MCTruth->Origin() == simb::kBeamNeutrino){
                simb::MCNeutrino neutrino = MCTruth->GetNeutrino();
                simb::MCParticle neutrinoParticle = neutrino.Nu();
                std::cout << "Neutrino: Track ID = " << neutrinoParticle.TrackId() << ", CCNC = " << neutrino.CCNC() << ", Interaction Type = " << neutrino.InteractionType() << "Vertex = (" << neutrinoParticle.Vx() << ", " << neutrinoParticle.Vy() << ", " << neutrinoParticle.Vz() << "), Number of Daughters = " << neutrinoParticle.NumberDaughters() << std::endl; 
                int neutrinoTrackID = neutrinoParticle.TrackId();

                if(neutrino.InteractionType() == 1098){
                    // This is a Nu+E elastic scattering event
                    for(int i = 0; i < MCTruth->NParticles(); i++){
                        simb::MCParticle particle = MCTruth->GetParticle(i);
                        if(particle.TrackId() != neutrinoTrackID){
                            // This is not the neutrino
                            if(particle.StatusCode() == 1 && (particle.PdgCode() == 11 || particle.PdgCode() == -11)){
                                // StatusCode = 1: G4 is tracking this particle, it is outgoing
                                // Particle has to be an electron or positron
                                // This is the recoil electron!
                                
                                double P_mag = std::sqrt((particle.Px() * particle.Px()) + (particle.Py() * particle.Py()) + (particle.Pz() * particle.Pz()));
                                double energy = particle.E() * 1000; // Converts from GeV to MeV
                                double theta = std::acos(particle.Pz() / P_mag);
                                double EThetaSquared = (energy * theta * theta); 
                                
                                std::cout << "MCParticle: Track ID = " << particle.TrackId() << ", Vertex = (" << particle.Vx() << ", " << particle.Vy() << ", " << particle.Vz() << "), PDG = " << particle.PdgCode() << ", Mother = " << particle.Mother() << ", Status Code = " <<  particle.StatusCode() << ", Energy = " << energy << ", Theta = " << theta << ", EThetaSquared = " << EThetaSquared << std::endl;

                            }

                        }
                    }
                } else{
                    // This is a BNB or Cosmic event
                }
            }
        }
    }

    std::cout << "_________________________" << std::endl;

}

void sbnd::NuE::Slices(art::Event const& e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    std::cout << "_________ Slices _________" << std::endl;
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
                std::cout << "Slice " << sliceID << std::endl;
                if(sliceID == std::numeric_limits<int>::max()) continue;

                art::FindManyP<recob::Hit> sliceHitAssns(sliceVec, e, sliceLabel);
                const std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssns.at(slice.key()));

                // Gets the true particle ID of the truth particle who owns the most hits in the slice. True is for rollup.
                const int sliceTrueParticleID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sliceHits, true);
                std::cout << "The slice has most hits coming from true particle with ID = " << sliceTrueParticleID << std::endl;
                std::cout << "Number of hits in slice = " << sliceHits.size() << std::endl;

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
                    // Interaction Type = 1098 is Nu+E elastic scattering
                } else if(sliceMCTruth->Origin() == simb::kCosmicRay){
                    // This is a cosmic ray
                }

                // Loop through all the MCParticles in the MCTruth of the slice
                //for(int j = 0; j < sliceMCTruth->NParticles(); j++){
                    //simb::MCParticle sliceMCParticles = sliceMCTruth->GetParticle(j);
                    //std::cout << "MCParticle: PDG = " << sliceMCParticles.PdgCode() << ", Status Code = " << sliceMCParticles.StatusCode() << std::endl;
                    // Status Code = 1 means it is outgoing, this is tracked by g4
                //}

                double numTruthMatchedHits = 0;
                double numTruthMatchedHitsSlice = 0;
                double numHitsSlice = sliceHits.size();

                // Looping through the hits in the slice
                for(const art::Ptr<recob::Hit> &sliceHit : sliceHits){
                    // Gets the true particle ID of the truth particle which owns the hit. True is for rollup.
                    const int sliceHitTrueParticleID = TruthMatchUtils::TrueParticleID(clockData, sliceHit, true);
                    //std::cout << "The slice hit comes from true particle with ID = " << sliceHitTrueParticleID << std::endl;
                    
                    if(sliceHitTrueParticleID == std::numeric_limits<int>::min()) continue;
                    
                    // Gets the MCParticle that the hit comes from
                    //const simb::MCParticle* sliceHitMCParticle = particleInv->TrackIdToParticle_P(sliceHitTrueParticleID);
                    
                    // Gets the MCTruth that the hit comes from   
                    const art::Ptr<simb::MCTruth> sliceHitMCTruth = particleInv->TrackIdToMCTruth_P(sliceHitTrueParticleID);

                    // If the MCTruth from the hit is the same as the MCTruth of the slice, add to the counter
                    if(sliceMCTruth == sliceHitMCTruth) numTruthMatchedHitsSlice++;
                }

                // Looping through all the hits in the event
                for(const art::Ptr<recob::Hit> &hit : hitVec){
                    // Gets the true particle ID of the truth particle which owns the hit. True is for rollup.
                    const int hitTrueParticleID = TruthMatchUtils::TrueParticleID(clockData, hit, true);
                    //std::cout << "The hit comes from true particle with ID = " << hitTrueParticleID << std::endl;
                    
                    if(hitTrueParticleID == std::numeric_limits<int>::min()) continue;
                    
                    // Gets the MCParticle that the hit comes from 
                    //const simb::MCParticle* hitMCParticle = particleInv->TrackIdToParticle_P(hitTrueParticleID);
                    
                    // Gets the MCTruth that the hit comes from
                    const art::Ptr<simb::MCTruth> hitMCTruth = particleInv->TrackIdToMCTruth_P(hitTrueParticleID);

                    // If the MCTruth from the hit is the same as the MCTruth of the slice, add to the counter
                    if(sliceMCTruth == hitMCTruth) numTruthMatchedHits++;
                }
               
                std::cout << "Total number of hits = " << hitVec.size() << std::endl; 
                std::cout << "Number of hits truth matched to MCTruth = " << numTruthMatchedHits << std::endl;
                std::cout << "Number of hits truth matched to MCTruth in slice = " << numTruthMatchedHitsSlice << std::endl;
                std::cout << "Number of hits in slice = " << numHitsSlice << std::endl;

                double slicePurity = (numTruthMatchedHitsSlice/numHitsSlice);
                double sliceCompleteness = (numTruthMatchedHitsSlice/numTruthMatchedHits);

                std::cout << "\nSlice Completeness = " << sliceCompleteness << ", Slice Purity = " << slicePurity << std::endl;

                int sliceCategory = -999999;
                // Slice Categories: Cosmic = 0, Signal = 1, Fuzzy Signal = 2, BNB = 3, Fuzzy BNB = 4
                if(sliceMCTruth->Origin() == simb::kCosmicRay) sliceCategory = 0;
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() == 1098 && sliceCompleteness > 0.5 && slicePurity > 0.3) sliceCategory = 1; 
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() == 1098 && sliceCompleteness < 0.5 && sliceCompleteness > 0.1) sliceCategory = 2;
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() != 1098 && sliceCompleteness > 0.5 && slicePurity > 0.3) sliceCategory = 3;
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() != 1098 && sliceCompleteness < 0.5 && sliceCompleteness > 0.1) sliceCategory = 4;

                std::cout << "Slice Category = " << sliceCategory << std::endl;
            
            }

        }
    }

    std::cout << "_________________________" << std::endl;

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
