////////////////////////////////////////////////////////////////////////
// Class:       NuEWithSysts
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuEWithSysts_module.cc
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
#include "sbnobj/Common/Reco/MVAPID.h"

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
#include <iostream>
#include <Eigen/Dense>
#include <array>

// ROOT
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <TFile.h>

constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NuEWithSysts;
}


class sbnd::NuEWithSysts : public art::EDAnalyzer {
public:
  explicit NuEWithSysts(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuEWithSysts(NuEWithSysts const&) = delete;
  NuEWithSysts(NuEWithSysts&&) = delete;
  NuEWithSysts& operator=(NuEWithSysts const&) = delete;
  NuEWithSysts& operator=(NuEWithSysts&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void Slices(art::Event const& e);
  void clearVectors();

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  unsigned int eventID; // Event num
  unsigned int runID; // Run num  
  unsigned int subRunID; // Subrun num

  TTree* NuEWeightsTree;
  std::vector<double>   reco_sliceID;
  std::vector<double>   reco_sliceInteraction;
  std::vector<double>   reco_sliceTrueVX;
  std::vector<double>   reco_sliceTrueVY;
  std::vector<double>   reco_sliceTrueVZ;
  std::vector<double>   reco_sliceOrigin;
  std::vector<double>   reco_sliceTrueCCNC;
  std::vector<double>   reco_sliceTrueNeutrinoType;

  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_horncurrent;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_expskin;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_pioninexsec;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_pionqexsec;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_piontotxsec;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_nucleoninexsec;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_nucleonqexsec;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_nucleontotxsec;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_kplus;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_kmin;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_kzero;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_piplus;
  std::vector<std::vector<double>>   reco_sliceMCTruthFlux_weight_piminus;

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<geo::Geometry> theGeometry;
  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  const std::string PFParticleLabel;
  const std::string sliceLabel;
  const std::string sliceSCELabel;
  const std::string vertexLabel;
  const std::string nuGenModuleLabel;
  const std::string hitLabel;
  const std::string crumbsLabel;
  const std::string razzledLabel;
  const std::string trackLabel;
  const std::string showerLabel;
  const std::string MCTruthLabel;
  const std::string TruthLabel;
  const std::string spacePointLabel;
  double DLCurrent;
  const std::string clusterLabel;
  const std::string POTModuleLabel;
  double signal;
  const std::string fluxWeightLabel;
  const std::string genieWeightLabel;
  const std::string mcTruthFluxLabel;

  TFile *outputFile = TFile::Open("NuEWeightsAnalyserOutput.root","RECREATE");

};


sbnd::NuEWithSysts::NuEWithSysts(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  PFParticleLabel(p.get<std::string>("PFParticleLabel")),
  sliceLabel(p.get<std::string>("SliceLabel")),
  sliceSCELabel(p.get<std::string>("SliceSCELabel")),
  vertexLabel(p.get<std::string>("VertexLabel")),
  nuGenModuleLabel(p.get<std::string>("NuGenModuleLabel")),
  hitLabel(p.get<std::string>("HitLabel")),
  crumbsLabel(p.get<std::string>("CRUMBSLabel")),
  razzledLabel(p.get<std::string>("RazzledLabel")),
  trackLabel(p.get<std::string>("TrackLabel")),
  showerLabel(p.get<std::string>("ShowerLabel")),
  MCTruthLabel(p.get<std::string>("MCTruthLabel")),
  TruthLabel(p.get<std::string>("TruthLabel")),
  spacePointLabel(p.get<std::string>("SpacePointLabel")),
  DLCurrent(p.get<double>("DLCurrent")),
  clusterLabel(p.get<std::string>("ClusterLabel")),
  POTModuleLabel(p.get<std::string>("POTLabel")),
  signal(p.get<double>("Signal")),
  fluxWeightLabel(p.get<std::string>("FluxWeightLabel")),
  genieWeightLabel(p.get<std::string>("GenieWeightLabel")),
  mcTruthFluxLabel(p.get<std::string>("MCTruthFluxLabel"))
{

  art::ServiceHandle<art::TFileService> fs;  

  NuEWeightsTree = fs->make<TTree>("NuEWeights","");
  NuEWeightsTree->Branch("eventID", &eventID);
  NuEWeightsTree->Branch("runID", &runID);
  NuEWeightsTree->Branch("subRunID", &subRunID);
  NuEWeightsTree->Branch("DLCurrent", &DLCurrent);
  NuEWeightsTree->Branch("signal", &signal);
  
  NuEWeightsTree->Branch("reco_sliceID", &reco_sliceID);
  NuEWeightsTree->Branch("reco_sliceInteraction", &reco_sliceInteraction);
  NuEWeightsTree->Branch("reco_sliceTrueVX", &reco_sliceTrueVX);
  NuEWeightsTree->Branch("reco_sliceTrueVY", &reco_sliceTrueVY);
  NuEWeightsTree->Branch("reco_sliceTrueVZ", &reco_sliceTrueVZ);
  NuEWeightsTree->Branch("reco_sliceOrigin", &reco_sliceOrigin);
  NuEWeightsTree->Branch("reco_sliceTrueCCNC", &reco_sliceTrueCCNC);
  NuEWeightsTree->Branch("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);
  
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_horncurrent", &reco_sliceMCTruthFlux_weight_horncurrent);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_expskin", &reco_sliceMCTruthFlux_weight_expskin);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_pioninexsec", &reco_sliceMCTruthFlux_weight_pioninexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_pionqexsec", &reco_sliceMCTruthFlux_weight_pionqexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_piontotxsec", &reco_sliceMCTruthFlux_weight_piontotxsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_nucleoninexsec", &reco_sliceMCTruthFlux_weight_nucleoninexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_nucleonqexsec", &reco_sliceMCTruthFlux_weight_nucleonqexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_nucleontotxsec", &reco_sliceMCTruthFlux_weight_nucleontotxsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_kplus", &reco_sliceMCTruthFlux_weight_kplus);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_kmin", &reco_sliceMCTruthFlux_weight_kmin);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_kzero", &reco_sliceMCTruthFlux_weight_kzero);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_piplus", &reco_sliceMCTruthFlux_weight_piplus);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_piminus", &reco_sliceMCTruthFlux_weight_piminus);
  
}

void sbnd::NuEWithSysts::analyze(art::Event const& e)
{
    clearVectors();

    eventID = e.id().event();
    runID = e.id().run();
    subRunID = e.id().subRun();

    std::cout << "" << std::endl;
    std::cout << "========================================================================================================" << std::endl; 
    std::cout << "Run: " << runID << ", Subrun: " << subRunID << ", Event: " << eventID << ", DL/Current: " << DLCurrent << std::endl;
    Slices(e);
    
    NuEWeightsTree->Fill();
}

void sbnd::NuEWithSysts::Slices(art::Event const& e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    std::cout << "_________ Slices _________" << std::endl;

    int counter = 0;

    art::Handle<std::vector<recob::Slice>>  sliceHandle;
    std::vector<art::Ptr<recob::Slice>>     sliceVec;
    if(e.getByLabel(sliceLabel, sliceHandle))
        art::fill_ptr_vector(sliceVec, sliceHandle);

    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>>    mcTruthVec;
    if(e.getByLabel(mcTruthFluxLabel, mcTruthHandle))
        art::fill_ptr_vector(mcTruthVec, mcTruthHandle);

    std::cout << "------------ Number of slices = " << sliceVec.size() << std::endl;

    if(!sliceVec.empty()){
        art::Handle<std::vector<recob::Hit>>    hitHandle;
        std::vector<art::Ptr<recob::Hit>>       hitVec;

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
                //const simb::MCParticle* sliceMCParticle = particleInv->TrackIdToParticle_P(sliceTrueParticleID);
                
                // Get the MCTruth associated with the majority of the hits in the slice
                const art::Ptr<simb::MCTruth> sliceMCTruth = particleInv->TrackIdToMCTruth_P(sliceTrueParticleID);

                // Get the MCNeutrino associated with the majority of the hits in the slice
                const simb::MCNeutrino sliceMCNeutrino = sliceMCTruth->GetNeutrino();
                const simb::MCParticle sliceMCNeutrinoParticle = sliceMCNeutrino.Nu();

                // Looking at flux systematics associated with MCTruth of slice
                art::FindManyP<std::map<std::string,std::vector<float>>> mcTruthFluxWeightMapAssns(mcTruthVec, e, fluxWeightLabel);
                const std::vector<art::Ptr<std::map<std::string,std::vector<float>>>> mcTruthFluxWeightMaps(mcTruthFluxWeightMapAssns.at(sliceMCTruth.key()));
        
                std::vector<double> reco_sliceMCTruthFlux_weight_horncurrent_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_expskin_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_pioninexsec_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_pionqexsec_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_piontotxsec_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_nucleoninexsec_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_nucleonqexsec_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_nucleontotxsec_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_kplus_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_kmin_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_kzero_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_piplus_vector;
                std::vector<double> reco_sliceMCTruthFlux_weight_piminus_vector;

                std::cout << "============================" << std::endl;
                std::cout << "Size of mcTruthFluxWeightMaps vector = " << mcTruthFluxWeightMaps.size() << std::endl;
                const art::Ptr<std::map<std::string, std::vector<float>>> &weightFluxMapPtr = mcTruthFluxWeightMaps.at(0);
                const std::map<std::string, std::vector<float>> &weightFluxMap = *weightFluxMapPtr;

                std::cout << "Number of map entries = " << weightFluxMap.size() << std::endl;

                std::unordered_map<std::string, std::vector<double>*> targetMap = {
                    {"horncurrent_Flux",      &reco_sliceMCTruthFlux_weight_horncurrent_vector},
                    {"expskin_Flux",          &reco_sliceMCTruthFlux_weight_expskin_vector},
                    {"pioninexsec_Flux",      &reco_sliceMCTruthFlux_weight_pioninexsec_vector},
                    {"pionqexsec_Flux",       &reco_sliceMCTruthFlux_weight_pionqexsec_vector},
                    {"piontotxsec_Flux",      &reco_sliceMCTruthFlux_weight_piontotxsec_vector},
                    {"nucleoninexsec_Flux",   &reco_sliceMCTruthFlux_weight_nucleoninexsec_vector},
                    {"nucleonqexsec_Flux",    &reco_sliceMCTruthFlux_weight_nucleonqexsec_vector},
                    {"nucleontotxsec_Flux",   &reco_sliceMCTruthFlux_weight_nucleontotxsec_vector},
                    {"kplus_Flux",            &reco_sliceMCTruthFlux_weight_kplus_vector},
                    {"kminus_Flux",           &reco_sliceMCTruthFlux_weight_kmin_vector},
                    {"kzero_Flux",            &reco_sliceMCTruthFlux_weight_kzero_vector},
                    {"piplus_Flux",           &reco_sliceMCTruthFlux_weight_piplus_vector},
                    {"piminus_Flux",          &reco_sliceMCTruthFlux_weight_piminus_vector} 
                };

                for(const auto &[name, values] : weightFluxMap){
                    std::cout << "Map key = " << name << std::endl;
                    std::cout << "Number of universes? = " << values.size() << std::endl;

                    for(float v : values){
                        std::cout << "  Weight = " << v << std::endl;
                    }

                    auto mapFill = targetMap.find(name);
                    if(mapFill == targetMap.end()) continue;

                    auto &target = *(mapFill->second);

                    target.assign(values.begin(), values.end());
                
                }
                
                std::cout << "============================" << std::endl;

                // Looking at GENIE systematics associated with MCTruth of slice
                /*
                art::FindManyP<std::map<std::string, std::vector<float>>> mcTruthGENIEWeightMapAssns(mcTruthVec, e, genieWeightLabel);
                const std::vector<art::Ptr<std::map<std::string, std::vector<float>>>> mcTruthGENIEWeightMaps(mcTruthGENIEWeightMapAssns.at(sliceMCTruth.key()));
                
                std::cout << "============================" << std::endl;
                std::cout << "Size of mcTruthGENIEWeightMaps vector = " << mcTruthGENIEWeightMaps.size() << std::endl;
                
                for(size_t i = 0; i < mcTruthGENIEWeightMaps.size(); ++i){
                    std::cout << "Vector " << i << ":" << std::endl;
                    const art::Ptr<std::map<std::string, std::vector<float>>> &weightGENIEMapPtr = mcTruthGENIEWeightMaps.at(i);
                    const std::map<std::string, std::vector<float>> &weightGENIEMap = *weightGENIEMapPtr;

                    std::cout << "  Number of map entries = " << weightGENIEMap.size() << std::endl;

                    for(const auto &[name, values] : weightGENIEMap){
                        std::cout << "  Map key = " << name << std::endl;
                        std::cout << "  Number of universes? = " << values.size() << std::endl;
                    }
                }

                std::cout << "============================" << std::endl;
                */

                counter++;
                reco_sliceID.push_back(sliceID);
                reco_sliceMCTruthFlux_weight_horncurrent.push_back(reco_sliceMCTruthFlux_weight_horncurrent_vector);
                reco_sliceMCTruthFlux_weight_expskin.push_back(reco_sliceMCTruthFlux_weight_expskin_vector);
                reco_sliceMCTruthFlux_weight_pioninexsec.push_back(reco_sliceMCTruthFlux_weight_pioninexsec_vector);
                reco_sliceMCTruthFlux_weight_pionqexsec.push_back(reco_sliceMCTruthFlux_weight_pionqexsec_vector);
                reco_sliceMCTruthFlux_weight_piontotxsec.push_back(reco_sliceMCTruthFlux_weight_piontotxsec_vector);
                reco_sliceMCTruthFlux_weight_nucleoninexsec.push_back(reco_sliceMCTruthFlux_weight_nucleoninexsec_vector);
                reco_sliceMCTruthFlux_weight_nucleonqexsec.push_back(reco_sliceMCTruthFlux_weight_nucleonqexsec_vector);
                reco_sliceMCTruthFlux_weight_nucleontotxsec.push_back(reco_sliceMCTruthFlux_weight_nucleontotxsec_vector);
                reco_sliceMCTruthFlux_weight_kplus.push_back(reco_sliceMCTruthFlux_weight_kplus_vector);
                reco_sliceMCTruthFlux_weight_kmin.push_back(reco_sliceMCTruthFlux_weight_kmin_vector);
                reco_sliceMCTruthFlux_weight_kzero.push_back(reco_sliceMCTruthFlux_weight_kzero_vector);
                reco_sliceMCTruthFlux_weight_piplus.push_back(reco_sliceMCTruthFlux_weight_piplus_vector);
                reco_sliceMCTruthFlux_weight_piminus.push_back(reco_sliceMCTruthFlux_weight_piminus_vector);

                reco_sliceMCTruthFlux_weight_horncurrent_vector.clear();
                reco_sliceMCTruthFlux_weight_expskin_vector.clear();
                reco_sliceMCTruthFlux_weight_pioninexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_pionqexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_piontotxsec_vector.clear();
                reco_sliceMCTruthFlux_weight_nucleoninexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_nucleonqexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_nucleontotxsec_vector.clear();
                reco_sliceMCTruthFlux_weight_kplus_vector.clear();
                reco_sliceMCTruthFlux_weight_kmin_vector.clear();
                reco_sliceMCTruthFlux_weight_kzero_vector.clear();
                reco_sliceMCTruthFlux_weight_piplus_vector.clear();
                reco_sliceMCTruthFlux_weight_piminus_vector.clear();

                if(sliceMCTruth->Origin() == simb::kBeamNeutrino){
                    reco_sliceInteraction.push_back(sliceMCNeutrino.InteractionType());
                    reco_sliceTrueVX.push_back(sliceMCNeutrinoParticle.Vx());
                    reco_sliceTrueVY.push_back(sliceMCNeutrinoParticle.Vy());
                    reco_sliceTrueVZ.push_back(sliceMCNeutrinoParticle.Vz());
                    reco_sliceTrueCCNC.push_back(sliceMCNeutrino.CCNC());
                    reco_sliceTrueNeutrinoType.push_back(sliceMCNeutrinoParticle.PdgCode());
                    if(sliceMCNeutrino.InteractionType() == 1098){
                        reco_sliceOrigin.push_back(1);
                    } else{
                        reco_sliceOrigin.push_back(3);
                    }
                } else if(sliceMCTruth->Origin() == simb::kCosmicRay || sliceMCTruth->Origin() == 0){
                    reco_sliceInteraction.push_back(-100);        
                    reco_sliceTrueVX.push_back(-999999);
                    reco_sliceTrueVY.push_back(-999999);
                    reco_sliceTrueVZ.push_back(-999999);
                    reco_sliceOrigin.push_back(0);
                    reco_sliceTrueCCNC.push_back(-999999);
                    reco_sliceTrueNeutrinoType.push_back(-999999);
                } else {
                    reco_sliceInteraction.push_back(-999999);        
                    reco_sliceTrueVX.push_back(-999999);
                    reco_sliceTrueVY.push_back(-999999);
                    reco_sliceTrueVZ.push_back(-999999);
                    reco_sliceOrigin.push_back(-999999);
                    reco_sliceTrueCCNC.push_back(-999999);
                    reco_sliceTrueNeutrinoType.push_back(-999999);
                }
            }
        }
    }

    if(counter == 0){
        reco_sliceID.push_back(-999999);        
        reco_sliceInteraction.push_back(-999999);       
        reco_sliceTrueVX.push_back(-999999); 
        reco_sliceTrueVY.push_back(-999999); 
        reco_sliceTrueVZ.push_back(-999999); 
        reco_sliceOrigin.push_back(-999999);
        reco_sliceTrueCCNC.push_back(-999999);
        reco_sliceTrueNeutrinoType.push_back(-999999);
       
        std::vector<double> reco_sliceMCTruthFlux_weight_horncurrent_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_expskin_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_pioninexsec_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_pionqexsec_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_piontotxsec_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_nucleoninexsec_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_nucleonqexsec_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_nucleontotxsec_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_kplus_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_kmin_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_kzero_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_piplus_vector;
        std::vector<double> reco_sliceMCTruthFlux_weight_piminus_vector;

        reco_sliceMCTruthFlux_weight_horncurrent_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_expskin_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_pioninexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_pionqexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_piontotxsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_nucleoninexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_nucleonqexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_nucleontotxsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_kplus_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_kmin_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_kzero_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_piplus_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_piminus_vector.push_back(-999999);

        reco_sliceMCTruthFlux_weight_horncurrent.push_back(reco_sliceMCTruthFlux_weight_horncurrent_vector);
        reco_sliceMCTruthFlux_weight_expskin.push_back(reco_sliceMCTruthFlux_weight_expskin_vector);
        reco_sliceMCTruthFlux_weight_pioninexsec.push_back(reco_sliceMCTruthFlux_weight_pioninexsec_vector);
        reco_sliceMCTruthFlux_weight_pionqexsec.push_back(reco_sliceMCTruthFlux_weight_pionqexsec_vector);
        reco_sliceMCTruthFlux_weight_piontotxsec.push_back(reco_sliceMCTruthFlux_weight_piontotxsec_vector);
        reco_sliceMCTruthFlux_weight_nucleoninexsec.push_back(reco_sliceMCTruthFlux_weight_nucleoninexsec_vector);
        reco_sliceMCTruthFlux_weight_nucleonqexsec.push_back(reco_sliceMCTruthFlux_weight_nucleonqexsec_vector);
        reco_sliceMCTruthFlux_weight_nucleontotxsec.push_back(reco_sliceMCTruthFlux_weight_nucleontotxsec_vector);
        reco_sliceMCTruthFlux_weight_kplus.push_back(reco_sliceMCTruthFlux_weight_kplus_vector);
        reco_sliceMCTruthFlux_weight_kmin.push_back(reco_sliceMCTruthFlux_weight_kmin_vector);
        reco_sliceMCTruthFlux_weight_kzero.push_back(reco_sliceMCTruthFlux_weight_kzero_vector);
        reco_sliceMCTruthFlux_weight_piplus.push_back(reco_sliceMCTruthFlux_weight_piplus_vector);
        reco_sliceMCTruthFlux_weight_piminus.push_back(reco_sliceMCTruthFlux_weight_piminus_vector);
    }

    std::cout << "_________________________" << std::endl;
}

void sbnd::NuEWithSysts::clearVectors(){

    reco_sliceID.clear();
    reco_sliceInteraction.clear();
    reco_sliceTrueVX.clear();
    reco_sliceTrueVY.clear();
    reco_sliceTrueVZ.clear();
    reco_sliceOrigin.clear();
    reco_sliceTrueCCNC.clear();
    reco_sliceTrueNeutrinoType.clear();
   
    reco_sliceMCTruthFlux_weight_horncurrent.clear();
    reco_sliceMCTruthFlux_weight_expskin.clear();
    reco_sliceMCTruthFlux_weight_pioninexsec.clear();
    reco_sliceMCTruthFlux_weight_pionqexsec.clear();
    reco_sliceMCTruthFlux_weight_piontotxsec.clear();
    reco_sliceMCTruthFlux_weight_nucleoninexsec.clear();
    reco_sliceMCTruthFlux_weight_nucleonqexsec.clear();
    reco_sliceMCTruthFlux_weight_nucleontotxsec.clear();
    reco_sliceMCTruthFlux_weight_kplus.clear();
    reco_sliceMCTruthFlux_weight_kmin.clear();
    reco_sliceMCTruthFlux_weight_kzero.clear();
    reco_sliceMCTruthFlux_weight_piplus.clear();
    reco_sliceMCTruthFlux_weight_piminus.clear();
}

void sbnd::NuEWithSysts::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::NuEWithSysts::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::NuEWithSysts)
