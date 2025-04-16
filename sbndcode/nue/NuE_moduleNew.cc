////////////////////////////////////////////////////////////////////////
// Class:       NuE
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuE_module.cc
//
// Generated at Fri Feb 14 10:01:33 2025 by Rachel Coackley using cetskelgen
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

// C++
#include <string>
#include <vector>
#include <map>
#include <cmath>

// ROOT
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <TFile.h>

constexpr double def_double = -std::numeric_limits<double>::max();

double dist_between(double x0, double y0, double z0, double x1, double y1, double z1) {
  return pow((pow(x0-x1, 2.0) + pow(y0-y1, 2.0) + pow(z0-z1, 2.0)), 0.5);
}

namespace sbnd {
  class NuE;
}


class sbnd::NuE : public art::EDAnalyzer {
public:
    explicit NuE(fhicl::ParameterSet const& p);

    NuE(NuE const&) = delete;
    NuE(NuE&&) = delete;
    NuE& operator=(NuE const&) = delete;
    NuE& operator=(NuE&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
    double Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID);
    double Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID);
    void ClearMaps(const art::Event &e);
    void SetupMaps(const art::Event &e);
    void showerEnergy(const art::Event &e);
    void recoNeutrino(const art::Event &e);
    void trueNeutrino(const art::Event &e);
    void slices(const art::Event &e);
    void PFPs(const art::Event &e);
    void MCParticles(const art::Event &e);
    void hits(const art::Event &e);
    void clearVectors();

private:

    art::ServiceHandle<cheat::BackTrackerService>       backTracker;

    TTree* NuETree;

    // Tree Variables
    unsigned int eventID;                    // Event num
    unsigned int runID;                      // Run num  
    unsigned int subRunID;                   // Subrun num

    std::vector<double>     truth_neutrinoVX;
    std::vector<double>     truth_neutrinoVY;
    std::vector<double>     truth_neutrinoVZ;
    std::vector<double>     truth_CCNC;                 // MCTruth CC or NC, CC = 0, NC = 1
    std::vector<double>     truth_neutrinoType;         // MCTruth neutrino flavour (PDG code)
    std::vector<double>     truth_chargedLepton;        // MCTruth charged lepton flavour (PDG code)
    std::vector<double>     truth_neutrinoTPCID;        // MCTruth neutrino in TPC ID
    std::vector<double>     truth_neutrinoTPCValid;     // MCTruth neutrino in TPC?
    
    std::vector<double>     truth_particlePDG;
    std::vector<double>     truth_particleVX;
    std::vector<double>     truth_particleVY;
    std::vector<double>     truth_particleVZ;
    std::vector<double>     truth_particlePX;
    std::vector<double>     truth_particlePY;
    std::vector<double>     truth_particlePZ;
    std::vector<double>     truth_particleEnergy;
    std::vector<double>     truth_particleAngle;
    std::vector<double>     truth_particleETheta2;
    std::vector<double>     truth_particleDirectionX;
    std::vector<double>     truth_particleDirectionY;
    std::vector<double>     truth_particleDirectionZ;

    std::vector<double>     reco_neutrinoPDG;
    std::vector<double>     reco_neutrinoIsPrimary;
    std::vector<double>     reco_neutrinoVX;
    std::vector<double>     reco_neutrinoVY;
    std::vector<double>     reco_neutrinoVZ;
    std::vector<double>     reco_neutrinoSliceID;

    std::vector<double>     reco_particlePDG;
    std::vector<double>     reco_particleIsPrimary;
    std::vector<double>     reco_particleVX;
    std::vector<double>     reco_particleVY;
    std::vector<double>     reco_particleVZ;
    std::vector<double>     reco_particleDirectionX;
    std::vector<double>     reco_particleDirectionY;
    std::vector<double>     reco_particleDirectionZ;
    std::vector<double>     reco_particleSliceID;
    std::vector<double>     reco_particleBestPlaneEnergy;
    std::vector<double>     reco_particleTheta;
    std::vector<double>     reco_particleTrackScore;
    std::vector<double>     reco_particleCompleteness;
    std::vector<double>     reco_particlePurity;

    std::vector<double>     reco_sliceID;
    std::vector<double>     reco_sliceCompleteness;
    std::vector<double>     reco_slicePurity;
    std::vector<double>     reco_sliceScore;

    std::vector<double>     reco_hitPlane;
    std::vector<double>     reco_hitX;
    std::vector<double>     reco_hitUVZ;

    std::map<int,int> fHitsMap;

    // Geometry information
    detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    detinfo::DetectorPropertiesData propData = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
    art::ServiceHandle<geo::Geometry> theGeometry;
    geo::GeometryCore const* geom = art::ServiceHandle<geo::Geometry>()->provider();
 
    // Label strings
    const std::string PFParticleLabel;
    const std::string sliceLabel;
    const std::string vertexLabel;
    const std::string nuGenModuleLabel;
    const std::string hitLabel;
    const std::string crumbsLabel;
    const std::string trackLabel;
    const std::string showerLabel;
    const std::string MCTruthLabel;

    double DLCurrent;

    // Output file
    TFile *outputFile = TFile::Open("NuEAnalyserOutput.root","RECREATE");
};


sbnd::NuE::NuE(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  PFParticleLabel(p.get<std::string>("PFParticleLabel")),
  sliceLabel(p.get<std::string>("SliceLabel")),
  vertexLabel(p.get<std::string>("VertexLabel")),
  nuGenModuleLabel(p.get<std::string>("NuGenModuleLabel")),
  hitLabel(p.get<std::string>("HitLabel")),
  crumbsLabel(p.get<std::string>("CRUMBSLabel")),
  trackLabel(p.get<std::string>("TrackLabel")),
  showerLabel(p.get<std::string>("ShowerLabel")),
  MCTruthLabel(p.get<std::string>("MCTruthLabel")),
  DLCurrent(p.get<double>("DLCurrent"))
{
    art::ServiceHandle<art::TFileService> fs;
    NuETree = fs->make<TTree>("NuE","");
    
    NuETree->Branch("eventID", &eventID);
    NuETree->Branch("runID", &runID);
    NuETree->Branch("subRunID", &subRunID);
    NuETree->Branch("DLCurrent", &DLCurrent);

    NuETree->Branch("truth_neutrinoVX", "std::vector<double>", &truth_neutrinoVX);
    NuETree->Branch("truth_neutrinoVY", "std::vector<double>", &truth_neutrinoVY);
    NuETree->Branch("truth_neutrinoVZ", "std::vector<double>", &truth_neutrinoVZ);
    NuETree->Branch("truth_CCNC", "std::vector<double>", &truth_CCNC);
    NuETree->Branch("truth_neutrinoType", "std::vector<double>", &truth_neutrinoType);
    NuETree->Branch("truth_chargedLepton", "std::vector<double>", &truth_chargedLepton);
    NuETree->Branch("truth_neutrinoTPCID", "std::vector<double>", &truth_neutrinoTPCID);
    NuETree->Branch("truth_neutrinoTPCValid", "std::vector<double>", &truth_neutrinoTPCValid);

    NuETree->Branch("truth_particlePDG", "std::vector<double>", &truth_particlePDG);
    NuETree->Branch("truth_particleVX", "std::vector<double>", &truth_particleVX);
    NuETree->Branch("truth_particleVY", "std::vector<double>", &truth_particleVY);
    NuETree->Branch("truth_particleVZ", "std::vector<double>", &truth_particleVZ);
    NuETree->Branch("truth_particlePX", "std::vector<double>", &truth_particlePX);
    NuETree->Branch("truth_particlePY", "std::vector<double>", &truth_particlePY);
    NuETree->Branch("truth_particlePZ", "std::vector<double>", &truth_particlePZ);
    NuETree->Branch("truth_particleEnergy", "std::vector<double>", &truth_particleEnergy);
    NuETree->Branch("truth_particleAngle", "std::vector<double>", &truth_particleAngle);
    NuETree->Branch("truth_particleETheta2", "std::vector<double>", &truth_particleETheta2);
    NuETree->Branch("truth_particleDirectionX", "std::vector<double>", &truth_particleDirectionX);
    NuETree->Branch("truth_particleDirectionY", "std::vector<double>", &truth_particleDirectionY);
    NuETree->Branch("truth_particleDirectionZ", "std::vector<double>", &truth_particleDirectionZ);

    NuETree->Branch("reco_neutrinoPDG", "std::vector<double>", &reco_neutrinoPDG);
    NuETree->Branch("reco_neutrinoIsPrimary", "std::vector<double>", &reco_neutrinoIsPrimary);
    NuETree->Branch("reco_neutrinoVX", "std::vector<double>", &reco_neutrinoVX);
    NuETree->Branch("reco_neutrinoVY", "std::vector<double>", &reco_neutrinoVY);
    NuETree->Branch("reco_neutrinoVZ", "std::vector<double>", &reco_neutrinoVZ);
    NuETree->Branch("reco_neutrinoSliceID", "std::vector<double>", &reco_neutrinoSliceID);

    NuETree->Branch("reco_particlePDG", "std::vector<double>", &reco_particlePDG);
    NuETree->Branch("reco_particleIsPrimary", "std::vector<double>", &reco_particleIsPrimary);
    NuETree->Branch("reco_particleVX", "std::vector<double>", &reco_particleVX);
    NuETree->Branch("reco_particleVY", "std::vector<double>", &reco_particleVY);
    NuETree->Branch("reco_particleVZ", "std::vector<double>", &reco_particleVZ);
    NuETree->Branch("reco_particleDirectionX", "std::vector<double>", &reco_particleDirectionX);
    NuETree->Branch("reco_particleDirectionY", "std::vector<double>", &reco_particleDirectionY);
    NuETree->Branch("reco_particleDirectionZ", "std::vector<double>", &reco_particleDirectionZ);
    NuETree->Branch("reco_particleSliceID", "std::vector<double>", &reco_particleSliceID);
    NuETree->Branch("reco_particleBestPlaneEnergy", "std::vector<double>", &reco_particleBestPlaneEnergy);
    NuETree->Branch("reco_particleTheta", "std::vector<double>", &reco_particleTheta);
    NuETree->Branch("reco_particleTrackScore", "std::vector<double>", &reco_particleTrackScore);
    NuETree->Branch("reco_particleCompleteness", "std::vector<double>", &reco_particleCompleteness);
    NuETree->Branch("reco_particlePurity", "std::vector<double>", &reco_particlePurity);

    NuETree->Branch("reco_sliceID", "std::vector<double>", &reco_sliceID);
    NuETree->Branch("reco_sliceCompleteness", "std::vector<double>", &reco_sliceCompleteness);
    NuETree->Branch("reco_slicePurity", "std::vector<double>", &reco_slicePurity);
    NuETree->Branch("reco_sliceScore", "std::vector<double>", &reco_sliceScore);

    NuETree->Branch("reco_hitPlane", "std::vector<double>", &reco_hitPlane);
    NuETree->Branch("reco_hitX", "std::vector<double>", &reco_hitX);
    NuETree->Branch("reco_hitUVZ", "std::vector<double>", &reco_hitUVZ);
}

void sbnd::NuE::analyze(art::Event const& e){

  // const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  ClearMaps(e);
  SetupMaps(e);
  clearVectors();

  eventID = e.id().event();
  runID = e.id().run();
  subRunID = e.id().subRun();
  // DLCurrent: 0 = uboone dl, 1 = dune dl, 2 = current, 3 = cheated

  std::cout << "" << std::endl;
  std::cout << "________________________________________________________________________________________" << std::endl;
  std::cout << "Run: " << runID << ", Subrun: " << subRunID << ", Event: " << eventID << ", DL/Current: " << DLCurrent << std::endl;
       
  trueNeutrino(e);
  recoNeutrino(e);
  slices(e);
  PFPs(e);
  showerEnergy(e);
  MCParticles(e);

  NuETree->Fill();
}

double sbnd::NuE::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    std::map<int, int> objectHitsMap;

    for(unsigned int i = 0; i < objectHits.size(); ++i)
        ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

    //std::cout << "Number in object: " << objectHitsMap[ID] << " Number of hits in object in truth: " << fHitsMap[ID] << std::endl;
    return (fHitsMap[ID] == 0) ? def_double : objectHitsMap[ID]/static_cast<double>(fHitsMap[ID]);
}

double sbnd::NuE::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    std::map<int, int> objectHitsMap;

    for(unsigned int i = 0; i < objectHits.size(); ++i)
        ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

    return (objectHits.size() == 0) ? def_double : objectHitsMap[ID]/static_cast<double>(objectHits.size());
}

void sbnd::NuE::ClearMaps(const art::Event &e){
    fHitsMap.clear();
}

void sbnd::NuE::SetupMaps(const art::Event &e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    art::Handle<std::vector<recob::Hit> > hitsHandle;
    e.getByLabel(hitLabel,hitsHandle);

    for(unsigned hit_i = 0; hit_i < hitsHandle->size(); ++hit_i) {
        const art::Ptr<recob::Hit> hit(hitsHandle,hit_i);
        fHitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
    }
}

void sbnd::NuE::showerEnergy(const art::Event &e){
    printf("\n");

    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    art::Handle<std::vector<recob::Shower>> showerHandle;
    std::vector<art::Ptr<recob::Shower>> showerVec;
    if(e.getByLabel(showerLabel, showerHandle))
        art::fill_ptr_vector(showerVec, showerHandle);

    //art::FindManyP<recob::Hit> showerHitAssns(showerVec, e, showerLabel);
   
    //std::map<unsigned int, std::vector<art::Ptr<recob::Hit>>> showerHitsByID;

    //for(const art::Ptr<recob::Shower> &shower : showerVec){
        //const std::vector<art::Ptr<recob::Hit>> showerHits(showerHitAssns.at(shower.key()));
        //std::cout << "Shower ID: " << shower->ID() << ", Number of hits in shower: " << showerHits.size() << std::endl;
        //showerHitsByID[shower->ID()] = showerHits;
    //}

    int counter = 0;

    if(!pfpVec.empty()){
        // Get associations between pfparticles and showers
        art::FindOneP<recob::Shower> pfpShowerAssns(pfpVec, e, showerLabel);
        art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpVec, e, PFParticleLabel);
        art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);
        art::FindOneP<recob::Slice> pfpSliceAssns(pfpVec, e, sliceLabel);    
        art::FindManyP<recob::Hit> showerHitAssns(showerVec, e, showerLabel);
        
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
            if(!(pfp->PdgCode() == std::numeric_limits<int>::max())){
               
                const art::Ptr<recob::Shower> pfpShower(pfpShowerAssns.at(pfp.key()));
                const art::Ptr<recob::Slice> pfpSlice(pfpSliceAssns.at(pfp.key()));
                const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
      
                const auto meta  = pfpMetadataAssns.at(pfp.key());
                const auto props = meta->GetPropertiesMap();
                const auto trackscoreobj = props.find("TrackScore");

                if(pfpVertexs.size() > 0){
                    if(!(pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
                        if(trackscoreobj->second <= 1 && trackscoreobj->second >= 0){
                                                                               
                            const std::vector<art::Ptr<recob::Hit>> showerHits(showerHitAssns.at(pfpShower.key()));

                            const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
                            counter++;

                            //const std::vector<art::Ptr<recob::Hit>> &hits = showerHitsByID[pfpShower->ID()];

                            const int showerID_truth = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, true);
                            
                            reco_particlePDG.push_back(pfp->PdgCode());
                            reco_particleIsPrimary.push_back(pfp->IsPrimary());
                            reco_particleVX.push_back(pfpVertex->position().X());
                            reco_particleVY.push_back(pfpVertex->position().Y());
                            reco_particleVZ.push_back(pfpVertex->position().Z());
                            reco_particleDirectionX.push_back(pfpShower->Direction().X());
                            reco_particleDirectionY.push_back(pfpShower->Direction().Y());
                            reco_particleDirectionZ.push_back(pfpShower->Direction().Z());
                            reco_particleSliceID.push_back(pfpSlice->ID());
                            reco_particleBestPlaneEnergy.push_back(pfpShower->Energy()[pfpShower->best_plane()]);
                            reco_particleTheta.push_back(pfpShower->Direction().Theta());
                            reco_particleTrackScore.push_back(trackscoreobj->second);
                            reco_particleCompleteness.push_back(Completeness(e, showerHits, showerID_truth));
                            reco_particlePurity.push_back(Purity(e, showerHits, showerID_truth));
                            
                            printf("Reco Particle %d: PDG Code = %d, Is Primary = %d, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %d, Best Plane Energy = %f, Theta = %f, Trackscore = %f, Completeness = %f, Purity = %f, Shower ID = %d, Number of Hits in Shower = %ld\n", counter, pfp->PdgCode(), pfp->IsPrimary(), pfpVertex->position().X(), pfpVertex->position().Y(), pfpVertex->position().Z(), pfpShower->Direction().X(), pfpShower->Direction().Y(), pfpShower->Direction().Z(), pfpSlice->ID(), pfpShower->Energy()[pfpShower->best_plane()], pfpShower->Direction().Theta(), trackscoreobj->second, Completeness(e, showerHits, showerID_truth), Purity(e, showerHits, showerID_truth), pfpShower->ID(), showerHits.size());
                        }
                    }
                }
            }
        }
    }
    
    if(counter == 0){
        reco_particlePDG.push_back(-999999);
        reco_particleIsPrimary.push_back(-999999);
        reco_particleVX.push_back(-999999);
        reco_particleVY.push_back(-999999);
        reco_particleVZ.push_back(-999999);
        reco_particleDirectionX.push_back(-999999);
        reco_particleDirectionY.push_back(-999999);
        reco_particleDirectionZ.push_back(-999999);
        reco_particleSliceID.push_back(-999999);
        reco_particleBestPlaneEnergy.push_back(-999999);
        reco_particleTheta.push_back(-999999);
        reco_particleTrackScore.push_back(-999999);
        reco_particleCompleteness.push_back(-999999);
        reco_particlePurity.push_back(-999999);
        printf("Reco Particle %d: PDG Code = %d, Is Primary = %d, Vertex = (%d, %d, %d), Direction = (%d, %d, %d), Slice ID = %d, Best Plane Energy = %d, Theta = %d, Trackscore = %d, Completeness = %d, Purity = %d\n", counter, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999);
    }
}

void sbnd::NuE::recoNeutrino(const art::Event &e){
    // Function for dealing with the reconstructed primary neutrino

    printf("\n");
    // Gets all the PFPs in the event
    art::Handle<std::vector<recob::PFParticle>>   pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>>      pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    art::FindManyP<recob::Slice> pfpSliceAssns(pfpVec, e, sliceLabel);      // Gets association between slices and PFPs
    art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);   // Gets association between PFPs and vertices
    art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpHandle, e, PFParticleLabel);   // Gets association between PFPs and Metadata

    int counter = 0;

    if(!pfpVec.empty()){
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
            if(!(pfp->PdgCode() == std::numeric_limits<int>::max())){
    
                const std::vector<art::Ptr<recob::Slice>> pfpSlices(pfpSliceAssns.at(pfp.key()));

                if(pfpSlices.size() == 1){
                    const art::Ptr<recob::Slice> &pfpSlice(pfpSlices.front());
        
                    bool pfpPrimary = pfp->IsPrimary();

                    if(pfpPrimary && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
                        // This is a primary neutrino
                        const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
                        if(pfpVertexs.size() == 1){
                            counter++;
                            const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
                            reco_neutrinoPDG.push_back(pfp->PdgCode());
                            reco_neutrinoIsPrimary.push_back(pfpPrimary);
                            reco_neutrinoVX.push_back(pfpVertex->position().X());
                            reco_neutrinoVY.push_back(pfpVertex->position().Y());
                            reco_neutrinoVZ.push_back(pfpVertex->position().Z());
                            reco_neutrinoSliceID.push_back(pfpSlice->ID()); 
                            printf("Reco Neutrino %d: PDG Code = %d, Is Primary = %d, Vertex = (%f, %f, %f), Slice ID: %d\n", counter, pfp->PdgCode(), pfpPrimary, pfpVertex->position().X(), pfpVertex->position().Y(), pfpVertex->position().Z(), pfpSlice->ID());
                        }
                    }
                }
            }
        }
    }

    if(counter == 0){
        reco_neutrinoPDG.push_back(-999999);
        reco_neutrinoIsPrimary.push_back(-999999);
        reco_neutrinoVX.push_back(-999999);
        reco_neutrinoVY.push_back(-999999);
        reco_neutrinoVZ.push_back(-999999);
        reco_neutrinoSliceID.push_back(-999999); 
        printf("Reco Neutrino %d: PDG Code = %d, Is Primary = %d, Vertex = (%d, %d, %d), Slice ID: %d\n", counter, -999999, -999999, -999999, -999999, -999999, -999999);
    }
}

void sbnd::NuE::trueNeutrino(const art::Event &e){
    // Function for dealing with the true neutrino

    printf("\n");

    art::Handle<std::vector<simb::MCTruth>> neutrinoHandle;
    e.getByLabel(nuGenModuleLabel, neutrinoHandle);    

    int counter = 0;

    for(unsigned int nu_i = 0; nu_i < neutrinoHandle->size(); ++nu_i){
        const art::Ptr<simb::MCTruth> trueNeutrino(neutrinoHandle, nu_i);
        if(trueNeutrino.isNull()) continue;

        if(trueNeutrino->Origin() == 1){                                    // Neutrino must be a beam neutrino (1 = kBeamNeutrino)
            counter++;
            const simb::MCNeutrino neutrino = trueNeutrino->GetNeutrino();  // References to the neutrino
            const simb::MCParticle neutrinoParticle = neutrino.Nu();        // The incoming neutrino
            const simb::MCParticle lepton = neutrino.Lepton();              // The outgoing lepton
        
            geo::Point_t pos = {neutrinoParticle.Vx(), neutrinoParticle.Vy(), neutrinoParticle.Vz()};
            geo::TPCID tpcID = theGeometry->FindTPCAtPosition(pos);

            std::cout << "TPC ID Num: " << tpcID.TPC << ", is Valid: " << tpcID.isValid << std::endl;
            truth_neutrinoVX.push_back(neutrinoParticle.Vx());
            truth_neutrinoVY.push_back(neutrinoParticle.Vy());
            truth_neutrinoVZ.push_back(neutrinoParticle.Vz());
            truth_CCNC.push_back(neutrino.CCNC());
            truth_neutrinoType.push_back(neutrinoParticle.PdgCode());    
            truth_chargedLepton.push_back(lepton.PdgCode());
            truth_neutrinoTPCID.push_back(tpcID.TPC);
            truth_neutrinoTPCValid.push_back(static_cast<double>(tpcID.isValid));
            printf("True Neutrino %d: Vertex = (%f, %f, %f), CCNC = %d, Neutrino Type = %d, Charged Lepton = %d, TPC ID = %d, TPC Valid = %f\n", counter, neutrinoParticle.Vx(), neutrinoParticle.Vy(), neutrinoParticle.Vz(), neutrino.CCNC(), neutrinoParticle.PdgCode(), lepton.PdgCode(), tpcID.TPC, static_cast<double>(tpcID.isValid));
        }
    }

    if(counter == 0){
        truth_neutrinoVX.push_back(-999999);
        truth_neutrinoVY.push_back(-999999);
        truth_neutrinoVZ.push_back(-999999);
        truth_CCNC.push_back(-999999);
        truth_neutrinoType.push_back(-999999);    
        truth_chargedLepton.push_back(-999999);
        truth_neutrinoTPCID.push_back(-999999);
        truth_neutrinoTPCValid.push_back(-999999);
        printf("True Neutrino %d: Vertex = (%d, %d, %d), CCNC = %d, Neutrino Type = %d, Charged Lepton = %d, TPC ID = %d, TPC Valid = %d\n", counter, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999);
    }
}

void sbnd::NuE::slices(const art::Event &e){
    printf("\n");
    
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
   
    // Gets all the slices in the event
    art::Handle<std::vector<recob::Slice>>  sliceHandle;
    std::vector<art::Ptr<recob::Slice>>     sliceVec;
    if(e.getByLabel(sliceLabel, sliceHandle))
        art::fill_ptr_vector(sliceVec, sliceHandle);
   
    int counter = 0;

    if(!sliceVec.empty()){;
        // Gets all the hits in the event
        art::Handle<std::vector<recob::Hit>> hitHandle;
        std::vector<art::Ptr<recob::Hit>> hitVec;
        if(e.getByLabel(hitLabel, hitHandle))
            art::fill_ptr_vector(hitVec, hitHandle);
        if(!hitVec.empty()){

            int sliceID(std::numeric_limits<int>::max());

            for(const art::Ptr<recob::Slice> &slice : sliceVec){
                counter++;
                sliceID = slice->ID();
                if(sliceID == std::numeric_limits<int>::max()) continue;

                art::FindManyP<recob::Hit> sliceHitAssns(sliceVec, e, sliceLabel);
                const std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssns.at(slice.key()));

                art::FindManyP<recob::PFParticle> slicePFPAssns(sliceVec, e, sliceLabel);
                const std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssns.at(slice.key()));
        
                const int sliceID_truth = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hitVec, true);
                double sliceCompleteness = Completeness(e, sliceHits, sliceID_truth);

                art::FindManyP<sbn::CRUMBSResult> sliceCrumbsAssns(sliceVec, e, crumbsLabel);
                const std::vector<art::Ptr<sbn::CRUMBSResult>> sliceCrumbsResults = sliceCrumbsAssns.at(slice.key());
       
                reco_sliceID.push_back(sliceID);
                reco_sliceCompleteness.push_back(sliceCompleteness);
                reco_slicePurity.push_back(Purity(e, sliceHits, sliceID_truth));

                double sliceScoreVar = 0;

                if(sliceCrumbsResults.size() == 1){
                    const art::Ptr<sbn::CRUMBSResult> sliceCrumbsResult(sliceCrumbsResults.front());
                    sliceScoreVar = sliceCrumbsResult->score;
                } else{
                    sliceScoreVar = -999999;
                }
                reco_sliceScore.push_back(sliceScoreVar);
                printf("Slice %d: ID = %d, Completeness = %f, Purity = %f, Score = %f\n", counter, sliceID, sliceCompleteness, Purity(e, sliceHits, sliceID_truth), sliceScoreVar);
            }
        }
    }

    if(counter == 0){
        reco_sliceID.push_back(-999999);
        reco_sliceCompleteness.push_back(-999999);
        reco_slicePurity.push_back(-999999);
        reco_sliceScore.push_back(-999999);
        printf("Slice %d: ID = %d, Completeness = %d, Purity = %d, Score = %d\n", counter, -999999, -999999, -999999, -999999);
    }
}

void sbnd::NuE::PFPs(const art::Event &e){
    printf("\n");
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    int counter = 0;
    for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
        if(!(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14))){
            counter++;
        }
    }

    std::cout << "Number of PFPs: " << counter << std::endl;
}

void sbnd::NuE::MCParticles(const art::Event &e){
    printf("\n");
    art::Handle<std::vector<simb::MCParticle>> particleHandle;
    std::vector<art::Ptr<simb::MCParticle>> particleVec;
    if(e.getByLabel(MCTruthLabel, particleHandle))
        art::fill_ptr_vector(particleVec, particleHandle);

    int counter = 0;

    if(!particleVec.empty()){
        for(auto &particle : particleVec){
            if(particle->PdgCode() == 11 && particle->Mother() == 0){
                counter++;
            
                double P_mag = std::sqrt((particle->Px() * particle->Px()) + (particle->Py() * particle->Py()) + (particle->Pz() * particle->Pz()));
                double energy = particle->E() * 1000; // Converts from GeV to MeV
                double theta = std::acos(particle->Pz() / P_mag);
                double true_ETheta2 = (energy * (theta * theta));
                double DX = (particle->Px() / P_mag);
                double DY = (particle->Py() / P_mag);
                double DZ = (particle->Pz() / P_mag);

                truth_particlePDG.push_back(particle->PdgCode());
                truth_particleVX.push_back(particle->Vx());
                truth_particleVY.push_back(particle->Vy());
                truth_particleVZ.push_back(particle->Vz());
                truth_particlePX.push_back(particle->Px());
                truth_particlePY.push_back(particle->Py());
                truth_particlePZ.push_back(particle->Pz());
                truth_particleEnergy.push_back(energy);
                truth_particleAngle.push_back(theta);
                truth_particleETheta2.push_back(true_ETheta2);
                truth_particleDirectionX.push_back(DX);
                truth_particleDirectionY.push_back(DY);
                truth_particleDirectionZ.push_back(DZ);
                printf("True Particle %d: PDG Code = %d, Vertex = (%f, %f, %f), Momentum = (%f, %f, %f), Energy = %f, Theta = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", counter, particle->PdgCode(), particle->Vx(), particle->Vy(), particle->Vz(), particle->Px(), particle->Py(), particle->Pz(), energy, theta, true_ETheta2, DX, DY, DZ);
            }
        }
    }

    if(counter == 0){
        truth_particlePDG.push_back(-999999);
        truth_particleVX.push_back(-999999);
        truth_particleVY.push_back(-999999);
        truth_particleVZ.push_back(-999999);
        truth_particlePX.push_back(-999999);
        truth_particlePY.push_back(-999999);
        truth_particlePZ.push_back(-999999);
        truth_particleEnergy.push_back(-999999);
        truth_particleAngle.push_back(-999999);
        truth_particleETheta2.push_back(-999999);
        truth_particleDirectionX.push_back(-999999);
        truth_particleDirectionY.push_back(-999999);
        truth_particleDirectionZ.push_back(-999999);
        printf("True Particle %d: PDG Code = %d, Vertex = (%d, %d, %d), Momentum = (%d, %d, %d), Energy = %d, Theta = %d, ETheta2 = %d, Direction = (%d, %d, %d)\n", counter, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999);
    }
}

void sbnd::NuE::hits(const art::Event &e){

    detinfo::DetectorPropertiesData propD = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);

    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVec;
    if(e.getByLabel(hitLabel, hitHandle))
        art::fill_ptr_vector(hitVec, hitHandle);
 
    int counter = 0;

    if(!hitVec.empty()){

        for(auto &hit : hitVec){
            if(hit->View() == 0 || hit->View() == 1 || hit->View() == 2){
                counter++;
                const geo::WireID& wireID(hit->WireID());
                const double xpos = propD.ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);
           
                auto const& channelMapAlg =  art::ServiceHandle<geo::WireReadout const>()->Get();

                auto const wire = channelMapAlg.Plane(geo::PlaneID(wireID.Cryostat, wireID.TPC, wireID.Plane)).Wire(wireID);
            
                geo::Point_t start = wire.GetStart();
                geo::Point_t end = wire.GetEnd();

                const double ay(start.Y());
                const double az(start.Z());
                const double by(end.Y());
                const double bz(end.Z());

                const double ny(by - ay);
                const double nz(bz - az);
                const double n2(ny * ny + nz * nz);

                const double ry(ay - (ay * ny + az * nz) * ny / n2);
                const double rz(az - (ay * ny + az * nz) * nz / n2);
                const double sign((rz > 0.0) ? +1.0 : -1.0);
                const double uvz = sign * std::sqrt(ry * ry + rz * rz); // this is the U/V/Z coordinate
                
                reco_hitPlane.push_back(hit->View());
                reco_hitX.push_back(xpos);
                reco_hitUVZ.push_back(uvz);
            }
        }
    }

    if(counter == 0){
        reco_hitPlane.push_back(-999999);
        reco_hitX.push_back(-999999);
        reco_hitUVZ.push_back(-999999);
    }
}

void sbnd::NuE::clearVectors(){
    truth_neutrinoVX.clear();
    truth_neutrinoVY.clear();
    truth_neutrinoVZ.clear();
    truth_CCNC.clear();
    truth_neutrinoType.clear();
    truth_chargedLepton.clear();
    truth_neutrinoTPCID.clear();
    truth_neutrinoTPCValid.clear();

    // Truth particle info
    truth_particlePDG.clear();
    truth_particleVX.clear();
    truth_particleVY.clear();
    truth_particleVZ.clear();
    truth_particlePX.clear();
    truth_particlePY.clear();
    truth_particlePZ.clear();
    truth_particleEnergy.clear();
    truth_particleAngle.clear();
    truth_particleETheta2.clear();
    truth_particleDirectionX.clear();
    truth_particleDirectionY.clear();
    truth_particleDirectionZ.clear();

    // Reco neutrino info
    reco_neutrinoPDG.clear();
    reco_neutrinoIsPrimary.clear();
    reco_neutrinoVX.clear();
    reco_neutrinoVY.clear();
    reco_neutrinoVZ.clear();
    reco_neutrinoSliceID.clear();

    // Reco particle info
    reco_particlePDG.clear();
    reco_particleIsPrimary.clear();
    reco_particleVX.clear();
    reco_particleVY.clear();
    reco_particleVZ.clear();
    reco_particleDirectionX.clear();
    reco_particleDirectionY.clear();
    reco_particleDirectionZ.clear();
    reco_particleSliceID.clear();
    reco_particleBestPlaneEnergy.clear();
    reco_particleTheta.clear();
    reco_particleTrackScore.clear();
    reco_particleCompleteness.clear();
    reco_particlePurity.clear();

    // Reco slice info
    reco_sliceID.clear();
    reco_sliceCompleteness.clear();
    reco_slicePurity.clear();
    reco_sliceScore.clear();

    // Reco hit info
    reco_hitPlane.clear();
    reco_hitX.clear();
    reco_hitUVZ.clear();
}

void sbnd::NuE::beginJob()
{
}
void sbnd::NuE::endJob()
{
}

DEFINE_ART_MODULE(sbnd::NuE)
