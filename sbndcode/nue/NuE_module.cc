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
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NuE(NuE const&) = delete;
    NuE(NuE&&) = delete;
    NuE& operator=(NuE const&) = delete;
    NuE& operator=(NuE&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    double Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID);
    double Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID);
    void ClearMaps(const art::Event &e);
    void SetupMaps(const art::Event &e);
    void ShowerEnergy(const art::Event &e, const std::vector<art::Ptr<recob::PFParticle>> &pfpVec, const double VX, const double VY, const double VZ);
    void recoNeutrino(const art::Event &e);
    void trueNeutrino(const art::Event &e);

private:

    art::ServiceHandle<cheat::BackTrackerService>       backTracker;

    TTree *NuETree = new TTree("NuETree", "NuE");
    TTree *NuEHitTree = new TTree("NuEHitTree", "NuEHit");

    unsigned int eventID;                    // Event num
    unsigned int runID;                      // Run num  
    unsigned int subRunID;                   // Subrun num
    unsigned int nPFParticles;               // Number of PFParticles
    unsigned int nSlices;                    // Number of Slices
    unsigned int nNeutrinos;                 // Number of MCTruth Neutrinos
    unsigned int nTracks;                    // Number of reconstructed Tracks
    unsigned int nShowers;                   // Number of reconstructed Showers

    // Truth variables
    double trueNeutrinoVX;                    // MCTruth vertex x coord
    double trueNeutrinoVY;                    // MCTruth vertex y coord
    double trueNeutrinoVZ;                    // MCTruth vertex z coord
    int trueCCNC;                             // MCTruth CC or NC, CC = 0, NC = 1 
    int trueNeutrinoType;                     // MCTruth neutrino flavour (PDG code)
    int trueChargedLepton;                    // MCTruth charged lepton flavour (PDG code)
  
    // Reco variables
    double recoNeutrinoVX;                    // Reco vertex x coord
    double recoNeutrinoVY;                    // Reco vertex y coord
    double recoNeutrinoVZ;                    // Reco vertex z coord

    double ETheta2;
    double totalShowerEnergy;                 // Total energy of the shower in the event
    double smallestDeltaR;
    double largestEnergy;
    double chosenShowerTheta;
    double smallestTheta;

    double chosenSliceCompleteness;           // The completeness of the chosen slice

    double numEventsTotal = 0;
    double numEventsSelected = 0;

    float highestSliceScore = -10000;
    int highestSliceScoreIndex;

    std::map<int,int> fHitsMap;
    std::vector<std::vector<double>> sliceInfo;

    geo::Point_t pos = {0, 0, 0};             // Geometry position (0, 0, 0), gets overwritten with actual coords

    // Geometry information
    detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    detinfo::DetectorPropertiesData propData = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
    art::ServiceHandle<geo::Geometry> theGeometry;
    geo::GeometryCore const* geom = art::ServiceHandle<geo::Geometry>()->provider();
 
    double tpcID_num = 100;

    // Vectors for the NuE Tree
    std::vector<double> fullyReco_tree = std::vector<double>(0);
    std::vector<int> DLCurrent_tree = std::vector<int>(0);
    std::vector<double> event_tree = std::vector<double>(0);
    std::vector<double> run_tree = std::vector<double>(0);
    std::vector<double> subrun_tree = std::vector<double>(0);
    std::vector<double> trueNeutrinoVX_tree = std::vector<double>(0);
    std::vector<double> trueNeutrinoVY_tree = std::vector<double>(0);
    std::vector<double> trueNeutrinoVZ_tree = std::vector<double>(0);
    std::vector<int> trueCCNC_tree = std::vector<int>(0);
    std::vector<int> trueNeutrinoType_tree = std::vector<int>(0);
    std::vector<int> trueLeptonType_tree = std::vector<int>(0);
    std::vector<double> recoNeutrinoVX_tree = std::vector<double>(0);
    std::vector<double> recoNeutrinoVY_tree = std::vector<double>(0);
    std::vector<double> recoNeutrinoVZ_tree = std::vector<double>(0);
    std::vector<int> numSlices_tree = std::vector<int>(0);
    std::vector<double> sliceCompleteness_tree = std::vector<double>(0);
    std::vector<double> sliceNumPFPs_tree = std::vector<double>(0);
    std::vector<int> numPFPs_tree = std::vector<int>(0); 
    std::vector<double> showerEnergy_tree = std::vector<double>(0);
    std::vector<double> showerTheta_tree = std::vector<double>(0);
    std::vector<double> showerTrackScore_tree = std::vector<double>(0);
    std::vector<double> numShowers_tree = std::vector<double>(0);
    std::vector<double> showerETheta2_tree = std::vector<double>(0);

    // Vectors for the NuEHit Tree
    std::vector<double> event_hitTree = std::vector<double>(0);
    std::vector<double> run_hitTree = std::vector<double>(0);
    std::vector<double> subrun_hitTree = std::vector<double>(0);
    std::vector<int> plane_hitTree = std::vector<int>(0);
    std::vector<double> x_hitTree = std::vector<double>(0);
    std::vector<double> uvz_hitTree = std::vector<double>(0);

    // Label strings
    const std::string PFParticleLabel;
    const std::string sliceLabel;
    const std::string vertexLabel;
    const std::string nuGenModuleLabel;
    const std::string hitLabel;
    const std::string crumbsLabel;
    const std::string trackLabel;
    const std::string showerLabel;

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
  showerLabel(p.get<std::string>("ShowerLabel"))
{
}

void sbnd::NuE::analyze(art::Event const& e){

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  ClearMaps(e);
  SetupMaps(e);

  eventID = e.id().event();
  runID = e.id().run();
  subRunID = e.id().subRun();

  std::cout << "" << std::endl;
  std::cout << "________________________________________________________________________________________" << std::endl;
  std::cout << "Run: " << runID << ", Subrun: " << subRunID << ", Event: " << eventID << std::endl;

  for(auto &hit : hitVec){
    if(hit->View() == 0 || hit->View() == 1 || hit->View() == 2){
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
 
        event_hitTree.push_back(eventID);
        run_hitTree.push_back(runID);
        subrun_hitTree.push_back(subRunID);
        plane_hitTree.push_back(hit->View());
        x_hitTree.push_back(xpos);
        uvz_hitTree.push_back(uvz);
    }
  }

  if(!std::isnan(recoNeutrinoVX)){
    // Event is fully reconstructed
    fullyReco_tree.push_back(1);
  } else{
      std::cout << "There is no reco neutrino, skipping event" << std::endl;
      // Event is not fully reconstructed
      fullyReco_tree.push_back(0);
  }
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

void sbnd::NuE::ShowerEnergy(const art::Event &e){
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);
    if(pfpVec.empty()) return;
    
    // Get associations between pfparticles and showers
    art::FindOneP<recob::Shower> pfpShowerAssns(pfpVec, e, showerLabel);
    art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpVec, e, PFParticleLabel);
    art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);

    double numShowers = 0;
    double totalShowerEnergy;
    double chosenShowerTheta;
    double chosenShowerTrackScore;
    double largestEnergy = 0;

    for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
        if(pfp->PdgCode() == std::numeric_limits<int>::max()) return;
  
        const art::Ptr<recob::Shower> pfpShower(pfpShowerAssns.at(pfp.key()));
        const auto meta  = pfpMetadataAssns.at(pfp.key());
        const auto props = meta->GetPropertiesMap();
        const auto trackscoreobj = props.find("TrackScore");
        const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
        if(pfpVertexs.size() > 0){
            if(!(pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
                if(trackscoreobj->second <= 1 && trackscoreobj->second >= 0){
                    numShowers++;
                    totalShowerEnergy += pfpShower->Energy()[pfpShower->best_plane()];
        
                    if(pfpShower->Energy()[pfpShower->best_plane()] > largestEnergy){
                        chosenShowerTheta = pfpShower->Direction().Theta();
                        chosenShowerTrackScore = trackscoreobj->second;
                        largestEnergy = pfpShower->Energy()[pfpShower->best_plane()];
                    }
                }
            }
        }
    }

    double ETheta2 = (totalShowerEnergy * std::pow(chosenShowerTheta, 2));

    std::cout << "Number of showers: " << numShowers << std::endl;
    std::cout << "Chosen Shower: Energy = " << largestEnergy << ", Theta = " << chosenShowerTheta << ", Trackscore = " << chosenShowerTrackScore << std::endl;
    std::cout << "Shower Energy = " << totalShowerEnergy << ", Shower Theta = " << chosenShowerTheta << ", ETheta^2 = " << ETheta2 << std::endl;
    numShowers_tree.push_back(numShowers);
    showerEnergy_tree.push_back(totalShowerEnergy);
    showerTheta_tree.push_back(chosenShowerTheta);
    showerTrackScore_tree.push_back(chosenShowerTrackScore);
    showerETheta2_tree.push_back(ETheta2);
}

void sbnd::NuE::recoNeutrino(const art::Event &e){
    // Function for dealing with the reconstructed primary neutrino
    
    double numRecoNeutrinos = 0;
    std::vector<std::tuple<art::Ptr<recob::PFParticle>, art::Ptr<recob::Vertex>, art::Ptr<recob::Slice>, art::Ptr<sbn::CRUMBSResult>>> v;   // Initialise the vector that holds the primary reco neutrino info

    // Gets all the PFPs in the event
    art::Handle<std::vector<recob::PFParticle>>   pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>>      pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    art::FindManyP<recob::Slice> pfpSliceAssns(pfpVec, e, sliceLabel);      // Gets association between slices and PFPs
    art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);   // Gets association between PFPs and vertices
    art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpHandle, e, PFParticleLabel);   // Gets association between PFPs and Metadata

    if(pfpVec.empty()) return;      // Skips the event if there are no PFPs.
    if(sliceVec.empty()) return;    // Skips the event if there are no slices.

    for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
        if(pfp->PdgCode() == std::numeric_limits<int>::max()) return;
    
        const std::vector<art::Ptr<recob::Slice>> pfpSlices(pfpSliceAssns.at(pfp.key()));

        if(pfpSlices.size() == 1){
             const art::Ptr<recob::Slice> &pfpSlice(pfpSlices.front());
        
             bool pfpPrimary = pfp->IsPrimary();

             if(pfpPrimary && (PFParticleID == 12 || PFParticleID == 14)){
                 // This is a primary neutrino
                 numRecoNeutrinos++;

                 const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
                 if(pfpVertexs.size() == 1){
                     const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());

                     art::FindManyP<sbn::CRUMBSResult> sliceCrumbsAssns(pfpSlices, e, crumbsLabel);
                     std::tuple<art::Ptr<recob::PFParticle>, art::Ptr<recob::Vertex>, art::Ptr<recob::Slice>, art::Ptr<sbn::CRUMBSResult>> tuple;

                     const std::vector<art::Ptr<sbn::CRUMBSResult>> sliceCrumbsResults = sliceCrumbsAssns.at(0);
                     
                     if(sliceCrumbsResults.size() != 1){
                         std::cout << "Error: Slice has multiple CRUMBS results" << std::endl;
                     }

                     const art::Ptr<sbn::CRUMBSResult> sliceCrumbsResult(sliceCrumbsResults.front());

                     tuple = std::make_tuple(pfp, pfpVertex, pfpSlice, sliceCrumbsResult); // Creates a tuple that contains the pfp, and the vertex coords
                 
                     v.push_back(tuple);
                 }
             }
        }
    }

    if(v.size() != 0){
        float highestSliceScore = -10000;
        int highestSliceScoreIndex;

        for(long unsigned int neutrinoCand = 0; neutrinoCand < v.size(); neutrinoCand++){
            auto crumbsPtr = std::get<3>(v[neutrinoCand]);
            if(crumbsPtr->score > highestSliceScore){
                highestSliceScore = crumbsPtr->score;
                highestSliceScoreIndex = neutrinoCand;
            }
        }

        recoNeutrinoVX = std::get<1>(v[highestSliceScoreIndex])->position().X();
        recoNeutrinoVY = std::get<1>(v[highestSliceScoreIndex])->position().Y();
        recoNeutrinoVZ = std::get<1>(v[highestSliceScoreIndex])->position().Z();

    } else{
        // There is no reconstructed neutrino in the event
        recoNeutrinoVX = NAN;
        recoNeutrinoVY = NAN;
        recoNeutrinoVZ = NAN;
    }

    std::cout << "Reconstructed Primary Neutrino Vertex: (" << recoNeutrinoVX << ", " << recoNeutrinoVY << ", " << recoNeutrinoVZ << ")" << std::endl;

    recoNeutrinoVX_tree.push_back(recoNeutrinoVX);
    recoNeutrinoVY_tree.push_back(recoNeutrinoVY);
    recoNeutrinoVZ_tree.push_back(recoNeutrinoVZ);
}

void sbnd::NuE::trueNeutrino(const art::Event &e){
    // Function for dealing with the true neutrino
 
    art::Handle<std::vector<simb::MCTruth>> neutrinoHandle;
    e.getByLabel(nuGenModuleLabel, neutrinoHandle);    

    int numTrueNeutrinos = 0;
 
    for(unsigned int nu_i = 0; nu_i < neutrinoHandle->size(); ++nu_i){
        const art::Ptr<simb::MCTruth> trueNeutrino(neutrinoHandle, nu_i);
        if(trueNeutrino.isNull()){
            continue;
        }

        if(trueNeutrino->Origin() == 1){                                    // Neutrino must be a beam neutrino (1 = kBeamNeutrino)
            numTrueNeutrinos++;

            const simb::MCNeutrino neutrino = trueNeutrino->GetNeutrino();  // References to the neutrino
            const simb::MCParticle neutrinoParticle = neutrino.Nu();        // The incoming neutrino
            const simb::MCParticle lepton = neutrino.Lepton();              // The outgoing lepton
        
            geo::Point_t pos = {neutrinoParticle.Vx(), neutrinoParticle.Vy(), neutrinoParticle.Vz()}
            geo::TPCID tpcID = theGeometry->FindTPCAtPosition(pos);

            if(tpcID.isValid){
                trueNeutrinoVX = neutrinoParticle.Vx();
                trueNeutrinoVY = neutrinoParticle.Vy();
                trueNeutrinoVZ = neutrinoParticle.Vz();
                trueCCNC = neutrino.CCNC();
                trueNeutrinoType = neutrinoParticle.PdgCode();
                trueChargedLepton = lepton.PdgCode();
                
            }
        }
    }
    
    std::cout << "Number of True Neutrinos = " << numTrueNeutrinos << std::endl;
    
    if(numTrueNeutrinos == 1){
        trueNeutrinoVX_tree.push_back(trueNeutrinoVX);
        trueNeutrinoVY_tree.push_back(trueNeutrinoVY);
        trueNeutrinoVZ_tree.push_back(trueNeutrinoVZ);
        trueCCNC_tree.push_back(trueCCNC);
        trueNeutrinoType_tree.push_back(trueNeutrinoType);
        trueLeptonType_tree.push_back(trueChargedLepton);
    } else{
        std::cout << "Unsure which true neutrino to put in tree" << std::endl;
    }
}

void sbnd::NuE::slices(const art::Event &e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
   
    std::vector<<std::tuple<art::Ptr<recob::Slice>, double, art::Ptr<sbn::CRUMBSResult>, double>> v; 

    // Gets all the slices in the event
    art::Handle<std::vector<recob::Slice>>  sliceHandle;
    std::vector<art::Ptr<recob::Slice>>     sliceVec;
    if(e.getByLabel(sliceLabel, sliceHandle))
        art::fill_ptr_vector(sliceVec, sliceHandle);
    if(sliceVec.empty()) return;

    // Gets all the hits in the event
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVec;
    if(e.getByLabel(hitLabel, hitHandle))
        art::fill_ptr_vector(hitVec, hitHandle);
    if(hitVec.empty()) return;

    int sliceID(std::numeric_limits<int>::max());
    int numSlices = 0;

    for(const art::Ptr<recob::Slice> &slice : sliceVec){
        numSlices++;
        sliceID = slice->ID();
        if(sliceID == std::numeric_limits<int>::max()) return;

        art::FindManyP<recob::Hit> sliceHitAssns(sliceVec, e, sliceLabel);
        const std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssns.at(slice.key()));

        art::FindManyP<sbn::CRUMBSResult> sliceCrumbsAssns(pfpSlices, e, crumbsLabel);
        const std::vector<art::Ptr<sbn::CRUMBSResult>> sliceCrumbsResults = sliceCrumbsAssns.at(0);
        if(sliceCrumbsResults.size() != 1) std::cout << "Error: slice has multiple CRUMBS results" << std::endl;

        const art::Ptr<sbn::CRUMBSResult> sliceCrumbsResult(sliceCrumbsResults.front());

        art::FindManyP<recob::PFParticle> slicePFPAssns(sliceVec, e, sliceLabel);
        const std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssns.at(slice.key()));
        double numPFPs = slicePFPs.size();
        std::cout << "Number of PFPs in slice: " << numPFPs << std::endl;

        const int sliceID_truth = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hitVec, true);
        
        double sliceCompleteness = Completeness(e, sliceHits, sliceID_truth);
        std::cout << "Slice Completeness: " << sliceCompleteness << std::endl;

        std::tuple<art::Ptr<recob::Slice>, double, art::Ptr<sbn::CRUMBSResult>, double> tuple = std::make_tuple(slice, sliceCompleteness, sliceCrumbsResult, numPFPs);
        v.push_back(tuple);
    }

    float highestSliceScore = -10000;
    int highestSliceScoreIndex;

    std::cout << "Number of slices in event: " << numSlices << std::endl;
    if(v.size() != 0){
        for(long unsigned int sliceCand = 0; sliceCand < v.size(); sliceCand++){
            auto crumbsPtr = std::get<2>(v[sliceCand]);
            if(crumbsPtr->score > highestSliceScore){
                highestSliceScore = crumbsPtr->score;
                highestSliceScoreIndex = sliceCand;
            }
        }

        double chosenSliceCompleteness = std::get<1>(v[highestSliceScoreIndex]);
        double chosenSliceNumPFPs = std::get<3>(v[highestSliceScoreIndex]);


        numSlices_tree.push_back(numSlices);
        sliceCompleteness_tree.push_back(chosenSliceCompleteness);
        sliceNumPFPs_tree.push_back(chosenSliceNumPFPs);
    }
}

void sbnd::NuE::PFPs(const art::Event &e){
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    double numPFPs = 0;
    for(const art::Ptr<recob::PFParticle &pfp : pfpVec){
        if(!(pfpPrimary && (PFParticleID == 12 || PFParticleID == 14))){
            numPFPs++;
        }
    }

    std::cout << "Number of PFPs: " << numPFPs << std::endl;
}

void sbnd::NuE::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    NuETree->Branch("fullyReco_tree", &fullyReco_tree);
    NuETree->Branch("DLCurrent_tree", &DLCurrent_tree);
    NuETree->Branch("event_tree", &event_tree);
    NuETree->Branch("run_tree", &run_tree);
    NuETree->Branch("subrun_tree", &subrun_tree);
    NuETree->Branch("trueNeutrinoVX_tree", &trueNeutrinoVX_tree);
    NuETree->Branch("trueNeutrinoVY_tree", &trueNeutrinoVY_tree);
    NuETree->Branch("trueNeutrinoVZ_tree", &trueNeutrinoVZ_tree);
    NuETree->Branch("trueCCNC_tree", &trueCCNC_tree);
    NuETree->Branch("trueNeutrinoType_tree", &trueNeutrinoType_tree);
    NuETree->Branch("trueLeptonType_tree", &trueLeptonType_tree);
    NuETree->Branch("recoNeutrinoVX_tree", &recoNeutrinoVX_tree);
    NuETree->Branch("recoNeutrinoVY_tree", &recoNeutrinoVY_tree);
    NuETree->Branch("recoNeutrinoVZ_tree", &recoNeutrinoVZ_tree); 
    NuETree->Branch("numSlices_tree", &numSlices_tree);
    NuETree->Branch("sliceCompleteness_tree", &sliceCompleteness_tree);
    NuETree->Branch("sliceNumPFPs_tree", &sliceNumPFPs_tree);
    NuETree->Branch("numPFPs_tree", &numPFPs_tree);
    NuETree->Branch("showerEnergy_tree", &showerEnergy_tree);
    NuETree->Branch("showerTheta_tree", &showerTheta_tree);
    NuETree->Branch("showerTrackScore_tree", &showerTrackScore_tree)
    NuETree->Branch("numShowers_tree", &numShowers_tree);
    NuETree->Branch("showerETheta2_tree", &showerETheta2_tree);

    NuEHitTree->Branch("event_hitTree", &event_hitTree);
    NuEHitTree->Branch("run_hitTree", &run_hitTree);
    NuEHitTree->Branch("subrun_hitTree", &subrun_hitTree);
    NuEHitTree->Branch("plane_hitTree", &plane_hitTree);
    NuEHitTree->Branch("x_hitTree", &x_hitTree);
    NuEHitTree->Branch("uvz_hitTree", &uvz_hitTree);
}

void sbnd::NuE::endJob()
{
    NuETree->Fill();
    NuEHitTree->Fill();

    //std::cout << "events selected: " << numEventsSelected << "/" << numEventsTotal << std::endl;
}

DEFINE_ART_MODULE(sbnd::NuE)
