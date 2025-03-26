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
  int numTrueNeutrinos = 0;                 // Counter for number of truth neutrinos 
  
  // Reco variables
  double recoNeutrinoVX;                    // Reco vertex x coord
  double recoNeutrinoVY;                    // Reco vertex y coord
  double recoNeutrinoVZ;                    // Reco vertex z coord
  size_t numRecoNeutrinos = 0;              // Counter for number of reco neutrinos

  double totalShowerEnergy;                 // Total energy of the shower in the event
  double smallestDeltaR;
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
  std::vector<double> event_tree = std::vector<double>(0);
  std::vector<double> run_tree = std::vector<double>(0);
  std::vector<double> subrun_tree = std::vector<double>(0);
  std::vector<double> trueNeutrinoVX_tree = std::vector<double>(0);
  std::vector<double> trueNeutrinoVY_tree = std::vector<double>(0);
  std::vector<double> trueNeutrinoVZ_tree = std::vector<double>(0);
  std::vector<double> recoNeutrinoVX_tree = std::vector<double>(0);
  std::vector<double> recoNeutrinoVY_tree = std::vector<double>(0);
  std::vector<double> recoNeutrinoVZ_tree = std::vector<double>(0);
  std::vector<int> trueCCNC_tree = std::vector<int>(0);
  std::vector<int> trueNeutrinoType_tree = std::vector<int>(0);
  std::vector<int> trueLeptonType_tree = std::vector<int>(0);
  std::vector<int> numSlices_tree = std::vector<int>(0);
  std::vector<int> numPfps_tree = std::vector<int>(0); 
  std::vector<int> numTracks_tree = std::vector<int>(0);
  std::vector<int> numShowers_tree = std::vector<int>(0);
  std::vector<double> tpcID_tree = std::vector<double>(0);
  std::vector<int> dl_current_tree = std::vector<int>(0);
  std::vector<double> sliceCompleteness_tree = std::vector<double>(0);
  std::vector<double> trackScore_tree = std::vector<double>(0);
  std::vector<double> showerEnergy_tree = std::vector<double>(0);
  std::vector<double> showerTheta_tree = std::vector<double>(0);
  std::vector<double> showerSmallestTheta_tree = std::vector<double>(0);
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

  double trackscore;
  totalShowerEnergy = 0;

  eventID = e.id().event();
  runID = e.id().run();
  subRunID = e.id().subRun();

  std::cout << "" << std::endl;
  std::cout << "________________________________________________________________________________________" << std::endl;
  std::cout << "Run: " << runID << ", Subrun: " << subRunID << ", Event: " << eventID << std::endl;

  numTrueNeutrinos = 0;
  numRecoNeutrinos = 0;
  nSlices = 0;
  nPFParticles = 0;

  std::map<int, std::vector<std::tuple<art::Ptr<recob::PFParticle>, double, double, double>>> slicePFParticleMap;

  art::Handle<std::vector<recob::PFParticle>>   pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>>      pfpVec;
  art::Handle<std::vector<simb::MCTruth>>       neutrinoHandle;
  art::Handle<std::vector<recob::Track>>        trackHandle;
  std::vector<art::Ptr<recob::Track>>           trackVec;
  art::Handle<std::vector<recob::Shower>>       showerHandle;
  std::vector<art::Ptr<recob::Shower>>          showerVec;  

  e.getByLabel(nuGenModuleLabel, neutrinoHandle);

  detinfo::DetectorPropertiesData propD = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
  art::ServiceHandle<geo::Geometry const> theGeometry;

  nNeutrinos = neutrinoHandle->size();

  for(unsigned int nu_i = 0; nu_i < nNeutrinos; ++nu_i){
    const art::Ptr<simb::MCTruth> trueNeutrino(neutrinoHandle, nu_i);
    if(trueNeutrino.isNull()){
        continue;
    }

    if(trueNeutrino->Origin() == 1){                                   // Neutrino must be a beam neutrino (1 = kBeamNeutrino)
        const simb::MCNeutrino neutrino = trueNeutrino->GetNeutrino(); // Reference to the neutrino
        const simb::MCParticle neutrinoParticle = neutrino.Nu();        // The incoming neutrino
        const simb::MCParticle lepton = neutrino.Lepton();              // The outgoing lepton
   
        pos = {neutrinoParticle.Vx(), neutrinoParticle.Vy(), neutrinoParticle.Vz()};
        geo::TPCID tpcID = theGeometry->FindTPCAtPosition(pos);
        tpcID_num = tpcID.TPC;

        if(tpcID.isValid){
            trueNeutrinoVX = neutrinoParticle.Vx();
            trueNeutrinoVY = neutrinoParticle.Vy();
            trueNeutrinoVZ = neutrinoParticle.Vz();
            trueCCNC = neutrino.CCNC();
            trueNeutrinoType = neutrinoParticle.PdgCode();
            trueChargedLepton = lepton.PdgCode();
            numTrueNeutrinos++;
        }
    }
  }

  art::Handle<std::vector<recob::Slice>>    sliceHandle;
  std::vector<art::Ptr<recob::Slice>>       sliceVec;

  art::Handle<std::vector<recob::Vertex>>   vertexHandle;
  std::vector<art::Ptr<recob::Vertex>>      vertexVec;

  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hitVec;

  if(e.getByLabel(hitLabel, hitHandle))
    art::fill_ptr_vector(hitVec, hitHandle);

  if(e.getByLabel(sliceLabel, sliceHandle))
    art::fill_ptr_vector(sliceVec, sliceHandle);

  if(sliceVec.empty()) return;


  int slice_ID(std::numeric_limits<int>::max());

  for(const art::Ptr<recob::Slice> &slice : sliceVec){
    nSlices++;
    slice_ID = slice->ID();

    std::vector<double> vec4SliceInfo = std::vector<double>(0);
    if(slice_ID == std::numeric_limits<int>::max()) return; 
 
    art::FindManyP<recob::Hit> sliceHitAssns(sliceVec, e, sliceLabel);
    const std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssns.at(slice.key())); 
    
    const int sliceID_truth = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hitVec, true);
    // Calculating slice completeness
    double sliceCompleteness = Completeness(e, sliceHits, sliceID_truth);
    std::cout << "Slice Completeness: " << sliceCompleteness << std::endl;
    vec4SliceInfo.push_back(slice_ID);
    vec4SliceInfo.push_back(sliceCompleteness);
    sliceInfo.push_back(vec4SliceInfo);
    vec4SliceInfo.clear();
  }

  if(e.getByLabel(PFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  if(pfpVec.empty()) return;

  if(e.getByLabel(trackLabel, trackHandle))
    art::fill_ptr_vector(trackVec, trackHandle);

  if(e.getByLabel(showerLabel, showerHandle))
    art::fill_ptr_vector(showerVec, showerHandle);

  nTracks = trackVec.size();
  nShowers = showerVec.size();

  for(const art::Ptr<recob::Shower> &shower : showerVec){
    art::FindManyP<recob::Hit> showerHitAssns(showerVec, e, showerLabel);
    const std::vector<art::Ptr<recob::Hit>> showerHits(showerHitAssns.at(shower.key()));
    
    //const int showerID_truth = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, true);

    // Calculating the completeness of the shower
    //double showerCompleteness = Completeness(e, showerHits, showerID_truth);
    //std::cout << "Shower Completeness: " << showerCompleteness << std::endl;

    // Calculating the purity of the shower
    //double showerPurity = Purity(e, showerHits, showerID_truth);
    //std::cout << "Shower Purity: " << showerPurity << std::endl;
  }

  std::cout << "num tracks: " << nTracks << " num showers: " << nShowers << std::endl;

  int PFParticleID(std::numeric_limits<int>::max());

  // Get associations between pfparticles and slices
  art::FindManyP<recob::Slice> pfpSliceAssns(pfpVec, e, sliceLabel);
  // Get associations between pfparticles and vertex
  art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);
  // Get associations between pfparticles and metadata
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpHandle, e, PFParticleLabel);
  // Get associations between pfparticles and showers
  art::FindOneP<recob::Shower> pfpShowerAssns(pfpVec, e, showerLabel);

  // initialises the vector that holds the primary reco neutrino info
  std::vector<std::tuple<art::Ptr<recob::PFParticle>, art::Ptr<recob::Vertex>, art::Ptr<recob::Slice>, art::Ptr<sbn::CRUMBSResult>>> v;

  for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
    nPFParticles++;
    PFParticleID = pfp->PdgCode();
  
    if(PFParticleID == std::numeric_limits<int>::max()) return;

    const std::vector<art::Ptr<recob::Slice>> pfpSlices(pfpSliceAssns.at(pfp.key()));
    
    // Condition on the number of slices and showers associated with the pfparticle. pfp has > 1 slice or no slice then it is skipped
    if(pfpSlices.size() == 1){
        const art::Ptr<recob::Slice> &pfpSlice(pfpSlices.front());
        //int SliceID = pfpSlice->ID();
    
        const auto meta  = pfpMetadataAssns.at(pfp.key());
        const auto props = meta->GetPropertiesMap();
        const auto trackscoreobj = props.find("TrackScore");
       
        if(!(PFParticleID == 12 || PFParticleID == 14)){
            trackscore = trackscoreobj->second;
        }

        const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));

        if(pfpVertexs.size() == 1){
            const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
            bool pfpPrimary = pfp->IsPrimary();

            // Print statements when looking at the cheated vertexing
            //std::cout << "PFP Vertex: " << pfpVertex->position().X() << ", " << pfpVertex->position().Y() << ", " << pfpVertex->position().Z() << std::endl;
            //std::cout << "Is Primary: " << pfpPrimary << " PFP ID: " << PFParticleID << std::endl;

            if(pfpPrimary && (PFParticleID == 12 || PFParticleID == 14)){
                art::FindManyP<sbn::CRUMBSResult> sliceCrumbsAssns(pfpSlices, e, crumbsLabel);
                // Initialise tuple to contain the pfp primary neutrino, its vertex, its slice, and the crumbs score associated to the slice
                std::tuple<art::Ptr<recob::PFParticle>, art::Ptr<recob::Vertex>, art::Ptr<recob::Slice>, art::Ptr<sbn::CRUMBSResult>> tuple;
                
                numRecoNeutrinos++;

                const std::vector<art::Ptr<sbn::CRUMBSResult>> sliceCrumbsResults = sliceCrumbsAssns.at(0);
                if(sliceCrumbsResults.size() != 1){
                    std::cout << "Errors: slice has multiple CRUMBS results" << std::endl;
                }

                const art::Ptr<sbn::CRUMBSResult> sliceCrumbsResult(sliceCrumbsResults.front());
                tuple = std::make_tuple(pfp, pfpVertex, pfpSlice, sliceCrumbsResult); // Creates a tuple that contains the pfp, and the vertex coords
                // Pushes the tuple containing the reco primary neutrino info to a vector called v
                v.push_back(tuple);
            }
        }
    } 
  }

  if(v.size() != 0){
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

    double chosenSliceID = std::get<2>(v[highestSliceScoreIndex])->ID();
    for(double sliceIndex = 0; sliceIndex < sliceInfo.size(); sliceIndex++){
        if(chosenSliceID == sliceInfo[sliceIndex][0]){
            chosenSliceCompleteness = sliceInfo[sliceIndex][1];
        } 
    }

  } else{
    recoNeutrinoVX = NAN;
    recoNeutrinoVY = NAN;
    recoNeutrinoVZ = NAN;
  }

  std::cout << "True vertex: (" << trueNeutrinoVX << ", " << trueNeutrinoVY << ", " << trueNeutrinoVZ << ")   Reco vertex: (" << recoNeutrinoVX << ", " << recoNeutrinoVY << ", " << recoNeutrinoVZ << ")" << std::endl;
  std::cout << "dx: " << recoNeutrinoVX - trueNeutrinoVX << ", dy: " << recoNeutrinoVY - trueNeutrinoVY << ", dz: " << recoNeutrinoVZ - trueNeutrinoVZ << std::endl; 
  std::cout << "Total shower energy: " << totalShowerEnergy << std::endl;

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
 
  std::cout << "completeness of chosen slice: " << chosenSliceCompleteness << std::endl;
  std::cout << "num slices: " << nSlices << ", num pfps: " << nPFParticles << std::endl;

  if(!std::isnan(recoNeutrinoVX)){
    numEventsTotal++;

    ShowerEnergy(e, pfpVec, recoNeutrinoVX, recoNeutrinoVY, recoNeutrinoVZ);
    std::cout << "total shower energy: " << totalShowerEnergy << ", theta: " << chosenShowerTheta << ", smallest theta: " << smallestTheta << std::endl;
   
    double largestThetaAllowed = std::sqrt(1.022/totalShowerEnergy); 
    double ETheta2 = totalShowerEnergy * std::pow(chosenShowerTheta, 2);
    std::cout << "largest theta value allowed: " << TMath::RadToDeg() * largestThetaAllowed << ", theta of shower: " << TMath::RadToDeg() * chosenShowerTheta << std::endl;
    if(largestThetaAllowed > chosenShowerTheta){
        std::cout << "This event would be selected!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        numEventsSelected++;
    }

    if(largestThetaAllowed > smallestTheta){
        std::cout << "event would be selected based on smallest theta!!!!!!!!!" << std::endl;
    }

    event_tree.push_back(eventID);
    run_tree.push_back(runID);
    subrun_tree.push_back(subRunID);
    trueNeutrinoVX_tree.push_back(trueNeutrinoVX);
    trueNeutrinoVY_tree.push_back(trueNeutrinoVY);
    trueNeutrinoVZ_tree.push_back(trueNeutrinoVZ);
    recoNeutrinoVX_tree.push_back(recoNeutrinoVX);
    recoNeutrinoVY_tree.push_back(recoNeutrinoVY);
    recoNeutrinoVZ_tree.push_back(recoNeutrinoVZ);
    trueCCNC_tree.push_back(trueCCNC);
    trueNeutrinoType_tree.push_back(trueNeutrinoType);
    trueLeptonType_tree.push_back(trueChargedLepton);
    numSlices_tree.push_back(nSlices);
    numPfps_tree.push_back(nPFParticles);
    numTracks_tree.push_back(nTracks);
    numShowers_tree.push_back(nShowers);
    tpcID_tree.push_back(tpcID_num);
    dl_current_tree.push_back(1); // 0 = uboone dl, 1 = dune dl, 2 = current
    sliceCompleteness_tree.push_back(chosenSliceCompleteness);
    trackScore_tree.push_back(trackscore);
    showerEnergy_tree.push_back(totalShowerEnergy);
    showerTheta_tree.push_back(chosenShowerTheta);
    showerSmallestTheta_tree.push_back(smallestTheta);
    showerETheta2_tree.push_back(ETheta2);
  } else{
    std::cout << "There is no reco neutrino, skipping event" << std::endl;
  }

 sliceInfo.clear(); 
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

void sbnd::NuE::ShowerEnergy(const art::Event &e, const std::vector<art::Ptr<recob::PFParticle>> &pfpVec, const double VX, const double VY, const double VZ){
  // Get associations between pfparticles and showers
  art::FindOneP<recob::Shower> pfpShowerAssns(pfpVec, e, showerLabel);
  art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpVec, e, PFParticleLabel);
  art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);

  smallestDeltaR = 1000000;
  smallestTheta = 1000000;
  std::cout << "reco vertex: " << VX << ", " << VY << ", " << VZ << std::endl;

  for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
    if(pfp->PdgCode() == std::numeric_limits<int>::max()) return;
  
    const art::Ptr<recob::Shower> pfpShower(pfpShowerAssns.at(pfp.key()));
    const auto meta  = pfpMetadataAssns.at(pfp.key());
    const auto props = meta->GetPropertiesMap();
    const auto trackscoreobj = props.find("TrackScore");
    const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
    const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
    
    if(!(pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
        std::cout << "2nd iteration PFP PDG: " << pfp.key() << ", Track Score: " << trackscoreobj->second << std::endl;
        if(trackscoreobj->second <= 1 && trackscoreobj->second >= 0){
            std::cout << "counting this as a shower" << std::endl;
            std::cout << "energy of this shower: " << pfpShower->Energy()[pfpShower->best_plane()] << ", best plane: " << pfpShower->best_plane() << std::endl;
            totalShowerEnergy += pfpShower->Energy()[pfpShower->best_plane()];
            std::cout << "theta: " << TMath::RadToDeg() * pfpShower->Direction().Theta() << std::endl;
            std::cout << "shower vertex: " << pfpVertex->position().X() << ", " << pfpVertex->position().Y() << ", " << pfpVertex->position().Z() << std::endl;
        
            double score = (std::sqrt(std::pow(VX - pfpVertex->position().X(), 2) + std::pow(VY - pfpVertex->position().Y(), 2) + std::pow(VZ - pfpVertex->position().Z(), 2)) + (std::pow(trackscoreobj->second, 2) * 100)); // + (std::pow(pfpShower->Energy()[pfpShower->best_plane()], -1) * 10000));
            std::cout << "score: " << score << std::endl;
            if(score < smallestDeltaR){ 
                chosenShowerTheta = pfpShower->Direction().Theta();
                smallestDeltaR = score;
            }

            if(pfpShower->Direction().Theta() < smallestTheta) smallestTheta = pfpShower->Direction().Theta();
        }
    }
  }
}

void sbnd::NuE::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  NuETree->Branch("event_tree", &event_tree);
  NuETree->Branch("run_tree", &run_tree);
  NuETree->Branch("subrun_tree", &subrun_tree);
  NuETree->Branch("trueNeutrinoVX_tree", &trueNeutrinoVX_tree);
  NuETree->Branch("trueNeutrinoVY_tree", &trueNeutrinoVY_tree);
  NuETree->Branch("trueNeutrinoVZ_tree", &trueNeutrinoVZ_tree);
  NuETree->Branch("recoNeutrinoVX_tree", &recoNeutrinoVX_tree);
  NuETree->Branch("recoNeutrinoVY_tree", &recoNeutrinoVY_tree);
  NuETree->Branch("recoNeutrinoVZ_tree", &recoNeutrinoVZ_tree);
  NuETree->Branch("trueCCNC_tree", &trueCCNC_tree);
  NuETree->Branch("trueNeutrinoType_tree", &trueNeutrinoType_tree);
  NuETree->Branch("trueLeptonType_tree", &trueLeptonType_tree);
  NuETree->Branch("numSlices_tree", &numSlices_tree);
  NuETree->Branch("numPfps_tree", &numPfps_tree);
  NuETree->Branch("numTracks_tree", &numTracks_tree);
  NuETree->Branch("numShowers_tree", &numShowers_tree);
  NuETree->Branch("tpcID_tree", &tpcID_tree);
  NuETree->Branch("dl_current_tree", &dl_current_tree);  
  NuETree->Branch("sliceCompleteness_tree", &sliceCompleteness_tree);
  NuETree->Branch("trackScore_tree", &trackScore_tree);
  NuETree->Branch("showerEnergy_tree", &showerEnergy_tree);
  NuETree->Branch("showerTheta_tree", &showerTheta_tree);
  NuETree->Branch("showerSmallestTheta_tree", &showerSmallestTheta_tree);
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

  std::cout << "events selected: " << numEventsSelected << "/" << numEventsTotal << std::endl;
}

DEFINE_ART_MODULE(sbnd::NuE)
