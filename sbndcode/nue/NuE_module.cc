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

constexpr double def_double = -std::numeric_limits<double>::max();

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
  void clearVectors();
  void ClearMaps(const art::Event &e);
  void SetupMaps(const art::Event &e);
  double Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID);
  double Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID);
  void PFPs(art::Event const& e);
  int GetNumGenEvents(const art::Event &e); 
 
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  unsigned int eventID; // Event num
  unsigned int runID; // Run num  
  unsigned int subRunID; // Subrun num

  TTree* SubRunTree;
  double pot;
  int spills, numGenEvents;
  unsigned int subRunNumber;
  unsigned int subRunRun;


  TTree* NuETree;
  std::vector<double>   truth_recoilElectronPDG;
  std::vector<double>   truth_recoilElectronVX;
  std::vector<double>   truth_recoilElectronVY;
  std::vector<double>   truth_recoilElectronVZ;
  std::vector<double>   truth_recoilElectronPX;
  std::vector<double>   truth_recoilElectronPY;
  std::vector<double>   truth_recoilElectronPZ;
  std::vector<double>   truth_recoilElectronEnergy;
  std::vector<double>   truth_recoilElectronAngle;
  std::vector<double>   truth_recoilElectronETheta2;
  std::vector<double>   truth_recoilElectronDX;
  std::vector<double>   truth_recoilElectronDY;
  std::vector<double>   truth_recoilElectronDZ;
 
  std::vector<double>   truth_neutrinoVX;
  std::vector<double>   truth_neutrinoVY;
  std::vector<double>   truth_neutrinoVZ;
  std::vector<double>   truth_neutrinoCCNC;
  std::vector<double>   truth_neutrinoTPCID;
  std::vector<double>   truth_neutrinoTPCValid;
  std::vector<double>   truth_neutrinoType;

  std::vector<double>   reco_sliceID;
  std::vector<double>   reco_sliceCompleteness;
  std::vector<double>   reco_slicePurity;
  std::vector<double>   reco_sliceScore;
  std::vector<double>   reco_sliceCategory;
  std::vector<double>   reco_sliceInteraction;
  std::vector<double>   reco_sliceTrueVX;
  std::vector<double>   reco_sliceTrueVY;
  std::vector<double>   reco_sliceTrueVZ;
  std::vector<double>   reco_sliceNumHits;
  std::vector<double>   reco_sliceNumHitsTruthMatched;
  std::vector<double>   reco_sliceNumTruthHits;
  std::vector<double>   reco_sliceOrigin;
  std::vector<double>   reco_sliceTrueCCNC;
  std::vector<double>   reco_sliceTrueNeutrinoType;

  std::vector<double>   truth_particleSliceID;
  std::vector<double>   truth_particlePrimary;
  std::vector<double>   truth_particleVX;
  std::vector<double>   truth_particleVY;
  std::vector<double>   truth_particleVZ;
  std::vector<double>   truth_particlePDG;
  std::vector<double>   truth_particleTrackID;
  std::vector<double>   truth_particleMother;
  std::vector<double>   truth_particleStatusCode;

  std::vector<double>   reco_particlePDG;
  std::vector<double>   reco_particleIsPrimary;
  std::vector<double>   reco_particleVX;
  std::vector<double>   reco_particleVY;
  std::vector<double>   reco_particleVZ;
  std::vector<double>   reco_particleDX;
  std::vector<double>   reco_particleDY;
  std::vector<double>   reco_particleDZ;
  std::vector<double>   reco_particleSliceID;
  std::vector<double>   reco_particleBestPlaneEnergy;
  std::vector<double>   reco_particleTheta;
  std::vector<double>   reco_particleTrackScore;
  std::vector<double>   reco_particleCompleteness;
  std::vector<double>   reco_particlePurity;
  std::vector<double>   reco_particleID;
  std::vector<double>   reco_particleTruePDG;
  std::vector<double>   reco_particleTrueOrigin;
  std::vector<double>   reco_particleTrueInteractionType;
  std::vector<double>   reco_particleNumHits;
  std::vector<double>   reco_particleNumHitsTruthMatched;
  std::vector<double>   reco_particleNumTruthHits;
  std::vector<double>   reco_particleClearCosmic;
  
  std::vector<double>   reco_neutrinoID;
  std::vector<double>   reco_neutrinoPDG;
  std::vector<double>   reco_neutrinoVX;
  std::vector<double>   reco_neutrinoVY;
  std::vector<double>   reco_neutrinoVZ;
  std::vector<double>   reco_neutrinoSliceID;

  std::map<int,int> fHitsMap;

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

  NuETree = fs->make<TTree>("NuE","");
  NuETree->Branch("eventID", &eventID);
  NuETree->Branch("runID", &runID);
  NuETree->Branch("subRunID", &subRunID);
  NuETree->Branch("DLCurrent", &DLCurrent);
  NuETree->Branch("signal", &signal);

  NuETree->Branch("truth_neutrinoVX", &truth_neutrinoVX);
  NuETree->Branch("truth_neutrinoVY", &truth_neutrinoVY);
  NuETree->Branch("truth_neutrinoVZ", &truth_neutrinoVZ);
  NuETree->Branch("truth_neutrinoCCNC", &truth_neutrinoCCNC);
  NuETree->Branch("truth_neutrinoTPCID", &truth_neutrinoTPCID);
  NuETree->Branch("truth_neutrinoTPCValid", &truth_neutrinoTPCValid);
  NuETree->Branch("truth_neutrinoType", &truth_neutrinoType);
  
  NuETree->Branch("truth_recoilElectronPDG", &truth_recoilElectronPDG);
  NuETree->Branch("truth_recoilElectronVX", &truth_recoilElectronVX);
  NuETree->Branch("truth_recoilElectronVY", &truth_recoilElectronVY);
  NuETree->Branch("truth_recoilElectronVZ", &truth_recoilElectronVZ);
  NuETree->Branch("truth_recoilElectronPX", &truth_recoilElectronPX);
  NuETree->Branch("truth_recoilElectronPY", &truth_recoilElectronPY);
  NuETree->Branch("truth_recoilElectronPZ", &truth_recoilElectronPZ);
  NuETree->Branch("truth_recoilElectronEnergy", &truth_recoilElectronEnergy);
  NuETree->Branch("truth_recoilElectronAngle", &truth_recoilElectronAngle);
  NuETree->Branch("truth_recoilElectronETheta2", &truth_recoilElectronETheta2);
  NuETree->Branch("truth_recoilElectronDX", &truth_recoilElectronDX);
  NuETree->Branch("truth_recoilElectronDY", &truth_recoilElectronDY);
  NuETree->Branch("truth_recoilElectronDZ", &truth_recoilElectronDZ);

  NuETree->Branch("reco_sliceID", &reco_sliceID);
  NuETree->Branch("reco_sliceCompleteness", &reco_sliceCompleteness);
  NuETree->Branch("reco_slicePurity", &reco_slicePurity);
  NuETree->Branch("reco_sliceScore", &reco_sliceScore);
  NuETree->Branch("reco_sliceCategory", &reco_sliceCategory);
  NuETree->Branch("reco_sliceInteraction", &reco_sliceInteraction);
  NuETree->Branch("reco_sliceTrueVX", &reco_sliceTrueVX);
  NuETree->Branch("reco_sliceTrueVY", &reco_sliceTrueVY);
  NuETree->Branch("reco_sliceTrueVZ", &reco_sliceTrueVZ);
  NuETree->Branch("reco_sliceNumHits", &reco_sliceNumHits);
  NuETree->Branch("reco_sliceNumHitsTruthMatched", &reco_sliceNumHitsTruthMatched);
  NuETree->Branch("reco_sliceNumTruthHits", &reco_sliceNumTruthHits);
  NuETree->Branch("reco_sliceOrigin", &reco_sliceOrigin);
  NuETree->Branch("reco_sliceTrueCCNC", &reco_sliceTrueCCNC);
  NuETree->Branch("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);
  
  NuETree->Branch("truth_particleSliceID", &truth_particleSliceID);
  NuETree->Branch("truth_particlePrimary", &truth_particlePrimary);
  NuETree->Branch("truth_particleVX", &truth_particleVX);
  NuETree->Branch("truth_particleVY", &truth_particleVY);
  NuETree->Branch("truth_particleVZ", &truth_particleVZ);
  NuETree->Branch("truth_particlePDG", &truth_particlePDG);
  NuETree->Branch("truth_particleTrackID", &truth_particleTrackID);
  NuETree->Branch("truth_particleMother", &truth_particleMother);
  NuETree->Branch("truth_particleStatusCode", &truth_particleStatusCode);

  NuETree->Branch("reco_particlePDG", &reco_particlePDG);
  NuETree->Branch("reco_particleIsPrimary", &reco_particleIsPrimary);
  NuETree->Branch("reco_particleVX", &reco_particleVX);
  NuETree->Branch("reco_particleVY", &reco_particleVY);
  NuETree->Branch("reco_particleVZ", &reco_particleVZ);
  NuETree->Branch("reco_particleDX", &reco_particleDX);
  NuETree->Branch("reco_particleDY", &reco_particleDY);
  NuETree->Branch("reco_particleDZ", &reco_particleDZ);
  NuETree->Branch("reco_particleSliceID", &reco_particleSliceID);
  NuETree->Branch("reco_particleBestPlaneEnergy", &reco_particleBestPlaneEnergy);
  NuETree->Branch("reco_particleTheta", &reco_particleTheta);
  NuETree->Branch("reco_particleTrackScore", &reco_particleTrackScore);
  NuETree->Branch("reco_particleCompleteness", &reco_particleCompleteness);
  NuETree->Branch("reco_particlePurity", &reco_particlePurity);
  NuETree->Branch("reco_particleID", &reco_particleID);
  NuETree->Branch("reco_particleTruePDG", &reco_particleTruePDG);
  NuETree->Branch("reco_particleTrueOrigin", &reco_particleTrueOrigin);
  NuETree->Branch("reco_particleTrueInteractionType", &reco_particleTrueInteractionType);
  NuETree->Branch("reco_particleNumHits", &reco_particleNumHits);
  NuETree->Branch("reco_particleNumHitsTruthMatched", &reco_particleNumHitsTruthMatched);
  NuETree->Branch("reco_particleNumTruthHits", &reco_particleNumTruthHits);
  NuETree->Branch("reco_particleClearCosmic", &reco_particleClearCosmic);
  
  NuETree->Branch("reco_neutrinoID", &reco_neutrinoID);
  NuETree->Branch("reco_neutrinoPDG", &reco_neutrinoPDG);
  NuETree->Branch("reco_neutrinoVX", &reco_neutrinoVX);
  NuETree->Branch("reco_neutrinoVY", &reco_neutrinoVY);
  NuETree->Branch("reco_neutrinoVZ", &reco_neutrinoVZ);
  NuETree->Branch("reco_neutrinoSliceID", &reco_neutrinoSliceID);
 
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
    clearVectors();
    ClearMaps(e);
    SetupMaps(e);

    eventID = e.id().event();
    runID = e.id().run();
    subRunID = e.id().subRun();

    numGenEvents = GetNumGenEvents(e);

    std::cout << "" << std::endl;
    std::cout << "========================================================================================================" << std::endl; 
    std::cout << "Run: " << runID << ", Subrun: " << subRunID << ", Event: " << eventID << ", DL/Current: " << DLCurrent << std::endl;
    Slices(e);
    MCParticles(e);
    PFPs(e);

    NuETree->Fill();
}

int sbnd::NuE::GetNumGenEvents(const art::Event &e){
    int nGenEvents = 0;
    for(const art::ProcessConfiguration &process: e.processHistory()){
        std::optional<fhicl::ParameterSet> genConfig = e.getProcessParameterSet(process.processName());
        if (genConfig && genConfig->has_key("source") && genConfig->has_key("source.maxEvents") && genConfig->has_key("source.module_type") ) {
            int maxEvents = genConfig->get<int>("source.maxEvents");
            std::string moduleType = genConfig->get<std::string>("source.module_type");
            if (moduleType == "EmptyEvent") {
                nGenEvents += maxEvents;
            }
        }
    }
    printf("Num Gen Events = %i\n", nGenEvents);
    return nGenEvents;
}

double sbnd::NuE::Completeness(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    std::map<int, int> objectHitsMap;

    for(unsigned int i = 0; i < objectHits.size(); ++i)
        ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

    //std::cout << "Completeness calculator:" << std::endl;
    //std::cout << "Number of truth-matched hits in reco object: " << objectHitsMap[ID] << " Number of hits in object in truth: " << fHitsMap[ID] << std::endl;
    return (fHitsMap[ID] == 0) ? def_double : objectHitsMap[ID]/static_cast<double>(fHitsMap[ID]);
}

double sbnd::NuE::Purity(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &objectHits, const int ID){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    std::map<int, int> objectHitsMap;

    for(unsigned int i = 0; i < objectHits.size(); ++i)
        ++objectHitsMap[TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true)];

    double skippedHitsPFP = 0;
    for(unsigned int i = 0; i < objectHits.size(); ++i){
        double trueID_particle = TruthMatchUtils::TrueParticleID(clockData,objectHits[i],true);
        if(trueID_particle != ID){
            skippedHitsPFP++;
        }
    }

    /*
    std::cout << "Purity calculator:" << std::endl;
    std::cout << "Total number of PFP hits = " << static_cast<double>(objectHits.size()) << std::endl;
    std::cout << "Total number of truth matched hits in PFP = " << objectHitsMap[ID] << std::endl;
    std::cout << "Number of hits not truth matched to PFP = " << skippedHitsPFP << std::endl;
    */
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

void sbnd::NuE::PFPs(art::Event const& e){
    std::cout << "_________ PFParticles _________" << std::endl;
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    art::Handle<std::vector<recob::Shower>> showerHandle;
    std::vector<art::Ptr<recob::Shower>> showerVec;
    if(e.getByLabel(showerLabel, showerHandle))
        art::fill_ptr_vector(showerVec, showerHandle);

    int counter = 0;
    int neutrinoCounter = 0;

    if(!pfpVec.empty()){
        art::FindOneP<recob::Shower> pfpShowerAssns(pfpVec, e, showerLabel);
        art::FindOneP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(pfpVec, e, PFParticleLabel);
        art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);
        art::FindOneP<recob::Slice> pfpSliceAssns(pfpVec, e, sliceSCELabel);    
        art::FindManyP<recob::Hit> showerHitAssns(showerVec, e, showerLabel);
                
        std::vector<double> recoNeutrinoIDs;        
        std::vector<double> recoNeutrinoSliceIDs;

        for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
            if(!(pfp->PdgCode() == std::numeric_limits<int>::max())){
                const art::Ptr<recob::Shower> pfpShower(pfpShowerAssns.at(pfp.key()));
                const art::Ptr<recob::Slice> pfpSlice(pfpSliceAssns.at(pfp.key()));
                const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
                
                if(pfpVertexs.size() > 0){
                    if(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
                        // This is the reco neutrino
                        // Add the reco neutrino to the vectors
                        recoNeutrinoIDs.push_back(pfp->Self());
                        recoNeutrinoSliceIDs.push_back(pfpSlice->ID());
                        std::cout << "Reco Neutrino found in slice " << pfpSlice->ID() << ", with ID = " << pfp->Self() << std::endl;
                    }
                }
            }
        }

        for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
            if(!(pfp->PdgCode() == std::numeric_limits<int>::max())){
                std::cout << "----" << std::endl;
                const art::Ptr<recob::Shower> pfpShower(pfpShowerAssns.at(pfp.key()));
                const art::Ptr<recob::Slice> pfpSlice(pfpSliceAssns.at(pfp.key()));
                const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
                
                const auto meta  = pfpMetadataAssns.at(pfp.key());
                const auto props = meta->GetPropertiesMap();
                const auto trackscoreobj = props.find("TrackScore");
                const auto isClearCosmicObj = props.find("IsClearCosmic");
                           
                /* 
                if(isClearCosmicObj != props.end()){
                    std::cout << "isClearCosmic = " << isClearCosmicObj->second << std::endl;
                } else{
                    std::cout << "isClearCosmic = props.end()" << std::endl;
                }
                */

                if(pfpVertexs.size() > 0){
                    if(!(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14))){
                        // PFP is not the reco neutrino
                        //if(trackscoreobj->second <= 1 && trackscoreobj->second >= 0){
                            const std::vector<art::Ptr<recob::Hit>> showerHits(showerHitAssns.at(pfpShower.key()));
                            const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
                            counter++;

                            const int showerID_truth = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, showerHits, true);
                            double pfpCompleteness = Completeness(e, showerHits, showerID_truth);
                            double pfpPurity = Purity(e, showerHits, showerID_truth);

                            const simb::MCParticle* pfpMCParticle = particleInv->TrackIdToParticle_P(showerID_truth);

                            art::Handle<std::vector<simb::MCParticle>> truthParticleHandle;
                            std::vector<art::Ptr<simb::MCParticle>> truthParticleVec;
                            if(e.getByLabel(MCTruthLabel, truthParticleHandle))
                                art::fill_ptr_vector(truthParticleVec, truthParticleHandle);

                            double primary = 0;
                            std::cout << "pfp->IsPrimary() = " << pfp->IsPrimary() << std::endl;
                            if(pfp->IsPrimary() == 1){
                                // This is a cosmic -> particle is primary
                                primary = 1;
                            } else{
                                // This could still be a direct daughter of the neutrino = primary
                                // Get the mother of the particle, if it equals the ID of the reco neutrino in the slice -> particle is primary
                                for(std::size_t i = 0; i < recoNeutrinoIDs.size(); ++i){
                                    if(recoNeutrinoSliceIDs[i] == pfpSlice->ID()){
                                        // PFP is in the same slice as this reco neutrino
                                        std::cout << "reco neutrino ID = " << recoNeutrinoIDs[i] << ", parent ID = " << pfp->Parent() << std::endl;
                                        std::cout << "matched slice" << std::endl;
                                        if(recoNeutrinoIDs[i] == pfp->Parent()){
                                            // The PFPs parent is the reco neutrino -> particle is primary
                                            primary = 1;
                                        }
                                    }
                                }
                            }
                            std::cout << "is PFP Primary = " << primary << std::endl;

                            reco_particlePDG.push_back(pfp->PdgCode());
                            reco_particleIsPrimary.push_back(primary);
                            reco_particleVX.push_back(pfpVertex->position().X());
                            reco_particleVY.push_back(pfpVertex->position().Y());
                            reco_particleVZ.push_back(pfpVertex->position().Z());
                            reco_particleDX.push_back(pfpShower->Direction().X());
                            reco_particleDY.push_back(pfpShower->Direction().Y());
                            reco_particleDZ.push_back(pfpShower->Direction().Z());
                            reco_particleSliceID.push_back(pfpSlice->ID());
                            reco_particleBestPlaneEnergy.push_back(pfpShower->Energy()[pfpShower->best_plane()]);
                            reco_particleTheta.push_back(pfpShower->Direction().Theta());
                            reco_particleCompleteness.push_back(pfpCompleteness);
                            reco_particlePurity.push_back(pfpPurity);
                            reco_particleID.push_back(pfp->Self());
                            reco_particleNumHits.push_back(showerHits.size());
                            reco_particleNumHitsTruthMatched.push_back(pfpPurity*showerHits.size());
                            reco_particleNumTruthHits.push_back(fHitsMap[showerID_truth]);

                           
                            if(trackscoreobj != props.end()){
                                reco_particleTrackScore.push_back(trackscoreobj->second);
                                if(trackscoreobj->second == 0) std::cout << "TRACKSCORE = 0" << std::endl;
                                else std::cout << "Trackscore = " << trackscoreobj->second << std::endl;
                            } else{
                                std::cout << "Trackscore == props.end()" << std::endl;
                                reco_particleTrackScore.push_back(-999999);
                            }
                           
                            if(isClearCosmicObj != props.end()){
                                std::cout << "isClearCosmic = " << isClearCosmicObj->second << std::endl;
                                reco_particleClearCosmic.push_back(isClearCosmicObj->second);
                            } else{
                                std::cout << "isClearCosmic = props.end()" << std::endl;
                                reco_particleClearCosmic.push_back(0);
                            }

                            if(isClearCosmicObj != props.end() && trackscoreobj != props.end()) std::cout << "PFP has Trackscore and ClearCosmic Score!! Trackscore = " << trackscoreobj->second << ", ClearCosmic = " << isClearCosmicObj->second << std::endl;
                            
                            //printf("Reco Particle %d: ID = %li, PDG Code = %d, Is Primary = %d, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %d, Best Plane Energy = %f, Theta = %f, Trackscore = %f, Completeness = %f, Purity = %f\n", counter, pfp->Self(), pfp->PdgCode(), pfp->IsPrimary(), pfpVertex->position().X(), pfpVertex->position().Y(), pfpVertex->position().Z(), pfpShower->Direction().X(), pfpShower->Direction().Y(), pfpShower->Direction().Z(), pfpSlice->ID(), pfpShower->Energy()[pfpShower->best_plane()], pfpShower->Direction().Theta(), trackscoreobj->second, pfpCompleteness, pfpPurity);
                            printf("Reco Particle %d: ID = %li, PDG Code = %d, Is Primary = %d, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %d, Best Plane Energy = %f, Theta = %f, Completeness = %f, Purity = %f, Num Hits in PFP = %zu\n", counter, pfp->Self(), pfp->PdgCode(), pfp->IsPrimary(), pfpVertex->position().X(), pfpVertex->position().Y(), pfpVertex->position().Z(), pfpShower->Direction().X(), pfpShower->Direction().Y(), pfpShower->Direction().Z(), pfpSlice->ID(), pfpShower->Energy()[pfpShower->best_plane()], pfpShower->Direction().Theta(), pfpCompleteness, pfpPurity, showerHits.size());
                            
                            if(pfpMCParticle){
                                // The PFP has truth info associated with it
                            
                                const art::Ptr<simb::MCTruth> pfpMCTruth = particleInv->TrackIdToMCTruth_P(showerID_truth);
                                const simb::MCNeutrino pfpMCNeutrino = pfpMCTruth->GetNeutrino();
                                
                                reco_particleTruePDG.push_back(pfpMCParticle->PdgCode());
                                reco_particleTrueOrigin.push_back(pfpMCTruth->Origin());
                                if(pfpMCTruth->Origin() == 1){
                                    // From a beam neutrino
                                    reco_particleTrueInteractionType.push_back(pfpMCNeutrino.InteractionType());
                                } else{
                                    // From a cosmic ray
                                    reco_particleTrueInteractionType.push_back(-999999);
                                }
                                
                                
                                std::cout << "Truth info associated with PFP: Interaction Type = " << pfpMCNeutrino.InteractionType() << ", Origin of Neutrino = " << pfpMCTruth->Origin() << ", True PDG Code = " << pfpMCParticle->PdgCode() << std::endl;
                            } else{
                                // No truth info associated with the PFP
                                reco_particleTruePDG.push_back(-999999);
                                reco_particleTrueOrigin.push_back(-999999);
                                reco_particleTrueInteractionType.push_back(-999999);
                                std::cout << "PFP has no truth info associated with it, showerID_truth = " << showerID_truth << std::endl;
                            }
                        //}
                    } else{
                        // This is the reco neutrino
                        neutrinoCounter++;
                        const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
                        printf("Reco Neutrino %d: ID = %li, PDG Code = %d, Is Primary = %d, Vertex = (%f, %f, %f), Slice ID = %d\n", counter, pfp->Self(), pfp->PdgCode(), pfp->IsPrimary(), pfpVertex->position().X(), pfpVertex->position().Y(), pfpVertex->position().Z(), pfpSlice->ID());
                    
                        reco_neutrinoID.push_back(pfp->Self());
                        reco_neutrinoPDG.push_back(pfp->PdgCode());
                        reco_neutrinoVX.push_back(pfpVertex->position().X());
                        reco_neutrinoVY.push_back(pfpVertex->position().Y());
                        reco_neutrinoVZ.push_back(pfpVertex->position().Z());
                        reco_neutrinoSliceID.push_back(pfpSlice->ID());
                    }
                }
            }
        } 
    }

    if(counter == 0){
        // There are no PFPs in the event
        reco_particlePDG.push_back(-999999);
        reco_particleIsPrimary.push_back(-999999);
        reco_particleVX.push_back(-999999);
        reco_particleVY.push_back(-999999);
        reco_particleVZ.push_back(-999999);
        reco_particleDX.push_back(-999999);
        reco_particleDY.push_back(-999999);
        reco_particleDZ.push_back(-999999);
        reco_particleSliceID.push_back(-999999);
        reco_particleBestPlaneEnergy.push_back(-999999);
        reco_particleTheta.push_back(-999999);
        reco_particleTrackScore.push_back(-999999);
        reco_particleCompleteness.push_back(-999999);
        reco_particlePurity.push_back(-999999);
        reco_particleID.push_back(-999999);
        reco_particleTruePDG.push_back(-999999);
        reco_particleTrueOrigin.push_back(-999999);
        reco_particleTrueInteractionType.push_back(-999999);
        reco_particleNumHits.push_back(-999999);
        reco_particleNumHitsTruthMatched.push_back(-999999);
        reco_particleNumTruthHits.push_back(-999999);
        reco_particleClearCosmic.push_back(-999999);
    }

    if(neutrinoCounter == 0){
        reco_neutrinoID.push_back(-999999);
        reco_neutrinoPDG.push_back(-999999);
        reco_neutrinoVX.push_back(-999999);
        reco_neutrinoVY.push_back(-999999);
        reco_neutrinoVZ.push_back(-999999);
        reco_neutrinoSliceID.push_back(-999999);
    }
    std::cout << "_________________________" << std::endl;
}

void sbnd::NuE::MCParticles(art::Event const& e){
    art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
    if(e.getByLabel(TruthLabel, MCTruthHandle))
        art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

    std::cout << "_________ MCTruth _________" << std::endl;
    int counter = 0;
    if(!MCTruthVec.empty()){
        for(auto &MCTruth : MCTruthVec){
            std::cout << "MCTruth: Origin = " << MCTruth->Origin() << ", Number of Particles = " << MCTruth->NParticles() << std::endl;
            if(MCTruth->Origin() == simb::kBeamNeutrino){
                simb::MCNeutrino neutrino = MCTruth->GetNeutrino();
                simb::MCParticle neutrinoParticle = neutrino.Nu();
                std::cout << "Neutrino: Track ID = " << neutrinoParticle.TrackId() << ", CCNC = " << neutrino.CCNC() << ", Interaction Type = " << neutrino.InteractionType() << "Vertex = (" << neutrinoParticle.Vx() << ", " << neutrinoParticle.Vy() << ", " << neutrinoParticle.Vz() << "), Number of Daughters = " << neutrinoParticle.NumberDaughters() << std::endl; 
                int neutrinoTrackID = neutrinoParticle.TrackId();

                geo::Point_t pos = {neutrinoParticle.Vx(), neutrinoParticle.Vy(), neutrinoParticle.Vz()};
                geo::TPCID tpcID = theGeometry->FindTPCAtPosition(pos);

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
                                
                                double DX = (particle.Px() / P_mag);
                                double DY = (particle.Py() / P_mag);
                                double DZ = (particle.Pz() / P_mag);

                                std::cout << "MCParticle: Track ID = " << particle.TrackId() << ", Vertex = (" << particle.Vx() << ", " << particle.Vy() << ", " << particle.Vz() << "), PDG = " << particle.PdgCode() << ", Mother = " << particle.Mother() << ", Status Code = " <<  particle.StatusCode() << ", Energy = " << energy << ", Theta = " << theta << ", EThetaSquared = " << EThetaSquared << std::endl;

                                counter++;
                                truth_recoilElectronPDG.push_back(particle.PdgCode());
                                truth_recoilElectronVX.push_back(particle.Vx());
                                truth_recoilElectronVY.push_back(particle.Vy());
                                truth_recoilElectronVZ.push_back(particle.Vz());
                                truth_recoilElectronPX.push_back(particle.Px());
                                truth_recoilElectronPY.push_back(particle.Py());
                                truth_recoilElectronPZ.push_back(particle.Pz());
                                truth_recoilElectronEnergy.push_back(energy);
                                truth_recoilElectronAngle.push_back(theta);
                                truth_recoilElectronETheta2.push_back(EThetaSquared);
                                truth_recoilElectronDX.push_back(DX);
                                truth_recoilElectronDY.push_back(DY);
                                truth_recoilElectronDZ.push_back(DZ);
                
                                truth_neutrinoVX.push_back(neutrinoParticle.Vx());
                                truth_neutrinoVY.push_back(neutrinoParticle.Vy());
                                truth_neutrinoVZ.push_back(neutrinoParticle.Vz());
                                truth_neutrinoCCNC.push_back(neutrino.CCNC());
                                truth_neutrinoTPCID.push_back(tpcID.TPC);
                                truth_neutrinoTPCValid.push_back(static_cast<double>(tpcID.isValid));
                                truth_neutrinoType.push_back(neutrinoParticle.PdgCode());
                            }

                        }
                    }
                } else{
                    // This is a BNB or Cosmic event
                }
            }
        }
    }
    
    // There is no true recoil electron from a Nu+E elastic scattering interaction in the event. 
    if(counter == 0){
        truth_recoilElectronPDG.push_back(-999999);
        truth_recoilElectronVX.push_back(-999999);
        truth_recoilElectronVY.push_back(-999999);
        truth_recoilElectronVZ.push_back(-999999);
        truth_recoilElectronPX.push_back(-999999);
        truth_recoilElectronPY.push_back(-999999);
        truth_recoilElectronPZ.push_back(-999999);
        truth_recoilElectronEnergy.push_back(-999999);
        truth_recoilElectronAngle.push_back(-999999);
        truth_recoilElectronETheta2.push_back(-999999);
        truth_recoilElectronDX.push_back(-999999);
        truth_recoilElectronDY.push_back(-999999);
        truth_recoilElectronDZ.push_back(-999999);
    
        truth_neutrinoVX.push_back(-999999);
        truth_neutrinoVY.push_back(-999999);
        truth_neutrinoVZ.push_back(-999999);
        truth_neutrinoCCNC.push_back(-999999);
        truth_neutrinoTPCID.push_back(-999999);
        truth_neutrinoTPCValid.push_back(-999999);
        truth_neutrinoType.push_back(-999999);
    }

    std::cout << "_________________________" << std::endl;

}

void sbnd::NuE::Slices(art::Event const& e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    std::cout << "_________ Slices _________" << std::endl;

    int counter = 0;
    int counter_trueParticles = 0;

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

                // Getting the CRUMBS score of the slice
                art::FindManyP<sbn::CRUMBSResult> sliceCrumbsAssns(sliceVec, e, crumbsLabel);
                const std::vector<art::Ptr<sbn::CRUMBSResult>> sliceCrumbsResults = sliceCrumbsAssns.at(slice.key()); 

                double sliceScoreVar = 0;

                if(sliceCrumbsResults.size() == 1){
                    const art::Ptr<sbn::CRUMBSResult> sliceCrumbsResult(sliceCrumbsResults.front());
                    sliceScoreVar = sliceCrumbsResult->score;
                } else{
                    sliceScoreVar = -999999;
                }

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

                int numTrueParticlesSlice = sliceMCTruth->NParticles();
                std::cout << "Number of particles in the MC Truth matched to the slice = " << numTrueParticlesSlice << std::endl;
                std::cout << "True Neutrino Track ID = " << sliceMCNeutrinoParticle.TrackId() << std::endl;
                if(!(sliceMCTruth->Origin() == simb::kCosmicRay || sliceMCTruth->Origin() == 0)){
                    for(int b = 0; b < numTrueParticlesSlice; b++){
                        // Loop through all of the particles in the MC Truth
                        counter_trueParticles++;
                        const simb::MCParticle truthParticle_slice = sliceMCTruth->GetParticle(b);
                        int primaryParticle = 0;
                        if(sliceMCNeutrinoParticle.TrackId() == truthParticle_slice.Mother()){
                            std::cout << "This particle has its mother as the neutrino!" << std::endl; 
                            primaryParticle = 1;
                        }
                        std::cout << "Truth Particle " << b << ": Slice ID = " << sliceID << ", Primary = " << primaryParticle << ", Vertex = (" << truthParticle_slice.Vx() << ", " << truthParticle_slice.Vy() << ", " << truthParticle_slice.Vz() << "), PDG = " << truthParticle_slice.PdgCode() << ", Track ID = " << truthParticle_slice.TrackId() << ", Mother Track ID = " << truthParticle_slice.Mother() << ", Status Code = " << truthParticle_slice.StatusCode() << std::endl;
                    
                        truth_particleSliceID.push_back(sliceID);
                        truth_particlePrimary.push_back(primaryParticle);
                        truth_particleVX.push_back(truthParticle_slice.Vx());
                        truth_particleVY.push_back(truthParticle_slice.Vy());
                        truth_particleVZ.push_back(truthParticle_slice.Vz());
                        truth_particlePDG.push_back(truthParticle_slice.PdgCode());
                        truth_particleTrackID.push_back(truthParticle_slice.TrackId());
                        truth_particleMother.push_back(truthParticle_slice.Mother());
                        truth_particleStatusCode.push_back(truthParticle_slice.StatusCode());
                    }
                }else{
                    std::cout << "slice is truth matched to cosmic origin, not looking at the true particles" << std::endl;
                }

                std::cout << "MCParticle that owns most hits in slice: PDG = " << sliceMCParticle->PdgCode() << std::endl;
                std::cout << "MCTruth: Origin = " << sliceMCTruth->Origin() << ", Num Particles = " << sliceMCTruth->NParticles() << std::endl;
                
                if(sliceMCTruth->Origin() == simb::kBeamNeutrino){
                    // This is a beam neutrino
                    std::cout << "MCNeutrino: Interaction Type = " << sliceMCNeutrino.InteractionType() << ", CCNC = " << sliceMCNeutrino.CCNC() << ", PDG = " << sliceMCNeutrinoParticle.PdgCode() << ", Vertex = (" << sliceMCNeutrinoParticle.Vx() << ", " << sliceMCNeutrinoParticle.Vy() << ", " << sliceMCNeutrinoParticle.Vz() << ")" << std::endl;
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

                if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() == 1098) std::cout << "SLICE MCTRUTH = NU+E SCATTER, LOOK HERE!" << std::endl;

                // Looping through the hits in the slice
                for(const art::Ptr<recob::Hit> &sliceHit : sliceHits){
                    // Gets the true particle ID of the truth particle which owns the hit. True is for rollup.
                    const int sliceHitTrueParticleID = TruthMatchUtils::TrueParticleID(clockData, sliceHit, true);
                    //std::cout << "The slice hit comes from true particle with ID = " << sliceHitTrueParticleID << std::endl;
                    
                    if(sliceHitTrueParticleID == std::numeric_limits<int>::min()){
                        //std::cout << "sliceHitTrueParticleID = " << sliceHitTrueParticleID << std::endl;
                        //geo::WireID wid = sliceHit->WireID();
                        //std::cout << "Peak Time = " << sliceHit->PeakTime() << ", Wire ID: Cryostat = " << wid.Cryostat << ", TPC = " << wid.TPC << ", Plane = " << wid.Plane << ", Wire = " << wid.Wire << std::endl;
                        continue;
                    }
                    // Gets the MCParticle that the hit comes from
                    //const simb::MCParticle* sliceHitMCParticle = particleInv->TrackIdToParticle_P(sliceHitTrueParticleID);
                    
                    // Gets the MCTruth that the hit comes from   
                    const art::Ptr<simb::MCTruth> sliceHitMCTruth = particleInv->TrackIdToMCTruth_P(sliceHitTrueParticleID);

                    // If the MCTruth from the hit is the same as the MCTruth of the slice, add to the counter
                    if(sliceMCTruth == sliceHitMCTruth){
                        numTruthMatchedHitsSlice++;
                    } else{
                        const simb::MCNeutrino sliceHitMCNeutrino = sliceHitMCTruth->GetNeutrino();
                        //std::cout << "This hit is in the slice but it isn't truth matched to the slice:" << std::endl;
                        //std::cout << "sliceHitMCTruth: Origin = " << sliceHitMCTruth->Origin() << ", Interaction Type = " << sliceHitMCNeutrino.InteractionType() << std::endl;
                    }
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
                if(sliceMCTruth->Origin() == simb::kCosmicRay || sliceMCTruth->Origin() == 0) sliceCategory = 0;
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() == 1098 && sliceCompleteness > 0.5 && slicePurity > 0.3) sliceCategory = 1; 
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() == 1098 && sliceCompleteness < 0.5 && sliceCompleteness > 0.1) sliceCategory = 2;
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() != 1098 && sliceCompleteness > 0.5 && slicePurity > 0.3) sliceCategory = 3;
                else if(sliceMCTruth->Origin() == simb::kBeamNeutrino && sliceMCNeutrino.InteractionType() != 1098 && sliceCompleteness < 0.5 && sliceCompleteness > 0.1) sliceCategory = 4;

                std::cout << "Slice Category = " << sliceCategory << std::endl;

                counter++;
                reco_sliceID.push_back(sliceID);        
                reco_sliceCompleteness.push_back(sliceCompleteness);        
                reco_slicePurity.push_back(slicePurity);        
                reco_sliceScore.push_back(sliceScoreVar);        
                reco_sliceCategory.push_back(sliceCategory);       
                reco_sliceNumHits.push_back(numHitsSlice);
                reco_sliceNumHitsTruthMatched.push_back(numTruthMatchedHitsSlice);
                reco_sliceNumTruthHits.push_back(numTruthMatchedHits);

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
                    std::cout << "True Neutrino Vertex = (" << sliceMCNeutrinoParticle.Vx() << ", " << sliceMCNeutrinoParticle.Vy() << ", " << sliceMCNeutrinoParticle.Vz() << "), Type of Interaction = " << sliceMCNeutrino.CCNC() << ", PDG Code of Neutrino = " << sliceMCNeutrinoParticle.PdgCode() << std::endl; 
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
        reco_sliceCompleteness.push_back(-999999);        
        reco_slicePurity.push_back(-999999);        
        reco_sliceScore.push_back(-999999);        
        reco_sliceCategory.push_back(-999999);        
        reco_sliceInteraction.push_back(-999999);       
        reco_sliceTrueVX.push_back(-999999); 
        reco_sliceTrueVY.push_back(-999999); 
        reco_sliceTrueVZ.push_back(-999999); 
        reco_sliceNumHits.push_back(-999999);
        reco_sliceNumHitsTruthMatched.push_back(-999999);
        reco_sliceNumTruthHits.push_back(-999999);
        reco_sliceOrigin.push_back(-999999);
        reco_sliceTrueCCNC.push_back(-999999);
        reco_sliceTrueNeutrinoType.push_back(-999999);
    }
   
    if(counter_trueParticles == 0){
        // There are no truth particles in the any of the slices in the event
        truth_particleSliceID.push_back(-999999);
        truth_particlePrimary.push_back(-999999);
        truth_particleVX.push_back(-999999);
        truth_particleVY.push_back(-999999);
        truth_particleVZ.push_back(-999999);
        truth_particlePDG.push_back(-999999);
        truth_particleTrackID.push_back(-999999);
        truth_particleMother.push_back(-999999);
        truth_particleStatusCode.push_back(-999999);
    }

    std::cout << "_________________________" << std::endl;

}

void sbnd::NuE::clearVectors(){
    truth_recoilElectronPDG.clear();
    truth_recoilElectronVX.clear();
    truth_recoilElectronVY.clear();
    truth_recoilElectronVZ.clear();
    truth_recoilElectronPX.clear();
    truth_recoilElectronPY.clear();
    truth_recoilElectronPZ.clear();
    truth_recoilElectronEnergy.clear();
    truth_recoilElectronAngle.clear();
    truth_recoilElectronETheta2.clear();
    truth_recoilElectronDX.clear();
    truth_recoilElectronDY.clear();
    truth_recoilElectronDZ.clear();

    truth_neutrinoVX.clear();
    truth_neutrinoVY.clear();
    truth_neutrinoVZ.clear();
    truth_neutrinoCCNC.clear();
    truth_neutrinoTPCID.clear();
    truth_neutrinoTPCValid.clear();
    truth_neutrinoType.clear();

    reco_sliceID.clear();
    reco_sliceCompleteness.clear();
    reco_slicePurity.clear();
    reco_sliceScore.clear();
    reco_sliceCategory.clear();
    reco_sliceInteraction.clear();
    reco_sliceTrueVX.clear();
    reco_sliceTrueVY.clear();
    reco_sliceTrueVZ.clear();
    reco_sliceNumHits.clear();
    reco_sliceNumHitsTruthMatched.clear();
    reco_sliceNumTruthHits.clear();
    reco_sliceOrigin.clear();
    reco_sliceTrueCCNC.clear();
    reco_sliceTrueNeutrinoType.clear();
    
    truth_particleSliceID.clear();
    truth_particlePrimary.clear();
    truth_particleVX.clear();
    truth_particleVY.clear();
    truth_particleVZ.clear();
    truth_particlePDG.clear();
    truth_particleTrackID.clear();
    truth_particleMother.clear();
    truth_particleStatusCode.clear();

    reco_particlePDG.clear();
    reco_particleIsPrimary.clear();
    reco_particleVX.clear();
    reco_particleVY.clear();
    reco_particleVZ.clear();
    reco_particleDX.clear();
    reco_particleDY.clear();
    reco_particleDZ.clear();
    reco_particleSliceID.clear();
    reco_particleBestPlaneEnergy.clear();
    reco_particleTheta.clear();
    reco_particleTrackScore.clear();
    reco_particleCompleteness.clear();
    reco_particlePurity.clear();
    reco_particleID.clear();
    reco_particleTruePDG.clear();
    reco_particleTrueOrigin.clear();
    reco_particleTrueInteractionType.clear();
    reco_particleNumHits.clear();
    reco_particleNumHitsTruthMatched.clear();
    reco_particleNumTruthHits.clear();
    reco_particleClearCosmic.clear();
    
    reco_neutrinoID.clear();
    reco_neutrinoPDG.clear();
    reco_neutrinoVX.clear();
    reco_neutrinoVY.clear();
    reco_neutrinoVZ.clear();
    reco_neutrinoSliceID.clear();

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
