#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "sbndcode/RecoUtils/RecoUtils.h"
#include "lardataobj/MCBase/MCShower.h"

// sbndcode includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"

#include "lardataobj/RecoBase/OpHit.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include <map>

using namespace std;

namespace sbnd{
	
class SingleMuonInfo : public art::EDAnalyzer {
public:
    using CRTHit = sbn::crt::CRTHit;
    explicit SingleMuonInfo(fhicl::ParameterSet const& pset);
    virtual ~SingleMuonInfo();

    void beginJob();
    void analyze(const art::Event& evt);

 private:
    void ClearVecs();
		 
    TTree* fTruthTree; // Let's save both genie and G4 level information
    TTree* fRecHitTree; // Let's save reconstructed hit level informaion
    TTree* fRecTrackTree; // Let's save reconstructed track level information
    TTree* fRecShowerTree; // Let's save reconstructed shower information
    
    // Common Variables
    
    Int_t    frun;                  
    Int_t    fsubrun;               
    Int_t    fevent;
    
    art::InputTag fGenLabel;
    art::InputTag fSimLabel;
    art::InputTag fOpHitModuleLabel;
    art::InputTag fOpFlashModuleLabel0;
    art::InputTag fOpFlashModuleLabel1;
    art::InputTag fCrtHitModuleLabel;
    art::InputTag fCrtTrackModuleLabel;
    
    map<int,art::InputTag> fFlashLabels;
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
 }; // class SingleMuonInfo

//========================================================================
SingleMuonInfo::SingleMuonInfo(fhicl::ParameterSet const& pset) :
EDAnalyzer(pset),
fGenLabel(pset.get<art::InputTag>("GenLabel","generator")),
fSimLabel(pset.get<art::InputTag>("SimLabel","largeant")),
fOpHitModuleLabel(pset.get<art::InputTag>("OpHitModuleLabel","ophitpmt")),				
fOpFlashModuleLabel0(pset.get<art::InputTag>("OpFlashModuleLabel0","opflashtpc0")),
fOpFlashModuleLabel1(pset.get<art::InputTag>("OpFlashModuleLabel1","opflashtpc1")),
fCrtHitModuleLabel(pset.get<art::InputTag>("CrtHitModuleLabel","crthit")),
fCrtTrackModuleLabel(pset.get<art::InputTag>("CrtTrackModuleLabel","crttrack"))	
{
 fFlashLabels[0] = fOpFlashModuleLabel0;
 fFlashLabels[1] = fOpFlashModuleLabel1;	
}
 
//========================================================================
SingleMuonInfo::~SingleMuonInfo(){
  //destructor
}
//========================================================================

//========================================================================
void SingleMuonInfo::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  
  // Get a pointer to the geometry service provider.
  fGeometryService = lar::providerFrom<geo::Geometry>();  
  
  art::ServiceHandle<art::TFileService> tfs;
  
  fTruthTree = tfs->make<TTree>("TruthTree","");
  fTruthTree->Branch("run", &frun, "run/I");
  fTruthTree->Branch("subrun", &fsubrun, "subrun/I");
  fTruthTree->Branch("event", &fevent, "event/I");
  
  fRecHitTree = tfs->make<TTree>("RecHitTree","");
  fRecHitTree->Branch("run", &frun, "run/I");
  fRecHitTree->Branch("subrun", &fsubrun, "subrun/I");
  fRecHitTree->Branch("event", &fevent, "event/I");
  
  fRecTrackTree = tfs->make<TTree>("RecTrackTree","");
  fRecTrackTree->Branch("run", &frun, "run/I");
  fRecTrackTree->Branch("subrun", &fsubrun, "subrun/I");
  fRecTrackTree->Branch("event", &fevent, "event/I");
  
  fRecShowerTree = tfs->make<TTree>("fRecShowerTree","");
  fRecShowerTree->Branch("run", &frun, "run/I");
  fRecShowerTree->Branch("subrun", &fsubrun, "subrun/I");
  fRecShowerTree->Branch("event", &fevent, "event/I");
}

void SingleMuonInfo::analyze( const art::Event& evt){
     ClearVecs();
     frun = evt.run();
     fsubrun = evt.subRun();
     fevent = evt.id().event(); 
     
     fTruthTree->Fill();
     fRecHitTree->Fill();
     fRecTrackTree->Fill();
     fRecShowerTree->Fill();
} // end of analyze function
 
/////////////////////////////////////////// ClearVecs ///////////////////////////////
 
void SingleMuonInfo::ClearVecs()
{
     frun = -9999;
     fsubrun = -9999;
     fevent = -9999;
}
////////////////////////////////////////////////////////////////////////////////////
DEFINE_ART_MODULE(SingleMuonInfo)
}		
		
		
