////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Tue Apr  2 04:53:46 2024 by Jorge Romeo araujo using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// SBN/SBND includes
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbnobj/SBND/Commissioning/MuonTrack.hh"
#include "sbnobj/SBND/Trigger/pmtTrigger.hh"
// #include "sbnobj/SBND/Trigger/pmtSoftwareTrigger.hh"
// #include "sbnobj/SBND/Trigger/CRTmetric.hh"

// Truth includes
//#include "larsim/MCCheater/BackTrackerService.h"
//#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h" 
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

const int MAX_INT = std::numeric_limits<int>::max();
const long int TIME_CORRECTION = (long int) std::numeric_limits<int>::max() * 2;
const int DEFAULT_VALUE = -9999;

namespace test {
  class AnalyzeEvents;
}


class test::AnalyzeEvents : public art::EDAnalyzer {
public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
 // void endJob() override;

private:
  void ResetVars();
  void ResizeMCParticle(int nParticles);
  /// Resize the data structure for MCTracks
  void ResizeMCTrack(int nTracks);
  /// Resize the data structure for MCShowers
  void ResizeMCShower(int nShowers);


  // Declare member data here.

  TTree* fTree;
    //run information
  int _run;        ///< The run number
  int _subrun;     ///< The subrun number
  int _event;      ///< The event number
  double _evttime; ///< The event time
  int _t0;         ///< The t0

  unsigned int fEventID;
  
   //MCParticle Info
  size_t MaxMCParticles = 0;
  Int_t     mcpart_no_primaries;                 
  std::vector<Int_t>    mcpart_pdg;              
  std::vector<Int_t>    mcpart_status;           
  std::vector<std::string>    mcpart_process;
  std::vector<std::string>    mcpart_endprocess;
  std::vector<Float_t>  mcpart_Eng;              
  std::vector<Float_t>  mcpart_EndE;
  std::vector<Float_t>  mcpart_Mass;
  std::vector<Float_t>  mcpart_Px;
  std::vector<Float_t>  mcpart_Py;
  std::vector<Float_t>  mcpart_Pz;
  std::vector<Float_t>  mcpart_P;
  std::vector<Float_t>  mcpart_StartPointx;
  std::vector<Float_t>  mcpart_StartPointy;
  std::vector<Float_t>  mcpart_StartPointz;
  std::vector<Float_t>  mcpart_StartT;  
  std::vector<Float_t>  mcpart_EndT;          
  std::vector<Float_t>  mcpart_EndPointx;
  std::vector<Float_t>  mcpart_EndPointy;
  std::vector<Float_t>  mcpart_EndPointz;
  std::vector<Float_t>  mcpart_theta_xz;    
  std::vector<Float_t>  mcpart_theta_yz;    
  std::vector<Int_t>    mcpart_NumberDaughters;
  std::vector<Int_t>    mcpart_TrackId;
  std::vector<Int_t>    mcpart_Mother;

  //MCTrack info
  size_t MaxMCTracks = 0;
  Int_t mctrack_no_primaries;
  std::vector<Int_t>    mctrack_pdg;                      
  std::vector<Int_t>    mctrack_TrackId;

  //MCShower info
  size_t MaxMCShowers = 0;
  Int_t mcshower_no_primaries;
  std::vector<Int_t>    mcshower_pdg;                       
  std::vector<Int_t>    mcshower_TrackId;

  // Event Tree: MeVPrtl Truth
  int n_hnl;
  std::vector<double> mevprtl_decay_pos_x, mevprtl_decay_pos_y, mevprtl_decay_pos_z, mevprtl_decay_pos_t;
  std::vector<double> meson_dmom_e, meson_dmom_x, meson_dmom_y, meson_dmom_z;
  std::vector<double> mevprtl_mom_x, mevprtl_mom_y, mevprtl_mom_z, mevprtl_mom_e;
  std::vector<double> mevprtl_mass;
  std::vector<double> meson_pdg;
  // Event Tree: MC Truth info (HNL only)
  std::vector<double> truth_decay_pos_x, 
                      truth_decay_pos_y, 
                      truth_decay_pos_z, 
                      truth_decay_pos_t;
  std::vector<double> truth_mom_x, 
                      truth_mom_y, 
                      truth_mom_z, 
                      truth_mom_e;
  std::vector<double> truth_pdg; 

 
  bool freadTruth;         ///< Add Truth info to output (to be set via fcl)
  bool fReadMeVPrtl;       ///< Add MeVPrtl info to output (to be set via fcl)
  bool freadMCParticle;    ///< Add MCParticle info to output (to be set via fcl)

  std::string fGenieGenModuleLabel; ///< Label for Genie dataproduct (to be set via fcl)
  std::string fMCParticleModuleLabel; ///< Label for MCParticle dataproduct (to be set via fcl)
  std::string fMCTrackModuleLabel; ///< Label for MCTrack dataproduct (to be set via fcl)
  std::string fMCShowerModuleLabel; ///< Label for MCShower dataproduct (to be set via fcl)
};


test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p) : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel","generator");
  fMCParticleModuleLabel    = p.get<std::string>("MCParticleModuleLabel","largeant");
  fMCTrackModuleLabel    = p.get<std::string>("MCTrackModuleLabel ", "mcreco");
  fMCShowerModuleLabel    = p.get<std::string>("MCShowerModuleLabel ", "mcreco");
  freadTruth         = p.get<bool>("readTruth",true);
  fReadMeVPrtl       = p.get<bool>("ReadMeVPrtl",true);
  freadMCParticle    = p.get<bool>("readMCParticle",true);
  // Call appropriate consumes<>() for any products to be retrieved by this module.
   art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  fTree->Branch("eventID", &fEventID);

if (freadMCParticle){
    //MCParticle
    fTree->Branch("mcpart_pdg",&mcpart_pdg);
    fTree->Branch("mcpart_status",&mcpart_status);    
    fTree->Branch("mcpart_process",&mcpart_process);  
    fTree->Branch("mcpart_endprocess",&mcpart_endprocess);  
    fTree->Branch("mcpart_Eng",&mcpart_Eng);
    fTree->Branch("mcpart_EndE",&mcpart_EndE);
    fTree->Branch("mcpart_Mass",&mcpart_Mass);
    fTree->Branch("mcpart_Px",&mcpart_Px);
    fTree->Branch("mcpart_Py",&mcpart_Py);
    fTree->Branch("mcpart_Pz",&mcpart_Pz);
    fTree->Branch("mcpart_P",&mcpart_P);
    fTree->Branch("mcpart_StartPointx",&mcpart_StartPointx);
    fTree->Branch("mcpart_StartPointy",&mcpart_StartPointy);
    fTree->Branch("mcpart_StartPointz",&mcpart_StartPointz);
    fTree->Branch("mcpart_StartT",&mcpart_StartT);
    fTree->Branch("mcpart_EndT",&mcpart_EndT);
    fTree->Branch("mcpart_EndPointx",&mcpart_EndPointx);
    fTree->Branch("mcpart_EndPointy",&mcpart_EndPointy);
    fTree->Branch("mcpart_EndPointz",&mcpart_EndPointz);         
    fTree->Branch("mcpart_theta_xz",&mcpart_theta_xz);
    fTree->Branch("mcpart_theta_yz",&mcpart_theta_yz);   
    fTree->Branch("mcpart_NumberDaughters",&mcpart_NumberDaughters);
    fTree->Branch("mcpart_TrackId",&mcpart_TrackId);
    fTree->Branch("mcpart_Mother",&mcpart_Mother);

    //MCTrack info
    fTree->Branch("mctrack_no_primaries",&mctrack_no_primaries);
    fTree->Branch("mctrack_pdg",&mctrack_pdg);                        
    fTree->Branch("mctrack_TrackId",&mctrack_TrackId);

    //MCShower info
    fTree->Branch("mcshower_no_primaries",&mcshower_no_primaries);
    fTree->Branch("mcshower_pdg",&mcshower_pdg);                        
    fTree->Branch("mcshower_TrackId",&mcshower_TrackId);
  }
  if (fReadMeVPrtl)
  {
    fTree->Branch("n_hnl",&n_hnl,"n_hnl/I");
    fTree->Branch("mevprtl_decay_pos_x",&mevprtl_decay_pos_x);
    fTree->Branch("mevprtl_decay_pos_y",&mevprtl_decay_pos_y);
    fTree->Branch("mevprtl_decay_pos_z",&mevprtl_decay_pos_z);
    fTree->Branch("mevprtl_decay_pos_t",&mevprtl_decay_pos_t);
    fTree->Branch("mevprtl_mom_x",&mevprtl_mom_x);
    fTree->Branch("mevprtl_mom_y",&mevprtl_mom_y);
    fTree->Branch("mevprtl_mom_z",&mevprtl_mom_z);
    fTree->Branch("mevprtl_mom_e",&mevprtl_mom_e);
    fTree->Branch("mevprtl_mass",&mevprtl_mass);
    fTree->Branch("truth_decay_pos_x",&truth_decay_pos_x);
    fTree->Branch("truth_decay_pos_y",&truth_decay_pos_y);
    fTree->Branch("truth_decay_pos_z",&truth_decay_pos_z);
    fTree->Branch("truth_decay_pos_t",&truth_decay_pos_t);
    fTree->Branch("truth_mom_x",&truth_mom_x);
    fTree->Branch("truth_mom_y",&truth_mom_y);
    fTree->Branch("truth_mom_z",&truth_mom_z);
    fTree->Branch("truth_mom_e",&truth_mom_e);
    fTree->Branch("truth_pdg",&truth_pdg);
    fTree->Branch("meson_dmom_x",&meson_dmom_x);
    fTree->Branch("meson_dmom_y",&meson_dmom_y);
    fTree->Branch("meson_dmom_z",&meson_dmom_z);
    fTree->Branch("meson_dmom_e",&meson_dmom_e);
    fTree->Branch("mmeson_pdg",&meson_pdg);
  }

}

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
 

}

void test::AnalyzeEvents::analyze(art::Event const& evt)
{
  ResetVars();
  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();
  _t0 = 0.;
  // Implementation of required member function here.
  fEventID = evt.id().event();
  
  if (freadMCParticle){
    //MCParticle
    art::Handle<std::vector<simb::MCParticle>> MCParticleListHandle;
    std::vector<art::Ptr<simb::MCParticle>> MCParticleList;
    if (evt.getByLabel(fMCParticleModuleLabel,MCParticleListHandle)){
      art::fill_ptr_vector(MCParticleList,MCParticleListHandle);
      mcpart_no_primaries = MCParticleList.size();
      ResizeMCParticle(MCParticleList.size()); //Set vectors
      for (size_t iMCPart = 0; iMCPart < MCParticleList.size(); iMCPart++){
        art::Ptr<simb::MCParticle> pPart = MCParticleList[iMCPart]; //get particle pointer
        //Geant info
        mcpart_Mother[iMCPart] = pPart->Mother();
        mcpart_TrackId[iMCPart] = pPart->TrackId();
        mcpart_pdg[iMCPart] = pPart->PdgCode();
        mcpart_status[iMCPart] =  pPart->StatusCode();
        mcpart_process[iMCPart] =  pPart->Process();
        mcpart_endprocess[iMCPart] =  pPart->EndProcess();
        mcpart_Eng[iMCPart] = pPart->E();
        mcpart_EndE[iMCPart] = pPart->EndE();
        mcpart_Mass[iMCPart] = pPart->Mass();
        mcpart_Px[iMCPart] = pPart->Px();
        mcpart_Py[iMCPart] = pPart->Py();
        mcpart_Pz[iMCPart] = pPart->Pz();
        mcpart_P[iMCPart] = pPart->Momentum().Vect().Mag();
        mcpart_StartPointx[iMCPart] = pPart->Vx();
        mcpart_StartPointy[iMCPart] = pPart->Vy();
        mcpart_StartPointz[iMCPart] = pPart->Vz();
        mcpart_StartT[iMCPart] = pPart->T();
        mcpart_EndPointx[iMCPart] = pPart->EndPosition()[0];
        mcpart_EndPointy[iMCPart] = pPart->EndPosition()[1];
        mcpart_EndPointz[iMCPart] = pPart->EndPosition()[2];
        mcpart_EndT[iMCPart] = pPart->EndT();
        mcpart_theta_xz[iMCPart] =  std::atan2(pPart->Px(), pPart->Pz());
        mcpart_theta_yz[iMCPart] =  std::atan2(pPart->Py(), pPart->Pz());
        mcpart_NumberDaughters[iMCPart] = pPart->NumberDaughters();
      }
    }//endif get label
    else {
      std::cout << "Failed to get MCParticle data product." << std::endl;
    }
    //MCTracks
    art::Handle<std::vector<sim::MCTrack>> MCTrackListHandle;
    std::vector<art::Ptr<sim::MCTrack>> MCTrackList;
    if (evt.getByLabel(fMCTrackModuleLabel,MCTrackListHandle)){
      art::fill_ptr_vector(MCTrackList,MCTrackListHandle);
      mctrack_no_primaries = MCTrackList.size();
      ResizeMCTrack(MCTrackList.size());
      for (size_t iMCTrack = 0; iMCTrack < MCTrackList.size(); iMCTrack++){
        art::Ptr<sim::MCTrack> pTrack = MCTrackList[iMCTrack]; //get particle pointer
        mctrack_pdg[iMCTrack] = pTrack->PdgCode();
        mctrack_TrackId[iMCTrack] = pTrack->TrackID();
      }
    }
    else{
      //std::cout << "Failed to get MCTrack data product." << std::endl;
    }

    //MCShowers
    art::Handle<std::vector<sim::MCShower>> MCShowerListHandle;
    std::vector<art::Ptr<sim::MCShower>> MCShowerList;
    if (evt.getByLabel(fMCShowerModuleLabel,MCShowerListHandle)){
      art::fill_ptr_vector(MCShowerList,MCShowerListHandle);
      mcshower_no_primaries = MCShowerList.size();
      ResizeMCShower(MCShowerList.size());
      for (size_t iMCShower = 0; iMCShower < MCShowerList.size(); iMCShower++){
        art::Ptr<sim::MCShower> pShower = MCShowerList[iMCShower]; //get particle pointer
        mcshower_pdg[iMCShower] = pShower->PdgCode();
        mcshower_TrackId[iMCShower] = pShower->TrackID();
      }
    }
    else{
     // std::cout << "Failed to get MCShower data product." << std::endl;
    }

  } // end read mcparticle

  if (fReadMeVPrtl)
  {

    // Clean up vectors before filling
    mevprtl_decay_pos_x.clear();
    mevprtl_decay_pos_y.clear();
    mevprtl_decay_pos_z.clear();
    mevprtl_decay_pos_t.clear();
    mevprtl_mom_x.clear();
    mevprtl_mom_y.clear();
    mevprtl_mom_z.clear();
    mevprtl_mom_e.clear();
    mevprtl_mass.clear();
    truth_decay_pos_x.clear();
    truth_decay_pos_y.clear();
    truth_decay_pos_z.clear();
    truth_decay_pos_t.clear();
    truth_mom_x.clear();
    truth_mom_y.clear();
    truth_mom_z.clear();
    truth_mom_e.clear();
    truth_pdg.clear();
    n_hnl = 0;
    meson_dmom_x.clear();
    meson_dmom_y.clear();
    meson_dmom_z.clear();
    meson_dmom_e.clear();
    meson_pdg.clear();

    // Get MeVPrtl (HNL)
    art::Handle<std::vector<evgen::ldm::MeVPrtlTruth>> mevptHandle;
    std::string fMeVPrtlTruthModuleLabel="generator";
    evt.getByLabel(fMeVPrtlTruthModuleLabel, mevptHandle);
    std::vector<art::Ptr<evgen::ldm::MeVPrtlTruth>> mevptVec;
    if (mevptHandle.isValid()) art::fill_ptr_vector(mevptVec, mevptHandle);
    // std::cout<<"MeVPrtl: "<<mevptVec.size()<<std::endl;
    // std::cout<<"MeVPrtl: "<<mevptHandle.isValid()<<std::endl;
    // std::cout<<"MeVPrtl: here!"<<std::endl;
    for (auto const &mevpt: mevptVec)
      {
        // std::cout<<"MeVPrtl: "<<mevpt->mass<<std::endl;
        mevprtl_decay_pos_x.push_back(mevpt->decay_pos.X()); 
        mevprtl_decay_pos_y.push_back(mevpt->decay_pos.Y()); 
        mevprtl_decay_pos_z.push_back(mevpt->decay_pos.Z()); 
        mevprtl_decay_pos_t.push_back(mevpt->decay_pos.T());

        meson_dmom_x.push_back(mevpt->meson_dmom.X());
        meson_dmom_y.push_back(mevpt->meson_dmom.Y());
        meson_dmom_z.push_back(mevpt->meson_dmom.Z());
        meson_dmom_e.push_back(mevpt->meson_dmom.E());
        meson_pdg.push_back(mevpt->meson_pdg);

        mevprtl_mom_x.push_back(mevpt->mevprtl_mom.X());
        mevprtl_mom_y.push_back(mevpt->mevprtl_mom.Y());
        mevprtl_mom_z.push_back(mevpt->mevprtl_mom.Z());
        mevprtl_mom_e.push_back(mevpt->mevprtl_mom.E());
        mevprtl_mass.push_back(mevpt->mass);
        n_hnl++;
      }
      
    // Get All MCTruth handles
    std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = evt.getMany<std::vector<simb::MCTruth>>();

      //Loop over handles
    for(auto const& MCTruthHandle : MCTruthHandles)
      {
        std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
        art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

        //Loop over MCTruths in handle
        for(auto const& truth : MCTruthVec)
        {
          const bool _is_hnl = truth->Origin() == 0; // 0: HNL, everything else: cosmics
          if (!_is_hnl) continue;
          
          int N = truth->NParticles();
        
          for (int i = 0; i < N; ++i) // Loop over particles in simb::MCtruth
          {
            const simb::MCParticle &nu = truth->GetParticle(i);
            float E = nu.E();
            const int pdg = nu.PdgCode();

            const TLorentzVector &v4_f = nu.EndPosition();
            const TLorentzVector &p4_f = nu.EndMomentum();

            auto x_f = v4_f.X();
            auto y_f = v4_f.Y();
            auto z_f = v4_f.Z();
            auto t_f = v4_f.T();

            auto px_f = p4_f.X();
            auto py_f = p4_f.Y();
            auto pz_f = p4_f.Z();

            truth_decay_pos_x.push_back(x_f);
            truth_decay_pos_y.push_back(y_f);
            truth_decay_pos_z.push_back(z_f);
            truth_decay_pos_t.push_back(t_f);
            truth_mom_x      .push_back(px_f);
            truth_mom_y      .push_back(py_f);
            truth_mom_z      .push_back(pz_f);
            truth_mom_e      .push_back(E);
            truth_pdg        .push_back(pdg);
            
          }
        }

      }
  }



  fTree->Fill();

}




void test::AnalyzeEvents::ResetVars() {


  _run = -99999;
  _subrun = -99999;
  _event = -99999;
  _evttime = -99999;
  _t0 = -99999;

  mcpart_no_primaries = 0;
  mctrack_no_primaries = 0;
  mcshower_no_primaries = 0;

}
void test::AnalyzeEvents::ResizeMCParticle(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxMCParticles = (size_t) std::max(nParticles, 1);
              
  mcpart_pdg.assign(MaxMCParticles,DEFAULT_VALUE);              
  mcpart_status.assign(MaxMCParticles,DEFAULT_VALUE); 
  mcpart_process.assign(MaxMCParticles,"Dummy"); 
  mcpart_endprocess.assign(MaxMCParticles,"Dummy");           
  mcpart_Eng.assign(MaxMCParticles,DEFAULT_VALUE);              
  mcpart_EndE.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Mass.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Px.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Py.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Pz.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_P.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartPointx.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartPointy.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartPointz.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_StartT.assign(MaxMCParticles,DEFAULT_VALUE);  
  mcpart_EndT.assign(MaxMCParticles,DEFAULT_VALUE);          
  mcpart_EndPointx.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_EndPointy.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_EndPointz.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_theta_xz.assign(MaxMCParticles,DEFAULT_VALUE);    
  mcpart_theta_yz.assign(MaxMCParticles,DEFAULT_VALUE);    
  mcpart_NumberDaughters.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_TrackId.assign(MaxMCParticles,DEFAULT_VALUE);
  mcpart_Mother.assign(MaxMCParticles,DEFAULT_VALUE);

}

void test::AnalyzeEvents::ResizeMCTrack(int nTracks) {
  MaxMCTracks = (size_t) std::max(nTracks,1);

  mctrack_pdg.assign(MaxMCTracks,DEFAULT_VALUE);
  mctrack_TrackId.assign(MaxMCTracks,DEFAULT_VALUE);
}

void test::AnalyzeEvents::ResizeMCShower(int nShowers) {
  MaxMCShowers = (size_t) std::max(nShowers,1);

  mcshower_pdg.assign(MaxMCShowers,DEFAULT_VALUE);
  mcshower_TrackId.assign(MaxMCShowers,DEFAULT_VALUE);
}


// void test::AnalyzeEvents::endJob()
// {
//   // Implementation of optional member function here.
// }

DEFINE_ART_MODULE(test::AnalyzeEvents)
