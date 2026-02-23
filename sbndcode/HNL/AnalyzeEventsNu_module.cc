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
  class AnalyzeEventsNu;
}


class test::AnalyzeEventsNu : public art::EDAnalyzer {
public:
  explicit AnalyzeEventsNu(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEventsNu(AnalyzeEventsNu const&) = delete;
  AnalyzeEventsNu(AnalyzeEventsNu&&) = delete;
  AnalyzeEventsNu& operator=(AnalyzeEventsNu const&) = delete;
  AnalyzeEventsNu& operator=(AnalyzeEventsNu&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
 // void endJob() override;

private:
  void ResetVars();
  void ResizeMCNeutrino(int nNeutrinos);
  /// Resize the data structure for Genie primaries
  void ResizeGenie(int nPrimaries);
    /// Resize the data structure for MCParticles
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
  //mctruth information
  size_t MaxMCNeutrinos;     ///! The number of MCNeutrinos there is currently room for
  Int_t     mcevts_truth;                     ///< number of neutrino Int_teractions in the spill
  std::vector<Int_t>     nuScatterCode_truth; ///< Scattering code given by Genie for each neutrino
  std::vector<Int_t>     nuID_truth;          ///< Unique ID of each true neutrino
  std::vector<Int_t>     nuPDG_truth;         ///< neutrino PDG code
  std::vector<Int_t>     ccnc_truth;          ///< 0=CC 1=NC
  std::vector<Int_t>     mode_truth;          ///< 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  std::vector<Float_t>   enu_truth;           ///< true neutrino energy
  std::vector<Float_t>   Q2_truth;            ///< Momentum transfer squared
  std::vector<Float_t>   W_truth;             ///< hadronic invariant mass
  std::vector<Int_t>     hitnuc_truth;        ///< hit nucleon
  std::vector<Float_t>   nuvtxx_truth;        ///< neutrino vertex x
  std::vector<Float_t>   nuvtxy_truth;        ///< neutrino vertex y
  std::vector<Float_t>   nuvtxz_truth;        ///< neutrino vertex z
  std::vector<Float_t>   nu_dcosx_truth;      ///< neutrino dcos x
  std::vector<Float_t>   nu_dcosy_truth;      ///< neutrino dcos y
  std::vector<Float_t>   nu_dcosz_truth;      ///< neutrino dcos z
  std::vector<Float_t>   lep_mom_truth;       ///< lepton momentum
  std::vector<Float_t>   lep_dcosx_truth;     ///< lepton dcos x
  std::vector<Float_t>   lep_dcosy_truth;     ///< lepton dcos y
  std::vector<Float_t>   lep_dcosz_truth;     ///< lepton dcos z

  // flux information
  std::vector<Float_t>  tpx_flux;             ///< Px of parent particle leaving BNB target
  std::vector<Float_t>  tpy_flux;             ///< Py of parent particle leaving BNB target
  std::vector<Float_t>  tpz_flux;             ///< Pz of parent particle leaving BNB target
  std::vector<Int_t>    tptype_flux;         ///< Type of parent particle leaving BNB target

  // //genie information
  size_t MaxGeniePrimaries = 0;
  Int_t     genie_no_primaries;
  std::vector<Int_t>     genie_primaries_pdg;
  std::vector<Float_t>   genie_Eng;
  std::vector<Float_t>   genie_Px;
  std::vector<Float_t>   genie_Py;
  std::vector<Float_t>   genie_Pz;
  std::vector<Float_t>   genie_P;
  std::vector<Float_t>   genie_Pt_em;
  std::vector<Float_t>   genie_Pt_ep;
  std::vector<Float_t>   genie_Pt_pion;
  std::vector<Int_t>     genie_status_code;
  std::vector<Float_t>   genie_mass;
  std::vector<Int_t>     genie_trackID;
  std::vector<Int_t>     genie_ND;
  std::vector<Int_t>     genie_mother;

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

  std::string fGenieGenModuleLabel; ///< Label for Genie dataproduct (to be set via fcl)
  std::string fMCParticleModuleLabel; ///< Label for MCParticle dataproduct (to be set via fcl)
  std::string fMCTrackModuleLabel; ///< Label for MCTrack dataproduct (to be set via fcl)
  std::string fMCShowerModuleLabel; ///< Label for MCShower dataproduct (to be set via fcl)
 
  bool freadTruth;         ///< Add Truth info to output (to be set via fcl)
  bool freadMCParticle;    ///< Add MCParticle info to output (to be set via fcl)

};


test::AnalyzeEventsNu::AnalyzeEventsNu(fhicl::ParameterSet const& p) : EDAnalyzer{p}  // ,
  // More initializers here.
{
  fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel","generator");
  fMCParticleModuleLabel    = p.get<std::string>("MCParticleModuleLabel ", "largeant");
  fMCTrackModuleLabel    = p.get<std::string>("MCTrackModuleLabel ", "mcreco");
  fMCShowerModuleLabel    = p.get<std::string>("MCShowerModuleLabel ", "mcreco");

  freadTruth         = p.get<bool>("readTruth",true);
  freadMCParticle         = p.get<bool>("readMCParticle",true);
  // Call appropriate consumes<>() for any products to be retrieved by this module.
   art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  fTree->Branch("eventID", &fEventID);

if (freadTruth) {
    fTree->Branch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
    fTree->Branch("nuScatterCode_truth",&nuScatterCode_truth);
    fTree->Branch("nuID_truth",&nuID_truth);
    fTree->Branch("nuPDG_truth",&nuPDG_truth);
    fTree->Branch("ccnc_truth",&ccnc_truth);
    fTree->Branch("mode_truth",&mode_truth);
    fTree->Branch("enu_truth",&enu_truth);
    fTree->Branch("Q2_truth",&Q2_truth);
    fTree->Branch("W_truth",&W_truth);
    fTree->Branch("hitnuc_truth",&hitnuc_truth);
    fTree->Branch("nuvtxx_truth",&nuvtxx_truth);
    fTree->Branch("nuvtxy_truth",&nuvtxy_truth);
    fTree->Branch("nuvtxz_truth",&nuvtxz_truth);
    fTree->Branch("nu_dcosx_truth",&nu_dcosx_truth);
    fTree->Branch("nu_dcosy_truth",&nu_dcosy_truth);
    fTree->Branch("nu_dcosz_truth",&nu_dcosz_truth);
    fTree->Branch("lep_mom_truth",&lep_mom_truth);
    fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth);
    fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth);
    fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth);

    fTree->Branch("tpx_flux",&tpx_flux);
    fTree->Branch("tpy_flux",&tpy_flux);
    fTree->Branch("tpz_flux",&tpz_flux);
    fTree->Branch("tptype_flux",&tptype_flux);

    fTree->Branch("genie_no_primaries",&genie_no_primaries);
    fTree->Branch("genie_primaries_pdg",&genie_primaries_pdg);
    fTree->Branch("genie_Eng",&genie_Eng);
    fTree->Branch("genie_Px",&genie_Px);
    fTree->Branch("genie_Py",&genie_Py);
    fTree->Branch("genie_Pz",&genie_Pz);
    fTree->Branch("genie_Pt_em",&genie_Pt_em);
    fTree->Branch("genie_Pt_ep",&genie_Pt_ep);
    fTree->Branch("genie_Pt_pion",&genie_Pt_pion);
    fTree->Branch("genie_P",&genie_P);
    fTree->Branch("genie_status_code",&genie_status_code);
    fTree->Branch("genie_mass",&genie_mass);
    fTree->Branch("genie_trackID",&genie_trackID);
    fTree->Branch("genie_ND",&genie_ND);
    fTree->Branch("genie_mother",&genie_mother);
  }
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
 
}

void test::AnalyzeEventsNu::beginJob()
{
  // Implementation of optional member function here.
 

}

void test::AnalyzeEventsNu::analyze(art::Event const& evt)
{
  ResetVars();
  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();
  _t0 = 0.;
  // Implementation of required member function here.
  fEventID = evt.id().event();
  
  if (freadTruth){
    //Genie
    int nGeniePrimaries = 0, nMCNeutrinos = 0;
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle)){
      art::fill_ptr_vector(mclist, mctruthListHandle);
    }else {
      std::cout << "Failed to get Genie data product." << std::endl;
    }

    art::Ptr<simb::MCTruth> mctruth;

      if (!mclist.empty()) {//at least one mc record

        mctruth = mclist[0];

        if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();

    } // if have MC truth
      MF_LOG_DEBUG("HitDumper") << "Expected " << nGeniePrimaries << " GENIE particles";

    //Initially call the number of neutrinos to be stored the number of MCTruth objects.  This is not strictly true i.e. BNB + cosmic overlay but we will count the number of neutrinos later
    nMCNeutrinos = mclist.size();

    ResizeGenie(nGeniePrimaries);
    ResizeMCNeutrino(nMCNeutrinos);

    mcevts_truth = mclist.size();
    //Brailsford 2017/10/16
    //Issue 17917
    //To keep a 1:1 between neutrinos and 'flux' we need the assns
    art::FindManyP<simb::MCFlux> fmFluxNeutrino(mctruthListHandle, evt, fGenieGenModuleLabel);
    // Get GTruth information for scattering code
    art::FindManyP< simb::GTruth > fmgt( mctruthListHandle, evt, fGenieGenModuleLabel );

    if (mcevts_truth > 0){//at least one mc record

      //Brailsford 2017/10/16
      //Issue 17917
      //Loop over every truth in the spill rather than just looking at the first one.
      //Because MCTruth could be a neutrino OR something else (e.g. cosmics) we are going to have to count up how many neutrinos there are
      mcevts_truth = 0;
      for (unsigned int i_mctruth = 0; i_mctruth < mclist.size(); i_mctruth++){
        //fetch an mctruth
        art::Ptr<simb::MCTruth> curr_mctruth = mclist[i_mctruth];
        //Check if it's a neutrino
        if (!curr_mctruth->NeutrinoSet()) continue;

        // Genie Truth association only for the neutrino
        if (fmgt.size()>i_mctruth) {
          std::vector< art::Ptr<simb::GTruth> > mcgtAssn = fmgt.at(i_mctruth);

          nuScatterCode_truth[i_mctruth] = mcgtAssn[0]->fGscatter;
        } else {
          nuScatterCode_truth[i_mctruth] = -1.;
        }

        nuPDG_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().PdgCode();
        ccnc_truth[i_mctruth] = curr_mctruth->GetNeutrino().CCNC();
        mode_truth[i_mctruth] = curr_mctruth->GetNeutrino().Mode();
        Q2_truth[i_mctruth] = curr_mctruth->GetNeutrino().QSqr();
        W_truth[i_mctruth] = curr_mctruth->GetNeutrino().W();
        hitnuc_truth[i_mctruth] = curr_mctruth->GetNeutrino().HitNuc();
        enu_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().E();
        nuvtxx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vx();
        nuvtxy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vy();
        nuvtxz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vz();
        if (curr_mctruth->GetNeutrino().Nu().P()){
          nu_dcosx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Px()/curr_mctruth->GetNeutrino().Nu().P();
          nu_dcosy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Py()/curr_mctruth->GetNeutrino().Nu().P();
          nu_dcosz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Pz()/curr_mctruth->GetNeutrino().Nu().P();
        }
        lep_mom_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().P();
        if (curr_mctruth->GetNeutrino().Lepton().P()){
          lep_dcosx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Px()/curr_mctruth->GetNeutrino().Lepton().P();
          lep_dcosy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Py()/curr_mctruth->GetNeutrino().Lepton().P();
          lep_dcosz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Pz()/curr_mctruth->GetNeutrino().Lepton().P();
        }
        //Brailsford
        //2017/10/17
        //Issue 12918
        //Use the art::Ptr key as the neutrino's unique ID
        nuID_truth[i_mctruth] = curr_mctruth.key();
        //We need to also store N 'flux' neutrinos per event so now check that the FindOneP is valid and, if so, use it!
        if (fmFluxNeutrino.isValid()){
          if (fmFluxNeutrino.at(0).size()>i_mctruth){
          art::Ptr<simb::MCFlux> curr_mcflux = fmFluxNeutrino.at(0).at(i_mctruth);
          tpx_flux[i_mctruth] = curr_mcflux->ftpx;
          tpy_flux[i_mctruth] = curr_mcflux->ftpy;
          tpz_flux[i_mctruth] = curr_mcflux->ftpz;
          tptype_flux[i_mctruth] = curr_mcflux->ftptype;
          }
        }

        //Let's increase the neutrino count
        mcevts_truth++;
      }

      if (mctruth->NeutrinoSet()){
        //genie particles information
        genie_no_primaries = mctruth->NParticles();

        size_t StoreParticles = std::min((size_t) genie_no_primaries, MaxGeniePrimaries);
        if (genie_no_primaries > (int) StoreParticles) {
          // got this error? it might be a bug,
          // since the structure should have enough room for everything
          mf::LogError("HitDumper") << "event has "
            << genie_no_primaries << " MC particles, only "
            << StoreParticles << " stored in tree";
        }
        for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
          const simb::MCParticle& part(mctruth->GetParticle(iPart));
          genie_primaries_pdg[iPart]=part.PdgCode();
          genie_Eng[iPart]=part.E();
          genie_Px[iPart]=part.Px();
          genie_Py[iPart]=part.Py();
          genie_Pz[iPart]=part.Pz();
          genie_P[iPart]=part.P();
          genie_status_code[iPart]=part.StatusCode();
          genie_mass[iPart]=part.Mass();
          genie_trackID[iPart]=part.TrackId();
          genie_ND[iPart]=part.NumberDaughters();
          genie_mother[iPart]=part.Mother();
          // if (part.PdgCode() == 111 && part.StatusCode() == 1)
          // {
          //   genie_Pt_pion[iPart] = sqrt(pow(part.Px(),2) + pow(part.Py(),2));
          //   std::cout << "status " << part.StatusCode() << std::endl;
          //   std::cout << "event " << _event << std::endl;
            //  }

          
          //Seleccioname los electrones con pdg 11 y positrones con pdg -11 que provengan de la interacción de un neutrino
          // if ((part.PdgCode() == 11 || part.PdgCode() == -11) && part.Mother() == 0) {
          //   genie_Pt[iPart] = sqrt(pow(part.Px(), 2) + pow(part.Py(), 2));
          // }
          // if (part.PdgCode() == abs(11) && part.Mother() == 0)
          //  {
          //    genie_Pt[iPart] = sqrt(pow(part.Px(),2) + pow(part.Py(),2));
          //  }
          
        } // for particle
        //Hazme un for quwe recorra todas las particulas del evento y seleccioname los eventos que tengan un alectrón y positrṕon con la misma madre
        // for (size_t iPart = 0; iPart < StoreParticles; ++iPart)
        // {
        //   const simb::MCParticle& part(mctruth->GetParticle(iPart));
        //   if (part.PdgCode() == 13|| part.PdgCode() == -13)
        //   {
        //     std::cout << "pdg: " << part.PdgCode() << std::endl;
        //     std::cout << "E: " << part.E() << std::endl;
        //     std::cout << "Mother: " << part.Mother() << std::endl;
        //     std::cout << "event: " << _event  << std::endl;
            
        //     // for (size_t jPart = 0; jPart < StoreParticles; ++jPart)
        //     // {
        //     //   const simb::MCParticle& part2(mctruth->GetParticle(jPart));
        //     //   if (part2.PdgCode() == -11 || part2.PdgCode() == 11) 
        //     //   {
        //     //     if (part.Mother() == part2.Mother())
        //     //      {
        //     //       std::cout << "ee " << part.Mother() << std::endl;
        //     //       std::cout << "pdg " << part.PdgCode() << std::endl;
        //     //       std::cout << "pdg " << part2.PdgCode() << std::endl;
        //     //       genie_Pt_em[iPart] = sqrt(pow(part.Px(),2) + pow(part.Py(),2));
        //     //       genie_Pt_ep[iPart] = sqrt(pow(part2.Px(),2) + pow(part2.Py(),2));
        //     //      }
        //     //   }
        //     // }
        //   }
       //}
      } //if neutrino set
    }//if (mcevts_truth)



  }//if (fReadTruth){

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
      //std::cout << "Failed to get MCShower data product." << std::endl;
    }

  } // end read mcparticle


  fTree->Fill();

}




void test::AnalyzeEventsNu::ResetVars() {


  _run = -99999;
  _subrun = -99999;
  _event = -99999;
  _evttime = -99999;
  _t0 = -99999;

  mcevts_truth = 0;
  genie_no_primaries = 0;
}

void test::AnalyzeEventsNu::ResizeMCNeutrino(int nNeutrinos) {

  //min size is 1, to guarantee an address
  MaxMCNeutrinos = (size_t) std::max(nNeutrinos, 1);

  nuScatterCode_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuID_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuPDG_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  ccnc_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  mode_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  enu_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  Q2_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  W_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  hitnuc_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuvtxx_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuvtxy_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nuvtxz_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nu_dcosx_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nu_dcosy_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  nu_dcosz_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_mom_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_dcosx_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_dcosy_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  lep_dcosz_truth.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  //Also resize the flux information here as it's a 1:1 with the MCNeutrino
  tpx_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  tpy_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  tpz_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);
  tptype_flux.assign(MaxMCNeutrinos, DEFAULT_VALUE);

}

void test::AnalyzeEventsNu::ResizeGenie(int nPrimaries) {

  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);

  genie_primaries_pdg.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Eng.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Px.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Py.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Pz.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_P.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_Pt_em.assign(MaxGeniePrimaries, -0.2);
  genie_Pt_ep.assign(MaxGeniePrimaries, -0.2);
  genie_Pt_pion.assign(MaxGeniePrimaries, -0.2);
  genie_status_code.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_mass.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_trackID.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_ND.assign(MaxGeniePrimaries, DEFAULT_VALUE);
  genie_mother.assign(MaxGeniePrimaries, DEFAULT_VALUE);

}
void test::AnalyzeEventsNu::ResizeMCParticle(int nParticles) {

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

void test::AnalyzeEventsNu::ResizeMCTrack(int nTracks) {
  MaxMCTracks = (size_t) std::max(nTracks,1);

  mctrack_pdg.assign(MaxMCTracks,DEFAULT_VALUE);
  mctrack_TrackId.assign(MaxMCTracks,DEFAULT_VALUE);
}

void test::AnalyzeEventsNu::ResizeMCShower(int nShowers) {
  MaxMCShowers = (size_t) std::max(nShowers,1);

  mcshower_pdg.assign(MaxMCShowers,DEFAULT_VALUE);
  mcshower_TrackId.assign(MaxMCShowers,DEFAULT_VALUE);
}




/*
void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}
*/
DEFINE_ART_MODULE(test::AnalyzeEventsNu)
