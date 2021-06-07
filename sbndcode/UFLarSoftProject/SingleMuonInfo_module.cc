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
double length(const simb::MCParticle& part, TVector3& start, TVector3& end)
{
  art::ServiceHandle<geo::Geometry> geom;
  double xmin = -2.0 * geom->DetHalfWidth();
  double xmax = 2.0 * geom->DetHalfWidth();
  double ymin = -geom->DetHalfHeight();
  double ymax = geom->DetHalfHeight();
  double zmin = 0.;
  double zmax = geom->DetLength();
  
  double result = 0.;
  TVector3 disp;
  int n = part.NumberTrajectoryPoints();
  bool first = true;
  
 for(int i = 0; i < n; ++i) {
    double mypos[3] = {part.Vx(i), part.Vy(i), part.Vz(i)};
    if(mypos[0] >= xmin && mypos[0] <= xmax && mypos[1] >= ymin && mypos[1] <= ymax && mypos[2] >= zmin && mypos[2] <= zmax){
       double xGen   = part.Vx(i);
       double newX = xGen;
       
       TVector3 pos(newX,part.Vy(i),part.Vz(i));
       if(first){
          start = pos;
       }
       else {
          disp -= pos;
          result += disp.Mag();
       }
       first=false;
       disp = pos;
       end = pos;
   }
  }
  return result;
}	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

vector<double> G4_eng_dep(int P_id){
       art::ServiceHandle<cheat::BackTrackerService> bt_serv;
       vector<double> eng_vec;
       const geo::View_t views[3]={geo::kU, geo::kV, geo::kW};
       double totalE_particle = 0.;
       for(auto const &view : views){
           std::vector< const sim::IDE * >  ides = bt_serv->TrackIdToSimIDEs_Ps(P_id, view);
	   double plne_eng = 0.;
	   if( ides.size() ){
	      for (auto const &ide: ides){ 
		   totalE_particle += ide->energy;
		   plne_eng += ide->energy;
	      }
	   }
	   
	   eng_vec.push_back(plne_eng);
       }
       
       totalE_particle /= 3.0;
       eng_vec.push_back(totalE_particle);
       return eng_vec;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
class SingleMuonInfo : public art::EDAnalyzer {
public:
    using CRTHit = sbn::crt::CRTHit;
    explicit SingleMuonInfo(fhicl::ParameterSet const& pset);
    virtual ~SingleMuonInfo();

    void beginJob();
    void analyze(const art::Event& evt);

 private:
    void ClearVecs();
		 
    TTree* fTruthTree; // Let's save G4 level information
    TTree* fRecHitTree; // Let's save reconstructed hit level informaion
    TTree* fRecTrackTree; // Let's save reconstructed track level information
    TTree* fRecShowerTree; // Let's save reconstructed shower information
    
    // Common Variables
    
    Int_t    frun;                  
    Int_t    fsubrun;               
    Int_t    fevent;
    
    // Genie level Variables
    
    Int_t fn_nu;
    vector<bool> fgen_nu_set;
    vector<int>  fgen_n_part;
    vector<int>  fgen_nu_pdg;
 
    Int_t fn_g4;
    vector<int> fg4_pdg;
    vector<string> fg4_st_process;
    vector<double> fg4_plen;
    vector<bool> fg4_nu_origin;
    vector<int> fg4_nu_ccnc;
    vector<int> fg4_nu_mode;
    vector<double> fg4_nu_stx;
    vector<double> fg4_nu_sty;
    vector<double> fg4_nu_stz;
    
    // Hit level Varialble;
    
    Int_t fn_hits;
    vector<int> fhit_channel;
    vector<int> fhit_wire;
    vector<int> fhit_plane;
    vector<int> fhit_tpc;
    vector<int> fhit_cryostat;
    
    art::InputTag fGenLabel;
    art::InputTag fHitsModuleLabel;
    art::InputTag fSimLabel;
    art::InputTag fOpHitModuleLabel;
    art::InputTag fOpFlashModuleLabel0;
    art::InputTag fOpFlashModuleLabel1;
    art::InputTag fCrtHitModuleLabel;
    art::InputTag fCrtTrackModuleLabel;
    bool fuse_run_genie;
    
    map<int,art::InputTag> fFlashLabels;
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    const cheat::ParticleInventory *inventory_service;
 }; // class SingleMuonInfo

//========================================================================
SingleMuonInfo::SingleMuonInfo(fhicl::ParameterSet const& pset) :
EDAnalyzer(pset),
fGenLabel(pset.get<art::InputTag>("GenLabel","generator")),
fHitsModuleLabel(pset.get<art::InputTag>("HitsModuleLabel","gaushit")),		
fSimLabel(pset.get<art::InputTag>("SimLabel","largeant")),
fOpHitModuleLabel(pset.get<art::InputTag>("OpHitModuleLabel","ophitpmt")),				
fOpFlashModuleLabel0(pset.get<art::InputTag>("OpFlashModuleLabel0","opflashtpc0")),
fOpFlashModuleLabel1(pset.get<art::InputTag>("OpFlashModuleLabel1","opflashtpc1")),
fCrtHitModuleLabel(pset.get<art::InputTag>("CrtHitModuleLabel","crthit")),
fCrtTrackModuleLabel(pset.get<art::InputTag>("CrtTrackModuleLabel","crttrack")),
fuse_run_genie(pset.get<bool>("use_run_genie",true))	
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
  
  fTruthTree->Branch("n_nu", &fn_nu, "n_nu/I");
  fTruthTree->Branch("gen_nu_set", &fgen_nu_set);
  fTruthTree->Branch("gen_n_part", &fgen_n_part);
  fTruthTree->Branch("gen_nu_pdg", &fgen_nu_pdg);
  
  fTruthTree->Branch("n_g4", &fn_g4, "n_g4/I");
  fTruthTree->Branch("g4_pdg", &fg4_pdg);
  fTruthTree->Branch("g4_st_process", &fg4_st_process);
  fTruthTree->Branch("g4_plen", &fg4_plen);
  fTruthTree->Branch("g4_nu_origin", &fg4_nu_origin);
  fTruthTree->Branch("g4_nu_ccnc", &fg4_nu_ccnc);
  fTruthTree->Branch("g4_nu_mode", &fg4_nu_mode);
  fTruthTree->Branch("g4_nu_stx", &fg4_nu_stx);
  fTruthTree->Branch("g4_nu_sty", &fg4_nu_sty);
  fTruthTree->Branch("g4_nu_stz", &fg4_nu_stz);
  
  fRecHitTree = tfs->make<TTree>("RecHitTree","");
  fRecHitTree->Branch("run", &frun, "run/I");
  fRecHitTree->Branch("subrun", &fsubrun, "subrun/I");
  fRecHitTree->Branch("event", &fevent, "event/I");
  fRecHitTree->Branch("n_hits", &fn_hits, "n_hits/I");
  fTruthTree->Branch("hit_channel", &fhit_channel);
  fTruthTree->Branch("hit_wire", &fhit_wire);
  fTruthTree->Branch("hit_tpc", &fhit_tpc);
  fTruthTree->Branch("hit_plane", &fhit_plane);
  fTruthTree->Branch("hit_cryostat", &fhit_cryostat);
  
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
     
     art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
     std::vector<art::Ptr<simb::MCTruth> > mclist;
     if (evt.getByLabel(fGenLabel,mctruthListHandle))
         art::fill_ptr_vector(mclist, mctruthListHandle);
     
     art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
     std::vector< art::Ptr<simb::MCParticle> > ptList;
     if (evt.getByLabel(fSimLabel, mcParticleHandle))
         art::fill_ptr_vector(ptList, mcParticleHandle); 
     
     art::Handle< std::vector<recob::Hit> > hitListHandle;
     std::vector<art::Ptr<recob::Hit> > hitlist;
     if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
     art::fill_ptr_vector(hitlist, hitListHandle);
     
     inventory_service=lar::providerFrom<cheat::ParticleInventoryService>();
     
     if(fuse_run_genie){
        fn_nu = mclist.size();
	for(auto const& my_Nu : mclist){
	    fgen_nu_set.push_back(my_Nu->NeutrinoSet());
	    fgen_n_part.push_back(my_Nu->NParticles());
	    const simb::MCNeutrino & my_neutrino = my_Nu->GetNeutrino();
	    fgen_nu_pdg.push_back(my_neutrino.Nu().PdgCode());
	    // save nu energery
	    // save nu st x,
	    // save nu st y,
	    // save nu st z,
	    // save interaction type
	    // save interaction mode
        }
     } // save nu info
     
     fn_g4=ptList.size();
     
     for(auto const& pPart : ptList){
         fg4_pdg.push_back(pPart->PdgCode());
	 fg4_st_process.push_back(pPart->Process());
	 
	 TVector3 mcstart, mcend;
         double plen = length(*pPart, mcstart, mcend);
	 fg4_plen.push_back(plen);
	 // Save trackID
	 // Save initial Process name
	 // Save end process name
	 // Save number of daughters
	 // Save Start X,Y,Z coordinates
	 // Save End X,Y,Z coordinates
	 // Save pathlength
	 // Save Momentum and Px,Py, and Pz commponents
	 // Save origin of the particle
	 // Save the statuscode
	 // Save start energy of the particle
	 // Save end energy of the particle
	 // Save the energy deposition inside detector as seen by different planes
	 
	 
	 art::Ptr<simb::MCTruth> truth=inventory_service->TrackIdToMCTruth_P(pPart->TrackId());   
	 fg4_nu_origin.push_back(truth->NeutrinoSet());
	 
	 if(fuse_run_genie){
	    int ccnc=-1;
	    int mode=-1;
	    double nu_x=-9999; double nu_y=-9999; double nu_z=-9999;
	    
	    if(fg4_nu_origin.back() && pPart->Process()=="primary"){
	       const simb::MCNeutrino & my_neutrino = truth->GetNeutrino();
	       ccnc = my_neutrino.CCNC(); //0=CC 1=NC
	       mode = my_neutrino.Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
	       nu_x = my_neutrino.Nu().Vx();
	       nu_y = my_neutrino.Nu().Vy();
	       nu_z = my_neutrino.Nu().Vz();
	       // save neutrino pdg;
	       // save neutrnio momentums PX, Py, and Pz
	       // Save energy of the neutrino
	    }
	    
	    fg4_nu_ccnc.push_back(ccnc);
	    fg4_nu_mode.push_back(mode);
	    fg4_nu_stx.push_back(nu_x);
	    fg4_nu_sty.push_back(nu_y);
	    fg4_nu_stz.push_back(nu_z);
	    
        } // saving genie information
     } // loop over ptlist
     
     ////////////////////////////////////// Reconstructed hits ///////////////////////////
     
     fn_hits = hitlist.size();
     
     for(auto const& hit : hitlist){
         fhit_channel.push_back(hit->Channel());
	 fhit_wire.push_back(hit->WireID().Wire);
	 fhit_plane.push_back(hit->WireID().Plane);
	 fhit_tpc.push_back(hit->WireID().TPC);
	 fhit_cryostat.push_back(hit->WireID().Cryostat);
     } // loop over recob hits
     
     
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
     fn_nu=-9999;
     fgen_nu_set.clear();
     fgen_n_part.clear();
     fgen_nu_pdg.clear();
     fn_g4=-9999;
     fg4_pdg.clear();
     fg4_st_process.clear();
     fg4_plen.clear();
     fg4_nu_origin.clear();
     fg4_nu_ccnc.clear();
     fg4_nu_mode.clear();
     fg4_nu_stx.clear();
     fg4_nu_sty.clear();
     fg4_nu_stz.clear();
     fn_hits=-9999;
     fhit_channel.clear();
     fhit_wire.clear();
     fhit_plane.clear();
     fhit_tpc.clear();
     fhit_cryostat.clear();
}
////////////////////////////////////////////////////////////////////////////////////
DEFINE_ART_MODULE(SingleMuonInfo)
}		
		
		
