////////////////////////////////////////////////////////////////////////
// Class:       ImageMaker
// Plugin Type: analyzer (Unknown Unknown)
// File:        ImageMaker_module.cc
//
// Generated at Fri Mar 25 16:08:24 2022 by Varuna Crishan Meddage using cetskelgen
// from  version .
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
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TPolyMarker.h"
#include "TStyle.h"
#include "TTree.h"

#include <string>

class ImageMaker;


class ImageMaker : public art::EDAnalyzer {
public:
  explicit ImageMaker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ImageMaker(ImageMaker const&) = delete;
  ImageMaker(ImageMaker&&) = delete;
  ImageMaker& operator=(ImageMaker const&) = delete;
  ImageMaker& operator=(ImageMaker&&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const& e) override;
  void DrawEVD(art::Event const& evt, int plID, int tpcID, int SLID, const std::vector<art::Ptr<recob::Hit> >& cat_1, const std::vector<art::Ptr<recob::Hit> >& cat_2, const std::vector<art::Ptr<recob::Hit> >& cat_3, const std::vector<art::Ptr<recob::Hit>>& cat_4, const std::vector<art::Ptr<recob::Hit> >& cat_5);
  void DrawEVD(art::Event const& evt, int tpcID, int SLID, const std::vector<art::Ptr<recob::Hit> >& hit_vec_pl0, const std::vector<art::Ptr<recob::Hit> >& hit_vec_pl1, const std::vector<art::Ptr<recob::Hit> >& hit_vec_pl2);
  void DrawEVD(art::Event const& evt, int tpcID, int SLID, const char* name, const std::vector<std::vector<art::Ptr<recob::Hit>>>& hit_vec_pl0, const std::vector<std::vector<art::Ptr<recob::Hit>>>&
		  hit_vec_pl1, const std::vector<std::vector<art::Ptr<recob::Hit>>>& hit_vec_pl2);
  bool Is_NuLike(const std::vector<art::Ptr<recob::PFParticle>>& pfp_vec);
  std::vector<art::Ptr<recob::Track>> Get_trks_frm_PFP(const std::vector<art::Ptr<recob::PFParticle>>& pfp_vector, const art::FindManyP<recob::Track>& trk_pfp_assn);
  std::vector<art::Ptr<recob::Shower>> Get_shws_frm_PFP(const std::vector<art::Ptr<recob::PFParticle>>& pfp_vector, const art::FindManyP<recob::Shower>& shw_pfp_assn);
  const simb::MCParticle* Back_tracking_trks(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Track> trk, const art::FindManyP<recob::Hit>& hit_trk_assn, std::map<int, const simb::MCParticle*>& particle_map);
  std::pair<const simb::MCParticle*, float> Back_tracking_shws(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Shower> shw, const art::FindManyP<recob::Hit>& hit_shw_assn, std::map<int,
		  std::vector<int> >& mom_id_map, std::map<int, const simb::MCParticle*>& particle_map);
  
  void Clear();

private:
  // Declare member data here.
  TTree *fEventTree;
  std::vector<int> frun;
  std::vector<int> fsubrun;
  std::vector<int> fevent;
  std::vector<int> fsliceID;	
  std::vector<bool> fIs_nu;
  
  bool fDrawEVD;	
		
  art::InputTag fSliceLabel;
  art::InputTag fPFPLabel;
  art::InputTag fTrackLabel;
  art::InputTag fShowerLabel;		
};


ImageMaker::ImageMaker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
fDrawEVD(p.get<bool>("DrawEVD")),
fSliceLabel(p.get<art::InputTag>("SliceLabel")),
fPFPLabel(p.get<art::InputTag>("PFPLabel")),
fTrackLabel(p.get<art::InputTag>("TrackLabel")),
fShowerLabel(p.get<art::InputTag>("ShowerLabel"))				
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ImageMaker::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("SliceTree", "Analyser Output Tree");
  fEventTree->Branch("run", &frun);
  fEventTree->Branch("subrun", &fsubrun);
  fEventTree->Branch("event", &fevent);
  fEventTree->Branch("sliceID", &fsliceID);
  fEventTree->Branch("Is_nu", &fIs_nu);
}

void ImageMaker::analyze(art::Event const& e)
{
 Clear();
 std::cout << "********* " << e.run() << "  " << e.subRun() << "  " << e.id().event() << "  **************\n";	
	
  // Implementation of required member function here.
 art::Handle< std::vector<recob::Slice> > SliceListHandle;
 std::vector< art::Ptr<recob::Slice> > SliceList;
 if( e.getByLabel(fSliceLabel,SliceListHandle))
     art::fill_ptr_vector(SliceList,SliceListHandle);
 
 art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
 std::vector< art::Ptr<recob::PFParticle> > PFPList;
 if( e.getByLabel(fPFPLabel,PFPListHandle))
     art::fill_ptr_vector(PFPList,PFPListHandle);
 
 art::Handle< std::vector<recob::Track> > TrackListHandle;
 std::vector< art::Ptr<recob::Track> > TrackList;
 if( e.getByLabel(fTrackLabel,TrackListHandle))
     art::fill_ptr_vector(TrackList,TrackListHandle);
 
 art::Handle< std::vector<recob::Shower> > ShowerListHandle;
 std::vector< art::Ptr<recob::Shower> > ShowerList;
 if( e.getByLabel(fShowerLabel,ShowerListHandle))
     art::fill_ptr_vector(ShowerList,ShowerListHandle);
 
 art::FindManyP<recob::Hit> findManyHits(SliceListHandle, e, fSliceLabel);
 art::FindManyP<recob::PFParticle> findManyPFPs(SliceListHandle, e, fSliceLabel);
 art::FindManyP<recob::Track> findManyTrks(PFPListHandle, e, fTrackLabel);
 art::FindManyP<recob::Shower> findManyShwrs(PFPListHandle, e, fShowerLabel);
 art::FindManyP<recob::Hit> findManyTrkHits(TrackListHandle, e, fTrackLabel);
 art::FindManyP<recob::Hit> findManyShwHits(ShowerListHandle, e, fShowerLabel);
 auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
 art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
 
 std::map<int, const simb::MCParticle*> trueparticlesMap;
 const sim::ParticleList& particles = particleInventory->ParticleList();
 for(auto const& particleIter : particles){
     const simb::MCParticle *particle = particleIter.second;
     trueparticlesMap[particle->TrackId()] = particle;
 }
 
 std::map<int, std::vector<int> > MotherTrackIdMap;
 for(auto const& mapiter : trueparticlesMap){
     const simb::MCParticle *particle = mapiter.second;
     unsigned int particlePdg = std::abs(particle->PdgCode());
     if(particlePdg != 11 && particlePdg != 22) continue;
     if(trueparticlesMap.find(particle->Mother()) != trueparticlesMap.end()){
        if(std::abs(trueparticlesMap[particle->Mother()]->PdgCode()) == 11 || std::abs(trueparticlesMap[particle->Mother()]->PdgCode()) == 22){
	   int MotherId = particle->Mother();
           int part_temp = MotherId;
	   while(MotherId != 0){
                 part_temp = trueparticlesMap[MotherId]->Mother();
                 if(trueparticlesMap.find(part_temp) == trueparticlesMap.end()) break;
                 if(std::abs(trueparticlesMap[part_temp]->PdgCode()) != 11 && std::abs(trueparticlesMap[part_temp]->PdgCode()) != 22) break;
                 MotherId = part_temp;
           }//close the while loop
	   MotherTrackIdMap[MotherId].push_back(particle->TrackId());
           continue;
        } // select pdg == 11 && pdg == 22
     } // mother is in the map
     MotherTrackIdMap[particle->TrackId()].push_back(particle->TrackId());
 }
 
 //std::cout << "******************** Just before the slice loop **************\n";
 
 for(auto const& slice : SliceList){
     frun.push_back(e.run());
     fsubrun.push_back(e.subRun());;
     fevent.push_back(e.id().event());
     fsliceID.push_back(slice.key());
     
     std::vector<art::Ptr<recob::PFParticle>> slicepfps = findManyPFPs.at(slice.key());
     fIs_nu.push_back(Is_NuLike(slicepfps));
     
     /*if(findManyTrks.isValid()){
	std::vector<art::Ptr<recob::Track>> slice_trks_vector = Get_trks_frm_PFP(slicepfps, findManyTrks);
	std::cout << "@@@@@@@@@@@@ Number of tracks associated with the slice " << slice.key() << "  " << slice_trks_vector.size() << " @@@@@@@@@@@@@@@@@@@@\n";
	if(findManyTrkHits.isValid()){
	   for(auto const& trk : slice_trks_vector){
	       const simb::MCParticle* paticle = Back_tracking_trks(clockData, trk, findManyTrkHits, trueparticlesMap);
	       if(paticle == nullptr) {std::cout << "Reco track ID : " << trk.key() << " has a null pointer associated\n";}
	       else std::cout << "Reco track ID : " << trk.key() << "  MC paricle PDG : " << paticle->PdgCode() << "\n";
	   }
	}
     }*/
      // trk-pfp association is valid
     
     if(findManyShwrs.isValid()){
	std::vector<art::Ptr<recob::Shower>> slice_shwrs_vector = Get_shws_frm_PFP(slicepfps, findManyShwrs);
	std::cout << "@@@@@@@@@@@@@@ Number of showers associated with the slice " << slice.key() << "  " << slice_shwrs_vector.size() << " @@@@@@@@@@@@@@@@@@@@\n";
	if(findManyShwHits.isValid()){
	   for(auto const& shw : slice_shwrs_vector){
	       auto const& my_pair = Back_tracking_shws(clockData, shw, findManyShwHits, MotherTrackIdMap, trueparticlesMap);
	       if(my_pair.first != nullptr)std::cout << "Reco shower ID : " << shw.key() << "  True trk ID : " << my_pair.first->TrackId() << "  True PDG : " << my_pair.first->PdgCode() << "Energy frac : " << my_pair.second << "\n";
	       else std::cout << "*********** Null ptr in the list " << " Energy frac : " << my_pair.second << "\n";
	   }
        }
     } // shw-pfp association is valid
     
     int n_total_tpc_0 = 0; int n_total_tpc_1 = 0;		     		     
     
     ////////////////////////////////////////////////////// Draw the event display ///////////////////////////////////////////////////////////////////////////////////////////		     
     
     if(fDrawEVD){
        std::vector<art::Ptr<recob::Hit>> slicehits = findManyHits.at(slice.key());
        for(auto const& hit : slicehits){
	    if(hit->WireID().TPC == 0)  n_total_tpc_0++;
	    else if(hit->WireID().TPC == 1) n_total_tpc_1++;
        }
     
        int TPC_ID = -9999;
     
        if(n_total_tpc_0 > n_total_tpc_1) TPC_ID = 0;
        else if(n_total_tpc_0 < n_total_tpc_1) TPC_ID = 1;
        else{
           if(n_total_tpc_0 != 0) TPC_ID = 0;
        }
     
        std::vector<art::Ptr<recob::Hit> > plane_0_hit_charge_less_200,plane_0_hit_charge_200_400,plane_0_hit_charge_400_600,plane_0_hit_charge_600_800,plane_0_hit_charge_over_800;
        std::vector<art::Ptr<recob::Hit> > plane_1_hit_charge_less_200,plane_1_hit_charge_200_400,plane_1_hit_charge_400_600,plane_1_hit_charge_600_800,plane_1_hit_charge_over_800; 
        std::vector<art::Ptr<recob::Hit> > plane_2_hit_charge_less_200,plane_2_hit_charge_200_400,plane_2_hit_charge_400_600,plane_2_hit_charge_600_800,plane_2_hit_charge_over_800;  
        std::vector<art::Ptr<recob::Hit> > plane_0_hits,plane_1_hits,plane_2_hits; 
        std::vector<std::vector<art::Ptr<recob::Hit>>> pl0_vectors,pl1_vectors,pl2_vectors;
        pl0_vectors.resize(5);pl1_vectors.resize(5);pl2_vectors.resize(5);  
     
       for(auto const& hit : slicehits){
            if(int(hit->WireID().TPC) == TPC_ID){
	    
	       if(hit->WireID().Plane == 0){
	          plane_0_hits.push_back(hit);
	          if(hit->Integral() < 200) plane_0_hit_charge_less_200.push_back(hit);
	          else if(hit->Integral() >= 200 && hit->Integral() < 400) plane_0_hit_charge_200_400.push_back(hit);  
	          else if(hit->Integral() >= 400 && hit->Integral() < 600) plane_0_hit_charge_400_600.push_back(hit);
	          else if(hit->Integral() >= 600 && hit->Integral() < 800) plane_0_hit_charge_600_800.push_back(hit);
	          else if(hit->Integral() >= 800) plane_0_hit_charge_over_800.push_back(hit);
	       }
	    
	       else if(hit->WireID().Plane == 1){
		    plane_1_hits.push_back(hit);
	            if(hit->Integral() < 200) plane_1_hit_charge_less_200.push_back(hit); 
	            else if(hit->Integral() >= 200 && hit->Integral() < 400) plane_1_hit_charge_200_400.push_back(hit);  
	            else if(hit->Integral() >= 400 && hit->Integral() < 600) plane_1_hit_charge_400_600.push_back(hit);
	            else if(hit->Integral() >= 600 && hit->Integral() < 800) plane_1_hit_charge_600_800.push_back(hit);
	            else if(hit->Integral() >= 800) plane_1_hit_charge_over_800.push_back(hit);  
	       }
	    
	       else if(hit->WireID().Plane == 2){
		    plane_2_hits.push_back(hit);
	            if(hit->Integral() < 200) plane_2_hit_charge_less_200.push_back(hit);
	            else if(hit->Integral() >= 200 && hit->Integral() < 400) plane_2_hit_charge_200_400.push_back(hit); 
	            else if(hit->Integral() >= 400 && hit->Integral() < 600) plane_2_hit_charge_400_600.push_back(hit);
	            else if(hit->Integral() >= 600 && hit->Integral() < 800) plane_2_hit_charge_600_800.push_back(hit);
	            else if(hit->Integral() >= 800) plane_2_hit_charge_over_800.push_back(hit);   
	       }
	    }
	 
	    pl0_vectors[0] = plane_0_hit_charge_less_200;
	    pl0_vectors[1] = plane_0_hit_charge_200_400; 
	    pl0_vectors[2] = plane_0_hit_charge_400_600;
	    pl0_vectors[3] = plane_0_hit_charge_600_800;
	    pl0_vectors[4] = plane_0_hit_charge_over_800;
	 
	    pl1_vectors[0] = plane_1_hit_charge_less_200;
	    pl1_vectors[1] = plane_1_hit_charge_200_400; 
	    pl1_vectors[2] = plane_1_hit_charge_400_600;
	    pl1_vectors[3] = plane_1_hit_charge_600_800;
	    pl1_vectors[4] = plane_1_hit_charge_over_800;
	 
	    pl2_vectors[0] = plane_2_hit_charge_less_200;
	    pl2_vectors[1] = plane_2_hit_charge_200_400; 
	    pl2_vectors[2] = plane_2_hit_charge_400_600;
	    pl2_vectors[3] = plane_2_hit_charge_600_800;
	    pl2_vectors[4] = plane_2_hit_charge_over_800;
        } // loop over slice hits
        std::string name = "cos-noise";
        if(Is_NuLike(slicepfps)) name = "nu";
        DrawEVD(e, TPC_ID, slice.key(), name.c_str(), pl0_vectors, pl1_vectors, pl2_vectors);
      }
      
      //////////////////////////////////////////////////////////////////// End of event display //////////////////////////////////////////////
     	     		     	     
 } // loop over slice list
 fEventTree->Fill();
 //std::cout << "************* End of the analyze function ****************\n";
} // End of analyzer module

//////////////////////////////////////////////////////////////////////////////////////////////

void ImageMaker::DrawEVD(art::Event const& evt, int plID, int tpcID, int SLID, const std::vector<art::Ptr<recob::Hit> >& cat_1, const std::vector<art::Ptr<recob::Hit> >& cat_2, const
		std::vector<art::Ptr<recob::Hit> >& cat_3, const std::vector<art::Ptr<recob::Hit>> & cat_4, const std::vector<art::Ptr<recob::Hit> >& cat_5){
     art::ServiceHandle<art::TFileService> tfs;
     TString canvasName = Form("Display_%i_%i_%i_%i_%i_%i",evt.run(), evt.subRun(), evt.event(), tpcID, plID, SLID);
     TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);
     //canvas->Divide(1,3,0,0.1);
     int WireID = 0;
     double PeakTime = 0;
     int x_min = std::numeric_limits<int>::max();
     int x_max = -std::numeric_limits<int>::max();
     double y_min = std::numeric_limits<double>::max();
     double y_max = -std::numeric_limits<double>::max();
     
     {
       int point = 0;
       auto cat_1_hits_Poly = std::make_unique<TPolyMarker>(cat_1.size());
       if(cat_1.size() !=0){
          for(auto const& cat_1_hit : cat_1){
	      int PlaneID = cat_1_hit->WireID().Plane;
	      if(PlaneID != plID) continue;
	      WireID   = cat_1_hit->WireID().Wire;
              PeakTime = cat_1_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      cat_1_hits_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over catogory 1 hit vector
       } // cat 1 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto cat_2_hits_Poly = std::make_unique<TPolyMarker>(cat_2.size());
       if(cat_2.size() !=0){
          for(auto const& cat_2_hit : cat_2){
	      int PlaneID = cat_2_hit->WireID().Plane;
	      if(PlaneID != plID) continue;
	      WireID = cat_2_hit->WireID().Wire;
              PeakTime = cat_2_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      cat_2_hits_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over catagory 2 hit vector
       } // cat 2 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0; 
       auto cat_3_hits_Poly = std::make_unique<TPolyMarker>(cat_3.size());
       if(cat_3.size() !=0){
	  for(auto const& cat_3_hit : cat_3){
	      int PlaneID = cat_3_hit->WireID().Plane;
	      if(PlaneID != plID) continue;
	      WireID = cat_3_hit->WireID().Wire;
              PeakTime = cat_3_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      cat_3_hits_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over catagory 3 hit vector
       } // cat 3 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0; 
       auto cat_4_hits_Poly = std::make_unique<TPolyMarker>(cat_4.size());
       if(cat_4.size() !=0){
	  for(auto const& cat_4_hit : cat_4){
	      int PlaneID = cat_4_hit->WireID().Plane;
	      if(PlaneID != plID) continue;
	      WireID = cat_4_hit->WireID().Wire;
              PeakTime = cat_4_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      cat_4_hits_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over catagory 4 hit vector
       } // cat 4 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0; 
       auto cat_5_hits_Poly = std::make_unique<TPolyMarker>(cat_5.size());
       if(cat_5.size() !=0){
	  for(auto const& cat_5_hit : cat_5){
	      int PlaneID = cat_5_hit->WireID().Plane;
	      if(PlaneID != plID) continue;
	      WireID = cat_5_hit->WireID().Wire;
              PeakTime = cat_5_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      cat_5_hits_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over catagory 5 hit vector
       } // cat 5 hit vector size is non zero
       
       //************************************************************************
       //canvas->cd(1);
       gStyle->SetOptStat(0);
       int new_xmin = std::max(1,x_min-25);
       double new_ymin = std::max(0.,y_min-25);
       double new_ymax = std::min(7680.,y_max+25);
       TH2F axes("axes","",1,new_xmin,x_max+25,1,new_ymin,new_ymax);
       axes.SetDirectory(0);
       axes.GetXaxis()->SetTitle(" WireID");
       axes.GetYaxis()->SetTitle(" time (tick)");
       axes.Draw();
       if(cat_1.size() !=0){
          cat_1_hits_Poly->SetMarkerStyle(20);
          cat_1_hits_Poly->SetMarkerSize(0.75);
          cat_1_hits_Poly->SetMarkerColor(8);
          cat_1_hits_Poly->Draw();
       }
       if(cat_2.size() !=0){
          cat_2_hits_Poly->SetMarkerStyle(20);
          cat_2_hits_Poly->SetMarkerSize(0.75);
          cat_2_hits_Poly->SetMarkerColor(4);
          cat_2_hits_Poly->Draw();
       }
       if(cat_3.size() !=0){
          cat_3_hits_Poly->SetMarkerStyle(20);
          cat_3_hits_Poly->SetMarkerSize(0.75);
          cat_3_hits_Poly->SetMarkerColor(6);
          cat_3_hits_Poly->Draw();
       }
       if(cat_4.size() !=0){
          cat_4_hits_Poly->SetMarkerStyle(20);
          cat_4_hits_Poly->SetMarkerSize(0.75);
          cat_4_hits_Poly->SetMarkerColor(5);
          cat_4_hits_Poly->Draw();
       }
       if(cat_5.size() !=0){
          cat_5_hits_Poly->SetMarkerStyle(20);
          cat_5_hits_Poly->SetMarkerSize(0.75);
          cat_5_hits_Poly->SetMarkerColor(2);
          cat_5_hits_Poly->Draw();
       }
       
       canvas->Write();
       canvas->SaveAs(canvasName + "_MC_reco_hits.png");
    } // close the block
  }

/////////////////////////////////////////////////////////////////////////////////////////////
  
void ImageMaker::DrawEVD(art::Event const& evt, int tpcID, int SLID, const std::vector<art::Ptr<recob::Hit> >& hit_vec_pl0, const std::vector<art::Ptr<recob::Hit> >& hit_vec_pl1, const std::vector<art::Ptr<recob::Hit> >& hit_vec_pl2){
     art::ServiceHandle<art::TFileService> tfs;
     TString canvasName = Form("AllHitsDisplay_%i_%i_%i_%i_%i",evt.run(), evt.subRun(), evt.event(), tpcID, SLID);
     TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);
     canvas->Divide(1,3,0,0.1);
     int WireID = 0;
     double PeakTime = 0;
     int x_min = std::numeric_limits<int>::max();
     int x_max = -std::numeric_limits<int>::max();
     double y_min = std::numeric_limits<double>::max();
     double y_max = -std::numeric_limits<double>::max();
     
     {
      int point = 0;
      auto all_hit_Poly2 = std::make_unique<TPolyMarker>(hit_vec_pl2.size());
      if(hit_vec_pl2.size() !=0){
         for(auto const& hit : hit_vec_pl2){
	      int PlaneID = hit->WireID().Plane;
	      if(PlaneID != 2) continue;
	      WireID   = hit->WireID().Wire;
              PeakTime = hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      all_hit_Poly2->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over  hit vector
       } // hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto all_hit_Poly1 = std::make_unique<TPolyMarker>(hit_vec_pl1.size());
       if(hit_vec_pl1.size() !=0){
         for(auto const& hit : hit_vec_pl1){
	      int PlaneID = hit->WireID().Plane;
	      if(PlaneID != 1) continue;
	      WireID   = hit->WireID().Wire;
              PeakTime = hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      all_hit_Poly1->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over  hit vector
       } // hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto all_hit_Poly0 = std::make_unique<TPolyMarker>(hit_vec_pl0.size());
       if(hit_vec_pl0.size() !=0){
         for(auto const& hit : hit_vec_pl0){
	      int PlaneID = hit->WireID().Plane;
	      if(PlaneID != 0) continue;
	      WireID   = hit->WireID().Wire;
              PeakTime = hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      all_hit_Poly0->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over  hit vector
       } // hit vector size is non zero
       
       //****************************************************************
       
       gStyle->SetOptStat(0);
       int new_xmin = std::max(1,x_min-25);
       double new_ymin = std::max(0.,y_min-25);
       double new_ymax = std::min(7680.,y_max+25);
       TH2F axes("axes","",1,new_xmin,x_max+25,1,new_ymin,new_ymax);
       axes.SetDirectory(0);
       axes.GetXaxis()->SetLabelSize(0.045);
       axes.GetYaxis()->SetLabelSize(0.045);
       canvas->cd(1);
       axes.Draw();
       if(hit_vec_pl2.size() !=0){
          all_hit_Poly2->SetMarkerStyle(20);
          all_hit_Poly2->SetMarkerSize(0.50);
          all_hit_Poly2->SetMarkerColor(2);
          all_hit_Poly2->Draw();
       }
       canvas->cd(2);
       axes.Draw();
       if(hit_vec_pl1.size() !=0){
          all_hit_Poly1->SetMarkerStyle(20);
          all_hit_Poly1->SetMarkerSize(0.50);
          all_hit_Poly1->SetMarkerColor(2);
          all_hit_Poly1->Draw();
       }
       canvas->cd(3);
       axes.Draw();
       if(hit_vec_pl0.size() !=0){
          all_hit_Poly0->SetMarkerStyle(20);
          all_hit_Poly0->SetMarkerSize(0.50);
          all_hit_Poly0->SetMarkerColor(2);
          all_hit_Poly0->Draw();
       }
       canvas->Write();
       canvas->SaveAs(canvasName + "_MC_reco_hits.png");
     }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ImageMaker::DrawEVD(art::Event const& evt, int tpcID, int SLID, const char* name, const std::vector<std::vector<art::Ptr<recob::Hit>>>& hit_vec_pl0, const std::vector<std::vector<art::Ptr<recob::Hit>>>& hit_vec_pl1, const std::vector<std::vector<art::Ptr<recob::Hit>>>& hit_vec_pl2){
     art::ServiceHandle<art::TFileService> tfs;
     TString canvasName = Form("Display_%i_%i_%i_%i_%i_%s",evt.run(), evt.subRun(), evt.event(), tpcID, SLID, name);
     TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);
     canvas->Divide(1,3,0,0.1);
     int WireID = 0;
     double PeakTime = 0;
     int x_min = std::numeric_limits<int>::max();
     int x_max = -std::numeric_limits<int>::max();
     double y_min = std::numeric_limits<double>::max();
     double y_max = -std::numeric_limits<double>::max();
     {
      int point = 0;
      auto pl_2_cat_1_Poly = std::make_unique<TPolyMarker>(hit_vec_pl2[0].size());
      if(hit_vec_pl2[0].size() !=0){
         for(auto const& cat_1_hit : hit_vec_pl2[0]){
	     int PlaneID = cat_1_hit->WireID().Plane;
	     if(PlaneID != 2) continue;
	     WireID   = cat_1_hit->WireID().Wire;
             PeakTime = cat_1_hit->PeakTime();
	     x_min = std::min(WireID,x_min);
             x_max = std::max(WireID,x_max);
             y_min = std::min(PeakTime,y_min);
             y_max = std::max(PeakTime,y_max);
	     pl_2_cat_1_Poly->SetPoint(point,WireID,PeakTime);
             ++point;
	  } // loop over plane 2 catogory 1 hit vector
      } // plane 2 cat 1 hit vector size is non zero
       
      point = 0; 
      WireID = 0;
      PeakTime = 0;
      auto pl_2_cat_2_Poly = std::make_unique<TPolyMarker>(hit_vec_pl2[1].size());
      if(hit_vec_pl2[1].size() !=0){
         for(auto const& cat_2_hit : hit_vec_pl2[1]){
	     int PlaneID = cat_2_hit->WireID().Plane;
	     if(PlaneID != 2) continue;
	     WireID = cat_2_hit->WireID().Wire;
             PeakTime = cat_2_hit->PeakTime();
	     x_min = std::min(WireID,x_min);
             x_max = std::max(WireID,x_max);
             y_min = std::min(PeakTime,y_min);
             y_max = std::max(PeakTime,y_max);
	     pl_2_cat_2_Poly->SetPoint(point,WireID,PeakTime);
             ++point;
	  } // loop over plane 2 catogory 2 hit vector
       } // plane 2 cat 2 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto pl_2_cat_3_Poly = std::make_unique<TPolyMarker>(hit_vec_pl2[2].size());
       if(hit_vec_pl2[2].size() !=0){
         for(auto const& cat_3_hit : hit_vec_pl2[2]){
	     int PlaneID = cat_3_hit->WireID().Plane;
	     if(PlaneID != 2) continue;
	     WireID = cat_3_hit->WireID().Wire;
             PeakTime = cat_3_hit->PeakTime();
	     x_min = std::min(WireID,x_min);
             x_max = std::max(WireID,x_max);
             y_min = std::min(PeakTime,y_min);
             y_max = std::max(PeakTime,y_max);
	     pl_2_cat_3_Poly->SetPoint(point,WireID,PeakTime);
             ++point;
	  } // loop over plane 2 catogory 3 hit vector
       } // plane 2 cat 3 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto pl_2_cat_4_Poly = std::make_unique<TPolyMarker>(hit_vec_pl2[3].size());
       if(hit_vec_pl2[3].size() !=0){
         for(auto const& cat_4_hit : hit_vec_pl2[3]){
	     int PlaneID = cat_4_hit->WireID().Plane;
	     if(PlaneID != 2) continue;
	     WireID = cat_4_hit->WireID().Wire;
             PeakTime = cat_4_hit->PeakTime();
	     x_min = std::min(WireID,x_min);
             x_max = std::max(WireID,x_max);
             y_min = std::min(PeakTime,y_min);
             y_max = std::max(PeakTime,y_max);
	     pl_2_cat_4_Poly->SetPoint(point,WireID,PeakTime);
             ++point;
	  } // loop over plane 2 catogory 4 hit vector
       } // plane 2 cat 4 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto pl_2_cat_5_Poly = std::make_unique<TPolyMarker>(hit_vec_pl2[4].size());
       if(hit_vec_pl2[4].size() !=0){
         for(auto const& cat_5_hit : hit_vec_pl2[4]){
	     int PlaneID = cat_5_hit->WireID().Plane;
	     if(PlaneID != 2) continue;
	     WireID = cat_5_hit->WireID().Wire;
             PeakTime = cat_5_hit->PeakTime();
	     x_min = std::min(WireID,x_min);
             x_max = std::max(WireID,x_max);
             y_min = std::min(PeakTime,y_min);
             y_max = std::max(PeakTime,y_max);
	     pl_2_cat_5_Poly->SetPoint(point,WireID,PeakTime);
             ++point;
	  } // loop over plane 2 catogory 5 hit vector
       } // plane 2 cat 5 hit vector size is non zero
       
       //********************* Plotting plane 1 information ***************************
       
       point = 0;
       WireID = 0;
       PeakTime = 0;
       auto pl_1_cat_1_Poly = std::make_unique<TPolyMarker>(hit_vec_pl1[0].size());
       if(hit_vec_pl1[0].size() !=0){
          for(auto const& cat_1_hit : hit_vec_pl1[0]){
	      int PlaneID = cat_1_hit->WireID().Plane;
	      if(PlaneID != 1) continue;
	      WireID   = cat_1_hit->WireID().Wire;
              PeakTime = cat_1_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      pl_1_cat_1_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over plane 1 catogory 1 hit vector
       } // plane 1 cat 1 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto pl_1_cat_2_Poly = std::make_unique<TPolyMarker>(hit_vec_pl1[1].size());
       if(hit_vec_pl1[1].size() !=0){
          for(auto const& cat_2_hit : hit_vec_pl1[1]){
	      int PlaneID = cat_2_hit->WireID().Plane;
	      if(PlaneID != 1) continue;
	      WireID = cat_2_hit->WireID().Wire;
              PeakTime = cat_2_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      pl_1_cat_2_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	   } // loop over plane 1 catogory 2 hit vector
        } // plane 1 cat 2 hit vector size is non zero
       
        point = 0; 
        WireID = 0;
        PeakTime = 0;
        auto pl_1_cat_3_Poly = std::make_unique<TPolyMarker>(hit_vec_pl1[2].size());
        if(hit_vec_pl1[2].size() !=0){
           for(auto const& cat_3_hit : hit_vec_pl1[2]){
	       int PlaneID = cat_3_hit->WireID().Plane;
	       if(PlaneID != 1) continue;
	       WireID = cat_3_hit->WireID().Wire;
               PeakTime = cat_3_hit->PeakTime();
	       x_min = std::min(WireID,x_min);
               x_max = std::max(WireID,x_max);
               y_min = std::min(PeakTime,y_min);
               y_max = std::max(PeakTime,y_max);
	       pl_1_cat_3_Poly->SetPoint(point,WireID,PeakTime);
               ++point;
	   } // loop over plane 1 catogory 3 hit vector
        } // plane 1 cat 3 hit vector size is non zero
       
        point = 0; 
        WireID = 0;
        PeakTime = 0;
        auto pl_1_cat_4_Poly = std::make_unique<TPolyMarker>(hit_vec_pl1[3].size());
        if(hit_vec_pl1[3].size() !=0){
           for(auto const& cat_4_hit : hit_vec_pl1[3]){
	       int PlaneID = cat_4_hit->WireID().Plane;
	       if(PlaneID != 1) continue;
	       WireID = cat_4_hit->WireID().Wire;
               PeakTime = cat_4_hit->PeakTime();
	       x_min = std::min(WireID,x_min);
               x_max = std::max(WireID,x_max);
               y_min = std::min(PeakTime,y_min);
               y_max = std::max(PeakTime,y_max);
	       pl_1_cat_4_Poly->SetPoint(point,WireID,PeakTime);
               ++point;
	   } // loop over plane 1 catogory 4 hit vector
        } // plane 1 cat 4 hit vector size is non zero
       
        point = 0; 
        WireID = 0;
        PeakTime = 0;
        auto pl_1_cat_5_Poly = std::make_unique<TPolyMarker>(hit_vec_pl1[4].size());
        if(hit_vec_pl1[4].size() !=0){
          for(auto const& cat_5_hit : hit_vec_pl1[4]){
	      int PlaneID = cat_5_hit->WireID().Plane;
	      if(PlaneID != 1) continue;
	      WireID = cat_5_hit->WireID().Wire;
              PeakTime = cat_5_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      pl_1_cat_5_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	   } // loop over plane 1 catogory 5 hit vector
        } // plane 1 cat 5 hit vector size is non zero
       
       //******************** Closing plotting palne 1 information ********************
	
      //********************* Plotting plane 0 information ****************************
       point = 0;
       WireID = 0;
       PeakTime = 0;
       auto pl_0_cat_1_Poly = std::make_unique<TPolyMarker>(hit_vec_pl0[0].size());
       if(hit_vec_pl0[0].size() !=0){
          for(auto const& cat_1_hit : hit_vec_pl0[0]){
	      int PlaneID = cat_1_hit->WireID().Plane;
	      if(PlaneID != 0) continue;
	      WireID   = cat_1_hit->WireID().Wire;
              PeakTime = cat_1_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      pl_0_cat_1_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	  } // loop over plane 0 catogory 1 hit vector
       } // plane 0 cat 1 hit vector size is non zero
       
       point = 0; 
       WireID = 0;
       PeakTime = 0;
       auto pl_0_cat_2_Poly = std::make_unique<TPolyMarker>(hit_vec_pl0[1].size());
       if(hit_vec_pl0[1].size() !=0){
          for(auto const& cat_2_hit : hit_vec_pl0[1]){
	      int PlaneID = cat_2_hit->WireID().Plane;
	      if(PlaneID != 0) continue;
	      WireID = cat_2_hit->WireID().Wire;
              PeakTime = cat_2_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      pl_0_cat_2_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	   } // loop over plane 0 catogory 2 hit vector
        } // plane 0 cat 2 hit vector size is non zero
       
        point = 0; 
        WireID = 0;
        PeakTime = 0;
        auto pl_0_cat_3_Poly = std::make_unique<TPolyMarker>(hit_vec_pl0[2].size());
        if(hit_vec_pl0[2].size() !=0){
           for(auto const& cat_3_hit : hit_vec_pl0[2]){
	       int PlaneID = cat_3_hit->WireID().Plane;
	       if(PlaneID != 0) continue;
	       WireID = cat_3_hit->WireID().Wire;
               PeakTime = cat_3_hit->PeakTime();
	       x_min = std::min(WireID,x_min);
               x_max = std::max(WireID,x_max);
               y_min = std::min(PeakTime,y_min);
               y_max = std::max(PeakTime,y_max);
	       pl_0_cat_3_Poly->SetPoint(point,WireID,PeakTime);
               ++point;
	   } // loop over plane 0 catogory 3 hit vector
        } // plane 0 cat 3 hit vector size is non zero
       
        point = 0; 
        WireID = 0;
        PeakTime = 0;
        auto pl_0_cat_4_Poly = std::make_unique<TPolyMarker>(hit_vec_pl0[3].size());
        if(hit_vec_pl0[3].size() !=0){
           for(auto const& cat_4_hit : hit_vec_pl0[3]){
	       int PlaneID = cat_4_hit->WireID().Plane;
	       if(PlaneID != 0) continue;
	       WireID = cat_4_hit->WireID().Wire;
               PeakTime = cat_4_hit->PeakTime();
	       x_min = std::min(WireID,x_min);
               x_max = std::max(WireID,x_max);
               y_min = std::min(PeakTime,y_min);
               y_max = std::max(PeakTime,y_max);
	       pl_0_cat_4_Poly->SetPoint(point,WireID,PeakTime);
               ++point;
	   } // loop over plane 0 catogory 4 hit vector
        } // plane 0 cat 4 hit vector size is non zero
       
        point = 0; 
        WireID = 0;
        PeakTime = 0;
        auto pl_0_cat_5_Poly = std::make_unique<TPolyMarker>(hit_vec_pl0[4].size());
        if(hit_vec_pl0[4].size() !=0){
          for(auto const& cat_5_hit : hit_vec_pl0[4]){
	      int PlaneID = cat_5_hit->WireID().Plane;
	      if(PlaneID != 0) continue;
	      WireID = cat_5_hit->WireID().Wire;
              PeakTime = cat_5_hit->PeakTime();
	      x_min = std::min(WireID,x_min);
              x_max = std::max(WireID,x_max);
              y_min = std::min(PeakTime,y_min);
              y_max = std::max(PeakTime,y_max);
	      pl_0_cat_5_Poly->SetPoint(point,WireID,PeakTime);
              ++point;
	   } // loop over plane 0 catogory 5 hit vector
        } // plane 0 cat 5 hit vector size is non zero	
	
	
      //********************* Closing plotting plane 0 information ********************	
       
       gStyle->SetOptStat(0);
       int new_xmin = std::max(1,x_min-25);
       double new_ymin = std::max(0.,y_min-25);
       double new_ymax = std::min(7680.,y_max+25);
       TH2F axes("axes","",1,new_xmin,x_max+25,1,new_ymin,new_ymax);
       axes.SetDirectory(0);
       axes.GetXaxis()->SetLabelSize(0.065);
       axes.GetYaxis()->SetLabelSize(0.065);
       canvas->cd(1);
       axes.Draw();
       if(hit_vec_pl2[0].size() !=0){
          pl_2_cat_1_Poly->SetMarkerStyle(20);
          pl_2_cat_1_Poly->SetMarkerSize(0.30);
          pl_2_cat_1_Poly->SetMarkerColor(8);
          pl_2_cat_1_Poly->Draw();
       }
       if(hit_vec_pl2[1].size() !=0){
          pl_2_cat_2_Poly->SetMarkerStyle(20);
          pl_2_cat_2_Poly->SetMarkerSize(0.30);
          pl_2_cat_2_Poly->SetMarkerColor(4);
          pl_2_cat_2_Poly->Draw();
       }
       if(hit_vec_pl2[2].size() !=0){
          pl_2_cat_3_Poly->SetMarkerStyle(20);
          pl_2_cat_3_Poly->SetMarkerSize(0.30);
          pl_2_cat_3_Poly->SetMarkerColor(6);
          pl_2_cat_3_Poly->Draw();
       }
       if(hit_vec_pl2[3].size() !=0){
          pl_2_cat_4_Poly->SetMarkerStyle(20);
          pl_2_cat_4_Poly->SetMarkerSize(0.30);
          pl_2_cat_4_Poly->SetMarkerColor(5);
          pl_2_cat_4_Poly->Draw();
       }
       if(hit_vec_pl2[4].size() !=0){
          pl_2_cat_5_Poly->SetMarkerStyle(20);
          pl_2_cat_5_Poly->SetMarkerSize(0.30);
          pl_2_cat_5_Poly->SetMarkerColor(2);
          pl_2_cat_5_Poly->Draw();
       }
       canvas->cd(2);
       axes.Draw();
       if(hit_vec_pl1[0].size() !=0){
          pl_1_cat_1_Poly->SetMarkerStyle(20);
          pl_1_cat_1_Poly->SetMarkerSize(0.30);
          pl_1_cat_1_Poly->SetMarkerColor(8);
          pl_1_cat_1_Poly->Draw();
       }
       if(hit_vec_pl1[1].size() !=0){
          pl_1_cat_2_Poly->SetMarkerStyle(20);
          pl_1_cat_2_Poly->SetMarkerSize(0.30);
          pl_1_cat_2_Poly->SetMarkerColor(4);
          pl_1_cat_2_Poly->Draw();
       }
       if(hit_vec_pl1[2].size() !=0){
          pl_1_cat_3_Poly->SetMarkerStyle(20);
          pl_1_cat_3_Poly->SetMarkerSize(0.30);
          pl_1_cat_3_Poly->SetMarkerColor(6);
          pl_1_cat_3_Poly->Draw();
       }
       if(hit_vec_pl1[3].size() !=0){
          pl_1_cat_4_Poly->SetMarkerStyle(20);
          pl_1_cat_4_Poly->SetMarkerSize(0.30);
          pl_1_cat_4_Poly->SetMarkerColor(5);
          pl_1_cat_4_Poly->Draw();
       }
       if(hit_vec_pl1[4].size() !=0){
          pl_1_cat_5_Poly->SetMarkerStyle(20);
          pl_1_cat_5_Poly->SetMarkerSize(0.30);
          pl_1_cat_5_Poly->SetMarkerColor(2);
          pl_1_cat_5_Poly->Draw();
       }
       canvas->cd(3);
       axes.Draw();
       if(hit_vec_pl0[0].size() !=0){
          pl_0_cat_1_Poly->SetMarkerStyle(20);
          pl_0_cat_1_Poly->SetMarkerSize(0.30);
          pl_0_cat_1_Poly->SetMarkerColor(8);
          pl_0_cat_1_Poly->Draw();
       }
       if(hit_vec_pl0[1].size() !=0){
          pl_0_cat_2_Poly->SetMarkerStyle(20);
          pl_0_cat_2_Poly->SetMarkerSize(0.30);
          pl_0_cat_2_Poly->SetMarkerColor(4);
          pl_0_cat_2_Poly->Draw();
       }
       if(hit_vec_pl0[2].size() !=0){
          pl_0_cat_3_Poly->SetMarkerStyle(20);
          pl_0_cat_3_Poly->SetMarkerSize(0.30);
          pl_0_cat_3_Poly->SetMarkerColor(6);
          pl_0_cat_3_Poly->Draw();
       }
       if(hit_vec_pl0[3].size() !=0){
          pl_0_cat_4_Poly->SetMarkerStyle(20);
          pl_0_cat_4_Poly->SetMarkerSize(0.30);
          pl_0_cat_4_Poly->SetMarkerColor(5);
          pl_0_cat_4_Poly->Draw();
       }
       if(hit_vec_pl0[4].size() !=0){
          pl_0_cat_5_Poly->SetMarkerStyle(20);
          pl_0_cat_5_Poly->SetMarkerSize(0.30);
          pl_0_cat_5_Poly->SetMarkerColor(2);
          pl_0_cat_5_Poly->Draw();
       }
       //canvas->Write();
       canvas->SaveAs(canvasName + ".png");
     } // close internal block
}

//////////////////////////////// Is_NuLike //////////////////////////////////////////
bool ImageMaker::Is_NuLike(const std::vector<art::Ptr<recob::PFParticle>>& pfp_vec){
     for(auto const& pfp : pfp_vec){
	 //std::cout << "************ key-self " << pfp.key() << "   " << pfp->Self() << "\n";
         if(pfp->PdgCode() == 12 || pfp->PdgCode() == 14) return true;
     }
     return false;
}
////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Clear /////////////////////////////////////////////
void ImageMaker::Clear(){
     frun.clear();
     fsubrun.clear();
     fevent.clear();
     fsliceID.clear();
     fIs_nu.clear();
}
////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Get_trks_frm_PFP //////////////////////////////////
std::vector<art::Ptr<recob::Track>> ImageMaker::Get_trks_frm_PFP(const std::vector<art::Ptr<recob::PFParticle>>& pfp_vector, const art::FindManyP<recob::Track>& trk_pfp_assn){
     std::vector<art::Ptr<recob::Track>> trks_frm_pfp_vec;
     for(auto const& pfp : pfp_vector){
         std::vector<art::Ptr<recob::Track>> trks = trk_pfp_assn.at(pfp.key());
	 if(trks.size() > 1) std::cout << "********* ImageMaker::Get_trks_frm_PFP::More than one track is associated with the pfparticle " << pfp.key() << "\n";
	 if(!trks.empty() && (std::find(trks_frm_pfp_vec.begin(), trks_frm_pfp_vec.end(), trks[0]) == trks_frm_pfp_vec.end()))trks_frm_pfp_vec.push_back(trks[0]);
     }
     return trks_frm_pfp_vec;
}
////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Get_shws_frm_PFP ///////////////////////////////////
std::vector<art::Ptr<recob::Shower>> ImageMaker::Get_shws_frm_PFP(const std::vector<art::Ptr<recob::PFParticle>>& pfp_vector, const art::FindManyP<recob::Shower>& shw_pfp_assn){
     std::vector<art::Ptr<recob::Shower>> shwrs_frm_pfp_vec;
     for(auto const& pfp : pfp_vector){
         std::vector<art::Ptr<recob::Shower>> shwrs = shw_pfp_assn.at(pfp.key());
	 if(shwrs.size() > 1) std::cout << "************ ImageMaker::Get_trks_frm_PFP::More than one shower is associated with the pfparticle " << pfp.key() << "\n"; 
	 if(!shwrs.empty() && (std::find(shwrs_frm_pfp_vec.begin(), shwrs_frm_pfp_vec.end(), shwrs[0]) == shwrs_frm_pfp_vec.end()))shwrs_frm_pfp_vec.push_back(shwrs[0]);
     }
     return shwrs_frm_pfp_vec;
}
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Back_tracking_trks /////////////////////////////////
const simb::MCParticle* ImageMaker::Back_tracking_trks(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Track> trk, const art::FindManyP<recob::Hit>& hit_trk_assn, std::map<int, const simb::MCParticle*>& particle_map){
     std::vector<art::Ptr<recob::Hit>> hits = hit_trk_assn.at(trk.key());
     std::map<int,double> id_to_energy_map;
     double particleEnergy = -99999;
     int likelyTrackID = 0;
     art::ServiceHandle<cheat::BackTrackerService> bt_serv;
     for(auto const& hit : hits){
	 std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
	 for(unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt){
             double energy = trackIDs.at(idIt).energy;
             int TrackID = trackIDs.at(idIt).trackID;
             id_to_energy_map[TrackID]+=energy;
         }
     }
     
     for(auto const& [id, erg] : id_to_energy_map){
         double particle_contrib_energy = erg;
         if(particle_contrib_energy > particleEnergy){
            particleEnergy = particle_contrib_energy;
            likelyTrackID = id;
         }
     }
     
     const simb::MCParticle* particle(nullptr);
     if(particle_map.find(likelyTrackID) != particle_map.end()){ 
	std::cout << "************* ImageMaker::Back_tracking_trks::Retruning a valid pointer a simb::MCParticle **************\n";
	return particle_map[likelyTrackID];
     }
     else{ 
	std::cout << "********** ImageMaker::Back_tracking_trks::Retruning a null simb::MCParticle pointer **********\n";
	return particle;
     }	 
}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Back_tracking_shws //////////////////////////////////
std::pair<const simb::MCParticle*, float> ImageMaker::Back_tracking_shws(detinfo::DetectorClocksData const& clockData, const art::Ptr<recob::Shower> shw, const art::FindManyP<recob::Hit>& hit_shw_assn, std::map<int,
		std::vector<int> >& mom_id_map, std::map<int, const simb::MCParticle*>& particle_map){
      std::vector<art::Ptr<recob::Hit>> hits = hit_shw_assn.at(shw.key());
      art::ServiceHandle<cheat::BackTrackerService> bt_serv;
      std::map<int, float> mapTrackIdToShowerhitEnergy;
      float tot_shw_eng = 0;
      //std::map<int, float> mapTrackIdToEDep;
      for(auto const& hit : hits){
          std::vector<sim::TrackIDE> showertrackides = bt_serv->HitToTrackIDEs(clockData,hit); 
	  for(auto const& trackIde : showertrackides){
	      mapTrackIdToShowerhitEnergy[(std::abs(trackIde.trackID))] += trackIde.energy;
	      tot_shw_eng += trackIde.energy;
	      //mapTrackIdToEDep[(trackIde.trackID)] += trackIde.energy;
	  }
      }
      
      std::map<int, float> mapMotherIdToEnergyDep;
      //std::map<int, float> mapMotherIdToEDep;
      for(auto const& showermotherIter : mom_id_map){
          for(auto const& ShDaughter : showermotherIter.second){
              mapMotherIdToEnergyDep[showermotherIter.first] += mapTrackIdToShowerhitEnergy[ShDaughter];
          }
      }
      
      float maxEnergy = -99999;
      int ShowerTrackId  = -99999;
      for(auto const& mapMEIter : mapMotherIdToEnergyDep){
          float MEnergy = mapMEIter.second;
          int MTrackId  = mapMEIter.first;
	  if(MEnergy > maxEnergy){
             maxEnergy = MEnergy;
             ShowerTrackId = MTrackId;
          }
      }
      std::pair<const simb::MCParticle*, float> shw_particle;
      const simb::MCParticle* particle(nullptr);
      if(particle_map.find(ShowerTrackId) != particle_map.end()){
         std::cout << "************* ImageMaker::Back_tracking_shws::Retruning a valid pointer a simb::MCParticle **************\n";
	 shw_particle = std::make_pair(particle_map[ShowerTrackId], float(maxEnergy)/tot_shw_eng);
      } 
      else{
          std::cout << "********** ImageMaker::Back_tracking_shws::Empty list **********\n";
	  shw_particle = std::make_pair(particle, 0);
      }
      return shw_particle;
}
//////////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(ImageMaker)
