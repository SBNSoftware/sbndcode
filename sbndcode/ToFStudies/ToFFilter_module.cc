////////////////////////////////////////////////////////////////////////
// Class:       ToFFilter
// Plugin Type: filter (art v3_05_01)
// File:        ToFFilter_module.cc
//
// Generated at Fri Jan 28 17:07:23 2022 by Varuna Crishan Meddage using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
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
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"

#include "lardataobj/RecoBase/OpHit.h"

#include "sbnobj/SBND/ToF/ToF.hh"

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
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include <map>
#include <memory>

using namespace std;

namespace sbnd{


class ToFFilter : public art::EDFilter {
public:
  explicit ToFFilter(fhicl::ParameterSet const& p);
  
  ToFFilter(ToFFilter const&) = delete;
  ToFFilter(ToFFilter&&) = delete;
  ToFFilter& operator=(ToFFilter const&) = delete;
  ToFFilter& operator=(ToFFilter&&) = delete;

  bool filter(art::Event& evt) override;

private:
  art::InputTag ftofLhitLabel;
  art::InputTag ftofChitLabel;
  art::InputTag ftofLflashLabel;
  art::InputTag ftofCflashLabel;
  art::InputTag ftofLflashhitLabel;
  art::InputTag ftofCflashhitLabel;
  bool fuse_Lhit;
  bool fuse_Chit;
  bool fuse_Lflsh;
  bool fuse_Cflsh;
  bool fuse_Lflsh_hit;
  bool fuse_Cflsh_hit;
  float ftof_Lhit_cut;
  float ftof_Chit_cut;
  float ftof_Lflsh_cut;
  float ftof_Cflsh_cut;
  float ftof_Lflshhit_cut;
  float ftof_Cflshhit_cut;
};


ToFFilter::ToFFilter(fhicl::ParameterSet const& pset): 
EDFilter{pset},
ftofLhitLabel(pset.get<art::InputTag>("tofLhitLabel")),
ftofChitLabel(pset.get<art::InputTag>("tofChitLabel")),
ftofLflashLabel(pset.get<art::InputTag>("tofLflashLabel")),
ftofCflashLabel(pset.get<art::InputTag>("tofCflashLabel")),
ftofLflashhitLabel(pset.get<art::InputTag>("tofLflashhitLabel")),
ftofCflashhitLabel(pset.get<art::InputTag>("tofCflashhitLabel")),
fuse_Lhit(pset.get<bool>("use_Lhit")),
fuse_Chit(pset.get<bool>("use_Chit")),
fuse_Lflsh(pset.get<bool>("use_Lflsh")),
fuse_Cflsh(pset.get<bool>("use_Cflsh")),
fuse_Lflsh_hit(pset.get<bool>("use_Lflsh_hit")),
fuse_Cflsh_hit(pset.get<bool>("use_Cflsh_hit")),		
ftof_Lhit_cut(pset.get<float>("tof_Lhit_cut")),
ftof_Chit_cut(pset.get<float>("tof_Chit_cut")),
ftof_Lflsh_cut(pset.get<float>("tof_Lflsh_cut")),
ftof_Cflsh_cut(pset.get<float>("tof_Cflsh_cut")),
ftof_Lflshhit_cut(pset.get<float>("tof_Lflshhit_cut")),
ftof_Cflshhit_cut(pset.get<float>("tof_Cflshhit_cut"))
{
}

bool ToFFilter::filter(art::Event& evt)
{
  bool keep_event = true;
  
  // ================================== Calculatin ToF values using Largest optical hit method =========================
  	
  if(fuse_Lhit){
     art::Handle< std::vector<sbnd::ToF::ToF> > tofLhitListHandle;
     std::vector< art::Ptr<sbnd::ToF::ToF> >    tofLhitList;
     if( evt.getByLabel(ftofLhitLabel,tofLhitListHandle))
         art::fill_ptr_vector(tofLhitList, tofLhitListHandle);
     
     std::vector<bool> is_cos_like;
     
     for(auto const& tf_lhit : tofLhitList){
         if(tf_lhit->tof >= ftof_Lhit_cut) is_cos_like.push_back(false);
	 else is_cos_like.push_back(true);
     }
     
     if(is_cos_like.size()){
	bool found_nu=false;
        for(auto itr: is_cos_like){
	    if(!itr){ 
	       found_nu=true; 
	       break;
	    } 
	}
	
	if(!found_nu) keep_event=false;
     }
  } // Using largest optical hit to calculate tof
  
  //======================================================================================================================
  
  // ================================== Calculatin ToF values using Closest optical hit method ===========================
  
  if(fuse_Chit){
     art::Handle< std::vector<sbnd::ToF::ToF> > tofChitListHandle;
     std::vector< art::Ptr<sbnd::ToF::ToF> >    tofChitList;
     if( evt.getByLabel(ftofChitLabel,tofChitListHandle))
         art::fill_ptr_vector(tofChitList, tofChitListHandle);
     
     std::vector<bool> is_cos_like;
     
     for(auto const& tf_chit : tofChitList){
         if(tf_chit->tof >= ftof_Chit_cut) is_cos_like.push_back(false);
	 else is_cos_like.push_back(true);
     }
     
     if(is_cos_like.size()){
	bool found_nu=false;
        for(auto itr: is_cos_like){
	    if(!itr){ 
	       found_nu=true; 
	       break;
	    } 
	}
	
	if(!found_nu) keep_event=false;
     }
  } // Using closest optical hit to calculate tof
  
  //======================================================================================================================
  
  //==================================== Calculation ToF values using Largest optical flash ==============================
  
  if(fuse_Lflsh){
     art::Handle< std::vector<sbnd::ToF::ToF> > tofLflashListHandle;
     std::vector< art::Ptr<sbnd::ToF::ToF> >    tofLflashList;
     if( evt.getByLabel(ftofLflashLabel,tofLflashListHandle))
         art::fill_ptr_vector(tofLflashList, tofLflashListHandle);
     
     std::vector<bool> is_cos_like;
     
     for(auto const& tf_lflash : tofLflashList){
         if(tf_lflash->tof >= ftof_Lflsh_cut) is_cos_like.push_back(false);
	 else is_cos_like.push_back(true);
     }
     
     if(is_cos_like.size()){
	bool found_nu=false;
        for(auto itr: is_cos_like){
	    if(!itr){ 
	       found_nu=true; 
	       break;
	    } 
	}
	
	if(!found_nu) keep_event=false;
     }
  } // Using largest optical flash to calculate tof
  
  //======================================================================================================================
  
  //==================================== Calculation ToF values using Largest optical flash ==============================
  
  if(fuse_Cflsh){
     art::Handle< std::vector<sbnd::ToF::ToF> > tofCflashListHandle;
     std::vector< art::Ptr<sbnd::ToF::ToF> >    tofCflashList;
     if( evt.getByLabel(ftofCflashLabel,tofCflashListHandle))
         art::fill_ptr_vector(tofCflashList, tofCflashListHandle);
     
     std::vector<bool> is_cos_like;
     
     for(auto const& tf_cflash : tofCflashList){
         if(tf_cflash->tof >= ftof_Cflsh_cut) is_cos_like.push_back(false);
	 else is_cos_like.push_back(true);
     }
     
     if(is_cos_like.size()){
	bool found_nu=false;
        for(auto itr: is_cos_like){
	    if(!itr){ 
	       found_nu=true; 
	       break;
	    } 
	}
	
	if(!found_nu) keep_event=false;
     }
  } // Using closest optical flash to calculate tof
  
  //======================================================================================================================
  
  //=========================Calculation ToF values using Earliest hit of the Largest flash ==============================
  
  if(fuse_Lflsh_hit){
     art::Handle< std::vector<sbnd::ToF::ToF> > tofLflashhitListHandle;
     std::vector< art::Ptr<sbnd::ToF::ToF> >    tofLflashhitList;
     if( evt.getByLabel(ftofLflashhitLabel,tofLflashhitListHandle))
         art::fill_ptr_vector(tofLflashhitList, tofLflashhitListHandle);
     
     std::vector<bool> is_cos_like;
     
     for(auto const& tf_lflashhit : tofLflashhitList){
         if(tf_lflashhit->tof >= ftof_Lflshhit_cut) is_cos_like.push_back(false);
	 else is_cos_like.push_back(true);
     }
     
     if(is_cos_like.size()){
	bool found_nu=false;
        for(auto itr: is_cos_like){
	    if(!itr){ 
	       found_nu=true; 
	       break;
	    } 
	}
	
	if(!found_nu) keep_event=false;
     }
  } // Using earliest hit of the largest optical flash to calculate tof
  
  //======================================================================================================================
  
  //=========================Calculation ToF values using Earliest hit of the Closest flash ==============================
  
  if(fuse_Cflsh_hit){
     art::Handle< std::vector<sbnd::ToF::ToF> > tofCflashhitListHandle;
     std::vector< art::Ptr<sbnd::ToF::ToF> >    tofCflashhitList;
     if( evt.getByLabel(ftofCflashhitLabel,tofCflashhitListHandle))
         art::fill_ptr_vector(tofCflashhitList, tofCflashhitListHandle);
     
     std::vector<bool> is_cos_like;
     
     for(auto const& tf_cflashhit : tofCflashhitList){
         if(tf_cflashhit->tof >= ftof_Cflshhit_cut) is_cos_like.push_back(false);
	 else is_cos_like.push_back(true);
     }
     
     if(is_cos_like.size()){
	bool found_nu=false;
        for(auto itr: is_cos_like){
	    if(!itr){ 
	       found_nu=true; 
	       break;
	    } 
	}
	
	if(!found_nu) keep_event=false;
     }
  } // Using earliest hit of the closest optical flash to calculate tof
  
  //======================================================================================================================
  
  return keep_event;
}

DEFINE_ART_MODULE(ToFFilter)
}
