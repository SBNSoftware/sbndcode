////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////
#ifndef BLIPANA_H
#define BLIPANA_H

//##########################
//### Framework includes ###
//##########################
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//########################
//### LArSoft includes ###
//########################
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

//########################
//### SBNDCode includes###
//########################
#include "sbndcode/BlipReco/Utils/BlipUtils.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

//########################
//### C++ includes     ###
//########################
#include <cstring>
#include <vector>
#include <map>
#include <utility>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cmath>

//########################
//### ROOT includes    ###
//########################
#include "TTree.h"
#include "TH1F.h"

// Set global constants and max array sizes
const int kNplanes  = 3;  
const int kMaxHits  = 25000;

namespace sbnd { 
  
  //====================================================
  //  Data storage structure
  //====================================================
  class BlipAnaTreeDataStruct 
  {
    public:
 
    // --- Main TTree object ---
    TTree* tree;

    // --- Event information ---   
    int   event;                    // event number
    int   run;                      // run number

    // --- Hit information ---
    int	  nhits;                    // number of hits
    int	  hit_tpc[kMaxHits];        // tpc number
    int	  hit_plane[kMaxHits];      // plane number
    int	  hit_wire[kMaxHits];       // wire number
    int	  hit_channel[kMaxHits];    // channel ID
    float	hit_peakT[kMaxHits];      // peak time (tick)
    float hit_rms[kMaxHits];        // shape RMS
    float	hit_ph[kMaxHits];         // amplitude
    float	hit_area[kMaxHits];       // charge (area) in ADC units
    float hit_charge[kMaxHits];     // reconstructed number of electrons
    int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
    int	  hit_mcid[kMaxHits];       // G4 TrackID of leading MCParticle that created the hit
    float hit_mcfrac[kMaxHits];     // fraction of hit energy from leading MCParticle
    float hit_mcenergy[kMaxHits];   // true energy
    float hit_mccharge[kMaxHits];   // true number of electrons (drift-attenuated)

    // === Function for resetting data ===
    void Clear(){ 
      event     = -9;
      run       = -9;
      nhits     = 0;
      FillWith(hit_tpc,     -999);
      FillWith(hit_plane,   -999);
      FillWith(hit_wire,    -999);
      FillWith(hit_channel, -999);
      FillWith(hit_peakT,   -999);
      FillWith(hit_rms,     -999);
      FillWith(hit_ph,      -999);
      FillWith(hit_area,    -999);
      FillWith(hit_trkid,   -999);
      FillWith(hit_mcid,    -999);
      FillWith(hit_mcfrac,  -999);
      FillWith(hit_mcenergy,-999);
      FillWith(hit_mccharge,-999);
      FillWith(hit_charge,  -999);
    }
        
    // === Function for initializing tree branches ===
    void MakeTree(){
      art::ServiceHandle<art::TFileService> tfs;
      tree = tfs->make<TTree>("anatree","analysis tree");
      tree->Branch("event",&event,"event/I");
      tree->Branch("run",&run,"run/I");
      tree->Branch("nhits",nhits,"nhits/I");
      tree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I"); 
      tree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I"); 
      tree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I"); 
      tree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I"); 
      tree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/F"); 
      tree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/F"); 
      tree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/F"); 
      tree->Branch("hit_area",hit_area,"hit_area[nhits]/F"); 
      tree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I"); 
      tree->Branch("hit_mcid",hit_mcid,"hit_mcid[nhits]/I");
      tree->Branch("hit_mcfrac",hit_mcfrac,"hit_mcfrac[nhits]/F"); 
      tree->Branch("hit_mcenergy",hit_mcenergy,"hit_mcenergy[nhits]/F"); 
      tree->Branch("hit_mccharge",hit_mccharge,"hit_mccharge[nhits]/F"); 
      tree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
    }
    
    // === Function for filling tree ===
    void FillTree(){ tree->Fill(); }
    

  };//BlipAnaTreeDataStruct class


  //====================================================
  //  BlipAna class definition
  //====================================================
  class BlipAna : public art::EDAnalyzer 
  { 
    public:
    // --- constructor/destructor ---
    explicit BlipAna(fhicl::ParameterSet const& pset);
    virtual ~BlipAna();
    // --- art functions ---
    void analyze(const art::Event& evt);

    private:

    // --- Data and calo alg objects ---
    BlipAnaTreeDataStruct*  fData;
    calo::CalorimetryAlg    fCaloAlg;
    
    // --- ParticleInventoryService ---
    //art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // --- BackTrackerService ---
    //art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // --- FCL configs ---
    std::string   fHitModuleLabel;
    std::string   fLArG4ModuleLabel;
    std::string   fTrackModuleLabel;

    // --- Histograms ---
    TH1F*         h_nhits;
  
    // === Function for initializing histograms ===
    void InitializeHistograms(){
      art::ServiceHandle<art::TFileService> tfs;
      h_nhits = tfs->make<TH1F>("nhits","Number of hits per event",1000,0,1000);
    }

  };//class BlipAna

}//namespace sbnd



//====================================================
//  BlipAna constructor and destructor
//====================================================
sbnd::BlipAna::BlipAna(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData(nullptr)
  ,fCaloAlg(pset.get< fhicl::ParameterSet >("CaloAlg"))
  ,fHitModuleLabel(pset.get<std::string>("HitModuleLabel","gaushit"))
  ,fLArG4ModuleLabel(pset.get<std::string>("LArG4ModuleLabel","largeant"))
  ,fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel","pandora"))
{
  fData = new BlipAnaTreeDataStruct();
  fData ->Clear();
  fData ->MakeTree();
  InitializeHistograms();
}
sbnd::BlipAna::~BlipAna(){}



//====================================================
//  Main event-by-event analysis
//====================================================
void sbnd::BlipAna::analyze(const art::Event& evt)
{
  std::cout
  <<"=========== BlipAna =========================\n"
  <<"Processing event "<<evt.id().event()<<"\n";
  fData->Clear();

  // Event information
  fData->event  = evt.id().event();
  fData->run    = evt.id().run();
  bool isMC     = !evt.isRealData();
 
  //----------------------------------
  // Get data products for this event
  //---------------------------------- 
  // * SimChannels
  std::vector<const sim::SimChannel*> simChannels;
  if (isMC) evt.getView(fLArG4ModuleLabel, simChannels);
  // * Hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);
  // * Tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
  // * Hit <--> Track associations
  art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel);

  //-------------------------------
  // Save hit information
  //------------------------------ 
  fData->nhits = (int)hitlist.size();
  for(size_t i=0; i<hitlist.size(); i++){
    fData->hit_tpc[i]     = hitlist[i]->WireID().TPC;
    fData->hit_channel[i] = hitlist[i]->Channel();
    fData->hit_plane[i]   = hitlist[i]->WireID().Plane;
    fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
    fData->hit_peakT[i]   = hitlist[i]->PeakTime();
    fData->hit_rms[i]     = hitlist[i]->RMS();
    fData->hit_ph[i]	    = hitlist[i]->PeakAmplitude();
    fData->hit_area[i]    = hitlist[i]->Integral();
    fData->hit_charge[i]  = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),hitlist[i]->WireID().Plane);
    // mcid, mcfrac, mcenergy, mccharge

    // Associated track
    if (fmtk.isValid()){
      if (fmtk.at(i).size())  fData->hit_trkid[i] = fmtk.at(i)[0]->ID();
      else                    fData->hit_trkid[i] = -1;
    }
    
    // Find G4 particle ID for leading contributor to this hit
    if( BlipUtils::DoesHitHaveSimChannel(simChannels,hitlist[i]) ) {
      BlipUtils::HitTruth( hitlist[i], fData->hit_mcid[i], fData->hit_mcfrac[i], fData->hit_mcenergy[i], fData->hit_charge[i]);
    }

  }//endloop hits

  // Procedure
  //  - Look for hits that were not included in a track
  //  - Plane-to-plane time matching
  //  - Wire intersection check to get XYZ
  //  - Create "blip" object and save to tree (nblips, blip_xyz, blip_charge, blip_mcenergy)

  // Fill TTree
  fData->FillTree();
}





DEFINE_ART_MODULE(sbnd::BlipAna)

#endif
