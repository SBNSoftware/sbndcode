////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////
#ifndef BLIPANA_H
#define BLIPANA_H

// Framework includes
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

// LArSoft includes
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

// SBNDCode includes
#include "sbndcode/BlipReco/Utils/BlipUtils.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

// C++ includes
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

// ROOT includes
#include "TTree.h"
#include "TH1F.h"

// Set global constants and max array sizes
const int kNplanes  = 3;  
const int kMaxHits  = 20000;
const int kMaxTrks  = 1000;
const int kMaxG4    = 10000;

namespace sbnd { 
  
  //###################################################
  //  Data storage structure
  //###################################################
  class BlipAnaTreeDataStruct 
  {
    public:
 
    // --- Main TTree object ---
    TTree* tree;

    // --- Event information ---   
    int   event;                    // event number
    int   run;                      // run number
    float lifetime;                 // electron lifetime (us)

    // --- G4 information ---
    float total_depEnergy;          // total deposited energy in AV
    float total_numElectrons;       // total electrons reaching anode wires
    int   nparticles;               // number of G4 particles
    int   isPrimary[kMaxG4];        // is primary particle
    int   trackID[kMaxG4];          // G4 track ID
    int   pdg[kMaxG4];              // PDG
    int   nDaughters[kMaxG4];       // number of daughters
    int   mother[kMaxG4];           // mother particle
    float E[kMaxG4];                // initial energy (MeV)
    float endE[kMaxG4];             // final energy (MeV)
    float mass[kMaxG4];             // mass (MeV)
    float P[kMaxG4];                // momentum (MeV)
    float Px[kMaxG4];               // momentum x (MeV)
    float Py[kMaxG4];               // momentum y (MeV)
    float Pz[kMaxG4];               // momentum z (MeV)
    float X[kMaxG4];                // starting x (cm)
    float Y[kMaxG4];                // starting y (cm)
    float Z[kMaxG4];                // starting y (cm)
    float endX[kMaxG4];             // ending x (cm)
    float endY[kMaxG4];             // ending y (cm)
    float endZ[kMaxG4];             // ending y (cm)
    float startT[kMaxG4];           // starting time (us)
    float endT[kMaxG4];             // ending time (us)
    float pathlen[kMaxG4];          // path length (cm)
    //int   inTPCActive[kMaxG4];      // in active volume
    //float startpointX_tpcAV[kMaxG4];// first x in AV (cm)
    //float startpointY_tpcAV[kMaxG4];// first y in AV (cm)
    //float startpointZ_tpcAV[kMaxG4];// first z in AV (cm)
    //float endpointX_tpcAV[kMaxG4];  // last x in AV (cm)
    //float endpointY_tpcAV[kMaxG4];  // last y in AV (cm)
    //float endpointZ_tpcAV[kMaxG4];  // last z in AV (cm)
    int   numElectrons[kMaxG4];     // electrons reaching anode wires
    float depEnergy[kMaxG4];        // energy deposited in AV
    std::vector<std::string> process;// process name
  
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
      event                 = -999; // event-wide info
      run                   = -999;
      lifetime              = -999;
      total_depEnergy       = -999;
      total_numElectrons    = -999;
      nparticles            = 0;    // G4 particles
      FillWith(isPrimary,   -9);
      FillWith(trackID,     -999);
      FillWith(pdg,         -99999);
      FillWith(nDaughters,  -999);
      FillWith(mother,      -999);
      FillWith(E,           -999.);
      FillWith(endE,        -999.);
      FillWith(mass,        -999.);
      FillWith(P,           -999.);
      FillWith(Px,          -999.);
      FillWith(Py,          -999.);
      FillWith(Pz,          -999.);
      FillWith(X, -99999.);
      FillWith(Y, -99999.);
      FillWith(Z, -99999.);
      FillWith(endX,   -99999.);
      FillWith(endY,   -99999.);
      FillWith(endZ,   -99999.);
      FillWith(startT,      -99999.);
      FillWith(endT,        -99999.);
      FillWith(pathlen,     -999.);
      //FillWith(inTPCActive, -9);
      //FillWith(startpointX_tpcAV, -99999.);
      //FillWith(startpointY_tpcAV, -99999.);
      //FillWith(startpointZ_tpcAV, -99999.);
      //FillWith(endpointX_tpcAV,   -99999.);
      //FillWith(endpointY_tpcAV,   -99999.);
      //FillWith(endpointZ_tpcAV,   -99999.);
      FillWith(numElectrons,-999);
      FillWith(depEnergy,   -999.);
      FillWith(process,     "");
      nhits                 = 0;    // TPC hits
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

    // === Function for resizing vectors (if necessary) ===
    // To be called after numbers of hits/tracks/particles
    // in the event has been determined
    void Resize() {
      if( nparticles >= 0 ) process.assign(nparticles,"");
    }
        
    // === Function for initializing tree branches ===
    void MakeTree(){
      art::ServiceHandle<art::TFileService> tfs;
      tree = tfs->make<TTree>("anatree","analysis tree");
      tree->Branch("event",&event,"event/I");
      tree->Branch("run",&run,"run/I");
      tree->Branch("lifetime",&lifetime,"lifetime/F");
      tree->Branch("total_depEnergy",&total_depEnergy,"total_depEnergy/F");
      tree->Branch("total_numElectrons",&total_numElectrons,"total_numElectrons/F");
      tree->Branch("nparticles",&nparticles,"nparticles/I");
      tree->Branch("isPrimary",isPrimary,"isPrimary[nparticles]/I");
      tree->Branch("trackID",trackID,"trackID[nparticles]/I");
      tree->Branch("pdg",pdg,"pdg[nparticles]/I");
      tree->Branch("nDaughters",nDaughters,"nDaughters[nparticles]/I");
      tree->Branch("mother",mother,"mother[nparticles]/I");
      tree->Branch("E",E,"E[nparticles]/F");
      tree->Branch("endE",endE,"endE[nparticles]/F");
      tree->Branch("mass",mass,"mass[nparticles]/F");
      tree->Branch("P",P,"P[nparticles]/F");
      tree->Branch("Px",Px,"Px[nparticles]/F");
      tree->Branch("Py",Py,"Py[nparticles]/F");
      tree->Branch("Pz",Pz,"Pz[nparticles]/F");
      tree->Branch("X",X,"X[nparticles]/F");
      tree->Branch("Y",Y,"Y[nparticles]/F");
      tree->Branch("Z",Z,"Z[nparticles]/F");
      tree->Branch("endX",endX,"endX[nparticles]/F");
      tree->Branch("endY",endY,"endY[nparticles]/F");
      tree->Branch("endZ",endZ,"endZ[nparticles]/F");
      tree->Branch("startT",startT,"startT[nparticles]/F");
      tree->Branch("endT",endT,"endT[nparticles]/F");
      tree->Branch("pathlen",pathlen,"pathlen[nparticles]/F");
      //tree->Branch("inTPCActive",inTPCActive,"inTPCActive[nparticles]/I");
      //tree->Branch("startpointX_tpcAV",startpointX_tpcAV,"startpointX_tpcAV[nparticles]/F");
      //tree->Branch("startpointY_tpcAV",startpointY_tpcAV,"startpointY_tpcAV[nparticles]/F");
      //tree->Branch("startpointZ_tpcAV",startpointZ_tpcAV,"startpointZ_tpcAV[nparticles]/F");
      //tree->Branch("endpointX_tpcAV",endpointX_tpcAV,"endpointX_tpcAV[nparticles]/F");
      //tree->Branch("endpointY_tpcAV",endpointY_tpcAV,"endpointY_tpcAV[nparticles]/F");
      //tree->Branch("endpointZ_tpcAV",endpointZ_tpcAV,"endpointZ_tpcAV[nparticles]/F");
      //tree->Branch("numElectrons",numElectrons,"numElectrons[nparticles]/I");
      tree->Branch("depEnergy",depEnergy,"depEnergy[nparticles]/F");
      tree->Branch("process",&process);
      tree->Branch("nhits",&nhits,"nhits/I");
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


  //###################################################
  //  BlipAna class definition
  //###################################################
  class BlipAna : public art::EDAnalyzer 
  { 
    public:
    explicit BlipAna(fhicl::ParameterSet const& pset);
    virtual ~BlipAna();
    void beginJob();
    void analyze(const art::Event& evt);

    private:
    
    // --- Services --- 
    art::ServiceHandle<geo::Geometry> fGeo;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    
    // --- Detector and clock data ---
    float clkdata_BeamGateTime;
    
    // --- Data and calo objects ---
    BlipAnaTreeDataStruct*  fData;
    calo::CalorimetryAlg    fCaloAlg;

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



//###################################################
//  BlipAna constructor and destructor
//###################################################
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


//###################################################
//  beginJob: retrieve relevant detector and clock
//  data here, and save them into class variables
//###################################################
void sbnd::BlipAna::beginJob() {
  auto const clockData  = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  //auto const detProp    = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
  clkdata_BeamGateTime = clockData.BeamGateTime();
}

//###################################################
//  Main event-by-event analysis
//###################################################
void sbnd::BlipAna::analyze(const art::Event& evt)
{
  std::cout
  <<"=========== BlipAna =========================\n"
  <<"Processing event "<<evt.id().event()<<"\n";
  fData->Clear();

  std::cout<<"Beam time "<<clkdata_BeamGateTime<<"\n";

  //=========================================
  // Event information
  //=========================================
  fData->event  = evt.id().event();
  fData->run    = evt.id().run();
  bool isMC     = !evt.isRealData();


  //=========================================
  // Get data products for this event
  //=========================================
  // -- simchannels
  std::vector<const sim::SimChannel*> simChannels;
  if (isMC) evt.getView(fLArG4ModuleLabel, simChannels);
  // -- G4 particles
  const sim::ParticleList& plist = pi_serv->ParticleList();
  // -- hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);
  // -- tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
  // -- hit<->track associations
  art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel);

  // Resize data struct objects
  fData->nhits      = (int)hitlist.size();
  fData->nparticles = (int)plist.size();
  fData->Resize();
  
  std::cout
  <<"Found "
  <<fData->nparticles<<" G4 particles, "
  <<fData->nhits<<" hits from "<<fHitModuleLabel
  <<"\n";


  //====================================
  // Save G4 particle information
  //====================================
  
  // -------------------------------------------------------
  // Find total visible energy deposited in the LAr AV 
  // by looking at the SimChannels for all three planes
  fData->total_depEnergy = 0;
  fData->total_numElectrons = 0;
  for(auto const &chan : simChannels ) {
    for(auto const &tdcide : chan->TDCIDEMap() ) {
      for(const auto& ide : tdcide.second) {
        fData->total_depEnergy += ide.energy/kNplanes;
        fData->total_numElectrons += ide.numElectrons/kNplanes;
      }
    }
  }
  
  // ------------------------------------------------------
  // Loop through the MCParticles
  sim::ParticleList::const_iterator itPart = plist.begin();
  for(size_t i = 0; (i<plist.size())&&(itPart!=plist.end())&&(i<kMaxG4); i++){
    const simb::MCParticle* pPart = (itPart++)->second;

    // Get track ID and status
    int trackID   = pPart->TrackId();
    int isPrimary = (int)(pPart->Process() == "primary");

    // Convert units and make time reference the
    // start of the main drift window (BeamTriggerTime)
    float mass    = /*GeV->MeV*/1e3 * pPart->Mass();
    float E       = /*GeV->MeV*/1e3 * pPart->E();
    float endE    = /*GeV->MeV*/1e3 * pPart->EndE();
    float P       = /*GeV->MeV*/1e3 * pPart->Momentum().Vect().Mag();
    float Px      = /*GeV->MeV*/1e3 * pPart->Px();
    float Py      = /*GeV->MeV*/1e3 * pPart->Py();
    float Pz      = /*GeV->MeV*/1e3 * pPart->Pz();
    float T       = pPart->T() - clkdata_BeamGateTime;
    float endT    = pPart->EndT() - clkdata_BeamGateTime;
 
    // Get the total energy deposited by this particle by looking at IDEs from 3 planes.
    // The value should be the same on all planes, but better safe than sorry...
    geo::View_t views[3]={geo::kU, geo::kV, geo::kW};
    int nPlanes = 0;
    float totalE_particle = 0, totalne_particle = 0;
    for(const geo::View_t view : views ) {
      std::vector<const sim::IDE* > ides = bt_serv->TrackIdToSimIDEs_Ps(trackID, view);
      if( ides.size() ) nPlanes++;
      for (auto ide: ides) {
        totalE_particle += ide->energy;
        totalne_particle += ide->numElectrons;
      }
    }
    if(nPlanes) {
      fData->depEnergy[i]     = totalE_particle/nPlanes;
      fData->numElectrons[i]  = totalne_particle/nPlanes;
    }
      
    fData->trackID[i]         = trackID;
    fData->isPrimary[i]       = isPrimary;
    fData->pdg[i]             = pPart->PdgCode();
    fData->nDaughters[i]      = pPart->NumberDaughters();
    fData->mother[i]          = pPart->Mother();
    fData->E[i]               = E;
    fData->endE[i]            = endE;
    fData->mass[i]            = mass;
    fData->P[i]               = P;
    fData->Px[i]              = Px;
    fData->Py[i]              = Py;
    fData->Pz[i]              = Pz;
    fData->X[i]               = pPart->Vx();
    fData->Y[i]               = pPart->Vy();
    fData->Z[i]               = pPart->Vz();
    fData->endX[i]            = pPart->EndPosition()[0];
    fData->endY[i]            = pPart->EndPosition()[1];
    fData->endZ[i]            = pPart->EndPosition()[2];
    fData->startT[i]          = T;
    fData->endT[i]            = endT;
    fData->pathlen[i]         = BlipUtils::PathLength(*pPart);
    fData->process[i]         = pPart->Process();
    
  }//endloop over G4 particles


  //====================================
  // Save hit information
  //====================================
  for(size_t i=0; (i<hitlist.size()) && (i<kMaxHits); i++){
    fData->hit_tpc[i]     = hitlist[i]->WireID().TPC;
    fData->hit_channel[i] = hitlist[i]->Channel();
    fData->hit_plane[i]   = hitlist[i]->WireID().Plane;
    fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
    fData->hit_peakT[i]   = hitlist[i]->PeakTime();
    fData->hit_rms[i]     = hitlist[i]->RMS();
    fData->hit_ph[i]	    = hitlist[i]->PeakAmplitude();
    fData->hit_area[i]    = hitlist[i]->Integral();
    fData->hit_charge[i]  = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),hitlist[i]->WireID().Plane);
    // Find associated track
    if (fmtk.isValid()){
      if (fmtk.at(i).size())  fData->hit_trkid[i] = fmtk.at(i)[0]->ID();
      else                    fData->hit_trkid[i] = -1;
    }
    // Find G4 particle ID for leading contributor
    if( BlipUtils::DoesHitHaveSimChannel(simChannels,hitlist[i]) )
      BlipUtils::HitTruth( hitlist[i], fData->hit_mcid[i], fData->hit_mcfrac[i], fData->hit_mcenergy[i], fData->hit_mccharge[i]);
  }//endloop hits

  // Procedure
  //  - Look for hits that were not included in a track
  //  - Plane-to-plane time matching
  //  - Wire intersection check to get XYZ
  //  - Create "blip" object and save to tree (nblips, blip_xyz, blip_charge, blip_mcenergy)

  //====================================
  // Fill TTree
  //====================================
  fData->FillTree();
}





DEFINE_ART_MODULE(sbnd::BlipAna)

#endif
