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
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/BlipReco/Utils/BlipUtils.h"

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

// Helper templates for initializing arrays
namespace{  
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, ITER to, TYPE value) 
    { std::fill(from, to, value); }
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }
  template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }
}

// Set global constants and max array sizes
const int kMaxHits  = 20000;
const int kMaxTrks  = 1000;
const int kMaxG4    = 10000;
const int kMaxEDeps = 10000;
const int kNplanes  = 3;  
    
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
    float startPointx[kMaxG4];      // starting x (cm)
    float startPointy[kMaxG4];      // starting y (cm)
    float startPointz[kMaxG4];      // starting y (cm)
    float endPointx[kMaxG4];        // ending x (cm)
    float endPointy[kMaxG4];        // ending y (cm)
    float endPointz[kMaxG4];        // ending y (cm)
    float startT[kMaxG4];           // starting time (us)
    float endT[kMaxG4];             // ending time (us)
    float pathlen[kMaxG4];          // path length (cm)
    int   numElectrons[kMaxG4];     // electrons reaching anode wires
    float depEnergy[kMaxG4];        // energy deposited in AV (MeV)
    std::vector<std::string> process;// process name
 
    // --- True energy deposit info (derived) ---
    int   nedeps;                   // number of true localized energy depositions
    int   edep_g4id[kMaxEDeps];     // leading G4 track ID
    float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
    float edep_x[kMaxEDeps];        // x (cm)
    float edep_y[kMaxEDeps];        // y (cm)
    float edep_z[kMaxEDeps];        // z (cm)
    float edep_ds[kMaxEDeps];       // extent (cm)

    // --- Hit information ---
    int	  nhits;                    // number of hits
    int	  hit_tpc[kMaxHits];        // tpc number
    int	  hit_plane[kMaxHits];      // plane number
    int	  hit_wire[kMaxHits];       // wire number
    int	  hit_channel[kMaxHits];    // channel ID
    float	hit_peakT[kMaxHits];      // raw peak time (tick)
    float	hit_time[kMaxHits];       // corrected peak time (tick)
    float hit_rms[kMaxHits];        // shape RMS
    float	hit_ph[kMaxHits];         // amplitude
    float	hit_area[kMaxHits];       // charge (area) in ADC units
    float hit_charge[kMaxHits];     // reconstructed number of electrons
    int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
    int   hit_isreal[kMaxHits];     // is this hit real?
    int	  hit_g4id[kMaxHits];       // G4 TrackID of leading particle
    float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
    float hit_g4energy[kMaxHits];   // true energy
    float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
    int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
    int   hit_blipid[kMaxHits];     // key of Blip in which hit was included


    // === Function for resetting data ===
    void Clear(){ 
      event                 = -999; // --- event-wide info ---
      run                   = -999;
      lifetime              = -999;
      total_depEnergy       = -999;
      total_numElectrons    = -999;
      nparticles            = 0;    // --- G4 particles ---
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
      FillWith(startPointx, -99999.);
      FillWith(startPointy, -99999.);
      FillWith(startPointz, -99999.);
      FillWith(endPointx,   -99999.);
      FillWith(endPointy,   -99999.);
      FillWith(endPointz,   -99999.);
      FillWith(startT,      -99999.);
      FillWith(endT,        -99999.);
      FillWith(pathlen,     -999.);
      FillWith(numElectrons,-999);
      FillWith(depEnergy,   -999.);
      FillWith(process,     "");
      nedeps                = 0;    // --- EDeps ---
      FillWith(edep_energy, -999);
      FillWith(edep_g4id,   -9);
      FillWith(edep_x,      -99999.);
      FillWith(edep_y,      -99999.);
      FillWith(edep_z,      -99999.);
      FillWith(edep_ds,     -999);
      nhits                 = 0;    // --- TPC hits ---
      FillWith(hit_tpc,     -999);
      FillWith(hit_plane,   -999);
      FillWith(hit_wire,    -999);
      FillWith(hit_channel, -999);
      FillWith(hit_peakT,   -999);
      FillWith(hit_time,    -999);
      FillWith(hit_rms,     -999);
      FillWith(hit_ph,      -999);
      FillWith(hit_area,    -999);
      FillWith(hit_charge,  -999);
      FillWith(hit_isreal,  -999);
      FillWith(hit_trkid,     -9);
      FillWith(hit_g4id,      -9);
      FillWith(hit_g4frac,  -999);
      FillWith(hit_g4energy,-999);
      FillWith(hit_g4charge,-999);
      FillWith(hit_clustid,   -9);
      FillWith(hit_blipid,     -9);
    }

    // === Function for resizing vectors (if necessary) ===
    // To be called after numbers of hits/tracks/particles
    // in the event has been determined
    void Resize() {
      if(nparticles) process.assign(nparticles,"");
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
      tree->Branch("startPointx",startPointx,"startPointx[nparticles]/F");
      tree->Branch("startPointy",startPointy,"startPointy[nparticles]/F");
      tree->Branch("startPointz",startPointz,"startPointz[nparticles]/F");
      tree->Branch("endPointx",endPointx,"endPointx[nparticles]/F");
      tree->Branch("endPointy",endPointy,"endPointy[nparticles]/F");
      tree->Branch("endPointz",endPointz,"endPointz[nparticles]/F");
      tree->Branch("startT",startT,"startT[nparticles]/F");
      tree->Branch("endT",endT,"endT[nparticles]/F");
      tree->Branch("pathlen",pathlen,"pathlen[nparticles]/F");
      tree->Branch("numElectrons",numElectrons,"numElectrons[nparticles]/I");
      tree->Branch("depEnergy",depEnergy,"depEnergy[nparticles]/F");
      tree->Branch("process",&process);
      tree->Branch("nedeps",&nedeps,"nedeps/I");
      tree->Branch("edep_g4id",edep_g4id,"edep_g4id[nedeps]/I"); 
      tree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
      tree->Branch("edep_x",edep_x,"edep_x[nedeps]/F"); 
      tree->Branch("edep_y",edep_y,"edep_y[nedeps]/F"); 
      tree->Branch("edep_z",edep_z,"edep_z[nedeps]/F"); 
      tree->Branch("edep_ds",edep_ds,"edep_ds[nedeps]/F"); 
      tree->Branch("nhits",&nhits,"nhits/I");
      tree->Branch("hit_tpc",hit_tpc,"hit_tpc[nhits]/I"); 
      tree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I"); 
      tree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I"); 
      tree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I"); 
      tree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/F"); 
      tree->Branch("hit_time",hit_time,"hit_time[nhits]/F"); 
      tree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/F"); 
      tree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/F"); 
      tree->Branch("hit_area",hit_area,"hit_area[nhits]/F"); 
      tree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
      tree->Branch("hit_isreal",hit_isreal,"hit_isreal[nhits]/I"); 
      tree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I"); 
      tree->Branch("hit_g4id",hit_g4id,"hit_g4id[nhits]/I");
      tree->Branch("hit_g4frac",hit_g4frac,"hit_g4frac[nhits]/F"); 
      tree->Branch("hit_g4energy",hit_g4energy,"hit_g4energy[nhits]/F"); 
      tree->Branch("hit_g4charge",hit_g4charge,"hit_g4charge[nhits]/F"); 
      tree->Branch("hit_clustid",hit_clustid,"hit_clustid[nhits]/I"); 
      tree->Branch("hit_blipid",hit_blipid,"hit_blipid[nhits]/I"); 
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
    void PrintParticleInfo(size_t);
  
    // --- Detector and clock data ---
    float detprop_XTicksOffset[kNplanes];

    // --- Data and calo objects ---
    BlipAnaTreeDataStruct*  fData;
    calo::CalorimetryAlg    fCaloAlg;

    // --- FCL configs ---
    std::string         fHitModuleLabel;
    std::string         fLArG4ModuleLabel;
    std::string         fTrackModuleLabel;
    bool                fSaveParticleList;
    float               fBlipMergeDist;
    std::vector<float>  fMinHitRMS;
    std::vector<float>  fMaxHitRMS;
//    std::vector<float>  fMinHitPhWidthRatio;

    // --- Histograms ---
    TH1F*   h_trueblip_ds;
    TH1F*   h_trueblip_energy;
    TH1F*   h_nhits[2][kNplanes];
    TH1F*   h_hitrms[2][kNplanes];
    TH1F*   h_hitph[2][kNplanes];

    // Initialize histograms
    void InitializeHistograms(){
      art::ServiceHandle<art::TFileService> tfs;
      //h_nhits = tfs->make<TH1F>("nhits","Number of hits per event",1000,0,1000);
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
  ,fSaveParticleList(pset.get<bool>("SaveParticleList",true))
  ,fBlipMergeDist(pset.get<float>("BlipMergeDist",0.2))
  ,fMinHitRMS(pset.get<std::vector<float>>("MinHitRMS",{1,1.5,1}))
  ,fMaxHitRMS(pset.get<std::vector<float>>("MaxHitRMS",{5,5,5}))
//  ,fMinHitPhWidthRatio(pset.get<std::vector<float>>("MinHitPhWidthRatio",{0., 5., 0.}
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
  //auto const clockData  = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  //clkdata_BeamGateTime = clockData.BeamGateTime();
  auto const detProp    = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  detprop_XTicksOffset[0] = detProp.GetXTicksOffset(0,0,0);
  detprop_XTicksOffset[1] = detProp.GetXTicksOffset(1,0,0);
  detprop_XTicksOffset[2] = detProp.GetXTicksOffset(2,0,0);
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
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
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
  std::cout
  <<"Total energy deposited: "<<fData->total_depEnergy<<" MeV \n";

  // Save all the "true" blips in the event
  std::vector<BlipUtils::TrueBlip> trueBlipsVec;

  // Loop through the MCParticles
  sim::ParticleList::const_iterator itPart = plist.begin();
  for(size_t i = 0; (i<plist.size())&&(itPart!=plist.end()); i++){
    const simb::MCParticle* pPart = (itPart++)->second;

    // Get important info about particle and convert units if necessary
    int trackID   = pPart->TrackId();
    int isPrimary = (int)(pPart->Process() == "primary");
    float mass    = /*GeV->MeV*/1e3 * pPart->Mass();
    float E       = /*GeV->MeV*/1e3 * pPart->E();
    float endE    = /*GeV->MeV*/1e3 * pPart->EndE();
    float P       = /*GeV->MeV*/1e3 * pPart->Momentum().Vect().Mag();
    float Px      = /*GeV->MeV*/1e3 * pPart->Px();
    float Py      = /*GeV->MeV*/1e3 * pPart->Py();
    float Pz      = /*GeV->MeV*/1e3 * pPart->Pz();
    float ne      = 0;
    float edep    = BlipUtils::PartEnergyDep(trackID,ne);
    float pathlen = BlipUtils::PathLength(*pPart);
   
    // Make true blips 
    BlipUtils::TrueBlip tb = BlipUtils::MakeTrueBlip(trackID);
    if( tb.isValid ) trueBlipsVec.push_back(tb);

    // Save to TTree object
    if(i<kMaxG4){
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
      fData->startPointx[i]     = pPart->Vx();
      fData->startPointy[i]     = pPart->Vy();
      fData->startPointz[i]     = pPart->Vz();
      fData->endPointx[i]       = pPart->EndPosition()[0];
      fData->endPointy[i]       = pPart->EndPosition()[1];
      fData->endPointz[i]       = pPart->EndPosition()[2];
      fData->startT[i]          = pPart->T();
      fData->endT[i]            = pPart->EndT();
      fData->pathlen[i]         = pathlen;
      fData->process[i]         = pPart->Process();
      fData->depEnergy[i]       = edep;
      fData->numElectrons[i]    = ne;
      PrintParticleInfo(i);
    }//endif saving particle list to TTree
  
  }//endloop over G4 particles



  //====================================
  // Merge and save true blip information
  //====================================
  MergeBlips(trueBlipsVec, fBlipMergeDist); 
  fData->nedeps = (int)trueBlipsVec.size();
  std::cout<<"Found "<<trueBlipsVec.size()<<" true blips:\n";
  for(size_t i=0; i<trueBlipsVec.size(); i++ ) {
    fData->edep_energy[i] = trueBlipsVec.at(i).Energy;
    fData->edep_ds[i]     = trueBlipsVec.at(i).Length;
    fData->edep_x[i]      = trueBlipsVec.at(i).Location.X();
    fData->edep_y[i]      = trueBlipsVec.at(i).Location.Y();
    fData->edep_z[i]      = trueBlipsVec.at(i).Location.Z();
    fData->edep_g4id[i]   = trueBlipsVec.at(i).LeadingG4TrackID;
    std::cout
    <<"   ~ "<<trueBlipsVec.at(i).Energy<<" MeV, "
    <<" ds= "<<trueBlipsVec.at(i).Length<<" cm, "
    <<" trkID= "<<trueBlipsVec.at(i).LeadingG4TrackID<<"\n";
  }



  //====================================
  // Save hit information
  //====================================
  
  std::vector<BlipUtils::HitInfo> hitinfo(hitlist.size());
  std::map<int,std::vector<int>> wirehitsMap;

  for(size_t i=0; i<hitlist.size(); i++){
    
    // Variables to be calculated
    int wire = hitlist[i]->WireID().Wire;
    int plane = hitlist[i]->WireID().Plane;
    int trkid = -1;
    int g4id = -1;
    std::set<int> g4ids;
    float g4frac, g4energy, g4charge;
    float charge = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),hitlist[i]->WireID().Plane);
    float time = hitlist[i]->PeakTime()-detprop_XTicksOffset[fData->hit_plane[i]];

    // Find associated track
    if (fmtk.isValid()){
      if (fmtk.at(i).size())  trkid = fmtk.at(i)[0]->ID();
    }
    
    // Find G4 particle ID for leading contributor
    if( BlipUtils::DoesHitHaveSimChannel(simChannels,hitlist[i]) ){
      BlipUtils::HitTruth( hitlist[i], g4id, g4frac, g4energy, g4charge);
      g4ids = BlipUtils::HitTruthIds(hitlist[i]);
    }

    // Add to the wire-by-wire hit map
    wirehitsMap[wire].push_back(i);

    // Bundle the info to access later
    hitinfo[i].hit      = hitlist[i];
    hitinfo[i].g4ids    = g4ids;
    hitinfo[i].trkid    = trkid;
    hitinfo[i].g4id     = g4id;
    hitinfo[i].g4frac   = g4frac;
    hitinfo[i].g4energy = g4energy;
    hitinfo[i].g4charge = g4charge;
    hitinfo[i].Charge   = charge;
    hitinfo[i].Time     = time;
    // Determine realness; group in overlapping hits
    bool isreal         = (g4id > 0);
    hitinfo[i].isreal   = isreal;
    if( isreal ) {
      for(size_t j=0; j<hitlist.size(); j++){
        if( BlipUtils::DoHitsOverlap(hitlist[i],hitlist[j]) )
          hitinfo[j].isreal = isreal;
      }
    }
   
    // Save TTree data
    if( i < kMaxHits ) {
      fData->hit_trkid[i]   = trkid;
      fData->hit_plane[i]   = plane;
      fData->hit_wire[i]    = wire;
      fData->hit_tpc[i]     = hitlist[i]->WireID().TPC;
      fData->hit_channel[i] = hitlist[i]->Channel();
      fData->hit_peakT[i]   = hitlist[i]->PeakTime();
      fData->hit_rms[i]     = hitlist[i]->RMS();
      fData->hit_ph[i]	    = hitlist[i]->PeakAmplitude();
      fData->hit_area[i]    = hitlist[i]->Integral();
      fData->hit_time[i]    = time;
      fData->hit_charge[i]  = charge;
      fData->hit_g4id[i]    = g4id;
      fData->hit_g4frac[i]  = g4frac;
      fData->hit_g4energy[i]= g4energy;
      fData->hit_g4charge[i]= g4charge;
    }
  
  }//endloop hits
  

  //=================================================================
  // Blip Reconstruction
  //================================================================
  //  
  //  Will eventually move these into separate alg class.
  //
  //  Procedure
  //  [x] Look for hits that were not included in a track 
  //  [x] Filter hits based on hit width, etc
  //  - Merge together closely-spaced hits on same wires, save average peakT +/- spread
  //  - Merge together clusters on adjacent wires (if they match up in time)
  //  - Plane-to-plane time matching
  //  - Wire intersection check to get XYZ
  //  - Create "blip" object and save to tree (nblips, blip_xyz, blip_charge, blip_g4energy)

  // Create a hit mask
  std::vector<bool> hitmask(hitlist.size(), false);
  for(size_t i=0; i<hitlist.size(); i++){
    int plane = hitlist[i]->WireID().Plane;
    if( hitlist[i]->RMS() > fMinHitRMS[plane] &&
        hitlist[i]->RMS() < fMaxHitRMS[plane] ) {
      hitmask[i] = true;
      // Any overlapping hits are probably good too
      for(size_t j=0; j<hitlist.size(); j++){
        if( BlipUtils::DoHitsOverlap(hitlist[i],hitlist[j]) ) {
          hitmask[j] = true;
        }
      }
    }
  }

  // Exclude any hits that were part of tracks
  for(size_t i=0; i<hitlist.size(); i++){
    if( hitinfo[i].trkid >= 0 ) hitmask[i] = false;
  }

  std::cout
  <<"We masked "<<std::count(hitmask.begin(),hitmask.end(),true)<<" hits based on width and track assn cuts\n";

  // Create collection of hit clusters on same wires
  std::vector<BlipUtils::HitClust> hitclust;
  std::vector<bool> hitIsClustered(hitlist.size(),false);
  for(auto const& wirehits : wirehitsMap){
    std::vector<int> hits = wirehits.second;
    for(auto const& hi : hits ){
      if( !hitmask[hi] ) continue;
      if( hitIsClustered[hi] ) continue;
      BlipUtils::HitClust hc = BlipUtils::MakeHitClust(hitlist[hi],hitinfo[hi]);
      hitIsClustered[hi] = true;
      // see if we can add other hits to it
      for(auto const& hj : hits ) {
        if( !hitmask[hj] ) continue;
        if( hitlist[hj] == hitlist[hi] ) continue;
        if( hitIsClustered[hj] ) continue;
        float rms = hitlist[hj]->RMS();
        float t = hitinfo[hj].Time;
        float t1 = t-rms;
        float t2 = t+rms;
        if( (t1 > hc.StartTime && t1 < hc.EndTime )
          ||(t2 > hc.StartTime && t2 < hc.EndTime ) ){
          BlipUtils::GrowHitClust(hitlist[hj],hitinfo[hj],hc);
          hitIsClustered[hj] = true;
        }
      }
      hitclust.push_back(hc);
    }
  }
  std::cout
  <<"Reconstructed "<<hitclust.size()<<" hit clusters:\n";
  for(auto const& hc : hitclust){
    std::cout<<"   * wire: "<<hc.LeadWire<<", charge: "<<hc.Charge<<", time: "<<hc.Time<<", fullwidth="<<hc.EndTime-hc.StartTime<<", hitMult="<<hc.HitIDs.size()<<"\n";
  }
   
   

  //====================================
  // Fill TTree
  //====================================
  fData->FillTree();
}



//###################################################
//  Printouts for debugging
//###################################################

void sbnd::BlipAna::PrintParticleInfo(size_t i){
  printf("  %5i  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, KE0=%8.1f, Edep=%8.3f, T=%8.1f --> %8.1f, moth=%5i, %12s, ND=%i\n",
   (int)i,
   fData->trackID[i],
   fData->pdg[i],
   fData->startPointx[i],
   fData->startPointy[i],
   fData->startPointz[i],
   fData->pathlen[i], 
   fData->E[i]-fData->mass[i],
   fData->depEnergy[i],
   fData->startT[i],
   fData->endT[i],
   fData->mother[i],
   fData->process[i].c_str(),
   fData->nDaughters[i]
  ); 
}



DEFINE_ART_MODULE(sbnd::BlipAna)

#endif
