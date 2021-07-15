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
//#include "sbndcode/RecoUtils/RecoUtils.h"
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
const int kMaxHits  = 100000;
const int kMaxBlips = 5000;
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

    // --- Configurations and switches
    std::string treeName            = "anatree";
    bool  saveParticleList          = true;

    // --- Event information ---   
    int   event;                    // event number
    int   run;                      // run number

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
    int   edep_clustid[kMaxEDeps];  // hitclust ID
    int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
    float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
    float edep_charge[kMaxEDeps];   // total electrons reaching anode wires
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

    // --- Hit cluster information ---
    int   nclusts;
    int   clust_tpc[kMaxHits];
    int   clust_plane[kMaxHits];
    int   clust_wire[kMaxHits];
    int   clust_nwires[kMaxHits];
    int   clust_nhits[kMaxHits];
    float clust_charge[kMaxHits];
    float clust_time[kMaxHits];
    float clust_startTime[kMaxHits];
    float clust_endTime[kMaxHits];
    float clust_g4energy[kMaxHits];
    float clust_g4charge[kMaxHits];
    int   clust_g4id[kMaxHits];
    int   clust_ismatched[kMaxHits];
    int   clust_blipid[kMaxHits];
    int   clust_edepid[kMaxHits];
    
    // --- Blip information ---
    int   nblips;
    int   blip_tpc[kMaxBlips];
    int   blip_bestplane[kMaxBlips];
    int   blip_ncrossings[kMaxBlips];
    float blip_x[kMaxBlips];
    float blip_y[kMaxBlips];
    float blip_z[kMaxBlips];
    float blip_prms[kMaxBlips];
    float blip_charge[kMaxBlips];
    int   blip_edepid[kMaxBlips];

    // === Function for resetting data ===
    void Clear(){ 
      event                 = -999; // --- event-wide info ---
      run                   = -999;
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
      FillWith(edep_charge, -999);
      FillWith(edep_x,      -99999.);
      FillWith(edep_y,      -99999.);
      FillWith(edep_z,      -99999.);
      FillWith(edep_ds,     -999);
      FillWith(edep_g4id,   -9);
      FillWith(edep_clustid,-9);
      FillWith(edep_blipid, -9);
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
      nclusts               = 0;    // --- Hit Clusters ---
      FillWith(clust_tpc,  -9);
      FillWith(clust_plane, -9);
      FillWith(clust_wire, -9);
      FillWith(clust_nwires, -9);
      FillWith(clust_nhits, -9);
      FillWith(clust_charge, -999);
      FillWith(clust_time, -999);
      FillWith(clust_startTime, -999);
      FillWith(clust_endTime, -999);
      FillWith(clust_g4id, -9);
      FillWith(clust_g4charge, -999);
      FillWith(clust_g4energy, -999);
      FillWith(clust_ismatched, 0);
      FillWith(clust_edepid, -9);
      FillWith(clust_blipid, -9);
      nblips                = 0;  // --- Blips ---
      FillWith(blip_tpc,  -9);
      FillWith(blip_bestplane, -9);
      FillWith(blip_x, -99999);
      FillWith(blip_y, -99999);
      FillWith(blip_z, -99999);
      FillWith(blip_prms, -999999);
      FillWith(blip_ncrossings, -9);
      FillWith(blip_charge, -999);
      FillWith(blip_edepid, -9);

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
      tree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
      tree->Branch("event",&event,"event/I");
      tree->Branch("run",&run,"run/I");
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
      tree->Branch("edep_blipid",edep_blipid,"edep_blipid[nedeps]/I"); 
      tree->Branch("edep_clustid",edep_clustid,"edep_clustid[nedeps]/I"); 
      tree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
      tree->Branch("edep_charge",edep_charge,"edep_charge[nedeps]/F"); 
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
      tree->Branch("nclusts",&nclusts,"nclusts/I");
      tree->Branch("clust_tpc",clust_tpc,"clust_tpc[nclusts]/I");
      tree->Branch("clust_plane",clust_plane,"clust_plane[nclusts]/I");
      tree->Branch("clust_wire",clust_wire,"clust_wire[nclusts]/I");
      tree->Branch("clust_nwires",clust_nwires,"clust_nwires[nclusts]/I");
      tree->Branch("clust_nhits",clust_nhits,"clust_nhits[nclusts]/I");
      tree->Branch("clust_charge",clust_charge,"clust_charge[nclusts]/F");
      tree->Branch("clust_time",clust_time,"clust_time[nclusts]/F");
      tree->Branch("clust_startTime",clust_startTime,"clust_startTime[nclusts]/F");
      tree->Branch("clust_endTime",clust_endTime,"clust_endTime[nclusts]/F");
      tree->Branch("clust_g4charge",clust_g4charge,"clust_g4charge[nclusts]/F");
      tree->Branch("clust_g4energy",clust_g4energy,"clust_g4energy[nclusts]/F");
      tree->Branch("clust_ismatched",clust_ismatched,"clust_ismatched[nclusts]/I");
      tree->Branch("clust_edepid",clust_edepid,"clust_edepid[nclusts]/I");
      tree->Branch("clust_blipid",clust_blipid,"clust_blipid[nclusts]/I");
      tree->Branch("nblips",&nblips,"nblips/I");
      tree->Branch("blip_tpc",blip_tpc,"blip_tpc[nblips]/I");
      tree->Branch("blip_bestplane",blip_bestplane,"blip_bestplane[nblips]/I");
      tree->Branch("blip_ncrossings",blip_ncrossings,"blip_ncrossings[nblips]/I");
      tree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
      tree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
      tree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
      tree->Branch("blip_prms",blip_prms,"blip_prms[nblips]/F");
      tree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/F");
      tree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
      
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
    float TickOffset[kNplanes];

    // --- Data and calo objects ---
    BlipAnaTreeDataStruct*  fData;
    calo::CalorimetryAlg    fCaloAlg;

    // --- FCL configs ---
    std::string         fAnaTreeName;
    std::string         fHitModuleLabel;
    std::string         fLArG4ModuleLabel;
    std::string         fTrackModuleLabel;
    bool                fSaveParticleList;
    float               fBlipMergeDist;
    std::vector<float>  fMinHitRMS;
    std::vector<float>  fMaxHitRMS;
    std::vector<int>    fPlanePriority;

//    std::vector<float>  fMinHitPhWidthRatio;

    // --- Histograms ---
    /*
    TH1F*   h_trueblip_ds;
    TH1F*   h_trueblip_energy;
    TH1F*   h_nhits[2][kNplanes];
    TH1F*   h_hitrms[2][kNplanes];
    TH1F*   h_hitph[2][kNplanes];
    */

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
  ,fData              (nullptr)
  ,fCaloAlg           (pset.get< fhicl::ParameterSet >    ("CaloAlg"))
  ,fAnaTreeName       (pset.get< std::string >            ("AnaTreeName",       "anatree"))
  ,fHitModuleLabel    (pset.get< std::string >            ("HitModuleLabel",    "gaushit"))
  ,fLArG4ModuleLabel  (pset.get< std::string >            ("LArG4ModuleLabel",  "largeant"))
  ,fTrackModuleLabel  (pset.get< std::string >            ("TrackModuleLabel",  "pandora"))
  ,fSaveParticleList  (pset.get< bool >                   ("SaveParticleList",  true))
  ,fBlipMergeDist     (pset.get< float >                  ("BlipMergeDist",     0.2))
  ,fMinHitRMS         (pset.get< std::vector< float > >   ("MinHitRMS",         {1.7, 1.8,  1.7}))
  ,fMaxHitRMS         (pset.get< std::vector< float > >   ("MaxHitRMS",         {5,   5,    5}))
  ,fPlanePriority     (pset.get< std::vector< int > >     ("PlanePriority",     {2,   0,    1}))
{
  fData = new BlipAnaTreeDataStruct();
  fData ->saveParticleList = fSaveParticleList;
  fData ->treeName = fAnaTreeName;
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
  // -- Detector Properties --
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  // First plane induction plane will be treated as reference. Offsets for subsequent
  // planes must be subtracted off those times to line up with plane 0.
  TickOffset[0] = 0;
  TickOffset[1] = detProp.GetXTicksOffset(1,0,0)-detProp.GetXTicksOffset(0,0,0);
  TickOffset[2] = detProp.GetXTicksOffset(2,0,0)-detProp.GetXTicksOffset(0,0,0);
  std::cout<<"Electron lifetime = "<<detProp.ElectronLifetime()<<"\n";
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
  
  //=========================================
  // Get data products for this event
  //=========================================
  
  
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
  <<"Found "<<fData->nparticles<<" G4 particles, "
            <<fData->nhits<<" hits from "<<fHitModuleLabel
  <<"\n";



  //====================================
  // Save G4 particle information
  //====================================
 
  // Find total visible energy and number electrons drifted to wires
  BlipUtils::CalcTotalDep(fData->total_depEnergy,fData->total_numElectrons);
  std::cout<<"Total energy deposited: "<<fData->total_depEnergy<<" MeV \n";

  // Create empty vector to save all the "true" blips in the event
  std::vector<BlipUtils::TrueBlip> trueblips;
  
  // Loop through the MCParticles
  sim::ParticleList::const_iterator itPart = plist.begin();
  for(size_t i = 0; (i<plist.size())&&(itPart!=plist.end()); i++){
    const simb::MCParticle* pPart = (itPart++)->second;
    
    // Get important info and do conversions
    int trackID   = pPart->TrackId();
    int isPrimary = (int)(pPart->Process() == "primary");
    float mass    = /*GeV->MeV*/1e3 * pPart->Mass();
    float E       = /*GeV->MeV*/1e3 * pPart->E();
    float endE    = /*GeV->MeV*/1e3 * pPart->EndE();
    float P       = /*GeV->MeV*/1e3 * pPart->Momentum().Vect().Mag();
    float Px      = /*GeV->MeV*/1e3 * pPart->Px();
    float Py      = /*GeV->MeV*/1e3 * pPart->Py();
    float Pz      = /*GeV->MeV*/1e3 * pPart->Pz();
    float pathlen = BlipUtils::PathLength(*pPart);

    float ne  = 0, edep = 0;
    BlipUtils::CalcPartDep(trackID,edep,ne);
    
    // Make true blips 
    BlipUtils::TrueBlip tb = BlipUtils::MakeTrueBlip(trackID);
    if( tb.isValid ) trueblips.push_back(tb);

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
    }
  
  }//endloop over G4 particles



  //====================================
  // Merge and save true blip information
  //====================================
  MergeTrueBlips(trueblips, fBlipMergeDist); 
  fData->nedeps = (int)trueblips.size();
  std::cout<<"Found "<<trueblips.size()<<" true blips:\n";
  for(size_t i=0; i<trueblips.size(); i++ ) {
    trueblips[i].ID       = i;
    fData->edep_energy[i] = trueblips.at(i).Energy;
    fData->edep_charge[i] = trueblips.at(i).NumElectrons;
    fData->edep_ds[i]     = trueblips.at(i).Length;
    fData->edep_x[i]      = trueblips.at(i).Position.X();
    fData->edep_y[i]      = trueblips.at(i).Position.Y();
    fData->edep_z[i]      = trueblips.at(i).Position.Z();
    fData->edep_g4id[i]   = trueblips.at(i).LeadG4ID;
    std::cout
    <<"   ~ "<<trueblips.at(i).Energy<<" MeV, "
    <<" ds= "<<trueblips.at(i).Length<<" cm, "
    <<" trkID= "<<trueblips.at(i).LeadG4ID<<"\n";
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
    int trkid = -9;
    int g4id = -9;
    std::set<int> g4ids;
    float g4frac, g4energy, g4charge;
    float charge = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),hitlist[i]->WireID().Plane);
    float time = hitlist[i]->PeakTime()-TickOffset[hitlist[i]->WireID().Plane];
  
    // Find associated track
    if (fmtk.isValid()){
      if (fmtk.at(i).size())  trkid = fmtk.at(i)[0]->ID();
    }
    
    // Find G4 particle ID for leading contributor
    if( BlipUtils::DoesHitHaveSimChannel(hitlist[i]) ){
      BlipUtils::HitTruth( hitlist[i], g4id, g4frac, g4energy, g4charge);
      g4ids = BlipUtils::HitTruthIds(hitlist[i]);
    }

    // Add to the wire-by-wire hit map
    wirehitsMap[wire].push_back(i);

    // Bundle the info to access later
    hitinfo[i].hitid    = i;
    hitinfo[i].hit      = hitlist[i];
    hitinfo[i].g4ids    = g4ids;
    hitinfo[i].trkid    = trkid;
    hitinfo[i].g4id     = g4id;
    hitinfo[i].g4frac   = g4frac;
    hitinfo[i].g4energy = g4energy;
    hitinfo[i].g4charge = g4charge;
    hitinfo[i].Charge   = charge;
    hitinfo[i].Time     = time;
    bool isreal         = (g4id > 0); // determine realness
    hitinfo[i].isreal   = isreal;
    if( isreal ) {                    // group in overlapping hits
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
 
  // Flag real hits
  for(size_t i=0; i<hitlist.size(); i++){
    fData->hit_isreal[i] = hitinfo[i].isreal;
  }


  //=================================================================
  // Blip Reconstruction
  //================================================================
  //  
  //  Will eventually move these into separate alg class.
  //
  //  Procedure
  //  [x] Look for hits that were not included in a track 
  //  [x] Filter hits based on hit width, etc
  //  [x] Merge together closely-spaced hits on same wires, save average peakT +/- spread
  //  - Merge together clusters on adjacent wires (if they match up in time)
  //  - Plane-to-plane time matching
  //  - Wire intersection check to get XYZ
  //  - Create "blip" object and save to tree (nblips, blip_xyz, blip_charge, blip_g4energy)

  // Create a series of masks that we'll update as we go along
  std::vector<bool> hitPassesCuts(hitlist.size(), false);
  std::vector<bool> hitIsClustered(hitlist.size(),false);
  std::vector<bool> hitIsBlipped(hitlist.size(),  false);
  
  // Check for hit quality via RMS and proximity to other good hits
  for(size_t i=0; i<hitlist.size(); i++){
    int plane = hitlist[i]->WireID().Plane;
    if( hitlist[i]->RMS() > fMinHitRMS[plane] &&
        hitlist[i]->RMS() < fMaxHitRMS[plane] ) {
      hitPassesCuts[i] = true;
      for(auto j : wirehitsMap[hitlist[i]->WireID().Wire] ) {
        if( BlipUtils::DoHitsOverlap(hitlist[i],hitlist[j]) ) 
          hitPassesCuts[j] = true;
      }
    }
  }
  // Exclude any hits that were part of tracks
  for(size_t i=0; i<hitlist.size(); i++){
    if( hitinfo[i].trkid >= 0 ) hitPassesCuts[i] = false;
  }

  std::cout
  <<"We selected "<<std::count(hitPassesCuts.begin(),hitPassesCuts.end(),true)<<" hits based on width and track assn cuts\n";

  // ---------------------------------------------------
  // Create collection of hit clusters on same wires
  // ---------------------------------------------------
  std::vector<BlipUtils::HitClust> hitclust;
  std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
  
  for(auto const& wirehits : wirehitsMap){
    for(auto const& hi : wirehits.second ){
      if( !hitPassesCuts[hi] || hitIsClustered[hi] ) continue;
      // cluster this hit
      BlipUtils::HitClust hc = BlipUtils::MakeHitClust(hitinfo[hi]);
      if( !hc.isValid ) continue;
      hitIsClustered[hi] = true;
      // see if we can add other hits to it; continue until 
      // no new hits can be lumped in with this clust
      int hitsAdded;
      do{
        hitsAdded = 0;  
        for(auto const& hj : wirehits.second ) {
          if( !hitPassesCuts[hj] || hitIsClustered[hj] ) continue;
          if( hitlist[hj] == hitlist[hi] ) continue;
          float rms = hitlist[hj]->RMS();
          float t1 = hitinfo[hj].Time - rms;
          float t2 = hitinfo[hj].Time + rms;
          if( (t1 > hc.StartTime && t1 < hc.EndTime )
            ||(t2 > hc.StartTime && t2 < hc.EndTime ) ){
            BlipUtils::GrowHitClust(hitinfo[hj],hc);
            hitIsClustered[hj] = true;
            hitsAdded++;
          }
        }
      } while ( hitsAdded!=0 );
      hitclust.push_back(hc);
      tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(hitclust.size()-1);
//      wireclustsMap[wire].push_back(hitclust.size()-1);
    }
  }
  //std::cout<<"Reconstructed "<<hitclust.size()<<" hit clusters:\n";
  //for(auto const& hc : hitclust)
  //  std::cout<<"   * wire: "<<hc.LeadHitWire<<", start/end wire: "<<hc.StartWire<<"/"<<hc.EndWire<<"  (span: "<<hc.Wires.size()<<"), plane: "<<hc.Plane<<", time: "<<hc.Time<<", half-width="<<(hc.EndTime-hc.StartTime)/2.<<"\n";

  // Look for clusters on adjacent wires (but same plane) and merge them together 
  // to account for scenarios where charge from a single blip is shared between wires
  std::vector<BlipUtils::HitClust> hitclust_merged;
//  std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap_merged;
  for(auto const& tpcMap : tpc_planeclustsMap ) {
    for(auto const& planeclusts : tpcMap.second ) {
      for( auto const& ci : planeclusts.second ){
        if( hitclust[ci].isMerged ) continue;
        BlipUtils::HitClust hc = hitclust[ci];
        hitclust[ci].isMerged = true;
        int clustsAdded = 0;
        do{
          int clustsAdded = 0;
          for(auto const& cj : planeclusts.second ){
            if( hitclust[cj].isMerged ) continue;
            int w1 = hitclust[cj].StartWire;
            int w2 = hitclust[cj].EndWire;
            // check if clusters are adjacent
            if( ((w1 - hc.EndWire) <= 1) || ((hc.StartWire-w2) <= 1) ){
              // check if there's a time match
              if( BlipUtils::DoHitClustsMatch(hc,hitclust[cj]) ) {
                hc = BlipUtils::MergeHitClusts(hc,hitclust[cj]);
                if( hitclust[cj].isMerged ) { clustsAdded++;}
              }
            }
          }
        } while ( clustsAdded!=0);
        hitclust_merged.push_back(hc);
  //      tpc_planeclustsMap_merged[hc.TPC][hc.Plane].push_back(hitclust_merged.size()-1);
      }
    }
  }
  hitclust = hitclust_merged;
  std::cout<<"After merging: "<<hitclust_merged.size()<<" hit clusts:\n";
  for(auto const& hc : hitclust_merged)
    std::cout<<"   * wire: "<<hc.LeadHitWire<<", start/end wire: "<<hc.StartWire<<"/"<<hc.EndWire<<"  (span: "<<hc.Wires.size()<<"), plane: "<<hc.Plane<<", time: "<<hc.Time<<", half-width="<<(hc.EndTime-hc.StartTime)/2.<<"\n";

  //--------------------------------------------
  // Save hit cluster info
  //--------------------------------------------
  fData->nclusts = (int)hitclust.size();
  for(size_t i=0; i<hitclust.size(); i++){
    hitclust[i].ID          = i;
    fData->clust_tpc[i]     = hitclust[i].TPC;
    fData->clust_plane[i]   = hitclust[i].Plane;
    fData->clust_wire[i]    = hitclust[i].LeadHitWire;
    fData->clust_nwires[i]  = (int)hitclust[i].Wires.size();
    fData->clust_nhits[i]   = (int)hitclust[i].HitIDs.size();
    fData->clust_charge[i]  = hitclust[i].Charge;
    fData->clust_time[i]    = hitclust[i].Time;
    fData->clust_startTime[i]=hitclust[i].StartTime;
    fData->clust_endTime[i] =hitclust[i].EndTime;
    // tag associated hits and get true G4 energy/charge
    for(auto const& hitID : hitclust[i].HitIDs)
      fData->hit_clustid[hitID] = i;
    // Find associated true energy dep
    //fData->clust_g4energy[i] = 0;
    //fData->clust_g4charge[i] = 0;
    for(size_t j=0; j< trueblips.size(); j++){
      int tbG4 = trueblips[j].LeadG4ID;
      if( tbG4 > 0 && tbG4 == hitclust[i].LeadHitG4ID ) {
        fData->edep_clustid[j] = hitclust[i].ID;
        fData->clust_edepid[i] = trueblips[j].ID;
        fData->clust_g4energy[i] = trueblips[j].Energy; 
        fData->clust_g4charge[i] = trueblips[j].NumElectrons;
        break;
      }
    }
  }

  // Now for the fun stuff! Match hit clusters in time between the different planes 
  // in order to make Blips. We will require only a match between at least 2 planes
  // for now (will make this fhicl-configurable eventually).
  // ---------------------------------------------------
  // Create collection of Blips
  // ---------------------------------------------------
  std::vector<BlipUtils::Blip> blips;
//  std::vector<bool> isClustMatched(hitclust.size(),false);
//  std::vector<int>  clustBlipID(hitclust.size(),-9);
  //std::cout<<"Looping over clusters to make blips...\n";
  for(size_t i=0; i<hitclust.size(); i++){
    int iPlane = hitclust[i].Plane;
    int iTPC = hitclust[i].TPC;
//    if( isClustMatched[i] ) continue;
    if( hitclust[i].isMatched ) continue;
    
    float t1 = hitclust[i].StartTime;
    float t2 = hitclust[i].EndTime;
    std::vector<BlipUtils::HitClust> hcgroup;
    //std::vector<int> hcindex;
    bool isPlaneMatched[3] = {false,false,false};
    //std::cout<<"    cl "<<i<<" on plane "<<iPlane<<"   t1-t2: "<<t1<<"   "<<t2<<"\n";
    for(size_t j=0; j<hitclust.size(); j++){
      int jPlane = hitclust[j].Plane;
      int jTPC = hitclust[j].TPC;
      if( hitclust[j].isMatched || isPlaneMatched[jPlane] ) continue;
      if( i==j || iPlane==jPlane || jTPC!=iTPC )        continue;
      //std::cout<<"     does cl "<<j<<" on plane "<<jPlane<<" match up? "<<hitclust[j].StartTime<<"  "<<hitclust[j].EndTime<<"\n";
      
      if( BlipUtils::DoHitClustsMatch(hitclust[j],t1,t2)) {
      //  std::cout<<"    yep!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        if( !hitclust[i].isMatched ) {
          hitclust[i].isMatched = true;
//          isClustMatched[i] = true;
          isPlaneMatched[iPlane] = true;
          hcgroup.push_back(hitclust[i]);
          //hcindex.push_back(i);
        }
        if( !hitclust[j].isMatched ) {
          hitclust[j].isMatched = true;
          //isClustMatched[j] = true;
          isPlaneMatched[jPlane] = true;
          hcgroup.push_back(hitclust[j]);
        //  hcindex.push_back(j);
        }
        t1 = std::min(hitclust[i].StartTime,hitclust[j].StartTime);
        t2 = std::max(hitclust[i].EndTime,hitclust[j].EndTime);
    //    std::cout<<" UPdating t1/t2 --> "<<t1<<"  "<<t2<<"\n";
      }
    }
    if(hcgroup.size()) {
  //    std::cout<<"Attempting to make blip from "<<hcgroup.size()<<" matches...\n";
      BlipUtils::Blip newBlip = BlipUtils::MakeBlip(hcgroup);
      if( newBlip.isValid ){
        blips.push_back(newBlip);
        for(auto hc : hcgroup ) hitclust[hc.ID].BlipID = blips.size()-1;
      }
    }

  }//end loop over hitclusts

  std::cout<<"Reconstructed "<<blips.size()<<" blips:\n";
  for(auto const& b : blips)
    std::cout
    <<"   -- TPC: "<<b.TPC
    <<"; charge: "<<b.Charge[0]<<" ,"<<b.Charge[1]<<" ,"<<b.Charge[2]
    <<"; Position: "<<b.Position.X()<<", "<<b.Position.Y()<<", "<<b.Position.Z()
    <<"; LocRMS: "<<b.PositionRMS
    <<"; NCross: "<<b.NCrossings
    <<"\n";
 
 
  // Save blip info to tree
  fData->nblips = blips.size();
  for(size_t i=0; i<blips.size(); i++){

    // determine "best" plane
    size_t bestplane = 2;
    for(auto prefPlane : fPlanePriority ) 
      if( blips[i].Charge[prefPlane] > 0 ) { bestplane = prefPlane; break;}

    fData->blip_tpc[i]        = blips[i].TPC;
    fData->blip_bestplane[i]  = bestplane;
    fData->blip_charge[i]     = blips[i].Charge[bestplane];
    fData->blip_x[i]          = (float)blips[i].Position.X();
    fData->blip_y[i]          = (float)blips[i].Position.Y();
    fData->blip_z[i]          = (float)blips[i].Position.Z();
    fData->blip_prms[i]       = blips[i].PositionRMS;
    fData->blip_ncrossings[i] = blips[i].NCrossings;
    
    // find associated true edep
    float max = 0;
    for(auto clustID : blips[i].ClustIDs ) {
      int edepid = fData->clust_edepid[clustID];
      float E = trueblips[edepid].Energy;
      if( E > max ) {
        fData->blip_edepid[i] = edepid;
        fData->edep_blipid[edepid] = i;
        max = E;
      }
    }

  }

  // Update clust data in Tree, so we can tie reconstructed hit clusters
  // to the blips they were grouped into
  for(size_t i=0; i<hitclust.size(); i++){
    fData->clust_ismatched[i] = hitclust[i].isMatched;
    fData->clust_blipid[i] = hitclust[i].BlipID;
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
