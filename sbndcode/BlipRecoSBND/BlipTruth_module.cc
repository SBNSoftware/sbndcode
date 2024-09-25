//#####################################################################
//###  BlipTruth analyzer module
//#####################################################################
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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "cetlib/search_path.h"

// SBND-specific includes
//#include "sbndcode/BlipRecoSBND/Alg/BlipRecoAlg.h"
#include "sbndcode/BlipRecoSBND/Utils/BlipUtils.h"

// C++ includes
#include <cstring>
#include <utility>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cmath>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

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
const int kMaxG4      = 100000;
const int kMaxEDeps   = 10000;

class BlipTruth;
  
//###################################################
//  Data storage structure
//###################################################
class BlipTruthTreeDataStruct 
{
  public:

  // --- TTrees
  TTree* evtTree;

  // --- Configurations and switches ---
  std::string treeName      = "tree";

  // --- Event information ---   
  int           event;                // event number
  int           run;                  // run number
  int           subrun;               // subrun number

  // --- G4 information ---
  int   nparticles;               // number of G4 particles
  bool  part_isPrimary[kMaxG4];        // is primary particle
  int   part_trackID[kMaxG4];          // G4 track ID
  int   part_pdg[kMaxG4];              // PDG
  int   part_nDaughters[kMaxG4];       // number of daughters
  int   part_mother[kMaxG4];           // mother particle
  float part_E[kMaxG4];                // initial energy (MeV)
  float part_KE[kMaxG4];               // initial kinetic energy (MeV)
  float part_endE[kMaxG4];             // final energy (MeV)
  float part_endKE[kMaxG4];             // final energy (MeV)
  float part_mass[kMaxG4];             // mass (MeV)
  float part_P[kMaxG4];                // momentum (MeV)
  float part_Px[kMaxG4];               // momentum x (MeV)
  float part_Py[kMaxG4];               // momentum y (MeV)
  float part_Pz[kMaxG4];               // momentum z (MeV)
  float part_startPointx[kMaxG4];      // starting x (cm)
  float part_startPointy[kMaxG4];      // starting y (cm)
  float part_startPointz[kMaxG4];      // starting y (cm)
  float part_endPointx[kMaxG4];        // ending x (cm)
  float part_endPointy[kMaxG4];        // ending y (cm)
  float part_endPointz[kMaxG4];        // ending y (cm)
  float part_startT[kMaxG4];           // starting time (us)
  float part_endT[kMaxG4];             // ending time (us)
  float part_pathlen[kMaxG4];          // path length (cm)
  int   part_numTrajPts[kMaxG4];       // number traj points
  float part_depEnergy[kMaxG4];        // energy deposited in AV (MeV)
  int   part_depElectrons[kMaxG4];     // electrons deposited
  std::vector<std::string> part_process;// process name

  // --- True energy deposit info (derived from SimChannels and SimEnergyDeposits) ---
  int   nedeps;                   // number of true localized energy depositions
  int   edep_tpc[kMaxEDeps];      // TPC
  int   edep_g4index[kMaxEDeps];     // leading G4 index ("part_variable[g4id]")
  int   edep_g4trkid[kMaxEDeps];  // leading G4 track ID ("part_trackID")
  float edep_g4qfrac[kMaxEDeps];  // fraction of total charge from lead particle
  bool  edep_isPrimary[kMaxEDeps];// matched to a primary generated particle?
  int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
  int   edep_blipid[kMaxEDeps];   // reconstructed 3D blip ID
  float edep_energy[kMaxEDeps];   // total energy deposited [MeV]
  int   edep_electrons[kMaxEDeps];// total ionization electrons deposited (e-)
  int   edep_charge[kMaxEDeps];   // total electrons reaching anode wires (e-)
  int   edep_tdrift[kMaxEDeps];   // drift time for this energy dep (us)
  float edep_x[kMaxEDeps];        // x (cm)
  float edep_y[kMaxEDeps];        // y (cm)
  float edep_z[kMaxEDeps];        // z (cm)
  float edep_pathlen[kMaxEDeps];  // cm
  int   edep_proc[kMaxEDeps];     // encodes particle process
                                  //  0 = primary
                                  //  1 = compton scatter ("compt")
                                  //  2 = photoelectric effect ("phot")
                                  //  3 = e+e- pair production ("conv")
                                  //  4 = other
  
  
  
  
  // === Function for resetting data ===
  void Clear(){ 
    event                 = -999; // --- event-wide info ---
    run                   = -999;
    subrun                = -999; 
    nparticles            = 0;    // --- G4 particles ---
    FillWith(part_isPrimary,   false);
    //FillWith(part_madeHitCol,        false);
    FillWith(part_trackID,     -999);
    FillWith(part_pdg,         -99999);
    FillWith(part_nDaughters,  -999);
    FillWith(part_mother,      -999);
    FillWith(part_E,           -999.);
    FillWith(part_endE,        -999.);
    FillWith(part_KE,           -999.);
    FillWith(part_endKE,        -999.);
    FillWith(part_mass,        -999.);
    FillWith(part_P,           -999.);
    FillWith(part_Px,          -999.);
    FillWith(part_Py,          -999.);
    FillWith(part_Pz,          -999.);
    FillWith(part_startPointx, -99999.);
    FillWith(part_startPointy, -99999.);
    FillWith(part_startPointz, -99999.);
    FillWith(part_endPointx,   -99999.);
    FillWith(part_endPointy,   -99999.);
    FillWith(part_endPointz,   -99999.);
    FillWith(part_startT,      -99999.);
    FillWith(part_endT,        -99999.);
    FillWith(part_pathlen,     -999.);
    FillWith(part_numTrajPts,   -9);
    FillWith(part_depElectrons,-999);
    FillWith(part_depEnergy,   -999.);
    FillWith(part_process,     "");
    nedeps                = 0;    // --- EDeps ---
    FillWith(edep_tpc,    -9);
    FillWith(edep_energy, -999);
    FillWith(edep_electrons,  -999);
    FillWith(edep_charge, -999);
    FillWith(edep_tdrift, -999);
    FillWith(edep_x,      -99999.);
    FillWith(edep_y,      -99999.);
    FillWith(edep_z,      -99999.);
    FillWith(edep_pathlen,      -99999.);
    FillWith(edep_g4trkid,  -9);
    FillWith(edep_g4index,     -9);
    FillWith(edep_g4qfrac,  -9);
    FillWith(edep_pdg,   -999);
    FillWith(edep_proc,   -9);
    FillWith(edep_isPrimary, false);
  }

  // === Function for resizing vectors (if necessary) ===
  // To be called after numbers of hits/tracks/particles
  // in the event has been determined
  void Resize() {
    if(nparticles) part_process.assign(nparticles,"");
  }
      
  // === Function for initializing tree branches ===
  void MakeTree(){
    
    art::ServiceHandle<art::TFileService> tfs;
   
    evtTree = tfs->make<TTree>(treeName.c_str(),"analysis tree");
    evtTree->Branch("event",&event,"event/I");
    evtTree->Branch("run",&run,"run/I");
    evtTree->Branch("subrun",&subrun,"subrun/I");

    evtTree->Branch("nparticles",&nparticles,"nparticles/I");
    evtTree->Branch("part_isPrimary",part_isPrimary,"part_isPrimary[nparticles]/O");
    evtTree->Branch("part_trackID",part_trackID,"part_trackID[nparticles]/I");
    evtTree->Branch("part_pdg",part_pdg,"part_pdg[nparticles]/I");
    evtTree->Branch("part_nDaughters",part_nDaughters,"part_nDaughters[nparticles]/I");
    evtTree->Branch("part_mother",part_mother,"part_mother[nparticles]/I");
    evtTree->Branch("part_KE",part_KE,"part_KE[nparticles]/F");
    evtTree->Branch("part_startT",part_startT,"part_startT[nparticles]/F");
    evtTree->Branch("part_endT",part_endT,"part_endT[nparticles]/F");
    evtTree->Branch("part_startPointx",part_startPointx,"part_startPointx[nparticles]/F");
    evtTree->Branch("part_startPointy",part_startPointy,"part_startPointy[nparticles]/F");
    evtTree->Branch("part_startPointz",part_startPointz,"part_startPointz[nparticles]/F");
    evtTree->Branch("part_pathlen",part_pathlen,"part_pathlen[nparticles]/F");
    evtTree->Branch("part_depEnergy",part_depEnergy,"part_depEnergy[nparticles]/F");
    evtTree->Branch("part_depElectrons",part_depElectrons,"part_depElectrons[nparticles]/I");
    evtTree->Branch("part_process",&part_process);
      evtTree->Branch("nedeps",&nedeps,"nedeps/I");
      evtTree->Branch("edep_g4index",edep_g4index,"edep_g4index[nedeps]/I"); 
      evtTree->Branch("edep_g4trkid",edep_g4trkid,"edep_g4trkid[nedeps]/I"); 
      evtTree->Branch("edep_g4qfrac",edep_g4qfrac,"edep_g4qfrac[nedeps]/F"); 
      evtTree->Branch("edep_isPrimary",edep_isPrimary,"edep_isPrimary[nedeps]/O"); 
      evtTree->Branch("edep_pdg",edep_pdg,"edep_pdg[nedeps]/I"); 
      evtTree->Branch("edep_proc",edep_proc,"edep_proc[nedeps]/I"); 
      evtTree->Branch("edep_blipid",edep_blipid,"edep_blipid[nedeps]/I"); 
      evtTree->Branch("edep_energy",edep_energy,"edep_energy[nedeps]/F"); 
      evtTree->Branch("edep_electrons",edep_electrons,"edep_electrons[nedeps]/I"); 
      evtTree->Branch("edep_charge",edep_charge,"edep_charge[nedeps]/I"); 
      evtTree->Branch("edep_tdrift",edep_tdrift,"edep_tdrift[nedeps]/I"); 
      evtTree->Branch("edep_x",edep_x,"edep_x[nedeps]/F"); 
      evtTree->Branch("edep_y",edep_y,"edep_y[nedeps]/F"); 
      evtTree->Branch("edep_z",edep_z,"edep_z[nedeps]/F");
      evtTree->Branch("edep_pathlen",edep_pathlen,"edep_pathlen[nedeps]/F");
  }


};//BlipTruthTreeDataStruct class



//###################################################
//  BlipTruth class definition
//###################################################
class BlipTruth : public art::EDAnalyzer 
{ 
  public:
  explicit BlipTruth(fhicl::ParameterSet const& pset);
  virtual ~BlipTruth();
  
  //void beginJob();                      // called once, at start of job
  void endJob();                        // called once, at end of job
  void analyze(const art::Event& evt);  // called per event

  private:
  void    PrintParticleInfo(size_t);
  void    PrintTrueBlipInfo(const blip::TrueBlip&);

  // --- Data and calo objects ---
  BlipTruthTreeDataStruct*  fData;
  //blip::BlipRecoAlg*      fBlipAlg;

  // --- FCL configs ---
  std::string         fGeantProducer;
  std::string         fSimDepProducer;
  int                 fCaloPlane;
  float               fTrueBlipMergeDist;
  
  int                 fMaxPrints;

  // --- Counters and such ---
  int   fNumEvents          = 0;
  int   fNumPrints          = 0;

  // --- Histograms ---
  TH1D*   h_part_process;

  // Initialize histograms
  void InitializeHistograms(){
    
    art::ServiceHandle<art::TFileService> tfs;

    // MC histograms related to truth
    h_part_process    = tfs->make<TH1D>("part_process","MCParticle->Process()",5,0,5);
    auto xa = h_part_process->GetXaxis();
      xa->SetBinLabel(1,"primary");
      xa->SetBinLabel(2,"compt");
      xa->SetBinLabel(3,"phot");
      xa->SetBinLabel(4,"conv");
      xa->SetBinLabel(5,"other");
    
    
 }

};//class BlipTruth



//###################################################
//  BlipTruth constructor and destructor
//###################################################
BlipTruth::BlipTruth(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset)
  ,fData  (nullptr)
{
  // blip reconstruction algorithm class
  fGeantProducer  = pset.get<std::string>            ("GeantProducer",   "largeant");
  fSimDepProducer = pset.get<std::string>           ("SimEDepProducer",  "ionandscint");
  fCaloPlane          = 2;
  fTrueBlipMergeDist  = 0.3;
  fMaxPrints          = 100;

  // data tree object
  fData = new BlipTruthTreeDataStruct();
  fData ->treeName = pset.get<std::string> ("tree", "tree");
  fData ->Clear();
  fData ->MakeTree();

  // initialize histograms
  InitializeHistograms();
    
}
BlipTruth::~BlipTruth(){}



//###################################################
//  Main event-by-event analysis
//###################################################
void BlipTruth::analyze(const art::Event& evt)
{ 
  //============================================
  // New event!
  //============================================
  fData            ->Clear();
  fData->event      = evt.id().event();
  fData->run        = evt.id().run();
  fData->subrun     = evt.id().subRun();
  fNumEvents++;

  auto const detProp    = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  
  // Tell us what's going on!
  if( fNumEvents < 200 || (fNumEvents % 100) == 0 ) {
    std::cout<<"\n"
    <<"=========== BlipTruth =========================\n"
    <<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"; total events processed: "<<fNumEvents<<"\n";
  }
  
  //=======================================
  // Get data products for this event
  //========================================
  
  // -- G4 particles
  art::Handle< std::vector<simb::MCParticle> > pHandle;
  std::vector<art::Ptr<simb::MCParticle> > plist;
  if (evt.getByLabel(fGeantProducer,pHandle))
    art::fill_ptr_vector(plist, pHandle);

  // -- SimEnergyDeposits
  art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
  std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
  if (evt.getByLabel(fSimDepProducer,sedHandle))
      art::fill_ptr_vector(sedlist, sedHandle); 
 
  // Resize data struct objects
  fData->nparticles = (int)plist.size();
  fData->Resize();

  if( plist.size()==0 ) return;
  
  bool printout = (fNumEvents < fMaxPrints);

  //====================================
  // Keep tabs on total energy, charge,
  // and electrons arriving at anode
  //====================================
  float total_depEnergy         = 0;
  float total_depElectrons      = 0;
  //float total_numElectrons      = 0;

  std::vector<blip::TrueBlip>     trueblips;
  std::vector<blip::ParticleInfo> pinfo;
  pinfo.resize(plist.size());

  //==================================================
  // Use G4 information to determine the "true" blips in this event.
  //==================================================
  for(size_t i = 0; i<plist.size(); i++){
    BlipUtils::FillParticleInfo( *plist[i], pinfo[i], sedlist, fCaloPlane);
    pinfo[i].index = i;
  }
  BlipUtils::MakeTrueBlips(pinfo, trueblips);
  BlipUtils::MergeTrueBlips(trueblips, fTrueBlipMergeDist);
  

  //====================================
  // Save MCParticle information
  //====================================
  std::map<int,int> map_g4trkid_index;
  if( plist.size() ) {
    
    // Loop through the MCParticles
    if( printout ){
      std::cout<<"\n";
      std::cout<<"-------------------------------------------------------------\n";
      std::cout<<"Looping over G4 MCParticles in event "<<evt.id().event()<<"\n";
    }
    for(size_t i = 0; i<plist.size(); i++){
      auto& pPart = plist[i];
      map_g4trkid_index[pPart->TrackId()] = i;
      total_depEnergy       += pinfo[i].depEnergy;
      total_depElectrons    += pinfo[i].depElectrons;
      
      // Save to TTree object
      if(i<kMaxG4){
        fData->part_trackID[i]         = pPart->TrackId();
        fData->part_pdg[i]             = pPart->PdgCode();
        fData->part_nDaughters[i]      = pPart->NumberDaughters();
        fData->part_mother[i]          = pPart->Mother();
        fData->part_E[i]               = pinfo[i].E;
        fData->part_endE[i]            = pinfo[i].endE;
        fData->part_mass[i]            = pinfo[i].mass;
        fData->part_KE[i]              = pinfo[i].KE;
        fData->part_endKE[i]           = pinfo[i].endKE;
        fData->part_P[i]               = pinfo[i].P;
        fData->part_Px[i]              = pinfo[i].Px;
        fData->part_Py[i]              = pinfo[i].Py;
        fData->part_Pz[i]              = pinfo[i].Pz;
        fData->part_startPointx[i]     = pPart->Vx();
        fData->part_startPointy[i]     = pPart->Vy();
        fData->part_startPointz[i]     = pPart->Vz();
        fData->part_endPointx[i]       = pPart->EndPosition()[0];
        fData->part_endPointy[i]       = pPart->EndPosition()[1];
        fData->part_endPointz[i]       = pPart->EndPosition()[2];
        fData->part_startT[i]          = pinfo[i].time;
        fData->part_endT[i]            = pinfo[i].endtime;
        fData->part_pathlen[i]         = pinfo[i].pathLength;
        fData->part_numTrajPts[i]      = pinfo[i].numTrajPts;
        fData->part_process[i]         = pPart->Process();
        fData->part_depEnergy[i]       = pinfo[i].depEnergy;
        fData->part_depElectrons[i]    = pinfo[i].depElectrons;
        fData->part_isPrimary[i]       = pinfo[i].isPrimary;
        
        
        if( printout ) PrintParticleInfo(i);
      }
    } // endloop over G4 particles
   
    if( printout ) std::cout<<"True total energy deposited: "<<total_depEnergy<<" MeV \n";
  
  }//endif particles found in event
  

  //====================================
  // Save TrueBlip information
  //====================================
  fData->nedeps = (int)trueblips.size();
  if( trueblips.size() ) {
    if( printout ) {
      std::cout<<"\n";
      std::cout<<"-------------------------------------------------------------\n";
      std::cout<<"Looping over true energy depositions ('trueblips'):\n";
    }
    for(auto& trueblip : trueblips ) {
      int i     = trueblip.ID;
      int ig4   = trueblip.LeadG4Index;
      auto& pPart = plist[ig4];
      fData->edep_tpc[i]      = trueblip.TPC;
      fData->edep_energy[i]   = trueblip.Energy;
      fData->edep_electrons[i]= trueblip.DepElectrons;
      fData->edep_charge[i]   = trueblip.NumElectrons;
      fData->edep_x[i]        = trueblip.Position.X();
      fData->edep_y[i]        = trueblip.Position.Y();
      fData->edep_z[i]        = trueblip.Position.Z();

      fData->edep_tdrift[i]   = trueblip.DriftTime;
      fData->edep_pdg[i]      = trueblip.LeadG4PDG;
      fData->edep_g4trkid[i]  = trueblip.LeadG4ID;
      fData->edep_g4index[i]  = trueblip.LeadG4Index;
      fData->edep_g4qfrac[i]  = trueblip.G4ChargeMap[trueblip.LeadG4ID] / trueblip.DepElectrons;
      fData->edep_isPrimary[i]= (pPart->Process() == "primary");
      fData->edep_pathlen[i]        = trueblip.PathLen;

      int proc_code = -9;
      std::string proc = fData->part_process[ig4];
      if(     proc == "primary") { h_part_process->Fill("primary",1); proc_code = 0; }
      else if(proc == "compt"  ) { h_part_process->Fill("compt",  1); proc_code = 1; }
      else if(proc == "phot"   ) { h_part_process->Fill("phot",   1); proc_code = 2; }
      else if(proc == "conv"   ) { h_part_process->Fill("conv",   1); proc_code = 3; }
      else                       { h_part_process->Fill("other",  1); proc_code = 4; }
      fData->edep_proc[i] = proc_code;
      //float ne_dep  = trueblip.DepElectrons;
      //float ne      = trueblip.NumElectrons;
      //total_numElectrons += ne;

      if( printout ) PrintTrueBlipInfo(trueblip);

    }
  }//endif trueblips were made
  

  //====================================
  // Fill TTree
  //====================================
  fData->evtTree->Fill();

}


//###################################################
//  endJob: output useful info to screen
//###################################################
void BlipTruth::endJob(){
  printf("\n***********************************************\n");
  printf("                      DONE                        ");
  printf("\n***********************************************\n");
}




//###################################################
//  Printouts for debugging
//###################################################

void BlipTruth::PrintParticleInfo(size_t i){
  printf("  %5i  trkID: %-6i PDG: %-10i XYZ= %7.1f %7.1f %7.1f, dL=%7.2f, Npts=%4i, KE0=%8.3f, Edep=%8.3f, T=%10.2f, moth=%5i, %12s, ND=%i\n",
   (int)i,
   fData->part_trackID[i],
   fData->part_pdg[i],
   fData->part_startPointx[i],
   fData->part_startPointy[i],
   fData->part_startPointz[i],
   fData->part_pathlen[i],
   fData->part_numTrajPts[i],
   fData->part_KE[i],
   fData->part_depEnergy[i],
   fData->part_startT[i]/1e3,
   fData->part_mother[i],
   fData->part_process[i].c_str(),
   fData->part_nDaughters[i]
  ); 
}

void BlipTruth::PrintTrueBlipInfo(const blip::TrueBlip& tb){
  printf("  edepID: %5i  G4ID: %-6i PDG: %-10i XYZ: %7.2f, %7.2f, %7.2f, %8.3f MeV, %8i e- deposited,  %12s\n",
   tb.ID,
   tb.LeadG4ID,
   tb.LeadG4PDG,
   tb.Position.X(),
   tb.Position.Y(),
   tb.Position.Z(),
   tb.Energy,
   tb.DepElectrons,
   fData->part_process[tb.LeadG4Index].c_str()
  ); 
}


DEFINE_ART_MODULE(BlipTruth)

#endif
