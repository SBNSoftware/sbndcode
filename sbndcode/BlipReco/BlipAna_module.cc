//#####################################################################
//###  BlipAna analyzer module
//###
//###  Contains algorithms for reconstructing isolated, MeV-scale energy
//###  depositions in the TPC, called "blips." A TTree is made for offline
//###  analysis and plot-making. Algs will eventually be migrated into 
//###  dedicated alg/tool classes as appropriate.
//###
//###  Author: Will Foreman (wforeman_at_iit.edu)
//###  Date:   Sept 2021
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
#include "lardataobj/AnalysisBase/CosmicTag.h"
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
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "cetlib/search_path.h"

// SBNDCode includes
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
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"

// Set global constants and max array sizes
const int kMaxHits  = 100000;
const int kMaxTrks  = 10000;
const int kMaxShwrs = 10000;
const int kMaxBlips = 5000;
const int kMaxG4    = 100000;
const int kMaxEDeps = 100000;
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

    // --- Configurations and switches ---
    std::string treeName            = "anatree";
    bool  saveParticleList          = true;
    
    // --- Event information ---   
    int   event;                    // event number
    int   run;                      // run number

    // --- G4 information ---
    float total_depEnergy;          // total deposited energy in AV
    float total_numElectrons;       // total electrons reaching anode wires
    float gamma_depEnergy;          // total gamma-induced energy deposited
    float gamma_numElectrons;       // total electrons from gamma-induced depositions
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
    int   edep_tpc[kMaxEDeps];      // TPC
    int   edep_g4id[kMaxEDeps];     // leading G4 track ID
    int   edep_pdg[kMaxEDeps];      // leading G4 track PDG
    int   edep_clustid[kMaxEDeps];  // hitclust ID
    int   edep_blipid[kMaxEDeps];   // reconstructed blip ID
    float edep_energy[kMaxEDeps];   // total energy deposited (MeV)
    float edep_energyESTAR[kMaxEDeps];   // total energy deposited (MeV)
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
    float	hit_sumadc[kMaxHits];     // summed ADC
    float hit_charge[kMaxHits];     // reconstructed number of electrons
    int   hit_mult[kMaxHits];       // multiplicity
    int	  hit_trkid[kMaxHits];      // is this hit associated with a reco track?
    int	  hit_shwrid[kMaxHits];      // is this hit associated with a reco shower?
    int   hit_ismatch[kMaxHits];    // does hit have time match on another plane?
    int   hit_isreal[kMaxHits];     // is this hit real?
    int	  hit_g4id[kMaxHits];       // G4 TrackID of leading particle
    float hit_g4frac[kMaxHits];     // fraction of hit energy from leading MCParticle
    float hit_g4energy[kMaxHits];   // true energy
    float hit_g4charge[kMaxHits];   // true number of electrons (drift-attenuated)
    int   hit_clustid[kMaxHits];    // key of HitClust in which hit was included
    int   hit_blipid[kMaxHits];     // key of Blip in which hit was included
    
    // --- Track information ---
    int   ntrks;                    // number tracks
    int   trk_id[kMaxTrks];         // trackID
    int   trk_npts[kMaxTrks];       // number 3D trajectory points
    float trk_length[kMaxTrks];     // track length [cm]
    float trk_startx[kMaxTrks];     // starting X coordinate
    float trk_starty[kMaxTrks];     // starting Y coordinate
    float trk_startz[kMaxTrks];     // starting Z coordinate
    float trk_startd[kMaxTrks];     // starting distance to boundary
    float trk_endx[kMaxTrks];       // ending X coordinate
    float trk_endy[kMaxTrks];       // ending Y coordinate
    float trk_endz[kMaxTrks];       // ending Z coordinate
    float trk_endd[kMaxTrks];       // ending distance to boundary
    //int   trk_origin[kMaxTrks];     // 0: unknown, 1: cosmic, 2: neutrino, 3: supernova, 4: singles
    //int   trk_g4id[kMaxTrks];       // G4 TrackID of leading contributing particle
    //int   trk_g4pdg[kMaxTrks];      // G4 PDG
    //float trk_purity[kMaxTrks];     // track hit purity
    //float trk_pitch[kMaxTrks];      // track pitch
    //float trk_ke[kMaxTrks];         // track kinetic energy
    //int   trk_cosmictag[kMaxTrks];  // cosmic tagg
    //float trk_cosmicscore[kMaxTrks];// cosmic score
    //int   trk_pidpdg[kMaxTrks];     // particle PID Pdg
    //float trk_pidchi[kMaxTrks];     // chisquared for PID
    //int   trk_bestplane[kMaxTrks];  // plane w/ most hits
    
    // --- Shower information ---
    int   nshwrs;                     // number showers
    int   shwr_id[kMaxShwrs];         // shower ID
    float shwr_dirx[kMaxShwrs];       // direction X
    float shwr_diry[kMaxShwrs];       // direction Y
    float shwr_dirz[kMaxShwrs];       // direction Z
    float shwr_startx[kMaxShwrs];     // starting X coordinate
    float shwr_starty[kMaxShwrs];     // starting X coordinate
    float shwr_startz[kMaxShwrs];     // starting X coordinate
    float shwr_length[kMaxShwrs];     // shower length
    float shwr_openangle[kMaxShwrs];  // opening angle

    // --- Hit cluster information ---
    int   nclusts;
    int   clust_tpc[kMaxHits];
    int   clust_plane[kMaxHits];
    int   clust_wire[kMaxHits];
    int   clust_nwires[kMaxHits];
    int   clust_nhits[kMaxHits];
    float clust_charge[kMaxHits];
    float clust_time[kMaxHits];
    float clust_time_w[kMaxHits];
    float clust_time_lh[kMaxHits];
    float clust_startTime[kMaxHits];
    float clust_endTime[kMaxHits];
    float clust_g4energy[kMaxHits];
    float clust_g4charge[kMaxHits];
    int   clust_g4id[kMaxHits];
    int   clust_ismatch[kMaxHits];
    int   clust_blipid[kMaxHits];
    int   clust_edepid[kMaxHits];
    
    // --- Blip information ---
    float total_blip_energy;
    int   nblips;
    int   blip_tpc[kMaxBlips];
    int   blip_nplanes[kMaxBlips];
    int   blip_caloplane[kMaxBlips];
    float blip_x[kMaxBlips];
    float blip_y[kMaxBlips];
    float blip_z[kMaxBlips];
    float blip_maxdiff[kMaxBlips];
    float blip_charge[kMaxBlips];
    float blip_energy[kMaxBlips];
    float blip_energyESTAR[kMaxBlips];
    int   blip_edepid[kMaxBlips];
    int   blip_clustid[kNplanes][kMaxBlips];

    // === Function for resetting data ===
    void Clear(){ 
      event                 = -999; // --- event-wide info ---
      run                   = -999;
      total_depEnergy       = -999;
      total_numElectrons    = -999;
      gamma_depEnergy       = -999;
      gamma_numElectrons    = -999;
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
      FillWith(edep_tpc,    -9);
      FillWith(edep_energy, -999);
      FillWith(edep_energyESTAR, -999);
      FillWith(edep_charge, -999);
      FillWith(edep_x,      -99999.);
      FillWith(edep_y,      -99999.);
      FillWith(edep_z,      -99999.);
      FillWith(edep_ds,     -999);
      FillWith(edep_g4id,   -9);
      FillWith(edep_pdg,   -999);
      FillWith(edep_clustid,-9);
      FillWith(edep_blipid, -9);
      nhits                 = 0;    // --- TPC hits ---
      FillWith(hit_tpc,     -9);
      FillWith(hit_plane,   -9);
      FillWith(hit_wire,    -999);
      FillWith(hit_channel, -999);
      FillWith(hit_peakT,   -999);
      FillWith(hit_time,    -999);
      FillWith(hit_rms,     -999);
      FillWith(hit_ph,      -999);
      FillWith(hit_area,    -999);
      FillWith(hit_sumadc,  -999);
      FillWith(hit_mult,    -999);
      FillWith(hit_charge,  -999);
      FillWith(hit_isreal,  -9);
      FillWith(hit_ismatch, -9);
      FillWith(hit_trkid,   -9);
      FillWith(hit_shwrid,  -9);
      FillWith(hit_g4id,    -999);
      FillWith(hit_g4frac,  -9);
      FillWith(hit_g4energy,-999);
      FillWith(hit_g4charge,-999);
      FillWith(hit_clustid, -9);
      FillWith(hit_blipid,  -9);
      ntrks                 = 0;    // --- Tracks --- 
      FillWith(trk_id,      -999); 
      FillWith(trk_npts,    -999); 
      FillWith(trk_length,  -999);    
      FillWith(trk_startx,  -999);    
      FillWith(trk_starty,  -999);    
      FillWith(trk_startz,  -999);    
      FillWith(trk_startd,  -999);   
      FillWith(trk_endx,    -999);      
      FillWith(trk_endy,    -999);      
      FillWith(trk_endz,    -999);      
      FillWith(trk_endd,    -999);      
      //FillWith(trk_origin,  -999);
      //FillWith(trk_g4id,    -999);      
      //FillWith(trk_g4pdg,   -999);     
      //FillWith(trk_purity,  -999);    
      //FillWith(trk_pitch,   -999);     
      //FillWith(trk_ke,      -999);        
      //FillWith(trk_cosmictag,-999); 
      //FillWith(trk_cosmicscore,-999);
      //FillWith(trk_pidpdg,  -999);     
      //FillWith(trk_pidchi,  -999);   
      //FillWith(trk_bestplane,-999);  
      nshwrs                = 0;    // --- Showers --- 
      FillWith(shwr_id,     -999);
      FillWith(shwr_dirx,   -999);
      FillWith(shwr_diry,   -999);
      FillWith(shwr_dirz,   -999);
      FillWith(shwr_startx, -999);
      FillWith(shwr_starty, -999);
      FillWith(shwr_startz, -999);
      FillWith(shwr_length, -999);
      FillWith(shwr_openangle, -999);
      nclusts               = 0;    // --- Hit Clusters ---
      FillWith(clust_tpc,  -9);
      FillWith(clust_plane, -9);
      FillWith(clust_wire, -9);
      FillWith(clust_nwires, -9);
      FillWith(clust_nhits, -9);
      FillWith(clust_charge, -999);
      FillWith(clust_time, -999);
      FillWith(clust_time_lh, -999);
      FillWith(clust_time_w, -999);
      FillWith(clust_startTime, -999);
      FillWith(clust_endTime, -999);
      FillWith(clust_g4id, -9);
      FillWith(clust_g4charge, -999);
      FillWith(clust_g4energy, -999);
      FillWith(clust_ismatch, 0);
      FillWith(clust_edepid, -9);
      FillWith(clust_blipid, -9);
      total_blip_energy     = -9;  // --- Blips ---
      nblips                = 0;
      FillWith(blip_tpc,  -9);
      FillWith(blip_nplanes,  -9);
      FillWith(blip_caloplane, -9);
      FillWith(blip_x, -99999);
      FillWith(blip_y, -99999);
      FillWith(blip_z, -99999);
      FillWith(blip_maxdiff, -9);
      FillWith(blip_charge, -999);
      FillWith(blip_energy, -999);
      FillWith(blip_energyESTAR, -999);
      FillWith(blip_edepid, -9);
      for(int i=0; i<kNplanes; i++)
        FillWith(blip_clustid[i], -9);

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
//      tree->Branch("gamma_depEnergy",&gamma_depEnergy,"gamma_depEnergy/F");
//      tree->Branch("gamma_numElectrons",&gamma_numElectrons,"gamma_numElectrons/F");
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
//      tree->Branch("edep_tpc",edep_tpc,"edep_tpc[nedeps]/I"); 
      tree->Branch("edep_g4id",edep_g4id,"edep_g4id[nedeps]/I"); 
      tree->Branch("edep_pdg",edep_pdg,"edep_pdg[nedeps]/I"); 
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
      tree->Branch("hit_sumadc",hit_sumadc,"hit_sumadc[nhits]/F"); 
      tree->Branch("hit_mult",hit_mult,"hit_mult[nhits]/I"); 
      tree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/F");
      tree->Branch("hit_ismatch",hit_ismatch,"hit_ismatch[nhits]/F");
      tree->Branch("hit_isreal",hit_isreal,"hit_isreal[nhits]/I"); 
      tree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I"); 
      tree->Branch("hit_shwrid",hit_shwrid,"hit_shwrid[nhits]/I"); 
      tree->Branch("hit_g4id",hit_g4id,"hit_g4id[nhits]/I");
      tree->Branch("hit_g4frac",hit_g4frac,"hit_g4frac[nhits]/F"); 
      tree->Branch("hit_g4energy",hit_g4energy,"hit_g4energy[nhits]/F"); 
      tree->Branch("hit_g4charge",hit_g4charge,"hit_g4charge[nhits]/F"); 
      tree->Branch("hit_clustid",hit_clustid,"hit_clustid[nhits]/I"); 
      tree->Branch("hit_blipid",hit_blipid,"hit_blipid[nhits]/I");
      tree->Branch("ntrks",&ntrks,"ntrks/I");          
      tree->Branch("trk_id",trk_id,"trk_id[ntrks]/I");       
      tree->Branch("trk_npts",trk_npts,"trk_npts[ntrks]/I");
      tree->Branch("trk_length",trk_length,"trk_length[ntrks]/F");
      tree->Branch("trk_startx",trk_startx,"trk_startx[ntrks]/F");
      tree->Branch("trk_starty",trk_starty,"trk_starty[ntrks]/F");
      tree->Branch("trk_startz",trk_startz,"trk_startz[ntrks]/F");
      tree->Branch("trk_startd",trk_startd,"trk_startd[ntrks]/F");
      tree->Branch("trk_endx",trk_endx,"trk_endx[ntrks]/F");
      tree->Branch("trk_endy",trk_endy,"trk_endy[ntrks]/F");
      tree->Branch("trk_endz",trk_endz,"trk_endz[ntrks]/F");
      tree->Branch("trk_endd",trk_endd,"trk_endd[ntrks]/F");
      //tree->Branch("trk_origin",trk_origin,"trk_origin[ntrks]/I");
      //tree->Branch("trk_g4id",trk_g4id,"trk_g4id[ntrks]/I");
      //tree->Branch("trk_g4pdg",trk_g4pdg,"trk_g4pdg[ntrks]/I");
      //tree->Branch("trk_purity",trk_purity,"trk_purity[ntrks]/F");
      //tree->Branch("trk_pitch",trk_pitch,"trk_pitch[ntrks]/F");
      //tree->Branch("trk_ke",trk_ke,"trk_ke[ntrks]/F");
      //tree->Branch("trk_cosmictag",trk_cosmictag,"trk_cosmictag[ntrks]/I"); 
      //tree->Branch("trk_cosmicscore",trk_cosmicscore,"trk_cosmicscore[ntrks]/F");
      //tree->Branch("trk_pidpdg",trk_pidpdg,"trk_pidpdg[ntrks]/I");
      //tree->Branch("trk_pidchi",trk_pidchi,"trk_pidchi[ntrks]/F");
      //tree->Branch("trk_bestplane",trk_bestplane,"trk_bestplane[ntrks]/I");
      tree->Branch("nshwrs",&nshwrs,"nshwrs/I");
      tree->Branch("shwr_id",shwr_id,"shwr_id[nshwrs]/I");
      tree->Branch("shwr_dirx",shwr_dirx,"shwr_dirx[nshwrs]/F");
      tree->Branch("shwr_diry",shwr_diry,"shwr_diry[nshwrs]/F");
      tree->Branch("shwr_dirz",shwr_dirz,"shwr_dirz[nshwrs]/F");
      tree->Branch("shwr_startx",shwr_startx,"shwr_startx[nshwrs]/F");
      tree->Branch("shwr_starty",shwr_starty,"shwr_starty[nshwrs]/F");
      tree->Branch("shwr_startz",shwr_startz,"shwr_startz[nshwrs]/F");
      tree->Branch("shwr_length",shwr_length,"shwr_length[nshwrs]/F");
      tree->Branch("shwr_openangle",shwr_openangle,"shwr_openangle[nshwrs]/F");
      tree->Branch("nclusts",&nclusts,"nclusts/I");
      tree->Branch("clust_tpc",clust_tpc,"clust_tpc[nclusts]/I");
      tree->Branch("clust_plane",clust_plane,"clust_plane[nclusts]/I");
      tree->Branch("clust_wire",clust_wire,"clust_wire[nclusts]/I");
      tree->Branch("clust_nwires",clust_nwires,"clust_nwires[nclusts]/I");
      tree->Branch("clust_nhits",clust_nhits,"clust_nhits[nclusts]/I");
      tree->Branch("clust_charge",clust_charge,"clust_charge[nclusts]/F");
      tree->Branch("clust_time",clust_time,"clust_time[nclusts]/F");
      tree->Branch("clust_time_w",clust_time_w,"clust_time_w[nclusts]/F");
      tree->Branch("clust_time_lh",clust_time_lh,"clust_time_lh[nclusts]/F");
      tree->Branch("clust_startTime",clust_startTime,"clust_startTime[nclusts]/F");
      tree->Branch("clust_endTime",clust_endTime,"clust_endTime[nclusts]/F");
      tree->Branch("clust_g4charge",clust_g4charge,"clust_g4charge[nclusts]/F");
      tree->Branch("clust_g4energy",clust_g4energy,"clust_g4energy[nclusts]/F");
      tree->Branch("clust_ismatch",clust_ismatch,"clust_ismatch[nclusts]/I");
      tree->Branch("clust_edepid",clust_edepid,"clust_edepid[nclusts]/I");
      tree->Branch("clust_blipid",clust_blipid,"clust_blipid[nclusts]/I");
      tree->Branch("total_blip_energy",&total_blip_energy,"total_blip_energy/F");
      tree->Branch("nblips",&nblips,"nblips/I");
      tree->Branch("blip_tpc",blip_tpc,"blip_tpc[nblips]/I");
      tree->Branch("blip_nplanes",blip_nplanes,"blip_nplanes[nblips]/I");
      tree->Branch("blip_caloplane",blip_caloplane,"blip_caloplane[nblips]/I");
      tree->Branch("blip_x",blip_x,"blip_x[nblips]/F");
      tree->Branch("blip_y",blip_y,"blip_y[nblips]/F");
      tree->Branch("blip_z",blip_z,"blip_z[nblips]/F");
      tree->Branch("blip_maxdiff",blip_maxdiff,"blip_maxdiff[nblips]/F");
      tree->Branch("blip_charge",blip_charge,"blip_charge[nblips]/F");
      tree->Branch("blip_energy",blip_energy,"blip_energy[nblips]/F");
      tree->Branch("blip_energyESTAR",blip_energyESTAR,"blip_energyESTAR[nblips]/F");
      tree->Branch("blip_edepid",blip_edepid,"blip_edepid[nblips]/I");
      for(int i=0; i<kNplanes; i++) 
        tree->Branch(Form("blip_clustid_pl%i",i),blip_clustid[i],Form("blip_clustid_pl%i[nblips]/I",i));

      
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
    
    void beginJob();                      // called once, at start of job
    void endJob();                        // called once, at end of job
    void analyze(const art::Event& evt);  // called per event

    private:
    void PrintParticleInfo(size_t);
    void PrintG4Info(const simb::MCParticle&);
  
    // --- Detector and clock data ---
    float TickOffset[kNplanes];

    // --- Data and calo objects ---
    BlipAnaTreeDataStruct*  fData;
    calo::CalorimetryAlg    fCaloAlg;
    TGraph2D*               ESTAR;

    // --- FCL configs ---
    std::string         fAnaTreeName;
    std::string         fHitModuleLabel;
    std::string         fLArG4ModuleLabel;
    std::string         fTrackModuleLabel;
    std::string         fShowerModuleLabel;
    std::string         fCaloModuleLabel;
    bool                fSaveParticleList;
    bool                fSaveDiagHistos;
    float               fBlipMergeDist;
    std::vector<float>  fMinHitRMS;
    std::vector<float>  fMaxHitRMS;
    std::vector<float>  fMinHitAmplitude;
    std::vector<int>    fPlanePriority;
    int                 fCaloPlane;         // use this plane for calorimetry
    int                 fMinMatchedPlanes;  // minimum number of planes needed to make a blip
    float               fDiffTolerance;     // 1cm in MicroBooNE
    float               fHitMatchTolerance;
  
    // --- Counters and such ---
    int   fNumEvents          = 0;
    int   fNumHits_pl[3][2]   = {};

    // --- Histograms ---
    TH1D*   h_nhits[kNplanes][2];       // --- hit diagnostics ---  
    TH1D*   h_hitph[kNplanes][2];
    TH1D*   h_hitrms[kNplanes][2];
    TH2D*   h_hitrms_vs_ph[kNplanes][2];
    TH1D*   h_sumadc_res[kNplanes][2];
    TH2D*   h_hitarea_vs_sumadc[kNplanes][2];
    TH2D*   h_nelec_TrueVsReco[kNplanes];
    TH1D*   h_nelec_Resolution[kNplanes];
    TH1D*   h_chargecomp[kNplanes];
    TH1D*   h_hitpur[kNplanes];
    TH1D*   h_hit_dT;

    // Initialize histograms
    void InitializeHistograms(){
      art::ServiceHandle<art::TFileService> tfs;
     
      h_hit_dT = tfs->make<TH1D>("hit_dT","Hit dT [ticks]",500,0,50);

      for(int i=0; i<kNplanes; i++) {
        h_chargecomp[i] = tfs->make<TH1D>(Form("pl%i_charge_completeness",i),Form("Charge completness, plane %i",i),101,0,1.01);
        h_hitpur[i]     = tfs->make<TH1D>(Form("pl%i_hit_purity",i),Form("Hit purity, plane %i",i),101,0,1.01);
      }

      if( fSaveDiagHistos ) {
        art::TFileDirectory diagDir = tfs->mkdir("Diagnostics");
        float hitMax  = 300;    int hitBins   = 300;
        float phMax = 100;      int phBins    = 200;
        float rmsMax = 10;      int rmsBins   = 200;
        float areaMax = 500;    int areaBins  = 100;
        for(int i=0; i<kNplanes; i++) {
          std::vector<std::string> label = { "fake", "real" };
          for(int j=0;j<2;j++) {
            h_nhits[i][j]= diagDir.make<TH1D>(Form("pl%i_%s_nhits",i,label[j].c_str()),  Form("Plane %i, %s hits;Number of hits",i,label[j].c_str()),hitBins,0,hitMax);
            h_hitph[i][j]= diagDir.make<TH1D>(Form("pl%i_%s_hit_ph",i,label[j].c_str()), Form("Plane %i, %s hits;Pulse height [ADC]",i,label[j].c_str()),phBins,0,phMax);
            h_hitrms[i][j]= diagDir.make<TH1D>(Form("pl%i_%s_hit_rms",i,label[j].c_str()),Form("Plane %i, %s hits;Hit RMS [ticks]",i,label[j].c_str()),rmsBins,0,rmsMax);
            h_hitrms_vs_ph[i][j]  = diagDir.make<TH2D>(Form("pl%i_%s_hit_rms_vs_ph",i,label[j].c_str()),
              Form("Plane %i, %s hits;Hit RMS [ticks];Pulse height [ADC]",i,label[j].c_str()),
              rmsBins/2,0,rmsMax,
              phBins/2,0,phMax);
            h_hitarea_vs_sumadc[i][j]  = diagDir.make<TH2D>(Form("pl%i_%s_hit_area_vs_sumadc",i,label[j].c_str()),
              Form("Plane %i, %s hits;Hit area [ticks*ADC];Summed ADC [ticks*ADC]",i,label[j].c_str()),
              areaBins,0,areaMax,
              areaBins,0,areaMax);
            h_sumadc_res[i][j] = diagDir.make<TH1D>(Form("pl%i_%s_sumadc_res",i,label[j].c_str()),
              Form("Plane %i, %s hits;(SummedADC - Integral) / Integral",i,label[j].c_str()),
              300,-3,3);
          }
            h_nelec_TrueVsReco[i] = diagDir.make<TH2D>( Form("pl%i_nelec_TrueVsReco",i),
              Form("Plane %i;True hit charge [ #times 10^{3} electrons ];Reconstructed hit charge [ #times 10^{3} electrons ]",i),60,0,30, 60,0,30);
            h_nelec_TrueVsReco[i] ->SetOption("colz");
            h_nelec_Resolution[i] = diagDir.make<TH1D>( Form("pl%i_nelec_res",i),Form("Plane %i;Hit charge resolution: (reco-true)/true",i),200,-2,2);
        }
      }
    
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
  ,fAnaTreeName       (pset.get< std::string >            ("AnaTreeName",     "anatree"))
  ,fHitModuleLabel    (pset.get< std::string >            ("HitModuleLabel",  "gaushit"))
  ,fLArG4ModuleLabel  (pset.get< std::string >            ("LArG4ModuleLabel","largeant"))
  ,fTrackModuleLabel  (pset.get< std::string >            ("TrackModuleLabel","pandoraTrack"))
  ,fShowerModuleLabel  (pset.get< std::string >           ("ShowerModuleLabel","pandoraShower"))
  ,fCaloModuleLabel  (pset.get< std::string >             ("CaloModuleLabel", "pandoraCalo"))
  ,fSaveParticleList  (pset.get< bool >                   ("SaveParticleList", true))
  ,fSaveDiagHistos    (pset.get< bool >                   ("SaveDiagnosticHistos",true))
  ,fBlipMergeDist     (pset.get< float >                  ("BlipMergeDist",     0.3))
  ,fMinHitRMS         (pset.get< std::vector< float > >   ("MinHitRMS",         {1.5, 1.5,  1.5}))
  ,fMaxHitRMS         (pset.get< std::vector< float > >   ("MaxHitRMS",         {5,   5,    5}))
  ,fMinHitAmplitude   (pset.get< std::vector< float > >   ("MinHitAmplitude",   {10, 14,    10}))
  ,fPlanePriority     (pset.get< std::vector< int > >     ("PlanePriority",     {2,   0,    1}))
  ,fCaloPlane         (pset.get< int >                    ("CaloPlane",         2))
  ,fMinMatchedPlanes  (pset.get< int >                    ("MinMatchedPlanes",  2))
  ,fDiffTolerance     (pset.get< float >                  ("DiffTolerance",     999.0))
  ,fHitMatchTolerance    (pset.get< float >                  ("HitMatchTolerance",    2.))
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

  // Load the ESTAR lookup table
  std::string fname;
  //std::string path="ShowerEnergyReco/ESTAREnergyLookupCurve.root";
  std::string path="ESTAREnergyLookupCurve.root";
  std::string gname="ESTAR_energy_lookup_curve";
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(path,fname);
  if( fname.empty() ) {
    throw cet::exception("BlipAna") << "Could not find ESTAR lookup curve.\n";
  } else {
    TFile fin(fname.c_str(),"READ");
    if(!fin.IsOpen()) {
      throw cet::exception("BlipAna") << "Could not open ESTAR file.\n";
    } else {
      ESTAR = dynamic_cast<TGraph2D*>(fin.Get(gname.c_str()));
      if(!ESTAR) throw cet::exception("BlipAna") << "Could not read the ESTAR TGraph";
    }
  }

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
  fNumEvents++;


  //=========================================
  // Get data products for this event
  //=========================================
  
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;
  
  // -- G4 particles
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();

  // -- hits
  art::Handle< std::vector<recob::Hit> > hitlistHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitlistHandle))
    art::fill_ptr_vector(hitlist, hitlistHandle);
  
  // -- tracks
  art::Handle< std::vector<recob::Track> > tracklistHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,tracklistHandle))
    art::fill_ptr_vector(tracklist, tracklistHandle);
  
  // -- showers
  art::Handle< std::vector<recob::Shower> > showerlistHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerlistHandle))
    art::fill_ptr_vector(showerlist, showerlistHandle);

  // -- hit<->track associations
  art::FindManyP<recob::Track> fmtrk(hitlistHandle,evt,fTrackModuleLabel);

  // -- hit<->shower associations
  art::FindManyP<recob::Shower> fmshwr(hitlistHandle,evt,fShowerModuleLabel);
  
  // -- track<->calorimetry
  art::FindMany<anab::Calorimetry> fmcal(tracklistHandle, evt, fCaloModuleLabel);

  // Resize data struct objects
  fData->nhits      = (int)hitlist.size();
  fData->nparticles = (int)plist.size();
  fData->ntrks      = (int)tracklist.size();
  fData->nshwrs     = (int)showerlist.size();
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

    float numElectrons  = 0, edep = 0;
    BlipUtils::CalcPartDep(trackID,edep,numElectrons);

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
      fData->numElectrons[i]    = numElectrons;
      //if( pathlen ) PrintParticleInfo(i);
      PrintParticleInfo(i);
    }
  
  }//endloop over G4 particles

  // Merge and save true blip information
  MergeTrueBlips(trueblips, fBlipMergeDist); 
  fData->nedeps = (int)trueblips.size();
  std::cout<<"Found "<<trueblips.size()<<" true blips:\n";
  for(size_t i=0; i<trueblips.size(); i++ ) {
    trueblips[i].ID       = i;
    fData->edep_tpc[i]    = trueblips.at(i).TPC;
    fData->edep_energy[i] = trueblips.at(i).Energy;
    fData->edep_charge[i] = trueblips.at(i).NumElectrons;
    fData->edep_ds[i]     = trueblips.at(i).Length;
    fData->edep_x[i]      = trueblips.at(i).Position.X();
    fData->edep_y[i]      = trueblips.at(i).Position.Y();
    fData->edep_z[i]      = trueblips.at(i).Position.Z();
    fData->edep_g4id[i]   = trueblips.at(i).LeadG4ID;
    fData->edep_pdg[i]    = trueblips.at(i).LeadG4PDG;
    std::cout
    <<"   ~ "<<trueblips.at(i).Energy<<" MeV, "
    <<" ds= "<<trueblips.at(i).Length<<" cm, "
    <<" trkID= "<<trueblips.at(i).LeadG4ID<<", pdg "<<trueblips.at(i).LeadG4PDG<<"\n";
  }


  //====================================
  // Save hit information
  //====================================
  std::vector<BlipUtils::HitInfo> hitinfo(hitlist.size());
  std::map<int,std::vector<int>> wirehitsMap;
  
  for(size_t i=0; i<hitlist.size(); i++){
    int   wire  = hitlist[i]->WireID().Wire;
    int   plane = hitlist[i]->WireID().Plane;
    int   tpc   = hitlist[i]->WireID().TPC;
    hitinfo[i].hitid      = i;
    hitinfo[i].hit        = hitlist[i];
    hitinfo[i].wire       = wire;
    hitinfo[i].tpc        = tpc;
    hitinfo[i].plane      = plane;
    hitinfo[i].driftTicks = hitlist[i]->PeakTime()-TickOffset[plane];
    hitinfo[i].qcoll      = fCaloAlg.ElectronsFromADCArea(hitlist[i]->Integral(),plane);
    
    // Find G4 particle ID for leading contributor
    if( BlipUtils::DoesHitHaveSimChannel(hitlist[i]) ){
      BlipUtils::HitTruth( hitlist[i], hitinfo[i].g4id, hitinfo[i].g4frac, hitinfo[i].g4energy, hitinfo[i].g4charge);
      hitinfo[i].g4ids = BlipUtils::HitTruthIds(hitlist[i]);
      hitinfo[i].isreal = (hitinfo[i].g4id > 0);
      //if( isreal ) {                    // group in overlapping hits
      //  for(size_t j=0; j<hitlist.size(); j++){
      //    if( BlipUtils::DoHitsOverlap(hitlist[i],hitlist[j]) ) hitinfo[j].isreal = isreal;
      //  }
      //}
    }
   
   // Find associated track
    if (fmtrk.isValid()){ 
      if (fmtrk.at(i).size())  hitinfo[i].trkid = fmtrk.at(i)[0]->ID();
    }
    
    // Find associated shower
    if (fmshwr.isValid()){
      if (fmshwr.at(i).size()) hitinfo[i].shwrid = fmshwr.at(i)[0]->ID();
    }
    
    // Add to the wire-by-wire hit map
    wirehitsMap[wire].push_back(i);

    // Save TTree data
    if( i < kMaxHits ) {
      fData->hit_plane[i]   = plane;
      fData->hit_wire[i]    = wire;
      fData->hit_tpc[i]     = tpc;
      fData->hit_trkid[i]   = hitinfo[i].trkid;
      fData->hit_shwrid[i]  = hitinfo[i].shwrid;
      fData->hit_channel[i] = hitlist[i]->Channel();
      fData->hit_peakT[i]   = hitlist[i]->PeakTime();
      fData->hit_rms[i]     = hitlist[i]->RMS();
      fData->hit_ph[i]	    = hitlist[i]->PeakAmplitude();
      fData->hit_area[i]    = hitlist[i]->Integral();
      fData->hit_sumadc[i]  = hitlist[i]->SummedADC();
      fData->hit_mult[i]    = hitlist[i]->Multiplicity();
      fData->hit_time[i]    = hitinfo[i].driftTicks;
      fData->hit_charge[i]  = hitinfo[i].qcoll;
      fData->hit_isreal[i]  = hitinfo[i].isreal;
      fData->hit_g4id[i]    = hitinfo[i].g4id;
      fData->hit_g4frac[i]  = hitinfo[i].g4frac;
      fData->hit_g4energy[i]= hitinfo[i].g4energy;
      fData->hit_g4charge[i]= hitinfo[i].g4charge;
    }
  }
 
  // Time matching
  for(size_t i=0; i<hitlist.size(); i++){
    for(size_t j=i+1; j<hitlist.size(); j++){
      if( hitinfo[i].plane == hitinfo[j].plane ) continue;
      if( hitinfo[i].tpc != hitinfo[j].tpc   ) continue;
      if( hitinfo[i].trkid >= 0 || hitinfo[j].trkid >= 0) continue;
      float dT = fabs(hitinfo[j].driftTicks - hitinfo[i].driftTicks);
      h_hit_dT->Fill( dT );
      if( dT <= fHitMatchTolerance ) {
        hitinfo[i].ismatch = true;
        hitinfo[j].ismatch = true;
        fData->hit_ismatch[i] = true;
        fData->hit_ismatch[j] = true;
      }
    }
  }
 
  // Flag real hits and fill hit diagnostic histograms
  int   nhits_pl[3][2]={};
  float total_chargeInHits[3] = {0};
  for(size_t i=0; i<hitlist.size(); i++){
    int   pl      = hitlist[i]->WireID().Plane; 
    int   ir      = hitinfo[i].isreal;
    float integral= hitlist[i]->Integral();
    float sumadc  = hitlist[i]->SummedADC();
    float q       = hitinfo[i].qcoll;
    float qtrue   = hitinfo[i].g4charge;
    nhits_pl[pl][ir]++;
    fNumHits_pl[pl][ir]++;
    if(ir) total_chargeInHits[pl] += qtrue;
    if( fSaveDiagHistos ) {
      h_hitph[pl][ir]         ->Fill(hitlist[i]->PeakAmplitude());
      h_hitrms[pl][ir]        ->Fill(hitlist[i]->RMS());
      h_hitrms_vs_ph[pl][ir]  ->Fill(hitlist[i]->RMS(),hitlist[i]->PeakAmplitude());
      h_hitarea_vs_sumadc[pl][ir]->Fill(integral, sumadc);
      if(integral>0) h_sumadc_res[pl][ir]    ->Fill( (sumadc-integral)/integral );
      if(ir){
        h_nelec_TrueVsReco[pl]->Fill(q/1e3,qtrue/1e3);
        h_nelec_Resolution[pl]->Fill((q-qtrue)/qtrue);
      }
    }
  }

  // Fill histograms
  for(size_t i=0; i<kNplanes; i++){
    int   totHits = nhits_pl[i][0]+nhits_pl[i][1];
    std::cout<<"charge in hits: "<<total_chargeInHits[i]<<", total charge: "<<fData->total_numElectrons<<"\n";
    float qcomp = -9, pur = -9;
    if( fData->total_numElectrons ) qcomp = total_chargeInHits[i]/fData->total_numElectrons;
    if( totHits ) pur   = nhits_pl[i][1] / float(totHits);
    h_chargecomp[i] ->Fill( qcomp );
    h_hitpur[i]     ->Fill( pur );
    std::cout<<"Hits on plane "<<i<<": "<<nhits_pl[i][1]<<" real, "<<nhits_pl[i][0]<<" fake\n";
    std::cout<<" - charge completeness: "<<qcomp<<"\n";
    std::cout<<" - hit purity: "<<pur<<"\n";
    if( fSaveDiagHistos ) for(size_t j=0; j<2; j++) h_nhits[i][j]->Fill(nhits_pl[i][j]);
  }
  
  
  
  //====================================
  // Save track information
  //====================================
  for(size_t i=0; i<tracklist.size(); i++){
    const auto& startPt = tracklist[i]->Vertex();
    const auto& endPt   = tracklist[i]->End();
    fData->trk_id[i]    = tracklist[i]->ID();
    fData->trk_npts[i]  = tracklist[i]->NumberTrajectoryPoints();
    fData->trk_length[i]= tracklist[i]->Length();
    fData->trk_startx[i]= startPt.X();
    fData->trk_starty[i]= startPt.Y();
    fData->trk_startz[i]= startPt.Z();
    fData->trk_endx[i]  = endPt.X();
    fData->trk_endy[i]  = endPt.Y();
    fData->trk_endz[i]  = endPt.Z();
    fData->trk_startd[i]= BlipUtils::DistToBoundary(startPt);
    fData->trk_endd[i]  = BlipUtils::DistToBoundary(endPt);
  }
  
  //====================================
  // Save shower information
  //====================================
  for(size_t i=0; i<showerlist.size(); i++){
    fData->shwr_id[i]     = showerlist[i]->ID();
    fData->shwr_startx[i] = showerlist[i]->ShowerStart()[0]; 
    fData->shwr_starty[i] = showerlist[i]->ShowerStart()[1]; 
    fData->shwr_startz[i] = showerlist[i]->ShowerStart()[2]; 
    fData->shwr_dirx[i]   = showerlist[i]->Direction()[0];
    fData->shwr_diry[i]   = showerlist[i]->Direction()[1];
    fData->shwr_dirz[i]   = showerlist[i]->Direction()[2];
    fData->shwr_length[i] = showerlist[i]->Length();
    fData->shwr_openangle[i] = showerlist[i]->OpenAngle();
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
  //  [x] Merge together clusters on adjacent wires (if they match up in time)
  //  [x] Plane-to-plane time matching
  //  [x] Wire intersection check to get XYZ
  //  [x] Create "blip" object and save to tree (nblips, blip_xyz, blip_charge, blip_g4energy)

  // Create a series of masks that we'll update as we go along
  std::vector<bool> hitPassesCuts(hitlist.size(), false);
  std::vector<bool> hitIsClustered(hitlist.size(),false);
  
  for(size_t i=0; i<hitlist.size(); i++){
    int plane = hitlist[i]->WireID().Plane;
    
    // Exclude any hits that were part of tracks
    if( hitinfo[i].trkid >= 0 ) continue;

    // Check hit widths and amplitudes
    if( hitlist[i]->RMS() > fMinHitRMS[plane] &&
        hitlist[i]->RMS() < fMaxHitRMS[plane] &&
        hitlist[i]->PeakAmplitude() > fMinHitAmplitude[plane] ) {
      hitPassesCuts[i] = true;
      // Since this hit looks good, also flag any overlapping hits
      // on the same wire as good too
      for(auto j : wirehitsMap[hitlist[i]->WireID().Wire] ) {
        if( BlipUtils::DoHitsOverlap(hitlist[i],hitlist[j]) ) 
          hitPassesCuts[j] = true;
      }
    }
    
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
          float t1 = hitinfo[hj].driftTicks - rms;
          float t2 = hitinfo[hj].driftTicks + rms;
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
  std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap_merged;
  for(auto const& tpcMap : tpc_planeclustsMap ) {
    for(auto const& planeclusts : tpcMap.second ) {
      for( auto const& ci : planeclusts.second ){
        if( hitclust[ci].isMerged ) continue;
//        std::cout<<"   looking for adjacent clusters to "<<ci<<", wire "<<hitclust[ci].StartWire<<" "<<hitclust[ci].EndWire<<"\n";
        BlipUtils::HitClust hc = hitclust[ci];
        hitclust[ci].isMerged = true;
        int clustsAdded;
        do{
          clustsAdded = 0;
          for(auto const& cj : planeclusts.second ){
            if( hitclust[cj].isMerged ) continue;
            int w1 = hitclust[cj].StartWire;
            int w2 = hitclust[cj].EndWire;
  //          std::cout<<"     comparing to clust "<<cj<<", which has w1: "<<w1<<" and w2: "<<w2<<", t: "<<hitclust[cj].Time<<"\n";
            // check if clusters are adjacent
            int dw1 = (w1-hc.EndWire);
            int dw2 = (hc.StartWire-w2);
            if( (dw1 >= 0 && dw1 < 2) || 
                (dw2 >= 0 && dw2 < 2) ){
              // check if there's a time match
    //          std::cout<<"          potential match!\n";
              if( BlipUtils::DoHitClustsOverlap(hc,hitclust[cj]) ) {
      //          std::cout<<"       confirmed match!\n";
                hc = BlipUtils::MergeHitClusts(hc,hitclust[cj]);
                if( hitclust[cj].isMerged ) { clustsAdded++;}
              }
            }
          }
        } while ( clustsAdded!=0);
        hitclust_merged.push_back(hc);
        tpc_planeclustsMap_merged[hc.TPC][hc.Plane].push_back(hitclust_merged.size()-1);
      }
    }
  }
  hitclust = hitclust_merged;
  tpc_planeclustsMap = tpc_planeclustsMap_merged;

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
    fData->clust_time_w[i]  = hitclust[i].WeightedTime;
    fData->clust_time_lh[i] = hitclust[i].LeadHitTime;
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
      if( tbG4 >= 0 && tbG4 == hitclust[i].LeadHitG4ID ) {
        fData->edep_clustid[j] = hitclust[i].ID;
        fData->clust_edepid[i] = trueblips[j].ID;
        fData->clust_g4energy[i] = trueblips[j].Energy; 
        fData->clust_g4charge[i] = trueblips[j].NumElectrons;
        break;
      }
    }
  }
  
  std::cout<<"After merging: "<<hitclust_merged.size()<<" hit clusts:\n";
  for(size_t i=0; i<hitclust.size(); i++){
    std::cout<<"   * "<<i<<", wire: "<<hitclust[i].LeadHitWire<<", start/end wire: "<<hitclust[i].StartWire<<"/"<<hitclust[i].EndWire<<"  (span: "<<hitclust[i].Wires.size()<<"), plane: "<<hitclust[i].Plane<<", time: "<<hitclust[i].Time<<", half-width="<<(hitclust[i].EndTime-hitclust[i].StartTime)/2.<<", leadG4: "<<hitclust[i].LeadHitG4ID<<", edepid: "<<fData->clust_edepid[i]<<", E: "<<fData->clust_g4energy[i]<<"\n";
  }

  // =============================================================================
  // Plane matching and 3D blip formation
  // =============================================================================
  std::vector<BlipUtils::Blip> blips;

  // --------------------------------------
  // Method 1:  Require match between the calo plane (typically collection)
  //            and one or two induction planes.
  for(auto const& tpcMap : tpc_planeclustsMap ) { // loop on TPCs
    std::cout<<"Performing cluster matching in TPC "<<tpcMap.first<<", which has "<<tpcMap.second.size()<<" planes\n";
    auto planeMap = tpcMap.second;
    if( planeMap.find(fCaloPlane) != planeMap.end() ){
      std::cout<<"  this TPC has clusters on our calo plane ("<<fCaloPlane<<")\n";
      int   planeA            = fCaloPlane;
      auto  hitclusts_planeA  = planeMap[planeA];
      for(auto const& i : hitclusts_planeA ) {
//        float t1 = hitclust[i].StartTime;
//        float t2 = hitclust[i].EndTime;
        std::vector<BlipUtils::HitClust> hcgroup;
        bool isPlaneMatched[3] = {false,false,false};
        for(auto  hitclusts_planeB : planeMap ) {
          int planeB = hitclusts_planeB.first;
          if( planeB == planeA ) continue;
          if( isPlaneMatched[planeB] ) continue;
          for(auto const& j : hitclusts_planeB.second ) {
//            if( BlipUtils::DoHitClustsOverlap(hitclust[j],t1,t2) ) {
            if( BlipUtils::DoHitClustsMatch(hitclust[j],hitclust[i],fHitMatchTolerance)) {
              if( !hitclust[i].isMatched ) {
                hitclust[i].isMatched = true;
                isPlaneMatched[planeA] = true;
                hcgroup.push_back(hitclust[i]);
              }
              if( !hitclust[j].isMatched ) {
                hitclust[j].isMatched = true;
                isPlaneMatched[planeB] = true;
                hcgroup.push_back(hitclust[j]);
              }
  //            t1 = std::min(hitclust[i].StartTime,hitclust[j].StartTime);
  //            t2 = std::max(hitclust[i].EndTime,hitclust[j].EndTime);
            }//endif match found
          }//endloop over clusters on this plane
        }//endloop over planes
        if(hcgroup.size()) {
          //std::cout<<"Attempting to make blip from "<<hcgroup.size()<<" matches...\n";
          BlipUtils::Blip newBlip = BlipUtils::MakeBlip(hcgroup);
          if( newBlip.isValid && newBlip.MaxIntersectDiff <= fDiffTolerance ){
            //std::cout<<"New blip created with max intersection difference of "<<newBlip.MaxIntersectDiff<<"\n";
            blips.push_back(newBlip);
            for(auto hc : hcgroup ) hitclust[hc.ID].BlipID = blips.size()-1;
          }
        }
      }//endloop over clusts on planeA
    }//endif planeA has clusts
  }//endloop over TPCs
  
  /*
  // =============================================================================
  // Now for the fun stuff! Match hit clusters in time between the different planes 
  // in order to make Blips. We will require only a match between at least 2 planes
  // for now (will make this fhicl-configurable eventually).
  // ---------------------------------------------------
  // Create collection of Blips
  // ---------------------------------------------------
  //std::cout<<"Looping over clusters to make blips...\n";
  for(size_t i=0; i<hitclust.size(); i++){
    int iPlane = hitclust[i].Plane;
    int iTPC = hitclust[i].TPC;
    if( hitclust[i].isMatched ) continue;
    
    float t1 = hitclust[i].StartTime;
    float t2 = hitclust[i].EndTime;
    std::vector<BlipUtils::HitClust> hcgroup;
    bool isPlaneMatched[3] = {false,false,false};
    //std::cout<<"    cl "<<i<<" on plane "<<iPlane<<"   t1-t2: "<<t1<<"   "<<t2<<"\n";
    for(size_t j=0; j<hitclust.size(); j++){
      int jPlane = hitclust[j].Plane;
      int jTPC = hitclust[j].TPC;
      if( hitclust[j].isMatched || isPlaneMatched[jPlane] ) continue;
      if( i==j || iPlane==jPlane || jTPC!=iTPC )        continue;
      //std::cout<<"     does cl "<<j<<" on plane "<<jPlane<<" match up? "<<hitclust[j].StartTime<<"  "<<hitclust[j].EndTime<<"\n";
      
      if( BlipUtils::DoHitClustsMatch(hitclust[j],t1,t2)) {
        if( !hitclust[i].isMatched ) {
          hitclust[i].isMatched = true;
          isPlaneMatched[iPlane] = true;
          hcgroup.push_back(hitclust[i]);
        }
        if( !hitclust[j].isMatched ) {
          hitclust[j].isMatched = true;
          isPlaneMatched[jPlane] = true;
          hcgroup.push_back(hitclust[j]);
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
  */


  // -----------------------------------------------------
  // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
  // like associating blip with some nearby track/shower and using its tagged T0)
  //    Method 1: Assume a dE/dx = 2 MeV/cm for electrons, use that + local E-field
  //              calculate recombination.
  //    Method 2: ESTAR lookup table method ala ArgoNeuT
  for(size_t i=0; i<blips.size(); i++){
    
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    float qColl = blips[i].Charge[fCaloPlane];
    float td    = blips[i].DriftTime;
    float depEl = qColl * exp( td / detProp.ElectronLifetime() ); 
    auto const blipPos = blips[i].Position;
    float Efield = detProp.Efield(0);
    if( SCE->EnableSimEfieldSCE() ) {
      geo::Point_t point = {double(blipPos.X()), double(blipPos.Y()), double(blipPos.Z())};
      auto const EfieldOffsets = SCE->GetEfieldOffsets(point);
      Efield *= std::hypot(1+EfieldOffsets.X(), EfieldOffsets.Y(), EfieldOffsets.Z());
    }
    
    // METHOD 1
    float dEdx = 2.0; // MeV/cm
    float recomb = BlipUtils::ModBoxRecomb(dEdx,Efield);
    blips[i].Energy = depEl * (1./recomb) * 23.6e-6;

    // METHOD 2
    blips[i].EnergyESTAR = ESTAR->Interpolate(depEl, Efield); 

  }
  
  fData->total_blip_energy = 0;
 

  // Save blip info to tree
  fData->nblips = blips.size();
  for(size_t i=0; i<blips.size(); i++){
    fData->blip_tpc[i]        = blips[i].TPC;
    fData->blip_nplanes[i]    = blips[i].NPlanes;
    fData->blip_caloplane[i]  = fCaloPlane;
    fData->blip_charge[i]     = blips[i].Charge[fCaloPlane];
    fData->blip_energy[i]     = blips[i].Energy;
    fData->blip_energyESTAR[i]= blips[i].EnergyESTAR;
    fData->blip_maxdiff[i]    = blips[i].MaxIntersectDiff;
    fData->blip_x[i]          = (float)blips[i].Position.X();
    fData->blip_y[i]          = (float)blips[i].Position.Y();
    fData->blip_z[i]          = (float)blips[i].Position.Z();

    // std::cout<<"Finding true energy dep for blip "<<i<<" (associated clusters: "<<blips[i].ClustIDs.size()<<")\n";
    // find associated true edep
    float max = 0;
    for(auto hitID : blips[i].HitIDs ) {
      if( hitID >= 0 ) {
        fData->hit_blipid[hitID] = i;
      }
    }
    for(auto clustID : blips[i].ClustIDs ) {
      int clustPlane  = fData->clust_plane[clustID];
      int edepid      = fData->clust_edepid[clustID];
      fData->blip_clustid[clustPlane][i] = clustID;
      if( edepid >= 0 && edepid < fData->nedeps ) {
        float E = trueblips[edepid].Energy;
        if( E > max ) {
          fData->blip_edepid[i] = edepid;
          fData->edep_blipid[edepid] = i;
          max = E;
        }
      }
    }
//    std::cout<<"  best-matched edepid for blip "<<i<<": "<<fData->blip_edepid[i]<<"\n";


  }
  
  std::cout<<"Reconstructed "<<blips.size()<<" blips:\n";
  for(size_t i=0; i<blips.size(); i++){
    auto const& b = blips[i];
    fData->total_blip_energy += b.Energy;
    std::cout
    <<"   -- "<<i<<", TPC: "<<b.TPC
    <<"; charge: "<<b.Charge[2]
    <<"; recoEnergy: "<<b.Energy<<" MeV"
    <<"; Position: "<<b.Position.X()<<", "<<b.Position.Y()<<", "<<b.Position.Z()
    <<"; MaxIntersectDiff: "<<b.MaxIntersectDiff
    <<"; EdepID: "<<fData->blip_edepid[i]
    <<"\n";
  }

  // Update clust data in Tree, so we can tie reconstructed hit clusters
  // to the blips they were grouped into
  for(size_t i=0; i<hitclust.size(); i++){
    fData->clust_ismatch[i] = hitclust[i].isMatched;
    fData->clust_blipid[i] = hitclust[i].BlipID;
  }


  //====================================
  // Fill TTree
  //====================================
  fData->FillTree();
}


//###################################################
//  endJob: output useful info to screen
//###################################################
void sbnd::BlipAna::endJob(){
  
  printf("\n=============================================\n");
  printf("BlipAna Statistics\n\n");
  printf("  Total events  :   %i\n",fNumEvents);
  for(size_t i=0; i<kNplanes; i++){
    printf("\n  Plane %lu \n",i);
    printf("   - ave total fake hits  : %.2f\n",fNumHits_pl[i][0]/(float)fNumEvents);
    printf("   - ave total real hits  : %.2f\n",fNumHits_pl[i][1]/(float)fNumEvents);
    printf("   - charge completeness  : %f\n",h_chargecomp[i]->GetMean());
    printf("   - hit purity           : %f\n",h_hitpur[i]->GetMean());
  } 
  printf("\n=============================================\n");

}




//###################################################
//  Printouts for debugging
//###################################################

void sbnd::BlipAna::PrintParticleInfo(size_t i){
  printf("  %5i  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, dL=%7.1f, KE0=%8.2f, Edep=%8.3f, T=%10.2f, moth=%5i, %12s, ND=%i\n",
   (int)i,
   fData->trackID[i],
   fData->pdg[i],
   fData->startPointx[i],
   fData->startPointy[i],
   fData->startPointz[i],
   fData->pathlen[i], 
   fData->E[i]-fData->mass[i],
   fData->depEnergy[i],
   fData->startT[i]/1e3,
   fData->mother[i],
   fData->process[i].c_str(),
   fData->nDaughters[i]
  ); 
}

void sbnd::BlipAna::PrintG4Info(const simb::MCParticle& part){
  printf("  trkID: %-6i PDG: %-8i XYZ= %7.1f %7.1f %7.1f, endXYZ: %7.1f %7.1f %7.1f, KE0=%8.1f, moth=%5i, %12s, ND=%i\n",
   (int)part.TrackId(),
   (int)part.PdgCode(),
   (float)part.Vx(),
   (float)part.Vy(),
   (float)part.Vz(),
   (float)(part.EndPosition()[0]),
   (float)(part.EndPosition()[1]),
   (float)(part.EndPosition()[2]),
   float(1e3*(part.E()-part.Mass())),
   (int)part.Mother(),
   part.Process().c_str(),
   (int)part.NumberDaughters()
  ); 
}


DEFINE_ART_MODULE(sbnd::BlipAna)

#endif
