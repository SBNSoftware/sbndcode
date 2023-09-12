/////////////////////////////////////////////////////////////////////
// 
//  BlipRecoAlg
//
//  To include in your module, add the header file at the top:
//    #include "sbndcode/BlipReco/Alg/BlipRecoAlg.h
//
//  then declare a private alg object in your class definition:
//    BlipRecoAlg fBlipAlg;
//    
//  W. Foreman, May 2022
//  wforeman @ iit.edu
//
/////////////////////////////////////////////////////////////////////
#ifndef BLIPRECOALG_H
#define BLIPRECOALG_H

// framework includes
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "art/Framework/Principal/Event.h"


// Microboone includes
//#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
//#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
//#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
//#include "ubevt/Database/UbooneElectronLifetimeService.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

// Blip-specific utils
#include "sbndcode/BlipRecoSBND/Utils/BlipUtils.h"

// ROOT stuff
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

// c++
#include <vector>
#include <iostream>
#include <memory>
#include <math.h>
#include <limits>


namespace blip {

  //--------------------------------------------
  class BlipRecoAlg{
   public:
    
    //Constructor/destructor
    BlipRecoAlg( fhicl::ParameterSet const& pset );
    BlipRecoAlg();
    ~BlipRecoAlg();
  
    void    reconfigure(fhicl::ParameterSet const& pset );
    void    RunBlipReco(const art::Event& evt);
    void    PrintConfig();
    
    // TO-DO: make these private and create getters instead
    std::vector<blip::HitInfo>      hitinfo;
    std::vector<blip::HitClust>     hitclust;
    std::vector<blip::Blip>         blips;  
    std::vector<blip::TrueBlip>     trueblips;
    std::vector<blip::ParticleInfo> pinfo;
    
    calo::CalorimetryAlg*   fCaloAlg;
    float   ModBoxRecomb(float,float);
    float   dQdx_to_dEdx(float,float);
    float   Q_to_E(float,float);

    float   fNominalRecombFactor;
   
    std::vector<bool>   fBadChanMask;
    std::vector<bool>   fBadChanMaskPerEvt;
    int                 EvtBadChanCount;
    
    TH1D*   h_recoWireEff_denom;
    TH1D*   h_recoWireEff_num;
    
    TH1D*   h_recoWireEffQ_denom;
    TH1D*   h_recoWireEffQ_num;


   private:
    
    //const detinfo::DetectorProperties& detProp;
    float kLArDensity;
    float kNominalEfield;
    float kDriftVelocity;
    
    float               mWion;

    // --- FCL configs ---
    std::string         fHitProducer;
    std::string         fTrkProducer;
    std::string         fGeantProducer;
    std::string         fSimDepProducer;
    std::string         fSimChanProducer;
    float               fSimGainFactor;
    bool                fDebugMode;
    float               fTrueBlipMergeDist;
    bool                fDoHitFiltering;
    float               fMaxHitTrkLength;
    float               fMaxHitAmp;
    std::vector<float>  fMinHitAmp;
    std::vector<float>  fMinHitRMS;
    std::vector<float>  fMaxHitRMS;
    std::vector<float>  fMinHitGOF;
    std::vector<float>  fMaxHitGOF;
    std::vector<float>  fMinHitRatio;
    std::vector<float>  fMaxHitRatio;
    int                 fMaxHitMult;
    float               fHitClustWidthFact;
    int                 fHitClustWireRange;
    float               fMatchQDiffLimit;
    float               fMatchMaxQRatio;
    //std::vector<float>  fTimeOffsets;
    float               fMatchMinOverlap;
    float               fMatchSigmaFact;
    float               fMatchMaxTicks;
    int                 fMaxWiresInCluster;
    float               fMinClusterCharge;
    float               fMaxClusterCharge;
    float               fMaxClusterSpan;
    int                 fMinMatchedPlanes;
    bool                fPickyBlips;
    bool                fApplyTrkCylinderCut;
    float               fCylinderRadius; 
    
    bool                fVetoBadChannels;
    std::string         fBadChanProducer;
    std::string         fBadChanFile;
    int                 fMinDeadWireGap;
    
    bool                keepAllClusts;
    bool                fKeepAllClusts[kNplanes];

    // --- Calorimetry configs ---
    int                 fCaloPlane;
    float               fCalodEdx;
    bool                fLifetimeCorr;
    bool                fSCECorr;
    bool                fYZUniformityCorr;
    float               fModBoxA;
    float               fModBoxB;
 

    // --- Histograms ---
    //TH1D*   h_chanstatus;
    //TH1D*   h_hit_chanstatus;
    TH1D*   h_chan_nhits;
    TH1D*   h_chan_nclusts;
    TH1D*   h_chan_bad;
    TH1D*   h_clust_nwires;
    TH1D*   h_clust_timespan;
    TH1D*   h_hit_maskfrac[kNplanes];
    TH1D*   h_hit_maskfrac_true[kNplanes];
    TH1D*   h_clust_overlap[kNplanes];
    TH1D*   h_clust_dt[kNplanes];
    TH1D*   h_clust_dtfrac[kNplanes];
    TH2D*   h_clust_q[kNplanes]; 
    TH2D*   h_clust_q_cut[kNplanes]; 
    TH1D*   h_clust_picky_overlap[kNplanes];
    TH1D*   h_clust_picky_dt[kNplanes];
    TH1D*   h_clust_picky_dtfrac[kNplanes];
    TH2D*   h_clust_picky_q[kNplanes]; 
    
    TH1D*   h_clust_truematch_overlap[kNplanes];
    TH1D*   h_clust_truematch_dt[kNplanes];
    TH1D*   h_clust_truematch_dtfrac[kNplanes];
    TH2D*   h_clust_truematch_q[kNplanes]; 

    TH1D*   h_nmatches[kNplanes];

    TH1D*   h_recomb;


  };

}

#endif
