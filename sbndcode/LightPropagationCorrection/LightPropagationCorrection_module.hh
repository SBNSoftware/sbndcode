
// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// Services
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT and C++ includes
#include <TTree.h>
#include <string.h>

// Services
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/Utils/TruthMatchUtils.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"


// G4 includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"


// Reco includes
// PDS
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
// TPC
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/MCSFitResult.h"


// CRT
#include "sbnobj/SBND/CRT/CRTTrack.hh"

// Cosmic rejection includes
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/Common/Reco/CorrectedOpFlashTiming.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/Common/Reco/RangeP.h"
#include "lardataobj/AnalysisBase/T0.h"

// Geometry and mapping
#include "larcore/Geometry/WireReadout.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include <iostream>
#include <fstream>

// Flash finder utilities
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderManager.h"
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderFMWKInterface.h"
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/PECalib.h"
#include "sbndcode/OpDetReco/OpFlash/FlashTools/FlashGeoBase.hh"
#include "sbndcode/OpDetReco/OpFlash/FlashTools/FlashT0Base.hh"
#include "sbndcode/OpDetReco/OpFlash/FlashTools/DriftEstimatorBase.hh"


#define fXFidCut1 1.5
#define fXFidCut2 190
#define fYFidCut 190
#define fZFidCut1 10
#define fZFidCut2 490

#define xdet_size 1000
#define ydet_size 1000
#define zmindet_size -500
#define zmaxdet_size 1800

#define fDefaulNeutrinoID 99999


namespace sbnd {
  class LightPropagationCorrection;

}

class sbnd::LightPropagationCorrection : public art::EDProducer {
public:
    explicit LightPropagationCorrection(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    LightPropagationCorrection(LightPropagationCorrection const&) = delete;
    LightPropagationCorrection(LightPropagationCorrection&&) = delete;
    LightPropagationCorrection& operator=(LightPropagationCorrection const&) = delete;
    LightPropagationCorrection& operator=(LightPropagationCorrection&&) = delete;


private:

    // Required functions.
    void produce(art::Event & e) override;
    void endJob() override;
    void beginJob() override;

    // Selected optional functions.
    void ResetEventVars();
    void ResetSliceInfo();
    size_t HighestOpT0ScoreIdx(const std::vector< art::Ptr<sbn::OpT0Finder> >);
    size_t HighestBFMScoreIdx(const std::vector< art::Ptr<sbn::TPCPMTBarycenterMatch> > );

    void GetPropagationTimeCorrectionPerChannel();
    void CorrectOpHitTime(std::vector<art::Ptr<recob::OpHit>> , std::vector<recob::OpHit> & );
    void FillLiteOpHit(std::vector<recob::OpHit> const& , std::vector<::lightana::LiteOpHit_t>& );
    void FillCorrectionTree(double & , recob::OpFlash const& , std::vector<recob::OpHit> const& , std::vector<recob::OpHit> const& );
    void CorrectOpFlash(art::Ptr<recob::OpFlash> const& flash, sbn::CorrectedOpFlashTiming &correctedOpFlashTiming, bool matched, int tpc);
    int GetTrackPID(art::FindManyP<anab::ParticleID> , art::Ptr<recob::Track>  );
    double GetTrackMomentum(int, art::FindManyP<sbn::RangeP>, art::FindManyP<sbn::RangeP>, art::FindManyP<recob::MCSFitResult> , art::FindManyP<recob::MCSFitResult>  , art::Ptr<recob::Track>);
    void SaveTrueTrajectory(art::Event const&);
    void GetMCNeutrino(art::Event const& );
    bool GetParticlePropagationTime(int, double, art::FindManyP<anab::Calorimetry>, art::Ptr<recob::Track>, int, int);
    void GetParticlePropagationTimeLite(int, std::vector<art::Ptr<recob::SpacePoint>>, art::FindManyP<recob::Hit>);
    void GetParticlePropagationTimeMC(art::Event const& , const std::vector<art::Ptr<recob::Hit> >& , int );
    void GetTruthMatchedID(art::Event const& , const std::vector<art::Ptr<recob::Hit> >& );

    ::lightana::LiteOpHitArray_t GetAssociatedLiteHits(::lightana::LiteOpFlash_t , ::lightana::LiteOpHitArray_t );

    trkf::TrackMomentumCalculator tmc;


    geo::WireReadoutGeom const& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    opdet::sbndPDMapAlg fPDSMap;
    //Flash finder manager
    ::lightana::FlashFinderManager _mgr;
    ::lightana::FlashFinderManager _mgr_tpc0;
    ::lightana::FlashFinderManager _mgr_tpc1;

    std::unique_ptr<art::FindManyP<recob::OpHit>> flashToOpHitAssns_tpc0;
    std::unique_ptr<art::FindManyP<recob::OpHit>> flashToOpHitAssns_tpc1;
    // Tool for calculating the OpFlash Y and Z centers
    std::unique_ptr<lightana::FlashGeoBase> _flashgeo;

    // Tool for calculating the OpFlash t0
    std::unique_ptr<lightana::FlashT0Base> _flasht0calculator;

    //Vector for PMT position
    std::vector<double> fOpDetID;
    std::vector<double> fOpDetX;
    std::vector<double> fOpDetY;
    std::vector<double> fOpDetZ;

    std::string fReco2Label;
    std::string fOpT0FinderModuleLabel;
    std::string fTPCPMTBarycenterFMModuleLabel;
    std::string fOpFlashLabel_tpc0;
    std::string fOpFlashLabel_tpc1;
    std::string fSpacePointLabel;
    std::string fTPCTrackLabel;
    std::string fParticleIDLabel;
    std::string fCalorimetryLabel;
    std::string fPandoraRangeLabel;
    std::string fPandoraMCSLabel;
    std::string fOpHitsModuleLabel;
    bool fUseCheatMC; 
    bool fUseMCVIS;
    std::string fMCModuleLabel;
    std::vector<std::string> fMCTruthModuleLabel;
    std::vector<std::string> fMCTruthInstanceLabel;
    std::vector<int> fMCTruthOrigin;
    std::vector<int> fMCTruthPDG;
    std::string fSimDepProducer;
    
    std::string fFlashMatchingTool;
    
    bool fSaveCorrectionTree;

    std::vector<double> fTimeCorrectionPerChannel;
    std::vector<double> fParticlePropagationTimePerChannel;
    std::vector<double> fPhotonPropagationTimePerChannel;
    std::vector<double> fPhotonPropagationDistancePerChannel;
    std::vector<size_t> fSelectedChannelList;

    double fRecoVx = 0.0;
    double fRecoVy = 0.0;
    double fRecoVz = 0.0;

    double fTrueVx = 0.0;
    double fTrueVy = 0.0;
    double fTrueVz = 0.0;
    double fTrueVt = 0.0;

    //Space Point Variables
    std::vector<double> fSpacePointX;
    std::vector<double> fSpacePointY;
    std::vector<double> fSpacePointZ;
    std::vector<double> fSpacePointPropagationTime;
    std::vector<double> fSpacePointIntegral;
    std::vector<double> fSpacePointMomentum;
    std::vector<double> fSpacePointEnergy;
    std::vector<double> fSpacePointBeta;
    std::vector<int> fSpacePointPFPID;
    std::vector<int> fSpacePointTMID;

    std::vector<double> fTruthMatchedX;
    std::vector<double> fTruthMatchedY;
    std::vector<double> fTruthMatchedZ;
    std::vector<double> fTruthMatchedT;
    std::vector<double> fTruthMatchedID;
    std::vector<double> fTruthMatchedE;

    //Charge Barycenter 
    std::vector<double> fChargeBarycenterX{0.,0.};
    std::vector<double> fChargeBarycenterY{0.,0.};
    std::vector<double> fChargeBarycenterZ{0.,0.};
    std::vector<double> fChargeWeightX{0.,0.};
    std::vector<double> fChargeWeightY{0.,0.};
    std::vector<double> fChargeWeightZ{0.,0.};
    std::vector<double> fChargeTotalWeight{0.,0.};

    bool fDoPIDCorrection;

    double fDriftDistance; // Total Drift Distance
    double fSpeedOfLight; // Speed of light in mm/ns
    double fVISLightPropTime;
    double fKinkDistance;
    double fVGroupVUV_I;
    double fVGroupVIS;
    double fVGroupVUV; // Speed of light in vacuum in mm/ns

    double fNuScoreThreshold;
    double fFMScoreThreshold;


    bool fDebug;

    art::ServiceHandle<art::TFileService> tfs;
    TTree *fTree;

    int fNeutrinoID;
    bool fAllPFPsCorrected;

    int fEvent;
    int fRun;
    int fSubrun;
    size_t fNOpChannels;
    double _fNuScore;
    double _fFMScore;
    std::vector<double> fNuScore;
    std::vector<double> fFMScore;
    std::vector<double> fOpFlashTimeOld;
    std::vector<double> fOpFlashTimeNew;
    std::vector<double> fOpFlashXCenter;
    std::vector<double> fOpFlashYCenter;
    std::vector<double> fOpFlashZCenter;
    std::vector<double> fOpFlashPE;
    std::vector<double> fSliceVx;
    std::vector<double> fSliceVy;
    std::vector<double> fSliceVz;
    std::vector<double> fSliceSPX;
    std::vector<double> fSliceSPY;
    std::vector<double> fSliceSPZ;
    std::vector<double> fSliceSPT;
    std::vector<double> fSliceSPMomentum;
    std::vector<double> fSliceSPEnergy;
    std::vector<double> fSliceSPBeta;
    std::vector<double> fSliceSPID;
    std::vector<int> fSliceSPTMID;
    std::vector<double> fRecoMomentum; 
    std::vector<int> fRecoPID;
    std::vector<double> fTMX;
    std::vector<double> fTMY;
    std::vector<double> fTMZ;
    std::vector<double> fTMT;
    std::vector<double> fTME;
    std::vector<double> fTMID;
    std::vector<double> fTMPartID;
    std::vector<double> truedE;
    std::vector<double> truedE_vecX;
    std::vector<double> truedE_vecY;
    std::vector<double> truedE_vecZ;
    std::vector<double> truedE_vecID;
    std::vector<double> truedE_vecT;

    std::vector<double> truedE_;
    std::vector<double> truedE_X;
    std::vector<double> truedE_Y;
    std::vector<double> truedE_Z;
    std::vector<double> truedE_ID;
    std::vector<double> truedE_T;

    std::vector<double> fTrueMomentum;
    std::vector<int> fTruePID;

    int fTruthMatchedTrackID;

    std::vector<std::vector<double>> fOpHitOldTime;
    std::vector<std::vector<double>> fOpHitNewTime;
    std::vector<std::vector<double>> fOpHitPE;
    std::vector<std::vector<int>> fOpHitOpCh;
};


DEFINE_ART_MODULE(sbnd::LightPropagationCorrection)
