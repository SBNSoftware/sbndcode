
#ifndef SBND_LIGHT_PROPAGATION_CORRECTION_MODULE_HH
#define SBND_LIGHT_PROPAGATION_CORRECTION_MODULE_HH

// -----------------------------------------------------------------------------
// C++ standard library
// -----------------------------------------------------------------------------

#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// ROOT
// -----------------------------------------------------------------------------

#include <TTree.h>

// -----------------------------------------------------------------------------
// art framework
// -----------------------------------------------------------------------------

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"

// -----------------------------------------------------------------------------
// canvas and FHiCL
// -----------------------------------------------------------------------------

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// -----------------------------------------------------------------------------
// LArSoft services and geometry
// -----------------------------------------------------------------------------

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// -----------------------------------------------------------------------------
// Reconstruction objects
// -----------------------------------------------------------------------------

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"

// -----------------------------------------------------------------------------
// Pandora and reconstruction algorithms
// -----------------------------------------------------------------------------

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// -----------------------------------------------------------------------------
// SBND data products and geometry
// -----------------------------------------------------------------------------

#include "sbnobj/Common/Reco/CorrectedOpFlashTiming.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

// -----------------------------------------------------------------------------
// Flash reconstruction and timing tools
// -----------------------------------------------------------------------------

#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderManager.h"
#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderFMWKInterface.h"
#include "sbndcode/OpDetReco/OpFlash/FlashTools/FlashGeoBase.hh"
#include "sbndcode/OpDetReco/OpFlash/FlashTools/FlashT0Base.hh"


namespace sbnd {

class LightPropagationCorrection : public art::EDProducer {

public:

    // -------------------------------------------------------------------------
    // Construction and framework interface
    // -------------------------------------------------------------------------

    explicit LightPropagationCorrection(fhicl::ParameterSet const& p);

    LightPropagationCorrection(LightPropagationCorrection const&) = delete;
    LightPropagationCorrection(LightPropagationCorrection&&) = delete;
    LightPropagationCorrection& operator=(LightPropagationCorrection const&) = delete;
    LightPropagationCorrection& operator=(LightPropagationCorrection&&) = delete;

private:

    // -------------------------------------------------------------------------
    // art module lifecycle
    // -------------------------------------------------------------------------

    void produce(art::Event& e) override;
    void beginJob() override;
    void endJob() override;

    // -------------------------------------------------------------------------
    // Event and slice bookkeeping
    // -------------------------------------------------------------------------

    void ResetEventVars();
    void ResetSliceInfo();

    size_t HighestOpT0ScoreIdx(std::vector<art::Ptr<sbn::OpT0Finder>> const& opt0);
    size_t HighestBFMScoreIdx(std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>> const& matches);

    // -------------------------------------------------------------------------
    // Particle identification and propagation time
    // -------------------------------------------------------------------------

    int GetTrackPID(art::FindManyP<anab::ParticleID> const& trackToPIDAssns, art::Ptr<recob::Track> const& track);

    double GetTrackMomentum(int pdg, art::FindManyP<sbn::RangeP> const& muonRangeAssns, art::FindManyP<sbn::RangeP> const& protonRangeAssns, art::FindManyP<recob::MCSFitResult> const& muonMCSAssns, art::FindManyP<recob::MCSFitResult> const& protonMCSAssns, art::Ptr<recob::Track> const& track);

    bool GetParticlePropagationTime(int pdg, double initialMomentum, art::FindManyP<anab::Calorimetry> const& trackToCaloAssns, art::Ptr<recob::Track> const& track, int pfpID, int parentPFPID);

    void GetParticlePropagationTimeLite(int pfpID, std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints, art::FindManyP<recob::Hit> const& spacePointToHitAssns);

    bool IsInFV(geo::Point_t const& start, geo::Point_t const& end);

    // -------------------------------------------------------------------------
    // Particle trajectory representation
    // -------------------------------------------------------------------------

    struct ParticleSpacePoint_t {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double propagationTime = 0.0;
        int pfpID = -1;
    };

    std::vector<ParticleSpacePoint_t> fParticleSpacePoints;

    // -------------------------------------------------------------------------
    // Light propagation corrections
    // -------------------------------------------------------------------------

    void GetPropagationTimeCorrectionPerChannelScintillation();
    void GetPropagationTimeCorrectionPerChannelCherenkov();

    // -------------------------------------------------------------------------
    // OpHit and OpFlash correction
    // -------------------------------------------------------------------------

    void CorrectOpHitTime(std::vector<art::Ptr<recob::OpHit>> const& oldOpHits,
                          std::vector<recob::OpHit>& newOpHits);

    void FillLiteOpHit(std::vector<recob::OpHit> const& opHits,
                       std::vector<::lightana::LiteOpHit_t>& liteOpHits);

    void FillFirstOpHit(double flash_time, 
                        std::vector<recob::OpHit> const& all_opHits, 
                        std::vector<recob::OpHit>& first_ophits);

    ::lightana::LiteOpHitArray_t GetAssociatedLiteHits(::lightana::LiteOpFlash_t lite_flash, 
                                                       ::lightana::LiteOpHitArray_t lite_hits_v);

    void GetBoxTimeDiffMap(std::vector<recob::OpHit>& first_ophits);

    bool CorrectOpFlash(art::Ptr<recob::OpFlash> const& flash,
                        sbn::CorrectedOpFlashTiming& correctedTiming,
                        art::FindManyP<recob::OpHit> const& flashToOpHitAssnsTPC0,
                        art::FindManyP<recob::OpHit> const& flashToOpHitAssnsTPC1);



    // -------------------------------------------------------------------------
    // Debugging and output
    // -------------------------------------------------------------------------

    void FillCorrectionTree(double& newFlashTime,
                            recob::OpFlash const& flash,
                            std::vector<recob::OpHit> const& oldOpHits,
                            std::vector<recob::OpHit> const& newOpHits);

    // -------------------------------------------------------------------------
    // Detector geometry and services
    // -------------------------------------------------------------------------

    geo::WireReadoutGeom const& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

    opdet::sbndPDMapAlg fPDSMap;
    std::map<size_t, size_t> fChannelBoxMap;
    std::map<size_t, size_t> fUncoatedChannelMap;
    std::map<size_t, double> fBoxTimeDiffMap;

    trkf::TrackMomentumCalculator tmc;

    // -------------------------------------------------------------------------
    // Flash reconstruction tools
    // -------------------------------------------------------------------------

    ::lightana::FlashFinderManager _mgr;
    ::lightana::FlashFinderManager _mgr_tpc0;
    ::lightana::FlashFinderManager _mgr_tpc1;

    std::unique_ptr<lightana::FlashGeoBase> _flashgeo;
    std::unique_ptr<lightana::FlashT0Base> _flasht0calculator;

    // -------------------------------------------------------------------------
    // Input labels
    // -------------------------------------------------------------------------

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
    std::string fFlashMatchingTool;

    // -------------------------------------------------------------------------
    // Detector optical channel information
    // -------------------------------------------------------------------------

    size_t fNOpChannels = 0;

    std::vector<size_t> fOpDetID;
    std::vector<double> fOpDetX;
    std::vector<double> fOpDetY;
    std::vector<double> fOpDetZ;

    std::vector<double> fTimeCorrectionPerChannelScintillation;
    std::vector<double> fTimeCorrectionPerChannelCherenkov;

    // -------------------------------------------------------------------------
    // Particle and light propagation configuration
    // -------------------------------------------------------------------------

    bool fDoPIDCorrection;
    bool fUseMCVIS;

    double fSpeedOfLight;
    double fVGroupVIS;
    double fVGroupVUV;

    double fMIPdEdX; // In MeV/cm, for muons
    double fHIPdEdX; // In MeV/cm, for protons

    double fCherenkovRefractiveIndex;
    double fCherenkovCosineTolerance;
    double fCherenkovPMTRadius;
    unsigned int fCherenkovPMTRings;
    unsigned int fCherenkovPMTAzimuthSteps;
    double fMaxCoatedUncoatedTimeDiff;

    double fMinHitPE;
    double fPreWindow;
    double fPostWindow;


    // -------------------------------------------------------------------------
    // Fiducial volume and event selection
    // -------------------------------------------------------------------------

    std::vector<double> fFVx;
    std::vector<double> fFVy;
    std::vector<double> fFVz;

    double fNuScoreThreshold;
    double fFMScoreThreshold;
    double fTrackScoreThreshold;

    // -------------------------------------------------------------------------
    // Reconstructed and true interaction vertices
    // -------------------------------------------------------------------------

    double fRecoVx = 0.0;
    double fRecoVy = 0.0;
    double fRecoVz = 0.0;

    double fTrueVx = 0.0;
    double fTrueVy = 0.0;
    double fTrueVz = 0.0;
    double fTrueVt = 0.0;

    int fNeutrinoID = -1;
    bool fAllPFPsCorrected = true;

    // -------------------------------------------------------------------------
    // Debugging configuration and ROOT output
    // -------------------------------------------------------------------------

    bool fSaveCorrectionTree;
    bool fDebug;

    art::ServiceHandle<art::TFileService> tfs;
    TTree* fTree = nullptr;

    // -------------------------------------------------------------------------
    // Event identifiers and scalar bookkeeping
    // -------------------------------------------------------------------------

    int fEvent = 0;
    int fRun = 0;
    int fSubrun = 0;

    double _fNuScore = 0.0;
    double _fFMScore = 0.0;

    // -------------------------------------------------------------------------
    // ROOT output vectors
    // -------------------------------------------------------------------------

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
    std::vector<double> fSliceSPID;

    std::vector<std::vector<double>> fOpHitOldTime;
    std::vector<std::vector<double>> fOpHitNewTime;
    std::vector<std::vector<double>> fOpHitPE;
    std::vector<std::vector<int>> fOpHitOpCh;
};

} // namespace sbnd

DEFINE_ART_MODULE(sbnd::LightPropagationCorrection)

#endif // SBND_LIGHT_PROPAGATION_CORRECTION_MODULE_HH
