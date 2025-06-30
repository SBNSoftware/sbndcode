
// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// ROOT and C++ includes
#include <TTree.h>
#include <string.h>

// Services
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/Utilities/AssociationUtil.h"


// G4 includes
#include "nusimdata/SimulationBase/MCParticle.h"
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

// CRT
#include "sbnobj/SBND/CRT/CRTTrack.hh"

// Cosmic rejection includes
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "lardataobj/AnalysisBase/T0.h"

// Geometry and mapping
#include "larcore/Geometry/WireReadout.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include <iostream>
#include <fstream>

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
    void GetPropagationTimeCorrectionPerChannel();
    double GetPropagationTime(double );
    double GetFlashT0(double , std::vector<recob::OpHit> );

    //---TREE PARAMETERS
    TTree *fTree;
    art::ServiceHandle<art::TFileService> tfs;

    geo::WireReadoutGeom const& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

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
    std::string fOpHitsModuleLabel;
    std::string fOpFlashNewLabel;

    std::vector<double> fTimeCorrectionPerChannel;
    double fRecoVx = 0.0;
    double fRecoVy = 0.0;
    double fRecoVz = 0.0;

    //Space Point Variables
    std::vector<double> fSpacePointX;
    std::vector<double> fSpacePointY;
    std::vector<double> fSpacePointZ;
    std::vector<double> fSpacePointIntegral;

    //Charge Barycenter 
    std::vector<double> fChargeBarycenterX{0.,0.};
    std::vector<double> fChargeBarycenterY{0.,0.};
    std::vector<double> fChargeBarycenterZ{0.,0.};
    std::vector<double> fChargeWeightX{0.,0.};
    std::vector<double> fChargeWeightY{0.,0.};
    std::vector<double> fChargeWeightZ{0.,0.};
    std::vector<double> fChargeTotalWeight{0.,0.};

    double fDriftDistance; // Total Drift Distance
    double fSpeedOfLight; // Speed of light in mm/ns
    double fVISLightPropTime;
    double fKinkDistance;
    double fVGroupVUV_I;
    double fVGroupVIS;
    double fVGroupVUV; // Speed of light in vacuum in mm/ns


    double fPDFraction;
    double fPreWindow;
    double fPostWindow;
    double fMinHitPE;


};


DEFINE_ART_MODULE(sbnd::LightPropagationCorrection)
