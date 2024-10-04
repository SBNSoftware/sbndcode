////////////////////////////////////////////////////////////////////////
// Class:       SBNDPDSAnalyzer
// Plugin Type: analyzer
// File:        SBNDPDSAnalyzer_module.hh
//
// Created by Francisco Nicolas-Arnaldos using cetskelgen
////////////////////////////////////////////////////////////////////////

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
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

// Cosmic rejection includes
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "lardataobj/AnalysisBase/T0.h"

// Geometry and mapping
#include "larcore/Geometry/Geometry.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "sbndcode/SERCalibration/Alg/SERPulseFinderBase.hh"

#include "TH1D.h"


namespace opdet {
  class SERCalibration;
}

class opdet::SERCalibration : public art::EDAnalyzer {

  public:
    explicit SERCalibration(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SERCalibration(SERCalibration const &) = delete;
    SERCalibration(SERCalibration &&) = delete;
    SERCalibration & operator = (SERCalibration const &) = delete;
    SERCalibration & operator = (SERCalibration &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    //Selected optional functions
    void beginJob() override;
    void endJob() override;

   private:

    // Declare member data here.
    std::string fInputLabel;
    std::vector<std::string> fPDTypes;
    std::vector<std::string> fElectronics;
    //OpDecoAlg tool
    std::unique_ptr<opdet::SERPulseFinderBase> fSERPulseFinderPtr;
    //PDS map
    opdet::sbndPDMapAlg pdsmap;
    std::vector<TH1D> calibratedSER_v;
    int fSERStart=-200; //TTicks before the peak to build the SER
    int fSEREnd=300; //TTicks after the peak to build the SER

    int _eventID;
    int _runID;
    int _subrunID;
    int fVerbosity=1;

    TTree * fTree;

};