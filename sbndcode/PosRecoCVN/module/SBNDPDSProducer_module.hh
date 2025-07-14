// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// ROOT and C++ includes
#include <string.h>
#include <vector>

struct PixelMapVars {
    std::vector<std::vector<float>> flash_ophit_pe;
    std::vector<std::vector<int>>   flash_ophit_ch;
    std::vector<std::vector<float>> flash_ophit_time;
    std::vector<float> nuvT;
    std::vector<float> dEpromx;
    std::vector<float> dEpromy;
    std::vector<float> dEpromz;
    std::vector<float> dEtpc;
    std::vector<float> nuvZ;
};

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
#include "larcore/Geometry/WireReadout.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"


#define xdet_size 1000
#define ydet_size 1000
#define zmindet_size -500
#define zmaxdet_size 1800

namespace opdet {
  class SBNDPDSProducer;

  enum SBNDPDSDetectorType {
    kPDUnknown = -1,   
    kPMTCoated = 0,       
    kPMTUncoated = 1,
    kXARAPUCAVUV,
    kXARAPUCAVIS
  };

}


class opdet::SBNDPDSProducer : public art::EDProducer {
public:
  explicit SBNDPDSProducer(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  SBNDPDSProducer(SBNDPDSProducer const&) = delete;
  SBNDPDSProducer(SBNDPDSProducer&&) = delete;
  SBNDPDSProducer& operator=(SBNDPDSProducer const&) = delete;
  SBNDPDSProducer& operator=(SBNDPDSProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  // void beginJob() override;

private:

  // Functions
  void FillMCTruth(art::Event const& e);

  void FillAverageDepositedEnergyVariables(std::vector<std::vector<double>> fenergydep, std::vector<std::vector<double>> fenergydepX,
  std::vector<std::vector<double>> fenergydepY, std::vector<std::vector<double>> fenergydepZ, std::vector<std::vector<double>> fstepT,
  std::vector<double> &dEtpc, std::vector<double> &dEpromx, std::vector<double> &dEpromy, std::vector<double> &dEpromz,
  std::vector<double> &dEspreadx, std::vector<double> &dEspready, std::vector<double> &dEspreadz,
  std::vector<std::vector<double>> &dElowedges, std::vector<std::vector<double>> &dEmaxedges);

  // Parámetros de configuración
  int fVerbosity;
  std::vector<std::string> fMCTruthModuleLabel;
  std::vector<std::string> fMCTruthInstanceLabel;
  std::string fMCModuleLabel;
  std::vector<std::string> fOpHitsModuleLabel;
  std::vector<std::string> fOpFlashesModuleLabel;
  std::vector<int> fG4BufferBoxX;
  std::vector<int> fG4BufferBoxY;
  std::vector<int> fG4BufferBoxZ;
  std::vector<int> fG4BeamWindow;

  // Variables internas necesarias
  std::vector<double> _nuvT;
  std::vector<double> _nuvZ;
  std::vector<double> _mc_dEpromx, _mc_dEpromy, _mc_dEpromz, _mc_dEtpc;
  std::vector<std::vector<float>> _flash_ophit_pe;
  std::vector<std::vector<int>> _flash_ophit_ch;
  std::vector<std::vector<float>> _flash_ophit_time;
};


DEFINE_ART_MODULE(opdet::SBNDPDSProducer)
