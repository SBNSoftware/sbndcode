////////////////////////////////////////////////////////////////////////
// Class:       LambdaPID
// Plugin Type: Producer 
// File:        LambdaPID_module.hh
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/Utils/TruthMatchUtils.h"

// truth includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// reco includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

// Lambda PID includes
#include "sbnobj/SBND/QELambdaReco/LambdaPID.hh"


// Define relevant invariant masses
#define ProtonMass 938.272
#define PionMass 139.57018

namespace qelambda {
  class LambdaPID;
}

class qelambda::LambdaPID : public art::EDProducer {
public:

  explicit LambdaPID(fhicl::ParameterSet const& p);

  // Required functions.
  void produce(art::Event &e) override;

  std::vector<art::Ptr<recob::PFParticle>> getSlicePrimPfp(art::Event const &e,std::vector<art::Ptr<recob::PFParticle>> slcPfpVec);

private:

  // Configuration parameters
  std::string fReco2Label;
  std::string fTrackLabel;
  std::string fCalorimetryLabel;
  std::string fParticleIDLabel;

  bool fRecoStatus;

  TVector3 fProtonMom;
  TVector3 fPionMom;

  // Lambda PID object
  qelambda::QELambdaPID fLambdaPID;

};

DEFINE_ART_MODULE(qelambda::LambdaPID)
