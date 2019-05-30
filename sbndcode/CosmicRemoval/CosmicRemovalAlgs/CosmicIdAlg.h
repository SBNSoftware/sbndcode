#ifndef COSMICIDALG_H_SEEN
#define COSMICIDALG_H_SEEN


///////////////////////////////////////////////
// CosmicIdAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/FiducialVolumeCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/StoppingParticleCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/GeometryCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CpaCrossCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/ApaCrossCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CrtHitCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/CrtTrackCosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalAlgs/PandoraT0CosmicTagAlg.h"
#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class CosmicIdAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> CrtHitModuleLabel {
        Name("CrtHitModuleLabel"),
        Comment("tag of CRT hit producer data product")
      };

      fhicl::Atom<art::InputTag> CrtTrackModuleLabel {
        Name("CrtTrackModuleLabel"),
        Comment("tag of CRT track producer data product")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<art::InputTag> CaloModuleLabel {
        Name("CaloModuleLabel"),
        Comment("tag of calorimetry producer data product")
      };

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("tag of pandora data product")
      };

      fhicl::Atom<bool> ApplyFiducialCut {
        Name("ApplyFiducialCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyStoppingCut {
        Name("ApplyStoppingCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyGeometryCut {
        Name("ApplyGeometryCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyCpaCrossCut {
        Name("ApplyCpaCrossCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyApaCrossCut {
        Name("ApplyApaCrossCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyCrtTrackCut {
        Name("ApplyCrtTrackCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyCrtHitCut {
        Name("ApplyCrtHitCut"),
        Comment("")
      };

      fhicl::Atom<bool> ApplyPandoraT0Cut {
        Name("ApplyPandoraT0Cut"),
        Comment("")
      };

      fhicl::Atom<bool> UseTrackAngleVeto {
        Name("UseTrackAngleVeto"),
        Comment("")
      };

      fhicl::Atom<double> MinSecondTrackLength {
        Name("MinSecondTrackLength"),
        Comment("")
      };

      fhicl::Atom<double> MinVertexDistance {
        Name("MinVertexDistance"),
        Comment("")
      };

      fhicl::Atom<double> MinMergeAngle {
        Name("MinMergeAngle"),
        Comment("")
      };

      fhicl::Table<FiducialVolumeCosmicTagAlg::Config> FVTagAlg {
        Name("FVTagAlg"),
      };

      fhicl::Table<StoppingParticleCosmicTagAlg::Config> SPTagAlg {
        Name("SPTagAlg"),
      };

      fhicl::Table<CpaCrossCosmicTagAlg::Config> CCTagAlg {
        Name("CCTagAlg"),
      };

      fhicl::Table<ApaCrossCosmicTagAlg::Config> ACTagAlg {
        Name("ACTagAlg"),
      };

      fhicl::Table<CrtHitCosmicTagAlg::Config> CHTagAlg {
        Name("CHTagAlg"),
      };

      fhicl::Table<CrtTrackCosmicTagAlg::Config> CTTagAlg {
        Name("CTTagAlg"),
      };

      fhicl::Table<PandoraT0CosmicTagAlg::Config> PTTagAlg {
        Name("PTTagAlg"),
      };

    };

    CosmicIdAlg(const Config& config);

    CosmicIdAlg(const fhicl::ParameterSet& pset) :
      CosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    CosmicIdAlg();

    ~CosmicIdAlg();

    void reconfigure(const Config& config);

    void SetCuts(bool FV, bool SP, bool Geo, bool CC, bool AC, bool CT, bool CH, bool PT);

    bool CosmicId(recob::Track track, const art::Event& event, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1);

    bool CosmicId(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1);

  private:

    art::InputTag fTpcTrackModuleLabel;
    art::InputTag fPandoraLabel;
    art::InputTag fCrtHitModuleLabel;
    art::InputTag fCrtTrackModuleLabel;
    art::InputTag fCaloModuleLabel;

    bool fApplyFiducialCut;
    bool fApplyStoppingCut;
    bool fApplyGeometryCut;
    bool fApplyCpaCrossCut;
    bool fApplyApaCrossCut;
    bool fApplyCrtTrackCut;
    bool fApplyCrtHitCut;
    bool fApplyPandoraT0Cut;

    bool fUseTrackAngleVeto;
    double fMinSecondTrackLength;
    double fMinVertexDistance;
    double fMinMergeAngle;

    FiducialVolumeCosmicTagAlg    fvTag;
    StoppingParticleCosmicTagAlg  spTag;
    GeometryCosmicTagAlg          geoTag;
    CpaCrossCosmicTagAlg          ccTag;
    ApaCrossCosmicTagAlg          acTag;
    CrtHitCosmicTagAlg            chTag;
    CrtTrackCosmicTagAlg          ctTag;
    PandoraT0CosmicTagAlg         ptTag;

  };

}

#endif
