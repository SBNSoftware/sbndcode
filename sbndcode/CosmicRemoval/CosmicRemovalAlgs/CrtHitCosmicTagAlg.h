#ifndef CRTHITCOSMICTAGALG_H_SEEN
#define CRTHITCOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// CrtHitCosmicTagAlg.h
//
// Functions for CRTHit match cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"

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

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class CrtHitCosmicTagAlg {
  public:
  
    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("Minimum t0 tagged time cut")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("Maximum t0 tagged time cut")
      };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<CRTT0MatchAlg::Config> T0Alg {
        Name("T0Alg"),
        Comment("Configuration for CRTHit matching")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("Configuration for t0 cut limits")
      };

    };

    CrtHitCosmicTagAlg(const Config& config);

    CrtHitCosmicTagAlg(const fhicl::ParameterSet& pset) :
      CrtHitCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    CrtHitCosmicTagAlg();

    ~CrtHitCosmicTagAlg();

    void reconfigure(const Config& config);

    // Returns true if matched to CRTHit outside beam time
    bool CrtHitCosmicTag(recob::Track track, std::vector<crt::CRTHit> crtHits, const art::Event& event);

  private:

    CRTT0MatchAlg t0Alg;
    double fBeamTimeMin;
    double fBeamTimeMax;

  };

}

#endif
