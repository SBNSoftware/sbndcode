#ifndef APACROSSCOSMICIDALG_H_SEEN
#define APACROSSCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// ApaCrossCosmicIdAlg.h
//
// Functions for APA crossing cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class ApaCrossCosmicIdAlg {
  public:

    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("")
      };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("")
      };

      fhicl::Atom<double> MaxApaDistance {
        Name("MaxApaDistance"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    ApaCrossCosmicIdAlg(const Config& config);

    ApaCrossCosmicIdAlg(const fhicl::ParameterSet& pset) :
      ApaCrossCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    ApaCrossCosmicIdAlg();

    ~ApaCrossCosmicIdAlg();

    void reconfigure(const Config& config);

    // Get time by matching tracks which cross the APA
    double T0FromApaCross(recob::Track track, std::vector<double> t0List, int tpc);

    // Tag tracks with times outside the beam
    bool ApaCrossCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1);

  private:

    double fDistanceLimit;
    double fMaxApaDistance;
    double fBeamTimeMin;
    double fBeamTimeMax;

    detinfo::DetectorProperties const* fDetectorProperties;
    TPCGeoAlg fTpcGeo;

  };

}

#endif
