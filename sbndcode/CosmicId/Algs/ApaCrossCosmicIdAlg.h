#ifndef APACROSSCOSMICIDALG_H_SEEN
#define APACROSSCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// ApaCrossCosmicIdAlg.h
//
// Functions for APA crossing cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
namespace detinfo {
  class DetectorPropertiesData;
}

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

    // Get the minimum distance from track to APA for different times
    std::pair<double, double> MinApaDistance(detinfo::DetectorPropertiesData const& detProp,
                                             recob::Track track, std::vector<double> t0List, int tpc);

    // Get time by matching tracks which cross the APA
    double T0FromApaCross(detinfo::DetectorPropertiesData const& detProp,
                          recob::Track track, std::vector<double> t0List, int tpc);

    // Get the distance from track to APA at fixed time
    double ApaDistance(detinfo::DetectorPropertiesData const& detProp,
                       recob::Track track, double t0, std::vector<art::Ptr<recob::Hit>> hits);

    // Work out what TPC track is in and get the minimum distance from track to APA for different times
    std::pair<double, double> MinApaDistance(detinfo::DetectorPropertiesData const& detProp,
                                             recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1);

    // Tag tracks with times outside the beam
    bool ApaCrossCosmicId(detinfo::DetectorPropertiesData const& detProp,
                          recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1);

  private:

    double fDistanceLimit;
    double fMaxApaDistance;
    double fBeamTimeMin;
    double fBeamTimeMax;

    TPCGeoAlg fTpcGeo;

  };

}

#endif
