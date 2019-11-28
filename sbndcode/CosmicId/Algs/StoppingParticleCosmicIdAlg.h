#ifndef STOPPINGPARTICLECOSMICIDALG_H_SEEN
#define STOPPINGPARTICLECOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// StoppingParticleCosmicIdAlg.h
//
// Functions for stopping particle cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// ROOT
#include "TGraph.h"
#include "TF1.h"

// c++
#include <vector>

namespace sbnd{

  class StoppingParticleCosmicIdAlg {
  public:

    struct Containment {
      using Name = fhicl::Name;

      fhicl::Atom<double> MinX { Name("MinX") };
      fhicl::Atom<double> MinY { Name("MinY") };
      fhicl::Atom<double> MinZ { Name("MinZ") };
      fhicl::Atom<double> MaxX { Name("MaxX") };
      fhicl::Atom<double> MaxY { Name("MaxY") };
      fhicl::Atom<double> MaxZ { Name("MaxZ") };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<Containment> ContainmentCuts {
        Name("ContainmentCuts"),
        Comment("Fiducial volume cut (cm) to decide if track exits")
      };

      fhicl::Atom<double> ResRangeMin {
        Name("ResRangeMin"),
        Comment("Minumum residual range (cm) to fit")
      };

      fhicl::Atom<double> ResRangeMax {
        Name("ResRangeMax"),
        Comment("Maximum residual range (cm) to fit")
      };

      fhicl::Atom<double> DEdxMax {
        Name("DEdxMax"),
        Comment("Maximum dE/dx (MeV/cm) to fit")
      };

      fhicl::Atom<double> StoppingChi2Limit {
        Name("StoppingChi2Limit"),
        Comment("Limit of pol/exp chi2 ratio to cut on to determine if stopping")
      };

    };

    StoppingParticleCosmicIdAlg(const Config& config);

    StoppingParticleCosmicIdAlg(const fhicl::ParameterSet& pset) :
      StoppingParticleCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    StoppingParticleCosmicIdAlg();

    ~StoppingParticleCosmicIdAlg();

    void reconfigure(const Config& config);

    // Calculate the chi2 ratio of pol0 and exp fit to dE/dx vs residual range
    double StoppingChiSq(geo::Point_t end, std::vector<art::Ptr<anab::Calorimetry>> calos);

    // Determine if the track end looks like it stops
    bool StoppingEnd(geo::Point_t end, std::vector<art::Ptr<anab::Calorimetry>> calos);

    // Determine if a track looks like a stopping cosmic
    bool StoppingParticleCosmicId(recob::Track track, std::vector<art::Ptr<anab::Calorimetry>> calos);

    // Determine if two tracks look like a stopping cosmic if they are merged
    bool StoppingParticleCosmicId(recob::Track track, recob::Track track2, std::vector<art::Ptr<anab::Calorimetry>> calos, std::vector<art::Ptr<anab::Calorimetry>> calos2);

    // Calculate chi2 for particles which do not stop
    std::vector<double> FlatChi2(std::vector<art::Ptr<anab::Calorimetry>> calos);

  private:

    double fMinX;
    double fMinY;
    double fMinZ;
    double fMaxX;
    double fMaxY;
    double fMaxZ;
    double fResRangeMin;
    double fResRangeMax;
    double fDEdxMax;
    double fStoppingChi2Limit;

    TPCGeoAlg fTpcGeo;

  };

}

#endif
