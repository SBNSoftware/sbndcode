#ifndef FIDUCIALVOLUMECOSMICIDALG_H_SEEN
#define FIDUCIALVOLUMECOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// FiducialVolumeCosmicIdAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"

// c++
#include <vector>


namespace sbnd{

  class FiducialVolumeCosmicIdAlg {
  public:

    struct Fiducial {
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

      fhicl::Table<Fiducial> FiducialCuts {
        Name("FiducialCuts"),
        Comment("Fiducial volume cuts (cm)")
      };

    };

    FiducialVolumeCosmicIdAlg(const Config& config);

    FiducialVolumeCosmicIdAlg(const fhicl::ParameterSet& pset) :
      FiducialVolumeCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    FiducialVolumeCosmicIdAlg();

    ~FiducialVolumeCosmicIdAlg();

    void reconfigure(const Config& config);

    // Check if point in fiducial volume used by this algorithm
    bool InFiducial(geo::Point_t point);

    // Check both start and end points of track are in fiducial volume
    bool FiducialVolumeCosmicId(recob::Track track);

  private:

    double fMinX;
    double fMinY;
    double fMinZ;
    double fMaxX;
    double fMaxY;
    double fMaxZ;

    TPCGeoAlg fTpcGeo;

  };

}

#endif
