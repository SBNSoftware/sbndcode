#ifndef GEOMETRYCOSMICIDALG_H_SEEN
#define GEOMETRYCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// GeometryCosmicIdAlg.h
//
// Functions for fiducial volume cosmic tagger
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

// c++
#include <vector>


namespace sbnd{

  class GeometryCosmicIdAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

    };

    GeometryCosmicIdAlg(const Config& config);

    GeometryCosmicIdAlg(const fhicl::ParameterSet& pset) :
      GeometryCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    GeometryCosmicIdAlg();

    ~GeometryCosmicIdAlg();

    void reconfigure(const Config& config);

    bool GeometryCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, bool tpc0Flash, bool tpc1Flash);

  private:

    TPCGeoAlg fTpcGeo;

  };

}

#endif
