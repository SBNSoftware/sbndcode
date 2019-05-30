#ifndef GEOMETRYCOSMICTAGALG_H_SEEN
#define GEOMETRYCOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// GeometryCosmicTagAlg.h
//
// Functions for fiducial volume cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

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

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>


namespace sbnd{

  class GeometryCosmicTagAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

    };

    GeometryCosmicTagAlg(const Config& config);

    GeometryCosmicTagAlg(const fhicl::ParameterSet& pset) :
      GeometryCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    GeometryCosmicTagAlg();

    ~GeometryCosmicTagAlg();

    void reconfigure(const Config& config);

    bool GeometryCosmicTag(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, bool tpc0Flash, bool tpc1Flash);

  private:

  };

}

#endif
