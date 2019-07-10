#ifndef CRTTRACKCOSMICIDALG_H_SEEN
#define CRTTRACKCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// CrtTrackCosmicIdAlg.h
//
// Functions for CRTTrack match cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

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
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class CrtTrackCosmicIdAlg {
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

      fhicl::Table<CRTTrackMatchAlg::Config> TrackMatchAlg {
        Name("TrackMatchAlg"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    CrtTrackCosmicIdAlg(const Config& config);

    CrtTrackCosmicIdAlg(const fhicl::ParameterSet& pset) :
      CrtTrackCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    CrtTrackCosmicIdAlg();

    ~CrtTrackCosmicIdAlg();

    void reconfigure(const Config& config);

    // Tags track as cosmic if it matches a CRTTrack
    bool CrtTrackCosmicId(recob::Track track, std::vector<crt::CRTTrack> crtTracks, const art::Event& event);

    // Getter
    CRTTrackMatchAlg TrackAlg() const {return trackMatchAlg;}

  private:

    CRTTrackMatchAlg trackMatchAlg;
    double fBeamTimeMin;
    double fBeamTimeMax;

  };

}

#endif
