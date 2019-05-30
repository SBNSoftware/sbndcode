#ifndef CPACROSSCOSMICTAGALG_H_SEEN
#define CPACROSSCOSMICTAGALG_H_SEEN


///////////////////////////////////////////////
// CpaCrossCosmicTagAlg.h
//
// Functions for CPA stitching cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

#include "sbndcode/CosmicRemoval/CosmicRemovalUtils/CosmicRemovalUtils.h"
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

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class CpaCrossCosmicTagAlg {
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

      fhicl::Atom<double> CpaStitchDistance {
        Name("CpaStitchDistance"),
        Comment("")
      };

      fhicl::Atom<double> CpaStitchAngle {
        Name("CpaStitchAngle"),
        Comment("")
      };

      fhicl::Atom<double> CpaXDifference {
        Name("CpaXDifference"),
        Comment("")
      };

      fhicl::Table<Fiducial> FiducialCuts {
        Name("FiducialCuts"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    CpaCrossCosmicTagAlg(const Config& config);

    CpaCrossCosmicTagAlg(const fhicl::ParameterSet& pset) :
      CpaCrossCosmicTagAlg(fhicl::Table<Config>(pset, {})()) {}

    CpaCrossCosmicTagAlg();

    ~CpaCrossCosmicTagAlg();

    void reconfigure(const Config& config);

    // Calculate the time by stitching tracks across the CPA
    std::pair<double, bool> T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks);

    // Tag tracks as cosmics from CPA stitching t0
    bool CpaCrossCosmicTag(recob::Track track, std::vector<recob::Track> tracks, art::FindManyP<recob::Hit> hitAssoc);

  private:

    double fCpaStitchDistance;
    double fCpaStitchAngle;
    double fCpaXDifference;
    double fMinX;
    double fMinY;
    double fMinZ;
    double fMaxX;
    double fMaxY;
    double fMaxZ;
    double fBeamTimeMin;
    double fBeamTimeMax;

    detinfo::DetectorProperties const* fDetectorProperties;
    TPCGeoAlg fTpcGeo;

  };

}

#endif
