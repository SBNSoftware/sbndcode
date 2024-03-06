#ifndef CRTEVENTDISPLAYALG_H_SEEN
#define CRTEVENTDISPLAYALG_H_SEEN

///////////////////////////////////////////////
// CRTEventDisplayAlg.h
//
// Quick and dirty event display for SBND CRT
// T Brooks (tbrooks@fnal.gov), November 2018
// Edited heavily - Henry Lay, November 2022
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

//larsoft
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

// sbnobj
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"

// ROOT
#include "TPolyLine3D.h"
#include "TCanvas.h"

namespace detinfo { class DetectorClocksData; }

namespace sbnd::crt {
  
  class CRTEventDisplayAlg {
  public:
    
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Table<fhicl::ParameterSet> GeoAlgConfig {
        Name("CRTGeoAlg"),
        Comment("Configuration parameters for the CRT geometry algorithm"),
        fhicl::ParameterSet()
      };

      fhicl::Table<CRTBackTrackerAlg::Config> BackTrackerAlgConfig {
        Name("CRTBackTrackerAlg"),
        Comment("Configuration parameters for the CRT back tracking algorithm")
      };
      
      fhicl::Atom<art::InputTag> SimLabel {
        Name("SimLabel")
          };
      fhicl::Atom<art::InputTag> SimDepositLabel {
        Name("SimDepositLabel")
          };
      fhicl::Atom<art::InputTag> StripHitLabel {
        Name("StripHitLabel")
          };
      fhicl::Atom<art::InputTag> ClusterLabel {
        Name("ClusterLabel")
          };
      fhicl::Atom<art::InputTag> SpacePointLabel {
        Name("SpacePointLabel")
          };
      fhicl::Atom<art::InputTag> TrackLabel {
        Name("TrackLabel")
          };

      fhicl::Atom<bool> DrawTaggers {
        Name("DrawTaggers")
          };
      fhicl::Atom<bool> DrawModules {
        Name("DrawModules")
          };
      fhicl::Atom<bool> DrawFEBs {
        Name("DrawFEBs")
          };
      fhicl::Atom<bool> DrawFEBEnds {
        Name("DrawFEBEnds")
          };
      fhicl::Atom<bool> DrawStrips {
        Name("DrawStrips")
          };
      fhicl::Atom<bool> DrawTpc {
        Name("DrawTpc")
          };
      fhicl::Atom<bool> DrawTrueTracks {
        Name("DrawTrueTracks")
          };
      fhicl::Atom<bool> DrawSimDeposits {
        Name("DrawSimDeposits")
          };
      fhicl::Atom<bool> DrawStripHits {
        Name("DrawStripHits")
          };
      fhicl::Atom<bool> DrawClusters {
        Name("DrawClusters")
          };
      fhicl::Atom<bool> DrawSpacePoints {
        Name("DrawSpacePoints")
          };
      fhicl::Atom<bool> DrawTracks {
        Name("DrawTracks")
          };

      fhicl::Atom<bool> ChoseTaggers {
        Name("ChoseTaggers")
          };
      fhicl::Sequence<int> ChosenTaggers {
        Name("ChosenTaggers")
          };

      fhicl::Atom<int> TaggerColour {
        Name("TaggerColour")
          };
      fhicl::Atom<int> FEBColour {
        Name("FEBColour")
          };
      fhicl::Atom<int> FEBEndColour {
        Name("FEBEndColour")
          };
      fhicl::Atom<int> TpcColour {
        Name("TpcColour")
          };
      fhicl::Atom<int> TrueTrackColour {
        Name("TrueTrackColour")
          };
      fhicl::Atom<int> SimDepositColour {
        Name("SimDepositColour")
          };
      fhicl::Atom<int> StripHitColour {
        Name("StripHitColour")
          };
      fhicl::Atom<int> ClusterStartingColour {
        Name("ClusterStartingColour")
          };
      fhicl::Atom<int> ClusterColourInterval {
        Name("ClusterColourInterval")
          };
      fhicl::Atom<int> SpacePointColour {
        Name("SpacePointColour")
          };
      fhicl::Atom<int> TrackColour {
        Name("TrackColour")
          };

      fhicl::Atom<double> MinTime {
        Name ("MinTime"),
        Comment ("Ignore truth & reco products before this time"),
        0.
      };
      fhicl::Atom<double> MaxTime {
        Name ("MaxTime"),
        Comment ("Ignore truth & reco products after this time"),
        3.2e6
      };

      fhicl::Atom<bool> Print {
        Name("Print")
          };

      fhicl::Atom<double> LineWidth {
        Name("LineWidth")
          };
    };
    
    CRTEventDisplayAlg(const Config& config);
    
    CRTEventDisplayAlg(const fhicl::ParameterSet& pset) :
      CRTEventDisplayAlg(fhicl::Table<Config>(pset, {})()) {}
    
    CRTEventDisplayAlg();

    ~CRTEventDisplayAlg();

    void reconfigure(const Config& config);

    void SetDrawTaggers(bool tf);
    void SetDrawTpc(bool tf);
    void SetDrawTrueTracks(bool tf);
    void SetDrawSimDeposits(bool tf);
    void SetDrawStripHits(bool tf);
    void SetDrawClusters(bool tf);

    void SetPrint(bool tf);

    void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour);

    void Draw(detinfo::DetectorClocksData const& clockData, const art::Event& event);

    bool IsPointInsideBox(const std::vector<double> &lims, const geo::Point_t &p);

  private:
    
    TPCGeoAlg         fTPCGeoAlg;
    CRTGeoAlg         fCRTGeoAlg;
    CRTBackTrackerAlg fCRTBackTrackerAlg;

    art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
    
    art::InputTag fSimLabel;
    art::InputTag fSimDepositLabel;
    art::InputTag fStripHitLabel;
    art::InputTag fClusterLabel;
    art::InputTag fSpacePointLabel;
    art::InputTag fTrackLabel;

    bool fDrawTaggers;
    bool fDrawModules;
    bool fDrawFEBs;
    bool fDrawFEBEnds;
    bool fDrawStrips;
    bool fDrawTpc;
    bool fDrawTrueTracks;
    bool fDrawSimDeposits;
    bool fDrawStripHits;
    bool fDrawClusters;
    bool fDrawSpacePoints;
    bool fDrawTracks;

    bool             fChoseTaggers;
    std::vector<int> fChosenTaggers;

    int fTaggerColour;
    int fFEBColour;
    int fFEBEndColour;
    int fTpcColour;
    int fTrueTrackColour;
    int fSimDepositColour;
    int fStripHitColour;
    int fClusterStartingColour;
    int fClusterColourInterval;
    int fSpacePointColour;
    int fTrackColour;

    double fMinTime;
    double fMaxTime;

    bool fPrint;

    double fLineWidth;
  };
}

#endif
