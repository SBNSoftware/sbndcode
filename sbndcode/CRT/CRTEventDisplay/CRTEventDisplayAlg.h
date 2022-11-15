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

// sbnobj
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// ROOT
#include "TPolyLine3D.h"
#include "TCanvas.h"

namespace detinfo { class DetectorClocksData; }

namespace sbnd
{
  
  class CRTEventDisplayAlg
  {
  public:
    
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
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

      fhicl::Atom<bool> DrawTaggers {
        Name("DrawTaggers")
          };
      fhicl::Atom<bool> DrawModules {
        Name("DrawModules")
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

      fhicl::Atom<int> TaggerColour {
        Name("TaggerColour")
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
      fhicl::Atom<int> ClusterColour {
        Name("ClusterColour")
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

  private:
    
    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;
    
    art::InputTag fSimLabel;
    art::InputTag fSimDepositLabel;
    art::InputTag fStripHitLabel;
    art::InputTag fClusterLabel;

    bool fDrawTaggers;
    bool fDrawModules;
    bool fDrawStrips;
    bool fDrawTpc;
    bool fDrawTrueTracks;
    bool fDrawSimDeposits;
    bool fDrawStripHits;
    bool fDrawClusters;

    int fTaggerColour;
    int fTpcColour;
    int fTrueTrackColour;
    int fSimDepositColour;
    int fStripHitColour;
    int fClusterColour;

    bool fPrint;

    double fLineWidth;
  };
}

#endif
