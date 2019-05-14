#ifndef CRTEVENTDISPLAY_H_SEEN
#define CRTEVENTDISPLAY_H_SEEN


///////////////////////////////////////////////
// CRTEventDisplay.h
//
// Quick and dirty event display for SBND CRT
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// c++
#include <vector>

// ROOT
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TMarker.h"
#include "TBox.h"
#include "TPad.h"


namespace sbnd{

  class CRTEventDisplay {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SimLabel {
        Name("SimLabel")
      };
      fhicl::Atom<art::InputTag> CRTDataLabel {
        Name("CRTDataLabel")
      };
      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel")
      };
      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel")
      };
      fhicl::Atom<art::InputTag> TPCTrackLabel {
        Name("TPCTrackLabel")
      };

      fhicl::Atom<bool> DrawTaggers {
        Name("DrawTaggers")
      };
      fhicl::Atom<bool> DrawModules {
        Name("DrawModules")
      };
      fhicl::Atom<bool> DrawTpc {
        Name("DrawTpc")
      };
      fhicl::Atom<bool> DrawCrtHits {
        Name("DrawCrtHits")
      };
      fhicl::Atom<bool> DrawCrtTracks {
        Name("DrawCrtTracks")
      };
      fhicl::Atom<bool> DrawIncompleteTracks {
        Name("DrawIncompleteTracks")
      };
      fhicl::Atom<bool> DrawTpcTracks {
        Name("DrawTpcTracks")
      };
      fhicl::Atom<bool> DrawTrueTracks {
        Name("DrawTrueTracks")
      };

      fhicl::Atom<int> TaggerColour {
        Name("TaggerColour")
      };
      fhicl::Atom<int> TpcColour {
        Name("TpcColour")
      };
      fhicl::Atom<int> CrtHitColour {
        Name("CrtHitColour")
      };
      fhicl::Atom<int> CrtTrackColour {
        Name("CrtTrackColour")
      };
      fhicl::Atom<int> TpcTrackColour {
        Name("TpcTrackColour")
      };
      fhicl::Atom<int> TrueTrackColour {
        Name("TrueTrackColour")
      };

      fhicl::Atom<bool> UseTrueID {
        Name("UseTrueID")
      };
      fhicl::Atom<int> TrueID {
        Name("TrueID")
      };

      fhicl::Atom<bool> Print {
        Name("Print")
      };

      fhicl::Atom<double> LineWidth {
        Name("LineWidth")
      };
      fhicl::Atom<double> IncompleteTrackLength {
        Name("IncompleteTrackLength")
      };
      fhicl::Atom<double> MinTime {
        Name("MinTime")
      };
      fhicl::Atom<double> MaxTime {
        Name("MaxTime")
      };

    };

    CRTEventDisplay(const Config& config);

    CRTEventDisplay(const fhicl::ParameterSet& pset) :
      CRTEventDisplay(fhicl::Table<Config>(pset, {})()) {}

    CRTEventDisplay();

    ~CRTEventDisplay();

    void reconfigure(const Config& config);

    void SetDrawTaggers(bool tf);
    void SetDrawTpc(bool tf);
    void SetDrawCrtHits(bool tf);
    void SetDrawCrtTracks(bool tf);
    void SetDrawTpcTracks(bool tf);
    void SetDrawTrueTracks(bool tf);

    bool IsVisible(const simb::MCParticle& particle);

    void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour);

    void Draw(const art::Event& event);

  private:

    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;

    art::InputTag fSimLabel;
    art::InputTag fCRTDataLabel;
    art::InputTag fCRTHitLabel;
    art::InputTag fCRTTrackLabel;
    art::InputTag fTPCTrackLabel;

    bool fDrawTaggers;
    bool fDrawModules;
    bool fDrawTpc;
    bool fDrawCrtHits;
    bool fDrawCrtTracks;
    bool fDrawIncompleteTracks;
    bool fDrawTpcTracks;
    bool fDrawTrueTracks;

    int fTaggerColour;
    int fTpcColour;
    int fCrtHitColour;
    int fCrtTrackColour;
    int fTpcTrackColour;
    int fTrueTrackColour;

    bool fUseTrueID;
    int fTrueID;

    bool fPrint;

    double fLineWidth;
    double fIncompleteTrackLength;
    double fMinTime;
    double fMaxTime;

  };

}

#endif
