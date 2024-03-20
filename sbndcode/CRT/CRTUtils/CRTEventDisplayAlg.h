#ifndef CRTEVENTDISPLAYALG_H_SEEN
#define CRTEVENTDISPLAYALG_H_SEEN


///////////////////////////////////////////////
// CRTEventDisplayAlg.h
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
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
namespace detinfo { class DetectorClocksData; }

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

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
//#include "TAttLine.h"


namespace sbnd{

  class CRTEventDisplayAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SimLabel {
        Name("SimLabel")
      };
      fhicl::Atom<art::InputTag> AuxDetIDEsLabel {
        Name("AuxDetIDEsLabel")
      };
      fhicl::Atom<art::InputTag> AuxDetHitsLabel {
        Name("AuxDetHitsLabel")
      };
      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel")
      };
      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel")
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
      fhicl::Atom<bool> DrawStrips {
        Name("DrawStrips")
      };
      fhicl::Atom<bool> DrawAuxDetIDEs {
        Name("DrawAuxDetIDEs")
      };
      fhicl::Atom<bool> DrawAuxDetHits {
        Name("DrawAuxDetHits")
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
      fhicl::Atom<bool> DrawInvisibleTracks {
        Name("DrawInvisibleTracks")
      };
      fhicl::Atom<bool> DrawTrueTracks {
        Name("DrawTrueTracks")
      };
      fhicl::Atom<bool> DataMode {
        Name("DataMode")
      };

      fhicl::Atom<int> TaggerColour {
        Name("TaggerColour")
      };
      fhicl::Atom<int> StripColour {
        Name("StripColour")
      };
      fhicl::Atom<int> FEBEndColour {
        Name("FEBEndColour")
      };
      fhicl::Atom<int> AuxDetIDEsColour {
        Name("AuxDetIDEsColour")
      };
      fhicl::Atom<int> AuxDetHitsColour {
        Name("AuxDetHitsColour")
      };
      fhicl::Atom<int> CrtHitColour {
        Name("CrtHitColour")
      };
      fhicl::Atom<int> CrtTrackColour {
        Name("CrtTrackColour")
      };
      fhicl::Atom<int> TrueTrackColour {
        Name("TrueTrackColour")
      };

      fhicl::Atom<bool> UseTrueID {
        Name("UseTrueID")
      };
      fhicl::Sequence<int> TrueID{
        Name("TrueID"),
        {1}
      };

      fhicl::Atom<bool> Print {
        Name("Print")
      };
      fhicl::Atom<double> LineWidth {
        Name("LineWidth")
      };
      fhicl::Atom<double> VolumeSizeOffset {
        Name("VolumeSizeOffset")
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

      /*fhicl::Table<CRTBackTracker::Config> CrtBackTrack {
        Name("CrtBackTrack"),
      };*/

    };

    CRTEventDisplayAlg(const Config& config);

    CRTEventDisplayAlg(const fhicl::ParameterSet& pset) :
      CRTEventDisplayAlg(fhicl::Table<Config>(pset, {})()) {}

    CRTEventDisplayAlg();

    ~CRTEventDisplayAlg();

    void reconfigure(const Config& config);

    void SetDrawTaggers(bool tf);
    void SetDrawFEBs(bool tf);
    void SetDrawStrips(bool tf);
    void SetDrawAuxDetIDEs(bool tf);
    void SetDrawAuxDetHits(bool tf);
    void SetDrawCrtHits(bool tf);
    void SetDrawCrtTracks(bool tf);
    void SetDrawTrueTracks(bool tf);
    void SetPrint(bool tf);

    //void SetTrueId(int id);

    bool IsVisible(const simb::MCParticle& particle);

    void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour);

    void Draw(detinfo::DetectorClocksData const& clockData, const art::Event& event);

  private:
    CRTGeoAlg fCrtGeo;

    //CRTBackTracker fCrtBackTrack;

    art::InputTag fSimLabel;
    art::InputTag fAuxDetHitsLabel;
    art::InputTag fAuxDetIDEsLabel;
    art::InputTag fCRTHitLabel;
    art::InputTag fCRTTrackLabel;

    bool fDrawTaggers;
    bool fDrawModules;
    bool fDrawStrips;
    bool fDrawFEBs;
    bool fDrawAuxDetIDEs;
    bool fDrawAuxDetHits;
    bool fDrawCrtHits;
    bool fDrawCrtTracks;
    bool fDrawIncompleteTracks;
    bool fDrawInvisibleTracks;
    bool fDrawTrueTracks;
    bool fDataMode;

    int fTaggerColour;
    int fStripColour;
    int fFEBEndColour;
    int fAuxDetIDEsColour;
    int fAuxDetHitsColour;
    int fCrtHitColour;
    int fCrtTrackColour;
    int fTrueTrackColour;

    bool fUseTrueID;
    std::vector<int> fTrueID;

    bool fPrint;

    double fLineWidth;
    double fVolumeSizeOffset;
    double fIncompleteTrackLength;
    double fMinTime;
    double fMaxTime;

  };

}

#endif
