////////////////////////////////////////////////////////////////////////
// Class:       EtaPMTAna
// Plugin Type: analyzer
// File:        EtaPMTAna_module.hh
//
//
////////////////////////////////////////////////////////////////////////

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"


// ROOT and C++ includes
#include <TTree.h>
#include <string.h>


// Reco includes
// PDS
#include "lardataobj/RecoBase/OpFlash.h"
#include "larsim/Simulation/LArG4Parameters.h"


// CRT
#include "sbnobj/SBND/CRT/CRTTrack.hh"


// Geometry and mapping
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"



namespace opdet {
  class EtaPMTAna;

    enum SBNDPDSDetectorType {
    kPDUnknown = -1,   
    kPMTCoated = 0,       
    kPMTUncoated = 1,
    kXARAPUCAVUV,
    kXARAPUCAVIS
  };
}


class opdet::EtaPMTAna : public art::EDAnalyzer {
public:
  explicit EtaPMTAna(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  EtaPMTAna(EtaPMTAna const&) = delete;
  EtaPMTAna(EtaPMTAna&&) = delete;
  EtaPMTAna& operator=(EtaPMTAna const&) = delete;
  EtaPMTAna& operator=(EtaPMTAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  
  bool CRTTrackCrossesAV(const int , const sbnd::crt::CRTTrack& , TVector3& , TVector3&);
  double GetPMTRatio(std::vector<double> PE_v);
  void FillAssociatedTracks(const recob::OpFlash& flash, const std::vector<art::Ptr<sbnd::crt::CRTTrack>>& matched_tracks);

  // Data product input labels
  std::vector<std::string> fOpFlashesModuleLabel;
  std::string fCRTTrackModuleLabel;
  double fMaxTimeDiff;
  double fMinTrackLength;
  bool fMakePDSGeoTree;
  int fVerbosity;
  std::vector<double> fAVmin_tpc0;
  std::vector<double> fAVmax_tpc0;
  std::vector<double> fAVmin_tpc1;
  std::vector<double> fAVmax_tpc1;
  double fIsochronousWidth;

  // CRT save window
  std::vector<double> fCRTSaveWindow={-1500000, 1500000};

  // PDS mapping and geometry
  opdet::sbndPDMapAlg fPDSMap;
  geo::WireReadoutGeom const& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  std::set<int> fPDSBoxIDs;
  geo::WireReadoutGeom const& channelMapAlg = art::ServiceHandle<geo::WireReadout const>()->Get();
  static constexpr double fDefaultSimIDE = -999.;
  art::ServiceHandle<geo::Geometry> fGeom;


  // TTrees
  TTree *fPDSMapTree;
  TTree *fTree;

  unsigned int _eventID, _runID, _subrunID;

  // Saving OpFlash
  std::vector<double> _flash_time;
  std::vector<double> _flash_total_pe;
  std::vector<int> _flash_tpc;
  std::vector<double> _flash_pmt_ratio;

  // Saving CRT information
  std::vector<double> _tr_midpoint_x;
  std::vector<double> _tr_length;

  // PDS geo tree
  std::vector <int> _opDetID;
  std::vector <double> _opDetX;
  std::vector <double> _opDetY;
  std::vector <double> _opDetZ;
  std::vector <int> _opDetType;
};


DEFINE_ART_MODULE(opdet::EtaPMTAna)
