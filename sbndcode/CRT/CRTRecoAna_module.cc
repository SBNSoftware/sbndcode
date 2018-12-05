////////////////////////////////////////////////////////////////////////
// Class:       CRTRecoAna
// Module Type: analyzer
// File:        CRTRecoAna_module.cc
//
// Analysis module for evaluating CRT reconstruction on through going
// muons.
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"


// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TBox.h"
#include "TPad.h"
#include "TString.h"
#include "TGeoManager.h"
#include "Math/Vector3D.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TRandom.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace {
  // Local namespace for local functions
  // Declare here, define later

}

namespace sbnd {

  struct CRTStrip {
    double t0; 
    uint32_t channel; 
    double x; 
    double ex; 
    int id1; 
    int id2;
    std::pair<std::string,unsigned> tagger;
  };

  struct CRTHit {
    std::string tagger;
    double xpos;
    double ypos;
    double zpos;
    double xerr;
    double yerr;
    double zerr;
    double ts0_s;
    double ts0_ns;
    double ts1_ns;
    std::vector<int> ids;
  }; // CRTHit

  struct CRTTrack {
    double xstart;
    double ystart;
    double zstart;
    double xend;
    double yend;
    double zend;
    double time;
    std::vector<int> ids;
    bool complete;
  };

  struct RecoTruth{
    std::vector<std::pair<std::pair<std::string,unsigned>,art::Ptr<crt::CRTData>>> sipmHits;
    std::vector<CRTStrip> stripHits;
    std::vector<CRTHit> crtHits;
    std::vector<CRTHit> crtAveHits;
    std::vector<CRTTrack> crtTracks;
  };

  class CRTRecoAna : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of detector simulation data product")
      };

      fhicl::Atom<art::InputTag> CRTModuleLabel {
        Name("CRTModuleLabel"),
        Comment("tag of CRT simulation data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Atom<bool> VeryVerbose {
        Name("VeryVerbose"),
        Comment("Print detailed information about every track")
      };

      fhicl::Atom<bool> Plot {
        Name("Plot"),
        Comment("Plot tracks")
      };

      fhicl::Atom<int> PlotTrackID {
        Name("PlotTrackID"),
        Comment("ID of track to plot")
      };

      fhicl::Atom<double> TimeCoincidenceLimit {
        Name("TimeCoincidenceLimit"),
        Comment("Minimum time between two overlapping hit crt strips")
      };

      fhicl::Atom<double> QPed {
        Name("QPed"),
        Comment("Pedestal offset [ADC]")
      };

      fhicl::Atom<double> QSlope {
        Name("QSlope"),
        Comment("Pedestal slope [ADC/photon]")
      };

      fhicl::Atom<bool> UseReadoutWindow {
        Name("UseReadoutWindow"),
        Comment("Only reconstruct hits within readout window")
      };

      fhicl::Atom<double> TimeLimit {
        Name("TimeLimit"),
        Comment("Maximum time difference to combine hits into a tzero")
      };

      fhicl::Atom<double> AverageHitDistance {
        Name("AverageHitDistance"),
        Comment("Maximum distance to avarage hits on the same tagger over ")
      };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("Maximum distance from CRTTrack to absorb a new CRTHit")
      };

      fhicl::Atom<bool> UseTopPlane {
        Name("UseTopPlane"),
        Comment("Use hits from the top plane (SBND specific)")
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTRecoAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    std::vector<double> ChannelToLimits(CRTStrip strip);

    std::vector<double> CrtOverlap(std::vector<double> strip1, std::vector<double> strip2);
  
    void DrawTrueTracks(const std::vector<simb::MCParticle>& particle, std::map<int, RecoTruth> truthMatch, bool truth, bool strips, bool hits, bool tracks, bool tpc, int ind);

    std::pair<std::string,unsigned> ChannelToTagger(uint32_t channel);

    bool CheckModuleOverlap(uint32_t channel);

    bool CrossesTPC(CRTTrack track);
    
    //Function to average hits within a certain distance of each other
    std::vector<CRTHit> AverageHits(std::vector<CRTHit> hits, std::map<int, RecoTruth>& truthMatch);

    //Take a list of hits and find average parameters
    CRTHit DoAverage(std::vector<CRTHit> hits);

    //Create CRTTracks from list of hits
    std::vector<CRTTrack> CreateTracks(std::vector<CRTHit> hits, std::map<int, RecoTruth>& truthMatch);

    //Calculate the tagger crossing point of CRTTrack candidate
    TVector3 CrossPoint(CRTHit hit, TVector3 start, TVector3 diff);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTModuleLabel;      ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on
    bool          fVeryVerbose;         ///< print more information about what's going on
    bool          fPlot;                ///< plot tracks
    int           fPlotTrackID;         ///< id of track to plot
    double        fTimeCoincidenceLimit;///< minimum time between two overlapping hit crt strips [ticks]
    double        fQPed;                ///< Pedestal offset of SiPMs [ADC]
    double        fQSlope;              ///< Pedestal slope of SiPMs [ADC/photon]
    bool          fUseReadoutWindow;    ///< Only reconstruct hits within readout window
    double        fTimeLimit;           // Maximum time difference to combine hits into a tzero [us]
    double        fAverageHitDistance;  // Maximum distance to avarage hits on the same tagger over [cm]
    double        fDistanceLimit;       // Maximum distance from CRTTrack to absorb a new CRTHit [cm]
    bool          fUseTopPlane;         // Use hits from the top plane (SBND specific)
    
    // n-tuples
    TH2D* fTagXYZResolution[7];
    TH1D* fTagXResolution[7];
    TH1D* fTagYResolution[7];

    TH1D* fSipmHitsPerTrack;
    TH1D* fStripHitsPerTrack;
    TH1D* fCrtHitsPerTrack;
    TH1D* fAveHitsPerTrack;
    TH1D* fTracksPerTrack;

    TH1D* fTrackDist;
    TH1D* fMuCharge;
    TH1D* fOtherCharge;

    TH1D* fSipmDist;

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider
    detinfo::ElecClock fTrigClock;
    art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
    const geo::AuxDetGeometry* fAuxDetGeo;
    const geo::AuxDetGeometryCore* fAuxDetGeoCore;

    CRTTruthRecoAlg truthAlg;

    // Positions of the CRT planes
    std::vector<double> crtPlanes = {-359.1, -357.3, 357.3, 359.1, -358.9, -357.1, 661.52, 663.32, 865.52, 867.32, -240.65, -238.85, 655.35, 657.15};
    std::vector<int> fixCoord   = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2}; // Fixed coordinate for each plane
    std::vector<int> widthCoord = {2, 1, 2, 1, 0, 2, 2, 0, 2, 0, 1, 0, 1, 0}; // Width direction for each plane
    std::vector<int> lenCoord   = {1, 2, 1, 2, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1}; // Length direction for each plane

    const size_t nTaggers = 7;

    // Performance Counters
    int nSipms = 0;
    int nStrips = 0;
    int nHits = 0;

    // Truth Counters
    int nParts = 0;
    int nSinglePlane = 0;
    int nMatchTrack = 0;
    int nTwoStripPlanes = 0;

    int nMissedHits = 0;
    int nMatchingHits = 0;
    int nNoTruth = 0;
    int nMissTag[7] = {0,0,0,0,0,0,0};
    int nMatchTag[7] = {0,0,0,0,0,0,0};
    int nNoTruthTag[7] = {0,0,0,0,0,0,0};

    std::map<std::string, int> nameToInd;
    std::map<int, std::string> indToName;

  }; // class CRTRecoAna

  // Constructor
  CRTRecoAna::CRTRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTModuleLabel       (config().CRTModuleLabel())
    , fVerbose              (config().Verbose())
    , fVeryVerbose          (config().VeryVerbose())
    , fPlot                 (config().Plot())
    , fPlotTrackID          (config().PlotTrackID())
    , fTimeCoincidenceLimit (config().TimeCoincidenceLimit())
    , fQPed                 (config().QPed())
    , fQSlope               (config().QSlope())
    , fUseReadoutWindow     (config().UseReadoutWindow())
    , fTimeLimit            (config().TimeLimit())
    , fAverageHitDistance   (config().AverageHitDistance())
    , fDistanceLimit        (config().DistanceLimit())
    , fUseTopPlane          (config().UseTopPlane())
    , truthAlg()
  {
    // Get a pointer to the fGeometryServiceetry service provider
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 
    fTrigClock = fDetectorClocks->TriggerClock();
    fAuxDetGeo = &(*fAuxDetGeoService);
    fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();
  }

  void CRTRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    // Define n-tuple
    fTagXYZResolution[0]    = tfs->make<TH2D>("tag0XYZresolution",  ";True Z - reco Z (cm);True Y - reco Y (cm)",  50, -25, 25,  50, -25, 25);
    fTagXYZResolution[1]    = tfs->make<TH2D>("tag1XYZresolution",  ";True Z - reco Z (cm);True Y - reco Y (cm)",  50, -25, 25,  50, -25, 25);
    fTagXYZResolution[2]    = tfs->make<TH2D>("tag2XYZresolution",  ";True X - reco X (cm);True Z - reco Z (cm)",  50, -25, 25,  50, -25, 25);
    fTagXYZResolution[3]    = tfs->make<TH2D>("tag3XYZresolution",  ";True Z - reco Z (cm);True X - reco X (cm)",  50, -25, 25,  50, -25, 25);
    fTagXYZResolution[4]    = tfs->make<TH2D>("tag4XYZresolution",  ";True Z - reco Z (cm);True X - reco X (cm)",  50, -25, 25,  50, -25, 25);
    fTagXYZResolution[5]    = tfs->make<TH2D>("tag5XYZresolution",  ";True Y - reco Y (cm);True X - reco X (cm)",  50, -25, 25,  50, -25, 25);
    fTagXYZResolution[6]    = tfs->make<TH2D>("tag6XYZresolution",  ";True Y - reco Y (cm);True X - reco X (cm)",  50, -25, 25,  50, -25, 25);

    fTagXResolution[0]    = tfs->make<TH1D>("tag0Xresolution",  ";True Z - reco Z (cm);",  50, -25, 25);
    fTagXResolution[1]    = tfs->make<TH1D>("tag1Xresolution",  ";True Z - reco Z (cm);",  50, -25, 25);
    fTagXResolution[2]    = tfs->make<TH1D>("tag2Xresolution",  ";True X - reco X (cm);",  50, -25, 25);
    fTagXResolution[3]    = tfs->make<TH1D>("tag3Xresolution",  ";True Z - reco Z (cm);",  50, -25, 25);
    fTagXResolution[4]    = tfs->make<TH1D>("tag4Xresolution",  ";True Z - reco Z (cm);",  50, -25, 25);
    fTagXResolution[5]    = tfs->make<TH1D>("tag5Xresolution",  ";True Y - reco Y (cm);",  50, -25, 25);
    fTagXResolution[6]    = tfs->make<TH1D>("tag6Xresolution",  ";True Y - reco Y (cm);",  50, -25, 25);

    fTagYResolution[0]    = tfs->make<TH1D>("tag0Yresolution",  ";True Y - reco Y (cm);",  50, -25, 25);
    fTagYResolution[1]    = tfs->make<TH1D>("tag1Yresolution",  ";True Y - reco Y (cm);",  50, -25, 25);
    fTagYResolution[2]    = tfs->make<TH1D>("tag2Yresolution",  ";True Z - reco Z (cm);",  50, -25, 25);
    fTagYResolution[3]    = tfs->make<TH1D>("tag3Yresolution",  ";True X - reco X (cm);",  50, -25, 25);
    fTagYResolution[4]    = tfs->make<TH1D>("tag4Yresolution",  ";True X - reco X (cm);",  50, -25, 25);
    fTagYResolution[5]    = tfs->make<TH1D>("tag5Yresolution",  ";True X - reco X (cm);",  50, -25, 25);
    fTagYResolution[6]    = tfs->make<TH1D>("tag6Yresolution",  ";True X - reco X (cm);",  50, -25, 25);

    fSipmHitsPerTrack  = tfs->make<TH1D>("sipmhitspertrack",   ";SiPM Hits;N Tracks",  20, 0, 20);
    fStripHitsPerTrack = tfs->make<TH1D>("striphitspertrack",  ";Strip Hits;N Tracks",  10, 0, 10);
    fCrtHitsPerTrack   = tfs->make<TH1D>("crthitspertrack",    ";CRT Hits;N Tracks",  10, 0, 10);
    fAveHitsPerTrack   = tfs->make<TH1D>("avehitspertrack",    ";Average Hits;N Tracks",  10, 0, 10);
    fTracksPerTrack    = tfs->make<TH1D>("trackshitspertrack", ";CRT Tracks;N Tracks",  10, 0, 10);

    fTrackDist    = tfs->make<TH1D>("trackdist", ";True reco dist^2 (cm^2);N Tracks",  50, 0, 200);
    fMuCharge    = tfs->make<TH1D>("mucharge", ";Muon charge (ADCs);N Tracks",  100, 0, 2000);
    fOtherCharge    = tfs->make<TH1D>("othercharge", ";Other charge (ADCs);N Tracks",  100, 0, 2000);
    fSipmDist    = tfs->make<TH1D>("sipmdist", ";Distance from SiPM (cm);N Tracks",  50, -3, 13);
    // Initial output
    std::cout<<"----------------- CRT Hit Reco Module -------------------"<<std::endl;

    // Take a position that is known to be the center of a strip

    nameToInd["volTaggerSideRight_0"] = 0;
    nameToInd["volTaggerSideLeft_0"] = 1;
    nameToInd["volTaggerBot_0"] = 2;
    nameToInd["volTaggerTopLow_0"] = 3;
    nameToInd["volTaggerTopHigh_0"] = 4;
    nameToInd["volTaggerFaceFront_0"] = 5;
    nameToInd["volTaggerFaceBack_0"] = 6;

    indToName[0] = "volTaggerSideRight_0";
    indToName[1] = "volTaggerSideLeft_0";
    indToName[2] = "volTaggerBot_0";
    indToName[3] = "volTaggerTopLow_0";
    indToName[4] = "volTaggerTopHigh_0";
    indToName[5] = "volTaggerFaceFront_0";
    indToName[6] = "volTaggerFaceBack_0";

  }// CRTRecoAna::beginJob()

  void CRTRecoAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Detector properties
    double readoutWindowMuS  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    double driftTimeMuS = (2.*fGeometryService->DetHalfWidth()+3.)/fDetectorProperties->DriftVelocity(); // [us]

    // Store the true x (CRT width direction) and y (CRT length direction) crossing points
    std::vector<TVector3> *taggerXYZ = new std::vector<TVector3>[nTaggers];
    std::map<int,TVector3> *partXYZ = new std::map<int,TVector3>[nTaggers];
    std::map<int,TVector3> *usedXYZ = new std::map<int,TVector3>[nTaggers];

    std::vector<int> muonIDs;

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // Fill a map of true particles
    std::vector<simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      double pt = particle.T() * 1e-3; // [us]

      // Check particle is a muon
      int pdg = std::abs(particle.PdgCode());
      if(pdg != 13){ muonIDs.push_back(partId); continue;}

      // Check if the particle is in the reconstructible time window
      if(fUseReadoutWindow){
        if(pt < -driftTimeMuS || pt > readoutWindowMuS) continue;
      }

      // Check particle is through-going
      if(!truthAlg.IsThroughGoing(particle)) continue;
      int ntag = 0;
      int nstrp = 0;
      // Loop over number of taggers
      for (size_t tag_i = 0; tag_i < nTaggers; tag_i++){
        if(!truthAlg.CrossesTagger(particle, tag_i)) continue;
        if(truthAlg.CrossesStrip(particle, tag_i) && tag_i!=4) nstrp++;
        TVector3 crossPoint = truthAlg.TaggerCrossPoint(particle, tag_i);
        taggerXYZ[tag_i].push_back(crossPoint);
        partXYZ[tag_i][partId] = crossPoint;
        ntag++;
      }

      particles.push_back(particle);
      nParts++;

      if(ntag==1) nSinglePlane++;
      if(nstrp>1) nTwoStripPlanes++;
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    // Retrieve list of CRT hits
    art::Handle< std::vector<crt::CRTData>> crtListHandle;
    std::vector<art::Ptr<crt::CRTData> > crtList;
    if (event.getByLabel(fCRTModuleLabel,crtListHandle))
      art::fill_ptr_vector(crtList, crtListHandle); 

    if(fVerbose) std::cout<<"Number of SiPM hits = "<<crtList.size()<<std::endl;

    // Fill a vector of pairs of time and width direction for each CRT plane
    // The y crossing point of z planes and z crossing point of y planes would be constant
    std::map<int, RecoTruth> truthMatch;
    std::map<std::pair<std::string,unsigned>, std::vector<CRTStrip>> taggerStrips;
    
    // Loop over all the SiPM hits in 2 (should be in pairs due to trigger)
    for (size_t i = 0; i < crtList.size(); i+=2){

      // Get the time, channel, center and width
      fTrigClock.SetTime(crtList[i]->T0());
      double t1 = fTrigClock.Time(); // [us]
      uint32_t channel = crtList[i]->Channel();
      int strip = (channel >> 1) & 15;
      int module = (channel >> 5);
      std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
      TVector3 center = fAuxDetGeoCore->AuxDetChannelToPosition(2*strip, name);
      const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);
      double width = 2*stripGeo.HalfWidth1();

      std::pair<std::string,unsigned> tagger = ChannelToTagger(channel);
      //std::cout<<"Tagger = "<<tagger.first<<" plane = "<<tagger.second<<" channel = "<<channel<<std::endl;

      // Check hit is inside the reconstructible window
      if(fUseReadoutWindow){
        if (!(t1 >= -driftTimeMuS && t1 <= readoutWindowMuS)) continue;
      }
      nSipms += 2;

      // Sort the hit SiPMs by what plane they're in
      nStrips++;

      int id1 = std::abs(crtList[i]->TrackID());
      int id2 = std::abs(crtList[i+1]->TrackID()); 
      truthMatch[id1].sipmHits.push_back(std::make_pair(tagger, crtList[i]));
      truthMatch[id2].sipmHits.push_back(std::make_pair(tagger, crtList[i+1]));

      // Get the time of hit on the second SiPM
      fTrigClock.SetTime(crtList[i+1]->T0());
      double t2 = fTrigClock.Time(); // [us]

      // Calculate the number of photoelectrons at each SiPM
      double npe1 = ((double)crtList[i]->ADC() - fQPed)/fQSlope;
      double npe2 = ((double)crtList[i+1]->ADC() - fQPed)/fQSlope;

      // Calculate the distance between the SiPMs
      double x = (width/2.)*atan(log(1.*npe2/npe1)) + (width/2.);
      if(tagger.first=="volTaggerTopLow_0"||tagger.first=="volTaggerTopHigh_0") fSipmDist->Fill(x);

      if(id1==id2 && tagger.first=="volTaggerTopLow_0"){
        if(std::find(muonIDs.begin(), muonIDs.end(), id1) != muonIDs.end()){
          fMuCharge->Fill(crtList[i]->ADC());
        }
        else {
          fOtherCharge->Fill(crtList[i]->ADC());
        }
      }

      // Calculate the error
      double normx = x + 0.344677*x - 1.92045;
      double ex = 1.92380e+00+1.47186e-02*normx-5.29446e-03*normx*normx;
      double time = (t1 + t2)/2.;

      CRTStrip stripHit = {time, channel, x, ex, id1, id2, tagger};
      taggerStrips[tagger].push_back(stripHit);

      truthMatch[id1].stripHits.push_back(stripHit);
      if(id1!=id2) truthMatch[id2].stripHits.push_back(stripHit);

    }

    // Remove any duplicate (same channel and time) hit strips
    for(auto &tagStrip : taggerStrips){
      std::sort(tagStrip.second.begin(), tagStrip.second.end(),
                [](const CRTStrip & a, const CRTStrip & b) -> bool{
                  return (a.t0 < b.t0) || 
                         ((a.t0 == b.t0) && (a.channel < b.channel));
                });
      // Remove hits with the same time and channel
      tagStrip.second.erase(std::unique(tagStrip.second.begin(), tagStrip.second.end(),
                                           [](const CRTStrip & a, const CRTStrip & b) -> bool{
                                             return a.t0 == b.t0 && a.channel == b.channel;
                                            }), tagStrip.second.end());
    }

    std::vector<CRTHit> crtHits;
    std::vector<std::string> usedTaggers;

    // TODO: Save any strips with no overlaps as hits (may increase noise hits too much)
    for (auto &tagStrip : taggerStrips){

      if (std::find(usedTaggers.begin(),usedTaggers.end(),tagStrip.first.first)!=usedTaggers.end()) continue;

      usedTaggers.push_back(tagStrip.first.first);

      unsigned planeID = 0;
      if(tagStrip.first.second==0) planeID = 1;
      std::pair<std::string,unsigned> otherPlane = std::make_pair(tagStrip.first.first, planeID);

      for (size_t hit_i = 0; hit_i < tagStrip.second.size(); hit_i++){
        // Get the position (in real space) of the 4 corners of the hit, taking charge sharing into account
        std::vector<double> limits1 =  ChannelToLimits(tagStrip.second[hit_i]);

        // Check for overlaps on the first plane
        if(CheckModuleOverlap(tagStrip.second[hit_i].channel)){

          // Loop over all the hits on the parallel (odd) plane
          for (size_t hit_j = 0; hit_j < taggerStrips[otherPlane].size(); hit_j++){

            // Get the limits in the two variable directions
            std::vector<double> limits2 = ChannelToLimits(taggerStrips[otherPlane][hit_j]);
            // If the time and position match then record the pair of hits
            std::vector<double> overlap = CrtOverlap(limits1, limits2);
            double t0_1 = tagStrip.second[hit_i].t0;
            double t0_2 = taggerStrips[otherPlane][hit_j].t0;

            if (overlap[0] != -99999 && std::abs(t0_1 - t0_2) < fTimeCoincidenceLimit){
              nHits++;

              // Calculate the mean and error in x, y, z
              TVector3 mean((overlap[0] + overlap[1])/2., 
                            (overlap[2] + overlap[3])/2., 
                            (overlap[4] + overlap[5])/2.);
              TVector3 error(std::abs((overlap[1] - overlap[0])/2.), 
                             std::abs((overlap[3] - overlap[2])/2.), 
                             std::abs((overlap[5] - overlap[4])/2.));

              // Average the time
              double time = (t0_1 + t0_2)/2; // [us]

              // If the PID matches one of the true crossing particles calculate the resolution
              std::vector<int> ids = {tagStrip.second[hit_i].id1, 
                                      tagStrip.second[hit_i].id2, 
                                      taggerStrips[otherPlane][hit_j].id1, 
                                      taggerStrips[otherPlane][hit_j].id2};

              std::sort(ids.begin(), ids.end());
              ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

              // Create a CRT hit
              CRTHit crtHit = {tagStrip.first.first, mean.X(), mean.Y(), mean.Z(), 
                               error.X(), error.Y(), error.Z(), 
                               time*1e-6, time*1e3, time*1e3, ids};
              crtHits.push_back(crtHit);

              int tag_i = nameToInd[tagStrip.first.first];
              for(int& ID : ids){
                if(partXYZ[tag_i].find(ID)!=partXYZ[tag_i].end()){
                  nMatchingHits++; 
                  nMatchTag[tag_i]++;
                  usedXYZ[tag_i][ID] = partXYZ[tag_i][ID];
                  TVector3 dist = partXYZ[tag_i][ID] - mean;
                  double distx = dist[widthCoord[tag_i*2]];
                  double disty = dist[lenCoord[tag_i*2]];
                  fTagXYZResolution[tag_i]->Fill(distx,disty); 
                  fTagXResolution[tag_i]->Fill(distx); 
                  fTagYResolution[tag_i]->Fill(disty);
                }
                else{ nNoTruth++; nNoTruthTag[tag_i]++; }
                truthMatch[ID].crtHits.push_back(crtHit);
              }

            }
          }
        }
        else{
          TVector3 mean((limits1[0] + limits1[1])/2., 
                        (limits1[2] + limits1[3])/2., 
                        (limits1[4] + limits1[5])/2.);
          TVector3 error(std::abs((limits1[1] - limits1[0])/2.), 
                         std::abs((limits1[3] - limits1[2])/2.), 
                         std::abs((limits1[5] - limits1[4])/2.));

          nHits++;

          double time = tagStrip.second[hit_i].t0; // [us]

          // If the PID matches calculate the resolution
          std::vector<int> ids = {tagStrip.second[hit_i].id1, tagStrip.second[hit_i].id2};

          std::sort(ids.begin(), ids.end());
          ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

          // Just use the single plane limits as the crt hit
          CRTHit crtHit = {tagStrip.first.first, mean.X(), mean.Y(), mean.Z(), 
                           error.X(), error.Y(), error.Z(), 
                           time*1e-6, time*1e3, time*1e3, ids};
          crtHits.push_back(crtHit);

          int tag_i = nameToInd[tagStrip.first.first];
          for(int& ID : ids){
            if(partXYZ[tag_i].find(ID)!=partXYZ[tag_i].end()){
              nMatchingHits++; 
              nMatchTag[tag_i]++;
              usedXYZ[tag_i][ID] = partXYZ[tag_i][ID];
              TVector3 dist = partXYZ[tag_i][ID] - mean;
              double distx = dist[widthCoord[tag_i*2]];
              double disty = dist[lenCoord[tag_i*2]];
              fTagXYZResolution[tag_i]->Fill(distx,disty); 
              fTagXResolution[tag_i]->Fill(distx); 
              fTagYResolution[tag_i]->Fill(disty);
            }
            else{ nNoTruth++; nNoTruthTag[tag_i]++; }
            truthMatch[ID].crtHits.push_back(crtHit);
          }

        }
      }
      for (size_t hit_j = 0; hit_j < taggerStrips[otherPlane].size(); hit_j++){
        // Get the limits in the two variable directions
        std::vector<double> limits1 = ChannelToLimits(taggerStrips[otherPlane][hit_j]);

        if(!CheckModuleOverlap(taggerStrips[otherPlane][hit_j].channel)){

          TVector3 mean((limits1[0] + limits1[1])/2., 
                        (limits1[2] + limits1[3])/2., 
                        (limits1[4] + limits1[5])/2.);
          TVector3 error(std::abs((limits1[1] - limits1[0])/2.), 
                         std::abs((limits1[3] - limits1[2])/2.), 
                         std::abs((limits1[5] - limits1[4])/2.));

          nHits++;
          double time = taggerStrips[otherPlane][hit_j].t0; // [us]

          // If the PID matches calculate the resolution
          std::vector<int> ids = {taggerStrips[otherPlane][hit_j].id1, taggerStrips[otherPlane][hit_j].id2};
          std::sort(ids.begin(), ids.end());
          ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

          // Just use the single plane limits as the crt hit
          CRTHit crtHit = {otherPlane.first, mean.X(), mean.Y(), mean.Z(), 
                           error.X(), error.Y(), error.Z(), 
                           time*1e-6, time*1e3, time*1e3, ids};
          crtHits.push_back(crtHit);

          int tag_i = nameToInd[tagStrip.first.first];
          for(int& ID : ids){
            if(partXYZ[tag_i].find(ID)!=partXYZ[tag_i].end()){
              nMatchingHits++; 
              nMatchTag[tag_i]++;
              usedXYZ[tag_i][ID] = partXYZ[tag_i][ID];
              TVector3 dist = partXYZ[tag_i][ID] - mean;
              double distx = dist[widthCoord[tag_i*2]];
              double disty = dist[lenCoord[tag_i*2]];
              fTagXYZResolution[tag_i]->Fill(distx,disty); 
              fTagXResolution[tag_i]->Fill(distx); 
              fTagYResolution[tag_i]->Fill(disty);
            }
            else{ nNoTruth++; nNoTruthTag[tag_i]++; }
            truthMatch[ID].crtHits.push_back(crtHit);
          }

        }
      }
    }

    // Loop over taggers
    for (size_t tag_i = 0; tag_i < nTaggers; tag_i++){
      // Loop over particle XYZ map
      for (auto & part : partXYZ[tag_i]) {
        // Find id in used list
        // If it's not there then count it
        if(usedXYZ[tag_i].find(part.first)==usedXYZ[tag_i].end()){
          nMissedHits++;
          nMissTag[tag_i]++;
        }
      }
    }

    // Sort CRTHits by time
    std::sort(crtHits.begin(), crtHits.end(), [](auto& left, auto& right)->bool{
              return left.ts1_ns < right.ts1_ns;});

    std::vector<std::vector<CRTHit>> CRTTzeroVect;
    std::vector<int> npvec;
    int iflag[2000] = {};

    // Loop over crt hits
    for(size_t i = 0; i<crtHits.size(); i++){

      if (iflag[i] == 0){
        std::vector<CRTHit> CRTTzero;
        std::map<std::string,int> nPlanes;

        double time_ns_A = crtHits[i].ts1_ns;

        iflag[i] = 1;
        CRTTzero.push_back(crtHits[i]);
        nPlanes[crtHits[i].tagger]++;

        // Get the t0
        // Loop over all the other CRT hits
        for(size_t j = i+1; j<crtHits.size(); j++){

          if(iflag[j] == 0){
            // If ts1_ns - ts1_ns < diff then put them in a vector
            double time_ns_B = crtHits[j].ts1_ns;
            double diff = std::abs(time_ns_B - time_ns_A) * 1e-3; // [us]

            if(diff < fTimeLimit){
              iflag[j]=1;
              CRTTzero.push_back(crtHits[j]);
              nPlanes[crtHits[j].tagger]++;
            }
          }
        }

        int np = 0;
        for(auto &nPlane : nPlanes){
          if(nPlane.second>0 && nPlane.first!="volTaggerTopHigh_0") np++;
        }

        CRTTzeroVect.push_back(CRTTzero);
        npvec.push_back(np);
      }
    }

    int nTracks = 0;
    int nCompTracks = 0;
    int nIncTracks = 0;

    // Loop over tzeros
    for(size_t i = 0; i<CRTTzeroVect.size(); i++){

      //loop over hits for this tzero, sort by tagger
      std::map<std::string, std::vector<CRTHit>> hits;
      for (size_t ah = 0; ah< CRTTzeroVect[i].size(); ++ah){        
        std::string ip = CRTTzeroVect[i][ah].tagger;       
        hits[ip].push_back(CRTTzeroVect[i][ah]);
      } // loop over hits
      
      //loop over planes and calculate average hits
      std::vector<CRTHit> allHits;
      for (auto &keyVal : hits){
        std::string ip = keyVal.first;
        std::vector<CRTHit> ahits = AverageHits(hits[ip], truthMatch);
        if(fUseTopPlane && ip == "volTaggerTopHigh_0"){ 
          allHits.insert(allHits.end(), ahits.begin(), ahits.end());
        }
        else if(ip != "volTaggerTopHigh_0"){ 
          allHits.insert(allHits.end(), ahits.begin(), ahits.end());
        }
      }

      std::vector<CRTTrack> rTracks = CreateTracks(allHits, truthMatch);
      nTracks += rTracks.size();
      for(size_t j = 0; j < rTracks.size(); j++){
        if(rTracks[j].complete) nCompTracks++;
        else nIncTracks++;
      }
    }

    if(fVerbose){
      std::cout<<"Number of hits = "<<crtHits.size()<<std::endl
               <<"Number of T zero = "<<CRTTzeroVect.size()<<std::endl
               <<"Number of tracks = "<<nTracks<<std::endl
               <<"Number of complete tracks = "<<nCompTracks<<std::endl
               <<"Number of incomplete tracks = "<<nIncTracks<<std::endl;
    }

    if(fPlot) DrawTrueTracks(particles, truthMatch, true, false, false, true, false, fPlotTrackID);

    for(auto const& particle : particles){
      int partId = particle.TrackId();

      RecoTruth rt = truthMatch[partId];
      fSipmHitsPerTrack->Fill(rt.sipmHits.size());
      fStripHitsPerTrack->Fill(rt.stripHits.size());
      fCrtHitsPerTrack->Fill(rt.crtHits.size());
      fAveHitsPerTrack->Fill(rt.crtAveHits.size());
      fTracksPerTrack->Fill(rt.crtTracks.size());

      if(rt.crtTracks.size()>0) nMatchTrack++;

      if(rt.crtTracks.size()>0){
        double minxres = 999999;
        double minipt = 1;

        for(size_t i = 0; i < rt.crtTracks.size(); i++){
          double sx = rt.crtTracks[i].xstart;
          double sy = rt.crtTracks[i].ystart;
          double sz = rt.crtTracks[i].zstart;
          double ex = rt.crtTracks[i].xend;
          double ey = rt.crtTracks[i].yend;
          double ez = rt.crtTracks[i].zend;
          double p0 = sx - (sz*(sx-ex)/(sz-ez));
          double p1 = (sx-ex)/(sz-ez);
          double p2 = sy - (sz*(sy-ey)/(sz-ez));
          double p3 = (sy-ey)/(sz-ez);
          ROOT::Math::XYZVector x0(p0, p2, 0.);
          ROOT::Math::XYZVector x1(p0+p1, p2+p3, 1.);
          ROOT::Math::XYZVector u = (x1-x0).Unit();

          int nTraj = particle.NumberTrajectoryPoints();
          int ipt = 0;
          double xres = 0;
          for(int j = 0; j < nTraj; j++){
            //if(abs(particle.Vx(j))<200 && particle.Vy(j)<200 && particle.Vy(j)>-200 && particle.Vz(j)>0 && particle.Vz(j)<500){
            if (particle.Vx(j) >= crtPlanes[0] && particle.Vx(j) <= crtPlanes[3] && particle.Vy(j) >= crtPlanes[4] && particle.Vy(j) <= -crtPlanes[4] && particle.Vz(j) >= crtPlanes[10] && particle.Vz(j) <= crtPlanes[13]){
              ROOT::Math::XYZVector xp(particle.Vx(j), particle.Vy(j), particle.Vz(j));
              double d2 = ((xp-x0).Cross(u)).Mag2();
              xres += d2;
              ipt++;
            }
          }

          if(xres < minxres) {minxres = xres; minipt = ipt;}
        }

        fTrackDist->Fill(minxres/minipt);
      }
    }

    // Detailed truth matching information for each true particle
    if(fVeryVerbose){
      for(auto const& particle : particles){
        int partId = particle.TrackId();
        //Print the crossing points
        std::cout<<"\nParticle "<<partId<<": Start("<<particle.Vx()<<","<<particle.Vy()<<","<<particle.Vz()<<") End("<<particle.EndX()<<","<<particle.EndY()<<","<<particle.EndZ()<<") time = "<<particle.T()*1e-3<<" us \nTrue crossing points:\n"; //FIXME to us
        for(size_t tag_i = 0; tag_i < nTaggers; tag_i++){
          if(partXYZ[tag_i].find(partId)!=partXYZ[tag_i].end()){
            std::cout<<"Tagger "<<tag_i<<": Coordinates = ("<<partXYZ[tag_i][partId].X()<<", "<<partXYZ[tag_i][partId].Y()<<", "<<partXYZ[tag_i][partId].Z()<<")\n";
          }
        }

        RecoTruth rt = truthMatch[partId];
        std::cout<<"->Reco tracks:\n";
        for(size_t trk_i = 0; trk_i < rt.crtTracks.size(); trk_i++){
          CRTTrack tr = rt.crtTracks[trk_i];
          std::cout<<"  Start = ("<<tr.xstart<<","<<tr.ystart<<","<<tr.zstart<<") End = ("<<tr.xend<<","<<tr.yend<<","<<tr.zend<<") time = "<<tr.time*1e-3<<" us, ids = ";
          for(int id : tr.ids) std::cout<<id<<" ";
          std::cout<<"\n";
        }

        std::cout<<"-->Average hits:\n";
        for(size_t hit_i = 0; hit_i < rt.crtAveHits.size(); hit_i++){
          CRTHit ah = rt.crtAveHits[hit_i];
          std::cout<<"   "<<ah.tagger<<": Coordinates = ("<<ah.xpos<<","<<ah.ypos<<","<<ah.zpos<<") time = "<<ah.ts1_ns*1e-3<<" us, ids = ";
          for(int id : ah.ids) std::cout<<id<<" ";
          std::cout<<"\n";
        }

        std::cout<<"--->Reco crossing points:\n";
        for(size_t hit_i = 0; hit_i < rt.crtHits.size(); hit_i++){
          CRTHit ht = rt.crtHits[hit_i];
          std::cout<<"    "<<ht.tagger<<": Coordinates = ("<<ht.xpos<<", "<<ht.ypos<<", "<<ht.zpos<<") time = "<<ht.ts1_ns*1e-3<<" us, ids = ";
          for(int id : ht.ids){ std::cout<<id<<" ";}
          std::cout<<"\n";
        }

        std::cout<<"---->Strip hits:\n";
        for(size_t hit_i = 0; hit_i < rt.stripHits.size(); hit_i++){
          CRTStrip sp = rt.stripHits[hit_i];
          std::cout<<"     "<<sp.tagger.first<<" ("<<sp.tagger.second<<"): time = "<<sp.t0<<" id1 = "<<sp.id1<<" id2 = "<<sp.id2<<"\n";
        }

        std::cout<<"----->SiPM hits:\n";
        for(size_t hit_i = 0; hit_i < rt.sipmHits.size(); hit_i++){
          std::pair<std::string,unsigned> tagger = rt.sipmHits[hit_i].first;
          art::Ptr<crt::CRTData> si = rt.sipmHits[hit_i].second;
          fTrigClock.SetTime(si->T0());
          double t1 = fTrigClock.Time(); // [us]
          std::cout<<"      "<<tagger.first<<" ("<<tagger.second<<"): Channel = "<<si->Channel()<<" time = "<<t1<<" us, ID = "<<si->TrackID()<<"\n";
        }
      }
    }

    delete[] taggerXYZ;
    delete[] partXYZ;
    delete[] usedXYZ;

  } // CRTRecoAna::analyze()

  void CRTRecoAna::endJob(){

    TF1 *fx[7];
    TF1 *fy[7];
    for(int i = 0; i < 7; i++){
      TString fxname = Form("f%ix", i);
      TString fyname = Form("f%iy", i);
      fx[i] = new TF1(fxname,"gaus");
      fTagXResolution[i]->Fit(fxname,"Q");
      fy[i] = new TF1(fyname,"gaus");
      fTagYResolution[i]->Fit(fyname,"Q");
    }

    std::string Xs[7] = {"Z","Z","X","Z","Z","Y","Y"};
    std::string Ys[7] = {"Y","Y","Z","X","X","X","X"};

    std::cout<<"=========================== GENERAL INFORMATION ===========================\n"
             <<"Number of through-going muons = "<<nParts
             <<"\nNumber of muons that pass through 1 plane = "<<nSinglePlane
             <<"\nNumber of muons with at least 1 matched track = "<<nMatchTrack
             <<"\nNumber of muons that cross 2 strips on 2 planes = "<<nTwoStripPlanes
             <<"\nNumber of hit SiPMs = "<<nSipms
             <<"\nNumber of hit strips = "<<nStrips
             <<"\nNumber of CRT hits = "<<nHits
             <<"\nTotal number of matched hits = "<<nMatchingHits
             <<"\nTotal number of missed hits = "<<nMissedHits
             <<"\nTotal number of hits with no truth = "<<nNoTruth<<"\n\n"
             <<"=========================== TAGGER INFORMATION ============================\n";
    for(int i = 0; i < 7; i++){
    std::cout<<"---->"<<indToName[i]<<":\n"
             <<Xs[i]<<" resolution = "<<fx[i]->GetParameter(2)<<" +/- "<<fx[i]->GetParError(2)<<" cm, bias = "<<fx[i]->GetParameter(1)<<" +/- "<<fx[i]->GetParError(1)<<" cm\n"
             <<Ys[i]<<" resolution = "<<fy[i]->GetParameter(2)<<" +/- "<<fy[i]->GetParError(2)<<" cm, bias = "<<fy[i]->GetParameter(1)<<" +/- "<<fy[i]->GetParError(1)<<" cm\n"
             <<"Efficiency = "<<(double)nMatchTag[i]/(nMatchTag[i]+nMissTag[i])<<"\n"
             <<nMatchTag[i]<<" matched hits, "<<nMissTag[i]<<" missed hits, "<<nNoTruthTag[i]<<" hits with no truth info\n\n";
    }

  } // CRTRecoAna::endJob()


  // Function to calculate the strip position limits in real space from channel
  std::vector<double> CRTRecoAna::ChannelToLimits(CRTStrip stripHit){
    int strip = (stripHit.channel >> 1) & 15;
    int module = (stripHit.channel >> 5);
    std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
    const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);
    double halfWidth = stripGeo.HalfWidth1();
    double halfHeight = stripGeo.HalfHeight();
    double halfLength = stripGeo.HalfLength();
    double l1[3] = {-halfWidth+stripHit.x+stripHit.ex, halfHeight, halfLength};
    double w1[3] = {0,0,0};
    stripGeo.LocalToWorld(l1, w1);
    double l2[3] = {-halfWidth+stripHit.x-stripHit.ex, -halfHeight, -halfLength};
    double w2[3] = {0,0,0};
    stripGeo.LocalToWorld(l2, w2);
    // Use this to get the limits in the two variable directions
    std::vector<double> limits = {std::min(w1[0],w2[0]), std::max(w1[0],w2[0]), 
                                  std::min(w1[1],w2[1]), std::max(w1[1],w2[1]), 
                                  std::min(w1[2],w2[2]), std::max(w1[2],w2[2])};
    return limits;
  } // CRTRecoAna::ChannelToLimits


  // Function to calculate the overlap between two crt strips
  std::vector<double> CRTRecoAna::CrtOverlap(std::vector<double> strip1, std::vector<double> strip2){
    double minX = std::max(strip1[0], strip2[0]);
    double maxX = std::min(strip1[1], strip2[1]);
    double minY = std::max(strip1[2], strip2[2]);
    double maxY = std::min(strip1[3], strip2[3]);
    double minZ = std::max(strip1[4], strip2[4]);
    double maxZ = std::min(strip1[5], strip2[5]);
    std::vector<double> null = {-99999, -99999, -99999, -99999, -99999, -99999};
    std::vector<double> overlap = {minX, maxX, minY, maxY, minZ, maxZ};
    if ((minX<maxX && minY<maxY) || (minX<maxX && minZ<maxZ) || (minY<maxY && minZ<maxZ)) return overlap;
    return null;
  } // CRTRecoAna::CRTOverlap()
    
  void CRTRecoAna::DrawTrueTracks(const std::vector<simb::MCParticle>& particle, std::map<int, RecoTruth> truthMatch, bool truth, bool strips, bool hits, bool tracks, bool tpc, int ind){
    // Create a canvas 
    TCanvas *c1 = new TCanvas("c1","",700,700);
    // Draw the tagger planes
    for(int tag_i = 0; tag_i < 7; tag_i++){
      double tagCenter[3] = {0, 0, 208.25};
      tagCenter[fixCoord[tag_i*2]] = (crtPlanes[tag_i*2]+crtPlanes[tag_i*2+1])/2;
      double tagDim[3] = {0, 0, 0};
      if(tag_i==0 || tag_i==1){ tagDim[0] = 1.8; tagDim[1] = 360; tagDim[2] = 450; }
      if(tag_i==2){ tagDim[0] = 399.5; tagDim[1] = 1.8; tagDim[2] = 478; }
      if(tag_i==3 || tag_i==4){ tagDim[0] = 450; tagDim[1] = 1.8; tagDim[2] = 450; }
      if(tag_i==5 || tag_i==6){ tagDim[0] = 360; tagDim[1] = 360; tagDim[2] = 1.8; }
      double rmin[3] = {tagCenter[0]-tagDim[0],tagCenter[1]-tagDim[1],tagCenter[2]-tagDim[2]};
      double rmax[3] = {tagCenter[0]+tagDim[0],tagCenter[1]+tagDim[1],tagCenter[2]+tagDim[2]};
      truthAlg.DrawCube(c1, rmin, rmax, 1);
    }

    if(tpc){
      double xmin = -2.0 * fGeometryService->DetHalfWidth();
      double xmax = 2.0 * fGeometryService->DetHalfWidth();
      double ymin = -fGeometryService->DetHalfHeight();
      double ymax = fGeometryService->DetHalfHeight();
      double zmin = 0.;
      double zmax = fGeometryService->DetLength();
      double rmin[3] = {xmin, ymin, zmin};
      double rmax[3] = {xmax, ymax, zmax};
      truthAlg.DrawCube(c1, rmin, rmax, 2);
    }

    // Draw the true particles
    TPolyLine3D *trajectories[100];
    TPolyLine3D *crttrack[100];
    int ncrtTracks = 0;
    size_t lim = particle.size();
    for(size_t i = 0; i < lim; i++){
      int id = particle[i].TrackId();
      bool plot = true;
      if(ind >= 0){ 
        plot = false;
        if(id == ind) plot = true;
      }
      if(truth && plot){
        int nTraj = particle[i].NumberTrajectoryPoints();
        trajectories[i] = new TPolyLine3D(nTraj);
        int ipt = 0;
        for(int j = 0; j < nTraj; j++){
          if(abs(particle[i].Vx(j))<500 && particle[i].Vy(j)<900 && particle[i].Vy(j)>-450 && particle[i].Vz(j)<700 && particle[i].Vz(j)>-400){
            trajectories[i]->SetPoint(ipt, particle[i].Vx(j), particle[i].Vy(j), particle[i].Vz(j));
            ipt++;
          }
        }
        trajectories[i]->SetLineColor(851+i);
        trajectories[i]->SetLineWidth(2);
        trajectories[i]->Draw();
      }
    }
    std::map<int, RecoTruth> tMatch = truthMatch;
    if(truthMatch.find(ind)!=truthMatch.end()){
      tMatch.clear();
      tMatch[ind] = truthMatch[ind];
    }
    for(auto &irt : tMatch){
      RecoTruth rt = irt.second;
      if(strips){
        // Plot the hit strips
        for(size_t j = 0; j < rt.stripHits.size(); j++){
          // Calculate the limits
          std::vector<double> limits = ChannelToLimits(rt.stripHits[j]);
          // Plot a rectangle
          double rmin[3] = {limits[0], limits[2], limits[4]};
          double rmax[3] = {limits[1], limits[3], limits[5]};
          truthAlg.DrawCube(c1, rmin, rmax, 2);
        }
      }
      if(hits){
        // Plot the hits
        for(size_t j = 0; j < rt.crtHits.size(); j++){
          CRTHit ht = rt.crtHits[j];
          // Get the limits
          double rmin[3] = {ht.xpos-ht.xerr,ht.ypos-ht.yerr,ht.zpos-ht.zerr};
          double rmax[3] = {ht.xpos+ht.xerr,ht.ypos+ht.yerr,ht.zpos+ht.zerr};
          truthAlg.DrawCube(c1, rmin, rmax, 3);
        }
      }
      if(tracks){
        // Plot the tracks
        for(size_t j = 0; j < rt.crtTracks.size(); j++){
          // Get the start and end points
          CRTTrack tr = rt.crtTracks[j];
          crttrack[ncrtTracks] = new TPolyLine3D(2);
          crttrack[ncrtTracks]->SetPoint(0, tr.xstart, tr.ystart, tr.zstart);
          crttrack[ncrtTracks]->SetPoint(1, tr.xend, tr.yend, tr.zend);
          // Draw a line between them
          crttrack[ncrtTracks]->SetLineColor(2);
          crttrack[ncrtTracks]->SetLineWidth(2);
          if (true/*tr.complete && CrossesTPC(tr)*/){
            crttrack[ncrtTracks]->Draw();
          }
          ncrtTracks++;
        }
      }
    }

    c1->SaveAs("crtTagger.root");
  }

  std::pair<std::string,unsigned> CRTRecoAna::ChannelToTagger(uint32_t channel){
    int strip = (channel >> 1) & 15;
    int module = (channel >> 5);
    std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
    TVector3 center = fAuxDetGeoCore->AuxDetChannelToPosition(2*strip, name);
    const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);

    std::set<std::string> volNames = {stripGeo.TotalVolume()->GetName()};
    std::vector<std::vector<TGeoNode const*> > paths = fGeometryService->FindAllVolumePaths(volNames);
    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
      path += paths.at(0).at(inode)->GetName();
      if (inode < paths.at(0).size() - 1) {
        path += "/";
      }
    }
    TGeoManager* manager = fGeometryService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeModule = manager->GetMother(2);
    TGeoNode* nodeTagger = manager->GetMother(3);
    // Module position in parent (tagger) frame
    double origin[3] = {0, 0, 0};
    double modulePosMother[3];
    nodeModule->LocalToMaster(origin, modulePosMother);
    unsigned planeID = (modulePosMother[2] > 0);
    std::string tagName = nodeTagger->GetName();
    std::pair<std::string, unsigned> output = std::make_pair(tagName, planeID);
    return output;
  }


  // WARNING: Relies on all modules in a tagger having the same dimensions
  bool CRTRecoAna::CheckModuleOverlap(uint32_t channel){
    bool hasOverlap = false;
    // Get the module ID
    int strip = (channel >> 1) & 15;
    int module = (channel >> 5);
    // Get the name of the module
    std::string name = fGeometryService->AuxDet(module).TotalVolume()->GetName();
    // Get the tagger TGeoNode
    const geo::AuxDetSensitiveGeo stripGeo = fAuxDetGeoCore->ChannelToAuxDetSensitive(name, 2*strip);
    std::set<std::string> volNames = {stripGeo.TotalVolume()->GetName()};
    std::vector<std::vector<TGeoNode const*> > paths = fGeometryService->FindAllVolumePaths(volNames);
    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
      path += paths.at(0).at(inode)->GetName();
      if (inode < paths.at(0).size() - 1) {
        path += "/";
      }
    }
    TGeoManager* manager = fGeometryService->ROOTGeoManager();
    manager->cd(path.c_str());
    TGeoNode* nodeModule = manager->GetMother(2);
    TGeoNode* nodeTagger = manager->GetMother(3);
    std::string modName = nodeModule->GetName();
    // Get the limits of the module in the tagger frame
    double height = fGeometryService->AuxDet(module).HalfHeight();
    double width = fGeometryService->AuxDet(module).HalfWidth1();
    double length = fGeometryService->AuxDet(module).Length()/2.;
    double pos1[3] = {width, height, length};
    double tagp1[3];
    nodeModule->LocalToMaster(pos1, tagp1);
    double pos2[3] = {-width, -height, -length};
    double tagp2[3];
    nodeModule->LocalToMaster(pos2, tagp2);
    std::vector<double> limits = {std::min(tagp1[0],tagp2[0]),
                                  std::max(tagp1[0],tagp2[0]),
                                  std::min(tagp1[1],tagp2[1]),
                                  std::max(tagp1[1],tagp2[1]),
                                  std::min(tagp1[2],tagp2[2]),
                                  std::max(tagp1[2],tagp2[2])};

    double origin[3] = {0, 0, 0};
    double modulePosMother[3];
    nodeModule->LocalToMaster(origin, modulePosMother);
    unsigned planeID = (modulePosMother[2] > 0);

    // Get the number of daughters from the tagger
    int nDaughters = nodeTagger->GetNdaughters();
    // Loop over the daughters
    for(int mod_i = 0; mod_i < nDaughters; mod_i++){
      // Check the name not the same as the current module
      TGeoNode* nodeDaughter = nodeTagger->GetDaughter(mod_i);
      std::string d_name = nodeDaughter->GetName();
      // Remove last two characters from name to match the AuxDet name
      if(d_name == modName) continue;
      // Get the limits of the module in the tagger frame
      double d_tagp1[3];
      nodeDaughter->LocalToMaster(pos1, d_tagp1);
      double d_tagp2[3];
      nodeDaughter->LocalToMaster(pos2, d_tagp2);
      std::vector<double> d_limits = {std::min(d_tagp1[0],d_tagp2[0]),
                                      std::max(d_tagp1[0],d_tagp2[0]),
                                      std::min(d_tagp1[1],d_tagp2[1]),
                                      std::max(d_tagp1[1],d_tagp2[1]),
                                      std::min(d_tagp1[2],d_tagp2[2]),
                                      std::max(d_tagp1[2],d_tagp2[2])};
      double d_modulePosMother[3];
      nodeDaughter->LocalToMaster(origin, d_modulePosMother);
      unsigned d_planeID = (d_modulePosMother[2] > 0);

      // Check the overlap of the two modules
      std::vector<double> overlap = CrtOverlap(limits, d_limits);
      // If there is an overlap set to true
      if(overlap[0]!=-99999 && d_planeID!=planeID) hasOverlap = true;
    }
    return hasOverlap;
  }


  bool CRTRecoAna::CrossesTPC(CRTTrack track){
    // Check if particle enters the TPC
    bool enters = false;
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();
    // Get track info
    TVector3 start(track.xstart, track.ystart, track.zstart);
    TVector3 end(track.xend, track.yend, track.zend);
    TVector3 diff = end - start;
    // Loop over trajectory points
    int npts = 100;
    for (int traj_i = 0; traj_i < npts; traj_i++){
      TVector3 trajPoint = start + ((traj_i+1)/100.)*diff;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
        enters = true;
      }
    }
    return enters;
  }


  std::vector<CRTHit> CRTRecoAna::AverageHits(std::vector<CRTHit> hits, std::map<int, RecoTruth>& truthMatch){
    std::vector<CRTHit> returnHits;
    std::vector<CRTHit> aveHits;
    std::vector<CRTHit> spareHits;
    if (hits.size()>0){
      // loop over size of tx
      bool first = true;
      TVector3 middle(0., 0., 0.);
      for (size_t i = 0; i < hits.size(); i++){
        // Get the position of the hit
        TVector3 pos(hits[i].xpos, hits[i].ypos, hits[i].zpos);
        // If first then set average = hit pos
        if(first){
          middle = pos;
          first = false;
        }
        // If distance from average < 10 cm then add to average
        if((pos-middle).Mag() < fAverageHitDistance){
          aveHits.push_back(hits[i]);
        }
        // Else add to another vector
        else{
          spareHits.push_back(hits[i]);
        }
      }

      CRTHit aveHit = DoAverage(aveHits);
      returnHits.push_back(aveHit);

      for(int& ID : aveHit.ids){
        truthMatch[ID].crtAveHits.push_back(aveHit);
      }

      std::vector<CRTHit> moreHits = AverageHits(spareHits, truthMatch);
      returnHits.insert(returnHits.end(), moreHits.begin(), moreHits.end());
      return returnHits;

    }
    else {
      return returnHits;
    }
  }


  CRTHit CRTRecoAna::DoAverage(std::vector<CRTHit> hits){
    std::string tagger = hits[0].tagger;
    double xpos = 0.; 
    double ypos = 0.;
    double zpos = 0.;
    double xmax = -99999; double xmin = 99999;
    double ymax = -99999; double ymin = 99999;
    double zmax = -99999; double zmin = 99999;
    double ts0_s = 0.;
    double ts0_ns = 0.;
    double ts1_ns = 0.;
    std::vector<int> ids;
    int nhits = 0;
    // Loop over hits
    for( auto& hit : hits ){
      // Get the mean x,y,z and times
      xpos += hit.xpos;
      ypos += hit.ypos;
      zpos += hit.zpos;
      ts0_s += hit.ts0_s;
      ts0_ns += hit.ts0_ns;
      ts1_ns += hit.ts1_ns;
      // For the errors get the maximum limits
      if(hit.xpos + hit.xerr > xmax) xmax = hit.xpos + hit.xerr;
      if(hit.xpos - hit.xerr < xmin) xmin = hit.xpos - hit.xerr;
      if(hit.ypos + hit.yerr > ymax) ymax = hit.ypos + hit.yerr;
      if(hit.ypos - hit.yerr < ymin) ymin = hit.ypos - hit.yerr;
      if(hit.zpos + hit.zerr > zmax) zmax = hit.zpos + hit.zerr;
      if(hit.zpos - hit.zerr < zmin) zmin = hit.zpos - hit.zerr;
      // Add all the unique IDs in the vector
      ids.insert(ids.end(), hit.ids.begin(), hit.ids.end());
      nhits++;
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    CRTHit crtHit = {tagger, xpos/nhits, ypos/nhits, zpos/nhits, 
                     (xmax-xmin)/2., (ymax-ymin)/2., (zmax-zmin)/2., 
                     ts0_s/nhits, ts0_ns/nhits, ts1_ns/nhits, ids};
    return crtHit;
  }
 

  //Function to create tracks from tzero hit collections
  std::vector<CRTTrack> CRTRecoAna::CreateTracks(std::vector<CRTHit> hits, std::map<int, RecoTruth>& truthMatch){
    std::vector<CRTTrack> returnTracks;

    //Store list of hit pairs with distance between them
    std::vector<std::pair<std::pair<size_t, size_t>, double>> hitPairDist;
    std::vector<std::pair<size_t, size_t>> usedPairs;
    //Calculate the distance between all hits on different planes
    for(size_t i = 0; i < hits.size(); i++){
      CRTHit hit1 = hits[i];
      for(size_t j = 0; j < hits.size(); j++){
        CRTHit hit2 = hits[j];
        std::pair<size_t, size_t> hitPair = std::make_pair(i, j);
        std::pair<size_t, size_t> rhitPair = std::make_pair(j, i);
        //Only compare hits on different taggers and don't reuse hits
        if(hit1.tagger!=hit2.tagger && std::find(usedPairs.begin(), usedPairs.end(), rhitPair)==usedPairs.end()){
          //Calculate the distance between hits and store
          TVector3 pos1(hit1.xpos, hit1.ypos, hit1.zpos);
          TVector3 pos2(hit2.xpos, hit2.ypos, hit2.zpos);
          double dist = (pos1 - pos2).Mag();
          usedPairs.push_back(hitPair);
          hitPairDist.push_back(std::make_pair(hitPair, dist));
        }
      }
    }

    //Sort map by distance
    std::sort(hitPairDist.begin(), hitPairDist.end(), [](auto& left, auto& right){
                return left.second > right.second;});

    //Store potential hit collections + distance along 1D hit
    std::vector<std::pair<std::vector<size_t>, double>> tracks;
    for(size_t i = 0; i < hitPairDist.size(); i++){
      size_t hit_i = hitPairDist[i].first.first;
      size_t hit_j = hitPairDist[i].first.second;
      //Make sure bottom plane hit is always hit_i
      if(hits[hit_j].tagger=="volTaggerBot_0") std::swap(hit_i, hit_j);
      CRTHit ihit = hits[hit_i];
      CRTHit jhit = hits[hit_j];

      //If the bottom plane hit is a 1D hit
      if(ihit.xerr>100. || ihit.zerr>100.){
        double facMax = 1;
        std::vector<size_t> nhitsMax;
        double minDist = 99999;

        //Loop over the length of the 1D hit
        for(int i = 0; i<21; i++){
          double fac = (i)/10.;
          std::vector<size_t> nhits;
          double totalDist = 0.;
           TVector3 start(ihit.xpos-(1.-fac)*ihit.xerr, ihit.ypos, ihit.zpos-(1.-fac)*ihit.zerr);
          TVector3 end(jhit.xpos, jhit.ypos, jhit.zpos);
          TVector3 diff = start - end;

          //Loop over the rest of the hits
          for(size_t k = 0; k < hits.size(); k++){
            if(k == hit_i || k == hit_j || hits[k].tagger == ihit.tagger || hits[k].tagger == jhit.tagger) continue;
            //Calculate the distance between the track crossing point and the true hit
            TVector3 mid(hits[k].xpos, hits[k].ypos, hits[k].zpos);
            TVector3 cross = CrossPoint(hits[k], start, diff);
            double dist = (cross-mid).Mag();
            //If the distance is less than some limit add the hit to the track and record the distance
            if(dist < fDistanceLimit){ 
              nhits.push_back(k);
              totalDist += dist;
            }
          }
          //If the distance down the 1D hit means more hits ore included and they are closer to the track record it
          if(nhits.size()>=nhitsMax.size() && totalDist/nhits.size() < minDist){
            nhitsMax = nhits;
            facMax = fac;
            minDist = totalDist/nhits.size();
          }
          nhits.clear();
        }
        //Record the track candidate
        std::vector<size_t> trackCand;
        trackCand.push_back(hit_i);
        trackCand.push_back(hit_j);
        trackCand.insert(trackCand.end(), nhitsMax.begin(), nhitsMax.end());
        tracks.push_back(std::make_pair(trackCand, facMax));
      }

      //If there is no 1D hit
      else{
        TVector3 start(ihit.xpos, ihit.ypos, ihit.zpos);
        TVector3 end(jhit.xpos, jhit.ypos, jhit.zpos);
        TVector3 diff = start - end;
        std::vector<size_t> trackCand;
        trackCand.push_back(hit_i);
        trackCand.push_back(hit_j);
        //Loop over all the other hits
        for(size_t k = 0; k < hits.size(); k++){
          if(k == hit_i || k == hit_j || hits[k].tagger == ihit.tagger || hits[k].tagger == jhit.tagger) continue;
          //Calculate distance to other hits not on the planes of the track hits
          TVector3 mid(hits[k].xpos, hits[k].ypos, hits[k].zpos);
          TVector3 cross = CrossPoint(hits[k], start, diff);
          double dist = (cross-mid).Mag();
          //Record any within a certain distance
          if(dist < fDistanceLimit){ 
            trackCand.push_back(k);
          }
        }
        tracks.push_back(std::make_pair(trackCand, 1));
      }
    }

    //Sort track candidates by number of hits
    std::sort(tracks.begin(), tracks.end(), [](auto& left, auto& right){
              return left.first.size() > right.first.size();});

    //Record used hits
    std::vector<size_t> usedHits;
    //Loop over candidates
    for(auto& track : tracks){
      size_t hit_i = track.first[0];
      size_t hit_j = track.first[1];

      // Make sure the first hit is the top high tagger if there are only two hits
      if(hits[hit_j].tagger=="volTaggerTopHigh_0") std::swap(hit_i, hit_j);
      CRTHit ihit = hits[track.first[0]];
      CRTHit jhit = hits[track.first[1]];

      //Calculate the time and record the ids
      double time = (ihit.ts1_ns + jhit.ts1_ns)/2.;
      std::vector<int> ids;
      //Check no hits in track have been used
      bool used = false;

      //Loop over hits in track candidate
      for(size_t i = 0; i < track.first.size(); i++){
        //Check if any of the hits have been used
        if(std::find(usedHits.begin(), usedHits.end(), track.first[i]) != usedHits.end()) used=true;
        CRTHit hit = hits[track.first[i]];
        ids.insert(ids.end(), hit.ids.begin(), hit.ids.end());
      }
      //If any of the hits have already been used skip this track
      if(used) continue;
      //Remove id duplicates
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      //Create track
      CRTTrack crtTrack = {ihit.xpos-(1.-track.second)*ihit.xerr, ihit.ypos, ihit.zpos-(1.-track.second)*ihit.zerr, jhit.xpos, jhit.ypos, jhit.zpos, time, ids, true};

      //If only the top two planes are hit create an incomplete/stopping track
      if(track.first.size()==2 && ihit.tagger == "volTaggerTopHigh_0" && jhit.tagger == "volTaggerTopLow_0"){ 
        TVector3 start(ihit.xpos, ihit.ypos, ihit.zpos);
        TVector3 end(jhit.xpos, jhit.ypos, jhit.zpos);
        TVector3 diff = (end - start).Unit();
        //Project the end of the track into the tpc
        TVector3 projEnd(start.X()+1000.*diff.X(), start.Y()+1000.*diff.Y(), start.Z()+1000.*diff.Z()); 
        CRTTrack t_crtTrack = {start.X(), start.Y(), start.Z(), projEnd.X(), projEnd.Y(), projEnd.Z(), time, ids, false};
        crtTrack = t_crtTrack;
      }
      returnTracks.push_back(crtTrack);
      for(auto& ID : ids){
        truthMatch[ID].crtTracks.push_back(crtTrack);
      }
      //Record which hits were used only if the track has more than two hits
      //If there are multiple 2 hit tracks there is no way to distinguish between them
      //TODO: Add charge matching to distinguish
      for(size_t i = 0; i < track.first.size(); i++){
        if(track.first.size()>2) usedHits.push_back(track.first[i]);
      }
    }
    return returnTracks;
  }


  //Function to calculate the crossing point of a track and tagger
  TVector3 CRTRecoAna::CrossPoint(CRTHit hit, TVector3 start, TVector3 diff){
    TVector3 cross;
    if(hit.xerr > 0.39 && hit.xerr < 0.41){
      double xc = hit.xpos;
      TVector3 crossp(xc, (xc-start.X())/(diff.X())*diff.Y()+start.Y(), (xc-start.X())/(diff.X())*diff.Z()+start.Z());
      cross = crossp;
    }
    else if(hit.yerr > 0.39 && hit.yerr < 0.41){
      double yc = hit.ypos;
      TVector3 crossp((yc-start.Y())/(diff.Y())*diff.X()+start.X(), yc, (yc-start.Y())/(diff.Y())*diff.Z()+start.Z());
      cross = crossp;
    }
    else if(hit.zerr > 0.39 && hit.zerr < 0.41){
      double zc = hit.zpos;
      TVector3 crossp((zc-start.Z())/(diff.Z())*diff.X()+start.X(), (zc-start.Z())/(diff.Z())*diff.Y()+start.Y(), zc);
      cross = crossp;
    }
    return cross;
  }

  DEFINE_ART_MODULE(CRTRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


