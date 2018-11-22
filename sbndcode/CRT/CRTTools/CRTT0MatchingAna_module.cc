////////////////////////////////////////////////////////////////////////
// Class:       CRTT0MatchingAnaAna
// Module Type: analyzer
// File:        CRTT0MatchingAnaAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

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
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace {
  // Local namespace for local functions
  // Declare here, define later
  // Print a TVector3
  void PrintVect(TVector3 vector, std::string name){
    std::cout<<name<<": ("<<vector.X()<<", "<<vector.Y()<<", "<<vector.Z()<<")\n";
  }

}

namespace sbnd {

  class CRTT0MatchingAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CrtHitModuleLabel {
        Name("CrtHitModuleLabel"),
        Comment("tag of CRT hit producer data product")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of tpc track producer data product")
      };

      fhicl::Atom<double> DistanceLimit {
        Name("DistanceLimit"),
        Comment("Maximum distance between projected crossing point and CRT hit for T0 to be considered")
      };

      fhicl::Atom<double> MinTrackLength {
        Name("MinTrackLength"),
        Comment("Minimum track length to perform matching on")
      };

      fhicl::Atom<double> TrackDirectionFrac {
        Name("TrackDirectionFrac"),
        Comment("Fraction of track to average direction over")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };
      
    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTT0MatchingAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Calculate the position that an MC particle crosses a crt plane
    TVector3 TaggerCrossPoint(const simb::MCParticle& particle, int tag_i);

    // Utility function that determines the possible x range of a track
    std::pair<double, double> TrackT0Range(double startX, double endX, int tpc);

    // Convert an xyz position to a CRT channel ID
    uint32_t PositionToChannel(TVector3 position);

    // Convert start time to CRT crossing point
    TVector3 T0ToXYZPosition(TVector3 position, TVector3 direction, std::string tagger, int tpc, double t0);

    // Calculate the distance from the track crossing point to CRT overlap coordinates
    double DistToCrtHit(TVector3 trackPos, crt::CRTHit crtHit);

    // Calculate the distance of closest approach between the end of a track and a crt hit
    double DistOfClosestApproach(TVector3 trackPos, TVector3 trackDir, crt::CRTHit crtHit, int tpc, double t0);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtHitModuleLabel;   ///< name of CRT producer
    art::InputTag fTpcTrackModuleLabel; ///< name of CRT producer
    double        fDistanceLimit;       ///< Maximum distance between projected crossing point and CRT hit
    double        fMinTrackLength;      ///< Minimum track length to perform T0 matching on
    double        fTrackDirectionFrac;  ///< Minimum track length to perform T0 matching on
    bool          fVerbose;             ///< print information about what's going on

    // n-tuples
    TH1D* fDistance;               ///< Distance between projected (at true T0) and true crossing points
    TH1D* fTrueTime;               ///< True time of track
    TH1D* fTrueRecoDist;           ///< Minimum distance between reco start and end points and matched true track trajectory
    TH1D* fMinDistCross;           ///< Minimum distance between projected crossing point and CRT hit for crossing track
    TH1D* fMinDistCont;            ///< Minimum distance between projected crossing point and CRT hit for a contained track
    TH1D* fTracksPerEvent;         ///< Number of tracks per event
    TH1D* fMatchesPerEvent;        ///< Number of T0 matched per event
    TH1D* fCrossLength;            ///< Length of tracks which cross CRTs
    TGraphAsymmErrors* fPurity;    ///< Correct matches / all matches
    TGraphAsymmErrors* fEffCross;  ///< Correct matches / all crossers
    TGraphAsymmErrors* fEffCont;   ///< Contained with no match / all contained
    TGraphAsymmErrors* fEffTotal;  ///< (Correct crossers + correct contained) / all tracks
    TGraphAsymmErrors* fEffLength; ///< (Correct crossers + correct contained) / all tracks

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider

    // Positions of the CRT planes
    std::vector<double> crtPlanes = {-359.1, -357.3, 357.3, 359.1, -358.9, -357.1, 661.52, 663.32, 865.52, 867.32, -240.65, -238.85, 655.35, 657.15};
    std::vector<int> fixCoord   = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2}; // Fixed coordinate for each plane
    std::vector<int> widthCoord = {2, 1, 2, 1, 0, 2, 2, 0, 2, 0, 1, 0, 1, 0}; // Width direction for each plane
    std::vector<int> lenCoord   = {1, 2, 1, 2, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1}; // Length direction for each plane

    // Performance Counters
    int nCorrectExit = 0; // True track crosses CRT and algorithm assigns the right T0
    int nIncorrectExit = 0; // True track crosses CRT and algorithm assigns wrong/no T0
    int nCorrectCont = 0; // True track contained and algorithm assigns no T0
    int nIncorrectCont = 0; // True track contained and algorithm assigns a T0

    // Counters for different limits
    int nLims = 10;
    double mindx = 6.;
    TH1D *hCorrectExit = new TH1D("hCorrectExit", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hTotalExit = new TH1D("hTotalExit", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hCorrectCont = new TH1D("hCorrectCont", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hTotalCont = new TH1D("hTotalCont", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hTotalCorrect = new TH1D("hTotalCorrect", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hTotalMatch = new TH1D("hTotalMatch", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hTotal = new TH1D("hTotal", "", nLims, mindx/2., nLims*mindx + mindx/2.);
    TH1D *hCorrectLen = new TH1D("hCorrectLen", "", 40, 0, 200);
    TH1D *hTotalLen = new TH1D("hTotalLen", "", 40, 0., 200);

    std::map<std::string, int> nameToInd;

  }; // class CRTT0MatchingAna


  // Constructor
  CRTT0MatchingAna::CRTT0MatchingAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel      (config().SimModuleLabel())
    , fCrtHitModuleLabel   (config().CrtHitModuleLabel())
    , fTpcTrackModuleLabel (config().TpcTrackModuleLabel())
    , fDistanceLimit       (config().DistanceLimit())
    , fMinTrackLength      (config().MinTrackLength())
    , fTrackDirectionFrac  (config().TrackDirectionFrac())
    , fVerbose             (config().Verbose())
  {

    // Get a pointer to the geometry service provider
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 

  } //CRTT0MatchingAna()


  void CRTT0MatchingAna::beginJob()
  {

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    fDistance        = tfs->make<TH1D>("distance",        ";Distance (cm)",                  100, 0,     500 );
    fTrueTime        = tfs->make<TH1D>("truetime",        ";True time (ticks)",              300, -3000, 3500);
    fTrueRecoDist    = tfs->make<TH1D>("truerecodist",    ";Distance (cm)",                  100, 0,     10  );
    fMinDistCross    = tfs->make<TH1D>("mindistcross",    ";Distance (cm)",                  100, 0,     1000);
    fMinDistCont     = tfs->make<TH1D>("mindistcont",     ";Distance (cm)",                  100, 0,     1000);
    fTracksPerEvent  = tfs->make<TH1D>("tracksperevent",  ";Reconstructed tracks per event", 40,  0,     40  );
    fMatchesPerEvent = tfs->make<TH1D>("matchesperevent", ";CRT matches per event",          40,  0,     40  );
    fCrossLength     = tfs->make<TH1D>("crosslength",     ";Reco length (cm)",               40,  0,     200 );
    fPurity          = tfs->makeAndRegister<TGraphAsymmErrors>("purity",   ";Distance Limit (cm);Purity"    );
    fEffCross        = tfs->makeAndRegister<TGraphAsymmErrors>("effcross", ";Distance Limit (cm);Efficiency");
    fEffCont         = tfs->makeAndRegister<TGraphAsymmErrors>("effcont",  ";Distance Limit (cm);Efficiency");
    fEffTotal        = tfs->makeAndRegister<TGraphAsymmErrors>("efftotal", ";Distance Limit (cm);Efficiency");
    fEffLength       = tfs->makeAndRegister<TGraphAsymmErrors>("efflength",";Reco length (cm);Efficiency");

    nameToInd["volTaggerSideRight_0"] = 0;
    nameToInd["volTaggerSideLeft_0"] = 1;
    nameToInd["volTaggerBot_0"] = 2;
    nameToInd["volTaggerTopLow_0"] = 3;
    nameToInd["volTaggerTopHigh_0"] = 4;
    nameToInd["volTaggerFaceFront_0"] = 5;
    nameToInd["volTaggerFaceBack_0"] = 6;

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT T0 Matching Ana Module -------------------"<<std::endl;

    // Take a position that is known to be the center of a strip

  } // CRTT0MatchingAna::beginJob()


  void CRTT0MatchingAna::analyze(const art::Event& event)
  {

    int nTracks = 0;
    int nMatches = 0;

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Detector properties
    double readoutWindow  = (double)fDetectorProperties->ReadOutWindowSize();
    double driftTimeTicks = fDetectorClocks->Time2Tick((2.*fGeometryService->DetHalfWidth()+3.)/fDetectorProperties->DriftVelocity());

    if(fVerbose) std::cout<<"Readout window = "<<readoutWindow<<" ticks, Drift time = "<<driftTimeTicks
                          <<" ticks, Drift velocity = "<<fDetectorProperties->DriftVelocity()<<" cm/us \n";

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // Fill a map of true particles
    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      particles[partId] = particle;
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    // Retrieve list of CRT hits
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
    std::map<std::string, std::vector<crt::CRTHit>> crtHitMap;
    for (auto const& crtHit : (*crtHitHandle)){
      crtHitMap[crtHit.tagger].push_back(crtHit);
    }

    // Retrieve the tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    if(fVerbose) std::cout<<"Number of CRT hits = "<<crtHitHandle->size()<<std::endl
                          <<"Number of TPC tracks = "<<tpcTrackHandle->size()<<std::endl;

    // Loop over all the tracks in the event
    for (auto const& track : (*tpcTrackHandle)){

      if (fVerbose) std::cout<<"\n---->Track "<<track.ID()<<":"<<std::endl;
      if (track.Length() < fMinTrackLength){ if(fVerbose) std::cout<<"Track too short!\n"; continue; }

      // Match the track to a true particle
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());
      int tpc = hits[0]->WireID().TPC;
      if (tpc != (int)hits[hits.size()-1]->WireID().TPC){ if(fVerbose) std::cout<<"Track has been stitched!\n"; continue; }

      // ===================== Get Truth Information ======================== //
      if (fVerbose) std::cout<<"Number of hits = "<<hits.size()<<", TPC = "<<tpc<<std::endl;
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits);
      if (particles.find(trueId) == particles.end()){ if (fVerbose) std::cout<<"No valid true track!\n"; continue; }

      // Get the true T0
      double trueTime = fDetectorClocks->TPCG4Time2Tick(particles[trueId].T());
      fTrueTime->Fill(trueTime);

      if(fVerbose) std::cout<<"True particle information:\n"<<"PDG = "<<particles[trueId].PdgCode()<<", length = "
                            <<particles[trueId].Trajectory().TotalLength()<<" cm, time = "<<trueTime<<" ticks\n\n";

      nTracks++;

      // Calculate direction as an average over directions
      size_t nTrackPoints = track.NumberTrajectoryPoints();
      int endPoint = (int)floor(nTrackPoints*fTrackDirectionFrac);
      double xTotStart = 0; double yTotStart = 0; double zTotStart = 0;
      double xTotEnd = 0; double yTotEnd = 0; double zTotEnd = 0;
      for(int i = 0; i < endPoint; i++){
        xTotStart -= track.DirectionAtPoint(i)[0];
        yTotStart -= track.DirectionAtPoint(i)[1];
        zTotStart -= track.DirectionAtPoint(i)[2];
        xTotEnd += track.DirectionAtPoint(nTrackPoints - (i+1))[0];
        yTotEnd += track.DirectionAtPoint(nTrackPoints - (i+1))[1];
        zTotEnd += track.DirectionAtPoint(nTrackPoints - (i+1))[2];
      } 
      TVector3 startDir = {xTotStart/endPoint, yTotStart/endPoint, zTotStart/endPoint};
      TVector3 endDir = {xTotEnd/endPoint, yTotEnd/endPoint, zTotEnd/endPoint};

      // Get the start and end points
      TVector3 start = track.Vertex();
      TVector3 end = track.End();

      if (fVerbose){
        PrintVect(start, "Track start");
        PrintVect(end, "Track end");
        PrintVect(startDir, "Start direction");
        PrintVect(endDir, "End direction");
      }

      //Shift start and end by true time
      TVector3 startShift = start;
      TVector3 endShift = end;
      double timeToShift = fDetectorClocks->TPCTick2Time(trueTime) * fDetectorProperties->DriftVelocity();
      if(tpc == 0){
        timeToShift = -timeToShift;
      }
      startShift[0] += timeToShift;
      endShift[0] += timeToShift;
      //Loop over the true tajectory and find minimum distance from the shifted reco track
      size_t npts = particles[trueId].NumberTrajectoryPoints();
      double minStartDist = 99999;
      double minEndDist = 99999;
      for (size_t i = 0; i < npts; i++){
        TVector3 trajpoint = {particles[trueId].Vx(i), particles[trueId].Vy(i), particles[trueId].Vz(i)};
        //Calculate the distance between the trajectory point and the start and end
        double startDist = (startShift - trajpoint).Mag();
        double endDist = (endShift - trajpoint).Mag();
         //Record the minimum distance for each of the start and end
        if(startDist<minStartDist) minStartDist = startDist;
        if(endDist<minEndDist) minEndDist = endDist;
      }
      fTrueRecoDist->Fill(minStartDist);
      fTrueRecoDist->Fill(minEndDist);

      // Determine whether the true particle crosses the CRT planes and if there's a corresponding CRT hit
      bool trueCrtCross = false;
      if(fVerbose) std::cout<<"\nTruth(ish) matched CRT hits:\n";
      std::map<std::string, TVector3> trueXYZ;
      std::map<std::string, bool> trueMatch;

      // Loop over the taggers
      for(auto &taggerHits : crtHitMap){
        std::string tagger = taggerHits.first;
        // Calculate the crossing point of the true particle
        int tag_i = nameToInd[tagger];
        TVector3 trueCross = TaggerCrossPoint(particles[trueId], tag_i);
        trueXYZ[tagger] = trueCross;
        bool isMatch = false;

        // Loop over the hits on the tagger
        for(auto &crtHit : taggerHits.second){
          double trueDist = DistToCrtHit(trueCross, crtHit);
          double crtTimeTicks = fDetectorClocks->TPCG4Time2Tick((double)(int)crtHit.ts1_ns);
          if(trueDist<20. && std::abs(crtTimeTicks-trueTime)<2.) {
            if(fVerbose) std::cout<<tagger<<": CRT pos = ("<<crtHit.x_pos<<", "<<crtHit.y_pos<<", "<<crtHit.z_pos
                                  <<"), time = "<<crtTimeTicks<<"\nTrue pos = ("<<trueCross.X()<<", "
                                  <<trueCross.Y()<<", "<<trueCross.Z()<<"), dist = "<<trueDist<<std::endl;
            trueCrtCross = true;
            isMatch = true;
          }
        }
        trueMatch[tagger] = isMatch;
      }
      if (trueCrtCross && fVerbose) std::cout<<"True particle crosses CRTs!\n";
      if (trueCrtCross) fCrossLength->Fill(track.Length());
      if (!trueCrtCross && fVerbose) std::cout<<"True particle doesn't cross CRTs!\n";

      // ====================== Matching Algorithm ========================== //
      // Get the allowed t0 range
      std::pair<double, double> t0MinMax = TrackT0Range(start.X(), end.X(), tpc);
      if(fVerbose) std::cout<<"T0 range: min = "<<t0MinMax.first<<", max = "<<t0MinMax.second<<", truth = "<<trueTime<<std::endl;

      std::vector<std::pair<double, double>> t0Candidates;
      // Fill vector for loop over distance limits
      std::vector<std::vector<std::pair<double,double>>> t0candLim;
      for(int i = 0; i < nLims; i++){
        std::vector<std::pair<double, double>> temp;
        t0candLim.push_back(temp);
      }
      double minDistance = 99999;

      // Loop over the taggers
      for(auto &taggerHits : crtHitMap){
        std::string tagger = taggerHits.first;
        if (fVerbose) std::cout<<"\nTagger "<<tagger<<"\n";

        // Calculate the start and end crossing points at the true T0
        TVector3 trueStartPos = T0ToXYZPosition(start, startDir, tagger, tpc, trueTime);
        TVector3 trueEndPos = T0ToXYZPosition(end, endDir, tagger, tpc, trueTime);
        if(fVerbose && trueStartPos[0] != -99999) PrintVect(trueStartPos, "Cross point at true time");
        if(fVerbose && trueEndPos[0] != -99999) PrintVect(trueEndPos, "Cross point at true time");

        // Calculate the distance between the true and reco crossing point
        double trueStartDist = (trueStartPos-trueXYZ[tagger]).Mag();
        double trueEndDist = (trueEndPos-trueXYZ[tagger]).Mag();
        if(trueStartPos[0] != -99999 && trueMatch[tagger]){ 
          fDistance->Fill(trueStartDist);
          if(fVerbose) std::cout<<"Distance to true cross point = "<<trueStartDist<<"\n";
        }
        if(trueEndPos[0] != -99999 && trueMatch[tagger]){ 
          fDistance->Fill(trueEndDist);
          if(fVerbose) std::cout<<"Distance to true cross point = "<<trueEndDist<<"\n";
        }

        // Loop over all the CRT hits
        for(auto &crtHit : taggerHits.second){
          // Check if hit is within the allowed t0 range
          double crtTimeTicks = fDetectorClocks->TPCG4Time2Tick((double)(int)crtHit.ts1_ns);
          if (!(crtTimeTicks >= t0MinMax.first-20. && crtTimeTicks <= t0MinMax.second+20.)) continue;
          TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);

          // Calculate the distance between the crossing point and the CRT hit
          TVector3 startTrackPos = T0ToXYZPosition(start, startDir, tagger, tpc, crtTimeTicks);
          //double startDist = DistToCrtHit(startTrackPos, crtHit);
          double startDist = DistOfClosestApproach(start, startDir, crtHit, tpc, crtTimeTicks);
          if (startDist < minDistance) minDistance = startDist;

          // If the distance is less than some limit record the time
          if (startDist < fDistanceLimit){ 
            if(fVerbose) std::cout<<"Match! Time = "<<crtTimeTicks<<", Position = ("<<startTrackPos.X()<<", "
                                  <<startTrackPos.Y()<<", "<<startTrackPos.Z()<<") dist = "<<startDist<<std::endl;
            t0Candidates.push_back(std::make_pair(startDist, crtTimeTicks));
          }

          // Calculate the distance between the crossing point and the CRT hit
          TVector3 endTrackPos = T0ToXYZPosition(end, endDir, tagger, tpc, crtTimeTicks);
          //double endDist = DistToCrtHit(endTrackPos, crtHit);
          double endDist = DistOfClosestApproach(end, endDir, crtHit, tpc, crtTimeTicks);
          if (endDist < minDistance) minDistance = endDist;

          // If the distance is less than some limit record the time
          if (endDist < fDistanceLimit){ 
            if(fVerbose) std::cout<<"Match! Time = "<<crtTimeTicks<<", Position = ("<<endTrackPos.X()<<", "
                                  <<endTrackPos.Y()<<", "<<endTrackPos.Z()<<") dist = "<<endDist<<std::endl;
            t0Candidates.push_back(std::make_pair(endDist, crtTimeTicks));
          }

          // Loop over different distance limits
          for(int i = 0; i < nLims; i++){
            if(startDist < (i+1.)*mindx) t0candLim[i].push_back(std::make_pair(startDist, crtTimeTicks));
            if(endDist < (i+1.)*mindx) t0candLim[i].push_back(std::make_pair(startDist, crtTimeTicks));
          }
        }
      }

      if(trueCrtCross) fMinDistCross->Fill(minDistance);
      else fMinDistCont->Fill(minDistance);

      // Use the t0 that appears the most often
      // Sort the candidates by distance
      std::sort(t0Candidates.begin(), t0Candidates.end(), [](auto& left, auto& right){
                return left.first < right.first;});
      if(fVerbose) std::cout<<"Number of t0 candidates = "<<t0Candidates.size()<<std::endl;
      // Loop over the candidates
      //double bestTime = -99999;
      //int bestOccur = 0;
      /*size_t t0_i = 0;
      bool containsTrue = false;
      while (t0_i < t0Candidates.size()){
        // Record the current time
        double currentTime = t0Candidates[t0_i].second;
        //int nOccur = 1;
        t0_i++;
        if(std::abs(currentTime-trueTime)<2.) containsTrue = true;
        // Loop over all times within 1 tick of this time
        while (std::abs(currentTime-t0Candidates[t0_i])<1.){
          nOccur++;
          t0_i++;
        }
        if (fVerbose) std::cout<<"t0 candidate = "<<currentTime<<": Number of occurences = "<<nOccur<<std::endl;
        if (nOccur>bestOccur){
          bestTime = currentTime;
          bestOccur = nOccur;
        }
      }*/
      double bestTime = -99999;
      if(t0Candidates.size()>0) bestTime = t0Candidates[0].second;

      if(bestTime != -99999) nMatches++;

      // Calculate performance
      if(trueCrtCross && std::abs(bestTime-trueTime)<2.){ if(fVerbose) std::cout<<"CORRECT CROSSER! <==\n"; nCorrectExit++; hCorrectLen->Fill(track.Length());}
      if(trueCrtCross && std::abs(bestTime-trueTime)>=2.){ if(fVerbose) std::cout<<"INCORRECT CROSSER! <==\n"; nIncorrectExit++; }
      if(!trueCrtCross && bestTime == -99999){ if(fVerbose) std::cout<<"CORRECT CONTAINED! <==\n"; nCorrectCont++; hCorrectLen->Fill(track.Length()); }
      if(!trueCrtCross && bestTime != -99999){ if(fVerbose) std::cout<<"INCORRECT CONTAINED! <==\n"; nIncorrectCont++; }

      hTotalLen->Fill(track.Length());

      for(int i = 0; i < nLims; i++){
        std::sort(t0candLim[i].begin(), t0candLim[i].end(), [](auto& left, auto& right){
                  return left.first < right.first;});
        double bestT = -99999;
        if(t0candLim[i].size()>0) bestT = t0candLim[i][0].second;
        if(trueCrtCross && std::abs(bestT-trueTime)<2.) {hCorrectExit->Fill((i+1.)*mindx); hTotalCorrect->Fill((i+1.)*mindx);}
        if(trueCrtCross) {hTotalExit->Fill((i+1.)*mindx); hTotal->Fill((i+1.)*mindx);}
        if(!trueCrtCross && bestT == -99999) {hCorrectCont->Fill((i+1.)*mindx); hTotalCorrect->Fill((i+1.)*mindx);}
        if(!trueCrtCross) {hTotalCont->Fill((i+1.)*mindx); hTotal->Fill((i+1.)*mindx);}
        if(bestT != -99999) hTotalMatch->Fill((i+1.)*mindx);
      }

    }

    fTracksPerEvent->Fill(nTracks);
    fMatchesPerEvent->Fill(nMatches);

  } // CRTT0MatchingAna::analyze()


  void CRTT0MatchingAna::endJob(){
  
    std::cout<<"--------------- Truth Matching Info: ----------------"<<std::endl
             <<"---------------- Performance Info: ------------------"<<std::endl
             <<"Number of correct exiters     = "<<nCorrectExit<<std::endl
             <<"Number of incorrect exiters   = "<<nIncorrectExit<<std::endl
             <<"Number of correct contained   = "<<nCorrectCont<<std::endl
             <<"Number of incorrect contained = "<<nIncorrectCont<<std::endl
             <<"Total correct = "<<(double)(nCorrectExit+nCorrectCont)/(nCorrectExit+nIncorrectExit+nCorrectCont+nIncorrectCont)<<std::endl;

    TCanvas *c1 = new TCanvas();
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);

    //Efficiency for crossing tracks = number of tracks with good match / total tracks with CRT hit
    fEffCross->SetMarkerStyle(8);
    fEffCross->SetMarkerColor(4);
    fEffCross->SetLineColor(4);
    fEffCross->SetLineWidth(3);
    fEffCross->BayesDivide(hCorrectExit, hTotalExit);
    fEffCross->Draw("ap");
    leg->AddEntry("fEffCross", "CRT crossers", "p");

    //Efficiency for contained tracks = number of contained tracks with no match / total tracks with no CRT hits
    fEffCont->SetMarkerStyle(8);
    fEffCont->SetMarkerColor(4);
    fEffCont->SetLineColor(4);
    fEffCont->SetLineWidth(3);
    // Divide the histograms to make a TGraphAsymmErrors
    fEffCont->BayesDivide(hCorrectCont, hTotalCont);
    fEffCont->Draw("p same");
    leg->AddEntry("fEffCont", "Contained", "p");    

    // Save to file
    leg->Draw("same");
    c1->SaveAs("effic.root");

    //Total efficiency = (tracks with good match + contained tracks with no match) / total tracks
    fEffTotal->SetTitle(";Distance Limit (cm);Correct/Total");
    fEffTotal->SetMarkerStyle(8);
    fEffTotal->SetMarkerColor(1);
    fEffTotal->SetLineColor(1);
    fEffTotal->SetLineWidth(3);
    fEffTotal->BayesDivide(hTotalCorrect, hTotal);
    fEffTotal->Draw("ap");

    //Purity = tracks with good match / all matched tracks
    fPurity->SetMarkerStyle(8);
    fPurity->SetMarkerColor(1);
    fPurity->SetLineColor(1);
    fPurity->SetLineWidth(3);
    fPurity->BayesDivide(hCorrectExit, hTotalMatch);
    fPurity->Draw("ap");

    fEffLength->SetMarkerStyle(8);
    fEffLength->SetMarkerColor(1);
    fEffLength->SetLineColor(1);
    fEffLength->SetLineWidth(3);
    fEffLength->BayesDivide(hCorrectLen, hTotalLen);
    fEffLength->Draw("ap");

  } // CRTT0MatchingAna::endJob()


  // Function to calculate the point where a true particle crosses a CRT tagger
  TVector3 CRTT0MatchingAna::TaggerCrossPoint(simb::MCParticle const& particle, int tag_i){

    double tagCenter[3] = {0, 0, 208.25};
    tagCenter[fixCoord[tag_i*2]] = (crtPlanes[tag_i*2]+crtPlanes[tag_i*2+1])/2;
    double tagDim[3] = {0, 0, 0};
    if(tag_i==0 || tag_i==1){ tagDim[0] = 1.8; tagDim[1] = 360; tagDim[2] = 450; }
    if(tag_i==2){ tagDim[0] = 399.5; tagDim[1] = 1.8; tagDim[2] = 478; }
    if(tag_i==3 || tag_i==4){ tagDim[0] = 450; tagDim[1] = 1.8; tagDim[2] = 450; }
    if(tag_i==5 || tag_i==6){ tagDim[0] = 360; tagDim[1] = 360; tagDim[2] = 1.8; }
    TVector3 start, end;
    bool first = true;
    // Get the trajectory of the true particle
    size_t npts = particle.NumberTrajectoryPoints();
    // Loop over particle trajectory
    for (size_t i = 0; i < npts; i++){
      TVector3 trajPoint(particle.Vx(i), particle.Vy(i), particle.Vz(i));
      // If the particle is inside the tagger volume then set to true.
      if(trajPoint[0]>tagCenter[0]-tagDim[0] && trajPoint[0]<tagCenter[0]+tagDim[0] &&
         trajPoint[1]>tagCenter[1]-tagDim[1] && trajPoint[1]<tagCenter[1]+tagDim[1] &&
         trajPoint[2]>tagCenter[2]-tagDim[2] && trajPoint[2]<tagCenter[2]+tagDim[2]){
        if(first) start = trajPoint;
        first = false;
        end = trajPoint;
      }
    }
    TVector3 crossPoint((start.X()+end.X())/2,(start.Y()+end.Y())/2,(start.Z()+end.Z())/2);

    return crossPoint;

  } // CRTT0MatchingAna::TaggerCrossPoint()


  // Utility function that determines the possible x range of a track
  std::pair<double, double> CRTT0MatchingAna::TrackT0Range(double startX, double endX, int tpc){

    std::pair<double, double> result;
    double Vd = fDetectorProperties->DriftVelocity();
    // Whole track must be within tpc
    // Find which TPC the track hits are detected in
    if(tpc == 0){
      // Lowest |X| is furthest from APA
      double lowX = std::max(startX, endX);
      // xmin is shift from furthest to 0 (the CPA)
      double xmax = 0 - lowX;
      // Highest |X| is closest to APA
      double highX = std::min(startX, endX);
      // xmax is shift from closest to APA
      double xmin = -(2.0*fGeometryService->DetHalfWidth()+3.) - highX;
      double t0max = -2.*xmin/Vd;
      double t0min = -2.*xmax/Vd;
      result = std::make_pair(t0min, t0max);
    }
    else{
      // Lowest |X| is furthest from APA
      double lowX = std::min(startX, endX);
      // xmin is shift from furthest to 0 (the CPA)
      double xmin = 0 - lowX;
      // Highest |X| is closest to APA
      double highX = std::max(startX, endX);
      // xmax is shift from closest to APA
      double xmax = (2.0*fGeometryService->DetHalfWidth()+3.) - highX;
      double t0min = 2.*xmin/Vd;
      double t0max = 2.*xmax/Vd;
      result = std::make_pair(t0min, t0max);
    }

    return result;

  } // CRTT0MatchingAna::TrackT0Range()


  // Function to project a track position on to a tagger
  TVector3 CRTT0MatchingAna::T0ToXYZPosition(TVector3 position, TVector3 direction, std::string tagger, int tpc, double t0){

    //Here crt_i is index of tagger, so runs from 0 to 6
    TVector3 returnVal(-99999, -99999, -99999);
    int crt_i = nameToInd[tagger];

    // Convert the t0 into an x shift
    double shift = fDetectorClocks->TPCTick2Time(t0) * fDetectorProperties->DriftVelocity();
    // Apply the shift depending on which TPC the track is in
    if (tpc == 1) position[0] += shift;
    if (tpc == 0) position[0] -= shift;

    // Calculate the step to the CRT plane
    double step = (crtPlanes[crt_i*2]- position[fixCoord[crt_i*2]])/direction[fixCoord[crt_i*2]];

    // If the step is < 0 return a null position
    if (step < 0) return returnVal;

    // Calculate the CRT crossing point of the output coordinate
    returnVal[lenCoord[crt_i*2]] = position[lenCoord[crt_i*2]] + step*direction[lenCoord[crt_i*2]];
    returnVal[widthCoord[crt_i*2]] = position[widthCoord[crt_i*2]] + step*direction[widthCoord[crt_i*2]];
    returnVal[fixCoord[crt_i*2]] = crtPlanes[crt_i*2];

    return returnVal;

  } // CRTT0MatchingAna::T0ToXYZPosition()


  // Function to calculate the distance between a projected cross point and a CRT hit
  double CRTT0MatchingAna::DistToCrtHit(TVector3 trackPos, crt::CRTHit crtHit){

    double minDistX = 99999;
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double xpos = crtHit.x_pos + ((i+1.)/10. - 1.)*crtHit.x_err;
      double distX = std::abs(trackPos.X() - xpos);
      if(distX < minDistX) minDistX = distX;
    }

    double minDistY = 99999;
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double ypos = crtHit.y_pos + ((i+1.)/10. - 1.)*crtHit.y_err;
      double distY = std::abs(trackPos.Y() - ypos);
      if(distY < minDistY) minDistY = distY;
    }

    double minDistZ = 99999;
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double zpos = crtHit.z_pos + ((i+1.)/10. - 1.)*crtHit.z_err;
      double distZ = std::abs(trackPos.Z() - zpos);
      if(distZ < minDistZ) minDistZ = distZ;
    }

    double dist = std::sqrt(std::pow(minDistX, 2) + std::pow(minDistY, 2) + std::pow(minDistZ, 2));
    return dist;

  } // CRTT0MatchingAna::DistToCrtHit()


  // Function to calculate the distance of closest approach to a CRT hit
  double CRTT0MatchingAna::DistOfClosestApproach(TVector3 trackPos, TVector3 trackDir, crt::CRTHit crtHit, int tpc, double t0){

    double minDist = 99999;

    // Convert the t0 into an x shift
    double shift = fDetectorClocks->TPCTick2Time(t0) * fDetectorProperties->DriftVelocity();
    // Apply the shift depending on which TPC the track is in
    if (tpc == 1) trackPos[0] += shift;
    if (tpc == 0) trackPos[0] -= shift;

    TVector3 endPos = trackPos + trackDir;
    double denominator = (endPos - trackPos).Mag();
    // 1D hits should only have a lot of variance in one direction
    if(crtHit.x_err > 50.){
      // Loop over size of hit to find the min dist
      for(int i = 0; i < 20.; i++){
        double xpos = crtHit.x_pos + ((i+1.)/10. - 1.)*crtHit.x_err;
        TVector3 crtPoint(xpos, crtHit.y_pos, crtHit.z_pos);
        double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
        double dca = numerator/denominator;
        if(dca < minDist) minDist = dca;
      }
    }
    else if(crtHit.y_err > 50.){
      // Loop over size of hit to find the min dist
      for(int i = 0; i < 20.; i++){
        double ypos = crtHit.y_pos + ((i+1.)/10. - 1.)*crtHit.y_err;
        TVector3 crtPoint(crtHit.x_pos, ypos, crtHit.z_pos);
        double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
        double dca = numerator/denominator;
        if(dca < minDist) minDist = dca;
      }
    }
    else if(crtHit.y_err > 50.){
      // Loop over size of hit to find the min dist
      for(int i = 0; i < 20.; i++){
        double zpos = crtHit.z_pos + ((i+1.)/10. - 1.)*crtHit.z_err;
        TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, zpos);
        double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
        double dca = numerator/denominator;
        if(dca < minDist) minDist = dca;
      }
    }
    else{
      TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }

    return minDist;

  } // CRTT0MatchingAna::DistToCrtHit()


  DEFINE_ART_MODULE(CRTT0MatchingAna)
} // namespace sbnd

// Back to local namespace.
namespace {


} // local namespace


