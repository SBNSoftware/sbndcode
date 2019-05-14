////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackMatchingAna
// Module Type: analyzer
// File:        CRTTrackMatchingAnaAna_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcorealg/CoreUtils/NumericUtils.h" // util::absDiff()

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/MCCheater/BackTrackerService.h"


// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

namespace sbnd {

  struct TrackMatch{
    int tpc;
    int trueID;
    std::vector<int> crtIDs;
    double trueTime;
  };

  struct RecoTruth{
    std::vector<simb::MCParticle> particles;
    std::vector<recob::Track> tpcTracks;
    std::vector<crt::CRTTrack> crtTracks;
    std::map<int, TrackMatch> matchingMap;
    std::vector<RecoCRTTrack> recoCrtTracks;
  };


  class CRTTrackMatchingAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CrtTrackModuleLabel {
        Name("CrtTrackModuleLabel"),
        Comment("tag of CRT track producer data product")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("tag of TPC track producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

      fhicl::Atom<double> MaxAngleDiff {
        Name("MaxAngleDiff"),
        Comment("Maximum difference between CRT and TPC track angles")
      };

      fhicl::Atom<double> MaxDistance {
        Name("MaxDistance"),
        Comment("Maximum distance between CRT and TPC track start/end positions")
      };

      fhicl::Atom<bool> Plot {
        Name("Plot"),
        Comment("Plot what's going on")
      };

      fhicl::Atom<int> PlotID {
        Name("PlotID"),
        Comment("TPC track ID to plot")
      };
      
    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTTrackMatchingAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Function to draw true and reco tracks
    void DrawTrueTracks(RecoTruth truthMatch, bool truth, bool tpctracks, bool crttracks, bool crtreco, int id);

    // Function to calculate the distance to a CRT track start/end point
    double DistToCrtTrack(TVector3 trackPos, TVector3 crtCross, TVector3 crtErr);

    // Function the calculate the average squared distance between a CRT and TPC track
    double TrackAveDCA(recob::Track track, TVector3 start, TVector3 end);

    // Function the calculate the average squared distance between a CRT and TPC track
    double TrueAveDCA(simb::MCParticle particle, TVector3 start, TVector3 end);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCrtTrackModuleLabel; ///< name of CRT track producer
    art::InputTag fTpcTrackModuleLabel; ///< name of TPC track producer
    bool          fVerbose;             ///< print information about what's going on
    double        fMaxAngleDiff;        ///< max difference between CRT and TPC angles
    double        fMaxDistance;         ///< max distance between CRT and TPC start/end positions
    bool          fPlot;                ///< plot what's going on
    int           fPlotID;              ///< tpc track ID to plot

    // histograms
    TH1D* fDeltaLengthMatch;
    TH1D* fDeltaThetaMatch;
    TH1D* fDeltaPhiMatch;
    TH1D* fDeltaDistMatch;
    TH1D* fChiSqMatch;
    TH2D* fDeltaAnglesMatch;
    TH1D* fDeltaLengthNoMatch;
    TH1D* fDeltaThetaNoMatch;
    TH1D* fDeltaPhiNoMatch;
    TH1D* fDeltaDistNoMatch;
    TH1D* fChiSqNoMatch;
    TH2D* fDeltaAnglesNoMatch;
    TH1D* fNCrtCandidates;
    TH1D* fNTpcCandidates;
    TH1D* fDCACrtTrue;
    TH1D* fDCACrtReco;
    TH1D* fDCARecoTrue;

    // Other variables shared between different methods.
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider
    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;

    CRTTrackMatchAlg trackAlg;
    CRTTruthRecoAlg truthAlg;

    // Performance Counters
    int nCorrectMatch = 0;
    int nIncorrectMatch = 0;
    int nCorrectNoMatch = 0;
    int nIncorrectNoMatch = 0;
    int nCorrectTime = 0;
    int nIncorrectTime = 0;

  }; // class CRTTrackMatchingAna


  // Constructor
  CRTTrackMatchingAna::CRTTrackMatchingAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCrtTrackModuleLabel  (config().CrtTrackModuleLabel())
    , fTpcTrackModuleLabel  (config().TpcTrackModuleLabel())
    , fVerbose              (config().Verbose())
    , fMaxAngleDiff         (config().MaxAngleDiff())
    , fMaxDistance          (config().MaxDistance())
    , fPlot                 (config().Plot())
    , fPlotID               (config().PlotID())
    , trackAlg()
    , truthAlg()
  {

    // Get a pointer to the geometry service provider
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks     = lar::providerFrom<detinfo::DetectorClocksService>(); 

  } // CRTTrackMatchingAna()


  void CRTTrackMatchingAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;
    fDeltaLengthMatch   = tfs->make<TH1D>("deltalengthmatch", ";#Delta length (cm)", 50, -700, 700);
    fDeltaThetaMatch    = tfs->make<TH1D>("deltathetamatch", ";#Delta #theta (deg)", 200, -180, 180);
    fDeltaPhiMatch      = tfs->make<TH1D>("deltaphimatch", ";#Delta #phi (deg)", 200, -180, 180);
    fDeltaDistMatch     = tfs->make<TH1D>("deltadistmatch", ";#Delta distance (cm)", 50, 0, 700);
    fChiSqMatch         = tfs->make<TH1D>("chisqmatch", ";#chi^2", 50, 0, 5000);
    fDeltaAnglesMatch   = tfs->make<TH2D>("deltaanglesmatch", "#Delta #theta (rad);#Delta #phi (rad)", 50, -4, 4, 50, -4, 4);
    fDeltaLengthNoMatch = tfs->make<TH1D>("deltalengthnomatch", ";#Delta length (cm)", 50, -700, 700);
    fDeltaThetaNoMatch  = tfs->make<TH1D>("deltathetanomatch", ";#Delta #theta (deg)", 200, -180, 180);
    fDeltaPhiNoMatch    = tfs->make<TH1D>("deltaphinomatch", ";#Delta #phi (deg)", 200, -180, 180);
    fDeltaDistNoMatch   = tfs->make<TH1D>("deltadistnomatch", ";#Delta distance (cm)", 50, 0, 700);
    fChiSqNoMatch       = tfs->make<TH1D>("chisqnomatch", ";#chi^2", 50, 0, 5000);
    fDeltaAnglesNoMatch = tfs->make<TH2D>("deltaanglesnomatch", "#Delta #theta (rad);#Delta #phi (rad)", 50, -4, 4, 50, -4, 4);
    fNCrtCandidates     = tfs->make<TH1D>("ncrtcandidates", ";N CRT candidates", 5, 0, 5);
    fNTpcCandidates     = tfs->make<TH1D>("ntpccandidates", ";N CRT candidates", 5, 0, 5);
    fDCACrtTrue         = tfs->make<TH1D>("distsqcrttrue", ";d^2", 50, 0, 100);
    fDCACrtReco         = tfs->make<TH1D>("distsqcrtreco", ";d^2", 50, 0, 100);
    fDCARecoTrue        = tfs->make<TH1D>("distsqrecotrue", ";d^2", 50, 0, 50);

    // Initial output
    if(fVerbose) std::cout<<"----------------- CRT Track Matching Ana Module -------------------"<<std::endl;

  } // CRTTrackMatchingAna::beginJob()


  void CRTTrackMatchingAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // Fill a map of true particles
    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    for (auto const& particle: (*particleHandle)){
      int partId = particle.TrackId();
      particles[partId] = particle;
      parts.push_back(particle);
    }

    if(fVerbose) std::cout<<"Number of true particles = "<<particles.size()<<std::endl;

    // Retrieve list of CRT tracks
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCrtTrackModuleLabel);

    // Retrieve the TPC tracks
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);

    // Get track to hit associations
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);

    if(fVerbose) std::cout<<"Number of CRT tracks = "<<crtTrackHandle->size()<<std::endl
                          <<"Number of TPC tracks = "<<tpcTrackHandle->size()<<std::endl;


    //----------------------------------   Truth Matching   -----------------------------------
    if(fVerbose) std::cout<<"\n----------------------------------   Truth Matching   -----------------------------------\n";

    std::map<int, TrackMatch> matchingMap;
    std::map<int, recob::Track> recoTracks;
    std::vector<recob::Track> tpcTracks;
    // Loop over TPC tracks
    for (auto const& tpcTrack : (*tpcTrackHandle)){
      TrackMatch trackMatch;
      trackMatch.trueTime = -99999;
      int trackID = tpcTrack.ID();
      if(fVerbose) std::cout<<"\n-->TPC Track:"<<trackID<<"\n";

      recoTracks[trackID] = tpcTrack;

      tpcTracks.push_back(tpcTrack);
      TVector3 tpcStart = tpcTrack.Vertex<TVector3>();
      TVector3 tpcEnd = tpcTrack.End<TVector3>();

      // Get the associated hits
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
      int tpc = fTpcGeo.DetectedInTPC(hits);
      trackMatch.tpc = tpc;

      // Find the true particle ID and check it exists
      if (fVerbose) std::cout<<"Number of hits = "<<hits.size()<<", TPC = "<<trackMatch.tpc<<std::endl;
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(hits);
      trackMatch.trueID = trueId;

      if (particles.find(trueId) == particles.end()){ 
        if (fVerbose) std::cout<<"No valid true tpcTrack!\n"; 
        matchingMap[trackID] = trackMatch; 
        continue; 
      }

      simb::MCParticle particle = particles[trueId];

      // Match angles to get rid of secondary particles
      std::pair<TVector3, TVector3> trueStartEnd = truthAlg.TpcCrossPoints(particle);
      TVector3 trueStart = trueStartEnd.first;
      TVector3 trueEnd = trueStartEnd.second;
      double tpcTheta = (tpcStart - tpcEnd).Theta();
      double tpcPhi = (tpcStart - tpcEnd).Phi();
      double trueTheta = (trueStart - trueEnd).Theta();
      double truePhi = (trueStart - trueEnd).Phi();
      if(std::max((trueStart-tpcStart).Mag(), (trueEnd-tpcEnd).Mag()) > std::max((trueStart-tpcEnd).Mag(), (trueEnd-tpcStart).Mag())){
        trueTheta = (trueEnd - trueStart).Theta();
        truePhi = (trueEnd - trueStart).Phi();
      }
      double dTheta = atan2(sin(tpcTheta - trueTheta), cos(tpcTheta - trueTheta));
      double dPhi = atan2(sin(tpcPhi - truePhi), cos(tpcPhi - truePhi));

      // Get the true t0
      double trueTime = particle.T() * 1e-3; // [us]
      trackMatch.trueTime = trueTime; 
      size_t nTrajPoints = particle.NumberTrajectoryPoints();

      if(fVerbose) std::cout<<"True particle information:\n"<<"PDG = "<<particle.PdgCode()<<", length = "
                            <<particle.Trajectory().TotalLength()<<" cm, time = "<<trueTime<<" us\n" 
                            <<"Reco start = ("<<tpcStart.X()<<", "<<tpcStart.Y()<<", "<<tpcStart.Z()
                            <<") reco end = ("<<tpcEnd.X()<<", "<<tpcEnd.Y()<<", "<<tpcEnd.Z()<<") length = "<<tpcTrack.Length()<<"\n";

      if(std::abs(dTheta) > 0.1 || std::abs(dPhi) > 0.1){ 
        if (fVerbose) std::cout<<"True-reco track angles don't match! dTheta = "<<dTheta<<", dPhi = "<<dPhi<<"\n";
        matchingMap[trackID] = trackMatch;
        continue;
      }

      // Loop over the CRT tracks
      int crtIndex = 0;
      for (auto const& crtTrack : (*crtTrackHandle)){
        // Only look at tracks which cross the TPC
        if(!trackAlg.CrossesTPC(crtTrack)){ crtIndex++; continue; }
        // Find the time difference with the true particle
        double crtTime = ((double)(int)crtTrack.ts1_ns) * 1e-3; // [us]
        double timeDiff = std::abs(crtTime - trueTime);

        // Find the distance between start and end points and the true particle
        TVector3 crtStart(crtTrack.x1_pos, crtTrack.y1_pos, crtTrack.z1_pos);
        TVector3 crtEnd(crtTrack.x2_pos, crtTrack.y2_pos, crtTrack.z2_pos);
        TVector3 crtStartErr(crtTrack.x1_err, crtTrack.y1_err, crtTrack.z1_err);
        TVector3 crtEndErr(crtTrack.x2_err, crtTrack.y2_err, crtTrack.z2_err);
        double minStartDist = 99999;
        double minEndDist = 99999;
        for (size_t i = 0; i < nTrajPoints; i++){
          TVector3 trajPoint(particle.Vx(i), particle.Vy(i), particle.Vz(i));
          double startDist = DistToCrtTrack(trajPoint, crtStart, crtStartErr);
          if(startDist < minStartDist) minStartDist = startDist;
          double endDist = DistToCrtTrack(trajPoint, crtEnd, crtEndErr);
          if(endDist < minEndDist) minEndDist = endDist;
        }
        // If there is a track at the same time and the start and end points are within a certain 
        //dist of the true trajectory then mark as "matchable"
        // Store in some kind of map with the track ID and crt track index
        if(timeDiff<10 && fVerbose){
          std::cout<<"CRT Track:"<<crtIndex<<"\n"
                   <<"Time = "<<crtTime<<" start = ("<<crtStart.X()<<", "<<crtStart.Y()<<", "<<crtStart.Z()
                   <<") end = ("<<crtEnd.X()<<", "<<crtEnd.Y()<<", "<<crtEnd.Z()<<")\n"
                   <<"Time diff = "<<timeDiff<<", min start dist = "<<minStartDist<<", min end dist = "<<minEndDist<<"\n";
        }

        if(timeDiff < 1. && minStartDist < 40. && minEndDist < 40.){ 
          trackMatch.crtIDs.push_back(crtIndex);
          if(fVerbose) std::cout<<"Matches CRT track, time diff = "<<timeDiff<<", min start dist = "<<minStartDist<<", min end dist = "<<minEndDist<<"\n";
        }

        crtIndex++;
      }
      matchingMap[trackID] = trackMatch;
    }

    //---------------------------------- Matching Algorithm -----------------------------------
    if(fVerbose) std::cout<<"\n---------------------------------- CRTTrack->RecoCRTTrack -----------------------------------\n";

    //TODO Account for crt track errors
    int crtIndex = 0;
    std::vector<crt::CRTTrack> crtTracks;
    std::vector<RecoCRTTrack> recoCrtTracks;
    // Loop over CRT tracks
    for (auto const& crtTrack : (*crtTrackHandle)){

      crtTracks.push_back(crtTrack);

      //Check that crt track crosses tpc volume, if not skip it
      if(!trackAlg.CrossesTPC(crtTrack)){ crtIndex++; continue; }

      std::vector<RecoCRTTrack> tempTracks = trackAlg.CrtToRecoTrack(crtTrack, crtIndex);
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());

      crtIndex++;
    }

    if(fVerbose) std::cout<<"\n---------------------------------- Matching Algorithm -----------------------------------\n";
    // Loop over the reco crt tracks
    std::map<int, int> nTpcCandidates;
    for (auto const& recoCrtTrack : recoCrtTracks){
      if(fVerbose) std::cout<<"\nRecoCRTTrack: "<<recoCrtTrack.crtID<<"\n";
      std::vector<std::pair<int, double>> crtTpcMatchCandidates;

      // Get the length, angle and start and end position of the TPC track
      TVector3 crtStart = recoCrtTrack.start;
      TVector3 crtEnd = recoCrtTrack.end;
      double crtLength = (crtStart - crtEnd).Mag();
      double crtTheta = (crtStart - crtEnd).Theta();
      double crtPhi = (crtStart - crtEnd).Phi();

      // Loop over the TPC tracks
      for (auto const& tpcTrack : (*tpcTrackHandle)){
        int trackID = tpcTrack.ID();
        //if (tpcTrack.Length() < 10.){ if(fVerbose) std::cout<<"Track too short!\n"; continue; }
        // If the tpcTrack has been stitched across the CPA it already has an associated t0
        std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(tpcTrack.ID());
        int tpc = hits[0]->WireID().TPC;
        //if (tpc != (int)hits[hits.size()-1]->WireID().TPC){ if(fVerbose) std::cout<<"Track has been stitched!\n"; continue;}

        // Get the length, angle and start and end position of the TPC track
        TVector3 tpcStart = tpcTrack.Vertex<TVector3>();
        TVector3 tpcEnd = tpcTrack.End<TVector3>();
        double tpcLength = (tpcStart - tpcEnd).Mag();
        double tpcTheta = (tpcStart - tpcEnd).Theta();
        double tpcPhi = (tpcStart - tpcEnd).Phi();

        // Find the difference with the CRT track
        double dDist1 = (crtStart-tpcStart).Mag();
        double dDist2 = (crtEnd-tpcEnd).Mag();
        if(std::max((crtStart-tpcStart).Mag(), (crtEnd-tpcEnd).Mag()) > std::max((crtStart-tpcEnd).Mag(), (crtEnd-tpcStart).Mag())){
          crtTheta = (crtEnd - crtStart).Theta();
          crtPhi = (crtEnd - crtStart).Phi();
          dDist1 = (crtEnd-tpcStart).Mag();
          dDist2 = (crtStart-tpcEnd).Mag();
        }
        double dLength = tpcLength - crtLength;
        double dTheta = atan2(sin(tpcTheta - crtTheta), cos(tpcTheta - crtTheta));
        double dPhi = atan2(sin(tpcPhi - crtPhi), cos(tpcPhi - crtPhi));
        double aveDCA = TrackAveDCA(tpcTrack, crtStart, crtEnd);

        // Plot difference for tracks with and without true match
        if(std::find(matchingMap[trackID].crtIDs.begin(), matchingMap[trackID].crtIDs.end(), recoCrtTrack.crtID) != matchingMap[trackID].crtIDs.end()
           && matchingMap[trackID].tpc == recoCrtTrack.tpc){

          if(fVerbose) std::cout<<"TPC track start = ("<<tpcStart.X()<<", "<<tpcStart.Y()<<", "<<tpcStart.Z()
                                <<"), end = ("<<tpcEnd.X()<<", "<<tpcEnd.Y()<<", "<<tpcEnd.Z()<<")\n"
                                <<"True match: ID = "<<trackID<<" tpc len = "<<tpcLength<<" crt len = "<<crtLength<<" diff = "<<dLength
                                <<" tpc theta = "<<tpcTheta<<" crt theta = "<<crtTheta<<" diff = "<<dTheta<<" tpc phi = "
                                <<tpcPhi<<" crt phi = "<<crtPhi<<" diff = "<<dPhi<<" dist1 = "<<dDist1<<" dist2 = "<<dDist2<<"\n";

          fDeltaLengthMatch->Fill(dLength);
          fDeltaThetaMatch->Fill(dTheta*180/(TMath::Pi()));
          fDeltaPhiMatch->Fill(dPhi*180/(TMath::Pi()));
          fDeltaDistMatch->Fill(std::min(dDist1,dDist2));
          fChiSqMatch->Fill(aveDCA);
          fDeltaAnglesMatch->Fill(dTheta, dPhi);
        }
        else {
          fDeltaLengthNoMatch->Fill(dLength);
          fDeltaThetaNoMatch->Fill(dTheta*180/(TMath::Pi()));
          fDeltaPhiNoMatch->Fill(dPhi*180/(TMath::Pi()));
          fDeltaDistNoMatch->Fill(std::min(dDist1,dDist2));
          fChiSqNoMatch->Fill(aveDCA);
          fDeltaAnglesNoMatch->Fill(dTheta, dPhi);
        }

        // Do the actual matching
        if((std::abs(dTheta) < fMaxAngleDiff && std::abs(dPhi) < fMaxAngleDiff && tpc == recoCrtTrack.tpc) && (dDist1<fMaxDistance||dDist2<fMaxDistance)){
          crtTpcMatchCandidates.push_back(std::make_pair(trackID, std::abs(dTheta)));
          nTpcCandidates[trackID]++;

          if(std::find(matchingMap[trackID].crtIDs.begin(), matchingMap[trackID].crtIDs.end(), recoCrtTrack.crtID) != matchingMap[trackID].crtIDs.end()){ 
            if(fVerbose) std::cout<<"Correct match! ID = "<<trackID<<" dTheta = "<<dTheta<<" dPhi = "<<dPhi<<" dist1 = "<<dDist1
                                  <<" dist2 = "<<dDist2<<" dLength = "<<dLength<<" aveDCA = "<<aveDCA<<"\n";
          }
          else{ 
            if(fVerbose) std::cout<<"Incorrect match! ID = "<<trackID<<" dTheta = "<<dTheta<<" dPhi = "<<dPhi<<" dist1 = "<<dDist1
                                  <<" dist2 = "<<dDist2<<" dLength = "<<dLength<<" aveDCA = "<<aveDCA<<"\n";
          }

        }
        else if(std::find(matchingMap[trackID].crtIDs.begin(), matchingMap[trackID].crtIDs.end(), recoCrtTrack.crtID) != matchingMap[trackID].crtIDs.end() && tpc == recoCrtTrack.tpc){
          if(fVerbose) std::cout<<"Incorrect non match! ID = "<<trackID<<" Theta = "<<dTheta<<" dPhi = "<<dPhi<<" dist1 = "<<dDist1
                                <<" dist2 = "<<dDist2<<" dLength = "<<dLength<<" aveDCA = "<<aveDCA<<"\n";
        } 
      }
      // Choose the track which matches the closest
      int matchedTrackID = -99999;
      if(crtTpcMatchCandidates.size() > 0){
        std::sort(crtTpcMatchCandidates.begin(), crtTpcMatchCandidates.end(), [](auto& left, auto& right){
                  return left.second < right.second;});
        matchedTrackID = crtTpcMatchCandidates[0].first;
      }
      fNCrtCandidates->Fill(crtTpcMatchCandidates.size());

      if(matchedTrackID != -99999){
        if(std::find(matchingMap[matchedTrackID].crtIDs.begin(), matchingMap[matchedTrackID].crtIDs.end(), recoCrtTrack.crtID) != matchingMap[matchedTrackID].crtIDs.end()){
          if(fVerbose) std::cout<<"Matched correctly: Tpc ID = "<<matchedTrackID<<" crt ID = "<<recoCrtTrack.crtID<<" true time = "<<matchingMap[matchedTrackID].trueTime
                                <<" crt time = "<<recoCrtTrack.trueTime<<" n cands = "<<crtTpcMatchCandidates.size()<<"\n";
          nCorrectMatch++;
        }
        else{ 
          if(fVerbose) std::cout<<"Matched incorrectly: Tpc ID = "<<matchedTrackID<<" crt ID = "<<recoCrtTrack.crtID<<" true time = "<<matchingMap[matchedTrackID].trueTime
                                <<" crt time = "<<recoCrtTrack.trueTime<<" n cands = "<<crtTpcMatchCandidates.size()<<"\n";
          nIncorrectMatch++;
        }
        if(std::abs(matchingMap[matchedTrackID].trueTime - recoCrtTrack.trueTime) < 5.) nCorrectTime++;
        else nIncorrectTime++;
        // Check the errors on the CRT track
        bool smallErrors = true;
        int ind = 0;
        TVector3 start, end;
        for (auto const& crtTrack : (*crtTrackHandle)){
          if(ind != recoCrtTrack.crtID){ind++; continue;}
          if(crtTrack.x1_err<10. && crtTrack.x2_err<10. && crtTrack.y1_err<10. && crtTrack.y2_err<10. && crtTrack.z1_err<10. && crtTrack.z2_err<10.){ 
            smallErrors = true;
            start[0] = crtTrack.x1_pos; start[1] = crtTrack.y1_pos; start[2] = crtTrack.z1_pos;
            end[0] = crtTrack.x2_pos; end[1] = crtTrack.y2_pos; end[2] = crtTrack.z2_pos;
          }
          ind++;
        }
        // If small then fill the dist^2 histograms
        if(smallErrors){
          double trackAveDca = TrackAveDCA(recoTracks[matchedTrackID], crtStart, crtEnd);
          double trueAveDca = TrueAveDCA(particles[matchingMap[matchedTrackID].trueID], start, end);
          fDCACrtReco->Fill(trackAveDca);
          fDCACrtTrue->Fill(trueAveDca);
          fDCARecoTrue->Fill(std::abs(trackAveDca-trueAveDca));
        }
      }
      else{
        bool matchable = false;
        for(auto const& matching : matchingMap){
          if(std::find(matching.second.crtIDs.begin(), matching.second.crtIDs.end(), recoCrtTrack.crtID) != matching.second.crtIDs.end() && matching.second.tpc == recoCrtTrack.tpc){ 
            if(fVerbose) std::cout<<"Should have matched: Tpc ID = "<<matching.first<<" crt ID = "<<recoCrtTrack.crtID
                                  <<" true time = "<<matching.second.trueTime<<" crt time = "<<recoCrtTrack.trueTime<<"\n";
            matchable = true;
          }
        }
        if(matchable){ 
          nIncorrectNoMatch++;
        }
        else{
          if(fVerbose) std::cout<<"Correctly didn't match anything: crt ID = "<<recoCrtTrack.crtID<<" crt time = "<<recoCrtTrack.trueTime<<"\n";
          nCorrectNoMatch++;
        }
      }
    }

    for(auto const& tpcCands : nTpcCandidates){
      fNTpcCandidates->Fill(tpcCands.second);
    }

    RecoTruth truthMatch;
    truthMatch.particles = parts;
    truthMatch.crtTracks = crtTracks;
    truthMatch.tpcTracks = tpcTracks;
    truthMatch.matchingMap = matchingMap;
    truthMatch.recoCrtTracks = recoCrtTracks;
    if(fPlot){
      DrawTrueTracks(truthMatch, true, true, false, false, fPlotID);
    }

    //--------------------------------- Performance Analysis ----------------------------------

  } // CRTTrackMatchingAna::analyze()


  void CRTTrackMatchingAna::endJob(){

    std::cout<<"----------------- RESULTS --------------------\n"
             <<"Number of correct matches       = "<<nCorrectMatch<<"\n"
             <<"Number of incorrect matches     = "<<nIncorrectMatch<<"\n"
             <<"Number of correct non matches   = "<<nCorrectNoMatch<<"\n"
             <<"Number of incorrect non matches = "<<nIncorrectNoMatch<<"\n"
             <<"Efficiency = "<<((double)(nCorrectMatch+nCorrectNoMatch))/((double)(nCorrectMatch+nCorrectNoMatch+nIncorrectMatch+nIncorrectNoMatch))<<"\n"
             <<"Purity = "<<((double)(nCorrectMatch))/((double)(nCorrectMatch+nIncorrectMatch))<<"\n"
             <<"Number of matches with correct time = "<<nCorrectTime<<", incorrect time = "<<nIncorrectTime<<" purity = "<<((double)nCorrectTime)/(nCorrectTime+nIncorrectTime)<<"\n";

    
    TCanvas *c2 = new TCanvas("c2","",700,700);
    fDeltaThetaNoMatch->Draw();
    fDeltaThetaMatch->Draw("SAME");
    c2->SaveAs("deltaTheta.root");
    TCanvas *c3 = new TCanvas("c3","",700,700);
    fDeltaPhiNoMatch->Draw();
    fDeltaPhiMatch->Draw("SAME");
    c3->SaveAs("deltaPhi.root");
    TCanvas *c4 = new TCanvas("c4","",700,700);
    fDeltaDistNoMatch->Draw();
    fDeltaDistMatch->Draw("SAME");
    c4->SaveAs("deltaDist.root");
  
  } // CRTTrackMatchingAna::endJob()

  
  // Function to draw true and reco tracks
  void CRTTrackMatchingAna::DrawTrueTracks(RecoTruth rt, bool truth, bool tpcTracks, bool crtTracks, bool crtreco, int id){

    // Create a canvas 
    TCanvas *c1 = new TCanvas("c1","",700,700);
    // Draw the tagger planes
    for(size_t tag_i = 0; tag_i < fCrtGeo.NumTaggers(); tag_i++){
      double rmin[3] = {fCrtGeo.GetTagger(tag_i).minX, fCrtGeo.GetTagger(tag_i).minY, fCrtGeo.GetTagger(tag_i).minZ};
      double rmax[3] = {fCrtGeo.GetTagger(tag_i).maxX, fCrtGeo.GetTagger(tag_i).maxY, fCrtGeo.GetTagger(tag_i).maxZ};
      truthAlg.DrawCube(c1, rmin, rmax, 1);
    }

    double xmin = fTpcGeo.MinX(); 
    double xmax = fTpcGeo.MaxX();
    double ymin = fTpcGeo.MinY();
    double ymax = fTpcGeo.MaxY();
    double zmin = fTpcGeo.MinZ();
    double zmax = fTpcGeo.MaxZ();
    double rmin[3] = {xmin, ymin, zmin};
    double rmax[3] = {0, ymax, zmax};
    truthAlg.DrawCube(c1, rmin, rmax, 1);
    double rmin1[3] = {0, ymin, zmin};
    double rmax1[3] = {xmax, ymax, zmax};
    truthAlg.DrawCube(c1, rmin1, rmax1, 1);

    // Draw the true particles
    TPolyLine3D *trajectories[100];
    TPolyLine3D *crttrack[100];
    TPolyLine3D *crtrecotrack[100];
    TPolyLine3D *tpctrack[100];
    int ncrtTracks = 0;
    int ncrt = 0;
    int nparts = 0;
    size_t lim = rt.particles.size();

    if(truth){
      for(size_t i = 0; i < lim; i++){
        int trueID = rt.particles[i].TrackId();
        bool plot = false;
        for(auto& matching : (rt.matchingMap)){
          if(matching.second.trueID == trueID) plot = true;
        }
        if(plot){
          int nTraj = rt.particles[i].NumberTrajectoryPoints();
          trajectories[nparts] = new TPolyLine3D(nTraj);
          int ipt = 0;
          for(int j = 0; j < nTraj; j++){
            double px = rt.particles[i].Vx(j);
            double py = rt.particles[i].Vy(j);
            double pz = rt.particles[i].Vz(j);
            if(abs(px) < 500 && py < 900 && py > -450 && pz < 700 && pz > -400){
              trajectories[nparts]->SetPoint(ipt, px, py, pz);
              ipt++;
            }
          }
          trajectories[nparts]->SetLineColor(4);
          trajectories[nparts]->SetLineWidth(2);
          if(id==-1||rt.matchingMap[id].trueID==trueID){ 
            trajectories[nparts]->Draw();
            nparts++;
          }
        }
      }
    }

    if(crtTracks){
      // Plot the tracks
      for(size_t i = 0; i < rt.crtTracks.size(); i++){
        // Get the start and end points
        crt::CRTTrack tr = rt.crtTracks[i];
        crttrack[ncrtTracks] = new TPolyLine3D(2);
        crttrack[ncrtTracks]->SetPoint(0, tr.x1_pos, tr.y1_pos, tr.z1_pos);
        crttrack[ncrtTracks]->SetPoint(1, tr.x2_pos, tr.y2_pos, tr.z2_pos);
        // Draw a line between them
        crttrack[ncrtTracks]->SetLineColor(2);
        crttrack[ncrtTracks]->SetLineWidth(2);
        if ((id==-1 && trackAlg.CrossesTPC(tr)) || std::find(rt.matchingMap[id].crtIDs.begin(), rt.matchingMap[id].crtIDs.end(), (int)i) != rt.matchingMap[id].crtIDs.end()){
          if(tr.complete){
            crttrack[ncrtTracks]->Draw();
            ncrtTracks++;
          }
          if(!tr.complete){
            TVector3 start(tr.x1_pos, tr.y1_pos, tr.z1_pos);
            TVector3 end(tr.x2_pos, tr.y2_pos, tr.z2_pos);
            if(start.Y() < end.Y()) std::swap(start, end);
            TVector3 diff = (end - start).Unit();
            TVector3 newEnd = start + 1000*diff;
            crttrack[ncrtTracks]->SetPoint(0, start.X(), start.Y(), start.Z());
            crttrack[ncrtTracks]->SetPoint(1, newEnd.X(), newEnd.Y(), newEnd.Z());
            crttrack[ncrtTracks]->Draw();
            ncrtTracks++;
          }
        }
      }
    }

    if(crtreco){
      // Plot the tracks
      for(size_t i = 0; i < rt.recoCrtTracks.size(); i++){
        // Get the start and end points
        TVector3 st = rt.recoCrtTracks[i].start;
        TVector3 ed = rt.recoCrtTracks[i].end;
        crtrecotrack[ncrt] = new TPolyLine3D(2);
        crtrecotrack[ncrt]->SetPoint(0, st.X(), st.Y(), st.Z());
        crtrecotrack[ncrt]->SetPoint(1, ed.X(), ed.Y(), ed.Z());
        // Draw a line between them
        crtrecotrack[ncrt]->SetLineColor(46);
        crtrecotrack[ncrt]->SetLineWidth(2);
        if (!(st.X()==0 && ed.X()==0) && (id==-1||std::find(rt.matchingMap[id].crtIDs.begin(), rt.matchingMap[id].crtIDs.end(), rt.recoCrtTracks[i].crtID) != rt.matchingMap[id].crtIDs.end())){
          crtrecotrack[ncrt]->Draw();
          ncrt++;
        }
      }
    }

    if(tpcTracks){
      // Plot the tracks
      for(size_t i = 0; i < rt.tpcTracks.size(); i++){
        // Get the start and end points
        recob::Track tr = rt.tpcTracks[i];
        size_t npts = tr.NumberTrajectoryPoints();
        tpctrack[i] = new TPolyLine3D(npts);
        for(size_t j = 0; j < npts; j++){
          auto& pos = tr.LocationAtPoint(j);
          tpctrack[i]->SetPoint(j, pos.X(), pos.Y(), pos.Z());
        }
        // Draw a line between them
        tpctrack[i]->SetLineColor(3);
        tpctrack[i]->SetLineWidth(2);
        //if(id == -1 || tr.ID() == id) tpctrack[i]->Draw();
        if(tr.ID() == 3 || tr.ID() == id) tpctrack[i]->Draw();
      }
    }

    c1->SaveAs("crtTagger.root");

  } // CRTTrackMatchingAna::DrawTrueTracks()

  // Function to calculate the distance to a CRT track start/end point
  double CRTTrackMatchingAna::DistToCrtTrack(TVector3 trackPos, TVector3 crtCross, TVector3 crtErr){

    double minDistX = 99999;
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double xpos = crtCross.X() + ((i+1.)/10. - 1.)*crtErr.X();
      double distX = std::abs(trackPos.X() - xpos);
      if(distX < minDistX) minDistX = distX;
    }

    double minDistY = 99999;
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double ypos = crtCross.Y() + ((i+1.)/10. - 1.)*crtErr.Y();
      double distY = std::abs(trackPos.Y() - ypos);
      if(distY < minDistY) minDistY = distY;
    }

    double minDistZ = 99999;
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double zpos = crtCross.Z() + ((i+1.)/10. - 1.)*crtErr.Z();
      double distZ = std::abs(trackPos.Z() - zpos);
      if(distZ < minDistZ) minDistZ = distZ;
    }

    double dist = std::sqrt(std::pow(minDistX, 2) + std::pow(minDistY, 2) + std::pow(minDistZ, 2));
    return dist;

  } // CRTTrackMatchingAna::DistToCrtTrack()


  // Function the calculate the average squared distance between a CRT and TPC track
  double CRTTrackMatchingAna::TrackAveDCA(recob::Track track, TVector3 start, TVector3 end){

    double denominator = (end - start).Mag();

    // Loop over track trajectory
    int nTraj = track.NumberTrajectoryPoints();
    int ipt = 0;
    double xres = 0;
    for(int j = 0; j < nTraj; j++){
      TVector3 xp = track.LocationAtPoint<TVector3>(j);
      double numerator = ((xp - start).Cross(xp - end)).Mag();
      if(xp[0]==-999) continue;
      double d2 = numerator/denominator;
      xres += d2;
      ipt++;
    }

    return xres/ipt;

  } // CRTTrackMatchingAna::TrackAveDCA()


  // Function the calculate the average squared distance between a CRT and TPC track
  double CRTTrackMatchingAna::TrueAveDCA(simb::MCParticle particle, TVector3 start, TVector3 end){

    double denominator = (end - start).Mag();

    double xmin = fTpcGeo.MinX(); 
    double xmax = fTpcGeo.MaxX();
    double ymin = fTpcGeo.MinY();
    double ymax = fTpcGeo.MaxY();
    double zmin = fTpcGeo.MinZ();
    double zmax = fTpcGeo.MaxZ();

    // Get the trajectory of the true particle
    size_t npts = particle.NumberTrajectoryPoints();
    int ipt = 0;
    double xres = 0;

    // Loop over particle trajectory
    for (size_t i = 0; i < npts; i++){
      TVector3 trajPoint(particle.Vx(i), particle.Vy(i), particle.Vz(i));
      // If the particle is inside the tagger volume then set to true.
      if(trajPoint[0]>xmin && trajPoint[0]<xmax &&
         trajPoint[1]>ymin && trajPoint[1]<ymax &&
         trajPoint[2]>zmin && trajPoint[2]<zmax){
        double numerator = ((trajPoint - start).Cross(trajPoint - end)).Mag();
        double d2 = numerator/denominator;
        xres += d2;
        ipt++;
      }
    }

    return xres/ipt;

  } // CRTTrackMatchingAna::TrueAveDCA()


  DEFINE_ART_MODULE(CRTTrackMatchingAna)
} // namespace sbnd

