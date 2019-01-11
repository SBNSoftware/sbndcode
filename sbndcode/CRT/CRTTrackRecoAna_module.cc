////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackRecoAna
// Module Type: analyzer
// File:        CRTTrackRecoAna_module.cc
//
// Analysis module for evaluating CRT track reconstruction on through going
// muons.
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTruthRecoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackRecoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTruthMatchUtils.h"

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

namespace sbnd {

  struct RecoTruth{
    std::vector<crt::CRTHit> crtHits;
    std::vector<crt::CRTHit> crtAveHits;
    std::vector<crt::CRTTrack> crtTracks;
  };

  class CRTTrackRecoAna : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> CRTHitModuleLabel {
        Name("CRTHitModuleLabel"),
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

      fhicl::Table<CRTTrackRecoAlg::Config> TrackAlg {
        Name("TrackAlg"),
        Comment("")
      };

      fhicl::Atom<bool> UseTopPlane {
        Name("UseTopPlane"),
        Comment("Use hits from the top plane (SBND specific)")
      };

      fhicl::Atom<bool> UseReadoutWindow {
        Name("UseReadoutWindow"),
        Comment("")
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTTrackRecoAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;
  
    void DrawTrueTracks(const std::vector<simb::MCParticle>& particle, std::map<int, RecoTruth> truthMatch, bool truth, bool hits, bool tracks, bool tpc, int ind);


  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTModuleLabel;      ///< name of CRT hit producer
    art::InputTag fCRTHitModuleLabel;      ///< name of CRT hit producer
    bool          fVerbose;             ///< print information about what's going on
    bool          fVeryVerbose;         ///< print more information about what's going on
    bool          fPlot;                ///< plot tracks
    int           fPlotTrackID;         ///< id of track to plot
    bool          fUseTopPlane;         // Use hits from the top plane (SBND specific)
    bool          fUseReadoutWindow;         // Use hits from the top plane (SBND specific)
    
    // n-tuples
    TH1D* fCrtHitsPerTrack;
    TH1D* fAveHitsPerTrack;
    TH1D* fTracksPerTrack;

    TH1D* fTrackDist;

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider

    CRTTruthRecoAlg truthAlg;
    CRTTrackRecoAlg trackAlg;

    // Positions of the CRT planes
    std::vector<double> crtPlanes = {-359.1, -357.3, 357.3, 359.1, -358.9, -357.1, 661.52, 663.32, 865.52, 867.32, -240.65, -238.85, 655.35, 657.15};
    std::vector<int> fixCoord   = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2}; // Fixed coordinate for each plane
    std::vector<int> widthCoord = {2, 1, 2, 1, 0, 2, 2, 0, 2, 0, 1, 0, 1, 0}; // Width direction for each plane
    std::vector<int> lenCoord   = {1, 2, 1, 2, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1}; // Length direction for each plane

    const size_t nTaggers = 7;

    // Performance Counters

    // Truth Counters
    int nParts = 0;
    int nSinglePlane = 0;
    int nMatchTrack = 0;
    int nTwoStripPlanes = 0;

  }; // class CRTTrackRecoAna

  // Constructor
  CRTTrackRecoAna::CRTTrackRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTModuleLabel       (config().CRTModuleLabel())
    , fCRTHitModuleLabel    (config().CRTHitModuleLabel())
    , fVerbose              (config().Verbose())
    , fVeryVerbose          (config().VeryVerbose())
    , fPlot                 (config().Plot())
    , fPlotTrackID          (config().PlotTrackID())
    , fUseTopPlane          (config().UseTopPlane())
    , fUseReadoutWindow     (config().UseReadoutWindow())
    , truthAlg()
    , trackAlg(config().TrackAlg())
  {
    // Get a pointer to the fGeometryServiceetry service provider
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 
  }

  void CRTTrackRecoAna::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fCrtHitsPerTrack   = tfs->make<TH1D>("crthitspertrack",    ";CRT Hits;N Tracks",  10, 0, 10);
    fAveHitsPerTrack   = tfs->make<TH1D>("avehitspertrack",    ";Average Hits;N Tracks",  10, 0, 10);
    fTracksPerTrack    = tfs->make<TH1D>("trackshitspertrack", ";CRT Tracks;N Tracks",  10, 0, 10);

    fTrackDist    = tfs->make<TH1D>("trackdist", ";True reco dist^2 (cm^2);N Tracks",  50, 0, 200);
    // Initial output
    std::cout<<"----------------- CRT Track Reco Ana Module -------------------"<<std::endl;

    // Take a position that is known to be the center of a strip


  }// CRTTrackRecoAna::beginJob()

  void CRTTrackRecoAna::analyze(const art::Event& event)
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
    std::map<int,TVector3> *partXYZ = new std::map<int,TVector3>[nTaggers];

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
    art::Handle< std::vector<crt::CRTHit>> crtListHandle;
    std::vector<art::Ptr<crt::CRTHit> > crtList;
    if (event.getByLabel(fCRTHitModuleLabel, crtListHandle))
      art::fill_ptr_vector(crtList, crtListHandle); 

    if(fVerbose) std::cout<<"Number of CRT hits = "<<crtList.size()<<std::endl;

    std::map<int, RecoTruth> truthMatch;

    std::map<art::Ptr<crt::CRTHit>, int> hitIds;
    for(size_t i = 0; i < crtList.size(); i++){
      hitIds[crtList[i]] = i;

      std::vector<int> ids = CRTTruthMatchUtils::AllTrueIds(crtListHandle, event, fCRTHitModuleLabel, fCRTModuleLabel, i);
      for(int& ID : ids){
        truthMatch[ID].crtHits.push_back(*crtList[i]);
      }
    }

    // Fill a vector of pairs of time and width direction for each CRT plane
    // The y crossing point of z planes and z crossing point of y planes would be constant

    std::vector<std::vector<art::Ptr<crt::CRTHit>>> crtTzeroVect = trackAlg.CreateCRTTzeros(crtList);
    
    int nTracks = 0;
    int nCompTracks = 0;
    int nIncTracks = 0;

    // Loop over tzeros
    for(size_t i = 0; i<crtTzeroVect.size(); i++){

      //loop over hits for this tzero, sort by tagger
      std::map<std::string, std::vector<art::Ptr<crt::CRTHit>>> hits;
      for (size_t ah = 0; ah< crtTzeroVect[i].size(); ++ah){        
        std::string ip = crtTzeroVect[i][ah]->tagger;       
        hits[ip].push_back(crtTzeroVect[i][ah]);
      } // loop over hits
      
      //loop over planes and calculate average hits
      std::vector<std::pair<crt::CRTHit, std::vector<int>>> allHits;
      for (auto &keyVal : hits){
        std::string ip = keyVal.first;
        std::vector<std::pair<crt::CRTHit, std::vector<int>>> ahits = trackAlg.AverageHits(hits[ip], hitIds);
        if(fUseTopPlane && ip == "volTaggerTopHigh_0"){ 
          allHits.insert(allHits.end(), ahits.begin(), ahits.end());
        }
        else if(ip != "volTaggerTopHigh_0"){ 
          allHits.insert(allHits.end(), ahits.begin(), ahits.end());
        }
      }

      // Do truth matching
      for(auto const& allHit : allHits){
        std::vector<int> ids;
        for(auto const& hit_i : allHit.second){
          std::vector<int> trueIds = CRTTruthMatchUtils::AllTrueIds(crtListHandle, event, fCRTHitModuleLabel, fCRTModuleLabel, hit_i);
          ids.insert(ids.end(), trueIds.begin(), trueIds.end());
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        for(int& ID : ids){
          truthMatch[ID].crtAveHits.push_back(allHit.first);
        }
      }
        

      std::vector<std::pair<crt::CRTTrack, std::vector<int>>> tracks = trackAlg.CreateTracks(allHits);
      nTracks += tracks.size();
      for(const auto& track : tracks){
        if(track.first.complete) nCompTracks++;
        else nIncTracks++;

        std::vector<int> ids;
        for(auto const& hit_i : track.second){
          std::vector<int> trueIds = CRTTruthMatchUtils::AllTrueIds(crtListHandle, event, fCRTHitModuleLabel, fCRTModuleLabel, hit_i);
          ids.insert(ids.end(), trueIds.begin(), trueIds.end());
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        for(int& ID : ids){
          truthMatch[ID].crtTracks.push_back(track.first);
        }

      }
    }

    if(fVerbose){
      std::cout<<"Number of hits = "<<crtList.size()<<std::endl
               <<"Number of T zero = "<<crtTzeroVect.size()<<std::endl
               <<"Number of tracks = "<<nTracks<<std::endl
               <<"Number of complete tracks = "<<nCompTracks<<std::endl
               <<"Number of incomplete tracks = "<<nIncTracks<<std::endl;
    }

    if(fPlot) DrawTrueTracks(particles, truthMatch, true, false, true, false, fPlotTrackID);

    for(auto const& particle : particles){
      int partId = particle.TrackId();

      RecoTruth rt = truthMatch[partId];
      fCrtHitsPerTrack->Fill(rt.crtHits.size());
      fAveHitsPerTrack->Fill(rt.crtAveHits.size());
      fTracksPerTrack->Fill(rt.crtTracks.size());

      if(rt.crtTracks.size()>0) nMatchTrack++;

      if(rt.crtTracks.size()>0){
        double minxres = 999999;
        double minipt = 1;

        for(size_t i = 0; i < rt.crtTracks.size(); i++){
          double sx = rt.crtTracks[i].x1_pos;
          double sy = rt.crtTracks[i].y1_pos;
          double sz = rt.crtTracks[i].z1_pos;
          double ex = rt.crtTracks[i].x2_pos;
          double ey = rt.crtTracks[i].y2_pos;
          double ez = rt.crtTracks[i].z2_pos;
          double p0 = sx - (sz*(sx-ex)/(sz-ez));
          double p1 = (sx-ex)/(sz-ez);
          double p2 = sy - (sz*(sy-ey)/(sz-ez));
          double p3 = (sy-ey)/(sz-ez);
          TVector3 x0(p0, p2, 0.);
          TVector3 x1(p0+p1, p2+p3, 1.);
          TVector3 u = (x1-x0).Unit();

          int nTraj = particle.NumberTrajectoryPoints();
          int ipt = 0;
          double xres = 0;
          for(int j = 0; j < nTraj; j++){
            if (particle.Vx(j) >= crtPlanes[0] && particle.Vx(j) <= crtPlanes[3] && particle.Vy(j) >= crtPlanes[4] && particle.Vy(j) <= -crtPlanes[4] && particle.Vz(j) >= crtPlanes[10] && particle.Vz(j) <= crtPlanes[13]){
              TVector3 xp(particle.Vx(j), particle.Vy(j), particle.Vz(j));
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
        std::cout<<"\nParticle "<<partId<<": Start("<<particle.Vx()<<","<<particle.Vy()<<","<<particle.Vz()<<") End("<<particle.EndX()<<","<<particle.EndY()<<","<<particle.EndZ()<<") time = "<<particle.T()*1e-3<<" us \nTrue crossing points:\n";
        for(size_t tag_i = 0; tag_i < nTaggers; tag_i++){
          if(partXYZ[tag_i].find(partId)!=partXYZ[tag_i].end()){
            std::cout<<"Tagger "<<tag_i<<": Coordinates = ("<<partXYZ[tag_i][partId].X()<<", "<<partXYZ[tag_i][partId].Y()<<", "<<partXYZ[tag_i][partId].Z()<<")\n";
          }
        }

        RecoTruth rt = truthMatch[partId];
        std::cout<<"->Reco tracks:\n";
        for(size_t trk_i = 0; trk_i < rt.crtTracks.size(); trk_i++){
          crt::CRTTrack tr = rt.crtTracks[trk_i];
          std::cout<<"  Start = ("<<tr.x1_pos<<","<<tr.y1_pos<<","<<tr.z1_pos<<") End = ("<<tr.x2_pos<<","<<tr.y2_pos<<","<<tr.z2_pos<<") time = "<<tr.ts1_ns*1e-3<<" us\n";
        }

        std::cout<<"-->Average hits:\n";
        for(size_t hit_i = 0; hit_i < rt.crtAveHits.size(); hit_i++){
          crt::CRTHit ah = rt.crtAveHits[hit_i];
          std::cout<<"   "<<ah.tagger<<": Coordinates = ("<<ah.x_pos<<","<<ah.y_pos<<","<<ah.z_pos<<") time = "<<ah.ts1_ns*1e-3<<" us\n";
        }

        std::cout<<"--->Reco crossing points:\n";
        for(size_t hit_i = 0; hit_i < rt.crtHits.size(); hit_i++){
          crt::CRTHit ht = rt.crtHits[hit_i];
          std::cout<<"    "<<ht.tagger<<": Coordinates = ("<<ht.x_pos<<", "<<ht.y_pos<<", "<<ht.z_pos<<") time = "<<ht.ts1_ns*1e-3<<" us\n";
        }
      }
    }

    delete[] partXYZ;

  } // CRTTrackRecoAna::analyze()

  void CRTTrackRecoAna::endJob(){

    std::cout<<"=========================== GENERAL INFORMATION ===========================\n"
             <<"Number of through-going muons = "<<nParts
             <<"\nNumber of muons that pass through 1 plane = "<<nSinglePlane
             <<"\nNumber of muons with at least 1 matched track = "<<nMatchTrack
             <<"\nNumber of muons that cross 2 strips on 2 planes = "<<nTwoStripPlanes<<"\n";


  } // CRTTrackRecoAna::endJob()


  void CRTTrackRecoAna::DrawTrueTracks(const std::vector<simb::MCParticle>& particle, std::map<int, RecoTruth> truthMatch, bool truth, bool hits, bool tracks, bool tpc, int ind){
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
      if(hits){
        // Plot the hits
        for(size_t j = 0; j < rt.crtHits.size(); j++){
          crt::CRTHit ht = rt.crtHits[j];
          // Get the limits
          double rmin[3] = {ht.x_pos-ht.x_err,ht.y_pos-ht.y_err,ht.z_pos-ht.z_err};
          double rmax[3] = {ht.x_pos+ht.x_err,ht.y_pos+ht.y_err,ht.z_pos+ht.z_err};
          truthAlg.DrawCube(c1, rmin, rmax, 3);
        }
      }
      if(tracks){
        // Plot the tracks
        for(size_t j = 0; j < rt.crtTracks.size(); j++){
          // Get the start and end points
          crt::CRTTrack tr = rt.crtTracks[j];
          crttrack[ncrtTracks] = new TPolyLine3D(2);
          crttrack[ncrtTracks]->SetPoint(0, tr.x1_pos, tr.y1_pos, tr.z1_pos);
          crttrack[ncrtTracks]->SetPoint(1, tr.x2_pos, tr.y2_pos, tr.z2_pos);
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


  DEFINE_ART_MODULE(CRTTrackRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


