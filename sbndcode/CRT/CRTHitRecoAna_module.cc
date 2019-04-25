////////////////////////////////////////////////////////////////////////
// Class:       CRTHitRecoAna
// Module Type: analyzer
// File:        CRTHitRecoAna_module.cc
//
// Analysis module for evaluating CRT hit reconstruction on through going
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
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"
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

  struct RecoTruth{
    std::vector<std::pair<std::pair<std::string,unsigned>,art::Ptr<crt::CRTData>>> sipmHits;
    std::vector<CRTStrip> stripHits;
    std::vector<crt::CRTHit> crtHits;
  };

  class CRTHitRecoAna : public art::EDAnalyzer {
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

      fhicl::Table<CRTHitRecoAlg::Config> HitAlg {
        Name("HitAlg"),
      };

    }; // Config
 
    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CRTHitRecoAna(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    void DrawTrueTracks(const std::vector<simb::MCParticle>& particle, std::map<int, RecoTruth> truthMatch, bool truth, bool strips, bool hits, bool tpc, int ind);

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fCRTModuleLabel;      ///< name of CRT producer
    bool          fVerbose;             ///< print information about what's going on
    bool          fVeryVerbose;         ///< print more information about what's going on
    bool          fPlot;                ///< plot tracks
    int           fPlotTrackID;         ///< id of track to plot
    bool          fUseReadoutWindow;    ///< Only reconstruct hits within readout window
    
    // n-tuples
    TH2D* fTagXYZResolution[7];
    TH1D* fTagXResolution[7];
    TH1D* fTagYResolution[7];

    TH1D* fSipmDist;

    // Other variables shared between different methods.
    geo::GeometryCore const* fGeometryService;                 ///< pointer to Geometry provider
    detinfo::DetectorProperties const* fDetectorProperties;    ///< pointer to detector properties provider
    detinfo::DetectorClocks const* fDetectorClocks;            ///< pointer to detector clocks provider
    detinfo::ElecClock fTrigClock;

    CRTTruthRecoAlg truthAlg;
    CRTHitRecoAlg hitAlg;

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

  }; // class CRTHitRecoAna

  // Constructor
  CRTHitRecoAna::CRTHitRecoAna(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fCRTModuleLabel       (config().CRTModuleLabel())
    , fVerbose              (config().Verbose())
    , fVeryVerbose          (config().VeryVerbose())
    , fPlot                 (config().Plot())
    , fPlotTrackID          (config().PlotTrackID())
    , fUseReadoutWindow     (config().HitAlg().UseReadoutWindow())
    , truthAlg()
    , hitAlg(config().HitAlg())
  {
    // Get a pointer to the fGeometryServiceetry service provider
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 
    fTrigClock = fDetectorClocks->TriggerClock();
  }

  void CRTHitRecoAna::beginJob()
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

    fSipmDist    = tfs->make<TH1D>("sipmdist", ";Distance from SiPM (cm);N Tracks",  50, -3, 13);
    // Initial output
    std::cout<<"----------------- CRT Hit Reco Ana Module -------------------"<<std::endl;

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

  }// CRTHitRecoAna::beginJob()

  void CRTHitRecoAna::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    std::vector<uint8_t> tfeb_id = {0};
    std::map<uint8_t, std::vector<std::pair<int,float>>> tpesmap;
    tpesmap[0] = {std::make_pair(0,0)};

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
    //std::map<std::pair<std::string,unsigned>, std::vector<CRTStrip>> taggerStrips;
    std::map<std::pair<std::string,unsigned>, std::vector<CRTStrip>> taggerStrips = hitAlg.CreateTaggerStrips(crtList);
    
    for(auto const& taggerStrip : taggerStrips){

      for(auto const& strip : taggerStrip.second){

        size_t ind = strip.dataID;

        int id1 = CRTTruthMatchUtils::TrueIdFromTotalEnergy(crtListHandle, event, fCRTModuleLabel, ind);
        int id2 = CRTTruthMatchUtils::TrueIdFromTotalEnergy(crtListHandle, event, fCRTModuleLabel, ind+1); 

        truthMatch[id1].sipmHits.push_back(std::make_pair(taggerStrip.first, crtList[ind]));
        truthMatch[id2].sipmHits.push_back(std::make_pair(taggerStrip.first, crtList[ind+1]));

        truthMatch[id1].stripHits.push_back(strip);
        if(id1!=id2) truthMatch[id2].stripHits.push_back(strip);

        nSipms += 2;
        nStrips ++;
      }

    }

    std::vector<std::pair<crt::CRTHit, std::vector<int>>> crtHitPairs = hitAlg.CreateCRTHits(taggerStrips);

    std::vector<crt::CRTHit> crtHits;
    for(auto const& crtHitPair : crtHitPairs){
      crtHits.push_back(crtHitPair.first);
      nHits++;
      // If the PID matches calculate the resolution
      std::vector<int> ids;
      for(auto const& data_i : crtHitPair.second){
        ids.push_back(CRTTruthMatchUtils::TrueIdFromTotalEnergy(crtListHandle, event, fCRTModuleLabel, data_i));
      }
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

      int tag_i = nameToInd[crtHitPair.first.tagger];
      for(int& ID : ids){
        if(partXYZ[tag_i].find(ID)!=partXYZ[tag_i].end()){
          nMatchingHits++; 
          nMatchTag[tag_i]++;
          usedXYZ[tag_i][ID] = partXYZ[tag_i][ID];
          TVector3 mean(crtHitPair.first.x_pos, crtHitPair.first.y_pos, crtHitPair.first.z_pos);
          TVector3 dist = partXYZ[tag_i][ID] - mean;
          double distx = dist[widthCoord[tag_i*2]];
          double disty = dist[lenCoord[tag_i*2]];
          fTagXYZResolution[tag_i]->Fill(distx,disty); 
          fTagXResolution[tag_i]->Fill(distx); 
          fTagYResolution[tag_i]->Fill(disty);
        }
        else{ nNoTruth++; nNoTruthTag[tag_i]++; }
        truthMatch[ID].crtHits.push_back(crtHitPair.first);
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

    if(fVerbose){
      std::cout<<"Number of hits = "<<crtHits.size()<<std::endl;
    }

    if(fPlot) DrawTrueTracks(particles, truthMatch, true, false, true, false, fPlotTrackID);
    
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
        std::cout<<"--->Reco crossing points:\n";
        for(size_t hit_i = 0; hit_i < rt.crtHits.size(); hit_i++){
          crt::CRTHit ht = rt.crtHits[hit_i];
          std::cout<<"    "<<ht.tagger<<": Coordinates = ("<<ht.x_pos<<", "<<ht.y_pos<<", "<<ht.z_pos<<") time = "<<ht.ts1_ns*1e-3<<" us\n";
        }

        std::cout<<"---->Strip hits:\n";
        for(size_t hit_i = 0; hit_i < rt.stripHits.size(); hit_i++){
          CRTStrip sp = rt.stripHits[hit_i];
          std::cout<<"     "<<sp.tagger.first<<" ("<<sp.tagger.second<<"): time = "<<sp.t0<<"\n";
        }

        std::cout<<"----->SiPM hits:\n";
        for(size_t hit_i = 0; hit_i < rt.sipmHits.size(); hit_i++){
          std::pair<std::string,unsigned> tagger = rt.sipmHits[hit_i].first;
          art::Ptr<crt::CRTData> si = rt.sipmHits[hit_i].second;
          fTrigClock.SetTime(si->T0());
          double t1 = fTrigClock.Time(); // [us]
          std::cout<<"      "<<tagger.first<<" ("<<tagger.second<<"): Channel = "<<si->Channel()<<" time = "<<t1<<" us\n";
        }
      }
    }

    delete[] taggerXYZ;
    delete[] partXYZ;
    delete[] usedXYZ;

  } // CRTHitRecoAna::analyze()

  void CRTHitRecoAna::endJob(){

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

  } // CRTHitRecoAna::endJob()

    
  void CRTHitRecoAna::DrawTrueTracks(const std::vector<simb::MCParticle>& particle, std::map<int, RecoTruth> truthMatch, bool truth, bool strips, bool hits, bool tpc, int ind){
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
    if(truthMatch.find(ind) != truthMatch.end()){
      tMatch.clear();
      tMatch[ind] = truthMatch[ind];
    }
    for(auto &irt : tMatch){
      RecoTruth rt = irt.second;
      if(strips){
        // Plot the hit strips
        for(size_t j = 0; j < rt.stripHits.size(); j++){
          // Calculate the limits
          std::vector<double> limits = hitAlg.ChannelToLimits(rt.stripHits[j]);
          // Plot a rectangle
          double rmin[3] = {limits[0], limits[2], limits[4]};
          double rmax[3] = {limits[1], limits[3], limits[5]};
          truthAlg.DrawCube(c1, rmin, rmax, 2);
        }
      }
      if(hits){
        // Plot the hits
        for(size_t j = 0; j < rt.crtHits.size(); j++){
          crt::CRTHit ht = rt.crtHits[j];
          // Get the limits
          double rmin[3] = {ht.x_pos-ht.x_err, ht.y_pos-ht.y_err, ht.z_pos-ht.z_err};
          double rmax[3] = {ht.x_pos+ht.x_err, ht.y_pos+ht.y_err, ht.z_pos+ht.z_err};
          truthAlg.DrawCube(c1, rmin, rmax, 3);
        }
      }
    }

    c1->SaveAs("crtTagger.root");
  }


  DEFINE_ART_MODULE(CRTHitRecoAna)
} // namespace sbnd

// Back to local namespace.
namespace {

} // local namespace


