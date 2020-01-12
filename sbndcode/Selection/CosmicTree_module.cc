////////////////////////////////////////////////////////////////////////
// Class:       CosmicTree
// Module Type: analyzer
// File:        CosmicTree_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <algorithm>

namespace sbnd {

  class CosmicTree : public art::EDAnalyzer {
  public:

    // Describes configuration parameters of the module
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
 
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> GenModuleLabel {
        Name("GenModuleLabel"),
        Comment("tag of generator data product")
      };

      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of g4 data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

    }; // Inputs

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit CosmicTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    std::vector<double> OpFlashes(std::vector<double> optimes);

    // Reset variables in each loop
    void ResetVars();
    void ResetEventVars();

  private:

    // fcl file parameters
    art::InputTag fGenModuleLabel;      ///< name of generator producer
    art::InputTag fSimModuleLabel;      ///< name of g4 producer
    bool          fVerbose;             ///< print information about what's going on

    TPCGeoAlg fTpcGeo;
    CRTGeoAlg fCrtGeo;

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorProperties const* fDetectorProperties;

    // Tree (One entry per primary muon)
    TTree *fParticleTree;

    //Particle tree parameters
    bool cross_apa;
    bool cross_cpa;
    int pdg;
    double time;
    double vtx_x;
    double vtx_y;
    double vtx_z;
    double end_x;
    double end_y;
    double end_z;
    double vtx_x_tpc;
    double vtx_y_tpc;
    double vtx_z_tpc;
    double end_x_tpc;
    double end_y_tpc;
    double end_z_tpc;
    double vtx_x_crt;
    double vtx_y_crt;
    double vtx_z_crt;
    double end_x_crt;
    double end_y_crt;
    double end_z_crt;
    double length;
    double contained_length;
    double momentum;
    double theta;
    double phi;

    TTree *fEventTree;

    int evt;

  }; // class CosmicTree


  // Constructor
  CosmicTree::CosmicTree(Parameters const& config)
    : EDAnalyzer(config)
    , fGenModuleLabel       (config().GenModuleLabel())
    , fSimModuleLabel       (config().SimModuleLabel())
    , fVerbose              (config().Verbose())
  {

  } // CosmicTree()


  void CosmicTree::beginJob()
  {

    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fParticleTree = tfs->make<TTree>("particles", "particles");

    fParticleTree->Branch("cross_apa",        &cross_apa);
    fParticleTree->Branch("cross_cpa",        &cross_cpa);
    fParticleTree->Branch("pdg",              &pdg);
    fParticleTree->Branch("time",             &time);
    fParticleTree->Branch("vtx_x",            &vtx_x);
    fParticleTree->Branch("vtx_y",            &vtx_y);
    fParticleTree->Branch("vtx_z",            &vtx_z);
    fParticleTree->Branch("end_x",            &end_x);
    fParticleTree->Branch("end_y",            &end_y);
    fParticleTree->Branch("end_z",            &end_z);
    fParticleTree->Branch("vtx_x_tpc",        &vtx_x_tpc);
    fParticleTree->Branch("vtx_y_tpc",        &vtx_y_tpc);
    fParticleTree->Branch("vtx_z_tpc",        &vtx_z_tpc);
    fParticleTree->Branch("end_x_tpc",        &end_x_tpc);
    fParticleTree->Branch("end_y_tpc",        &end_y_tpc);
    fParticleTree->Branch("end_z_tpc",        &end_z_tpc);
    fParticleTree->Branch("vtx_x_crt",        &vtx_x_crt);
    fParticleTree->Branch("vtx_y_crt",        &vtx_y_crt);
    fParticleTree->Branch("vtx_z_crt",        &vtx_z_crt);
    fParticleTree->Branch("end_x_crt",        &end_x_crt);
    fParticleTree->Branch("end_y_crt",        &end_y_crt);
    fParticleTree->Branch("end_z_crt",        &end_z_crt);
    fParticleTree->Branch("length",           &length);
    fParticleTree->Branch("contained_length", &contained_length);
    fParticleTree->Branch("momentum",         &momentum);
    fParticleTree->Branch("theta",            &theta);
    fParticleTree->Branch("phi",              &phi);

    fEventTree = tfs->make<TTree>("events", "events");
    fEventTree->Branch("evt",  &evt);

    // Initial output
    if(fVerbose) std::cout<<"----------------- PDS Ana Module -------------------"<<std::endl;

  }// CosmicTree::beginJob()


  void CosmicTree::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    evt = event.id().event();
    fEventTree->Fill();

    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------

    // Get truth info and matching
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    art::Handle<std::vector<simb::MCTruth>> genHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthList;
    if(event.getByLabel(fGenModuleLabel, genHandle)) art::fill_ptr_vector(mctruthList, genHandle);

    //----------------------------------------------------------------------------------------------------------
    //                                        COSMIC MUON ANALYSIS
    //----------------------------------------------------------------------------------------------------------

    for (auto const& particle: (*particleHandle)){
      // Only interested in muons
      if(!(std::abs(particle.PdgCode()) == 13)) continue;
      // Only want primary particles
      if(particle.Mother() != 0) continue;
      // Only want stable particles (post fsi)
      if(particle.StatusCode() != 1) continue;
      // Only want particles that are inside the TPC
      if(!fTpcGeo.InVolume(particle)) continue;

      ResetVars();

      int id = particle.TrackId();
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(id);
      if(truth->Origin() != simb::kCosmicRay) continue;

      pdg = particle.PdgCode();

      cross_apa = fTpcGeo.CrossesApa(particle);
      cross_cpa = fTpcGeo.CrossesCpa(particle);

      time = particle.T(); //[ns]

      vtx_x = particle.Vx();
      vtx_y = particle.Vy();
      vtx_z = particle.Vz();

      end_x = particle.EndX();
      end_y = particle.EndY();
      end_z = particle.EndZ();

      length = particle.Trajectory().TotalLength();
      contained_length = fTpcGeo.TpcLength(particle);
      momentum = particle.P();
      std::pair<TVector3, TVector3> se = fTpcGeo.CrossingPoints(particle);
      theta = (se.second-se.first).Theta();
      phi = (se.second-se.first).Phi();

      vtx_x_tpc = se.first.X();
      vtx_y_tpc = se.first.Y();
      vtx_z_tpc = se.first.Z();

      end_x_tpc = se.second.X();
      end_y_tpc = se.second.Y();
      end_z_tpc = se.second.Z();

      for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
        std::string name = fCrtGeo.GetTagger(i).name;
        if(name.find("Front") != std::string::npos){
          geo::Point_t cross = fCrtGeo.TaggerCrossingPoint(name, particle);
          vtx_x_crt = cross.X();
          vtx_y_crt = cross.Y();
          vtx_z_crt = cross.Z();
        }
        if(name.find("Back") != std::string::npos){
          geo::Point_t cross = fCrtGeo.TaggerCrossingPoint(name, particle);
          end_x_crt = cross.X();
          end_y_crt = cross.Y();
          end_z_crt = cross.Z();
        }
      }

      fParticleTree->Fill();
    }

  } // CosmicTree::analyze()


  void CosmicTree::endJob(){

  } // CosmicTree::endJob()

  // Reset the tree variables
  void CosmicTree::ResetVars(){
    cross_apa = false;
    cross_cpa = false;
    pdg = -99999;
    time = -99999;
    vtx_x = -99999;
    vtx_y = -99999;
    vtx_z = -99999;
    end_x = -99999;
    end_y = -99999;
    end_z = -99999;
    vtx_x_tpc = -99999;
    vtx_y_tpc = -99999;
    vtx_z_tpc = -99999;
    end_x_tpc = -99999;
    end_y_tpc = -99999;
    end_z_tpc = -99999;
    vtx_x_crt = -99999;
    vtx_y_crt = -99999;
    vtx_z_crt = -99999;
    end_x_crt = -99999;
    end_y_crt = -99999;
    end_z_crt = -99999;
    length = -99999;
    contained_length = -99999;
    momentum = -99999;
    theta = -99999;
    phi = -99999;
    
  }

  DEFINE_ART_MODULE(CosmicTree)
} // namespace sbnd

