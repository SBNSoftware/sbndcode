////////////////////////////////////////////////////////////////////////
// Class:       DirtTree
// Module Type: analyzer
// File:        DirtTree_module.cc
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
#include "larcoreobj/SummaryData/POTSummary.h"

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

  class DirtTree : public art::EDAnalyzer {
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
    explicit DirtTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;

    // Called for every sub run
    virtual void beginSubRun(const art::SubRun& subrun) override;
 
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

    // Tree (One entry per dirt particle)
    TTree *fNuTree;

    //Neutrino tree parameters
    int pdg;
    bool is_cc;
    bool in_tpc;
    bool in_cryo;
    bool in_crt;
    double vtx_x;
    double vtx_y;
    double vtx_z;
    double edep_tpc;
    int n_parts;

    TTree *fEventTree;

    int evt;
    bool nu_trigger;
    bool dirt_trigger;
    int n_nu;

    TTree *fPotTree;

    double pot;

  }; // class DirtTree


  // Constructor
  DirtTree::DirtTree(Parameters const& config)
    : EDAnalyzer(config)
    , fGenModuleLabel       (config().GenModuleLabel())
    , fSimModuleLabel       (config().SimModuleLabel())
    , fVerbose              (config().Verbose())
  {

  } // DirtTree()


  void DirtTree::beginJob()
  {

    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fNuTree = tfs->make<TTree>("nus", "nus");
    fNuTree->Branch("pdg",              &pdg);
    fNuTree->Branch("is_cc",            &is_cc);
    fNuTree->Branch("in_tpc",           &in_tpc);
    fNuTree->Branch("in_cryo",          &in_cryo);
    fNuTree->Branch("in_crt",           &in_crt);
    fNuTree->Branch("vtx_x",            &vtx_x);
    fNuTree->Branch("vtx_y",            &vtx_y);
    fNuTree->Branch("vtx_z",            &vtx_z);
    fNuTree->Branch("edep_tpc",         &edep_tpc);
    fNuTree->Branch("n_parts",         &n_parts);

    fEventTree = tfs->make<TTree>("events", "events");
    fEventTree->Branch("evt",  &evt);
    fEventTree->Branch("nu_trigger",  &nu_trigger);
    fEventTree->Branch("dirt_trigger",  &dirt_trigger);
    fEventTree->Branch("n_nu",  &n_nu);

    fPotTree = tfs->make<TTree>("pots", "pots");
    fPotTree->Branch("pot",  &pot);

    // Initial output
    if(fVerbose) std::cout<<"----------------- PDS Ana Module -------------------"<<std::endl;

  }// DirtTree::beginJob()


  // Called for every sub run
  void DirtTree::beginSubRun(const art::SubRun& subrun){

    art::Handle< sumdata::POTSummary > potHandle;
    subrun.getByLabel(fGenModuleLabel, potHandle);
    const sumdata::POTSummary& potSum = (*potHandle);
    pot = potSum.totpot;

    fPotTree->Fill();

    return;
  } // XSecTree::beginSubRun()


  void DirtTree::analyze(const art::Event& event)
  {
    ResetEventVars();

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    evt = event.id().event();

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

    n_nu = mctruthList.size();

    for (size_t i = 0; i < mctruthList.size(); i++){
      art::Ptr<simb::MCTruth> truth = mctruthList.at(i);

      ResetVars();
      if(truth->Origin() != simb::kBeamNeutrino) continue;

      pdg = truth->GetNeutrino().Nu().PdgCode();
      geo::Point_t vertex {truth->GetNeutrino().Nu().Vx(), 
                           truth->GetNeutrino().Nu().Vy(), 
                           truth->GetNeutrino().Nu().Vz()};

      vtx_x = vertex.X();
      vtx_y = vertex.Y();
      vtx_z = vertex.Z();
      is_cc = truth->GetNeutrino().CCNC() == simb::kCC;
      in_cryo = fTpcGeo.InCryo(vertex);
      in_crt = fCrtGeo.IsInsideCRT(vertex);
      in_tpc = fTpcGeo.InFiducial(vertex, 0.);

      std::vector<const simb::MCParticle*> particles = pi_serv->MCTruthToParticles_Ps(truth);
      n_parts = particles.size();
      for(size_t j = 0; j < particles.size(); j++){
        if(particles.at(j)->StatusCode() != 1) continue;
        int p_pdg = std::abs(particles.at(j)->PdgCode());
        if(!(p_pdg == 11 || p_pdg == 13 || p_pdg == 111 || p_pdg == 211 || p_pdg == 2212)) continue;
        edep_tpc += fTpcGeo.EDep(*particles.at(j));
      }

      fNuTree->Fill();
      if(in_tpc && edep_tpc > 0.02) nu_trigger = true;
      if(!in_tpc && edep_tpc > 0.02) dirt_trigger = true;

    }


    fEventTree->Fill();

  } // DirtTree::analyze()


  void DirtTree::endJob(){

  } // DirtTree::endJob()

  // Reset the tree variables
  void DirtTree::ResetVars(){
    pdg = -99999;
    is_cc = false;
    in_tpc = false;
    in_cryo = false;
    in_crt = false;
    vtx_x = -99999;
    vtx_y = -99999;
    vtx_z = -99999;
    edep_tpc = 0;
    n_parts = 0;
  }

  void DirtTree::ResetEventVars(){
    nu_trigger = false;
    dirt_trigger = false;
    n_nu = 0;
  }

  DEFINE_ART_MODULE(DirtTree)
} // namespace sbnd

