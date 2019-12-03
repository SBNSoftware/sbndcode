////////////////////////////////////////////////////////////////////////
// Class:       PdsTree
// Module Type: analyzer
// File:        PdsTree_module.cc
//
// Tom Brooks (tbrooks@fnal.gov)
////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbndcode/RecoUtils/RecoUtils.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.h"

// LArSoft includes
#include "lardataobj/RecoBase/OpHit.h"
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

  class PdsTree : public art::EDAnalyzer {
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

      fhicl::Atom<art::InputTag> PdsModuleLabel {
        Name("PdsModuleLabel"),
        Comment("tag of pds producer data product")
      };

      fhicl::Atom<bool> Verbose {
        Name("Verbose"),
        Comment("Print information about what's going on")
      };

    }; // Inputs

    using Parameters = art::EDAnalyzer::Table<Config>;
 
    // Constructor: configures module
    explicit PdsTree(Parameters const& config);
 
    // Called once, at start of the job
    virtual void beginJob() override;
 
    // Called once per event
    virtual void analyze(const art::Event& event) override;

    // Called once, at end of the job
    virtual void endJob() override;

    // Reset variables in each loop
    void ResetVars();

  private:

    // fcl file parameters
    art::InputTag fSimModuleLabel;      ///< name of detsim producer
    art::InputTag fPdsModuleLabel; ///< name of PDS track producer
    bool          fVerbose;             ///< print information about what's going on

    TPCGeoAlg fTpcGeo;

    opdet::sbndPDMapAlg fChannelMap; //map for photon detector types

    geo::GeometryCore const* fGeometryService;
    detinfo::DetectorProperties const* fDetectorProperties;

    std::vector<std::string> opdets {"pmt", "barepmt"};
    
    // Tree (One entry per primary muon)
    TTree *fParticleTree;

    //Particle tree parameters
    bool is_cosmic;         // True origin of PFP is cosmic
    bool is_nu;             // True origin of PFP is nu
    bool cross_apa;
    bool is_cc;             // True origin of PFP is CC nu
    int nu_pdg;
    int pdg;
    double time;
    double vtx_x;
    double vtx_y;
    double vtx_z;
    double end_x;
    double end_y;
    double end_z;
    double length;
    double contained_length;
    double momentum;
    double theta;
    double phi;
    double vtx_x_tpc;
    double vtx_y_tpc;
    double vtx_z_tpc;
    double end_x_tpc;
    double end_y_tpc;
    double end_z_tpc;
    double e_dep_tpc0;
    double e_dep_tpc1;

    std::map<std::string, int> n_ophits_tpc0;
    std::map<std::string, int> n_ophits_tpc1;
    std::map<std::string, double> ophit_pe_tpc0;
    std::map<std::string, double> ophit_pe_tpc1;
    std::map<std::string, double> ophit_area_tpc0;
    std::map<std::string, double> ophit_area_tpc1;
    std::map<std::string, double> ophit_amp_tpc0;
    std::map<std::string, double> ophit_amp_tpc1;
    std::map<std::string, double> ave_time_diff;
    std::map<std::string, double> time_std_dev;


  }; // class PdsTree


  // Constructor
  PdsTree::PdsTree(Parameters const& config)
    : EDAnalyzer(config)
    , fSimModuleLabel       (config().SimModuleLabel())
    , fPdsModuleLabel       (config().PdsModuleLabel())
    , fVerbose              (config().Verbose())
  {

  } // PdsTree()


  void PdsTree::beginJob()
  {

    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fParticleTree = tfs->make<TTree>("particles", "particles");

    fParticleTree->Branch("is_cosmic",        &is_cosmic);
    fParticleTree->Branch("is_nu",            &is_nu);
    fParticleTree->Branch("cross_apa",        &cross_apa);
    fParticleTree->Branch("is_cc",            &is_cc);
    fParticleTree->Branch("nu_pdg",           &nu_pdg);
    fParticleTree->Branch("pdg",              &pdg);
    fParticleTree->Branch("time",             &time);
    fParticleTree->Branch("vtx_x",            &vtx_x);
    fParticleTree->Branch("vtx_y",            &vtx_y);
    fParticleTree->Branch("vtx_z",            &vtx_z);
    fParticleTree->Branch("end_x",            &end_x);
    fParticleTree->Branch("end_y",            &end_y);
    fParticleTree->Branch("end_z",            &end_z);
    fParticleTree->Branch("length",           &length);
    fParticleTree->Branch("contained_length", &contained_length);
    fParticleTree->Branch("momentum",         &momentum);
    fParticleTree->Branch("theta",            &theta);
    fParticleTree->Branch("phi",              &phi);
    fParticleTree->Branch("vtx_x_tpc",        &vtx_x_tpc);
    fParticleTree->Branch("vtx_y_tpc",        &vtx_y_tpc);
    fParticleTree->Branch("vtx_z_tpc",        &vtx_z_tpc);
    fParticleTree->Branch("end_x_tpc",        &end_x_tpc);
    fParticleTree->Branch("end_y_tpc",        &end_y_tpc);
    fParticleTree->Branch("end_z_tpc",        &end_z_tpc);
    fParticleTree->Branch("e_dep_tpc0",       &e_dep_tpc0);
    fParticleTree->Branch("e_dep_tpc1",       &e_dep_tpc1);
    for(auto const& opdet : opdets){
      fParticleTree->Branch((opdet+"_n_ophits_tpc0").c_str(),    &n_ophits_tpc0[opdet]);
      fParticleTree->Branch((opdet+"_n_ophits_tpc1").c_str(),    &n_ophits_tpc1[opdet]);
      fParticleTree->Branch((opdet+"_ophit_pe_tpc0").c_str(),    &ophit_pe_tpc0[opdet]);
      fParticleTree->Branch((opdet+"_ophit_pe_tpc1").c_str(),    &ophit_pe_tpc1[opdet]);
      fParticleTree->Branch((opdet+"_ophit_area_tpc0").c_str(),  &ophit_area_tpc0[opdet]);
      fParticleTree->Branch((opdet+"_ophit_area_tpc1").c_str(),  &ophit_area_tpc1[opdet]);
      fParticleTree->Branch((opdet+"_ophit_amp_tpc0").c_str(),   &ophit_amp_tpc0[opdet]);
      fParticleTree->Branch((opdet+"_ophit_amp_tpc1").c_str(),   &ophit_amp_tpc1[opdet]);
      fParticleTree->Branch((opdet+"_ave_time_diff").c_str(),    &ave_time_diff[opdet]);
      fParticleTree->Branch((opdet+"_time_std_dev").c_str(),     &time_std_dev[opdet]);
    }

    // Initial output
    if(fVerbose) std::cout<<"----------------- PDS Ana Module -------------------"<<std::endl;

  }// PdsTree::beginJob()


  void PdsTree::analyze(const art::Event& event)
  {

    // Fetch basic event info
    if(fVerbose){
      std::cout<<"============================================"<<std::endl
               <<"Run = "<<event.run()<<", SubRun = "<<event.subRun()<<", Event = "<<event.id().event()<<std::endl
               <<"============================================"<<std::endl;
    }

    //----------------------------------------------------------------------------------------------------------
    //                                          GETTING PRODUCTS
    //----------------------------------------------------------------------------------------------------------

    // Get truth info and matching
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    // Retrieve all the truth info in the events
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimModuleLabel);

    // Get PDS handle
    auto pdsHandle = event.getValidHandle<std::vector<recob::OpHit>>(fPdsModuleLabel);

    //----------------------------------------------------------------------------------------------------------
    //                                        MUON PDS RECO ANALYSIS
    //----------------------------------------------------------------------------------------------------------


    std::map<int, simb::MCParticle> particles;
    for (auto const& particle: (*particleHandle)){
      // Only interested in muons
      if(!(std::abs(particle.PdgCode()) == 13)) continue;
      // Only want primary particles
      if(particle.Mother() != 0) continue;
      // Only want stable particles (post fsi)
      if(particle.StatusCode() != 1) continue;
      // Only want particles that are inside the TPC
      if(!fTpcGeo.InVolume(particle)) continue;
      double ptime = particle.T()/1e3;
      // PDS only simulated in this window
      if(ptime < -1250 || ptime > 2500) continue;
      int id = particle.TrackId();
      particles[id] = particle;
    }

    // Put them in a map for easier access
    for (auto const& part: particles){
      simb::MCParticle particle = part.second;

      ResetVars();

      pdg = particle.PdgCode();

      // True variables
      int id = particle.TrackId();
      art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(id);
      if(truth->Origin() == simb::kBeamNeutrino){ 
        is_nu = true;
        nu_pdg = truth->GetNeutrino().Nu().PdgCode();
        if(truth->GetNeutrino().CCNC() == simb::kCC) is_cc = true;
      }
      if(truth->Origin() == simb::kCosmicRay) is_cosmic = true;

      cross_apa = fTpcGeo.CrossesApa(particle);

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

      for(size_t i = 0; i < particle.NumberTrajectoryPoints() - 1; i++){
        geo::Point_t pos {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
        if(fTpcGeo.InFiducial(pos, 0)){ 
          if(pos.X() <= 0) e_dep_tpc0 += particle.E(i) - particle.E(i+1);
          else e_dep_tpc1 += particle.E(i) - particle.E(i+1);
        }
      }

      std::map<std::string, int> nhits;
      std::map<std::string, std::vector<double>> optimes;
      for(auto const& ophit : (*pdsHandle)){
        // Only look at PMTs
        std::string od = fChannelMap.pdName(ophit.OpChannel());
        if(fChannelMap.pdName(ophit.OpChannel()) != "pmt" && fChannelMap.pdName(ophit.OpChannel()) != "barepmt") continue;
        if(ophit.PeakTime() < (time/1e3 - 10) || ophit.PeakTime() > (time/1e3 + 10)) continue;
        ave_time_diff[od] += ophit.PeakTime() - time/1e3;
        optimes[od].push_back(ophit.PeakTime());
        nhits[od]++;
        // Only look at hits within 1 us of the true time, PeakTime() in [us]
        if(ophit.PeakTime() < (time/1e3) || ophit.PeakTime() > (time/1e3 + 5)) continue;
        // Work out what TPC detector is in odd = TPC1, even = TPC0
        if(ophit.OpChannel() % 2 == 0){
          n_ophits_tpc1[od]++;
          ophit_pe_tpc1[od] += ophit.PE();
          ophit_area_tpc1[od] += ophit.Area();
          ophit_amp_tpc1[od] += ophit.Amplitude();
        }
        else{
          n_ophits_tpc0[od]++;
          ophit_pe_tpc0[od] += ophit.PE();
          ophit_area_tpc0[od] += ophit.Area();
          ophit_amp_tpc0[od] += ophit.Amplitude();
        }
      }
      for(auto const& kv : nhits){
        // The mean time difference
        ave_time_diff[kv.first] /= nhits[kv.first];
        double time_mean = std::accumulate(optimes[kv.first].begin(), optimes[kv.first].end(), 0.)/optimes[kv.first].size();
        double std_dev = 0.;
        for(auto const& t : optimes[kv.first]){
          std_dev += std::pow(t - time_mean, 2.);
        }
        time_std_dev[kv.first] = std_dev/(optimes.size()-1);
      }

      fParticleTree->Fill();
    }

  } // PdsTree::analyze()


  void PdsTree::endJob(){

  } // PdsTree::endJob()

  // Reset the tree variables
  void PdsTree::ResetVars(){
    is_cosmic = false;
    is_nu = false;
    cross_apa = false;
    is_cc = false;
    nu_pdg = -99999;
    pdg = -99999;
    time = -99999;
    vtx_x = -99999;
    vtx_y = -99999;
    vtx_z = -99999;
    end_x = -99999;
    end_y = -99999;
    end_z = -99999;
    length = -99999;
    contained_length = -99999;
    momentum = -99999;
    theta = -99999;
    phi = -99999;
    e_dep_tpc0 = 0;
    e_dep_tpc1 = 0;
    for(auto const& opdet : opdets){
      n_ophits_tpc0[opdet] = 0;
      n_ophits_tpc1[opdet] = 0;
      ophit_pe_tpc0[opdet] = 0;
      ophit_pe_tpc1[opdet] = 0;
      ophit_area_tpc0[opdet] = 0;
      ophit_area_tpc1[opdet] = 0;
      ophit_amp_tpc0[opdet] = 0;
      ophit_amp_tpc1[opdet] = 0;
      ave_time_diff[opdet] = 0;
      time_std_dev[opdet] = 0;
    }
    
  }

  DEFINE_ART_MODULE(PdsTree)
} // namespace sbnd

