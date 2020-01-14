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
#include "sbndcode/OpDetSim/OpT0FinderTypes.h"
#include "sbndcode/CosmicId/Utils/CosmicIdUtils.h"

// LArSoft includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
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
#include "Pandora/PdgTable.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TVector3.h"
#include "TH1.h"
#include "TFile.h"

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
      fhicl::Atom<art::InputTag> GenModuleLabel {
        Name("GenModuleLabel"),
        Comment("tag of generator data product")
      };

      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
        Comment("tag of g4 data product")
      };

      fhicl::Atom<art::InputTag> PdsModuleLabel {
        Name("PdsModuleLabel"),
        Comment("tag of pds producer data product")
      };

      fhicl::Atom<art::InputTag> PandoraModuleLabel {
        Name("PandoraModuleLabel"),
        Comment("tag of pandora producer data product")
      };

      fhicl::Atom<art::InputTag> SpModuleLabel {
        Name("SpModuleLabel"),
        Comment("tag of spacepoint producer data product")
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

    std::vector<double> OpFlashes(std::vector<recob::OpHit> optimes);
    std::vector<double> OpVariables(std::vector<recob::OpHit> ophits, int tpc, double start_t, double end_t);
    double FlashScore(double x, double y, double z, std::vector<double> variables);

    // Reset variables in each loop
    void ResetVars();
    void ResetEventVars();
    void ResetPfpVars();

    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;

  private:

    // fcl file parameters
    art::InputTag fGenModuleLabel;      ///< name of generator producer
    art::InputTag fSimModuleLabel;      ///< name of g4 producer
    art::InputTag fPdsModuleLabel;      ///< name of PDS producer
    art::InputTag fPandoraModuleLabel;  ///< name of Pandora producer
    art::InputTag fSpModuleLabel;  ///< name of Pandora producer
    bool          fVerbose;             ///< print information about what's going on

    std::string fInputFilename;
    float fBeamWindowEnd, fBeamWindowStart;
    float fLightWindowEnd, fLightWindowStart;
    float fMinFlashPE;
    float fPEscale;
    float fTermThreshold;

    std::map<unsigned int, unsigned int> _pfpmap;

    std::vector<float> dysp, dzsp, rrsp, pesp, dymean, dzmean, rrmean, pemean;
    int rr_nbins, dy_nbins, dz_nbins, pe_nbins;

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
    double closest_flash_tpc0;
    double closest_flash_tpc1;

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
    std::map<std::string, double> ave_time_diff_pe;

    TTree *fEventTree;

    bool nu_tpc0;
    bool nu_tpc1;
    double beam_edep_tpc0;
    double beam_edep_tpc1;
    int n_flashes_tpc0;
    int n_flashes_tpc1;
    int n_fake_flashes_tpc0;
    int n_fake_flashes_tpc1;
    bool score_nu;
    std::map<std::string, int> n_beam_hits_tpc0;
    std::map<std::string, int> n_beam_hits_tpc1;
    std::map<std::string, double> n_beam_pe_tpc0;
    std::map<std::string, double> n_beam_pe_tpc1;

    TTree *fPfpTree;
    bool pfp_is_nu;
    double pfp_time;
    double pfp_score;
    double pfp_vtx_x;
    double pfp_dy;
    double pfp_dz;
    double pfp_rr;
    double pfp_pe;

    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

  }; // class PdsTree


  // Constructor
  PdsTree::PdsTree(Parameters const& config)
    : EDAnalyzer(config)
    , fGenModuleLabel       (config().GenModuleLabel())
    , fSimModuleLabel       (config().SimModuleLabel())
    , fPdsModuleLabel       (config().PdsModuleLabel())
    , fPandoraModuleLabel   (config().PandoraModuleLabel())
    , fSpModuleLabel        (config().SpModuleLabel())
    , fVerbose              (config().Verbose())
  {

    fBeamWindowStart =  -0.2;
    fBeamWindowEnd   =  2.0;  // in ns
    fInputFilename = "/sbnd/data/fm_metrics_sbnd.root";
    fLightWindowEnd = 0.09;
    fLightWindowStart = -0.01;
    fMinFlashPE = 0.;
    fPEscale = 100.;
    fTermThreshold = 30.;

    //read histograms and fill vectors for match score calculation
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename, fname);

    TFile *infile = new TFile(fInputFilename.c_str(), "READ");

    TH1 *temphisto = (TH1*)infile->Get("rrp1");
    rr_nbins = temphisto->GetNbinsX();
    if (rr_nbins<=0) {
      rr_nbins=1;
      rrmean.push_back(0);
      rrsp.push_back(0.001);
    }
    else {
      for (int ib = 1; ib <= rr_nbins; ++ib) {
        rrmean.push_back(temphisto->GetBinContent(ib));
        float tt = temphisto->GetBinError(ib);
        if(tt<=0) tt = 100.;
        rrsp.push_back(tt);
      }
    }

    temphisto = (TH1*)infile->Get("dyp1");
    dy_nbins = temphisto->GetNbinsX();
    if (dy_nbins<=0) {
      dy_nbins=1;
      dymean.push_back(0);
      dysp.push_back(0.001);
    }
    else {
      for (int ib = 1;ib <= dy_nbins; ++ib) {
        dymean.push_back(temphisto->GetBinContent(ib));
        float tt = temphisto->GetBinError(ib);
        if(tt<=0) tt = 100.;
        dysp.push_back(tt);
      }
    }
    
    temphisto = (TH1*)infile->Get("dzp1");
    dz_nbins = temphisto->GetNbinsX();
    if (dz_nbins<=0) {
      dz_nbins=1;
      dzmean.push_back(0);
      dzsp.push_back(0.001);
    }
    else {
      for (int ib = 1; ib <= dz_nbins; ++ib) {
        dzmean.push_back(temphisto->GetBinContent(ib));
        float tt = temphisto->GetBinError(ib);
        if(tt<=0) tt = 100.;
        dzsp.push_back(tt);
      }
    }

    temphisto = (TH1*)infile->Get("pep1");
    pe_nbins = temphisto->GetNbinsX();
    if (pe_nbins<=0) {
      pe_nbins=1;
      pemean.push_back(0);
      pesp.push_back(0.001);
    }
    else {
      for (int ib=1; ib <= pe_nbins; ++ib) {
        pemean.push_back(temphisto->GetBinContent(ib));
        float tt = temphisto->GetBinError(ib);
        if(tt<=0) tt = 100.;
        pesp.push_back(tt);
      }
    }
    
    infile->Close();

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
    fParticleTree->Branch("closest_flash_tpc0",       &closest_flash_tpc0);
    fParticleTree->Branch("closest_flash_tpc1",       &closest_flash_tpc1);
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
      fParticleTree->Branch((opdet+"_ave_time_diff_pe").c_str(), &ave_time_diff_pe[opdet]);
    }

    fEventTree = tfs->make<TTree>("events", "events");

    fEventTree->Branch("nu_tpc0",          &nu_tpc0);
    fEventTree->Branch("nu_tpc1",          &nu_tpc1);
    fEventTree->Branch("beam_edep_tpc0",   &beam_edep_tpc0);
    fEventTree->Branch("beam_edep_tpc1",   &beam_edep_tpc1);
    fEventTree->Branch("n_flashes_tpc0",        &n_flashes_tpc0);
    fEventTree->Branch("n_flashes_tpc1",        &n_flashes_tpc1);
    fEventTree->Branch("n_fake_flashes_tpc0",   &n_fake_flashes_tpc0);
    fEventTree->Branch("n_fake_flashes_tpc1",   &n_fake_flashes_tpc1);
    fEventTree->Branch("score_nu",          &score_nu);
    for(auto const& opdet : opdets){
      fEventTree->Branch((opdet+"_n_beam_hits_tpc0").c_str(),    &n_beam_hits_tpc0[opdet]);
      fEventTree->Branch((opdet+"_n_beam_hits_tpc1").c_str(),    &n_beam_hits_tpc1[opdet]);
      fEventTree->Branch((opdet+"_n_beam_pe_tpc0").c_str(),    &n_beam_pe_tpc0[opdet]);
      fEventTree->Branch((opdet+"_n_beam_pe_tpc1").c_str(),    &n_beam_pe_tpc1[opdet]);
    }

    fPfpTree = tfs->make<TTree>("pfps", "pfps");
    fPfpTree->Branch("is_nu",          &pfp_is_nu);
    fPfpTree->Branch("time",          &pfp_time);
    fPfpTree->Branch("score",          &pfp_score);
    fPfpTree->Branch("vtx_x",          &pfp_vtx_x);
    fPfpTree->Branch("dy",          &pfp_dy);
    fPfpTree->Branch("dz",           &pfp_dz);
    fPfpTree->Branch("rr",           &pfp_rr);
    fPfpTree->Branch("pe",           &pfp_pe);

    // Initial output
    if(fVerbose) std::cout<<"----------------- PDS Ana Module -------------------"<<std::endl;

  }// PdsTree::beginJob()


  void PdsTree::analyze(const art::Event& event)
  {

    ResetEventVars();

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

    art::Handle<std::vector<simb::MCTruth>> genHandle;
    std::vector<art::Ptr<simb::MCTruth>> mctruthList;
    if(event.getByLabel(fGenModuleLabel, genHandle)) art::fill_ptr_vector(mctruthList, genHandle);

    // Get PDS handle
    auto pdsHandle = event.getValidHandle<std::vector<recob::OpHit>>(fPdsModuleLabel);

    // grab PFParticles in event
    PFParticleHandle pfParticleHandle;
    event.getByLabel(fPandoraModuleLabel, pfParticleHandle);
    PFParticleIdMap pfParticleMap;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
    
    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfParticleHandle, event, fPandoraModuleLabel);
    art::FindManyP<recob::Track> pfp_track_assn_v(pfParticleHandle, event, fPandoraModuleLabel);

    // grab associated metadata
    auto const& spacepoint_h = event.getValidHandle<std::vector<recob::SpacePoint> >(fSpModuleLabel);
    art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, event, fSpModuleLabel);

    std::map<int, simb::MCParticle> particles;
    std::vector<simb::MCParticle> parts;
    for (auto const& particle: (*particleHandle)){
      double ptime = particle.T()/1e3;
      // PDS only simulated in this window
      if(ptime < -1250 || ptime > 2500) continue;
      parts.push_back(particle);
      // Only interested in muons
      //if(!(std::abs(particle.PdgCode()) == 13)) continue;
      // Only want primary particles
      if(particle.Mother() != 0) continue;
      // Only want stable particles (post fsi)
      if(particle.StatusCode() != 1) continue;
      // Only want particles that are inside the TPC
      if(!fTpcGeo.InVolume(particle)) continue;
      int id = particle.TrackId();
      particles[id] = particle;
    }

    //----------------------------------------------------------------------------------------------------------
    //                                        FLASH RECONSTRUCTION
    //----------------------------------------------------------------------------------------------------------

    // get flash time
    int nbins = 500 * (fBeamWindowEnd - fBeamWindowStart);
    TH1F *ophittime = new TH1F("ophittime", "ophittime", nbins, fBeamWindowStart, fBeamWindowEnd); // in us
    std::vector<recob::OpHit> ophs;
    for(auto const& oph : (*pdsHandle)){
      ophs.push_back(oph);
      if ( fChannelMap.pdType(oph.OpChannel(),"pmt")) {
        if ( (oph.PeakTime()>fBeamWindowStart) && (oph.PeakTime()< fBeamWindowEnd) ) {
	        ophittime->Fill(oph.PeakTime(), fPEscale*oph.PE());
        }
      }
    }
    auto ibin =  ophittime->GetMaximumBin();
    float flashtime = (ibin*0.002) + fBeamWindowStart;  // in us
    float lowedge = flashtime + fLightWindowStart;
    float highedge = flashtime + fLightWindowEnd;

/*
    std::vector<double> flash_x {0, 0};
    std::vector<double> flash_y {0, 0};
    std::vector<double> flash_z {0, 0};
    std::vector<double> flash_r {-99999, -99999};
    std::vector<double> flash_pe {0, 0};
*/
    std::vector<std::vector<double>> opvars;

    for (size_t it=0; it<fGeometryService->NTPC(); ++it) {
      opvars.push_back(OpVariables(ophs, it, lowedge, highedge));
      /*
      double PMTxyz[3];
	    double sum_Ay=0; double sum_Az=0;
	    double sum_Cy=0; double sum_Cz=0;
	    double sum_D=0;
      double unpe_tot = 0;
      double pe_tot = 0;
	    for(auto const& oph : (*pdsHandle)){
	      if ( !((oph.PeakTime()>lowedge) && (oph.PeakTime()< highedge)) ) continue;
	      fGeometryService->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
	      if ((it==0 && PMTxyz[0]>0) || (it==1 && PMTxyz[0]<0) ) continue;

	      if ( fChannelMap.pdType(oph.OpChannel(),"pmt")){
	 	      // Add up the position, weighting with PEs
          pe_tot  += oph.PE();
          flash_y[it] += oph.PE() * PMTxyz[1];
          flash_z[it] += oph.PE() * PMTxyz[2];
	 	      sum_Ay  += pow(oph.PE(),2.) * pow(PMTxyz[1],2.);
	 	      sum_Az  += pow(oph.PE(),2.) * pow(PMTxyz[2],2.);
	 	      sum_D   += pow(oph.PE(),2.);
	 	      sum_Cy  += pow(oph.PE(),2.) * PMTxyz[1];
	 	      sum_Cz  += pow(oph.PE(),2.) * PMTxyz[2];
        }
	      else if ( fChannelMap.pdType(oph.OpChannel(),"barepmt")){
          unpe_tot += oph.PE();
          flash_y[it] += oph.PE() * PMTxyz[1];
          flash_z[it] += oph.PE() * PMTxyz[2];
        }
	    }
	 
      if(pe_tot <= 0) continue;
	    flash_y[it] /= (pe_tot + unpe_tot);
	    flash_z[it] /= (pe_tot + unpe_tot);
	    flash_r[it] = sqrt((sum_Ay - pow(sum_Cy,2.)/sum_D + sum_Az - pow(sum_Cz,2.)/sum_D)/sum_D);
      flash_pe[it] = 100. * unpe_tot / pe_tot;

      std::vector<double> allowed_x_pe;
      for(int i = 0; i < pe_nbins; i++){
        if(flash_pe[it] >= pemean[i]-pesp[i] && flash_pe[it] <= pemean[i]+pesp[i]) allowed_x_pe.push_back(200. - i*200./nbins);
        if(i == 0 && flash_pe[it] < pemean[i]-pesp[i]) allowed_x_pe.push_back(0.);
        if(i == pe_nbins-1 && flash_pe[it] > pemean[i]+pesp[i]) allowed_x_pe.push_back(200.);
      }
      std::sort(allowed_x_pe.begin(), allowed_x_pe.end());

      std::vector<double> allowed_x_rr;
      for(int i = 0; i < rr_nbins; i++){
        if(flash_r[it] >= rrmean[i]-rrsp[i] && flash_r[it] <= rrmean[i]+rrsp[i]) allowed_x_rr.push_back(200. - i*200./nbins);
        if(i == 0 && flash_r[it] < rrmean[i]-rrsp[i]) allowed_x_rr.push_back(0.);
        if(i == rr_nbins-1 && flash_r[it] > rrmean[i]+rrsp[i]) allowed_x_rr.push_back(200.);
      }
      std::sort(allowed_x_rr.begin(), allowed_x_rr.end());

      if(allowed_x_pe.size()==0 || allowed_x_rr.size() ==0){ 
        std::cout<<"No valid X!\n";
        flash_x[it] = 100.;
      }
      else{
        double min_x = std::min(allowed_x_pe[0], allowed_x_rr[0]);
        double max_x = std::max(allowed_x_pe.back(), allowed_x_rr.back());
        flash_x[it] = (max_x + min_x)/2.;
      }
      */
    }
  
    double min_score = 99999;
    // Loop over pandora pfp particles, select primary particles identified as neutrinos
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it){
      const art::Ptr<recob::PFParticle> pParticle(it->second);
      if (!pParticle->IsPrimary()) continue;
      // Check if this particle is identified as the neutrino
      const int pdg(pParticle->PdgCode());
      const bool isNeutrino(std::abs(pdg) == pandora::NU_E || 
                            std::abs(pdg) == pandora::NU_MU || 
                            std::abs(pdg) == pandora::NU_TAU);
      //Find neutrino pfparticle
      if(!isNeutrino) continue;

      ResetPfpVars();

      // go through these pfparticles and fill info needed for matching
      int pfp_tpc = -1;
      double nuvtx_x = 0.;
      double nuvtx_y = 0.;
      double nuvtx_z = 0.;
      double norm = 0.;

      std::vector<art::Ptr<recob::Hit>> allHits;

      for (const size_t daughterId : pParticle->Daughters()) {
        art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);

        const std::vector< art::Ptr<recob::SpacePoint> >& SPs = pfp_spacepoint_assn_v.at(pDaughter.key());
        for (size_t sp=0; sp < SPs.size(); sp++) {
          const std::vector< art::Ptr<recob::Hit> > hits = spacepoint_hit_assn_v.at( SPs[sp].key() );

          for (size_t h=0; h < hits.size(); h++) {

            // Only use hits from the collection plane
            if (fGeometryService->SignalType(hits[h]->WireID()) != geo::kCollection) continue;

            // Add the charged point to the vector
            const auto &position(SPs[sp]->XYZ());
            const auto charge(hits[h]->Integral());

            int tpcindex = (hits[h]->WireID()).TPC;
            if(pfp_tpc == -1) pfp_tpc = tpcindex;
            if(tpcindex != pfp_tpc) pfp_tpc = -2;

            nuvtx_x += charge * position[0];
            nuvtx_y += charge * position[1];
            nuvtx_z += charge * position[2];
            norm += charge;

            allHits.push_back(hits[h]);

          } // for all hits associated to this spacepoint
        } // for all spacepoints
      } // for all pfp pointers

      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(allHits, false);
      if(particles.find(trueId) != particles.end()){
        art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(trueId);
        if(truth->Origin() == simb::kBeamNeutrino) pfp_is_nu = true;
        pfp_time = particles[trueId].T()*1e-3;
      }

      // No charge deposition in PFP, very strange
      if(norm <= 0) continue;
      // PFP in two TPCs, should have t0 tag already
      if(pfp_tpc <= -1) continue;
      // No flash in TPC
      if(opvars[pfp_tpc][2] == -99999) continue;

      nuvtx_x /= norm;
      nuvtx_y /= norm;
      nuvtx_z /= norm;
/*	
	    //      calculate match score here, put association on the event
	    float slice = 200. - abs(nuvtx_x);
      float drift_distance = 200.;
	    int y_bin = int(dy_nbins * (slice / drift_distance));
	    int z_bin = int(dz_nbins * (slice / drift_distance));
	    int r_bin = int(rr_nbins * (slice / drift_distance));
	    //int pe_bin = int(pe_nbins * (slice / drift_distance));
  
      double score = 0;
      if(dysp[y_bin] > 0) score += abs(abs(opvars[pfp_tpc][0] - nuvtx_y) - dymean[y_bin])/dysp[y_bin];
      if(dzsp[z_bin] > 0) score += abs(abs(opvars[pfp_tpc][1] - nuvtx_z) - dzmean[z_bin])/dzsp[z_bin];
      if(rrsp[r_bin] > 0) score += abs(opvars[pfp_tpc][2] - rrmean[r_bin])/rrsp[r_bin];
      //if(pesp[pe_bin] > 0) score += abs(flash_pe[pfp_tpc] - pemean[pe_bin])/pesp[pe_bin];

      pfp_score = score;*/
      pfp_score = FlashScore(nuvtx_x, nuvtx_y, nuvtx_z, opvars[pfp_tpc]);
/*
      pfp_dy = abs(abs(flash_y[pfp_tpc] - nuvtx_y)-dymean[y_bin])/dysp[y_bin];
      pfp_dz = abs(abs(flash_z[pfp_tpc] - nuvtx_z)-dzmean[z_bin])/dzsp[z_bin];
      pfp_rr = abs(flash_r[pfp_tpc] - rrmean[r_bin])/rrsp[r_bin];
      pfp_pe = abs(flash_pe[pfp_tpc] - pemean[pe_bin])/pesp[pe_bin];
*/
      if(pfp_score < min_score){
        min_score = pfp_score;
        if(pfp_is_nu) score_nu = true;
        else score_nu = false;
      }
      
      fPfpTree->Fill();

    } // over all PFparticles

    //----------------------------------------------------------------------------------------------------------
    //                                        MUON PDS RECO ANALYSIS
    //----------------------------------------------------------------------------------------------------------


    // Optical flash reconstruction for numuCC
    std::vector<recob::OpHit> ophits_tpc0;
    std::vector<recob::OpHit> ophits_tpc1;
    for(auto const& ophit : (*pdsHandle)){
      // Only look at PMTs
      std::string od = fChannelMap.pdName(ophit.OpChannel());
      if( od != "pmt" ) continue;
      // Work out what TPC detector is in odd = TPC1, even = TPC0
      double PMTxyz[3];
	    fGeometryService->OpDetGeoFromOpChannel(ophit.OpChannel()).GetCenter(PMTxyz);
      if(PMTxyz[0] > 0){ 
        ophits_tpc1.push_back(ophit);
        // Beam activity
        if(ophit.PeakTime() >= 0 && ophit.PeakTime() <= 1.7){ 
          n_beam_hits_tpc1[od]++;
          n_beam_pe_tpc1[od] += ophit.PE();
        }
      }
      else{ 
        ophits_tpc0.push_back(ophit);
        // Beam activity
        if(ophit.PeakTime() >= 0 && ophit.PeakTime() <= 1.7){ 
          n_beam_hits_tpc0[od]++;
          n_beam_pe_tpc0[od] += ophit.PE();
        }
      }
    }

    std::vector<double> opflashes_tpc0 = OpFlashes(ophits_tpc0);
    n_flashes_tpc0 = opflashes_tpc0.size();

    std::vector<double> opflashes_tpc1 = OpFlashes(ophits_tpc1);
    n_flashes_tpc1 = opflashes_tpc1.size();

    std::pair<std::vector<double>, std::vector<double>> fake_flashes = CosmicIdUtils::FakeTpcFlashes(parts);
    n_fake_flashes_tpc0 = fake_flashes.first.size();
    n_fake_flashes_tpc1 = fake_flashes.second.size();

    // Put them in a map for easier access
    for (auto const& part: particles){
      simb::MCParticle particle = part.second;

      ResetVars();

      pdg = particle.PdgCode();
      if(std::abs(pdg) != 13) continue;

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
          if(pos.X() <= 0){ 
            e_dep_tpc0 += particle.E(i) - particle.E(i+1);
            if(time >=0 && time <= 1600) beam_edep_tpc0 += particle.E(i) - particle.E(i+1);
          }
          else{ 
            e_dep_tpc1 += particle.E(i) - particle.E(i+1);
            if(time >=0 && time <= 1600) beam_edep_tpc1 += particle.E(i) - particle.E(i+1);
          }
        }
      }

      // Find the closest optical flash to the true time
      for(auto const& flash : opflashes_tpc0){
        if(std::abs(flash - (time/1e3)) < std::abs(closest_flash_tpc0))
          closest_flash_tpc0 = flash - (time/1e3);
      }
      for(auto const& flash : opflashes_tpc1){
        if(std::abs(flash - (time/1e3)) < std::abs(closest_flash_tpc1))
          closest_flash_tpc1 = flash - (time/1e3);
      }

      std::map<std::string, int> nhits;
      std::map<std::string, double> npe;
      std::map<std::string, std::vector<double>> optimes;
      for(auto const& ophit : (*pdsHandle)){
        // Only look at PMTs
        std::string od = fChannelMap.pdName(ophit.OpChannel());
        if( std::find(opdets.begin(), opdets.end(), od) == opdets.end()) continue;
        if(ophit.PeakTime() < (time/1e3 - 5) || ophit.PeakTime() > (time/1e3 + 5)) continue;
        ave_time_diff[od] += ophit.PeakTime() - time/1e3;
        ave_time_diff_pe[od] += (ophit.PeakTime() - time/1e3)*ophit.PE();
        optimes[od].push_back(ophit.PeakTime());
        nhits[od]++;
        npe[od] += ophit.PE();
        // Only look at hits within 1 us of the true time, PeakTime() in [us]
        if(ophit.PeakTime() < (time/1e3) || ophit.PeakTime() > (time/1e3 + 5)) continue;
        // Work out what TPC detector is in odd = TPC1, even = TPC0
        double PMTxyz[3];
	      fGeometryService->OpDetGeoFromOpChannel(ophit.OpChannel()).GetCenter(PMTxyz);
        if(PMTxyz[0] > 0){ 
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
        ave_time_diff_pe[kv.first] /= npe[kv.first];
        double time_mean = std::accumulate(optimes[kv.first].begin(), optimes[kv.first].end(), 0.)/optimes[kv.first].size();
        double std_dev = 0.;
        for(auto const& t : optimes[kv.first]){
          std_dev += std::pow(t - time_mean, 2.);
        }
        time_std_dev[kv.first] = std_dev/(optimes[kv.first].size()-1);
      }

      fParticleTree->Fill();
    }

    // Determine if there are neutrinos in the AV
    for (size_t i = 0; i < mctruthList.size(); i++){

      // Get the pointer to the MCTruth object
      art::Ptr<simb::MCTruth> truth = mctruthList.at(i);

      if(truth->Origin() != simb::kBeamNeutrino) continue;

      // Get truth info if numuCC in AV
      geo::Point_t vtx {truth->GetNeutrino().Nu().Vx(), 
                        truth->GetNeutrino().Nu().Vy(), 
                        truth->GetNeutrino().Nu().Vz()};
      if(!fTpcGeo.InFiducial(vtx, 0.)) continue;

      if(vtx.X() < 0) nu_tpc0 = true;
      if(vtx.X() > 0) nu_tpc1 = true;

    }

    fEventTree->Fill();

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
    closest_flash_tpc0 = -99999;
    closest_flash_tpc1 = -99999;
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
      ave_time_diff_pe[opdet] = 0;
    }
    
  }

  void PdsTree::ResetEventVars(){
    nu_tpc0 = false;
    nu_tpc1 = false;
    beam_edep_tpc0 = 0;
    beam_edep_tpc1 = 0;
    n_flashes_tpc0 = 0;
    n_flashes_tpc1 = 0;
    n_fake_flashes_tpc0 = 0;
    n_fake_flashes_tpc1 = 0;
    score_nu = false;
    for(auto const& opdet : opdets){
      n_beam_hits_tpc0[opdet] = 0;
      n_beam_hits_tpc1[opdet] = 0;
      n_beam_pe_tpc0[opdet] = 0;
      n_beam_pe_tpc1[opdet] = 0;
    }
  }

  void PdsTree::ResetPfpVars(){
    pfp_is_nu = false;
    pfp_time = -99999;
    pfp_score = -99999;
    pfp_dy = -99999;
    pfp_dz = -99999;
    pfp_rr = -99999;
    pfp_pe = -99999;
  }

  std::vector<double> PdsTree::OpFlashes(std::vector<recob::OpHit> optimes){
    std::vector<double> opflashes;

    // get flash time
    int nbins = 500 * (2500 - -1250);
    TH1F *ophittimes = new TH1F("ophittime", "ophittime", nbins, -1250, 2500); // in us
    for(size_t i = 0; i < optimes.size(); i++){
      if ( fChannelMap.pdType(optimes[i].OpChannel(), "pmt")) {
	      ophittimes->Fill(optimes[i].PeakTime(), optimes[i].PE());
      }
    }
    for(int i = 0; i < ophittimes->GetNbinsX(); i++){
      if(ophittimes->GetBinContent(i) > 0.2){
        opflashes.push_back(((double)i*0.002) - 1250 - 0.2);
        while(ophittimes->GetBinContent(i) > 0.2) i++;
      }
    }

    delete ophittimes;
    return opflashes;
  }

  std::vector<double> PdsTree::OpVariables(std::vector<recob::OpHit> ophits, int tpc, double start_t, double end_t){

    std::vector<double> variables {0, 0, -99999, -99999};

    double PMTxyz[3];
	  double sum_Ay=0; 
    double sum_Az=0;
	  double sum_Cy=0; 
    double sum_Cz=0;
	  double sum_D=0;
    double unpe_tot = 0;
    double pe_tot = 0;

	  for(auto const& oph : ophits){

      // Time window of flash
	    if ( !((oph.PeakTime() > start_t) && (oph.PeakTime() < end_t)) ) continue;

      // Check if in correct TPC
	    fGeometryService->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
	    if ((tpc==0 && PMTxyz[0]>0) || (tpc==1 && PMTxyz[0]<0) ) continue;

      // For TPB coated PMTs
	    if ( fChannelMap.pdType(oph.OpChannel(),"pmt")){
	 	    // Add up the position, weighting with PEs
        pe_tot  += oph.PE();
        variables[0] += oph.PE() * PMTxyz[1];
        variables[1] += oph.PE() * PMTxyz[2];
	 	    sum_Ay  += pow(oph.PE(),2.) * pow(PMTxyz[1],2.);
	 	    sum_Az  += pow(oph.PE(),2.) * pow(PMTxyz[2],2.);
	 	    sum_D   += pow(oph.PE(),2.);
	 	    sum_Cy  += pow(oph.PE(),2.) * PMTxyz[1];
	 	    sum_Cz  += pow(oph.PE(),2.) * PMTxyz[2];
      }
      // For uncoated PMTs
	    else if ( fChannelMap.pdType(oph.OpChannel(),"barepmt")){
        unpe_tot += oph.PE();
        variables[0] += oph.PE() * PMTxyz[1];
        variables[1] += oph.PE() * PMTxyz[2];
      }
	  }
	 
    // If no hits return null vector
    if(pe_tot <= 0) return variables;

    // Calculate charge weighted center
	  variables[0] /= (pe_tot + unpe_tot);
	  variables[1] /= (pe_tot + unpe_tot);
    // Calculate PE spread
	  variables[2] = sqrt((sum_Ay - pow(sum_Cy,2.)/sum_D + sum_Az - pow(sum_Cz,2.)/sum_D)/sum_D);
    // Calculate PE ratio
    variables[3] = 100. * unpe_tot / pe_tot;

    return variables;
  }

  double PdsTree::FlashScore(double x, double y, double z, std::vector<double> variables){

    //      calculate match score here, put association on the event
	  float slice = 200. - abs(x);
    float drift_distance = 200.;
	  int y_bin = int(dy_nbins * (slice / drift_distance));
	  int z_bin = int(dz_nbins * (slice / drift_distance));
	  int r_bin = int(rr_nbins * (slice / drift_distance));
	  int pe_bin = int(pe_nbins * (slice / drift_distance));
  
    double score = 0;
    if(dysp[y_bin] > 0) score += abs(abs(variables[0] - y) - dymean[y_bin])/dysp[y_bin];
    if(dzsp[z_bin] > 0) score += abs(abs(variables[1] - z) - dzmean[z_bin])/dzsp[z_bin];
    if(rrsp[r_bin] > 0) score += abs(variables[2] - rrmean[r_bin])/rrsp[r_bin];

    pfp_dy = abs(abs(variables[0] - y)-dymean[y_bin])/dysp[y_bin];
    pfp_dz = abs(abs(variables[1] - z)-dzmean[z_bin])/dzsp[z_bin];
    pfp_rr = abs(variables[2] - rrmean[r_bin])/rrsp[r_bin];
    pfp_pe = abs(variables[3] - pemean[pe_bin])/pesp[pe_bin];

	  return score;

  }

  void PdsTree::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap){
      for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
          const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
          if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
              std::cout << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!" <<"\n";
          }
      }
  }
  
  DEFINE_ART_MODULE(PdsTree)
} // namespace sbnd

