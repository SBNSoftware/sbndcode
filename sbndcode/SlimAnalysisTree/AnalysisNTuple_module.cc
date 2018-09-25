////////////////////////////////////////////////////////////////////////
// Class:       AnalysisNTuple
// Plugin Type: analyzer (art v2_08_04)
// File:        AnalysisNTuple_module.cc
//
// Generated at Tue Jan 23 18:32:51 2018 by Rhiannon Jones using 
// cetskelgen from cetlib version v3_01_01.
//
////////////////////////////////////////////////////////////////////////
//
// This module constructs a root file which contains 
//  MCParticle, 
//  Reco Track, 
//  Reco Shower,
//  True/Reco event information 
//  for use in topology-based analyses
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include <sstream>
#include <cmath>
#include <ctime>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"

#include <typeinfo>

namespace sbnd {
  class AnalysisNTuple;
}


class sbnd::AnalysisNTuple : public art::EDAnalyzer {

public:
 
  struct Inputs {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<float> detectorLengthX{
      Name("DetectorLengthX"),
      Comment("Detetector length in the x direction")
    };
    fhicl::Atom<float> detectorLengthY{
      Name("DetectorLengthY"),
      Comment("Detetector length in the y direction")
    };
    fhicl::Atom<float> detectorLengthZ{
      Name("DetectorLengthZ"),
      Comment("Detetector length in the z direction")
    };
    fhicl::Atom<float> coordinateOffsetX{
      Name("CoordinateOffsetX"),
      Comment("Shift in the x direction of the origin")
    };
    fhicl::Atom<float> coordinateOffsetY{
      Name("CoordinateOffsetY"),
      Comment("Shift in the y direction of the origin")
    };
    fhicl::Atom<float> coordinateOffsetZ{
      Name("CoordinateOffsetZ"),
      Comment("Shift in the z direction of the origin")
    };
    fhicl::Atom<float> selectedBorderX{
      Name("SelectedBorderX"),
      Comment("Shift in the x direction of the origin")
    };
    fhicl::Atom<float> selectedBorderY{
      Name("SelectedBorderY"),
      Comment("Shift in the y direction of the origin")
    };
    fhicl::Atom<float> selectedBorderZ{
      Name("SelectedBorderZ"),
      Comment("Shift in the z direction of the origin")
    };
    fhicl::Atom<art::InputTag> generator_label{
      Name("TruthLabel"),
      Comment("Generator label")
    };
    fhicl::Atom<art::InputTag> geant_label{
      Name("G4Label"),
      Comment("G4 label")
    };
    fhicl::Atom<art::InputTag> pandora_label{
      Name("PandoraLabel"),
      Comment("Pandora label")
    };
    fhicl::Atom<art::InputTag> reco_track_label{
      Name("RecoTrackLabel"),
      Comment("Pandora track label")
    };
    fhicl::Atom<art::InputTag> reco_shower_label{
      Name("RecoShowerLabel"),
      Comment("EM shower label")
    };
    fhicl::Atom<art::InputTag> reco_track_calorimetry_label{
      Name("RecoTrackCalorimetryLabel"),
      Comment("Track calorimetry label")
    };
    fhicl::Atom<art::InputTag> reco_track_particleid_label{
      Name("RecoTrackParticleIDLabel"),
      Comment("Track PID label")
    };
    fhicl::Atom<art::InputTag> hit_label{
      Name("HitLabel"),
      Comment("Hit label")
    };
  };
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
  
    fhicl::Table<AnalysisNTuple::Inputs> inputs {
      Name("inputs"), 
    };
    fhicl::Table<trkf::TrajectoryMCSFitter::Config> fitter {
      Name("fitter"),
    };
  };

  using Parameters = art::EDAnalyzer::Table<Config>;
  
  explicit AnalysisNTuple(Parameters const & config);

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalysisNTuple(AnalysisNTuple const &) = delete;
  AnalysisNTuple(AnalysisNTuple &&) = delete;
  AnalysisNTuple & operator = (AnalysisNTuple const &) = delete;
  AnalysisNTuple & operator = (AnalysisNTuple &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Geometry 
  float m_detectorLengthX;
  float m_detectorLengthY;
  float m_detectorLengthZ;
  float m_coordinateOffsetX;
  float m_coordinateOffsetY;
  float m_coordinateOffsetZ;
  float m_selectedBorderX;
  float m_selectedBorderY;
  float m_selectedBorderZ;

  // Handle labels 
  art::InputTag m_generator_label;
  art::InputTag m_geant_label;
  art::InputTag m_pandora_label;
  art::InputTag m_reco_track_label;
  art::InputTag m_reco_shower_label;
  art::InputTag m_reco_track_calorimetry_label;
  art::InputTag m_reco_track_particleid_label;
  art::InputTag m_hit_label;
  
  // Momentum fitters
  trkf::TrajectoryMCSFitter m_mcs_fitter;
  trkf::TrackMomentumCalculator m_range_fitter;

  // Counters
  unsigned int pfparticle, primary_pfparticle;

  // ROOT
  TTree *event_tree, *mcparticle_tree, *recotrack_tree, *recoshower_tree, *mcsstudy_tree;

  // Variables associated with event_tree
  // Universal
  int event_id;
  std::time_t time_now;

  // Truth
  int t_nu_pdgcode, t_interaction, t_scatter;
  int t_particles, t_protons, t_neutrons, t_muons, t_charged_pions, t_kaons, t_neutral_pions, t_photons, t_electrons;
  double t_inv_mass, t_nu_lepton_angle, t_vertex_energy, t_bjorkenx, t_inelasticity, t_qsqr, t_transverse_momentum;
  double t_vertex[3], t_momentum[3];
  bool t_iscc;

  // Reco
  int r_particles, r_tracks, r_showers;
  double r_vertex[3];

  // Variables associated with mcparticle_tree
  int p_pdgcode, p_id, p_n_hits, p_status;
  double p_vertex[3], p_end[3], p_momentum[3];
  double p_energy, p_mass, p_transverse_momentum, p_momentum_magnitude;

  // Variables associated with recotrack_tree
  unsigned int tr_dedx_size, tr_residual_range_size;
  int tr_n_hits, tr_id_energy, tr_id_charge, tr_id_hits;
  double tr_pida, tr_chi2_mu, tr_chi2_pi, tr_chi2_pr, tr_chi2_ka;
  double tr_range_mom_proton, tr_range_mom_muon, tr_mcs_mom_muon;
  double tr_vertex[3], tr_end[3], tr_dedx[100000], tr_residual_range[100000];
  double tr_length, tr_kinetic_energy, tr_range, tr_missing_energy;

  // Variables associated with recoshower_tree
  unsigned int sh_dedx_size;
  int sh_n_hits, sh_id_energy, sh_id_charge, sh_id_hits;
  double sh_start[3], sh_direction[3], sh_length, sh_open_angle, sh_energy, sh_dedx;

  // Variables for MCS momentum study
  double mcs_cluster_width;
  int mcs_pdgcode;
};

// Constructor
sbnd::AnalysisNTuple::AnalysisNTuple(Parameters const & config)
  :
  EDAnalyzer(config),
  m_detectorLengthX(config().inputs().detectorLengthX()),
  m_detectorLengthY(config().inputs().detectorLengthY()),
  m_detectorLengthZ(config().inputs().detectorLengthZ()),
  m_coordinateOffsetX(config().inputs().coordinateOffsetX()),
  m_coordinateOffsetY(config().inputs().coordinateOffsetY()),
  m_coordinateOffsetZ(config().inputs().coordinateOffsetZ()),
  m_selectedBorderX(config().inputs().selectedBorderX()),
  m_selectedBorderY(config().inputs().selectedBorderY()),
  m_selectedBorderZ(config().inputs().selectedBorderZ()),
  m_generator_label(config().inputs().generator_label()),
  m_geant_label(config().inputs().geant_label()),
  m_pandora_label(config().inputs().pandora_label()),
  m_reco_track_label(config().inputs().reco_track_label()),
  m_reco_shower_label(config().inputs().reco_shower_label()),
  m_reco_track_calorimetry_label(config().inputs().reco_track_calorimetry_label()),
  m_reco_track_particleid_label(config().inputs().reco_track_particleid_label()),
  m_hit_label(config().inputs().hit_label()),
  m_mcs_fitter(config().fitter)
{
}

void sbnd::AnalysisNTuple::analyze(art::Event const & e)
{

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  // Add one to the event counter
  event_id++;

  // Get current time stamp
  time_now = std::time(nullptr);
  
  // Proceed with the nTuple-filling if 
  //    True neutrino vertex contained
  //  & True and reconstructed neutrino vertex contained
  //  & At least 1 primary reconstructed track
  //  & All primary reconstructed tracks contained 
  
// ------------------------------------------------------------------------------
//                           GET THE HANDLES
// ------------------------------------------------------------------------------
  // Get the MCTruth handle
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel(m_generator_label, mct_handle);
  /*std::vector< art::Ptr<simb::MCTruth> > mct_vector;
  if (e.getByLabel(m_generator_label, mct_handle)) art::fill_ptr_vector(mct_vector, mct_handle);*/
  int mct_size = mct_handle->size();
 
  // Get the Track handle
  art::Handle< std::vector< recob::Track > > trk_handle;
  e.getByLabel(m_reco_track_label, trk_handle);
/*  std::vector< art::Ptr< recob::Track > > trk_vector;
  if (e.getByLabel(m_reco_track_label, trk_handle)) art::fill_ptr_vector(trk_vector, trk_handle);*/
  //int trk_size = trk_handle->size();
  
  // Get the Shower handle
  art::Handle< std::vector< recob::Shower > > shw_handle;
  e.getByLabel(m_reco_shower_label, shw_handle);
/*  std::vector< art::Ptr< recob::Shower > > shw_vector;
  if (e.getByLabel(m_reco_shower_label, shw_handle)) art::fill_ptr_vector(shw_vector, shw_handle);*/
  int shw_size = shw_handle->size();

  // Get the PFParticle handle
  art::Handle< std::vector< recob::PFParticle > > pfp_handle;
  e.getByLabel(m_pandora_label, pfp_handle );
  int pfp_size = pfp_handle->size();
  
  // Get the hit handle
  art::Handle< std::vector< recob::Hit > > hit_handle;
  e.getByLabel(m_hit_label, hit_handle );
  std::vector< art::Ptr< recob::Hit > > hit_vector;
  if (e.getByLabel(m_hit_label, hit_handle)) art::fill_ptr_vector(hit_vector, hit_handle);
  //int hit_size = hit_handle->size();

// ------------------------------------------------------------------------------
//                           EVENT-TREE INFORMATION
// ------------------------------------------------------------------------------
  // Loop over the truth handle and get everything we can
  if(mct_size != 1) return;
  
  if(mct_handle.isValid() &&
     shw_handle.isValid() && 
     pfp_handle.isValid() && pfp_size && 
     trk_handle.isValid() ) {

// ------------------------------------------------------------------------------
//                     EVENT-TREE RECONSTRUCTED INFORMATION
// ------------------------------------------------------------------------------
    
    // Initialise counter for the number of reconstructed primary particles
    unsigned int n_primaries          = 0;
    unsigned int n_primary_tracks     = 0;
    unsigned int n_primary_showers    = 0;
    unsigned int neutrino_id          = 0;
   
    bool neutrino_found = false;

    // Loop over PFParticles and find how many are primary and whether 
    // only 1 neutrino was found
    for(int i = 0; i < pfp_size; ++i) {
  
      art::Ptr< recob::PFParticle > pfp( pfp_handle, i );

      // Only look for a single muon neutrino event
      if(!neutrino_found && pfp->IsPrimary()){
        neutrino_found = true;
        neutrino_id = pfp->Self();
      }
      else if(neutrino_found && pfp->IsPrimary()) return;
    
    }

    if(!neutrino_found) return;

    // Get the GEANT information of the particles 
    art::FindManyP< simb::MCParticle > fmcp( mct_handle, e, m_geant_label );
    art::FindManyP< simb::GTruth > fmgt( mct_handle, e, m_generator_label );
    art::Ptr< simb::MCTruth > mct(mct_handle, 0);
    std::vector< art::Ptr<simb::MCParticle> > mcp_assn = fmcp.at(0);
    std::vector< art::Ptr<simb::GTruth> > mcgt_assn = fmgt.at(0);
    simb::MCNeutrino nu = mct->GetNeutrino();

    // Get vertex association
    art::FindManyP< recob::Vertex  > fvtx( pfp_handle, e, m_pandora_label );
    std::vector< art::Ptr<recob::Vertex> > vtx_assn = fvtx.at(neutrino_id);

    if(vtx_assn.size()  > 1) return;
    if(vtx_assn.size() == 0) return;

    // Set array to be current vertex position
    vtx_assn[0]->XYZ(r_vertex);
    float reco_vertex_x = r_vertex[0];
    float reco_vertex_y = r_vertex[1];
    float reco_vertex_z = r_vertex[2];
    
    // Make sure reconstructed vertex is within the TPC fiducial volume
    if (    (reco_vertex_x > (m_detectorLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
         || (reco_vertex_x < (-m_coordinateOffsetX + m_selectedBorderX)) 
         || (reco_vertex_y > (m_detectorLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
         || (reco_vertex_y < (-m_coordinateOffsetY + m_selectedBorderY)) 
         || (reco_vertex_z > (m_detectorLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
         || (reco_vertex_z < (-m_coordinateOffsetZ + m_selectedBorderZ))) return;
    
    // Get track associations with PFParticles from Pandora
    art::FindManyP< recob::Track   > fmtrk( pfp_handle, e, m_reco_track_label );
    art::FindManyP< recob::Cluster > fmclu( pfp_handle, e, m_pandora_label );

    // Find the number of reconstructed primary final state particles 
    for(int k = 0; k < pfp_size; ++k) {
  
      art::Ptr< recob::PFParticle > pfp( pfp_handle, k );

      pfparticle++;

      // If the PFParticle is no a primary, continue
      if(neutrino_id != pfp->Parent()) continue;
     
      // Counter for all jobs
      primary_pfparticle++;

      // For primary PFParticles get associated tracks and their calorimetry
      n_primaries++;
    
      std::vector< art::Ptr<recob::Track>   > trk_assn = fmtrk.at(pfp->Self());
      std::vector< art::Ptr<recob::Cluster> > clu_assn = fmclu.at(pfp->Self());

      if(trk_assn.size()) {
      
        art::FindManyP< anab::Calorimetry  > fmcal( trk_handle, e, m_reco_track_calorimetry_label );
        art::FindManyP< anab::ParticleID   > fmpid( trk_handle, e, m_reco_track_particleid_label );
        art::FindManyP< recob::Hit         > fmhit( trk_handle, e, m_reco_track_label );
        
        // Loop over tracks associated with primary PFParticles
        for(size_t i = 0; i < trk_assn.size(); ++i) {
        
          float track_vtx_x = trk_assn[i]->Vertex()[0];
          float track_vtx_y = trk_assn[i]->Vertex()[1];
          float track_vtx_z = trk_assn[i]->Vertex()[2];
          float track_end_x = trk_assn[i]->End()[0];
          float track_end_y = trk_assn[i]->End()[1];
          float track_end_z = trk_assn[i]->End()[2];

          // The border for contained tracks should be the edge of the active volume,
          // since this is where we can measure energy up to
          // Find out if one end of a track escapes (if so, MCS)
          bool does_vtx_escape =                                                                   
            (     (track_vtx_x > (m_detectorLengthX - m_coordinateOffsetX))                                     
               || (track_vtx_x < (-m_coordinateOffsetX))                                                    
               || (track_vtx_y > (m_detectorLengthY - m_coordinateOffsetY))                                     
               || (track_vtx_y < (-m_coordinateOffsetY))     
               || (track_vtx_z > (m_detectorLengthZ - m_coordinateOffsetZ))                                     
               || (track_vtx_z < (-m_coordinateOffsetZ)));                                                  
                                                                                                   
          bool does_end_escape =                                                                   
            (     (track_end_x > (m_detectorLengthX - m_coordinateOffsetX))                                     
               || (track_end_x < (-m_coordinateOffsetX))                                                    
               || (track_end_y > (m_detectorLengthY - m_coordinateOffsetY))                                     
               || (track_end_y < (-m_coordinateOffsetY))                                                    
               || (track_end_z > (m_detectorLengthZ - m_coordinateOffsetZ))                                     
               || (track_end_z < (-m_coordinateOffsetZ)));                                                  
                                                                                                   
          bool one_end_escapes = true;                                                             
          if(does_vtx_escape && does_end_escape)   one_end_escapes = false;                        
          if(!does_vtx_escape && !does_end_escape) one_end_escapes = false;      
          
          // Get the track-based variables
          std::vector< art::Ptr<anab::Calorimetry> > cal_assn = fmcal.at(trk_assn[i]->ID());
          std::vector< art::Ptr<anab::ParticleID> >  pid_assn = fmpid.at(trk_assn[i]->ID());
          std::vector< art::Ptr<recob::Hit> >        hit_assn = fmhit.at(trk_assn[i]->ID());
   
          // Loop over PID association
          for ( size_t j = 0; j < pid_assn.size(); ++j ){

            if (!pid_assn[j]) continue;
            if (!pid_assn[j]->PlaneID().isValid) continue;
              
            // Get the plane number
            int planenum = pid_assn[j]->PlaneID().Plane;

            // Only look at the collection plane, since this is where the dEdx
            // is acquired and we need this for the PIDA values
            if (planenum!=2) continue;
              
            // Loop over cal association
            for ( size_t k = 0; k < cal_assn.size(); ++k ){

              if (!cal_assn[k]) continue;
              if (!cal_assn[k]->PlaneID().isValid) continue;
                
              // Get the plane number
              int planenumcal = cal_assn[k]->PlaneID().Plane;

              // Only look at the collection plane, since this is where the dEdx
              // is acquired and we need this for the PIDA values
              if (planenumcal!=2) continue;

              // Add one to the counter for the event tree
              n_primary_tracks++;

              // Get associated MCParticle ID using 3 different methods:
              //    Which particle contributes the most energy to all the hits
              //    Which particle contributes the reco charge to all the hits
              //    Which particle is the biggest contributor to all the hits
              
              tr_id_energy      = RecoUtils::TrueParticleIDFromTotalTrueEnergy(hit_assn);
              tr_id_charge      = RecoUtils::TrueParticleIDFromTotalRecoCharge(hit_assn);
              tr_id_hits        = RecoUtils::TrueParticleIDFromTotalRecoHits(hit_assn);

              tr_n_hits         = hit_assn.size(); 
              tr_chi2_pr        = pid_assn[j]->Chi2Proton();
              tr_chi2_mu        = pid_assn[j]->Chi2Muon();
              tr_chi2_pi        = pid_assn[j]->Chi2Pion();
              tr_chi2_ka        = pid_assn[j]->Chi2Kaon();
              tr_pida           = pid_assn[j]->PIDA();
              tr_missing_energy = pid_assn[j]->MissingE();

              tr_kinetic_energy      = cal_assn[k]->KineticEnergy();
              tr_range               = cal_assn[k]->Range();
              tr_dedx_size           = cal_assn[k]->dEdx().size();
              tr_residual_range_size = cal_assn[k]->ResidualRange().size();
              for(unsigned int l = 0; l < tr_dedx_size; ++l) tr_dedx[l]                     = cal_assn[k]->dEdx()[l];
              for(unsigned int l = 0; l < tr_residual_range_size; ++l) tr_residual_range[l] = cal_assn[k]->ResidualRange()[l];

              tr_vertex[0] = track_vtx_x;
              tr_vertex[1] = track_vtx_y;
              tr_vertex[2] = track_vtx_z;
              
              tr_end[0] = track_end_x;
              tr_end[1] = track_end_y;
              tr_end[2] = track_end_z;

              tr_length = trk_assn[i]->Length();

              // Momentum
              //    Assign a momentum of 0 to every parameter and then 
              //    fill with the relevant value based on:
              //    Whether the track escapes
              //      - can have an MCS momentum
              //    If it is contained
              //      - proton range-based
              //      - muon range-based
              //
              // Further down the line, once PID has been performed properly, 
              // this should be strictly selected.
              tr_range_mom_proton = 0.;
              tr_range_mom_muon   = 0.;
              tr_mcs_mom_muon     = 0.;
              double reco_momentum_muon, reco_momentum_proton;
              if(one_end_escapes){
                recob::MCSFitResult mcs_result = m_mcs_fitter.fitMcs(*trk_assn[i]);
                reco_momentum_muon = mcs_result.bestMomentum();
                if(reco_momentum_muon < 0) tr_mcs_mom_muon = 0.;
                tr_mcs_mom_muon = reco_momentum_muon;
              }
              else if(!does_vtx_escape && !does_end_escape){
                reco_momentum_muon   = m_range_fitter.GetTrackMomentum(tr_length, 13); 
                reco_momentum_proton = m_range_fitter.GetTrackMomentum(tr_length, 2212); 
                if(reco_momentum_muon   == 0) tr_range_mom_muon   = 0.;
                if(reco_momentum_proton == 0) tr_range_mom_proton = 0.;
                tr_range_mom_muon   = reco_momentum_muon;
                tr_range_mom_proton = reco_momentum_proton;
              }
              
              // MCS Muon-Pion study
              //    Find a known muon
              //      Get the width of the cluster associated to the PFParticle
              //    Repeat for known charged pions
              //    Plot widths for both muons and pions
              if(clu_assn.size() == 0) continue;
              for(art::Ptr<simb::MCParticle> part : mcp_assn) {

                // Define a temporary histogram to fit
                // ATTENTION: Cut on PGD codes which refer to elements (Argon39 and above) 
                // Only interested in the final state PARTICLES
                std::cout << " MC ID : " << part->TrackId() << ", reco ID : " << tr_id_hits << std::endl;
                if(part->Process() != "primary" || part->PdgCode() >= 1000018039) continue;
                if(part->TrackId() != tr_id_hits) continue;
                if(part->PdgCode() == 13 || part->PdgCode() == 211 || part->PdgCode() == -211){
                  mcs_cluster_width = clu_assn[i]->Width();
                  mcs_pdgcode       = part->PdgCode();
                }
              }
              mcsstudy_tree->Fill();
              recotrack_tree->Fill();
            }
          }
        }
      }
    }
    if(shw_handle.isValid() && shw_size) {

      art::FindManyP< recob::Hit > fmhit( shw_handle, e, m_reco_shower_label );
      
      // Loop over PMA tracks and find any within 2 cm of the primary vertex
      // count them
      for(int i = 0; i < shw_size; ++i) {
      
        art::Ptr< recob::Shower > shw( shw_handle, i );
       
        // Distance s of the shower start from the reconstructed 
        // neutrino vertex
        double s_vtx = sqrt(pow(shw->ShowerStart()[0] - r_vertex[0],2) + pow(shw->ShowerStart()[1] - r_vertex[1],2) + pow(shw->ShowerStart()[2] - r_vertex[2],2));
       
        int bp = shw->best_plane();

        // Cut of 40 cm for showers to accommodate photon conversion after
        // neutral pion decay
        if(s_vtx < 40) {

          std::vector< art::Ptr<recob::Hit> > hit_assn = fmhit.at(i);
          
          // Add one to the counter for the event tree
          n_primary_showers++;

          // Get associated MCParticle ID using 3 different methods:
          //    Which particle contributes the most energy to all the hits
          //    Which particle contributes the reco charge to all the hits
          //    Which particle is the biggest contributor to all the hits
          sh_id_energy    = RecoUtils::TrueParticleIDFromTotalTrueEnergy(hit_assn);
          sh_id_charge    = RecoUtils::TrueParticleIDFromTotalRecoCharge(hit_assn);
          sh_id_hits      = RecoUtils::TrueParticleIDFromTotalRecoHits(hit_assn);
          
          sh_n_hits       = hit_assn.size();
          sh_dedx_size    = shw->dEdx().size();
          sh_start[0]     = shw->ShowerStart()[0];
          sh_start[1]     = shw->ShowerStart()[1];
          sh_start[2]     = shw->ShowerStart()[2];
          sh_direction[0] = shw->Direction()[0];
          sh_direction[1] = shw->Direction()[1];
          sh_direction[2] = shw->Direction()[2];
          sh_length       = shw->Length();
          sh_open_angle   = shw->OpenAngle();
          sh_energy       = shw->Energy()[bp];
          sh_dedx         = shw->dEdx()[bp];  

          recoshower_tree->Fill();
        }
      }
    }

    // Number of reconstructed primary particles, primary tracks and 
    // all reconstructed showers
    // Location of the primary reconstructed vertex
    r_tracks    = n_primary_tracks;
    r_showers   = n_primary_showers; 
    r_particles = n_primaries;
// ------------------------------------------------------------------------------
//                     EVENT-TREE TRUTH INFORMATION
// ------------------------------------------------------------------------------

    // Start defining truth variables
    t_nu_pdgcode          = nu.Nu().PdgCode();
    t_iscc                = nu.CCNC() == simb::curr_type_::kCC;
    t_interaction         = mcgt_assn[0]->fGint;
    t_scatter             = mcgt_assn[0]->fGscatter;
    t_vertex[0]           = nu.Nu().Vx();
    t_vertex[1]           = nu.Nu().Vy();
    t_vertex[2]           = nu.Nu().Vz();
    t_momentum[0]         = nu.Nu().Px();
    t_momentum[1]         = nu.Nu().Py();
    t_momentum[2]         = nu.Nu().Pz();
    t_vertex_energy       = nu.Nu().E();
    t_inv_mass            = nu.W();
    t_nu_lepton_angle     = nu.Theta();
    t_qsqr                = nu.QSqr();
    t_transverse_momentum = nu.Pt();
    t_bjorkenx            = nu.X();
    t_inelasticity        = nu.Y();

    t_particles       = 0;
    t_photons         = 0;
    t_electrons       = 0;
    t_neutral_pions   = 0;
    t_charged_pions   = 0;
    t_protons         = 0;
    t_neutrons        = 0;
    t_kaons           = 0;
    t_muons           = 0;

    // Counters and filling the mcparticle_tree
    for(art::Ptr<simb::MCParticle> part : mcp_assn) {

      //art::FindManyP< recob::Hit > fmhit( trk_handle, e, m_reco_track_label );
      
      // ATTENTION: Cut on PGD codes which refer to elements (Argon39 and above) 
      // Only interested in the final state PARTICLES
      if(part->Process() != "primary" || part->PdgCode() >= 1000018039) continue;
  
      std::vector< art::Ptr< recob::Hit > > hits = bt_serv->TrackIdToHits_Ps(part->TrackId(), hit_vector); 

      t_particles++;

      // Find the number of individual particle types
      if( part->PdgCode() == 22 )   t_photons++;
      if( part->PdgCode() == 111 )  t_neutral_pions++;
      if( part->PdgCode() == 2212 ) t_protons++;
      if( part->PdgCode() == 2112 ) t_neutrons++;
      if( part->PdgCode() == 11   || part->PdgCode() == -11 )  t_electrons++;
      if( part->PdgCode() == 13   || part->PdgCode() == -13 )  t_muons++;
      if( part->PdgCode() == 211  || part->PdgCode() == -211 ) t_charged_pions++;
      if( part->PdgCode() == 321  || part->PdgCode() == -321 || part->PdgCode() == 311 ) t_kaons++;

      p_id                  = part->TrackId();
      p_pdgcode             = part->PdgCode();
      p_status              = part->StatusCode();
      p_n_hits              = hits.size();
      p_vertex[0]           = part->Vx();
      p_vertex[1]           = part->Vy();
      p_vertex[2]           = part->Vz();
      p_end[0]              = part->EndX();
      p_end[1]              = part->EndY();
      p_end[2]              = part->EndZ();
      p_momentum[0]         = part->Px();
      p_momentum[1]         = part->Py();
      p_momentum[2]         = part->Pz();
      p_energy              = part->E();
      p_mass                = part->Mass();
      p_transverse_momentum = part->Pt();
      p_momentum_magnitude  = part->P();

      mcparticle_tree->Fill();

    }
    
    // Fill the event tree once everything has been set
    event_tree->Fill();

  }
}

void sbnd::AnalysisNTuple::beginJob()
{
  // Implementation of optional member function here.
  // Event id
  event_id = 0;

  // Initialise the counters
  // Particle based
  pfparticle         = 0;
  primary_pfparticle = 0;

  // Initiate trees
  event_tree      = new TTree("event_tree",      "Event tree: True and reconstructed SBND event information");
  mcparticle_tree = new TTree("particle_tree",   "MCParticle tree: True SBND initial and final state topology information");
  recotrack_tree  = new TTree("recotrack_tree",  "Track tree: Reconstructed final state track information");
  recoshower_tree = new TTree("recoshower_tree", "Shower tree: Reconstructed final state shower information");
  mcsstudy_tree   = new TTree("mcsstudy_tree",   "MCS tree: Cluster information for studying the MCS behaviour of muons and pions");

  // Event tree branches
  event_tree->Branch("event_id",                   &event_id,              "event_id/I");
  event_tree->Branch("time_now",                   &time_now,              "time_now/I");
  event_tree->Branch("t_nu_pdgcode",               &t_nu_pdgcode,          "t_nu_pdgcode/I");
  event_tree->Branch("t_iscc",                     &t_iscc,                "t_iscc/O");
  event_tree->Branch("t_interaction",              &t_interaction,         "t_interaction/I");
  event_tree->Branch("t_scatter",                  &t_scatter,             "t_scatter/I");
  event_tree->Branch("t_vertex",                   &t_vertex,              "t_vertex[3]/D");
  event_tree->Branch("t_momentum",                 &t_momentum,            "t_momentum[3]/D");
  event_tree->Branch("t_particles",                &t_particles,           "t_particles/I");
  event_tree->Branch("t_protons",                  &t_protons,             "t_protons/I");
  event_tree->Branch("t_neutrons",                 &t_neutrons,            "t_neutrons/I");
  event_tree->Branch("t_muons",                    &t_muons,               "t_muons/I");
  event_tree->Branch("t_charged_pions",            &t_charged_pions,       "t_charged_pions/I");
  event_tree->Branch("t_neutral_pions",            &t_neutral_pions,       "t_neutral_pions/I");
  event_tree->Branch("t_kaons",                    &t_kaons,               "t_kaons/I");
  event_tree->Branch("t_photons",                  &t_photons,             "t_photons/I");
  event_tree->Branch("t_electrons",                &t_electrons,           "t_electrons/I");
  event_tree->Branch("t_inv_mass",                 &t_inv_mass,            "t_inv_mass/D");
  event_tree->Branch("t_qsqr",                     &t_qsqr,                "t_qsqr/D");
  event_tree->Branch("t_transverse_momentum",      &t_transverse_momentum, "t_transverse_momentum/D");
  event_tree->Branch("t_bjorkenx",                 &t_bjorkenx,            "t_bjorkenx/D");
  event_tree->Branch("t_nu_lepton_angle",          &t_nu_lepton_angle,     "t_nu_lepton_angle/D");
  event_tree->Branch("t_vertex_energy",            &t_vertex_energy,       "t_vertex_energy/D");
  event_tree->Branch("r_particles",                &r_particles,           "r_particles/I");
  event_tree->Branch("r_tracks",                   &r_tracks,              "r_tracks/I");
  event_tree->Branch("r_showers",                  &r_showers,             "r_showers/I");
  event_tree->Branch("r_vertex",                   &r_vertex,              "r_vertex[3]/D");

  // MCParticle tree branches
  // Variables associated with mcparticle_tree
  mcparticle_tree->Branch("event_id",              &event_id,              "event_id/I");
  mcparticle_tree->Branch("time_now",              &time_now,              "time_now/I");
  mcparticle_tree->Branch("p_id",                  &p_id,                  "p_id/I");
  mcparticle_tree->Branch("p_pdgcode",             &p_pdgcode,             "p_pdgcode/I");
  mcparticle_tree->Branch("p_status",              &p_status,              "p_status/I");
  mcparticle_tree->Branch("p_n_hits",              &p_n_hits,              "p_n_hits/I");
  mcparticle_tree->Branch("p_vertex",              &p_vertex,              "p_vertex[3]/D");
  mcparticle_tree->Branch("p_end",                 &p_end,                 "p_end[3]/D");
  mcparticle_tree->Branch("p_momentum",            &p_momentum,            "p_momentum[3]/D");
  mcparticle_tree->Branch("p_energy",              &p_energy,              "p_energy/D");
  mcparticle_tree->Branch("p_mass",                &p_mass,                "p_mass/D");
  mcparticle_tree->Branch("p_transverse_momentum", &p_transverse_momentum, "p_transverse_momentum/D");
  mcparticle_tree->Branch("p_momentum_magnitude",  &p_momentum_magnitude,  "p_momentum_magnitude/D");
  
  // Reco Track tree branches
  recotrack_tree->Branch("event_id",               &event_id,               "event_id/I");
  recotrack_tree->Branch("time_now",               &time_now,               "time_now/I");
  recotrack_tree->Branch("tr_n_hits",              &tr_n_hits,              "tr_n_hits/I");
  recotrack_tree->Branch("tr_id_energy",           &tr_id_energy,           "tr_id_energy/I");
  recotrack_tree->Branch("tr_id_charge",           &tr_id_charge,           "tr_id_charge/I");
  recotrack_tree->Branch("tr_id_hits",             &tr_id_hits,             "tr_id_hits/I");
  recotrack_tree->Branch("tr_dedx_size",           &tr_dedx_size,           "tr_dedx_size/i");
  recotrack_tree->Branch("tr_residual_range_size", &tr_residual_range_size, "tr_residual_range_size/i");
  recotrack_tree->Branch("tr_pida",                &tr_pida,                "tr_pida/D");
  recotrack_tree->Branch("tr_chi2_mu",             &tr_chi2_mu,             "tr_chi2_mu/D");
  recotrack_tree->Branch("tr_chi2_pi",             &tr_chi2_pi,             "tr_chi2_pi/D");
  recotrack_tree->Branch("tr_chi2_pr",             &tr_chi2_pr,             "tr_chi2_pr/D");
  recotrack_tree->Branch("tr_chi2_ka",             &tr_chi2_ka,             "tr_chi2_ka/D");
  recotrack_tree->Branch("tr_range_mom_proton",    &tr_range_mom_proton,    "tr_range_mom_proton/D");
  recotrack_tree->Branch("tr_range_mom_muon",      &tr_range_mom_muon,      "tr_range_mom_muon/D");
  recotrack_tree->Branch("tr_mcs_mom_muon",        &tr_mcs_mom_muon,        "tr_mcs_mom_muon/D");
  recotrack_tree->Branch("tr_vertex",              &tr_vertex,              "tr_vertex[3]/D");
  recotrack_tree->Branch("tr_end",                 &tr_end,                 "tr_end[3]/D");
  recotrack_tree->Branch("tr_dedx",                &tr_dedx,                ("tr_dedx[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree->Branch("tr_residual_range",      &tr_residual_range,      ("tr_residual_range[" + std::to_string(100000)+"]/D").c_str());
  recotrack_tree->Branch("tr_length",              &tr_length,              "tr_length/D");
  recotrack_tree->Branch("tr_range",               &tr_range,               "tr_range/D");
  recotrack_tree->Branch("tr_kinetic_energy",      &tr_kinetic_energy,      "tr_kinetic_energy/D");
  recotrack_tree->Branch("tr_missing_energy",      &tr_missing_energy,      "tr_missing_energy/D");
  
  // Reco Shower tree branches
  recoshower_tree->Branch("event_id",              &event_id,               "event_id/I");
  recoshower_tree->Branch("time_now",              &time_now,               "time_now/I");
  recoshower_tree->Branch("sh_n_hits",             &sh_n_hits,              "sh_n_hits/I");
  recoshower_tree->Branch("sh_id_energy",          &sh_id_energy,           "sh_id_energy/I");
  recoshower_tree->Branch("sh_id_charge",          &sh_id_charge,           "sh_id_charge/I");
  recoshower_tree->Branch("sh_id_hits",            &sh_id_hits,             "sh_id_hits/I");
  recoshower_tree->Branch("sh_dedx_size",          &sh_dedx_size,           "sh_dedx_size/i");
  recoshower_tree->Branch("sh_start",              &sh_start,               "sh_start[3]/D");
  recoshower_tree->Branch("sh_direction",          &sh_direction,           "sh_direction[3]/D");
  recoshower_tree->Branch("sh_length",             &sh_length,              "sh_length/D");
  recoshower_tree->Branch("sh_open_angle",         &sh_open_angle,          "sh_open_angle/D");
  recoshower_tree->Branch("sh_energy",             &sh_energy,              "sh_energy/D");
  recoshower_tree->Branch("sh_dedx",               &sh_dedx,                "sh_dedx/D");

  // Muon, pion cluster info
  mcsstudy_tree->Branch("mcs_cluster_width", &mcs_cluster_width, "mcs_cluster_width/D");
  mcsstudy_tree->Branch("mcs_pdgcode",       &mcs_pdgcode,       "mcs_pdgcode/I");

  // Set directories
  event_tree->SetDirectory(0);
  mcparticle_tree->SetDirectory(0);
  recotrack_tree->SetDirectory(0);
  recoshower_tree->SetDirectory(0);
  mcsstudy_tree->SetDirectory(0);
}

void sbnd::AnalysisNTuple::endJob()
{
  // Print the tree, write the file, close
  // This relative path is needed for grid jobs
  TFile file("output_file.root", "RECREATE");
  //TFile file("/pnfs/sbnd/persistent/users/rsjones/test_analysis/output_file.root", "RECREATE");
  event_tree->Write();
  mcparticle_tree->Write();
  recotrack_tree->Write();
  recoshower_tree->Write();
  mcsstudy_tree->Write();
  file.Write();
  file.Close();

  delete event_tree;
  delete mcparticle_tree;
  delete recotrack_tree;
  delete recoshower_tree;
  delete mcsstudy_tree;

}
DEFINE_ART_MODULE(sbnd::AnalysisNTuple)
