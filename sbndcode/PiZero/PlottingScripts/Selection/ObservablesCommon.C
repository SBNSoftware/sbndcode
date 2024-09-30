#pragma once

#include "Common.C"

constexpr double kPiZeroMass = 134.9769;

std::vector<int> *slc_true_event_type_incl = 0;
std::vector<float> *slc_comp = 0;
std::vector<bool> *slc_sel_incl = 0;
std::vector<double> *slc_best_pzc_pizero_mom = 0, *slc_best_pzc_invariant_mass = 0, *slc_best_pzc_cos_theta_pizero = 0;

std::vector<size_t> *slc_best_pzc_photon_0_id = 0, *slc_best_pzc_photon_1_id = 0;
std::vector<int> *slc_best_pzc_photon_0_true_trackid = 0, *slc_best_pzc_photon_1_true_trackid = 0;
std::vector<float> *slc_best_pzc_photon_0_comp = 0, *slc_best_pzc_photon_1_comp = 0, *slc_best_pzc_photon_0_pur = 0, *slc_best_pzc_photon_1_pur = 0;

std::vector<std::vector<double>> *slc_true_pz_pizero_mom = 0, *slc_true_pz_cos_theta_pizero = 0, *slc_pfp_shower_dir_x = 0,
  *slc_pfp_shower_dir_y = 0, *slc_pfp_shower_dir_z = 0, *slc_pfp_track_dir_x = 0, *slc_pfp_track_dir_y = 0, *slc_pfp_track_dir_z = 0,
  *slc_pfp_shower_energy = 0, *slc_pfp_true_energy = 0, *slc_true_pz_open_angle = 0, *slc_true_pz_gamma0_energy = 0, *slc_true_pz_gamma1_energy = 0,
  *slc_pfp_true_p_x = 0, *slc_pfp_true_p_y = 0, *slc_pfp_true_p_z = 0;
std::vector<std::vector<int>> *slc_true_pz_gamma0_trackid = 0, *slc_true_pz_gamma1_trackid = 0;

TFile* file = TFile::Open("/exp/sbnd/app/users/hlay/ncpizero/srcs/sbndcode/sbndcode/PiZero/ShowerEnergyCorrection/shower_energy_correction_hist_NCPiZeroBv2.root");
TProfile *fShowerEnergyCorrectionHist = (TProfile*) file->Get("hShowerEnergy2DRecoFractionalResolution_pfx");

void InitialiseTree(TChain *tree);

double CorrectEnergy(const double &energy);

void InitialiseTree(TChain *tree)
{
  tree->SetBranchStatus("*", 0);

  tree->SetBranchAddress("slc_true_event_type_incl", &slc_true_event_type_incl);
  tree->SetBranchAddress("slc_comp", &slc_comp);
  tree->SetBranchAddress("slc_sel_incl", &slc_sel_incl);
  tree->SetBranchAddress("slc_best_pzc_pizero_mom", &slc_best_pzc_pizero_mom);
  tree->SetBranchAddress("slc_best_pzc_invariant_mass", &slc_best_pzc_invariant_mass);
  tree->SetBranchAddress("slc_best_pzc_cos_theta_pizero", &slc_best_pzc_cos_theta_pizero);
  tree->SetBranchAddress("slc_best_pzc_photon_0_id", &slc_best_pzc_photon_0_id);
  tree->SetBranchAddress("slc_best_pzc_photon_1_id", &slc_best_pzc_photon_1_id);
  tree->SetBranchAddress("slc_best_pzc_photon_0_true_trackid", &slc_best_pzc_photon_0_true_trackid);
  tree->SetBranchAddress("slc_best_pzc_photon_1_true_trackid", &slc_best_pzc_photon_1_true_trackid);
  tree->SetBranchAddress("slc_best_pzc_photon_0_comp", &slc_best_pzc_photon_0_comp);
  tree->SetBranchAddress("slc_best_pzc_photon_1_comp", &slc_best_pzc_photon_1_comp);
  tree->SetBranchAddress("slc_best_pzc_photon_0_pur", &slc_best_pzc_photon_0_pur);
  tree->SetBranchAddress("slc_best_pzc_photon_1_pur", &slc_best_pzc_photon_1_pur);
  tree->SetBranchAddress("slc_true_pz_pizero_mom", &slc_true_pz_pizero_mom);
  tree->SetBranchAddress("slc_true_pz_cos_theta_pizero", &slc_true_pz_cos_theta_pizero);
  tree->SetBranchAddress("slc_pfp_shower_dir_x", &slc_pfp_shower_dir_x);
  tree->SetBranchAddress("slc_pfp_shower_dir_y", &slc_pfp_shower_dir_y);
  tree->SetBranchAddress("slc_pfp_shower_dir_z", &slc_pfp_shower_dir_z);
  tree->SetBranchAddress("slc_pfp_track_dir_x", &slc_pfp_track_dir_x);
  tree->SetBranchAddress("slc_pfp_track_dir_y", &slc_pfp_track_dir_y);
  tree->SetBranchAddress("slc_pfp_track_dir_z", &slc_pfp_track_dir_z);
  tree->SetBranchAddress("slc_pfp_true_p_x", &slc_pfp_true_p_x);
  tree->SetBranchAddress("slc_pfp_true_p_y", &slc_pfp_true_p_y);
  tree->SetBranchAddress("slc_pfp_true_p_z", &slc_pfp_true_p_z);
  tree->SetBranchAddress("slc_pfp_shower_energy", &slc_pfp_shower_energy);
  tree->SetBranchAddress("slc_pfp_true_energy", &slc_pfp_true_energy);
  tree->SetBranchAddress("slc_true_pz_open_angle", &slc_true_pz_open_angle);
  tree->SetBranchAddress("slc_true_pz_gamma0_energy", &slc_true_pz_gamma0_energy);
  tree->SetBranchAddress("slc_true_pz_gamma1_energy", &slc_true_pz_gamma1_energy);
  tree->SetBranchAddress("slc_true_pz_gamma0_trackid", &slc_true_pz_gamma0_trackid);
  tree->SetBranchAddress("slc_true_pz_gamma1_trackid", &slc_true_pz_gamma1_trackid);
}

double CorrectEnergy(const double &energy)
{
  const int bin = fShowerEnergyCorrectionHist->FindBin(energy);

  return energy * (1 - fShowerEnergyCorrectionHist->GetBinContent(bin));
}
