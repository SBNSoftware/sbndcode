#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

std::vector<Cut> ncpizero_incl_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
};

std::vector<Cut> ncpizero_0p0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
  { "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" },
  { "no_razzled_protons", "slc_n_primary_razzled_protons_thresh==0", "No Razzled Protons" },
};

std::vector<Cut> ncpizero_1p0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
  { "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" },
  { "one_razzled_proton", "slc_n_primary_razzled_protons_thresh==1", "One Razzled Proton" },
};

std::vector<Cut> ncpizero_Np0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
  { "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" },
  { "has_one_razzled_proton", "slc_n_primary_razzled_protons_thresh>0", "Has One Razzled Proton" },
};

std::vector<Cut> ncpizero_Xp0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
  { "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" },
};
