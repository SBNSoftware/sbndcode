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

std::vector<Cut> ncpizero_incl_broad_cuts = {
  { "presel", "!slc_is_clear_cosmic && slc_is_fv", "Pre-Selection", kOrange-2 },
  { "cosmic_rejection", "slc_crumbs_score>-0.025", "Cosmic Rejection", kGreen+1 },
  { "muon_rejection", "slc_n_primary_razzled_muons==0", "Muon Rejection", kRed-9 },
  { "photons_selection", "slc_n_pfps>1 && slc_n_primary_razzled_photons>1 && slc_best_pzc_good_kinematics", "Photon Selection", kBlue-9 },
};

std::vector<Cut> ncpizero_0p0pi_cuts = ncpizero_incl_cuts;
ncpizero_0p0pi_cuts.push_back({ "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" });
ncpizero_0p0pi_cuts.push_back({ "no_razzled_protons", "slc_n_primary_razzled_protons_thresh==0", "No Razzled Protons" });

std::vector<Cut> ncpizero_0p0pi_broad_cuts = ncpizero_incl_broad_cuts;
ncpizero_0p0pi_broad_cuts.push_back({ "pion_rejection", "slc_n_primary_razzled_pions_thresh==0", "Pion Rejection" });
ncpizero_0p0pi_broad_cuts.push_back({ "proton_rejection", "slc_n_primary_razzled_protons_thresh==0", "Proton Rejection" });

std::vector<Cut> ncpizero_1p0pi_cuts = ncpizero_incl_cuts;
ncpizero_1p0pi_cuts.push_back({ "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" });
ncpizero_1p0pi_cuts.push_back({ "one_razzled_proton", "slc_n_primary_razzled_protons_thresh==1", "One Razzled Proton" });

std::vector<Cut> ncpizero_1p0pi_broad_cuts = ncpizero_incl_broad_cuts;
ncpizero_1p0pi_broad_cuts.push_back({ "pion_rejection", "slc_n_primary_razzled_pions_thresh==0", "Pion Rejection" });
ncpizero_1p0pi_broad_cuts.push_back({ "proton_selection", "slc_n_primary_razzled_protons_thresh==1", "Proton Selection" });

std::vector<Cut> ncpizero_Np0pi_cuts = ncpizero_incl_cuts;
ncpizero_Np0pi_cuts.push_back({ "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" });
ncpizero_Np0pi_cuts.push_back({ "at_least_one_razzled_proton", "slc_n_primary_razzled_protons_thresh>0", "At Least One Razzled Proton" });

std::vector<Cut> ncpizero_Np0pi_broad_cuts = ncpizero_incl_broad_cuts;
ncpizero_Np0pi_broad_cuts.push_back({ "pion_rejection", "slc_n_primary_razzled_pions_thresh==0", "Pion Rejection" });
ncpizero_Np0pi_broad_cuts.push_back({ "proton_selection", "slc_n_primary_razzled_protons_thresh>0", "Proton Selection" });

std::vector<Cut> ncpizero_Xp0pi_cuts = ncpizero_incl_cuts;
ncpizero_Xp0pi_cuts.push_back({ "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" });

std::vector<Cut> ncpizero_Xp0pi_broad_cuts = ncpizero_incl_broad_cuts;
ncpizero_Xp0pi_broad_cuts.push_back({ "pion_rejection", "slc_n_primary_razzled_pions_thresh==0", "Pion Rejection" });
