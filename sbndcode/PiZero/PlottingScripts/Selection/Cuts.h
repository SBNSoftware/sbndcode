#include "/sbnd/app/users/hlay/plotting_utils/Structs.h"

std::vector<Cut> razzle_dazzle_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_dazzle_muons", "slc_n_dazzle_muons==0", "No Dazzle Muons" },
  { "has_two_shws", "slc_n_shws>1", "Has Two Showers" },
  { "has_two_razzle_photons", "slc_n_razzle_photons>1", "Has Two Razzle Photons" },
};

std::vector<Cut> dazzle_muons_razzled_photons_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_dazzle_muons", "slc_n_dazzle_muons==0", "No Dazzle Muons" },
  { "has_two_pfps", "slc_n_pfps>2", "Has Two PFPs" }, // note the 'neutrino' counts as a PFP so really we're requiring 3
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> dazzle_muons_pions_razzled_photons_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_dazzle_muons", "slc_n_dazzle_muons==0", "No Dazzle Muons" },
  { "no_dazzle_pions", "slc_n_dazzle_pions==0", "No Dazzle Pions" },
  { "has_two_pfps", "slc_n_pfps>2", "Has Two PFPs" }, // note the 'neutrino' counts as a PFP so really we're requiring 3
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> razzled_muons_razzle_photons_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "has_two_shws", "slc_n_shws>1", "Has Two Showers" },
  { "has_two_razzle_photons", "slc_n_razzle_photons>1", "Has Two Razzle Photons" },
};

std::vector<Cut> razzled_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "has_two_pfps", "slc_n_pfps>2", "Has Two PFPs" }, // note the 'neutrino' counts as a PFP so really we're requiring 3
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> razzled_cuts_no_pion = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "no_razzled_pions", "slc_n_razzled_pions==0", "No Razzled Pions" },
  { "has_two_pfps", "slc_n_pfps>2", "Has Two PFPs" }, // note the 'neutrino' counts as a PFP so really we're requiring 3
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ncpizero_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
};

std::vector<Cut> ncpizero_cuts_dazzle_muons = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_dazzle_muons", "slc_n_dazzle_muons==0", "No Dazzle Muons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
  { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" },
};

std::vector<Cut> ncpizero_0p0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.125", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "no_razzled_pions", "slc_n_razzled_pions==0", "No Razzled Pions" },
  { "no_razzled_protons", "slc_n_razzled_protons==0", "No Razzled Protons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ncpizero_1p0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>0", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "no_razzled_pions", "slc_n_razzled_pions==0", "No Razzled Pions" },
  { "one_razzled_proton", "slc_n_razzled_protons==1", "One Razzled Proton" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ncpizero_Np0pi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>0", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "no_razzled_pions", "slc_n_razzled_pions==0", "No Razzled Pions" },
  { "has_one_razzled_proton", "slc_n_razzled_protons>0", "Has One Razzled Proton" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ncpizero_0pXpi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.125", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "no_razzled_protons", "slc_n_razzled_protons==0", "No Razzled Protons" },
  { "has_two_pfps", "slc_n_pfps>1", "Has Two PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ncpizero_1pXpi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>0", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "one_razzled_proton", "slc_n_razzled_protons==1", "One Razzled Proton" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ncpizero_NpXpi_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>0", "CRUMBS Cut" },
  { "no_razzled_muons", "slc_n_razzled_muons==0", "No Razzled Muons" },
  { "has_one_razzled_proton", "slc_n_razzled_protons>0", "Has One Razzled Proton" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};

std::vector<Cut> ccpizero_cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>0.325", "CRUMBS Cut" },
  { "has_one_razzled_muon", "slc_n_razzled_muons>0", "Has One Razzled Muon" },
  { "has_three_pfps", "slc_n_pfps>2", "Has Three PFPs" },
  { "has_two_razzled_photons", "slc_n_razzled_photons>1", "Has Two Razzled Photons" },
};
