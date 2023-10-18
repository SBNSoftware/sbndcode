#include "/sbnd/app/users/hlay/plotting_utils/Structs.h"

std::vector<Plot> selection_plots = {
  { "slc_is_clear_cosmic", "slc_is_clear_cosmic", ";Is Clear Cosmic?;Slices",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_is_fv", "slc_is_fv", ";IsFV?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_crumbs_score", "slc_crumbs_score", ";CRUMBS Score;Slices",
    50, -1.5, 1. },
  { "slc_crumbs_nc_score", "slc_crumbs_nc_score", ";CRUMBS NC Score;Slices",
    50, -1.5, 1. },
  { "slc_crumbs_ccnue_score", "slc_crumbs_ccnue_score", ";CRUMBS CC#nu_{e} Score;Slices",
    50, -1.5, 1. },
  { "slc_n_dazzle_muons", "slc_n_dazzle_muons", ";N Dazzle Muons;Slices",
    5, -0.5, 4.5 },
  { "slc_n_shws", "slc_n_shws", ";N Showers;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzle_photons", "slc_n_razzle_photons", ";N Razzle Photons;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzled_photons", "slc_n_razzled_photons", ";N Razzled Photons;Slices",
    5, -0.5, 4.5 },
};

std::vector<Plot> old_plots = {
  { "slc_is_clear_cosmic", "slc_is_clear_cosmic", ";Is Clear Cosmic?;Slices",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_is_fv", "slc_is_fv", ";IsFV?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_crumbs_score", "slc_crumbs_score", ";CRUMBS Score;Slices",
    50, -1.5, 1. },
  { "slc_crumbs_nc_score", "slc_crumbs_nc_score", ";CRUMBS NC Score;Slices",
    50, -1.5, 1. },
  { "slc_crumbs_ccnue_score", "slc_crumbs_ccnue_score", ";CRUMBS CC#nu_{e} Score;Slices",
    50, -1.5, 1. },
  { "slc_n_dazzle_muons", "slc_n_dazzle_muons", ";N Dazzle Muons;Slices",
    5, -0.5, 4.5 },
  { "slc_n_shws", "slc_n_shws", ";N Showers;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzle_photons", "slc_n_razzle_photons", ";N Razzle Photons;Slices",
    5, -0.5, 4.5 },
};

std::vector<Plot> true_observables = {
  { "nu_pz_pizero_mom", "nu_pz_pizero_mom", ";p_{#pi^{0}} (GeV/c);Events",
    25, 0, 2 },
  { "nu_pz_cos_theta_pizero", "nu_pz_cos_theta_pizero", ";cos(#theta_{#pi^{0}});Events",
    25, -1, 1 },
  { "nu_pz_cos_com", "nu_pz_cos_com", ";cos(#theta_{CoM});Events",
    25, 0, 1 },
  { "nu_pz_decay_asymmetry", "nu_pz_decay_asymmetry", ";Decay Asymmetry;Events",
    25, 0, 1 },
};

std::vector<Plot> observables = {
  { "slc_pzc_invariant_mass", "slc_pzc_invariant_mass[0]", ";M_{#gamma#gamma} (MeV/c^{2});Events",
    25, 0, 500 },
  { "slc_pzc_pizero_mom", "slc_pzc_pizero_mom[0]", ";p_{#pi^{0}} (GeV/c);Events",
    25, 0, 2 },
  { "slc_pzc_cos_theta_pizero", "slc_pzc_cos_theta_pizero[0]", ";cos(#theta_{#pi^{0}});Events",
    25, -1, 1 },
  { "slc_pzc_cos_com", "slc_pzc_cos_com[0]", ";cos(#theta_{CoM});Events",
    25, 0, 1 },
  { "slc_pzc_decay_asymmetry", "slc_pzc_decay_asymmetry[0]", ";Decay Asymmetry;Events",
    25, 0, 1 },
};

std::vector<Plot> no_plots = {};
