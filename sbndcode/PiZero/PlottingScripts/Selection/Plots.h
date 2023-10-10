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