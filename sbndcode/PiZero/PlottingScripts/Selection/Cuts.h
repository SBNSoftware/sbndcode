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
