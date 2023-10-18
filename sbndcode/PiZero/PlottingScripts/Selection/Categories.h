#include "/sbnd/app/users/hlay/plotting_utils/Structs.h"

const TCut nu_base = "slc_true_event_type!=6 && slc_true_event_type!=7";
const TCut nc_base = nu_base + "slc_true_fv && slc_true_ccnc==1";
const TCut ccnumu_base = nu_base + "slc_true_fv && slc_true_ccnc==0 && abs(slc_true_pdg)==14";
const TCut ccnue_base = nu_base + "slc_true_fv && slc_true_ccnc==0 && abs(slc_true_pdg)==12";
const TCut dirt_base = nu_base + "!slc_true_av";
const TCut nonfv_base = nu_base + "slc_true_av && !slc_true_fv";
const TCut cosmic_base = !nu_base;

const TCut true_nu_base = "nu_event_type!=6 && nu_event_type!=7";
const TCut true_nc_base = true_nu_base + "nu_fv && nu_ccnc==1";
const TCut true_ccnumu_base = true_nu_base + "nu_fv && nu_ccnc==0 && abs(nu_pdg)==14";

const TCut ncpizero_cut = nc_base + "slc_true_n_neutral_pions>=1";
const TCut true_ncpizero_cut = true_nc_base + "nu_n_neutral_pions>=1";
const TCut ncpizero_0p0pi_cut = ncpizero_cut + "slc_true_n_neutral_pions==1 && slc_true_n_protons==0 && slc_true_n_charged_pions==0";
const TCut true_ncpizero_0p0pi_cut = true_ncpizero_cut + "nu_n_neutral_pions==1 && nu_n_protons==0 && nu_n_charged_pions==0";
const TCut ncpizero_1p0pi_cut = ncpizero_cut + "slc_true_n_neutral_pions==1 && slc_true_n_protons==1 && slc_true_n_charged_pions==0";
const TCut true_ncpizero_1p0pi_cut = true_ncpizero_cut + "nu_n_neutral_pions==1 && nu_n_protons==1 && nu_n_charged_pions==0";
const TCut ncpizero_Np0pi_cut = ncpizero_cut + "slc_true_n_neutral_pions==1 && slc_true_n_protons>0 && slc_true_n_charged_pions==0";
const TCut true_ncpizero_Np0pi_cut = true_ncpizero_cut + "nu_n_neutral_pions==1 && nu_n_protons>0 && nu_n_charged_pions==0";
const TCut ncpizero_0pXpi_cut = ncpizero_cut + "slc_true_n_neutral_pions==1 && slc_true_n_protons==0";
const TCut true_ncpizero_0pXpi_cut = true_ncpizero_cut + "nu_n_neutral_pions==1 && nu_n_protons==0";
const TCut ncpizero_1pXpi_cut = ncpizero_cut + "slc_true_n_neutral_pions==1 && slc_true_n_protons==1";
const TCut true_ncpizero_1pXpi_cut = true_ncpizero_cut + "nu_n_neutral_pions==1 && nu_n_protons==1";
const TCut ncpizero_NpXpi_cut = ncpizero_cut + "slc_true_n_neutral_pions==1 && slc_true_n_protons>0";
const TCut true_ncpizero_NpXpi_cut = true_ncpizero_cut + "nu_n_neutral_pions==1 && nu_n_protons>0";

const TCut ccpizero_cut = ccnumu_base + "slc_true_n_neutral_pions>=1";
const TCut true_ccpizero_cut = true_ccnumu_base + "nu_n_neutral_pions>=1";

const std::vector<Cut> ncpizero_categories = {
  { "Signal", ncpizero_cut + "slc_comp>.5", "Signal (NC #pi^{0})", kMagenta+2 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_0p0pi_categories = {
  { "Signal", ncpizero_0p0pi_cut + "slc_comp>.5", "Signal (NC 1#pi^{0}0p0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", !ncpizero_0p0pi_cut + ncpizero_cut, "Other NC #pi^{0}", kMagenta-9 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_0p0pi_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_1p0pi_categories = {
  { "Signal", ncpizero_1p0pi_cut + "slc_comp>.5", "Signal (NC 1#pi^{0}1p0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", !ncpizero_1p0pi_cut + ncpizero_cut, "Other NC #pi^{0}", kMagenta-9 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_1p0pi_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_Np0pi_categories = {
  { "Signal", ncpizero_Np0pi_cut + "slc_comp>.5", "Signal (NC 1#pi^{0}Np0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", !ncpizero_Np0pi_cut + ncpizero_cut, "Other NC #pi^{0}", kMagenta-9 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_Np0pi_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_0pXpi_categories = {
  { "Signal", ncpizero_0pXpi_cut + "slc_comp>.5", "Signal (NC 1#pi^{0}0p)", kMagenta+2 },
  { "NCPiZero", !ncpizero_0pXpi_cut + ncpizero_cut, "Other NC #pi^{0}", kMagenta-9 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_0pXpi_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_1pXpi_categories = {
  { "Signal", ncpizero_1pXpi_cut + "slc_comp>.5", "Signal (NC 1#pi^{0}1p)", kMagenta+2 },
  { "NCPiZero", !ncpizero_1pXpi_cut + ncpizero_cut, "Other NC #pi^{0}", kMagenta-9 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_1pXpi_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_NpXpi_categories = {
  { "Signal", ncpizero_NpXpi_cut + "slc_comp>.5", "Signal (NC 1#pi^{0}Np)", kMagenta+2 },
  { "NCPiZero", !ncpizero_NpXpi_cut + ncpizero_cut, "Other NC #pi^{0}", kMagenta-9 },
  { "NC", nc_base + !ncpizero_cut, "Other NC", kOrange+2 },
  { "CCNuMu", ccnumu_base, "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ncpizero_NpXpi_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ccpizero_categories = {
  { "Signal", ccpizero_cut + "slc_comp>.5", "Signal (CC #nu_{#mu} #pi^{0})", kMagenta+2 },
  { "CCNuMu", ccnumu_base + !ccpizero_cut, "Other CC #nu_{#mu}", kGreen+2 },
  { "NC", nc_base, "NC", kOrange+2 },
  { "CCNuE", ccnue_base, "CC #nu_{e}", kCyan+2 },
  { "Dirt", dirt_base, "Dirt", kOrange+3 },
  { "NonFVNu", nonfv_base, "Non-FV #nu", kGray+2 },
  { "Cosmic", cosmic_base, "Cosmic", kRed+1 },
  { "BadRecoSignal", ccpizero_cut + "slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> selection_categories = ncpizero_categories;

const std::vector<Cut> event_modes = {
  { "Coherent", "nu_mode==3", "Coherent", kRed },
  { "QE", "nu_mode==0", "QE", kBlue },
  { "MEC", "nu_mode==10", "MEC", kMagenta+2 },
  { "Resonant", "nu_mode==1", "Resonant", kGreen+1 },
  { "DIS", "nu_mode==2", "DIS", kOrange+2 }
};
