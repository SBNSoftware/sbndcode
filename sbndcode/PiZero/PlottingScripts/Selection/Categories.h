#pragma once

#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

const TCut true_ncpizero_incl_cut = "nu_event_type_incl==0";
const TCut true_ncpizero_0p0pi_cut = "nu_event_type_0p0pi==0";
const TCut true_ncpizero_1p0pi_cut = "nu_event_type_1p0pi==0";
const TCut true_ncpizero_Np0pi_cut = "nu_event_type_Np0pi==0";
const TCut true_ncpizero_Xp0pi_cut = "nu_event_type_Xp0pi==0";

const std::vector<Cut> true_categories = {
  { "ncpizero_incl", true_ncpizero_incl_cut, "", kBlack },
  { "ncpizero_0p0pi", true_ncpizero_0p0pi_cut, "", kBlack },
  { "ncpizero_1p0pi", true_ncpizero_1p0pi_cut, "", kBlack },
  { "ncpizero_Xp0pi", true_ncpizero_Xp0pi_cut, "", kBlack },
  { "ncpizero_Np0pi", true_ncpizero_Np0pi_cut, "", kBlack },
};

const std::vector<Cut> ncpizero_incl_categories = {
  { "Signal", "slc_true_event_type_incl==0 && slc_comp>.5", "Signal (NC 1#pi^{0})", kMagenta+2 },
  { "NC", "slc_true_event_type_incl==2", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type_incl==3", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type_incl==4", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type_incl==5", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type_incl==6", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type_incl==7 || slc_true_event_type_incl==8", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type_incl==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_0p0pi_categories = {
  { "Signal", "slc_true_event_type_0p0pi==0 && slc_comp>.5", "Signal (NC 1#pi^{0}0p0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", "slc_true_event_type_0p0pi==1", "Other NC 1#pi^{0}", kMagenta-9 },
  { "NC", "slc_true_event_type_0p0pi==2", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type_0p0pi==3", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type_0p0pi==4", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type_0p0pi==5", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type_0p0pi==6", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type_0p0pi==7 || slc_true_event_type_0p0pi==8", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type_0p0pi==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_1p0pi_categories = {
  { "Signal", "slc_true_event_type_1p0pi==0 && slc_comp>.5", "Signal (NC 1#pi^{0}1p0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", "slc_true_event_type_1p0pi==1", "Other NC 1#pi^{0}", kMagenta-9 },
  { "NC", "slc_true_event_type_1p0pi==2", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type_1p0pi==3", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type_1p0pi==4", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type_1p0pi==5", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type_1p0pi==6", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type_1p0pi==7 || slc_true_event_type_1p0pi==8", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type_1p0pi==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_Np0pi_categories = {
  { "Signal", "slc_true_event_type_Np0pi==0 && slc_comp>.5", "Signal (NC 1#pi^{0}Np0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", "slc_true_event_type_Np0pi==1", "Other NC 1#pi^{0}", kMagenta-9 },
  { "NC", "slc_true_event_type_Np0pi==2", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type_Np0pi==3", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type_Np0pi==4", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type_Np0pi==5", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type_Np0pi==6", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type_Np0pi==7 || slc_true_event_type_Np0pi==8", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type_Np0pi==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> ncpizero_Xp0pi_categories = {
  { "Signal", "slc_true_event_type_Xp0pi==0 && slc_comp>.5", "Signal (NC 1#pi^{0}Xp0#pi^{#pm})", kMagenta+2 },
  { "NCPiZero", "slc_true_event_type_Xp0pi==1", "Other NC 1#pi^{0}", kMagenta-9 },
  { "NC", "slc_true_event_type_Xp0pi==2", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type_Xp0pi==3", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type_Xp0pi==4", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type_Xp0pi==5", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type_Xp0pi==6", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type_Xp0pi==7 || slc_true_event_type_Xp0pi==8", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type_Xp0pi==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> event_modes = {
  { "Coherent", "nu_mode==3", "Coherent", kRed },
  { "QE", "nu_mode==0", "QE", kBlue },
  { "MEC", "nu_mode==10", "MEC", kMagenta+2 },
  { "Resonant", "nu_mode==1", "Resonant", kGreen+1 },
  { "DIS", "nu_mode==2", "DIS", kOrange+2 }
};
