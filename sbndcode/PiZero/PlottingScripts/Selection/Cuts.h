#pragma once

#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

const Cut noCut                    = { "no_cut", "", "No Cut" };
const Cut notClearCosmic           = { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" };
const Cut fv                       = { "fv", "slc_is_fv", "FV" };
const Cut crumbs                   = { "crumbs", "slc_crumbs_score>-0.005", "CRUMBS Cut" };
const Cut noRazzledMuons           = { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" };
const Cut atLeastTwoPFPs           = { "at_least_two_pfps", "slc_n_pfps>1", "At Least Two PFPs" };
const Cut atLeastTwoRazzledPhotons = { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" };
const Cut goodPiZeroKinematics     = { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" };
const Cut noRazzledPions           = { "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" };
const Cut noRazzledProtons         = { "no_razzled_protons", "slc_n_primary_razzled_protons_thresh==0", "No Razzled Protons" };
const Cut oneRazzledProton         = { "one_razzled_proton", "slc_n_primary_razzled_protons_thresh==1", "One Razzled Proton" };
const Cut atLeastOneRazzledProton  = { "at_least_one_razzled_proton", "slc_n_primary_razzled_protons_thresh>0", "At Least One Razzled Proton" };
const Cut goodOpT0FracHigh         = { "good_opt0_frac_high", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE<0.756", "Good OpT0 Frac High" };
const Cut goodOpT0FracLow          = { "good_opt0_frac_low", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE>-0.700", "Good OpT0 Frac Low" };
const Cut opT0Score                = { "opt0_score", "slc_opt0_score>30", "OpT0 Score" };

const Cut preSel              = { "presel", "!slc_is_clear_cosmic && slc_is_fv", "Pre-Sel", kOrange-2 };
const Cut cosmicRej           = { "cosmic_rejection", "slc_crumbs_score>-0.005", "Cosmic Rej", kGreen+1 };
const Cut muonRej             = { "muon_rejection", "slc_n_primary_razzled_muons ==0", "Muon Rej", kRed-9 };
const Cut photonSel           = { "photons_selection", "slc_n_pfps>1 && slc_n_primary_razzled_photons>1 && slc_best_pzc_good_kinematics", "Photon Sel", kBlue-9 };
const Cut pionRej             = { "pion_rejection", "slc_n_primary_razzled_pions_thresh==0", "Pion Rej", kViolet - 5 };
const Cut protonRej           = { "proton_rejection", "slc_n_primary_razzled_protons_thresh==0", "Proton Rej", kViolet - 6 };
const Cut oneProtonSel        = { "proton_selection", "slc_n_primary_razzled_protons_thresh==1", "Proton Sel", kViolet - 6 };
const Cut atLeastOneProtonSel = { "proton_selection", "slc_n_primary_razzled_protons_thresh>0", "Proton Sel", kViolet - 6 };
const Cut opT0Cuts            = { "opt0_cuts",
                                  "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE<0.756 && (slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE>-0.700 && slc_opt0_score>30",
                                  "OpT0 Cuts" };

const std::vector<Cut> ncpizero_incl_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbs,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  goodOpT0FracHigh,
  goodOpT0FracLow,
  opT0Score
};

const std::vector<Cut> ncpizero_incl_broad_cuts = {
  noCut,
  preSel,
  cosmicRej,
  muonRej,
  photonSel,
  opT0Cuts
};

const std::vector<Cut> ncpizero_0p0pi_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbs,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  noRazzledPions,
  noRazzledProtons,
  goodOpT0FracHigh,
  goodOpT0FracLow,
  opT0Score
};

const std::vector<Cut> ncpizero_0p0pi_broad_cuts = {
  noCut,
  preSel,
  cosmicRej,
  muonRej,
  photonSel,
  pionRej,
  protonRej,
  opT0Cuts,
};

const std::vector<Cut> ncpizero_1p0pi_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbs,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  noRazzledPions,
  oneRazzledProton,
  goodOpT0FracHigh,
  goodOpT0FracLow,
  opT0Score
};

const std::vector<Cut> ncpizero_1p0pi_broad_cuts = {
  noCut,
  preSel,
  cosmicRej,
  muonRej,
  photonSel,
  pionRej,
  oneProtonSel,
  opT0Cuts,
};

const std::vector<Cut> ncpizero_Np0pi_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbs,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  noRazzledPions,
  atLeastOneRazzledProton,
  goodOpT0FracHigh,
  goodOpT0FracLow,
  opT0Score
};

const std::vector<Cut> ncpizero_Np0pi_broad_cuts = {
  noCut,
  preSel,
  cosmicRej,
  muonRej,
  photonSel,
  pionRej,
  atLeastOneProtonSel,
  opT0Cuts,
};

const std::vector<Cut> ncpizero_Xp0pi_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbs,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  noRazzledPions,
  goodOpT0FracHigh,
  goodOpT0FracLow,
  opT0Score
};

const std::vector<Cut> ncpizero_Xp0pi_broad_cuts = {
  noCut,
  preSel,
  cosmicRej,
  muonRej,
  photonSel,
  pionRej,
  opT0Cuts,
};
