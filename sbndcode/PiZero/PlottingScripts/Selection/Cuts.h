#pragma once

#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

const Cut noCut                    = { "no_cut", "", "No Cut" };
const Cut notClearCosmic           = { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" };
const Cut fv                       = { "fv", "slc_is_fv", "FV" };
const Cut crumbsIncl               = { "crumbs", "slc_crumbs_score>-0.1", "CRUMBS Cut" };
const Cut crumbs0p0pi              = { "crumbs", "slc_crumbs_score>-0.145", "CRUMBS Cut" };
const Cut crumbsNp0pi              = { "crumbs", "slc_crumbs_score>-0.085", "CRUMBS Cut" };
const Cut noRazzledMuons           = { "no_razzled_muons", "slc_n_primary_razzled_muons==0", "No Razzled Muons" };
const Cut atLeastTwoPFPs           = { "at_least_two_pfps", "slc_n_pfps>1", "At Least Two PFPs" };
const Cut atLeastTwoRazzledPhotons = { "has_two_razzled_photons", "slc_n_primary_razzled_photons>1", "Has Two Razzled Photons" };
const Cut goodPiZeroKinematics     = { "good_pizero_kinematics", "slc_best_pzc_good_kinematics", "Good PiZero Kinematics" };
const Cut noRazzledPions           = { "no_razzled_pions", "slc_n_primary_razzled_pions_thresh==0", "No Razzled Pions" };
const Cut noRazzledProtons         = { "no_razzled_protons", "slc_n_primary_razzled_protons_thresh==0", "No Razzled Protons" };
const Cut atLeastOneRazzledProton  = { "at_least_one_razzled_proton", "slc_n_primary_razzled_protons_thresh>0", "At Least One Razzled Proton" };
const Cut goodOpT0FracHighIncl     = { "good_opt0_frac_high", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE<0.756", "Good OpT0 Frac High" };
const Cut goodOpT0FracHigh0p0pi    = { "good_opt0_frac_high", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE<0.408", "Good OpT0 Frac High" };
const Cut goodOpT0FracHighNp0pi    = { "good_opt0_frac_high", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE<0.792", "Good OpT0 Frac High" };
const Cut goodOpT0FracLowIncl      = { "good_opt0_frac_low", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE>-0.7", "Good OpT0 Frac Low" };
const Cut goodOpT0FracLow0p0pi     = { "good_opt0_frac_low", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE>-0.704", "Good OpT0 Frac Low" };
const Cut goodOpT0FracLowNp0pi     = { "good_opt0_frac_low", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE>-0.424", "Good OpT0 Frac Low" };
const Cut opT0ScoreIncl            = { "opt0_score", "slc_opt0_score>70", "OpT0 Score" };
const Cut opT0Score0p0pi           = { "opt0_score", "slc_opt0_score>145", "OpT0 Score" };
const Cut opT0ScoreNp0pi           = { "opt0_score", "slc_opt0_score>175", "OpT0 Score" };

const Cut preSel              = { "presel", "!slc_is_clear_cosmic && slc_is_fv", "Pre-Sel", kOrange-2 };
const Cut cosmicRejIncl       = { "cosmic_rejection", "slc_crumbs_score>-0.1", "Cosmic Rej", kGreen+1 };
const Cut cosmicRej0p0pi      = { "cosmic_rejection", "slc_crumbs_score>-0.145", "Cosmic Rej", kGreen+1 };
const Cut cosmicRejNp0pi      = { "cosmic_rejection", "slc_crumbs_score>-0.085", "Cosmic Rej", kGreen+1 };
const Cut muonRej             = { "muon_rejection", "slc_n_primary_razzled_muons ==0", "Muon Rej", kRed-9 };
const Cut photonSel           = { "photons_selection", "slc_n_pfps>1 && slc_n_primary_razzled_photons>1 && slc_best_pzc_good_kinematics", "Photon Sel", kBlue-9 };
const Cut pionRej             = { "pion_rejection", "slc_n_primary_razzled_pions_thresh==0", "Pion Rej", kViolet - 5 };
const Cut protonRej           = { "proton_rejection", "slc_n_primary_razzled_protons_thresh==0", "Proton Rej", kViolet - 6 };
const Cut atLeastOneProtonSel = { "proton_selection", "slc_n_primary_razzled_protons_thresh>0", "Proton Sel", kViolet - 6 };
const Cut opT0CutsIncl        = { "opt0_cuts", goodOpT0FracHighIncl.cut + goodOpT0FracLowIncl.cut + opT0ScoreIncl.cut, "OpT0 Cuts" };
const Cut opT0Cuts0p0pi       = { "opt0_cuts", goodOpT0FracHigh0p0pi.cut + goodOpT0FracLow0p0pi.cut + opT0Score0p0pi.cut, "OpT0 Cuts" };
const Cut opT0CutsNp0pi       = { "opt0_cuts", goodOpT0FracHighNp0pi.cut + goodOpT0FracLowNp0pi.cut + opT0ScoreNp0pi.cut, "OpT0 Cuts" };

const std::vector<Cut> ncpizero_incl_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbsIncl,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  goodOpT0FracHighIncl,
  goodOpT0FracLowIncl,
  opT0ScoreIncl
};

const std::vector<Cut> ncpizero_incl_broad_cuts = {
  noCut,
  preSel,
  cosmicRejIncl,
  muonRej,
  photonSel,
  opT0CutsIncl
};

const std::vector<Cut> ncpizero_0p0pi_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbs0p0pi,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  noRazzledPions,
  noRazzledProtons,
  goodOpT0FracHigh0p0pi,
  goodOpT0FracLow0p0pi,
  opT0Score0p0pi
};

const std::vector<Cut> ncpizero_0p0pi_broad_cuts = {
  noCut,
  preSel,
  cosmicRej0p0pi,
  muonRej,
  photonSel,
  pionRej,
  protonRej,
  opT0Cuts0p0pi,
};

const std::vector<Cut> ncpizero_Np0pi_cuts = {
  noCut,
  notClearCosmic,
  fv,
  crumbsNp0pi,
  noRazzledMuons,
  atLeastTwoPFPs,
  atLeastTwoRazzledPhotons,
  goodPiZeroKinematics,
  noRazzledPions,
  atLeastOneRazzledProton,
  goodOpT0FracHighNp0pi,
  goodOpT0FracLowNp0pi,
  opT0ScoreNp0pi
};

const std::vector<Cut> ncpizero_Np0pi_broad_cuts = {
  noCut,
  preSel,
  cosmicRejNp0pi,
  muonRej,
  photonSel,
  pionRej,
  atLeastOneProtonSel,
  opT0CutsNp0pi,
};