#pragma once

#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

std::vector<Plot> selection_plots = {
  { "slc_is_clear_cosmic", "slc_is_clear_cosmic", ";Is Clear Cosmic?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_is_fv", "slc_is_fv", ";IsFV?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_crumbs_nc_score", "slc_crumbs_nc_score", ";CRUMBS NC Score;Candidates",
    50, -1.5, 1. },
  { "slc_n_razzled_muons", "slc_n_razzled_muons", ";N Razzled Muons;Candidates",
    5, -0.5, 4.5 },
  { "slc_n_razzled_pions_thresh", "slc_n_razzled_pions_thresh", ";N Razzled Pions;Candidates",
    5, -0.5, 4.5 },
  { "slc_n_razzled_protons_thresh", "slc_n_razzled_protons_thresh", ";N Razzled Protons;Candidates",
    5, -0.5, 4.5 },
  { "slc_n_shws", "slc_n_shws", ";N Showers;Candidates",
    5, -0.5, 4.5 },
  { "slc_n_trks", "slc_n_trks", ";N Tracks;Candidates",
    5, -0.5, 4.5 },
  { "slc_n_pfps", "slc_n_pfps", ";N PFPs;Candidates",
    5, -0.5, 4.5 },
  { "slc_n_razzled_photons", "slc_n_razzled_photons", ";N Razzled Photons;Candidates",
    5, -0.5, 4.5 },
  { "slc_best_corr_pzc_good_kinematics", "slc_best_corr_pzc_good_kinematics", ";Best #pi^{0} Candidates Good Kinematics?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_opt0_fracPE", "slc_opt0_fracPE", ";OpT0 Fractional PE Difference;Candidates",
    100, -2, 8 },
  { "slc_opt0_fracPE_log", "slc_opt0_fracPE", ";OpT0 Fractional PE Difference;Candidates",
    100, -2, 8, kBlack, true },
  { "slc_opt0_score", "slc_opt0_score", ";OpT0 Score;Candidates",
    100, 0, 2e5 }, // 100, 0, 5e3 for cut optimisation
};

std::vector<Plot> selection_plots_tmp = {
  { "slc_crumbs_nc_score", "slc_crumbs_nc_score", ";CRUMBS NC Score;Candidates",
    50, -1.5, 1. },
};

std::vector<Plot> optimisation_plots = {
  { "slc_crumbs_score", "slc_crumbs_score", ";CRUMBS Score;Candidates",
    50, -1.5, 1. },
  { "slc_crumbs_nc_score", "slc_crumbs_nc_score", ";CRUMBS NC Score;Candidates",
    50, -1.5, 1. },
  { "slc_opt0_fracPE", "slc_opt0_fracPE", ";OpT0 Fractional PE Difference;Candidates",
    100, -2, 2 },
  { "slc_opt0_score", "slc_opt0_score", ";OpT0 Score;Candidates",
    100, 0, 5e3 },
  { "slc_crumbs_ccnue_score", "slc_crumbs_ccnue_score", ";CRUMBS CCNuE Score;Candidates",
    50, -1.5, 1. },
};
std::vector<Plot> extension_plots = {
  { "slc_opt0_fracPE", "slc_opt0_fracPE", ";OpT0 Fractional PE Difference;Candidates",
    100, -2, 8 },
  { "slc_opt0_fracPE_log", "slc_opt0_fracPE", ";OpT0 Fractional PE Difference;Candidates",
    100, -2, 8, kBlack, true },
  { "slc_opt0_fracPE_log_narrow", "slc_opt0_fracPE", ";OpT0 Fractional PE Difference;Candidates",
    100, -2, 2, kBlack, true },
  { "slc_opt0_score", "slc_opt0_score", ";OpT0 Score;Candidates",
    100, 0, 2e5 },
  { "slc_opt0_score_log", "slc_opt0_score", ";OpT0 Score;Candidates",
    100, 0, 2e5, kBlack, true },
  { "slc_opt0_score_narrow", "slc_opt0_score", ";OpT0 Score;Candidates",
    100, 0, 2e3 },
};

const std::vector<double> pizeromombins      = { 0., 60., 100, 140., 180., 220., 260., 300., 350., 400., 500., 600., 1000. };
const std::vector<double> costhetapizerobins = { -1., -0.5, 0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };
const std::vector<double> coscombins         = { 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1. };
const std::vector<double> decayasymbins      = { 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1. };
const std::vector<double> openanglebins      = { 0., 20., 40., 60., 80., 100., 120., 140., 180. };
const std::vector<double> invmassbins        = { 0, 40, 60, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200, 250, 350, 500 };

const std::vector<double> pizeromombins2d      = { 0., 60., 120, 180., 240., 300., 400., 600., 1000. };
const std::vector<double> costhetapizerobins2d = { -1., -0.5, 0, 0.4, 0.65, 0.9, 1. };

std::vector<VarBinPlot> true_observables_var_bin = {
  { "nu_pz_pizero_mom", "nu_pz_pizero_mom*1e3", ";True p_{#pi^{0}} (MeV/c);Events / 60 MeV/c",
    pizeromombins.size() - 1, pizeromombins, 60. },
  { "nu_pz_cos_theta_pizero", "nu_pz_cos_theta_pizero", ";True cos(#theta_{#pi^{0}});Events / 0.1",
    costhetapizerobins.size() - 1, costhetapizerobins, 0.1 },
  { "nu_pz_cos_com", "nu_pz_cos_com", ";True cos(#theta_{CoM});Events / 0.1",
    coscombins.size() - 1, coscombins, 0.1 },
  { "nu_pz_decay_asymmetry", "nu_pz_decay_asymmetry", ";True |E_{#gamma_{1}} - E_{#gamma_{2}}| / E_{#gamma_{1}} + E_{#gamma_{2}};Events / 0.1",
    decayasymbins.size() - 1, decayasymbins, 0.1 },
  { "nu_pz_open_angle", "nu_pz_open_angle", ";True #theta_{#gamma#gamma} (#circ);Events / 20#circ",
    openanglebins.size() - 1, openanglebins, 20. },
};

std::vector<VarBinPlot> true_slc_observables = {
  { "slc_true_pz_pizero_mom", "slc_true_pz_pizero_mom*1e3", ";True p_{#pi^{0}} (MeV/c);Candidates / 60 MeV/c",
    pizeromombins.size() - 1, pizeromombins, 60. },
  { "slc_true_pz_cos_theta_pizero", "slc_true_pz_cos_theta_pizero", ";True cos(#theta_{#pi^{0}});Candidates / 0.1",
    costhetapizerobins.size() - 1, costhetapizerobins, 0.1 },
  { "slc_true_pz_cos_com", "slc_true_pz_cos_com", ";True cos(#theta_{CoM});Candidates / 0.1",
    coscombins.size() - 1, coscombins, 0.1 },
  { "slc_true_pz_decay_asymmetry", "slc_true_pz_decay_asymmetry", ";True |E_{#gamma_{1}} - E_{#gamma_{2}}| / E_{#gamma_{1}} + E_{#gamma_{2}};Candidates / 0.1",
    decayasymbins.size() - 1, decayasymbins, 0.1 },
  { "slc_true_pz_open_angle", "slc_true_pz_open_angle", ";True #theta_{#gamma#gamma} (#circ);Candidates / 20#circ",
    openanglebins.size() - 1, openanglebins, 20. },
};

std::vector<VarBinPlot> observables = {
  { "slc_best_corr_pzc_invariant_mass", "slc_best_corr_pzc_invariant_mass", ";M_{#gamma#gamma} (MeV/c^{2});Candidates / 10 MeV/c^{2}",
    invmassbins.size() - 1, invmassbins, 10 },
  { "slc_best_corr_pzc_pizero_mom", "slc_best_corr_pzc_pizero_mom", ";p_{#pi^{0}} (MeV/c);Candidates / 60 MeV/c",
    pizeromombins.size() - 1, pizeromombins, 60. },
  { "slc_best_corr_pzc_cos_theta_pizero", "slc_best_corr_pzc_cos_theta_pizero", ";cos(#theta_{#pi^{0}});Candidates / 0.1",
    costhetapizerobins.size() - 1, costhetapizerobins, 0.1 },
  { "slc_best_corr_pzc_cos_com", "slc_best_corr_pzc_cos_com", ";cos(#theta_{CoM});Candidates / 0.1",
    coscombins.size() - 1, coscombins, 0.1 },
  { "slc_best_corr_pzc_decay_asymmetry", "slc_best_corr_pzc_decay_asymmetry", ";Decay Asymmetry;Candidates / 0.1",
    decayasymbins.size() - 1, decayasymbins, 0.1 },
};

std::vector<VarBinPlot> observables_corr = {
  { "slc_best_corr_pzc_invariant_mass_corr", "slc_best_corr_pzc_invariant_mass_corr", ";M_{#gamma#gamma} (MeV/c^{2});Candidates / 10 MeV/c^{2}",
    invmassbins.size() - 1, invmassbins, 10 },
  { "slc_best_corr_pzc_pizero_mom_corr", "slc_best_corr_pzc_pizero_mom_corr", ";p_{#pi^{0}} (MeV/c);Candidates / 60 MeV/c",
    pizeromombins.size() - 1, pizeromombins, 60. },
  { "slc_best_corr_pzc_cos_theta_pizero_corr", "slc_best_corr_pzc_cos_theta_pizero_corr", ";cos(#theta_{#pi^{0}});Candidates / 0.1",
    costhetapizerobins.size() - 1, costhetapizerobins, 0.1 },
  { "slc_best_corr_pzc_cos_com_corr", "slc_best_corr_pzc_cos_com_corr", ";cos(#theta_{CoM});Candidates / 0.1",
    coscombins.size() - 1, coscombins, 0.1 },
  { "slc_best_corr_pzc_decay_asymmetry_corr", "slc_best_corr_pzc_decay_asymmetry_corr", ";Decay Asymmetry;Candidates / 0.1",
    decayasymbins.size() - 1, decayasymbins, 0.1 },
};

std::vector<Plot> no_plots = {};

std::vector<TwoDPlotSet> true_observables_twod_sets = {
  { "pizero_momentum_and_cos_theta", "1e3 * nu_pz_pizero_mom", "nu_pz_cos_theta_pizero",
    pizeromombins2d.size() - 1, pizeromombins2d, "True p_{#pi^{0}} (MeV/c)", "60 MeV/c", 60.,
    costhetapizerobins2d.size() - 1, costhetapizerobins2d, "True cos(#theta_{#pi^{0}} )", "0.2", 0.2 },
};

const TwoDPlotSet observable_set = { "pizero_momentum_and_cos_theta", "slc_best_corr_pzc_pizero_mom", "slc_best_corr_pzc_cos_theta_pizero",
                                     pizeromombins2d.size() - 1, pizeromombins2d, "p_{#pi^{0}} (MeV/c)", " MeV/c", 1.,
                                     costhetapizerobins2d.size() - 1, costhetapizerobins2d, "cos(#theta_{#pi^{0}} )" };

std::vector<TwoDPlotSet> observables_twod_sets = { observable_set };

std::vector<Plot> reco_eff_plots = {
  { "energy", "mc_energy0*1e3", ";E (MeV);", 0, 0, 0 },
  { "momentum", "mc_momentum*1e3", ";p (MeV/c);", 0, 0, 0 },
};
