#pragma once

#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

std::vector<Plot> selection_plots = {
  { "slc_is_clear_cosmic", "slc_is_clear_cosmic", ";Is Clear Cosmic?;Slices",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_is_fv", "slc_is_fv", ";IsFV?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_crumbs_score", "slc_crumbs_score", ";CRUMBS Score;Slices",
    50, -1.5, 1. },
  { "slc_n_razzled_muons", "slc_n_razzled_muons", ";N Razzled Muons;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzled_pions_thresh", "slc_n_razzled_pions_thresh", ";N Razzled Pions;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzled_protons_thresh", "slc_n_razzled_protons_thresh", ";N Razzled Protons;Slices",
    5, -0.5, 4.5 },
  { "slc_n_shws", "slc_n_shws", ";N Showers;Slices",
    5, -0.5, 4.5 },
  { "slc_n_trks", "slc_n_trks", ";N Tracks;Slices",
    5, -0.5, 4.5 },
  { "slc_n_pfps", "slc_n_pfps", ";N PFPs;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzled_photons", "slc_n_razzled_photons", ";N Razzled Photons;Slices",
    5, -0.5, 4.5 },
  { "slc_best_pzc_good_kinematics", "slc_best_pzc_good_kinematics", ";Best #pi^{0} Candidates Good Kinematics?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_opt0_frac", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE", ";OpT0 Fraction;Slices",
    100, -2, 8 }, // 100, -2, 2 for cut optimisation
  { "slc_opt0_score", "slc_opt0_score", ";OpT0 Score;Slices",
    100, 0, 2e5 }, // 100, 0, 5e3 for cut optimisation
};

std::vector<Plot> extension_plots = {
  { "slc_opt0_frac", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE", ";OpT0 Fraction;Slices",
    100, -2, 8 },
  { "slc_opt0_frac_log", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE", ";OpT0 Fraction;Slices",
    100, -2, 8, kBlack, true },
  { "slc_opt0_frac_log_narrow", "(slc_opt0_hypPE - slc_opt0_measPE)/slc_opt0_measPE", ";OpT0 Fraction;Slices",
    100, -2, 2, kBlack, true },
  { "slc_opt0_score", "slc_opt0_score", ";OpT0 Score;Slices",
    100, 0, 2e5 },
  { "slc_opt0_score_log", "slc_opt0_score", ";OpT0 Score;Slices",
    100, 0, 2e5, kBlack, true },
  { "slc_opt0_score_narrow", "slc_opt0_score", ";OpT0 Score;Slices",
    100, 0, 2e3 },
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

std::vector<VarBinPlot> true_slc_observables = {
  { "slc_true_pz_pizero_mom", "slc_true_pz_pizero_mom*1e3", ";p_{#pi^{0}} (MeV/c);Events",
    8, { 0., 60., 120., 180., 240., 300., 400., 600., 1000. } },
  { "slc_true_pz_cos_theta_pizero", "slc_true_pz_cos_theta_pizero", ";cos(#theta_{#pi^{0}});Events",
    9, { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. } },
  { "slc_true_pz_cos_com", "slc_true_pz_cos_com", ";cos(#theta_{CoM});Events",
    10, { 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1. } },
  { "slc_true_pz_decay_asymmetry", "slc_true_pz_decay_asymmetry", ";Decay Asymmetry;Events",
    10, { 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1. } },
};

std::vector<Plot> observables = {
  { "slc_best_pzc_invariant_mass", "slc_best_pzc_invariant_mass", ";M_{#gamma#gamma} (MeV/c^{2});Events",
    25, 0, 500 },
  { "slc_best_pzc_pizero_mom", "slc_best_pzc_pizero_mom", ";p_{#pi^{0}} (MeV/c);Events",
    25, 0, 1000 },
  { "slc_best_pzc_cos_theta_pizero", "slc_best_pzc_cos_theta_pizero", ";cos(#theta_{#pi^{0}});Events",
    25, -1, 1 },
  { "slc_best_pzc_cos_com", "slc_best_pzc_cos_com", ";cos(#theta_{CoM});Events",
    25, 0, 1 },
  { "slc_best_pzc_decay_asymmetry", "slc_best_pzc_decay_asymmetry", ";Decay Asymmetry;Events",
    25, 0, 1 },
};

std::vector<Plot> no_plots = {};

std::vector<TwoDPlotSet> true_observables_twod_sets = {
  { "pizero_momentum_and_cos_theta", "1e3 * nu_pz_pizero_mom", "nu_pz_cos_theta_pizero",
    8, { 0., 60., 120., 180., 240., 300., 400., 600., 1000. }, "p_{#pi^{0}} (MeV/c)", "60 MeV/c", 60.,
    9, { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. }, "cos(#theta_{#pi^{0}} )", "0.2", 0.2 },
  { "cos_com_and_cos_theta", "nu_pz_cos_com", "nu_pz_cos_theta_pizero",
    10, { 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1. }, "cos(#theta_{CoM} )", "0.1", 0.1,
    9, { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. }, "cos(#theta_{#pi^{0}} )", "0.2", 0.2 },
  { "cos_theta_and_pizero_momentum", "nu_pz_cos_theta_pizero", "1e3 * nu_pz_pizero_mom",
    8, { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 1. }, "cos(#theta_{#pi^{0}} )", "0.2", 0.2,
    9, { 0., 60., 120., 180., 240., 300., 400., 500., 600., 1000. }, "p_{#pi^{0}} (MeV/c)", "60 MeV/c", 60. },
};

const TwoDPlotSet observable_set = { "pizero_momentum_and_cos_theta", "slc_best_pzc_pizero_mom", "slc_best_pzc_cos_theta_pizero",
                                     8, { 0., 60., 120., 180., 240., 300., 400., 600., 1000. }, "p_{#pi^{0}} (MeV/c)", " MeV/c", 1.,
                                     9, { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. }, "cos(#theta_{#pi^{0}} )" };

std::vector<TwoDPlotSet> observables_twod_sets = { observable_set };

std::vector<Plot> reco_eff_plots = {
  { "energy", "mc_energy0*1e3", ";E (MeV);", 0, 0, 0 },
  { "momentum", "mc_momentum*1e3", ";p (MeV/c);", 0, 0, 0 },
};
