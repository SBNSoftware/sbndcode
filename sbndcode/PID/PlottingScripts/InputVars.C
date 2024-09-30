#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void InputVars()
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv14/razzled/training/inputvars";
  gSystem->Exec("mkdir -p " + saveDir);

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("razzled/pfpTree");
  tree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv14/NCPiZeroAv14_rockbox.root");
  tree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv14/NCPiZeroAv14_intrnue.root");
  tree->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv14/NCPiZeroAv14_intime.root");

  std::vector<Cut> particles = { {"electron", "abs(truePDG)==11 && !unambiguousSlice", "e^{#pm}", kMagenta+2},
                                 {"muon", "abs(truePDG)==13 && !unambiguousSlice", "#mu^{#pm}", kRed+2},
                                 {"photon", "abs(truePDG)==22 && !unambiguousSlice", "#gamma", kBlue+2},
                                 {"pion", "abs(truePDG)==211 && !unambiguousSlice", "#pi^{#pm}", kGreen+2},
                                 {"proton", "abs(truePDG)==2212 && !unambiguousSlice", "p", kOrange+2},
                                 //                              {"other", "abs(truePDG)!=11 && abs(truePDG)!=13 && abs(truePDG)!=22 && abs(truePDG)!=211 && abs(truePDG)!=2212 && !unambiguousSlice", "Other", kBlack}
  };

  TCut quality_cut = "trk_length > 5 && showerEnergy > 10 && trackContained && showerContained";

  for(auto& particle : particles)
    particle.cut += quality_cut;

  std::vector<Plot> plots = { {"pfp_numDaughters", "pfp_numDaughters", ";PFP N Daughters;PFPs",
                               19, -8.5, 8.5},
                              {"pfp_numDaughters_filled", "pfp_numDaughters", ";PFP N Daughters;PFPs",
                               9, -.5, 8.5},
                              {"pfp_maxDaughterHits", "pfp_maxDaughterHits", ";PFP Max Daughter Hits;PFPs",
                               100, -100, 400},
                              {"pfp_maxDaughterHits_filled", "pfp_maxDaughterHits", ";PFP Max Daughter Hits;PFPs",
                               80, 0, 400},
                              {"pfp_trackScore", "pfp_trackScore", ";PFP Track Score;PFPs",
                               300, -1.5, 1.5},
                              {"pfp_trackScore_filled", "pfp_trackScore", ";PFP Track Score;PFPs",
                               140, -.2, 1.2},
                              {"pfp_chargeEndFrac", "pfp_chargeEndFrac", ";PFP Charge End Fraction;PFPs",
                               300, -1.5, 1.5},
                              {"pfp_chargeEndFrac_filled", "pfp_chargeEndFrac", ";PFP Charge End Fraction;PFPs",
                               140, -.2, 1.2},
                              {"pfp_chargeFracSpread", "pfp_chargeFracSpread", ";PFP Charge Fraction Spread;PFPs",
                               200, -1.5, 2.5},
                              {"pfp_chargeFracSpread_filled", "pfp_chargeFracSpread", ";PFP Charge Fraction Spread;PFPs",
                               120, -.2, 2.2},
                              {"pfp_linearFitDiff", "pfp_linearFitDiff", ";PFP Linear Fit Difference;PFPs",
                               180, -1.2, .6},
                              {"pfp_linearFitDiff_filled", "pfp_linearFitDiff", ";PFP Linear Fit Difference;PFPs",
                               80, -.1, .3},
                              {"pfp_linearFitLength", "pfp_linearFitLength", ";PFP Linear Fit Length (cm);PFPs",
                               60, -100, 500},
                              {"pfp_linearFitLength_filled", "pfp_linearFitLength", ";PFP Linear Fit Length (cm);PFPs",
                               51, -10, 500},
                              {"pfp_linearFitGapLength", "pfp_linearFitGapLength", ";PFP Linear Fit Gap Length (cm);PFPs",
                               300, -1.5, 1.5},
                              {"pfp_linearFitGapLength_filled", "pfp_linearFitGapLength", ";PFP Linear Fit Gap Length (cm);PFPs",
                               140, -.2, 1.2},
                              {"pfp_linearFitRMS", "pfp_linearFitRMS", ";PFP Linear Fit RMS;PFPs",
                               250, -1.5, 3.5},
                              {"pfp_linearFitRMS_filled", "pfp_linearFitRMS", ";PFP Linear Fit RMS;PFPs",
                               170, -.2, 3.2},
                              {"pfp_openAngleDiff", "pfp_openAngleDiff", ";PFP Open Angle Difference (#circ);PFPs",
                               60, -20, 100},
                              {"pfp_openAngleDiff_filled", "pfp_openAngleDiff", ";PFP Open Angle Difference (#circ);PFPs",
                               50, 0, 100},
                              {"pfp_secondaryPCARatio", "pfp_secondaryPCARatio", ";PFP Secondary PCA Ratio;PFPs",
                               300, -1.5, 1.5},
                              {"pfp_secondaryPCARatio_log", "pfp_secondaryPCARatio", ";PFP Secondary PCA Ratio;PFPs",
                               300, -1.5, 1.5, kBlack, true},
                              {"pfp_secondaryPCARatio_filled", "pfp_secondaryPCARatio", ";PFP Secondary PCA Ratio;PFPs",
                               140, -.2, 1.2},
                              {"pfp_secondaryPCARatio_filled_log", "pfp_secondaryPCARatio", ";PFP Secondary PCA Ratio;PFPs",
                               140, -.2, 1.2, kBlack, true},
                              {"pfp_tertiaryPCARatio", "pfp_tertiaryPCARatio", ";PFP Tertiary PCA Ratio;PFPs",
                               300, -1.5, 1.5},
                              {"pfp_tertiaryPCARatio_log", "pfp_tertiaryPCARatio", ";PFP Tertiary PCA Ratio;PFPs",
                               300, -1.5, 1.5, kBlack, true},
                              {"pfp_tertiaryPCARatio_filled", "pfp_tertiaryPCARatio", ";PFP Tertiary PCA Ratio;PFPs",
                               140, -.2, 1.2},
                              {"pfp_tertiaryPCARatio_filled_log", "pfp_tertiaryPCARatio", ";PFP Tertiary PCA Ratio;PFPs",
                               140, -.2, 1.2, kBlack, true},
                              {"pfp_vertexDist", "pfp_vertexDist", ";PFP Vertex Distance (cm);PFPs",
                               250, -100, 400},
                              {"pfp_vertexDist_log", "pfp_vertexDist", ";PFP Vertex Distance (cm);PFPs",
                               250, -100, 400, kBlack, true},
                              {"pfp_vertexDist_filled", "pfp_vertexDist", ";PFP Vertex Distance (cm);PFPs",
                               205, -10, 400},
                              {"pfp_vertexDist_filled_log", "pfp_vertexDist", ";PFP Vertex Distance (cm);PFPs",
                               205, -10, 400, kBlack, true},
                              {"trk_length", "trk_length", ";Track Length (cm);PFPs",
                               100, -200, 800},
                              {"trk_length_filled", "trk_length", ";Track Length (cm);PFPs",
                               85, -50, 800},
                              {"trk_chi2PIDMuon", "trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;PFPs",
                               60, -20, 100},
                              {"trk_chi2PIDMuon_filled", "trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;PFPs",
                               52, -4, 100},
                              {"trk_chi2PIDPion", "trackChi2PIDPion", ";Track #chi^{2} PID Pion;PFPs",
                               60, -20, 100},
                              {"trk_chi2PIDKaon", "trackChi2PIDKaon", ";Track #chi^{2} PID Kaon;PFPs",
                               60, -20, 100},
                              {"trk_chi2PIDProton", "trk_chi2PIDProton", ";Track #chi^{2} PID Proton;PFPs",
                               84, -20, 400},
                              {"trk_chi2PIDProton_filled", "trk_chi2PIDProton", ";Track #chi^{2} PID Proton;PFPs",
                               81, -5, 400},
                              {"trk_chi2PIDMuonPionDiff", "trk_chi2PIDMuonPionDiff", ";Track #chi^{2} PID Muon Pion Difference;PFPs",
                               90, -120, 60},
                              {"trk_chi2PIDMuonPionDiff_filled", "trk_chi2PIDMuonPionDiff", ";Track #chi^{2} PID Muon Pion Difference;PFPs",
                               60, -60, 60},
                              {"trk_mcsScatterMean", "trk_mcsScatterMean", ";Track MCS Mean (mrad);PFPs",
                               60, -200, 1000},
                              {"trk_mcsScatterMean_log", "trk_mcsScatterMean", ";Track MCS Mean (mrad);PFPs",
                               60, -200, 1000, kBlack, true},
                              {"trk_mcsScatterMean_filled", "trk_mcsScatterMean", ";Track MCS Mean (mrad);PFPs",
                               50, 0, 1000},
                              {"trk_mcsScatterMean_filled_log", "trk_mcsScatterMean", ";Track MCS Mean (mrad);PFPs",
                               50, 0, 1000, kBlack, true},
                              {"trk_mcsScatterMaxRatio", "trk_mcsScatterMaxRatio", ";Track MCS Max Ratio;PFPs",
                               200, -1.5, 2.5},
                              {"trk_mcsScatterMaxRatio_filled", "trk_mcsScatterMaxRatio", ";Track MCS Max Ratio;PFPs",
                               120, -.2, 2.2},
                              {"trk_meanDCA", "trk_meanDCA", ";Track Mean DCA (cm);PFPs",
                               70, -10, 60},
                              {"trk_meanDCA_log", "trk_meanDCA", ";Track Mean DCA (cm);PFPs",
                               70, -10, 60, kBlack, true},
                              {"trk_meanDCA_filled", "trk_meanDCA", ";Track Mean DCA (cm);PFPs",
                               60, 0, 60},
                              {"trk_meanDCA_filled_log", "trk_meanDCA", ";Track Mean DCA (cm);PFPs",
                               60, 0, 60, kBlack, true},
                              {"trk_stoppingdEdxChi2Ratio", "trk_stoppingdEdxChi2Ratio", ";Track Stopping dE/dx Fit #chi^{2} Ratio;PFPs",
                               50, -10, 15},
                              {"trk_stoppingdEdxChi2Ratio_log", "trk_stoppingdEdxChi2Ratio", ";Track Stopping dE/dx Fit #chi^{2} Ratio;PFPs",
                               50, -10, 15, kBlack, true},
                              {"trk_stoppingdEdxChi2Ratio_filled", "trk_stoppingdEdxChi2Ratio", ";Track Stopping dE/dx Fit #chi^{2} Ratio;PFPs",
                               30, 0, 15},
                              {"trk_stoppingdEdxChi2Ratio_filled_log", "trk_stoppingdEdxChi2Ratio", ";Track Stopping dE/dx Fit #chi^{2} Ratio;PFPs",
                               30, 0, 15, kBlack, true},
                              {"trk_chi2Pol0dEdxFit", "trk_chi2Pol0dEdxFit", ";Track dE/dx Pol0 Fit Value (MeV/cm);PFPs",
                               50, -10, 15},
                              {"trk_chi2Pol0dEdxFit_filled", "trk_chi2Pol0dEdxFit", ";Track dE/dx Pol0 Fit Value (MeV/cm);PFPs",
                               30, 0, 15},
                              {"trk_momDiff", "trk_momDiff", ";Track (p_{MCS} - p_{Range}) / p_{Range};PFPs",
                               90, -20, 70},
                              {"trk_momDiff_filled", "trk_momDiff", ";Track (p_{MCS} - p_{Range}) / p_{Range};PFPs",
                               70, 0, 70},
                              {"shw_bestdEdx", "shw_bestdEdx", ";Shower Best Plane dE/dx (MeV/cm);PFPs",
                               60, -10, 20},
                              {"shw_bestdEdx_filled", "shw_bestdEdx", ";Shower Best Plane dE/dx (MeV/cm);PFPs",
                               40, 0, 20},
                              {"shw_convGap", "shw_convGap", ";Shower Conversion Gap (cm);PFPs",
                               60, -10, 50},
                              {"shw_convGap_log", "shw_convGap", ";Shower Conversion Gap (cm);PFPs",
                               60, -10, 50, kBlack, true},
                              {"shw_convGap_filled", "shw_convGap", ";Shower Conversion Gap (cm);PFPs",
                               52, -2, 50},
                              {"shw_convGap_filled_log", "shw_convGap", ";Shower Conversion Gap (cm);PFPs",
                               52, -2, 50, kBlack, true},
                              {"shw_openAngle", "shw_openAngle", ";Shower Opening Angle (#circ);PFPs",
                               60, -20, 100},
                              {"shw_openAngle_log", "shw_openAngle", ";Shower Opening Angle (#circ);PFPs",
                               60, -20, 100, kBlack, true},
                              {"shw_openAngle_filled", "shw_openAngle", ";Shower Opening Angle (#circ);PFPs",
                               52, -4, 100},
                              {"shw_openAngle_filled_log", "shw_openAngle", ";Shower Opening Angle (#circ);PFPs",
                               52, -4, 100, kBlack, true},
                              {"shw_modHitDensity", "shw_modHitDensity", ";Shower |Hit Density| (hits/wire);PFPs",
                               50, -10, 40},
                              {"shw_modHitDensity_filled", "shw_modHitDensity", ";Shower |Hit Density| (hits/wire);PFPs",
                               44, -4, 40},
                              {"shw_sqrtEnergyDensity", "shw_sqrtEnergyDensity", ";Shower #sqrt{Energy Density} (#sqrt{MeV}/cm);PFPs",
                               150, -6, 4},
                              {"shw_sqrtEnergyDensity_filled", "shw_sqrtEnergyDensity", ";Shower #sqrt{Energy Density} (#sqrt{MeV}/cm);PFPs",
                               75, -1, 4},
  };
  
  for(auto &plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name + "_byType",
                                    "c_" + plot.name + "_byType");
      canvas->cd();

      MakeComparisonPlot(canvas, tree, plot, particles, {.35, .77, .8, .87});

      canvas->SaveAs(saveDir + "/" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
          
      delete canvas;
    }
}
