#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

void InputVars()
{
  const TString save_dir = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv12/crumbs/training/inputvars";
  gSystem->Exec("mkdir -p " + save_dir);

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *slices = new TChain("crumbs/SliceTree");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  std::vector<Plot> plots = { { "tpc_CRFracHitsInLongestTrack", "tpc_CRFracHitsInLongestTrack",
                                ";CRFracHitsInLongestTrack;Slices (normalised)",
                                50, 0, 1 },
                              { "tpc_CRLongestTrackDeflection", "tpc_CRLongestTrackDeflection",
                                ";CRLongestTrackDeflection;Slices (normalised)",
                                50, 0, 1 },
                              { "tpc_CRLongestTrackDirY", "tpc_CRLongestTrackDirY",
                                ";CRLongestTrackDirY;Slices (normalised)",
                                50, -1, .5 },
                              { "tpc_CRNHitsMax", "tpc_CRNHitsMax", ";CRNHitsMax;Slices (normalised)",
                                50, 0, 3000 },
                              { "tpc_NuEigenRatioInSphere", "tpc_NuEigenRatioInSphere",
                                ";NuEigenRatioInSphere;Slices (normalised)",
                                50, 0, 1 },
                              { "tpc_NuNFinalStatePfos", "tpc_NuNFinalStatePfos",
                                ";NuNFinalStatePfos;Slices (normalised)",
                                10, -0.5, 9.5 },
                              { "tpc_NuNHitsTotal", "tpc_NuNHitsTotal", ";NuNHitsTotal;Slices (normalised)",
                                50, 0, 5000 },
                              { "tpc_NuNSpacePointsInSphere", "tpc_NuNSpacePointsInSphere",
                                ";NuNSpacePointsInSphere;Slices (normalised)",
                                50, 0, 300 },
                              { "tpc_NuVertexY", "tpc_NuVertexY", ";NuVertexY (cm);Slices (normalised)",
                                50, -500, 500 },
                              { "tpc_NuWeightedDirZ", "tpc_NuWeightedDirZ",
                                ";NuWeightedDirZ;Slices (normalised)",
                                50, -1, 1 },
                              { "tpc_StoppingChi2CosmicRatio", "tpc_StoppingChi2CosmicRatio",
                                ";StoppingChi2CosmicRatio;Slices (normalised)",
                                50, -5, 10 },
                              { "pds_FMTotalScore", "pds_FMTotalScore", ";FMTotalScore;Slices (normalised)",
                                50, -10, 60 },
                              { "pds_FMPE", "pds_FMPE", ";FMPE;Slices (normalised)",
                                50, -10000, 40000 },
                              { "pds_FMTime", "pds_FMTime", ";FMTime (#mus);Slices (normalised)",
                                50, -500, 200 },
                              { "crt_TrackScore", "crt_TrackScore", ";CRTTrackScore;Slices (normalised)",
                                50, -4, 200 },
                              { "crt_TrackScore_log", "crt_TrackScore", ";CRTTrackScore;Slices (normalised)",
                                50, -4, 200, kBlack, true },
                              { "crt_SPScore", "crt_SPScore", ";CRTSPScore (cm);Slices (normalised)",
                                50, -4, 100 },
                              { "crt_SPScore_log", "crt_SPScore", ";CRTSPScore (cm);Slices (normalised)",
                                50, -4, 100, kBlack, true },
                              { "crt_TrackTime", "crt_TrackTime",  ";CRTTrackTime (#mus);Slices (normalised)",
                                50, -3000, 1500 },
                              { "crt_TrackTime_log", "crt_TrackTime",  ";CRTTrackTime (#mus);Slices (normalised)",
                                50, -3000, 1500, kBlack, true },
                              { "crt_SPTime", "crt_SPTime", ";CRTSPTime (#mus);Slices (normalised)",
                                50, -3000, 1500 },
                              { "crt_SPTime_log", "crt_SPTime", ";CRTSPTime (#mus);Slices (normalised)",
                                50, -3000, 1500, kBlack, true },
                              { "pds_OpT0Score", "pds_OpT0Score", ";OpT0Score;Slices (normalised)",
                                50, -5000, 40000 },
                              { "pds_OpT0MeasuredPE", "pds_OpT0MeasuredPE",
                                ";OpT0MeasuredPE;Slices (normalised)",
                                50, -10000, 100000 },
  };

  std::vector<Cut> cuts = { { "signal", "strstr(matchedType,\"Nu\") && !strstr(matchedType,\"DirtNu\") && matchedPurity > 0.8 && matchedCompleteness > 0.8", "True Neutrino Slices", kBlue+2 },
                            { "background", "!strstr(matchedType,\"Nu\")", "Cosmic Slices", kRed+2 }
  };

  for(auto const& plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);

      MakeComparisonPlot(canvas, slices, plot, cuts, {.25, .82, .8, .87});

      canvas->SaveAs(save_dir + "/" + plot.name + ".png");
      canvas->SaveAs(save_dir + "/" + plot.name + ".pdf");
    }
}
