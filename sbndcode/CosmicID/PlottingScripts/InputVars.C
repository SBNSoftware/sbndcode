#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

void InputVars()
{
  const TString save_dir = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv12/crumbs/training/inputvars/thesis";
  gSystem->Exec("mkdir -p " + save_dir);

  gROOT->SetStyle("henrySBND");
  gStyle->SetTitleSize(0.055, "xyz");
  gStyle->SetLabelSize(0.055, "xyz");
  gStyle->SetTitleOffset(1.1, "x");
  gROOT->ForceStyle();

  TChain *slices = new TChain("crumbs/SliceTree");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  /*  std::vector<Plot> plots = { { "tpc_CRFracHitsInLongestTrack", "tpc_CRFracHitsInLongestTrack",
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
  */
  std::vector<Plot> plots = { { "tpc_CRFracHitsInLongestTrack", "tpc_CRFracHitsInLongestTrack",
                                ";Fraction of Hits in Longest Track (Cosmic Reco);Slices (normalised)",
                                50, 0, 1 },
                              { "tpc_CRLongestTrackDeflection", "tpc_CRLongestTrackDeflection",
                                ";Longest Track Deflection (Cosmic Reco);Slices (normalised)",
                                50, 0, 1 },
                              { "tpc_CRLongestTrackDirY", "tpc_CRLongestTrackDirY",
                                ";Longest Track Y Direction (Cosmic Reco);Slices (normalised)",
                                50, -1, .5 },
                              { "tpc_CRNHitsMax", "tpc_CRNHitsMax", ";Number of Space Points in Longest Track (Cosmic Reco);Slices (normalised)",
                                50, 0, 3000 },
                              { "tpc_NuEigenRatioInSphere", "tpc_NuEigenRatioInSphere",
                                ";Eigen Ratio in Sphere (Nu Reco);Slices (normalised)",
                                50, 0, 1 },
                              { "tpc_NuNFinalStatePfos", "tpc_NuNFinalStatePfos",
                                ";Number of Final State PFOs (Nu Reco);Slices (normalised)",
                                10, -0.5, 9.5 },
                              { "tpc_NuNHitsTotal", "tpc_NuNHitsTotal", ";Number of Space Points (Nu Reco);Slices (normalised)",
                                50, 0, 5000 },
                              { "tpc_NuNSpacePointsInSphere", "tpc_NuNSpacePointsInSphere",
                                ";Number of Space Points in Sphere (Nu Reco);Slices (normalised)",
                                50, 0, 300 },
                              { "tpc_NuVertexY", "tpc_NuVertexY", ";Vertex Y (Nu Reco) (cm);Slices (normalised)",
                                50, -500, 500 },
                              { "tpc_NuWeightedDirZ", "tpc_NuWeightedDirZ",
                                ";Weighted Z Direction (Nu Reco);Slices (normalised)",
                                50, -1, 1 },
                              { "tpc_StoppingChi2CosmicRatio", "tpc_StoppingChi2CosmicRatio",
                                ";Stopping Fits Chi2 Ratio;Slices (normalised)",
                                50, -5, 10 },
                              { "tpc_StoppingChi2CosmicRatio_central", "tpc_StoppingChi2CosmicRatio",
                                ";Stopping Fits Chi2 Ratio;Slices (normalised)",
                                50, -5, 10, kBlack, false, "tpc_NuVertexY>180 && tpc_CRLongestTrackDirY>-0.6" },
                              { "pds_FMTotalScore", "pds_FMTotalScore", ";SimpleFlash Total Score;Slices (normalised)",
                                50, -10, 60 },
                              { "pds_FMPE", "pds_FMPE", ";SimpleFlash Flash PE Count;Slices (normalised)",
                                50, -10000, 40000 },
                              { "pds_FMTime", "pds_FMTime", ";SimpleFlash Flash Time (#mus);Slices (normalised)",
                                50, -500, 200 },
                              { "crt_TrackScore", "crt_TrackScore", ";CRT Track Match Score;Slices (normalised)",
                                50, -4, 200 },
                              { "crt_TrackScore_log", "crt_TrackScore", ";CRT Track Match Score;Slices (normalised)",
                                50, -4, 200, kBlack, true },
                              { "crt_SPScore", "crt_SPScore", ";CRT SpacePoint Match Score (cm);Slices (normalised)",
                                50, -4, 100 },
                              { "crt_SPScore_log", "crt_SPScore", ";CRT SpacePoint Match Score (cm);Slices (normalised)",
                                50, -4, 100, kBlack, true },
                              { "crt_TrackTime", "crt_TrackTime",  ";CRT Track Match Time (#mus);Slices (normalised)",
                                50, -3000, 1500 },
                              { "crt_TrackTime_log", "crt_TrackTime",  ";CRT Track Match Time (#mus);Slices (normalised)",
                                50, -3000, 1500, kBlack, true },
                              { "crt_SPTime", "crt_SPTime", ";CRT SpacePoint Match Time (#mus);Slices (normalised)",
                                50, -3000, 1500 },
                              { "crt_SPTime_log", "crt_SPTime", ";CRT SpacePoint Match Time (#mus);Slices (normalised)",
                                50, -3000, 1500, kBlack, true },
                              { "pds_OpT0Score", "pds_OpT0Score", ";OpT0 Score;Slices (normalised)",
                                50, -5000, 40000 },
                              { "pds_OpT0MeasuredPE", "pds_OpT0MeasuredPE",
                                ";OpT0 Measured PE;Slices (normalised)",
                                50, -10000, 100000 },
  };

  std::vector<Cut> cuts = { { "signal", "strstr(matchedType,\"Nu\") && !strstr(matchedType,\"DirtNu\") && matchedPurity > 0.8 && matchedCompleteness > 0.8", "True Neutrino Slices", kBlue+2 },
                            { "background", "!strstr(matchedType,\"Nu\")", "Cosmic Slices", kRed+2 }
  };

  for(auto const& plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);

      if(plot.name == "tpc_CRNHitsMax")
	{
	  gStyle->SetTitleSize(0.04, "x");
	  gStyle->SetTitleOffset(1.6, "x");
	}
      else
	{
	  gStyle->SetTitleSize(0.055, "x");
	  gStyle->SetTitleOffset(1.1, "x");
	}
      
      MakeComparisonPlot(canvas, slices, plot, cuts, {.25, .82, .8, .87});

      canvas->SaveAs(save_dir + "/" + plot.name + ".png");
      canvas->SaveAs(save_dir + "/" + plot.name + ".pdf");
    }
}
