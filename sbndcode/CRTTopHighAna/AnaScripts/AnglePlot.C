#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void AnglePlot()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt_top_high/plots/angleplots";
  const bool save = true;

  if(save)
    gSystem->Exec("mkdir -p " + saveDir);

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt_top_high/production/crttophighana.root");

  const int nbins    = 45;
  const double xlow  = 0;
  const double xhigh = 180;

  std::vector<Plot> plots = { {"angle_to_vertical", "angleToVertical", ";#theta (#circ);#mu^{#pm}s", nbins, xlow, xhigh },
                              {"angle_to_vertical_weighted", "angleToVertical", ";#theta / #Omega (#circ);#mu^{#pm}s (A.U.)", nbins, xlow, xhigh },
                              {"cos_angle_to_vertical", "cos(TMath::DegToRad() * angleToVertical)", ";cos(#theta);#mu^{#pm}s", 40, 0, 1},
                              {"cossq_angle_to_vertical", "pow(cos(TMath::DegToRad() * angleToVertical), 2)", ";cos^{2}(#theta);#mu^{#pm}s", 40, 0, 1}
  };

  
  TH1F *weighting = new TH1F("weighting", "", nbins, xlow, xhigh);
  for(unsigned i = 1; i <= nbins; ++i)
    {
      double costheta0 = TMath::Cos(TMath::DegToRad() * (i - 1) * (xhigh / nbins));
      double costheta1 = TMath::Cos(TMath::DegToRad() * i * (xhigh / nbins));
      double area      = 2 * TMath::Pi() * (costheta0 - costheta1);
      weighting->SetBinContent(i, area);
      weighting->SetBinError(i, 0.);
    }
  
  TCut basecut = "abs(pdg)==13";

  std::vector<Cut> cuts = { {"all", basecut, "All", kMagenta+2},
                            {"enters_exits_tpc", basecut + "entersOrExits", "Enters or Exits TPC", kGreen+2},
                            {"enters_exits_crt", basecut + "entersOrExits && entersOrExitsCRT", "Enters or Exits CRT", kBlue+2},
                            {"has_tpc_track", basecut + "entersOrExits && entersOrExitsCRT && hasTPCTrack", "Has TPC Track", kOrange+2},
                            {"has_tpc_long_track", basecut + "entersOrExits && entersOrExitsCRT && hasTPCLongTrack", "Has Long TPC Track", kRed+2},
  };

  for(auto const& plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      if(plot.name.Contains("weight"))
        MakeComparisonPlot(canvas, tree, plot, cuts, {.23, .77, .83, .87}, false, weighting);
      else
        MakeComparisonPlot(canvas, tree, plot, cuts, {.23, .77, .83, .87}, false);

      if(save)
        {
          canvas->SaveAs(saveDir + "/" + plot.name + ".png");
          canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
        }
      delete canvas;
    }
}
