#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"

void ReconstructionEfficiency()
{
  const TString saveDir = "/sbnd/data/users/hlay/pizero/plots/pizeroreco/reconstructionefficiency";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  gSystem->Exec("mkdir -p " + saveDir);

  TChain *tree = new TChain("pizeroana/pizeros");
  tree->Add("/sbnd/data/users/hlay/pizero/pizeroana_sbnd_test.root");

  const Cut true_cut = {"neutrino_fv", "mct_origin==1 && mc_fv", "", kGreen+2};
  const Cut reco_cut = {"neutrino_fv_reco", "mct_origin==1 && mc_fv && reco_gamma1_npfps>0 && reco_gamma2_npfps>0 && reco_gamma1_bestpfp_pur>.5 && reco_gamma2_bestpfp_pur>.5",
			"", kGreen+2};
  const Cut reco_cut_gamma1 = {"neutrino_fv_reco_gamma1", "mct_origin==1 && mc_fv && reco_gamma1_npfps>0 && reco_gamma1_bestpfp_pur>.5", "", kRed+2};
  const Cut reco_cut_gamma2 = {"neutrino_fv_reco_gamma2", "mct_origin==1 && mc_fv && reco_gamma2_npfps>0 && reco_gamma2_bestpfp_pur>.5", "", kBlue+2};

  std::vector<Plot> plots = { {"pizero_decay_open_angle", "mc_open_angle",
			       ";#theta_{#gamma#gamma} (#circ);#pi^{0}s",
			       45, 0, 180},
  };

  for(Plot plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      TH1F* true_hist = new TH1F(plot.name + "_" + true_cut.name, plot.axes_labels, plot.nbins, plot.xlow, plot.xhigh);
      tree->Draw(plot.var + ">>" + plot.name + "_" + true_cut.name, plot.req + true_cut.cut, "histE");

      TH1F* reco_hist = new TH1F(plot.name + "_" + reco_cut.name, plot.axes_labels, plot.nbins, plot.xlow, plot.xhigh);
      tree->Draw(plot.var + ">>" + plot.name + "_" + reco_cut.name, plot.req + reco_cut.cut, "histE");

      TH1F* reco_hist_gamma1 = new TH1F(plot.name + "_" + reco_cut_gamma1.name, plot.axes_labels, plot.nbins, plot.xlow, plot.xhigh);
      tree->Draw(plot.var + ">>" + plot.name + "_" + reco_cut_gamma1.name, plot.req + reco_cut_gamma1.cut, "histE");

      TH1F* reco_hist_gamma2 = new TH1F(plot.name + "_" + reco_cut_gamma2.name, plot.axes_labels, plot.nbins, plot.xlow, plot.xhigh);
      tree->Draw(plot.var + ">>" + plot.name + "_" + reco_cut_gamma2.name, plot.req + reco_cut_gamma2.cut, "histE");

      MakeComparisonPlotEff(canvas,
			    {true_hist, true_hist, true_hist},
			    {reco_hist, reco_hist_gamma1, reco_hist_gamma2},
			    plot.axes_labels,
			    {"Reco #gamma_{1}#gamma_{2}", "Reco #gamma_{1}", "Reco #gamma_{2}"},
			    {reco_cut.colour, reco_cut_gamma1.colour, reco_cut_gamma2.colour},
			    {.39, .825, .85, .9},
			    2);

      if(save)
        {
          canvas->SaveAs(saveDir + "/" + plot.name + "_eff" + ".png");
          canvas->SaveAs(saveDir + "/" + plot.name + "_eff" + ".pdf");
        }
    }
}
