#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void InputVariablesChi2Only()
{
  const TString saveDir = "/sbnd/data/users/hlay/razzled/plots/investigations/inputvariables_chi2only";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + saveDir);

  const bool all = false, split = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/bnb/razzled*.root");
  tree->Add("/sbnd/data/users/hlay/razzled/intrnue/razzled*.root");

  std::vector<Cut> cuts = { {"all", "", "", kMagenta+2},
                            {"ambiguous", "!unambiguousSlice", "", kGreen+2},
  };

  std::vector<Cut> particles = { {"muon", "abs(truePDG)==13 && !unambiguousSlice", "#mu^{#pm}", kRed+2},
				 {"pion", "abs(truePDG)==211 && !unambiguousSlice", "#pi^{#pm}", kGreen+2},
				 //				 {"kaon", "abs(truePDG)==321 && !unambiguousSlice", "K^{#pm}", kBlue+2},
				 {"proton", "abs(truePDG)==2212 && !unambiguousSlice", "p", kOrange+2},
				 //				 {"other", "abs(truePDG)!=13 && abs(truePDG)!=211 && abs(truePDG)!=321 && abs(truePDG)!=2212 && !unambiguousSlice", "Other", kBlack}
  };

  TCut quality_cut = "trk_length > 5 && trackContained && pfp_trackScore > .5";

  for(auto& particle : particles)
    particle.cut += quality_cut;

  std::vector<Plot> plots = { {"trk_chi2PIDMuon", "trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;PFPs",
                               60, -20, 100},
                              {"trk_chi2PIDMuon_filled", "trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;PFPs",
                               52, -4, 100},
			      {"trk_chi2PIDPion", "trackChi2PIDPion", ";Track #chi^{2} PID Pion;PFPs",
			       60, -20, 100},
                              {"trk_chi2PIDPion_filled", "trackChi2PIDPion", ";Track #chi^{2} PID Pion;PFPs",
                               52, -4, 100},
			      {"trk_chi2PIDKaon", "trackChi2PIDKaon", ";Track #chi^{2} PID Kaon;PFPs",
			       60, -40, 200},
                              {"trk_chi2PIDKaon_filled", "trackChi2PIDKaon", ";Track #chi^{2} PID Kaon;PFPs",
                               52, -8, 200},
                              {"trk_chi2PIDProton", "trk_chi2PIDProton", ";Track #chi^{2} PID Proton;PFPs",
                               84, -20, 400},
                              {"trk_chi2PIDProton_filled", "trk_chi2PIDProton", ";Track #chi^{2} PID Proton;PFPs",
                               81, -5, 400},
  };

  std::vector<TwoDPlot> twodplots = {
    {"trk_chi2PIDMuonPion","trackChi2PIDPion:trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;Track #chi^{2} PID Pion;Tracks"
     , 52, -4, 100, 52, -4, 100},
    {"trk_chi2PIDMuonKaon","trackChi2PIDKaon:trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;Track #chi^{2} PID Kaon;Tracks"
     , 52, -4, 100, 52, -8, 200},
    {"trk_chi2PIDMuonProton","trk_chi2PIDProton:trk_chi2PIDMuon", ";Track #chi^{2} PID Muon;Track #chi^{2} PID Proton;Tracks"
     , 52, -4, 100, 81, -5, 400},
    {"trk_chi2PIDPionKaon","trackChi2PIDKaon:trackChi2PIDPion", ";Track #chi^{2} PID Pion;Track #chi^{2} PID Kaon;Tracks"
     , 52, -4, 100, 52, -8, 200},
    {"trk_chi2PIDPionProton","trk_chi2PIDProton:trackChi2PIDPion", ";Track #chi^{2} PID Pion;Track #chi^{2} PID Proton;Tracks"
     , 52, -4, 100, 81, -5, 400},
    {"trk_chi2PIDKaonProton","trk_chi2PIDProton:trackChi2PIDKaon", ";Track #chi^{2} PID Kaon;Track #chi^{2} PID Proton;Tracks"
     , 52, -8, 200, 81, -5, 400},
  };

  
  for(auto &plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name + "_byType",
				    "c_" + plot.name + "_byType");
      canvas->cd();

      MakeComparisonPlot(canvas, tree, plot, particles, {.35, .77, .8, .87});

      if(save)
	{
	  canvas->SaveAs(saveDir + "/" + plot.name + ".png");
	  canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
	}
      delete canvas;
    }

  for(auto const &particle : particles)
    {
      for(auto const &plot : twodplots)
	{
	  TString plot_name = plot.name + "_" + particle.name;

	  gStyle->SetNdivisions(505, "x");
	  TCanvas *canvas = new TCanvas("c_" + plot_name, "c_" + plot_name);
	  canvas->cd();
	  canvas->SetRightMargin(0.25);
	  gStyle->SetPalette(kBlueRedYellow);
	  TH2F* hist = new TH2F(plot_name, plot.axes_labels, plot.nbinsx, plot.xlow, plot.xhigh,
				plot.nbinsy, plot.ylow, plot.yhigh);
	  tree->Draw(plot.var + ">>" + plot_name, plot.req + particle.cut, "colz");
	  hist->GetYaxis()->SetTitleOffset(1.25);
	  hist->GetZaxis()->SetTitleOffset(1.3);

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + plot_name + ".png");
	      canvas->SaveAs(saveDir + "/" + plot_name + ".pdf");
	    }
	  delete canvas, hist;
	}
    }
}	  
