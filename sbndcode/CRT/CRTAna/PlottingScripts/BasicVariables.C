void BasicVariables()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/basicvariables";
  gSystem->Exec("mkdir -p " + saveDir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt/crt_clustering_bnb_cosmics/crtana_sbnd.root");

  struct datacut {
    TCut cut;
    TString name;
  };

  struct plt {
    TString name;
    TString var;
    TString axes_labels;
    int nbins;
    double xlow;
    double xhigh;
    int colour;
    std::vector<TString> binlabels = {};
    TCut req = "";
    bool yaxisfromzero = false;
  };

  struct plttwod {
    TString name;
    TString var;
    TString axes_labels;
    int nbinsx;
    double xlow;
    double xhigh;
    int nbinsy;
    double ylow;
    double yhigh;
    TCut req = "";
  };

  const int plotcolour = kMagenta+2;
  
  std::vector<plt> plots = { {"n_strip_hits", "@sh_channel.size()", ";nStripHits;Events",
                              50, 0, 500, plotcolour},
                             {"n_clusters", "@cl_ts0.size()", ";nClusters;Events",
                              50, 0, 500, plotcolour},
                             {"sh_channel","sh_channel", ";Channel ID;Strip Hits",
                              10000, 0, 10000, plotcolour},
                             {"sh_ts0","sh_ts0", ";Ts0 [ns];Strip Hits",
                              120, 0, 1.2e9, plotcolour},
                             {"sh_ts1","sh_ts1", ";Ts1 [ns];Strip Hits",
                              100, 0, 1e7, plotcolour},
                             {"sh_pos","sh_pos", ";Lateral Position [cm];Strip Hits",
                              30, 0, 15, plotcolour},
                             {"sh_err","sh_err", ";Lateral Position Error [cm];Strip Hits",
                              30, 0, 15, plotcolour},
                             {"sh_adc1","sh_adc1", ";ADC1;Strip Hits",
                              50, 0, 5000, plotcolour},
                             {"sh_adc2","sh_adc2", ";ADC2;Strip Hits",
                              50, 0, 5000, plotcolour},
                             {"sh_saturated","sh_saturated", ";;Strip Hits",
                              2, 0, 2, plotcolour, {"Not Saturated", "Saturated"}, "", true},
                             {"sh_completeness","sh_truth_completeness", ";Completeness;Strip Hits",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
                             {"sh_purity","sh_truth_purity", ";Purity;Strip Hits",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
                             {"cl_ts0","cl_ts0", ";Ts0 [ns];Clusters",
                              120, 0, 1.2e9, plotcolour},
                             {"cl_ts1","cl_ts1", ";Ts1 [ns];Clusters",
                              100, 0, 1e7, plotcolour},
                             {"cl_nhits","cl_nhits", ";Hits per cluster;Clusters",
                              10, 0, 10, plotcolour},
                             {"cl_threed","cl_threed", ";;Clusters",
                              2, 0, 2, plotcolour, {"TwoD","ThreeD"}, "", true},
                             {"cl_completeness","cl_truth_completeness", ";Completeness;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
                             {"cl_purity","cl_truth_purity", ";Purity;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
			     {"sp_x","cl_sp_x", ";x [cm];Space Points",
                              100, -500, 500, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_y","cl_sp_y", ";y [cm];Space Points",
                              100, -500, 1500, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_z","cl_sp_z", ";z [cm];Space Points",
                              100, -250, 850, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_pe","cl_sp_pe", ";PE;Space Points",
                              100, 0, 800, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_time","cl_sp_time", ";Time in readout [ns];Space Points",
                              100, 0, 3.5e6, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"cl_has_sp","cl_has_sp", ";Has SpacePoint?;Clusters",
                              2, 0, 2, plotcolour, {"No","Yes"}, "", true},
                             {"cl_sp_complete","cl_sp_complete", ";Is SpacePoint Complete?;Clusters",
                              2, 0, 2, plotcolour, {"No","Yes"}, "cl_has_sp"},
  };
  
  std::vector<plttwod> twodplots = {
    {"cl_completeness_nhits","cl_nhits:cl_truth_completeness", ";Completeness;Hits per cluster;Clusters",
     50, 0, 1 + std::numeric_limits<double>::epsilon(), 10, 0, 10},
    {"cl_purity_nhits","cl_nhits:cl_truth_purity", ";Purity;Hits per cluster;Clusters",
     50, 0, 1 + std::numeric_limits<double>::epsilon(), 10, 0, 10},
  };

  for(auto const &plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();
      TH1F* hist = new TH1F(plot.name, plot.axes_labels,
                            plot.nbins, plot.xlow, plot.xhigh);
      hist->SetLineColor(plot.colour);
      hist->GetYaxis()->SetTitleOffset(1.3);
      hist->GetYaxis()->SetNdivisions(507);
      if(plot.yaxisfromzero) hist->SetMinimum(0);
      tree->Draw(plot.var + ">>" + plot.name, plot.req,"histE");

      if(plot.binlabels.size() == plot.nbins)
        {
          for(int bin = 1; bin <= hist->GetXaxis()->GetNbins(); ++bin)
            hist->GetXaxis()->SetBinLabel(bin, plot.binlabels[bin-1]);
        }

      if(save)
        {
          canvas->SaveAs(saveDir + "/" + plot.name + ".png");
          canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
        }
      delete canvas, hist;
    }

  for(auto const &plot : twodplots)
    {
      gStyle->SetNdivisions(505, "x");
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();
      canvas->SetRightMargin(0.25);
      gStyle->SetPalette(kBlueRedYellow);
      TH2F* hist = new TH2F(plot.name, plot.axes_labels, plot.nbinsx, plot.xlow, plot.xhigh,
                            plot.nbinsy, plot.ylow, plot.yhigh);
      tree->Draw(plot.var + ">>" + plot.name, plot.req, "colz");
      hist->GetYaxis()->SetTitleOffset(1.25);
      hist->GetZaxis()->SetTitleOffset(1.3);

      if(save)
        {
          canvas->SaveAs(saveDir + "/" + plot.name + ".png");
          canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
        }
      delete canvas, hist;
    }
}
