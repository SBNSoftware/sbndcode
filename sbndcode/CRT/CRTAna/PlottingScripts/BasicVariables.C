void BasicVariables()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/v09_66_02/basicvariables";
  gSystem->Exec("mkdir -p " + saveDir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_66_02.root");

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
    bool logy = false;
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
                             {"sh_saturated1","sh_saturated1", ";SiPM 1;Strip Hits",
                              2, 0, 2, plotcolour, {"Not Saturated", "Saturated"}, "", true},
                             {"sh_saturated2","sh_saturated2", ";SiPM 2;Strip Hits",
                              2, 0, 2, plotcolour, {"Not Saturated", "Saturated"}, "", true},
                             {"sh_saturated","sh_saturated1 || sh_saturated2", ";Either SiPM;Strip Hits",
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
                             {"cl_composition","cl_composition", ";Composition;Clusters",
                              8, 0, 8, plotcolour, {"Undefined","X", "Y", "XY", "Z", "XZ", "YZ", "XYZ"}, "", true},
                             {"cl_composition_log","cl_composition", ";Composition;Clusters",
                              8, 0, 8, plotcolour, {"Undefined","X", "Y", "XY", "Z", "XZ", "YZ", "XYZ"}, "", true, true},
                             {"cl_completeness","cl_truth_completeness", ";Completeness;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
                             {"cl_purity","cl_truth_purity", ";Purity;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
                             {"cl_completeness_multi_strip_clusters","cl_truth_completeness", ";Completeness;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour, {}, "cl_nhits>2"},
                             {"cl_purity_multi_strip_clusters","cl_truth_purity", ";Purity;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour, {}, "cl_nhits>2"},
                             {"sp_x","cl_sp_x", ";x [cm];Space Points",
                              100, -500, 500, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_ex","cl_sp_ex", ";x error [cm];Space Points",
                              50, 0, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_y","cl_sp_y", ";y [cm];Space Points",
                              100, -500, 1500, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_ey","cl_sp_ey", ";y error [cm];Space Points",
                              50, 0, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_z","cl_sp_z", ";z [cm];Space Points",
                              100, -250, 850, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_ez","cl_sp_ez", ";z error [cm];Space Points",
                              50, 0, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_pe","cl_sp_pe", ";PE;Space Points",
                              100, 0, 800, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_time","cl_sp_time", ";Time in readout [ns];Space Points",
                              100, 0, 3.5e6, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_etime","cl_sp_etime", ";Time error [ns];Space Points",
                              50, 0, 100, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"cl_has_sp","cl_has_sp", ";Has SpacePoint?;Clusters",
                              2, 0, 2, plotcolour, {"No","Yes"}, "", true},
                             {"cl_sp_complete","cl_sp_complete", ";Is SpacePoint Complete?;Clusters",
                              2, 0, 2, plotcolour, {"No","Yes"}, "cl_has_sp"},
                             {"tr_start_x","tr_start_x", ";x start [cm];Tracks",
                              100, -500, 500, plotcolour},
                             {"tr_start_y","tr_start_y", ";y start [cm];Tracks",
                              100, -500, 1500, plotcolour},
                             {"tr_start_z","tr_start_z", ";z start [cm];Tracks",
                              100, -250, 850, plotcolour},
                             {"tr_end_x","tr_end_x", ";x end [cm];Tracks",
                              100, -500, 500, plotcolour},
                             {"tr_end_y","tr_end_y", ";y end [cm];Tracks",
                              100, -500, 1500, plotcolour},
                             {"tr_end_z","tr_end_z", ";z end [cm];Tracks",
                              100, -250, 850, plotcolour},
                             {"tr_dir_x","tr_dir_x", ";x direction [cm];Tracks",
                              50, -1, 1, plotcolour},
                             {"tr_dir_y","tr_dir_y", ";y direction [cm];Tracks",
                              50, -1, 1, plotcolour},
                             {"tr_dir_z","tr_dir_z", ";z direction [cm];Tracks",
                              50, -1, 1, plotcolour},
                             {"tr_time","tr_time", ";Time in readout [ns];Tracks",
                              100, 0, 3.5e6, plotcolour},
                             {"tr_pe","tr_pe", ";PE;Tracks",
                              100, 0, 2000, plotcolour},
                             {"tr_length", "tr_length", ";Length [cm];Tracks",
                              100, 0, 1500, plotcolour},
                             {"tr_tof", "tr_tof", ";ToF [ns];Tracks",
                              100, 0, 100, plotcolour},
                             {"tr_theta", "tr_theta", ";#theta [#circ];Tracks",
                              90, 0, 180, plotcolour},
                             {"tr_phi", "tr_phi", ";#phi [#circ];Tracks",
                              90, -180, 180, plotcolour},
                             {"tr_triple", "tr_triple", ";Has 3 SP?;Tracks",
                              2, 0, 2, plotcolour, {"No", "Yes"}, "", true},
                             {"tr_completeness","tr_truth_completeness", ";Completeness;Tracks",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
                             {"tr_purity","tr_truth_purity", ";Purity;Tracks",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour},
  };

  const auto copy_of_plots = plots;

  for(auto const &plot : copy_of_plots)
    {
      if(plot.name.BeginsWith("tr_"))
        {
          plt a = plot;
          plt b = plot;

          a.req = "!tr_triple";
          a.name.Append("_two_sp");

          b.req = "tr_triple";
          b.name.Append("_three_sp");

          plots.push_back(a);
          plots.push_back(b);
        }
    }
  
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

      if(plot.logy)
        {
          hist->SetMinimum(1);
          canvas->SetLogy();
        }

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
