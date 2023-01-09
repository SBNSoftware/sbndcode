void SpacePointVariables()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/spacepointvariables";
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
    bool yaxislog = false;
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

  std::vector<plt> plots = { {"spacepoint_x","cl_sp_x", ";x [cm];SpacePoints",
                              100, -500, 500, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"spacepoint_y","cl_sp_y", ";y [cm];SpacePoints",
                              100, -500, 1500, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"spacepoint_z","cl_sp_z", ";z [cm];SpacePoints",
                              100, -250, 850, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"spacepoint_pe","cl_sp_pe", ";PE;SpacePoints",
                              100, 0, 800, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"spacepoint_time","cl_sp_time", ";Time in readout [ns];SpacePoints",
                              100, 0, 3.5e6, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"cl_has_sp","cl_has_sp", ";Has SpacePoint?;Clusters",
                              2, 0, 2, plotcolour, {"No","Yes"}, "", true},
                             {"cl_sp_complete","cl_sp_complete", ";Is SpacePoint Complete?;Clusters",
                              2, 0, 2, plotcolour, {"No","Yes"}, "cl_has_sp", true},
                             {"cl_x_accuracy", "cl_sp_x - cl_truth_x", ";(reco - true) x [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"cl_x_accuracy_no_sides", "cl_sp_x - cl_truth_x", ";(reco - true) x [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=3 && cl_tagger!=4"},
                             {"cl_x_accuracy_no_sides_log", "cl_sp_x - cl_truth_x", ";(reco - true) x [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=3 && cl_tagger!=4", false, true},
                             {"cl_y_accuracy", "cl_sp_y - cl_truth_y", ";(reco - true) y [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"cl_y_accuracy_no_top_bottom", "cl_sp_y - cl_truth_y", ";(reco - true) y [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=0 && cl_tagger!=5 && cl_tagger!=6"},
                             {"cl_y_accuracy_no_top_bottom_log", "cl_sp_y - cl_truth_y", ";(reco - true) y [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=0 && cl_tagger!=5 && cl_tagger!=6", false, true},
                             {"cl_z_accuracy", "cl_sp_z - cl_truth_z", ";(reco - true) z [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"cl_z_accuracy_no_front_back", "cl_sp_z - cl_truth_z", ";(reco - true) z [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=1 && cl_tagger!=2"},
                             {"cl_z_accuracy_no_front_back_log", "cl_sp_z - cl_truth_z", ";(reco - true) z [cm];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=1 && cl_tagger!=2", false, true},
                             {"cl_time_accuracy", "cl_sp_time - (cl_truth_time + 1.7e6)", ";(reco - (true - G4RefTime)) t [ns];SpacePoints",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
  };

  std::vector<plttwod> twodplots = {
    {"cl_pe_energy_relation", "cl_truth_energy:cl_sp_pe", ";PE;True Energy (MeV);SpacePoints",
     50, 0, 700, 50, 0, 0.1},
    {"cl_pe_energy_relation_zoom", "cl_truth_energy:cl_sp_pe", ";PE;True Energy (MeV);SpacePoints",
     50, 0, 700, 50, 0, 0.02},
  };

  gSystem->Exec("mkdir -p " + saveDir);

  for(auto const &plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();
      if(plot.yaxislog) canvas->SetLogy();
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
