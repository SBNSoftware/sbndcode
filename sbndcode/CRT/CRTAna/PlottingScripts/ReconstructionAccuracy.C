void ReconstructionAccuracy()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/v09_67_00/reconstructionaccuracy";
  gSystem->Exec("mkdir -p " + saveDir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_67_00.root");

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
    bool zaxislog = false;
  };

  const int plotcolour = kMagenta+2;

  std::vector<plt> plots = { {"sp_x_accuracy", "cl_sp_x - cl_truth_x", ";(reco - true) x [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_x_accuracy_no_sides", "cl_sp_x - cl_truth_x", ";(reco - true) x [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=3 && cl_tagger!=4"},
                             {"sp_x_accuracy_no_sides_log", "cl_sp_x - cl_truth_x", ";(reco - true) x [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=3 && cl_tagger!=4", false, true},
                             {"sp_y_accuracy", "cl_sp_y - cl_truth_y", ";(reco - true) y [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_y_accuracy_no_top_bottom", "cl_sp_y - cl_truth_y", ";(reco - true) y [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=0 && cl_tagger!=5 && cl_tagger!=6"},
                             {"sp_y_accuracy_no_top_bottom_log", "cl_sp_y - cl_truth_y", ";(reco - true) y [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=0 && cl_tagger!=5 && cl_tagger!=6", false, true},
                             {"sp_z_accuracy", "cl_sp_z - cl_truth_z", ";(reco - true) z [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_z_accuracy_no_front_back", "cl_sp_z - cl_truth_z", ";(reco - true) z [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=1 && cl_tagger!=2"},
                             {"sp_z_accuracy_no_front_back_log", "cl_sp_z - cl_truth_z", ";(reco - true) z [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=1 && cl_tagger!=2", false, true},
                             {"sp_pos_accuracy", "sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Space Points",
                              55, -5, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_pos_accuracy_multi_strip_clusters", "sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Space Points",
                              55, -5, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_nhits>2"},
                             {"sp_core_x_accuracy", "cl_sp_x - cl_truth_core_x", ";(reco - true) x [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_core_x_accuracy_no_sides", "cl_sp_x - cl_truth_core_x", ";(reco - true) x [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=3 && cl_tagger!=4"},
                             {"sp_core_x_accuracy_no_sides_log", "cl_sp_x - cl_truth_core_x", ";(reco - true) x [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=3 && cl_tagger!=4", false, true},
                             {"sp_core_y_accuracy", "cl_sp_y - cl_truth_core_y", ";(reco - true) y [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_core_y_accuracy_no_top_bottom", "cl_sp_y - cl_truth_core_y", ";(reco - true) y [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=0 && cl_tagger!=5 && cl_tagger!=6"},
                             {"sp_core_y_accuracy_no_top_bottom_log", "cl_sp_y - cl_truth_core_y", ";(reco - true) y [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=0 && cl_tagger!=5 && cl_tagger!=6", false, true},
                             {"sp_core_z_accuracy", "cl_sp_z - cl_truth_core_z", ";(reco - true) z [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_core_z_accuracy_no_front_back", "cl_sp_z - cl_truth_core_z", ";(reco - true) z [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=1 && cl_tagger!=2"},
                             {"sp_core_z_accuracy_no_front_back_log", "cl_sp_z - cl_truth_core_z", ";(reco - true) z [cm];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_tagger!=1 && cl_tagger!=2", false, true},
                             {"sp_core_pos_accuracy", "sqrt((cl_sp_x - cl_truth_core_x)^2 + (cl_sp_y - cl_truth_core_y)^2 + (cl_sp_z - cl_truth_core_z)^2)", ";(reco - true) position [cm];Space Points",
                              55, -5, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_core_pos_accuracy_multi_strip_clusters", "sqrt((cl_sp_x - cl_truth_core_x)^2 + (cl_sp_y - cl_truth_core_y)^2 + (cl_sp_z - cl_truth_core_z)^2)", ";(reco - true) position [cm];Space Points",
                              55, -5, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete && cl_nhits>2"},
                             {"sp_time_accuracy", "cl_sp_time - (cl_truth_time + 1.7e6)", ";(reco - (true - G4RefTime)) time [ns];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sp_core_time_accuracy", "cl_sp_time - (cl_truth_core_time + 1.7e6)", ";(reco - (true - G4RefTime)) time [ns];Space Points",
                              100, -50, 50, plotcolour, {}, "cl_has_sp && cl_sp_complete"},
                             {"sh_pos_accuracy", "sh_pos - sh_truth_pos", ";(reco - true) position [cm];Strip Hits",
                              100, -50, 50, plotcolour},
                             {"sh_time_accuracy", "sh_ts1 - (sh_truth_time + 1.7e6)", ";(reco - (true - G4RefTime)) time [ns];Strip Hits",
                              100, -50, 50, plotcolour},
                             {"tr_start_x_accuracy", "tr_start_x - tr_truth_start_x", ";(reco - true) Start x [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_start_y_accuracy", "tr_start_y - tr_truth_start_y", ";(reco - true) Start y [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_start_z_accuracy", "tr_start_z - tr_truth_start_z", ";(reco - true) Start z [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_end_x_accuracy", "tr_end_x - tr_truth_end_x", ";(reco - true) End x [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_end_y_accuracy", "tr_end_y - tr_truth_end_y", ";(reco - true) End y [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_end_z_accuracy", "tr_end_z - tr_truth_end_z", ";(reco - true) End z [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_dir_x_accuracy", "tr_dir_x - tr_truth_dir_x", ";(reco - true) Direction x;Tracks",
                              100, -1, 1, plotcolour},
                             {"tr_dir_y_accuracy", "tr_dir_y - tr_truth_dir_y", ";(reco - true) Direction y;Tracks",
                              100, -1, 1, plotcolour},
                             {"tr_dir_z_accuracy", "tr_dir_z - tr_truth_dir_z", ";(reco - true) Direction z;Tracks",
                              100, -1, 1, plotcolour},
                             {"tr_time_accuracy", "tr_time - (tr_truth_time + 1.7e6)", ";(reco - (true - G4RefTime)) time [ns];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_length_accuracy", "tr_length - tr_truth_length", ";(reco - true) Length [cm];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_tof_accuracy", "tr_tof - tr_truth_tof", ";(reco - true) ToF [ns];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_theta_accuracy", "tr_theta - tr_truth_theta", ";(reco - true) #theta [#circ];Tracks",
                              100, -50, 50, plotcolour},
                             {"tr_phi_accuracy", "tr_phi - tr_truth_phi", ";(reco - true) #phi [#circ];Tracks",
                              100, -50, 50, plotcolour},
  };

  std::vector<plttwod> twodplots = {
    {"sp_pe_energy_relation", "cl_truth_energy*1e3:cl_sp_pe", ";PE;True Energy (MeV);SpacePoints",
     50, 0, 700, 50, 0, 100, "cl_has_sp && cl_sp_complete"},
    {"sp_pe_energy_relation_zoom", "cl_truth_energy*1e3:cl_sp_pe", ";PE;True Energy (MeV);SpacePoints",
     50, 0, 700, 50, 0, 20, "cl_has_sp && cl_sp_complete"},
    {"sp_pe_energy_relation_log", "cl_truth_energy*1e3:cl_sp_pe", ";PE;True Energy (MeV);SpacePoints",
     50, 0, 700, 50, 0, 100, "cl_has_sp && cl_sp_complete", true},
    {"sp_pe_energy_relation_zoom_log", "cl_truth_energy*1e3:cl_sp_pe", ";PE;True Energy (MeV);SpacePoints",
     50, 0, 700, 50, 0, 20, "cl_has_sp && cl_sp_complete", true},
    {"sp_pos_accuracy_nhits", "cl_nhits:sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];# Hits in Cluster;SpacePoints",
     55, -5, 50, 10, 0, 10, "cl_has_sp && cl_sp_complete", true},
    {"sp_pos_accuracy_completeness", "cl_truth_completeness:sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Cluster completeness;SpacePoints",
     55, -5, 50, 50, 0, 1 + std::numeric_limits<double>::epsilon(), "cl_has_sp && cl_sp_complete", true},
    {"sp_pos_accuracy_purity", "cl_truth_purity:sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Cluster purity;SpacePoints",
     55, -5, 50, 50, 0, 1 + std::numeric_limits<double>::epsilon(), "cl_has_sp && cl_sp_complete", true},
    {"sp_pos_accuracy_completeness_multi_strip_clusters", "cl_truth_completeness:sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Cluster completeness;SpacePoints",
     55, -5, 50, 50, 0, 1 + std::numeric_limits<double>::epsilon(), "cl_has_sp && cl_sp_complete && cl_nhits>2", true},
    {"sp_pos_accuracy_purity_multi_strip_clusters", "cl_truth_purity:sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Cluster purity;SpacePoints",
     55, -5, 50, 50, 0, 1 + std::numeric_limits<double>::epsilon(), "cl_has_sp && cl_sp_complete && cl_nhits>2", true},
    {"tr_pe_energy_relation", "tr_truth_energy*1e3:tr_pe", ";PE;True Energy (MeV);Tracks",
     50, 0, 1400, 50, 0, 100, "cl_has_sp && cl_sp_complete"},
  };

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
      if(plot.zaxislog) canvas->SetLogz();
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
