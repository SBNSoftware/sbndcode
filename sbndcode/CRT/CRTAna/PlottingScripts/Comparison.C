#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"

#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"

void Comparison()
{
  const TString save_dir = "/sbnd/data/users/hlay/crt/clustering/plots/comparison";
  gSystem->Exec("mkdir -p " + save_dir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree1 = new TChain("crtana/tree");
  TChain *tree2 = new TChain("crtana/tree");
  tree1->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_66_00_updated_clustering.root");
  tree2->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_66_00_fixed_truth_matching.root");

  const int plotcolour = kMagenta+2;
  
  std::vector<Cut> cuts = { { "", ""},
			    { "multi_strip_clusters", "cl_nhits>2"},
			    { "not_multi_strip_clusters", "cl_nhits<3"}
  };

  std::vector<Plot> plots = { {"cl_completeness","cl_truth_completeness", ";Completeness;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour, false},
                             {"cl_purity","cl_truth_purity", ";Purity;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour, false},
			     {"cl_completeness_log","cl_truth_completeness", ";Completeness;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour, true},
                             {"cl_purity_log","cl_truth_purity", ";Purity;Clusters",
                              50, 0, 1 + std::numeric_limits<double>::epsilon(), plotcolour, true},
                             {"sp_pos_accuracy", "sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Space Points",
                              105, -5, 100, plotcolour, false, "cl_has_sp && cl_sp_complete"},
                             {"sp_pos_accuracy_log", "sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Space Points",
                              105, -5, 100, plotcolour, true, "cl_has_sp && cl_sp_complete"},
                             {"sp_pos_accuracy_long_tail", "sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Space Points",
                              505, -5, 500, plotcolour, false, "cl_has_sp && cl_sp_complete"},
                             {"sp_pos_accuracy_long_tail_log", "sqrt((cl_sp_x - cl_truth_x)^2 + (cl_sp_y - cl_truth_y)^2 + (cl_sp_z - cl_truth_z)^2)", ";(reco - true) position [cm];Space Points",
                              505, -5, 500, plotcolour, true, "cl_has_sp && cl_sp_complete"},
  };
  
  for(auto const &cut : cuts)
    {
      for(auto const &plot : plots)
	{
	  TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
	  canvas->cd();

	  MakeComparisonPlot(canvas, tree1, tree2, plot, "Improved Clustering", "Fixed Truth Matching", cut);

	  if(save)
	    {
	      TString plot_name = plot.name;
	      if(cut.name != "")
		plot_name.Append("_" + cut.name);

	      canvas->SaveAs(save_dir + "/" + plot_name + "_compare_fixed_truth_matching.png");
	      canvas->SaveAs(save_dir + "/" + plot_name + "_compare_fixed_truth_matching.pdf");
	    }

	  delete canvas;
	}
    }
}
