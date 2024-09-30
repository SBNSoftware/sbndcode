#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void RecoVariables()
{
  const TString saveDir = "/sbnd/data/users/hlay/pizero/plots/pizeroreco/recovariables";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  gSystem->Exec("mkdir -p " + saveDir);

  TChain *tree = new TChain("pizeroana/pizeros");
  tree->Add("/sbnd/data/users/hlay/pizero/pizeroana_sbnd_test.root");

  const Cut cut = {"neutrino_fv", "mct_origin==1 && mc_fv", "", kGreen+2};

  const double one_plus_epsilon = 1 + std::numeric_limits<double>::epsilon();

  std::vector<Plot> plots = { {"reco_gamma1_nhits", "reco_gamma1_nhits", ";Hits;#pi^{0}s",
			       50, 0, 1000},
			      {"reco_gamma1_npfps", "reco_gamma1_npfps", ";PFPs;#pi^{0}s",
			       5, -.5, 4.5},
			      {"reco_gamma1_bestpfp_pfp_track_score", "reco_gamma1_bestpfp_pfp_track_score",
			       ";Track Score;#pi^{0}s",
			       50, -one_plus_epsilon, one_plus_epsilon, kBlack, false, "reco_gamma1_npfps!=0"},
			      {"reco_gamma1_bestpfp_comp", "reco_gamma1_bestpfp_comp",
			       ";Completeness;#pi^{0}s",
			       50, 0, one_plus_epsilon, kBlack, false, "reco_gamma1_npfps!=0"},
			      {"reco_gamma1_bestpfp_pur", "reco_gamma1_bestpfp_pur",
			       ";Purity;#pi^{0}s",
			       50, 0, one_plus_epsilon, kBlack, false, "reco_gamma1_npfps!=0"},
			      {"reco_gamma2_nhits", "reco_gamma2_nhits", ";Hits;#pi^{0}s",
			       50, 0, 1000},
			      {"reco_gamma2_npfps", "reco_gamma2_npfps", ";PFPs;#pi^{0}s",
			       5, -.5, 4.5},
			      {"reco_gamma2_bestpfp_pfp_track_score", "reco_gamma2_bestpfp_pfp_track_score",
			       ";Track Score;#pi^{0}s",
			       50, -one_plus_epsilon, one_plus_epsilon, kBlack, false, "reco_gamma2_npfps!=0"},
			      {"reco_gamma2_bestpfp_comp", "reco_gamma2_bestpfp_comp",
			       ";Completeness;#pi^{0}s",
			       50, 0, one_plus_epsilon, kBlack, false, "reco_gamma2_npfps!=0"},
			      {"reco_gamma2_bestpfp_pur", "reco_gamma2_bestpfp_pur",
			       ";Purity;#pi^{0}s",
			       50, 0, one_plus_epsilon, kBlack, false, "reco_gamma2_npfps!=0"},
  };
  
  for(auto &plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      plot.colour = cut.colour;

      MakePlot(canvas, tree, plot, cut);

      if(save)
	{
	  canvas->SaveAs(saveDir + "/" + plot.name + ".png");
	  canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
	}
      delete canvas;
    }
}
