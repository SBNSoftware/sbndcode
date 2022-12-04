void BasicFEBDataAna()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/basicfebdataana";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");

  gSystem->Exec("mkdir -p " + saveDir + "/" + run_name);

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

  const int plotcolour = kBlue-4;
  
  std::vector<datacut> cuts = { { "", "all" },
				{ "feb_flags == 3", "data"},
				{ "feb_flags == 7", "t0reset"},
				{ "feb_flags == 11", "t1reset"}
  };

  std::vector<plt> plots = { {"feb_mac5_id", "feb_mac5", ";FEB Geometry ID;FEBDatas",
			      10, 0, 10, plotcolour},
			     {"feb_flags", "feb_flags", ";FEB Flags;FEBDatas",
			      16, 0, 16, plotcolour},
			     {"feb_t0", "1e-6 * feb_ts0", ";FEB t0 (ms);FEBDatas",
			      40, 0, 1200, plotcolour},
			     {"feb_t1_wide", "1e-6 * feb_ts1", ";FEB t1 (ms);FEBDatas",
			      40, 0, 1200, plotcolour},
			     {"feb_t1_trigger_region", "1e-6 * feb_ts1", ";FEB t1 (ms);FEBDatas",
			      100, 0, 100, plotcolour},
			     {"feb_adc", "feb_adc", ";SiPM ADC;SiPMs",
			      50, 0, 5000, plotcolour},
  };

  for(auto const &cut : cuts)
    {
      for(auto const &plot : plots)
	{
	  TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name, 
					"c_" + plot.name + "_" + cut.name);
	  canvas->cd();
	  TH1F* hist = new TH1F(plot.name + "_" + cut.name, plot.axes_labels, 
				plot.nbins, plot.xlow, plot.xhigh);
	  hist->SetLineColor(plot.colour);
	  hist->GetYaxis()->SetTitleOffset(1.3);
	  if(plot.yaxisfromzero) hist->SetMinimum(0);
	  tree->Draw(plot.var + ">>" + plot.name + "_" + cut.name, cut.cut + plot.req,"histE");

	  if(plot.binlabels.size() == plot.nbins)
	    {
	      for(int bin = 1; bin <= hist->GetXaxis()->GetNbins(); ++bin)
		hist->GetXaxis()->SetBinLabel(bin, plot.binlabels[bin-1]);
	    }

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + run_name + "/" + plot.name + "_" + cut.name + ".png");
	      canvas->SaveAs(saveDir + "/" + run_name + "/" + plot.name + "_" + cut.name + ".pdf");
	    }
	  delete canvas, hist;
	}
    }
}

      
      
    
