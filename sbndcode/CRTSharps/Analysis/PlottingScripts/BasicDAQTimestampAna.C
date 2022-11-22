void BasicDAQTimestampAna()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/basicdaqtimestampana";
  const bool save = false;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4460";
  gSystem->Exec("mkdir -p " + saveDir + "/" + run_name);

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/ana/data*_ana.root");
  std::cout << "N Entries: " << tree->GetEntries() << std::endl;

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
				{ "tdc_channel == 0", "PPS"},
				{ "tdc_channel == 1", "BES"},
				{ "tdc_channel == 2", "RWM"},
				{ "tdc_channel == 3", "FTRIG"},
  };

  std::vector<plt> plots = { {"tdc_channel", "tdc_channel", ";SPEC TDC Channel;TDC Readouts",
			      4, 0, 4, plotcolour},
			     {"tdc_timestamp", "tdc_timestamp % 1e9", ";SPEC TDC Time;TDC Readouts",
			      60, 0, 1.2e9, plotcolour},
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

      
      
    
