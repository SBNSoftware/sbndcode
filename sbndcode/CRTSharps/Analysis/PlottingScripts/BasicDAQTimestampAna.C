void BasicDAQTimestampAna()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/basicdaqtimestampana";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";
  gSystem->Exec("mkdir -p " + saveDir + "/" + run_name);

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");
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
  
  std::vector<plt> plots = { {"tdc_channel", "tdc_channel", ";SPEC TDC Channel;TDC Readouts",
			      4, 0, 4, plotcolour},
			     {"tdc_timestamp_pps", "(tdc_timestamp % 1e9)", "Channel 0 PPS;SPEC TDC Time -1s [ns];TDC Readouts",
			      20, 246, 266, plotcolour, {}, "tdc_channel==0"},
			     {"tdc_timestamp_bes", "tdc_timestamp % 1e9", "Channel 1 BES;SPEC TDC Time [ns];TDC Readouts",
			      60, 0, 1.2e9, plotcolour, {}, "tdc_channel==1"},
			     {"tdc_timestamp_rwm", "tdc_timestamp % 1e9", "Channel 2 RWM;SPEC TDC Time [ns];TDC Readouts",
			      60, 0, 1.2e9, plotcolour, {}, "tdc_channel==2"},
			     {"tdc_timestamp_ftrig", "tdc_timestamp % 1e9", "Channel 3 FTRIG;SPEC TDC Time [ns];TDC Readouts",
			      60, 0, 1.2e9, plotcolour, {}, "tdc_channel==3"},
  };

  for(auto const &plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();
      canvas->SetRightMargin(.1);
      canvas->SetTopMargin(.1);
      TH1F* hist = new TH1F(plot.name, plot.axes_labels, 
			    plot.nbins, plot.xlow, plot.xhigh);
      hist->SetLineColor(plot.colour);
      hist->GetYaxis()->SetTitleOffset(1.3);
      if(plot.yaxisfromzero) hist->SetMinimum(0);
      tree->Draw(plot.var + ">>" + plot.name, plot.req,"histE");

      if(plot.binlabels.size() == plot.nbins)
	{
	  for(int bin = 1; bin <= hist->GetXaxis()->GetNbins(); ++bin)
	    hist->GetXaxis()->SetBinLabel(bin, plot.binlabels[bin-1]);
	}

      if(save)
	{
	  canvas->SaveAs(saveDir + "/" + run_name + "/" + plot.name + ".png");
	  canvas->SaveAs(saveDir + "/" + run_name + "/" + plot.name + ".pdf");
	}
      //	  delete canvas, hist;
    }
}
