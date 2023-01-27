void BasicVariables()
{
  const TString saveDir = "/sbnd/data/users/hlay/pizero/plots/basicvariables";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pizeroana/pizeros");
  tree->Add("/sbnd/data/users/hlay/pizero/pizeroana_sbnd.root");

  struct datacut {
    TCut cut;
    TString name;
    int colour;
  };

  struct plt {
    TString name;
    TString var;
    TString axes_labels;
    int nbins;
    double xlow;
    double xhigh;
    std::vector<TString> binlabels = {};
    TCut req = "";
    bool yaxisfromzero = true;
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

  std::vector<datacut> cuts = { {"", "all", kMagenta+2},
				{"mct_origin==1", "neutrino", kGreen+2},
				{"mct_origin==2", "cosmic", kRed-2},
				{"mct_origin==1 && mc_av", "neutrino_av", kGreen+2},
				{"mct_origin==1 && mc_fv", "neutrino_fv", kGreen+2},
				{"mct_origin==2 && mc_av", "cosmic_av", kRed-2},
				{"mct_origin==2 && mc_fv", "cosmic_fv", kRed-2},
  };

  std::vector<plt> plots = { {"mctruth_origin", "mct_origin", ";Origin;#pi^{0}s",
                              3, 0, 3, {"Unknown", "Neutrino Beam", "Cosmic Ray"}},
			     {"pizero_vx", "mc_vx", ";x position [cm];#pi^{0}s",
			      100, -400, 400},
			     {"pizero_vx_detector", "mc_vx", ";x position [cm];#pi^{0}s",
			      105, -210, 210},
			     {"pizero_vy", "mc_vy", ";y position [cm];#pi^{0}s",
			      100, -400, 400},
			     {"pizero_vy_detector", "mc_vy", ";y position [cm];#pi^{0}s",
			      105, -210, 210},
			     {"pizero_vz", "mc_vz", ";z position [cm];#pi^{0}s",
			      125, -250, 750},
			     {"pizero_vz_detector", "mc_vz", ";z position [cm];#pi^{0}s",
			      130, -10, 510},
			     {"pizero_vt", "mc_vt", ";time [ns];#pi^{0}s",
			      340, -1.8e6, 1.6e6},
			     {"pizero_vt_beam", "mc_vt", ";time [ns];#pi^{0}s",
			      90, -100, 1700},
			     {"pizero_cos_theta", "cos(TMath::DegToRad() * mc_vtheta)",
			      ";cos(#theta);#pi^{0}s",
			      50, 0, 1},
			     {"pizero_cos_phi", "cos(TMath::DegToRad() * mc_vphi)",
			      ";cos(#phi);#pi^{0}s",
			      100, -1, 1},
			     {"pizero_energy", "mc_ve", ";E (GeV);#pi^{0}s",
			      100, 0, 4},
			     {"pizero_decay_open_angle", "mc_open_angle",
			      ";#theta_{#gamma#gamma} (#circ);#pi^{0}s",
			      45, 0, 180},
			     {"pizero_decay_cos_open_angle", "cos(TMath::DegToRad() * mc_open_angle)",
			      ";cos(#theta_{#gamma#gamma});#pi^{0}s",
			      50, 0, 1},
  };
  
  std::vector<plttwod> twodplots = {
    {"pizero_energy_versus_decay_open_angle", "mc_ve:mc_open_angle",
     ";#theta_{#gamma#gamma} (#circ);E (GeV);#pi^{0}s",
     45, 0, 180, 100, 0, 4},
  };
  
  for(auto const & cut : cuts)
    {
      gSystem->Exec("mkdir -p " + saveDir + "/" + cut.name);

      for(auto const &plot : plots)
	{
	  TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name, 
					"c_" + plot.name + "_" + cut.name);
	  canvas->cd();
	  TH1F* hist = new TH1F(plot.name + "_" + cut.name, plot.axes_labels,
				plot.nbins, plot.xlow, plot.xhigh);
	  hist->SetLineColor(cut.colour);
	  hist->GetYaxis()->SetTitleOffset(1.3);
	  hist->GetYaxis()->SetNdivisions(507);
	  if(plot.yaxisfromzero) hist->SetMinimum(0);
	  tree->Draw(plot.var + ">>" + plot.name + "_" + cut.name, plot.req && cut.cut, "histE");

	  if(plot.binlabels.size() == plot.nbins)
	    {
	      for(int bin = 1; bin <= hist->GetXaxis()->GetNbins(); ++bin)
		hist->GetXaxis()->SetBinLabel(bin, plot.binlabels[bin-1]);
	    }

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
	      canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");
	    }
	  delete canvas, hist;
	}

      for(auto const &plot : twodplots)
	{
	  gStyle->SetNdivisions(505, "x");
	  TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name, 
					"c_" + plot.name + "_" + cut.name);
	  canvas->cd();
	  canvas->SetRightMargin(0.25);
	  gStyle->SetPalette(kBlueRedYellow);
	  TH2F* hist = new TH2F(plot.name + "_" + cut.name, plot.axes_labels, 
				plot.nbinsx, plot.xlow, plot.xhigh,
				plot.nbinsy, plot.ylow, plot.yhigh);
	  tree->Draw(plot.var + ">>" + plot.name + "_" + cut.name, plot.req && cut.cut, "colz");
	  hist->GetYaxis()->SetTitleOffset(1.25);
	  hist->GetZaxis()->SetTitleOffset(1.3);

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
	      canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");
	    }
	  delete canvas, hist;
	}
    }
}
