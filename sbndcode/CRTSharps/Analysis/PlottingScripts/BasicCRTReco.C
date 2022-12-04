void BasicCRTReco()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/basiccrtreco";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");

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
  
  std::vector<datacut> cuts = { { "", "all"},
				{ "chit_t1 < 3.32e5", "prebeam"},
				{ "chit_t1 > 3.32e5 && chit_t1 < 3.35e5", "beam"},
  };

  std::vector<plt> plots = { {"n_crthits", "@chit_x.size()", ";nCRTHits;Events",
			      50, 0, 50, plotcolour},
			     {"n_crttracks", "@ct_time.size()", ";nCRTTracks;Events",
			      5, 0, 5, plotcolour},
			     {"n_febdatas", "@feb_mac5.size()", ";nFEBDatas;Events",
			      80, 0, 80, plotcolour},
			     {"crthit_x", "-chit_x", ";CRTHit -x position (cm);CRTHits",
			      50, -50, 200, plotcolour},
			     {"crthit_y", "chit_y", ";CRTHit y position (cm);CRTHits",
			      50, -120, 130, plotcolour},
			     {"crthit_z", "chit_z", ";CRTHit z position (cm);CRTHits",
			      120, -300, 900, plotcolour},
			     {"crthit_t0", "1e-6 * chit_t0", ";CRTHit t0 (ms);CRTHits",
			      120, 0, 1200, plotcolour, {}, "", true},
			     {"crthit_t1", "1e-6 * chit_t1", ";CRTHit t1 (ms);CRTHits",
			      100, 0, 1200, plotcolour, {}, "", true},
			     {"crthit_t1_wide_beam", "1e-3 * chit_t1", ";CRTHit t1 (#mus);CRTHits",
			      100, 0, 500, plotcolour, {}, "", true},
			     {"crthit_t1_beam", "1e-3 * chit_t1", ";CRTHit t1 (#mus);CRTHits",
			      120, 320, 350, plotcolour, {}, "", true},
			     {"crthit_t1_fine_beam", "1e-3 * chit_t1", ";CRTHit t1 (#mus);CRTHits",
			      120, 332, 336, plotcolour, {}, "", true},
			     {"crthit_t1_diff", "chit_t1_diff", ";CRTHit #Delta t1 (ns);CRTHits",
			      100, -50, 50, plotcolour, {}, "", true},
    			     {"crthit_pes", "chit_pes", ";CRTHit PEs;CRTHits",
			      40, 0, 800, plotcolour},
			     {"crthit_panel", "chit_z > 0", ";CRTHit panel;CRTHits",
			      2, 0, 2, plotcolour, {"Upstream", "Downstream"}, "", true},
			     {"crthit_sipm_raw_adc", "chit_sipm_raw_adc", ";CRTHit Raw ADC;CRTHits",
			      100, 0, 5000, plotcolour},
			     {"crthit_sipm_adc", "chit_sipm_adc", ";CRTHit ADC;CRTHits",
			      100, 0, 5000, plotcolour},
			     {"crthit_sipm_corr_adc", "chit_sipm_corr_adc", ";CRTHit Corrected ADC;CRTHits",
			      160, 0, 8000, plotcolour},
			     {"feb_geo_id", "feb_mac5", ";FEB Geo ID;FEBDatas",
			      9, 0, 9, plotcolour},
			     {"feb_t0", "1e-6 * feb_ts0", ";FEB t0 (ms);FEBDatas",
			      40, 0, 1200, plotcolour},
			     {"feb_t1", "1e-6 * feb_ts1", ";FEB t1 (ms);FEBDatas",
			      40, 0, 1200, plotcolour},
			     {"feb_t1_beam", "1e-3 * feb_ts1", ";FEB t1 (#mus);FEBDatas",
			      40, 0, 1000, plotcolour},
			     {"feb_adc", "feb_adc", ";SiPM ADC;SiPMs",
			      50, 0, 5000, plotcolour},
			     {"track_length","ct_length", ";Track Length (cm);Tracks",
			      50, 900, 1100, plotcolour},
			     {"track_hit1_x","-ct_hit1_x", ";Track hit 1 -x position (cm);Tracks",
			      50, -50, 200, plotcolour},
			     {"track_hit1_y","ct_hit1_y", ";Track hit 1 y position (cm);Tracks",
			      50, -120, 130, plotcolour},
			     {"track_hit1_z","ct_hit1_z", ";Track hit 1 z position (cm);Tracks",
			      120, -300, 900, plotcolour},
			     {"track_hit2_x","-ct_hit2_x", ";Track hit 2 -x position (cm);Tracks",
			      50, -50, 200, plotcolour},
			     {"track_hit2_y","ct_hit2_y", ";Track hit 2 y position (cm);Tracks",
			      50, -120, 130, plotcolour},
			     {"track_hit2_z","ct_hit2_z", ";Track hit 2 z position (cm);Tracks",
			      120, -300, 900, plotcolour},
			     {"track_hit1_sipm_raw_adc", "ct_hit1_sipm_raw_adc", ";Track hit 1 Raw ADC;Tracks",
			      100, 0, 5000, plotcolour},
			     {"track_hit1_sipm_adc", "ct_hit1_sipm_adc", ";Track hit 1 ADC;Tracks",
			      100, 0, 5000, plotcolour},
			     {"track_hit1_sipm_corr_adc", "ct_hit1_sipm_corr_adc", ";Track hit 1 Corrected ADC;Tracks",
			      160, 0, 8000, plotcolour},
			     {"track_hit2_sipm_raw_adc", "ct_hit2_sipm_raw_adc", ";Track hit 2 Raw ADC;Tracks",
			      100, 0, 5000, plotcolour},
			     {"track_hit2_sipm_adc", "ct_hit2_sipm_adc", ";Track hit 2 ADC;Tracks",
			      100, 0, 5000, plotcolour},
			     {"track_hit2_sipm_corr_adc", "ct_hit2_sipm_corr_adc", ";Track hit 2 Corrected ADC;Tracks",
			      160, 0, 8000, plotcolour},
  };

  std::vector<plttwod> twodplots = { 
    {"crthit_xy", "chit_y:-chit_x", 
     ";CRTHit -x position (cm);CRTHit y position (cm);CRTHits",
     50, -50, 200, 50, -120, -130},
    {"crthit_xy_upstream", "chit_y:-chit_x", 
     ";CRTHit -x position (cm);CRTHit y position (cm);CRTHits",
     50, -50, 200, 50, -120, -130, "chit_z < 0"},
    {"crthit_xy_downstream", "chit_y:-chit_x", 
     ";CRTHit -x position (cm);CRTHit y position (cm);CRTHits",
     50, -50, 200, 50, -120, -130, "chit_z > 0"},
    {"crthit_xz", "chit_z:-chit_x", 
     ";CRTHit -x position (cm);CRTHit z position (cm);CRTHits",
     50, -50, 200, 120, -300, 900},
    {"crthit_zy", "chit_y:chit_z", 
     ";CRTHit z position (cm);CRTHit y position (cm);CRTHits",
     120, -300, 900, 50, -120, 130},
  };

  for(auto const &cut : cuts)
    {
      gSystem->Exec("mkdir -p " + saveDir + "/" + run_name + "/" + cut.name);

      for(auto const &plot : plots)
	{
	  if(cut.name != "all" && !plot.name.Contains("crthit"))
	    continue;
	  
	  TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name, 
					"c_" + plot.name + "_" + cut.name);
	  canvas->cd();
	  TH1F* hist = new TH1F(plot.name + "_" + cut.name, plot.axes_labels, 
				plot.nbins, plot.xlow, plot.xhigh);
	  hist->SetLineColor(plot.colour);
	  hist->GetYaxis()->SetTitleOffset(1.3);
	  hist->GetYaxis()->SetNdivisions(507);
	  if(plot.yaxisfromzero) hist->SetMinimum(0);
	  tree->Draw(plot.var + ">>" + plot.name + "_" + cut.name, cut.cut + plot.req,"histE");

	  if(plot.binlabels.size() == plot.nbins)
	    {
	      for(int bin = 1; bin <= hist->GetXaxis()->GetNbins(); ++bin)
		hist->GetXaxis()->SetBinLabel(bin, plot.binlabels[bin-1]);
	    }

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + run_name + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
	      canvas->SaveAs(saveDir + "/" + run_name + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");
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
	  TH2F* hist = new TH2F(plot.name + "_" + cut.name, plot.axes_labels, plot.nbinsx, plot.xlow, plot.xhigh,
				plot.nbinsy, plot.ylow, plot.yhigh);
	  tree->Draw(plot.var + ">>" + plot.name + "_" + cut.name, cut.cut + plot.req, "colz");
	  hist->GetYaxis()->SetTitleOffset(1.25);
	  hist->GetZaxis()->SetTitleOffset(1.3);

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + run_name + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
	      canvas->SaveAs(saveDir + "/" + run_name + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");
	    }
	  delete canvas, hist;
	}
    }
}
