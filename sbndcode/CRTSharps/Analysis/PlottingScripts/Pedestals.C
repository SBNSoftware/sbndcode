
void Pedestals()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/pedestals";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4460/ana/data*_ana.root");

  const TString run_name = "run4460";
  gSystem->Exec("mkdir -p /sbnd/data/users/hlay/crt/sharps/plots/pedestals/" + run_name);

  std::vector<std::vector<uint16_t>> *feb_adc = 0;
  std::vector<uint16_t>              *feb_mac5 = 0;
  std::vector<uint16_t>              *feb_flags = 0;

  tree->SetBranchAddress("feb_adc", &feb_adc);
  tree->SetBranchAddress("feb_mac5", &feb_mac5);
  tree->SetBranchAddress("feb_flags", &feb_flags);

  const int plotcolour = kBlue-4;

  std::vector<std::vector<TH1F*>> hists;
  hists.resize(8);
  
  for(unsigned geoID = 0; geoID < 8; ++geoID)
    {
      for(unsigned channel = 0; channel < 32; ++channel)
	{
	  TString boardname = std::to_string(geoID) + "_" + std::to_string(channel);
	  hists[geoID].emplace_back(new TH1F("hPedestal" + boardname, ";ADC;T1 Reset Events", 
					     100, 0, 500));
	}
    }

  const unsigned N = tree->GetEntries();
  std::cout << "==== TOTAL ENTRIES: " << N << std::endl;

  for(unsigned i = 0; i < N; ++i)
    {
      if(i%1000 == 0) std::cout << i << " / " << N << std::endl;
      tree->GetEntry(i);
      for(unsigned ii = 0; ii < feb_flags->size(); ++ii)
	{
	  if(feb_flags->at(ii) != 11)
	    continue;

	  for(unsigned ch = 0; ch < 32; ++ch)
	    hists[feb_mac5->at(ii)][ch]->Fill((feb_adc->at(ii)).at(ch));
	}
    }

  std::stringstream ss;

  for(unsigned geoID = 0; geoID < 8; ++geoID)
    {
      for(unsigned channel = 0; channel < 32; ++channel)
	{
	  TString boardname = std::to_string(geoID) + "_" + std::to_string(channel);
	  TCanvas *canvas = new TCanvas("cPedestal" + boardname,
					"cPedestal" + boardname);
	  canvas->cd();
	  hists[geoID][channel]->SetLineColor(plotcolour);
	  hists[geoID][channel]->GetYaxis()->SetTitleOffset(1.3);
	  hists[geoID][channel]->Draw("histE");
	  TF1 *gausFunc = new TF1("gausFunc", "gaus",0,500);
	  TFitResultPtr fit = hists[geoID][channel]->Fit("gausFunc","S");
	  gausFunc->SetLineColor(kRed+2);
	  gausFunc->SetLineWidth(5);
	  gausFunc->Draw("same");

	  TPaveText *pt = new TPaveText(.7,.78,.9,.9, "NDC");
	  pt->AddText(Form("GeoID: %i, Channel: %i", geoID, channel));
	  pt->AddText(Form("Mean:  %.3f", gausFunc->GetParameter("Mean")));
	  pt->AddText(Form("Sigma: %.3f", gausFunc->GetParameter("Sigma")));
	  pt->SetFillColor(kWhite);
	  pt->SetLineColor(kBlack);
	  pt->SetLineWidth(4);
	  pt->SetTextColor(kRed+2);
	  TText *t1 = pt->GetLineWith("GeoID");
	  t1->SetTextColor(kBlack);
	  pt->SetTextAlign(12);
	  pt->Draw("same");

	  if(save)
	    {
	      canvas->SaveAs(saveDir + "/" + run_name + "/pedestal" + boardname + "_fit.png");
	      canvas->SaveAs(saveDir + "/" + run_name + "/pedestal" + boardname + "_fit.pdf");
	    }
	  
	  delete canvas;
	  
	  ss << "[" << geoID * 32 + channel << ", " << gausFunc->GetParameter("Mean") 
	     << "]," << std::endl;
	}
    }
  std::cout << ss.str() << std::endl;
}

      
      
    
