void ChannelMapping()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/channel_mapping";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run2100/ana/data*_ana.root");

  const TString run_name = "run2100";
  gSystem->Exec("mkdir -p " + saveDir + "/" + run_name);

  std::vector<double> *chit_x = 0, *chit_y = 0, *chit_z = 0;
  std::vector<std::vector<uint16_t> > *chit_sipm_feb_mac5 = 0;

  tree->SetBranchAddress("chit_x", &chit_x);
  tree->SetBranchAddress("chit_y", &chit_y);
  tree->SetBranchAddress("chit_z", &chit_z);
  tree->SetBranchAddress("chit_sipm_feb_mac5", &chit_sipm_feb_mac5);

  const unsigned N = tree->GetEntries();

  TH2F *hists_xy[8];
  TH2F *hists_yz[8];
  
  for(unsigned int i = 0; i < 8; ++i)
    {
      hists_xy[i] = new TH2F(Form("hXY%d",i), ";CRTHit x position (cm);CRTHit y position (cm);CRTHits", 30, -220, 80, 30, -150, 150);
      hists_yz[i] = new TH2F(Form("hYZ%d",i), ";CRTHit y position (cm);CRTHit z position (cm);CRTHits", 30, -150, 150, 120, -300, 900);
    }

  for(unsigned ev = 0; ev < N; ++ev)
    {
      tree->GetEntry(ev);

      for(unsigned int ht = 0; ht < chit_x->size(); ++ht)
	{
	  const double x = chit_x->at(ht);
	  const double y = chit_y->at(ht);
	  const double z = chit_z->at(ht);

	  hists_xy[chit_sipm_feb_mac5->at(ht)[0]]->Fill(x,y);
	  hists_xy[chit_sipm_feb_mac5->at(ht)[1]]->Fill(x,y);
	  hists_yz[chit_sipm_feb_mac5->at(ht)[0]]->Fill(y,z);
	  hists_yz[chit_sipm_feb_mac5->at(ht)[1]]->Fill(y,z);
	}
    }

  gStyle->SetNdivisions(505, "x");
  gStyle->SetPalette(kBlueRedYellow);
  
  for(unsigned int i = 0; i < 8; ++i)
    {
      TCanvas *canvas_xy = new TCanvas(Form("c_crthit_xy_%d", i), Form("c_crthit_xy_%d", i));
      canvas_xy->cd();
      canvas_xy->SetRightMargin(0.25);
      hists_xy[i]->Draw("colz");
      hists_xy[i]->GetYaxis()->SetTitleOffset(1.25);
      hists_xy[i]->GetZaxis()->SetTitleOffset(1.3);
	  
      if(save)
	{
	  canvas_xy->SaveAs(saveDir + "/" + run_name + "/crthit_xy" + Form("_board_%d", i) + ".png");
	  canvas_xy->SaveAs(saveDir + "/" + run_name + "/crthit_xy" + Form("_board_%d", i) + ".pdf");
	}

      TCanvas *canvas_yz = new TCanvas(Form("c_crthit_yz_%d", i), Form("c_crthit_yz_%d", i));
      canvas_yz->cd();
      canvas_yz->SetRightMargin(0.25);
      hists_yz[i]->Draw("colz");
      hists_yz[i]->GetYaxis()->SetTitleOffset(1.25);
      hists_yz[i]->GetZaxis()->SetTitleOffset(1.3);
	  
      if(save)
	{
	  canvas_yz->SaveAs(saveDir + "/" + run_name + "/crthit_yz" + Form("_board_%d", i) + ".png");
	  canvas_yz->SaveAs(saveDir + "/" + run_name + "/crthit_yz" + Form("_board_%d", i) + ".pdf");
	}
    }
}

      
      
    
