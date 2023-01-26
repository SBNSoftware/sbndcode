void ReconstructionEfficiency()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/v09_64_01/reconstructionefficiency";
  gSystem->Exec("mkdir -p " + saveDir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt/crt_clustering_bnb_cosmics/crtana_sbnd.root");

  std::vector<double> *td_energy = 0;
  std::vector<bool> *td_reco_status = 0;
  std::vector<int> *td_pdg = 0;

  tree->SetBranchAddress("td_energy", &td_energy);
  tree->SetBranchAddress("td_reco_status", &td_reco_status);
  tree->SetBranchAddress("td_pdg", &td_pdg);

  double bins[23];
  for(int i = 0; i < 16; ++i)
    bins[i] = i * .5;
  for(int i = 16; i < 19; ++i)
    bins[i] = 8 + (i - 16) * 1;
  bins[19] = 11; bins[20] = 15; bins[21] = 20; bins[22] = 25;

  TH1F *hTrueDepositEnergy = new TH1F("hTrueDepositEnergy", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyReco = new TH1F("hTrueDepositEnergyReco", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyNoDropped = new TH1F("hTrueDepositEnergyNoDropped", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyNoDroppedReco = new TH1F("hTrueDepositEnergyRecoNoDropped", ";True energy (MeV);Deposits", 22, bins);

  const unsigned N = tree->GetEntries();

  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);

      for(unsigned ii = 0; ii < td_energy->size(); ++ii)
	{
	  hTrueDepositEnergy->Fill(1e3*td_energy->at(ii));
	  
	  if(td_reco_status->at(ii))
	    hTrueDepositEnergyReco->Fill(1e3*td_energy->at(ii));

	  if(td_pdg->at(ii) == -999999) 
	    continue;

	  hTrueDepositEnergyNoDropped->Fill(1e3*td_energy->at(ii));
	  
	  if(td_reco_status->at(ii))
	    hTrueDepositEnergyNoDroppedReco->Fill(1e3*td_energy->at(ii));
	}
    }

  for(unsigned i = 1; i <= hTrueDepositEnergy->GetNbinsX(); ++i)
    {
      hTrueDepositEnergy->SetBinContent(i, hTrueDepositEnergy->GetBinContent(i) / hTrueDepositEnergy->GetBinWidth(i));
      hTrueDepositEnergyReco->SetBinContent(i, hTrueDepositEnergyReco->GetBinContent(i) / hTrueDepositEnergyReco->GetBinWidth(i));
      hTrueDepositEnergyNoDropped->SetBinContent(i, hTrueDepositEnergyNoDropped->GetBinContent(i) / hTrueDepositEnergyNoDropped->GetBinWidth(i));
      hTrueDepositEnergyNoDroppedReco->SetBinContent(i, hTrueDepositEnergyNoDroppedReco->GetBinContent(i) / hTrueDepositEnergyNoDroppedReco->GetBinWidth(i));
    }

  TCanvas *cEnergyRecoEff = new TCanvas("cEnergyRecoEff", "cEnergyRecoEff");
  cEnergyRecoEff->cd();
  
  TH1F *hTrueDepositEnergyClone = (TH1F*) hTrueDepositEnergy->Clone("hTrueDepositEnergyClone");
  hTrueDepositEnergyClone->Scale(1 / hTrueDepositEnergyClone->GetMaximum());
  hTrueDepositEnergyClone->SetLineColor(kGray+2);

  TEfficiency *eEnergyRecoEff = new TEfficiency(*hTrueDepositEnergyReco, *hTrueDepositEnergy);
  eEnergyRecoEff->SetTitle(";True energy (MeV);Efficiency");
  eEnergyRecoEff->SetLineColor(kRed+2);
  eEnergyRecoEff->SetMarkerColor(kRed+2);
  eEnergyRecoEff->SetLineWidth(2);

  eEnergyRecoEff->Draw();
  hTrueDepositEnergyClone->Draw("samehist");

  TLegend *lEnergyRecoEff = new TLegend(.6,.45,.85,.55);
  lEnergyRecoEff->SetBorderSize(0);
  lEnergyRecoEff->AddEntry(eEnergyRecoEff, "Efficiency","ple");
  lEnergyRecoEff->AddEntry(hTrueDepositEnergyClone, "True Distribution","l");
  lEnergyRecoEff->Draw();

  if(save)
    {
      cEnergyRecoEff->SaveAs(saveDir + "/true_energy_deposit_reco_eff.png");
      cEnergyRecoEff->SaveAs(saveDir + "/true_energy_deposit_reco_eff.pdf");
    }

  TCanvas *cEnergyNoDroppedRecoEff = new TCanvas("cEnergyNoDroppedRecoEff", "cEnergyNoDroppedRecoEff");
  cEnergyNoDroppedRecoEff->cd();
  
  TH1F *hTrueDepositEnergyNoDroppedClone = (TH1F*) hTrueDepositEnergyNoDropped->Clone("hTrueDepositEnergyNoDroppedClone");
  hTrueDepositEnergyNoDroppedClone->Scale(1 / hTrueDepositEnergyNoDroppedClone->GetMaximum());
  hTrueDepositEnergyNoDroppedClone->SetLineColor(kGray+2);

  TEfficiency *eEnergyNoDroppedRecoEff = new TEfficiency(*hTrueDepositEnergyNoDroppedReco, *hTrueDepositEnergyNoDropped);
  eEnergyNoDroppedRecoEff->SetTitle(";True energy (MeV);Efficiency");
  eEnergyNoDroppedRecoEff->SetLineColor(kRed+2);
  eEnergyNoDroppedRecoEff->SetMarkerColor(kRed+2);
  eEnergyNoDroppedRecoEff->SetLineWidth(2);

  eEnergyNoDroppedRecoEff->Draw();
  hTrueDepositEnergyNoDroppedClone->Draw("samehist");

  TLegend *lEnergyNoDroppedRecoEff = new TLegend(.6,.45,.85,.55);
  lEnergyNoDroppedRecoEff->SetBorderSize(0);
  lEnergyNoDroppedRecoEff->AddEntry(eEnergyNoDroppedRecoEff, "Efficiency","ple");
  lEnergyNoDroppedRecoEff->AddEntry(hTrueDepositEnergyNoDroppedClone, "True Distribution","l");
  lEnergyNoDroppedRecoEff->Draw();

  if(save)
    {
      cEnergyNoDroppedRecoEff->SaveAs(saveDir + "/true_energy_deposit_no_dropped_reco_eff.png");
      cEnergyNoDroppedRecoEff->SaveAs(saveDir + "/true_energy_deposit_no_dropped_reco_eff.pdf");
    }
}
