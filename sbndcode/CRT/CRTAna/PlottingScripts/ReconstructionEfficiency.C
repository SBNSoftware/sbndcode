void ReconstructionEfficiency()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/v09_67_00/reconstructionefficiency";
  gSystem->Exec("mkdir -p " + saveDir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_67_00.root");

  std::vector<double> *td_tag_energy = 0;
  std::vector<bool> *td_tag_reco_status = 0;
  std::vector<double> *td_energy = 0;
  std::vector<bool> *td_reconstructable = 0;
  std::vector<bool> *td_reco_status = 0;
  std::vector<bool> *td_reco_triple = 0;

  tree->SetBranchAddress("td_tag_energy", &td_tag_energy);
  tree->SetBranchAddress("td_tag_reco_status", &td_tag_reco_status);

  tree->SetBranchAddress("td_energy", &td_energy);
  tree->SetBranchAddress("td_reconstructable", &td_reconstructable);
  tree->SetBranchAddress("td_reco_status", &td_reco_status);
  tree->SetBranchAddress("td_reco_triple", &td_reco_triple);

  double bins[23];
  for(int i = 0; i < 16; ++i)
    bins[i] = i * .5;
  for(int i = 16; i < 19; ++i)
    bins[i] = 8 + (i - 16) * 1;
  bins[19] = 11; bins[20] = 15; bins[21] = 20; bins[22] = 25;

  double binstracks[23];
  for(int i = 0; i < 23; ++i)
    binstracks[i] = 2 * bins[i];

  TH1F *hTrueDepositEnergy = new TH1F("hTrueDepositEnergy", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyReco = new TH1F("hTrueDepositEnergyReco", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueTrackDepositEnergy = new TH1F("hTrueTrackDepositEnergy", ";True energy (MeV);Deposits", 22, binstracks);
  TH1F *hTrueTrackDepositEnergyReco = new TH1F("hTrueTrackDepositEnergyReco", ";True energy (MeV);Deposits", 22, binstracks);
  TH1F *hTrue2TrackDepositEnergyReco = new TH1F("hTrue2TrackDepositEnergyReco", ";True energy (MeV);Deposits", 22, binstracks);
  TH1F *hTrue3TrackDepositEnergyReco = new TH1F("hTrue3TrackDepositEnergyReco", ";True energy (MeV);Deposits", 22, binstracks);

  const unsigned N = tree->GetEntries();

  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);

      for(unsigned ii = 0; ii < td_tag_energy->size(); ++ii)
        {
          hTrueDepositEnergy->Fill(1e3*td_tag_energy->at(ii));

          if(td_tag_reco_status->at(ii))
            hTrueDepositEnergyReco->Fill(1e3*td_tag_energy->at(ii));
        }

      for(unsigned ii = 0; ii < td_energy->size(); ++ii)
        {
          if(!td_reconstructable->at(ii))
            continue;

          hTrueTrackDepositEnergy->Fill(1e3*td_energy->at(ii));

          if(td_reco_status->at(ii))
            hTrueTrackDepositEnergyReco->Fill(1e3*td_energy->at(ii));

          if(td_reco_triple->at(ii) && td_reco_status->at(ii))
            hTrue3TrackDepositEnergyReco->Fill(1e3*td_energy->at(ii));

          if(!td_reco_triple->at(ii) && td_reco_status->at(ii))
            hTrue2TrackDepositEnergyReco->Fill(1e3*td_energy->at(ii));
        }
    }

  for(unsigned i = 1; i <= hTrueDepositEnergy->GetNbinsX(); ++i)
    {
      hTrueDepositEnergy->SetBinContent(i, hTrueDepositEnergy->GetBinContent(i) / hTrueDepositEnergy->GetBinWidth(i));
      hTrueDepositEnergyReco->SetBinContent(i, hTrueDepositEnergyReco->GetBinContent(i) / hTrueDepositEnergyReco->GetBinWidth(i));
      hTrueTrackDepositEnergy->SetBinContent(i, hTrueTrackDepositEnergy->GetBinContent(i) / hTrueTrackDepositEnergy->GetBinWidth(i));
      hTrueTrackDepositEnergyReco->SetBinContent(i, hTrueTrackDepositEnergyReco->GetBinContent(i) / hTrueTrackDepositEnergyReco->GetBinWidth(i));
      hTrue2TrackDepositEnergyReco->SetBinContent(i, hTrue2TrackDepositEnergyReco->GetBinContent(i) / hTrue2TrackDepositEnergyReco->GetBinWidth(i));
      hTrue3TrackDepositEnergyReco->SetBinContent(i, hTrue3TrackDepositEnergyReco->GetBinContent(i) / hTrue3TrackDepositEnergyReco->GetBinWidth(i));
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
  gPad->Update();
  eEnergyRecoEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
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

  TCanvas *cTrackEnergyRecoEff = new TCanvas("cTrackEnergyRecoEff", "cTrackEnergyRecoEff");
  cTrackEnergyRecoEff->cd();

  TH1F *hTrueTrackDepositEnergyClone = (TH1F*) hTrueTrackDepositEnergy->Clone("hTrueTrackDepositEnergyClone");
  hTrueTrackDepositEnergyClone->Scale(1 / hTrueTrackDepositEnergyClone->GetMaximum());
  hTrueTrackDepositEnergyClone->SetLineColor(kGray+2);

  TEfficiency *eTrackEnergyRecoEff = new TEfficiency(*hTrueTrackDepositEnergyReco, *hTrueTrackDepositEnergy);
  eTrackEnergyRecoEff->SetTitle(";True energy (MeV);Efficiency");
  eTrackEnergyRecoEff->SetLineColor(kRed+2);
  eTrackEnergyRecoEff->SetMarkerColor(kRed+2);
  eTrackEnergyRecoEff->SetLineWidth(2);

  eTrackEnergyRecoEff->Draw();
  gPad->Update();
  eTrackEnergyRecoEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTrueTrackDepositEnergyClone->Draw("samehist");

  TLegend *lTrackEnergyRecoEff = new TLegend(.6,.45,.85,.55);
  lTrackEnergyRecoEff->SetBorderSize(0);
  lTrackEnergyRecoEff->AddEntry(eTrackEnergyRecoEff, "Efficiency","ple");
  lTrackEnergyRecoEff->AddEntry(hTrueTrackDepositEnergyClone, "True Distribution","l");
  lTrackEnergyRecoEff->Draw();

  if(save)
    {
      cTrackEnergyRecoEff->SaveAs(saveDir + "/true_track_energy_deposit_reco_eff.png");
      cTrackEnergyRecoEff->SaveAs(saveDir + "/true_track_energy_deposit_reco_eff.pdf");
    }

  TCanvas *c2TrackEnergyRecoEff = new TCanvas("c2TrackEnergyRecoEff", "c2TrackEnergyRecoEff");
  c2TrackEnergyRecoEff->cd();

  TEfficiency *e2TrackEnergyRecoEff = new TEfficiency(*hTrue2TrackDepositEnergyReco, *hTrueTrackDepositEnergy);
  e2TrackEnergyRecoEff->SetTitle(";True energy (MeV);Efficiency");
  e2TrackEnergyRecoEff->SetLineColor(kRed+2);
  e2TrackEnergyRecoEff->SetMarkerColor(kRed+2);
  e2TrackEnergyRecoEff->SetLineWidth(2);

  e2TrackEnergyRecoEff->Draw();
  gPad->Update();
  e2TrackEnergyRecoEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTrueTrackDepositEnergyClone->Draw("samehist");

  TLegend *l2TrackEnergyRecoEff = new TLegend(.6,.45,.85,.55);
  l2TrackEnergyRecoEff->SetBorderSize(0);
  l2TrackEnergyRecoEff->AddEntry(e2TrackEnergyRecoEff, "Efficiency","ple");
  l2TrackEnergyRecoEff->AddEntry(hTrueTrackDepositEnergyClone, "True Distribution","l");
  l2TrackEnergyRecoEff->Draw();

  if(save)
    {
      c2TrackEnergyRecoEff->SaveAs(saveDir + "/true_track_two_sp_energy_deposit_reco_eff.png");
      c2TrackEnergyRecoEff->SaveAs(saveDir + "/true_track_two_sp_energy_deposit_reco_eff.pdf");
    }

  TCanvas *c3TrackEnergyRecoEff = new TCanvas("c3TrackEnergyRecoEff", "c3TrackEnergyRecoEff");
  c3TrackEnergyRecoEff->cd();

  TEfficiency *e3TrackEnergyRecoEff = new TEfficiency(*hTrue3TrackDepositEnergyReco, *hTrueTrackDepositEnergy);
  e3TrackEnergyRecoEff->SetTitle(";True energy (MeV);Efficiency");
  e3TrackEnergyRecoEff->SetLineColor(kRed+2);
  e3TrackEnergyRecoEff->SetMarkerColor(kRed+2);
  e3TrackEnergyRecoEff->SetLineWidth(2);

  e3TrackEnergyRecoEff->Draw();
  gPad->Update();
  e3TrackEnergyRecoEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTrueTrackDepositEnergyClone->Draw("samehist");

  TLegend *l3TrackEnergyRecoEff = new TLegend(.6,.45,.85,.55);
  l3TrackEnergyRecoEff->SetBorderSize(0);
  l3TrackEnergyRecoEff->AddEntry(e3TrackEnergyRecoEff, "Efficiency","ple");
  l3TrackEnergyRecoEff->AddEntry(hTrueTrackDepositEnergyClone, "True Distribution","l");
  l3TrackEnergyRecoEff->Draw();

  if(save)
    {
      c3TrackEnergyRecoEff->SaveAs(saveDir + "/true_track_three_sp_energy_deposit_reco_eff.png");
      c3TrackEnergyRecoEff->SaveAs(saveDir + "/true_track_three_sp_energy_deposit_reco_eff.pdf");
    }
}
