void MatchingPerformance()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/clustering/plots/v09_67_00_no_duplicates/matchingperformance";
  gSystem->Exec("mkdir -p " + saveDir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_67_00_no_duplicates.root");

  std::vector<double> *tpc_truth_energy = 0, *tpc_length = 0, *tpc_sp_score = 0, *tpc_tr_score = 0;
  std::vector<bool> *tpc_sp_matchable = 0, *tpc_sp_matched = 0, *tpc_sp_good_match = 0;
  std::vector<bool> *tpc_tr_matchable = 0, *tpc_tr_matched = 0, *tpc_tr_good_match = 0;

  tree->SetBranchAddress("tpc_truth_energy", &tpc_truth_energy);
  tree->SetBranchAddress("tpc_length", &tpc_length);
  tree->SetBranchAddress("tpc_sp_score", &tpc_sp_score);
  tree->SetBranchAddress("tpc_tr_score", &tpc_tr_score);
  tree->SetBranchAddress("tpc_sp_matchable", &tpc_sp_matchable);
  tree->SetBranchAddress("tpc_sp_matched", &tpc_sp_matched);
  tree->SetBranchAddress("tpc_sp_good_match", &tpc_sp_good_match);
  tree->SetBranchAddress("tpc_tr_matchable", &tpc_tr_matchable);
  tree->SetBranchAddress("tpc_tr_matched", &tpc_tr_matched);
  tree->SetBranchAddress("tpc_tr_good_match", &tpc_tr_good_match);

  double bins_energy[17];
  for(int i = 0; i < 10; ++i)
    bins_energy[i] = i * 2;
  for(int i = 10; i < 14; ++i)
    bins_energy[i] = 20 + (i - 10) * 5;
  bins_energy[14] = 40; bins_energy[15] = 60; bins_energy[16] = 100;

  for(unsigned i = 0; i < 17; ++i)
    std::cout << bins_energy[i] << " ";
  std::cout << std::endl;

  double bins_length[21];
  for(int i = 0; i < 10; ++i)
    bins_length[i] = i * 20;
  for(int i = 10; i < 14; ++i)
    bins_length[i] = 200 + (i - 10) * 50;
  for(int i = 14; i < 19; ++i)
    bins_length[i] = 400 + (i - 14) * 20;
  bins_length[19] = 500; bins_length[20] = 800;

  for(unsigned i = 0; i < 21; ++i)
    std::cout << bins_length[i] << " ";
  std::cout << std::endl;

  TH1F *hSPMatchableEnergy = new TH1F("hSPMatchableEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hSPMatchableMatchedEnergy = new TH1F("hSPMatchableMatchedEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hTRMatchableEnergy = new TH1F("hTRMatchableEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hTRMatchableMatchedEnergy = new TH1F("hTRMatchableMatchedEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hSPMatchedEnergy = new TH1F("hSPMatchedEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hSPWellMatchedEnergy = new TH1F("hSPWellMatchedEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hTRMatchedEnergy = new TH1F("hTRMatchedEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);
  TH1F *hTRWellMatchedEnergy = new TH1F("hTRWellMatchedEnergy", ";True energy (GeV);TPC Tracks", 16, bins_energy);

  TH1F *hSPMatchableLength = new TH1F("hSPMatchableLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hSPMatchableMatchedLength = new TH1F("hSPMatchableMatchedLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hTRMatchableLength = new TH1F("hTRMatchableLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hTRMatchableMatchedLength = new TH1F("hTRMatchableMatchedLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hSPMatchedLength = new TH1F("hSPMatchedLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hSPWellMatchedLength = new TH1F("hSPWellMatchedLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hTRMatchedLength = new TH1F("hTRMatchedLength", ";Track length (cm);TPC Tracks", 20, bins_length);
  TH1F *hTRWellMatchedLength = new TH1F("hTRWellMatchedLength", ";Track length (cm);TPC Tracks", 20, bins_length);

  TH1F *hSPMatchableShortLength = new TH1F("hSPMatchableShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hSPMatchableMatchedShortLength = new TH1F("hSPMatchableMatchedShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hTRMatchableShortLength = new TH1F("hTRMatchableShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hTRMatchableMatchedShortLength = new TH1F("hTRMatchableMatchedShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hSPMatchedShortLength = new TH1F("hSPMatchedShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hSPWellMatchedShortLength = new TH1F("hSPWellMatchedShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hTRMatchedShortLength = new TH1F("hTRMatchedShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);
  TH1F *hTRWellMatchedShortLength = new TH1F("hTRWellMatchedShortLength", ";Track length (cm);TPC Tracks", 100, 0, 100);

  TH1F *hSPMatchedScore = new TH1F("hSPMatchedScore", ";DCA (cm);TPC Tracks", 50, 0, 200);
  TH1F *hSPWellMatchedScore = new TH1F("hSPWellMatchedScore", ";DCA (cm);TPC Tracks", 50, 0, 200);
  TH1F *hTRMatchedScore = new TH1F("hTRMatchedScore", ";DCA (cm) + 4<#theta> (#circ);TPC Tracks", 50, 0, 200);
  TH1F *hTRWellMatchedScore = new TH1F("hTRWellMatchedScore", ";DCA (cm) + 4<#theta> (#circ);TPC Tracks", 50, 0, 200);

  const unsigned N = tree->GetEntries();

  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);

      for(unsigned ii = 0; ii < tpc_sp_matchable->size(); ++ii)
        {
          if(tpc_sp_matchable->at(ii))
            {
              hSPMatchableEnergy->Fill(tpc_truth_energy->at(ii));
              hSPMatchableLength->Fill(tpc_length->at(ii));
              hSPMatchableShortLength->Fill(tpc_length->at(ii));
              if(tpc_sp_matched->at(ii))
                {
                  hSPMatchableMatchedEnergy->Fill(tpc_truth_energy->at(ii));
                  hSPMatchableMatchedLength->Fill(tpc_length->at(ii));
                  hSPMatchableMatchedShortLength->Fill(tpc_length->at(ii));
                }
            }

          if(tpc_sp_matched->at(ii))
            {
              hSPMatchedEnergy->Fill(tpc_truth_energy->at(ii));
              hSPMatchedLength->Fill(tpc_length->at(ii));
              hSPMatchedShortLength->Fill(tpc_length->at(ii));
              hSPMatchedScore->Fill(tpc_sp_score->at(ii));
              if(tpc_sp_good_match->at(ii))
                {
                  hSPWellMatchedEnergy->Fill(tpc_truth_energy->at(ii));
                  hSPWellMatchedLength->Fill(tpc_length->at(ii));
                  hSPWellMatchedShortLength->Fill(tpc_length->at(ii));
                  hSPWellMatchedScore->Fill(tpc_sp_score->at(ii));
                }
            }

          if(tpc_tr_matchable->at(ii))
            {
              hTRMatchableEnergy->Fill(tpc_truth_energy->at(ii));
              hTRMatchableLength->Fill(tpc_length->at(ii));
              hTRMatchableShortLength->Fill(tpc_length->at(ii));
              if(tpc_tr_matched->at(ii))
                {
                  hTRMatchableMatchedEnergy->Fill(tpc_truth_energy->at(ii));
                  hTRMatchableMatchedLength->Fill(tpc_length->at(ii));
                  hTRMatchableMatchedShortLength->Fill(tpc_length->at(ii));
                }
            }

          if(tpc_tr_matched->at(ii))
            {
              hTRMatchedEnergy->Fill(tpc_truth_energy->at(ii));
              hTRMatchedLength->Fill(tpc_length->at(ii));
              hTRMatchedShortLength->Fill(tpc_length->at(ii));
              hTRMatchedScore->Fill(tpc_tr_score->at(ii));
              if(tpc_tr_good_match->at(ii))
                {
                  hTRWellMatchedEnergy->Fill(tpc_truth_energy->at(ii));
                  hTRWellMatchedLength->Fill(tpc_length->at(ii));
                  hTRWellMatchedShortLength->Fill(tpc_length->at(ii));
                  hTRWellMatchedScore->Fill(tpc_tr_score->at(ii));
                }
            }
        }
    }

  for(unsigned i = 1; i <= hSPMatchableEnergy->GetNbinsX(); ++i)
    {
      hSPMatchableEnergy->SetBinContent(i, hSPMatchableEnergy->GetBinContent(i) / hSPMatchableEnergy->GetBinWidth(i));
      hSPMatchableMatchedEnergy->SetBinContent(i, hSPMatchableMatchedEnergy->GetBinContent(i) / hSPMatchableMatchedEnergy->GetBinWidth(i));
      hSPMatchedEnergy->SetBinContent(i, hSPMatchedEnergy->GetBinContent(i) / hSPMatchedEnergy->GetBinWidth(i));
      hSPWellMatchedEnergy->SetBinContent(i, hSPWellMatchedEnergy->GetBinContent(i) / hSPWellMatchedEnergy->GetBinWidth(i));
      hTRMatchableEnergy->SetBinContent(i, hTRMatchableEnergy->GetBinContent(i) / hTRMatchableEnergy->GetBinWidth(i));
      hTRMatchableMatchedEnergy->SetBinContent(i, hTRMatchableMatchedEnergy->GetBinContent(i) / hTRMatchableMatchedEnergy->GetBinWidth(i));
      hTRMatchedEnergy->SetBinContent(i, hTRMatchedEnergy->GetBinContent(i) / hTRMatchedEnergy->GetBinWidth(i));
      hTRWellMatchedEnergy->SetBinContent(i, hTRWellMatchedEnergy->GetBinContent(i) / hTRWellMatchedEnergy->GetBinWidth(i));
    }

  for(unsigned i = 1; i <= hSPMatchableLength->GetNbinsX(); ++i)
    {
      hSPMatchableLength->SetBinContent(i, hSPMatchableLength->GetBinContent(i) / hSPMatchableLength->GetBinWidth(i));
      hSPMatchableMatchedLength->SetBinContent(i, hSPMatchableMatchedLength->GetBinContent(i) / hSPMatchableMatchedLength->GetBinWidth(i));
      hSPMatchedLength->SetBinContent(i, hSPMatchedLength->GetBinContent(i) / hSPMatchedLength->GetBinWidth(i));
      hSPWellMatchedLength->SetBinContent(i, hSPWellMatchedLength->GetBinContent(i) / hSPWellMatchedLength->GetBinWidth(i));
      hTRMatchableLength->SetBinContent(i, hTRMatchableLength->GetBinContent(i) / hTRMatchableLength->GetBinWidth(i));
      hTRMatchableMatchedLength->SetBinContent(i, hTRMatchableMatchedLength->GetBinContent(i) / hTRMatchableMatchedLength->GetBinWidth(i));
      hTRMatchedLength->SetBinContent(i, hTRMatchedLength->GetBinContent(i) / hTRMatchedLength->GetBinWidth(i));
      hTRWellMatchedLength->SetBinContent(i, hTRWellMatchedLength->GetBinContent(i) / hTRWellMatchedLength->GetBinWidth(i));
    }

  for(unsigned i = 1; i <= hSPMatchableShortLength->GetNbinsX(); ++i)
    {
      hSPMatchableShortLength->SetBinContent(i, hSPMatchableShortLength->GetBinContent(i) / hSPMatchableShortLength->GetBinWidth(i));
      hSPMatchableMatchedShortLength->SetBinContent(i, hSPMatchableMatchedShortLength->GetBinContent(i) / hSPMatchableMatchedShortLength->GetBinWidth(i));
      hSPMatchedShortLength->SetBinContent(i, hSPMatchedShortLength->GetBinContent(i) / hSPMatchedShortLength->GetBinWidth(i));
      hSPWellMatchedShortLength->SetBinContent(i, hSPWellMatchedShortLength->GetBinContent(i) / hSPWellMatchedShortLength->GetBinWidth(i));
      hTRMatchableShortLength->SetBinContent(i, hTRMatchableShortLength->GetBinContent(i) / hTRMatchableShortLength->GetBinWidth(i));
      hTRMatchableMatchedShortLength->SetBinContent(i, hTRMatchableMatchedShortLength->GetBinContent(i) / hTRMatchableMatchedShortLength->GetBinWidth(i));
      hTRMatchedShortLength->SetBinContent(i, hTRMatchedShortLength->GetBinContent(i) / hTRMatchedShortLength->GetBinWidth(i));
      hTRWellMatchedShortLength->SetBinContent(i, hTRWellMatchedShortLength->GetBinContent(i) / hTRWellMatchedShortLength->GetBinWidth(i));
    }

  TCanvas *cSPMatchEffEnergy = new TCanvas("cSPMatchEffEnergy", "cSPMatchEffEnergy");
  cSPMatchEffEnergy->cd();

  TH1F *hSPMatchableEnergyClone = (TH1F*) hSPMatchableEnergy->Clone("hSPMatchableEnergy");
  hSPMatchableEnergyClone->Scale(1 / hSPMatchableEnergyClone->GetMaximum());
  hSPMatchableEnergyClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchEffEnergy = new TEfficiency(*hSPMatchableMatchedEnergy, *hSPMatchableEnergy);
  eSPMatchEffEnergy->SetTitle(";True energy (GeV);Efficiency");
  eSPMatchEffEnergy->SetLineColor(kRed+2);
  eSPMatchEffEnergy->SetMarkerColor(kRed+2);
  eSPMatchEffEnergy->SetLineWidth(2);

  eSPMatchEffEnergy->Draw();
  gPad->Update();
  eSPMatchEffEnergy->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchableEnergyClone->Draw("samehist");

  TLegend *lSPMatchEffEnergy = new TLegend(.6,.45,.85,.55);
  lSPMatchEffEnergy->SetBorderSize(0);
  lSPMatchEffEnergy->AddEntry(eSPMatchEffEnergy, "Efficiency","ple");
  lSPMatchEffEnergy->AddEntry(hSPMatchableEnergyClone, "True Distribution","l");
  lSPMatchEffEnergy->Draw();

  if(save)
    {
      cSPMatchEffEnergy->SaveAs(saveDir + "/space_point_matching_efficiency_energy.png");
      cSPMatchEffEnergy->SaveAs(saveDir + "/space_point_matching_efficiency_energy.pdf");
    }

  TCanvas *cTRMatchEffEnergy = new TCanvas("cTRMatchEffEnergy", "cTRMatchEffEnergy");
  cTRMatchEffEnergy->cd();

  TH1F *hTRMatchableEnergyClone = (TH1F*) hTRMatchableEnergy->Clone("hTRMatchableEnergy");
  hTRMatchableEnergyClone->Scale(1 / hTRMatchableEnergyClone->GetMaximum());
  hTRMatchableEnergyClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchEffEnergy = new TEfficiency(*hTRMatchableMatchedEnergy, *hTRMatchableEnergy);
  eTRMatchEffEnergy->SetTitle(";True energy (GeV);Efficiency");
  eTRMatchEffEnergy->SetLineColor(kRed+2);
  eTRMatchEffEnergy->SetMarkerColor(kRed+2);
  eTRMatchEffEnergy->SetLineWidth(2);

  eTRMatchEffEnergy->Draw();
  gPad->Update();
  eTRMatchEffEnergy->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchableEnergyClone->Draw("samehist");

  TLegend *lTRMatchEffEnergy = new TLegend(.6,.45,.85,.55);
  lTRMatchEffEnergy->SetBorderSize(0);
  lTRMatchEffEnergy->AddEntry(eTRMatchEffEnergy, "Efficiency","ple");
  lTRMatchEffEnergy->AddEntry(hTRMatchableEnergyClone, "True Distribution","l");
  lTRMatchEffEnergy->Draw();

  if(save)
    {
      cTRMatchEffEnergy->SaveAs(saveDir + "/track_matching_efficiency_energy.png");
      cTRMatchEffEnergy->SaveAs(saveDir + "/track_matching_efficiency_energy.pdf");
    }

  TCanvas *cSPMatchEffLength = new TCanvas("cSPMatchEffLength", "cSPMatchEffLength");
  cSPMatchEffLength->cd();

  TH1F *hSPMatchableLengthClone = (TH1F*) hSPMatchableLength->Clone("hSPMatchableLength");
  hSPMatchableLengthClone->Scale(1 / hSPMatchableLengthClone->GetMaximum());
  hSPMatchableLengthClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchEffLength = new TEfficiency(*hSPMatchableMatchedLength, *hSPMatchableLength);
  eSPMatchEffLength->SetTitle(";Track length (cm);Efficiency");
  eSPMatchEffLength->SetLineColor(kRed+2);
  eSPMatchEffLength->SetMarkerColor(kRed+2);
  eSPMatchEffLength->SetLineWidth(2);

  eSPMatchEffLength->Draw();
  gPad->Update();
  eSPMatchEffLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchableLengthClone->Draw("samehist");

  TLegend *lSPMatchEffLength = new TLegend(.6,.45,.85,.55);
  lSPMatchEffLength->SetBorderSize(0);
  lSPMatchEffLength->AddEntry(eSPMatchEffLength, "Efficiency","ple");
  lSPMatchEffLength->AddEntry(hSPMatchableLengthClone, "True Distribution","l");
  lSPMatchEffLength->Draw();

  if(save)
    {
      cSPMatchEffLength->SaveAs(saveDir + "/space_point_matching_efficiency_length.png");
      cSPMatchEffLength->SaveAs(saveDir + "/space_point_matching_efficiency_length.pdf");
    }

  TCanvas *cTRMatchEffLength = new TCanvas("cTRMatchEffLength", "cTRMatchEffLength");
  cTRMatchEffLength->cd();

  TH1F *hTRMatchableLengthClone = (TH1F*) hTRMatchableLength->Clone("hTRMatchableLength");
  hTRMatchableLengthClone->Scale(1 / hTRMatchableLengthClone->GetMaximum());
  hTRMatchableLengthClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchEffLength = new TEfficiency(*hTRMatchableMatchedLength, *hTRMatchableLength);
  eTRMatchEffLength->SetTitle(";Track length (cm);Efficiency");
  eTRMatchEffLength->SetLineColor(kRed+2);
  eTRMatchEffLength->SetMarkerColor(kRed+2);
  eTRMatchEffLength->SetLineWidth(2);

  eTRMatchEffLength->Draw();
  gPad->Update();
  eTRMatchEffLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchableLengthClone->Draw("samehist");

  TLegend *lTRMatchEffLength = new TLegend(.6,.45,.85,.55);
  lTRMatchEffLength->SetBorderSize(0);
  lTRMatchEffLength->AddEntry(eTRMatchEffLength, "Efficiency","ple");
  lTRMatchEffLength->AddEntry(hTRMatchableLengthClone, "True Distribution","l");
  lTRMatchEffLength->Draw();

  if(save)
    {
      cTRMatchEffLength->SaveAs(saveDir + "/track_matching_efficiency_length.png");
      cTRMatchEffLength->SaveAs(saveDir + "/track_matching_efficiency_length.pdf");
    }

  TCanvas *cSPMatchEffShortLength = new TCanvas("cSPMatchEffShortLength", "cSPMatchEffShortLength");
  cSPMatchEffShortLength->cd();

  TH1F *hSPMatchableShortLengthClone = (TH1F*) hSPMatchableShortLength->Clone("hSPMatchableShortLength");
  hSPMatchableShortLengthClone->Scale(1 / hSPMatchableShortLengthClone->GetMaximum());
  hSPMatchableShortLengthClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchEffShortLength = new TEfficiency(*hSPMatchableMatchedShortLength, *hSPMatchableShortLength);
  eSPMatchEffShortLength->SetTitle(";Track length (cm);Efficiency");
  eSPMatchEffShortLength->SetLineColor(kRed+2);
  eSPMatchEffShortLength->SetMarkerColor(kRed+2);
  eSPMatchEffShortLength->SetLineWidth(2);

  eSPMatchEffShortLength->Draw();
  gPad->Update();
  eSPMatchEffShortLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchableShortLengthClone->Draw("samehist");

  TLegend *lSPMatchEffShortLength = new TLegend(.6,.45,.85,.55);
  lSPMatchEffShortLength->SetBorderSize(0);
  lSPMatchEffShortLength->AddEntry(eSPMatchEffShortLength, "Efficiency","ple");
  lSPMatchEffShortLength->AddEntry(hSPMatchableShortLengthClone, "True Distribution","l");
  lSPMatchEffShortLength->Draw();

  if(save)
    {
      cSPMatchEffShortLength->SaveAs(saveDir + "/space_point_matching_efficiency_length_short.png");
      cSPMatchEffShortLength->SaveAs(saveDir + "/space_point_matching_efficiency_length_short.pdf");
    }

  TCanvas *cTRMatchEffShortLength = new TCanvas("cTRMatchEffShortLength", "cTRMatchEffShortLength");
  cTRMatchEffShortLength->cd();

  TH1F *hTRMatchableShortLengthClone = (TH1F*) hTRMatchableShortLength->Clone("hTRMatchableShortLength");
  hTRMatchableShortLengthClone->Scale(1 / hTRMatchableShortLengthClone->GetMaximum());
  hTRMatchableShortLengthClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchEffShortLength = new TEfficiency(*hTRMatchableMatchedShortLength, *hTRMatchableShortLength);
  eTRMatchEffShortLength->SetTitle(";Track length (cm);Efficiency");
  eTRMatchEffShortLength->SetLineColor(kRed+2);
  eTRMatchEffShortLength->SetMarkerColor(kRed+2);
  eTRMatchEffShortLength->SetLineWidth(2);

  eTRMatchEffShortLength->Draw();
  gPad->Update();
  eTRMatchEffShortLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchableShortLengthClone->Draw("samehist");

  TLegend *lTRMatchEffShortLength = new TLegend(.6,.45,.85,.55);
  lTRMatchEffShortLength->SetBorderSize(0);
  lTRMatchEffShortLength->AddEntry(eTRMatchEffShortLength, "Efficiency","ple");
  lTRMatchEffShortLength->AddEntry(hTRMatchableShortLengthClone, "True Distribution","l");
  lTRMatchEffShortLength->Draw();

  if(save)
    {
      cTRMatchEffShortLength->SaveAs(saveDir + "/track_matching_efficiency_length_short.png");
      cTRMatchEffShortLength->SaveAs(saveDir + "/track_matching_efficiency_length_short.pdf");
    }

  TCanvas *cSPMatchQualEnergy = new TCanvas("cSPMatchQualEnergy", "cSPMatchQualEnergy");
  cSPMatchQualEnergy->cd();

  TH1F *hSPMatchedEnergyClone = (TH1F*) hSPMatchedEnergy->Clone("hSPMatchedEnergy");
  hSPMatchedEnergyClone->Scale(1 / hSPMatchedEnergyClone->GetMaximum());
  hSPMatchedEnergyClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchQualEnergy = new TEfficiency(*hSPWellMatchedEnergy, *hSPMatchedEnergy);
  eSPMatchQualEnergy->SetTitle(";True energy (GeV);Quality");
  eSPMatchQualEnergy->SetLineColor(kRed+2);
  eSPMatchQualEnergy->SetMarkerColor(kRed+2);
  eSPMatchQualEnergy->SetLineWidth(2);

  eSPMatchQualEnergy->Draw();
  gPad->Update();
  eSPMatchQualEnergy->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchedEnergyClone->Draw("samehist");

  TLegend *lSPMatchQualEnergy = new TLegend(.6,.45,.85,.55);
  lSPMatchQualEnergy->SetBorderSize(0);
  lSPMatchQualEnergy->AddEntry(eSPMatchQualEnergy, "Quality","ple");
  lSPMatchQualEnergy->AddEntry(hSPMatchedEnergyClone, "True Distribution","l");
  lSPMatchQualEnergy->Draw();

  if(save)
    {
      cSPMatchQualEnergy->SaveAs(saveDir + "/space_point_matching_quality_energy.png");
      cSPMatchQualEnergy->SaveAs(saveDir + "/space_point_matching_quality_energy.pdf");
    }

  TCanvas *cTRMatchQualEnergy = new TCanvas("cTRMatchQualEnergy", "cTRMatchQualEnergy");
  cTRMatchQualEnergy->cd();

  TH1F *hTRMatchedEnergyClone = (TH1F*) hTRMatchedEnergy->Clone("hTRMatchedEnergy");
  hTRMatchedEnergyClone->Scale(1 / hTRMatchedEnergyClone->GetMaximum());
  hTRMatchedEnergyClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchQualEnergy = new TEfficiency(*hTRWellMatchedEnergy, *hTRMatchedEnergy);
  eTRMatchQualEnergy->SetTitle(";True energy (GeV);Quality");
  eTRMatchQualEnergy->SetLineColor(kRed+2);
  eTRMatchQualEnergy->SetMarkerColor(kRed+2);
  eTRMatchQualEnergy->SetLineWidth(2);

  eTRMatchQualEnergy->Draw();
  gPad->Update();
  eTRMatchQualEnergy->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchedEnergyClone->Draw("samehist");

  TLegend *lTRMatchQualEnergy = new TLegend(.6,.45,.85,.55);
  lTRMatchQualEnergy->SetBorderSize(0);
  lTRMatchQualEnergy->AddEntry(eTRMatchQualEnergy, "Quality","ple");
  lTRMatchQualEnergy->AddEntry(hTRMatchedEnergyClone, "True Distribution","l");
  lTRMatchQualEnergy->Draw();

  if(save)
    {
      cTRMatchQualEnergy->SaveAs(saveDir + "/track_matching_quality_energy.png");
      cTRMatchQualEnergy->SaveAs(saveDir + "/track_matching_quality_energy.pdf");
    }

  TCanvas *cSPMatchQualLength = new TCanvas("cSPMatchQualLength", "cSPMatchQualLength");
  cSPMatchQualLength->cd();

  TH1F *hSPMatchedLengthClone = (TH1F*) hSPMatchedLength->Clone("hSPMatchedLength");
  hSPMatchedLengthClone->Scale(1 / hSPMatchedLengthClone->GetMaximum());
  hSPMatchedLengthClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchQualLength = new TEfficiency(*hSPWellMatchedLength, *hSPMatchedLength);
  eSPMatchQualLength->SetTitle(";Track length (cm);Quality");
  eSPMatchQualLength->SetLineColor(kRed+2);
  eSPMatchQualLength->SetMarkerColor(kRed+2);
  eSPMatchQualLength->SetLineWidth(2);

  eSPMatchQualLength->Draw();
  gPad->Update();
  eSPMatchQualLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchedLengthClone->Draw("samehist");

  TLegend *lSPMatchQualLength = new TLegend(.6,.45,.85,.55);
  lSPMatchQualLength->SetBorderSize(0);
  lSPMatchQualLength->AddEntry(eSPMatchQualLength, "Quality","ple");
  lSPMatchQualLength->AddEntry(hSPMatchedLengthClone, "True Distribution","l");
  lSPMatchQualLength->Draw();

  if(save)
    {
      cSPMatchQualLength->SaveAs(saveDir + "/space_point_matching_quality_length.png");
      cSPMatchQualLength->SaveAs(saveDir + "/space_point_matching_quality_length.pdf");
    }

  TCanvas *cTRMatchQualLength = new TCanvas("cTRMatchQualLength", "cTRMatchQualLength");
  cTRMatchQualLength->cd();

  TH1F *hTRMatchedLengthClone = (TH1F*) hTRMatchedLength->Clone("hTRMatchedLength");
  hTRMatchedLengthClone->Scale(1 / hTRMatchedLengthClone->GetMaximum());
  hTRMatchedLengthClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchQualLength = new TEfficiency(*hTRWellMatchedLength, *hTRMatchedLength);
  eTRMatchQualLength->SetTitle(";Track length (cm);Quality");
  eTRMatchQualLength->SetLineColor(kRed+2);
  eTRMatchQualLength->SetMarkerColor(kRed+2);
  eTRMatchQualLength->SetLineWidth(2);

  eTRMatchQualLength->Draw();
  gPad->Update();
  eTRMatchQualLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchedLengthClone->Draw("samehist");

  TLegend *lTRMatchQualLength = new TLegend(.6,.45,.85,.55);
  lTRMatchQualLength->SetBorderSize(0);
  lTRMatchQualLength->AddEntry(eTRMatchQualLength, "Quality","ple");
  lTRMatchQualLength->AddEntry(hTRMatchedLengthClone, "True Distribution","l");
  lTRMatchQualLength->Draw();

  if(save)
    {
      cTRMatchQualLength->SaveAs(saveDir + "/track_matching_quality_length.png");
      cTRMatchQualLength->SaveAs(saveDir + "/track_matching_quality_length.pdf");
    }

  TCanvas *cSPMatchQualShortLength = new TCanvas("cSPMatchQualShortLength", "cSPMatchQualShortLength");
  cSPMatchQualShortLength->cd();

  TH1F *hSPMatchedShortLengthClone = (TH1F*) hSPMatchedShortLength->Clone("hSPMatchedShortLength");
  hSPMatchedShortLengthClone->Scale(1 / hSPMatchedShortLengthClone->GetMaximum());
  hSPMatchedShortLengthClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchQualShortLength = new TEfficiency(*hSPWellMatchedShortLength, *hSPMatchedShortLength);
  eSPMatchQualShortLength->SetTitle(";Track length (cm);Quality");
  eSPMatchQualShortLength->SetLineColor(kRed+2);
  eSPMatchQualShortLength->SetMarkerColor(kRed+2);
  eSPMatchQualShortLength->SetLineWidth(2);

  eSPMatchQualShortLength->Draw();
  gPad->Update();
  eSPMatchQualShortLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchedShortLengthClone->Draw("samehist");

  TLegend *lSPMatchQualShortLength = new TLegend(.6,.45,.85,.55);
  lSPMatchQualShortLength->SetBorderSize(0);
  lSPMatchQualShortLength->AddEntry(eSPMatchQualShortLength, "Quality","ple");
  lSPMatchQualShortLength->AddEntry(hSPMatchedShortLengthClone, "True Distribution","l");
  lSPMatchQualShortLength->Draw();

  if(save)
    {
      cSPMatchQualShortLength->SaveAs(saveDir + "/space_point_matching_quality_length_short.png");
      cSPMatchQualShortLength->SaveAs(saveDir + "/space_point_matching_quality_length_short.pdf");
    }

  TCanvas *cTRMatchQualShortLength = new TCanvas("cTRMatchQualShortLength", "cTRMatchQualShortLength");
  cTRMatchQualShortLength->cd();

  TH1F *hTRMatchedShortLengthClone = (TH1F*) hTRMatchedShortLength->Clone("hTRMatchedShortLength");
  hTRMatchedShortLengthClone->Scale(1 / hTRMatchedShortLengthClone->GetMaximum());
  hTRMatchedShortLengthClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchQualShortLength = new TEfficiency(*hTRWellMatchedShortLength, *hTRMatchedShortLength);
  eTRMatchQualShortLength->SetTitle(";Track length (cm);Quality");
  eTRMatchQualShortLength->SetLineColor(kRed+2);
  eTRMatchQualShortLength->SetMarkerColor(kRed+2);
  eTRMatchQualShortLength->SetLineWidth(2);

  eTRMatchQualShortLength->Draw();
  gPad->Update();
  eTRMatchQualShortLength->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchedShortLengthClone->Draw("samehist");

  TLegend *lTRMatchQualShortLength = new TLegend(.6,.45,.85,.55);
  lTRMatchQualShortLength->SetBorderSize(0);
  lTRMatchQualShortLength->AddEntry(eTRMatchQualShortLength, "Quality","ple");
  lTRMatchQualShortLength->AddEntry(hTRMatchedShortLengthClone, "True Distribution","l");
  lTRMatchQualShortLength->Draw();

  if(save)
    {
      cTRMatchQualShortLength->SaveAs(saveDir + "/track_matching_quality_length_short.png");
      cTRMatchQualShortLength->SaveAs(saveDir + "/track_matching_quality_length_short.pdf");
    }

  TCanvas *cSPMatchQualScore = new TCanvas("cSPMatchQualScore", "cSPMatchQualScore");
  cSPMatchQualScore->cd();

  TH1F *hSPMatchedScoreClone = (TH1F*) hSPMatchedScore->Clone("hSPMatchedScore");
  hSPMatchedScoreClone->Scale(1 / hSPMatchedScoreClone->GetMaximum());
  hSPMatchedScoreClone->SetLineColor(kGray+2);

  TEfficiency *eSPMatchQualScore = new TEfficiency(*hSPWellMatchedScore, *hSPMatchedScore);
  eSPMatchQualScore->SetTitle(";DCA (cm);Quality");
  eSPMatchQualScore->SetLineColor(kRed+2);
  eSPMatchQualScore->SetMarkerColor(kRed+2);
  eSPMatchQualScore->SetLineWidth(2);

  eSPMatchQualScore->Draw();
  gPad->Update();
  eSPMatchQualScore->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hSPMatchedScoreClone->Draw("samehist");

  TLegend *lSPMatchQualScore = new TLegend(.6,.45,.85,.55);
  lSPMatchQualScore->SetBorderSize(0);
  lSPMatchQualScore->AddEntry(eSPMatchQualScore, "Quality","ple");
  lSPMatchQualScore->AddEntry(hSPMatchedScoreClone, "True Distribution","l");
  lSPMatchQualScore->Draw();

  if(save)
    {
      cSPMatchQualScore->SaveAs(saveDir + "/space_point_matching_quality_score.png");
      cSPMatchQualScore->SaveAs(saveDir + "/space_point_matching_quality_score.pdf");
    }

  TCanvas *cTRMatchQualScore = new TCanvas("cTRMatchQualScore", "cTRMatchQualScore");
  cTRMatchQualScore->cd();

  TH1F *hTRMatchedScoreClone = (TH1F*) hTRMatchedScore->Clone("hTRMatchedScore");
  hTRMatchedScoreClone->Scale(1 / hTRMatchedScoreClone->GetMaximum());
  hTRMatchedScoreClone->SetLineColor(kGray+2);

  TEfficiency *eTRMatchQualScore = new TEfficiency(*hTRWellMatchedScore, *hTRMatchedScore);
  eTRMatchQualScore->SetTitle(";DCA (cm) + 4<#theta> (#circ);Quality");
  eTRMatchQualScore->SetLineColor(kRed+2);
  eTRMatchQualScore->SetMarkerColor(kRed+2);
  eTRMatchQualScore->SetLineWidth(2);

  eTRMatchQualScore->Draw();
  gPad->Update();
  eTRMatchQualScore->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
  hTRMatchedScoreClone->Draw("samehist");

  TLegend *lTRMatchQualScore = new TLegend(.6,.45,.85,.55);
  lTRMatchQualScore->SetBorderSize(0);
  lTRMatchQualScore->AddEntry(eTRMatchQualScore, "Quality","ple");
  lTRMatchQualScore->AddEntry(hTRMatchedScoreClone, "True Distribution","l");
  lTRMatchQualScore->Draw();

  if(save)
    {
      cTRMatchQualScore->SaveAs(saveDir + "/track_matching_quality_score.png");
      cTRMatchQualScore->SaveAs(saveDir + "/track_matching_quality_score.pdf");
    }
}
