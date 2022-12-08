void ExtrapolatingDataRates()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/plots/extrapolatingdatarates";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");

  std::vector<uint64_t> *chit_unix_s = 0;
  std::vector<double> *chit_z = 0, *chit_t1 = 0;

  tree->SetBranchAddress("chit_unix_s", &chit_unix_s);
  tree->SetBranchAddress("chit_z", &chit_z);
  tree->SetBranchAddress("chit_t1", &chit_t1);

  TH1F *hEventRate = new TH1F("hEventRate", ";Hours into run;Events", 180, 0, 180);
  TH1F *hEventRateStart = new TH1F("hEventRateStart", ";Hours into run;Events", 100, 70, 80);
  TH1F *hEventRateRun = new TH1F("hEventRateRun", ";Hours into run;Events", 66, 74, 140);
  TH1F *hHitsPerEvent = new TH1F("hHitsPerEvent", ";Hits per event;Events", 30, 0, 30);
  TH1F *hBeamWindow = new TH1F("hBeamWindow", ";Hit t1;Events", 100, 332000, 336000);
  TH1F *hBeamHitsPerEvent = new TH1F("hBeamHitsPerEvent", ";Beam hits per event;Events", 30, 0, 30);

  const int N = tree->GetEntries();

  uint64_t run_start = std::numeric_limits<uint64_t>::max();

  for(int i = 0; i < N; ++i)
    {
      tree->GetEntry(i);

      for(int ii = 0; ii < chit_unix_s->size(); ++ii)
        {
          uint64_t time = chit_unix_s->at(ii);

          if(time < run_start)
            run_start = time;
        }
    }

  for(int i = 0; i < N; ++i)
    {
      tree->GetEntry(i);

      int upstream_hits = 0, upstream_beam_hits = 0;

      for(int ii = 0; ii < chit_unix_s->size(); ++ii)
        {
          if(chit_z->at(ii) > 0) continue;
          ++upstream_hits;
          uint64_t time = chit_unix_s->at(ii);
          hEventRate->Fill((time - run_start) / 3600.);
          hEventRateStart->Fill((time - run_start) / 3600.);
          hEventRateRun->Fill((time - run_start) / 3600.);
          double t1 = chit_t1->at(ii);
          hBeamWindow->Fill(t1);

          if(t1 > 332900 && t1 < 334400)
            ++upstream_beam_hits;
        }

      hHitsPerEvent->Fill(upstream_hits);
      hBeamHitsPerEvent->Fill(upstream_beam_hits);
    }

  TCanvas *cEventRate = new TCanvas("cEventRate", "cEventRate");
  cEventRate->cd();

  hEventRate->Draw("hist");

  TCanvas *cEventRateStart = new TCanvas("cEventRateStart", "cEventRateStart");
  cEventRateStart->cd();

  hEventRateStart->Draw("hist");

  TCanvas *cEventRateRun = new TCanvas("cEventRateRun", "cEventRateRun");
  cEventRateRun->cd();

  hEventRateRun->Draw("hist");
  TF1 *uniformFunc = new TF1("uniformFunc", "pol0", 80,130);
  TFitResultPtr fEventRateRun = hEventRateRun->Fit("uniformFunc","S", "", 80, 130);
  uniformFunc->SetLineColor(kRed+1);
  uniformFunc->SetLineWidth(4);
  uniformFunc->Draw("same");

  TPaveText *ptUnif = new TPaveText(.6,.38,.8,.5, "NDC");
  ptUnif->AddText(Form("Fitted Rate: %.2f", uniformFunc->GetParameter(0)));
  ptUnif->SetFillColor(kWhite);
  ptUnif->SetLineColor(kBlack);
  ptUnif->SetLineWidth(4);
  ptUnif->SetTextColor(kRed+1);
  ptUnif->Draw("same");

  if(save)
    {
      cEventRateRun->SaveAs(saveDir + "/hits_per_hour_4525.png");
      cEventRateRun->SaveAs(saveDir + "/hits_per_hour_4525.pdf");
    }

  TCanvas *cHitsPerEvent = new TCanvas("cHitsPerEvent", "cHitsPerEvent");
  cHitsPerEvent->cd();

  hHitsPerEvent->Draw("hist");
  TF1 *poissonFunc = new TF1("poissonFunc", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", 0, 30);
  //  TF1 *poissonFunc = new TF1("poissonFunc", "[0] * TMath::Poisson(x, [1])", 0, 30);
  //  TF1 *poissonFunc = new TF1("poissonFunc", "gaus", 0, 30);
  poissonFunc->SetParameters(600,10, 1);
  hHitsPerEvent->Fit("poissonFunc");
  poissonFunc->SetLineColor(kRed+1);
  poissonFunc->SetLineWidth(4);
  poissonFunc->Draw("same");

  TPaveText *ptPoisson = new TPaveText(.7,.78,.9,.9, "NDC");
  ptPoisson->AddText(Form("Fitted Mean: %.2f", poissonFunc->GetParameter(1)));
  ptPoisson->SetFillColor(kWhite);
  ptPoisson->SetLineColor(kBlack);
  ptPoisson->SetLineWidth(4);
  ptPoisson->SetTextColor(kRed+1);
  ptPoisson->Draw("same");

  if(save)
    {
      cHitsPerEvent->SaveAs(saveDir + "/hits_per_event_4525.png");
      cHitsPerEvent->SaveAs(saveDir + "/hits_per_event_4525.pdf");
    }

  TCanvas *cBeamWindow = new TCanvas("cBeamWindow", "cBeamWindow");
  cBeamWindow->cd();

  hBeamWindow->Draw("hist");

  TCanvas *cBeamHitsPerEvent = new TCanvas("cBeamHitsPerEvent", "cBeamHitsPerEvent");
  cBeamHitsPerEvent->cd();

  hBeamHitsPerEvent->Draw("hist");
  TF1 *poissonFunc2 = new TF1("poissonFunc2", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", 0, 30);
  poissonFunc2->SetParameters(5000, 1, 1);
  hBeamHitsPerEvent->Fit("poissonFunc2");
  poissonFunc2->SetLineColor(kRed+1);
  poissonFunc2->SetLineWidth(4);
  poissonFunc2->Draw("same");

  TPaveText *ptPoisson2 = new TPaveText(.7,.78,.9,.9, "NDC");
  ptPoisson2->AddText(Form("Fitted Mean: %.2f", poissonFunc2->GetParameter(1)));
  ptPoisson2->SetFillColor(kWhite);
  ptPoisson2->SetLineColor(kBlack);
  ptPoisson2->SetLineWidth(4);
  ptPoisson2->SetTextColor(kRed+1);
  ptPoisson2->Draw("same");

  if(save)
    {
      cBeamHitsPerEvent->SaveAs(saveDir + "/beam_hits_per_event_4525.png");
      cBeamHitsPerEvent->SaveAs(saveDir + "/beam_hits_per_event_4525.pdf");
    }
}
