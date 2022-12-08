void BucketStructureT1(const double period = 18.94)
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/bucketstructuret1";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";
  const double beam_start = 332899.1;
  const double beam_end   = beam_start + 81 * period;

  std::stringstream ss; ss.precision(4); ss << period;
  const TString period_s = ss.str();

  const int period_int = period;
  const double period_frac = period - period_int;
  const TString period_int_s = std::to_string(period_int);
  ss.str(""); ss.precision(2); ss << period_frac * 100;
  const TString period_frac_s = ss.str();

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");

  gSystem->Exec("mkdir -p " + saveDir + "/" + run_name);

  std::vector<double> *chit_t1 = 0;
  std::vector<std::vector<uint16_t>> *chit_sipm_raw_adc = 0, *chit_sipm_feb_mac5 = 0;

  tree->SetBranchAddress("chit_t1", &chit_t1);
  tree->SetBranchAddress("chit_sipm_raw_adc", &chit_sipm_raw_adc);
  tree->SetBranchAddress("chit_sipm_feb_mac5", &chit_sipm_feb_mac5);

  const unsigned N = tree->GetEntries();
  std::cout << "==== TOTAL ENTRIES: " << N << std::endl;

  TH1F *ts1_hist = new TH1F("ts1_hist", ";ts1 (#mus);CRTHits", 50,332.5,334.5);
  TH1F *ts1_bucket_hist = new TH1F("ts1_bucket_hist", ";CRT - Beam | Mod (" + period_s + ") [ns];CRTHits", 20, 0, 20);

  for(unsigned i = 0; i < N; ++i)
    {
      if(i%10000 == 0) std::cout << i << " / " << N << std::endl;
      tree->GetEntry(i);
      
      for(unsigned ii = 0; ii < chit_t1->size(); ++ii)
        {
          bool saturation = false;
          for(auto adc : chit_sipm_raw_adc->at(ii))
            if(adc > 4080) saturation = true;

          if(saturation)
            continue;

          double time = chit_t1->at(ii);
          ts1_hist->Fill(1e-3 * time);
          if(time > beam_start && time < beam_end)
            {
              double time_wrt_beam_start = time - beam_start;
              int bucketnumber = time_wrt_beam_start / period;
              double remainder = time_wrt_beam_start - period * bucketnumber;
              
              ts1_bucket_hist->Fill(remainder);
            }
        }
    }

  TCanvas *ts1_canvas = new TCanvas("ts1_canvas", "ts1_canvas");
  ts1_canvas->cd();
  ts1_canvas->SetRightMargin(.17);
  ts1_canvas->SetLeftMargin(.18);
  ts1_canvas->SetTopMargin(.13);
  ts1_hist->SetLineColor(kMagenta+2);
  ts1_hist->Draw("hist");

  if(save)
    {
      ts1_canvas->SaveAs(saveDir + "/" + run_name + "/ts1_beam_window.png");
      ts1_canvas->SaveAs(saveDir + "/" + run_name + "/ts1_beam_window.pdf");
    }

  TCanvas *ts1_bucket_canvas = new TCanvas("ts1_bucket_canvas", "ts1_bucket_canvas");
  ts1_bucket_canvas->cd();
  ts1_bucket_canvas->SetRightMargin(.15);
  ts1_bucket_canvas->SetLeftMargin(.15);
  ts1_bucket_canvas->SetTopMargin(.1);
  ts1_bucket_hist->SetLineColor(kCyan+2);
  ts1_bucket_hist->SetMinimum(0);
  ts1_bucket_hist->Draw("histE");

  TPaveText *pt = new TPaveText(.18,.79,.31,.87, "NDC");
  pt->AddText("SBND CRT##");
  pt->AddText("Run4500 Data");
  pt->SetFillColor(kWhite);
  pt->SetTextColor(kGray+2);
  pt->SetTextSize(0.035);
  pt->Draw("same");

  if(save)
    {
      ts1_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts1_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + ".png");
      ts1_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts1_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + ".pdf");
    }

  TF1 *gausFunc = new TF1("gausFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", 0, period);
  gausFunc->SetParameters(500,10,4,100);
  TFitResultPtr fit = ts1_bucket_hist->Fit("gausFunc", "S");
  gausFunc->SetLineColor(kRed+2);
  gausFunc->SetLineWidth(4);
  gausFunc->Draw("same");

  TPaveText *ptFit = new TPaveText(.65,.71,.8,.85, "NDC");
  ptFit->AddText("Gaussian Fit");
  ptFit->AddText(Form("Mean:      %.3f", gausFunc->GetParameter(1)));
  ptFit->AddText(Form("Sigma:     %.3f", gausFunc->GetParameter(2)));
  ptFit->AddText(Form("Bkgd Rate: %.3f", gausFunc->GetParameter(3)));
  ptFit->AddText(Form("#chi^{2} / dof:   %.3f / %i", fit->Chi2(), fit->Ndf()));
  ptFit->SetFillColor(kWhite);
  ptFit->SetTextColor(kBlack);
  ptFit->SetTextAlign(11);
  ptFit->SetTextSize(0.025);
  TText *t1 = ptFit->GetLineWith("Gaus");
  t1->SetTextColor(kRed+2);
  ptFit->Draw("same");

  if(save)
    {
      ts1_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts1_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + "_fit.png");
      ts1_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts1_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + "_fit.pdf");
    }
}
