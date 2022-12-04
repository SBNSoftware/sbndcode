void BucketStructureT0TDC(const double period = 18.94)
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/bucketstructuret0tdc";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";
  const double beam_start = 332899.25;
  const double beam_end   = beam_start + 1500;
  const double beam_start_rwm = -1045.25;
  const double beam_end_rwm   = beam_start_rwm + 1500;

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

  std::vector<double> *chit_t0 = 0, *chit_unix_s = 0;
  std::vector<std::vector<uint16_t>> *chit_sipm_raw_adc = 0;
  std::vector<uint32_t> *tdc_channel = 0;
  std::vector<uint64_t> *tdc_timestamp = 0;

  tree->SetBranchAddress("chit_t0", &chit_t0);
  tree->SetBranchAddress("chit_unix_s", &chit_unix_s);
  tree->SetBranchAddress("chit_sipm_raw_adc", &chit_sipm_raw_adc);
  tree->SetBranchAddress("tdc_channel", &tdc_channel);
  tree->SetBranchAddress("tdc_timestamp", &tdc_timestamp);

  const unsigned N = tree->GetEntries();

  std::cout << "==== TOTAL ENTRIES: " << N << std::endl;

  TH1F *ts0_to_bes_hist = new TH1F("ts0_to_bes_hist", ";CRT t0 - TDC BES [#kern[0.4]{#mus}];CRTHits", 40,332.5,334.5);
  TH1F *ts0_to_bes_bucket_hist = new TH1F("ts0_to_bes_bucket_hist", ";CRT t0 - TDC BES | Mod (" + period_s + ") [ns];CRTHits", 20, 0, 20);

  TH1F *ts0_to_rwm_hist = new TH1F("ts0_to_rwm_hist", ";CRT t0 - TDC RWM [#kern[0.4]{#mus}];CRTHits", 40,-1.45,.55);
  TH1F *ts0_to_rwm_bucket_hist = new TH1F("ts0_to_rwm_bucket_hist", ";CRT t0 - TDC RWM | Mod (" + period_s + ") [ns];CRTHits", 20, 0, 20);

  for(unsigned i = 0; i < N; ++i)
    {
      if(i%10000 == 0) std::cout << i << " / " << N << std::endl;
      tree->GetEntry(i);

      bool found_bes = false, found_rwm = false, found_ftrig = false;
      uint64_t bes = 0, rwm = 0, ftrig = 0;

      for(unsigned ii = 0; ii < tdc_channel->size(); ++ii)
	{
	  if(tdc_channel->at(ii) == 1)
	    {
	      found_bes = true;
	      bes = tdc_timestamp->at(ii);
	    }
	  else if(tdc_channel->at(ii) == 2)
	    {
	      found_rwm = true;
	      rwm = tdc_timestamp->at(ii);
	    }
	  else if(tdc_channel->at(ii) == 3)
	    {
	      found_ftrig = true;
	      ftrig = tdc_timestamp->at(ii);
	    }

	}

      if(!found_bes)
	{
	  std::cout << "Failed to find BES for event: " << i << std::endl;
	  continue;
	}
      if(!found_rwm)
	{
	  std::cout << "Failed to find RWM for event: " << i << std::endl;
	  continue;
	}
      if(!found_ftrig)
	{
	  std::cout << "Failed to find FTRIG for event: " << i << std::endl;
	  continue;
	}

      // if(ftrig < rwm)
      // 	continue;

      for(unsigned ii = 0; ii < chit_t0->size(); ++ii)
	{
	  bool saturation = false;
	  for(auto adc : chit_sipm_raw_adc->at(ii))
	    if(adc > 4080) saturation = true;

	  if(saturation)
	    continue;

	  double t0 = chit_t0->at(ii);
	  double unixs = chit_unix_s->at(ii);

	  if(std::abs(bes / 1e9 - unixs) < std::numeric_limits<double>::epsilon())
	    std::cout << std::setprecision(9) << "Disagreement between unit part of BES: " << bes / 1e9 << " & UnixS: " << unixs << std::endl;

	  double bes_time = t0 - (bes % 1000000000 - 193);
	  ts0_to_bes_hist->Fill(1e-3 * bes_time);

	  if(bes_time > beam_start && bes_time < beam_end)
	    {
	      double bes_time_wrt_beam_start = bes_time - beam_start;
	      int bes_bucketnumber = bes_time_wrt_beam_start / period;
	      double bes_remainder = bes_time_wrt_beam_start - period * bes_bucketnumber;
	      ts0_to_bes_bucket_hist->Fill(bes_remainder);
	    }

	  if(std::abs(rwm / 1e9 - unixs) < std::numeric_limits<double>::epsilon())
	    std::cout << std::setprecision(9) << "Disagreement between unit part of RWM: " << rwm / 1e9 << " & UnixS: " << unixs << std::endl;

	  double rwm_time = t0 - (rwm % 1000000000 - 193);
	  ts0_to_rwm_hist->Fill(1e-3 * rwm_time);

	  if(rwm_time > beam_start_rwm && rwm_time < beam_end_rwm)
	    {
	      double rwm_time_wrt_beam_start = rwm_time - beam_start_rwm;
	      int rwm_bucketnumber = rwm_time_wrt_beam_start / period;
	      double rwm_remainder = rwm_time_wrt_beam_start - period * rwm_bucketnumber;
	      ts0_to_rwm_bucket_hist->Fill(rwm_remainder);
	    }
	}
    }

  TCanvas *ts0_to_bes_canvas = new TCanvas("ts0_to_bes_canvas", "ts0_to_bes_canvas");
  ts0_to_bes_canvas->cd();
  ts0_to_bes_canvas->SetRightMargin(.17);
  ts0_to_bes_canvas->SetLeftMargin(.18);
  ts0_to_bes_canvas->SetTopMargin(.13);
  ts0_to_bes_hist->SetLineColor(kMagenta+2);
  ts0_to_bes_hist->Draw("hist");

  if(save)
    {
      ts0_to_bes_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_bes_beam_window.png");
      ts0_to_bes_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_bes_beam_window.pdf");
    }

  TCanvas *ts0_to_bes_bucket_canvas = new TCanvas("ts0_to_bes_bucket_canvas", "ts0_to_bes_bucket_canvas");
  ts0_to_bes_bucket_canvas->cd();
  ts0_to_bes_bucket_canvas->SetRightMargin(.15);
  ts0_to_bes_bucket_canvas->SetLeftMargin(.15);
  ts0_to_bes_bucket_canvas->SetTopMargin(.1);
  ts0_to_bes_bucket_hist->SetLineColor(kCyan+2);
  ts0_to_bes_bucket_hist->SetMinimum(0);
  ts0_to_bes_bucket_hist->Draw("hist");

  TPaveText *pt = new TPaveText(.18,.79,.31,.87, "NDC");
  pt->AddText("SBND CRT##");
  pt->AddText("Run4500 Data");
  pt->SetFillColor(kWhite);
  pt->SetTextColor(kGray+2);
  pt->SetTextSize(0.035);
  pt->Draw("same");

  if(save)
    {
      ts0_to_bes_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_bes_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + ".png");
      ts0_to_bes_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_bes_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + ".pdf");
    }

  TF1 *gausFunc = new TF1("gausFunc", "gaus", 0, period);
  TFitResultPtr fit = ts0_to_bes_bucket_hist->Fit("gausFunc", "S");
  gausFunc->SetLineColor(kRed+2);
  gausFunc->SetLineWidth(4);
  gausFunc->Draw("same");

  TPaveText *ptFit = new TPaveText(.65,.73,.8,.85, "NDC");
  ptFit->AddText("Gaussian Fit");
  ptFit->AddText(Form("Mean:     %.3f", gausFunc->GetParameter("Mean")));
  ptFit->AddText(Form("Sigma:    %.3f", gausFunc->GetParameter("Sigma")));
  ptFit->AddText(Form("#chi^{2} / dof:  %.3f / %i", fit->Chi2(), fit->Ndf()));
  ptFit->SetFillColor(kWhite);
  ptFit->SetTextColor(kBlack);
  ptFit->SetTextAlign(11);
  ptFit->SetTextSize(0.025);
  TText *t1 = ptFit->GetLineWith("Gaus");
  t1->SetTextColor(kRed+2);
  ptFit->Draw("same");

  if(save)
    {
      ts0_to_bes_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_bes_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + "_fit.png");
      ts0_to_bes_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_bes_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + "_fit.pdf");
    }

  TCanvas *ts0_to_rwm_canvas = new TCanvas("ts0_to_rwm_canvas", "ts0_to_rwm_canvas");
  ts0_to_rwm_canvas->cd();
  ts0_to_rwm_canvas->SetRightMargin(.17);
  ts0_to_rwm_canvas->SetLeftMargin(.18);
  ts0_to_rwm_canvas->SetTopMargin(.13);
  ts0_to_rwm_hist->SetLineColor(kMagenta+2);
  ts0_to_rwm_hist->Draw("hist");

  if(save)
    {
      ts0_to_rwm_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_rwm_beam_window.png");
      ts0_to_rwm_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_rwm_beam_window.pdf");
    }

  TCanvas *ts0_to_rwm_bucket_canvas = new TCanvas("ts0_to_rwm_bucket_canvas", "ts0_to_rwm_bucket_canvas");
  ts0_to_rwm_bucket_canvas->cd();
  ts0_to_rwm_bucket_canvas->SetRightMargin(.15);
  ts0_to_rwm_bucket_canvas->SetLeftMargin(.15);
  ts0_to_rwm_bucket_canvas->SetTopMargin(.1);
  ts0_to_rwm_bucket_hist->SetLineColor(kCyan+2);
  ts0_to_rwm_bucket_hist->SetMinimum(0);
  ts0_to_rwm_bucket_hist->Draw("hist");

  TPaveText *pt_rwm = new TPaveText(.18,.79,.31,.87, "NDC");
  pt_rwm->AddText("SBND CRT##");
  pt_rwm->AddText("Run4500 Data");
  pt_rwm->SetFillColor(kWhite);
  pt_rwm->SetTextColor(kGray+2);
  pt_rwm->SetTextSize(0.035);
  pt_rwm->Draw("same");

  if(save)
    {
      ts0_to_rwm_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_rwm_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + ".png");
      ts0_to_rwm_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_rwm_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + ".pdf");
    }

  TF1 *gausFunc_rwm = new TF1("gausFunc_rwm", "gaus", 0, period);
  TFitResultPtr fit_rwm = ts0_to_rwm_bucket_hist->Fit("gausFunc_rwm", "S");
  gausFunc_rwm->SetLineColor(kRed+2);
  gausFunc_rwm->SetLineWidth(4);
  gausFunc_rwm->Draw("same");

  TPaveText *ptFit_rwm = new TPaveText(.65,.73,.8,.85, "NDC");
  ptFit_rwm->AddText("Gaussian Fit");
  ptFit_rwm->AddText(Form("Mean:     %.3f", gausFunc_rwm->GetParameter("Mean")));
  ptFit_rwm->AddText(Form("Sigma:    %.3f", gausFunc_rwm->GetParameter("Sigma")));
  ptFit_rwm->AddText(Form("#chi^{2} / dof:  %.3f / %i", fit_rwm->Chi2(), fit_rwm->Ndf()));
  ptFit_rwm->SetFillColor(kWhite);
  ptFit_rwm->SetTextColor(kBlack);
  ptFit_rwm->SetTextAlign(11);
  ptFit_rwm->SetTextSize(0.025);
  TText *t1_rwm = ptFit_rwm->GetLineWith("Gaus");
  t1_rwm->SetTextColor(kRed+2);
  ptFit_rwm->Draw("same");

  if(save)
    {
      ts0_to_rwm_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_rwm_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + "_fit.png");
      ts0_to_rwm_bucket_canvas->SaveAs(saveDir + "/" + run_name + "/ts0_to_rwm_overlapped_buckets_period_" + period_int_s + "p" + period_frac_s + "_fit.pdf");
    }
}
