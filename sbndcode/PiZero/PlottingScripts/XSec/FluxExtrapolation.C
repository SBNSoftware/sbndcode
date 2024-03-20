#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

void Extrapolate(float &nu_x, float &nu_y, const float nu_z, const float &nu_other_x,
                 const float &nu_other_y, const float nu_other_z, const float extrap_z);

void FluxExtrapolation(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/flux_extrapolation";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_flux_configL_*.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *nus = new TChain("fluxana/tree");
  nus->Add(file);

  TChain *pot = new TChain("fluxana/pottree");
  pot->Add(file);

  const double nuPOT = GetPOT(pot);
  const double scaling = goalPOT / nuPOT;
  const double fv_face_area = 175 * 360 * 2.;

  const int N = nus->GetEntries();

  int nu_pdg;
  float nu_x, nu_y, nu_other_x, nu_other_y, nu_e;
  double ray_fv_length;

  nus->SetBranchStatus("*", 0);
  nus->SetBranchStatus("nu_x", 1);
  nus->SetBranchStatus("nu_y", 1);
  nus->SetBranchStatus("nu_other_x", 1);
  nus->SetBranchStatus("nu_other_y", 1);
  nus->SetBranchStatus("nu_pdg", 1);
  nus->SetBranchStatus("nu_e", 1);
  nus->SetBranchStatus("ray_fv_length", 1);

  nus->SetBranchAddress("nu_x", &nu_x);
  nus->SetBranchAddress("nu_y", &nu_y);
  nus->SetBranchAddress("nu_other_x", &nu_other_x);
  nus->SetBranchAddress("nu_other_y", &nu_other_y);
  nus->SetBranchAddress("nu_pdg", &nu_pdg);
  nus->SetBranchAddress("nu_e", &nu_e);
  nus->SetBranchAddress("ray_fv_length", &ray_fv_length);

  const int n_samples = 101;
  double baseline[n_samples], flux[n_samples];
  double ebaseline[n_samples], eflux[n_samples];

  double ray_traced_flux = 0;

  for(int s = 0; s < n_samples; ++s)
    {
      baseline[s]  = 11000 + 5 * s;
      flux[s]      = 0;
      ebaseline[s] = 0;
    }

  for(int i = 0; i < N; ++i)
    {
      if(!(i%1000000))
        std::cout << i << " / " << N << " (" << (100. * i) / N << "%)" << std::endl;

      nus->GetEntry(i);

      if(ray_fv_length > 0)
        ray_traced_flux += (ray_fv_length / 440);

      const float orig_x = nu_x;
      const float orig_y = nu_y;

      for(int s = 0; s < n_samples; ++s)
        {
          nu_x = orig_x;
          nu_y = orig_y;

          Extrapolate(nu_x, nu_y, 0, nu_other_x, nu_other_y, 500, baseline[s] - 11000);

          if(abs(nu_x) > 180 || abs(nu_x) < 5 || abs(nu_y) > 180)
            continue;

          flux[s] += 1;
        }
    }

  for(int s = 0; s < n_samples; ++s)
    {
      eflux[s]  = sqrt(flux[s]) / flux[s];
      flux[s]  *= (scaling / fv_face_area);
      eflux[s] *= flux[s];
    }

  ray_traced_flux *= (scaling / fv_face_area);

  TCanvas *c = new TCanvas("c", "c");
  c->SetTopMargin(.1);
  c->cd();

  TGraphErrors *g = new TGraphErrors(101, baseline, flux, ebaseline, eflux);
  g->Draw("AP");
  g->GetYaxis()->SetTitleOffset(1.15);
  g->GetXaxis()->SetLabelSize(0.06);
  g->SetTitle(";Baseline (cm);Integrated #nu Flux (cm^{-2})");

  c->SaveAs(saveDir + "/integrated_flux_baseline_dependence.png");
  c->SaveAs(saveDir + "/integrated_flux_baseline_dependence.pdf");

  TF1 *fLin = new TF1("lin", "[0] + [1] * x", baseline[0], baseline[100]);
  fLin->SetParameter(0, 5e13);
  fLin->SetParameter(1, -3e9);
  fLin->SetLineColor(kRed+2);

  TF1 *fRSq = new TF1("rSq", "[0] + [1] / (x * x)", baseline[0], baseline[100]);
  fRSq->SetParameter(1, 1e20);
  fRSq->SetLineColor(kBlue+2);
  fRSq->SetRange(baseline[0], baseline[100]);

  TFitResultPtr linFitRes = g->Fit(fLin);
  TFitResultPtr rsqFitRes = g->Fit(fRSq);

  fLin->Draw("same");
  fRSq->Draw("same");

  TPaveText *text = new TPaveText(.6, .6, .8, .85, "NDC");
  TText *title1 =   text->AddText("Linear Fit (a + bx)");
  title1->SetTextColor(kRed+2);
  TText *titleGap1 = text->AddText("");
  titleGap1->SetTextSize(0.005);
  text->AddText(Form("a = %.3e #pm %.3e", fLin->GetParameter(0), fLin->GetParError(0)));
  text->AddText(Form("b = %.3e #pm %.3e", fLin->GetParameter(1), fLin->GetParError(1)));
  text->AddText(Form("#chi^{2}/DoF = %.3f", fLin->GetChisquare() / fLin->GetNDF()));
  text->AddText("");
  text->AddText("");
  TText *title2 = text->AddText("R^{-2} Fit (a + bx^{-2})");
  title2->SetTextColor(kBlue+2);
  TText *titleGap2 = text->AddText("");
  titleGap2->SetTextSize(0.005);
  text->AddText(Form("a = %.3e #pm %.3e", fRSq->GetParameter(0), fRSq->GetParError(0)));
  text->AddText(Form("b = %.3e #pm %.3e", fRSq->GetParameter(1), fRSq->GetParError(1)));
  text->AddText(Form("#chi^{2}/DoF = %.3f", fRSq->GetChisquare() / fRSq->GetNDF()));
  text->SetTextAlign(12);
  text->SetTextSize(0.02);
  title1->SetTextSize(0.03);
  title2->SetTextSize(0.03);
  text->SetBorderSize(0);
  text->SetFillColor(kWhite);
  text->Draw();

  c->SaveAs(saveDir + "/integrated_flux_baseline_dependence_fit.png");
  c->SaveAs(saveDir + "/integrated_flux_baseline_dependence_fit.pdf");

  g->Draw("AP");
  fRSq->Draw("same");

  TLine *lowLine = new TLine(11010, g->GetYaxis()->GetXmin(), 11010, fRSq->Eval(11010));
  lowLine->SetNDC(false);
  lowLine->SetLineStyle(9);
  lowLine->SetLineWidth(3);
  lowLine->SetLineColor(kGray+2);

  lowLine->Draw();

  TLine *highLine = new TLine(11450, g->GetYaxis()->GetXmin(), 11450, fRSq->Eval(11450));
  highLine->SetNDC(false);
  highLine->SetLineStyle(9);
  highLine->SetLineWidth(3);
  highLine->SetLineColor(kGray+2);

  highLine->Draw();

  const double integral    = fRSq->Integral(11010, 11450);
  const double expectation = integral / 440;
  const double effbaseline = std::sqrt(fRSq->GetParameter(1) / (expectation - fRSq->GetParameter(0)));

  TLine *upLine = new TLine(effbaseline, g->GetYaxis()->GetXmin(), effbaseline, fRSq->Eval(effbaseline));
  upLine->SetNDC(false);
  upLine->SetLineStyle(7);
  upLine->SetLineWidth(4);
  upLine->SetLineColor(kRed+2);
  upLine->Draw();

  TLine *upLineLinear = new TLine(11230, g->GetYaxis()->GetXmin(), 11230, fRSq->Eval(11230));
  upLineLinear->SetNDC(false);
  upLineLinear->SetLineStyle(7);
  upLineLinear->SetLineWidth(4);
  upLineLinear->SetLineColor(kPink+2);
  upLineLinear->Draw();

  TLine *acrossLine = new TLine(g->GetXaxis()->GetXmin(), fRSq->Eval(effbaseline), effbaseline, fRSq->Eval(effbaseline));
  acrossLine->SetNDC(false);
  acrossLine->SetLineStyle(7);
  acrossLine->SetLineWidth(4);
  acrossLine->SetLineColor(kRed+2);
  acrossLine->Draw();

  TPaveText *text2 = new TPaveText(.4, .74, .8, .85, "NDC");
  text2->AddText(Form("Expected Value (11010 -> 11450): %.4e cm^{-2}", expectation));
  text2->AddText(Form("Effective Baseline: %.0f cm", effbaseline));
  text2->SetTextAlign(12);
  text2->SetTextSize(0.03);
  text2->SetBorderSize(0);
  text2->SetFillColor(kWhite);
  text2->Draw();

  std::cout << "Average from curve: " << expectation << std::endl;
  std::cout << "Effective baseline: " << effbaseline << std::endl;
  std::cout << "Curve evaluated at 1/2 FV: " << fRSq->Eval(11230) << std::endl;
  std::cout << "Ray-traced flux: " << ray_traced_flux << std::endl;

  c->SaveAs(saveDir + "/integrated_flux_baseline_dependence_expectation.png");
  c->SaveAs(saveDir + "/integrated_flux_baseline_dependence_expectation.pdf");
}

void Extrapolate(float &nu_x, float &nu_y, const float nu_z, const float &nu_other_x,
                 const float &nu_other_y, const float nu_other_z, const float extrap_z)
{
  const TVector3 start(nu_x, nu_y, nu_z);
  const TVector3 end(nu_other_x, nu_other_y, nu_other_z);

  const float k = extrap_z / (nu_other_z - nu_z);

  const TVector3 extrap = start + k * (end - start);

  nu_x = extrap.X();
  nu_y = extrap.Y();
}
