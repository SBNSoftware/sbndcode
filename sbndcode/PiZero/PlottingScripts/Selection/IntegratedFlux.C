#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

const double goalPOT = 10e20;

double GetPOT(TChain *subruns);

void IntegratedFlux(const TString productionVersion)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/integrated_flux";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = "/exp/sbnd/data/users/hlay/ncpizero/production/NCPiZeroAv7/tmpflux/flux_hist.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *nus = new TChain("fluxana/tree");
  nus->Add(file);

  TChain *pot = new TChain("fluxana/pottree");
  pot->Add(file);

  const double nuPOT = GetPOT(pot);
  const double scaling = goalPOT / nuPOT;
  const double fv_face_area = 175 * 360 * 2.;

  TString potString = Form(" %g POT", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const int N = nus->GetEntries();

  float nu_x, nu_y, nu_e;
  std::vector<float> *weights = 0;

  nus->SetBranchAddress("nu_x", &nu_x);
  nus->SetBranchAddress("nu_y", &nu_y);
  nus->SetBranchAddress("nu_e", &nu_e);
  nus->SetBranchAddress("evtwgt_flux_oneweight", &weights);

  const double true_nu_e_bins[15] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1., 1.2,
                                      1.5, 2., 3., 5. };

  std::vector<TH1F*> hNuEnergyFluxUniverses;
  std::vector<float> counts(1000, 0.);

  TH1F *hNuEnergyNominal = new TH1F("hNuEnergyNominal", ";True E_{#nu} (GeV);#nu s", 14, true_nu_e_bins);
  float nominalCount = 0.;

  for(int j = 0; j < 1000; ++j)
    hNuEnergyFluxUniverses.push_back(new TH1F(Form("hNuEnergyFluxUniverses%d", j), ";True E_{#nu} (GeV);#nu s", 14, true_nu_e_bins));

  for(int i = 0; i < N; ++i)
    {
      nus->GetEntry(i);

      if(abs(nu_x) > 180 || abs(nu_x) < 5 || abs(nu_y) > 180)
        continue;

      hNuEnergyNominal->Fill(nu_e);
      nominalCount += 1.;
      
      for(int j = 0; j < 1000; ++j)
        {
          counts[j] += weights->at(j);
          hNuEnergyFluxUniverses[j]->Fill(nu_e, weights->at(j));
        }
    }

  TH1F *hFluxUniverses = new TH1F("hFluxUniverses", ";Integrated #nu Flux (cm^{-2});Universes", 26, 1.2e13, 2.5e13);
  const float nominalFlux  = nominalCount * scaling / fv_face_area;

  for(int j = 0; j < 1000; ++j)
    {
      float flux = counts[j] * scaling / fv_face_area;

      hFluxUniverses->Fill(flux);
    }

  TCanvas *cFluxUniverses = new TCanvas("cFluxUniverses", "cFluxUniverses");
  cFluxUniverses->cd();
  cFluxUniverses->SetRightMargin(0.15);

  hFluxUniverses->SetLineColor(kBlue-3);
  hFluxUniverses->Draw("histe");

  TLine *nominalLine = new TLine();
  nominalLine->SetLineColor(kMagenta+2);
  nominalLine->SetLineWidth(5);
  nominalLine->DrawLine(nominalFlux, 0., nominalFlux, 1.12 * hFluxUniverses->GetMaximum());

  TLatex *potLatex = new TLatex(hFluxUniverses->GetBinLowEdge(22), 1.14 * hFluxUniverses->GetMaximum(), potString);
  potLatex->SetTextColor(kGray+2);
  potLatex->SetTextSize(0.035);
  potLatex->Draw();

  TString nominalString = Form("#splitline{Nominal:}{%.2ecm^{-2}}", nominalFlux);
  nominalString.ReplaceAll("e+","x10^{");
  nominalString.ReplaceAll("cm^{-2}","}cm^{-2}");

  TLatex *nominalLatex = new TLatex(0.8 * nominalFlux, 0.8 * hFluxUniverses->GetMaximum(), nominalString);
  nominalLatex->SetTextColor(kGray+2);
  nominalLatex->SetTextSize(0.035);
  nominalLatex->SetTextColor(kMagenta+2);
  nominalLatex->Draw();

  cFluxUniverses->SaveAs(saveDir + "/integrated_flux_universes.png");
  cFluxUniverses->SaveAs(saveDir + "/integrated_flux_universes.pdf");

  TCanvas *cNuEnergyFluxUniverses = new TCanvas("cNuEnergyFluxUniverses", "cNuEnergyFluxUniverses");
  cNuEnergyFluxUniverses->cd();
  cNuEnergyFluxUniverses->SetTopMargin(.1);

  NormaliseEntriesByBinWidth(hNuEnergyNominal, .1);
  hNuEnergyNominal->Scale(scaling);
  hNuEnergyNominal->SetMaximum(1.4 * hNuEnergyNominal->GetMaximum());
  hNuEnergyNominal->Draw("hist");

  for(int j = 0; j < 1000; ++j)
    {
      NormaliseEntriesByBinWidth(hNuEnergyFluxUniverses[j], .1);
      hNuEnergyFluxUniverses[j]->Scale(scaling);
      hNuEnergyFluxUniverses[j]->SetLineColor(kMagenta-10);
      hNuEnergyFluxUniverses[j]->SetLineWidth(1);
      hNuEnergyFluxUniverses[j]->Draw("histsame");
    }

  hNuEnergyNominal->Draw("histsame");

  TLegend *lNuEnergyFluxUniverses = new TLegend(.5, .65, .7, .8);
  lNuEnergyFluxUniverses->AddEntry(hNuEnergyNominal, "Nominal", "l");
  lNuEnergyFluxUniverses->AddEntry(hNuEnergyFluxUniverses[0], "Universes", "l");
  lNuEnergyFluxUniverses->Draw();

  cNuEnergyFluxUniverses->SaveAs(saveDir + "/nu_energy_flux_universes.png");
  cNuEnergyFluxUniverses->SaveAs(saveDir + "/nu_energy_flux_universes.pdf");
}

double GetPOT(TChain *subruns)
{
  double sum = 0., pot = 0;

  subruns->SetBranchAddress("pot", &pot);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += pot;
    }

  return sum;
}

