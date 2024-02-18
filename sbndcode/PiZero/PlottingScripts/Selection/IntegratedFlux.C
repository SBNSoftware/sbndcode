#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"
#include "WeightNames.h"

void IntegratedFlux(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/integrated_flux";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_flux_configI.root";

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

  std::vector<std::string> weight_names = flux_weight_names;
  const int n_weights = weight_names.size() + 1;

  float nu_x, nu_y, nu_e;
  std::vector<std::vector<float>*> weights = std::vector<std::vector<float>*>(n_weights, 0);

  nus->SetBranchAddress("nu_x", &nu_x);
  nus->SetBranchAddress("nu_y", &nu_y);
  nus->SetBranchAddress("nu_e", &nu_e);

  for(auto&& [ weight_i, name ] : enumerate(weight_names))
    nus->SetBranchAddress(Form("evtwgt_flux_weight_%s", name.c_str()), &weights[weight_i]);

  nus->SetBranchAddress("evtwgt_flux_oneweight", &weights[weight_names.size()]);

  weight_names.push_back("all");

  const double true_nu_e_bins[15] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1., 1.2,
                                      1.5, 2., 3., 5. };

  std::vector<std::vector<TH1F*>> hNuEnergyFluxUniverses = std::vector<std::vector<TH1F*>>(n_weights, std::vector<TH1F*>());
  std::vector<std::vector<float>> counts(n_weights, std::vector<float>(1000, 0.));

  TH1F *hNuEnergyNominal = new TH1F("hNuEnergyNominal", ";True E_{#nu} (GeV);#nu", 14, true_nu_e_bins);
  float nominalCount = 0.;

  ofstream outFile;
  outFile.open(saveDir + "/FluxMap.h");
  outFile << "const std::map<std::string, std::map<int, float>> univsIntegratedFluxMap = {" << std::endl;

  for(auto&& [ weight_i, name ] : enumerate(weight_names))
    {
      for(int j = 0; j < 1000; ++j)
        hNuEnergyFluxUniverses[weight_i].push_back(new TH1F(Form("hNuEnergyFluxUniverses%s%d", name.c_str(), j), ";True E_{#nu} (GeV);#nu", 14, true_nu_e_bins));
    }

  for(int i = 0; i < N; ++i)
    {
      nus->GetEntry(i);

      if(abs(nu_x) > 180 || abs(nu_x) < 5 || abs(nu_y) > 180)
        continue;

      hNuEnergyNominal->Fill(nu_e);
      nominalCount += 1.;

      for(auto&& [ weight_i, name ] : enumerate(weight_names))
        {
          for(int j = 0; j < 1000; ++j)
            {
              counts[weight_i][j] += weights[weight_i]->at(j);
              hNuEnergyFluxUniverses[weight_i][j]->Fill(nu_e, weights[weight_i]->at(j));
            }
        }
    }

  for(auto&& [ weight_i, name ] : enumerate(weight_names))
    {
      gSystem->Exec("mkdir -p " + saveDir + Form("/%s", name.c_str()));

      TH1F *hFluxUniverses = new TH1F(Form("hFluxUniverses%s", name.c_str()), ";Integrated #nu Flux (cm^{-2});Universes", 26, 1.2e13, 2.5e13);
      const float nominalFlux  = nominalCount * scaling / fv_face_area;

      outFile << Form("{ \"%s\", {", name.c_str()) << std::endl;

      for(int j = 0; j < 1000; ++j)
        {
          float flux = counts[weight_i][j] * scaling / fv_face_area;

          hFluxUniverses->Fill(flux);

          outFile << Form("{ %i, %f },", j, flux) << std::endl;
        }

      outFile << "}," << std::endl;

      TCanvas *cFluxUniverses = new TCanvas(Form("cFluxUniverses%s", name.c_str()), Form("cFluxUniverses%s", name.c_str()));
      cFluxUniverses->cd();
      cFluxUniverses->SetRightMargin(0.15);

      hFluxUniverses->SetLineColor(kBlue-3);
      hFluxUniverses->Draw("histe");

      TLine *nominalLine = new TLine();
      nominalLine->SetLineColor(kMagenta+2);
      nominalLine->SetLineWidth(5);
      nominalLine->DrawLine(nominalFlux, 0., nominalFlux, 1.12 * hFluxUniverses->GetMaximum());

      TLatex *potLatex = new TLatex(hFluxUniverses->GetBinLowEdge(22), 1.14 * hFluxUniverses->GetMaximum(), POTString());
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

      cFluxUniverses->SaveAs(saveDir + Form("/%s/integrated_flux_universes_%s.png", name.c_str(), name.c_str()));
      cFluxUniverses->SaveAs(saveDir + Form("/%s/integrated_flux_universes_%s.pdf", name.c_str(), name.c_str()));

      TCanvas *cNuEnergyFluxUniverses = new TCanvas(Form("cNuEnergyFluxUniverses%s", name.c_str()), Form("cNuEnergyFluxUniverses%s", name.c_str()));
      cNuEnergyFluxUniverses->cd();
      cNuEnergyFluxUniverses->SetTopMargin(.1);

      NormaliseEntriesByBinWidth(hNuEnergyNominal, .1);
      hNuEnergyNominal->Scale(scaling);
      hNuEnergyNominal->SetMaximum(1.5 * hNuEnergyNominal->GetMaximum());
      hNuEnergyNominal->Draw("hist");

      for(int j = 0; j < 1000; ++j)
        {
          NormaliseEntriesByBinWidth(hNuEnergyFluxUniverses[weight_i][j], .1);
          hNuEnergyFluxUniverses[weight_i][j]->Scale(scaling);
          hNuEnergyFluxUniverses[weight_i][j]->SetLineColor(kMagenta-10);
          hNuEnergyFluxUniverses[weight_i][j]->SetLineWidth(1);
          hNuEnergyFluxUniverses[weight_i][j]->Draw("histsame");
        }

      hNuEnergyNominal->Draw("histsame");

      TLegend *lNuEnergyFluxUniverses = new TLegend(.5, .65, .7, .8);
      lNuEnergyFluxUniverses->AddEntry(hNuEnergyNominal, "Nominal", "l");
      lNuEnergyFluxUniverses->AddEntry(hNuEnergyFluxUniverses[weight_i][0], "Universes", "l");
      lNuEnergyFluxUniverses->Draw();

      cNuEnergyFluxUniverses->SaveAs(saveDir + Form("/%s/nu_energy_flux_universes_%s.png", name.c_str(), name.c_str()));
      cNuEnergyFluxUniverses->SaveAs(saveDir + Form("/%s/nu_energy_flux_universes_%s.pdf", name.c_str(), name.c_str()));

      std::vector<TH1F*> hNuEnergyFluxUniversesPerBin;
      TH1F* hNuEnergyFluxUniversesCV = new TH1F(Form("hNuEnergyFluxUniversesCV%s", name.c_str()), ";True E_{#nu} (GeV);#nu", 14, true_nu_e_bins);

      for(int bin = 0; bin < hNuEnergyNominal->GetNbinsX(); ++bin)
        {
          hNuEnergyFluxUniversesPerBin.push_back(new TH1F(Form("hNuEnergyFluxUniverses%sBin%i", name.c_str(), bin+1),
                                                          Form(";%.1f < True E_{#nu} (GeV) < %.1f;Universes", hNuEnergyNominal->GetBinLowEdge(bin+1),
                                                               hNuEnergyNominal->GetBinLowEdge(bin+1) + hNuEnergyNominal->GetBinWidth(bin+1)),
                                                          25, .5 * hNuEnergyNominal->GetBinContent(bin+1), 1.5 * hNuEnergyNominal->GetBinContent(bin+1)));

          for(int j = 0; j < 1000; ++j)
            hNuEnergyFluxUniversesPerBin[bin]->Fill(hNuEnergyFluxUniverses[weight_i][j]->GetBinContent(bin+1));

          TCanvas *cNuEnergyFluxUniversesPerBin = new TCanvas(Form("cNuEnergyFluxUniverses%sBin%i", name.c_str(), bin+1),
                                                              Form("cNuEnergyFluxUniverses%sBin%i", name.c_str(), bin+1));
          cNuEnergyFluxUniversesPerBin->cd();

          hNuEnergyFluxUniversesPerBin[bin]->SetLineColor(kBlue+2);

          TF1 *fGaus = new TF1("fGaus", "gaus", .5 * hNuEnergyNominal->GetBinContent(bin+1), 1.5 * hNuEnergyNominal->GetBinContent(bin+1));
          hNuEnergyFluxUniversesPerBin[bin]->Fit(fGaus);
          fGaus->SetLineColor(kRed+2);

          hNuEnergyFluxUniversesPerBin[bin]->Draw("hist");
          fGaus->Draw("same");

          cNuEnergyFluxUniversesPerBin->SaveAs(saveDir + Form("/%s/nu_energy_bin%i_flux_variations_%s.png", name.c_str(), bin+1, name.c_str()));
          cNuEnergyFluxUniversesPerBin->SaveAs(saveDir + Form("/%s/nu_energy_bin%i_flux_variations_%s.pdf", name.c_str(), bin+1, name.c_str()));

          hNuEnergyFluxUniversesCV->SetBinContent(bin+1, fGaus->GetParameter("Mean"));
          hNuEnergyFluxUniversesCV->SetBinError(bin+1, fGaus->GetParameter("Sigma"));
        }

      TCanvas *cNuEnergyFluxUniversesOneSigma = new TCanvas(Form("cNuEnergyFluxUniversesOneSigma%s", name.c_str()),
                                                            Form("cNuEnergyFluxUniversesOneSigma%s", name.c_str()));
      cNuEnergyFluxUniversesOneSigma->cd();
      cNuEnergyFluxUniversesOneSigma->SetTopMargin(.1);

      hNuEnergyNominal->Draw("hist");
      hNuEnergyFluxUniversesCV->SetLineColor(kViolet-6);
      hNuEnergyFluxUniversesCV->Draw("histesame");

      TLegend *lNuEnergyFluxUniversesOneSigma = new TLegend(.5, .65, .7, .8);
      lNuEnergyFluxUniversesOneSigma->AddEntry(hNuEnergyNominal, "Nominal", "l");
      lNuEnergyFluxUniversesOneSigma->AddEntry(hNuEnergyFluxUniversesCV, "CV #pm 1 #sigma", "l");
      lNuEnergyFluxUniversesOneSigma->Draw();

      potLatex->Draw();

      cNuEnergyFluxUniversesOneSigma->SaveAs(saveDir + Form("/%s/nu_energy_flux_universes_one_sigma_%s.png", name.c_str(), name.c_str()));
      cNuEnergyFluxUniversesOneSigma->SaveAs(saveDir + Form("/%s/nu_energy_flux_universes_one_sigma_%s.pdf", name.c_str(), name.c_str()));

      TH2F *hNuEnergyFluxUniversesCovariance = new TH2F(Form("hNuEnergyFluxUniversesCovariance%s", name.c_str()), ";Reco Bin Index;Reco Bin Index",
                                                        14, 0.5, 14.5, 14, 0.5, 14.5);

      TH2F *hNuEnergyFluxUniversesCorrelation = new TH2F(Form("hNuEnergyFluxUniversesCorrelation%s", name.c_str()), ";Reco Bin Index;Reco Bin Index",
                                                         14, 0.5, 14.5, 14, 0.5, 14.5);

      for(int binX = 1; binX <= hNuEnergyFluxUniversesCovariance->GetNbinsX(); ++binX)
        {
          for(int binY = 1; binY <= hNuEnergyFluxUniversesCovariance->GetNbinsY(); ++binY)
            {
              double covXY = 0.;
              double covXX = 0.;
              double covYY = 0.;

              for(int univ = 0; univ < 1000; ++univ)
                {
                  covXY += (hNuEnergyFluxUniverses[weight_i][univ]->GetBinContent(binX) - hNuEnergyNominal->GetBinContent(binX))
                    * (hNuEnergyFluxUniverses[weight_i][univ]->GetBinContent(binY) - hNuEnergyNominal->GetBinContent(binY));

                  covXX += (hNuEnergyFluxUniverses[weight_i][univ]->GetBinContent(binX) - hNuEnergyNominal->GetBinContent(binX))
                    * (hNuEnergyFluxUniverses[weight_i][univ]->GetBinContent(binX) - hNuEnergyNominal->GetBinContent(binX));

                  covYY += (hNuEnergyFluxUniverses[weight_i][univ]->GetBinContent(binY) - hNuEnergyNominal->GetBinContent(binY))
                    * (hNuEnergyFluxUniverses[weight_i][univ]->GetBinContent(binY) - hNuEnergyNominal->GetBinContent(binY));
                }

              covXY /= 1000;
              covXX /= 1000;
              covYY /= 1000;

              const double corrXY = covXY / (TMath::Sqrt(covXX) * TMath::Sqrt(covYY));

              hNuEnergyFluxUniversesCovariance->SetBinContent(binX, binY, covXY);
              hNuEnergyFluxUniversesCorrelation->SetBinContent(binX, binY, corrXY);
            }
        }

      TCanvas *cNuEnergyFluxUniversesCovariance = new TCanvas(Form("cNuEnergyFluxUniversesCovariance%s", name.c_str()),
                                                              Form("cNuEnergyFluxUniversesCovariance%s", name.c_str()));
      cNuEnergyFluxUniversesCovariance->cd();
      cNuEnergyFluxUniversesCovariance->SetRightMargin(.2);

      hNuEnergyFluxUniversesCovariance->Draw("colz");

      cNuEnergyFluxUniversesCovariance->SaveAs(saveDir + Form("/%s/nu_energy_flux_covariance_%s.png", name.c_str(), name.c_str()));
      cNuEnergyFluxUniversesCovariance->SaveAs(saveDir + Form("/%s/nu_energy_flux_covariance_%s.pdf", name.c_str(), name.c_str()));

      TCanvas *cNuEnergyFluxUniversesCorrelation = new TCanvas("cNuEnergyFluxUniversesCorrelation", "cNuEnergyFluxUniversesCorrelation");
      cNuEnergyFluxUniversesCorrelation->cd();
      cNuEnergyFluxUniversesCorrelation->SetRightMargin(.2);

      hNuEnergyFluxUniversesCorrelation->Draw("colz");

      cNuEnergyFluxUniversesCorrelation->SaveAs(saveDir + Form("/%s/nu_energy_flux_correlation_%s.png", name.c_str(), name.c_str()));
      cNuEnergyFluxUniversesCorrelation->SaveAs(saveDir + Form("/%s/nu_energy_flux_correlation_%s.pdf", name.c_str(), name.c_str()));
    }

  outFile << "};";
  outFile.close();
}
