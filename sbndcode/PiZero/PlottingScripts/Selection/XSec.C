#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

const int n_flux_univs = 1000;

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void XSec(const TString productionVersion, const TString saveDirExt)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/xsec/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/xsec_trees/" + productionVersion + "_rockbox_xsec_trees.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/xsec_trees/" + productionVersion + "_intime_xsec_trees.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxslices = new TChain("slices");
  rockboxslices->Add(rockboxFile);
  TChain *intimeslices = new TChain("slices");
  intimeslices->Add(intimeFile);

  TChain *rockboxsubruns = new TChain("subruns");
  rockboxsubruns->Add(rockboxFile);
  TChain *intimesubruns = new TChain("subruns");
  intimesubruns->Add(intimeFile);

  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const double rockboxPOT = GetPOT(rockboxsubruns);
  const int rockboxSpills = GetGenEvents(rockboxsubruns);
  const int intimeSpills  = GetGenEvents(intimesubruns);

  const double rockboxScaling      = goalPOT / rockboxPOT;
  const double scaledRockboxSpills = rockboxScaling * rockboxSpills;
  const double intimeScaling       = (goalSpills - scaledRockboxSpills) / intimeSpills;

  const double pizero_momentum_bins[15] = { 0., 40., 80., 120., 160., 200., 240., 280., 320., 360., 400.,
                                            480., 560., 700., 1000. };

  TH1F *hNominalPiZeroMomentumRockbox = new TH1F("hNominalPiZeroMomentumRockbox", ";p_{#pi^{0}} (MeV/c);Slices / MeV/c", 14, pizero_momentum_bins);
  TH1F *hNominalPiZeroMomentumIntime = new TH1F("hNominalPiZeroMomentumIntime", ";p_{#pi^{0}} (MeV/c);Slices / MeV/c", 14, pizero_momentum_bins);
  rockboxslices->Draw("pzc_pizero_mom>>hNominalPiZeroMomentumRockbox", "selected");
  intimeslices->Draw("pzc_pizero_mom>>hNominalPiZeroMomentumIntime", "selected");
  hNominalPiZeroMomentumRockbox->Scale(rockboxScaling);
  hNominalPiZeroMomentumIntime->Scale(intimeScaling);
  TH1F* hNominalPiZeroMomentum = (TH1F*) hNominalPiZeroMomentumRockbox->Clone("hNominalPiZeroMomentum");
  hNominalPiZeroMomentum->Add(hNominalPiZeroMomentumIntime);
  NormaliseEntriesByBinWidth(hNominalPiZeroMomentum, .1);

  std::vector<TH1F*> hFluxUniversesPiZeroMomentum;
  for(int univ = 0; univ < n_flux_univs; ++univ)
    {
      TH1F *hFluxUniversesPiZeroMomentumRockbox = new TH1F(Form("hFluxUniverse%iPiZeroMomentumRockbox", univ), ";p_{#pi^{0}} (MeV/c);Slices / MeV/c", 14, pizero_momentum_bins);
      TH1F *hFluxUniversesPiZeroMomentumIntime = new TH1F(Form("hFluxUniverse%iPiZeroMomentumIntime", univ), ";p_{#pi^{0}} (MeV/c);Slices / MeV/c", 14, pizero_momentum_bins);
      rockboxslices->Draw(Form("pzc_pizero_mom>>hFluxUniverse%iPiZeroMomentumRockbox", univ), Form("flux_weights[%i] * selected",univ));
      intimeslices->Draw(Form("pzc_pizero_mom>>hFluxUniverse%iPiZeroMomentumIntime", univ), Form("flux_weights[%i] * selected",univ));
      hFluxUniversesPiZeroMomentumRockbox->Scale(rockboxScaling);
      hFluxUniversesPiZeroMomentumIntime->Scale(intimeScaling);
      TH1F *hFluxUniversesPiZeroMomentumCombined = (TH1F*) hFluxUniversesPiZeroMomentumRockbox->Clone("hFluxUniversesPiZeroMomentum");
      hFluxUniversesPiZeroMomentumCombined->Add(hFluxUniversesPiZeroMomentumIntime);

      hFluxUniversesPiZeroMomentum.push_back(hFluxUniversesPiZeroMomentumCombined);
      hFluxUniversesPiZeroMomentum[univ]->SetLineColor(kMagenta-10);
      hFluxUniversesPiZeroMomentum[univ]->SetLineWidth(1);
      NormaliseEntriesByBinWidth(hFluxUniversesPiZeroMomentum[univ], .1);
    }

  TCanvas *cFluxUniversesPiZeroMomentum = new TCanvas("cFluxUniversesPiZeroMomentum", "cFluxUniversesPiZeroMomentum");
  cFluxUniversesPiZeroMomentum->cd();

  hNominalPiZeroMomentum->SetMaximum(1.4*hNominalPiZeroMomentum->GetMaximum());
  hNominalPiZeroMomentum->Draw("hist");

  for(int univ = 0; univ < n_flux_univs; ++univ)
    hFluxUniversesPiZeroMomentum[univ]->Draw("histsame");
  
  cFluxUniversesPiZeroMomentum->Modified();
  hNominalPiZeroMomentum->Draw("histsame");

  TLegend *lFluxUniversesPiZeroMomentum = new TLegend(.5, .65, .7, .8);
  lFluxUniversesPiZeroMomentum->AddEntry(hNominalPiZeroMomentum, "Nominal", "l");
  lFluxUniversesPiZeroMomentum->AddEntry(hFluxUniversesPiZeroMomentum[0], "Universes", "l");
  lFluxUniversesPiZeroMomentum->Draw();

  cFluxUniversesPiZeroMomentum->SaveAs(saveDir + "/pizero_momentum_flux_universes.png");
  cFluxUniversesPiZeroMomentum->SaveAs(saveDir + "/pizero_momentum_flux_universes.pdf");

  std::vector<TH1F*> hFluxUniversesPiZeroMomentumPerBin;
  TH1F *hFluxCVPiZeroMomentum = new TH1F("hFluxCVPiZeroMomentum", ";p_{#pi^{0}} (MeV/c);Slices / MeV/c", 14, pizero_momentum_bins);

  for(int bin = 0; bin < hNominalPiZeroMomentum->GetNbinsX(); ++bin)
    {
      hFluxUniversesPiZeroMomentumPerBin.push_back(new TH1F(Form("hFluxUniversesPiZeroMomentumBin%i",bin+1),
                                                            Form(";%.0f < p_{#pi^{0}} (MeV/c) < %.0f;Universes", hNominalPiZeroMomentum->GetBinLowEdge(bin+1),
                                                                 hNominalPiZeroMomentum->GetBinLowEdge(bin+1) + hNominalPiZeroMomentum->GetBinWidth(bin+1)),
                                                            25, .75 * hNominalPiZeroMomentum->GetBinContent(bin+1), 1.25 * hNominalPiZeroMomentum->GetBinContent(bin+1)));
      
      for(int univ = 0; univ < n_flux_univs; ++univ)
        hFluxUniversesPiZeroMomentumPerBin[bin]->Fill(hFluxUniversesPiZeroMomentum[univ]->GetBinContent(bin+1));

      TCanvas *cFluxUniversesPiZeroMomentumPerBin = new TCanvas(Form("cFluxUniversesPiZeroMomentumBin%i", bin+1), Form("cFluxUniversesPiZeroMomentumBin%i", bin+1));
      cFluxUniversesPiZeroMomentumPerBin->cd();

      hFluxUniversesPiZeroMomentumPerBin[bin]->SetLineColor(kBlue+2);

      TF1 *fGaus = new TF1("fGaus", "gaus", .75 * hNominalPiZeroMomentum->GetBinContent(bin+1), 1.25 * hNominalPiZeroMomentum->GetBinContent(bin+1));
      hFluxUniversesPiZeroMomentumPerBin[bin]->Fit(fGaus);
      fGaus->SetLineColor(kRed+2);

      hFluxUniversesPiZeroMomentumPerBin[bin]->Draw("hist");
      fGaus->Draw("same");

      cFluxUniversesPiZeroMomentumPerBin->SaveAs(saveDir + Form("/pizero_momentum_bin%i_flux_variations.png", bin+1));
      cFluxUniversesPiZeroMomentumPerBin->SaveAs(saveDir + Form("/pizero_momentum_bin%i_flux_variations.pdf", bin+1));

      hFluxCVPiZeroMomentum->SetBinContent(bin+1, fGaus->GetParameter("Mean"));
      hFluxCVPiZeroMomentum->SetBinError(bin+1, fGaus->GetParameter("Sigma"));
    }

  TCanvas *cFluxOneSigmaPiZeroMomentum = new TCanvas("cFluxOneSigmaPiZeroMomentum", "cFluxOneSigmaPiZeroMomentum");
  cFluxOneSigmaPiZeroMomentum->cd();

  hNominalPiZeroMomentum->Draw("hist");
  hFluxCVPiZeroMomentum->SetLineColor(kViolet-6);
  hFluxCVPiZeroMomentum->Draw("histesame");

  TLegend *lFluxOneSigmaPiZeroMomentum = new TLegend(.5, .65, .7, .8);
  lFluxOneSigmaPiZeroMomentum->AddEntry(hNominalPiZeroMomentum, "Nominal", "l");
  lFluxOneSigmaPiZeroMomentum->AddEntry(hFluxCVPiZeroMomentum, "CV #pm 1 #sigma", "l");
  lFluxOneSigmaPiZeroMomentum->Draw();

  cFluxOneSigmaPiZeroMomentum->SaveAs(saveDir + "/pizero_momentum_flux_one_sigma.png");
  cFluxOneSigmaPiZeroMomentum->SaveAs(saveDir + "/pizero_momentum_flux_one_sigma.pdf");

  TH2F *hFluxCovariancePiZeroMomentum = new TH2F("hFluxCovariancePiZeroMomentum", ";Reco Bin Index;Reco Bin Index", 14, 0.5, 14.5, 14, 0.5, 14.5);

  for(int binX = 1; binX <= hFluxCovariancePiZeroMomentum->GetNbinsX(); ++binX)
    {
      for(int binY = 1; binY <= hFluxCovariancePiZeroMomentum->GetNbinsY(); ++binY)
        {
          double covXY = 0.;

          for(int univ = 0; univ < n_flux_univs; ++univ)
            covXY += (hFluxUniversesPiZeroMomentum[univ]->GetBinContent(binX) - hNominalPiZeroMomentum->GetBinContent(binX))
              * (hFluxUniversesPiZeroMomentum[univ]->GetBinContent(binY) - hNominalPiZeroMomentum->GetBinContent(binY));

          covXY /= n_flux_univs;

          hFluxCovariancePiZeroMomentum->SetBinContent(binX, binY, covXY);
        }
    }

  TCanvas *cFluxConvariancePiZeroMomentum = new TCanvas("cFluxConvariancePiZeroMomentum", "cFluxConvariancePiZeroMomentum");
  cFluxConvariancePiZeroMomentum->cd();
  cFluxConvariancePiZeroMomentum->SetRightMargin(.2);

  hFluxCovariancePiZeroMomentum->Draw("colz");
  
  cFluxConvariancePiZeroMomentum->SaveAs(saveDir + "/pizero_momentum_flux_covariance.png");
  cFluxConvariancePiZeroMomentum->SaveAs(saveDir + "/pizero_momentum_flux_covariance.pdf");
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

int GetGenEvents(TChain *subruns)
{
  int sum = 0., ngenevts = 0;

  subruns->SetBranchAddress("ngenevts", &ngenevts);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += ngenevts;
    }

  return sum;
}
