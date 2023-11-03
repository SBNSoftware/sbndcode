#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Categories.h"

const double goalPOT     = 10e20;

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void TrueEventMode2DPlots(const TString productionVersion)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/true_event_modes_two_d";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxevents = new TChain("ncpizeroana/events");
  rockboxevents->Add(rockboxFile);

  TChain *rockboxsubruns = new TChain("ncpizeroana/subruns");
  rockboxsubruns->Add(rockboxFile);

  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const double rockboxPOT = GetPOT(rockboxsubruns);
  const double rockboxScaling = goalPOT / rockboxPOT;

  const double pizero_momentum_bins[9] = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double pizero_direction_bins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  std::vector<THStack*> sPiZeroMomentumInDirectionBins;
  for(int i = 0; i < 9; ++i)
    sPiZeroMomentumInDirectionBins.push_back(new THStack(Form("sPiZeroMomentumDirectionBins%i", i),
                                                         Form("%.2f < cos(#theta_{#pi^{0}} ) < %.2f;p_{#pi^{0}} (MeV/c);Slices / MeV/c", pizero_direction_bins[i], pizero_direction_bins[i+1])));
    
  TLegend *lEventModes = new TLegend(.55, .28, .9, .8);

  for(auto const& category : event_modes)
    {
      for(int i = 0; i < 9; ++i)
        {
          TH1F *hTemp = new TH1F(Form("hPiZeroMomentumDirectionBin%iCategory%s", i, category.name.Data()), "", 8, pizero_momentum_bins);

          rockboxevents->Draw(Form("1e3 * nu_pz_pizero_mom>>hPiZeroMomentumDirectionBin%iCategory%s", i, category.name.Data()),
                              Form("nu_pz_cos_theta_pizero>%f && nu_pz_cos_theta_pizero<%f && nu_signal", pizero_direction_bins[i], pizero_direction_bins[i+1]) + category.cut);

          hTemp->Scale(rockboxScaling);
          NormaliseEntriesByBinWidth(hTemp);
          hTemp->SetFillColorAlpha(category.colour, 0.4);
          hTemp->SetLineColor(category.colour);
          hTemp->SetLineWidth(2);

          if(i == 0)
            lEventModes->AddEntry(hTemp, category.name, "f");

          sPiZeroMomentumInDirectionBins[i]->Add(hTemp);
        }
    }

  gStyle->SetPaperSize(45,80);

  TCanvas *cPiZeroMomentumInDirectionBins = new TCanvas("cPiZeroMomentumInDirectionBins","cPiZeroMomentumInDirectionBins");
  cPiZeroMomentumInDirectionBins->cd();
  cPiZeroMomentumInDirectionBins->Draw();
  cPiZeroMomentumInDirectionBins->Divide(3, 3);

  for(int i = 0; i < 9; ++i)
    {
      cPiZeroMomentumInDirectionBins->cd(i+1);
      gPad->SetTopMargin(0.15);
      sPiZeroMomentumInDirectionBins[i]->Draw("hist");
      if(i==0)
        lEventModes->Draw();
    }

  cPiZeroMomentumInDirectionBins->SaveAs(saveDir + "/pizero_momentum_and_cos_theta_by_event_mode.png");
  cPiZeroMomentumInDirectionBins->SaveAs(saveDir + "/pizero_momentum_and_cos_theta_by_event_mode.pdf");
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
