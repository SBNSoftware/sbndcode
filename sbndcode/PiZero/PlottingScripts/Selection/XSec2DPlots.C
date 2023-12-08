#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Categories.h"
#include "Plots.h"
#include "Cuts.h"

const double goalPOT = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void XSec2DPlots(const TString productionVersion, const TwoDPlotSet &plotSet, const std::vector<Cut> &categories,
                 const Cut &cut, const TCut trueCategory);

void XSec2DPlots(const TString productionVersion, const TwoDPlotSet &plotSet, const std::vector<Cut> &categories,
                 const std::vector<Cut> cuts, const TCut trueCategory, const TString name)
{
  TCut totalCut = "";

  for(auto const& cut : cuts)
    totalCut += cut.cut;

  Cut cut = { name, totalCut, "Full Selection" };

  XSec2DPlots(productionVersion, plotSet, categories, cut, trueCategory);
}

void XSec2DPlots(const TString productionVersion, const TwoDPlotSet &plotSet, const std::vector<Cut> &categories,
                 const Cut &cut, const TCut trueCategory)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/xsec_twod";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxevents = new TChain("ncpizeroana/events");
  rockboxevents->Add(rockboxFile);
  TChain *intimeevents = new TChain("ncpizeroana/events");
  intimeevents->Add(intimeFile);

  TChain *rockboxsubruns = new TChain("ncpizeroana/subruns");
  rockboxsubruns->Add(rockboxFile);
  TChain *intimesubruns = new TChain("ncpizeroana/subruns");
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

  const double nTargets = 4.36e31;
  const double intFlux  = 1.73e13;
  const double nonBinnedFactor = 1 / (nTargets * intFlux);

  double bins1[plotSet.nbins1 + 1];
  std::copy(plotSet.bins1.begin(), plotSet.bins1.end(), bins1);

  TString trueUnSelVar1 = plotSet.var1;
  TString trueUnSelVar2 = plotSet.var2;
  trueUnSelVar1.ReplaceAll("slc_best_pzc", "nu_pz");
  trueUnSelVar2.ReplaceAll("slc_best_pzc", "nu_pz");

  TString trueSelVar1 = plotSet.var1;
  TString trueSelVar2 = plotSet.var2;
  trueSelVar1.ReplaceAll("slc_best_pzc", "slc_true_pz");
  trueSelVar2.ReplaceAll("slc_best_pzc", "slc_true_pz");

  std::vector<TH1*> trueHistsUnselected;
  std::vector<TH1*> trueHistsSelected;

  std::vector<TH1*> recoHistsSig;
  std::vector<TH1*> recoHistsBack;

  for(int i = 0; i < plotSet.nbins2; ++i)
    {
      TH1D *tmpRockboxTrueUnSel = new TH1D(Form("tmpRockboxTrueUnSel%i", i),
                                           Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
                                                plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
                                                plotSet.normalisationUnit.Data()),
                                           plotSet.nbins1, bins1);

      rockboxevents->Draw(Form("%s*1e3>>tmpRockboxTrueUnSel%i", trueUnSelVar1.Data(), i),
                          Form("%s>%f && %s<%f", trueUnSelVar2.Data(), plotSet.bins2[i], trueUnSelVar2.Data(), plotSet.bins2[i+1]) + trueCategory);
      
      trueHistsUnselected.push_back(tmpRockboxTrueUnSel);

      TH1D *tmpRockboxTrueSel = new TH1D(Form("tmpRockboxTrueSel%i", i),
                                         Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
                                              plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
                                              plotSet.normalisationUnit.Data()),
                                         plotSet.nbins1, bins1);

      rockboxevents->Draw(Form("%s*1e3>>tmpRockboxTrueSel%i", trueSelVar1.Data(), i),
                          Form("%s>%f && %s<%f", trueSelVar2.Data(), plotSet.bins2[i], trueSelVar2.Data(), plotSet.bins2[i+1]) + trueCategory + cut.cut);
      
      trueHistsSelected.push_back(tmpRockboxTrueSel);

      TH1D *tmpRockboxRecoSig = new TH1D(Form("tmpRockboxRecoSig%i", i),
                                         Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
                                              plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
                                              plotSet.normalisationUnit.Data()),
                                         plotSet.nbins1, bins1);

      rockboxevents->Draw(Form("%s>>tmpRockboxRecoSig%i", plotSet.var1.Data(), i),
                          Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + categories[0].cut + cut.cut);
      
      tmpRockboxRecoSig->Scale(rockboxScaling);
      NormaliseEntriesByBinWidth(tmpRockboxRecoSig, plotSet.scale);
      recoHistsSig.push_back(tmpRockboxRecoSig);

      TH1D *tmpRockboxRecoBack = new TH1D(Form("tmpRockboxRecoBack%i", i),
                                          Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
                                               plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
                                               plotSet.normalisationUnit.Data()),
                                          plotSet.nbins1, bins1);
      TH1D *tmpIntimeRecoBack = new TH1D(Form("tmpIntimeRecoBack%i", i),
                                         Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
                                              plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
                                              plotSet.normalisationUnit.Data()),
                                         plotSet.nbins1, bins1);

      rockboxevents->Draw(Form("%s>>tmpRockboxRecoBack%i", plotSet.var1.Data(), i),
                          Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + !categories[0].cut + cut.cut);
      intimeevents->Draw(Form("%s>>tmpIntimeRecoBack%i", plotSet.var1.Data(), i),
                         Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + !categories[0].cut + cut.cut);
      
      tmpRockboxRecoBack->Scale(rockboxScaling);
      NormaliseEntriesByBinWidth(tmpRockboxRecoBack, plotSet.scale);
      NormaliseEntriesByBinWidth(tmpIntimeRecoBack, plotSet.scale);
      tmpIntimeRecoBack->Scale(intimeScaling);

      TH1D *hRecoBackClone = (TH1D*) tmpRockboxRecoBack->Clone();
      hRecoBackClone->Add(tmpIntimeRecoBack);

      recoHistsBack.push_back(hRecoBackClone);

      for(int j = 1; j < plotSet.nbins1 + 1; ++j)
        {
          const double eff_ij = trueHistsUnselected[i]->GetBinContent(j) == 0 ? 0. : trueHistsSelected[i]->GetBinContent(j) / trueHistsUnselected[i]->GetBinContent(j);

          const double binWidthI = plotSet.bins2[i+1] - plotSet.bins2[i];
          //const double binWidthJ = bins1[j] - bins1[j-1];      // Already corrected for this when we ran the NormaliseEntriesByBinWidth function

          const double multiplicativeFactor = eff_ij == 0 ? 0. : nonBinnedFactor / (eff_ij * binWidthI);

          recoHistsSig[i]->SetBinContent(j, multiplicativeFactor * recoHistsSig[i]->GetBinContent(j));
        }
    }

  TCanvas *canvas = new TCanvas("canvas","canvas");
  canvas->cd();
  canvas->Draw();
  canvas->Divide(3, 3); // note this hardcodes for 9 bins or less for var2

  for(int i = 0; i < plotSet.nbins2; ++i)
    {
      canvas->cd(i+1);
      gPad->SetTopMargin(0.18);
      recoHistsSig[i]->SetLineColor(kMagenta+2);
      recoHistsSig[i]->SetMaximum(1.2 * recoHistsSig[i]->GetMaximum());
      recoHistsSig[i]->Draw("hist");
      recoHistsSig[i]->GetYaxis()->SetTitleOffset(1.32);
      recoHistsSig[i]->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{#pi^{0}}dcos(#theta_{#pi^{0}})} (#frac{cm^{2}}{MeV/c})");
    }

  canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + cut.name + "_xsec.png");
  canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + cut.name + "_xsec.pdf");
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
