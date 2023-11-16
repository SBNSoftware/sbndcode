#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Categories.h"
#include "Plots.h"

const double goalPOT = 10e20;

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void TrueEventMode2DPlots(const TString productionVersion, const std::vector<TwoDPlotSet> &plotSets, const std::vector<Cut> &signals)
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

  for(auto const& signal : signals)
    {
      for(auto const& plotSet : plotSets)
        {
          double bins1[plotSet.nbins1 + 1];
          std::copy(plotSet.bins1.begin(), plotSet.bins1.end(), bins1);

          std::vector<THStack*> stacks;
          for(int i = 0; i < plotSet.nbins2; ++i)
            stacks.push_back(new THStack(Form("stack%i", i),
                                         Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
                                              plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
                                              plotSet.normalisationUnit.Data())));

          TLegend *lEventModes = new TLegend(.55, .28, .9, .8);

          for(auto const& category : event_modes)
            {
              for(int i = 0; i < plotSet.nbins2; ++i)
                {
                  TH1F *hTemp = new TH1F(Form("h%s%s%i%s", plotSet.name.Data(), signal.name.Data(), i, category.name.Data()), "", plotSet.nbins1, bins1);

                  rockboxevents->Draw(Form("%s>>h%s%s%i%s", plotSet.var1.Data(), plotSet.name.Data(), signal.name.Data(), i, category.name.Data()),
                                      Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + category.cut + signal.cut);

                  hTemp->Scale(rockboxScaling);
                  NormaliseEntriesByBinWidth(hTemp, plotSet.scale);
                  hTemp->SetFillColorAlpha(category.colour, 0.4);
                  hTemp->SetLineColor(category.colour);
                  hTemp->SetLineWidth(2);

                  if(i == 0)
                    lEventModes->AddEntry(hTemp, category.name, "f");

                  stacks[i]->Add(hTemp);
                }
            }

          gStyle->SetPaperSize(45,80);

          TCanvas *canvas = new TCanvas("canvas","canvas");
          canvas->cd();
          canvas->Draw();
          canvas->Divide(3, 3); // note this hardcodes for 9 bins or less for var2

          for(int i = 0; i < plotSet.nbins2; ++i)
            {
              canvas->cd(i+1);
              gPad->SetTopMargin(0.15);
              stacks[i]->Draw("hist");
              stacks[i]->GetYaxis()->SetTitleOffset(1.3);
              if(i==0)
                lEventModes->Draw();
            }

          canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + signal.name + "_by_event_mode.png");
          canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + signal.name + "_by_event_mode.pdf");

          delete canvas;
        }
    }
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
