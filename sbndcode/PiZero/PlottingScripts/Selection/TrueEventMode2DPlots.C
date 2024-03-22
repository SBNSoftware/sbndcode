#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void TrueEventMode2DPlots(const TString productionVersion, const std::vector<TwoDPlotSet> &plotSets, const std::vector<Cut> &signals)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_event_modes_two_d";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  TChain *ncpizeroSubruns = new TChain("ncpizeroana/subruns");
  ncpizeroSubruns->Add(ncpizeroFile);

  const double ncpizeroPOT     = GetPOT(ncpizeroSubruns);
  const double ncpizeroScaling = goalPOT / ncpizeroPOT;

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
                                              plotSet.normalisationUnit1.Data())));

          TLegend *lEventModes = new TLegend(.55, .28, .9, .8);

          for(auto const& category : event_modes)
            {
              for(int i = 0; i < plotSet.nbins2; ++i)
                {
                  TH1F *hTemp = new TH1F(Form("h%s%s%i%s", plotSet.name.Data(), signal.name.Data(), i, category.name.Data()), "", plotSet.nbins1, bins1);

                  ncpizeroEvents->Draw(Form("%s>>h%s%s%i%s", plotSet.var1.Data(), plotSet.name.Data(), signal.name.Data(), i, category.name.Data()),
                                       Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + category.cut + signal.cut);

                  hTemp->Scale(ncpizeroScaling);
                  NormaliseEntriesByBinWidth(hTemp, plotSet.scale1);
                  hTemp->SetFillColorAlpha(category.colour, 0.4);
                  hTemp->SetLineColor(category.colour);
                  hTemp->SetLineWidth(2);
                  hTemp->SetMarkerStyle(1);

                  if(i == 0)
                    lEventModes->AddEntry(hTemp, category.name, "f");

                  stacks[i]->Add(hTemp);
                }
            }

          gStyle->SetPaperSize(45,80);

          TCanvas *canvas = new TCanvas("canvas","canvas");
          canvas->cd();
          canvas->Draw();

          const TString wip = "SBND Work-in-progress";
          AddText(canvas, wip, kGray+2, {0., .97, .2, .99}, 0.025, 12);
          const TString potString = POTString(false);
          AddText(canvas, potString, kGray+2, {.8, .97, 1., .99}, 0.025, 32);

          for(int i = 0; i < plotSet.nbins2; ++i)
            {
              canvas->cd();
              int x = i % 3;
              int y = i / 3;

              TPad *pad = new TPad(Form("pad%i", i), "",
                                   0.01 + x * 0.33,
                                   0.97 - (y+1) * 0.32,
                                   (x+1) * 0.33,
                                   0.96 - y * 0.32);

              pad->Draw();
              pad->cd();
              pad->SetTopMargin(0.15);

              stacks[i]->SetMaximum(1.2 * stacks[i]->GetMaximum());
              stacks[i]->Draw("hist");
              stacks[i]->GetYaxis()->SetTitleOffset(1.4);
              stacks[i]->GetYaxis()->SetNdivisions(405);

              if(i==0)
                lEventModes->Draw();
            }

          canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + signal.name + "_by_event_mode.png");
          canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + signal.name + "_by_event_mode.pdf");

          delete canvas;
        }
    }
}
