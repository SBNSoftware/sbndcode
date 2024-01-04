#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void Observables2DPlots(const TString productionVersion, const std::vector<TwoDPlotSet> &plotSets, const SelectionParams &selectionParams)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_twod";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile  = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);
  TChain *intimeSubruns = new TChain("ncpizeroana/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  Cut cut = TotalCut(selectionParams.cuts);

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

      TLegend *lSelCategories = new TLegend(.55, .28, .9, .8);

      for(auto const& category : selectionParams.categories)
        {
          for(int i = 0; i < plotSet.nbins2; ++i)
            {
              TH1F *hTemp = new TH1F(Form("h%s%s%i%s", plotSet.name.Data(), cut.name.Data(), i, category.name.Data()), "", plotSet.nbins1, bins1);

              rockboxEvents->Draw(Form("%s>>h%s%s%i%s", plotSet.var1.Data(), plotSet.name.Data(), cut.name.Data(), i, category.name.Data()),
                                  Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + category.cut + cut.cut);

              TH1F *hTempIntime = new TH1F(Form("hIntime%s%s%i%s", plotSet.name.Data(), cut.name.Data(), i, category.name.Data()), "", plotSet.nbins1, bins1);

              intimeEvents->Draw(Form("%s>>hIntime%s%s%i%s", plotSet.var1.Data(), plotSet.name.Data(), cut.name.Data(), i, category.name.Data()),
                                 Form("%s>%f && %s<%f", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1]) + category.cut + cut.cut);

              hTemp->Scale(rockboxScaling);
              hTempIntime->Scale(intimeScaling);
              hTemp->Add(hTempIntime);
              NormaliseEntriesByBinWidth(hTemp, plotSet.scale);
              hTemp->SetFillColorAlpha(category.colour, 0.4);
              hTemp->SetLineColor(category.colour);
              hTemp->SetLineWidth(2);
              hTemp->SetMarkerStyle(1);

              if(i == 0)
                lSelCategories->AddEntry(hTemp, category.printed_name, "f");

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
          stacks[i]->SetMaximum(1.2 * stacks[i]->GetMaximum());
          stacks[i]->Draw("histe");
          stacks[i]->GetYaxis()->SetTitleOffset(1.3);
          if(i==0)
            lSelCategories->Draw();
        }

      canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + selectionParams.name + "_by_event_mode.png");
      canvas->SaveAs(saveDir + "/" + plotSet.name + "_" + selectionParams.name + "_by_event_mode.pdf");

      delete canvas;
    }
}
