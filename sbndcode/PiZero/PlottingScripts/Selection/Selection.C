#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"
#include "Selections.h"
#include "Plots.h"
#include "Common.C"

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples, const SelectionParams &selectionParams);

void Selection(const TString productionVersion, const SelectionParams &selectionParams, std::vector<Plot> &plots)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/selection/" + selectionParams.name;
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

  std::vector<Sample<TChain>> samples = { { "rockbox", rockboxEvents, rockboxScaling },
                                          { "intime", intimeEvents, intimeScaling }
  };

  ProduceCutTable(saveDir, samples, selectionParams);

  TCut currentCut = "";

  for(auto cut : selectionParams.cuts)
    {
      currentCut += cut.cut;
      cut.cut = currentCut;
      
      if(plots.size() != 0)
        gSystem->Exec("mkdir -p " + saveDir + "/" + cut.name);
      for(auto plot : plots)
        {
          TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name,
                                        "c_" + plot.name + "_" + cut.name);
          canvas->cd();

          plot.axes_labels += POTString();

          MakeStackedPlot(canvas, samples, plot, cut, selectionParams.categories, {.25, .8, .8, .87}, 4);

          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");

          delete canvas;
        }
    }

  gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/cut_table.tex");
}

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples, const SelectionParams &selectionParams)
{
  ofstream texFile;
  texFile.open(saveDir + "/cut_table.tex");

  double totalSignal = 0, totalSignalSlices = 0, totalBackSlices = 0;

  for(auto const& sample : samples)
    {
      totalSignal       += sample.scaling * sample.tree->Draw("", selectionParams.true_category);
      totalSignalSlices += sample.scaling * sample.tree->Draw("", selectionParams.categories[0].cut);
      totalBackSlices   += sample.scaling * sample.tree->Draw("", !(selectionParams.categories[0].cut));
    }

  texFile << docStart;

  texFile << '\n'
          << "Total Signal: " << totalSignal << "\\\\ \n"
          << "Total Signal Slices: " << totalSignalSlices << "\\\\ \n"
          << "Total Background Slices: " << totalBackSlices << "\\\\ \n" << std::endl;

  texFile << tableStart 
          << "\\hline\n"
          << "Cut Name & $\\epsilon$ (\\%) & $\\rho$ (\\%) & $\\epsilon\\rho$ & Selection $\\epsilon$ (\\%) & Selection $\\epsilon\\rho$ & BR (\\%) \\\\ \\hline" 
          << std::endl;

  TCut currentCut = "";

  for(unsigned i = 0; i < selectionParams.cuts.size(); ++i)
    {
      Cut cut = selectionParams.cuts[i];
      currentCut += cut.cut;
      cut.cut = currentCut;

      double sigSlices = 0., backSlices = 0.;
      for(unsigned j = 0; j < selectionParams.categories.size(); ++j)
        {
          for(auto const& sample : samples)
            {
              if(j == 0)
                sigSlices += sample.scaling * sample.tree->Draw("", cut.cut + selectionParams.categories[j].cut);
              else
                backSlices += sample.scaling * sample.tree->Draw("", cut.cut + selectionParams.categories[j].cut);
            }
        }

      const double eff     = sigSlices * 100. / totalSignal;
      const double selEff  = sigSlices * 100. / totalSignalSlices;
      const double pur     = sigSlices * 100./ (sigSlices + backSlices);
      const double backRej = 100. - 100. * (backSlices / totalBackSlices);
      
      texFile << cut.printed_name << " & " << Form("%.2f", eff) << " & " << Form("%.2f", pur)
              << " & " << Form("%.2f", (eff * pur) / 100.)
              << " & " << Form("%.2f", selEff) << " & " << Form("%.2f", (selEff * pur) / 100.)
              << " & " << Form("%.2f", backRej) << "\\\\ \\hline" << std::endl;
    }

  texFile << tableEnd << docEnd;
}
