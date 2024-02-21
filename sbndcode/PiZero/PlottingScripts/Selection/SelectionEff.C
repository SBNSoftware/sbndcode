#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Particles.h"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void SelectionEff(const TString productionVersion, const SelectionParams &selectionParams, const std::vector<VarBinPlot> &plots)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/selection_eff/" + selectionParams.name;
  gSystem->Exec("mkdir -p " + saveDir);

  const std::array<float, 4> legend_position = { .23, .82, .87, .91 };
  const int ncolumns                         = 3;

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  for(auto const &plot : plots)
    {
      TCanvas *canvas = new TCanvas("canvas" + plot.name, "canvas" + plot.name);
      canvas->cd();

      TString trueVarName = plot.var;
      trueVarName.ReplaceAll("slc_true", "nu");

      float bins[plot.bins.size()];
      for(int i = 0; i < plot.bins.size(); ++i)
        bins[i] = plot.bins[i];

      TH1F *trueHist = new TH1F("trueHist" + plot.name, plot.axes_labels, plot.nbins, bins);
      rockboxEvents->Draw(trueVarName + ">>trueHist" + plot.name, selectionParams.true_category);
      NormaliseEntriesByBinWidth(trueHist);

      TCut currentCut = "";

      std::vector<TH1*> recoHists;
      std::vector<int> colours;
      std::vector<TString> names;

      for(auto const& cut : selectionParams.broad_cuts)
        {
          currentCut += cut.cut;

          recoHists.emplace_back(new TH1F("recoHist" + plot.name + cut.name, plot.axes_labels, plot.nbins, bins));
          rockboxEvents->Draw(plot.var + ">>recoHist" + plot.name + cut.name, selectionParams.categories[0].cut + currentCut);
          NormaliseEntriesByBinWidth(recoHists.back());

          colours.push_back(cut.colour);
          names.push_back(cut.printed_name);
        }

      MakePlotMultiEff(canvas, trueHist, recoHists, plot.axes_labels, colours, names, legend_position, ncolumns);

      canvas->SaveAs(saveDir + "/selection_efficiency_" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/selection_efficiency_" + plot.name + ".pdf");
    }
}
