#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Particles.h"
#include "Plots.h"
#include "Cuts.h"
#include "Categories.h"

void SelectionEff(const TString productionVersion, const TString saveDirExt, std::vector<Cut> &cuts, const std::vector<Cut> &categories, const TCut &trueCategory)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/selection_eff/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const std::array<float, 4> legend_position = { .23, .82, .87, .91 };
  const int ncolumns                         = 3;

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  for(auto const &plot : true_slc_observables)
    {
      TCanvas *canvas = new TCanvas("canvas" + plot.name, "canvas" + plot.name);
      canvas->cd();

      TString trueVarName = plot.var;
      trueVarName.ReplaceAll("slc_true", "nu");

      float bins[plot.bins.size()];
      for(int i = 0; i < plot.bins.size(); ++i)
        bins[i] = plot.bins[i];

      TH1F *trueHist = new TH1F("trueHist" + plot.name, plot.axes_labels, plot.nbins, bins);
      rockboxEvents->Draw(trueVarName + ">>trueHist" + plot.name, trueCategory);
      NormaliseEntriesByBinWidth(trueHist);

      TCut currentCut = "";

      std::vector<TH1*> recoHists;
      std::vector<int> colours;
      std::vector<TString> names;

      for(auto const& cut : cuts)
        {
          currentCut += cut.cut;

          recoHists.emplace_back(new TH1F("recoHist" + plot.name + cut.name, plot.axes_labels, plot.nbins, bins));
          rockboxEvents->Draw(plot.var + ">>recoHist" + plot.name + cut.name, categories[0].cut + currentCut);
          NormaliseEntriesByBinWidth(recoHists.back());

          colours.push_back(cut.colour);
          names.push_back(cut.printed_name);
        }

      MakePlotMultiEff(canvas, trueHist, recoHists, plot.axes_labels, colours, names, legend_position, ncolumns);

      canvas->SaveAs(saveDir + "/selection_efficiency_" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/selection_efficiency_" + plot.name + ".pdf");
    }
}
