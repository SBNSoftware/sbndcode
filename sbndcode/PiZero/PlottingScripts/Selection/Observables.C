#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void Observables(const TString productionVersion, const SelectionParams &selectionParams, std::vector<Plot> &plots)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables/" + selectionParams.name;
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

  Cut cut = TotalCut(selectionParams.cuts);

  for(auto plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      plot.axes_labels += POTString();

      const int ncolumns = selectionParams.name == "ncpizero_incl" ? 4 : 3;
      const float xlow   = selectionParams.name == "ncpizero_incl" ? .25 : .24;
      const float xhigh  = selectionParams.name == "ncpizero_incl" ? .8 : .83;
      const float ylow   = selectionParams.name == "ncpizero_incl" ? .8 : .78;

      MakeStackedPlot(canvas, samples, plot, cut, selectionParams.categories, {xlow, ylow, xhigh, .87}, ncolumns);

      const TString wip = "SBND Work-in-progress";
      AddText(canvas, wip, kGray+2, {.8, .92, .9, .93}, 0.025, 32);

      canvas->SaveAs(saveDir + "/" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");

      delete canvas;
    }
}
