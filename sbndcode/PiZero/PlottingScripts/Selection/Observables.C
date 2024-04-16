#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void Observables(const TString productionVersion, const SelectionParams &selectionParams, std::vector<VarBinPlot> &plots)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables/" + selectionParams.name;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile  = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";
  const TString intimeFile   = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);
  TChain *ncpizeroSubruns = new TChain("ncpizeroana/subruns");
  ncpizeroSubruns->Add(ncpizeroFile);
  TChain *intimeSubruns = new TChain("ncpizeroana/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, ncpizeroScaling, intimeScaling;
  GetScaling(rockboxSubruns, ncpizeroSubruns, intimeSubruns, rockboxScaling, ncpizeroScaling, intimeScaling);

  std::vector<Sample<TChain>> samples = { { "rockbox", rockboxEvents, rockboxScaling, selectionParams.rockbox_mask },
                                          { "ncpizero", ncpizeroEvents, ncpizeroScaling, selectionParams.ncpizero_mask },
                                          { "intime", intimeEvents, intimeScaling }
  };

  Cut cut = TotalCut(selectionParams.cuts);

  for(auto plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      const int ncolumns = 3;
      const float xlow   = .24;
      const float xhigh  = .83;
      const float ylow   = .78;
      const float yhigh  = .87;

      MakeStackedPlot(canvas, samples, plot, cut, selectionParams.categories, {xlow, ylow, xhigh, yhigh}, ncolumns);

      const TString wip = "SBND Work-in-progress";
      AddText(canvas, wip, kGray+2, {.2, .92, .25, .93}, 0.035, 12);
      const TString sim = "SBND Simulation";
      AddText(canvas, sim, kGray+2, {.2, .96, .25, .97}, 0.035, 12);
      const TString potString = POTString(false);
      AddText(canvas, potString, kGray+2, {.8, .92, .9, .93}, 0.035, 32);

      canvas->SaveAs(saveDir + "/" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");

      delete canvas;
    }
}
