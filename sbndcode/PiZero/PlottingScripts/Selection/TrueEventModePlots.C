#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void TrueEventModePlots(const TString productionVersion, const std::vector<Plot> &plots, const std::vector<Cut> &signals)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_event_modes";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);

  const double rockboxPOT     = GetPOT(rockboxSubruns);
  const double rockboxScaling = goalPOT / rockboxPOT;

  std::vector<Sample<TChain>> samples = { { "rockbox", rockboxEvents, rockboxScaling } };

  for(auto const& signal : signals)
    {
      gSystem->Exec("mkdir -p " + saveDir + "/" + signal.name);

      for(auto plot : plots)
        {
          TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
          canvas->cd();

          plot.axes_labels += POTString();

          MakeStackedPlot(canvas, samples, plot, signal, event_modes, {.25, .8, .8, .87}, 5);

          canvas->SaveAs(saveDir + "/" + signal.name + "/" + plot.name + ".png");
          canvas->SaveAs(saveDir + "/" + signal.name + "/" + plot.name + ".pdf");

          delete canvas;
        }
    }
}
