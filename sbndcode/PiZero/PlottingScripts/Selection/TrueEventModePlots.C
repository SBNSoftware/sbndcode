#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void TrueEventModePlots(const TString productionVersion, const std::vector<VarBinPlot> &plots, const std::vector<Cut> &signals)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_event_modes";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(rockboxFile);
  events->Add(ncpizeroFile);
  TChain *subruns = new TChain("ncpizeroana/subruns");
  subruns->Add(rockboxFile);
  subruns->Add(ncpizeroFile);

  const double pot     = GetPOT(subruns);
  const double scaling = goalPOT / pot;

  std::vector<Sample<TChain>> samples = { { "rockbox", events, scaling } };

  for(auto const& signal : signals)
    {
      gSystem->Exec("mkdir -p " + saveDir + "/" + signal.name);

      for(auto plot : plots)
        {
          TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
          canvas->cd();

          gStyle->SetLabelSize(0.06, "y");
          gStyle->SetTitleSize(0.06, "y");
          gStyle->SetLabelSize(0.06, "x");
          gStyle->SetTitleSize(0.06, "x");

          MakeStackedPlot(canvas, samples, plot, signal, event_modes_simple, {.25, .8, .8, .87}, 4);

          const TString wip = "SBND Work-in-progress";
          AddText(canvas, wip, kGray+2, {.2, .92, .25, .93}, 0.03, 12);
          const TString sim = "SBND Simulation";
          AddText(canvas, sim, kGray+2, {.2, .95, .25, .96}, 0.03, 12);
          const TString potString = POTString(false);
          AddText(canvas, potString, kGray+2, {.8, .92, .9, .93}, 0.03, 32);

          canvas->SaveAs(saveDir + "/" + signal.name + "/" + plot.name + ".png");
          canvas->SaveAs(saveDir + "/" + signal.name + "/" + plot.name + ".pdf");

          delete canvas;
        }
    }
}
