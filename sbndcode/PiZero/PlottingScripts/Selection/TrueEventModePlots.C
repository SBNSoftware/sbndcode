#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"
#include "Categories.h"
#include "Plots.h"

const double goalPOT = 10e20;

void TrueEventModePlots(const TString productionVersion, const TString saveDirExt, const std::vector<Cut> &categories,
                        std::vector<Plot> &plots, const Cut &cut);

double GetPOT(TChain *subruns);

void RunMulti()
{
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero", event_modes, true_observables, {"ncpizero", true_ncpizero_cut});
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero_0p0pi", event_modes, true_observables, {"ncpizero_0p0pi", true_ncpizero_0p0pi_cut});
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero_1p0pi", event_modes, true_observables, {"ncpizero_1p0pi", true_ncpizero_1p0pi_cut});
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero_Np0pi", event_modes, true_observables, {"ncpizero_Np0pi", true_ncpizero_Np0pi_cut});
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero_0pXpi", event_modes, true_observables, {"ncpizero_0pXpi", true_ncpizero_0pXpi_cut});
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero_1pXpi", event_modes, true_observables, {"ncpizero_1pXpi", true_ncpizero_1pXpi_cut});
  TrueEventModePlots("NCPiZeroAv4_1", "ncpizero_NpXpi", event_modes, true_observables, {"ncpizero_NpXpi", true_ncpizero_NpXpi_cut});
}

void TrueEventModePlots(const TString productionVersion, const TString saveDirExt, const std::vector<Cut> &categories,
                        std::vector<Plot> &plots, const Cut &cut)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/true_event_modes/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *rockboxsubruns = new TChain("ncpizeroana/subruns");
  rockboxsubruns->Add(rockboxFile);

  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const double rockboxPOT = GetPOT(rockboxsubruns);
  const double rockboxScaling = goalPOT / rockboxPOT;

  std::vector<Sample<TChain>> samples = { { "rockbox", rockboxEvents, rockboxScaling },
  };

  for(auto plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      plot.axes_labels += potString;

      MakeStackedPlot(canvas, samples, plot, cut, categories, {.25, .8, .8, .87}, 5);

      canvas->SaveAs(saveDir + "/" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");

      delete canvas;
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
