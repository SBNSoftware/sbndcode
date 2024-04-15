#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"

void SlicingPerformance(const TString productionVersion)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/thesis/slicingperformance/" + productionVersion;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  std::vector<Plot> plots = { { "slc_comp", "slc_comp * 100.", ";Completeness (%);Slices (Normalised)", 25, 0, 100, kMagenta+2 },
                              { "slc_pur", "slc_pur * 100.", ";Purity (%);Slices (Normalised)", 25, 0, 100, kMagenta+2 }
  };

  const Cut cut = { "neutrinoSlicesMatchedToNeutrinos", "slc_true_event_type_incl<5 && !slc_is_clear_cosmic" };

  for(auto plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      MakePlot(canvas, rockboxEvents, plot, cut, 1., true);

      canvas->SaveAs(saveDir + "/" + plot.name + ".png");
      canvas->SaveAs(saveDir + "/" + plot.name + ".pdf");
      canvas->SaveAs(saveDir + "/" + plot.name + ".C");

      delete canvas;
    }
}
