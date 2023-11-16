#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Cuts.h"
#include "Categories.h"
#include "Plots.h"

const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

void Observables(const TString productionVersion, const TString saveDirExt = "tmp", const std::vector<Cut> cuts = ncpizero_cuts,
                 const std::vector<Cut> &categories = selection_categories, std::vector<Plot> &plots = selection_plots);

void Observables(const TString productionVersion, const TString saveDirExt = "tmp", const Cut &selection = ncpizero_cuts[0],
                 const std::vector<Cut> &categories = selection_categories, std::vector<Plot> &plots = selection_plots);

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void RunMultiObservables()
{
  Observables("NCPiZeroAv6", "ncpizero", ncpizero_cuts, ncpizero_categories, observables);
}

void Observables(const TString productionVersion, const TString saveDirExt = "tmp", const std::vector<Cut> cuts = ncpizero_cuts,
                 const std::vector<Cut> &categories = selection_categories, std::vector<Plot> &plots = selection_plots)
{
  TCut totalCut = "";

  for(auto const& cut : cuts)
    totalCut += cut.cut;

  Cut cut = { "full_selection", totalCut, "Full Selection" };

  Observables(productionVersion, saveDirExt, cut, categories, plots);
}
void Observables(const TString productionVersion, const TString saveDirExt, const Cut &selection,
                 const std::vector<Cut> &categories, std::vector<Plot> &plots)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/observables/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxsubruns = new TChain("ncpizeroana/subruns");
  rockboxsubruns->Add(rockboxFile);
  TChain *intimesubruns = new TChain("ncpizeroana/subruns");
  intimesubruns->Add(intimeFile);

  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const double rockboxPOT = GetPOT(rockboxsubruns);
  const int rockboxSpills = GetGenEvents(rockboxsubruns);
  const int intimeSpills  = GetGenEvents(intimesubruns);

  const double rockboxScaling      = goalPOT / rockboxPOT;
  const double scaledRockboxSpills = rockboxScaling * rockboxSpills;
  const double intimeScaling       = (goalSpills - scaledRockboxSpills) / intimeSpills;

  std::vector<Sample<TChain>> samples = { { "rockbox", rockboxEvents, rockboxScaling },
                                          { "intime", intimeEvents, intimeScaling }
  };

  for(auto plot : plots)
    {
      TCanvas *canvas = new TCanvas("c_" + plot.name, "c_" + plot.name);
      canvas->cd();

      plot.axes_labels += potString;

      MakeStackedPlot(canvas, samples, plot, selection, categories, {.25, .8, .8, .87}, 4);

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

int GetGenEvents(TChain *subruns)
{
  int sum = 0., ngenevts = 0;

  subruns->SetBranchAddress("ngenevts", &ngenevts);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += ngenevts;
    }

  return sum;
}
