#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"
#include "Cuts.h"
#include "Categories.h"
#include "Plots.h"

const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

void Selection(const TString productionVersion, const TString saveDirExt = "tmp", std::vector<Cut> &cuts = ncpizero_incl_cuts,
               const std::vector<Cut> &categories = ncpizero_incl_categories, std::vector<Plot> &plots = selection_plots,
               const TCut &trueCategory = "");

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples, std::vector<Cut> &cuts,
                     const std::vector<Cut> &categories, const TCut &trueCategory);

void RunMultiSelection()
{
  Selection("NCPiZeroAv10", "ncpizero_incl", ncpizero_incl_cuts, ncpizero_incl_categories, no_plots, true_ncpizero_incl_cut);
  Selection("NCPiZeroAv10", "ncpizero_0p0pi", ncpizero_0p0pi_cuts, ncpizero_0p0pi_categories, no_plots, true_ncpizero_0p0pi_cut);
  Selection("NCPiZeroAv10", "ncpizero_1p0pi", ncpizero_1p0pi_cuts, ncpizero_1p0pi_categories, no_plots, true_ncpizero_1p0pi_cut);
  Selection("NCPiZeroAv10", "ncpizero_Np0pi", ncpizero_Np0pi_cuts, ncpizero_Np0pi_categories, no_plots, true_ncpizero_Np0pi_cut);
  Selection("NCPiZeroAv10", "ncpizero_Xp0pi", ncpizero_Xp0pi_cuts, ncpizero_Xp0pi_categories, no_plots, true_ncpizero_Xp0pi_cut);

  Selection("NCPiZeroAv10", "ncpizero_incl", ncpizero_incl_cuts, ncpizero_incl_categories, selection_plots, true_ncpizero_incl_cut);
  Selection("NCPiZeroAv10", "ncpizero_0p0pi", ncpizero_0p0pi_cuts, ncpizero_0p0pi_categories, selection_plots, true_ncpizero_0p0pi_cut);
  Selection("NCPiZeroAv10", "ncpizero_1p0pi", ncpizero_1p0pi_cuts, ncpizero_1p0pi_categories, selection_plots, true_ncpizero_1p0pi_cut);
  Selection("NCPiZeroAv10", "ncpizero_Np0pi", ncpizero_Np0pi_cuts, ncpizero_Np0pi_categories, selection_plots, true_ncpizero_Np0pi_cut);
  Selection("NCPiZeroAv10", "ncpizero_Xp0pi", ncpizero_Xp0pi_cuts, ncpizero_Xp0pi_categories, selection_plots, true_ncpizero_Xp0pi_cut);
}

void Selection(const TString productionVersion, const TString saveDirExt, std::vector<Cut> &cuts, const std::vector<Cut> &categories,
               std::vector<Plot> &plots, const TCut &trueCategory)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/selection/" + saveDirExt;
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

  ProduceCutTable(saveDir, samples, cuts, categories, trueCategory);

  TCut currentCut = "";

  for(auto cut : cuts)
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

          plot.axes_labels += potString;

          MakeStackedPlot(canvas, samples, plot, cut, categories, {.25, .8, .8, .87}, 4);

          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");

          delete canvas;
        }
    }

  gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/cut_table.tex");
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

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples, std::vector<Cut> &cuts,
                     const std::vector<Cut> &categories, const TCut &trueCategory)
{
  ofstream texFile;
  texFile.open(saveDir + "/cut_table.tex");

  double totalSignal = 0, totalSignalSlices = 0, totalBackSlices = 0;

  for(auto const& sample : samples)
    {
      totalSignal       += sample.scaling * sample.tree->Draw("", trueCategory);
      totalSignalSlices += sample.scaling * sample.tree->Draw("", categories[0].cut);
      totalBackSlices   += sample.scaling * sample.tree->Draw("", !categories[0].cut);
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

  for(unsigned i = 0; i < cuts.size(); ++i)
    {
      Cut cut = cuts[i];
      currentCut += cut.cut;
      cut.cut = currentCut;

      double sigSlices = 0., backSlices = 0.;
      for(unsigned j = 0; j < categories.size(); ++j)
        {
          for(auto const& sample : samples)
            {
              if(j == 0)
                sigSlices += sample.scaling * sample.tree->Draw("", cut.cut + categories[j].cut);
              else
                backSlices += sample.scaling * sample.tree->Draw("", cut.cut + categories[j].cut);
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
