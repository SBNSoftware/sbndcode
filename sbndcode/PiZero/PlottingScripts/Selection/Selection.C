#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"
#include "Cuts.h"
#include "Categories.h"
#include "Plots.h"

const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

void Selection(const TString productionVersion, const TString saveDirExt = "tmp", std::vector<Cut> &cuts = razzled_cuts,
               const std::vector<Cut> &categories = selection_categories, std::vector<Plot> &plots = selection_plots);

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples, std::vector<Cut> &cuts,
                     const std::vector<Cut> &categories);

void RunMultiSelection()
{
  Selection("NCPiZeroAv3", "razzled_muons_razzle_photons_cuts", razzled_muons_razzle_photons_cuts);
  Selection("NCPiZeroAv3", "razzled_cuts", razzled_cuts);
}

void Selection(const TString productionVersion, const TString saveDirExt, std::vector<Cut> &cuts, const std::vector<Cut> &categories,
               std::vector<Plot> &plots)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/selection/" + saveDirExt;
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

  ProduceCutTable(saveDir, samples, cuts, categories);

  TCut currentCut = "";

  for(auto cut : cuts)
    {
      currentCut += cut.cut;
      cut.cut = currentCut;

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
                     const std::vector<Cut> &categories)
{
  ofstream texFile;
  texFile.open(saveDir + "/cut_table.tex");

  double totalSignal = 0, totalSignalSlices = 0, totalBackSlices = 0, totalNuFVBackSlices = 0;

  for(auto const& sample : samples)
    {
      totalSignal         += sample.tree->Draw("", "nu_signal");
      totalSignalSlices   += sample.tree->Draw("", "slc_true_signal && slc_comp>.5");
      totalBackSlices     += sample.tree->Draw("", "!slc_true_signal");
      totalNuFVBackSlices += sample.tree->Draw("", "!slc_true_signal && slc_true_event_type!=4 && slc_true_event_type!=5 && slc_true_event_type!=6");
    }

  texFile << docStart;

  texFile << '\n'
          << "Total Signal: " << totalSignal << "\\\\ \n"
          << "Total Signal Slices: " << totalSignalSlices << "\\\\ \n"
          << "Total Background Slices: " << totalBackSlices << "\\\\ \n"
          << "Total Nu FV Background Slices: " << totalNuFVBackSlices << "\\\\ \n" << std::endl;

  texFile << tableStart 
          << "\\hline\n"
          << "Cut Name & $\\epsilon$ (\\%) & $\\rho$ (\\%) & $\\epsilon\\rho$ & Selection $\\epsilon$ (\\%) & BR (\\%) & FV $\\nu$ BR (\\%)\\\\ \\hline" 
          << std::endl;

  TCut currentCut = "";

  for(unsigned i = 0; i < cuts.size(); ++i)
    {
      Cut cut = cuts[i];
      currentCut += cut.cut;
      cut.cut = currentCut;

      double sigSlices = 0., nuFVBackSlices = 0., otherBackSlices = 0.;
      for(unsigned j = 0; j < categories.size(); ++j)
        {
          for(auto const& sample : samples)
            {
              if(j == 0)
                sigSlices += sample.tree->Draw("", cut.cut + categories[j].cut);
              else if(j==4 || j == 5 || j == 6)
                otherBackSlices += sample.tree->Draw("", cut.cut + categories[j].cut);
              else
                nuFVBackSlices += sample.tree->Draw("", cut.cut + categories[j].cut);
            }
        }

      const double eff         = sigSlices * 100. / totalSignal;
      const double selEff      = sigSlices * 100. / totalSignalSlices;
      const double pur         = sigSlices * 100./ (sigSlices + nuFVBackSlices + otherBackSlices);
      const double backRej     = 100. - 100. * ((nuFVBackSlices + otherBackSlices) / totalBackSlices);
      const double nuFVBackRej = 100. - 100. * (nuFVBackSlices / totalNuFVBackSlices);

      texFile << cut.printed_name << " & " << Form("%.2f", eff) << " & " << Form("%.2f", pur)
              << " & " << Form("%.2f", (eff * pur) / 100.)
              << " & " << Form("%.2f", selEff) << " & " << Form("%.2f", backRej)
              << " & " << Form("%.2f", nuFVBackRej) << "\\\\ \\hline" << std::endl;
    }

  texFile << tableEnd << docEnd;
}
