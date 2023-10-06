#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"

const std::vector<Cut> categories = {
  { "Signal", "signal_cand && matching_flash_pe", "Signal (x100)", kMagenta+2, 100 },
  { "SignalDiffFlash", "signal_cand && !matching_flash_pe", "Signal (Diff Flash) (x100)", kMagenta-7, 100 },
  { "Background", "!signal_cand && matching_flash_pe", "Background", kRed+1 },
  { "BackgroundDiffFlash", "!signal_cand && !matching_flash_pe", "Background (Diff Flash)", kRed-7 },
};

/*
  const std::vector<Cut> categories = {
  { "Signal", "good_cand && matching_flash_pe", "Signal", kMagenta+2 },
  { "SignalDiffFlash", "good_cand && !matching_flash_pe", "Signal (Diff Flash)", kMagenta-7 },
  { "Background", "!good_cand && matching_flash_pe", "Background", kRed+1 },
  { "BackgroundDiffFlash", "!good_cand && !matching_flash_pe", "Background (Diff Flash)", kRed-7 },
  };
*/

std::vector<Cut> cuts = {
  { "no_cut", "", "No Cut" },
  { "good_opT0s_cut", "goodOpT00 && goodOpT01", "Good OpT0s" },
  { "opT0Frac0_cut", "opT0Frac0<0", "OpT0Frac0 \\textless 0" },
  { "opT0Frac1_cut", "opT0Frac1<0", "OpT0Frac1 \\textless 0" },
  { "sumOpT0Frac_cut", "!matching_flash_pe || (sumOpT0Frac<0.2 && sumOpT0Frac>-0.6)", "-0.6 \\textless sumOpT0Frac \\textless 0.2" },
  //  { "either_opT0Frac_cut", "matching_flash_pe || (opT0Frac0>-0.6 && opT0Frac1>-0.6)", "-0.6 \\textless OpT0Frac0,1 \\textless 0" },
  { "sep_cut", "sep<200", "Separation \\textless 200cm" },
  { "crumbs_cut", "crumbs0>-0.2 || crumbs1>-0.2", "Higher CRUMBS Score \\textgreater -0.2"},
};

std::vector<Plot> plots = {
  { "opT0Frac0", "opT0Frac0", ";OpT0Frac0;Candidates",
    100, -2, 8 },
  { "opT0Frac0_close", "opT0Frac0", ";OpT0Frac0;Candidates",
    48, -1.2, 1 },
  { "opT0Frac1", "opT0Frac1", ";OpT0Frac1;Candidates",
    100, -2, 8 },
  { "opT0Frac1_close", "opT0Frac1", ";OpT0Frac1;Candidates",
    48, -1.2, 1 },
  { "sumOpT0Frac", "sumOpT0Frac", ";SumOpT0Frac;Candidates",
    100, -2, 8 },
  { "sumOpT0Frac_close", "sumOpT0Frac", ";SumOpT0Frac;Candidates",
    48, -1.2, 1 },
  { "sep", "sep", ";Separation (cm);Candidates",
    50, 0, 600 },
  { "matching_flash", "matching_flash_pe", ";Same Flash?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, { "No", "Yes" } },
  { "max_crumbs", "max(crumbs0, crumbs1)", ";Higher CRUMBS Score;Candidates",
    28, -0.8, 0.6 },
  { "crumbs0", "crumbs0", ";CRUMBS Score 0;Candidates",
    28, -0.8, 0.6 },
  { "crumbs1", "crumbs1", ";CRUMBS Score 1;Candidates",
    28, -0.8, 0.6 },
  { "same_tpcs", "(vtx_x0 > 0 && vtx_x1 > 0) + (vtx_x0 < 0 && vtx_x1 < 0)", ";Same TPC?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, { "No", "Yes" } },
  { "good_opT0s", "goodOpT00 && goodOpT01", ";Good OpT0 Matches?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, { "No", "Yes" } },
};

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples);

void MultiSliceSelection()
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/split_slice_selection/selection";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = "/sbnd/data/users/hlay/ncpizero/split_slice_selection/NCPiZeroAv2_rockbox_multislice_tree.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("events");
  events->Add(file);

  std::vector<Sample<TChain>> samples = { { "all", events, 1. }
  };

  ProduceCutTable(saveDir, samples);

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

          MakeStackedPlot(canvas, samples, plot, cut, categories, {.25, .8, .8, .87}, 2);

          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
          canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");

          delete canvas;
        }
    }

  gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/cut_table.tex");
}

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples)
{
  ofstream texFile;
  texFile.open(saveDir + "/cut_table.tex");

  double totalSignal = 0, totalSignalCandidates = 0, totalBackCandidates = 0;

  for(auto const& sample : samples)
    {
      totalSignal           += sample.tree->Draw("", "signal");
      totalSignalCandidates += sample.tree->Draw("", "signal_cand");
      totalBackCandidates   += sample.tree->Draw("", "!signal_cand");
    }

  texFile << docStart;

  texFile << '\n'
          << "Total Signal: " << totalSignal << "\\\\ \n"
          << "Total Signal Candidates: " << totalSignalCandidates << "\\\\ \n"
          << "Total Background Candidates: " << totalBackCandidates << "\\\\ \n"
          << std::endl;

  texFile << tableStart 
          << "\\hline\n"
          << "Cut Name & $\\epsilon$ (\\%) & $\\rho$ (\\%) & $\\epsilon\\rho$ & BR (\\%)\\\\ \\hline"
          << std::endl;

  TCut currentCut = "";

  for(unsigned i = 0; i < cuts.size(); ++i)
    {
      Cut cut = cuts[i];
      currentCut += cut.cut;
      cut.cut = currentCut;

      double sigCandidates = 0., backCandidates = 0;

      for(unsigned j = 0; j < categories.size(); ++j)
        {
          for(auto const& sample : samples)
            {
              if(j == 0 || j == 1)
                sigCandidates += sample.tree->Draw("", cut.cut + categories[j].cut);
              else
                backCandidates += sample.tree->Draw("", cut.cut + categories[j].cut);
            }
        }

      const double eff         = sigCandidates * 100. / totalSignal;
      const double pur         = sigCandidates * 100./ (sigCandidates + backCandidates);
      const double backRej     = 100. - 100. * backCandidates / totalBackCandidates;

      texFile << cut.printed_name << " & " << Form("%.2f", eff) << " & " << Form("%.2f", pur)
              << " & " << Form("%.2f", (eff * pur) / 100.)
              << " & " << Form("%.2f", backRej)
              << "\\\\ \\hline" << std::endl;
    }

  texFile << tableEnd << docEnd;
}
