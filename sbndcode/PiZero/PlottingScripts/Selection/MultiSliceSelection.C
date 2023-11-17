#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"

const std::vector<Cut> categories = {
  { "Signal", "msc_signal && msc_matching_flash_pe", "Signal (x100)", kMagenta+2, 100 },
  { "SignalDiffFlash", "msc_signal && !msc_matching_flash_pe", "Signal (Diff Flash) (x100)", kMagenta-7, 100 },
  { "Background", "!msc_signal && msc_matching_flash_pe", "Background", kRed+1 },
  { "BackgroundDiffFlash", "!msc_signal && !msc_matching_flash_pe", "Background (Diff Flash)", kRed-7 },
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
  { "good_opt0s_cut", "msc_0_good_opt0 && msc_1_good_opt0", "Good OpT0s" },
  { "opt0_frac_0_cut", "msc_0_opt0_frac<0", "OpT0Frac0 \\textless 0" },
  { "opt0_frac_1_cut", "msc_1_opt0_frac<0", "OpT0Frac1 \\textless 0" },
  { "sum_opt0_frac_cut", "!msc_matching_flash_pe || (msc_sum_opt0_frac<0.2 && msc_sum_opt0_frac>-0.6)", "-0.6 \\textless sumOpT0Frac \\textless 0.2" },
  //  { "sep_cut", "msc_sep<200", "Separation \\textless 200cm" },
  { "dca_cut", "msc_dca<80", "DCA \\textless 80cm" },
  { "crumbs_cut", "msc_0_crumbs_score>-0.2 || msc_1_crumbs_score>-0.2", "Higher CRUMBS Score \\textgreater -0.2"},
  { "razzled_muons_cut", "msc_0_n_razzled_muons==0 && msc_1_n_razzled_muons==0", "No Razzled Muons" },
  { "razzled_photons_cut", "(msc_0_n_razzled_photons + msc_1_n_razzled_photons) > 1", "Has Two Razzled Photons" },
  { "razzled_photons_exact_cut", "(msc_0_n_razzled_photons + msc_1_n_razzled_photons) == 2", "Has Exactly Two Razzled Photons" },
};

std::vector<Plot> plots = {
  { "msc_0_opt0_frac", "msc_0_opt0_frac", ";OpT0Frac0;Candidates",
    100, -2, 8 },
  { "msc_0_opt0_frac_close", "msc_0_opt0_frac", ";OpT0Frac0;Candidates",
    48, -1.2, 1 },
  { "msc_1_opt0_frac", "msc_1_opt0_frac", ";OpT0Frac1;Candidates",
    100, -2, 8 },
  { "msc_1_opt0_frac_close", "msc_1_opt0_frac", ";OpT0Frac1;Candidates",
    48, -1.2, 1 },
  { "msc_sum_opt0_frac", "msc_sum_opt0_frac", ";SumOpT0Frac;Candidates",
    100, -2, 8 },
  { "msc_sum_opt0_frac_close", "msc_sum_opt0_frac", ";SumOpT0Frac;Candidates",
    48, -1.2, 1 },
  { "msc_sep", "msc_sep", ";Separation (cm);Candidates",
    50, 0, 600 },
  { "msc_dca", "msc_dca", ";DCA (cm);Candidates",
    50, 0, 600 },
  { "msc_matching_flash", "msc_matching_flash_pe", ";Same Flash?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, { "No", "Yes" } },
  { "msc_max_crumbs_score", "max(msc_0_crumbs_score, msc_1_crumbs_score)", ";Higher CRUMBS Score;Candidates",
    28, -0.8, 0.6 },
  { "msc_0_crumbs_score", "msc_0_crumbs_score", ";CRUMBS Score 0;Candidates",
    28, -0.8, 0.6 },
  { "msc_1_crumbs_score", "msc_1_crumbs_score", ";CRUMBS Score 1;Candidates",
    28, -0.8, 0.6 },
  { "msc_same_tpc", "msc_same_tpc", ";Same TPC?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, { "No", "Yes" } },
  { "msc_good_opt0s", "msc_0_good_opt0 && msc_1_good_opt0", ";Good OpT0 Matches?;Candidates",
    2, -0.5, 1.5, kBlack, false, "", true, { "No", "Yes" } },
  { "msc_0_n_razzled_muons", "msc_0_n_razzled_muons", ";N Razzled Muons 0;Candidate",
    5, -0.5, 4.5 },
  { "msc_1_n_razzled_muons", "msc_1_n_razzled_muons", ";N Razzled Muons 1;Candidate",
    5, -0.5, 4.5 },
  { "msc_0_n_razzled_photons", "msc_0_n_razzled_photons", ";N Razzled Photons 0;Candidate",
    5, -0.5, 4.5 },
  { "msc_1_n_razzled_photons", "msc_1_n_razzled_photons", ";N Razzled Photons 1;Candidate",
    5, -0.5, 4.5 },
  { "msc_n_razzled_photons", "msc_0_n_razzled_photons + msc_1_n_razzled_photons", ";#Sigma N Razzled Photons;Candidate",
    5, -0.5, 4.5 },
};

template<class T>
void ProduceCutTable(const TString &saveDir, std::vector<Sample<T>> &samples);

void MultiSliceSelection()
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv4/split_slice_selection/selection";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv4/NCPiZeroAv4_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
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
      totalSignalCandidates += sample.tree->Draw("", "msc_signal");
      totalBackCandidates   += sample.tree->Draw("", "!msc_signal");
    }

  texFile << docStart;

  texFile << '\n'
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

      const double eff         = sigCandidates * 100. / totalSignalCandidates;
      const double pur         = sigCandidates * 100./ (sigCandidates + backCandidates);
      const double backRej     = 100. - 100. * backCandidates / totalBackCandidates;

      texFile << cut.printed_name << " & " << Form("%.2f", eff) << " & " << Form("%.2f", pur)
              << " & " << Form("%.2f", (eff * pur) / 100.)
              << " & " << Form("%.2f", backRej)
              << "\\\\ \\hline" << std::endl;
    }

  texFile << tableEnd << docEnd;
}
