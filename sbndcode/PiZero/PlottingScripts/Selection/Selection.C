#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "LatexHeaders.h"

const std::vector<Cut> categories = {
  { "Signal", "slc_true_event_type==0 && slc_comp>.5", "Signal (NC #pi^{0})", kMagenta+2 },
  { "NC", "slc_true_event_type==1", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type==2", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type==3", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type==4", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type==5", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type==6", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

std::vector<Cut> cuts = {
  { "no_cut", "", "No Cut" },
  { "not_clear_cosmic", "!slc_is_clear_cosmic", "Not Clear Cosmic" },
  { "fv", "slc_is_fv", "FV" },
  { "crumbs", "slc_crumbs_score>-0.025", "CRUMBS Cut" },
  { "no_dazzle_muons", "slc_n_dazzle_muons_cut_based==0", "No Dazzle Muons" },
  { "has_two_shws", "slc_n_shws>1", "Has Two Showers" },
  { "has_two_razzle_photons", "slc_n_razzle_photons>1", "Has Two Razzle Photons" },
};

std::vector<Plot> plots = {
  { "slc_is_clear_cosmic", "slc_is_clear_cosmic", ";Is Clear Cosmic?;Slices",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_is_fv", "slc_is_fv", ";IsFV?;Slice",
    2, -0.5, 1.5, kBlack, false, "", true, {"No", "Yes"} },
  { "slc_crumbs_score", "slc_crumbs_score", ";CRUMBS Score;Slices",
    42, -1.5, .6 },
  { "slc_crumbs_nc_score", "slc_crumbs_nc_score", ";CRUMBS NC Score;Slices",
    42, -1.5, .6 },
  { "slc_crumbs_ccnue_score", "slc_crumbs_ccnue_score", ";CRUMBS CC#nu_{e} Score;Slices",
    42, -1.5, .6 },
  { "slc_n_dazzle_muons_cut_based", "slc_n_dazzle_muons_cut_based", ";N Dazzle Muons;Slices",
    5, -0.5, 4.5 },
  { "slc_n_shws", "slc_n_shws", ";N Showers;Slices",
    5, -0.5, 4.5 },
  { "slc_n_razzle_photons", "slc_n_razzle_photons", ";N Razzle Photons;Slices",
    5, -0.5, 4.5 },
};

double GetPOT(TChain *subruns);
void ProduceCutTable(const TString &saveDir, TChain *events);

void Selection(const TString saveDirExt = "tmp")
{
  const TString saveDir = "/sbnd/data/users/hlay/pizero/plots/selection/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);
  const TString file    = "/sbnd/data/users/hlay/pizero/ncpizeroana_sbnd.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add("/sbnd/data/users/hlay/pizero/ncpizeroana_sbnd.root");

  TChain *subruns = new TChain("ncpizeroana/subruns");
  subruns->Add("/sbnd/data/users/hlay/pizero/ncpizeroana_sbnd.root");

  const double goalPOT = 10e20;
  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const double totalPOT = GetPOT(subruns);
  const double scaling  = goalPOT / totalPOT;

  ProduceCutTable(saveDir, events);

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

          MakeStackedPlot(canvas, events, plot, cut, categories, scaling, {.25, .8, .8, .87}, 4);

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

void ProduceCutTable(const TString &saveDir, TChain *events)
{
  ofstream texFile;
  texFile.open(saveDir + "/cut_table.tex");

  const double totalSignal         = events->Draw("", "nu_signal");
  const double totalSignalSlices   = events->Draw("", "slc_true_signal && slc_comp>.5");
  const double totalBackSlices     = events->Draw("", "!slc_true_signal");
  const double totalNuFVBackSlices = events->Draw("", "!slc_true_signal && slc_true_event_type!=5 && slc_true_event_type!=6");

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
          if(j == 0)
            sigSlices = events->Draw("", cut.cut + categories[j].cut);
          else if(j == 5 || j == 6)
            otherBackSlices += events->Draw("", cut.cut + categories[j].cut);
          else
            nuFVBackSlices += events->Draw("", cut.cut + categories[j].cut);
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
