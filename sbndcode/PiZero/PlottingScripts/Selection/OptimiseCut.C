#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Cuts.h"
#include "Categories.h"
#include "Plots.h"

const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

void OptimiseCut(const TString productionVersion, const TString saveDirExt, Plot plot, const Cut signal_def, const Cut base_cut, const bool rej_cut);

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void RunMultiOptimiseCut()
{
  OptimiseCut("NCPiZeroAv4_1", "ncpizero", selection_plots[2], ncpizero_categories[0], ncpizero_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ncpizero_0p0pi", selection_plots[2], ncpizero_0p0pi_categories[0], ncpizero_0p0pi_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ncpizero_1p0pi", selection_plots[2], ncpizero_1p0pi_categories[0], ncpizero_1p0pi_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ncpizero_Np0pi", selection_plots[2], ncpizero_Np0pi_categories[0], ncpizero_Np0pi_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ncpizero_0pXpi", selection_plots[2], ncpizero_0pXpi_categories[0], ncpizero_0pXpi_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ncpizero_1pXpi", selection_plots[2], ncpizero_1pXpi_categories[0], ncpizero_1pXpi_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ncpizero_NpXpi", selection_plots[2], ncpizero_NpXpi_categories[0], ncpizero_NpXpi_cuts[2], false);
  OptimiseCut("NCPiZeroAv4_1", "ccpizero", selection_plots[2], ccpizero_categories[0], ccpizero_cuts[2], false);
}

void OptimiseCut(const TString productionVersion, const TString saveDirExt, Plot plot, const Cut signal_def, const Cut base_cut, const bool rej_cut)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/cut_optimisation/" + saveDirExt;
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

  const double rockboxPOT = GetPOT(rockboxsubruns);
  const int rockboxSpills = GetGenEvents(rockboxsubruns);
  const int intimeSpills  = GetGenEvents(intimesubruns);

  const double rockboxScaling      = goalPOT / rockboxPOT;
  const double scaledRockboxSpills = rockboxScaling * rockboxSpills;
  const double intimeScaling       = (goalSpills - scaledRockboxSpills) / intimeSpills;

  TH1F* hSig     = new TH1F("Sig", plot.axes_labels, plot.nbins * 2, plot.xlow, plot.xhigh);
  TH1F* hBackNu  = new TH1F("BackNu", plot.axes_labels, plot.nbins * 2, plot.xlow, plot.xhigh);
  TH1F* hBackCos = new TH1F("BackCos", plot.axes_labels, plot.nbins * 2, plot.xlow, plot.xhigh);

  rockboxEvents->Draw(plot.var + ">>Sig", signal_def.cut + base_cut.cut);
  rockboxEvents->Draw(plot.var + ">>BackNu", !signal_def.cut + base_cut.cut);
  intimeEvents->Draw(plot.var + ">>BackCos", !signal_def.cut + base_cut.cut);

  hSig->Scale(rockboxScaling);
  hBackNu->Scale(rockboxScaling);
  hBackCos->Scale(intimeScaling);

  TH1F* hBack = (TH1F*) hBackNu->Clone();
  hBack->Add(hBackCos);
  std::cout << signal_def.name << std::endl;
  std::cout << base_cut.name << std::endl;

  const double totalSignal = rockboxScaling * rockboxEvents->Draw("", signal_def.cut);
  const double initSignal  = rockboxScaling * rockboxEvents->Draw("", signal_def.cut + base_cut.cut);
  const double initBack    = rockboxScaling * rockboxEvents->Draw("", !signal_def.cut + base_cut.cut) + intimeScaling * intimeEvents->Draw("", !signal_def.cut + base_cut.cut);

  TGraph *gSelEff = new TGraph();
  gSelEff->SetMarkerColor(kRed+2);
  TGraph *gPur = new TGraph();
  gPur->SetMarkerColor(kBlue+2);
  TGraph *gEP  = new TGraph();
  gEP->SetMarkerColor(kMagenta+2);

  gSelEff->SetTitle(TString(hSig->GetXaxis()->GetTitle()) + ";Cut Value;Fraction");

  float accumSignal = initSignal, accumBack = initBack, maxEP = -std::numeric_limits<float>::max(), optimal_cut = -std::numeric_limits<float>::max(),
    max = -std::numeric_limits<float>::max();

  int i = rej_cut ? hSig->GetNbinsX()+1 : 0;
  int N = rej_cut ? 0 : hSig->GetNbinsX()+1;

  while((rej_cut && i > N) || (!rej_cut && i < N))
    {
      accumSignal -= hSig->GetBinContent(i);
      accumBack   -= hBack->GetBinContent(i);

      if(accumSignal < 0)
        break;

      const float eff = (accumSignal / totalSignal) * 100;
      const float pur = (accumSignal / (accumSignal + accumBack)) * 100;
      const float ep  = eff * pur / 100;
      const float cut = rej_cut ? hSig->GetBinLowEdge(i) : hSig->GetBinLowEdge(i) + hSig->GetBinWidth(i);

      std::cout << "Cut: " << cut << " Signal: " << accumSignal << " Back: " << accumBack << std::endl;
      std::cout << "\tEff: " << eff << " Pur: " << pur << " Eff*Pur: " << ep << std::endl;

      if(ep > maxEP)
        {
          maxEP = ep;
          optimal_cut = cut;
        }

      if(eff > max)
        max = eff;
      if(pur > max)
        max = pur;
      if(ep > max)
        max = ep;

      gSelEff->SetPoint(i, cut, eff);
      gPur->SetPoint(i, cut, pur);
      gEP->SetPoint(i, cut, eff * pur / 100);

      if(rej_cut)
        --i;
      else
        ++i;
    }

  TLine *line = new TLine(optimal_cut, 0, optimal_cut, 110);
  line->SetLineStyle(9);
  line->SetLineWidth(3);

  TPaveText *pt = new TPaveText(.26, .38, .5, .43, "NDC");
  pt->AddText(Form("Optimal Cut: %.4f", optimal_cut));
  pt->AddText(Form("Efficiency x Purity: %.2f %%", maxEP));
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.03);
  pt->SetFillColor(kWhite);
  pt->SetBorderSize(0);
  pt->SetLineWidth(0);
  pt->SetTextAlign(11);

  TLegend *leg = new TLegend(.23, .45, .5, .55);
  leg->AddEntry(gSelEff, "Efficiency", "p");
  leg->AddEntry(gPur, "Purity", "p");
  leg->AddEntry(gEP, "Efficiency x Purity", "p");

  TCanvas *c = new TCanvas("c", "c");
  c->SetTopMargin(.15);
  gSelEff->Draw("AP");
  gSelEff->SetMinimum(0);
  gSelEff->SetMaximum(110);
  gPur->Draw("Psame");
  gEP->Draw("Psame");
  line->Draw();
  pt->Draw();
  leg->Draw();

  c->SaveAs(saveDir + "/" + plot.name + ".png");
  c->SaveAs(saveDir + "/" + plot.name + ".pdf");
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
