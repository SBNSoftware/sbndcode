#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Plots.h"
#include "Common.C"
#include "Selections.h"

void OptimiseCut(const TString productionVersion, const TString saveDirExt, Plot plot, const Cut signal_def, const Cut base_cut, const bool rej_cut, const bool optimiseEPP);

void OptimiseCut(const TString productionVersion, const SelectionParams &selectionParams, const int excludingCut, Plot plot, const bool rej_cut, const bool optimiseEPP = false)
{
  Cut cut = TotalCutExcluding(selectionParams.cuts, excludingCut);

  OptimiseCut(productionVersion, selectionParams.name, plot, selectionParams.categories[0], cut, rej_cut, optimiseEPP);
}

void OptimiseCut(const TString productionVersion)
{
  OptimiseCut(productionVersion, ncpizero_incl, 3, optimisation_plots[1], false);
  OptimiseCut(productionVersion, ncpizero_incl, 7, optimisation_plots[2], true);
  OptimiseCut(productionVersion, ncpizero_incl, 8, optimisation_plots[2], false);
  OptimiseCut(productionVersion, ncpizero_incl, 9, optimisation_plots[3], false);

  OptimiseCut(productionVersion, ncpizero_0p0pi, 3, optimisation_plots[1], false);
  OptimiseCut(productionVersion, ncpizero_0p0pi, 9, optimisation_plots[2], true);
  OptimiseCut(productionVersion, ncpizero_0p0pi, 10, optimisation_plots[2], false);
  OptimiseCut(productionVersion, ncpizero_0p0pi, 11, optimisation_plots[3], false);

  OptimiseCut(productionVersion, ncpizero_Np0pi, 3, optimisation_plots[1], false);
  OptimiseCut(productionVersion, ncpizero_Np0pi, 9, optimisation_plots[2], true);
  OptimiseCut(productionVersion, ncpizero_Np0pi, 10, optimisation_plots[2], false);
  OptimiseCut(productionVersion, ncpizero_Np0pi, 11, optimisation_plots[3], false);
}

void OptimiseCut(const TString productionVersion, const TString saveDirExt, Plot plot, const Cut signal_def, const Cut base_cut, const bool rej_cut, const bool optimiseEPP)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/cut_optimisation/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile  = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);
  TChain *intimeSubruns = new TChain("ncpizeroana/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  TH1F* hSig     = new TH1F(Form("Sig_%s_%s", plot.name.Data(), saveDirExt.Data()), plot.axes_labels, plot.nbins * 10, plot.xlow, plot.xhigh);
  TH1F* hBackNu  = new TH1F(Form("BackNu_%s_%s", plot.name.Data(), saveDirExt.Data()), plot.axes_labels, plot.nbins * 10, plot.xlow, plot.xhigh);
  TH1F* hBackCos = new TH1F(Form("BackCos_%s_%s", plot.name.Data(), saveDirExt.Data()), plot.axes_labels, plot.nbins * 10, plot.xlow, plot.xhigh);

  rockboxEvents->Draw(Form(plot.var + ">>Sig_%s_%s", plot.name.Data(), saveDirExt.Data()), signal_def.cut + base_cut.cut);
  rockboxEvents->Draw(Form(plot.var + ">>BackNu_%s_%s", plot.name.Data(), saveDirExt.Data()), !signal_def.cut + base_cut.cut);
  intimeEvents->Draw(Form(plot.var + ">>BackCos_%s_%s", plot.name.Data(), saveDirExt.Data()), !signal_def.cut + base_cut.cut);

  hSig->Scale(rockboxScaling);
  hBackNu->Scale(rockboxScaling);
  hBackCos->Scale(intimeScaling);

  TH1F* hBack = (TH1F*) hBackNu->Clone();
  hBack->Add(hBackCos);

  const double totalSignal = rockboxScaling * rockboxEvents->Draw("", signal_def.cut);
  const double initSignal  = rockboxScaling * rockboxEvents->Draw("", signal_def.cut + base_cut.cut);
  const double initBack    = rockboxScaling * rockboxEvents->Draw("", !signal_def.cut + base_cut.cut) + intimeScaling * intimeEvents->Draw("", !signal_def.cut + base_cut.cut);

  TGraph *gSelEff = new TGraph();
  gSelEff->SetMarkerColor(kRed+2);
  TGraph *gPur = new TGraph();
  gPur->SetMarkerColor(kBlue+2);
  TGraph *gEP  = new TGraph();
  gEP->SetMarkerColor(kMagenta+2);
  TGraph *gEPP = new TGraph();
  gEPP->SetMarkerColor(kGreen+2);

  gSelEff->SetTitle(TString(hSig->GetXaxis()->GetTitle()) + ";Cut Value;Fraction");

  float accumSignal = initSignal, accumBack = initBack,
    maxEP = -std::numeric_limits<float>::max(),
    maxEPP = -std::numeric_limits<float>::max(),
    maxE = -std::numeric_limits<float>::max(),
    maxP = -std::numeric_limits<float>::max(),
    optimal_cut = -std::numeric_limits<float>::max(),
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
      const float epp = eff * pur * pur / 10000;
      const float cut = rej_cut ? hSig->GetBinLowEdge(i) : hSig->GetBinLowEdge(i) + hSig->GetBinWidth(i);

      if((!optimiseEPP && ep > maxEP) || (optimiseEPP && epp > maxEPP))
        {
          maxEP  = ep;
          maxEPP = epp;
          maxE   = eff;
          maxP   = pur;

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
      gEPP->SetPoint(i, cut, eff * pur * pur / 10000);

      if(rej_cut)
        --i;
      else
        ++i;
    }

  std::cout << "Optimal Cut: " << optimal_cut << " (Eff x Pur: " << maxEP << ")" << std::endl;

  TLine *line = new TLine(optimal_cut, 0, optimal_cut, 110);
  line->SetLineStyle(9);
  line->SetLineWidth(3);

  TPaveText *pt = new TPaveText(.25, .62, .5, .75, "NDC");
  pt->AddText(Form("Optimal Cut: %.4f", optimal_cut));
  pt->AddText(Form("Efficiency x Purity: %.2f %%", maxEP));
  if(optimiseEPP)
    pt->AddText(Form("Efficiency x Purity^2: %.2f %%", maxEPP));
  pt->AddText(Form("Efficiency: %.2f %%", maxE));
  pt->AddText(Form("Purity: %.2f %%", maxP));
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.03);
  pt->SetFillColor(kWhite);
  pt->SetBorderSize(0);
  pt->SetLineWidth(0);
  pt->SetTextAlign(11);

  TLegend *leg = new TLegend(.25, .77, .85, .82);
  leg->AddEntry(gSelEff, "Efficiency", "p");
  leg->AddEntry(gPur, "Purity", "p");
  leg->AddEntry(gEP, "Efficiency x Purity", "p");
  if(optimiseEPP)
    leg->AddEntry(gEPP, "Efficiency x Purity^2", "p");
  leg->SetNColumns(4);
  leg->SetTextSize(0.03);

  TCanvas *c = new TCanvas("c", "c");
  c->SetTopMargin(.15);
  gSelEff->Draw("AP");
  gSelEff->SetMinimum(0);
  if(max>80)
    gSelEff->SetMaximum(130);
  else
    gSelEff->SetMaximum(110);
  gPur->Draw("Psame");
  gEP->Draw("Psame");
  if(optimiseEPP)
    gEPP->Draw("Psame");
  line->Draw();
  pt->Draw();
  leg->Draw();

  TString name = plot.name;
  if(rej_cut)
    name += "_reverse";

  if(optimiseEPP)
    name += "_max_epp";

  c->SaveAs(saveDir + "/" + name + ".png");
  c->SaveAs(saveDir + "/" + name + ".pdf");
}
