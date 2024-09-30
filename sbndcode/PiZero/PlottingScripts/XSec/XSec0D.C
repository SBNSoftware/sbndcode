#include "XSecCommon.C"
#include "LatexHeaders.h"

void MakePlot(const int type, const Selections &selections, const TString &saveDir,
              const std::string &weightName = "", const int &nunivs = 0,
              const std::vector<std::string> &weightNames = std::vector<std::string>());

void MakeTables(Selections &selections, WeightSets &weightSets, const TString &saveDir);

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir);

void MakeSystSummaryPlot(const Selections &selections, const TString &saveDir, const std::string &weightName,
                         const std::vector<Syst> &systs);

void XSec0D(const TString &productionVersion, const TString &saveDirExt)
{
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_zero_d/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  XSecSamples samples = SetupSamples(productionVersion);

  std::vector<double> piZeroMomBins      = { -10., 1e10 };
  std::vector<double> cosThetaPiZeroBins = { -10., 1e10 };

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", "NC 1#pi^{0};;#sigma (cm^{2}/nucleon)",
                                      1, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", "NC 1#pi^{0}0p0#pi^{#pm};;#sigma (cm^{2}/nucleon)",
                                      1, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);

  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi", "NC 1#pi^{0}Np0#pi^{#pm};;#sigma (cm^{2}/nucleon)",
                                      1, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);

  selections[0].plot = xsec_incl;
  selections[1].plot = xsec_0p0pi;
  selections[2].plot = xsec_Np0pi;

  FillPlots(samples, selections, weightSets);

  MakePlot(0, selections, saveDir);
  MakeSummaryPlot(0, selections, saveDir);

  for(WeightSet &weightSet : weightSets)
    {
      for(std::string &name : weightSet.list)
        {
          MakePlot(1, selections, saveDir, name, weightSet.nunivs);
          MakePlot(2, selections, saveDir, name, weightSet.nunivs);
        }

      if(weightSet.name == "genie")
        MakePlot(3, selections, saveDir, weightSet.name + "_all", 0, weightSet.list);
    }

  const std::vector<std::string> all_systs_list = SystSetToWeightList(all_systs);

  MakePlot(3, selections, saveDir, "all", 0, all_systs_list);
  MakeSystSummaryPlot(selections, saveDir, "all", all_systs);

  MakeTables(selections, weightSets, saveDir);
}

void MakePlot(const int type, const Selections &selections, const TString &saveDir,
              const std::string &weightName, const int &nunivs, const std::vector<std::string> &weightNames)
{
  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();
  canvas->Divide(3, 1);

  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      canvas->cd(selection_i + 1);
      gPad->SetBottomMargin(0.05);
      gPad->SetTopMargin(0.1);
      gPad->SetLeftMargin(0.25);
      gPad->SetRightMargin(0.05);

      TH1F *hist = selection.plot->GetNominalHist0D(type == 0);
      hist->GetYaxis()->SetTitleOffset(1.9);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(5e-40);
      gPad->Update();

      if(type == 1)
        {
          for(int univ = 0; univ < nunivs; ++univ)
            {
              TH1F *unihist = selection.plot->GetUniverseHist0D(weightName, univ);
              unihist->Draw("hist][same");
            }

          hist->Draw("histe][same");
        }

      if(type == 2)
        {
          TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightName);
          graph->Draw("PEsame");
        }

      if(type == 3)
        {
          selection.plot->CombineErrorsInQuaderature(weightNames, weightName);
          TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightName);
          graph->Draw("PEsame");
        }

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();
    }

  if(type == 0)
    {
      canvas->SaveAs(saveDir + "/nominal.png");
      canvas->SaveAs(saveDir + "/nominal.pdf");
    }
  else if(type == 1)
    {
      canvas->SaveAs(saveDir + "/" + weightName.c_str() + "_univs.png");
      canvas->SaveAs(saveDir + "/" + weightName.c_str() + "_univs.pdf");
    }
  else if(type == 2 || type == 3)
    {
      canvas->SaveAs(saveDir + "/" + weightName.c_str() + "_cv_err.png");
      canvas->SaveAs(saveDir + "/" + weightName.c_str() + "_cv_err.pdf");
    }
}

void MakeTables(Selections &selections, WeightSets &weightSets, const TString &saveDir)
{
  for(WeightSet &weightSet : weightSets)
    {
      ofstream texFile;
      texFile.open(saveDir + "/" + weightSet.name.c_str() + "_systematic_fractional_errors.tex");
      texFile << docStart
              << tableStart
              << "\\hline\n"
              << "Systematic & \\multicolumn{2}{|c|}{NC1$\\pi^{0}$} & \\multicolumn{2}{|c|}{NC1$\\pi^{0}$0p0$\\pi^{\\p\
m}$} & \\multicolumn{2}{|c|}{NC1$\\pi^{0}$Np0$\\pi^{\\pm}$} \\\\ \\hline"
              << "& Resolution (\\%) & Bias (\\%) & Resolution (\\%) & Bias (\\%) & Resolution (\\%) & Bias (\\%)"
              << "\\\\ \\hline" << std::endl;

      for(std::string &name : weightSet.list)
        {
          texFile << "\\url{" << name << "}";

          for(Selection &selection : selections)
            {
              XSecPlot *plot = selection.plot;
              Bin *bin = plot->GetBins()[0];

              texFile << Form(" & %.2f & %.2f", 100. * bin->GetFracSystResAve(name),
                              100. * bin->GetFracSystBias(name));
            }

          texFile << "\\\\ \\hline" << std::endl;
        }

      texFile << "all";

      for(Selection &selection : selections)
        {
          XSecPlot *plot = selection.plot;
          Bin *bin = plot->GetBins()[0];

          texFile << Form(" & %.2f & %.2f", 100. * bin->GetFracSystResAveQuadSum(weightSet.list),
                          100. * bin->GetFracSystBiasQuadSum(weightSet.list));
        }

      texFile << "\\\\ \\hline";

      texFile << tableEnd << docEnd;
      texFile.close();
      gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/" + weightSet.name.c_str() + "_systematic_fractional_errors.tex");
    }
}

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir)
{
  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();
  canvas->Divide(3, 1);

  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      canvas->cd(selection_i + 1);
      gPad->SetBottomMargin(0.05);
      gPad->SetTopMargin(0.1);
      gPad->SetLeftMargin(0.25);
      gPad->SetRightMargin(0.05);

      TH1F *hist = selection.plot->GetNominalHist0D(type == 0);
      hist->GetYaxis()->SetTitleOffset(1.9);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(5e-40);
      gPad->Update();

      TH1F *geniePred = selection.plot->GetPredictedHist0D(selection.name, "genie", 1e-38);
      geniePred->SetLineColor(kOrange+2);
      geniePred->Draw("hist][same");

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      if(selection_i == selections.size() - 1)
        {
          TLegend *leg = new TLegend(.35, .6, .7, .8);
          leg->AddEntry(hist, "Nominal + Stat", "le");
          leg->AddEntry(geniePred, "GENIEv3 AR23_20i_00_000", "l");
          leg->Draw();
        }
    }

  canvas->SaveAs(saveDir + "/nominal_genie_compare.png");
  canvas->SaveAs(saveDir + "/nominal_genie_compare.pdf");
}

void MakeSystSummaryPlot(const Selections &selections, const TString &saveDir, const std::string &weightName,
                         const std::vector<Syst> &systs)
{
  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();
  canvas->Divide(3, 1);

  TLegend *leg = new TLegend(.77, .65, .9, .85);
  leg->SetTextSize(0.02);

  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      canvas->cd(selection_i + 1);

      TPad *p1 = new TPad(Form("p1%s", selection.name.Data()), Form("p1%s", selection.name.Data()),
                          0., .3, 1., 1.);
      p1->Draw();
      p1->SetLeftMargin(.25);
      p1->SetBottomMargin(0.05);
      p1->SetTopMargin(0.1);
      p1->SetRightMargin(0.05);
      p1->cd();

      TH1F *hist = selection.plot->GetNominalHist0D();
      hist->GetYaxis()->SetTitleOffset(1.9);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(5e-40);
      gPad->Update();

      const std::vector<std::string> systs_list = SystSetToWeightList(systs);
      selection.plot->CombineErrorsInQuaderature(systs_list, weightName);
      TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightName);
      graph->Draw("PEsame");

      if(selection_i == 0)
        leg->AddEntry(graph, "#splitline{MC XSec +}{Combined Systematics}", "le");

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      canvas->cd(selection_i + 1);

      TPad *p2 = new TPad(Form("p2%s", selection.name.Data()), Form("p2%s", selection.name.Data()),
                          0., 0., 1., .31);
      p2->Draw();
      p2->SetLeftMargin(.25);
      p2->SetBottomMargin(0.05);
      p2->SetTopMargin(0.04);
      p2->SetRightMargin(0.05);
      p2->cd();

      for(auto const& syst : systs)
        {
          TH1F *systHist = selection.plot->GetFracErrorHist0D(syst.name);
          systHist->SetTitle(";;Fractional Uncertainty");
          systHist->GetYaxis()->SetTitleOffset(1.2);
          systHist->GetXaxis()->SetLabelSize(0);
          systHist->GetXaxis()->SetLabelOffset(999);
          systHist->SetLineColor(syst.colour);
          systHist->Draw("hist][same");
          systHist->SetMinimum(0);
          systHist->SetMaximum(0.25);
          gPad->Update();

          if(selection_i == 0)
            leg->AddEntry(systHist, syst.printName.c_str(), "l");
        }
    }

  canvas->cd(0);
  leg->Draw();

  canvas->SaveAs(saveDir + "/all_systs.png");
  canvas->SaveAs(saveDir + "/all_systs.pdf");
}
