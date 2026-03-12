#include "XSecCommon.C"
#include "LatexHeaders.h"

void MakePlot(const int type, const std::vector<int> &plot_types, const Selections &selections,
              const TString &saveDir, const std::string &weightName = "", const int &nunivs = 0,
              const std::vector<std::string> &weightNames = {});

void MakePlot(const int type, const Selections &selections, const TString &saveDir, const std::string &weightName,
	      const int &nunivs = 0, const std::vector<std::string> &weightNames = {});

void MakeTables(Selections &selections, WeightSets &weightSets, const TString &saveDir);

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir, const std::string &weightName = "");

void MakeSystSummaryPlot(const Selections &selections, const TString &saveDir, const std::string &weightName,
                         const std::vector<Syst> &systs);

void XSec0D(const TString &productionVersion, const TString &saveDirExt)
{
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const std::vector<int> plot_types = { 0, 1, 2, 3 };

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_zero_d/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  for(auto const& plot_type : plot_types)
    gSystem->Exec("mkdir -p " + saveDir + "/" + PlotTypeMap.at(plot_type));

  XSecSamples samples = SetupSamples(productionVersion);

  std::vector<double> bins = { -10., 1e10 };

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", "NC 1#pi^{0};;#sigma (cm^{2}/nucleon)",
                                      1, bins, "pizero_mom", nTargets, intFlux);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", "NC 1#pi^{0}0p0#pi^{#pm};;#sigma (cm^{2}/nucleon)",
                                      1, bins, "pizero_mom", nTargets, intFlux);
  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi", "NC 1#pi^{0}Np0#pi^{#pm};;#sigma (cm^{2}/nucleon)",
                                      1, bins, "pizero_mom", nTargets, intFlux);

  selections[0].plot = xsec_incl;
  selections[1].plot = xsec_0p0pi;
  selections[2].plot = xsec_Np0pi;

  FillPlots(samples, selections, weightSets);

  MakePlot(0, plot_types, selections, saveDir);
  MakeSummaryPlot(0, selections, saveDir);

  for(WeightSet &weightSet : weightSets)
    {
      for(std::string &name : weightSet.list)
        {
          MakePlot(1, plot_types, selections, saveDir, name, weightSet.nunivs);
          MakePlot(2, plot_types, selections, saveDir, name, weightSet.nunivs);
          MakePlot(0, selections, saveDir, name, weightSet.nunivs);
        }

      if(weightSet.name == "genie")
	{
	  MakePlot(3, plot_types, selections, saveDir, weightSet.name + "_all", 0, weightSet.list);
	  MakePlot(1, selections, saveDir, weightSet.name + "_all", 0, weightSet.list);
	}
    }

  const std::vector<std::string> all_systs_list = SystSetToWeightList(all_systs);

  MakePlot(3, plot_types, selections, saveDir, "all", 0, all_systs_list);
  MakePlot(1, selections, saveDir, "all", 0, all_systs_list);
  MakeSystSummaryPlot(selections, saveDir, "all", all_systs);
  MakeSummaryPlot(1, selections, saveDir, "all");

  MakeTables(selections, weightSets, saveDir);
}

void MakePlot(const int type, const std::vector<int> &plot_types, const Selections &selections, const TString &saveDir,
              const std::string &weightName, const int &nunivs, const std::vector<std::string> &weightNames)
{
  for(auto const& plot_type : plot_types)
    {
      TCanvas *canvas = new TCanvas(Form("canvas_%d", plot_type),
                                    Form("canvas_%d", plot_type));
      canvas->cd();
      canvas->Divide(3, 1);

      for(auto&& [ selection_i, selection ] : enumerate(selections))
        {
          canvas->cd(selection_i + 1);
          gPad->SetBottomMargin(0.05);
          gPad->SetTopMargin(0.1);
          gPad->SetLeftMargin(0.25);
          gPad->SetRightMargin(0.05);

          TH1F *hist;

          switch(plot_type)
            {
            case 0:
              hist = selection.plot->GetNominalXSecHist(type == 0);
              hist->SetMaximum(5e-40);
              break;
            case 1:
              hist = selection.plot->GetNominalEfficiencyHist(type == 0);
              hist->SetMaximum(0.5);
              break;
            case 2:
              hist = selection.plot->GetNominalPurityHist(type == 0);
              hist->SetMaximum(0.6);
              break;
            case 3:
              hist = selection.plot->GetNominalBkgdCountHist(type == 0, true);
              hist->SetMaximum(2e5);
              break;
            default:
              std::runtime_error(Form("What is variable %d ???", plot_type));
            }

          hist->GetYaxis()->SetTitleOffset(1.9);
          hist->GetXaxis()->SetLabelSize(0);
          hist->GetXaxis()->SetLabelOffset(999);

          if(plot_type != 0)
            hist->GetYaxis()->SetTitle(PlotAxisMap.at(plot_type));

          hist->Draw("histe][");
          hist->SetMinimum(0);
          gPad->Update();

          if(type == 1)
            {
              for(int univ = 0; univ < nunivs; ++univ)
                {
                  TH1F *unihist;

                  switch(plot_type)
                    {
                    case 0:
                      unihist = selection.plot->GetUniverseXSecHist(weightName, univ);
                      break;
                    case 1:
                      unihist = selection.plot->GetUniverseEfficiencyHist(weightName, univ);
                      break;
                    case 2:
                      unihist = selection.plot->GetUniversePurityHist(weightName, univ);
                      break;
                    case 3:
                      unihist = selection.plot->GetUniverseBkgdCountHist(weightName, univ, true);
                      break;
                    default:
                      std::runtime_error(Form("What is variable %d ???", plot_type));
                    }

                  unihist->Draw("hist][same");
                }

              hist->Draw("histe][same");
            }

          if(type == 2)
            {
              TGraphAsymmErrors *graph;

              switch(plot_type)
                {
                case 0:
                  graph = selection.plot->GetCVErrXSecGraph(weightName);
                  break;
                case 1:
                  graph = selection.plot->GetCVErrEfficiencyGraph(weightName);
                  break;
                case 2:
                  graph = selection.plot->GetCVErrPurityGraph(weightName);
                  break;
                case 3:
                  graph = selection.plot->GetCVErrBkgdCountGraph(weightName, true);
                  break;
                default:
                  std::runtime_error(Form("What is variable %d ???", plot_type));
                }

              graph->Draw("PEsame");
            }

          if(type == 3)
            {
              selection.plot->CombineErrorsInQuaderature(weightNames, weightName);
              TGraphAsymmErrors *graph;

              switch(plot_type)
                {
                case 0:
                  graph = selection.plot->GetCVErrXSecGraph(weightName);
                  break;
                case 1:
                  graph = selection.plot->GetCVErrEfficiencyGraph(weightName);
                  break;
                case 2:
                  graph = selection.plot->GetCVErrPurityGraph(weightName);
                  break;
                case 3:
                  graph = selection.plot->GetCVErrBkgdCountGraph(weightName, true);
                  break;
                default:
                  std::runtime_error(Form("What is variable %d ???", plot_type));
                }

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
          canvas->SaveAs(saveDir + "/" + PlotTypeMap.at(plot_type) + "/" + PlotTypeMap.at(plot_type) + "_nominal.png");
          canvas->SaveAs(saveDir + "/" + PlotTypeMap.at(plot_type) + "/" + PlotTypeMap.at(plot_type) + "_nominal.pdf");
        }
      else if(type == 1)
        {
          canvas->SaveAs(saveDir + "/" + PlotTypeMap.at(plot_type) + "/" + PlotTypeMap.at(plot_type) + "_" + weightName.c_str() + "_univs.png");
          canvas->SaveAs(saveDir + "/" + PlotTypeMap.at(plot_type) + "/" + PlotTypeMap.at(plot_type) + "_" + weightName.c_str() + "_univs.pdf");
        }
      else if(type == 2 || type == 3)
        {
          canvas->SaveAs(saveDir + "/" + PlotTypeMap.at(plot_type) + "/" + PlotTypeMap.at(plot_type) +  "_" + weightName.c_str() + "_cv_err.png");
          canvas->SaveAs(saveDir + "/" + PlotTypeMap.at(plot_type) + "/" + PlotTypeMap.at(plot_type) +  "_" + weightName.c_str() + "_cv_err.pdf");
        }
    }
}

void MakePlot(const int type, const Selections &selections, const TString &saveDir, const std::string &weightName,
	      const int &nunivs, const std::vector<std::string> &weightNames)
{
  TCanvas *canvas = new TCanvas("canvas_vary_all", "canvas_vary_all");
  canvas->cd();
  canvas->Divide(3, 1);

  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      canvas->cd(selection_i + 1);
      gPad->SetBottomMargin(0.05);
      gPad->SetTopMargin(0.1);
      gPad->SetLeftMargin(0.25);
      gPad->SetRightMargin(0.05);

      TH1F *hist = selection.plot->GetNominalXSecHist(false);
      hist->Draw("histe][");
      gPad->Modified();
      gPad->Update();

      if(type == 1)
        selection.plot->CombineErrorsInQuaderature(weightNames, weightName);

      TGraphAsymmErrors *effGraph  = selection.plot->GetCVErrEfficiencyGraph(weightName);
      TGraphAsymmErrors *purGraph  = selection.plot->GetCVErrPurityGraph(weightName);
      TGraphAsymmErrors *backGraph = selection.plot->GetCVErrBkgdCountGraph(weightName);

      backGraph->Draw();
      const float backMax = backGraph->GetYaxis()->GetXmax();
      purGraph->Draw();
      const float purMax = purGraph->GetYaxis()->GetXmax();

      purGraph->SetFillColor(kRed+2);
      purGraph->SetFillStyle(3010);
      purGraph->SetLineWidth(1);
      effGraph->SetFillColor(kBlue+2);
      effGraph->SetFillStyle(3010);
      effGraph->SetLineWidth(1);
      backGraph->SetFillColor(kGreen+2);
      backGraph->SetFillStyle(3010);
      backGraph->SetLineWidth(1);

      TMultiGraph *mg = new TMultiGraph();
      mg->Add(purGraph);
      mg->Add(effGraph);
      mg->Add(backGraph);
      mg->Draw("PE3");

      gPad->Modified();
      gPad->Update();

      const float scale = purMax / backMax;
      backGraph->Scale(scale);


      mg->Draw("PE3same");
      mg->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
      mg->GetXaxis()->SetTitleOffset(1.1);
      mg->GetYaxis()->SetTitle("Fraction");
      mg->GetXaxis()->SetLabelSize(0.05);
      mg->GetYaxis()->SetLabelSize(0.05);

      gPad->Modified();
      gPad->Update();
      gPad->SetTicky(0);

      TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                gPad->GetUxmax(), gPad->GetUymax(),
                                0, backMax, 505, "+L");
      axis->SetTitleOffset(1.1);
      axis->SetLabelSize(0.05);
      axis->SetTitleSize(0.04);

      if(selection.plot->GetVar() == "pizero_mom")
        axis->SetTitle(PlotAxisMap.at(3) + " / MeV");
      else
        axis->SetTitle(PlotAxisMap.at(3) + " / 1");
      axis->Draw();

      AddText(canvas, wip, kGray+2, { .8, .895, .85, .905 }, 0.025, 32);
    }

  gSystem->Exec("mkdir -p " + saveDir + "/all");
  
  canvas->SaveAs(saveDir + "/all/varying_all_" + weightName.c_str() + "_cv_err.png");
  canvas->SaveAs(saveDir + "/all/varying_all_" + weightName.c_str() + "_cv_err.pdf");
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

              texFile << Form(" & %.2f & %.2f", 100. * bin->GetFracSystResAveXSec(name),
                              100. * bin->GetFracSystBiasXSec(name));
            }

          texFile << "\\\\ \\hline" << std::endl;
        }

      texFile << "all";

      for(Selection &selection : selections)
        {
          XSecPlot *plot = selection.plot;
          Bin *bin = plot->GetBins()[0];

          texFile << Form(" & %.2f & %.2f", 100. * bin->GetFracSystResAveQuadSumXSec(weightSet.list),
                          100. * bin->GetFracSystBiasQuadSumXSec(weightSet.list));
        }

      texFile << "\\\\ \\hline";

      texFile << tableEnd << docEnd;
      texFile.close();
      gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/" + weightSet.name.c_str() + "_systematic_fractional_errors.tex");
    }
}

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir, const std::string &weightName)
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

      TLegend *leg = new TLegend(.3, .65, .7, .8);

      TH1F *hist = selection.plot->GetNominalXSecHist(type == 0);
      hist->GetYaxis()->SetTitleOffset(1.9);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("hist][");
      hist->SetMinimum(0);
      hist->SetMaximum(5e-40);
      gPad->Update();

      if(type == 0)
	leg->AddEntry(hist, "Nominal + Stat", "le");
      else if(type == 1)
	{
	  TGraphAsymmErrors *graph = selection.plot->GetCVErrXSecGraph(weightName);
	  graph->Draw("PEsame");

	  leg->AddEntry(graph, "Nominal + Systs", "pe");
	}

      TH1F *geniePred = selection.plot->GetPredictedXSecHist(selection.name, "genie", 1e-38);
      geniePred->SetLineColor(kOrange+2);
      geniePred->SetMarkerStyle(1);
      geniePred->Draw("histe][same");

      TH1F *nuwroPred = selection.plot->GetPredictedXSecHist(selection.name, "nuwro", 1e-38);
      nuwroPred->SetLineColor(kGreen+2);
      nuwroPred->SetMarkerStyle(1);
      nuwroPred->Draw("histe][same");

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      std::cout << selection.name << " " << geniePred->GetBinContent(1) / hist->GetBinContent(1) << std::endl;
      std::cout << "\t" << geniePred->GetBinContent(1) << std::endl;
      std::cout << "\t" << nuwroPred->GetBinContent(1) << std::endl;

      if(selection_i == selections.size() - 1)
        {
          leg->AddEntry(geniePred, "GENIEv3 AR23_20i_00_000", "l");
          leg->AddEntry(nuwroPred, "NuWro v21.09.2", "l");
	  leg->SetTextSize(0.045);
          leg->Draw();
        }
    }

  if(type == 0)
    {
      canvas->SaveAs(saveDir + "/nominal_generator_compare.png");
      canvas->SaveAs(saveDir + "/nominal_generator_compare.pdf");
    }
  else if(type == 1)
    {
      canvas->SaveAs(saveDir + "/nominal_generator_compare_systs.png");
      canvas->SaveAs(saveDir + "/nominal_generator_compare_systs.pdf");
      canvas->SaveAs(saveDir + "/nominal_generator_compare_systs.C");
    }
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

      TH1F *hist = selection.plot->GetNominalXSecHist();
      hist->GetYaxis()->SetTitleOffset(1.9);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(5e-40);
      gPad->Update();

      const std::vector<std::string> systs_list = SystSetToWeightList(systs);
      selection.plot->CombineErrorsInQuaderature(systs_list, weightName);
      TGraphAsymmErrors *graph = selection.plot->GetCVErrXSecGraph(weightName);
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
          TH1F *systHist = selection.plot->GetFracErrorXSecHist(syst.name);
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
