#include "XSecCommon.C"
#include "LatexHeaders.h"

void MakePlot(const int type, const std::vector<int> &plot_types, const Selections &selections, const TString &saveDir,
              const std::string &weightName = "", const int &nunivs = 0,
              const std::vector<std::string> &weightNames = std::vector<std::string>());

void MakeCorrelationMatrix(const Selections &selections, const TString &saveDir, const std::string &weightName = "");

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir);

void MakeSystSummaryPlot(const Selections &selections, const TString &saveDir, const std::string &weightName,
                         const std::vector<Syst> &systs);

void XSec1D(const TString &productionVersion, const TString &saveDirExt, const int var)
{
  gROOT->SetStyle("henrySBND");
  gStyle->SetPaintTextFormat("1.2g");
  gROOT->ForceStyle();

  const std::vector<int> plot_types = { 0, 1, 2, 3 };

  TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_one_d/" + saveDirExt;

  if(var == 0)
    saveDir += "/pizero_mom";
  else if(var == 1)
    saveDir += "/cos_theta_pizero";
  else
    throw std::runtime_error(Form("Variable %i not available", var));

  gSystem->Exec("mkdir -p " + saveDir);

  XSecSamples samples = SetupSamples(productionVersion);

  const std::vector<double> piZeroMomBins = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const std::vector<double> cosThetaPiZeroBins = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  const TString pizeroMomAxes      = "p_{#pi^{0}} (MeV/c);#frac{d#sigma}{dp_{#pi^{0}}} (#frac{cm^{2}}{MeV/c nucleon})";
  const TString cosThetaPiZeroAxes = "cos(#theta_{#pi^{0}});#frac{d#sigma}{dcos(#theta_{#pi^{0}})} (#frac{cm^{2}}{nucleon})";

  const std::vector<double> bins = var == 0 ? piZeroMomBins : cosThetaPiZeroBins;
  const TString axes             = var == 0 ? pizeroMomAxes : cosThetaPiZeroAxes;
  const std::string varName      = var == 0 ? "pizero_mom" : "cos_theta_pizero";

  const int N = bins.size() - 1;

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", ("NC 1#pi^{0};" + axes).Data(),
                                      N, bins, varName, nTargets, intFlux);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", ("NC 1#pi^{0}0p0#pi^{#pm};" + axes).Data(),
                                      N, bins, varName, nTargets, intFlux);
  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi", ("NC 1#pi^{0}Np0#pi^{#pm};" + axes).Data(),
                                      N, bins, varName, nTargets, intFlux);

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

          MakeCorrelationMatrix(selections, saveDir, name);
        }

      if(weightSet.name == "genie")
        MakePlot(3, plot_types, selections, saveDir, weightSet.name + "_all", 0, weightSet.list);
    }

  const std::vector<std::string> all_systs_list = SystSetToWeightList(all_systs);

  MakePlot(3, plot_types, selections, saveDir, "all", 0, all_systs_list);
  MakeSystSummaryPlot(selections, saveDir, "all", all_systs);
}

void MakePlot(const int type, const std::vector<int> &plot_types, const Selections &selections, const TString &saveDir,
              const std::string &weightName, const int &nunivs, const std::vector<std::string> &weightNames)
{
  for(auto const& plot_type : plot_types)
    {
      for(auto&& [ selection_i, selection ] : enumerate(selections))
        {
          TCanvas *canvas = new TCanvas(Form("canvas_%s_%d", selection.name.Data(), plot_type),
                                        Form("canvas_%s_%d", selection.name.Data(), plot_type));
          canvas->cd();

          const TString saveSubDir = saveDir + "/" + selection.name + "/" + PlotTypeMap.at(plot_type);
          gSystem->Exec("mkdir -p " + saveSubDir);

          gPad->SetTopMargin(0.12);
          gPad->SetLeftMargin(0.2);
          gPad->SetRightMargin(0.1);

          TH1F *hist;

          switch(plot_type)
            {
            case 0:
              hist = selection.plot->GetNominalXSecHist(type == 0, false, selection.name);
              break;
            case 1:
              hist = selection.plot->GetNominalEfficiencyHist(type == 0);
              break;
            case 2:
              hist = selection.plot->GetNominalPurityHist(type == 0);
              break;
            case 3:
              hist = selection.plot->GetNominalBkgdCountHist(type == 0);
              break;
            default:
              std::runtime_error(Form("What is variable %d ???", plot_type));
            }

          hist->GetYaxis()->SetTitleOffset(1.5);
          hist->GetYaxis()->SetTitleSize(0.05);

          if(plot_type != 0)
            hist->GetYaxis()->SetTitle(PlotAxisMap.at(plot_type));

          if(plot_type == 3 && selection.plot->GetVar() == "pizero_mom")
            hist->GetYaxis()->SetTitle(PlotAxisMap.at(plot_type) + " / MeV");
          else if(plot_type == 3)
            hist->GetYaxis()->SetTitle(PlotAxisMap.at(plot_type) + " / 1");

          hist->Draw("histe][");
          hist->SetMinimum(0);
          hist->SetMaximum(1.25 * hist->GetMaximum());
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
                      unihist = selection.plot->GetUniverseBkgdCountHist(weightName, univ);
                      break;
                    default:
                      std::runtime_error(Form("What is variable %d ???", plot_type));
                    }

                  unihist->Draw("hist][same");
                }

              hist->Draw("histesame][");
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
                  graph = selection.plot->GetCVErrBkgdCountGraph(weightName);
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
                  graph = selection.plot->GetCVErrBkgdCountGraph(weightName);
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
    
          AddText(canvas, wip, kGray+2, {.8, .895, .91, .905}, 0.025, 32);

          if(type == 0)
            {
              canvas->SaveAs(saveSubDir + "/" + PlotTypeMap.at(plot_type) + "_" + selection.name + "_nominal.png");
              canvas->SaveAs(saveSubDir + "/" + PlotTypeMap.at(plot_type) + "_" + selection.name + "_nominal.pdf");
            }
          else if(type == 1)
            {
              canvas->SaveAs(saveSubDir + "/" + PlotTypeMap.at(plot_type) + "_" + selection.name + "_" + weightName.c_str() + "_univs.png");
              canvas->SaveAs(saveSubDir + "/" + PlotTypeMap.at(plot_type) + "_" + selection.name + "_" + weightName.c_str() + "_univs.pdf");
            }
          else if(type == 2 || type == 3)
            {
              canvas->SaveAs(saveSubDir + "/" + PlotTypeMap.at(plot_type) + "_" + selection.name + "_" + weightName.c_str() + "_cv_err.png");
              canvas->SaveAs(saveSubDir + "/" + PlotTypeMap.at(plot_type) + "_" + selection.name + "_" + weightName.c_str() + "_cv_err.pdf");
            }

          delete hist;
          delete canvas;
        }
    }
}

void MakeCorrelationMatrix(const Selections &selections, const TString &saveDir, const std::string &weightName)
{
  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      TCanvas *canvas = new TCanvas(Form("matrix_canvas%s", selection.name.Data()),
                                    Form("matrix_canvas%s", selection.name.Data()));
      canvas->cd();

      const TString saveSubDir = saveDir + "/" + selection.name + "/correlation_matrices";
      gSystem->Exec("mkdir -p " + saveSubDir);

      gPad->SetTopMargin(0.12);
      gPad->SetRightMargin(0.2);

      TH2D *frac_cov_matrix = selection.plot->CreateFractionalCovarianceMatrix(weightName);
      frac_cov_matrix->Draw("colztext");

      canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_frac_cov_matrix.png");
      canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_frac_cov_matrix.pdf");

      delete frac_cov_matrix;

      TH2D *corr_matrix = selection.plot->CreateCorrelationMatrix(weightName);
      corr_matrix->Draw("colztext");

      canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_corr_matrix.png");
      canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_corr_matrix.pdf");

      delete corr_matrix;
    }
}

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir)
{
  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      TCanvas *canvas = new TCanvas(Form("canvas%s", selection.name.Data()),
                                    Form("canvas%s", selection.name.Data()));
      canvas->cd();

      const TString saveSubDir = saveDir + "/" + selection.name;
      gSystem->Exec("mkdir -p " + saveSubDir);

      gPad->SetTopMargin(0.12);
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.1);

      TH1F *hist = selection.plot->GetNominalXSecHist(type == 0, false, selection.name);
      hist->GetYaxis()->SetTitleOffset(1.5);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(1.25 * hist->GetMaximum());
      gPad->Update();

      TH1F *geniePred = selection.plot->GetPredictedXSecHist(selection.name, "genie", 1e-38, true);
      geniePred->SetLineColor(kOrange+2);
      geniePred->SetMarkerStyle(1);
      geniePred->Draw("histeqsame");

      TH1F *nuwroPred = selection.plot->GetPredictedXSecHist(selection.name, "nuwro", 1e-38, true);
      nuwroPred->SetLineColor(kGreen+2);
      nuwroPred->SetMarkerStyle(1);
      nuwroPred->Draw("histeqsame");

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      AddText(canvas, wip, kGray+2, {.8, .895, .91, .905}, 0.025, 32);

      TLegend *leg = new TLegend(.35, .6, .7, .8);
      leg->AddEntry(hist, "Nominal + Stat", "le");
      leg->AddEntry(geniePred, "GENIEv3 AR23_20i_00_000", "l");
      leg->AddEntry(nuwroPred, "NuWro", "l");
      leg->Draw();

      if(type == 0)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare.pdf");
        }

      delete hist;
      delete canvas;
    }
}

void MakeSystSummaryPlot(const Selections &selections, const TString &saveDir, const std::string &weightName,
                         const std::vector<Syst> &systs)
{
  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      TCanvas *canvas = new TCanvas(Form("canvas%s", selection.name.Data()),
                                    Form("canvas%s", selection.name.Data()));
      canvas->cd();

      TLegend *leg = new TLegend(.6, .65, .8, .85);
      leg->SetTextSize(0.02);

      const TString saveSubDir = saveDir + "/" + selection.name;
      gSystem->Exec("mkdir -p " + saveSubDir);

      TPad *p1 = new TPad(Form("p1%s", selection.name.Data()), Form("p1%s", selection.name.Data()),
                          0., .3, 1., 1.);
      p1->Draw();
      p1->SetLeftMargin(.2);
      p1->SetBottomMargin(0.05);
      p1->SetTopMargin(0.12);
      p1->SetRightMargin(0.05);
      p1->cd();

      TH1F *hist = selection.plot->GetNominalXSecHist(false, false, selection.name);
      hist->GetYaxis()->SetTitleOffset(1.5);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(1.25 * hist->GetMaximum());
      gPad->Update();

      const std::vector<std::string> systs_list = SystSetToWeightList(systs);
      selection.plot->CombineErrorsInQuaderature(systs_list, weightName);
      TGraphAsymmErrors *graph = selection.plot->GetCVErrXSecGraph(weightName);
      graph->Draw("PEsame");

      leg->AddEntry(graph, "#splitline{MC XSec +}{Combined Systematics}", "le");

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      canvas->cd(0);

      TPad *p2 = new TPad(Form("p2%s", selection.name.Data()), Form("p2%s", selection.name.Data()),
                          0., 0., 1., .31);
      p2->Draw();
      p2->SetLeftMargin(.2);
      p2->SetBottomMargin(0.25);
      p2->SetTopMargin(0.04);
      p2->SetRightMargin(0.05);
      p2->cd();

      for(auto const& syst : systs)
        {
          TH1F *systHist = selection.plot->GetFracErrorXSecHist(syst.name);
          systHist->SetBit(TH1::kNoTitle);
          systHist->GetYaxis()->SetTitle("#splitline{Fractional}{Uncertainty}");
          systHist->GetYaxis()->SetTitleOffset(0.5);
          systHist->GetYaxis()->SetTitleSize(0.1);
          systHist->GetYaxis()->SetLabelSize(0.09);
          systHist->GetXaxis()->SetLabelSize(0.11);
          systHist->GetXaxis()->SetTitleOffset(1);
          systHist->GetXaxis()->SetTitleSize(0.11);
          systHist->SetLineColor(syst.colour);
          systHist->Draw("hist][same");
          systHist->SetMinimum(0);
          systHist->SetMaximum(0.28);
          gPad->Update();

          leg->AddEntry(systHist, syst.printName.c_str(), "l");
        }

      canvas->cd(0);
      leg->Draw();

      canvas->SaveAs(saveSubDir + "/" + selection.name + "_all_systs.png");
      canvas->SaveAs(saveSubDir + "/" + selection.name + "_all_systs.pdf");

      delete hist;
      delete canvas;
    }
}
