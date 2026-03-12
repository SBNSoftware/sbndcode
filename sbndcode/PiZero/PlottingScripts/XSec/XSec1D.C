#include "XSecCommon.C"
#include "LatexHeaders.h"

void MakePlot(const int type, const std::vector<int> &plot_types, const Selections &selections, const TString &saveDir,
              const std::string &weightName = "", const int &nunivs = 0,
              const std::vector<std::string> &weightNames = {});

void MakePlot(const int &type, const Selections &selections, const TString &saveDir, const std::string &weightName,
	      const int &nunivs = 0, const std::vector<std::string> &weightNames = {});

void MakeCorrelationMatrix(const Selections &selections, const TString &saveDir, const std::string &weightName = "");

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir, const std::string &weightName = "");

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
          MakePlot(0, selections, saveDir, name, weightSet.nunivs);

          MakeCorrelationMatrix(selections, saveDir, name);
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
  MakeSummaryPlot(2, selections, saveDir, "all");
  MakeSummaryPlot(3, selections, saveDir, "all");
  MakeSummaryPlot(4, selections, saveDir, "all");
  MakeSummaryPlot(5, selections, saveDir, "all");
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

void MakePlot(const int &type, const Selections &selections, const TString &saveDir, const std::string &weightName,
	      const int &nunivs, const std::vector<std::string> &weightNames)
{
  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      TCanvas *canvas = new TCanvas(Form("canvas_%s_vary_all", selection.name.Data()),
				    Form("canvas_%s_vary_all", selection.name.Data()));
      canvas->cd();

      const TString saveSubDir = saveDir + "/" + selection.name + "/all";
      gSystem->Exec("mkdir -p " + saveSubDir);

      gPad->SetTopMargin(0.12);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.15);

      TH1F *hist = selection.plot->GetNominalXSecHist(false, false, selection.name);
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
      purGraph->SetLineColor(kRed+2);
      purGraph->SetFillStyle(3010);
      purGraph->SetLineWidth(1);
      effGraph->SetFillColor(kBlue+2);
      effGraph->SetLineColor(kBlue+2);
      effGraph->SetFillStyle(3010);
      effGraph->SetLineWidth(1);
      backGraph->SetFillColor(kGreen+2);
      backGraph->SetLineColor(kGreen+2);
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
      const float gap = 0.6 / purMax;

      mg->SetMinimum(0);
      mg->SetMaximum(0.6);

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
				0, gap * backMax, 505, "+L");
      axis->SetTitleSize(0.04);

      std::vector<double> legPos;

      if(selection.plot->GetVar() == "pizero_mom")
	{
	  axis->SetLabelSize(0.05);
	  axis->SetTitleOffset(1.1);
	  axis->SetTitle(PlotAxisMap.at(3) + " / MeV");
	  legPos = { .56, .29, .76, .4 };
	}
      else
	{
	  axis->SetLabelSize(0.045);
	  axis->SetTitleOffset(1.45);
	  axis->SetTitle(PlotAxisMap.at(3) + " / 1");
	  legPos = { .22, .39, .42, .5 };
	}

      axis->Draw();
      
      TLegend *legend = new TLegend(legPos[0], legPos[1], legPos[2], legPos[3]);
      legend->AddEntry(effGraph, "Efficiency", "f");
      legend->AddEntry(purGraph, "Purity", "f");
      legend->AddEntry(backGraph, "Background Count", "f");
      legend->Draw();

      AddText(canvas, wip, kGray+2, { .8, .895, .85, .905 }, 0.025, 32);

      canvas->SaveAs(saveSubDir + "/varying_all_" + selection.name + "_" + weightName.c_str() + "_cv_err.png");
      canvas->SaveAs(saveSubDir + "/varying_all_" + selection.name + "_" + weightName.c_str() + "_cv_err.pdf");

      delete canvas;
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

void MakeSummaryPlot(const int type, const Selections &selections, const TString &saveDir, const std::string &weightName)
{
  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      TCanvas *canvas = new TCanvas(Form("canvas%s", selection.name.Data()),
                                    Form("canvas%s", selection.name.Data()));
      canvas->cd();

      const TString saveSubDir = saveDir + "/" + selection.name;
      gSystem->Exec("mkdir -p " + saveSubDir);

      const TString covSubDir = saveSubDir + "/covariance_matrices";
      gSystem->Exec("mkdir -p " + covSubDir);

      gPad->SetTopMargin(0.12);
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.1);

      std::vector<float> legPos = { .38, .67, .68, .82 };
      if(selection.plot->GetVar() != "pizero_mom")
	legPos = { .3, .65, .6, .8 };

      TLegend *leg = new TLegend(legPos[0], legPos[1], legPos[2], legPos[3]);

      TH1F *hist = selection.plot->GetNominalXSecHist(type == 0, false, selection.name);
      hist->GetYaxis()->SetTitleOffset(1.5);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->Draw("hist][");
      hist->SetMinimum(0);
      hist->SetMaximum(1.25 * hist->GetMaximum());
      gPad->Update();

      TH2D *cov_matrix, *genie_cov_matrix;

      unsigned syst_i = 0, genie_syst_i = 0;

      if(type != 0)
	{
	  for(WeightSet &weightSet : weightSets)
	    {
	      for(std::string &name : weightSet.list)
		{
		  TH2D *temp_cov_matrix = selection.plot->CreateCovarianceMatrix(name);
		  if(syst_i == 0)
		    cov_matrix = temp_cov_matrix;
		  else
		    cov_matrix->Add(temp_cov_matrix);

		  if(name.find("GENIEReWeight") != std::string::npos)
		    {
		      if(genie_syst_i == 0)
			genie_cov_matrix = temp_cov_matrix;
		      else
			genie_cov_matrix->Add(temp_cov_matrix);

		      ++genie_syst_i;
		    }

		  if(type == 1)
		    {
		      TCanvas *canvas_cov = new TCanvas(Form("canvas_cov%s", selection.name.Data()),
							Form("canvas_cov%s", selection.name.Data()));
		      canvas_cov->cd();

		      gPad->SetTopMargin(0.12);
		      gPad->SetRightMargin(0.2);

		      temp_cov_matrix->Draw("colztext");

		      canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_" + name.c_str() + "_cov_matrix.png");
		      canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_" + name.c_str() + "_cov_matrix.pdf");

		      TH2D *norm = NormMatrix(hist, temp_cov_matrix);
		      norm->Draw("colztext");

		      canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_" + name.c_str() + "_cov_norm_matrix.png");
		      canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_" + name.c_str() + "_cov_norm_matrix.pdf");

		      TH2D *shape = ShapeMatrix(hist, temp_cov_matrix);
		      shape->Draw("colztext");

		      canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_" + name.c_str() + "_cov_shape_matrix.png");
		      canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_" + name.c_str() + "_cov_shape_matrix.pdf");

		      delete norm;
		      delete shape;
		      delete canvas_cov;
		    }

		  ++syst_i;
		}
	    }

	  for(int i = 1; i < cov_matrix->GetNbinsX() + 1; ++i)
	    {
	      std::cout << cov_matrix->GetBinContent(i, i) << " " << std::pow(0.02 * hist->GetBinContent(i), 2) << " " << std::pow(0.01 * hist->GetBinContent(i), 2) << std::endl;
	      cov_matrix->SetBinContent(i, i, cov_matrix->GetBinContent(i, i) + std::pow(0.02 * hist->GetBinContent(i), 2)
					+ std::pow(0.01 * hist->GetBinContent(i), 2));
	    }
	}

      if(type == 1)
	{
	  TCanvas *canvas_cov = new TCanvas(Form("canvas_cov%s", selection.name.Data()),
					    Form("canvas_cov%s", selection.name.Data()));
	  canvas_cov->cd();

	  gPad->SetTopMargin(0.12);
	  gPad->SetRightMargin(0.2);

	  cov_matrix->Draw("colztext");

	  canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_all_systs_cov_matrix.png");
	  canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_all_systs_cov_matrix.pdf");

	  TH2D *norm = NormMatrix(hist, cov_matrix);
	  norm->Draw("colztext");

	  canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_all_systs_cov_norm_matrix.png");
	  canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_all_systs_cov_norm_matrix.pdf");

	  TH2D *shape = ShapeMatrix(hist, cov_matrix);
	  shape->Draw("colztext");

	  canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_all_systs_cov_shape_matrix.png");
	  canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_all_systs_cov_shape_matrix.pdf");

	  delete norm;
	  delete shape;
	  delete canvas_cov;

	  TCanvas *genie_canvas_cov = new TCanvas(Form("genie_canvas_cov%s", selection.name.Data()),
						  Form("genie_canvas_cov%s", selection.name.Data()));
	  genie_canvas_cov->cd();

	  gPad->SetTopMargin(0.12);
	  gPad->SetRightMargin(0.2);

	  genie_cov_matrix->Draw("colztext");

	  genie_canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_genie_all_cov_matrix.png");
	  genie_canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_genie_all_cov_matrix.pdf");

	  TH2D *genie_norm = NormMatrix(hist, genie_cov_matrix);
	  genie_norm->Draw("colztext");

	  genie_canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_genie_all_cov_norm_matrix.png");
	  genie_canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_genie_all_cov_norm_matrix.pdf");

	  TH2D *genie_shape = ShapeMatrix(hist, genie_cov_matrix);
	  genie_shape->Draw("colztext");

	  genie_canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_genie_all_cov_shape_matrix.png");
	  genie_canvas_cov->SaveAs(covSubDir + "/" + selection.name + "_genie_all_cov_shape_matrix.pdf");

	  delete genie_norm;
	  delete genie_shape;
	  delete genie_canvas_cov;
	}

      canvas->cd();

      if(type == 0)
	leg->AddEntry(hist, "Nominal + Stat", "le");
      else if(type == 1)
	{
	  TGraphAsymmErrors *graph = selection.plot->GetCVErrXSecGraph(weightName);
	  graph->Draw("PEsame");

	  leg->AddEntry(graph, "Nominal + Systs", "pe");
	}
      else if(type == 2 || type == 4)
	{
	  TH2D *norm = NormMatrix(hist, cov_matrix);
	  TGraphAsymmErrors *graph = selection.plot->GetXSecGraphWithMatrixErrors(hist, norm);
	  graph->SetLineColor(kRed+2);
	  graph->Draw("PEsame");

	  leg->AddEntry(graph, "Nominal + Norm Systs", "pe");
	}
      else if(type == 3 || type == 5)
	{
	  TH2D *shape = ShapeMatrix(hist, cov_matrix);
	  TGraphAsymmErrors *graph = selection.plot->GetXSecGraphWithMatrixErrors(hist, shape);
	  graph->SetLineColor(kMagenta-3);
	  graph->Draw("PEsame");

	  leg->AddEntry(graph, "Nominal + Shape Systs", "pe");
	}

      TH1F *geniePred = selection.plot->GetPredictedXSecHist(selection.name, "genie", 1e-38, true);
      geniePred->SetLineColor(kOrange+2);
      geniePred->SetMarkerStyle(1);

      if(type != 4 && type != 5)
	geniePred->Draw("histeqsame");

      TH1F *nuwroPred = selection.plot->GetPredictedXSecHist(selection.name, "nuwro", 1e-38, true);
      nuwroPred->SetLineColor(kGreen+2);
      nuwroPred->SetMarkerStyle(1);

      if(type != 4 && type != 5)
	nuwroPred->Draw("histeqsame");

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      AddText(canvas, wip, kGray+2, {.8, .895, .91, .905}, 0.025, 32);

      if(type == 1)
	{
	  const double genieChi2 = Chi2Comp(hist, geniePred, cov_matrix);
	  const double nuwroChi2 = Chi2Comp(hist, nuwroPred, cov_matrix);
	  leg->AddEntry(geniePred, Form("GENIEv3 AR23_20i_00_000 (%.2f/%d)", genieChi2, selection.plot->GetNBins()), "l");
	  leg->AddEntry(nuwroPred, Form("NuWro v21.09.2 (%.2f/%d)", nuwroChi2, selection.plot->GetNBins()), "l");
	}
      else if(type != 4 && type != 5)
	{
	  leg->AddEntry(geniePred, "GENIEv3 AR23_20i_00_000", "l");
	  leg->AddEntry(nuwroPred, "NuWro v21.09.2", "l");
	}

      leg->Draw();

      if(type == 0)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare.pdf");
        }
      else if(type == 1)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare_systs.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare_systs.pdf");
        }
      else if(type == 2)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare_systs_norm.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare_systs_norm.pdf");
        }
      else if(type == 3)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare_systs_shape.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_generator_compare_systs_shape.pdf");
        }
      else if(type == 4)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_systs_norm.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_systs_norm.pdf");
        }
      else if(type == 5)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_systs_shape.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal_systs_shape.pdf");
        }

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
