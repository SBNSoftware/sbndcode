#include "XSecCommon.C"
#include "LatexHeaders.h"

void MakePlot(const int type, const Selections selections, const TString saveDir,
              const std::string weightName = "", const int nunivs = 0,
              const std::vector<std::string> weightNames = std::vector<std::string>());

void MakeCorrelationMatrix(const Selections selections, const TString saveDir, const std::string weightName = "");

void XSec1D(const TString productionVersion, const TString saveDirExt, const int var)
{
  gROOT->SetStyle("henrySBND");
  gStyle->SetPaintTextFormat("1.2g");
  gROOT->ForceStyle();

  TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_one_d/" + saveDirExt;

  if(var == 0)
    saveDir += "/pizero_mom";
  else if(var == 1)
    saveDir += "/cos_theta_pizero";
  else
    throw std::runtime_error(Form("Variable %i not available", var));

  gSystem->Exec("mkdir -p " + saveDir);

  XSecSamples samples = SetupSamples(productionVersion);

  std::vector<double> piZeroMomBins;
  std::vector<double> cosThetaPiZeroBins;
  TString varAxis;

  if(var == 0)
    {
      piZeroMomBins      = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
      cosThetaPiZeroBins = { def_double, def_double_high };
      varAxis            = "p_{#pi^{0}};#frac{d#sigma}{dp_{#pi^{0}}} (#frac{cm^{2}}{MeV/c nucleon})";
    }
  else if(var == 1)
    {
      piZeroMomBins      = { def_double, def_double_high };
      cosThetaPiZeroBins = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };
      varAxis            = "cos(#theta_{#pi^{0}});#frac{d#sigma}{dcos(#theta_{#pi^{0}})} (#frac{cm^{2}}{nucleon})";
    }

  const int N = (piZeroMomBins.size() - 1) * (cosThetaPiZeroBins.size() - 1);

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", ("NC 1#pi^{0};" + varAxis).Data(),
                                      N, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", ("NC 1#pi^{0}0p0#pi^{#pm};" + varAxis).Data(),
                                      N, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);

  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi",
                                      ("NC 1#pi^{0}Np0#pi^{#pm};" + varAxis).Data(),
                                      N, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);

  selections[0].plot = xsec_incl;
  selections[1].plot = xsec_0p0pi;
  selections[2].plot = xsec_Np0pi;

  FillPlots(samples, selections, weightSets);

  MakePlot(0, selections, saveDir);

  std::vector<std::string> all_weights;

  for(WeightSet &weightSet : weightSets)
    {
      all_weights.insert(all_weights.end(), weightSet.list.begin(), weightSet.list.end());

      for(std::string &name : weightSet.list)
        {
          MakePlot(1, selections, saveDir, name, weightSet.nunivs);
          MakePlot(2, selections, saveDir, name, weightSet.nunivs);

          MakeCorrelationMatrix(selections, saveDir, name);
        }

      MakePlot(3, selections, saveDir, weightSet.name + "_all", 0, weightSet.list);
    }

  MakePlot(3, selections, saveDir, "all", 0, all_weights);
}

void MakePlot(const int type, const Selections selections, const TString saveDir,
              const std::string weightName = "", const int nunivs = 0,
              const std::vector<std::string> weightNames = std::vector<std::string>())
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

      TH1F *hist = selection.plot->GetNominalHist1D(type == 0);
      hist->GetYaxis()->SetTitleOffset(1.5);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->Draw("histe][");
      hist->SetMinimum(0);
      hist->SetMaximum(1.25 * hist->GetMaximum());
      gPad->Update();

      if(type == 1)
        {
          for(int univ = 0; univ < nunivs; ++univ)
            {
              TH1F *unihist = selection.plot->GetUniverseHist1D(weightName, univ);
              unihist->Draw("hist][same");
            }

          hist->Draw("histesame][");
        }

      if(type == 2)
        {
          TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightName);
          graph->Draw("PEsame");
        }

      if(type == 3)
        {
          TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightNames);
          graph->Draw("PEsame");
        }

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();
    
      if(type == 0)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_nominal.pdf");
        }
      else if(type == 1)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_univs.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_univs.pdf");
        }
      else if(type == 2 || type == 3)
        {
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_cv_err.png");
          canvas->SaveAs(saveSubDir + "/" + selection.name + "_" + weightName.c_str() + "_cv_err.pdf");
        }

      delete hist;
      delete canvas;
    }
}

void MakeCorrelationMatrix(const Selections selections, const TString saveDir, const std::string weightName = "")
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
