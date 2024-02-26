#include "XSecCommon.C"
#include "LatexHeaders.h"

void MakePlot(const int type, const Selections selections, const TString saveDir,
              const std::string weightName = "", const int nunivs = 0,
              const std::vector<std::string> weightNames = std::vector<std::string>());

void XSec2D(const TString productionVersion, const TString saveDirExt)
{
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_two_d/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  XSecSamples samples = SetupSamples(productionVersion);

  std::vector<double> piZeroMomBins      = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  std::vector<double> cosThetaPiZeroBins = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };
  TString varAxis                        = "p_{#pi^{0}};#frac{d#sigma}{dp_{#pi^{0}}dcos(#theta_{#pi^{0}})} (#frac{cm^{2}}{MeV/c nucleon})";

  const int N = (piZeroMomBins.size() - 1) * (cosThetaPiZeroBins.size() - 1);

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", (";" + varAxis).Data(),
                                      N, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", (";" + varAxis).Data(),
                                      N, piZeroMomBins, cosThetaPiZeroBins,
                                      "pizero_mom", "cos_theta_pizero", nTargets, intFlux);

  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi",
                                      (";" + varAxis).Data(),
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
      canvas->Divide(3, 3);

      const TString saveSubDir = saveDir + "/" + selection.name;
      gSystem->Exec("mkdir -p " + saveSubDir);

      for(int i = 1; i < 10; ++i)
        {
          canvas->cd(i);

          gPad->SetTopMargin(0.12);
          gPad->SetLeftMargin(0.2);
          gPad->SetRightMargin(0.1);

          TH1F *hist = selection.plot->GetNominalHist2D(type == 0, i);
          TString title = selection.plot->GetVar1BinString("cos(#theta_{#pi^{0}})", i);
          hist->SetTitle(title);
          hist->GetYaxis()->SetTitleOffset(1.5);
          hist->GetYaxis()->SetTitleSize(0.05);
          hist->Draw("histe][");
          hist->SetMinimum(0);
          hist->SetMaximum(1.25 * hist->GetMaximum());
          gPad->Update();

          TPaveText* titleTxt = (TPaveText*)gPad->FindObject("title");
          titleTxt->SetY1NDC(0.9);
          titleTxt->SetY2NDC(1);
          titleTxt->SetX1NDC(0.35);
          titleTxt->SetX2NDC(.85);
          gPad->Modified();
          gPad->Update();

          if(type == 1)
            {
              for(int univ = 0; univ < nunivs; ++univ)
                {
                  TH1F *unihist = selection.plot->GetUniverseHist2D(weightName, univ, i);
                  unihist->Draw("hist][same");
                }

              hist->Draw("histesame][");
            }

          if(type == 2)
            {
              TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightName, i);
              graph->Draw("PEsame");
            }

          if(type == 3)
            {
              TGraphAsymmErrors *graph = selection.plot->GetCVErrGraph(weightNames, i);
              graph->Draw("PEsame");
            }
        }

      canvas->cd(0);

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

      delete canvas;
    }
}
