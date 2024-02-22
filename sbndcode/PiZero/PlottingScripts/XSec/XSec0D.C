#include "XSecCommon.C"

void MakePlot(const int type, const Selections selections, const TString saveDir,
              const std::string weightName = "", const int nunivs = 0);

void XSec0D(const TString productionVersion, const TString saveDirExt)
{
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_zero_d/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  XSecSamples samples = SetupSamples(productionVersion);

  std::vector<double> piZeroMomBins      = { def_double, def_double_high };
  std::vector<double> cosThetaPiZeroBins = { def_double, def_double_high };

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

  for(WeightSet &weightSet : weightSets)
    {
      for(std::string &name : weightSet.list)
        {
          MakePlot(1, selections, saveDir, name, weightSet.nunivs);
          MakePlot(2, selections, saveDir, name, weightSet.nunivs);
        }
    }
}

void MakePlot(const int type, const Selections selections, const TString saveDir,
              const std::string weightName = "", const int nunivs = 0)
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

      TH1F *hist = selection.plot->GetNominalHist(type == 0);
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
              TH1F *unihist = selection.plot->GetUniverseHist(weightName, univ);
              unihist->Draw("hist][same");
            }

          hist->Draw("histe][same");
        }

      if(type == 2)
        {
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
  else if(type == 2)
    {
      canvas->SaveAs(saveDir + "/" + weightName.c_str() + "_cv_err.png");
      canvas->SaveAs(saveDir + "/" + weightName.c_str() + "_cv_err.pdf");
    }
}
