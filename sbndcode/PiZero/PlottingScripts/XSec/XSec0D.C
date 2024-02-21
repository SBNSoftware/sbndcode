#include "XSecCommon.C"
#include "Common.C"
#include "WeightNames.h"

void MakePlot(const int type, const Selections selections, const TString saveDir,
              const std::string weightName = "", const int nunivs = 0);

void XSec0D(const TString productionVersion, const TString saveDirExt)
{
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_zero_d/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  TChain *rockboxNus = new TChain("ncpizeroxsectrees/neutrinos");
  rockboxNus->Add(rockboxFile);
  TChain *rockboxSlices = new TChain("ncpizeroxsectrees/slices");
  rockboxSlices->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroxsectrees/subruns");
  rockboxSubruns->Add(rockboxFile);

  TChain *intimeNus = new TChain("ncpizeroxsectrees/neutrinos");
  intimeNus->Add(intimeFile);
  TChain *intimeSlices = new TChain("ncpizeroxsectrees/slices");
  intimeSlices->Add(intimeFile);
  TChain *intimeSubruns = new TChain("ncpizeroxsectrees/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  XSecSamples samples = { { "Rockbox", rockboxNus, rockboxSlices, rockboxScaling },
                          { "Intime", intimeNus, intimeSlices, intimeScaling }
  };

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", "NC 1#pi^{0};;#sigma (cm^{2}/nucleon)", 1, nTargets, intFlux);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", "NC 1#pi^{0}0p0#pi^{#pm};;#sigma (cm^{2}/nucleon)", 1, nTargets, intFlux);
  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi", "NC 1#pi^{0}Np0#pi^{#pm};;#sigma (cm^{2}/nucleon)", 1, nTargets, intFlux);

  Selection ncpizero_incl  = { "ncpizero_incl", "sel_incl", "event_type_incl", xsec_incl };
  Selection ncpizero_0p0pi = { "ncpizero_0p0pi", "sel_0p0pi", "event_type_0p0pi", xsec_0p0pi };
  Selection ncpizero_Np0pi = { "ncpizero_Np0pi", "sel_Np0pi", "event_type_Np0pi", xsec_Np0pi };

  Selections selections = { ncpizero_incl, ncpizero_0p0pi, ncpizero_Np0pi };

  WeightSets weightSets = { { "flux", tmp_flux_weight_names, 10 } };

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
