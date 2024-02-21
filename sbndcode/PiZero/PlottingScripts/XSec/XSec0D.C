#include "XSecCommon.C"
#include "Common.C"

void MakePlot(const Selections selections, const TString saveDir);

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

  XSecSamples samples = { { rockboxNus, rockboxSlices, rockboxScaling },
                          { intimeNus, intimeSlices, intimeScaling }
  };

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", "NC 1#pi^{0};;#sigma (cm^{2}/nucleon)", 1);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", "NC 1#pi^{0}0p0#pi^{#pm};;#sigma (cm^{2}/nucleon)", 1);
  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi", "NC 1#pi^{0}Np0#pi^{#pm};;#sigma (cm^{2}/nucleon)", 1);

  Selection ncpizero_incl  = { "ncpizero_incl", "sel_incl", "event_type_incl==0", xsec_incl };
  Selection ncpizero_0p0pi = { "ncpizero_0p0pi", "sel_0p0pi", "event_type_0p0pi==0", xsec_0p0pi };
  Selection ncpizero_Np0pi = { "ncpizero_Np0pi", "sel_Np0pi", "event_type_Np0pi==0", xsec_Np0pi };

  Selections selections = { ncpizero_incl, ncpizero_0p0pi, ncpizero_Np0pi };

  for(auto const& selection : selections)
    {
      CalculateNominal(samples, selection);
      CalculateNominalXSecPurity(selection, nTargets, intFlux);
    }

  MakePlot(selections, saveDir);
}

void MakePlot(const Selections selections, const TString saveDir)
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

      TH1F *hist = selection.plot->GetNominalHist();
      hist->GetYaxis()->SetTitleOffset(1.9);
      hist->GetXaxis()->SetLabelSize(0);
      hist->GetXaxis()->SetLabelOffset(999);
      hist->Draw("histe][");
      gPad->Update();

      TPaveText* title = (TPaveText*)gPad->FindObject("title");
      title->SetY1NDC(0.92);
      title->SetY2NDC(1);
      title->SetX1NDC(0.45);
      title->SetX2NDC(.8);
      gPad->Modified();
      gPad->Update();

      hist->SetMinimum(0);
      hist->SetMaximum(5e-40);
    }

  canvas->SaveAs(saveDir + "/nominal.png");
  canvas->SaveAs(saveDir + "/nominal.pdf");
}
