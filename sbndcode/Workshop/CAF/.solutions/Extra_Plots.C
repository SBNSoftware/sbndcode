#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#include "ExtraVars.h"
#include "ExtraCuts.h"
#include "Cuts.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"

void SetupExtraStyle()
{
  gStyle->SetLineWidth(4);
  gStyle->SetHistLineWidth(3);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(.5);
  gStyle->SetTitleOffset(.8, "y");
  gROOT->ForceStyle();
}

void Extra_Plots()
{
  SetupExtraStyle();
  const std::string inputName = "/Workshop/CAF/reco2_withmuon_10evts.flat.caf.root";
  const TString saveDir = "./plots";

  SpectrumLoader loader(inputName);
  
  Binning sliceBins = Binning::Simple(25, 0, 25);
  Spectrum *sTotalSlices = new Spectrum("sTotalSlices", sliceBins, loader, kTotalSlices, kNoSpillCut);
  Spectrum *sNuSlices = new Spectrum("sNuSlices", sliceBins, loader, kNuSlices, kNoSpillCut);

  Binning compBins = Binning::Simple(20, 0, 1.05);
  Spectrum *sCompleteness = new Spectrum("sCompleteness", compBins, loader, kCompleteness, kExactlyOneNuSlice, kIsNeutrinoCandidateSlice);

  loader.Go();
  
  sTotalSlices->OverridePOT(1e20);
  sNuSlices->OverridePOT(1e20);
  sCompleteness->OverridePOT(1e20);

  double target_pot = 1e21;

  TCanvas *cSlices = new TCanvas("cSlices", "cSlices");
  cSlices->cd();
  TH1D *hSlices = sTotalSlices->ToTH1(target_pot);
  hSlices->SetTitle(";N Slices;Events (1x10^{21} POT)");
  hSlices->SetLineColor(kGreen+1);
  hSlices->Draw("histe");
  cSlices->SaveAs(saveDir + "/total_slices.png");

  TCanvas *cSub = new TCanvas("cSub", "cSub");
  cSub->cd();
  Spectrum sClearCosmic = (*sTotalSlices) - (*sNuSlices);
  TH1D *hNuSlices = sNuSlices->ToTH1(target_pot);
  TH1D *hClearCosmic = sClearCosmic.ToTH1(target_pot);
  hNuSlices->SetTitle(";N Slices;Events (1x10^{21} POT)");
  hNuSlices->SetLineColor(kBlue+2);
  hClearCosmic->SetLineColor(kRed+2);
  hNuSlices->Draw("histe");
  hClearCosmic->Draw("histe same");
  TLegend *leg2 = new TLegend(0.5, 0.7, 0.85, 0.85);
  leg2->AddEntry(hNuSlices, "Neutrino Candidates", "le");
  leg2->AddEntry(hClearCosmic, "Clear Cosmics", "le");
  leg2->Draw();
  cSub->SaveAs(saveDir + "/slices_subtraction.png");

  TCanvas *cComp = new TCanvas("cComp", "cComp");
  cComp->cd();
  TH1D *hComp = sCompleteness->ToTH1(target_pot);
  hComp->SetTitle(";Completeness;Neutrino Candidate Slices (1x10^{21} POT)");
  hComp->SetLineColor(kGreen+2);
  hComp->Draw("histe");
  cComp->SaveAs(saveDir + "/slice_completeness.png");
}
