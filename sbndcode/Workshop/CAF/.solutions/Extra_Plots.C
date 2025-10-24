#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#include "ExtraVars.h"
#include "ExtraCuts.h"
#include "HenryCuts.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"

TString POTString(const double &pot, const bool brackets = true)
{
  TString potString = Form("%g POT", pot);

  if(brackets)
    {
      potString.Prepend(" (");
      potString.Append(")");
    }

  potString.ReplaceAll("e+", "x10^{");
  potString.ReplaceAll(" POT", "} POT");

  return potString;
}

void SetupStyle()
{
  gStyle->SetLineWidth(4);
  gStyle->SetHistLineWidth(3);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(.5);
  gStyle->SetPadLeftMargin(.12);
  gStyle->SetTitleOffset(1, "y");

  gROOT->ForceStyle();
}

void Extra_Plots()
{
  SetupStyle();

  // This is the input dataset we're going to use
  const std::string inputName = "/scratch/LAR25/caf/reco2_neutrino.flat.caf.root";

  // The directory we want to save plots in
  const TString saveDir = "./plots";

  // The amount of POT to scale the plots to
  const double pot = 1e21;

  // This object takes the dataset and fills all the 'spectra' you request from it
  SpectrumLoader loader(inputName);

  // Binning schemes we want to use for our plots
  Binning neutrinoEnergyBins = Binning::Custom({0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 3, 5});
  Binning nSliceBins         = Binning::Simple(25, -0.5, 24.5);
  Binning completenessBins   = Binning::Simple(22, -0.05, 1.05);

  // The 'spectra' we want to create from our dataset
  Spectrum *sNeutrinoEnergy = new Spectrum("sNeutrinoEnergy", neutrinoEnergyBins, loader,
                                           kNeutrinoEnergy, kNoSpillCut);

  Spectrum *sNSlices = new Spectrum("sNSlices", nSliceBins, loader, kNSlices, kNoSpillCut);

  Spectrum *sNNeutrinoCandidateSlices = new Spectrum("sNNeutrinoCandidateSlices", nSliceBins, loader,
                                                     kNNeutrinoCandidateSlices, kNoSpillCut);

  Spectrum *sNeutrinoCandidateCompleteness = new Spectrum("sNeutrinoCandidateCompleteness", completenessBins,
                                                          loader, kSliceCompleteness, kSingleNeutrinoCandidateSlice,
                                                          kIsNeutrinoCandidateSlice);

  // Tell the loader we've declared all our spectra, time to fill them!
  loader.Go();

  // Neutrino energy plot
  TCanvas *cNeutrinoEnergy = new TCanvas("cNeutrinoEnergy", "cNeutrinoEnergy");
  cNeutrinoEnergy->cd();

  TH1D *hNeutrinoEnergy = sNeutrinoEnergy->ToTH1(pot);
  hNeutrinoEnergy->SetTitle(";E_{#nu} (GeV);Neutrinos" + POTString(pot));
  hNeutrinoEnergy->SetLineColor(kCyan+2);
  hNeutrinoEnergy->Draw("histe");

  cNeutrinoEnergy->SaveAs(saveDir + "/neutrino_energy.pdf");
  cNeutrinoEnergy->SaveAs(saveDir + "/neutrino_energy.png");

  // Number of slices
  TCanvas *cNSlices = new TCanvas("cNSlices", "cNSlices");
  cNSlices->cd();

  TH1D *hNSlices = sNSlices->ToTH1(pot);
  hNSlices->SetTitle(";N Slices;Events" + POTString(pot));
  hNSlices->SetLineColor(kSpring+2);
  hNSlices->Draw("histe");

  cNSlices->SaveAs(saveDir + "/n_slices.pdf");
  cNSlices->SaveAs(saveDir + "/n_slices.png");

  // Number of slices by type
  TCanvas *cNSlicesByType = new TCanvas("cNSlicesByType", "cNSlicesByType");
  cNSlicesByType->cd();

  TH1D *hNNeutrinoCandidateSlices = sNNeutrinoCandidateSlices->ToTH1(pot);
  hNNeutrinoCandidateSlices->SetTitle(";N Slices;Events" + POTString(pot));
  hNNeutrinoCandidateSlices->SetLineColor(kBlue+2);
  hNNeutrinoCandidateSlices->Draw("histe");

  Spectrum sNClearCosmicSlices = *sNSlices - *sNNeutrinoCandidateSlices;
  TH1D *hNClearCosmicSlices = sNClearCosmicSlices.ToTH1(pot);
  hNClearCosmicSlices->SetLineColor(kRed+2);
  hNClearCosmicSlices->Draw("histesame");

  TLegend *lNSlicesByType = new TLegend(.6, .65, .88, .8);
  lNSlicesByType->AddEntry(hNNeutrinoCandidateSlices, "Neutrino Candidates", "le");
  lNSlicesByType->AddEntry(hNClearCosmicSlices, "Clear Cosmics", "le");
  lNSlicesByType->Draw();

  cNSlicesByType->SaveAs(saveDir + "/n_slices_by_type.pdf");
  cNSlicesByType->SaveAs(saveDir + "/n_slices_by_type.png");

  // Neutrino candidate slice completeness plot
  TCanvas *cNeutrinoCandidateCompleteness = new TCanvas("cNeutrinoCandidateCompleteness", "cNeutrinoCandidateCompleteness");
  cNeutrinoCandidateCompleteness->cd();

  TH1D *hNeutrinoCandidateCompleteness = sNeutrinoCandidateCompleteness->ToTH1(pot);
  hNeutrinoCandidateCompleteness->SetTitle(";Completeness;Neutrino Candidate Slices" + POTString(pot));
  hNeutrinoCandidateCompleteness->SetLineColor(kGreen+2);
  hNeutrinoCandidateCompleteness->Draw("histe");

  cNeutrinoCandidateCompleteness->SaveAs(saveDir + "/neutrino_candidate_slices.pdf");
  cNeutrinoCandidateCompleteness->SaveAs(saveDir + "/neutrino_candidate_slices.png");
}
