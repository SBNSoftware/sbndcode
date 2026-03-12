#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Common.C"

void ObservablesResolutionTmp(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gStyle->SetTitleOffset(1.35, "y");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(rockboxFile);
  events->Add(ncpizeroFile);

  const double pizeroMomBins[9] = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double showerEnBins[13] = { 0., 30., 60., 90., 120., 150., 180., 220., 260., 300., 400., 600., 1000. };
  const double cosThetaBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  TCanvas *cShowerEnergy2DRecoFractionalResolutionSimple = new TCanvas("cShowerEnergy2DRecoFractionalResolutionSimple", "cShowerEnergy2DRecoFractionalResolutionSimple");
  cShowerEnergy2DRecoFractionalResolutionSimple->cd();
  cShowerEnergy2DRecoFractionalResolutionSimple->SetRightMargin(.2);

  TH2F *hShowerEnergy2DRecoFractionalResolutionSimple = new TH2F("hShowerEnergy2DRecoFractionalResolutionSimple", ";Shower E (MeV);E (#frac{True}{Reco});#gamma", 12, showerEnBins, 40, 0, 2);
  events->Draw("slc_pfp_true_energy*1e3/slc_pfp_shower_energy:slc_pfp_shower_energy>>hShowerEnergy2DRecoFractionalResolutionSimple",
               "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergy2DRecoFractionalResolutionSimple->Draw("colz");
  hShowerEnergy2DRecoFractionalResolutionSimple->GetYaxis()->SetTitleOffset(1.25);

  cShowerEnergy2DRecoFractionalResolutionSimple->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_simple.pdf");
  cShowerEnergy2DRecoFractionalResolutionSimple->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_simple.png");

  TCanvas *cShowerEnergy2DRecoFractionalResolutionSimpleProfile = new TCanvas("cShowerEnergy2DRecoFractionalResolutionSimpleProfile", "cShowerEnergy2DRecoFractionalResolutionSimpleProfile");
  cShowerEnergy2DRecoFractionalResolutionSimpleProfile->cd();
  cShowerEnergy2DRecoFractionalResolutionSimpleProfile->SetRightMargin(.2);

  TH1F *hShowerEnergy2DRecoFractionalResolutionSimpleProfile = (TH1F*) hShowerEnergy2DRecoFractionalResolutionSimple->ProfileX();

  hShowerEnergy2DRecoFractionalResolutionSimpleProfile->SetLineColor(kMagenta+2);
  hShowerEnergy2DRecoFractionalResolutionSimpleProfile->GetYaxis()->SetTitle("E (#frac{True}{Reco})");
  hShowerEnergy2DRecoFractionalResolutionSimpleProfile->Draw("histe][");

  cShowerEnergy2DRecoFractionalResolutionSimpleProfile->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_simple.profile.pdf");
  cShowerEnergy2DRecoFractionalResolutionSimpleProfile->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_simple.profile.png");
}
