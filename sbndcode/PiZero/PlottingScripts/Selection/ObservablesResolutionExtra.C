#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Common.C"

void ObservablesResolutionExtra(const TString productionVersion)
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

  TString dotProductShw = "(slc_pfp_shower_dir_x * slc_pfp_true_p_x + slc_pfp_shower_dir_y * slc_pfp_true_p_y + slc_pfp_shower_dir_z * slc_pfp_true_p_z)";
  TString dotProductTrk = "(slc_pfp_track_dir_x * slc_pfp_true_p_x + slc_pfp_track_dir_y * slc_pfp_true_p_y + slc_pfp_track_dir_z * slc_pfp_true_p_z)";
  TString magRecoShw    = "sqrt(slc_pfp_shower_dir_x * slc_pfp_shower_dir_x + slc_pfp_shower_dir_y * slc_pfp_shower_dir_y + slc_pfp_shower_dir_z * slc_pfp_shower_dir_z)";
  TString magRecoTrk    = "sqrt(slc_pfp_track_dir_x * slc_pfp_track_dir_x + slc_pfp_track_dir_y * slc_pfp_track_dir_y + slc_pfp_track_dir_z * slc_pfp_track_dir_z)";
  TString magTrue       = "sqrt(slc_pfp_true_p_x * slc_pfp_true_p_x + slc_pfp_true_p_y * slc_pfp_true_p_y + slc_pfp_true_p_z * slc_pfp_true_p_z)";
  TString angleShw      = "acos(" + dotProductShw + "/(" + magRecoShw + "*" + magTrue + ")) * TMath::RadToDeg()";
  TString angleTrk      = "acos(" + dotProductTrk + "/(" + magRecoTrk + "*" + magTrue + ")) * TMath::RadToDeg()";

  TCanvas *cPhotonShowerDirection2DResolutionProfile = new TCanvas("cPhotonShowerDirection2DResolutionProfile", "cPhotonShowerDirection2DResolutionProfile");
  cPhotonShowerDirection2DResolutionProfile->cd();

  TH2F *hPhotonShowerDirection2DResolution = new TH2F("hPhotonShowerDirection2DResolution", ";True E (MeV);#theta_{True Reco} (#circ);#gamma", 12, showerEnBins, 36, 0, 180);
  events->Draw(angleShw + ":slc_pfp_true_energy*1e3>>hPhotonShowerDirection2DResolution",
               "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  TH2F *hPhotonTrackDirection2DResolution = new TH2F("hPhotonTrackDirection2DResolution", ";True E (MeV);#theta_{True Reco} (#circ);#gamma", 12, showerEnBins, 36, 0, 180);
  events->Draw(angleTrk + ":slc_pfp_true_energy*1e3>>hPhotonTrackDirection2DResolution",
               "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  TH1F *hPhotonShowerDirection2DResolutionProfile = (TH1F*) hPhotonShowerDirection2DResolution->ProfileX();
  hPhotonShowerDirection2DResolutionProfile->SetMarkerStyle(0);
  hPhotonShowerDirection2DResolutionProfile->SetLineColor(kMagenta+2);
  hPhotonShowerDirection2DResolutionProfile->GetYaxis()->SetTitle("#theta_{True Reco} (#circ)");
  hPhotonShowerDirection2DResolutionProfile->SetMinimum(0);
  hPhotonShowerDirection2DResolutionProfile->Draw("hist][e");

  TH1F *hPhotonTrackDirection2DResolutionProfile = (TH1F*) hPhotonTrackDirection2DResolution->ProfileX();
  hPhotonTrackDirection2DResolutionProfile->SetMarkerStyle(0);
  hPhotonTrackDirection2DResolutionProfile->SetLineColor(kCyan+2);
  hPhotonTrackDirection2DResolutionProfile->GetYaxis()->SetTitle("#theta_{True Reco} (#circ)");
  hPhotonTrackDirection2DResolutionProfile->Draw("hist][esame");

  TLegend *lPhotonShowerDirection2DResolutionProfile = new TLegend(.55, .5, .85, .65);
  lPhotonShowerDirection2DResolutionProfile->AddEntry(hPhotonShowerDirection2DResolutionProfile, "Shower Characterisation", "l");
  lPhotonShowerDirection2DResolutionProfile->AddEntry(hPhotonTrackDirection2DResolutionProfile, "Track Characterisation", "l");
  lPhotonShowerDirection2DResolutionProfile->Draw();

  cPhotonShowerDirection2DResolutionProfile->SaveAs(saveDir + "/photon_shower_direction_twod_resolution_profile.pdf");
  cPhotonShowerDirection2DResolutionProfile->SaveAs(saveDir + "/photon_shower_direction_twod_resolution_profile.png");

  TCanvas *cPhotonShowerDirection2DResolutionProfileReco = new TCanvas("cPhotonShowerDirection2DResolutionProfileReco", "cPhotonShowerDirection2DResolutionProfileReco");
  cPhotonShowerDirection2DResolutionProfileReco->cd();

  TH2F *hPhotonShowerDirection2DResolutionReco = new TH2F("hPhotonShowerDirection2DResolutionReco", ";Shower E (MeV);#theta_{True Reco} (#circ);#gamma", 12, showerEnBins, 36, 0, 180);
  events->Draw(angleShw + ":slc_pfp_shower_energy>>hPhotonShowerDirection2DResolutionReco",
               "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  TH2F *hPhotonTrackDirection2DResolutionReco = new TH2F("hPhotonTrackDirection2DResolutionReco", ";Shower E (MeV);#theta_{True Reco} (#circ);#gamma", 12, showerEnBins, 36, 0, 180);
  events->Draw(angleTrk + ":slc_pfp_shower_energy>>hPhotonTrackDirection2DResolutionReco",
               "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  TH1F *hPhotonShowerDirection2DResolutionProfileReco = (TH1F*) hPhotonShowerDirection2DResolutionReco->ProfileX();
  hPhotonShowerDirection2DResolutionProfileReco->SetMarkerStyle(0);
  hPhotonShowerDirection2DResolutionProfileReco->SetLineColor(kMagenta+2);
  hPhotonShowerDirection2DResolutionProfileReco->GetYaxis()->SetTitle("#theta_{True Reco} (#circ)");
  hPhotonShowerDirection2DResolutionProfileReco->SetMinimum(0);
  hPhotonShowerDirection2DResolutionProfileReco->Draw("hist][e");

  TH1F *hPhotonTrackDirection2DResolutionProfileReco = (TH1F*) hPhotonTrackDirection2DResolutionReco->ProfileX();
  hPhotonTrackDirection2DResolutionProfileReco->SetMarkerStyle(0);
  hPhotonTrackDirection2DResolutionProfileReco->SetLineColor(kCyan+2);
  hPhotonTrackDirection2DResolutionProfileReco->GetYaxis()->SetTitle("#theta_{True Reco} (#circ)");
  hPhotonTrackDirection2DResolutionProfileReco->Draw("hist][esame");

  lPhotonShowerDirection2DResolutionProfile->Draw();

  cPhotonShowerDirection2DResolutionProfileReco->Update();

  TLine *cutLine = new TLine(150, gPad->GetUymin(), 150, gPad->GetUymax());
  cutLine->SetLineWidth(3);
  cutLine->SetLineColor(kRed+2);
  cutLine->SetLineStyle(9);
  cutLine->Draw();

  cPhotonShowerDirection2DResolutionProfileReco->SaveAs(saveDir + "/photon_shower_direction_twod_resolution_profile_reco.pdf");
  cPhotonShowerDirection2DResolutionProfileReco->SaveAs(saveDir + "/photon_shower_direction_twod_resolution_profile_reco.png");
}
