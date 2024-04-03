#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "Common.C"

void ObservablesResolution(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  const double pizeroMomBins[9] = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double showerEnBins[13] = { 0., 30., 60., 90., 120., 150., 180., 220., 260., 300., 400., 600., 1000. };
  const double cosThetaBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  TCanvas *cPiZeroMomFractionalResolution = new TCanvas("cPiZeroMomFractionalResolution", "cPiZeroMomFractionalResolution");
  cPiZeroMomFractionalResolution->cd();

  TH1F *hPiZeroMomFractionalResolution = new TH1F("hPiZeroMomFractionalResolution", ";p_{#pi^{0}} (#frac{Reco - True}{True});#pi^{0}", 25, -1, 1);
  ncpizeroEvents->Draw("(slc_best_pzc_pizero_mom-slc_true_pz_pizero_mom*1e3)/(slc_true_pz_pizero_mom*1e3)>>hPiZeroMomFractionalResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hPiZeroMomFractionalResolution->SetLineColor(kMagenta+2);
  hPiZeroMomFractionalResolution->Draw();

  cPiZeroMomFractionalResolution->SaveAs(saveDir + "/pizero_momentum_fractional_resolution.pdf");
  cPiZeroMomFractionalResolution->SaveAs(saveDir + "/pizero_momentum_fractional_resolution.png");

  TCanvas *cPiZeroMomResolution = new TCanvas("cPiZeroMomResolution", "cPiZeroMomResolution");
  cPiZeroMomResolution->cd();

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (MeV) (Reco - True);#pi^{0}", 25, -300, 200);
  ncpizeroEvents->Draw("(slc_best_pzc_pizero_mom-slc_true_pz_pizero_mom*1e3)>>hPiZeroMomResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hPiZeroMomResolution->SetLineColor(kMagenta+2);
  hPiZeroMomResolution->Draw();

  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_momentum_resolution.pdf");
  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_momentum_resolution.png");

  TCanvas *cPiZeroMom2DFractionalResolution = new TCanvas("cPiZeroMom2DFractionalResolution", "cPiZeroMom2DFractionalResolution");
  cPiZeroMom2DFractionalResolution->cd();
  cPiZeroMom2DFractionalResolution->SetRightMargin(.2);

  TH2F *hPiZeroMom2DFractionalResolution = new TH2F("hPiZeroMom2DFractionalResolution", ";True p_{#pi^{0}} (MeV);p_{#pi^{0}} (#frac{Reco - True}{True});#pi^{0}", 8, pizeroMomBins, 25, -1, 1);
  ncpizeroEvents->Draw("(slc_best_pzc_pizero_mom-slc_true_pz_pizero_mom*1e3)/(slc_true_pz_pizero_mom*1e3):slc_true_pz_pizero_mom*1e3>>hPiZeroMom2DFractionalResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hPiZeroMom2DFractionalResolution->Draw("colz");
  hPiZeroMom2DFractionalResolution->GetYaxis()->SetTitleOffset(1.25);

  cPiZeroMom2DFractionalResolution->SaveAs(saveDir + "/pizero_momentum_twod_fractional_resolution.pdf");
  cPiZeroMom2DFractionalResolution->SaveAs(saveDir + "/pizero_momentum_twod_fractional_resolution.png");

  TCanvas *cPiZeroMom2DResolution = new TCanvas("cPiZeroMom2DResolution", "cPiZeroMom2DResolution");
  cPiZeroMom2DResolution->cd();
  cPiZeroMom2DResolution->SetRightMargin(.2);

  TH2F *hPiZeroMom2DResolution = new TH2F("hPiZeroMom2DResolution", ";True p_{#pi^{0}} (MeV);p_{#pi^{0}} (MeV) (Reco - True);#pi^{0}", 8, pizeroMomBins, 25, -300, 200);
  ncpizeroEvents->Draw("(slc_best_pzc_pizero_mom-slc_true_pz_pizero_mom*1e3):slc_true_pz_pizero_mom*1e3>>hPiZeroMom2DResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hPiZeroMom2DResolution->Draw("colz");
  hPiZeroMom2DResolution->GetYaxis()->SetTitleOffset(1.25);

  cPiZeroMom2DResolution->SaveAs(saveDir + "/pizero_momentum_twod_resolution.pdf");
  cPiZeroMom2DResolution->SaveAs(saveDir + "/pizero_momentum_twod_resolution.png");

  TCanvas *cCosineThetaPiZeroFractionalResolution = new TCanvas("cCosineThetaPiZeroFractionalResolution", "cCosineThetaPiZeroFractionalResolution");
  cCosineThetaPiZeroFractionalResolution->cd();

  TH1F *hCosineThetaPiZeroFractionalResolution = new TH1F("hCosineThetaPiZeroFractionalResolution", ";cos(#theta_{#pi^{0}}) (#frac{Reco - True}{True});#pi^{0}", 40, -1, 1);
  ncpizeroEvents->Draw("(slc_best_pzc_cos_theta_pizero-slc_true_pz_cos_theta_pizero)/slc_true_pz_cos_theta_pizero>>hCosineThetaPiZeroFractionalResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hCosineThetaPiZeroFractionalResolution->SetLineColor(kMagenta+2);
  hCosineThetaPiZeroFractionalResolution->Draw();

  cCosineThetaPiZeroFractionalResolution->SaveAs(saveDir + "/cosine_theta_pizero_fractional_resolution.pdf");
  cCosineThetaPiZeroFractionalResolution->SaveAs(saveDir + "/cosine_theta_pizero_fractional_resolution.png");

  TCanvas *cCosineThetaPiZeroResolution = new TCanvas("cCosineThetaPiZeroResolution", "cCosineThetaPiZeroResolution");
  cCosineThetaPiZeroResolution->cd();

  TH1F *hCosineThetaPiZeroResolution = new TH1F("hCosineThetaPiZeroResolution", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 40, -2, 2);
  ncpizeroEvents->Draw("(slc_best_pzc_cos_theta_pizero-slc_true_pz_cos_theta_pizero)>>hCosineThetaPiZeroResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hCosineThetaPiZeroResolution->SetLineColor(kMagenta+2);
  hCosineThetaPiZeroResolution->Draw();

  cCosineThetaPiZeroResolution->SaveAs(saveDir + "/cosine_theta_pizero_resolution.pdf");
  cCosineThetaPiZeroResolution->SaveAs(saveDir + "/cosine_theta_pizero_resolution.png");

  TCanvas *cCosineThetaPiZero2DFractionalResolution = new TCanvas("cCosineThetaPiZero2DFractionalResolution", "cCosineThetaPiZero2DFractionalResolution");
  cCosineThetaPiZero2DFractionalResolution->cd();
  cCosineThetaPiZero2DFractionalResolution->SetRightMargin(.2);

  TH2F *hCosineThetaPiZero2DFractionalResolution = new TH2F("hCosineThetaPiZero2DFractionalResolution", ";True cos(#theta_{#pi^{0}});cos(#theta_{#pi^{0}}) (#frac{Reco - True}{True});#pi^{0}", 9, cosThetaBins, 40, -1, 1);
  ncpizeroEvents->Draw("(slc_best_pzc_cos_theta_pizero-slc_true_pz_cos_theta_pizero)/(slc_true_pz_cos_theta_pizero):slc_true_pz_cos_theta_pizero>>hCosineThetaPiZero2DFractionalResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hCosineThetaPiZero2DFractionalResolution->Draw("colz");
  hCosineThetaPiZero2DFractionalResolution->GetYaxis()->SetTitleOffset(1.25);

  cCosineThetaPiZero2DFractionalResolution->SaveAs(saveDir + "/cosine_theta_pizero_twod_fractional_resolution.pdf");
  cCosineThetaPiZero2DFractionalResolution->SaveAs(saveDir + "/cosine_theta_pizero_twod_fractional_resolution.png");

  TCanvas *cCosineThetaPiZero2DResolution = new TCanvas("cCosineThetaPiZero2DResolution", "cCosineThetaPiZero2DResolution");
  cCosineThetaPiZero2DResolution->cd();
  cCosineThetaPiZero2DResolution->SetRightMargin(.2);

  TH2F *hCosineThetaPiZero2DResolution = new TH2F("hCosineThetaPiZero2DResolution", ";True cos(#theta_{#pi^{0}});cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 9, cosThetaBins, 40, -2, 2);
  ncpizeroEvents->Draw("(slc_best_pzc_cos_theta_pizero-slc_true_pz_cos_theta_pizero):slc_true_pz_cos_theta_pizero>>hCosineThetaPiZero2DResolution",
                       "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hCosineThetaPiZero2DResolution->Draw("colz");
  hCosineThetaPiZero2DResolution->GetYaxis()->SetTitleOffset(1.25);

  cCosineThetaPiZero2DResolution->SaveAs(saveDir + "/cosine_theta_pizero_twod_resolution.pdf");
  cCosineThetaPiZero2DResolution->SaveAs(saveDir + "/cosine_theta_pizero_twod_resolution.png");

  TCanvas *cShowerEnergy = new TCanvas("cShowerEnergy", "cShowerEnergy");
  cShowerEnergy->cd();

  TH1F *hShowerEnergy = new TH1F("hShowerEnergy", ";True E (MeV));#gamma/ 30 MeV", 12, showerEnBins);
  ncpizeroEvents->Draw("1e3*slc_pfp_true_energy>>hShowerEnergy",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  NormaliseEntriesByBinWidth(hShowerEnergy, 30);
  hShowerEnergy->SetLineColor(kMagenta+2);
  hShowerEnergy->GetYaxis()->SetTitleOffset(1.25);
  hShowerEnergy->Draw();

  cShowerEnergy->SaveAs(saveDir + "/shower_energy.pdf");
  cShowerEnergy->SaveAs(saveDir + "/shower_energy.png");

  TCanvas *cShowerEnergyReco = new TCanvas("cShowerEnergyReco", "cShowerEnergyReco");
  cShowerEnergyReco->cd();

  TH1F *hShowerEnergyReco = new TH1F("hShowerEnergyReco", ";Reco E (MeV));#gamma / 30 MeV", 12, showerEnBins);
  ncpizeroEvents->Draw("slc_pfp_shower_energy>>hShowerEnergyReco",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  NormaliseEntriesByBinWidth(hShowerEnergyReco, 30);
  hShowerEnergyReco->SetLineColor(kMagenta+2);
  hShowerEnergyReco->GetYaxis()->SetTitleOffset(1.25);
  hShowerEnergyReco->Draw();

  cShowerEnergyReco->SaveAs(saveDir + "/shower_energy_reco.pdf");
  cShowerEnergyReco->SaveAs(saveDir + "/shower_energy_reco.png");

  TCanvas *cShowerEnergyFractionalResolution = new TCanvas("cShowerEnergyFractionalResolution", "cShowerEnergyFractionalResolution");
  cShowerEnergyFractionalResolution->cd();

  TH1F *hShowerEnergyFractionalResolution = new TH1F("hShowerEnergyFractionalResolution", ";E (#frac{Reco - True}{True});#gamma", 40, -1, 1);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-1e3*slc_pfp_true_energy)/(1e3*slc_pfp_true_energy)>>hShowerEnergyFractionalResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergyFractionalResolution->SetLineColor(kMagenta+2);
  hShowerEnergyFractionalResolution->GetYaxis()->SetTitleOffset(1.25);
  hShowerEnergyFractionalResolution->Draw();

  cShowerEnergyFractionalResolution->SaveAs(saveDir + "/shower_energy_fractional_resolution.pdf");
  cShowerEnergyFractionalResolution->SaveAs(saveDir + "/shower_energy_fractional_resolution.png");

  TCanvas *cShowerEnergyResolution = new TCanvas("cShowerEnergyResolution", "cShowerEnergyResolution");
  cShowerEnergyResolution->cd();

  TH1F *hShowerEnergyResolution = new TH1F("hShowerEnergyResolution", ";E (MeV) (Reco - True);#gamma", 50, -300, 200);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy*1e3)>>hShowerEnergyResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergyResolution->SetLineColor(kMagenta+2);
  hShowerEnergyResolution->Draw();

  cShowerEnergyResolution->SaveAs(saveDir + "/shower_energy_resolution.pdf");
  cShowerEnergyResolution->SaveAs(saveDir + "/shower_energy_resolution.png");

  TCanvas *cShowerEnergy2DFractionalResolution = new TCanvas("cShowerEnergy2DFractionalResolution", "cShowerEnergy2DFractionalResolution");
  cShowerEnergy2DFractionalResolution->cd();
  cShowerEnergy2DFractionalResolution->SetRightMargin(.2);

  TH2F *hShowerEnergy2DFractionalResolution = new TH2F("hShowerEnergy2DFractionalResolution", ";True E (MeV);E (#frac{Reco - True}{True});#gamma", 12, showerEnBins, 40, -1, 1);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy*1e3)/(slc_pfp_true_energy*1e3):slc_pfp_true_energy*1e3>>hShowerEnergy2DFractionalResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergy2DFractionalResolution->Draw("colz");
  hShowerEnergy2DFractionalResolution->GetYaxis()->SetTitleOffset(1.25);

  cShowerEnergy2DFractionalResolution->SaveAs(saveDir + "/shower_energy_twod_fractional_resolution.pdf");
  cShowerEnergy2DFractionalResolution->SaveAs(saveDir + "/shower_energy_twod_fractional_resolution.png");

  TCanvas *cShowerEnergy2DFractionalResolutionProfile = new TCanvas("cShowerEnergy2DFractionalResolutionProfile", "cShowerEnergy2DFractionalResolutionProfile");
  cShowerEnergy2DFractionalResolutionProfile->cd();
  cShowerEnergy2DFractionalResolutionProfile->SetRightMargin(.2);

  TH1F *hShowerEnergy2DFractionalResolutionProfile = (TH1F*) hShowerEnergy2DFractionalResolution->ProfileX();

  hShowerEnergy2DFractionalResolutionProfile->SetLineColor(kMagenta+2);
  hShowerEnergy2DFractionalResolutionProfile->GetYaxis()->SetTitle("E (#frac{Reco - True}{True})");
  hShowerEnergy2DFractionalResolutionProfile->Draw("histe][");

  cShowerEnergy2DFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_fractional_resolution_profile.pdf");
  cShowerEnergy2DFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_fractional_resolution_profile.png");

  TCanvas *cShowerEnergy2DRecoFractionalResolution = new TCanvas("cShowerEnergy2DRecoFractionalResolution", "cShowerEnergy2DRecoFractionalResolution");
  cShowerEnergy2DRecoFractionalResolution->cd();
  cShowerEnergy2DRecoFractionalResolution->SetRightMargin(.2);

  TH2F *hShowerEnergy2DRecoFractionalResolution = new TH2F("hShowerEnergy2DRecoFractionalResolution", ";Reco E (MeV);E (#frac{Reco - True}{True});#gamma", 12, showerEnBins, 40, -1, 1);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy*1e3)/(slc_pfp_true_energy*1e3):slc_pfp_shower_energy>>hShowerEnergy2DRecoFractionalResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergy2DRecoFractionalResolution->Draw("colz");
  hShowerEnergy2DRecoFractionalResolution->GetYaxis()->SetTitleOffset(1.25);

  cShowerEnergy2DRecoFractionalResolution->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution.pdf");
  cShowerEnergy2DRecoFractionalResolution->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution.png");

  TCanvas *cShowerEnergy2DRecoFractionalResolutionProfile = new TCanvas("cShowerEnergy2DRecoFractionalResolutionProfile", "cShowerEnergy2DRecoFractionalResolutionProfile");
  cShowerEnergy2DRecoFractionalResolutionProfile->cd();
  cShowerEnergy2DRecoFractionalResolutionProfile->SetRightMargin(.2);

  const TString outFileName = saveDir + "/shower_energy_correction_hist.root";
  TFile *outFile = new TFile(outFileName, "recreate");

  TH1F *hShowerEnergy2DRecoFractionalResolutionProfile = (TH1F*) hShowerEnergy2DRecoFractionalResolution->ProfileX();

  hShowerEnergy2DRecoFractionalResolutionProfile->SetLineColor(kMagenta+2);
  hShowerEnergy2DRecoFractionalResolutionProfile->GetYaxis()->SetTitle("E (#frac{Reco - True}{True})");
  hShowerEnergy2DRecoFractionalResolutionProfile->Draw("histe][");

  cShowerEnergy2DRecoFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_profile.pdf");
  cShowerEnergy2DRecoFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_profile.png");

  hShowerEnergy2DRecoFractionalResolutionProfile->Write();

  TCanvas *cShowerEnergy2DResolution = new TCanvas("cShowerEnergy2DResolution", "cShowerEnergy2DResolution");
  cShowerEnergy2DResolution->cd();
  cShowerEnergy2DResolution->SetRightMargin(.2);

  TH2F *hShowerEnergy2DResolution = new TH2F("hShowerEnergy2DResolution", ";True E (MeV);E (MeV) (Reco - True);#gamma", 12, showerEnBins, 50, -300, 200);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy*1e3):slc_pfp_true_energy*1e3>>hShowerEnergy2DResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergy2DResolution->Draw("colz");
  hShowerEnergy2DResolution->GetYaxis()->SetTitleOffset(1.25);

  cShowerEnergy2DResolution->SaveAs(saveDir + "/shower_energy_twod_resolution.pdf");
  cShowerEnergy2DResolution->SaveAs(saveDir + "/shower_energy_twod_resolution.png");

  TCanvas *cPhotonShowerTrackDirectionResolution = new TCanvas("cPhotonShowerTrackDirectionResolution", "cPhotonShowerTrackDirectionResolution");
  cPhotonShowerTrackDirectionResolution->cd();

  TString dotProductShw = "(slc_pfp_shower_dir_x * slc_pfp_true_p_x + slc_pfp_shower_dir_y * slc_pfp_true_p_y + slc_pfp_shower_dir_z * slc_pfp_true_p_z)";
  TString dotProductTrk = "(slc_pfp_track_dir_x * slc_pfp_true_p_x + slc_pfp_track_dir_y * slc_pfp_true_p_y + slc_pfp_track_dir_z * slc_pfp_true_p_z)";
  TString magRecoShw    = "sqrt(slc_pfp_shower_dir_x * slc_pfp_shower_dir_x + slc_pfp_shower_dir_y * slc_pfp_shower_dir_y + slc_pfp_shower_dir_z * slc_pfp_shower_dir_z)";
  TString magRecoTrk    = "sqrt(slc_pfp_track_dir_x * slc_pfp_track_dir_x + slc_pfp_track_dir_y * slc_pfp_track_dir_y + slc_pfp_track_dir_z * slc_pfp_track_dir_z)";
  TString magTrue       = "sqrt(slc_pfp_true_p_x * slc_pfp_true_p_x + slc_pfp_true_p_y * slc_pfp_true_p_y + slc_pfp_true_p_z * slc_pfp_true_p_z)";
  TString angleShw      = "acos(" + dotProductShw + "/(" + magRecoShw + "*" + magTrue + ")) * TMath::RadToDeg()";
  TString angleTrk      = "acos(" + dotProductTrk + "/(" + magRecoTrk + "*" + magTrue + ")) * TMath::RadToDeg()";

  TH1F *hPhotonShowerDirectionResolution = new TH1F("hPhotonShowerDirectionResolution", ";#theta_{True Reco} (#circ);#gamma", 36, 0, 180);
  ncpizeroEvents->Draw(angleShw + ">>hPhotonShowerDirectionResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  TH1F *hPhotonTrackDirectionResolution = new TH1F("hPhotonTrackDirectionResolution", ";#theta_{True Reco} (#circ);#gamma", 36, 0, 180);
  ncpizeroEvents->Draw(angleTrk + ">>hPhotonTrackDirectionResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hPhotonShowerDirectionResolution->SetLineColor(kMagenta+2);
  hPhotonTrackDirectionResolution->SetLineColor(kCyan+2);
  hPhotonTrackDirectionResolution->Draw();
  hPhotonShowerDirectionResolution->Draw("same");

  TLegend *lPhotonShowerDirectionResolution = new TLegend(.55, .4, .85, .55);
  lPhotonShowerDirectionResolution->AddEntry(hPhotonShowerDirectionResolution, "Shower Characterisation", "l");
  lPhotonShowerDirectionResolution->AddEntry(hPhotonTrackDirectionResolution, "Track Characterisation", "l");
  lPhotonShowerDirectionResolution->Draw();

  cPhotonShowerTrackDirectionResolution->SaveAs(saveDir + "/photon_shower_track_direction_resolution.pdf");
  cPhotonShowerTrackDirectionResolution->SaveAs(saveDir + "/photon_shower_track_direction_resolution.png");

  TCanvas *cPhotonShowerTrackDirectionResolutionCumulative = new TCanvas("cPhotonShowerTrackDirectionResolutionCumulative", "cPhotonShowerTrackDirectionResolutionCumulative");
  cPhotonShowerTrackDirectionResolutionCumulative->cd();

  TH1 *hPhotonTrackDirectionResolutionCumulative = hPhotonTrackDirectionResolution->GetCumulative();
  hPhotonTrackDirectionResolutionCumulative->Scale(100./hPhotonTrackDirectionResolution->GetEntries());
  hPhotonTrackDirectionResolutionCumulative->SetMinimum(0);
  hPhotonTrackDirectionResolutionCumulative->GetYaxis()->SetTitle("Cumulative Proportion (%)");
  hPhotonTrackDirectionResolutionCumulative->SetLineColor(kCyan+2);
  hPhotonTrackDirectionResolutionCumulative->Draw();

  TH1 *hPhotonShowerDirectionResolutionCumulative = hPhotonShowerDirectionResolution->GetCumulative();
  hPhotonShowerDirectionResolutionCumulative->Scale(100./hPhotonShowerDirectionResolution->GetEntries());
  hPhotonShowerDirectionResolutionCumulative->SetMinimum(0);
  hPhotonShowerDirectionResolutionCumulative->SetLineColor(kMagenta+2);
  hPhotonShowerDirectionResolutionCumulative->Draw("same");

  lPhotonShowerDirectionResolution->Draw();

  cPhotonShowerTrackDirectionResolutionCumulative->SaveAs(saveDir + "/photon_shower_track_direction_resolution_cumulative.pdf");
  cPhotonShowerTrackDirectionResolutionCumulative->SaveAs(saveDir + "/photon_shower_track_direction_resolution_cumulative.png");

  TCanvas *cPhotonShowerDirectionResolution = new TCanvas("cPhotonShowerDirectionResolution", "cPhotonShowerDirectionResolution");
  cPhotonShowerDirectionResolution->cd();
  hPhotonShowerDirectionResolution->Draw();
  cPhotonShowerDirectionResolution->SaveAs(saveDir + "/photon_shower_direction_resolution.pdf");
  cPhotonShowerDirectionResolution->SaveAs(saveDir + "/photon_shower_direction_resolution.png");

  TCanvas *cPhotonShowerDirectionResolutionCumulative = new TCanvas("cPhotonShowerDirectionResolutionCumulative", "cPhotonShowerDirectionResolutionCumulative");
  cPhotonShowerDirectionResolutionCumulative->cd();
  hPhotonShowerDirectionResolutionCumulative->Draw();
  cPhotonShowerDirectionResolutionCumulative->SaveAs(saveDir + "/photon_shower_direction_resolution_cumulative.pdf");
  cPhotonShowerDirectionResolutionCumulative->SaveAs(saveDir + "/photon_shower_direction_resolution_cumulative.png");

  TCanvas *cPhotonShowerDirection2DResolution = new TCanvas("cPhotonShowerDirection2DResolution", "cPhotonShowerDirection2DResolution");
  cPhotonShowerDirection2DResolution->cd();
  cPhotonShowerDirection2DResolution->SetRightMargin(.2);

  TH2F *hPhotonShowerDirection2DResolution = new TH2F("hPhotonShowerDirection2DResolution", ";True E (MeV);#theta_{True Reco} (#circ);#gamma", 12, showerEnBins, 36, 0, 180);
  ncpizeroEvents->Draw(angleShw + ":slc_pfp_true_energy*1e3>>hPhotonShowerDirection2DResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hPhotonShowerDirection2DResolution->Draw("colz");
  hPhotonShowerDirection2DResolution->GetYaxis()->SetTitleOffset(1.25);

  cPhotonShowerDirection2DResolution->SaveAs(saveDir + "/photon_shower_direction_twod_resolution.pdf");
  cPhotonShowerDirection2DResolution->SaveAs(saveDir + "/photon_shower_direction_twod_resolution.png");

  TCanvas *cPhotonTrackDirection2DResolution = new TCanvas("cPhotonTrackDirection2DResolution", "cPhotonTrackDirection2DResolution");
  cPhotonTrackDirection2DResolution->cd();
  cPhotonTrackDirection2DResolution->SetRightMargin(.2);

  TH2F *hPhotonTrackDirection2DResolution = new TH2F("hPhotonTrackDirection2DResolution", ";True E (MeV);#theta_{True Reco} (#circ);#gamma", 12, showerEnBins, 36, 0, 180);
  ncpizeroEvents->Draw(angleTrk + ":slc_pfp_true_energy*1e3>>hPhotonTrackDirection2DResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hPhotonTrackDirection2DResolution->Draw("colz");
  hPhotonTrackDirection2DResolution->GetYaxis()->SetTitleOffset(1.25);

  cPhotonTrackDirection2DResolution->SaveAs(saveDir + "/photon_track_direction_twod_resolution.pdf");
  cPhotonTrackDirection2DResolution->SaveAs(saveDir + "/photon_track_direction_twod_resolution.png");
}
