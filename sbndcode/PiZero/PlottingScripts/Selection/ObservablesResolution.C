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

  TH2F *hShowerEnergy2DFractionalResolution = new TH2F("hShowerEnergy2DFractionalResolution", ";True E (MeV);E (#frac{Reco - True}{True});#gamma", 8, pizeroMomBins, 40, -1, 1);
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
  hShowerEnergy2DFractionalResolutionProfile->Draw("histe][");

  cShowerEnergy2DFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_fractional_resolution_profile.pdf");
  cShowerEnergy2DFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_fractional_resolution_profile.png");

  TCanvas *cShowerEnergy2DRecoFractionalResolution = new TCanvas("cShowerEnergy2DRecoFractionalResolution", "cShowerEnergy2DRecoFractionalResolution");
  cShowerEnergy2DRecoFractionalResolution->cd();
  cShowerEnergy2DRecoFractionalResolution->SetRightMargin(.2);

  TH2F *hShowerEnergy2DRecoFractionalResolution = new TH2F("hShowerEnergy2DRecoFractionalResolution", ";Reco E (MeV);E (#frac{Reco - True}{True});#gamma", 8, pizeroMomBins, 40, -1, 1);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy*1e3)/(slc_pfp_true_energy*1e3):slc_pfp_shower_energy>>hShowerEnergy2DRecoFractionalResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergy2DRecoFractionalResolution->Draw("colz");
  hShowerEnergy2DRecoFractionalResolution->GetYaxis()->SetTitleOffset(1.25);

  cShowerEnergy2DRecoFractionalResolution->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution.pdf");
  cShowerEnergy2DRecoFractionalResolution->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution.png");

  TCanvas *cShowerEnergy2DRecoFractionalResolutionProfile = new TCanvas("cShowerEnergy2DRecoFractionalResolutionProfile", "cShowerEnergy2DRecoFractionalResolutionProfile");
  cShowerEnergy2DRecoFractionalResolutionProfile->cd();
  cShowerEnergy2DRecoFractionalResolutionProfile->SetRightMargin(.2);

  TH1F *hShowerEnergy2DRecoFractionalResolutionProfile = (TH1F*) hShowerEnergy2DRecoFractionalResolution->ProfileX();

  hShowerEnergy2DRecoFractionalResolutionProfile->SetLineColor(kMagenta+2);
  hShowerEnergy2DRecoFractionalResolutionProfile->Draw("histe][");

  cShowerEnergy2DRecoFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_profile.pdf");
  cShowerEnergy2DRecoFractionalResolutionProfile->SaveAs(saveDir + "/shower_energy_twod_reco_fractional_resolution_profile.png");

  TCanvas *cShowerEnergy2DResolution = new TCanvas("cShowerEnergy2DResolution", "cShowerEnergy2DResolution");
  cShowerEnergy2DResolution->cd();
  cShowerEnergy2DResolution->SetRightMargin(.2);

  TH2F *hShowerEnergy2DResolution = new TH2F("hShowerEnergy2DResolution", ";True E (MeV);E (MeV) (Reco - True);#gamma", 8, pizeroMomBins, 50, -300, 200);
  ncpizeroEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy*1e3):slc_pfp_true_energy*1e3>>hShowerEnergy2DResolution",
                       "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5 && slc_pfp_shower_contained && slc_pfp_shower_energy>=0 && (slc_pfp_true_pdg==22)");

  hShowerEnergy2DResolution->Draw("colz");
  hShowerEnergy2DResolution->GetYaxis()->SetTitleOffset(1.25);

  cShowerEnergy2DResolution->SaveAs(saveDir + "/shower_energy_twod_resolution.pdf");
  cShowerEnergy2DResolution->SaveAs(saveDir + "/shower_energy_twod_resolution.png");
}
