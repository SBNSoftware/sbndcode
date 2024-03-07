#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

std::vector<int> *slc_true_event_type_incl = 0;
std::vector<float> *slc_comp = 0;
std::vector<bool> *slc_sel_incl = 0;
std::vector<double> *slc_best_pzc_pizero_mom = 0, *slc_best_pzc_invariant_mass = 0, *slc_best_pzc_cos_theta_pizero = 0;

std::vector<size_t> *slc_best_pzc_photon_0_id = 0, *slc_best_pzc_photon_1_id = 0;

std::vector<std::vector<double>> *slc_true_pz_pizero_mom = 0, *slc_true_pz_cos_theta_pizero = 0, *slc_pfp_shower_dir_x = 0,
  *slc_pfp_shower_dir_y = 0, *slc_pfp_shower_dir_z = 0, *slc_pfp_track_dir_x = 0, *slc_pfp_track_dir_y = 0, *slc_pfp_track_dir_z = 0,
  *slc_pfp_shower_energy = 0;

TFile* file = TFile::Open("/exp/sbnd/app/users/hlay/ncpizero/srcs/sbndcode/sbndcode/PiZero/ShowerEnergyCorrection/shower_energy_correction_hist_NCPiZeroAv12.root");
TProfile *fShowerEnergyCorrectionHist = (TProfile*) file->Get("hShowerEnergy2DRecoFractionalResolution_pfx");

void InitialiseTree(TChain *tree);

double CorrectEnergy(const double &energy);

void ObservablesResolution2(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution_2";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  InitialiseTree(ncpizeroEvents);

  const double pizeroMomBins[9] = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double cosThetaBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  TH1F *hInvariantMass = new TH1F("hInvariantMass", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);
  TH1F *hInvariantMassTrackDir = new TH1F("hInvariantMassTrackDir", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);
  TH1F *hInvariantMassEnergyCorr = new TH1F("hInvariantMassEnergyCorr", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);
  TH1F *hInvariantMassBoth = new TH1F("hInvariantMassBoth", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);

  TH1F *hPiZeroMom = new TH1F("hPiZeroMom", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomTrackDir = new TH1F("hPiZeroMomTrackDir", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomEnergyCorr = new TH1F("hPiZeroMomEnergyCorr", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomBoth = new TH1F("hPiZeroMomBoth", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomEn0Constr = new TH1F("hPiZeroMomEn0Constr", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomEn1Constr = new TH1F("hPiZeroMomEn1Constr", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomBothConstr = new TH1F("hPiZeroMomBothConstr", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomBestConstr = new TH1F("hPiZeroMomBestConstr", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionTrackDir = new TH1F("hPiZeroMomResolutionTrackDir", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionEnergyCorr = new TH1F("hPiZeroMomResolutionEnergyCorr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionBoth = new TH1F("hPiZeroMomResolutionBoth", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionEn0Constr = new TH1F("hPiZeroMomResolutionEn0Constr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionEn1Constr = new TH1F("hPiZeroMomResolutionEn1Constr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionBothConstr = new TH1F("hPiZeroMomResolutionBothConstr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionBestConstr = new TH1F("hPiZeroMomResolutionBestConstr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);

  TH1F *hCosThetaPiZero = new TH1F("hCosThetaPiZero", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroTrackDir = new TH1F("hCosThetaPiZeroTrackDir", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroEnergyCorr = new TH1F("hCosThetaPiZeroEnergyCorr", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroBoth = new TH1F("hCosThetaPiZeroBoth", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);

  TH1F *hCosThetaPiZeroResolution = new TH1F("hCosThetaPiZeroResolution", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionTrackDir = new TH1F("hCosThetaPiZeroResolutionTrackDir", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionEnergyCorr = new TH1F("hCosThetaPiZeroResolutionEnergyCorr", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionBoth = new TH1F("hCosThetaPiZeroResolutionBoth", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);

  const int N = ncpizeroEvents->GetEntries();
  
  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      ncpizeroEvents->GetEntry(ev_i);
      
      for(int slc_i = 0; slc_i < slc_true_event_type_incl->size(); ++slc_i)
        {
          if(slc_true_event_type_incl->at(slc_i) == 0 && slc_comp->at(slc_i) > .5 && slc_sel_incl->at(slc_i))
            {
              const double shwEn0 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              const double shwEn1 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));

              const double corrEn0 = CorrectEnergy(shwEn0);
              const double corrEn1 = CorrectEnergy(shwEn1);

              const TVector3 shwDir0 = TVector3(slc_pfp_shower_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_shower_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_shower_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              const TVector3 shwDir1 = TVector3(slc_pfp_shower_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_shower_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_shower_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              const TVector3 trkDir0 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              const TVector3 trkDir1 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              const double cosineThetaGammaGammaShw = shwDir0.Dot(shwDir1) / (shwDir0.Mag() * shwDir1.Mag());
              const double cosineThetaGammaGammaTrk = trkDir0.Dot(trkDir1) / (trkDir0.Mag() * trkDir1.Mag());

              const double invariantMass           = sqrt(2 * shwEn0 * shwEn1 * (1 - cosineThetaGammaGammaShw));
              const double invariantMassTrackDir   = sqrt(2 * shwEn0 * shwEn1 * (1 - cosineThetaGammaGammaTrk));
              const double invariantMassEnergyCorr = sqrt(2 * corrEn0 * corrEn1 * (1 - cosineThetaGammaGammaShw));
              const double invariantMassBoth       = sqrt(2 * corrEn0 * corrEn1 * (1 - cosineThetaGammaGammaTrk));

              const TVector3 pizeroMom           = shwEn0 * shwDir0 + shwEn1 * shwDir1;
              const TVector3 pizeroMomTrackDir   = shwEn0 * trkDir0 + shwEn1 * trkDir1;
              const TVector3 pizeroMomEnergyCorr = corrEn0 * shwDir0 + corrEn1 * shwDir1;
              const TVector3 pizeroMomBoth       = corrEn0 * trkDir0 + corrEn1 * trkDir1;

              const double corrEn0Constr = TMath::Power(134.9769, 2) / (2 * shwEn1 * (1 - cosineThetaGammaGammaTrk));
              const double corrEn1Constr = TMath::Power(134.9769, 2) / (2 * shwEn0 * (1 - cosineThetaGammaGammaTrk));

              /*
                const double pizeroMomSimp = corrEn0 + corrEn1;
                const double alpha         = std::abs(corrEn0 - corrEn1) / (corrEn0 + corrEn1);
                const double pizeroMomUB   = 134.9769 * TMath::Sqrt(2 / ((1 - TMath::Power(alpha, 2)) * (1 - cosineThetaGammaGammaTrk)));
              */

              const TVector3 pizeroMomEn0Constr  = corrEn0Constr * trkDir0 + corrEn1 * trkDir1;
              const TVector3 pizeroMomEn1Constr  = corrEn0 * trkDir0 + corrEn1Constr * trkDir1;
              const TVector3 pizeroMomBothConstr = corrEn0Constr * trkDir0 + corrEn1Constr * trkDir1;
              const TVector3 pizeroMomBestConstr = abs(corrEn0Constr - corrEn0) < abs(corrEn1Constr - corrEn1) ?
                                                                                  corrEn0Constr * trkDir0 + corrEn1 * trkDir1
                                                                                  : corrEn0 * trkDir0 + corrEn1Constr * trkDir1;

              const TVector3 zaxis = TVector3(0., 0., 1.);

              hInvariantMass->Fill(slc_best_pzc_invariant_mass->at(slc_i));
              hInvariantMassTrackDir->Fill(invariantMassTrackDir);
              hInvariantMassEnergyCorr->Fill(invariantMassEnergyCorr);
              hInvariantMassBoth->Fill(invariantMassBoth);

              hPiZeroMom->Fill(pizeroMom.Mag());
              hPiZeroMomTrackDir->Fill(pizeroMomTrackDir.Mag());
              hPiZeroMomEnergyCorr->Fill(pizeroMomEnergyCorr.Mag());
              hPiZeroMomBoth->Fill(pizeroMomBoth.Mag());
              hPiZeroMomEn0Constr->Fill(pizeroMomEn0Constr.Mag());
              hPiZeroMomEn1Constr->Fill(pizeroMomEn1Constr.Mag());
              hPiZeroMomBothConstr->Fill(pizeroMomBothConstr.Mag());
              hPiZeroMomBestConstr->Fill(pizeroMomBestConstr.Mag());

              hPiZeroMomResolution->Fill(pizeroMom.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionTrackDir->Fill(pizeroMomTrackDir.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionEnergyCorr->Fill(pizeroMomEnergyCorr.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionBoth->Fill(pizeroMomBoth.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionEn0Constr->Fill(pizeroMomEn0Constr.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionEn1Constr->Fill(pizeroMomEn1Constr.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionBothConstr->Fill(pizeroMomBothConstr.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionBestConstr->Fill(pizeroMomBestConstr.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));

              hCosThetaPiZero->Fill(cos(pizeroMom.Angle(zaxis)));
              hCosThetaPiZeroTrackDir->Fill(cos(pizeroMomTrackDir.Angle(zaxis)));
              hCosThetaPiZeroEnergyCorr->Fill(cos(pizeroMomEnergyCorr.Angle(zaxis)));
              hCosThetaPiZeroBoth->Fill(cos(pizeroMomBoth.Angle(zaxis)));

              hCosThetaPiZeroResolution->Fill(cos(pizeroMom.Angle(zaxis)) - slc_true_pz_cos_theta_pizero->at(slc_i).at(0));
              hCosThetaPiZeroResolutionTrackDir->Fill(cos(pizeroMomTrackDir.Angle(zaxis)) - slc_true_pz_cos_theta_pizero->at(slc_i).at(0));
              hCosThetaPiZeroResolutionEnergyCorr->Fill(cos(pizeroMomEnergyCorr.Angle(zaxis)) - slc_true_pz_cos_theta_pizero->at(slc_i).at(0));
              hCosThetaPiZeroResolutionBoth->Fill(cos(pizeroMomBoth.Angle(zaxis)) - slc_true_pz_cos_theta_pizero->at(slc_i).at(0));
            }
        }
    }

  TCanvas *cInvariantMass = new TCanvas("cInvariantMass", "cInvariantMass");
  cInvariantMass->cd();

  hInvariantMass->SetLineColor(kMagenta+2);
  hInvariantMassTrackDir->SetLineColor(kCyan+2);
  hInvariantMassEnergyCorr->SetLineColor(kGreen+2);
  hInvariantMassBoth->SetLineColor(kOrange+2);

  hInvariantMassTrackDir->Draw();
  hInvariantMass->Draw("same");
  hInvariantMassEnergyCorr->Draw("same");
  hInvariantMassBoth->Draw("same");

  TLine *nominalLine = new TLine();
  nominalLine->SetLineColor(kRed+2);
  nominalLine->SetLineWidth(5);
  nominalLine->DrawLine(134.9769, 0., 134.9769, 1.05 * hInvariantMass->GetMaximum());

  TLegend *lInvariantMass = new TLegend(.6, .4, .85, .6);
  lInvariantMass->AddEntry(hInvariantMass, "Standard", "l");
  lInvariantMass->AddEntry(hInvariantMassTrackDir, "Track Direction", "l");
  lInvariantMass->AddEntry(hInvariantMassEnergyCorr, "Energy Correction", "l");
  lInvariantMass->AddEntry(hInvariantMassBoth, "Both", "l");
  lInvariantMass->Draw();

  cInvariantMass->SaveAs(saveDir + "/invariant_mass.png");
  cInvariantMass->SaveAs(saveDir + "/invariant_mass.pdf");

  TCanvas *cPiZeroMom = new TCanvas("cPiZeroMom", "cPiZeroMom");
  cPiZeroMom->cd();

  NormaliseEntriesByBinWidth(hPiZeroMom);
  NormaliseEntriesByBinWidth(hPiZeroMomTrackDir);
  NormaliseEntriesByBinWidth(hPiZeroMomEnergyCorr);
  NormaliseEntriesByBinWidth(hPiZeroMomBoth);
  NormaliseEntriesByBinWidth(hPiZeroMomEn0Constr);
  NormaliseEntriesByBinWidth(hPiZeroMomEn1Constr);
  NormaliseEntriesByBinWidth(hPiZeroMomBothConstr);
  NormaliseEntriesByBinWidth(hPiZeroMomBestConstr);

  hPiZeroMom->SetLineColor(kMagenta+2);
  hPiZeroMomTrackDir->SetLineColor(kCyan+2);
  hPiZeroMomEnergyCorr->SetLineColor(kGreen+2);
  hPiZeroMomBoth->SetLineColor(kOrange+2);
  hPiZeroMomEn0Constr->SetLineColor(kRed+2);
  hPiZeroMomEn1Constr->SetLineColor(kBlue+2);
  hPiZeroMomBothConstr->SetLineColor(kBlack);
  hPiZeroMomBestConstr->SetLineColor(kGray+2);

  hPiZeroMomTrackDir->Draw("hist");
  hPiZeroMom->Draw("histsame");
  hPiZeroMomEnergyCorr->Draw("histsame");
  hPiZeroMomBoth->Draw("histsame");
  hPiZeroMomEn0Constr->Draw("histsame");
  hPiZeroMomEn1Constr->Draw("histsame");
  hPiZeroMomBothConstr->Draw("histsame");
  hPiZeroMomBestConstr->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.png");
  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.pdf");

  TCanvas *cPiZeroMomResolution = new TCanvas("cPiZeroMomResolution", "cPiZeroMomResolution");
  cPiZeroMomResolution->cd();

  hPiZeroMomResolution->SetLineColor(kMagenta+2);
  hPiZeroMomResolutionTrackDir->SetLineColor(kCyan+2);
  hPiZeroMomResolutionEnergyCorr->SetLineColor(kGreen+2);
  hPiZeroMomResolutionBoth->SetLineColor(kOrange+2);
  hPiZeroMomResolutionEn0Constr->SetLineColor(kRed+2);
  hPiZeroMomResolutionEn1Constr->SetLineColor(kBlue+2);
  hPiZeroMomResolutionBothConstr->SetLineColor(kBlack);
  hPiZeroMomResolutionBestConstr->SetLineColor(kGray+2);

  hPiZeroMomResolutionTrackDir->Draw("hist");
  hPiZeroMomResolution->Draw("histsame");
  hPiZeroMomResolutionEnergyCorr->Draw("histsame");
  hPiZeroMomResolutionBoth->Draw("histsame");
  hPiZeroMomResolutionEn0Constr->Draw("histsame");
  hPiZeroMomResolutionEn1Constr->Draw("histsame");
  hPiZeroMomResolutionBothConstr->Draw("histsame");
  hPiZeroMomResolutionBestConstr->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_mom_resolution.png");
  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_mom_resolution.pdf");

  TCanvas *cCosThetaPiZero = new TCanvas("cCosThetaPiZero", "cCosThetaPiZero");
  cCosThetaPiZero->cd();

  NormaliseEntriesByBinWidth(hCosThetaPiZero);
  NormaliseEntriesByBinWidth(hCosThetaPiZeroTrackDir);
  NormaliseEntriesByBinWidth(hCosThetaPiZeroEnergyCorr);
  NormaliseEntriesByBinWidth(hCosThetaPiZeroBoth);

  hCosThetaPiZero->SetLineColor(kMagenta+2);
  hCosThetaPiZeroTrackDir->SetLineColor(kCyan+2);
  hCosThetaPiZeroEnergyCorr->SetLineColor(kGreen+2);
  hCosThetaPiZeroBoth->SetLineColor(kOrange+2);

  hCosThetaPiZeroTrackDir->Draw("hist");
  hCosThetaPiZero->Draw("histsame");
  hCosThetaPiZeroEnergyCorr->Draw("histsame");
  hCosThetaPiZeroBoth->Draw("histsame");

  lInvariantMass->Draw();

  cCosThetaPiZero->SaveAs(saveDir + "/cos_theta_pizero.png");
  cCosThetaPiZero->SaveAs(saveDir + "/cos_theta_pizero.pdf");

  TCanvas *cCosThetaPiZeroResolution = new TCanvas("cCosThetaPiZeroResolution", "cCosThetaPiZeroResolution");
  cCosThetaPiZeroResolution->cd();

  hCosThetaPiZeroResolution->SetLineColor(kMagenta+2);
  hCosThetaPiZeroResolutionTrackDir->SetLineColor(kCyan+2);
  hCosThetaPiZeroResolutionEnergyCorr->SetLineColor(kGreen+2);
  hCosThetaPiZeroResolutionBoth->SetLineColor(kOrange+2);

  hCosThetaPiZeroResolutionTrackDir->Draw("hist");
  hCosThetaPiZeroResolution->Draw("histsame");
  hCosThetaPiZeroResolutionEnergyCorr->Draw("histsame");
  hCosThetaPiZeroResolutionBoth->Draw("histsame");

  lInvariantMass->Draw();

  cCosThetaPiZeroResolution->SaveAs(saveDir + "/cos_theta_pizero_resolution.png");
  cCosThetaPiZeroResolution->SaveAs(saveDir + "/cos_theta_pizero_resolution.pdf");
}

void InitialiseTree(TChain *tree)
{
  tree->SetBranchStatus("*", 0);

  tree->SetBranchAddress("slc_true_event_type_incl", &slc_true_event_type_incl);
  tree->SetBranchAddress("slc_comp", &slc_comp);
  tree->SetBranchAddress("slc_sel_incl", &slc_sel_incl);
  tree->SetBranchAddress("slc_best_pzc_pizero_mom", &slc_best_pzc_pizero_mom);
  tree->SetBranchAddress("slc_best_pzc_invariant_mass", &slc_best_pzc_invariant_mass);
  tree->SetBranchAddress("slc_best_pzc_cos_theta_pizero", &slc_best_pzc_cos_theta_pizero);
  tree->SetBranchAddress("slc_best_pzc_photon_0_id", &slc_best_pzc_photon_0_id);
  tree->SetBranchAddress("slc_best_pzc_photon_1_id", &slc_best_pzc_photon_1_id);
  tree->SetBranchAddress("slc_true_pz_pizero_mom", &slc_true_pz_pizero_mom);
  tree->SetBranchAddress("slc_true_pz_cos_theta_pizero", &slc_true_pz_cos_theta_pizero);
  tree->SetBranchAddress("slc_pfp_shower_dir_x", &slc_pfp_shower_dir_x);
  tree->SetBranchAddress("slc_pfp_shower_dir_y", &slc_pfp_shower_dir_y);
  tree->SetBranchAddress("slc_pfp_shower_dir_z", &slc_pfp_shower_dir_z);
  tree->SetBranchAddress("slc_pfp_track_dir_x", &slc_pfp_track_dir_x);
  tree->SetBranchAddress("slc_pfp_track_dir_y", &slc_pfp_track_dir_y);
  tree->SetBranchAddress("slc_pfp_track_dir_z", &slc_pfp_track_dir_z);
  tree->SetBranchAddress("slc_pfp_shower_energy", &slc_pfp_shower_energy);
}

double CorrectEnergy(const double &energy)
{
  const int bin = fShowerEnergyCorrectionHist->FindBin(energy);

  return energy * (1 - fShowerEnergyCorrectionHist->GetBinContent(bin));
}
