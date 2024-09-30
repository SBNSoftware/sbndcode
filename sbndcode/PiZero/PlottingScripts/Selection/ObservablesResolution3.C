#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "FittedObservables.C"

void ObservablesResolution3(const TString productionVersion)
{
  std::vector<double> covMatrix = GetCovMatrix(productionVersion, true);

  std::cout << "std::vector<double> covMatrix = {";
  for(int i = 0; i < covMatrix.size(); ++i)
    {
      std::cout << covMatrix[i];

      if(i < covMatrix.size() - 1)
        std::cout << ", ";
    }
  std::cout << "};" << std::endl;

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution_3";
  gSystem->Exec("mkdir -p " + saveDir);

  //  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";
  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  InitialiseTree(ncpizeroEvents);

  const double pizeroMomBins[9] = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double showerEnBins[13] = { 0., 30., 60., 90., 120., 150., 180., 220., 260., 300., 400., 600., 1000. };
  const double cosThetaBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  TH1F *hInvariantMass = new TH1F("hInvariantMass", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);
  TH1F *hInvariantMassCorr = new TH1F("hInvariantMassCorr", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);
  TH1F *hInvariantMassFitted = new TH1F("hInvariantMassFitted", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);

  TH1F *hPiZeroMom = new TH1F("hPiZeroMom", ";p_{#pi^{0}} (MeV/c);#pi^{0} / 60 MeV/c", 8, pizeroMomBins);
  TH1F *hPiZeroMomCorr = new TH1F("hPiZeroMomCorr", ";p_{#pi^{0}} (MeV/c);#pi^{0} / 60 MeV/c", 8, pizeroMomBins);
  TH1F *hPiZeroMomFitted = new TH1F("hPiZeroMomFitted", ";p_{#pi^{0}} (MeV/c);#pi^{0} / 60 MeV/c", 8, pizeroMomBins);

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionCorr = new TH1F("hPiZeroMomResolutionCorr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionFitted = new TH1F("hPiZeroMomResolutionFitted", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);

  TH1F *hLeadingShowerEn = new TH1F("hLeadingShowerEn", ";E_{#gamma_{1}} (MeV);#pi^{0} / 30 MeV", 8, pizeroMomBins);
  TH1F *hLeadingShowerEnCorr = new TH1F("hLeadingShowerEnCorr", ";E_{#gamma_{1}} (MeV);#pi^{0} / 30 MeV", 8, pizeroMomBins);
  TH1F *hLeadingShowerEnFitted = new TH1F("hLeadingShowerEnFitted", ";E_{#gamma_{1}} (MeV);#pi^{0} / 30 MeV", 8, pizeroMomBins);

  TH1F *hLeadingShowerEnResolution = new TH1F("hLeadingShowerEnResolution", ";E_{#gamma_{1}} (Reco - True) (MeV);#pi^{0}", 35, -310, 390);
  TH1F *hLeadingShowerEnResolutionCorr = new TH1F("hLeadingShowerEnResolutionCorr", ";E_{#gamma_{1}} (Reco - True) (MeV);#pi^{0}", 35, -310, 390);
  TH1F *hLeadingShowerEnResolutionFitted = new TH1F("hLeadingShowerEnResolutionFitted", ";E_{#gamma_{1}} (Reco - True) (MeV);#pi^{0}", 35, -310, 390);

  TH1F *hSubLeadingShowerEn = new TH1F("hSubLeadingShowerEn", ";E_{#gamma_{2}} (MeV);#pi^{0} / 30 MeV", 8, pizeroMomBins);
  TH1F *hSubLeadingShowerEnCorr = new TH1F("hSubLeadingShowerEnCorr", ";E_{#gamma_{2}} (MeV);#pi^{0} / 30 MeV", 8, pizeroMomBins);
  TH1F *hSubLeadingShowerEnFitted = new TH1F("hSubLeadingShowerEnFitted", ";E_{#gamma_{2}} (MeV);#pi^{0} / 30 MeV", 8, pizeroMomBins);

  TH1F *hSubLeadingShowerEnResolution = new TH1F("hSubLeadingShowerEnResolution", ";E_{#gamma_{2}} (Reco - True) (MeV);#pi^{0}", 35, -310, 390);
  TH1F *hSubLeadingShowerEnResolutionCorr = new TH1F("hSubLeadingShowerEnResolutionCorr", ";E_{#gamma_{2}} (Reco - True) (MeV);#pi^{0}", 35, -310, 390);
  TH1F *hSubLeadingShowerEnResolutionFitted = new TH1F("hSubLeadingShowerEnResolutionFitted", ";E_{#gamma_{2}} (Reco - True) (MeV);#pi^{0}", 35, -310, 390);

  TH1F *hLeadingShowerDirResolution = new TH1F("hLeadingShowerDirResolution", ";#theta_{True_{1} Reco_{1}} (#circ);#pi^{0}", 36, 0, 180);
  TH1F *hLeadingShowerDirResolutionCorr = new TH1F("hLeadingShowerDirResolutionCorr", ";#theta_{True_{1} Reco_{1}} (#circ);#pi^{0}", 36, 0, 180);
  TH1F *hLeadingShowerDirResolutionFitted = new TH1F("hLeadingShowerDirResolutionFitted", ";#theta_{True_{1} Reco_{1}} (#circ);#pi^{0}", 36, 0, 180);

  TH1F *hSubLeadingShowerDirResolution = new TH1F("hSubLeadingShowerDirResolution", ";#theta_{True_{2} Reco_{2}} (#circ);#pi^{0}", 36, 0, 180);
  TH1F *hSubLeadingShowerDirResolutionCorr = new TH1F("hSubLeadingShowerDirResolutionCorr", ";#theta_{True_{2} Reco_{2}} (#circ);#pi^{0}", 36, 0, 180);
  TH1F *hSubLeadingShowerDirResolutionFitted = new TH1F("hSubLeadingShowerDirResolutionFitted", ";#theta_{True_{2} Reco_{2}} (#circ);#pi^{0}", 36, 0, 180);

  TH1F *hCosThetaPiZero = new TH1F("hCosThetaPiZero", ";cos(#theta_{#pi^{0}});#pi^{0} / 60 MeV/c", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroCorr = new TH1F("hCosThetaPiZeroCorr", ";cos(#theta_{#pi^{0}});#pi^{0} / 60 MeV/c", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroFitted = new TH1F("hCosThetaPiZeroFitted", ";cos(#theta_{#pi^{0}});#pi^{0} / 60 MeV/c", 9, cosThetaBins);

  TH1F *hCosThetaPiZeroResolution = new TH1F("hCosThetaPiZeroResolution", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 40, -1, 1);
  TH1F *hCosThetaPiZeroResolutionCorr = new TH1F("hCosThetaPiZeroResolutionCorr", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 40, -1, 1);
  TH1F *hCosThetaPiZeroResolutionFitted = new TH1F("hCosThetaPiZeroResolutionFitted", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 40, -1, 1);

  TH1F *hPiZeroOpenAngleResolution = new TH1F("hPiZeroOpenAngleResolution", ";#theta_{#gamma_{1}#gamma_{2}} (Reco - True) (#circ);#pi^{0}", 33, -150, 180);
  TH1F *hPiZeroOpenAngleResolutionCorr = new TH1F("hPiZeroOpenAngleResolutionCorr", ";#theta_{#gamma_{1}#gamma_{2}} (Reco - True) (#circ);#pi^{0}", 33, -150, 180);
  TH1F *hPiZeroOpenAngleResolutionFitted = new TH1F("hPiZeroOpenAngleResolutionFitted", ";#theta_{#gamma_{1}#gamma_{2}} (Reco - True) (#circ);#pi^{0}", 33, -150, 180);

  const int N = ncpizeroEvents->GetEntries();

  int nsig = 0, nsigGood = 0;
  
  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      ncpizeroEvents->GetEntry(ev_i);
      
      for(int slc_i = 0; slc_i < slc_true_event_type_incl->size(); ++slc_i)
        {
          if(slc_true_event_type_incl->at(slc_i) == 0 && slc_comp->at(slc_i) > .5 && slc_sel_incl->at(slc_i))
            {
              const bool goodReco = ((slc_best_pzc_photon_0_true_trackid->at(slc_i) == slc_true_pz_gamma0_trackid->at(slc_i).at(0) &&
                                      slc_best_pzc_photon_1_true_trackid->at(slc_i) == slc_true_pz_gamma1_trackid->at(slc_i).at(0)) ||
                                     (slc_best_pzc_photon_0_true_trackid->at(slc_i) == slc_true_pz_gamma1_trackid->at(slc_i).at(0) &&
                                      slc_best_pzc_photon_1_true_trackid->at(slc_i) == slc_true_pz_gamma0_trackid->at(slc_i).at(0)))
                && slc_best_pzc_photon_0_comp->at(slc_i) > .8 && slc_best_pzc_photon_1_comp->at(slc_i) > .8
                && slc_best_pzc_photon_0_pur->at(slc_i) > .8 && slc_best_pzc_photon_1_pur->at(slc_i) > .8;

              /*
                if(!goodReco)
                continue;
              */

              ++nsig;

              double en0 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              double en1 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));

              const bool invert = en0 < en1;

              double corrEn0 = CorrectEnergy(en0);
              double corrEn1 = CorrectEnergy(en1);

              TVector3 dir0 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                       slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                       slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              TVector3 dir1 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                       slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                       slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              TVector3 shwDir0 = TVector3(slc_pfp_shower_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                          slc_pfp_shower_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                          slc_pfp_shower_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              TVector3 shwDir1 = TVector3(slc_pfp_shower_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                          slc_pfp_shower_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                          slc_pfp_shower_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              TVector3 trueDir0 = TVector3(slc_pfp_true_p_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                           slc_pfp_true_p_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                           slc_pfp_true_p_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              TVector3 trueDir1 = TVector3(slc_pfp_true_p_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                           slc_pfp_true_p_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                           slc_pfp_true_p_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              if(invert)
                {
                  std::swap(en0, en1);
                  std::swap(corrEn0, corrEn1);
                  std::swap(dir0, dir1);
                  std::swap(shwDir0, shwDir1);
                }

              const double cosineThetaGammaGamma = dir0.Dot(dir1) / (dir0.Mag() * dir1.Mag());
              const double thetaGammaGamma       = acos(cosineThetaGammaGamma);
              const double shwThetaGammaGamma    = shwDir0.Angle(shwDir1);
              const double trueThetaGammaGamma   = trueDir0.Angle(trueDir1);

              const double invariantMass = sqrt(2 * corrEn0 * corrEn1 * (1 - cosineThetaGammaGamma));
              const TVector3 pizeroMom   = corrEn0 * dir0 + corrEn1 * dir1;

              const TVector3 zaxis = TVector3(0., 0., 1.);

              bool good = false;
              std::vector<double> updated = DoKF(corrEn0, corrEn1, thetaGammaGamma, covMatrix, good);

              if(good)
                ++nsigGood;

              const double invariantMassFitted = TMath::Sqrt(2 * updated[0] * updated[1] * (1 - cos(updated[2])));
              const TVector3 pizeroMomFitted   = updated[0] * dir0 + updated[1] * dir1;

              const double pizeroMomFitted2    = TMath::Sqrt(TMath::Power(updated[0] + updated[1], 2) - TMath::Power(kPiZeroMass, 2));

              const float updatedEn0Proportion   = updated[0] / corrEn0;
              const float updatedEn1Proportion   = updated[1] / corrEn1;
              const float updatedAngleProportion = updated[2] / thetaGammaGamma;

              const float fitFactor0             = (updatedEn1Proportion / (updatedEn0Proportion + updatedEn1Proportion)) * (updatedAngleProportion - 1) + 1;
              const float fitFactor1             = (updatedEn0Proportion / (updatedEn0Proportion + updatedEn1Proportion)) * (updatedAngleProportion - 1) + 1;
              TVector3 fittedDir0 = fitFactor0 * (dir0 - dir1) + dir1;
              TVector3 fittedDir1 = fitFactor1 * (dir1 - dir0) + dir0;

              if(updated[2] - fittedDir1.Angle(fittedDir0) > 0.01)
                std::cout << updated[2] << " " << fittedDir1.Angle(fittedDir0) << std::endl;

              hInvariantMass->Fill(slc_best_pzc_invariant_mass->at(slc_i));
              hInvariantMassCorr->Fill(invariantMass);
              hInvariantMassFitted->Fill(invariantMassFitted);

              hPiZeroMom->Fill(slc_best_pzc_pizero_mom->at(slc_i));
              hPiZeroMomCorr->Fill(pizeroMom.Mag());
              hPiZeroMomFitted->Fill(pizeroMomFitted2);

              hPiZeroMomResolution->Fill(slc_best_pzc_pizero_mom->at(slc_i) - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionCorr->Fill(pizeroMom.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionFitted->Fill(pizeroMomFitted2 - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));

              hLeadingShowerEn->Fill(en0);
              hLeadingShowerEnCorr->Fill(corrEn0);
              hLeadingShowerEnFitted->Fill(updated[0]);

              hLeadingShowerEnResolution->Fill(en0 - 1e3 * slc_true_pz_gamma0_energy->at(slc_i).at(0));
              hLeadingShowerEnResolutionCorr->Fill(corrEn0 - 1e3 * slc_true_pz_gamma0_energy->at(slc_i).at(0));
              hLeadingShowerEnResolutionFitted->Fill(updated[0] - 1e3 * slc_true_pz_gamma0_energy->at(slc_i).at(0));

              hSubLeadingShowerEn->Fill(en1);
              hSubLeadingShowerEnCorr->Fill(corrEn1);
              hSubLeadingShowerEnFitted->Fill(updated[1]);

              hSubLeadingShowerEnResolution->Fill(en1 - 1e3 * slc_true_pz_gamma1_energy->at(slc_i).at(0));
              hSubLeadingShowerEnResolutionCorr->Fill(corrEn1 - 1e3 * slc_true_pz_gamma1_energy->at(slc_i).at(0));
              hSubLeadingShowerEnResolutionFitted->Fill(updated[1] - 1e3 * slc_true_pz_gamma1_energy->at(slc_i).at(0));

              hLeadingShowerDirResolution->Fill(TMath::RadToDeg() * shwDir0.Angle(trueDir0));
              hLeadingShowerDirResolutionCorr->Fill(TMath::RadToDeg() * dir0.Angle(trueDir0));
              hLeadingShowerDirResolutionFitted->Fill(TMath::RadToDeg() * fittedDir0.Angle(trueDir0));

              hSubLeadingShowerDirResolution->Fill(TMath::RadToDeg() * shwDir1.Angle(trueDir1));
              hSubLeadingShowerDirResolutionCorr->Fill(TMath::RadToDeg() * dir1.Angle(trueDir1));
              hSubLeadingShowerDirResolutionFitted->Fill(TMath::RadToDeg() * fittedDir1.Angle(trueDir1));

              hPiZeroOpenAngleResolution->Fill(TMath::RadToDeg() * (shwThetaGammaGamma - trueThetaGammaGamma));
              hPiZeroOpenAngleResolutionCorr->Fill(TMath::RadToDeg() * (thetaGammaGamma - trueThetaGammaGamma));
              hPiZeroOpenAngleResolutionFitted->Fill(TMath::RadToDeg() * (updated[2] - trueThetaGammaGamma));
            }
        }
    }

  std::cout << "Fitting good: " << nsigGood << " / " << nsig << " (" << (100. * nsigGood) / nsig << "%)" << std::endl;

  TCanvas *cInvariantMass = new TCanvas("cInvariantMass", "cInvariantMass");
  cInvariantMass->cd();

  hInvariantMass->SetLineColor(kMagenta+2);
  hInvariantMassCorr->SetLineColor(kOrange+2);
  hInvariantMassFitted->SetLineColor(kCyan+2);

  hInvariantMassFitted->Draw();
  hInvariantMassCorr->Draw("same");
  hInvariantMass->Draw("same");

  TLine *nominalLine = new TLine();
  nominalLine->SetLineColor(kRed+2);
  nominalLine->SetLineWidth(5);
  nominalLine->DrawLine(kPiZeroMass, 0., kPiZeroMass, 1.05 * hInvariantMassFitted->GetMaximum());

  TLegend *lInvariantMass = new TLegend(.6, .4, .85, .6);
  lInvariantMass->AddEntry(hInvariantMass, "Original", "l");
  lInvariantMass->AddEntry(hInvariantMassCorr, "Corrected", "l");
  lInvariantMass->AddEntry(hInvariantMassFitted, "Kinematic Fitting", "l");
  lInvariantMass->Draw();

  cInvariantMass->SaveAs(saveDir + "/invariant_mass.png");
  cInvariantMass->SaveAs(saveDir + "/invariant_mass.pdf");

  TCanvas *cPiZeroMom = new TCanvas("cPiZeroMom", "cPiZeroMom");
  cPiZeroMom->cd();

  NormaliseEntriesByBinWidth(hPiZeroMom, 60);
  NormaliseEntriesByBinWidth(hPiZeroMomCorr, 60);
  NormaliseEntriesByBinWidth(hPiZeroMomFitted, 60);

  hPiZeroMom->SetLineColor(kMagenta+2);
  hPiZeroMomCorr->SetLineColor(kOrange+2);
  hPiZeroMomFitted->SetLineColor(kCyan+2);

  hPiZeroMom->Draw("hist");
  hPiZeroMomCorr->Draw("histsame");
  hPiZeroMomFitted->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.png");
  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.pdf");

  TCanvas *cPiZeroMomResolution = new TCanvas("cPiZeroMomResolution", "cPiZeroMomResolution");
  cPiZeroMomResolution->cd();

  hPiZeroMomResolution->SetLineColor(kMagenta+2);
  hPiZeroMomResolutionCorr->SetLineColor(kOrange+2);
  hPiZeroMomResolutionFitted->SetLineColor(kCyan+2);

  hPiZeroMomResolutionCorr->Draw("histsame");
  hPiZeroMomResolution->Draw("histsame");
  hPiZeroMomResolutionFitted->Draw("histsame");

  TPaveText *text = new TPaveText(.6, .65, .68, .75, "NDC");
  text->AddText("Standard");
  text->AddText("Corrected");
  text->AddText("Fitted");
  text->SetTextAlign(12);
  text->SetTextSize(0.02);
  text->SetBorderSize(0);
  text->SetFillColor(kWhite);
  text->Draw();

  TPaveText *text2 = new TPaveText(.68, .65, .85, .75, "NDC");
  text2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolution->GetMean(), hPiZeroMomResolution->GetStdDev()));
  text2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolutionCorr->GetMean(), hPiZeroMomResolutionCorr->GetStdDev()));
  text2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolutionFitted->GetMean(), hPiZeroMomResolutionFitted->GetStdDev()));
  text2->SetTextAlign(12);
  text2->SetTextSize(0.02);
  text2->SetBorderSize(0);
  text2->SetFillColor(kWhite);
  text2->Draw();

  lInvariantMass->Draw();

  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_mom_resolution.png");
  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_mom_resolution.pdf");

  TCanvas *cLeadingShowerEn = new TCanvas("cLeadingShowerEn", "cLeadingShowerEn");
  cLeadingShowerEn->cd();

  NormaliseEntriesByBinWidth(hLeadingShowerEn, 30);
  NormaliseEntriesByBinWidth(hLeadingShowerEnCorr, 30);
  NormaliseEntriesByBinWidth(hLeadingShowerEnFitted, 30);

  hLeadingShowerEn->SetLineColor(kMagenta+2);
  hLeadingShowerEnCorr->SetLineColor(kOrange+2);
  hLeadingShowerEnFitted->SetLineColor(kCyan+2);

  hLeadingShowerEnFitted->Draw("histsame");
  hLeadingShowerEn->Draw("histsame");
  hLeadingShowerEnCorr->Draw("histsame");

  lInvariantMass->Draw();

  cLeadingShowerEn->SaveAs(saveDir + "/leading_shower_en.png");
  cLeadingShowerEn->SaveAs(saveDir + "/leading_shower_en.pdf");

  TCanvas *cLeadingShowerEnResolution = new TCanvas("cLeadingShowerEnResolution", "cLeadingShowerEnResolution");
  cLeadingShowerEnResolution->cd();

  hLeadingShowerEnResolution->SetLineColor(kMagenta+2);
  hLeadingShowerEnResolutionCorr->SetLineColor(kOrange+2);
  hLeadingShowerEnResolutionFitted->SetLineColor(kCyan+2);

  hLeadingShowerEnResolutionCorr->Draw("hist");
  hLeadingShowerEnResolution->Draw("histsame");
  hLeadingShowerEnResolutionFitted->Draw("histsame");

  lInvariantMass->Draw();

  cLeadingShowerEnResolution->SaveAs(saveDir + "/leading_shower_en_resolution.png");
  cLeadingShowerEnResolution->SaveAs(saveDir + "/leading_shower_en_resolution.pdf");

  TCanvas *cSubLeadingShowerEn = new TCanvas("cSubLeadingShowerEn", "cSubLeadingShowerEn");
  cSubLeadingShowerEn->cd();

  NormaliseEntriesByBinWidth(hSubLeadingShowerEn, 30);
  NormaliseEntriesByBinWidth(hSubLeadingShowerEnCorr, 30);
  NormaliseEntriesByBinWidth(hSubLeadingShowerEnFitted, 30);

  hSubLeadingShowerEn->SetLineColor(kMagenta+2);
  hSubLeadingShowerEnCorr->SetLineColor(kOrange+2);
  hSubLeadingShowerEnFitted->SetLineColor(kCyan+2);

  hSubLeadingShowerEn->Draw("hist");
  hSubLeadingShowerEnCorr->Draw("histsame");
  hSubLeadingShowerEnFitted->Draw("histsame");

  lInvariantMass->Draw();

  cSubLeadingShowerEn->SaveAs(saveDir + "/subleading_shower_en.png");
  cSubLeadingShowerEn->SaveAs(saveDir + "/subleading_shower_en.pdf");

  TCanvas *cSubLeadingShowerEnResolution = new TCanvas("cSubLeadingShowerEnResolution", "cSubLeadingShowerEnResolution");
  cSubLeadingShowerEnResolution->cd();

  hSubLeadingShowerEnResolution->SetLineColor(kMagenta+2);
  hSubLeadingShowerEnResolutionCorr->SetLineColor(kOrange+2);
  hSubLeadingShowerEnResolutionFitted->SetLineColor(kCyan+2);

  hSubLeadingShowerEnResolutionCorr->Draw("hist");
  hSubLeadingShowerEnResolution->Draw("histsame");
  hSubLeadingShowerEnResolutionFitted->Draw("histsame");

  lInvariantMass->Draw();

  cSubLeadingShowerEnResolution->SaveAs(saveDir + "/subleading_shower_en_resolution.png");
  cSubLeadingShowerEnResolution->SaveAs(saveDir + "/subleading_shower_en_resolution.pdf");

  TCanvas *cLeadingShowerDirResolution = new TCanvas("cLeadingShowerDirResolution", "cLeadingShowerDirResolution");
  cLeadingShowerDirResolution->cd();

  hLeadingShowerDirResolution->SetLineColor(kMagenta+2);
  hLeadingShowerDirResolutionCorr->SetLineColor(kOrange+2);
  hLeadingShowerDirResolutionFitted->SetLineColor(kCyan+2);

  hLeadingShowerDirResolutionCorr->Draw("hist");
  hLeadingShowerDirResolution->Draw("histsame");
  hLeadingShowerDirResolutionFitted->Draw("histsame");

  lInvariantMass->Draw();

  cLeadingShowerDirResolution->SaveAs(saveDir + "/leading_shower_dir_resolution.png");
  cLeadingShowerDirResolution->SaveAs(saveDir + "/leading_shower_dir_resolution.pdf");

  TCanvas *cSubLeadingShowerDirResolution = new TCanvas("cSubLeadingShowerDirResolution", "cSubLeadingShowerDirResolution");
  cSubLeadingShowerDirResolution->cd();

  hSubLeadingShowerDirResolution->SetLineColor(kMagenta+2);
  hSubLeadingShowerDirResolutionCorr->SetLineColor(kOrange+2);
  hSubLeadingShowerDirResolutionFitted->SetLineColor(kCyan+2);

  hSubLeadingShowerDirResolutionFitted->Draw("histsame");
  hSubLeadingShowerDirResolutionCorr->Draw("histsame");
  hSubLeadingShowerDirResolution->Draw("histsame");

  lInvariantMass->Draw();

  cSubLeadingShowerDirResolution->SaveAs(saveDir + "/subleading_shower_dir_resolution.png");
  cSubLeadingShowerDirResolution->SaveAs(saveDir + "/subleading_shower_dir_resolution.pdf");

  TCanvas *cPiZeroOpenAngleResolution = new TCanvas("cPiZeroOpenAngleResolution", "cPiZeroOpenAngleResolution");
  cPiZeroOpenAngleResolution->cd();

  hPiZeroOpenAngleResolution->SetLineColor(kMagenta+2);
  hPiZeroOpenAngleResolutionCorr->SetLineColor(kOrange+2);
  hPiZeroOpenAngleResolutionFitted->SetLineColor(kCyan+2);

  hPiZeroOpenAngleResolutionFitted->Draw("histsame");
  hPiZeroOpenAngleResolutionCorr->Draw("histsame");
  hPiZeroOpenAngleResolution->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroOpenAngleResolution->SaveAs(saveDir + "/pizero_opening_angle_resolution.png");
  cPiZeroOpenAngleResolution->SaveAs(saveDir + "/pizero_opening_angle_resolution.pdf");
}
