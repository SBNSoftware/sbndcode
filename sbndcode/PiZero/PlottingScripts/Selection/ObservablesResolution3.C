#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "FittedObservables.C"

void ObservablesResolution3(const TString productionVersion)
{
  std::vector<double> covMatrix = GetCovMatrix(productionVersion, true);

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution_3";
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
  TH1F *hInvariantMassCorr = new TH1F("hInvariantMassCorr", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);
  TH1F *hInvariantMassFitted = new TH1F("hInvariantMassFitted", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);

  TH1F *hPiZeroMom = new TH1F("hPiZeroMom", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomCorr = new TH1F("hPiZeroMomCorr", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomFitted = new TH1F("hPiZeroMomFitted", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionCorr = new TH1F("hPiZeroMomResolutionCorr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionFitted = new TH1F("hPiZeroMomResolutionFitted", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);

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

              const double en0 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              const double en1 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));

              double corrEn0 = CorrectEnergy(en0);
              double corrEn1 = CorrectEnergy(en1);

              const TVector3 dir0 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                             slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                             slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              const TVector3 dir1 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                             slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                             slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              const double cosineThetaGammaGamma = dir0.Dot(dir1) / (dir0.Mag() * dir1.Mag());
              const double thetaGammaGamma       = acos(cosineThetaGammaGamma);

              const double invariantMass = sqrt(2 * corrEn0 * corrEn1 * (1 - cosineThetaGammaGamma));
              const TVector3 pizeroMom   = corrEn0 * dir0 + corrEn1 * dir1;

              const TVector3 zaxis = TVector3(0., 0., 1.);

              bool good = false;
              std::vector<double> updated = DoKF(corrEn0, corrEn1, thetaGammaGamma, covMatrix, good);

              if(good)
                ++nsigGood;

              const double invariantMassFitted = TMath::Sqrt(2 * updated[0] * updated[1] * (1 - cos(updated[2])));
              const double pizeroMomFitted     = TMath::Sqrt(TMath::Power(updated[0] + updated[1], 2) - TMath::Power(kPiZeroMass, 2));

              hInvariantMass->Fill(slc_best_pzc_invariant_mass->at(slc_i));
              hInvariantMassCorr->Fill(invariantMass);
              hInvariantMassFitted->Fill(invariantMassFitted);

              hPiZeroMom->Fill(slc_best_pzc_pizero_mom->at(slc_i));
              hPiZeroMomCorr->Fill(pizeroMom.Mag());
              hPiZeroMomFitted->Fill(pizeroMomFitted);

              hPiZeroMomResolution->Fill(slc_best_pzc_pizero_mom->at(slc_i) - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionCorr->Fill(pizeroMom.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionFitted->Fill(pizeroMomFitted - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
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
  nominalLine->DrawLine(134.9769, 0., 134.9769, 1.05 * hInvariantMassFitted->GetMaximum());

  TLegend *lInvariantMass = new TLegend(.6, .4, .85, .6);
  lInvariantMass->AddEntry(hInvariantMass, "Original", "l");
  lInvariantMass->AddEntry(hInvariantMassCorr, "Corrected", "l");
  lInvariantMass->AddEntry(hInvariantMassFitted, "Kinematic Fitting", "l");
  lInvariantMass->Draw();

  cInvariantMass->SaveAs(saveDir + "/invariant_mass.png");
  cInvariantMass->SaveAs(saveDir + "/invariant_mass.pdf");

  TCanvas *cPiZeroMom = new TCanvas("cPiZeroMom", "cPiZeroMom");
  cPiZeroMom->cd();

  NormaliseEntriesByBinWidth(hPiZeroMom);
  NormaliseEntriesByBinWidth(hPiZeroMomCorr);
  NormaliseEntriesByBinWidth(hPiZeroMomFitted);

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
}
