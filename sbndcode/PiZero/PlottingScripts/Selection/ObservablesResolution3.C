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
  TH1F *hInvariantMassFitted = new TH1F("hInvariantMassFitted", ";M_{#gamma#gamma} (MeV/c^{2});#pi^{0}", 50, 0, 500);

  TH1F *hPiZeroMom = new TH1F("hPiZeroMom", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);
  TH1F *hPiZeroMomFitted = new TH1F("hPiZeroMomFitted", ";p_{#pi^{0}} (MeV/c);#pi^{0}", 8, pizeroMomBins);

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);
  TH1F *hPiZeroMomResolutionFitted = new TH1F("hPiZeroMomResolutionFitted", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -300, 400);

  TH1F *hCosThetaPiZero = new TH1F("hCosThetaPiZero", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroFitted = new TH1F("hCosThetaPiZeroFitted", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);

  TH1F *hCosThetaPiZeroResolution = new TH1F("hCosThetaPiZeroResolution", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionFitted = new TH1F("hCosThetaPiZeroResolutionFitted", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);

  const int N = ncpizeroEvents->GetEntries();

  int nsig = 0, nsigGood = 0;
  
  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      ncpizeroEvents->GetEntry(ev_i);
      
      for(int slc_i = 0; slc_i < slc_true_event_type_incl->size(); ++slc_i)
        {
          if(slc_true_event_type_incl->at(slc_i) == 0 && slc_comp->at(slc_i) > .5 && slc_sel_incl->at(slc_i))
            {
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

              hInvariantMass->Fill(invariantMass);
              hInvariantMassFitted->Fill(invariantMassFitted);

              hPiZeroMom->Fill(pizeroMom.Mag());
              hPiZeroMomFitted->Fill(pizeroMomFitted);

              hPiZeroMomResolution->Fill(pizeroMom.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionFitted->Fill(pizeroMomFitted - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));

              hCosThetaPiZero->Fill(cos(pizeroMom.Angle(zaxis)));
              //              hCosThetaPiZeroFitted->Fill();

              hCosThetaPiZeroResolution->Fill(cos(pizeroMom.Angle(zaxis)) - slc_true_pz_cos_theta_pizero->at(slc_i).at(0));
              //              hCosThetaPiZeroResolutionFitted->Fill();
            }
        }
    }

  std::cout << "Fitting good: " << nsigGood << " / " << nsig << " (" << (100. * nsigGood) / nsig << "%)" << std::endl;

  TCanvas *cInvariantMass = new TCanvas("cInvariantMass", "cInvariantMass");
  cInvariantMass->cd();

  hInvariantMass->SetLineColor(kMagenta+2);
  hInvariantMassFitted->SetLineColor(kCyan+2);

  hInvariantMassFitted->Draw();
  hInvariantMass->Draw("same");

  TLine *nominalLine = new TLine();
  nominalLine->SetLineColor(kRed+2);
  nominalLine->SetLineWidth(5);
  nominalLine->DrawLine(134.9769, 0., 134.9769, 1.05 * hInvariantMass->GetMaximum());

  TLegend *lInvariantMass = new TLegend(.6, .4, .85, .6);
  lInvariantMass->AddEntry(hInvariantMass, "Corrections", "l");
  lInvariantMass->AddEntry(hInvariantMassFitted, "Kinematic Fitting", "l");
  lInvariantMass->Draw();

  cInvariantMass->SaveAs(saveDir + "/invariant_mass.png");
  cInvariantMass->SaveAs(saveDir + "/invariant_mass.pdf");

  TCanvas *cPiZeroMom = new TCanvas("cPiZeroMom", "cPiZeroMom");
  cPiZeroMom->cd();

  NormaliseEntriesByBinWidth(hPiZeroMom);
  NormaliseEntriesByBinWidth(hPiZeroMomFitted);

  hPiZeroMom->SetLineColor(kMagenta+2);
  hPiZeroMomFitted->SetLineColor(kCyan+2);

  hPiZeroMomFitted->Draw("hist");
  hPiZeroMom->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.png");
  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.pdf");

  TCanvas *cPiZeroMomResolution = new TCanvas("cPiZeroMomResolution", "cPiZeroMomResolution");
  cPiZeroMomResolution->cd();

  hPiZeroMomResolution->SetLineColor(kMagenta+2);
  hPiZeroMomResolutionFitted->SetLineColor(kCyan+2);

  std::cout << "\nStandard === Mean: " << hPiZeroMomResolution->GetMean() << " SD: " << hPiZeroMomResolution->GetStdDev() << std::endl;
  std::cout << "Fitted   === Mean: " << hPiZeroMomResolutionFitted->GetMean() << " SD: " << hPiZeroMomResolutionFitted->GetStdDev() << '\n' << std::endl;

  hPiZeroMomResolutionFitted->Draw("hist");
  hPiZeroMomResolution->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_mom_resolution.png");
  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_mom_resolution.pdf");

  TCanvas *cCosThetaPiZero = new TCanvas("cCosThetaPiZero", "cCosThetaPiZero");
  cCosThetaPiZero->cd();

  NormaliseEntriesByBinWidth(hCosThetaPiZero);
  NormaliseEntriesByBinWidth(hCosThetaPiZeroFitted);

  hCosThetaPiZero->SetLineColor(kMagenta+2);
  hCosThetaPiZeroFitted->SetLineColor(kCyan+2);

  hCosThetaPiZeroFitted->Draw("hist");
  hCosThetaPiZero->Draw("histsame");

  lInvariantMass->Draw();

  cCosThetaPiZero->SaveAs(saveDir + "/cos_theta_pizero.png");
  cCosThetaPiZero->SaveAs(saveDir + "/cos_theta_pizero.pdf");

  TCanvas *cCosThetaPiZeroResolution = new TCanvas("cCosThetaPiZeroResolution", "cCosThetaPiZeroResolution");
  cCosThetaPiZeroResolution->cd();

  hCosThetaPiZeroResolution->SetLineColor(kMagenta+2);
  hCosThetaPiZeroResolutionFitted->SetLineColor(kCyan+2);

  hCosThetaPiZeroResolutionFitted->Draw("hist");
  hCosThetaPiZeroResolution->Draw("histsame");

  lInvariantMass->Draw();

  cCosThetaPiZeroResolution->SaveAs(saveDir + "/cos_theta_pizero_resolution.png");
  cCosThetaPiZeroResolution->SaveAs(saveDir + "/cos_theta_pizero_resolution.pdf");
}
