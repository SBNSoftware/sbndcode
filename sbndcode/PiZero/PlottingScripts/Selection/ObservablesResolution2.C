#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "ObservablesCommon.C"

void ObservablesResolution2(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution_2";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gStyle->SetTitleOffset(1.3, "y");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(rockboxFile);
  events->Add(ncpizeroFile);

  InitialiseTree(events);

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

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionTrackDir = new TH1F("hPiZeroMomResolutionTrackDir", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionEnergyCorr = new TH1F("hPiZeroMomResolutionEnergyCorr", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);
  TH1F *hPiZeroMomResolutionBoth = new TH1F("hPiZeroMomResolutionBoth", ";p_{#pi^{0}} (Reco - True) (MeV/c);#pi^{0}", 35, -310, 390);

  TH1F *hCosThetaPiZero = new TH1F("hCosThetaPiZero", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroTrackDir = new TH1F("hCosThetaPiZeroTrackDir", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroEnergyCorr = new TH1F("hCosThetaPiZeroEnergyCorr", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);
  TH1F *hCosThetaPiZeroBoth = new TH1F("hCosThetaPiZeroBoth", ";cos(#theta_{#pi^{0}});#pi^{0}", 9, cosThetaBins);

  TH1F *hCosThetaPiZeroResolution = new TH1F("hCosThetaPiZeroResolution", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionTrackDir = new TH1F("hCosThetaPiZeroResolutionTrackDir", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionEnergyCorr = new TH1F("hCosThetaPiZeroResolutionEnergyCorr", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);
  TH1F *hCosThetaPiZeroResolutionBoth = new TH1F("hCosThetaPiZeroResolutionBoth", ";cos(#theta_{#pi^{0}}) (Reco - True);#pi^{0}", 50, -2, 3);

  const int N = events->GetEntries();
  
  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      events->GetEntry(ev_i);
      
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

              const TVector3 dir0 = shwEn0 > 150 ? shwDir0 : trkDir0;
              const TVector3 dir1 = shwEn1 > 150 ? shwDir1 : trkDir1;

              const double cosineThetaGammaGammaShw = shwDir0.Dot(shwDir1) / (shwDir0.Mag() * shwDir1.Mag());
              const double cosineThetaGammaGammaCor = dir0.Dot(dir1) / (dir0.Mag() * dir1.Mag());

              const double invariantMass           = sqrt(2 * shwEn0 * shwEn1 * (1 - cosineThetaGammaGammaShw));
              const double invariantMassTrackDir   = sqrt(2 * shwEn0 * shwEn1 * (1 - cosineThetaGammaGammaCor));
              const double invariantMassEnergyCorr = sqrt(2 * corrEn0 * corrEn1 * (1 - cosineThetaGammaGammaShw));
              const double invariantMassBoth       = sqrt(2 * corrEn0 * corrEn1 * (1 - cosineThetaGammaGammaCor));

              const TVector3 pizeroMom           = shwEn0 * shwDir0 + shwEn1 * shwDir1;
              const TVector3 pizeroMomTrackDir   = shwEn0 * dir0 + shwEn1 * dir1;
              const TVector3 pizeroMomEnergyCorr = corrEn0 * shwDir0 + corrEn1 * shwDir1;
              const TVector3 pizeroMomBoth       = corrEn0 * dir0 + corrEn1 * dir1;

              const TVector3 zaxis = TVector3(0., 0., 1.);

              hInvariantMass->Fill(slc_best_pzc_invariant_mass->at(slc_i));
              hInvariantMassTrackDir->Fill(invariantMassTrackDir);
              hInvariantMassEnergyCorr->Fill(invariantMassEnergyCorr);
              hInvariantMassBoth->Fill(invariantMassBoth);

              hPiZeroMom->Fill(pizeroMom.Mag());
              hPiZeroMomTrackDir->Fill(pizeroMomTrackDir.Mag());
              hPiZeroMomEnergyCorr->Fill(pizeroMomEnergyCorr.Mag());
              hPiZeroMomBoth->Fill(pizeroMomBoth.Mag());

              hPiZeroMomResolution->Fill(pizeroMom.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionTrackDir->Fill(pizeroMomTrackDir.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionEnergyCorr->Fill(pizeroMomEnergyCorr.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));
              hPiZeroMomResolutionBoth->Fill(pizeroMomBoth.Mag() - 1e3 * slc_true_pz_pizero_mom->at(slc_i).at(0));

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
  nominalLine->DrawLine(kPiZeroMass, 0., kPiZeroMass, 1.05 * hInvariantMass->GetMaximum());

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

  hPiZeroMom->SetLineColor(kMagenta+2);
  hPiZeroMomTrackDir->SetLineColor(kCyan+2);
  hPiZeroMomEnergyCorr->SetLineColor(kGreen+2);
  hPiZeroMomBoth->SetLineColor(kOrange+2);

  hPiZeroMomTrackDir->Draw("hist");
  hPiZeroMom->Draw("histsame");
  hPiZeroMomEnergyCorr->Draw("histsame");
  hPiZeroMomBoth->Draw("histsame");

  lInvariantMass->Draw();

  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.png");
  cPiZeroMom->SaveAs(saveDir + "/pizero_mom.pdf");

  TCanvas *cPiZeroMomResolution = new TCanvas("cPiZeroMomResolution", "cPiZeroMomResolution");
  cPiZeroMomResolution->cd();

  hPiZeroMomResolution->SetLineColor(kMagenta+2);
  hPiZeroMomResolutionTrackDir->SetLineColor(kCyan+2);
  hPiZeroMomResolutionEnergyCorr->SetLineColor(kGreen+2);
  hPiZeroMomResolutionBoth->SetLineColor(kOrange+2);

  hPiZeroMomResolutionTrackDir->Draw("hist");
  hPiZeroMomResolution->Draw("histsame");
  hPiZeroMomResolutionEnergyCorr->Draw("histsame");
  hPiZeroMomResolutionBoth->Draw("histsame");

  lInvariantMass->Draw();

  TPaveText *tPiZeroMomResolution = new TPaveText(.56, .65, .68, .75, "NDC");
  tPiZeroMomResolution->AddText("Standard");
  tPiZeroMomResolution->AddText("Track Direction");
  tPiZeroMomResolution->AddText("Energy Correction");
  tPiZeroMomResolution->AddText("Both");
  tPiZeroMomResolution->SetTextAlign(12);
  tPiZeroMomResolution->SetTextSize(0.02);
  tPiZeroMomResolution->SetBorderSize(0);
  tPiZeroMomResolution->SetFillColor(kWhite);
  tPiZeroMomResolution->Draw();

  TPaveText *tPiZeroMomResolution2 = new TPaveText(.68, .65, .85, .75, "NDC");
  tPiZeroMomResolution2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolution->GetMean(), hPiZeroMomResolution->GetStdDev()));
  tPiZeroMomResolution2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolutionTrackDir->GetMean(), hPiZeroMomResolutionTrackDir->GetStdDev()));
  tPiZeroMomResolution2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolutionEnergyCorr->GetMean(), hPiZeroMomResolutionEnergyCorr->GetStdDev()));
  tPiZeroMomResolution2->AddText(Form("Mean: %.2f MeV/c #sigma: %.2f MeV/c", hPiZeroMomResolutionBoth->GetMean(), hPiZeroMomResolutionBoth->GetStdDev()));
  tPiZeroMomResolution2->SetTextAlign(12);
  tPiZeroMomResolution2->SetTextSize(0.02);
  tPiZeroMomResolution2->SetBorderSize(0);
  tPiZeroMomResolution2->SetFillColor(kWhite);
  tPiZeroMomResolution2->Draw();

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

  TPaveText *tCosThetaPiZeroResolution = new TPaveText(.56, .65, .68, .75, "NDC");
  tCosThetaPiZeroResolution->AddText("Standard");
  tCosThetaPiZeroResolution->AddText("Track Direction");
  tCosThetaPiZeroResolution->AddText("Energy Correction");
  tCosThetaPiZeroResolution->AddText("Both");
  tCosThetaPiZeroResolution->SetTextAlign(12);
  tCosThetaPiZeroResolution->SetTextSize(0.02);
  tCosThetaPiZeroResolution->SetBorderSize(0);
  tCosThetaPiZeroResolution->SetFillColor(kWhite);
  tCosThetaPiZeroResolution->Draw();

  TPaveText *tCosThetaPiZeroResolution2 = new TPaveText(.68, .65, .85, .75, "NDC");
  tCosThetaPiZeroResolution2->AddText(Form("Mean: %.4f", hCosThetaPiZeroResolution->GetMean()));
  tCosThetaPiZeroResolution2->AddText(Form("Mean: %.4f", hCosThetaPiZeroResolutionTrackDir->GetMean()));
  tCosThetaPiZeroResolution2->AddText(Form("Mean: %.4f", hCosThetaPiZeroResolutionEnergyCorr->GetMean()));
  tCosThetaPiZeroResolution2->AddText(Form("Mean: %.4f", hCosThetaPiZeroResolutionBoth->GetMean()));
  tCosThetaPiZeroResolution2->SetTextAlign(12);
  tCosThetaPiZeroResolution2->SetTextSize(0.02);
  tCosThetaPiZeroResolution2->SetBorderSize(0);
  tCosThetaPiZeroResolution2->SetFillColor(kWhite);
  tCosThetaPiZeroResolution2->Draw();

  TPaveText *tCosThetaPiZeroResolution3 = new TPaveText(.78, .65, .85, .75, "NDC");
  tCosThetaPiZeroResolution3->AddText(Form("#sigma: %.4f", hCosThetaPiZeroResolution->GetStdDev()));
  tCosThetaPiZeroResolution3->AddText(Form("#sigma: %.4f", hCosThetaPiZeroResolutionTrackDir->GetStdDev()));
  tCosThetaPiZeroResolution3->AddText(Form("#sigma: %.4f", hCosThetaPiZeroResolutionEnergyCorr->GetStdDev()));
  tCosThetaPiZeroResolution3->AddText(Form("#sigma: %.4f", hCosThetaPiZeroResolutionBoth->GetStdDev()));
  tCosThetaPiZeroResolution3->SetTextAlign(12);
  tCosThetaPiZeroResolution3->SetTextSize(0.02);
  tCosThetaPiZeroResolution3->SetBorderSize(0);
  tCosThetaPiZeroResolution3->SetFillColor(kWhite);
  tCosThetaPiZeroResolution3->Draw();

  cCosThetaPiZeroResolution->SaveAs(saveDir + "/cos_theta_pizero_resolution.png");
  cCosThetaPiZeroResolution->SaveAs(saveDir + "/cos_theta_pizero_resolution.pdf");
}
