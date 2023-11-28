#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

void PhotonRecoEff(const TString productionVersion)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/photon_reco_eff";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  const double gammaEnergyBins[11] = { 0., 50, 100, 150, 200, 250, 300, 400, 500, 750, 1e3 };
  const double pizeroMomBins[9]    = { 0., 60., 120., 180., 240., 300., 400., 600., 1e3 };
  const double cosThetaBins[10]    = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };
  const double cosComBins[11]      = { 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1. };

  TCanvas *cGamma0Energy = new TCanvas("cGamma0Energy", "cGamma0Energy");
  cGamma0Energy->cd();

  TH1F *hGamma0Energy = new TH1F("hGamma0Energy", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma0_energy*1e3>>hGamma0Energy", "nu_event_type_incl==0");
  NormaliseEntriesByBinWidth(hGamma0Energy);

  TH1F *hGamma0EnergyReco = new TH1F("hGamma0EnergyReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma0_energy*1e3>>hGamma0EnergyReco", "nu_event_type_incl==0 && nu_pz_gamma0_best_pfp_comp>.5 && nu_pz_gamma0_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hGamma0EnergyReco);

  TH1F *hGamma0EnergyGoodReco = new TH1F("hGamma0EnergyGoodReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma0_energy*1e3>>hGamma0EnergyGoodReco", "nu_event_type_incl==0 && nu_pz_gamma0_best_pfp_comp>.5 && nu_pz_gamma0_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hGamma0EnergyGoodReco);

  MakePlotEff(cGamma0Energy, hGamma0Energy, hGamma0EnergyReco, ";E_{#gamma_{0}} (MeV);Efficiency", kMagenta+2);

  cGamma0Energy->SaveAs(saveDir + "/gamma0_energy_reco_eff.png");
  cGamma0Energy->SaveAs(saveDir + "/gamma0_energy_reco_eff.pdf");

  TCanvas *cGamma1Energy = new TCanvas("cGamma1Energy", "cGamma1Energy");
  cGamma1Energy->cd();

  TH1F *hGamma1Energy = new TH1F("hGamma1Energy", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma1_energy*1e3>>hGamma1Energy", "nu_event_type_incl==0");
  NormaliseEntriesByBinWidth(hGamma1Energy);

  TH1F *hGamma1EnergyReco = new TH1F("hGamma1EnergyReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma1_energy*1e3>>hGamma1EnergyReco", "nu_event_type_incl==0 && nu_pz_gamma1_best_pfp_comp>.5 && nu_pz_gamma1_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hGamma1EnergyReco);

  MakePlotEff(cGamma1Energy, hGamma1Energy, hGamma1EnergyReco, ";E_{#gamma_{1}} (MeV);Efficiency", kMagenta+2);

  cGamma1Energy->SaveAs(saveDir + "/gamma1_energy_reco_eff.png");
  cGamma1Energy->SaveAs(saveDir + "/gamma1_energy_reco_eff.pdf");

  TCanvas *cPiZeroDecayAsym = new TCanvas("cPiZeroDecayAsym", "cPiZeroDecayAsym");
  cPiZeroDecayAsym->cd();

  TH1F *hPiZeroDecayAsym = new TH1F("hPiZeroDecayAsym", ";Decay Asymmetry;", 10, 0, 1);
  rockboxEvents->Draw("nu_pz_decay_asymmetry>>hPiZeroDecayAsym", "nu_event_type_incl==0");
  NormaliseEntriesByBinWidth(hPiZeroDecayAsym);

  TH1F *hPiZeroDecayAsymReco = new TH1F("hPiZeroDecayAsymReco", ";Decay Asymmetry;", 10, 0, 1);
  rockboxEvents->Draw("nu_pz_decay_asymmetry>>hPiZeroDecayAsymReco",
                      "nu_event_type_incl==0 && nu_pz_gamma0_best_pfp_comp>.5 && nu_pz_gamma0_best_pfp_pur>.5 && nu_pz_gamma1_best_pfp_comp>.5 && nu_pz_gamma1_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hPiZeroDecayAsymReco);

  MakePlotEff(cPiZeroDecayAsym, hPiZeroDecayAsym, hPiZeroDecayAsymReco, ";Decay Asymmetry;Efficiency", kMagenta+2, {.6, .77, .8, .92});

  cPiZeroDecayAsym->SaveAs(saveDir + "/decay_asymmetry_reco_eff.png");
  cPiZeroDecayAsym->SaveAs(saveDir + "/decay_asymmetry_reco_eff.pdf");

  TCanvas *cPiZeroMom = new TCanvas("cPiZeroMom", "cPiZeroMom");
  cPiZeroMom->cd();

  TH1F *hPiZeroMom = new TH1F("hPiZeroMom", ";p_{#pi^{0}} (MeV/c);", 8, pizeroMomBins);
  rockboxEvents->Draw("nu_pz_pizero_mom*1e3>>hPiZeroMom", "nu_event_type_incl==0");
  NormaliseEntriesByBinWidth(hPiZeroMom);

  TH1F *hPiZeroMomReco = new TH1F("hPiZeroMomReco", ";p_{#pi^{0}} (MeV/c);", 8, pizeroMomBins);
  rockboxEvents->Draw("nu_pz_pizero_mom*1e3>>hPiZeroMomReco",
                      "nu_event_type_incl==0 && nu_pz_gamma0_best_pfp_comp>.5 && nu_pz_gamma0_best_pfp_pur>.5 && nu_pz_gamma1_best_pfp_comp>.5 && nu_pz_gamma1_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hPiZeroMomReco);

  MakePlotEff(cPiZeroMom, hPiZeroMom, hPiZeroMomReco, ";p_{#pi^{0}} (MeV/c);Efficiency", kMagenta+2, {.6, .72, .8, .87});

  cPiZeroMom->SaveAs(saveDir + "/pizero_mom_reco_eff.png");
  cPiZeroMom->SaveAs(saveDir + "/pizero_mom_reco_eff.pdf");

  TCanvas *cPiZeroCosCom = new TCanvas("cPiZeroCosCom", "cPiZeroCosCom");
  cPiZeroCosCom->cd();

  TH1F *hPiZeroCosCom = new TH1F("hPiZeroCosCom", ";cos(#theta_{CoM} );", 10, cosComBins);
  rockboxEvents->Draw("nu_pz_cos_com>>hPiZeroCosCom", "nu_event_type_incl==0");
  NormaliseEntriesByBinWidth(hPiZeroCosCom);

  TH1F *hPiZeroCosComReco = new TH1F("hPiZeroCosComReco", ";cos(#theta_{CoM} );", 10, cosComBins);
  rockboxEvents->Draw("nu_pz_cos_com>>hPiZeroCosComReco",
                      "nu_event_type_incl==0 && nu_pz_gamma0_best_pfp_comp>.5 && nu_pz_gamma0_best_pfp_pur>.5 && nu_pz_gamma1_best_pfp_comp>.5 && nu_pz_gamma1_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hPiZeroCosComReco);

  MakePlotEff(cPiZeroCosCom, hPiZeroCosCom, hPiZeroCosComReco, ";cos(#theta_{CoM} );Efficiency", kMagenta+2, {.4, .27, .6, .42});

  cPiZeroCosCom->SaveAs(saveDir + "/cos_com_reco_eff.png");
  cPiZeroCosCom->SaveAs(saveDir + "/cos_com_reco_eff.pdf");

  TCanvas *cPiZeroCosTheta = new TCanvas("cPiZeroCosTheta", "cPiZeroCosTheta");
  cPiZeroCosTheta->cd();

  TH1F *hPiZeroCosTheta = new TH1F("hPiZeroCosTheta", ";cos(#theta_{#pi^{0}} );", 9, cosThetaBins);
  rockboxEvents->Draw("nu_pz_cos_theta_pizero>>hPiZeroCosTheta", "nu_event_type_incl==0");
  NormaliseEntriesByBinWidth(hPiZeroCosTheta);

  TH1F *hPiZeroCosThetaReco = new TH1F("hPiZeroCosThetaReco", ";cos(#theta_{#pi^{0}} );", 9, cosThetaBins);
  rockboxEvents->Draw("nu_pz_cos_theta_pizero>>hPiZeroCosThetaReco",
                      "nu_event_type_incl==0 && nu_pz_gamma0_best_pfp_comp>.5 && nu_pz_gamma0_best_pfp_pur>.5 && nu_pz_gamma1_best_pfp_comp>.5 && nu_pz_gamma1_best_pfp_pur>.5");
  NormaliseEntriesByBinWidth(hPiZeroCosThetaReco);

  MakePlotEff(cPiZeroCosTheta, hPiZeroCosTheta, hPiZeroCosThetaReco, ";cos(#theta_{#pi^{0}} );Efficiency", kMagenta+2, {.4, .72, .6, .87});

  cPiZeroCosTheta->SaveAs(saveDir + "/cos_theta_pizero_reco_eff.png");
  cPiZeroCosTheta->SaveAs(saveDir + "/cos_theta_pizero_reco_eff.pdf");
}
