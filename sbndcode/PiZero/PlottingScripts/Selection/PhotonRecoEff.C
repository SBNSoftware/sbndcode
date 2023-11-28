#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

void PhotonRecoEff(const TString productionVersion)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/photon_reco_eff";
  gSystem->Exec("mkdir -p " + saveDir);

  const TCut base_cut = "nu_event_type_incl == 0";

  const TCut gamma0_reco_cut = base_cut + "nu_pz_gamma0_best_pfp_comp > .5 && nu_pz_gamma0_best_pfp_pur > .5";
  const TCut gamma1_reco_cut = base_cut + "nu_pz_gamma1_best_pfp_comp > .5 && nu_pz_gamma1_best_pfp_pur > .5";
  const TCut pizero_reco_cut = gamma0_reco_cut + gamma1_reco_cut;

  const TCut gamma0_good_reco_cut = base_cut + "nu_pz_gamma0_best_pfp_comp > .8 && nu_pz_gamma0_best_pfp_pur > .8";
  const TCut gamma1_good_reco_cut = base_cut + "nu_pz_gamma1_best_pfp_comp > .8 && nu_pz_gamma1_best_pfp_pur > .8";
  const TCut pizero_good_reco_cut = gamma0_good_reco_cut + gamma1_good_reco_cut;

  const std::vector<int> colours             = { kMagenta + 2, kRed - 4 };
  const std::vector<TString> names           = { "Reco Eff", "Good Reco Eff" };
  const std::array<float, 4> legend_position = { .25, .86, .87, .91 };
  const int ncolumns                         = 3;

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
  rockboxEvents->Draw("nu_pz_gamma0_energy*1e3>>hGamma0Energy", base_cut);
  NormaliseEntriesByBinWidth(hGamma0Energy);

  TH1F *hGamma0EnergyReco = new TH1F("hGamma0EnergyReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma0_energy*1e3>>hGamma0EnergyReco", gamma0_reco_cut);
  NormaliseEntriesByBinWidth(hGamma0EnergyReco);

  TH1F *hGamma0EnergyGoodReco = new TH1F("hGamma0EnergyGoodReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma0_energy*1e3>>hGamma0EnergyGoodReco", gamma0_good_reco_cut);
  NormaliseEntriesByBinWidth(hGamma0EnergyGoodReco);

  MakePlotMultiEff(cGamma0Energy, hGamma0Energy, { hGamma0EnergyReco, hGamma0EnergyGoodReco }, ";E_{#gamma_{0}} (MeV);Efficiency", colours, names, legend_position, ncolumns);

  cGamma0Energy->SaveAs(saveDir + "/gamma0_energy_reco_eff.png");
  cGamma0Energy->SaveAs(saveDir + "/gamma0_energy_reco_eff.pdf");

  TCanvas *cGamma1Energy = new TCanvas("cGamma1Energy", "cGamma1Energy");
  cGamma1Energy->cd();

  TH1F *hGamma1Energy = new TH1F("hGamma1Energy", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma1_energy*1e3>>hGamma1Energy", base_cut);
  NormaliseEntriesByBinWidth(hGamma1Energy);

  TH1F *hGamma1EnergyReco = new TH1F("hGamma1EnergyReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma1_energy*1e3>>hGamma1EnergyReco", gamma1_reco_cut);
  NormaliseEntriesByBinWidth(hGamma1EnergyReco);

  TH1F *hGamma1EnergyGoodReco = new TH1F("hGamma1EnergyGoodReco", ";E_{#gamma} (MeV);", 10, gammaEnergyBins);
  rockboxEvents->Draw("nu_pz_gamma1_energy*1e3>>hGamma1EnergyGoodReco", gamma1_good_reco_cut);
  NormaliseEntriesByBinWidth(hGamma1EnergyGoodReco);

  MakePlotMultiEff(cGamma1Energy, hGamma1Energy, { hGamma1EnergyReco, hGamma1EnergyGoodReco}, ";E_{#gamma_{1}} (MeV);Efficiency", colours, names, legend_position, ncolumns);

  cGamma1Energy->SaveAs(saveDir + "/gamma1_energy_reco_eff.png");
  cGamma1Energy->SaveAs(saveDir + "/gamma1_energy_reco_eff.pdf");

  TCanvas *cPiZeroDecayAsym = new TCanvas("cPiZeroDecayAsym", "cPiZeroDecayAsym");
  cPiZeroDecayAsym->cd();

  TH1F *hPiZeroDecayAsym = new TH1F("hPiZeroDecayAsym", ";Decay Asymmetry;", 10, 0, 1);
  rockboxEvents->Draw("nu_pz_decay_asymmetry>>hPiZeroDecayAsym", base_cut);
  NormaliseEntriesByBinWidth(hPiZeroDecayAsym);

  TH1F *hPiZeroDecayAsymReco = new TH1F("hPiZeroDecayAsymReco", ";Decay Asymmetry;", 10, 0, 1);
  rockboxEvents->Draw("nu_pz_decay_asymmetry>>hPiZeroDecayAsymReco",
                      pizero_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroDecayAsymReco);

  TH1F *hPiZeroDecayAsymGoodReco = new TH1F("hPiZeroDecayAsymGoodReco", ";Decay Asymmetry;", 10, 0, 1);
  rockboxEvents->Draw("nu_pz_decay_asymmetry>>hPiZeroDecayAsymGoodReco",
                      pizero_good_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroDecayAsymGoodReco);

  MakePlotMultiEff(cPiZeroDecayAsym, hPiZeroDecayAsym, { hPiZeroDecayAsymReco, hPiZeroDecayAsymGoodReco}, ";Decay Asymmetry;Efficiency", colours, names, legend_position, ncolumns);

  cPiZeroDecayAsym->SaveAs(saveDir + "/decay_asymmetry_reco_eff.png");
  cPiZeroDecayAsym->SaveAs(saveDir + "/decay_asymmetry_reco_eff.pdf");

  TCanvas *cPiZeroMom = new TCanvas("cPiZeroMom", "cPiZeroMom");
  cPiZeroMom->cd();

  TH1F *hPiZeroMom = new TH1F("hPiZeroMom", ";p_{#pi^{0}} (MeV/c);", 8, pizeroMomBins);
  rockboxEvents->Draw("nu_pz_pizero_mom*1e3>>hPiZeroMom", base_cut);
  NormaliseEntriesByBinWidth(hPiZeroMom);

  TH1F *hPiZeroMomReco = new TH1F("hPiZeroMomReco", ";p_{#pi^{0}} (MeV/c);", 8, pizeroMomBins);
  rockboxEvents->Draw("nu_pz_pizero_mom*1e3>>hPiZeroMomReco",
                      pizero_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroMomReco);

  TH1F *hPiZeroMomGoodReco = new TH1F("hPiZeroMomGoodReco", ";p_{#pi^{0}} (MeV/c);", 8, pizeroMomBins);
  rockboxEvents->Draw("nu_pz_pizero_mom*1e3>>hPiZeroMomGoodReco",
                      pizero_good_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroMomGoodReco);

  MakePlotMultiEff(cPiZeroMom, hPiZeroMom, { hPiZeroMomReco, hPiZeroMomGoodReco}, ";p_{#pi^{0}} (MeV/c);Efficiency", colours, names, legend_position, ncolumns);

  cPiZeroMom->SaveAs(saveDir + "/pizero_mom_reco_eff.png");
  cPiZeroMom->SaveAs(saveDir + "/pizero_mom_reco_eff.pdf");

  TCanvas *cPiZeroCosCom = new TCanvas("cPiZeroCosCom", "cPiZeroCosCom");
  cPiZeroCosCom->cd();

  TH1F *hPiZeroCosCom = new TH1F("hPiZeroCosCom", ";cos(#theta_{CoM} );", 10, cosComBins);
  rockboxEvents->Draw("nu_pz_cos_com>>hPiZeroCosCom", base_cut);
  NormaliseEntriesByBinWidth(hPiZeroCosCom);

  TH1F *hPiZeroCosComReco = new TH1F("hPiZeroCosComReco", ";cos(#theta_{CoM} );", 10, cosComBins);
  rockboxEvents->Draw("nu_pz_cos_com>>hPiZeroCosComReco",
                      pizero_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroCosComReco);

  TH1F *hPiZeroCosComGoodReco = new TH1F("hPiZeroCosComGoodReco", ";cos(#theta_{CoM} );", 10, cosComBins);
  rockboxEvents->Draw("nu_pz_cos_com>>hPiZeroCosComGoodReco",
                      pizero_good_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroCosComGoodReco);

  MakePlotMultiEff(cPiZeroCosCom, hPiZeroCosCom, { hPiZeroCosComReco, hPiZeroCosComGoodReco }, ";cos(#theta_{CoM} );Efficiency", colours, names, legend_position, ncolumns);

  cPiZeroCosCom->SaveAs(saveDir + "/cos_com_reco_eff.png");
  cPiZeroCosCom->SaveAs(saveDir + "/cos_com_reco_eff.pdf");

  TCanvas *cPiZeroCosTheta = new TCanvas("cPiZeroCosTheta", "cPiZeroCosTheta");
  cPiZeroCosTheta->cd();

  TH1F *hPiZeroCosTheta = new TH1F("hPiZeroCosTheta", ";cos(#theta_{#pi^{0}} );", 9, cosThetaBins);
  rockboxEvents->Draw("nu_pz_cos_theta_pizero>>hPiZeroCosTheta", base_cut);
  NormaliseEntriesByBinWidth(hPiZeroCosTheta);

  TH1F *hPiZeroCosThetaReco = new TH1F("hPiZeroCosThetaReco", ";cos(#theta_{#pi^{0}} );", 9, cosThetaBins);
  rockboxEvents->Draw("nu_pz_cos_theta_pizero>>hPiZeroCosThetaReco",
                      pizero_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroCosThetaReco);

  TH1F *hPiZeroCosThetaGoodReco = new TH1F("hPiZeroCosThetaGoodReco", ";cos(#theta_{#pi^{0}} );", 9, cosThetaBins);
  rockboxEvents->Draw("nu_pz_cos_theta_pizero>>hPiZeroCosThetaGoodReco",
                      pizero_good_reco_cut);
  NormaliseEntriesByBinWidth(hPiZeroCosThetaGoodReco);

  MakePlotMultiEff(cPiZeroCosTheta, hPiZeroCosTheta, { hPiZeroCosThetaReco, hPiZeroCosThetaGoodReco }, ";cos(#theta_{#pi^{0}} );Efficiency", colours, names, legend_position, ncolumns);

  cPiZeroCosTheta->SaveAs(saveDir + "/cos_theta_pizero_reco_eff.png");
  cPiZeroCosTheta->SaveAs(saveDir + "/cos_theta_pizero_reco_eff.pdf");
}
