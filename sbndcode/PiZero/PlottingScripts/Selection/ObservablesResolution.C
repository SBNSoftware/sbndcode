void ObservablesResolution(const TString productionVersion)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/observables_resolution";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  TCanvas *cPiZeroMomResolution = new TCanvas("cPiZeroMomResolution", "cPiZeroMomResolution");
  cPiZeroMomResolution->cd();

  TH1F *hPiZeroMomResolution = new TH1F("hPiZeroMomResolution", ";p_{#pi^{0}} (#frac{Reco - True}{True});#pi^{0}", 40, -1, 1);
  rockboxEvents->Draw("(slc_best_pzc_pizero_mom-slc_true_pz_pizero_mom*1e3)/(slc_true_pz_pizero_mom*1e3)>>hPiZeroMomResolution",
                      "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hPiZeroMomResolution->SetLineColor(kMagenta+2);
  hPiZeroMomResolution->Draw();

  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_momentum_resolution.pdf");
  cPiZeroMomResolution->SaveAs(saveDir + "/pizero_momentum_resolution.png");

  TCanvas *cCosineThetaPiZeroResolution = new TCanvas("cCosineThetaPiZeroResolution", "cCosineThetaPiZeroResolution");
  cCosineThetaPiZeroResolution->cd();

  TH1F *hCosineThetaPiZeroResolution = new TH1F("hCosineThetaPiZeroResolution", ";cos(#theta_{#pi^{0}}) (#frac{Reco - True}{True});#pi^{0}", 40, -1, 1);
  rockboxEvents->Draw("(slc_best_pzc_cos_theta_pizero-slc_true_pz_cos_theta_pizero)/slc_true_pz_cos_theta_pizero>>hCosineThetaPiZeroResolution",
                      "slc_sel_incl && slc_true_event_type_incl==0 && slc_comp>.5");

  hCosineThetaPiZeroResolution->SetLineColor(kMagenta+2);
  hCosineThetaPiZeroResolution->Draw();

  cCosineThetaPiZeroResolution->SaveAs(saveDir + "/cosine_theta_pizero_resolution.pdf");
  cCosineThetaPiZeroResolution->SaveAs(saveDir + "/cosine_theta_pizero_resolution.png");

  TCanvas *cShowerEnergyResolution = new TCanvas("cShowerEnergyResolution", "cShowerEnergyResolution");
  cShowerEnergyResolution->cd();

  TH1F *hShowerEnergyResolution = new TH1F("hShowerEnergyResolution", ";E (#frac{Reco - True}{True});Showers", 40, -1, 1);
  rockboxEvents->Draw("(slc_pfp_shower_energy-slc_pfp_true_energy)/(slc_pfp_true_energy)>>hShowerEnergyResolution",
                      "slc_pfp_pdg==11 && slc_pfp_comp>.5 && slc_pfp_pur>.5");

  hShowerEnergyResolution->SetLineColor(kMagenta+2);
  hShowerEnergyResolution->Draw();

  cShowerEnergyResolution->SaveAs(saveDir + "/shower_energy_resolution.pdf");
  cShowerEnergyResolution->SaveAs(saveDir + "/shower_energy_resolution.png");
}
