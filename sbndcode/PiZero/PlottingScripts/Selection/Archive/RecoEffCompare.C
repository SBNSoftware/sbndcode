#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Particles.h"

void RecoEffCompare(const TString productionVersionA = "NCPiZeroAv10", const TString productionVersionB = "NCPiZeroAv8")
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersionA + "/reco_eff_compare";
  gSystem->Exec("mkdir -p " + saveDir);

  const std::vector<int> colours             = { kMagenta + 2, kMagenta - 7, kRed + 2, kRed - 7 };
  const std::vector<TString> names           = { "Fixed", "Fixed - HQ", "Unfixed", "Unfixed - HQ" };
  const std::array<float, 4> legend_position = { .25, .82, .87, .93 };
  const int ncolumns                         = 3;
  const double energyBins[14]                = { 0., 50, 100, 150, 200, 250, 300, 400, 500, 750, 1e3, 1.25e3, 1.5e3, 2e3 };

  const TCut reco_cut      = "(reco_nTracks>0 || reco_nShowers>0) && ((reco_track_purity>.5 && reco_track_completeness>.5) || (reco_shower_purity>.5 && reco_shower_completeness>.5))";
  const TCut good_reco_cut = "(reco_nTracks>0 || reco_nShowers>0) && ((reco_track_purity>.8 && reco_track_completeness>.8) || (reco_shower_purity>.8 && reco_shower_completeness>.8))";

  const TString rockboxFileA = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersionA + "/" + productionVersionA + "_rockbox_recoeff.root";
  const TString rockboxFileB = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersionB + "/" + productionVersionB + "_rockbox_recoeff.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *particleTreeA = new TChain("recoeff/ParticleTree");
  particleTreeA->Add(rockboxFileA);

  TChain *particleTreeB = new TChain("recoeff/ParticleTree");
  particleTreeB->Add(rockboxFileB);

  for(auto const& particle : particles)
    {
      TCanvas *canvas = new TCanvas("canvas" + particle.name, "canvas" + particle.name);
      canvas->cd();

      const TCut base_cut = Form("abs(mc_PDG)==%d", particle.pdg);

      double bins[particle.energybins.size()];
      for(int i = 0; i < particle.energybins.size(); ++i)
        bins[i] = particle.energybins[i];

      TH1F *trueHist = new TH1F("trueHist" + particle.name, ";E (MeV);" + particle.latex_name + " Efficiency", particle.energybins.size() - 1, bins);
      particleTreeA->Draw("mc_energy0 * 1e3>>trueHist" + particle.name, base_cut);
      NormaliseEntriesByBinWidth(trueHist);

      TH1F *recoHistA = new TH1F("recoHistA" + particle.name, ";E (MeV);" + particle.latex_name + " Efficiency", particle.energybins.size() - 1, bins);
      particleTreeA->Draw("mc_energy0 * 1e3>>recoHistA" + particle.name, base_cut + reco_cut);
      NormaliseEntriesByBinWidth(recoHistA);

      TH1F *recoHistB = new TH1F("recoHistB" + particle.name, ";E (MeV);" + particle.latex_name + " Efficiency", particle.energybins.size() - 1, bins);
      particleTreeB->Draw("mc_energy0 * 1e3>>recoHistB" + particle.name, base_cut + reco_cut);
      NormaliseEntriesByBinWidth(recoHistB);

      TH1F *goodRecoHistA = new TH1F("goodRecoHistA" + particle.name, ";E (MeV);" + particle.latex_name + " Efficiency", particle.energybins.size() - 1, bins);
      particleTreeA->Draw("mc_energy0 * 1e3>>goodRecoHistA" + particle.name, base_cut + good_reco_cut);
      NormaliseEntriesByBinWidth(goodRecoHistA);

      TH1F *goodRecoHistB = new TH1F("goodRecoHistB" + particle.name, ";E (MeV);" + particle.latex_name + " Efficiency", particle.energybins.size() - 1, bins);
      particleTreeB->Draw("mc_energy0 * 1e3>>goodRecoHistB" + particle.name, base_cut + good_reco_cut);
      NormaliseEntriesByBinWidth(goodRecoHistB);

      MakePlotMultiEff(canvas, trueHist, { recoHistA, goodRecoHistA, recoHistB, goodRecoHistB }, ";E (MeV);" + particle.latex_name + " Efficiency", colours, names, legend_position, ncolumns);

      canvas->SaveAs(saveDir + "/" + particle.name + "_energy_reco_eff.png");
      canvas->SaveAs(saveDir + "/" + particle.name + "_energy_reco_eff.pdf");
    }
}
