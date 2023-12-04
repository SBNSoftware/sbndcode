#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Particles.h"

void RecoEff(const TString productionVersion)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/reco_eff";
  gSystem->Exec("mkdir -p " + saveDir);

  const std::vector<int> colours             = { kMagenta + 2, kRed - 4 };
  const std::vector<TString> names           = { "Reco Eff", "Good Reco Eff" };
  const std::array<float, 4> legend_position = { .25, .86, .87, .91 };
  const int ncolumns                         = 3;
  const double energyBins[14]                = { 0., 50, 100, 150, 200, 250, 300, 400, 500, 750, 1e3, 1.25e3, 1.5e3, 2e3 };

  const TCut reco_cut      = "(reco_nTracks>0 || reco_nShowers>0) && ((reco_track_purity>.5 && reco_track_completeness>.5) || (reco_shower_purity>.5 && reco_shower_completeness>.5))";
  const TCut good_reco_cut = "(reco_nTracks>0 || reco_nShowers>0) && ((reco_track_purity>.8 && reco_track_completeness>.8) || (reco_shower_purity>.8 && reco_shower_completeness>.8))";

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox_recoeff.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *particleTree = new TChain("recoeff/ParticleTree");
  particleTree->Add(rockboxFile);

  for(auto const& particle : particles)
    {
      TCanvas *canvas = new TCanvas("canvas" + particle.name, "canvas" + particle.name);
      canvas->cd();

      const TCut base_cut = Form("abs(mc_PDG)==%d", particle.pdg);

      double bins[particle.energybins.size()];
      for(int i = 0; i < particle.energybins.size(); ++i)
        bins[i] = particle.energybins[i];

      TH1F *trueHist = new TH1F("trueHist" + particle.name, ";E (MeV);" + particle.latex_name, particle.energybins.size() - 1, bins);
      particleTree->Draw("mc_energy0 * 1e3>>trueHist" + particle.name, base_cut);
      NormaliseEntriesByBinWidth(trueHist);

      TH1F *recoHist = new TH1F("recoHist" + particle.name, ";E (MeV);" + particle.latex_name, particle.energybins.size() - 1, bins);
      particleTree->Draw("mc_energy0 * 1e3>>recoHist" + particle.name, base_cut + reco_cut);
      NormaliseEntriesByBinWidth(recoHist);

      TH1F *goodRecoHist = new TH1F("goodRecoHist" + particle.name, ";E (MeV);" + particle.latex_name, particle.energybins.size() - 1, bins);
      particleTree->Draw("mc_energy0 * 1e3>>goodRecoHist" + particle.name, base_cut + good_reco_cut);
      NormaliseEntriesByBinWidth(goodRecoHist);  

      MakePlotMultiEff(canvas, trueHist, { recoHist, goodRecoHist }, ";E (MeV);" + particle.latex_name, colours, names, legend_position, ncolumns);

      canvas->SaveAs(saveDir + "/" + particle.name + "_energy_reco_eff.png");
      canvas->SaveAs(saveDir + "/" + particle.name + "_energy_reco_eff.pdf");
    }
}
