#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"

#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"

template<class T>
void FillHists(T *t, TH1* hist, TH1* recoHist, const bool drop_bad_pdg = false);

void ComparisonEfficiency(const TString &file1, const TString &file2, const TString &save_name_extra,
			  const int &colour1 = kRed+2, const int &colour2 = kBlue+2,
			  std::array<float, 4> legend_position = {.6, .45, .85, .55})
{
  const TString save_dir = "/sbnd/data/users/hlay/crt/clustering/plots/comparison";
  gSystem->Exec("mkdir -p " + save_dir);
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree1 = new TChain("crtana/tree");
  TChain *tree2 = new TChain("crtana/tree");
  tree1->Add("/sbnd/data/users/hlay/crt/clustering/" + file1);
  tree2->Add("/sbnd/data/users/hlay/crt/clustering/" + file2);

  double bins[23];
  for(int i = 0; i < 16; ++i)
    bins[i] = i * .5;
  for(int i = 16; i < 19; ++i)
    bins[i] = 8 + (i - 16) * 1;
  bins[19] = 11; bins[20] = 15; bins[21] = 20; bins[22] = 25;

  TH1F *hTrueDepositEnergy1 = new TH1F("hTrueDepositEnergy1", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyReco1 = new TH1F("hTrueDepositEnergyReco1", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyNoDropped1 = new TH1F("hTrueDepositEnergyNoDropped1", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyNoDroppedReco1 = new TH1F("hTrueDepositEnergyRecoNoDropped1", ";True energy (MeV);Deposits", 22, bins);

  TH1F *hTrueDepositEnergy2 = new TH1F("hTrueDepositEnergy2", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyReco2 = new TH1F("hTrueDepositEnergyReco2", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyNoDropped2 = new TH1F("hTrueDepositEnergyNoDropped2", ";True energy (MeV);Deposits", 22, bins);
  TH1F *hTrueDepositEnergyNoDroppedReco2 = new TH1F("hTrueDepositEnergyRecoNoDropped2", ";True energy (MeV);Deposits", 22, bins);

  FillHists(tree1, hTrueDepositEnergy1, hTrueDepositEnergyReco1);
  FillHists(tree1, hTrueDepositEnergyNoDropped1, hTrueDepositEnergyNoDroppedReco1, true);
  FillHists(tree2, hTrueDepositEnergy2, hTrueDepositEnergyReco2);
  FillHists(tree2, hTrueDepositEnergyNoDropped2, hTrueDepositEnergyNoDroppedReco2, true);

  TCanvas *canvasEnergy = new TCanvas("canvasEnergy", "canvasEnergy");
  canvasEnergy->cd();

  MakeComparisonPlotEff(canvasEnergy, hTrueDepositEnergy1, hTrueDepositEnergyReco1,
			hTrueDepositEnergy2, hTrueDepositEnergyReco2,
			";True energy (MeV);Efficiency", "Improved Clustering", "Fixed Truth Matching",
			colour1, colour2, legend_position);

  if(save)
    {
      canvasEnergy->SaveAs(save_dir + "/true_energy_deposit_reco_eff_compare_" + save_name_extra + ".png");
      canvasEnergy->SaveAs(save_dir + "/true_energy_deposit_reco_eff_compare_" + save_name_extra + ".pdf");
    }

  TCanvas *canvasEnergyNoDropped = new TCanvas("canvasEnergyNoDropped", "canvasEnergyNoDropped");
  canvasEnergyNoDropped->cd();

  MakeComparisonPlotEff(canvasEnergyNoDropped, hTrueDepositEnergyNoDropped1, hTrueDepositEnergyNoDroppedReco1,
			hTrueDepositEnergyNoDropped2, hTrueDepositEnergyNoDroppedReco2,
			";True energy (MeV);Efficiency", "Improved Clustering", "Fixed Truth Matching",
			colour1, colour2, legend_position);
  if(save)
    {
      canvasEnergyNoDropped->SaveAs(save_dir + "/true_energy_deposit_no_dropped_reco_eff_compare_" + save_name_extra + ".png");
      canvasEnergyNoDropped->SaveAs(save_dir + "/true_energy_deposit_no_dropped_reco_eff_compare_" + save_name_extra + ".pdf");
    }
}

template<class T>
void FillHists(T *t, TH1* hist, TH1* recoHist, const bool drop_bad_pdg)
{
  std::vector<double> *td_energy = 0;
  std::vector<bool> *td_reco_status = 0;
  std::vector<int> *td_pdg = 0;

  t->SetBranchAddress("td_energy", &td_energy);
  t->SetBranchAddress("td_reco_status", &td_reco_status);
  t->SetBranchAddress("td_pdg", &td_pdg);

  const unsigned N = t->GetEntries();

  for(unsigned i = 0; i < N; ++i)
    {
      t->GetEntry(i);

      for(unsigned ii = 0; ii < td_energy->size(); ++ii)
        {
          if(drop_bad_pdg && td_pdg->at(ii) == -999999)
            continue;

          hist->Fill(1e3*td_energy->at(ii));

          if(td_reco_status->at(ii))
            recoHist->Fill(1e3*td_energy->at(ii));
	}
    }

  NormaliseEntriesByBinWidth(hist);
  NormaliseEntriesByBinWidth(recoHist);
}
