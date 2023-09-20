#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void ConfusionMatrixChi2(const bool efficiency_mode = true, const bool purity_mode = false,
			 const bool raw_mode = false, const bool require_primary = false,
			 const double comp_thresh = -1, const double pur_thresh = -1)
{
  if(efficiency_mode && purity_mode)
    {
      std::cout << "Cannot run both efficiency and purity mode" << std::endl;
      return;
    }

  const TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/confusion_matrices_chi2_pure_stopping";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/bnb/razzled*.root");
  tree->Add("/sbnd/data/users/hlay/razzled/intrnue/razzled*.root");

  const std::map<int, int> razzledMap = { { 13, 0 }, { 211, 1 }, { 321, 2 }, { 2212, 3 }, { 0, 4 } };
  std::vector<TString> axisLabels  = { "", "#mu^{#pm}", "#pi^{#pm}","K^{#pm}", "p", "other" };

  int truePDG;
  std::string *trueEndProcess = 0;
  float trueEndMomentum, energyComp, energyPurity, trackStartX, trackStartY, trackStartZ;
  bool recoPrimary, unambiguousSlice, trackContained;

  float pfp_trackScore, trk_length, trk_chi2PIDMuon, trk_chi2PIDProton, trackChi2PIDPion, trackChi2PIDKaon;
  int chi2PDG;

  tree->SetBranchAddress("pfp_trackScore", &pfp_trackScore);
  tree->SetBranchAddress("trk_length", &trk_length);
  tree->SetBranchAddress("trk_chi2PIDMuon", &trk_chi2PIDMuon);
  tree->SetBranchAddress("trk_chi2PIDProton", &trk_chi2PIDProton);
  tree->SetBranchAddress("trackChi2PIDPion", &trackChi2PIDPion);
  tree->SetBranchAddress("trackChi2PIDKaon", &trackChi2PIDKaon);
  tree->SetBranchAddress("chi2PDG", &chi2PDG);

  tree->SetBranchAddress("truePDG", &truePDG);
  tree->SetBranchAddress("trueEndProcess", &trueEndProcess);
  tree->SetBranchAddress("trueEndMomentum", &trueEndMomentum);
  tree->SetBranchAddress("energyComp", &energyComp);
  tree->SetBranchAddress("energyPurity", &energyPurity);
  tree->SetBranchAddress("recoPrimary", &recoPrimary);
  tree->SetBranchAddress("unambiguousSlice", &unambiguousSlice);

  tree->SetBranchAddress("trackStartX", &trackStartX);
  tree->SetBranchAddress("trackStartY", &trackStartY);
  tree->SetBranchAddress("trackStartZ", &trackStartZ);
  tree->SetBranchAddress("trackContained", &trackContained);

  const unsigned N = tree->GetEntries();

  const unsigned n_categories = 5;
  TH2F *hConfusionMatrix = new TH2F("hConfusionMatrix", ";True Class;Chi2 Class;",
				    n_categories, 0, n_categories, n_categories, 0, n_categories);
  
  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);
      
      if(unambiguousSlice)
	continue;
      
      if(abs(truePDG) != 13 && abs(truePDG) != 211 && abs(truePDG) != 321 && abs(truePDG) != 2212)
	truePDG = 0;

      if(abs(truePDG) == 211 && trueEndProcess->find("Inel") != std::string::npos)
        continue;

      if(abs(truePDG) == 321 && trueEndProcess->find("Inel") != std::string::npos)
        continue;

      if(abs(truePDG) == 13 && (trueEndProcess->find("Transp") != std::string::npos
                                || (trueEndProcess->find("Decay") != std::string::npos && trueEndMomentum > 0)))
        continue;

      if(abs(truePDG) == 2212 && trueEndProcess->find("Inel") != std::string::npos)
        continue;
      if(chi2PDG == -1)
	chi2PDG = 0;

      if(require_primary && !recoPrimary)
	continue;
      
      if(energyComp < comp_thresh || energyPurity < pur_thresh)
	continue;

      if(pfp_trackScore < .5)
	continue;

      if(abs(trackStartX) > 175 || abs(trackStartY) > 175 || trackStartZ < 25 || trackStartZ > 450)
        continue;

      if(!trackContained)
        continue;

      if(trk_length < 5.)
        continue;

      const int trueClass = razzledMap.at(abs(truePDG));
      const int recoClass = razzledMap.at(abs(chi2PDG));

      hConfusionMatrix->Fill(trueClass, recoClass);
    }

  TCanvas *c = new TCanvas("c", "c");
  c->cd();
  
  if(efficiency_mode)
    NormaliseEntriesByXTotal(hConfusionMatrix);
  if(purity_mode)
    NormaliseEntriesByYTotal(hConfusionMatrix);

  for(int i = 1; i <= hConfusionMatrix->GetNbinsX(); ++i)
    {
      hConfusionMatrix->GetXaxis()->SetBinLabel(i, axisLabels[i]);
      hConfusionMatrix->GetYaxis()->SetBinLabel(i, axisLabels[i]);
    }

  gStyle->SetPaintTextFormat("1.2g"),
  hConfusionMatrix->SetMarkerSize(3);
  hConfusionMatrix->Draw("col text");

  TString file_name = "chi2";

  if(efficiency_mode)
    file_name += "_efficiency";
  if(purity_mode)
    file_name += "_purity";
  if(raw_mode)
    file_name += "_raw";
  if(!require_primary)
    file_name += "_incl_non_reco_primaries";
  if(comp_thresh != -1)
    file_name += Form("_comp_thresh_%f", comp_thresh);
  if(pur_thresh != -1)
    file_name += Form("_pur_thresh_%f", pur_thresh);

  if(save)
    {
      c->SaveAs(save_dir + "/" + file_name + ".pdf");
      c->SaveAs(save_dir + "/" + file_name + ".png");
    }
}
