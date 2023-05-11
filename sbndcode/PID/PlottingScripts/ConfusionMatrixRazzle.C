#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"

void ConfusionMatrixRazzle(const bool efficiency_mode = true, const bool purity_mode = false,
			   const bool require_primary = false, const double comp_thresh = -1,
			   const double pur_thresh = -1)
{
  const TString method_name = "BDT::BDTG";

  if(efficiency_mode && purity_mode)
    {
      std::cout << "Cannot run both efficiency and purity mode" << std::endl;
      return;
    }

  const TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/confusion_matrices_razzle";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  const TString weights_file = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_19_00/PID/Razzle.weights.xml";

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/razzled_trees.root");

  const std::map<int, int> razzledMap = { { 11, 0 }, { 22, 1 }, { 0, 2 } };
  std::vector<TString> axisLabels  = { "", "e^{#pm}", "#gamma", "other" };

  int truePdg;
  float energyComp, energyPurity, showerStartX, showerStartY, showerStartZ;
  bool recoPrimary, unambiguousSlice, showerContained;

  float pfp_trackScore;
  float shw_bestdEdx, shw_convGap, shw_openAngle, shw_modHitDensity, shw_sqrtEnergyDensity;

  tree->SetBranchAddress("pfp_trackScore", &pfp_trackScore);

  tree->SetBranchAddress("shw_bestdEdx", &shw_bestdEdx);
  tree->SetBranchAddress("shw_convGap", &shw_convGap);
  tree->SetBranchAddress("shw_openAngle", &shw_openAngle);
  tree->SetBranchAddress("shw_modHitDensity", &shw_modHitDensity);
  tree->SetBranchAddress("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);
  
  tree->SetBranchAddress("truePdg", &truePdg);
  tree->SetBranchAddress("energyComp", &energyComp);
  tree->SetBranchAddress("energyPurity", &energyPurity);
  tree->SetBranchAddress("recoPrimary", &recoPrimary);
  tree->SetBranchAddress("unambiguousSlice", &unambiguousSlice);

  tree->SetBranchAddress("showerStartX", &showerStartX);
  tree->SetBranchAddress("showerStartY", &showerStartY);
  tree->SetBranchAddress("showerStartZ", &showerStartZ);
  tree->SetBranchAddress("showerContained", &showerContained);

  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  reader->AddVariable("bestdEdx", &shw_bestdEdx);
  reader->AddVariable("convGap", &shw_convGap);
  reader->AddVariable("openAngle", &shw_openAngle);
  reader->AddVariable("modHitDensity", &shw_modHitDensity);
  reader->AddVariable("sqrtEnergyDensity", &shw_sqrtEnergyDensity);
  
  reader->BookMVA(method_name, weights_file);

  const unsigned N = tree->GetEntries();

  const unsigned n_categories = 3;
  TH2F *hConfusionMatrix = new TH2F("hConfusionMatrix", ";True Class;Razzle Class;",
				    n_categories, 0, n_categories, n_categories, 0, n_categories);
  
  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);
      
      if(unambiguousSlice)
       	continue;
      
      if(abs(truePdg) != 11 && abs(truePdg) != 22)
	truePdg = 0;

      if(require_primary && !recoPrimary)
	continue;
      
      if(energyComp < comp_thresh || energyPurity < pur_thresh)
	continue;

      if(pfp_trackScore > .5)
        continue;

      if(abs(showerStartX) > 175 || abs(showerStartY) > 175 || showerStartZ < 25 || showerStartZ > 450)
        continue;

      if(!showerContained)
        continue;

      shw_openAngle *= TMath::DegToRad();
      const std::vector<float> bdtScores = reader->EvaluateMulticlass(method_name);

      float bestScore = -std::numeric_limits<float>::max();
      int recoClass = -1;

      for(unsigned j = 0; j < bdtScores.size(); ++j)
	{
	  if(bdtScores[j] > bestScore)
	    {
	      bestScore = bdtScores[j];
	      recoClass = j;
	    }
	}
      
      const int trueClass = razzledMap.at(abs(truePdg));

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

  TString file_name = "razzle";

  if(efficiency_mode)
    file_name += "_efficiency";
  if(purity_mode)
    file_name += "_purity";
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
