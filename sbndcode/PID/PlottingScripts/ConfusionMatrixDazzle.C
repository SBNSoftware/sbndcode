#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"

void ConfusionMatrixDazzle(const bool efficiency_mode = true, const bool purity_mode = false,
			   const bool require_primary = false, const double comp_thresh = -1,
			   const double pur_thresh = -1)
{
  const TString method_name = "BDT::BDTG";

  if(efficiency_mode && purity_mode)
    {
      std::cout << "Cannot run both efficiency and purity mode" << std::endl;
      return;
    }

  const TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/confusion_matrices_dazzle";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  const TString weights_file = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_19_00/PID/Dazzle.weights.xml";

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/razzled_trees.root");

  const std::map<int, int> razzledMap = { { 13, 0 }, { 211, 1 }, { 2212, 2 }, { 0, 3 } };
  std::vector<TString> axisLabels  = { "", "#mu^{#pm}", "#pi^{#pm}", "p", "other" };

  int truePdg;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ;
  bool recoPrimary, unambiguousSlice, trackContained;

  float pfp_numDaughters, pfp_maxDaughterHits, pfp_trackScore;

  float trk_length, trk_chi2PIDMuon, trk_chi2PIDProton, trk_chi2PIDMuonPionDiff, trk_mcsScatterMean,
    trk_mcsScatterMaxRatio, trk_momDiff, trk_meanDCA, trk_stoppingdEdxChi2Ratio, trk_chi2Pol0dEdxFit;

  tree->SetBranchAddress("pfp_numDaughters", &pfp_numDaughters);
  tree->SetBranchAddress("pfp_maxDaughterHits", &pfp_maxDaughterHits);
  tree->SetBranchAddress("pfp_trackScore", &pfp_trackScore);

  tree->SetBranchAddress("trk_length", &trk_length);
  tree->SetBranchAddress("trk_chi2PIDMuon", &trk_chi2PIDMuon);
  tree->SetBranchAddress("trk_chi2PIDProton", &trk_chi2PIDProton);
  tree->SetBranchAddress("trk_chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
  tree->SetBranchAddress("trk_mcsScatterMean", &trk_mcsScatterMean);
  tree->SetBranchAddress("trk_mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
  tree->SetBranchAddress("trk_meanDCA", &trk_meanDCA);
  tree->SetBranchAddress("trk_stoppingdEdxChi2Ratio", &trk_stoppingdEdxChi2Ratio);
  tree->SetBranchAddress("trk_chi2Pol0dEdxFit", &trk_chi2Pol0dEdxFit);
  tree->SetBranchAddress("trk_momDiff", &trk_momDiff);

  tree->SetBranchAddress("truePdg", &truePdg);
  tree->SetBranchAddress("energyComp", &energyComp);
  tree->SetBranchAddress("energyPurity", &energyPurity);
  tree->SetBranchAddress("recoPrimary", &recoPrimary);
  tree->SetBranchAddress("unambiguousSlice", &unambiguousSlice);

  tree->SetBranchAddress("trackStartX", &trackStartX);
  tree->SetBranchAddress("trackStartY", &trackStartY);
  tree->SetBranchAddress("trackStartZ", &trackStartZ);
  tree->SetBranchAddress("trackContained", &trackContained);

  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  reader->AddVariable("recoLen", &trk_length);
  reader->AddVariable("chi2PIDMuon", &trk_chi2PIDMuon);
  reader->AddVariable("chi2PIDProton", &trk_chi2PIDProton);
  reader->AddVariable("chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
  reader->AddVariable("mcsScatterMean", &trk_mcsScatterMean);
  reader->AddVariable("mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
  reader->AddVariable("meanDCA", &trk_meanDCA);
  reader->AddVariable("stoppingChi2Ratio", &trk_stoppingdEdxChi2Ratio);
  reader->AddVariable("chi2Pol0Fit", &trk_chi2Pol0dEdxFit);
  reader->AddVariable("pDiff", &trk_momDiff);
  reader->AddVariable("numDaughters", &pfp_numDaughters);
  reader->AddVariable("maxDaughterHits", &pfp_maxDaughterHits);


  reader->BookMVA(method_name, weights_file);

  const unsigned N = tree->GetEntries();

  const unsigned n_categories = 4;
  TH2F *hConfusionMatrix = new TH2F("hConfusionMatrix", ";True Class;Dazzle Class;",
				    n_categories, 0, n_categories, n_categories, 0, n_categories);
  
  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);
      
      if(unambiguousSlice)
	continue;
      
      if(abs(truePdg) != 13 && abs(truePdg) != 211 && abs(truePdg) != 2212)
	truePdg = 0;

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

  TString file_name = "dazzle";

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
