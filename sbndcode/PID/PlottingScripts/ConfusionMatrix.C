#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"

void ConfusionMatrix(const TString weights_name = "Razzled_standard",
		     const TString method_name = "BDTG", const bool other_category = false,
		     const bool efficiency_mode = true, const bool purity_mode = false,
		     const bool require_primary = false, const double comp_thresh = -1,
		     const double pur_thresh = -1, const bool trad_tracks = false,
		     const bool trad_showers = false)
{
  if(efficiency_mode && purity_mode)
    {
      std::cout << "Cannot run both efficiency and purity mode" << std::endl;
      return;
    }

  if(trad_tracks && trad_showers)
    {
      std::cout << "Cannot run both traditional tracks and traditional showers simultaneously" << std::endl;
      return;
    }

  TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/confusion_matrices";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  const TString weights_file = "/sbnd/data/users/hlay/razzled/training/second_pass/" + weights_name + "/weights/" + weights_name + "_BDTG.weights.xml";

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/razzled_trees.root");

  const std::map<int, int> razzledMap = { { 11, 0 }, { 13, 1 }, { 22, 2 }, { 211, 3 }, { 2212, 4 }, { 0, 5 } };
  std::vector<TString> axisLabels  = { "", "e^{#pm}", "#mu^{#pm}", "#gamma", "#pi^{#pm}", "p" };
  if(other_category)
    axisLabels.push_back("Other");

  int truePDG;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ, showerStartX, showerStartY, showerStartZ, showerEnergy;
  bool recoPrimary, unambiguousSlice, trackContained, showerContained;

  float pfp_numDaughters, pfp_maxDaughterHits, pfp_trackScore, pfp_chargeEndFrac, pfp_chargeFracSpread,
    pfp_linearFitDiff, pfp_linearFitLength, pfp_linearFitGapLength, pfp_linearFitRMS, pfp_openAngleDiff,
    pfp_secondaryPCARatio, pfp_tertiaryPCARatio, pfp_vertexDist;

  float trk_length, trk_chi2PIDMuon, trk_chi2PIDProton, trk_chi2PIDMuonPionDiff, trk_mcsScatterMean,
    trk_mcsScatterMaxRatio, trk_momDiff, trk_meanDCA, trk_stoppingdEdxChi2Ratio, trk_chi2Pol0dEdxFit;

  float shw_bestdEdx, shw_convGap, shw_openAngle, shw_modHitDensity, shw_sqrtEnergyDensity;

  tree->SetBranchAddress("pfp_numDaughters", &pfp_numDaughters);
  tree->SetBranchAddress("pfp_maxDaughterHits", &pfp_maxDaughterHits);
  tree->SetBranchAddress("pfp_trackScore", &pfp_trackScore);
  tree->SetBranchAddress("pfp_chargeEndFrac", &pfp_chargeEndFrac);
  tree->SetBranchAddress("pfp_chargeFracSpread", &pfp_chargeFracSpread);
  tree->SetBranchAddress("pfp_linearFitDiff", &pfp_linearFitDiff);
  tree->SetBranchAddress("pfp_linearFitLength", &pfp_linearFitLength);
  tree->SetBranchAddress("pfp_linearFitGapLength", &pfp_linearFitGapLength);
  tree->SetBranchAddress("pfp_linearFitRMS", &pfp_linearFitRMS);
  tree->SetBranchAddress("pfp_openAngleDiff", &pfp_openAngleDiff);
  tree->SetBranchAddress("pfp_secondaryPCARatio", &pfp_secondaryPCARatio);
  tree->SetBranchAddress("pfp_tertiaryPCARatio", &pfp_tertiaryPCARatio);
  tree->SetBranchAddress("pfp_vertexDist", &pfp_vertexDist);

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

  tree->SetBranchAddress("shw_bestdEdx", &shw_bestdEdx);
  tree->SetBranchAddress("shw_convGap", &shw_convGap);
  tree->SetBranchAddress("shw_openAngle", &shw_openAngle);
  tree->SetBranchAddress("shw_modHitDensity", &shw_modHitDensity);
  tree->SetBranchAddress("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);
  
  tree->SetBranchAddress("truePDG", &truePDG);
  tree->SetBranchAddress("energyComp", &energyComp);
  tree->SetBranchAddress("energyPurity", &energyPurity);
  tree->SetBranchAddress("recoPrimary", &recoPrimary);
  tree->SetBranchAddress("unambiguousSlice", &unambiguousSlice);

  tree->SetBranchAddress("trackStartX", &trackStartX);
  tree->SetBranchAddress("trackStartY", &trackStartY);
  tree->SetBranchAddress("trackStartZ", &trackStartZ);
  tree->SetBranchAddress("trackContained", &trackContained);

  tree->SetBranchAddress("showerStartX", &showerStartX);
  tree->SetBranchAddress("showerStartY", &showerStartY);
  tree->SetBranchAddress("showerStartZ", &showerStartZ);
  tree->SetBranchAddress("showerContained", &showerContained);
  tree->SetBranchAddress("showerEnergy", &showerEnergy);

  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  reader->AddVariable("pfp_numDaughters", &pfp_numDaughters);
  reader->AddVariable("pfp_maxDaughterHits", &pfp_maxDaughterHits);
  if(weights_name != "Razzled_no_track_score")
    {
      reader->AddVariable("pfp_trackScore", &pfp_trackScore);
      if(weights_name != "Razzled_no_track_score_inputs")
	{
	  reader->AddVariable("pfp_chargeEndFrac", &pfp_chargeEndFrac);
	  reader->AddVariable("pfp_chargeFracSpread", &pfp_chargeFracSpread);
	  reader->AddVariable("pfp_linearFitDiff", &pfp_linearFitDiff);
	  reader->AddVariable("pfp_linearFitLength", &pfp_linearFitLength);
	  reader->AddVariable("pfp_linearFitGapLength", &pfp_linearFitGapLength);
	  reader->AddVariable("pfp_linearFitRMS", &pfp_linearFitRMS);
	  reader->AddVariable("pfp_openAngleDiff", &pfp_openAngleDiff);
	  reader->AddVariable("pfp_secondaryPCARatio", &pfp_secondaryPCARatio);
	  reader->AddVariable("pfp_tertiaryPCARatio", &pfp_tertiaryPCARatio);
	  reader->AddVariable("pfp_vertexDist", &pfp_vertexDist);
	}
    }

  reader->AddVariable("trk_length", &trk_length);
  reader->AddVariable("trk_chi2PIDMuon", &trk_chi2PIDMuon);
  reader->AddVariable("trk_chi2PIDProton", &trk_chi2PIDProton);
  reader->AddVariable("trk_chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
  reader->AddVariable("trk_mcsScatterMean", &trk_mcsScatterMean);
  reader->AddVariable("trk_mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
  reader->AddVariable("trk_meanDCA", &trk_meanDCA);
  reader->AddVariable("trk_stoppingdEdxChi2Ratio", &trk_stoppingdEdxChi2Ratio);
  reader->AddVariable("trk_chi2Pol0dEdxFit", &trk_chi2Pol0dEdxFit);
  reader->AddVariable("trk_momDiff", &trk_momDiff);

  reader->AddVariable("shw_bestdEdx", &shw_bestdEdx);
  reader->AddVariable("shw_convGap", &shw_convGap);
  reader->AddVariable("shw_openAngle", &shw_openAngle);
  reader->AddVariable("shw_modHitDensity", &shw_modHitDensity);
  reader->AddVariable("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);
  
  reader->BookMVA(method_name, weights_file);

  const unsigned N = tree->GetEntries();

  const unsigned n_categories = other_category ? 6 : 5;
  TH2F *hConfusionMatrix = new TH2F("hConfusionMatrix", ";True Class;Razzled Class;",
				    n_categories, 0, n_categories, n_categories, 0, n_categories);
  
  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);
      
      if(unambiguousSlice)
	continue;
      
      if(abs(truePDG) != 11 && abs(truePDG) != 13 && abs(truePDG) != 22 
	 && abs(truePDG) != 211 && abs(truePDG) != 2212)
	{
	  if(other_category)
	    truePDG = 0;
	  else
	    continue;
	}

      if(require_primary && !recoPrimary)
	continue;
      
      if(energyComp < comp_thresh || energyPurity < pur_thresh)
	continue;

      if(trad_tracks && pfp_trackScore < .5)
	continue;

      if(trad_showers && pfp_trackScore > .5)
	continue;

      if(abs(trackStartX) > 175 || abs(trackStartY) > 175 || trackStartZ < 25 || trackStartZ > 450 ||
	 abs(showerStartX) > 175 || abs(showerStartY) > 175 || showerStartZ < 25 || showerStartZ > 450)
	continue;

      if(!trackContained || !showerContained)
	continue;

      if(trk_length < 5 || showerEnergy < 10)
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
      
      const int trueClass = razzledMap.at(abs(truePDG));

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

  TString file_name = weights_name;
  file_name.ToLower();

  save_dir += "/" + file_name;

  if(save)
    gSystem->Exec("mkdir -p " + save_dir);

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
  if(trad_tracks)
    file_name += "_traditional_tracks";
  if(trad_showers)
    file_name += "_traditional_showers";

  if(save)
    {
      c->SaveAs(save_dir + "/" + file_name + ".pdf");
      c->SaveAs(save_dir + "/" + file_name + ".png");
    }
}
