#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"

void ScoreDistributions(const TString weights_name = "Razzled_standard",
			const TString method_name = "BDT::BDTG", const bool other_category = false,
			const bool require_primary = false, const double comp_thresh = -1,
			const double pur_thresh = -1)
{
  const TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/score_distributions";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  const TString weights_file = "/sbnd/data/users/hlay/razzled/training/second_pass/" + weights_name + "/weights/" + weights_name + "_BDTG.weights.xml";

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/razzled_trees.root");

  struct RazzledParticle
  {
    int pdg;
    int id;
    int colour;
    TString label;
    TString name;
  };

  std::vector<RazzledParticle> particles = { { 11, 0, kMagenta+2, "e^{#pm}", "Electron" }, 
					     { 13, 1, kRed+2, "#mu^{#pm}", "Muon" },
					     { 22, 2, kBlue+2, "#gamma", "Photon" },
					     { 211, 3, kGreen+2, "#pi^{#pm}", "Pion" },
					     { 2212, 4, kOrange+2, "p", "Proton" }
  };

  std::map<int, int> razzledMap = { { 11, 0 }, { 13, 1 }, { 22, 2 }, { 211, 3 }, { 2212, 4 }, { 0, 5 } };
  
  if(other_category)
    particles.push_back({ 0, 5, kBlack, "Other", "Other" });

  int truePDG;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ, showerStartX, showerStartY, showerStartZ;
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

  std::vector<std::vector<TH1F*>> hScores;

  for(unsigned i = 0; i < particles.size(); ++i)
    {
      hScores.push_back(std::vector<TH1F*>());
      for(unsigned j = 0; j < particles.size(); ++j)
	hScores[i].push_back(new TH1F("h" + particles[i].name + "Score" + particles[j].name, ";" + particles[i].name + " Score;PFPs",
				      50, 0, 1));
    }
  
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

      if(abs(trackStartX) > 175 || abs(trackStartY) > 175 || trackStartZ < 25 || trackStartZ > 450 ||
	 abs(showerStartX) > 175 || abs(showerStartY) > 175 || showerStartZ < 25 || showerStartZ > 450)
	continue;

      if(!trackContained || !showerContained)
	continue;

      const std::vector<float> bdtScores = reader->EvaluateMulticlass(method_name);
      const int trueClass = razzledMap.at(abs(truePDG));

      for(unsigned j = 0; j < bdtScores.size(); ++j)
	hScores[j][trueClass]->Fill(bdtScores[j]);
    }

  for(unsigned i = 0; i < particles.size(); ++i)
    {
      TString lower_case_name = particles[i].name;
      lower_case_name.ToLower();

      TCanvas *c = new TCanvas("c_" + lower_case_name, "c_" + lower_case_name);
      c->cd();

      TLegend *legend = new TLegend(.35, .77, .8, .87);
      legend->SetBorderSize(0);
      legend->SetNColumns(2);
	  
      double max = -std::numeric_limits<double>::max();
      
      for(unsigned j = 0; j < particles.size(); ++j)
	{
	  hScores[i][j]->Scale(1./hScores[i][j]->GetEntries());

	  if(hScores[i][j]->GetMaximum() > max)
	    max = hScores[i][j]->GetMaximum();
	}

      for(unsigned j = 0; j < particles.size(); ++j)
	{
	  if(j == 0)
	    hScores[i][j]->SetMaximum(1.2 * max);

	  hScores[i][j]->SetLineColor(particles[j].colour);
	  hScores[i][j]->Draw("histEsame");

	  legend->AddEntry(hScores[i][j], particles[j].label, "lep");
	}

      legend->Draw();

      TString file_name = weights_name;
      file_name.ToLower();
      file_name += ("_" + lower_case_name);
      
      if(!require_primary)
	file_name += "_incl_non_reco_primaries";
      if(comp_thresh != -1)
	file_name += Form("_comp_thresh_%f", comp_thresh);
      if(pur_thresh != -1)
	file_name += Form("_pur_thresh_%f", pur_thresh);

      if(save)
	{
	  c->SaveAs(save_dir + "/" + file_name + ".png");
	  c->SaveAs(save_dir + "/" + file_name + ".pdf");
	}
    }
}
