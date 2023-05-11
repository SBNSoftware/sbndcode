#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"

void ScoreDistributionsDazzle(const bool require_primary = false, const double comp_thresh = -1,
			      const double pur_thresh = -1)
{
  const TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/score_distributions_dazzle";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  const TString method_name = "BDT::BDTG";
  const TString weights_file = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_19_00/PID/Dazzle.weights.xml";
  
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

  std::vector<RazzledParticle> particles = { { 13, 0, kRed+2, "#mu^{#pm}", "Muon" },
					     { 211, 1, kGreen+2, "#pi^{#pm}", "Pion" },
					     { 2212, 2, kOrange+2, "p", "Proton" },
					     { 0, 3, kBlack, "Other", "Other" }
  };

  std::map<int, int> razzledMap = { { 13, 0 }, { 211, 1 }, { 2212, 2 }, { 0, 3 } };
  
  int truePdg;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ;
  bool recoPrimary, unambiguousSlice, trackContained;

  float pfp_trackScore, pfp_numDaughters, pfp_maxDaughterHits;

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
      
      if(abs(truePdg) != 13 && abs(truePdg) != 211 && abs(truePdg) != 2212)
	truePdg = 0;

      if(require_primary && !recoPrimary)
	continue;
      
      if(energyComp < comp_thresh || energyPurity < pur_thresh)
	continue;

      if(abs(trackStartX) > 175 || abs(trackStartY) > 175 || trackStartZ < 25 || trackStartZ > 450)
	continue;

      if(!trackContained)
	continue;

      if(pfp_trackScore < .5)
	continue;

      const std::vector<float> bdtScores = reader->EvaluateMulticlass(method_name);
      const int trueClass = razzledMap.at(abs(truePdg));

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
	  hScores[i][j]->DrawNormalized("histEsame");

	  legend->AddEntry(hScores[i][j], particles[j].label, "lep");
	}

      legend->Draw();

      TString file_name = lower_case_name;
      
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
