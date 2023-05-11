#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"

void ScoreDistributionsRazzle(const bool require_primary = false, const double comp_thresh = -1,
			      const double pur_thresh = -1)
{
  const TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/score_distributions_razzle";
  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);
  
  const TString method_name = "BDT::BDTG";
  const TString weights_file = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_19_00/PID/Razzle.weights.xml";

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
					     { 22, 1, kBlue+2, "#gamma", "Photon" },
					     { 0, 2, kBlack, "Other", "Other" }
  };

  std::map<int, int> razzledMap = { { 11, 0 }, { 22, 1 }, { 0, 2 } };
  
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
      
      if(abs(truePdg) != 11 && abs(truePdg) != 22)
	truePdg = 0;

      if(require_primary && !recoPrimary)
	continue;
      
      if(energyComp < comp_thresh || energyPurity < pur_thresh)
	continue;

      if(abs(showerStartX) > 175 || abs(showerStartY) > 175 || showerStartZ < 25 || showerStartZ > 450)
	continue;

      if(!showerContained)
	continue;

      if(pfp_trackScore > .5)
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
