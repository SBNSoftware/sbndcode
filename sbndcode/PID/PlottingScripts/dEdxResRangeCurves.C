#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void dEdxResRangeCurves(const int determined_id = -1)
{
  const std::map<int, TString> pdgNameMap = { { 13, "Muon" },
					      { 211, "Pion" },
					      { 321, "Kaon" },
					      { 2212, "Proton" }
  };

  const std::map<int, TString> pdgLabelMap = { { 13, "#mu^{#pm}" },
					       { 211, "#pi^{#pm}" },
					       { 321, "K^{#pm}" },
					       { 2212, "p" }
  };
  
  TString save_dir = "/sbnd/data/users/hlay/razzled/plots/investigations/dEdxresrangecurves_test";
  TString determined_name = "";
  if(determined_id != -1){
    determined_name = pdgNameMap.at(determined_id);
    determined_name.ToLower();
    save_dir += "/" + determined_name + "_id";
  }

  const bool save = true;
  if(save)
    gSystem->Exec("mkdir -p " + save_dir);

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pandoraRazzled/pfpTree");
  tree->Add("/sbnd/data/users/hlay/razzled/bnb/razzled*.root");
  tree->Add("/sbnd/data/users/hlay/razzled/intrnue/razzled*.root");

  TFile* template_file = TFile::Open("/cvmfs/larsoft.opensciencegrid.org/products/larsoft_data/v1_02_02/ParticleIdentification/dEdxrestemplates.root");

  TProfile* muon_template   = (TProfile*) template_file->Get("dedx_range_mu");
  muon_template->SetLineColor(kRed);
  TProfile* pion_template   = (TProfile*) template_file->Get("dedx_range_pi");
  pion_template->SetLineColor(kOrange-3);
  TProfile* kaon_template   = (TProfile*) template_file->Get("dedx_range_ka");
  kaon_template->SetLineColor(kGreen+1);
  TProfile* proton_template = (TProfile*) template_file->Get("dedx_range_pro");
  proton_template->SetLineColor(kPink-9);

  std::map<int, TProfile*> profiles = { { 13, muon_template },
					{ 211, pion_template },
					{ 321, kaon_template },
					{ 2212, proton_template },
  };
  
  TPaveText *simText = new TPaveText(.75,.91,.85,.95, "NB NDC");
  simText->SetTextAlign(32);
  simText->SetTextSize(0.03);
  simText->SetTextColor(kGray+2);
  simText->SetFillColor(kWhite);
  simText->AddText("SBND Simulation");

  TLegend *l = new TLegend(.62, .76, .79, .88, "Theory Prediction", "NDC C");
  l->SetNColumns(2);
  l->SetFillColor(kWhite);
    
  for(auto const& [pdg, profile] : profiles)
    {
      profile->Rebin(8);
      l->AddEntry(profile, pdgLabelMap.at(pdg), "l");
    }

  TLegendEntry *header = (TLegendEntry*)l->GetListOfPrimitives()->First();
  header->SetTextAlign(22);
  header->SetTextSize(.035);

  int truePDG, chi2PDG;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ, trueEndMomentum;
  bool recoPrimary, unambiguousSlice, trackContained;

  float pfp_trackScore, trk_length;
  std::string *trueEndProcess = 0;

  std::vector<float> *trackdEdx = 0, *trackResRange = 0;

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

  tree->SetBranchAddress("pfp_trackScore", &pfp_trackScore);
  tree->SetBranchAddress("trk_length", &trk_length);

  tree->SetBranchAddress("trackdEdx", &trackdEdx);
  tree->SetBranchAddress("trackResRange", &trackResRange);
  tree->SetBranchAddress("chi2PDG", &chi2PDG);

  std::map<int, TH2F*> hists;
  hists[13]    = new TH2F("dEdxResRangeMuon", "#mu^{#pm};Residual Range (cm);dE/dx (MeV/cm);",
			 200, 0, 50, 150, 0, 15);
  hists[211]   = new TH2F("dEdxResRangePion", "#pi^{#pm};Residual Range (cm);dE/dx (MeV/cm);",
			 200, 0, 50, 150, 0, 15);
  hists[321]  = new TH2F("dEdxResRangeKaon", "K^{#pm};Residual Range (cm);dE/dx (MeV/cm);",
			 200, 0, 50, 150, 0, 15);
  hists[2212] = new TH2F("dEdxResRangeProton", "p;Residual Range (cm);dE/dx (MeV/cm);",
			 200, 0, 50, 150, 0, 15);

  const unsigned N = tree->GetEntries();

  for(unsigned i = 0; i < N; ++i)
    {
      tree->GetEntry(i);

      if(unambiguousSlice)
	continue;

      if(determined_id != -1 && chi2PDG != determined_id)
	continue;

      if(abs(truePDG) != 13 && abs(truePDG) != 211 && abs(truePDG) != 321 && abs(truePDG) != 2212)
	continue;

      if(abs(truePDG) == 211 && trueEndProcess->find("Inel") != std::string::npos)
      	continue;

      if(abs(truePDG) == 321 && trueEndProcess->find("Inel") != std::string::npos)
      	continue;

      if(abs(truePDG) == 13 && (trueEndProcess->find("Transp") != std::string::npos 
				|| (trueEndProcess->find("Decay") != std::string::npos && trueEndMomentum > 0)))
      	continue;

      if(abs(truePDG) == 2212 && trueEndProcess->find("Inel") != std::string::npos)
	continue;

      if(!recoPrimary)
	continue;

      if(energyComp < -1. || energyPurity < -1.)
        continue;

      if(pfp_trackScore < .5)
	continue;

      if(abs(trackStartX) > 175 || abs(trackStartY) > 175 || trackStartZ < 25 || trackStartZ > 450)
        continue;

      if(!trackContained)
        continue;

      if(trk_length < 5)
        continue;

      if(trackdEdx->size() != trackResRange->size())
	continue;

      if(trackdEdx->size() == 0)
	continue;

      for(unsigned j = 1; j < trackdEdx->size() - 1; ++j)
	hists[std::abs(truePDG)]->Fill(trackResRange->at(j), trackdEdx->at(j));
    }

  for(auto const& [pdg, hist] : hists)
    {
      TString name = pdgNameMap.at(pdg);
      
      TCanvas *c = new TCanvas("c" + name, "c" + name);
      c->cd();

      c->SetTopMargin(.1);
      c->SetRightMargin(.15);

      if(determined_id != -1)
	hist->SetTitle(pdgLabelMap.at(pdg) + " ID'd as " + pdgLabelMap.at(determined_id) + ";Residual Range (cm);dE/dx (MeV/cm);");

      hist->Draw("colz");

      simText->Draw();
      name.ToLower();

      if(save)
	{
	  c->SaveAs(save_dir + "/" + name + ".png");
	  c->SaveAs(save_dir + "/" + name + ".pdf");
	}

      profiles[pdg]->Draw("histLsame");

      if(save)
	{
	  c->SaveAs(save_dir + "/" + name + "_with_predicted.png");
	  c->SaveAs(save_dir + "/" + name + "_with_predicted.pdf");
	}

      for(auto const& [pdg_prof, profile] : profiles)
	profile->Draw("histLsame");

      l->Draw();

      if(save)
	{
	  c->SaveAs(save_dir + "/" + name + "_with_all_predicted.png");
	  c->SaveAs(save_dir + "/" + name + "_with_all_predicted.pdf");
	}
    }
}
