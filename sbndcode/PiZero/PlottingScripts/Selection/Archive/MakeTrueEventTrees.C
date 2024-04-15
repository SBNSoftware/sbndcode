#include "Enumerate.h"
#include "WeightNames.h"

void MakeTrueEventTrees(const TString productionVersion, const TString sampleName);

void MakeTrueEventTrees(const TString productionVersion)
{
  MakeTrueEventTrees(productionVersion, "rockbox");
  MakeTrueEventTrees(productionVersion, "intime");
}

void MakeTrueEventTrees(const TString productionVersion, const TString sampleName)
{
  const TString infile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_" + sampleName + ".root";

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(infile);

  int run, subrun, event;
  std::vector<int> *nu_event_type_incl = 0, *nu_event_type_0p0pi = 0, *nu_event_type_Np0pi = 0, *nu_mode = 0;
  std::vector<std::vector<double>> *nu_pz_pizero_mom = 0, *nu_pz_cos_theta_pizero = 0;
  std::vector<std::vector<std::vector<double>>*> flux_parameter_weights(flux_weight_names.size(), 0);

  events->SetBranchStatus("*", 0);

  events->SetBranchAddress("run", &run);
  events->SetBranchAddress("subrun", &subrun);
  events->SetBranchAddress("event", &event);

  events->SetBranchAddress("nu_event_type_incl", &nu_event_type_incl);
  events->SetBranchAddress("nu_event_type_0p0pi", &nu_event_type_0p0pi);
  events->SetBranchAddress("nu_event_type_Np0pi", &nu_event_type_Np0pi);
  events->SetBranchAddress("nu_mode", &nu_mode);
  events->SetBranchAddress("nu_pz_pizero_mom", &nu_pz_pizero_mom);
  events->SetBranchAddress("nu_pz_cos_theta_pizero", &nu_pz_cos_theta_pizero);

  for(auto&& [i, name] : enumerate(flux_weight_names))
    events->SetBranchAddress(("nu_weight_" + name).c_str(), &flux_parameter_weights[i]);

  TChain *subruns = new TChain("ncpizeroana/subruns");
  subruns->Add(infile);

  const TString outfilename = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/true_event_trees/" + productionVersion + "_" + sampleName + "_true_event_trees.root";
  TFile *outfile = new TFile(outfilename, "recreate");
  TTree *outnus = new TTree("nus","nus");
  TTree *outsubruns = (TTree*) subruns->Clone("subruns");

  int nuID, event_type_incl, event_type_0p0pi, event_type_Np0pi, mode;
  std::vector<double> pz_pizero_mom, pz_cos_theta_pizero;
  std::vector<std::vector<double>> flux_weights(flux_weight_names.size(), std::vector<double>());
  std::vector<double> flux_weights_all;

  outnus->Branch("run", &run);
  outnus->Branch("subrun", &subrun);
  outnus->Branch("event", &event);
  outnus->Branch("nuID", &nuID);
  outnus->Branch("nu_event_type_incl", &event_type_incl);
  outnus->Branch("nu_event_type_0p0pi", &event_type_0p0pi);
  outnus->Branch("nu_event_type_Np0pi", &event_type_Np0pi);
  outnus->Branch("nu_mode", &mode);
  outnus->Branch("nu_pz_pizero_mom", &pz_pizero_mom);
  outnus->Branch("nu_pz_cos_theta_pizero", &pz_cos_theta_pizero);

  for(auto&& [i, name] : enumerate(flux_weight_names))
    outnus->Branch(("nu_weight_" + name).c_str(), &flux_weights[i]);

  outnus->Branch("nu_weight_flux_all", &flux_weights_all);

  const int N = events->GetEntries();

  for(int i = 0; i < N; ++i)
    {
      if(!(i%10000))
        std::cout << "Event: " << i << " / " << N << std::endl;
      events->GetEntry(i);

      for(int j = 0; j < nu_event_type_incl->size(); ++j)
        {
          flux_weights_all.resize(1000, 1);

          nuID = j;

	  event_type_incl = nu_event_type_incl->at(j);
	  event_type_0p0pi = nu_event_type_0p0pi->at(j);
	  event_type_Np0pi = nu_event_type_Np0pi->at(j);

	  mode = nu_mode->at(j);

	  pz_pizero_mom = nu_pz_pizero_mom->at(j);
	  pz_cos_theta_pizero = nu_pz_cos_theta_pizero->at(j);

	  for(int l = 0; l < flux_weight_names.size(); ++l)
	    flux_weights[l] = flux_parameter_weights[l]->at(j);

	  if(nu_event_type_incl->at(j) != 6 && nu_event_type_incl->at(j) != 7)
	    {
	      for(int k = 0; k < 1000; ++k)
		{
		  double univ_weight = 1.;
		  for(int l = 0; l < flux_weight_names.size(); ++l)
		    univ_weight *= flux_parameter_weights.at(l)->at(j).at(k);

		  flux_weights_all[k] = univ_weight;
		}
	    }

	  outnus->Fill();
	}
    }

  outnus->Write();
  outsubruns->Write();
}
