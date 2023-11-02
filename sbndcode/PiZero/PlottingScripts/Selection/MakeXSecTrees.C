#include "Enumerate.h"
#include "WeightNames.h"

void MakeXSecTrees(const TString productionVersion, const TString sampleName);

void MakeXSecTrees(const TString productionVersion)
{
  MakeXSecTrees(productionVersion, "rockbox");
  MakeXSecTrees(productionVersion, "intime");
}

void MakeXSecTrees(const TString productionVersion, const TString sampleName)
{
  const TString infile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_" + sampleName + ".root";
  TChain *events = new TChain("ncpizeroana/events");
  events->Add(infile);

  int run, subrun, event;
  std::vector<int> *slc_true_event_type = 0, *slc_true_pdg = 0, *slc_true_ccnc = 0, *slc_true_n_neutral_pions = 0,
    *slc_n_razzled_muons = 0, *slc_n_razzled_photons = 0;
  std::vector<bool> *slc_true_fv = 0, *slc_true_av = 0, *slc_is_clear_cosmic = 0, *slc_is_fv = 0;
  std::vector<float> *slc_comp = 0, *slc_crumbs_score = 0;
  std::vector<size_t> *slc_n_pfps = 0;
  std::vector<double> *slc_best_pzc_pizero_mom = 0;
  std::vector<std::vector<std::vector<double>>*> flux_parameter_weights(flux_weight_names.size(), 0);

  events->SetBranchStatus("*", 0);

  events->SetBranchAddress("run", &run);
  events->SetBranchAddress("subrun", &subrun);
  events->SetBranchAddress("event", &event);

  events->SetBranchAddress("slc_true_event_type", &slc_true_event_type);
  events->SetBranchAddress("slc_true_fv", &slc_true_fv);
  events->SetBranchAddress("slc_true_av", &slc_true_av);
  events->SetBranchAddress("slc_true_pdg", &slc_true_pdg);
  events->SetBranchAddress("slc_true_ccnc", &slc_true_ccnc);
  events->SetBranchAddress("slc_true_n_neutral_pions", &slc_true_n_neutral_pions);
  events->SetBranchAddress("slc_comp", &slc_comp);
  events->SetBranchAddress("slc_best_pzc_pizero_mom", &slc_best_pzc_pizero_mom);

  events->SetBranchAddress("slc_is_clear_cosmic", &slc_is_clear_cosmic);
  events->SetBranchAddress("slc_is_fv", &slc_is_fv);
  events->SetBranchAddress("slc_crumbs_score", &slc_crumbs_score);
  events->SetBranchAddress("slc_n_razzled_muons", &slc_n_razzled_muons);
  events->SetBranchAddress("slc_n_pfps", &slc_n_pfps);
  events->SetBranchAddress("slc_n_razzled_photons", &slc_n_razzled_photons);

  for(auto&& [i, name] : enumerate(flux_weight_names))
    events->SetBranchAddress(("slc_true_weight_" + name).c_str(), &flux_parameter_weights[i]);

  TChain *subruns = new TChain("ncpizeroana/subruns");
  subruns->Add(infile);

  const TString outfilename = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/xsec_trees/" + productionVersion + "_" + sampleName + "_xsec_trees.root";
  TFile *outfile = new TFile(outfilename, "recreate");
  TTree *outslices = new TTree("slices","slices");
  TTree *outsubruns = (TTree*) subruns->Clone("subruns");

  int sliceID, category;
  double pzc_pizero_mom;
  bool selected;
  std::vector<double> flux_weights;

  outslices->Branch("run", &run);
  outslices->Branch("subrun", &subrun);
  outslices->Branch("event", &event);
  outslices->Branch("sliceID", &sliceID);
  outslices->Branch("category", &category);
  outslices->Branch("pzc_pizero_mom", &pzc_pizero_mom);
  outslices->Branch("flux_weights", &flux_weights);
  outslices->Branch("selected", &selected);

  const int N = events->GetEntries();

  for(int i = 0; i < N; ++i)
    {
      if(!(i%10000))
        std::cout << "Event: " << i << " / " << N << std::endl;
      events->GetEntry(i);

      for(int j = 0; j < slc_is_clear_cosmic->size(); ++j)
        {
          flux_weights.resize(1000, 1);

          sliceID = j;

          if(slc_true_event_type->at(j) == 6 || slc_true_event_type->at(j) == 7)
            category = 6;
          else if(slc_true_fv->at(j) && slc_true_ccnc->at(j) == 1 && slc_true_n_neutral_pions->at(j) > 0 && slc_comp->at(j) > .5)
            category = 0;
          else if(slc_true_fv->at(j) && slc_true_ccnc->at(j) == 1 && slc_true_n_neutral_pions->at(j) == 0)
            category = 1;
          else if(slc_true_fv->at(j) && slc_true_ccnc->at(j) == 0 && abs(slc_true_pdg->at(j)) == 14)
            category = 2;
          else if(slc_true_fv->at(j) && slc_true_ccnc->at(j) == 0 && abs(slc_true_pdg->at(j)) == 12)
            category = 3;
          else if(!slc_true_fv->at(j) && !slc_true_av->at(j))
            category = 4;
          else if(!slc_true_fv->at(j) && slc_true_av->at(j))
            category = 5;
          else if(slc_true_fv->at(j) && slc_true_ccnc->at(j) == 1 && slc_true_n_neutral_pions->at(j) > 0)
            category = 7;

          selected = !slc_is_clear_cosmic->at(j) && slc_is_fv->at(j) && slc_crumbs_score->at(j) > -0.025
            && slc_n_razzled_muons->at(j) == 0 && slc_n_pfps->at(j) > 1 && slc_n_razzled_photons->at(j) > 1;

          if(category == 0 || selected)
            {
              pzc_pizero_mom = slc_best_pzc_pizero_mom->at(j);

              if(category != 6)
                {
                  for(int k = 0; k < 1000; ++k)
                    {
                      double univ_weight = 1.;
                      for(int l = 0; l < flux_weight_names.size(); ++l)
                        univ_weight *= flux_parameter_weights.at(l)->at(j).at(k);

                      flux_weights[k] = univ_weight;
                    }
                }

              outslices->Fill();
            }
        }
    }

  outslices->Write();
  outsubruns->Write();
}
