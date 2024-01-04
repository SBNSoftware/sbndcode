void ListEvents(const TString productionVersion)
{
  const TString rockboxfile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  TChain *events = new TChain("ncpizeroana/events");
  events->Add(rockboxfile);

  int run, subrun, event;
  std::vector<int> *slc_true_event_type = 0, *slc_true_pdg = 0, *slc_true_ccnc = 0, *slc_true_n_neutral_pions = 0,
    *slc_n_razzled_muons = 0, *slc_n_razzled_photons = 0;
  std::vector<bool> *slc_true_fv = 0, *slc_true_av = 0, *slc_is_clear_cosmic = 0, *slc_is_fv = 0;
  std::vector<float> *slc_comp = 0, *slc_crumbs_score = 0;
  std::vector<size_t> *slc_n_pfps = 0;

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

  events->SetBranchAddress("slc_is_clear_cosmic", &slc_is_clear_cosmic);
  events->SetBranchAddress("slc_is_fv", &slc_is_fv);
  events->SetBranchAddress("slc_crumbs_score", &slc_crumbs_score);
  events->SetBranchAddress("slc_n_razzled_muons", &slc_n_razzled_muons);
  events->SetBranchAddress("slc_n_pfps", &slc_n_pfps);
  events->SetBranchAddress("slc_n_razzled_photons", &slc_n_razzled_photons);

  const int N = events->GetEntries();

  for(int i = 0; i < N; ++i)
    {
      events->GetEntry(i);

      for(int j = 0; j < slc_is_clear_cosmic->size(); ++j)
        {
          int category = -1;

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

          if(!slc_is_clear_cosmic->at(j) && slc_is_fv->at(j) && slc_crumbs_score->at(j) > -0.025 && slc_n_razzled_muons->at(j) == 0 && slc_n_pfps->at(j) > 1 && slc_n_razzled_photons->at(j) < 2 && category == 0)
            std::cout << run << " " << subrun << " " << event << std::endl;
        }
    }
}
