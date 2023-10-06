void MakeMultiSliceTrees()
{
  const TString infile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root";
  TChain *events = new TChain("ncpizeroana/events");
  events->Add(infile);

  int run, subrun, event;
  std::vector<bool> *slc_is_clear_cosmic = 0;
  std::vector<size_t> *slc_true_mctruth_id = 0;
  std::vector<int> *nu_event_type = 0, *slc_true_event_type = 0, *slc_n_dazzle_muons_cut_based = 0;
  std::vector<double> *slc_comp = 0, *slc_pur = 0, *slc_vtx_x = 0, *slc_vtx_y = 0, *slc_vtx_z = 0,
    *slc_opt0_measPE = 0, *slc_opt0_hypPE = 0, *slc_crumbs_score = 0, *slc_crumbs_nc_score = 0,
    *slc_crumbs_ccnue_score = 0;

  events->SetBranchAddress("run", &run);
  events->SetBranchAddress("subrun", &subrun);
  events->SetBranchAddress("event", &event);
  events->SetBranchAddress("nu_event_type", &nu_event_type);
  events->SetBranchAddress("slc_is_clear_cosmic", &slc_is_clear_cosmic);
  events->SetBranchAddress("slc_true_event_type", &slc_true_event_type);
  events->SetBranchAddress("slc_true_mctruth_id", &slc_true_mctruth_id);
  events->SetBranchAddress("slc_comp", &slc_comp);
  events->SetBranchAddress("slc_pur", &slc_pur);
  events->SetBranchAddress("slc_n_dazzle_muons_cut_based", &slc_n_dazzle_muons_cut_based);
  events->SetBranchAddress("slc_vtx_x", &slc_vtx_x);
  events->SetBranchAddress("slc_vtx_y", &slc_vtx_y);
  events->SetBranchAddress("slc_vtx_z", &slc_vtx_z);
  events->SetBranchAddress("slc_opt0_measPE", &slc_opt0_measPE);
  events->SetBranchAddress("slc_opt0_hypPE", &slc_opt0_hypPE);
  events->SetBranchAddress("slc_crumbs_score", &slc_crumbs_score);
  events->SetBranchAddress("slc_crumbs_nc_score", &slc_crumbs_nc_score);
  events->SetBranchAddress("slc_crumbs_ccnue_score", &slc_crumbs_ccnue_score);

  const TString outfilename = "/sbnd/data/users/hlay/ncpizero/split_slice_selection/NCPiZeroAv2_rockbox_multislice_tree.root";
  TFile *outfile = new TFile(outfilename, "recreate");
  TTree *outevents = new TTree("events", "events");

  bool   signal;
  size_t true_slc0, true_slc1, signus;

  std::vector<size_t> slc0, slc1;
  std::vector<double> vtx_x0, vtx_y0, vtx_z0, crumbs0, crumbs_nc0, crumbs_ccnue0, opT0Frac0,
    vtx_x1, vtx_y1, vtx_z1, crumbs1, crumbs_nc1, crumbs_ccnue1, opT0Frac1, sep, sumOpT0Frac;
  std::vector<int> n_dazzle_muons0, n_dazzle_muons1;
  std::vector<bool> matching_flash_pe, signal_cand, good_cand, goodOpT00, goodOpT01;

  outevents->Branch("run", &run);
  outevents->Branch("subrun", &subrun);
  outevents->Branch("event", &event);
  outevents->Branch("signal", &signal);
  outevents->Branch("signus", &signus);
  outevents->Branch("true_slc0", &true_slc0);
  outevents->Branch("true_slc1", &true_slc1);

  outevents->Branch("slc0", &slc0);
  outevents->Branch("vtx_x0", &vtx_x0);
  outevents->Branch("vtx_y0", &vtx_y0);
  outevents->Branch("vtx_z0", &vtx_z0);
  outevents->Branch("crumbs0", &crumbs0);
  outevents->Branch("crumbs_nc0", &crumbs_nc0);
  outevents->Branch("crumbs_ccnue0", &crumbs_ccnue0);
  outevents->Branch("n_dazzle_muons0", &n_dazzle_muons0);
  outevents->Branch("opT0Frac0", &opT0Frac0);
  outevents->Branch("goodOpT00", &goodOpT00);

  outevents->Branch("slc1", &slc1);
  outevents->Branch("vtx_x1", &vtx_x1);
  outevents->Branch("vtx_y1", &vtx_y1);
  outevents->Branch("vtx_z1", &vtx_z1);
  outevents->Branch("crumbs1", &crumbs1);
  outevents->Branch("crumbs_nc1", &crumbs_nc1);
  outevents->Branch("crumbs_ccnue1", &crumbs_ccnue1);
  outevents->Branch("n_dazzle_muons1", &n_dazzle_muons1);
  outevents->Branch("opT0Frac1", &opT0Frac1);
  outevents->Branch("goodOpT01", &goodOpT01);

  outevents->Branch("sep", &sep);
  outevents->Branch("sumOpT0Frac", &sumOpT0Frac);
  outevents->Branch("matching_flash_pe", &matching_flash_pe);
  outevents->Branch("signal_cand", &signal_cand);
  outevents->Branch("good_cand", &good_cand);

  for(size_t i = 0; i < events->GetEntries(); ++i)
    {
      events->GetEntry(i);

      slc0.clear(); vtx_x0.clear(); vtx_y0.clear(); vtx_z0.clear();
      crumbs0.clear(); crumbs_nc0.clear(); crumbs_ccnue0.clear();
      n_dazzle_muons0.clear(); opT0Frac0.clear(); goodOpT00.clear();

      slc1.clear(); vtx_x1.clear(); vtx_y1.clear(); vtx_z1.clear();
      crumbs1.clear(); crumbs_nc1.clear(); crumbs_ccnue1.clear();
      n_dazzle_muons1.clear(); opT0Frac1.clear(); goodOpT01.clear();

      sep.clear(); matching_flash_pe.clear(); sumOpT0Frac.clear();
      signal_cand.clear(); good_cand.clear();

      signus = 0;

      for(size_t j = 0; j < nu_event_type->size(); ++j)
        {
          if(nu_event_type->at(j) == 0)
            ++signus;
        }

      size_t sigslices = 0;
      double totalComp = 0., totalPur = 0.;

      true_slc0 = -1; true_slc1 = -1;

      for(size_t j = 0; j < slc_true_event_type->size(); ++j)
        {
          if(slc_is_clear_cosmic->at(j))
            continue;

          if(slc_true_event_type->at(j) == 0 && slc_comp->at(j) > .15)
            {
              ++sigslices;

              if(sigslices == 1)
                true_slc0 = j;
              else if(sigslices == 2)
                true_slc1 = j;

              totalComp += slc_comp->at(j);
              totalPur  += slc_pur->at(j);
            }
        }

      signal = totalComp > .8 && totalComp <= 1.0001 && sigslices == 2
        && (totalPur / 2) > .8;

      if(!signal)
        {
          true_slc0 = -1;
          true_slc1 = -1;
        }
    
      for(size_t j = 0; j < slc_true_event_type->size(); ++j)
        {
          if(slc_is_clear_cosmic->at(j))
            continue;

          TVector3 vertex0     = TVector3(slc_vtx_x->at(j), slc_vtx_y->at(j), slc_vtx_z->at(j));
          double hypPE0        = slc_opt0_hypPE->at(j);
          double measPE0       = slc_opt0_measPE->at(j);

          for(size_t k = j+1; k < slc_true_event_type->size(); ++k)
            {
              if(slc_is_clear_cosmic->at(k))
                continue;

              TVector3 vertex1     = TVector3(slc_vtx_x->at(k), slc_vtx_y->at(k), slc_vtx_z->at(k));
              double hypPE1        = slc_opt0_hypPE->at(k);
              double measPE1       = slc_opt0_measPE->at(k);

              slc0.push_back(j);
              vtx_x0.push_back(vertex0.X());
              vtx_y0.push_back(vertex0.Y());
              vtx_z0.push_back(vertex0.Z());
              crumbs0.push_back(slc_crumbs_score->at(j));
              crumbs_nc0.push_back(slc_crumbs_nc_score->at(j));
              crumbs_ccnue0.push_back(slc_crumbs_ccnue_score->at(j));
              n_dazzle_muons0.push_back(slc_n_dazzle_muons_cut_based->at(j));
              opT0Frac0.push_back((hypPE0 - measPE0) / measPE0);
              goodOpT00.push_back((hypPE0 > 0 && measPE0 > 0));

              slc1.push_back(k);
              vtx_x1.push_back(vertex1.X());
              vtx_y1.push_back(vertex1.Y());
              vtx_z1.push_back(vertex1.Z());
              crumbs1.push_back(slc_crumbs_score->at(k));
              crumbs_nc1.push_back(slc_crumbs_nc_score->at(k));
              crumbs_ccnue1.push_back(slc_crumbs_ccnue_score->at(k));
              n_dazzle_muons1.push_back(slc_n_dazzle_muons_cut_based->at(k));
              opT0Frac1.push_back((hypPE1 - measPE1) / measPE1);
              goodOpT01.push_back((hypPE1 > 0 && measPE1 > 0));

              sep.push_back((vertex0 - vertex1).Mag());
              matching_flash_pe.push_back(measPE0 == measPE1);
              
              if(measPE0 == measPE1)
                sumOpT0Frac.push_back((hypPE0 + hypPE1 - measPE0) / measPE0);
              else
                sumOpT0Frac.push_back((hypPE0 + hypPE1 - (measPE0+measPE1)) / (measPE0 + measPE1));

              signal_cand.push_back(signal && ((true_slc0 == j && true_slc1==k)
                                               || (true_slc0 == k && true_slc1 == k)));
              
              good_cand.push_back(slc_true_mctruth_id->at(j) == slc_true_mctruth_id->at(k)
                                  && slc_true_event_type->at(j) != 6
                                  && slc_true_event_type->at(k) != 6
                                  && (slc_comp->at(j) + slc_comp->at(k)) > 0.8
                                  && (slc_comp->at(j) + slc_comp->at(k)) < 1.0001);
            }
        }

      outevents->Fill();
    }

  outevents->Write();
}
