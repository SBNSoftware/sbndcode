void MultiSlices()
{
  const TString file = "/sbnd/data/users/hlay/pizero/ncpizeroana_sbnd.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add("/sbnd/data/users/hlay/pizero/ncpizeroana_sbnd.root");

  std::vector<int> *nu_event_type = 0, *slc_true_event_type = 0, *slc_n_trks = 0, *slc_n_shws = 0,
    *slc_n_dazzle_muons_cut_based = 0;
  std::vector<size_t> *slc_true_mctruth_id = 0;
  std::vector<double> *slc_comp = 0, *slc_pur = 0, *slc_vtx_x = 0, *slc_vtx_y = 0, *slc_vtx_z = 0;

  events->SetBranchAddress("nu_event_type", &nu_event_type);
  events->SetBranchAddress("slc_true_event_type", &slc_true_event_type);
  events->SetBranchAddress("slc_true_mctruth_id", &slc_true_mctruth_id);
  events->SetBranchAddress("slc_comp", &slc_comp);
  events->SetBranchAddress("slc_pur", &slc_pur);
  events->SetBranchAddress("slc_n_trks", &slc_n_trks);
  events->SetBranchAddress("slc_n_shws", &slc_n_shws);
  events->SetBranchAddress("slc_n_dazzle_muons_cut_based", &slc_n_dazzle_muons_cut_based);
  events->SetBranchAddress("slc_vtx_x", &slc_vtx_x);
  events->SetBranchAddress("slc_vtx_y", &slc_vtx_y);
  events->SetBranchAddress("slc_vtx_z", &slc_vtx_z);

  size_t recoverables = 0;

  for(size_t i = 0; i < events->GetEntries(); ++i)
    {
      events->GetEntry(i);

      size_t signus = 0;

      for(size_t j = 0; j < nu_event_type->size(); ++j)
        {
          if(nu_event_type->at(j) == 0)
            ++signus;
        }

      size_t sigslices = 0;
      double totalComp = 0., totalPur = 0.;
      std::vector<TVector3> vertices;

      for(size_t j = 0; j < slc_true_event_type->size(); ++j)
        {
          if(slc_true_event_type->at(j) == 0 && slc_comp->at(j) > .15)
            {
              ++sigslices;

              totalComp += slc_comp->at(j);
              totalPur  += slc_pur->at(j);
              vertices.emplace_back(slc_vtx_x->at(j), slc_vtx_y->at(j), slc_vtx_z->at(j));
            }
        }

      if(totalComp > .8 && totalComp <= 1 + 0.0001 && sigslices == 2 && (totalPur / 2) > .8)
        {
          ++recoverables;

          std::cout << "Recoverable event" << std::endl;
          if(signus > 1)
            std::cout << "WARNING - PILEUP SIGNAL" << std::endl;

          for(size_t j = 0; j < slc_true_event_type->size(); ++j)
            {
              if(slc_true_event_type->at(j) == 0)
                {
                  std::cout << "Slice with completeness " << slc_comp->at(j)
                            << " and purity " << slc_pur->at(j) << "\n\t"
                            << "Trks: " << slc_n_trks->at(j) << "\n\t"
                            << "Shws: " << slc_n_shws->at(j) << "\n\t"
                            << "DazzleMuons: " << slc_n_dazzle_muons_cut_based->at(j)
                            << std::endl;
                }
            }
          std::cout << "Separation: " << (vertices[0] - vertices[1]).Mag() << '\n'
                    << std::endl;
        }
    }

  std::cout << '\n' << recoverables << " potentially recoverable events\n" << std::endl;
}
