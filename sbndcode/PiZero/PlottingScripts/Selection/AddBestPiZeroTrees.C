constexpr double def_double = std::numeric_limits<double>::lowest();

void AddBestPiZeroTrees()
{
  const TString infilename = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv4_1/NCPiZeroAv4_1_rockbox.root.bck";
  TFile *infile = TFile::Open(infilename, "update");

  TTree *events = (TTree*) infile->Get("ncpizeroana/events");
  TTree *subruns = (TTree*) infile->Get("ncpizeroana/subruns");

  // const TString outfilename = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv4_1/NCPiZeroAv4_1_rockbox_best_pzcs.root";
  // TFile *outfile = new TFile(outfilename, "recreate");
  // outfile->mkdir("ncpizeroana");
  // outfile->cd("ncpizeroana");

  // TTree *outevents = (TTree*) events->Clone();
  // TTree *outsubruns = (TTree*) subruns->Clone();

  std::vector<size_t> *slc_n_pzcs = 0;
  std::vector<std::vector<double>> *slc_pzc_invariant_mass = 0, *slc_pzc_pizero_mom = 0,
    *slc_pzc_cos_theta_pizero = 0, *slc_pzc_cos_com = 0, *slc_pzc_decay_asymmetry = 0;
  std::vector<double> *slc_best_pzc_invariant_mass = 0, *slc_best_pzc_pizero_mom = 0,
    *slc_best_pzc_cos_theta_pizero = 0, *slc_best_pzc_cos_com = 0, *slc_best_pzc_decay_asymmetry = 0;

  events->SetBranchAddress("slc_n_pzcs", &slc_n_pzcs);
  events->SetBranchAddress("slc_pzc_invariant_mass", &slc_pzc_invariant_mass);
  events->SetBranchAddress("slc_pzc_pizero_mom", &slc_pzc_pizero_mom);
  events->SetBranchAddress("slc_pzc_cos_theta_pizero", &slc_pzc_cos_theta_pizero);
  events->SetBranchAddress("slc_pzc_cos_com", &slc_pzc_cos_com);
  events->SetBranchAddress("slc_pzc_decay_asymmetry", &slc_pzc_decay_asymmetry);

  TBranch *b_best_pzc_invariant_mass = events->Branch("slc_best_pzc_invariant_mass", &slc_best_pzc_invariant_mass);
  TBranch *b_best_pzc_pizero_mom = events->Branch("slc_best_pzc_pizero_mom", &slc_best_pzc_pizero_mom);
  TBranch *b_best_pzc_cos_theta_pizero = events->Branch("slc_best_pzc_cos_theta_pizero", &slc_best_pzc_cos_theta_pizero);
  TBranch *b_best_pzc_cos_com = events->Branch("slc_best_pzc_cos_com", &slc_best_pzc_cos_com);
  TBranch *b_best_pzc_decay_asymmetry = events->Branch("slc_best_pzc_decay_asymmetry", &slc_best_pzc_decay_asymmetry);

  for(size_t i = 0; i < events->GetEntries(); ++i)
    {
      if(!i%1000) std::cout << i << " / " << events->GetEntries() << std::endl;
      events->GetEntry(i);

      slc_best_pzc_invariant_mass->clear(); slc_best_pzc_pizero_mom->clear();
      slc_best_pzc_cos_theta_pizero->clear(); slc_best_pzc_cos_com->clear(); slc_best_pzc_decay_asymmetry->clear();
      
      for(size_t j = 0; j < slc_n_pzcs->size(); ++j)
        {
          if(slc_pzc_invariant_mass->at(j).size() == 0)
            {
              slc_best_pzc_invariant_mass->push_back(def_double);
              slc_best_pzc_pizero_mom->push_back(def_double);
              slc_best_pzc_cos_theta_pizero->push_back(def_double);
              slc_best_pzc_cos_com->push_back(def_double);
              slc_best_pzc_decay_asymmetry->push_back(def_double);
            }
          else
            {
              double bestInvMass = def_double; int bestK = -1;
              for(size_t k = 0; k < slc_pzc_invariant_mass->at(j).size(); ++k)
                {
                  if(abs(slc_pzc_invariant_mass->at(j).at(k) - 139.570) < abs(bestInvMass))
                    {
                      bestInvMass = abs(slc_pzc_invariant_mass->at(j).at(k) - 139.570);
                      bestK = k;
                    }
                }

              slc_best_pzc_invariant_mass->push_back(slc_pzc_invariant_mass->at(j).at(bestK));
              slc_best_pzc_pizero_mom->push_back(slc_pzc_pizero_mom->at(j).at(bestK));
              slc_best_pzc_cos_theta_pizero->push_back(slc_pzc_cos_theta_pizero->at(j).at(bestK));
              slc_best_pzc_cos_com->push_back(slc_pzc_cos_com->at(j).at(bestK));
              slc_best_pzc_decay_asymmetry->push_back(slc_pzc_decay_asymmetry->at(j).at(bestK));
            }
        }

      b_best_pzc_invariant_mass->Fill();
      b_best_pzc_pizero_mom->Fill();
      b_best_pzc_cos_theta_pizero->Fill();
      b_best_pzc_cos_com->Fill();
      b_best_pzc_decay_asymmetry->Fill();
    }

  events->Write();
  subruns->Write();
}
