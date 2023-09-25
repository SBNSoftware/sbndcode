void MultiSlices()
{
  const TString file = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(file);

  std::vector<bool> *slc_is_clear_cosmic = 0;
  std::vector<int> *nu_event_type = 0, *slc_true_event_type = 0, *slc_n_trks = 0, *slc_n_shws = 0,
    *slc_n_dazzle_muons_cut_based = 0;
  std::vector<size_t> *slc_true_mctruth_id = 0;
  std::vector<double> *slc_comp = 0, *slc_pur = 0, *slc_vtx_x = 0, *slc_vtx_y = 0, *slc_vtx_z = 0,
    *slc_opt0_measPE = 0, *slc_opt0_hypPE = 0;

  events->SetBranchAddress("nu_event_type", &nu_event_type);
  events->SetBranchAddress("slc_is_clear_cosmic", &slc_is_clear_cosmic);
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
  events->SetBranchAddress("slc_opt0_measPE", &slc_opt0_measPE);
  events->SetBranchAddress("slc_opt0_hypPE", &slc_opt0_hypPE);

  size_t recoverables = 0;

  TH1F *hSumOpT0Frac = new TH1F("hSumOpT0Frac", ";#frac{#Sigma PE_{hyp} - PE_{meas}}{PE_{meas}};Events", 40, -1, 1);

  TH2F *hOpT0Fracs = new TH2F("hOpT0Fracs", ";#frac{PE_{hyp, 1} - PE_{meas}}{PE_{meas}};;#frac{PE_{hyp, 2} - PE_{meas}}{PE_{meas}}", 40, -1, 1, 40, -1, 1);

  TH2F *hSumOpT0FracSep = new TH2F("hOpT0FracsSep", ";#frac{#Sigma PE_{hyp} - PE_{meas}}{PE_{meas}};Separation (cm)", 40, -1, 1, 60, 0, 600);

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
          if(slc_is_clear_cosmic->at(j))
            continue;

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

          double hypPE1 = 0., hypPE2 = 0., measPE = 0;
          int nslc = 0;

          for(size_t j = 0; j < slc_true_event_type->size(); ++j)
            {
              if(slc_true_event_type->at(j) == 0 && slc_comp->at(j) > .15)
                {
                  std::cout << "Slice with completeness " << slc_comp->at(j)
                            << " and purity " << slc_pur->at(j) << "\n\t"
                            << "Trks: " << slc_n_trks->at(j) << "\n\t"
                            << "Shws: " << slc_n_shws->at(j) << "\n\t"
                            << "DazzleMuons: " << slc_n_dazzle_muons_cut_based->at(j) << "\n\t"
                            << "OpT0Frac: "
                            << (slc_opt0_hypPE->at(j) - slc_opt0_measPE->at(j)) / slc_opt0_measPE->at(j)
                            << std::endl;

                  if(nslc == 0)
                    {
                      hypPE1 = slc_opt0_hypPE->at(j);
                      measPE = slc_opt0_measPE->at(j);
                    }
                  else if(nslc == 1)
                    {
                      if(measPE != slc_opt0_measPE->at(j))
                        std::cout << "Whoa... " << measPE << " vs. " 
                                  << slc_opt0_measPE->at(j) << std::endl;
                      hypPE2 = slc_opt0_hypPE->at(j);
                    }
                  else
                    std::cout << "Whoa... " << nslc << std::endl;

                  ++nslc;
                }
            }
          double sumOpT0Frac = (hypPE1 + hypPE2 - measPE) / measPE;
          double opT0Frac1   = (hypPE1 - measPE) / measPE;
          double opT0Frac2   = (hypPE2 - measPE) / measPE;
          double sep         = (vertices[0] - vertices[1]).Mag();

          hSumOpT0Frac->Fill(sumOpT0Frac);
          hOpT0Fracs->Fill(opT0Frac1, opT0Frac2);
          hSumOpT0FracSep->Fill(sumOpT0Frac, sep);

          std::cout << "Separation: " << sep << '\n' << std::endl;
        }
    }

  std::cout << '\n' << recoverables << " potentially recoverable events\n" << std::endl;

  TCanvas *cSumOpT0Frac = new TCanvas("cSumOpT0Frac", "cSumOpT0Frac");
  cSumOpT0Frac->cd();

  hSumOpT0Frac->Draw();

  TCanvas *cOpT0Fracs = new TCanvas("cOpT0Fracs", "cOpT0Fracs");
  cOpT0Fracs->cd();

  hOpT0Fracs->Draw("colz");

  TCanvas *cSumOpT0FracSep = new TCanvas("cSumOpT0FracSep", "cSumOpT0FracSep");
  cSumOpT0FracSep->cd();

  hSumOpT0FracSep->Draw("colz");
}
