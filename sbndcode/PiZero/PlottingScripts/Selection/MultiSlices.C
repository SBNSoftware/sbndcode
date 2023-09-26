struct Candidate
{
  size_t id0;
  size_t id1;
  double sep;
  double opT0Frac0;
  double opT0Frac1;
  double sumOpT0Frac;
  double crumbs0;
  double crumbs1;
  int    ndazzle0;
  int    ndazzle1;
};

void MultiSlices()
{
  const TString file = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(file);

  int run, subrun, event;
  std::vector<bool> *slc_is_clear_cosmic = 0;
  std::vector<int> *nu_event_type = 0, *slc_true_event_type = 0, *slc_n_trks = 0, *slc_n_shws = 0,
    *slc_n_dazzle_muons_cut_based = 0;
  std::vector<size_t> *slc_true_mctruth_id = 0;
  std::vector<double> *slc_comp = 0, *slc_pur = 0, *slc_vtx_x = 0, *slc_vtx_y = 0, *slc_vtx_z = 0,
    *slc_opt0_measPE = 0, *slc_opt0_hypPE = 0, *slc_crumbs_score = 0;

  events->SetBranchAddress("run", &run);
  events->SetBranchAddress("subrun", &subrun);
  events->SetBranchAddress("event", &event);
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
  events->SetBranchAddress("slc_crumbs_ccnue_score", &slc_crumbs_score);

  size_t recoverables = 0, selected = 0, signalSelected = 0;

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

      size_t sigslices = 0, slc0 = -1, slc1 = -1;
      double totalComp = 0., totalPur = 0.;
      std::vector<TVector3> vertices;

      for(size_t j = 0; j < slc_true_event_type->size(); ++j)
        {
          if(slc_is_clear_cosmic->at(j))
            continue;

          if(slc_true_event_type->at(j) == 0 && slc_comp->at(j) > .15)
            {
              ++sigslices;

              if(sigslices == 1)
                slc0 = j;
              else if(sigslices == 2)
                slc1 = j;

              totalComp += slc_comp->at(j);
              totalPur  += slc_pur->at(j);
              vertices.emplace_back(slc_vtx_x->at(j), slc_vtx_y->at(j), slc_vtx_z->at(j));
            }
        }

      const bool signal    = totalComp > .8 && totalComp <= 1 + 0.0001 && sigslices == 2
        && (totalPur / 2) > .8;
      const size_t sigSlc0 = signal ? slc0 : -1;
      const size_t sigSlc1 = signal ? slc1 : -1;

      if(signal)
        {
          ++recoverables;

          if(signus > 1)
            std::cout << "WARNING - PILEUP SIGNAL" << std::endl;

          double hypPE1 = 0., hypPE2 = 0., measPE = 0;
          int nslc = 0;

          for(size_t j = 0; j < slc_true_event_type->size(); ++j)
            {
              if(slc_true_event_type->at(j) == 0 && slc_comp->at(j) > .15)
                {
                  if(nslc == 0)
                    {
                      hypPE1 = slc_opt0_hypPE->at(j);
                      measPE = slc_opt0_measPE->at(j);
                      std::cout << "CRUMBS1: " << slc_crumbs_score->at(j) << std::endl;
                    }
                  else if(nslc == 1)
                    {
                      if(measPE != slc_opt0_measPE->at(j))
                        std::cout << "Whoa... " << measPE << " vs. " 
                                  << slc_opt0_measPE->at(j) << '\n'
                                  << '\t' << vertices[0].X() << " vs. "
                                  << vertices[1].X() << std::endl;
                      hypPE2 = slc_opt0_hypPE->at(j);
                      std::cout << "CRUMBS2: " << slc_crumbs_score->at(j) << std::endl;
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
        }

      std::vector<Candidate> candidates;

      for(size_t j = 0; j < slc_true_event_type->size(); ++j)
        {
          if(slc_is_clear_cosmic->at(j))
            continue;

          TVector3 vertex0 = TVector3(slc_vtx_x->at(j), slc_vtx_y->at(j), slc_vtx_z->at(j));
          double hypPE0    = slc_opt0_hypPE->at(j);
          double measPE    = slc_opt0_measPE->at(j);
          double opT0Frac0 = (hypPE0 - measPE) / measPE;
          double crumbs0   = slc_crumbs_score->at(j);
          int ndazzle0     = slc_n_dazzle_muons_cut_based->at(j);

          for(size_t k = j+1; k < slc_true_event_type->size(); ++k)
            {
              if(slc_is_clear_cosmic->at(k))
                continue;

              TVector3 vertex1   = TVector3(slc_vtx_x->at(k), slc_vtx_y->at(k), slc_vtx_z->at(k));
              double hypPE1      = slc_opt0_hypPE->at(k);
              double opT0Frac1   = (hypPE1 - measPE) / measPE;
              double sumOpT0Frac = (hypPE0 + hypPE1 - measPE) / measPE;
              double sep         = (vertex0 - vertex1).Mag();
              double crumbs1   = slc_crumbs_score->at(k);
              int ndazzle1     = slc_n_dazzle_muons_cut_based->at(k);

              candidates.push_back({j, k, sep, opT0Frac0, opT0Frac1, sumOpT0Frac,
                    crumbs0, crumbs1, ndazzle0, ndazzle1});
            }
        }

      for(auto it = candidates.begin(); it != candidates.end();)
        {
          if(it->sep < 200 && it->opT0Frac0 < -0.2 && it->opT0Frac1 < -0.2
             && it->sumOpT0Frac < 0.0 && it->sumOpT0Frac > -0.6 && it->crumbs0 > -0.1
             && it->crumbs1 > -0.1)// && it->ndazzle0 == 0 && it->ndazzle1 == 0)
            ++it;
          else
            it = candidates.erase(it);
        }

      if(candidates.size() > 0)
        {
          ++selected;

          std::sort(candidates.begin(), candidates.end(),
                    [](const auto &a, const auto &b)
                    { return abs(a.sumOpT0Frac) < abs(b.sumOpT0Frac); });

          if(signal && ((sigSlc0 == candidates[0].id0 && sigSlc1 == candidates[0].id1)
                        || (sigSlc0 == candidates[0].id1 && sigSlc1 == candidates[0].id0)))
            ++signalSelected;
        }

    }

  std::cout << '\n' << recoverables << " potentially recoverable events" << std::endl;
  std::cout << '\t' << selected << " selected events" << std::endl;
  std::cout << '\t' << "of which " << signalSelected << " signal\n" << std::endl;


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
