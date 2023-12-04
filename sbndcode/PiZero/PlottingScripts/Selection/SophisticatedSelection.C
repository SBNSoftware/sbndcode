const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

std::vector<int> *nu_event_type_incl = 0, *nu_event_type_0p0pi, *nu_event_type_Xp0pi = 0, *nu_event_type_1p0pi = 0,
  *nu_event_type_Np0pi = 0, *nu_pdg = 0, *nu_ccnc = 0, *nu_mode = 0, *nu_n_protons, *nu_n_charged_pions = 0,
  *nu_n_neutral_pions = 0, *nu_n_dalitz_neutral_pions = 0, *slc_n_trks = 0, *slc_n_shws = 0, *slc_n_dazzle_muons = 0,
  *slc_n_dazzle_pions = 0, *slc_n_dazzle_pions_thresh = 0, *slc_n_dazzle_protons = 0, *slc_n_dazzle_protons_thresh = 0,
  *slc_n_dazzle_other = 0, *slc_n_razzle_electrons = 0, *slc_n_razzle_photons = 0, *slc_n_razzle_other = 0,
  *slc_n_razzled_electrons = 0, *slc_n_razzled_muons = 0, *slc_n_razzled_photons = 0, *slc_n_razzled_pions = 0,
  *slc_n_razzled_pions_thresh = 0, *slc_n_razzled_protons = 0, *slc_n_razzled_protons_thresh = 0, *slc_true_event_type_incl = 0,
  *slc_true_event_type_0p0pi = 0, *slc_true_event_type_Xp0pi = 0, *slc_true_event_type_1p0pi = 0, *slc_true_event_type_Np0pi = 0,
  *slc_true_pdg = 0, *slc_true_ccnc = 0, *slc_true_mode = 0, *slc_true_n_protons = 0, *slc_true_n_charged_pions = 0, 
  *slc_true_n_neutral_pions = 0, *slc_true_n_dalitz_neutral_pions = 0, *slc_best_pzc_photon_0_id = 0, *slc_best_pzc_photon_1_id = 0;

std::vector<bool> *nu_signal = 0, *nu_av  = 0, *nu_fv = 0, *slc_is_clear_cosmic = 0, *slc_true_signal = 0,
  *slc_true_av = 0, *slc_true_fv = 0, *slc_is_fv = 0, *slc_best_pzc_good_kinematics = 0;
std::vector<size_t> *slc_n_pfps = 0;
std::vector<float> *slc_comp = 0, *slc_pur = 0, *slc_crumbs_score = 0;

std::vector<std::vector<double>> *slc_pfp_shower_length = 0, *slc_pfp_shower_energy = 0, *slc_pfp_shower_dir_x = 0, *slc_pfp_shower_dir_y = 0,
  *slc_pfp_shower_dir_z = 0;
std::vector<std::vector<int>> *slc_pfp_razzled_pdg = 0;
std::vector<std::vector<bool>> *slc_pzc_good_kinematics = 0;

std::vector<float> true_signal(5, 0);

std::vector<std::vector<float>> total_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> presel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> photonsel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> elecsel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> invmasssel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> pionsel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> protonsel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> protonantisel_slc(5, std::vector<float>(10,0));
std::vector<std::vector<float>> oneprotonsel_slc(5, std::vector<float>(10,0));

void InitialiseTree(TChain *tree);

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void ProcessTree(TChain* chain, const float weight);

void Categorise(int slc_true_event_type_incl, int slc_true_event_type_0p0pi, int slc_true_event_type_Xp0pi, int slc_true_event_type_1p0pi,
                int slc_true_event_type_Np0pi, float slc_comp, std::vector<std::vector<float>> &counts, const float weight);

void Categorise(int slc_true_event_type, float slc_comp, std::vector<std::vector<float>> &counts, const float weight, const int sel_type);

void Evaluate(const TString name, std::vector<float> &counts, const float total_signal, const float total_true_signal);

void SophisticatedSelection(const TString productionVersion)
{
  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxsubruns = new TChain("ncpizeroana/subruns");
  rockboxsubruns->Add(rockboxFile);
  TChain *intimesubruns = new TChain("ncpizeroana/subruns");
  intimesubruns->Add(intimeFile);

  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+","x10^{");
  potString.ReplaceAll(" POT","} POT");

  const double rockboxPOT = GetPOT(rockboxsubruns);
  const int rockboxSpills = GetGenEvents(rockboxsubruns);
  const int intimeSpills  = GetGenEvents(intimesubruns);

  const double rockboxScaling      = goalPOT / rockboxPOT;
  const double scaledRockboxSpills = rockboxScaling * rockboxSpills;
  const double intimeScaling       = (goalSpills - scaledRockboxSpills) / intimeSpills;

  InitialiseTree(rockboxEvents);
  ProcessTree(rockboxEvents, 1.);

  InitialiseTree(intimeEvents);
  ProcessTree(intimeEvents, intimeScaling / rockboxScaling);
  
  std::cout << "Total POT: " << rockboxPOT << " (" << Form("%.2f times smaller than 1e21)", rockboxScaling) << '\n'
            << "...also " << intimeEvents->GetEntries() << " intime cosmic triggers (" << Form("%.2f times smaller than 1e21)", intimeScaling) << '\n'
            << "NC1PiZero events: " << true_signal[0] << '\n'
            << "\t\t 0p0pi: " << true_signal[1] << '\n'
            << "\t\t Xp0pi: " << true_signal[2] << '\n'
            << "\t\t Np0pi: " << true_signal[3] << '\n'
            << "\t\t 1p0pi: " << true_signal[4] << '\n' << std::endl;

  Evaluate("Signal One PiZero - No Sel", total_slc[0], total_slc[0][0], true_signal[0]);
  Evaluate("Signal One PiZero - Standard", photonsel_slc[0], total_slc[0][0], true_signal[0]);
  Evaluate("Signal One PiZero - Allow Electron", elecsel_slc[0], total_slc[0][0], true_signal[0]);
  Evaluate("Signal One PiZero - InvMass", invmasssel_slc[0], total_slc[0][0], true_signal[0]);
  Evaluate("Signal 0p0pi", protonantisel_slc[1], total_slc[1][0], true_signal[1]);
  Evaluate("Signal Xp0pi", pionsel_slc[2], total_slc[2][0], true_signal[2]);
  Evaluate("Signal 1p0pi", oneprotonsel_slc[3], total_slc[3][0], true_signal[3]);
  Evaluate("Signal Np0pi", protonsel_slc[4], total_slc[4][0], true_signal[4]);

}

void InitialiseTree(TChain *tree)
{
  tree->SetBranchStatus("*", 0);

  tree->SetBranchAddress("nu_event_type_incl", &nu_event_type_incl);
  tree->SetBranchAddress("nu_event_type_0p0pi", &nu_event_type_0p0pi);
  tree->SetBranchAddress("nu_event_type_Xp0pi", &nu_event_type_Xp0pi);
  tree->SetBranchAddress("nu_event_type_1p0pi", &nu_event_type_1p0pi);
  tree->SetBranchAddress("nu_event_type_Np0pi", &nu_event_type_Np0pi);
  tree->SetBranchAddress("nu_pdg", &nu_pdg);
  tree->SetBranchAddress("nu_ccnc", &nu_ccnc);
  tree->SetBranchAddress("nu_mode", &nu_mode);
  tree->SetBranchAddress("nu_n_protons", &nu_n_protons);
  tree->SetBranchAddress("nu_n_charged_pions", &nu_n_charged_pions);
  tree->SetBranchAddress("nu_n_neutral_pions", &nu_n_neutral_pions);
  tree->SetBranchAddress("nu_n_dalitz_neutral_pions", &nu_n_dalitz_neutral_pions);
  tree->SetBranchAddress("nu_signal", &nu_signal);
  tree->SetBranchAddress("nu_av", &nu_av);
  tree->SetBranchAddress("nu_fv", &nu_fv);
  tree->SetBranchAddress("slc_n_trks", &slc_n_trks);
  tree->SetBranchAddress("slc_n_shws", &slc_n_shws);
  tree->SetBranchAddress("slc_n_dazzle_muons", &slc_n_dazzle_muons);
  tree->SetBranchAddress("slc_n_dazzle_pions", &slc_n_dazzle_pions);
  tree->SetBranchAddress("slc_n_dazzle_pions_thresh", &slc_n_dazzle_pions_thresh);
  tree->SetBranchAddress("slc_n_dazzle_protons", &slc_n_dazzle_protons);
  tree->SetBranchAddress("slc_n_dazzle_protons_thresh", &slc_n_dazzle_protons_thresh);
  tree->SetBranchAddress("slc_n_razzle_electrons", &slc_n_razzle_electrons);
  tree->SetBranchAddress("slc_n_razzle_photons", &slc_n_razzle_photons);
  tree->SetBranchAddress("slc_n_razzle_other", &slc_n_razzle_other);
  tree->SetBranchAddress("slc_n_razzled_electrons", &slc_n_razzled_electrons);
  tree->SetBranchAddress("slc_n_razzled_muons", &slc_n_razzled_muons);
  tree->SetBranchAddress("slc_n_razzled_photons", &slc_n_razzled_photons);
  tree->SetBranchAddress("slc_n_razzled_pions", &slc_n_razzled_pions);
  tree->SetBranchAddress("slc_n_razzled_pions_thresh", &slc_n_razzled_pions_thresh);
  tree->SetBranchAddress("slc_n_razzled_protons", &slc_n_razzled_protons);
  tree->SetBranchAddress("slc_n_razzled_protons_thresh", &slc_n_razzled_protons_thresh);
  tree->SetBranchAddress("slc_true_event_type_incl", &slc_true_event_type_incl);
  tree->SetBranchAddress("slc_true_event_type_0p0pi", &slc_true_event_type_0p0pi);
  tree->SetBranchAddress("slc_true_event_type_Xp0pi", &slc_true_event_type_Xp0pi);
  tree->SetBranchAddress("slc_true_event_type_1p0pi", &slc_true_event_type_1p0pi);
  tree->SetBranchAddress("slc_true_event_type_Np0pi", &slc_true_event_type_Np0pi);
  tree->SetBranchAddress("slc_true_pdg", &slc_true_pdg);
  tree->SetBranchAddress("slc_true_ccnc", &slc_true_ccnc);
  tree->SetBranchAddress("slc_true_mode", &slc_true_mode);
  tree->SetBranchAddress("slc_true_n_protons", &slc_true_n_protons);
  tree->SetBranchAddress("slc_true_n_charged_pions", &slc_true_n_charged_pions);
  tree->SetBranchAddress("slc_true_n_neutral_pions", &slc_true_n_neutral_pions);
  tree->SetBranchAddress("slc_true_n_dalitz_neutral_pions", &slc_true_n_dalitz_neutral_pions);
  tree->SetBranchAddress("slc_is_clear_cosmic", &slc_is_clear_cosmic);
  tree->SetBranchAddress("slc_true_signal", &slc_true_signal);
  tree->SetBranchAddress("slc_true_av", &slc_true_av);
  tree->SetBranchAddress("slc_true_fv", &slc_true_fv);
  tree->SetBranchAddress("slc_is_fv", &slc_is_fv);
  tree->SetBranchAddress("slc_n_pfps", &slc_n_pfps);
  tree->SetBranchAddress("slc_comp", &slc_comp);
  tree->SetBranchAddress("slc_pur", &slc_pur);
  tree->SetBranchAddress("slc_crumbs_score", &slc_crumbs_score);
  tree->SetBranchAddress("slc_best_pzc_good_kinematics", &slc_best_pzc_good_kinematics);
  tree->SetBranchAddress("slc_best_pzc_photon_0_id", &slc_best_pzc_photon_0_id);
  tree->SetBranchAddress("slc_best_pzc_photon_1_id", &slc_best_pzc_photon_1_id);
  tree->SetBranchAddress("slc_pfp_shower_length", &slc_pfp_shower_length);
  tree->SetBranchAddress("slc_pzc_good_kinematics", &slc_pzc_good_kinematics);
  tree->SetBranchAddress("slc_pfp_razzled_pdg", &slc_pfp_razzled_pdg);
  tree->SetBranchAddress("slc_pfp_shower_energy", &slc_pfp_shower_energy);
  tree->SetBranchAddress("slc_pfp_shower_dir_x", &slc_pfp_shower_dir_x);
  tree->SetBranchAddress("slc_pfp_shower_dir_y", &slc_pfp_shower_dir_y);
  tree->SetBranchAddress("slc_pfp_shower_dir_z", &slc_pfp_shower_dir_z);
}

double GetPOT(TChain *subruns)
{
  double sum = 0., pot = 0;

  subruns->SetBranchAddress("pot", &pot);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += pot;
    }

  return sum;
}

int GetGenEvents(TChain *subruns)
{
  int sum = 0., ngenevts = 0;

  subruns->SetBranchAddress("ngenevts", &ngenevts);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += ngenevts;
    }

  return sum;
}

void ProcessTree(TChain* chain, const float weight)
{
  const int N = chain->GetEntries();

  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      chain->GetEntry(ev_i);
      
      for(int nu_i = 0; nu_i < nu_event_type_incl->size(); ++nu_i)
        {
          if(nu_event_type_incl->at(nu_i) == 0)
            true_signal[0] += weight;
          if(nu_event_type_0p0pi->at(nu_i) == 0)
            true_signal[1] += weight;
          if(nu_event_type_Xp0pi->at(nu_i) == 0)
            true_signal[2] += weight;
          if(nu_event_type_1p0pi->at(nu_i) == 0)
            true_signal[3] += weight;
          if(nu_event_type_Np0pi->at(nu_i) == 0)
            true_signal[4] += weight;
        }

      for(int slc_i = 0; slc_i < slc_true_event_type_incl->size(); ++slc_i)
        {
          Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                     slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), total_slc, weight);

          if(!slc_is_clear_cosmic->at(slc_i) && slc_is_fv->at(slc_i) && slc_crumbs_score->at(slc_i) > -0.025 && slc_n_razzled_muons->at(slc_i) == 0)
            {
              Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                         slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), presel_slc, weight);

              if((slc_n_razzled_electrons->at(slc_i) + slc_n_razzled_photons->at(slc_i)) > 1 && slc_n_razzled_photons->at(slc_i) > 0)
                {
                  Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                             slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), elecsel_slc, weight);

                  double bestInvMassDiff = std::numeric_limits<double>::max();

                  int best_i = -1, best_j = -1;

                  for(int pfp_i = 0; pfp_i < slc_pfp_razzled_pdg->at(slc_i).size(); ++pfp_i)
                    {
                      if(slc_pfp_razzled_pdg->at(slc_i).at(pfp_i) != 22)
                        continue;

                      double e1 = slc_pfp_shower_energy->at(slc_i).at(pfp_i);
                      TVector3 dir1(slc_pfp_shower_dir_x->at(slc_i).at(pfp_i), slc_pfp_shower_dir_y->at(slc_i).at(pfp_i), slc_pfp_shower_dir_z->at(slc_i).at(pfp_i));

                      for(int pfp_j = 0; pfp_j < slc_pfp_razzled_pdg->at(slc_i).size(); ++pfp_j)
                        {
                          if(pfp_i == pfp_j)
                            continue;

                          if(slc_pfp_razzled_pdg->at(slc_i).at(pfp_i) != 22 && slc_pfp_razzled_pdg->at(slc_i).at(pfp_i) != 11)
                            continue;

                          double e2 = slc_pfp_shower_energy->at(slc_i).at(pfp_j);
                          TVector3 dir2(slc_pfp_shower_dir_x->at(slc_i).at(pfp_j), slc_pfp_shower_dir_y->at(slc_i).at(pfp_j), slc_pfp_shower_dir_z->at(slc_i).at(pfp_j));

                          const double costheta = dir1.Dot(dir2) / (dir1.Mag() * dir2.Mag());
                          const double invMass = sqrt(2 * e1 * e2 * (1 - costheta));

                          if(abs(134.9769 - invMass) < bestInvMassDiff)
                            {
                              bestInvMassDiff = abs(134.9769 - invMass);
                              best_i = pfp_i;
                              best_j = pfp_j;
                            }
                        }
                    }

                  if(best_i != -1 && best_j != -1 && bestInvMassDiff < 25.)
                    {
                      Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                                 slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), invmasssel_slc, weight);
                    }
                }

              bool anygood = false;
              for(auto const& good : slc_pzc_good_kinematics->at(slc_i))
                anygood |= good;

              if(slc_n_razzled_photons->at(slc_i) > 1 && anygood)
                {
                  Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                             slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), photonsel_slc, weight);

                  if(slc_n_razzled_pions_thresh->at(slc_i) == 0)
                    {
                      Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                                 slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), pionsel_slc, weight);

                      if(slc_n_razzled_protons_thresh->at(slc_i) > 0)
                        {
                          Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                                     slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), protonsel_slc, weight);
                        }
                      else
                        {
                          Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                                     slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), protonantisel_slc, weight);
                        }

                      if(slc_n_razzled_protons_thresh->at(slc_i) == 1)
                        {
                          Categorise(slc_true_event_type_incl->at(slc_i), slc_true_event_type_0p0pi->at(slc_i), slc_true_event_type_Xp0pi->at(slc_i), slc_true_event_type_1p0pi->at(slc_i),
                                     slc_true_event_type_Np0pi->at(slc_i), slc_comp->at(slc_i), oneprotonsel_slc, weight);
                        }

                    }
                }
            }
        }
    }
}

void Categorise(int slc_true_event_type_incl, int slc_true_event_type_0p0pi, int slc_true_event_type_Xp0pi, int slc_true_event_type_1p0pi,
                int slc_true_event_type_Np0pi, float slc_comp, std::vector<std::vector<float>> &counts, const float weight)
{
  Categorise(slc_true_event_type_incl, slc_comp, counts, weight, 0);
  Categorise(slc_true_event_type_0p0pi, slc_comp, counts, weight, 1);
  Categorise(slc_true_event_type_Xp0pi, slc_comp, counts, weight, 2);
  Categorise(slc_true_event_type_1p0pi, slc_comp, counts, weight, 3);
  Categorise(slc_true_event_type_Np0pi, slc_comp, counts, weight, 4);
}

void Categorise(int slc_true_event_type, float slc_comp, std::vector<std::vector<float>> &counts, const float weight, const int sel_type)
{
  if(slc_true_event_type == 0 && slc_comp < .5)
    slc_true_event_type = 9;

  if(slc_true_event_type < 0)
    slc_true_event_type = 8;

  counts[sel_type][slc_true_event_type] += weight;
}

void Evaluate(const TString name, std::vector<float> &counts, const float total_signal, const float total_signal_events)
{
  const float total = std::accumulate(counts.begin(), counts.end(), 0);

  std::cout << "\nEvaluating... " << name << '\n'
            << "TOTAL:           " << total << '\n'
            << "Signal:          " << counts[0] << Form(" (%.2f%%)", 100. * counts[0] / total) << '\n'
            << "Other NCPiZero:  " << counts[1] << Form(" (%.2f%%)", 100. * counts[1] / total) << '\n'
            << "Other NC:        " << counts[2] << Form(" (%.2f%%)", 100. * counts[2] / total) << '\n'
            << "CCNuMu:          " << counts[3] << Form(" (%.2f%%)", 100. * counts[3] / total) << '\n'
            << "CCNuE:           " << counts[4] << Form(" (%.2f%%)", 100. * counts[4] / total) << '\n'
            << "Dirt:            " << counts[5] << Form(" (%.2f%%)", 100. * counts[5] / total) << '\n'
            << "NonFV:           " << counts[6] << Form(" (%.2f%%)", 100. * counts[6] / total) << '\n'
            << "Cosmic:          " << counts[7] << Form(" (%.2f%%)", 100. * counts[7] / total) << '\n'
            << "No Truth Match:  " << counts[8] << Form(" (%.2f%%)", 100. * counts[8] / total) << '\n'
            << "Bad Reco Signal: " << counts[9] << Form(" (%.2f%%)", 100. * counts[9] / total) << '\n'
            << "EFF:             " << Form("%.2f%%", 100. * counts[0] / total_signal_events) << '\n'
            << "EFF * PUR:       " << Form("%.2f%%", 100. * counts[0] / total * counts[0] / total_signal_events) << '\n'
            << "SELEFF:          " << Form("%.2f%%", 100. * counts[0] / total_signal) << '\n'
            << "SELEFF * PUR:    " << Form("%.2f%%", 100. * counts[0] / total * counts[0] / total_signal) << '\n'
            << std::endl;
}
