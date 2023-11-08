const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

std::vector<int> *nu_event_type = 0, *nu_pdg = 0, *nu_ccnc = 0, *nu_mode = 0,
  *nu_n_protons, *nu_n_charged_pions = 0, *nu_n_neutral_pions = 0, *nu_n_dalitz_neutral_pions = 0,
  *slc_n_trks = 0, *slc_n_shws = 0, *slc_n_dazzle_muons = 0, *slc_n_dazzle_pions = 0, *slc_n_dazzle_pions_thresh = 0,
  *slc_n_dazzle_protons = 0, *slc_n_dazzle_protons_thresh = 0, *slc_n_dazzle_other = 0, *slc_n_razzle_electrons = 0,
  *slc_n_razzle_photons = 0, *slc_n_razzle_other = 0, *slc_n_razzled_electrons = 0, *slc_n_razzled_muons = 0,
  *slc_n_razzled_photons = 0, *slc_n_razzled_pions = 0, *slc_n_razzled_pions_thresh = 0, *slc_n_razzled_protons = 0,
  *slc_n_razzled_protons_thresh = 0, *slc_true_event_type = 0, *slc_true_pdg = 0, *slc_true_ccnc = 0, *slc_true_mode = 0,
  *slc_true_n_protons = 0, *slc_true_n_charged_pions = 0, *slc_true_n_neutral_pions = 0, *slc_true_n_dalitz_neutral_pions = 0,
  *slc_best_pzc_photon_0_id = 0, *slc_best_pzc_photon_1_id = 0;

std::vector<bool> *nu_signal = 0, *nu_av  = 0, *nu_fv = 0, *slc_is_clear_cosmic = 0, *slc_true_signal = 0,
  *slc_true_av = 0, *slc_true_fv = 0, *slc_is_fv = 0, *slc_best_pzc_good_kinematics = 0;
std::vector<size_t> *slc_n_pfps = 0;
std::vector<float> *slc_comp = 0, *slc_pur = 0, *slc_crumbs_score = 0;

std::vector<std::vector<double>> *slc_pfp_shower_length = 0;
std::vector<std::vector<bool>> *slc_pzc_good_kinematics = 0;

std::vector<float> true_signal(6, 0);

std::vector<std::vector<float>> total_slc(6, std::vector<float>(9,0));
std::vector<std::vector<float>> presel_slc(6, std::vector<float>(9,0));
std::vector<std::vector<float>> photonsel_slc(6, std::vector<float>(9,0));
std::vector<std::vector<float>> pionsel_slc(6, std::vector<float>(9,0));
std::vector<std::vector<float>> protonsel_slc(6, std::vector<float>(9,0));
std::vector<std::vector<float>> protonantisel_slc(6, std::vector<float>(9,0));
std::vector<std::vector<float>> oneprotonsel_slc(6, std::vector<float>(9,0));

void InitialiseTree(TChain *tree);

double GetPOT(TChain *subruns);
int GetGenEvents(TChain *subruns);

void ProcessTree(TChain* chain, const float weight);

void Categorise(int slc_true_event_type, int slc_true_n_neutral_pions, int slc_true_n_charged_pions,
                int slc_true_n_protons, float slc_comp, std::vector<std::vector<float>> &counts, const float weight);

void Evaluate(const TString name, std::vector<float> &counts, const float total_signal);

void SophisticatedSelection(const TString productionVersion, const TString saveDirExt)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/selection/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

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
            << "NCPiZero events: " << true_signal[0] << '\n'
            << "\t 1pi0: " << true_signal[1] << '\n'
            << "\t\t 0p0pi: " << true_signal[2] << '\n'
            << "\t\t Xp0pi: " << true_signal[3] << '\n'
            << "\t\t Np0pi: " << true_signal[4] << '\n'
            << "\t\t 1p0pi: " << true_signal[5] << '\n' << std::endl;

  Evaluate("Signal Inclusive - No Sel", total_slc[0], total_slc[0][0]);  
  Evaluate("Signal One PiZero", photonsel_slc[1], total_slc[1][0]);  
  Evaluate("Signal Xp0pi", pionsel_slc[3], total_slc[3][0]);
  Evaluate("Signal 0p0pi", protonantisel_slc[2], total_slc[2][0]);
  Evaluate("Signal Np0pi", protonsel_slc[4], total_slc[4][0]);
  Evaluate("Signal 1p0pi", oneprotonsel_slc[5], total_slc[5][0]);
}

void InitialiseTree(TChain *tree)
{
  tree->SetBranchStatus("*", 0);

  tree->SetBranchAddress("nu_event_type", &nu_event_type);
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
  tree->SetBranchAddress("slc_true_event_type", &slc_true_event_type);
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
      
      for(int nu_i = 0; nu_i < nu_event_type->size(); ++nu_i)
        {
          if(nu_event_type->at(nu_i) == 0)
            {
              true_signal[0] += weight;

              if(nu_n_neutral_pions->at(nu_i) == 1)
                {
                  true_signal[1] += weight;

                  if(nu_n_charged_pions->at(nu_i) == 0)
                    {
                      true_signal[3] += weight;

                      if(nu_n_protons->at(nu_i) == 0)
                        true_signal[2] += weight;
                      else
                        true_signal[4] += weight;
                      if(nu_n_protons->at(nu_i) == 1)
                        true_signal[5] += weight;

                    }
                }
            }
        }

      for(int slc_i = 0; slc_i < slc_true_event_type->size(); ++slc_i)
        {
          Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                     slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), total_slc, weight);

          if(!slc_is_clear_cosmic->at(slc_i) && slc_is_fv->at(slc_i) && slc_crumbs_score->at(slc_i) > -0.025 && slc_n_razzled_muons->at(slc_i) == 0)
            {
              Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                         slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), presel_slc, weight);

              bool anygood = false;
              for(auto const& good : slc_pzc_good_kinematics->at(slc_i))
                anygood |= good;

              if(slc_n_razzled_photons->at(slc_i) > 1 && anygood)
                {
                  Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                             slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), photonsel_slc, weight);

                  if(slc_n_razzled_pions_thresh->at(slc_i) == 0)
                    {
                      Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                                 slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), pionsel_slc, weight);

                      if(slc_n_razzled_protons_thresh->at(slc_i) > 0)
                        {
                          Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                                     slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), protonsel_slc, weight);
                        }
                      else
                        {
                          Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                                     slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), protonantisel_slc, weight);
                        }

                      if(slc_n_razzled_protons_thresh->at(slc_i) == 1)
                        {
                          Categorise(slc_true_event_type->at(slc_i), slc_true_n_neutral_pions->at(slc_i), slc_true_n_charged_pions->at(slc_i),
                                     slc_true_n_protons->at(slc_i), slc_comp->at(slc_i), oneprotonsel_slc, weight);
                        }

                    }
                }
            }
        }
    }
}

void Categorise(int slc_true_event_type, int slc_true_n_neutral_pions, int slc_true_n_charged_pions,
                int slc_true_n_protons, float slc_comp, std::vector<std::vector<float>> &counts, const float weight)
{
  if(slc_true_event_type == 0)
    {
      if(slc_comp > .5)
        counts[0][0] += weight;
      else
        counts[0][8] += weight;

      if(slc_true_n_neutral_pions == 1)
        {
          if(slc_comp > .5)
            counts[1][0] += weight;
          else
            counts[1][8] += weight;

          if(slc_true_n_charged_pions == 0)
            {
              if(slc_comp > .5)
                counts[3][0] += weight;
              else
                counts[3][8] += weight;

              if(slc_true_n_protons == 0)
                {
                  if(slc_comp > .5)
                    counts[2][0] += weight;
                  else
                    counts[2][8] += weight;
                }
              else
                counts[2][1] += weight;

              if(slc_true_n_protons != 0)
                {
                  if(slc_comp > .5)
                    counts[4][0] += weight;
                  else
                    counts[4][8] += weight;
                }
              else
                counts[4][1] += weight;

              if(slc_true_n_protons == 1)
                {
                  if(slc_comp > .5)
                    counts[5][0] += weight;
                  else
                    counts[5][8] += weight;
                }
              else
                counts[5][1] += weight;
            }
          else
            {
              counts[2][1] += weight;
              counts[3][1] += weight;
              counts[4][1] += weight;
              counts[5][1] += weight;
            }
        }
      else
        {
          counts[1][1] += weight;
          counts[2][1] += weight;
          counts[3][1] += weight;
          counts[4][1] += weight;
          counts[5][1] += weight;
        }
    }
  else if(slc_true_event_type < 0)
    for(int i = 0; i < 6; ++i)
      counts[i][7] += weight;
  else
    for(int i = 0; i < 6; ++i)
      counts[i][slc_true_event_type] += weight;
}

void Evaluate(const TString name, std::vector<float> &counts, const float total_signal)
{
  const float total = std::accumulate(counts.begin(), counts.end(), 0);

  std::cout << "\nEvaluating... " << name << '\n'
            << "TOTAL:           " << total << '\n'
            << "Signal:          " << counts[0] << Form(" (%.2f%%)", 100. * counts[0] / total) << '\n'
            << "Other NC:        " << counts[1] << Form(" (%.2f%%)", 100. * counts[1] / total) << '\n'
            << "CCNuMu:          " << counts[2] << Form(" (%.2f%%)", 100. * counts[2] / total) << '\n'
            << "CCNuE:           " << counts[3] << Form(" (%.2f%%)", 100. * counts[3] / total) << '\n'
            << "Dirt:            " << counts[4] << Form(" (%.2f%%)", 100. * counts[4] / total) << '\n'
            << "NonFV:           " << counts[5] << Form(" (%.2f%%)", 100. * counts[5] / total) << '\n'
            << "Cosmic:          " << counts[6] << Form(" (%.2f%%)", 100. * counts[6] / total) << '\n'
            << "No Truth Match:  " << counts[7] << Form(" (%.2f%%)", 100. * counts[7] / total) << '\n'
            << "Bad Reco Signal: " << counts[8] << Form(" (%.2f%%)", 100. * counts[8] / total) << '\n'
            << "EFF:             " << Form("%.2f%%", 100. * counts[0] / total_signal) << '\n'
            << "EFF * PUR:       " << Form("%.2f%%", 100. * counts[0] / total * counts[0] / total_signal) << '\n'
            << std::endl;
}
