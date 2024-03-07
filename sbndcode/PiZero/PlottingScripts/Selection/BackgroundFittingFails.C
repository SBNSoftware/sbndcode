#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

#include "FittedObservables.C"

void ProcessTree(TChain *tree, const double &scale, const std::vector<double> &covMatrix);

void BackgroundFittingFails(const TString productionVersion)
{
  std::vector<double> covMatrix = GetCovMatrix(productionVersion, true);

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile  = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);

  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);
  TChain *intimeSubruns = new TChain("ncpizeroana/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  InitialiseTree(rockboxEvents);
  ProcessTree(rockboxEvents, rockboxScaling, covMatrix);

  InitialiseTree(intimeEvents);
  ProcessTree(intimeEvents, intimeScaling, covMatrix);
}

void ProcessTree(TChain *tree, const double &scale, const std::vector<double> &covMatrix)
{
  const int N = tree->GetEntries();

  double nback = 0, nbackGood = 0;
  
  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      tree->GetEntry(ev_i);
      
      for(int slc_i = 0; slc_i < slc_true_event_type_incl->size(); ++slc_i)
        {
          if(!(slc_true_event_type_incl->at(slc_i) == 0 && slc_comp->at(slc_i) > .5) && slc_sel_incl->at(slc_i))
            {
              ++nback;

              const double en0 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              const double en1 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));

              double corrEn0 = CorrectEnergy(en0);
              double corrEn1 = CorrectEnergy(en1);

              const TVector3 dir0 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                             slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                             slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              const TVector3 dir1 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                             slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                             slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              const double cosineThetaGammaGamma = dir0.Dot(dir1) / (dir0.Mag() * dir1.Mag());
              const double thetaGammaGamma       = acos(cosineThetaGammaGamma);

              bool good = false;
              std::vector<double> updated = DoKF(corrEn0, corrEn1, thetaGammaGamma, covMatrix, good);

              if(good)
                ++nbackGood;
            }
        }
    }

  std::cout << "Fitting good: " << nbackGood * scale << " / " << nback * scale << " (" << (100. * nbackGood) / nback << "%)" << std::endl;
}
