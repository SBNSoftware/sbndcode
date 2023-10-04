#include "/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "RazzledHeaders.h"

void ConfusionMatrix(const PIDTraining &training, const bool efficiency_mode = true,
                     const bool purity_mode = false, const double comp_thresh = 0.5,
                     const double pur_thresh = 0.5);

void ConfusionMatrix()
{
  ConfusionMatrix(razzle_trainings[0], true, false);
  ConfusionMatrix(razzle_trainings[0], false, true);
  ConfusionMatrix(razzle_trainings[1], true, false);
  ConfusionMatrix(razzle_trainings[1], false, true);
  ConfusionMatrix(razzle_trainings[2], true, false);
  ConfusionMatrix(razzle_trainings[2], false, true);
  ConfusionMatrix(dazzle_trainings[0], true, false);
  ConfusionMatrix(dazzle_trainings[0], false, true);
  ConfusionMatrix(dazzle_trainings[1], true, false);
  ConfusionMatrix(dazzle_trainings[1], false, true);

  ConfusionMatrix(razzle_trainings[0], true, false, 0.8, 0.8);
  ConfusionMatrix(razzle_trainings[0], false, true, 0.8, 0.8);
  ConfusionMatrix(razzle_trainings[1], true, false, 0.8, 0.8);
  ConfusionMatrix(razzle_trainings[1], false, true, 0.8, 0.8);
  ConfusionMatrix(razzle_trainings[2], true, false, 0.8, 0.8);
  ConfusionMatrix(razzle_trainings[2], false, true, 0.8, 0.8);
  ConfusionMatrix(dazzle_trainings[0], true, false, 0.8, 0.8);
  ConfusionMatrix(dazzle_trainings[0], false, true, 0.8, 0.8);
  ConfusionMatrix(dazzle_trainings[1], true, false, 0.8, 0.8);
  ConfusionMatrix(dazzle_trainings[1], false, true, 0.8, 0.8);
}

void ConfusionMatrix(const PIDTraining &training, const bool efficiency_mode,
                     const bool purity_mode, const double comp_thresh,
                     const double pur_thresh)
{
  if(efficiency_mode && purity_mode)
    {
      std::cout << "Cannot run both efficiency and purity mode" << std::endl;
      return;
    }

  TString save_dir = "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/razzled/training/confusionmatrices";
  gSystem->Exec("mkdir -p " + save_dir);
  
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *pfps = new TChain("razzled/pfpTree");
  pfps->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root");
  pfps->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intrnue.root");
  pfps->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intime.root");

  int truePDG, recoPDG;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ,
    showerStartX, showerStartY, showerStartZ, showerEnergy;
  bool recoPrimary, unambiguousSlice, trackContained, showerContained;

  float pfp_numDaughters, pfp_maxDaughterHits, pfp_trackScore,
    trk_length, trk_chi2PIDMuon, trk_chi2PIDProton, trk_chi2PIDMuonPionDiff, trk_mcsScatterMean,
    trk_mcsScatterMaxRatio, trk_momDiff, trk_meanDCA, trk_stoppingdEdxChi2Ratio,
    trk_chi2Pol0dEdxFit, shw_bestdEdx, shw_convGap, shw_openAngle, shw_modHitDensity,
    shw_sqrtEnergyDensity;

  pfps->SetBranchAddress("pfp_numDaughters", &pfp_numDaughters);
  pfps->SetBranchAddress("pfp_maxDaughterHits", &pfp_maxDaughterHits);
  pfps->SetBranchAddress("pfp_trackScore", &pfp_trackScore);

  pfps->SetBranchAddress("trk_length", &trk_length);
  pfps->SetBranchAddress("trk_chi2PIDMuon", &trk_chi2PIDMuon);
  pfps->SetBranchAddress("trk_chi2PIDProton", &trk_chi2PIDProton);
  pfps->SetBranchAddress("trk_chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
  pfps->SetBranchAddress("trk_mcsScatterMean", &trk_mcsScatterMean);
  pfps->SetBranchAddress("trk_mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
  pfps->SetBranchAddress("trk_meanDCA", &trk_meanDCA);
  pfps->SetBranchAddress("trk_stoppingdEdxChi2Ratio", &trk_stoppingdEdxChi2Ratio);
  pfps->SetBranchAddress("trk_chi2Pol0dEdxFit", &trk_chi2Pol0dEdxFit);
  pfps->SetBranchAddress("trk_momDiff", &trk_momDiff);

  pfps->SetBranchAddress("shw_bestdEdx", &shw_bestdEdx);
  pfps->SetBranchAddress("shw_convGap", &shw_convGap);
  pfps->SetBranchAddress("shw_openAngle", &shw_openAngle);
  pfps->SetBranchAddress("shw_modHitDensity", &shw_modHitDensity);
  pfps->SetBranchAddress("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);

  pfps->SetBranchAddress("truePDG", &truePDG);
  pfps->SetBranchAddress("energyComp", &energyComp);
  pfps->SetBranchAddress("energyPurity", &energyPurity);
  pfps->SetBranchAddress("recoPrimary", &recoPrimary);
  pfps->SetBranchAddress("unambiguousSlice", &unambiguousSlice);

  pfps->SetBranchAddress("trackStartX", &trackStartX);
  pfps->SetBranchAddress("trackStartY", &trackStartY);
  pfps->SetBranchAddress("trackStartZ", &trackStartZ);
  pfps->SetBranchAddress("trackContained", &trackContained);

  pfps->SetBranchAddress("showerStartX", &showerStartX);
  pfps->SetBranchAddress("showerStartY", &showerStartY);
  pfps->SetBranchAddress("showerStartZ", &showerStartZ);
  pfps->SetBranchAddress("showerContained", &showerContained);
  pfps->SetBranchAddress("showerEnergy", &showerEnergy);

  pfps->SetBranchAddress("recoPDG", &recoPDG);

  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

  if(training.razzled)
    {
      reader->AddVariable("pfp_numDaughters", &pfp_numDaughters);
      reader->AddVariable("pfp_maxDaughterHits", &pfp_maxDaughterHits);
      reader->AddVariable("pfp_trackScore", &pfp_trackScore);
      reader->AddVariable("trk_length", &trk_length);
      reader->AddVariable("trk_chi2PIDMuon", &trk_chi2PIDMuon);
      reader->AddVariable("trk_chi2PIDProton", &trk_chi2PIDProton);
      reader->AddVariable("trk_chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
      reader->AddVariable("trk_mcsScatterMean", &trk_mcsScatterMean);
      reader->AddVariable("trk_mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
      reader->AddVariable("trk_meanDCA", &trk_meanDCA);
      reader->AddVariable("trk_stoppingdEdxChi2Ratio", &trk_stoppingdEdxChi2Ratio);
      reader->AddVariable("trk_chi2Pol0dEdxFit", &trk_chi2Pol0dEdxFit);
      reader->AddVariable("trk_momDiff", &trk_momDiff);
      reader->AddVariable("shw_bestdEdx", &shw_bestdEdx);
      reader->AddVariable("shw_convGap", &shw_convGap);
      reader->AddVariable("shw_openAngle", &shw_openAngle);
      reader->AddVariable("shw_modHitDensity", &shw_modHitDensity);
      reader->AddVariable("shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);
    }
  else if(training.razzle)
    {
      reader->AddVariable("bestdEdx", &shw_bestdEdx);
      reader->AddVariable("convGap", &shw_convGap);
      reader->AddVariable("openAngle", &shw_openAngle);
      reader->AddVariable("modHitDensity", &shw_modHitDensity);
      reader->AddVariable("sqrtEnergyDensity", &shw_sqrtEnergyDensity);
    }
  else if(training.dazzle)
    {
      reader->AddVariable("recoLen", &trk_length);
      reader->AddVariable("chi2PIDMuon", &trk_chi2PIDMuon);
      reader->AddVariable("chi2PIDProton", &trk_chi2PIDProton);
      reader->AddVariable("chi2PIDMuonPionDiff", &trk_chi2PIDMuonPionDiff);
      reader->AddVariable("mcsScatterMean", &trk_mcsScatterMean);
      reader->AddVariable("mcsScatterMaxRatio", &trk_mcsScatterMaxRatio);
      reader->AddVariable("meanDCA", &trk_meanDCA);
      reader->AddVariable("stoppingChi2Ratio", &trk_stoppingdEdxChi2Ratio);
      reader->AddVariable("chi2Pol0Fit", &trk_chi2Pol0dEdxFit);
      reader->AddVariable("pDiff", &trk_momDiff);
      reader->AddVariable("numDaughters", &pfp_numDaughters); 
      reader->AddVariable("maxDaughterHits", &pfp_maxDaughterHits);
    }
  
  reader->BookMVA("BDTG Method", training.path);

  const unsigned N = pfps->GetEntries();

  unsigned n_categories_true = 1, n_categories_reco = 1;

  if(training.razzled)
    { n_categories_true = 5; n_categories_reco = 5; }
  else if(training.razzle)
    { n_categories_true = 3; n_categories_reco = 4; }
  else if(training.dazzle)
    { n_categories_true = 4; n_categories_reco = 5; }

  TH2F *hConfusionMatrix = new TH2F("hConfusionMatrix",
                                    ";True Class;" + training.printed_name + " Class;",
                                    n_categories_true, 0, n_categories_true,
                                    n_categories_reco, 0, n_categories_reco);

  for(unsigned i = 0; i < N; ++i)
    {
      pfps->GetEntry(i);

      if(unambiguousSlice || !recoPrimary)
        continue;
      
      if(abs(truePDG) != 11 && abs(truePDG) != 13 && abs(truePDG) != 22 
         && abs(truePDG) != 211 && abs(truePDG) != 2212)
        continue;
      
      if(training.razzle && (abs(truePDG) == 13 || abs(truePDG) == 211 || abs(truePDG) == 2212))
        truePDG = -1;

      if(training.dazzle && (abs(truePDG) == 11 || abs(truePDG) == 22))
        truePDG = -1;

      if(energyComp < comp_thresh || energyPurity < pur_thresh)
        continue;

      if(abs(trackStartX) > 180 || abs(trackStartY) > 180 || trackStartZ < 10 || trackStartZ > 450 ||
         abs(showerStartX) > 180 || abs(showerStartY) > 180 || showerStartZ < 10 || showerStartZ > 450)
        continue;

      const std::vector<float> bdtScores = reader->EvaluateMulticlass("BDTG Method");

      float bestScore = -std::numeric_limits<float>::max();
      int recoClass = -1, trueClass = -1;

      for(unsigned j = 0; j < bdtScores.size(); ++j)
        {
          if(bdtScores[j] > bestScore)
            {
              bestScore = bdtScores[j];
              recoClass = j;
            }
        }
      
      if(training.razzled)
        {
          trueClass = razzledMap.at(abs(truePDG));
        }
      else if(training.razzle)
        {
          trueClass = razzleMap.at(abs(truePDG));

          if(recoClass == 2)
            recoClass = 3;
          if(recoPDG == 13)
            recoClass = 2;
        }
      else if(training.dazzle)
        {
          trueClass = dazzleMap.at(abs(truePDG));

          if(recoClass == 3)
            recoClass = 4;
          if(recoPDG == 11)
            recoClass = 3;
        }

      hConfusionMatrix->Fill(trueClass, recoClass);
    }

  TCanvas *c = new TCanvas("c", "c");
  c->cd();
  
  if(efficiency_mode)
    NormaliseEntriesByXTotal(hConfusionMatrix);
  if(purity_mode)
    NormaliseEntriesByYTotal(hConfusionMatrix);

  std::vector<TString> axisLabels;
  if(training.razzled)
    axisLabels = razzledAxisLabels;
  else if(training.razzle)
    axisLabels = razzleAxisLabels;
  else if(training.dazzle)
    axisLabels = dazzleAxisLabels;

  for(int i = 1; i <= hConfusionMatrix->GetNbinsX(); ++i)
    hConfusionMatrix->GetXaxis()->SetBinLabel(i, axisLabels[i]);
  for(int i = 1; i <= hConfusionMatrix->GetNbinsY(); ++i)
    hConfusionMatrix->GetYaxis()->SetBinLabel(i, axisLabels[i]);

  gStyle->SetPaintTextFormat("1.2g");
  hConfusionMatrix->SetMarkerSize(3);
  hConfusionMatrix->Draw("col text");

  TString file_name = training.name;

  if(efficiency_mode)
    file_name += "_efficiency";
  if(purity_mode)
    file_name += "_purity";
  if(comp_thresh > 0.5 && pur_thresh > 0.5)
    file_name += "_high_quality";

  c->SaveAs(save_dir + "/" + file_name + ".pdf");
  c->SaveAs(save_dir + "/" + file_name + ".png");
}
