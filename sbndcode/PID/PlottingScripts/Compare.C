#include "RazzledHeaders.h"

void Compare(const std::vector<PIDTraining> &trainings, const TString save_file_suffix,
             const int pdg, const int motherPDG = -1, const double purThresh = 0.5, const double compThresh = 0.5);

void Compare()
{
  Compare(razzled_v12_trainings, "electron", 11);
  Compare(razzled_v12_trainings, "muon", 13);
  Compare(razzled_v12_trainings, "photon", 22);
  Compare(razzled_v12_trainings, "pizero_photon", 22, 111);
  Compare(razzled_v12_trainings, "pion", 211);
  Compare(razzled_v12_trainings, "proton", 2212);

  Compare(razzled_v12_trainings, "electron_high_quality", 11, -1, 0.8, 0.8);
  Compare(razzled_v12_trainings, "muon_high_quality", 13, -1, 0.8, 0.8);
  Compare(razzled_v12_trainings, "photon_high_quality", 22, -1, 0.8, 0.8);
  Compare(razzled_v12_trainings, "pizero_photon_high_quality", 22, 111, 0.8, 0.8);
  Compare(razzled_v12_trainings, "pion_high_quality", 211, -1, 0.8, 0.8);
  Compare(razzled_v12_trainings, "proton_high_quality", 2212, -1, 0.8, 0.8);
}

void Compare(const std::vector<PIDTraining> &trainings, const TString save_file_suffix,
             const int pdg, const int motherPDG, const double purThresh, const double compThresh)
{
  const TString save_dir = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv12/razzled/comparisons";
  gSystem->Exec("mkdir -p " + save_dir);
  gSystem->Exec("mkdir -p " + save_dir + "/sep");

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *pfps = new TChain("razzled/pfpTree");
  pfps->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  pfps->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  pfps->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  int truePDG, trueMotherPDG, recoPDG;
  float energyComp, energyPurity, trackStartX, trackStartY, trackStartZ,
    showerStartX, showerStartY, showerStartZ, showerEnergy;
  bool recoPrimary, unambiguousSlice, trackContained, showerContained;

  float pfp_numDaughters, pfp_maxDaughterHits, pfp_trackScore,
    trk_length, trk_chi2PIDMuon, trk_chi2PIDProton, trk_chi2PIDMuonPionDiff, trk_mcsScatterMean,
    trk_mcsScatterMaxRatio, trk_momDiff, trk_meanDCA, trk_stoppingdEdxChi2Ratio, trk_chi2Pol0dEdxFit,
    shw_bestdEdx, shw_convGap, shw_openAngle, shw_modHitDensity, shw_sqrtEnergyDensity;

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
  pfps->SetBranchAddress("trueMotherPDG", &trueMotherPDG);
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

  TLegend *rocLeg = new TLegend(.25,.25,.55,.45);
  TMultiGraph *multi = new TMultiGraph();

  for(auto const& training : trainings)
    {
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
          if(training.name == "razzled_v12")
            reader->AddVariable("shw_sqrtEnergyDensity>2.5?2.5:shw_sqrtEnergyDensity", &shw_sqrtEnergyDensity);
          else
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

      int N_pfps = pfps->GetEntries();

      TH1F *hSignalPFP     = new TH1F("hSignalPFP" + training.name,
                                      ";" + training.printed_name + " Score;PFPs",
                                      1200, -0.1, 1.1);
      TH1F *hBackgroundPFP = new TH1F("hBackgroundPFP" + training.name,
                                      ";" + training.printed_name + " Score;PFPs",
                                      1200, -0.1, 1.1);

      for(int i = 0; i < N_pfps; ++i)
        {
          pfps->GetEntry(i);

          if(unambiguousSlice || !recoPrimary)
            continue;

          if(abs(trackStartX) > 180 || abs(trackStartY) > 180 || trackStartZ < 10
             || trackStartZ > 450)
            continue;

          if(abs(showerStartX) > 180 || abs(showerStartY) > 180 || showerStartZ < 10
             || showerStartZ > 450)
            continue;

          std::vector<float> bdtscores = reader->EvaluateMulticlass("BDTG Method");
          float bdtscore = -.05;

          if(training.razzled)
            bdtscore = bdtscores[razzledMap.at(pdg)];
          else if(training.razzle)
            bdtscore = bdtscores[razzleMap.at(pdg)];
          else if(training.dazzle)
            bdtscore = bdtscores[dazzleMap.at(pdg)];

          if((training.razzle && recoPDG!=11) || (training.dazzle && recoPDG!=13))
            bdtscore = -.05;

          if(truePDG == pdg && (motherPDG == -1 || trueMotherPDG == motherPDG) && energyPurity > purThresh && energyComp > compThresh)
            hSignalPFP->Fill(bdtscore);
          else if(truePDG != pdg)
            hBackgroundPFP->Fill(bdtscore); 
        }

      float eff[1100], rej[1100];
      int signalSum = 0, backSum = 0;
      const int signalTotal = hSignalPFP->GetEntries();
      const int backTotal   = hBackgroundPFP->GetEntries();

      for(int i = 1; i < 1101; ++i)
        {
          signalSum += hSignalPFP->GetBinContent(i);
          backSum   += hBackgroundPFP->GetBinContent(i);

          eff[i-1] = (float) (signalTotal - signalSum) / (float) signalTotal;
          rej[i-1] = (float) backSum / (float) backTotal;
        }

      TGraph *ROC = new TGraph(1100, eff, rej);
      ROC->SetLineColor(training.colour);
      ROC->SetLineWidth(5);
      multi->Add(ROC);
      rocLeg->AddEntry(ROC, training.printed_name, "l");

      TCanvas *cSep = new TCanvas("cSep", "cSep");
      cSep->cd();
      cSep->SetLogy();

      hSignalPFP->SetLineColor(kBlue+2);
      hBackgroundPFP->SetLineColor(kRed+2);
      hSignalPFP->SetFillColorAlpha(kBlue+2, 0.4);
      hBackgroundPFP->SetFillColorAlpha(kRed+2, 0.4);
      
      hBackgroundPFP->Rebin(20);
      hSignalPFP->Rebin(20);

      hBackgroundPFP->DrawNormalized("][hist");
      hSignalPFP->DrawNormalized("][histsame");

      TLegend* sepLeg = new TLegend(.3, .85, .85, .92);
      sepLeg->SetNColumns(2);
      sepLeg->AddEntry(hSignalPFP, pdgStrings.at(pdg), "lf");
      sepLeg->AddEntry(hBackgroundPFP, "Backgrounds", "lf");
      sepLeg->Draw();

      cSep->SaveAs(save_dir + "/sep/sep_" + save_file_suffix + "_" + training.name + ".png");
      cSep->SaveAs(save_dir + "/sep/sep_" + save_file_suffix + "_" + training.name + ".pdf");

      delete hSignalPFP;
      delete hBackgroundPFP;
    }

  TCanvas *cROC = new TCanvas("cROC","cROC");
  cROC->SetGrid();
  cROC->cd();

  gStyle->SetLabelSize(0.06,"xy");

  multi->SetTitle(";Signal Efficiency;Background Rejection");
  multi->Draw("AC");
  multi->GetXaxis()->SetNdivisions(20);
  multi->GetYaxis()->SetNdivisions(20);
  rocLeg->Draw();

  cROC->SaveAs(save_dir + "/roc_" + save_file_suffix + ".png");
  cROC->SaveAs(save_dir + "/roc_" + save_file_suffix + ".pdf");
}
