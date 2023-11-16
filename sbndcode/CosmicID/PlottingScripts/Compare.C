#include "/sbnd/app/users/hlay/plotting_utils/Structs.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TMVA/Reader.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <set>

const std::vector<CRUMBSTraining> inclusive_trainings = { { "Nominal", "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_22_00/CRUMBS/CRUMBS_SBND.weights.xml", false, true, false, "Nominal", kCyan-3 },
                                                          { "Updated", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/SimpleFlash/CRUMBS_Inclusive/weights/CrumbsTMVAClassification_BDT.weights.xml", false, true, true, "Updated", kMagenta+2 },
                                                          { "OpT0", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/OpT0Finder/CRUMBS_Inclusive/weights/CrumbsTMVAClassification_BDT.weights.xml", true, false, true, "OpT0", kViolet-5 },
};

const std::vector<CRUMBSTraining> ccnumu_trainings = { { "Nominal", "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_22_00/CRUMBS/CRUMBS_CCNuMu_SBND.weights.xml", false, true, false, "Nominal", kCyan-3 },
                                                       { "Updated", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/SimpleFlash/CRUMBS_CCNuMu/weights/CrumbsTMVAClassification_BDT.weights.xml", false, true, true, "Updated", kMagenta+2 },
                                                       { "OpT0", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/OpT0Finder/CRUMBS_CCNuMu/weights/CrumbsTMVAClassification_BDT.weights.xml", true, false, true, "OpT0", kViolet-5 },
};

const std::vector<CRUMBSTraining> ccnue_trainings = { { "Nominal", "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_22_00/CRUMBS/CRUMBS_CCNuE_SBND.weights.xml", false, true, false, "Nominal", kCyan-3 },
                                                      { "Updated", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/SimpleFlash/CRUMBS_CCNuE/weights/CrumbsTMVAClassification_BDT.weights.xml", false, true, true, "Updated", kMagenta+2 },
                                                      { "OpT0", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/OpT0Finder/CRUMBS_CCNuE/weights/CrumbsTMVAClassification_BDT.weights.xml", true, false, true, "OpT0", kViolet-5 },
};

const std::vector<CRUMBSTraining> nc_trainings = { { "Nominal", "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_22_00/CRUMBS/CRUMBS_NC_SBND.weights.xml", false, true, false, "Nominal", kCyan-3 },
                                                   { "Updated", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/SimpleFlash/CRUMBS_NC/weights/CrumbsTMVAClassification_BDT.weights.xml", false, true, true, "Updated", kMagenta+2 },
                                                   { "OpT0", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/OpT0Finder/CRUMBS_NC/weights/CrumbsTMVAClassification_BDT.weights.xml", true, false, true, "OpT0", kViolet-5 },
};

void Compare(const std::vector<CRUMBSTraining> &trainings, const TString save_file_suffix,
             const std::set<int> set_ccnc, const std::set<int> set_nutype,
             const double purThresh = 0.5, const double compThresh = 0.5);

void Compare()
{
  Compare(inclusive_trainings, "inclusive", {0, 1}, {12, 14});
  Compare(ccnumu_trainings, "ccnumu", {0}, {14});
  Compare(ccnue_trainings, "ccnue", {0}, {12});
  Compare(nc_trainings, "nc", {1}, {12, 14});
  Compare(inclusive_trainings, "inclusive_high_quality", {0, 1}, {12, 14}, 0.8, 0.8);
  Compare(ccnumu_trainings, "ccnumu_high_quality", {0}, {14}, 0.8, 0.8);
  Compare(ccnue_trainings, "ccnue_high_quality", {0}, {12}, 0.8, 0.8);
  Compare(nc_trainings, "nc_high_quality", {1}, {12, 14}, 0.8, 0.8);
}

void Compare(const std::vector<CRUMBSTraining> &trainings, const TString save_file_suffix,
             const std::set<int> set_ccnc, const std::set<int> set_nutype,
             const double purThresh, const double compThresh)
{
  const TString save_dir = "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/crumbs/training/comparisons";
  gSystem->Exec("mkdir -p " + save_dir);

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *slices = new TChain("crumbs/SliceTree");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_rockbox.root");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intrnue.root");
  slices->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv2/NCPiZeroAv2_intime.root");

  float reader_tpc_CRFracHitsInLongestTrack, reader_tpc_CRLongestTrackDeflection,
    reader_tpc_CRLongestTrackDirY, reader_tpc_CRNHitsMax,
    reader_tpc_NuEigenRatioInSphere, reader_tpc_NuNFinalStatePfos, reader_tpc_NuNHitsTotal,
    reader_tpc_NuNSpacePointsInSphere, reader_tpc_NuVertexY, reader_tpc_NuWeightedDirZ,
    reader_tpc_StoppingChi2CosmicRatio, reader_pds_FMTotalScore, reader_pds_FMPE,
    reader_pds_FMTime, reader_crt_TrackScore, reader_crt_SPScore, reader_crt_TrackTime,
    reader_crt_SPTime, reader_pds_OpT0Score, reader_pds_OpT0MeasuredPE;

  float tpc_CRFracHitsInLongestTrack, tpc_CRLongestTrackDeflection,
    tpc_CRLongestTrackDirY, tpc_CRNHitsMax,
    tpc_NuEigenRatioInSphere, tpc_NuNFinalStatePfos, tpc_NuNHitsTotal,
    tpc_NuNSpacePointsInSphere, tpc_NuVertexY, tpc_NuWeightedDirZ,
    tpc_StoppingChi2CosmicRatio, pds_FMTotalScore, pds_FMPE,
    pds_FMTime, crt_TrackScore, crt_SPScore, crt_TrackTime,
    crt_SPTime, pds_OpT0Score, pds_OpT0MeasuredPE;

  unsigned eventID, subRunID, runID, slicePDG, sliceIndex;
  int ccnc, nutype;
  std::string *matchedType = 0;
  double matchedPurity, matchedCompleteness;
 
  slices->SetBranchAddress("tpc_CRFracHitsInLongestTrack", &tpc_CRFracHitsInLongestTrack);
  slices->SetBranchAddress("tpc_CRLongestTrackDeflection", &tpc_CRLongestTrackDeflection);
  slices->SetBranchAddress("tpc_CRLongestTrackDirY", &tpc_CRLongestTrackDirY);
  slices->SetBranchAddress("tpc_CRNHitsMax", &tpc_CRNHitsMax);
  slices->SetBranchAddress("tpc_NuEigenRatioInSphere", &tpc_NuEigenRatioInSphere);
  slices->SetBranchAddress("tpc_NuNFinalStatePfos", &tpc_NuNFinalStatePfos);
  slices->SetBranchAddress("tpc_NuNHitsTotal", &tpc_NuNHitsTotal);
  slices->SetBranchAddress("tpc_NuNSpacePointsInSphere", &tpc_NuNSpacePointsInSphere);
  slices->SetBranchAddress("tpc_NuVertexY", &tpc_NuVertexY);
  slices->SetBranchAddress("tpc_NuWeightedDirZ", &tpc_NuWeightedDirZ);
  slices->SetBranchAddress("tpc_StoppingChi2CosmicRatio", &tpc_StoppingChi2CosmicRatio);

  slices->SetBranchAddress("pds_FMTotalScore", &pds_FMTotalScore);
  slices->SetBranchAddress("pds_FMPE", &pds_FMPE);
  slices->SetBranchAddress("pds_FMTime", &pds_FMTime);

  slices->SetBranchAddress("pds_OpT0Score", &pds_OpT0Score);
  slices->SetBranchAddress("pds_OpT0MeasuredPE", &pds_OpT0MeasuredPE);

  slices->SetBranchAddress("crt_TrackScore", &crt_TrackScore);
  slices->SetBranchAddress("crt_SPScore", &crt_SPScore);
  slices->SetBranchAddress("crt_TrackTime", &crt_TrackTime);
  slices->SetBranchAddress("crt_SPTime", &crt_SPTime);

  slices->SetBranchAddress("eventID", &eventID);
  slices->SetBranchAddress("subRunID", &subRunID);
  slices->SetBranchAddress("runID", &runID);
  slices->SetBranchAddress("slicePDG", &slicePDG);
  slices->SetBranchAddress("matchedType", &matchedType);
  slices->SetBranchAddress("matchedPurity", &matchedPurity);
  slices->SetBranchAddress("matchedCompleteness", &matchedCompleteness);
  slices->SetBranchAddress("ccnc", &ccnc);
  slices->SetBranchAddress("nutype", &nutype);

  TLegend *rocLeg = new TLegend(.25,.25,.55,.45);
  TMultiGraph *multi = new TMultiGraph();

  for(auto const& training : trainings)
    {
      TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
      reader->AddVariable("tpc_CRFracHitsInLongestTrack", &reader_tpc_CRFracHitsInLongestTrack);
      reader->AddVariable("tpc_CRLongestTrackDeflection", &reader_tpc_CRLongestTrackDeflection);
      reader->AddVariable("tpc_CRLongestTrackDirY", &reader_tpc_CRLongestTrackDirY);
      reader->AddVariable("tpc_CRNHitsMax", &reader_tpc_CRNHitsMax);
      reader->AddVariable("tpc_NuEigenRatioInSphere", &reader_tpc_NuEigenRatioInSphere);
      reader->AddVariable("tpc_NuNFinalStatePfos", &reader_tpc_NuNFinalStatePfos);
      reader->AddVariable("tpc_NuNHitsTotal", &reader_tpc_NuNHitsTotal);
      reader->AddVariable("tpc_NuNSpacePointsInSphere", &reader_tpc_NuNSpacePointsInSphere);
      reader->AddVariable("tpc_NuVertexY", &reader_tpc_NuVertexY);
      reader->AddVariable("tpc_NuWeightedDirZ", &reader_tpc_NuWeightedDirZ);
      reader->AddVariable("tpc_StoppingChi2CosmicRatio", &reader_tpc_StoppingChi2CosmicRatio);
  
      if(training.simpleFlash)
        {
          reader->AddVariable("pds_FMTotalScore", &reader_pds_FMTotalScore);
          reader->AddVariable("pds_FMPE", &reader_pds_FMPE);
          reader->AddVariable("pds_FMTime", &reader_pds_FMTime);
        }
      
      if(training.opT0)
        {
          reader->AddVariable("pds_OpT0Score", &reader_pds_OpT0Score);
          reader->AddVariable("isinf(pds_OpT0MeasuredPE)?-10000:pds_OpT0MeasuredPE", &reader_pds_OpT0MeasuredPE);
        }

      if(training.newCRT)
        {
          reader->AddVariable("crt_TrackScore", &reader_crt_TrackScore);
          reader->AddVariable("crt_SPScore", &reader_crt_SPScore);
          reader->AddVariable("crt_TrackTime", &reader_crt_TrackTime);
          reader->AddVariable("crt_SPTime", &reader_crt_SPTime);
        }
      else
        {
          reader->AddVariable("crt_TrackScore", &reader_crt_TrackScore);
          reader->AddVariable("crt_HitScore", &reader_crt_SPScore);
          reader->AddVariable("crt_TrackTime", &reader_crt_TrackTime);
          reader->AddVariable("crt_HitTime", &reader_crt_SPTime);
        }
      
      reader->BookMVA("BDT Method", training.path);

      int N_slices = slices->GetEntries();

      TH1F *hNuSlice    = new TH1F("hNuSlice" + training.name, ";Score;Slices", 100, -1, 1);
      TH1F *hNotNuSlice = new TH1F("hNotNuSlice" + training.name, ";Score;Slices", 100, -1, 1);

      for(int i = 0; i < N_slices; ++i)
        {
          slices->GetEntry(i);

          reader_tpc_CRFracHitsInLongestTrack = tpc_CRFracHitsInLongestTrack;
          reader_tpc_CRLongestTrackDeflection = tpc_CRLongestTrackDeflection;
          reader_tpc_CRLongestTrackDirY = tpc_CRLongestTrackDirY;
          reader_tpc_CRNHitsMax = tpc_CRNHitsMax;
          reader_tpc_NuEigenRatioInSphere = tpc_NuEigenRatioInSphere;
          reader_tpc_NuNFinalStatePfos = tpc_NuNFinalStatePfos;
          reader_tpc_NuNHitsTotal = tpc_NuNHitsTotal;
          reader_tpc_NuNSpacePointsInSphere = tpc_NuNSpacePointsInSphere;
          reader_tpc_NuVertexY = tpc_NuVertexY; 
          reader_tpc_NuWeightedDirZ = tpc_NuWeightedDirZ;
          reader_tpc_StoppingChi2CosmicRatio = tpc_StoppingChi2CosmicRatio;

          if(training.simpleFlash)
            {
              reader_pds_FMTotalScore = pds_FMTotalScore;
              reader_pds_FMPE = pds_FMPE;
              reader_pds_FMTime = pds_FMTime;
            }

          if(training.opT0)
            {
              reader_pds_OpT0Score = pds_OpT0Score;
              reader_pds_OpT0MeasuredPE = isinf(pds_OpT0MeasuredPE)?-10000:pds_OpT0MeasuredPE;
            }

          reader_crt_TrackScore = crt_TrackScore;
          reader_crt_SPScore = crt_SPScore;
          reader_crt_TrackTime = crt_TrackTime;
          reader_crt_SPTime = crt_SPTime;

          float bdtscore = reader->EvaluateMVA("BDT Method");

          if(*matchedType == "Nu" && matchedPurity > purThresh && matchedCompleteness > compThresh
             && set_ccnc.count(ccnc) && set_nutype.count(abs(nutype)))
            hNuSlice->Fill(bdtscore);
          else if(*matchedType != "Nu")
            hNotNuSlice->Fill(bdtscore); 
        }

      float eff[100], rej[100];
      int signalSum = 0, backSum = 0;
      const int signalTotal = hNuSlice->GetEntries();
      const int backTotal   = hNotNuSlice->GetEntries();

      for(int i = 1; i < 101; ++i)
        {
          signalSum += hNuSlice->GetBinContent(i);
          backSum   += hNotNuSlice->GetBinContent(i);

          eff[i-1] = (float) (signalTotal - signalSum) / (float) signalTotal;
          rej[i-1] = (float) backSum / (float) backTotal;
        }

      TGraph *ROC = new TGraph(100, rej, eff);
      ROC->SetLineColor(training.colour);
      ROC->SetLineWidth(5);
      multi->Add(ROC);
      rocLeg->AddEntry(ROC, training.printed_name, "l");
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
