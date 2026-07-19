#include <vector>
#include <map>
#include <tuple>
#include <utility>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TH1F.h>
#include <string>
#include <sstream>
#include <TBenchmark.h>
#include <TRandom.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <fstream>
#include <TLegend.h>
#include <THStack.h>
#include <set>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <iomanip>
#include <TH2D.h>
#include <TProfile.h>
#include <TGaxis.h>
#include "TRandom3.h"
#include <algorithm>
#include <functional>
#include <array>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TMatrixDSym.h>
#include <TChain.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

struct weights_struct{
    double signalNuE = 0;
    double BNBNuE = 0;
    double cosmicsNuE = 0;
    double NuENuE = 0;
};

struct recoilElectron_struct{
    double energy;
    double angle;
    double dx;
    double dy;
    double dz;
};

struct highestEnergyPFP_struct{
    double PFPID = -999999;
    double energy = -999999;
    double theta = -999999;
    double dx = -999999;
    double dy = -999999;
    double dz = -999999;
    double vx = -999999;
    double vy = -999999;
    double vz = -999999;
    double completeness = -999999;
    double purity = -999999;
    double trackscore = -999999;
    double primary = -999999;
    double truePDG = -999999;
    double trueOrigin = -999999;
    double trueInt = -999999;
    double bestPlanedEdx = -999999;
    double razzledPDG11 = -999999;
    double razzledPDG13 = -999999;
    double razzledPDG22 = -999999;
    double razzledPDG211 = -999999;
    double razzledPDG2212 = -999999;
    double razzledBestPDG = -999999;
    double trueVX = -999999;
    double trueVY = -999999;
    double trueVZ = -999999;
    double trueEndX = -999999;
    double trueEndY = -999999;
    double trueEndZ = -999999;
    double numHits = -999999;
    double clearCosmic = -999999;
};

struct GenieParam_struct{
    std::string mapKey;
    std::string shortName;
    int nUniv;
    bool isMultisim;
    int origNUniv;   // Raw universe count before any expansion (6, 1, 100, or 2/4/7/10)
    bool skipForNow; // True for the 5 knobs not understood yet
};

std::vector<std::string> listRootFiles(const std::string& dirPath){
    std::vector<std::string> fileList;
    TSystemDirectory dir("inputDir", dirPath.c_str());
    TList *files = dir.GetListOfFiles();
    if(files){
        TSystemFile *file;
        TIter next(files);
        while((file = (TSystemFile*)next())){
            std::string fname = file->GetName();
            if(!file->IsDirectory() && fname.size() > 5 && fname.compare(fname.size()-5, 5, ".root") == 0){
                fileList.push_back(dirPath + "/" + fname);
            }
        }
    }
    std::sort(fileList.begin(), fileList.end());
    return fileList;
}

void nuECovarainceMatrixXSec_macro(){

    std::string cutsApplied = "allCuts";
    std::string base_path = "/nashome/c/coackley/systPlotsGenieXSecCov_" + cutsApplied + "/";

    int clearCosmicCut = 1;
    int numPFPs0Cut = 1;
    int numRecoNeutrinosCut = 1;
    int CRUMBSCut = 1;
    int FVCut = 1;
    int primaryPFPCut = 1;
    int ETheta2Cut = 1;
    int razzledPDG11Cut = 1;
    int razzledPDG211Cut = 1;
    int dEdxCut = 1;

    double crumbsScoreCut_low = 0.2;
    double crumbsScoreCut_high = 0.76;

    double FVCut_xHigh = 192;
    double FVCut_xLow = -192;
    double FVCut_xCentre = 10;

    double FVCut_yHigh = 194;
    double FVCut_yLow = -194;

    double FVCut_zHigh = 450;
    double FVCut_zLow = 6;

    double primaryPFPCutValue = 1;

    double razzled211High_highestEnergyPFP = 0.0125;
    double razzled211Low_highestEnergyPFP = 0.00;
    double razzled11High_highestEnergyPFP = 1;
    double razzled11Low_highestEnergyPFP = 0.875;

    double dEdxHigh_highestEnergyPFP = 3.25;
    double dEdxLow_highestEnergyPFP = 0;

    double ETheta2High_highestEnergyPFP = 3.066;
    double ETheta2Low_highestEnergyPFP = 0;

    if (gSystem->AccessPathName(base_path.c_str())) {
        gSystem->mkdir(base_path.c_str(), kTRUE);
    }

    std::string inputDir = "/exp/sbnd/data/users/coackley/analysisFiles_14Jul/";
    std::vector<std::string> inputFiles = listRootFiles(inputDir);
    std::cout << "Found " << inputFiles.size() << " input files in " << inputDir << std::endl;
    if(inputFiles.empty()){
        std::cerr << "No input files found in " << inputDir << std::endl;
        return;
    }

    TChain *tree = new TChain("ana/NuE");
    TChain *subRunTree = new TChain("ana/SubRun");
    TChain *weightsTree = new TChain("ana/NuEWeights");

    for(const auto& f : inputFiles){
        tree->Add(f.c_str());
        weightsTree->Add(f.c_str());
        subRunTree->Add(f.c_str());
    }

    if(tree->GetEntries() != weightsTree->GetEntries()){
        std::cerr << "FATAL: NuE has " << tree->GetEntries() << " entries but NuEWeights has " << weightsTree->GetEntries() << " entries — they must be 1:1" << std::endl;
        return;
    }

    std::cout << "Chained " << tree->GetEntries() << " events across " << inputFiles.size() << " files (" << subRunTree->GetEntries() << " subrun entries)" << std::endl;

    TFile *fOut = new TFile("/exp/sbnd/data/users/coackley/selectionCovarianceMatrixGenieXSec.root", "RECREATE");
    if(!fOut || fOut->IsZombie()){
        std::cerr << "Error creating output ROOT file" << std::endl;
        return;
    }

    double subRunSignal, subRunDLCurrent, subRunPOT;
    int subRunSpills, subRunNumGenEvents;
    unsigned int subRunNumber, subRunRun;

    subRunTree->SetBranchAddress("signal", &subRunSignal);
    subRunTree->SetBranchAddress("DLCurrent", &subRunDLCurrent);
    subRunTree->SetBranchAddress("pot", &subRunPOT);
    subRunTree->SetBranchAddress("spills", &subRunSpills);
    subRunTree->SetBranchAddress("numGenEvents", &subRunNumGenEvents);
    subRunTree->SetBranchAddress("subRun", &subRunNumber);
    subRunTree->SetBranchAddress("run", &subRunRun);

    Long64_t numEntriesSubRun = subRunTree->GetEntries();

    double cosmicSpillsSumNuE = 0;
    double BNBSpillsSumNuE = 0;
    double SignalSpillsSumNuE = 0;

    double POTSignalNuE_notMissing = 0;
    double POTBNBNuE_notMissing = 0;
    double POTNuENuE_notMissing = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 2 && subRunDLCurrent == 5) BNBSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 1 && subRunDLCurrent == 5) SignalSpillsSumNuE += subRunNumGenEvents;

        if(subRunSignal == 1){
            if(subRunDLCurrent == 5) POTSignalNuE_notMissing += subRunPOT;
        } else if(subRunSignal == 2){
            if(subRunDLCurrent == 5) POTBNBNuE_notMissing += subRunPOT;
        } else if(subRunSignal == 4){
            if(subRunDLCurrent == 5) POTNuENuE_notMissing += subRunPOT;
        }
    }

    double targetPOT = 1e21;
    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);
    double totalPOTNuENuE_notMissing = (POTNuENuE_notMissing + POTBNBNuE_notMissing);

    weights_struct weights;
    weights.signalNuE = targetPOT / POTSignalNuE_notMissing;
    weights.BNBNuE = targetPOT / POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;
    weights.NuENuE = targetPOT / totalPOTNuENuE_notMissing;

    std::cout << "Weights DLNu+E: BNB = " << weights.BNBNuE << ", Signal = " << weights.signalNuE << ", Intime Cosmics = " << weights.cosmicsNuE << ", nue = " << weights.NuENuE << std::endl;

    UInt_t eventID, runID, subRunID;
    int nuEScatter;
    double nuEScatterTrueVX, nuEScatterTrueVY, nuEScatterTrueVZ;
    double DLCurrent, signal;

    std::vector<double> *truth_recoilElectronPDG = nullptr;
    std::vector<double> *truth_recoilElectronEnergy = nullptr;
    std::vector<double> *truth_recoilElectronAngle = nullptr;
    std::vector<double> *truth_recoilElectronDX = nullptr;
    std::vector<double> *truth_recoilElectronDY = nullptr;
    std::vector<double> *truth_recoilElectronDZ = nullptr;

    std::vector<double> *reco_sliceID = nullptr;
    std::vector<double> *reco_sliceCompleteness = nullptr;
    std::vector<double> *reco_sliceScore = nullptr;
    std::vector<double> *reco_sliceTrueVX = nullptr;
    std::vector<double> *reco_sliceTrueVY = nullptr;
    std::vector<double> *reco_sliceTrueVZ = nullptr;
    std::vector<double> *reco_sliceOrigin = nullptr;
    std::vector<double> *reco_sliceTrueCCNC = nullptr;
    std::vector<double> *reco_sliceTrueNeutrinoType = nullptr;

    std::vector<double> *reco_particlePDG = nullptr;
    std::vector<double> *reco_particleIsPrimary = nullptr;
    std::vector<double> *reco_particleVX = nullptr;
    std::vector<double> *reco_particleVY = nullptr;
    std::vector<double> *reco_particleVZ = nullptr;
    std::vector<double> *reco_particleDX = nullptr;
    std::vector<double> *reco_particleDY = nullptr;
    std::vector<double> *reco_particleDZ = nullptr;
    std::vector<double> *reco_particleSliceID = nullptr;
    std::vector<double> *reco_particleBestPlaneEnergy = nullptr;
    std::vector<double> *reco_particleTheta = nullptr;
    std::vector<double> *reco_particleTrackScore = nullptr;
    std::vector<double> *reco_particleCompleteness = nullptr;
    std::vector<double> *reco_particlePurity = nullptr;
    std::vector<double> *reco_particleID = nullptr;
    std::vector<double> *reco_particleTruePDG = nullptr;
    std::vector<double> *reco_particleTrueOrigin = nullptr;
    std::vector<double> *reco_particleTrueInteractionType = nullptr;
    std::vector<double> *reco_particleNumHits = nullptr;
    std::vector<double> *reco_particleClearCosmic = nullptr;
    std::vector<double> *reco_particleBestPlanedEdx = nullptr;
    std::vector<double> *reco_particleRazzledPDG11 = nullptr;
    std::vector<double> *reco_particleRazzledPDG13 = nullptr;
    std::vector<double> *reco_particleRazzledPDG22 = nullptr;
    std::vector<double> *reco_particleRazzledPDG211 = nullptr;
    std::vector<double> *reco_particleRazzledPDG2212 = nullptr;
    std::vector<double> *reco_particleRazzledBestPDG = nullptr;
    std::vector<double> *reco_particleTrueVX = nullptr;
    std::vector<double> *reco_particleTrueVY = nullptr;
    std::vector<double> *reco_particleTrueVZ = nullptr;
    std::vector<double> *reco_particleTrueEndX = nullptr;
    std::vector<double> *reco_particleTrueEndY = nullptr;
    std::vector<double> *reco_particleTrueEndZ = nullptr;

    std::vector<double> *reco_neutrinoID = nullptr;
    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;

    std::vector<double> *angleRecalculationPCAPFP10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_pfpID = nullptr;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("nuEScatter", &nuEScatter);
    tree->SetBranchAddress("nuEScatterTrueVX", &nuEScatterTrueVX);
    tree->SetBranchAddress("nuEScatterTrueVY", &nuEScatterTrueVY);
    tree->SetBranchAddress("nuEScatterTrueVZ", &nuEScatterTrueVZ);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);
    tree->SetBranchAddress("signal", &signal);

    tree->SetBranchAddress("truth_recoilElectronPDG", &truth_recoilElectronPDG);
    tree->SetBranchAddress("truth_recoilElectronEnergy", &truth_recoilElectronEnergy);
    tree->SetBranchAddress("truth_recoilElectronAngle", &truth_recoilElectronAngle);
    tree->SetBranchAddress("truth_recoilElectronDX", &truth_recoilElectronDX);
    tree->SetBranchAddress("truth_recoilElectronDY", &truth_recoilElectronDY);
    tree->SetBranchAddress("truth_recoilElectronDZ", &truth_recoilElectronDZ);

    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness", &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_sliceScore", &reco_sliceScore);
    tree->SetBranchAddress("reco_sliceTrueVX", &reco_sliceTrueVX);
    tree->SetBranchAddress("reco_sliceTrueVY", &reco_sliceTrueVY);
    tree->SetBranchAddress("reco_sliceTrueVZ", &reco_sliceTrueVZ);
    tree->SetBranchAddress("reco_sliceOrigin", &reco_sliceOrigin);
    tree->SetBranchAddress("reco_sliceTrueCCNC", &reco_sliceTrueCCNC);
    tree->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);

    tree->SetBranchAddress("reco_particlePDG", &reco_particlePDG);
    tree->SetBranchAddress("reco_particleIsPrimary", &reco_particleIsPrimary);
    tree->SetBranchAddress("reco_particleVX", &reco_particleVX);
    tree->SetBranchAddress("reco_particleVY", &reco_particleVY);
    tree->SetBranchAddress("reco_particleVZ", &reco_particleVZ);
    tree->SetBranchAddress("reco_particleDX", &reco_particleDX);
    tree->SetBranchAddress("reco_particleDY", &reco_particleDY);
    tree->SetBranchAddress("reco_particleDZ", &reco_particleDZ);
    tree->SetBranchAddress("reco_particleSliceID", &reco_particleSliceID);
    tree->SetBranchAddress("reco_particleBestPlaneEnergy", &reco_particleBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTheta", &reco_particleTheta);
    tree->SetBranchAddress("reco_particleTrackScore", &reco_particleTrackScore);
    tree->SetBranchAddress("reco_particleCompleteness", &reco_particleCompleteness);
    tree->SetBranchAddress("reco_particlePurity", &reco_particlePurity);
    tree->SetBranchAddress("reco_particleID", &reco_particleID);
    tree->SetBranchAddress("reco_particleTruePDG", &reco_particleTruePDG);
    tree->SetBranchAddress("reco_particleTrueOrigin", &reco_particleTrueOrigin);
    tree->SetBranchAddress("reco_particleTrueInteractionType", &reco_particleTrueInteractionType);
    tree->SetBranchAddress("reco_particleNumHits", &reco_particleNumHits);
    tree->SetBranchAddress("reco_particleClearCosmic", &reco_particleClearCosmic);
    tree->SetBranchAddress("reco_particleBestPlanedEdx", &reco_particleBestPlanedEdx);
    tree->SetBranchAddress("reco_particleRazzledPDG11", &reco_particleRazzledPDG11);
    tree->SetBranchAddress("reco_particleRazzledPDG13", &reco_particleRazzledPDG13);
    tree->SetBranchAddress("reco_particleRazzledPDG22", &reco_particleRazzledPDG22);
    tree->SetBranchAddress("reco_particleRazzledPDG211", &reco_particleRazzledPDG211);
    tree->SetBranchAddress("reco_particleRazzledPDG2212", &reco_particleRazzledPDG2212);
    tree->SetBranchAddress("reco_particleRazzledBestPDG", &reco_particleRazzledBestPDG);
    tree->SetBranchAddress("reco_particleTrueVX", &reco_particleTrueVX);
    tree->SetBranchAddress("reco_particleTrueVY", &reco_particleTrueVY);
    tree->SetBranchAddress("reco_particleTrueVZ", &reco_particleTrueVZ);
    tree->SetBranchAddress("reco_particleTrueEndX", &reco_particleTrueEndX);
    tree->SetBranchAddress("reco_particleTrueEndY", &reco_particleTrueEndY);
    tree->SetBranchAddress("reco_particleTrueEndZ", &reco_particleTrueEndZ);

    tree->SetBranchAddress("reco_neutrinoID", &reco_neutrinoID);
    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);

    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_angle", &angleRecalculationPCAPFP10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_pfpID", &angleRecalculationPCAPFP10cm_pfpID);

    Long64_t numEntries = tree->GetEntries();

    double nuEScatterTrueVX_weights, nuEScatterTrueVY_weights, nuEScatterTrueVZ_weights;
    std::vector<double> *reco_sliceID_weights = nullptr;

    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE = nullptr;

    weightsTree->SetBranchAddress("nuEScatterTrueVX_weightTree", &nuEScatterTrueVX_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVY_weightTree", &nuEScatterTrueVY_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVZ_weightTree", &nuEScatterTrueVZ_weights);
    weightsTree->SetBranchAddress("reco_sliceID_weightTree", &reco_sliceID_weights);

    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE);
    weightsTree->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE);

    std::vector<GenieParam_struct> genieParams = {
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi","NonResPionNorm_NR_nu_n_CC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi","NonResPionNorm_NR_nu_n_CC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi","NonResPionNorm_NR_nu_n_NC_1Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi","NonResPionNorm_NR_nu_n_NC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi","NonResPionNorm_NR_nu_n_NC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi","NonResPionNorm_NR_nu_np_CC_1Pi",7,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi","NonResPionNorm_NR_nu_p_CC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi","NonResPionNorm_NR_nu_p_CC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi","NonResPionNorm_NR_nu_p_NC_1Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi","NonResPionNorm_NR_nu_p_NC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi","NonResPionNorm_NR_nu_p_NC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi","NonResPionNorm_NR_nubar_n_CC_1Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi","NonResPionNorm_NR_nubar_n_CC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi","NonResPionNorm_NR_nubar_n_CC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi","NonResPionNorm_NR_nubar_n_NC_1Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi","NonResPionNorm_NR_nubar_n_NC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi","NonResPionNorm_NR_nubar_n_NC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi","NonResPionNorm_NR_nubar_p_CC_1Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi","NonResPionNorm_NR_nubar_p_CC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi","NonResPionNorm_NR_nubar_p_CC_3Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi","NonResPionNorm_NR_nubar_p_NC_1Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi","NonResPionNorm_NR_nubar_p_NC_2Pi",6,false},
        {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi","NonResPionNorm_NR_nubar_p_NC_3Pi",6,false},

        {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu","Misc_C12ToAr40_2p2hScaling_nu",6,false},
        {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar","Misc_C12ToAr40_2p2hScaling_nubar",6,false},
        {"MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression","Misc_SPPLowQ2Suppression",10,false},
        {"MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio","Misc_nuenuebar_xsec_ratio",2,false},
        {"MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio","Misc_nuenumu_xsec_ratio",2,false},

        {"MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement","MINERvAq0q3_Mnv2p2hGaussEnhancement",4,false},

        {"MINERvAE2p2h_SBN_v1_E2p2h_A_nu","MINERvAE2p2h_E2p2h_A_nu",6,false},
        {"MINERvAE2p2h_SBN_v1_E2p2h_A_nubar","MINERvAE2p2h_E2p2h_A_nubar",6,false},
        {"MINERvAE2p2h_SBN_v1_E2p2h_B_nu","MINERvAE2p2h_E2p2h_B_nu",6,false},
        {"MINERvAE2p2h_SBN_v1_E2p2h_B_nubar","MINERvAE2p2h_E2p2h_B_nubar",6,false},

        {"GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse","multisim_CCRESVariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_COHVariationResponse","multisim_COHVariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_CoulombCCQE","multisim_CoulombCCQE",100,true},
        {"GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse","multisim_DISBYVariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse","multisim_FSI_N_VariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse","multisim_FSI_pi_VariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NCELVariationResponse","multisim_NCELVariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse","multisim_NCRESVariationResponse",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi","multisim_NonRESBGvbarnCC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi","multisim_NonRESBGvbarnCC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi","multisim_NonRESBGvbarnNC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi","multisim_NonRESBGvbarnNC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi","multisim_NonRESBGvbarpCC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi","multisim_NonRESBGvbarpCC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi","multisim_NonRESBGvbarpNC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi","multisim_NonRESBGvbarpNC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi","multisim_NonRESBGvnCC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi","multisim_NonRESBGvnCC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi","multisim_NonRESBGvnNC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi","multisim_NonRESBGvnNC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi","multisim_NonRESBGvpCC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi","multisim_NonRESBGvpCC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi","multisim_NonRESBGvpNC1pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi","multisim_NonRESBGvpNC2pi",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NormCCMEC","multisim_NormCCMEC",100,true},
        {"GENIEReWeight_SBN_v1_multisim_NormNCMEC","multisim_NormNCMEC",100,true},
        {"GENIEReWeight_SBN_v1_multisim_RDecBR1eta","multisim_RDecBR1eta",100,true},
        {"GENIEReWeight_SBN_v1_multisim_RDecBR1gamma","multisim_RDecBR1gamma",100,true},
        {"GENIEReWeight_SBN_v1_multisim_RPA_CCQE","multisim_RPA_CCQE",100,true},
        {"GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse","multisim_ZExpAVariationResponse",100,true},

        {"GENIEReWeight_SBN_v1_multisigma_AhtBY","multisigma_AhtBY",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_BhtBY","multisigma_BhtBY",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_CV1uBY","multisigma_CV1uBY",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_CV2uBY","multisigma_CV2uBY",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_CoulombCCQE","multisigma_CoulombCCQE",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_DecayAngMEC","multisigma_DecayAngMEC",1,false},
        {"GENIEReWeight_SBN_v1_multisigma_EtaNCEL","multisigma_EtaNCEL",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrAbs_N","multisigma_FrAbs_N",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrAbs_pi","multisigma_FrAbs_pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrCEx_N","multisigma_FrCEx_N",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrCEx_pi","multisigma_FrCEx_pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrInel_N","multisigma_FrInel_N",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrInel_pi","multisigma_FrInel_pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrPiProd_N","multisigma_FrPiProd_N",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi","multisigma_FrPiProd_pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MFP_N","multisigma_MFP_N",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MFP_pi","multisigma_MFP_pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MaCCRES","multisigma_MaCCRES",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MaNCEL","multisigma_MaNCEL",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MaNCRES","multisigma_MaNCRES",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MvCCRES","multisigma_MvCCRES",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_MvNCRES","multisigma_MvNCRES",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi","multisigma_NonRESBGvbarnCC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi","multisigma_NonRESBGvbarnCC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi","multisigma_NonRESBGvbarnNC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi","multisigma_NonRESBGvbarnNC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi","multisigma_NonRESBGvbarpCC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi","multisigma_NonRESBGvbarpCC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi","multisigma_NonRESBGvbarpNC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi","multisigma_NonRESBGvbarpNC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi","multisigma_NonRESBGvnCC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi","multisigma_NonRESBGvnCC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi","multisigma_NonRESBGvnNC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi","multisigma_NonRESBGvnNC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi","multisigma_NonRESBGvpCC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi","multisigma_NonRESBGvpCC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi","multisigma_NonRESBGvpNC1pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi","multisigma_NonRESBGvpNC2pi",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NormCCCOH","multisigma_NormCCCOH",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NormCCMEC","multisigma_NormCCMEC",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NormNCCOH","multisigma_NormNCCOH",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_NormNCMEC","multisigma_NormNCMEC",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_RDecBR1eta","multisigma_RDecBR1eta",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma","multisigma_RDecBR1gamma",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_RPA_CCQE","multisigma_RPA_CCQE",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad","multisigma_ThetaDelta2NRad",1,false},
        {"GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi","multisigma_Theta_Delta2Npi",1,false},
        {"GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape","multisigma_VecFFCCQEshape",1,false},
        {"GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE","multisigma_ZExpA1CCQE",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE","multisigma_ZExpA2CCQE",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE","multisigma_ZExpA3CCQE",6,false},
        {"GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE","multisigma_ZExpA4CCQE",6,false},
    };

    const int NPARAMS_GENIE = (int)genieParams.size();
    std::cout << "Number of GENIE parameters loaded: " << NPARAMS_GENIE << " out of 115" << std::endl;

    for(auto& gp : genieParams){
        gp.origNUniv = gp.nUniv;
        gp.skipForNow = false;

        if(gp.isMultisim) continue;

        if(gp.nUniv == 6 || gp.nUniv == 1){
            gp.nUniv = 100;
            gp.isMultisim = true;
        } else {
            gp.skipForNow = true;
        }
    }

    int nActiveGenieParams = 0;
    for(const auto& gp : genieParams) if(!gp.skipForNow) nActiveGenieParams++;
    std::cout << "Active GENIE parameters after skipping not-yet-understood knobs: " << nActiveGenieParams << " / " << NPARAMS_GENIE << std::endl;

    std::map<std::string, std::vector<float>> sigmaUMap;
    {
        TFile* sigmaUFile = TFile::Open("/exp/sbnd/data/users/coackley/sigma_u_universes.root", "READ");
        if(!sigmaUFile || sigmaUFile->IsZombie()){
            std::cerr << "FATAL: could not open sigma_u_universes.root" << std::endl;
            return;
        }
        TTree* sigmaUTree = (TTree*)sigmaUFile->Get("sigma_u_tree");
        if(!sigmaUTree){
            std::cerr << "FATAL: sigma_u_tree not found in sigma_u_universes.root" << std::endl;
            return;
        }

        std::vector<std::string> knobsNeedingSigmaU;
        for(const auto& gp : genieParams){
            if(!gp.skipForNow && (gp.origNUniv == 6 || gp.origNUniv == 1)){
                knobsNeedingSigmaU.push_back(gp.mapKey);
            }
        }

        std::map<std::string, float> sigmaUBranchVal;
        for(const auto& key : knobsNeedingSigmaU){
            sigmaUBranchVal[key] = 0.0f;
            sigmaUTree->SetBranchAddress((key + "_sigmau").c_str(), &sigmaUBranchVal[key]);
            sigmaUMap[key].assign(100, 0.0f);
        }

        Long64_t nSigmaURows = sigmaUTree->GetEntries();
        for(Long64_t row = 0; row < nSigmaURows && row < 100; ++row){
            sigmaUTree->GetEntry(row);
            for(const auto& key : knobsNeedingSigmaU){
                sigmaUMap[key][row] = sigmaUBranchVal[key];
            }
        }

        sigmaUFile->Close();
        std::cout << "Loaded sigma_u draws for " << knobsNeedingSigmaU.size() << " knobs from sigma_u_universes.root" << std::endl;
    }

    auto expandToPseudoMultisim = [&](double wRaw, int rawN, double sigmaU) -> double {
        double wu;
        if(rawN == 6) wu = 1.0 + (wRaw - 1.0) * sigmaU;
        else wu = 1.0 + (wRaw - 1.0) * 2.0 * std::fabs(sigmaU);
        return std::clamp(wu, 0.0, 10.0);
    };

    std::vector<std::vector<std::vector<float>>*> reco_slice_GENIEWeightVecs = {
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi,
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi,
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu,
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar,
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression,
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio,
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio,
        reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement,
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu,
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar,
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu,
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE,
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE
    };

    if((int)reco_slice_GENIEWeightVecs.size() != NPARAMS_GENIE){
        std::cerr << "ERROR: reco_slice_GENIEWeightVecs size doesn't match genieParams (" << NPARAMS_GENIE << ")" << std::endl;
        return;
    }

    long long sliceWeightFallbacks = 0, sliceWeightCalls = 0;

    auto getGenieSliceWeight = [&](std::vector<std::vector<float>>* vec, size_t sliceIdx, int u, bool wFound, int expectedN) -> double {
        sliceWeightCalls++;
        if(!wFound || !vec || sliceIdx >= vec->size()){
            sliceWeightFallbacks++;
            return 1.0;
        }
        if((int)vec->at(sliceIdx).size() != expectedN){
            sliceWeightFallbacks++;
            return 1.0;
        }
        if(u < 0 || u >= expectedN){
            sliceWeightFallbacks++;
            return 1.0;
        }
        return vec->at(sliceIdx).at(u);
    };


    const int NSTAGES = 10; // 0 = no cuts, 1..9 = cumulative cuts applied in order below
    const std::vector<std::string> stageNames = {
        "s0_noCuts","s1_numPFPs0","s2_numRecoNeutrinos","s3_CRUMBS","s4_FV",
        "s5_primaryPFP","s6_razzledPDG11","s7_razzledPDG211","s8_ETheta2","s9_dEdx"
    };

    const int NCATS = 3; // 0 = all, 1 = signal, 2 = background
    const std::vector<std::string> catNames = {"all","signal","background"};

    const int NUNIV_GENIE = 100; // every ACTIVE genie param has been expanded/validated to 100 universes

    const int NBINS_STUDY  = 15;
    double angleLow_study  = 0.0;
    double angleHigh_study = 20.0;

    std::vector<double> CVnominal((size_t)NSTAGES*NCATS*NBINS_STUDY, 0.0);

    std::vector<double> UNIVsum((size_t)NSTAGES*NCATS*NPARAMS_GENIE*NUNIV_GENIE*NBINS_STUDY, 0.0);

    auto CVidx = [&](int stage,int cat,int bin)->size_t{
        return (size_t)((stage*NCATS + cat)*NBINS_STUDY + bin);
    };
    auto UNIVidx = [&](int stage,int cat,int param,int u,int bin)->size_t{
        return ((((size_t)stage*NCATS + cat)*NPARAMS_GENIE + param)*NUNIV_GENIE + u)*NBINS_STUDY + bin;
    };
    auto angleToBin_study = [&](double angleDeg)->int{
        if(angleDeg < angleLow_study || angleDeg >= angleHigh_study) return -1;
        int b = (int)((angleDeg - angleLow_study)/(angleHigh_study - angleLow_study)*NBINS_STUDY);
        if(b < 0) b = 0;
        if(b >= NBINS_STUDY) b = NBINS_STUDY - 1;
        return b;
    };

    std::vector<std::vector<double>> sliceUnivWeight_genie(NPARAMS_GENIE, std::vector<double>(NUNIV_GENIE, 1.0));

    Long64_t numEntries2 = numEntries;

    for(Long64_t e = 0; e < numEntries2; ++e){
        tree->GetEntry(e);
        weightsTree->GetEntry(e);

        const double epsCheck = 1e-6;
        if(std::abs(nuEScatterTrueVX - nuEScatterTrueVX_weights) > epsCheck || std::abs(nuEScatterTrueVY - nuEScatterTrueVY_weights) > epsCheck || std::abs(nuEScatterTrueVZ - nuEScatterTrueVZ_weights) > epsCheck){
            std::cerr << "ERROR: entry " << e << " misaligned between NuE and NuEWeights trees!" << std::endl;
            continue;
        }

        recoilElectron_struct recoilElectron;
        recoilElectron.energy = -999999; recoilElectron.dx = -999999; recoilElectron.dy = -999999; recoilElectron.dz = -999999;
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->at(i) != -999999){
                recoilElectron.energy = truth_recoilElectronEnergy->at(i);
                recoilElectron.angle = truth_recoilElectronAngle->at(i);
                recoilElectron.dx = truth_recoilElectronDX->at(i);
                recoilElectron.dy = truth_recoilElectronDY->at(i);
                recoilElectron.dz = truth_recoilElectronDZ->at(i);
            } else if(truth_recoilElectronPDG->size() == 1 && truth_recoilElectronPDG->at(i) == -999999){
                recoilElectron.energy = -999999;
                recoilElectron.angle = -999999;
                recoilElectron.dx = -999999;
                recoilElectron.dy = -999999;
                recoilElectron.dz = -999999;
            }
        }

        if(reco_sliceID->size() == 0) continue;

        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) == -999999) continue;

            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                sliceCategoryPlottingMacro = 0;
            } else if(reco_sliceOrigin->at(slice) == 1){
                if(reco_sliceCompleteness->at(slice) > 0.5 && recoilElectron.energy > 150){
                    if(FVCut == 0 && (reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVX->at(slice) > -201.3) && (reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVY->at(slice) > -203.8) && (reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5)){
                        sliceCategoryPlottingMacro = 1;
                    } else if(FVCut == 1 && (reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow && std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre) && (reco_sliceTrueVY->at(slice) < FVCut_yHigh && reco_sliceTrueVY->at(slice) > FVCut_yLow) && (reco_sliceTrueVZ->at(slice) < FVCut_zHigh && reco_sliceTrueVZ->at(slice) > FVCut_zLow)){
                        sliceCategoryPlottingMacro = 1;
                    } else{
                        sliceCategoryPlottingMacro = 2;
                    }
                } else{
                    sliceCategoryPlottingMacro = 2;
                }
            } else if(reco_sliceOrigin->at(slice) == 3){
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    if(reco_sliceTrueCCNC->at(slice) == 0 && reco_sliceTrueNeutrinoType->at(slice) == 12){
                        sliceCategoryPlottingMacro = 5;
                    } else{
                        sliceCategoryPlottingMacro = 3;
                    }
                } else{
                    if(reco_sliceTrueCCNC->at(slice) == 0 && reco_sliceTrueNeutrinoType->at(slice) == 12){
                        sliceCategoryPlottingMacro = 6;
                    } else{
                        sliceCategoryPlottingMacro = 4;
                    }
                }
            }

            double weight = 0;
            if(signal == 1 && DLCurrent == 5){
                weight = weights.signalNuE;
            } else if(signal == 2 && DLCurrent == 5 && sliceCategoryPlottingMacro != 5 && sliceCategoryPlottingMacro != 6){
                weight = weights.BNBNuE;
            } else if((signal == 2 || signal == 4) && DLCurrent == 5 && (sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6)){
                weight = weights.NuENuE;
            } else if(signal == 3 && DLCurrent == 5){
                weight = weights.cosmicsNuE;
            }

            double numPFPsSlice = 0;
            double numPrimaryPFPs10Slice = 0;
            highestEnergyPFP_struct highestEnergyPFP;

            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                    if((clearCosmicCut == 1 && reco_particleClearCosmic->at(pfp) == 0) || clearCosmicCut == 0){
                        numPFPsSlice++;
                        if(reco_particleIsPrimary->at(pfp) == 1){
                            if(reco_particleNumHits->at(pfp) >= 10) numPrimaryPFPs10Slice++;
                        }

                        if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergyPFP.energy){
                            highestEnergyPFP.energy = reco_particleBestPlaneEnergy->at(pfp);
                            highestEnergyPFP.PFPID = reco_particleID->at(pfp);
                            highestEnergyPFP.bestPlanedEdx = reco_particleBestPlanedEdx->at(pfp);
                            highestEnergyPFP.razzledPDG11 = reco_particleRazzledPDG11->at(pfp);
                            highestEnergyPFP.razzledPDG211 = reco_particleRazzledPDG211->at(pfp);
                        }
                    }
                }
            }

            double pfp10cm_PCAAngle = -999999;
            for(size_t pfpAngle = 0; pfpAngle < angleRecalculationPCAPFP10cm_pfpID->size(); ++pfpAngle){
                if(angleRecalculationPCAPFP10cm_pfpID->at(pfpAngle) == highestEnergyPFP.PFPID){
                    pfp10cm_PCAAngle = angleRecalculationPCAPFP10cm_angle->at(pfpAngle);
                }
            }

            double recoVX = -999999, recoVY = -999999, recoVZ = -999999;
            int numRecoNeutrinos = 0;
            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    numRecoNeutrinos++;
                    recoVX = reco_neutrinoVX->at(recoNeut);
                    recoVY = reco_neutrinoVY->at(recoNeut);
                    recoVZ = reco_neutrinoVZ->at(recoNeut);
                }
            }

            size_t wSliceIdx_cached = 999999;
            bool sliceWeightValid_cached = false;

            if(DLCurrent == 5 && signal != 3){
                if(slice < reco_sliceID_weights->size() && reco_sliceID_weights->at(slice) == reco_sliceID->at(slice)){
                    wSliceIdx_cached = slice;
                    sliceWeightValid_cached = true;
                } else {
                    for(size_t ws = 0; ws < reco_sliceID_weights->size(); ++ws){
                        if(reco_sliceID_weights->at(ws) == reco_sliceID->at(slice)){
                            wSliceIdx_cached = ws;
                            sliceWeightValid_cached = true;
                            break;
                        }
                    }
                }

                for(int p = 0; p < NPARAMS_GENIE; p++){
                    if(genieParams[p].skipForNow) continue;

                    auto* vec2D = reco_slice_GENIEWeightVecs[p];

                    if(genieParams[p].origNUniv == 6 || genieParams[p].origNUniv == 1){
                        int rawN = genieParams[p].origNUniv;
                        sliceWeightCalls++;
                        bool valid = sliceWeightValid_cached && vec2D
                                     && wSliceIdx_cached < vec2D->size()
                                     && (int)vec2D->at(wSliceIdx_cached).size() == rawN;
                        if(!valid){
                            sliceWeightFallbacks++;
                            std::fill(sliceUnivWeight_genie[p].begin(), sliceUnivWeight_genie[p].end(), 1.0);
                            continue;
                        }

                        double wRaw = vec2D->at(wSliceIdx_cached).at(0);
                        const std::vector<float>& sigmaU = sigmaUMap[genieParams[p].mapKey];

                        for(int u = 0; u < NUNIV_GENIE; u++){
                            sliceUnivWeight_genie[p][u] = expandToPseudoMultisim(wRaw, rawN, sigmaU[u]);
                        }
                    } else {
                        int N = genieParams[p].nUniv;
                        for(int u = 0; u < N; u++){
                            sliceUnivWeight_genie[p][u] = getGenieSliceWeight(vec2D, wSliceIdx_cached, u, sliceWeightValid_cached, N);
                        }
                    }
                }
            } else {
                for(int p = 0; p < NPARAMS_GENIE; p++){
                    if(genieParams[p].skipForNow) continue;
                    std::fill(sliceUnivWeight_genie[p].begin(), sliceUnivWeight_genie[p].end(), 1.0);
                }
            }

            if(DLCurrent == 5 && pfp10cm_PCAAngle != -999999){

                int binStudy = angleToBin_study(pfp10cm_PCAAngle*TMath::RadToDeg());

                if(binStudy != -1){

                    bool cutPass[9];
                    cutPass[0] = !(numPFPs0Cut == 1 && numPFPsSlice == 0);
                    cutPass[1] = !(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0);
                    cutPass[2] = !(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high));
                    cutPass[3] = !(FVCut == 1 && !(recoVX < FVCut_xHigh && recoVX > FVCut_xLow && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow && recoVZ < FVCut_zHigh));
                    cutPass[4] = !(primaryPFPCut == 1 && numPrimaryPFPs10Slice != primaryPFPCutValue);
                    cutPass[5] = !(razzledPDG11Cut == 1 && ((highestEnergyPFP.razzledPDG11 > razzled11High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG11 < razzled11Low_highestEnergyPFP)));
                    cutPass[6] = !(razzledPDG211Cut == 1 && ((highestEnergyPFP.razzledPDG211 > razzled211High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG211 < razzled211Low_highestEnergyPFP)));
                    cutPass[7] = !(ETheta2Cut == 1 && ((highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) > ETheta2High_highestEnergyPFP || (highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) < ETheta2Low_highestEnergyPFP));
                    cutPass[8] = !(dEdxCut == 1 && (highestEnergyPFP.bestPlanedEdx > dEdxHigh_highestEnergyPFP || highestEnergyPFP.bestPlanedEdx < dEdxLow_highestEnergyPFP));

                    int stageReached = 0;
                    for(int c = 0; c < 9; c++){
                        if(cutPass[c]) stageReached = c+1;
                        else break;
                    }

                    bool isSignal = (sliceCategoryPlottingMacro == 1 && signal == 1);
                    bool isBackground = ((sliceCategoryPlottingMacro == 2 && signal == 1) ||
                                          (sliceCategoryPlottingMacro == 0 && signal != 4) ||
                                          sliceCategoryPlottingMacro == 3 || sliceCategoryPlottingMacro == 4 ||
                                          sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6);

                    for(int s = 0; s <= stageReached; s++){
                        for(int catIdx = 0; catIdx < NCATS; catIdx++){
                            if(catIdx == 1 && !isSignal) continue;
                            if(catIdx == 2 && !isBackground) continue;

                            CVnominal[CVidx(s,catIdx,binStudy)] += weight;

                            for(int p = 0; p < NPARAMS_GENIE; p++){
                                if(genieParams[p].skipForNow) continue;
                                for(int u = 0; u < NUNIV_GENIE; u++){
                                    UNIVsum[UNIVidx(s,catIdx,p,u,binStudy)] += weight * sliceUnivWeight_genie[p][u];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "GENIE slice weight fallback: " << sliceWeightFallbacks << "/" << sliceWeightCalls << " (" << (sliceWeightCalls ? 100.0*sliceWeightFallbacks/sliceWeightCalls : 0.) << "%)" << std::endl;

    std::cout << "Building " << NSTAGES*NCATS*nActiveGenieParams << " cut-stage/GENIE-parameter covariance matrices (plus " << NSTAGES*NCATS << " summed totals)..." << std::endl;

    fOut->cd();
    TDirectory *dirStage = fOut->mkdir("covariance_angle_GENIE_byCutStage");

    std::vector<TMatrixDSym> finalCovSum(NCATS, TMatrixDSym(NBINS_STUDY));
    std::vector<TMatrixDSym> finalCorrSum(NCATS, TMatrixDSym(NBINS_STUDY));

    for(int s = 0; s < NSTAGES; s++){
        TDirectory *dS = dirStage->mkdir(stageNames[s].c_str());
        for(int catIdx = 0; catIdx < NCATS; catIdx++){
            TDirectory *dC = dS->mkdir(catNames[catIdx].c_str());

            std::vector<double> cv(NBINS_STUDY);
            for(int b = 0; b < NBINS_STUDY; b++) cv[b] = CVnominal[CVidx(s,catIdx,b)];

            TMatrixDSym cov_sumIndividual(NBINS_STUDY);
            cov_sumIndividual.Zero();

            for(int p = 0; p < NPARAMS_GENIE; p++){
                if(genieParams[p].skipForNow) continue;

                TDirectory *dP = dC->mkdir(genieParams[p].shortName.c_str());
                dP->cd();

                TMatrixDSym cov(NBINS_STUDY), corr(NBINS_STUDY);
                for(int i = 0; i < NBINS_STUDY; i++){
                    for(int j = 0; j < NBINS_STUDY; j++){
                        double sum = 0.0;
                        for(int u = 0; u < NUNIV_GENIE; u++){
                            double di = UNIVsum[UNIVidx(s,catIdx,p,u,i)] - cv[i];
                            double dj = UNIVsum[UNIVidx(s,catIdx,p,u,j)] - cv[j];
                            sum += di*dj;
                        }
                        cov(i,j) = sum / NUNIV_GENIE;
                    }
                }
                for(int i = 0; i < NBINS_STUDY; i++){
                    for(int j = 0; j < NBINS_STUDY; j++){
                        double denom = std::sqrt(cov(i,i)*cov(j,j));
                        corr(i,j) = (denom > 0) ? cov(i,j)/denom : 0.0;
                    }
                }

                for(int i = 0; i < NBINS_STUDY; i++)
                    for(int j = 0; j < NBINS_STUDY; j++)
                        cov_sumIndividual(i,j) += cov(i,j);

                std::string tag = stageNames[s] + "_" + catNames[catIdx] + "_" + genieParams[p].shortName;
                TH1D *hCV   = new TH1D(("h_CV_"+tag).c_str(), "", NBINS_STUDY, angleLow_study, angleHigh_study);
                TH2D *hCov  = new TH2D(("h_cov_"+tag).c_str(), "", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
                TH2D *hCorr = new TH2D(("h_corr_"+tag).c_str(), "", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);

                for(int b = 0; b < NBINS_STUDY; b++) hCV->SetBinContent(b+1, cv[b]);
                for(int i = 0; i < NBINS_STUDY; i++){
                    for(int j = 0; j < NBINS_STUDY; j++){
                        hCov->SetBinContent(i+1, j+1, cov(i,j));
                        hCorr->SetBinContent(i+1, j+1, corr(i,j));
                    }
                }

                hCV->Write();
                hCov->Write();
                hCorr->Write();

                delete hCV;
                delete hCov;
                delete hCorr;
            }

            TMatrixDSym corr_sumIndividual(NBINS_STUDY);
            for(int i = 0; i < NBINS_STUDY; i++){
                for(int j = 0; j < NBINS_STUDY; j++){
                    double denom = std::sqrt(cov_sumIndividual(i,i) * cov_sumIndividual(j,j));
                    corr_sumIndividual(i,j) = (denom > 0) ? cov_sumIndividual(i,j) / denom : 0.0;
                }
            }

            dC->cd();
            std::string tagSum = stageNames[s] + "_" + catNames[catIdx] + "_sumOfIndividualParams_GENIE";
            TH1D *hCVsum = new TH1D(("h_CV_"+tagSum).c_str(), "", NBINS_STUDY, angleLow_study, angleHigh_study);
            TH2D *hCovSum  = new TH2D(("h_cov_"+tagSum).c_str(), "", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
            TH2D *hCorrSum = new TH2D(("h_corr_"+tagSum).c_str(), "", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
            for(int b = 0; b < NBINS_STUDY; b++) hCVsum->SetBinContent(b+1, cv[b]);
            for(int i = 0; i < NBINS_STUDY; i++){
                for(int j = 0; j < NBINS_STUDY; j++){
                    hCovSum->SetBinContent(i+1, j+1, cov_sumIndividual(i,j));
                    hCorrSum->SetBinContent(i+1, j+1, corr_sumIndividual(i,j));
                }
            }
            hCVsum->Write();
            hCovSum->Write();
            hCorrSum->Write();
            delete hCVsum;
            delete hCovSum;
            delete hCorrSum;

            if(s == NSTAGES-1){
                finalCovSum[catIdx]  = cov_sumIndividual;
                finalCorrSum[catIdx] = corr_sumIndividual;
            }
        }
    }
    fOut->cd();
    std::cout << "Done building cut-stage/GENIE-parameter matrices." << std::endl;

    std::vector<double> cvFinal_all(NBINS_STUDY), cvFinal_signal(NBINS_STUDY), cvFinal_back(NBINS_STUDY);
    for(int b = 0; b < NBINS_STUDY; b++){
        cvFinal_all[b]    = CVnominal[CVidx(NSTAGES-1, 0, b)];
        cvFinal_signal[b] = CVnominal[CVidx(NSTAGES-1, 1, b)];
        cvFinal_back[b]   = CVnominal[CVidx(NSTAGES-1, 2, b)];
    }

    TH1D* h_angle_CV_genie        = new TH1D("h_angle_CV_genie", "Reconstructed Recoil Angle (Nominal, After All Cuts);Angle [deg];Weighted Events", NBINS_STUDY, angleLow_study, angleHigh_study);
    TH1D* h_angle_signal_CV_genie = new TH1D("h_angle_signal_CV_genie", "Reconstructed Recoil Angle (Nominal, After All Cuts, Signal Only);Angle [deg];Weighted Events", NBINS_STUDY, angleLow_study, angleHigh_study);
    TH1D* h_angle_back_CV_genie   = new TH1D("h_angle_back_CV_genie", "Reconstructed Recoil Angle (Nominal, After All Cuts, Background Only);Angle [deg];Weighted Events", NBINS_STUDY, angleLow_study, angleHigh_study);
    for(int b = 0; b < NBINS_STUDY; b++){
        h_angle_CV_genie->SetBinContent(b+1, cvFinal_all[b]);
        h_angle_signal_CV_genie->SetBinContent(b+1, cvFinal_signal[b]);
        h_angle_back_CV_genie->SetBinContent(b+1, cvFinal_back[b]);
    }

    TH2D* h_cov_angle_genie        = new TH2D("h_covMatrix_angle_genie", "Total GENIE Systematic Covariance Matrix (Recoil Angle);Angle bin;Angle bin", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
    TH2D* h_corr_angle_genie       = new TH2D("h_corrMatrix_angle_genie", "Total GENIE Systematic Correlation Matrix (Recoil Angle);Angle bin;Angle bin", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
    TH2D* h_cov_angle_signal_genie = new TH2D("h_covMatrix_angle_signal_genie", "Total GENIE Systematic Covariance Matrix (Recoil Angle, Signal Only);Angle bin;Angle bin", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
    TH2D* h_corr_angle_signal_genie= new TH2D("h_corrMatrix_angle_signal_genie", "Total GENIE Systematic Correlation Matrix (Recoil Angle, Signal Only);Angle bin;Angle bin", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
    TH2D* h_cov_angle_back_genie   = new TH2D("h_covMatrix_angle_back_genie", "Total GENIE Systematic Covariance Matrix (Recoil Angle, Background Only);Angle bin;Angle bin", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);
    TH2D* h_corr_angle_back_genie  = new TH2D("h_corrMatrix_angle_back_genie", "Total GENIE Systematic Correlation Matrix (Recoil Angle, Background Only);Angle bin;Angle bin", NBINS_STUDY, angleLow_study, angleHigh_study, NBINS_STUDY, angleLow_study, angleHigh_study);

    for(int i = 0; i < NBINS_STUDY; i++){
        for(int j = 0; j < NBINS_STUDY; j++){
            h_cov_angle_genie->SetBinContent(i+1, j+1, finalCovSum[0](i,j));
            h_corr_angle_genie->SetBinContent(i+1, j+1, finalCorrSum[0](i,j));
            h_cov_angle_signal_genie->SetBinContent(i+1, j+1, finalCovSum[1](i,j));
            h_corr_angle_signal_genie->SetBinContent(i+1, j+1, finalCorrSum[1](i,j));
            h_cov_angle_back_genie->SetBinContent(i+1, j+1, finalCovSum[2](i,j));
            h_corr_angle_back_genie->SetBinContent(i+1, j+1, finalCorrSum[2](i,j));
        }
    }

    TCanvas *c_cov = new TCanvas("c_cov_angle_genie", "", 800, 700);
    c_cov->SetRightMargin(0.15);
    h_cov_angle_genie->Draw("COLZ");
    c_cov->SaveAs((base_path + "covMatrix_angle_genie_total.pdf").c_str());

    TCanvas *c_corr = new TCanvas("c_corr_angle_genie", "", 800, 700);
    c_corr->SetRightMargin(0.15);
    h_corr_angle_genie->SetMinimum(-1.0);
    h_corr_angle_genie->SetMaximum(1.0);
    h_corr_angle_genie->Draw("COLZ");
    c_corr->SaveAs((base_path + "corrMatrix_angle_genie_total.pdf").c_str());

    TCanvas *c_cov_signal = new TCanvas("c_cov_angle_signal_genie", "", 800, 700);
    c_cov_signal->SetRightMargin(0.15);
    h_cov_angle_signal_genie->Draw("COLZ");
    c_cov_signal->SaveAs((base_path + "covMatrix_angle_signal_genie_total.pdf").c_str());

    TCanvas *c_corr_signal = new TCanvas("c_corr_angle_signal_genie", "", 800, 700);
    c_corr_signal->SetRightMargin(0.15);
    h_corr_angle_signal_genie->SetMinimum(-1.0);
    h_corr_angle_signal_genie->SetMaximum(1.0);
    h_corr_angle_signal_genie->Draw("COLZ");
    c_corr_signal->SaveAs((base_path + "corrMatrix_angle_signal_genie_total.pdf").c_str());

    TCanvas *c_cov_back = new TCanvas("c_cov_angle_back_genie", "", 800, 700);
    c_cov_back->SetRightMargin(0.15);
    h_cov_angle_back_genie->Draw("COLZ");
    c_cov_back->SaveAs((base_path + "covMatrix_angle_back_genie_total.pdf").c_str());

    TCanvas *c_corr_back = new TCanvas("c_corr_angle_back_genie", "", 800, 700);
    c_corr_back->SetRightMargin(0.15);
    h_corr_angle_back_genie->SetMinimum(-1.0);
    h_corr_angle_back_genie->SetMaximum(1.0);
    h_corr_angle_back_genie->Draw("COLZ");
    c_corr_back->SaveAs((base_path + "corrMatrix_angle_back_genie_total.pdf").c_str());

    fOut->cd();
    TDirectory *dirCov = fOut->GetDirectory("covariance_angle_GENIE");
    if(!dirCov) dirCov = fOut->mkdir("covariance_angle_GENIE");
    dirCov->cd();
    h_angle_CV_genie->Write();
    h_cov_angle_genie->Write();
    h_corr_angle_genie->Write();
    h_angle_signal_CV_genie->Write();
    h_cov_angle_signal_genie->Write();
    h_corr_angle_signal_genie->Write();
    h_angle_back_CV_genie->Write();
    h_cov_angle_back_genie->Write();
    h_corr_angle_back_genie->Write();
    fOut->cd();

    delete h_angle_CV_genie;
    delete h_cov_angle_genie;
    delete h_corr_angle_genie;
    delete c_cov;
    delete c_corr;

    delete h_angle_signal_CV_genie;
    delete h_cov_angle_signal_genie;
    delete h_corr_angle_signal_genie;
    delete c_cov_signal;
    delete c_corr_signal;

    delete h_angle_back_CV_genie;
    delete h_cov_angle_back_genie;
    delete h_corr_angle_back_genie;
    delete c_cov_back;
    delete c_corr_back;

    fOut->Write();
    fOut->Close();
    delete fOut;
}
