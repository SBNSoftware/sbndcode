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
    double signalCurrent = 0;
    double BNBCurrent = 0;
    double cosmicsCurrent = 0;
    double signalUboone = 0;
    double BNBUboone = 0;
    double cosmicsUboone = 0;
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

void nuECovarainceMatrixFlux_macro(){

    std::string cutsApplied = "allCuts";
    std::string base_path = "/nashome/c/coackley/systPlotsFluxCov18July_" + cutsApplied + "/";

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

    TFile *fOut = new TFile("/exp/sbnd/data/users/coackley/selectionCovarianceMatrixFlux_18July.root", "RECREATE");
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
    weights.BNBNuE = targetPOT /POTBNBNuE_notMissing;
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

    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_horncurrent = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_expskin = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_pioninexsec = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_pionqexsec = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_piontotxsec = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_nucleoninexsec = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_nucleonqexsec = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_nucleontotxsec = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_kplus = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_kmin = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_kzero = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_piplus = nullptr;
    std::vector<std::vector<float>> *reco_sliceMCTruthFlux_weight_piminus = nullptr;

    weightsTree->SetBranchAddress("nuEScatterTrueVX_weightTree", &nuEScatterTrueVX_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVY_weightTree", &nuEScatterTrueVY_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVZ_weightTree", &nuEScatterTrueVZ_weights);
    weightsTree->SetBranchAddress("reco_sliceID_weightTree", &reco_sliceID_weights);

    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_horncurrent", &reco_sliceMCTruthFlux_weight_horncurrent);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_expskin", &reco_sliceMCTruthFlux_weight_expskin);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_pioninexsec", &reco_sliceMCTruthFlux_weight_pioninexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_pionqexsec", &reco_sliceMCTruthFlux_weight_pionqexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_piontotxsec", &reco_sliceMCTruthFlux_weight_piontotxsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_nucleoninexsec", &reco_sliceMCTruthFlux_weight_nucleoninexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_nucleonqexsec", &reco_sliceMCTruthFlux_weight_nucleonqexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_nucleontotxsec", &reco_sliceMCTruthFlux_weight_nucleontotxsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_kplus", &reco_sliceMCTruthFlux_weight_kplus);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_kmin", &reco_sliceMCTruthFlux_weight_kmin);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_kzero", &reco_sliceMCTruthFlux_weight_kzero);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_piplus", &reco_sliceMCTruthFlux_weight_piplus);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_piminus", &reco_sliceMCTruthFlux_weight_piminus);

    const int NUNIV = 1000;
    long long sliceWeightFallbacks = 0, sliceWeightCalls = 0;

    auto getSliceWeight = [&](std::vector<std::vector<float>>* vec, size_t sliceIdx, int u, bool wFound) -> double {
        sliceWeightCalls++;
        if(!wFound || !vec || sliceIdx >= vec->size() || (int)vec->at(sliceIdx).size() != NUNIV){
            sliceWeightFallbacks++;
            return 1.0;
        }
        return vec->at(sliceIdx).at(u);
    };

// ================= ADDITIONAL: cut-stage x parameter covariance study =================

    const int NSTAGES = 10; // 0 = no cuts, 1..9 = cumulative cuts applied in order below
    const std::vector<std::string> stageNames = {
        "s0_noCuts","s1_numPFPs0","s2_numRecoNeutrinos","s3_CRUMBS","s4_FV",
        "s5_primaryPFP","s6_razzledPDG11","s7_razzledPDG211","s8_ETheta2","s9_dEdx"
    };

    const int NPARAMS = 14; // 0 = all 13 combined, 1-13 = individual flux params
    const std::vector<std::string> paramNames = {
        "combined","horncurrent","expskin","pioninexsec","pionqexsec","piontotxsec",
        "nucleoninexsec","nucleonqexsec","nucleontotxsec","kplus","kmin","kzero",
        "piplus","piminus"
    };
    // order must match paramNames[1..13] above
    std::vector<std::vector<std::vector<float>>*> paramPtrs = {
        reco_sliceMCTruthFlux_weight_horncurrent, reco_sliceMCTruthFlux_weight_expskin,
        reco_sliceMCTruthFlux_weight_pioninexsec, reco_sliceMCTruthFlux_weight_pionqexsec,
        reco_sliceMCTruthFlux_weight_piontotxsec, reco_sliceMCTruthFlux_weight_nucleoninexsec,
        reco_sliceMCTruthFlux_weight_nucleonqexsec, reco_sliceMCTruthFlux_weight_nucleontotxsec,
        reco_sliceMCTruthFlux_weight_kplus, reco_sliceMCTruthFlux_weight_kmin,
        reco_sliceMCTruthFlux_weight_kzero, reco_sliceMCTruthFlux_weight_piplus,
        reco_sliceMCTruthFlux_weight_piminus
    };

    const int NCATS = 3; // 0 = all, 1 = signal, 2 = background
    const std::vector<std::string> catNames = {"all","signal","background"};

    const int NUNIV_STUDY  = 1000;
    const int NBINS_STUDY  = 15;
    double angleLow_study  = 0.0;
    double angleHigh_study = 20.0;

    std::vector<double> CVsum((size_t)NSTAGES*NCATS*NPARAMS*NBINS_STUDY, 0.0);
    std::vector<double> UNIVsum((size_t)NSTAGES*NCATS*NPARAMS*NUNIV_STUDY*NBINS_STUDY, 0.0);

    auto CVidx = [&](int stage,int cat,int param,int bin)->size_t{
        return (size_t)(((stage*NCATS + cat)*NPARAMS + param)*NBINS_STUDY + bin);
    };
    auto UNIVidx = [&](int stage,int cat,int param,int u,int bin)->size_t{
        return ((((size_t)stage*NCATS + cat)*NPARAMS + param)*NUNIV_STUDY + u)*NBINS_STUDY + bin;
    };
    auto angleToBin_study = [&](double angleDeg)->int{
        if(angleDeg < angleLow_study || angleDeg >= angleHigh_study) return -1;
        int b = (int)((angleDeg - angleLow_study)/(angleHigh_study - angleLow_study)*NBINS_STUDY);
        if(b < 0) b = 0;
        if(b >= NBINS_STUDY) b = NBINS_STUDY - 1;
        return b;
    };

    double paramUnivWeight_study[13][1000];

    const int NANGLEBINS = 60;
    double angleBinLow  = 0.0;
    double angleBinHigh = 30.0;

    TH1D* h_angle_CV = new TH1D("h_angle_CV", "Reconstructed Recoil Angle (Nominal, After All Cuts);Angle [deg];Weighted Events", NANGLEBINS, angleBinLow, angleBinHigh);
    TH1D* h_angle_signal_CV = new TH1D("h_angle_signal_CV", "Reconstructed Recoil Angle (Nominal, After All Cuts, Signal Only);Angle [deg];Weighted Events", NANGLEBINS, angleBinLow, angleBinHigh);
    TH1D* h_angle_back_CV = new TH1D("h_angle_back_CV", "Reconstructed Recoil Angle (Nominal, After All Cuts, Background Only);Angle [deg];Weighted Events", NANGLEBINS, angleBinLow, angleBinHigh);

    std::vector<TH1D*> h_angle_univ(NUNIV, nullptr);
    std::vector<TH1D*> h_angle_signal_univ(NUNIV, nullptr);
    std::vector<TH1D*> h_angle_back_univ(NUNIV, nullptr);

    for(int u = 0; u < NUNIV; u++){
        h_angle_univ[u] = new TH1D(Form("h_angle_univ_%d", u), "", NANGLEBINS, angleBinLow, angleBinHigh);
        h_angle_signal_univ[u] = new TH1D(Form("h_angle_signal_univ_%d", u), "", NANGLEBINS, angleBinLow, angleBinHigh);
        h_angle_back_univ[u] = new TH1D(Form("h_angle_back_univ_%d", u), "", NANGLEBINS, angleBinLow, angleBinHigh);
    }

    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);
        weightsTree->GetEntry(e);

        const double epsCheck = 1e-6;
        if(std::abs(nuEScatterTrueVX - nuEScatterTrueVX_weights) > epsCheck || std::abs(nuEScatterTrueVY - nuEScatterTrueVY_weights) > epsCheck || std::abs(nuEScatterTrueVZ - nuEScatterTrueVZ_weights) > epsCheck){
            std::cerr << "ERROR: entry " << e << " misaligned between NuE and NuEWeights trees!" << std::endl;
            continue;
        }

        recoilElectron_struct recoilElectron;
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
            std::vector<double> sliceUnivWeight_combined(NUNIV, 1.0);

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

                for(int u = 0; u < NUNIV; u++){
                    double wHorn = getSliceWeight(reco_sliceMCTruthFlux_weight_horncurrent, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wExp = getSliceWeight(reco_sliceMCTruthFlux_weight_expskin, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wKplus = getSliceWeight(reco_sliceMCTruthFlux_weight_kplus, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wKmin = getSliceWeight(reco_sliceMCTruthFlux_weight_kmin, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wKzero = getSliceWeight(reco_sliceMCTruthFlux_weight_kzero, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wNinex = getSliceWeight(reco_sliceMCTruthFlux_weight_nucleoninexsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wNqex = getSliceWeight(reco_sliceMCTruthFlux_weight_nucleonqexsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wNtotx = getSliceWeight(reco_sliceMCTruthFlux_weight_nucleontotxsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPiminus = getSliceWeight(reco_sliceMCTruthFlux_weight_piminus, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPinex = getSliceWeight(reco_sliceMCTruthFlux_weight_pioninexsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPiqex = getSliceWeight(reco_sliceMCTruthFlux_weight_pionqexsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPitotx = getSliceWeight(reco_sliceMCTruthFlux_weight_piontotxsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPiplus = getSliceWeight(reco_sliceMCTruthFlux_weight_piplus, wSliceIdx_cached, u, sliceWeightValid_cached);
                    sliceUnivWeight_combined[u] = wHorn*wExp*wKplus*wKmin*wKzero*wNinex*wNqex*wNtotx*wPiminus*wPinex*wPiqex*wPitotx*wPiplus;
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

                    for(int pi = 0; pi < 13; pi++){
                        for(int u = 0; u < NUNIV_STUDY; u++){
                            paramUnivWeight_study[pi][u] = getSliceWeight(paramPtrs[pi], wSliceIdx_cached, u, sliceWeightValid_cached);
                        }
                    }

                    for(int s = 0; s <= stageReached; s++){
                        for(int catIdx = 0; catIdx < NCATS; catIdx++){
                            if(catIdx == 1 && !isSignal) continue;
                            if(catIdx == 2 && !isBackground) continue;

                            CVsum[CVidx(s,catIdx,0,binStudy)] += weight;
                            for(int u = 0; u < NUNIV_STUDY; u++){
                                UNIVsum[UNIVidx(s,catIdx,0,u,binStudy)] += weight * sliceUnivWeight_combined[u];
                            }

                            for(int pi = 0; pi < 13; pi++){
                                CVsum[CVidx(s,catIdx,pi+1,binStudy)] += weight;
                                for(int u = 0; u < NUNIV_STUDY; u++){
                                    UNIVsum[UNIVidx(s,catIdx,pi+1,u,binStudy)] += weight * paramUnivWeight_study[pi][u];
                                }
                            }
                        }
                    }
                }
            }

            if(numPFPs0Cut == 1 && numPFPsSlice == 0) continue;

            if(numPFPs0Cut == 1 && numPFPsSlice == 0) continue;
            if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0) continue;
            if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high)) continue;
            if(FVCut == 1){
                if(!(recoVX < FVCut_xHigh && recoVX > FVCut_xLow && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow && recoVZ < FVCut_zHigh)) continue;
            }
            if(primaryPFPCut == 1 && numPrimaryPFPs10Slice != primaryPFPCutValue) continue;
            if(razzledPDG11Cut == 1 && ((highestEnergyPFP.razzledPDG11 > razzled11High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG11 < razzled11Low_highestEnergyPFP))) continue;
            if(razzledPDG211Cut == 1 && ((highestEnergyPFP.razzledPDG211 > razzled211High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG211 < razzled211Low_highestEnergyPFP))) continue;
            if(ETheta2Cut == 1 && ((highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) > ETheta2High_highestEnergyPFP || (highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) < ETheta2Low_highestEnergyPFP)) continue;
            if(dEdxCut == 1 && (highestEnergyPFP.bestPlanedEdx > dEdxHigh_highestEnergyPFP || highestEnergyPFP.bestPlanedEdx < dEdxLow_highestEnergyPFP)) continue;

            if(DLCurrent == 5 && pfp10cm_PCAAngle != -999999){
                h_angle_CV->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);
                if(sliceCategoryPlottingMacro == 1 && signal == 1) h_angle_signal_CV->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);
                else if((sliceCategoryPlottingMacro == 2 && signal == 1) || (sliceCategoryPlottingMacro == 0 && signal != 4) || sliceCategoryPlottingMacro == 3 || sliceCategoryPlottingMacro == 4 || sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6) h_angle_back_CV->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);

                for(int u = 0; u < NUNIV; u++){
                    double wComb = sliceUnivWeight_combined[u];
                    h_angle_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight * wComb);
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) h_angle_signal_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight * wComb);
                    else if((sliceCategoryPlottingMacro == 2 && signal == 1) || sliceCategoryPlottingMacro == 3 || sliceCategoryPlottingMacro == 4 || sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6) h_angle_back_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight * wComb);
                    else if(sliceCategoryPlottingMacro == 0 && signal != 4) h_angle_back_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);
                }
            }
        }
    }

    std::cout << "slice weight fallback: " << sliceWeightFallbacks << "/" << sliceWeightCalls << " (" << (sliceWeightCalls? 100.0*sliceWeightFallbacks/sliceWeightCalls : 0.) << "%)" << std::endl;

    std::cout << "Building " << NSTAGES*NCATS*NPARAMS << " cut-stage/parameter covariance matrices..." << std::endl;

    fOut->cd();
    TDirectory *dirStage = fOut->mkdir("covariance_angle_byCutStage");

    for(int s = 0; s < NSTAGES; s++){
        TDirectory *dS = dirStage->mkdir(stageNames[s].c_str());
        for(int catIdx = 0; catIdx < NCATS; catIdx++){
            TDirectory *dC = dS->mkdir(catNames[catIdx].c_str());
            for(int p = 0; p < NPARAMS; p++){
                TDirectory *dP = dC->mkdir(paramNames[p].c_str());
                dP->cd();

                std::vector<double> cv(NBINS_STUDY);
                for(int b = 0; b < NBINS_STUDY; b++) cv[b] = CVsum[CVidx(s,catIdx,p,b)];

                TMatrixDSym cov(NBINS_STUDY), corr(NBINS_STUDY);
                for(int i = 0; i < NBINS_STUDY; i++){
                    for(int j = 0; j < NBINS_STUDY; j++){
                        double sum = 0.0;
                        for(int u = 0; u < NUNIV_STUDY; u++){
                            double di = UNIVsum[UNIVidx(s,catIdx,p,u,i)] - cv[i];
                            double dj = UNIVsum[UNIVidx(s,catIdx,p,u,j)] - cv[j];
                            sum += di*dj;
                        }
                        cov(i,j) = sum / NUNIV_STUDY;
                    }
                }
                for(int i = 0; i < NBINS_STUDY; i++){
                    for(int j = 0; j < NBINS_STUDY; j++){
                        double denom = std::sqrt(cov(i,i)*cov(j,j));
                        corr(i,j) = (denom > 0) ? cov(i,j)/denom : 0.0;
                    }
                }

                std::string tag = stageNames[s] + "_" + catNames[catIdx] + "_" + paramNames[p];
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
        }
    }
    fOut->cd();
    std::cout << "Done building cut-stage/parameter matrices." << std::endl;

    TMatrixDSym covMatrix_angle(NANGLEBINS), corrMatrix_angle(NANGLEBINS);

    TMatrixDSym covMatrix_angle(NANGLEBINS), corrMatrix_angle(NANGLEBINS);
    TMatrixDSym covMatrix_angle_signal(NANGLEBINS), corrMatrix_angle_signal(NANGLEBINS);
    TMatrixDSym covMatrix_angle_back(NANGLEBINS), corrMatrix_angle_back(NANGLEBINS);

    std::vector<double> CV_binContent(NANGLEBINS), CV_signal_binContent(NANGLEBINS), CV_back_binContent(NANGLEBINS);
    for(int i = 0; i < NANGLEBINS; i++){
        CV_binContent[i] = h_angle_CV->GetBinContent(i+1);
        CV_signal_binContent[i] = h_angle_signal_CV->GetBinContent(i+1);
        CV_back_binContent[i] = h_angle_back_CV->GetBinContent(i+1);
    }

    for(int i = 0; i < NANGLEBINS; i++){
        for(int j = 0; j < NANGLEBINS; j++){
            double sum = 0.0, sum_signal = 0.0, sum_back = 0.0;
            for(int u = 0; u < NUNIV; u++){
                double di = h_angle_univ[u]->GetBinContent(i+1) - CV_binContent[i];
                double dj = h_angle_univ[u]->GetBinContent(j+1) - CV_binContent[j];
                sum += di * dj;

                double di_signal = h_angle_signal_univ[u]->GetBinContent(i+1) - CV_signal_binContent[i];
                double dj_signal = h_angle_signal_univ[u]->GetBinContent(j+1) - CV_signal_binContent[j];
                sum_signal += di_signal * dj_signal;

                double di_back = h_angle_back_univ[u]->GetBinContent(i+1) - CV_back_binContent[i];
                double dj_back = h_angle_back_univ[u]->GetBinContent(j+1) - CV_back_binContent[j];
                sum_back += di_back * dj_back;
            }
            covMatrix_angle(i,j) = sum / NUNIV;
            covMatrix_angle_signal(i,j) = sum_signal / NUNIV;
            covMatrix_angle_back(i,j) = sum_back / NUNIV;
        }
    }

    for(int i = 0; i < NANGLEBINS; i++){
        for(int j = 0; j < NANGLEBINS; j++){
            double denom = std::sqrt(covMatrix_angle(i,i) * covMatrix_angle(j,j));
            corrMatrix_angle(i,j) = (denom > 0) ? covMatrix_angle(i,j) / denom : 0.0;
            double denom_signal = std::sqrt(covMatrix_angle_signal(i,i) * covMatrix_angle_signal(j,j));
            corrMatrix_angle_signal(i,j) = (denom_signal > 0) ? covMatrix_angle_signal(i,j) / denom_signal : 0.0;
            double denom_back = std::sqrt(covMatrix_angle_back(i,i) * covMatrix_angle_back(j,j));
            corrMatrix_angle_back(i,j) = (denom_back > 0) ? covMatrix_angle_back(i,j) / denom_back : 0.0;
        }
    }

    TH2D* h_cov_angle = new TH2D("h_covMatrix_angle_flux", "Flux Systematic Covariance Matrix (Recoil Angle);Angle bin;Angle bin", NANGLEBINS, angleBinLow, angleBinHigh, NANGLEBINS, angleBinLow, angleBinHigh);
    TH2D* h_corr_angle = new TH2D("h_corrMatrix_angle_flux", "Flux Systematic Correlation Matrix (Recoil Angle);Angle bin;Angle bin", NANGLEBINS, angleBinLow, angleBinHigh, NANGLEBINS, angleBinLow, angleBinHigh);
    TH2D* h_cov_angle_signal = new TH2D("h_covMatrix_angle_signal_flux", "Flux Systematic Covariance Matrix (Recoil Angle, Signal Only);Angle bin;Angle bin", NANGLEBINS, angleBinLow, angleBinHigh, NANGLEBINS, angleBinLow, angleBinHigh);
    TH2D* h_corr_angle_signal = new TH2D("h_corrMatrix_angle_signal_flux", "Flux Systematic Correlation Matrix (Recoil Angle, Signal Only);Angle bin;Angle bin", NANGLEBINS, angleBinLow, angleBinHigh, NANGLEBINS, angleBinLow, angleBinHigh);
    TH2D* h_cov_angle_back = new TH2D("h_covMatrix_angle_back_flux", "Flux Systematic Covariance Matrix (Recoil Angle, Background Only);Angle bin;Angle bin", NANGLEBINS, angleBinLow, angleBinHigh, NANGLEBINS, angleBinLow, angleBinHigh);
    TH2D* h_corr_angle_back = new TH2D("h_corrMatrix_angle_back_flux", "Flux Systematic Correlation Matrix (Recoil Angle, Background Only);Angle bin;Angle bin", NANGLEBINS, angleBinLow, angleBinHigh, NANGLEBINS, angleBinLow, angleBinHigh);

    for(int i = 0; i < NANGLEBINS; i++){
        for(int j = 0; j < NANGLEBINS; j++){
            h_cov_angle->SetBinContent(i+1, j+1, covMatrix_angle(i,j));
            h_corr_angle->SetBinContent(i+1, j+1, corrMatrix_angle(i,j));
            h_cov_angle_signal->SetBinContent(i+1, j+1, covMatrix_angle_signal(i,j));
            h_corr_angle_signal->SetBinContent(i+1, j+1, corrMatrix_angle_signal(i,j));
            h_cov_angle_back->SetBinContent(i+1, j+1, covMatrix_angle_back(i,j));
            h_corr_angle_back->SetBinContent(i+1, j+1, corrMatrix_angle_back(i,j));
        }
    }

    TCanvas *c_cov = new TCanvas("c_cov_angle_flux", "", 800, 700);
    c_cov->SetRightMargin(0.15);
    h_cov_angle->Draw("COLZ");
    c_cov->SaveAs((base_path + "covMatrix_angle_flux.pdf").c_str());

    TCanvas *c_corr = new TCanvas("c_corr_angle_flux", "", 800, 700);
    c_corr->SetRightMargin(0.15);
    h_corr_angle->SetMinimum(-1.0);
    h_corr_angle->SetMaximum(1.0);
    h_corr_angle->Draw("COLZ");
    c_corr->SaveAs((base_path + "corrMatrix_angle_flux.pdf").c_str());

    TCanvas *c_cov_signal = new TCanvas("c_cov_angle_signal_flux", "", 800, 700);
    c_cov_signal->SetRightMargin(0.15);
    h_cov_angle_signal->Draw("COLZ");
    c_cov_signal->SaveAs((base_path + "covMatrix_angle_signal_flux.pdf").c_str());

    TCanvas *c_corr_signal = new TCanvas("c_corr_angle_signal_flux", "", 800, 700);
    c_corr_signal->SetRightMargin(0.15);
    h_corr_angle_signal->SetMinimum(-1.0);
    h_corr_angle_signal->SetMaximum(1.0);
    h_corr_angle_signal->Draw("COLZ");
    c_corr_signal->SaveAs((base_path + "corrMatrix_angle_signal_flux.pdf").c_str());

    TCanvas *c_cov_back = new TCanvas("c_cov_angle_back_flux", "", 800, 700);
    c_cov_back->SetRightMargin(0.15);
    h_cov_angle_back->Draw("COLZ");
    c_cov_back->SaveAs((base_path + "covMatrix_angle_back_flux.pdf").c_str());

    TCanvas *c_corr_back = new TCanvas("c_corr_angle_back_flux", "", 800, 700);
    c_corr_back->SetRightMargin(0.15);
    h_corr_angle_back->SetMinimum(-1.0);
    h_corr_angle_back->SetMaximum(1.0);
    h_corr_angle_back->Draw("COLZ");
    c_corr_back->SaveAs((base_path + "corrMatrix_angle_back_flux.pdf").c_str());

    fOut->cd();
    TDirectory *dirCov = fOut->GetDirectory("covariance_angle");
    if(!dirCov) dirCov = fOut->mkdir("covariance_angle");
    dirCov->cd();
    h_angle_CV->Write();
    h_cov_angle->Write();
    h_corr_angle->Write();
    h_angle_signal_CV->Write();
    h_cov_angle_signal->Write();
    h_corr_angle_signal->Write();
    h_angle_back_CV->Write();
    h_cov_angle_back->Write();
    h_corr_angle_back->Write();
    fOut->cd();

    delete h_angle_CV;
    delete h_cov_angle;
    delete h_corr_angle;
    delete c_cov;
    delete c_corr;
    for(int u = 0; u < NUNIV; u++) delete h_angle_univ[u];

    delete h_angle_signal_CV;
    delete h_cov_angle_signal;
    delete h_corr_angle_signal;
    delete c_cov_signal;
    delete c_corr_signal;
    for(int u = 0; u < NUNIV; u++) delete h_angle_signal_univ[u];

    delete h_angle_back_CV;
    delete h_cov_angle_back;
    delete h_corr_angle_back;
    delete c_cov_back;
    delete c_corr_back;
    for(int u = 0; u < NUNIV; u++) delete h_angle_back_univ[u];

    fOut->Write();
    fOut->Close();
    delete fOut;
}
