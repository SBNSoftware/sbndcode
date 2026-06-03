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
#include <utility>
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

struct eventKey_struct{
    UInt_t runID;
    UInt_t subRunID;
    UInt_t eventID;
    int signal;
    int DLCurrent;

    bool operator == (const eventKey_struct& other) const {
        return runID == other.runID && subRunID == other.subRunID && eventID == other.eventID && signal == other.signal && DLCurrent == other.DLCurrent;
    }
};

struct eventKeyHash_struct{
    std::size_t operator()(const eventKey_struct& k) const{
        std::size_t h = 0;

        h ^= std::hash<UInt_t>{}(k.runID);
        h ^= std::hash<UInt_t>{}(k.subRunID) << 1;
        h ^= std::hash<UInt_t>{}(k.eventID) << 2;
        h ^= std::hash<int>{}(k.signal) << 3;
        h ^= std::hash<int>{}(k.DLCurrent) << 4;

        return h;
    }
};

struct weights_struct{
    double signalNuE = 0;
    double BNBNuE = 0;
    double cosmicsNuE = 0;
    double signalCurrent = 0;
    double BNBCurrent = 0;
    double cosmicsCurrent = 0;
    double signalUboone = 0;
    double BNBUboone = 0;
    double cosmicsUboone = 0;
};

void nuESystWeightMatching_macro(){
    // Load in the NuE and SubRun TTrees
    TFile *fNuE = TFile::Open("/exp/sbnd/app/users/coackley/nue/srcs/sbndcode/sbndcode/nue/merged_noWeights.root");
    if(!fNuE){
        std::cerr << "Error opening the NuE TTree file" << std::endl;
        return;
    }

    TDirectory *dirNuE = (TDirectory*)fNuE->Get("ana");
    if(!dirNuE){
        std::cerr << "Directory 'ana' not found in NuE file" << std::endl;
        return;
    }

    TTree *tree = (TTree*)dirNuE->Get("NuE");
    if(!tree){
        std::cerr << "NuE TTree not found" << std::endl;
        return;
    }

    TTree *subRunTree = (TTree*)dirNuE->Get("SubRun");
    if(!subRunTree){
        std::cerr << "SubRun TTree not found" << std::endl;
        return;
    }

    // Load in the NuEWeights TTree
    TFile *fNuEWeights = TFile::Open("/exp/sbnd/app/users/coackley/nue/srcs/sbndcode/sbndcode/nue/merged_weights.root");
    
    if(!fNuEWeights){
        std::cerr << "Error opening the NuEWeights TTree file" << std::endl;
        return;
    }

    TDirectory *dirNuEWeights = (TDirectory*)fNuEWeights->Get("ana");
    if(!dirNuEWeights){
        std::cerr << "Directory 'ana' not found in the NuEWeights file" << std::endl;
        return;
    }

    TTree *weightsTree = (TTree*)dirNuEWeights->Get("NuEWeights");
    if(!weightsTree){
        std::cerr << "NuEWeights TTree not found" << std::endl;
        return;
    }

    // SubRun Tree branch variables
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

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBUboone;

    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;

    double cosmicSpillsSumNuE = 0;
    double BNBSpillsSumNuE = 0;
    double NuESpillsSumNuE = 0;

    double POTSignalNuE_notMissing = 0;
    double POTBNBNuE_notMissing = 0;
    
    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 2 && subRunDLCurrent == 5) BNBSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 1 && subRunDLCurrent == 5) NuESpillsSumNuE += subRunNumGenEvents;

        std::pair<unsigned int, unsigned int> key = std::make_pair(subRunRun, subRunNumber);

        if(subRunSignal == 1){
            if(subRunDLCurrent == 5 && seenSubRunsSignalNuE.find(key) == seenSubRunsSignalNuE.end()){
                totalPOTSignalNuE += subRunPOT;
                seenSubRunsSignalNuE.insert(key);
            }
            
            if(subRunDLCurrent == 5) POTSignalNuE_notMissing += subRunPOT;
                
        } else if(subRunSignal == 2){
            if(subRunDLCurrent == 5 && seenSubRunsBNBNuE.find(key) == seenSubRunsBNBNuE.end()){
                totalPOTBNBNuE += subRunPOT;
                seenSubRunsBNBNuE.insert(key);
            }
            
            if(subRunDLCurrent == 5) POTBNBNuE_notMissing += subRunPOT;
        }
    }

    double targetPOT = 1e21;
    double targetSpills = (targetPOT/(5e12));

    double BNBScaledSpills_NuE = ((targetPOT/POTBNBNuE_notMissing) * BNBSpillsSumNuE);
    double SignalScaledSpills_NuE = ((targetPOT/POTSignalNuE_notMissing) * NuESpillsSumNuE);

    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    weights_struct weights;
    weights.signalNuE = targetPOT / POTSignalNuE_notMissing;
    weights.BNBNuE = targetPOT /POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;

    std::cout << "Weights DLNu+E: BNB = " << weights.BNBNuE << ", Signal = " << weights.signalNuE << ", Intime Cosmics = " << weights.cosmicsNuE << std::endl;


    // NuE Tree branch variables
    UInt_t eventID, runID, subRunID;
    int nuEScatter;
    double nuEScatterTrueVX, nuEScatterTrueVY, nuEScatterTrueVZ;
    double DLCurrent, signal;

    std::vector<double> *truth_recoilElectronPDG = nullptr;  
    std::vector<double> *truth_recoilElectronVX = nullptr;  
    std::vector<double> *truth_recoilElectronVY = nullptr;  
    std::vector<double> *truth_recoilElectronVZ = nullptr;  
    std::vector<double> *truth_recoilElectronPX = nullptr;  
    std::vector<double> *truth_recoilElectronPY = nullptr;  
    std::vector<double> *truth_recoilElectronPZ = nullptr;  
    std::vector<double> *truth_recoilElectronEnergy = nullptr;  
    std::vector<double> *truth_recoilElectronAngle = nullptr;  
    std::vector<double> *truth_recoilElectronETheta2 = nullptr;  
    std::vector<double> *truth_recoilElectronDX = nullptr;  
    std::vector<double> *truth_recoilElectronDY = nullptr;  
    std::vector<double> *truth_recoilElectronDZ = nullptr;  
    
    std::vector<double> *reco_sliceID = nullptr;  
    std::vector<double> *reco_sliceCompleteness = nullptr;  
    std::vector<double> *reco_slicePurity = nullptr;  
    std::vector<double> *reco_sliceScore = nullptr;  
    std::vector<double> *reco_sliceCategory = nullptr;  
    std::vector<double> *reco_sliceInteraction = nullptr;  
    std::vector<double> *reco_sliceTrueVX = nullptr;  
    std::vector<double> *reco_sliceTrueVY = nullptr;  
    std::vector<double> *reco_sliceTrueVZ = nullptr;  
    std::vector<double> *reco_sliceNumHits = nullptr;  
    std::vector<double> *reco_sliceNumHitsTruthMatched = nullptr;  
    std::vector<double> *reco_sliceNumTruthHits = nullptr;  
    std::vector<double> *reco_sliceOrigin = nullptr; 
    std::vector<double> *reco_sliceTrueCCNC = nullptr;  
    std::vector<double> *reco_sliceTrueNeutrinoType = nullptr;

    std::vector<double> *truth_particleSliceID = nullptr;  
    std::vector<double> *truth_particlePrimary = nullptr;  
    std::vector<double> *truth_particleVX = nullptr;  
    std::vector<double> *truth_particleVY = nullptr;  
    std::vector<double> *truth_particleVZ = nullptr;  
    std::vector<double> *truth_particlePDG = nullptr;  
    std::vector<double> *truth_particleTrackID = nullptr;  
    std::vector<double> *truth_particleMother = nullptr;  
    std::vector<double> *truth_particleStatusCode = nullptr;  

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
    std::vector<double> *reco_particleNumHitsTruthMatched = nullptr;
    std::vector<double> *reco_particleNumTruthHits = nullptr;
    std::vector<double> *reco_particleClearCosmic = nullptr;  
    std::vector<double> *reco_particlePlane0dEdx = nullptr;  
    std::vector<double> *reco_particlePlane1dEdx = nullptr;  
    std::vector<double> *reco_particlePlane2dEdx = nullptr;  
    std::vector<double> *reco_particleBestPlanedEdx = nullptr;  
    std::vector<double> *reco_particleRazzledPDG11 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG13 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG22 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG211 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG2212 = nullptr;  
    std::vector<double> *reco_particleRazzledBestPDG = nullptr;  
    std::vector<double> *reco_particleShowerLength = nullptr;  
    std::vector<double> *reco_particleShowerOpenAngle = nullptr;  
    std::vector<double> *reco_particleShowerBestPlaneEnergy = nullptr;  
    std::vector<double> *reco_particleTrueVX = nullptr;  
    std::vector<double> *reco_particleTrueVY = nullptr;  
    std::vector<double> *reco_particleTrueVZ = nullptr;  
    std::vector<double> *reco_particleTrueEndX = nullptr;  
    std::vector<double> *reco_particleTrueEndY = nullptr;  
    std::vector<double> *reco_particleTrueEndZ = nullptr;  
  
    std::vector<double> *reco_neutrinoID = nullptr;
    std::vector<double> *reco_neutrinoPDG = nullptr;
    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;
    
    std::vector<double> *angleRecalculationPCASlice_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dz = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dz = nullptr;
    
    std::vector<double> *angleRecalculationPCAPFP_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dz = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dz = nullptr;

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
    tree->SetBranchAddress("truth_recoilElectronVX", &truth_recoilElectronVX);
    tree->SetBranchAddress("truth_recoilElectronVY", &truth_recoilElectronVY);
    tree->SetBranchAddress("truth_recoilElectronVZ", &truth_recoilElectronVZ);
    tree->SetBranchAddress("truth_recoilElectronPX", &truth_recoilElectronPX);
    tree->SetBranchAddress("truth_recoilElectronPY", &truth_recoilElectronPY);
    tree->SetBranchAddress("truth_recoilElectronPZ", &truth_recoilElectronPZ);
    tree->SetBranchAddress("truth_recoilElectronEnergy", &truth_recoilElectronEnergy);
    tree->SetBranchAddress("truth_recoilElectronAngle", &truth_recoilElectronAngle);
    tree->SetBranchAddress("truth_recoilElectronETheta2", &truth_recoilElectronETheta2);
    tree->SetBranchAddress("truth_recoilElectronDX", &truth_recoilElectronDX);
    tree->SetBranchAddress("truth_recoilElectronDY", &truth_recoilElectronDY);
    tree->SetBranchAddress("truth_recoilElectronDZ", &truth_recoilElectronDZ);
    
    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness", &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_slicePurity", &reco_slicePurity);
    tree->SetBranchAddress("reco_sliceScore", &reco_sliceScore);
    tree->SetBranchAddress("reco_sliceCategory", &reco_sliceCategory);
    tree->SetBranchAddress("reco_sliceInteraction", &reco_sliceInteraction);
    tree->SetBranchAddress("reco_sliceTrueVX", &reco_sliceTrueVX);
    tree->SetBranchAddress("reco_sliceTrueVY", &reco_sliceTrueVY);
    tree->SetBranchAddress("reco_sliceTrueVZ", &reco_sliceTrueVZ);
    tree->SetBranchAddress("reco_sliceNumHits", &reco_sliceNumHits);
    tree->SetBranchAddress("reco_sliceNumHitsTruthMatched", &reco_sliceNumHitsTruthMatched);
    tree->SetBranchAddress("reco_sliceNumTruthHits", &reco_sliceNumTruthHits);
    tree->SetBranchAddress("reco_sliceOrigin", &reco_sliceOrigin);
    tree->SetBranchAddress("reco_sliceTrueCCNC", &reco_sliceTrueCCNC);
    tree->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);

    tree->SetBranchAddress("truth_particleSliceID", &truth_particleSliceID);
    tree->SetBranchAddress("truth_particlePrimary", &truth_particlePrimary);
    tree->SetBranchAddress("truth_particleVX", &truth_particleVX);
    tree->SetBranchAddress("truth_particleVY", &truth_particleVY);
    tree->SetBranchAddress("truth_particleVZ", &truth_particleVZ);
    tree->SetBranchAddress("truth_particlePDG", &truth_particlePDG);
    tree->SetBranchAddress("truth_particleTrackID", &truth_particleTrackID);
    tree->SetBranchAddress("truth_particleMother", &truth_particleMother);
    tree->SetBranchAddress("truth_particleStatusCode", &truth_particleStatusCode);

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
    tree->SetBranchAddress("reco_particleNumHitsTruthMatched", &reco_particleNumHitsTruthMatched);
    tree->SetBranchAddress("reco_particleNumTruthHits", &reco_particleNumTruthHits);
    tree->SetBranchAddress("reco_particleClearCosmic", &reco_particleClearCosmic);
    tree->SetBranchAddress("reco_particlePlane0dEdx", &reco_particlePlane0dEdx);
    tree->SetBranchAddress("reco_particlePlane1dEdx", &reco_particlePlane1dEdx);
    tree->SetBranchAddress("reco_particlePlane2dEdx", &reco_particlePlane2dEdx);
    tree->SetBranchAddress("reco_particleBestPlanedEdx", &reco_particleBestPlanedEdx);
    tree->SetBranchAddress("reco_particleRazzledPDG11", &reco_particleRazzledPDG11);
    tree->SetBranchAddress("reco_particleRazzledPDG13", &reco_particleRazzledPDG13);
    tree->SetBranchAddress("reco_particleRazzledPDG22", &reco_particleRazzledPDG22);
    tree->SetBranchAddress("reco_particleRazzledPDG211", &reco_particleRazzledPDG211);
    tree->SetBranchAddress("reco_particleRazzledPDG2212", &reco_particleRazzledPDG2212);
    tree->SetBranchAddress("reco_particleRazzledBestPDG", &reco_particleRazzledBestPDG);
    tree->SetBranchAddress("reco_particleShowerLength", &reco_particleShowerLength);
    tree->SetBranchAddress("reco_particleShowerOpenAngle", &reco_particleShowerOpenAngle);
    tree->SetBranchAddress("reco_particleShowerBestPlaneEnergy", &reco_particleShowerBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTrueVX", &reco_particleTrueVX);
    tree->SetBranchAddress("reco_particleTrueVY", &reco_particleTrueVY);
    tree->SetBranchAddress("reco_particleTrueVZ", &reco_particleTrueVZ);
    tree->SetBranchAddress("reco_particleTrueEndX", &reco_particleTrueEndX);
    tree->SetBranchAddress("reco_particleTrueEndY", &reco_particleTrueEndY);
    tree->SetBranchAddress("reco_particleTrueEndZ", &reco_particleTrueEndZ);
    
    tree->SetBranchAddress("reco_neutrinoID", &reco_neutrinoID);
    tree->SetBranchAddress("reco_neutrinoPDG", &reco_neutrinoPDG);
    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);
    
    tree->SetBranchAddress("angleRecalculationPCASlice_angle", &angleRecalculationPCASlice_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice_sliceID", &angleRecalculationPCASlice_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice_dx", &angleRecalculationPCASlice_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice_dy", &angleRecalculationPCASlice_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice_dz", &angleRecalculationPCASlice_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_angle", &angleRecalculationPCASlice5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_sliceID", &angleRecalculationPCASlice5cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dx", &angleRecalculationPCASlice5cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dy", &angleRecalculationPCASlice5cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dz", &angleRecalculationPCASlice5cm_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_angle", &angleRecalculationPCASlice10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_sliceID", &angleRecalculationPCASlice10cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dx", &angleRecalculationPCASlice10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dy", &angleRecalculationPCASlice10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dz", &angleRecalculationPCASlice10cm_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_angle", &angleRecalculationPCASlice15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_sliceID", &angleRecalculationPCASlice15cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dx", &angleRecalculationPCASlice15cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dy", &angleRecalculationPCASlice15cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dz", &angleRecalculationPCASlice15cm_dz);
    
    tree->SetBranchAddress("angleRecalculationPCAPFP_angle", &angleRecalculationPCAPFP_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP_pfpID", &angleRecalculationPCAPFP_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dx", &angleRecalculationPCAPFP_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dy", &angleRecalculationPCAPFP_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dz", &angleRecalculationPCAPFP_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_angle", &angleRecalculationPCAPFP5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_pfpID", &angleRecalculationPCAPFP5cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dx", &angleRecalculationPCAPFP5cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dy", &angleRecalculationPCAPFP5cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dz", &angleRecalculationPCAPFP5cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_angle", &angleRecalculationPCAPFP10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_pfpID", &angleRecalculationPCAPFP10cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dx", &angleRecalculationPCAPFP10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dy", &angleRecalculationPCAPFP10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dz", &angleRecalculationPCAPFP10cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_angle", &angleRecalculationPCAPFP15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_pfpID", &angleRecalculationPCAPFP15cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dx", &angleRecalculationPCAPFP15cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dy", &angleRecalculationPCAPFP15cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dz", &angleRecalculationPCAPFP15cm_dz);

    Long64_t numEntries = tree->GetEntries();    
    
    // NuEWeights Tree branch variable
    UInt_t eventID_weights, runID_weights, subRunID_weights;
    int nuEScatter_weights;
    double nuEScatterTrueVX_weights, nuEScatterTrueVY_weights, nuEScatterTrueVZ_weights;
    double DLCurrent_weights, signal_weights;

    std::vector<double> *nuEScatter_MCTruthFlux_weight_horncurrent = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_expskin = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_pioninexsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_pionqexsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_piontotxsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_nucleoninexsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_nucleonqexsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_nucleontotxsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_kplus = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_kmin = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_kzero = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_piplus = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_piminus = nullptr;

    std::vector<double> *reco_sliceID_weights = nullptr;  
    std::vector<double> *reco_sliceInteraction_weights = nullptr;  
    std::vector<double> *reco_sliceTrueVX_weights = nullptr;  
    std::vector<double> *reco_sliceTrueVY_weights = nullptr;  
    std::vector<double> *reco_sliceTrueVZ_weights = nullptr;  
    std::vector<double> *reco_sliceOrigin_weights = nullptr;  
    std::vector<double> *reco_sliceTrueCCNC_weights = nullptr;  
    std::vector<double> *reco_sliceTrueNeutrinoType_weights = nullptr;  
    
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_horncurrent = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_expskin = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_pioninexsec = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_pionqexsec = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_piontotxsec = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_nucleoninexsec = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_nucleonqexsec = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_nucleontotxsec = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_kplus = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_kmin = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_kzero = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_piplus = nullptr;  
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_piminus = nullptr;  

    weightsTree->SetBranchAddress("eventID", &eventID_weights);
    weightsTree->SetBranchAddress("runID", &runID_weights);
    weightsTree->SetBranchAddress("subRunID", &subRunID_weights);
    weightsTree->SetBranchAddress("DLCurrent", &DLCurrent_weights);
    weightsTree->SetBranchAddress("signal", &signal_weights);
    
    weightsTree->SetBranchAddress("nuEScatter", &nuEScatter_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVX", &nuEScatterTrueVX_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVY", &nuEScatterTrueVY_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVZ", &nuEScatterTrueVZ_weights);
    
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_horncurrent", &nuEScatter_MCTruthFlux_weight_horncurrent);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_expskin", &nuEScatter_MCTruthFlux_weight_expskin);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_pioninexsec", &nuEScatter_MCTruthFlux_weight_pioninexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_pionqexsec", &nuEScatter_MCTruthFlux_weight_pionqexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_piontotxsec", &nuEScatter_MCTruthFlux_weight_piontotxsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_nucleoninexsec", &nuEScatter_MCTruthFlux_weight_nucleoninexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_nucleonqexsec", &nuEScatter_MCTruthFlux_weight_nucleonqexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_nucleontotxsec", &nuEScatter_MCTruthFlux_weight_nucleontotxsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_kplus", &nuEScatter_MCTruthFlux_weight_kplus);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_kmin", &nuEScatter_MCTruthFlux_weight_kmin);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_kzero", &nuEScatter_MCTruthFlux_weight_kzero);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_piplus", &nuEScatter_MCTruthFlux_weight_piplus);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_piminus", &nuEScatter_MCTruthFlux_weight_piminus);
    
    weightsTree->SetBranchAddress("reco_sliceID", &reco_sliceID_weights);
    weightsTree->SetBranchAddress("reco_sliceInteraction", &reco_sliceInteraction_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueVX", &reco_sliceTrueVX_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueVY", &reco_sliceTrueVY_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueVZ", &reco_sliceTrueVZ_weights);
    weightsTree->SetBranchAddress("reco_sliceOrigin", &reco_sliceOrigin_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueCCNC", &reco_sliceTrueCCNC_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType_weights);
    
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

    Long64_t numEntries_weights = weightsTree->GetEntries();    

    std::unordered_map<eventKey_struct, Long64_t, eventKeyHash_struct> weightEntryMap;
    for(Long64_t i = 0; i < numEntries_weights; ++i){
        weightsTree->GetEntry(i);

        eventKey_struct key{runID_weights, subRunID_weights, eventID_weights, static_cast<int>(signal_weights), static_cast<int>(DLCurrent_weights)};

        weightEntryMap[key] = i;
    }

    const int NUNIV = 1000;

    // Per-parameter histograms: each shows how total signal count varies across 1000 universes
    TH1D* h_horncurrent  = new TH1D("h_horncurrent",  "Horn Current;Total nu+e count;Universes",     60, 0, 600);
    TH1D* h_expskin      = new TH1D("h_expskin",      "Exp Skin;Total nu+e count;Universes",         60, 0, 600);
    TH1D* h_kplus        = new TH1D("h_kplus",        "K+;Total nu+e count;Universes",               60, 0, 600);
    TH1D* h_kmin         = new TH1D("h_kmin",         "K-;Total nu+e count;Universes",               60, 0, 600);
    TH1D* h_kzero        = new TH1D("h_kzero",        "K0;Total nu+e count;Universes",               60, 0, 600);
    TH1D* h_nucleoninex  = new TH1D("h_nucleoninex",  "Nucleon InelXsec;Total nu+e count;Universes", 60, 0, 600);
    TH1D* h_nucleonqex   = new TH1D("h_nucleonqex",   "Nucleon QelXsec;Total nu+e count;Universes",  60, 0, 600);
    TH1D* h_nucleontotx  = new TH1D("h_nucleontotx",  "Nucleon TotXsec;Total nu+e count;Universes",  60, 0, 600);
    TH1D* h_piminus      = new TH1D("h_piminus",      "Pi-;Total nu+e count;Universes",              60, 0, 600);
    TH1D* h_pioninex     = new TH1D("h_pioninex",     "Pion InelXsec;Total nu+e count;Universes",    60, 0, 600);
    TH1D* h_pionqex      = new TH1D("h_pionqex",      "Pion QelXsec;Total nu+e count;Universes",     60, 0, 600);
    TH1D* h_piontotx     = new TH1D("h_piontotx",     "Pion TotXsec;Total nu+e count;Universes",     60, 0, 600);
    TH1D* h_piplus       = new TH1D("h_piplus",       "Pi+;Total nu+e count;Universes",              60, 0, 600);

    // Running totals across events, one entry per universe
    std::vector<double> count_horncurrent(NUNIV, 0.0);
    std::vector<double> count_expskin(NUNIV, 0.0);
    std::vector<double> count_kplus(NUNIV, 0.0);
    std::vector<double> count_kmin(NUNIV, 0.0);
    std::vector<double> count_kzero(NUNIV, 0.0);
    std::vector<double> count_nucleoninex(NUNIV, 0.0);
    std::vector<double> count_nucleonqex(NUNIV, 0.0);
    std::vector<double> count_nucleontotx(NUNIV, 0.0);
    std::vector<double> count_piminus(NUNIV, 0.0);
    std::vector<double> count_pioninex(NUNIV, 0.0);
    std::vector<double> count_pionqex(NUNIV, 0.0);
    std::vector<double> count_piontotx(NUNIV, 0.0);
    std::vector<double> count_piplus(NUNIV, 0.0);

    double actualSignalCount = 0.0;


    for(Long64_t e = 0; e < numEntries; ++e){
        //std::cout << "============= New Event =============" << std::endl;
        tree->GetEntry(e);

        //std::cout << "DLCurrent = " << DLCurrent << ", signal = " << signal << ", eventID = " << eventID << ", subRunID = " << subRunID << ", runID = " << runID << std::endl;
        //std::cout << "True nu+e elastic scatter in event = " << nuEScatter << ", True vertex = (" << nuEScatterTrueVX << ", " << nuEScatterTrueVY << ", " << nuEScatterTrueVZ << ")" << std::endl;

        if(reco_sliceID->size() == 0) continue;

        //std::cout << "--- Slices for event in NuE tree ---" << std::endl;
        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) == -999999) continue;
            //std::cout << "Slice " << slice << ": ID = " << reco_sliceID->at(slice) << ", Interaction = " << reco_sliceInteraction->at(slice) << ", True Vertex = (" << reco_sliceTrueVX->at(slice) << ", " << reco_sliceTrueVY->at(slice) << ", " << reco_sliceTrueVZ->at(slice) << "), Origin = " << reco_sliceOrigin->at(slice) << ", True CCNC = " << reco_sliceTrueCCNC->at(slice) << ", True Neutrino Type = " << reco_sliceTrueNeutrinoType->at(slice) << std::endl;
        }
        //std::cout << "------------------------------------" << std::endl;

        if(signal == 2 || signal == 1){
            // This is either a BNB or signal event
            //std::cout << "This is a BNB or signal event -> Look for weights" << std::endl;

            eventKey_struct key{runID, subRunID, eventID, static_cast<int>(signal), static_cast<int>(DLCurrent)};

            auto it = weightEntryMap.find(key);

            if(it == weightEntryMap.end()){
                //std::cout << "No matching weights event found" << std::endl;
            } else{
                Long64_t weightEntry = it->second;
                weightsTree->GetEntry(weightEntry);
                //std::cout << "Found matching weights event at entry " << weightEntry << std::endl;

                //std::cout << "DLCurrent = " << DLCurrent_weights << ", signal = " << signal_weights << ", eventID = " << eventID_weights << ", subRunID = " << subRunID_weights << ", runID = " << runID_weights << std::endl;
                //std::cout << "True nu+e elastic scatter in event = " << nuEScatter_weights << ", True vertex = (" << nuEScatterTrueVX_weights << ", " << nuEScatterTrueVY_weights << ", " << nuEScatterTrueVZ_weights << ")" << std::endl;
                //std::cout << "  weight of horncurrent universe 1 = " << nuEScatter_MCTruthFlux_weight_horncurrent->at(0) << std::endl;
                
                if(nuEScatter == 1 && signal == 1 && DLCurrent == 5){
                    // This is an event with a nu+e elastic scatter in it (from the signal files)
                    actualSignalCount += weights.signalNuE; // nominal value

                    bool weightsValid = (nuEScatter_MCTruthFlux_weight_horncurrent->size() == NUNIV);

                    if(weightsValid){
                        for(int u = 0; u < NUNIV; u++){
                            count_horncurrent[u] += weights.signalNuE * nuEScatter_MCTruthFlux_weight_horncurrent->at(u);
                            count_expskin[u]     += weights.signalNuE * nuEScatter_MCTruthFlux_weight_expskin->at(u);
                            count_kplus[u]       += weights.signalNuE * nuEScatter_MCTruthFlux_weight_kplus->at(u);
                            count_kmin[u]        += weights.signalNuE * nuEScatter_MCTruthFlux_weight_kmin->at(u);
                            count_kzero[u]       += weights.signalNuE * nuEScatter_MCTruthFlux_weight_kzero->at(u);
                            count_nucleoninex[u] += weights.signalNuE * nuEScatter_MCTruthFlux_weight_nucleoninexsec->at(u);
                            count_nucleonqex[u]  += weights.signalNuE * nuEScatter_MCTruthFlux_weight_nucleonqexsec->at(u);
                            count_nucleontotx[u] += weights.signalNuE * nuEScatter_MCTruthFlux_weight_nucleontotxsec->at(u);
                            count_piminus[u]     += weights.signalNuE * nuEScatter_MCTruthFlux_weight_piminus->at(u);
                            count_pioninex[u]    += weights.signalNuE * nuEScatter_MCTruthFlux_weight_pioninexsec->at(u);
                            count_pionqex[u]     += weights.signalNuE * nuEScatter_MCTruthFlux_weight_pionqexsec->at(u);
                            count_piontotx[u]    += weights.signalNuE * nuEScatter_MCTruthFlux_weight_piontotxsec->at(u);
                            count_piplus[u]      += weights.signalNuE * nuEScatter_MCTruthFlux_weight_piplus->at(u);
                        }
                    }
                }
                
                //std::cout << "Number of slices = " << reco_sliceID_weights->size() << std::endl;
                //std::cout << "--- Slices for event in NuEWeights tree ---" << std::endl;

                //std::cout << "reco_sliceID_weights->size() = " << reco_sliceID_weights->size() << ", reco_sliceMCTruthFlux_weight_horncurrent->size() = " << reco_sliceMCTruthFlux_weight_horncurrent->size() << std::endl;
                for(size_t sliceWeight = 0; sliceWeight < reco_sliceID_weights->size(); ++sliceWeight){
                    if(reco_sliceID_weights->at(sliceWeight) == -999999) continue;
                    //std::cout << "Slice " << sliceWeight << ": ID = " << reco_sliceID_weights->at(sliceWeight) << ", Interaction = " << reco_sliceInteraction_weights->at(sliceWeight) << ", True Vertex = (" << reco_sliceTrueVX_weights->at(sliceWeight) << ", " << reco_sliceTrueVY_weights->at(sliceWeight) << ", " << reco_sliceTrueVZ_weights->at(sliceWeight) << "), Origin = " << reco_sliceOrigin_weights->at(sliceWeight) << ", True CCNC = " << reco_sliceTrueCCNC_weights->at(sliceWeight) << ", True Neutrino Type = " << reco_sliceTrueNeutrinoType_weights->at(sliceWeight) << ", Number of entries in horn current vector = " << reco_sliceMCTruthFlux_weight_horncurrent->at(sliceWeight).size() << std::endl;
                    //std::cout << "  weight of horncurrent universe 1 = " << reco_sliceMCTruthFlux_weight_horncurrent->at(sliceWeight).at(0) << std::endl;
                    //std::cout << "  weight of horncurrent universe 2 = " << reco_sliceMCTruthFlux_weight_horncurrent->at(sliceWeight).at(1) << std::endl;
                    std::cout << "Slice " << sliceWeight << ": ID = " << reco_sliceID_weights->at(sliceWeight) << ", CRUMBS Score = " << reco_sliceScore->at(sliceWeight) << std::endl; 
                    std::cout << "  Universe 1 Horn Current Weight for Slice = " << reco_sliceMCTruthFlux_weight_horncurrent->at(sliceWeight).at(0) << std::endl; 
                }
                
                //std::cout << "-------------------------------------------" << std::endl;
            }


        } else{
            //std::cout << "Signal = " << signal << " -> cosmic slice, no weights" << std::endl;
        }
        
        //std::cout << "=====================================" << std::endl;
    }

    // Fill each histogram with 1000 universe totals
    for(int u = 0; u < NUNIV; u++){
        h_horncurrent->Fill(count_horncurrent[u]);
        h_expskin->Fill(count_expskin[u]);
        h_kplus->Fill(count_kplus[u]);
        h_kmin->Fill(count_kmin[u]);
        h_kzero->Fill(count_kzero[u]);
        h_nucleoninex->Fill(count_nucleoninex[u]);
        h_nucleonqex->Fill(count_nucleonqex[u]);
        h_nucleontotx->Fill(count_nucleontotx[u]);
        h_piminus->Fill(count_piminus[u]);
        h_pioninex->Fill(count_pioninex[u]);
        h_pionqex->Fill(count_pionqex[u]);
        h_piontotx->Fill(count_piontotx[u]);
        h_piplus->Fill(count_piplus[u]);
    }

    // Systematic uncertainty per parameter
    std::cout << "\n=== Systematic Uncertainties on nu+e Signal Count ===" << std::endl;
    std::cout << Form("Nominal: %.2f", actualSignalCount) << std::endl;

    std::vector<double> systValues;

    auto computeSyst = [&](const std::string& name, TH1D* h, double nominal){
        double mean   = h->GetMean();
        double stddev = h->GetStdDev();
        double shift  = mean - nominal;
        std::cout << Form("%-20s  mean=%.2f  shift=%.2f (%+.1f%%)  syst=%.2f (%.1f%%)", name.c_str(), mean, shift, 100.*shift/nominal, stddev, 100.*stddev/nominal) << std::endl;
        systValues.push_back(stddev);
    };

    computeSyst("horncurrent",  h_horncurrent,  actualSignalCount);
    computeSyst("expskin",      h_expskin,      actualSignalCount);
    computeSyst("kplus",        h_kplus,        actualSignalCount);
    computeSyst("kmin",         h_kmin,         actualSignalCount);
    computeSyst("kzero",        h_kzero,        actualSignalCount);
    computeSyst("nucleoninex",  h_nucleoninex,  actualSignalCount);
    computeSyst("nucleonqex",   h_nucleonqex,   actualSignalCount);
    computeSyst("nucleontotx",  h_nucleontotx,  actualSignalCount);
    computeSyst("piminus",      h_piminus,      actualSignalCount);
    computeSyst("pioninex",     h_pioninex,     actualSignalCount);
    computeSyst("pionqex",      h_pionqex,      actualSignalCount);
    computeSyst("piontotx",     h_piontotx,     actualSignalCount);
    computeSyst("piplus",       h_piplus,       actualSignalCount);

    // Total systematic in quadrature
    double totalSystSq = 0.0;
    for(double s : systValues) totalSystSq += s * s;
    double totalSyst = sqrt(totalSystSq);

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << Form("%-20s  syst=%.2f (%.1f%%)", "TOTAL (quadrature)", totalSyst, 100.*totalSyst/actualSignalCount) << std::endl;
    std::cout << Form("%-20s  %.2f +/- %.2f (syst)", "Signal count", actualSignalCount, totalSyst) << std::endl;

    auto plotUniverseDist = [&](const std::string& paramName, TH1D* h, double nominal){

        TCanvas *c = new TCanvas(("c_" + paramName).c_str(), "", 800, 600);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.08);

        // Style the histogram
        h->SetLineColor(kBlue+1);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitle("Total nu+e Signal Count");
        h->GetYaxis()->SetTitle("Universes");
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.1);
        h->SetStats(0);  // hide the stats box

        h->Draw("HIST E");  // E draws error bars on the histogram bins

        // Nominal vertical line
        TLine *nomLine = new TLine(nominal, 0, nominal, h->GetMaximum() * 1.05);
        nomLine->SetLineColor(kMagenta+1);
        nomLine->SetLineWidth(2);
        nomLine->Draw("SAME");

        // Nominal label
        TLatex latex;
        latex.SetTextColor(kMagenta+1);
        latex.SetTextSize(0.04);
        latex.SetTextFont(62);  // bold
        // Position label to the left of the nominal line
        double labelX = nominal - (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) * 0.25;
        double labelY = h->GetMaximum() * 0.65;
        latex.DrawLatex(labelX, labelY, Form("Nominal: %.2f", nominal));

        // POT label top right (grey, like reference)
        TLatex potLabel;
        potLabel.SetTextColor(kGray+1);
        potLabel.SetTextSize(0.035);
        potLabel.SetNDC();  // use normalised coordinates
        potLabel.DrawLatex(0.70, 0.93, "1#times10^{21} POT");

        // Parameter name label top left
        TLatex paramLabel;
        paramLabel.SetTextSize(0.04);
        paramLabel.SetNDC();
        paramLabel.DrawLatex(0.15, 0.85, paramName.c_str());

        c->Update();

        std::string outPath = "/nashome/c/coackley/systPlots/nuE_signalCount_" + paramName + ".pdf";
        c->SaveAs(outPath.c_str());

        delete nomLine;
        delete c;
    };

    // Plot all 13
    plotUniverseDist("horncurrent", h_horncurrent, actualSignalCount);
    plotUniverseDist("expskin",     h_expskin,     actualSignalCount);
    plotUniverseDist("kplus",       h_kplus,       actualSignalCount);
    plotUniverseDist("kmin",        h_kmin,        actualSignalCount);
    plotUniverseDist("kzero",       h_kzero,       actualSignalCount);
    plotUniverseDist("nucleoninex", h_nucleoninex, actualSignalCount);
    plotUniverseDist("nucleonqex",  h_nucleonqex,  actualSignalCount);
    plotUniverseDist("nucleontotx", h_nucleontotx, actualSignalCount);
    plotUniverseDist("piminus",     h_piminus,     actualSignalCount);
    plotUniverseDist("pioninex",    h_pioninex,    actualSignalCount);
    plotUniverseDist("pionqex",     h_pionqex,     actualSignalCount);
    plotUniverseDist("piontotx",    h_piontotx,    actualSignalCount);
    plotUniverseDist("piplus",      h_piplus,      actualSignalCount);

}
