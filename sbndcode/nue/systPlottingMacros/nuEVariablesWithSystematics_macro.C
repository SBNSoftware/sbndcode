// Branch naming convention
//   {varName}_{catName}_{paramName}
//   -> vector<vector<double>> of shape [1001][nBins]
//      index 0          = nominal histogram bin contents
//      indices 1..1000  = universe 0..999 histogram bin contents
//
// One TTree per cut level, named after the cut (e.g. "clearCosmic").
// A second small TTree "sigCount" holds the per-universe signal
// counts (shape [nParams][1001], index 0 = nominal).
//
// Usage:
//   root -l -b -q 'computeSystematics_variables.C(0)' -> for no cuts
//   root -l -b -q 'computeSystematics_variables.C(1)' -> for clear cosmic cut

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
    double showerLength = -999999;
    double showerOpenAngle = -999999;
    double showerBestPlaneEnergy = -999999;
    double trueVX = -999999;
    double trueVY = -999999;
    double trueVZ = -999999;
    double trueEndX = -999999;
    double trueEndY = -999999;
    double trueEndZ = -999999;
    double trueLength = -999999;
    double numHits = -999999;
    double clearCosmic = -999999;
};

struct VarConfig {
    std::string name;
    int nBins;
    double min;
    double max;
    std::string xAxisLabel;
    std::string titleLabel;
};

// Returns a list fo 25 variables to plot. Inlcudes the number of bins, xmin, xmax, x axis label and title
std::vector<VarConfig> getAllVariableConfigs(){
    return {
        {"sliceCompleteness", 102, 0, 1.02, "Completeness", "Slice Completeness"},
        {"sliceCRUMBS", 25, -1, 1, "CRUMBS Score", "Slice CRUMBS Score"},
        {"slicePurity", 100, 0, 1, "Purity", "Slice Purity"},
        {"sliceNumRecoNeut", 10, 0, 10, "Number of Reco Neutrinos", "Number of Reco Neutrinos in Slice"},
        {"sliceNumPFPs", 20, 0, 20, "Number of PFPs", "Number of PFPs in Slice"},
        {"sliceNumPrimaryPFPs", 20, 0, 20, "Number of Primary PFPs", "Number of Primary PFPs in Slice"},
        {"sliceNumPrimaryPFPs10", 20, 0, 20, "Number of Primary PFPs", "Number of Primary PFPs with Number of Hits > 10"},
        {"sliceFracHitsInPFPs", 20, 0, 1, "Fraction", "Fraction of Hits in Slice Contained in PFPs"},
        {"sliceFracHitsInHighestEnergyPFPs", 20, 0, 1, "Fraction", "Fraction of Hits in Slice Contained in Highest Energy PFP in Slice"},
        {"ERecoSumThetaReco", 27, 0, 13.797, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice"},
        {"ERecoHighestThetaReco", 27, 0, 13.797, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice"},
        {"ERecoHighestThetaReco_pfp10cmPoints", 27, 0, 13.797, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice"},
        {"dEdx", 40, 0, 10, "dE/dx (MeV/cm)", "dE/dx of the PFP in the Slice with the Highest Energy"},
        {"razzledPDG11", 80, 0, 1, "Score", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy"},
        {"razzledPDG13", 80, 0, 1, "Score", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy"},
        {"razzledPDG22", 80, 0, 1, "Score", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy"},
        {"razzledPDG211", 80, 0, 1, "Score", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy"},
        {"razzledPDG2212", 80, 0, 1, "Score", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy"},
        {"pfpCompleteness", 50, 0, 1, "Completeness", "Completeness of the PFP in the Slice with the Highest Energy"},
        {"pfpPurity", 50, 0, 1, "Purity", "Purity of the PFP in the Slice with the Highest Energy"},
        {"pfpNumHits", 60, 0, 3000, "Number of Hits", "Number of Hits in the PFP with the Highest Energy in the Slice"},
        {"sliceNumHits", 60, 0, 3000, "Number of Hits", "Number of Hits in the Slice"},
        {"recoVXSmallerBins", 808, -202, 202, "x_{Reco} (cm)", "X Coordinate of Reco Neutrino in Slice"},
        {"recoVYSmallerBins", 816, -204, 204, "y_{Reco} (cm)", "Y Coordinate of Reco Neutrino in Slice"},
        {"recoVZSmallerBins", 1020, 0, 510, "z_{Reco} (cm)", "Z Coordinate of Reco Neutrino in Slice"},
    };
}

// Struct which indicates the cuts that are active
struct CutConfig {
    std::string name;
    int clearCosmicCut;
    int numPFPs0Cut;
    int numRecoNeutrinosCut;
    int CRUMBSCut;
    int FVCut;
    int primaryPFPCut;
    int ETheta2Cut;
    int razzledPDG11Cut;
    int razzledPDG211Cut;
    int dEdxCut;
};

std::vector<CutConfig> getAllCutConfigs(){
    return {
        // name                                                                                 CC  nPFP nRN  CRUMBS FV  prim ETheta2 r11  r211 dEdx
        {"clearCosmic",                                                                         1,  0,   0,   0,     0,  0,   0,      0,   0,   0},
        {"clearCosmic_numPFPs",                                                                 1,  1,   0,   0,     0,  0,   0,      0,   0,   0},
        {"clearCosmic_numPFPs_recoNeut",                                                        1,  1,   1,   0,     0,  0,   0,      0,   0,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS",                                                 1,  1,   1,   1,     0,  0,   0,      0,   0,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS_FV",                                              1,  1,   1,   1,     1,  0,   0,      0,   0,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS_FV_primary",                                      1,  1,   1,   1,     1,  1,   0,      0,   0,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS_FV_primary_ETheta2",                              1,  1,   1,   1,     1,  1,   1,      0,   0,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS_FV_primary_ETheta2_razzled11",                    1,  1,   1,   1,     1,  1,   1,      1,   0,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS_FV_primary_ETheta2_razzled11_razzled211",         1,  1,   1,   1,     1,  1,   1,      1,   1,   0},
        {"clearCosmic_numPFPs_recoNeut_CRUMBS_FV_primary_ETheta2_razzled11_razzled211_dEdx",    1,  1,   1,   1,     1,  1,   1,      1,   1,   1},
    };
}

// Helper function that given a variable nname and slice quanities it will return the value to be put in the histogram
// i.e. "sliceCompleteness" will return reco_sliceCompleteness_val
double getVariableValue(const std::string& varName, double reco_sliceCompleteness_val, double reco_sliceScore_val, double reco_slicePurity_val, int numRecoNeutrinos, double numPFPsSlice, double numPrimaryPFPsSlice, double numPrimaryPFPs10Slice, double numHitsInPFPs, double reco_sliceNumHits_val, double summedEnergy, const highestEnergyPFP_struct& hePFP, double pfp10cm_PCAAngle, double recoVX, double recoVY, double recoVZ){
    if(varName == "sliceCompleteness") return reco_sliceCompleteness_val;
    if(varName == "sliceCRUMBS") return reco_sliceScore_val;
    if(varName == "slicePurity") return reco_slicePurity_val;
    if(varName == "sliceNumRecoNeut") return numRecoNeutrinos;
    if(varName == "sliceNumPFPs") return numPFPsSlice;
    if(varName == "sliceNumPrimaryPFPs") return numPrimaryPFPsSlice;
    if(varName == "sliceNumPrimaryPFPs10") return numPrimaryPFPs10Slice;
    if(varName == "sliceFracHitsInPFPs") return (numHitsInPFPs / reco_sliceNumHits_val);
    if(varName == "sliceFracHitsInHighestEnergyPFPs") return (hePFP.numHits / reco_sliceNumHits_val);
    if(varName == "ERecoSumThetaReco") return (summedEnergy * hePFP.theta * hePFP.theta);
    if(varName == "ERecoHighestThetaReco") return (hePFP.energy * hePFP.theta * hePFP.theta);
    if(varName == "ERecoHighestThetaReco_pfp10cmPoints") return (hePFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle);
    if(varName == "dEdx") return hePFP.bestPlanedEdx;
    if(varName == "razzledPDG11") return hePFP.razzledPDG11;
    if(varName == "razzledPDG13") return hePFP.razzledPDG13;
    if(varName == "razzledPDG22") return hePFP.razzledPDG22;
    if(varName == "razzledPDG211") return hePFP.razzledPDG211;
    if(varName == "razzledPDG2212") return hePFP.razzledPDG2212;
    if(varName == "pfpCompleteness") return hePFP.completeness;
    if(varName == "pfpPurity") return hePFP.purity;
    if(varName == "pfpNumHits") return hePFP.numHits;
    if(varName == "sliceNumHits") return reco_sliceNumHits_val;
    if(varName == "recoVXSmallerBins") return recoVX;
    if(varName == "recoVYSmallerBins") return recoVY;
    if(varName == "recoVZSmallerBins") return recoVZ;
    return -999999;
}

// Helper function to find the bin index the variable belongs to
// if > hi then goes into overflow bin, if < low then goes into bin 0
inline int findBin(double val, int nBins, double lo, double hi){
    if(val < lo) return 0;
    if(val >= hi) return nBins - 1;
    return (int)((val - lo) / (hi - lo) * nBins);
}

// CutIndex specifies which cuts to apply
void computeSystematics_variables(int cutIndex = 0){

    // Given the cutIndex, returns the cut names
    std::vector<CutConfig> allCutConfigs = getAllCutConfigs();
    if(cutIndex < 0 || cutIndex >= (int)allCutConfigs.size()){
        std::cerr << "Invalid cutIndex " << cutIndex << ". Must be 0-" << (int)allCutConfigs.size()-1 << std::endl;
        return;
    }

    CutConfig thisCut = allCutConfigs[cutIndex];
    std::cout << "Processing cut " << cutIndex << ": " << thisCut.name << std::endl;

    // Name of the output root file
    std::string outputFile = "/exp/sbnd/data/users/coackley/systHistograms_" + thisCut.name + ".root";

    // Cut values
    double crumbsScoreCut_low  = 0.2;
    double crumbsScoreCut_high = 0.68;

    double FVCut_xHigh   = 194;
    double FVCut_xLow    = -196;
    double FVCut_xCentre = 10;
    double FVCut_yHigh   = 196;
    double FVCut_yLow    = -196;
    double FVCut_zHigh   = 450;
    double FVCut_zLow    = 6.5;

    double primaryPFPCutValue = 1;

    double razzled211High_highestEnergyPFP = 0.05;
    double razzled211Low_highestEnergyPFP  = 0.00;
    double razzled11High_highestEnergyPFP  = 1;
    double razzled11Low_highestEnergyPFP   = 0.75;

    double dEdxHigh_highestEnergyPFP = 3;
    double dEdxLow_highestEnergyPFP  = 0.25;

    double ETheta2High_highestEnergyPFP = 1.533;
    double ETheta2Low_highestEnergyPFP  = 0;

    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0;      double zMax = 509.4;

    // Loads all the variables to be plotted 
    std::vector<VarConfig>  varConfigs = getAllVariableConfigs();
    int nVars = varConfigs.size();

    // Slice categories
    std::vector<std::string> catNames = {"cosmic", "signal", "signal_fuzzy", "BNB", "BNB_fuzzy"};
    
    // Flux parameters
    std::vector<std::string> paramNames = {"horncurrent", "expskin", "kplus", "kmin", "kzero", "nucleoninex", "nucleonqex", "nucleontotx", "piminus", "pioninex", "pionqex", "piontotx", "piplus", "combined_allParams"};
    
    int nCats   = catNames.size();
    int nParams = paramNames.size();
    const int NUNIV = 1000;

    // Input root files
    TFile *fNuE = TFile::Open("/exp/sbnd/data/users/coackley/signalBNBIntimeCosmic14June_withoutWeights.root");
    if(!fNuE){ std::cerr << "Error opening NuE file" << std::endl; return; }
    TDirectory *dirNuE = (TDirectory*)fNuE->Get("ana");
    if(!dirNuE){ std::cerr << "Directory 'ana' not found in NuE file" << std::endl; return; }
    TTree *tree = (TTree*)dirNuE->Get("NuE");
    if(!tree){ std::cerr << "NuE TTree not found" << std::endl; return; }
    TTree *subRunTree = (TTree*)dirNuE->Get("SubRun");
    if(!subRunTree){ std::cerr << "SubRun TTree not found" << std::endl; return; }

    TFile *fNuEWeights = TFile::Open("/exp/sbnd/data/users/coackley/signalBNBIntimeCosmic14June_withWeights.root");
    if(!fNuEWeights){ std::cerr << "Error opening NuEWeights file" << std::endl; return; }
    TDirectory *dirNuEWeights = (TDirectory*)fNuEWeights->Get("ana");
    if(!dirNuEWeights){ std::cerr << "Directory 'ana' not found in NuEWeights file" << std::endl; return; }
    TTree *weightsTree = (TTree*)dirNuEWeights->Get("NuEWeights");
    if(!weightsTree){ std::cerr << "NuEWeights TTree not found" << std::endl; return; }

    double subRunSignal, subRunDLCurrent, subRunPOT;
    int subRunSpills, subRunNumGenEvents;
    unsigned int subRunNumber, subRunRun;
    subRunTree->SetBranchAddress("signal",       &subRunSignal);
    subRunTree->SetBranchAddress("DLCurrent",    &subRunDLCurrent);
    subRunTree->SetBranchAddress("pot",          &subRunPOT);
    subRunTree->SetBranchAddress("spills",       &subRunSpills);
    subRunTree->SetBranchAddress("numGenEvents", &subRunNumGenEvents);
    subRunTree->SetBranchAddress("subRun",       &subRunNumber);
    subRunTree->SetBranchAddress("run",          &subRunRun);

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBNuE;
    double totalPOTSignalNuE = 0, totalPOTBNBNuE = 0;
    double cosmicSpillsSumNuE = 0, BNBSpillsSumNuE = 0, NuESpillsSumNuE = 0;
    double POTSignalNuE_notMissing = 0, POTBNBNuE_notMissing = 0;

    Long64_t numEntriesSubRun = subRunTree->GetEntries();
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

    double targetPOT    = 1e21;
    double targetGates  = ((1333568 / 6.293443e+18) * targetPOT);
    double cosmicsWeights_NuE = (((1 - 0.0754) * targetGates) / cosmicSpillsSumNuE);

    weights_struct weights;
    weights.signalNuE  = targetPOT / POTSignalNuE_notMissing;
    weights.BNBNuE     = targetPOT / POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;
    std::cout << "Weights DLNu+E: BNB = " << weights.BNBNuE << ", Signal = " << weights.signalNuE << ", Intime Cosmics = " << weights.cosmicsNuE << std::endl;

    UInt_t eventID, runID, subRunID;
    int nuEScatter;
    double nuEScatterTrueVX, nuEScatterTrueVY, nuEScatterTrueVZ;
    double DLCurrent, signal;

    std::vector<double> *truth_recoilElectronPDG    = nullptr;
    std::vector<double> *truth_recoilElectronEnergy = nullptr;
    std::vector<double> *truth_recoilElectronAngle  = nullptr;
    std::vector<double> *truth_recoilElectronDX     = nullptr;
    std::vector<double> *truth_recoilElectronDY     = nullptr;
    std::vector<double> *truth_recoilElectronDZ     = nullptr;
    std::vector<double> *truth_recoilElectronVX     = nullptr;
    std::vector<double> *truth_recoilElectronVY     = nullptr;
    std::vector<double> *truth_recoilElectronVZ     = nullptr;
    std::vector<double> *truth_recoilElectronPX     = nullptr;
    std::vector<double> *truth_recoilElectronPY     = nullptr;
    std::vector<double> *truth_recoilElectronPZ     = nullptr;
    std::vector<double> *truth_recoilElectronETheta2 = nullptr;

    std::vector<double> *reco_sliceID              = nullptr;
    std::vector<double> *reco_sliceCompleteness    = nullptr;
    std::vector<double> *reco_slicePurity          = nullptr;
    std::vector<double> *reco_sliceScore           = nullptr;
    std::vector<double> *reco_sliceCategory        = nullptr;
    std::vector<double> *reco_sliceInteraction     = nullptr;
    std::vector<double> *reco_sliceTrueVX          = nullptr;
    std::vector<double> *reco_sliceTrueVY          = nullptr;
    std::vector<double> *reco_sliceTrueVZ          = nullptr;
    std::vector<double> *reco_sliceNumHits         = nullptr;
    std::vector<double> *reco_sliceNumHitsTruthMatched = nullptr;
    std::vector<double> *reco_sliceNumTruthHits    = nullptr;
    std::vector<double> *reco_sliceOrigin          = nullptr;
    std::vector<double> *reco_sliceTrueCCNC        = nullptr;
    std::vector<double> *reco_sliceTrueNeutrinoType = nullptr;

    std::vector<double> *truth_particleSliceID    = nullptr;
    std::vector<double> *truth_particlePrimary    = nullptr;
    std::vector<double> *truth_particleVX         = nullptr;
    std::vector<double> *truth_particleVY         = nullptr;
    std::vector<double> *truth_particleVZ         = nullptr;
    std::vector<double> *truth_particlePDG        = nullptr;
    std::vector<double> *truth_particleTrackID    = nullptr;
    std::vector<double> *truth_particleMother     = nullptr;
    std::vector<double> *truth_particleStatusCode = nullptr;

    std::vector<double> *reco_particlePDG                  = nullptr;
    std::vector<double> *reco_particleIsPrimary             = nullptr;
    std::vector<double> *reco_particleVX                   = nullptr;
    std::vector<double> *reco_particleVY                   = nullptr;
    std::vector<double> *reco_particleVZ                   = nullptr;
    std::vector<double> *reco_particleDX                   = nullptr;
    std::vector<double> *reco_particleDY                   = nullptr;
    std::vector<double> *reco_particleDZ                   = nullptr;
    std::vector<double> *reco_particleSliceID              = nullptr;
    std::vector<double> *reco_particleBestPlaneEnergy      = nullptr;
    std::vector<double> *reco_particleTheta                = nullptr;
    std::vector<double> *reco_particleTrackScore           = nullptr;
    std::vector<double> *reco_particleCompleteness         = nullptr;
    std::vector<double> *reco_particlePurity               = nullptr;
    std::vector<double> *reco_particleID                   = nullptr;
    std::vector<double> *reco_particleTruePDG              = nullptr;
    std::vector<double> *reco_particleTrueOrigin           = nullptr;
    std::vector<double> *reco_particleTrueInteractionType  = nullptr;
    std::vector<double> *reco_particleNumHits              = nullptr;
    std::vector<double> *reco_particleNumHitsTruthMatched  = nullptr;
    std::vector<double> *reco_particleNumTruthHits         = nullptr;
    std::vector<double> *reco_particleClearCosmic          = nullptr;
    std::vector<double> *reco_particlePlane0dEdx           = nullptr;
    std::vector<double> *reco_particlePlane1dEdx           = nullptr;
    std::vector<double> *reco_particlePlane2dEdx           = nullptr;
    std::vector<double> *reco_particleBestPlanedEdx        = nullptr;
    std::vector<double> *reco_particleRazzledPDG11         = nullptr;
    std::vector<double> *reco_particleRazzledPDG13         = nullptr;
    std::vector<double> *reco_particleRazzledPDG22         = nullptr;
    std::vector<double> *reco_particleRazzledPDG211        = nullptr;
    std::vector<double> *reco_particleRazzledPDG2212       = nullptr;
    std::vector<double> *reco_particleRazzledBestPDG       = nullptr;
    std::vector<double> *reco_particleShowerLength         = nullptr;
    std::vector<double> *reco_particleShowerOpenAngle      = nullptr;
    std::vector<double> *reco_particleShowerBestPlaneEnergy = nullptr;
    std::vector<double> *reco_particleTrueVX               = nullptr;
    std::vector<double> *reco_particleTrueVY               = nullptr;
    std::vector<double> *reco_particleTrueVZ               = nullptr;
    std::vector<double> *reco_particleTrueEndX             = nullptr;
    std::vector<double> *reco_particleTrueEndY             = nullptr;
    std::vector<double> *reco_particleTrueEndZ             = nullptr;

    std::vector<double> *reco_neutrinoID      = nullptr;
    std::vector<double> *reco_neutrinoPDG     = nullptr;
    std::vector<double> *reco_neutrinoVX      = nullptr;
    std::vector<double> *reco_neutrinoVY      = nullptr;
    std::vector<double> *reco_neutrinoVZ      = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;

    std::vector<double> *angleRecalculationPCASlice_angle    = nullptr;
    std::vector<double> *angleRecalculationPCASlice_sliceID  = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dx       = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dy       = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dz       = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_angle   = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dx      = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dy      = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dz      = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_angle   = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dx      = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dy      = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dz      = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_angle   = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dx      = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dy      = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dz      = nullptr;

    std::vector<double> *angleRecalculationPCAPFP_angle     = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_pfpID     = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dx        = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dy        = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dz        = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_angle  = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_pfpID  = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dx     = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dy     = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dz     = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dx    = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dy    = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dz    = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dx    = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dy    = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dz    = nullptr;

    tree->SetBranchAddress("eventID",           &eventID);
    tree->SetBranchAddress("runID",             &runID);
    tree->SetBranchAddress("subRunID",          &subRunID);
    tree->SetBranchAddress("nuEScatter",        &nuEScatter);
    tree->SetBranchAddress("nuEScatterTrueVX",  &nuEScatterTrueVX);
    tree->SetBranchAddress("nuEScatterTrueVY",  &nuEScatterTrueVY);
    tree->SetBranchAddress("nuEScatterTrueVZ",  &nuEScatterTrueVZ);
    tree->SetBranchAddress("DLCurrent",         &DLCurrent);
    tree->SetBranchAddress("signal",            &signal);
    tree->SetBranchAddress("truth_recoilElectronPDG",     &truth_recoilElectronPDG);
    tree->SetBranchAddress("truth_recoilElectronVX",      &truth_recoilElectronVX);
    tree->SetBranchAddress("truth_recoilElectronVY",      &truth_recoilElectronVY);
    tree->SetBranchAddress("truth_recoilElectronVZ",      &truth_recoilElectronVZ);
    tree->SetBranchAddress("truth_recoilElectronPX",      &truth_recoilElectronPX);
    tree->SetBranchAddress("truth_recoilElectronPY",      &truth_recoilElectronPY);
    tree->SetBranchAddress("truth_recoilElectronPZ",      &truth_recoilElectronPZ);
    tree->SetBranchAddress("truth_recoilElectronEnergy",  &truth_recoilElectronEnergy);
    tree->SetBranchAddress("truth_recoilElectronAngle",   &truth_recoilElectronAngle);
    tree->SetBranchAddress("truth_recoilElectronETheta2", &truth_recoilElectronETheta2);
    tree->SetBranchAddress("truth_recoilElectronDX",      &truth_recoilElectronDX);
    tree->SetBranchAddress("truth_recoilElectronDY",      &truth_recoilElectronDY);
    tree->SetBranchAddress("truth_recoilElectronDZ",      &truth_recoilElectronDZ);
    tree->SetBranchAddress("reco_sliceID",               &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness",     &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_slicePurity",           &reco_slicePurity);
    tree->SetBranchAddress("reco_sliceScore",            &reco_sliceScore);
    tree->SetBranchAddress("reco_sliceCategory",         &reco_sliceCategory);
    tree->SetBranchAddress("reco_sliceInteraction",      &reco_sliceInteraction);
    tree->SetBranchAddress("reco_sliceTrueVX",           &reco_sliceTrueVX);
    tree->SetBranchAddress("reco_sliceTrueVY",           &reco_sliceTrueVY);
    tree->SetBranchAddress("reco_sliceTrueVZ",           &reco_sliceTrueVZ);
    tree->SetBranchAddress("reco_sliceNumHits",          &reco_sliceNumHits);
    tree->SetBranchAddress("reco_sliceNumHitsTruthMatched", &reco_sliceNumHitsTruthMatched);
    tree->SetBranchAddress("reco_sliceNumTruthHits",     &reco_sliceNumTruthHits);
    tree->SetBranchAddress("reco_sliceOrigin",           &reco_sliceOrigin);
    tree->SetBranchAddress("reco_sliceTrueCCNC",         &reco_sliceTrueCCNC);
    tree->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);
    tree->SetBranchAddress("truth_particleSliceID",    &truth_particleSliceID);
    tree->SetBranchAddress("truth_particlePrimary",    &truth_particlePrimary);
    tree->SetBranchAddress("truth_particleVX",         &truth_particleVX);
    tree->SetBranchAddress("truth_particleVY",         &truth_particleVY);
    tree->SetBranchAddress("truth_particleVZ",         &truth_particleVZ);
    tree->SetBranchAddress("truth_particlePDG",        &truth_particlePDG);
    tree->SetBranchAddress("truth_particleTrackID",    &truth_particleTrackID);
    tree->SetBranchAddress("truth_particleMother",     &truth_particleMother);
    tree->SetBranchAddress("truth_particleStatusCode", &truth_particleStatusCode);
    tree->SetBranchAddress("reco_particlePDG",                 &reco_particlePDG);
    tree->SetBranchAddress("reco_particleIsPrimary",           &reco_particleIsPrimary);
    tree->SetBranchAddress("reco_particleVX",                  &reco_particleVX);
    tree->SetBranchAddress("reco_particleVY",                  &reco_particleVY);
    tree->SetBranchAddress("reco_particleVZ",                  &reco_particleVZ);
    tree->SetBranchAddress("reco_particleDX",                  &reco_particleDX);
    tree->SetBranchAddress("reco_particleDY",                  &reco_particleDY);
    tree->SetBranchAddress("reco_particleDZ",                  &reco_particleDZ);
    tree->SetBranchAddress("reco_particleSliceID",             &reco_particleSliceID);
    tree->SetBranchAddress("reco_particleBestPlaneEnergy",     &reco_particleBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTheta",               &reco_particleTheta);
    tree->SetBranchAddress("reco_particleTrackScore",          &reco_particleTrackScore);
    tree->SetBranchAddress("reco_particleCompleteness",        &reco_particleCompleteness);
    tree->SetBranchAddress("reco_particlePurity",              &reco_particlePurity);
    tree->SetBranchAddress("reco_particleID",                  &reco_particleID);
    tree->SetBranchAddress("reco_particleTruePDG",             &reco_particleTruePDG);
    tree->SetBranchAddress("reco_particleTrueOrigin",          &reco_particleTrueOrigin);
    tree->SetBranchAddress("reco_particleTrueInteractionType", &reco_particleTrueInteractionType);
    tree->SetBranchAddress("reco_particleNumHits",             &reco_particleNumHits);
    tree->SetBranchAddress("reco_particleNumHitsTruthMatched", &reco_particleNumHitsTruthMatched);
    tree->SetBranchAddress("reco_particleNumTruthHits",        &reco_particleNumTruthHits);
    tree->SetBranchAddress("reco_particleClearCosmic",         &reco_particleClearCosmic);
    tree->SetBranchAddress("reco_particlePlane0dEdx",          &reco_particlePlane0dEdx);
    tree->SetBranchAddress("reco_particlePlane1dEdx",          &reco_particlePlane1dEdx);
    tree->SetBranchAddress("reco_particlePlane2dEdx",          &reco_particlePlane2dEdx);
    tree->SetBranchAddress("reco_particleBestPlanedEdx",       &reco_particleBestPlanedEdx);
    tree->SetBranchAddress("reco_particleRazzledPDG11",        &reco_particleRazzledPDG11);
    tree->SetBranchAddress("reco_particleRazzledPDG13",        &reco_particleRazzledPDG13);
    tree->SetBranchAddress("reco_particleRazzledPDG22",        &reco_particleRazzledPDG22);
    tree->SetBranchAddress("reco_particleRazzledPDG211",       &reco_particleRazzledPDG211);
    tree->SetBranchAddress("reco_particleRazzledPDG2212",      &reco_particleRazzledPDG2212);
    tree->SetBranchAddress("reco_particleRazzledBestPDG",      &reco_particleRazzledBestPDG);
    tree->SetBranchAddress("reco_particleShowerLength",        &reco_particleShowerLength);
    tree->SetBranchAddress("reco_particleShowerOpenAngle",     &reco_particleShowerOpenAngle);
    tree->SetBranchAddress("reco_particleShowerBestPlaneEnergy", &reco_particleShowerBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTrueVX",              &reco_particleTrueVX);
    tree->SetBranchAddress("reco_particleTrueVY",              &reco_particleTrueVY);
    tree->SetBranchAddress("reco_particleTrueVZ",              &reco_particleTrueVZ);
    tree->SetBranchAddress("reco_particleTrueEndX",            &reco_particleTrueEndX);
    tree->SetBranchAddress("reco_particleTrueEndY",            &reco_particleTrueEndY);
    tree->SetBranchAddress("reco_particleTrueEndZ",            &reco_particleTrueEndZ);
    tree->SetBranchAddress("reco_neutrinoID",      &reco_neutrinoID);
    tree->SetBranchAddress("reco_neutrinoPDG",     &reco_neutrinoPDG);
    tree->SetBranchAddress("reco_neutrinoVX",      &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY",      &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ",      &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice_angle",   &angleRecalculationPCASlice_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice_sliceID", &angleRecalculationPCASlice_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice_dx",      &angleRecalculationPCASlice_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice_dy",      &angleRecalculationPCASlice_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice_dz",      &angleRecalculationPCASlice_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_angle",   &angleRecalculationPCASlice5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_sliceID", &angleRecalculationPCASlice5cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dx",      &angleRecalculationPCASlice5cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dy",      &angleRecalculationPCASlice5cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dz",      &angleRecalculationPCASlice5cm_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_angle",   &angleRecalculationPCASlice10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_sliceID", &angleRecalculationPCASlice10cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dx",      &angleRecalculationPCASlice10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dy",      &angleRecalculationPCASlice10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dz",      &angleRecalculationPCASlice10cm_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_angle",   &angleRecalculationPCASlice15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_sliceID", &angleRecalculationPCASlice15cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dx",      &angleRecalculationPCASlice15cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dy",      &angleRecalculationPCASlice15cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dz",      &angleRecalculationPCASlice15cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP_angle",   &angleRecalculationPCAPFP_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP_pfpID",   &angleRecalculationPCAPFP_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dx",      &angleRecalculationPCAPFP_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dy",      &angleRecalculationPCAPFP_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dz",      &angleRecalculationPCAPFP_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_angle",  &angleRecalculationPCAPFP5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_pfpID",  &angleRecalculationPCAPFP5cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dx",     &angleRecalculationPCAPFP5cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dy",     &angleRecalculationPCAPFP5cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dz",     &angleRecalculationPCAPFP5cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_angle", &angleRecalculationPCAPFP10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_pfpID", &angleRecalculationPCAPFP10cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dx",    &angleRecalculationPCAPFP10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dy",    &angleRecalculationPCAPFP10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dz",    &angleRecalculationPCAPFP10cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_angle", &angleRecalculationPCAPFP15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_pfpID", &angleRecalculationPCAPFP15cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dx",    &angleRecalculationPCAPFP15cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dy",    &angleRecalculationPCAPFP15cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dz",    &angleRecalculationPCAPFP15cm_dz);

    Long64_t numEntries = tree->GetEntries();

    UInt_t eventID_weights, runID_weights, subRunID_weights;
    int nuEScatter_weights;
    double nuEScatterTrueVX_weights, nuEScatterTrueVY_weights, nuEScatterTrueVZ_weights;
    double DLCurrent_weights, signal_weights;

    std::vector<double> *nuEScatter_MCTruthFlux_weight_horncurrent  = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_expskin      = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_pioninexsec  = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_pionqexsec   = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_piontotxsec  = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_nucleoninexsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_nucleonqexsec  = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_nucleontotxsec = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_kplus        = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_kmin         = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_kzero        = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_piplus       = nullptr;
    std::vector<double> *nuEScatter_MCTruthFlux_weight_piminus      = nullptr;

    std::vector<double> *reco_sliceID_weights              = nullptr;
    std::vector<double> *reco_sliceInteraction_weights     = nullptr;
    std::vector<double> *reco_sliceTrueVX_weights          = nullptr;
    std::vector<double> *reco_sliceTrueVY_weights          = nullptr;
    std::vector<double> *reco_sliceTrueVZ_weights          = nullptr;
    std::vector<double> *reco_sliceOrigin_weights          = nullptr;
    std::vector<double> *reco_sliceTrueCCNC_weights        = nullptr;
    std::vector<double> *reco_sliceTrueNeutrinoType_weights = nullptr;

    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_horncurrent   = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_expskin       = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_pioninexsec   = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_pionqexsec    = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_piontotxsec   = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_nucleoninexsec = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_nucleonqexsec  = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_nucleontotxsec = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_kplus         = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_kmin          = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_kzero         = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_piplus        = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthFlux_weight_piminus       = nullptr;

    weightsTree->SetBranchAddress("eventID",   &eventID_weights);
    weightsTree->SetBranchAddress("runID",     &runID_weights);
    weightsTree->SetBranchAddress("subRunID",  &subRunID_weights);
    weightsTree->SetBranchAddress("DLCurrent", &DLCurrent_weights);
    weightsTree->SetBranchAddress("signal",    &signal_weights);
    weightsTree->SetBranchAddress("nuEScatter",        &nuEScatter_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVX",  &nuEScatterTrueVX_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVY",  &nuEScatterTrueVY_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVZ",  &nuEScatterTrueVZ_weights);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_horncurrent",   &nuEScatter_MCTruthFlux_weight_horncurrent);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_expskin",       &nuEScatter_MCTruthFlux_weight_expskin);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_pioninexsec",   &nuEScatter_MCTruthFlux_weight_pioninexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_pionqexsec",    &nuEScatter_MCTruthFlux_weight_pionqexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_piontotxsec",   &nuEScatter_MCTruthFlux_weight_piontotxsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_nucleoninexsec",&nuEScatter_MCTruthFlux_weight_nucleoninexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_nucleonqexsec", &nuEScatter_MCTruthFlux_weight_nucleonqexsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_nucleontotxsec",&nuEScatter_MCTruthFlux_weight_nucleontotxsec);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_kplus",         &nuEScatter_MCTruthFlux_weight_kplus);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_kmin",          &nuEScatter_MCTruthFlux_weight_kmin);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_kzero",         &nuEScatter_MCTruthFlux_weight_kzero);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_piplus",        &nuEScatter_MCTruthFlux_weight_piplus);
    weightsTree->SetBranchAddress("nuEScatter_MCTruthFlux_weight_piminus",       &nuEScatter_MCTruthFlux_weight_piminus);
    weightsTree->SetBranchAddress("reco_sliceID",               &reco_sliceID_weights);
    weightsTree->SetBranchAddress("reco_sliceInteraction",      &reco_sliceInteraction_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueVX",           &reco_sliceTrueVX_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueVY",           &reco_sliceTrueVY_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueVZ",           &reco_sliceTrueVZ_weights);
    weightsTree->SetBranchAddress("reco_sliceOrigin",           &reco_sliceOrigin_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueCCNC",         &reco_sliceTrueCCNC_weights);
    weightsTree->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType_weights);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_horncurrent",    &reco_sliceMCTruthFlux_weight_horncurrent);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_expskin",        &reco_sliceMCTruthFlux_weight_expskin);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_pioninexsec",    &reco_sliceMCTruthFlux_weight_pioninexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_pionqexsec",     &reco_sliceMCTruthFlux_weight_pionqexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_piontotxsec",    &reco_sliceMCTruthFlux_weight_piontotxsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_nucleoninexsec", &reco_sliceMCTruthFlux_weight_nucleoninexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_nucleonqexsec",  &reco_sliceMCTruthFlux_weight_nucleonqexsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_nucleontotxsec", &reco_sliceMCTruthFlux_weight_nucleontotxsec);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_kplus",          &reco_sliceMCTruthFlux_weight_kplus);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_kmin",           &reco_sliceMCTruthFlux_weight_kmin);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_kzero",          &reco_sliceMCTruthFlux_weight_kzero);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_piplus",         &reco_sliceMCTruthFlux_weight_piplus);
    weightsTree->SetBranchAddress("reco_sliceMCTruthFlux_weight_piminus",        &reco_sliceMCTruthFlux_weight_piminus);

    Long64_t numEntries_weights = weightsTree->GetEntries();

    // Builds a map which maps every events key (eventID, subRunID, runID, signal and DLCurrent) to its entry number in the tree
    std::unordered_map<eventKey_struct, Long64_t, eventKeyHash_struct> weightEntryMap;
    for(Long64_t i = 0; i < numEntries_weights; ++i){
        weightsTree->GetEntry(i);
        eventKey_struct key{runID_weights, subRunID_weights, eventID_weights, static_cast<int>(signal_weights), static_cast<int>(DLCurrent_weights)};
        weightEntryMap[key] = i;
    }

    // Helper function to get universe weights for true nu+e elastic scattering events
    auto getNuEWeight = [&](std::vector<double>* vec, int u) -> double {
        if(!vec || (int)vec->size() != NUNIV) return 1.0;
        return vec->at(u);
    };

    // Helper function to get universe weights for slices
    auto getSliceWeight = [&](std::vector<std::vector<double>>* vec, size_t sliceIdx, int u, bool wFound) -> double {
        if(!wFound || !vec || sliceIdx >= vec->size() || (int)vec->at(sliceIdx).size() != NUNIV) return 1.0;
        return vec->at(sliceIdx).at(u);
    };

    // Allocates a 5D nested vector: bins[variable][category][parameter][universe_index][bin]
    // Universe index 0 = nominal, 1-1000 = universe 1-1000
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> bins(nVars, std::vector<std::vector<std::vector<std::vector<double>>>>(nCats, std::vector<std::vector<std::vector<double>>>(nParams, std::vector<std::vector<double>>(NUNIV + 1, std::vector<double>()))));

    // Loop over variables
    for(int v = 0; v < nVars; v++){
        // Get the number of bins for this variable
        int nb = varConfigs[v].nBins;
        // Loop through the slice categories
        for(int cat = 0; cat < nCats; cat++){
            // Loop though the flux parameters
            for(int p = 0; p < nParams; p++){
                // Loop through the universes
                for(int ui = 0; ui <= NUNIV; ui++){
                    // Allocate bin contents. This resizes the innermost vector to length = number of bins and sets every element to 0
                    bins[v][cat][p][ui].assign(nb, 0.0);
                }
            }
        }
    }

    // A 2D array with [parameter][universe], everything initialised to 0, will be used to collect the weighted signal count for each syst universe
    std::vector<std::vector<double>> sigCountUniv(nParams, std::vector<double>(NUNIV, 0.0));
    
    // Nominal number of true nu+e elastic scattering events
    double actualSignalCount = 0.0;

    std::cout << "Starting event loop over " << numEntries << " entries..." << std::endl;

    // Loop throguh the events
    for(Long64_t e = 0; e < numEntries; ++e){
        // Progress bar every 100000 events so I can check it hasn't frozen
        if(e % 100000 == 0) std::cout << "  Entry " << e << " / " << numEntries << std::endl;

        tree->GetEntry(e);

        // Find a match between the NuE and NuEWeights Tree
        bool weightsFound = false;
        if(signal == 1 || signal == 2){
            eventKey_struct key{runID, subRunID, eventID, static_cast<int>(signal), static_cast<int>(DLCurrent)};
            auto it = weightEntryMap.find(key);
            if(it != weightEntryMap.end()){
                weightsTree->GetEntry(it->second);
                weightsFound = true;
            }
        }

        // Count the number of true nu+e elastic scattering events
        if(nuEScatter == 1 && signal == 1 && DLCurrent == 5){
            bool nuEWeightsValid = weightsFound && nuEScatter_MCTruthFlux_weight_horncurrent && (int)nuEScatter_MCTruthFlux_weight_horncurrent->size() == NUNIV;

            int FVCutFlag = thisCut.FVCut;
            bool inFV = false;
            // If FV cut not being applied then check vertex is within AV
            if(FVCutFlag == 0) inFV = (nuEScatterTrueVX > xMin && nuEScatterTrueVX < xMax && nuEScatterTrueVY > yMin && nuEScatterTrueVY < yMax && nuEScatterTrueVZ > zMin && nuEScatterTrueVZ < zMax);

            // If FV is applied then check vertex is within FV
            else inFV = (nuEScatterTrueVX > FVCut_xLow && nuEScatterTrueVX < FVCut_xHigh && std::abs(nuEScatterTrueVX) > FVCut_xCentre && nuEScatterTrueVY > FVCut_yLow && nuEScatterTrueVY < FVCut_yHigh && nuEScatterTrueVZ > FVCut_zLow && nuEScatterTrueVZ < FVCut_zHigh);
            
            // If its in FV/AV then add to the count
            if(inFV){
                actualSignalCount += weights.signalNuE;

                if(nuEWeightsValid){
                    for(int u = 0; u < NUNIV; u++){
                        double w_horn = getNuEWeight(nuEScatter_MCTruthFlux_weight_horncurrent, u);
                        double w_exp = getNuEWeight(nuEScatter_MCTruthFlux_weight_expskin, u);
                        double w_kp = getNuEWeight(nuEScatter_MCTruthFlux_weight_kplus, u);
                        double w_km = getNuEWeight(nuEScatter_MCTruthFlux_weight_kmin, u);
                        double w_k0 = getNuEWeight(nuEScatter_MCTruthFlux_weight_kzero, u);
                        double w_nie = getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleoninexsec, u);
                        double w_nqe = getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleonqexsec, u);
                        double w_ntx = getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleontotxsec, u);
                        double w_pim = getNuEWeight(nuEScatter_MCTruthFlux_weight_piminus, u);
                        double w_pie = getNuEWeight(nuEScatter_MCTruthFlux_weight_pioninexsec, u);
                        double w_piq = getNuEWeight(nuEScatter_MCTruthFlux_weight_pionqexsec, u);
                        double w_pit = getNuEWeight(nuEScatter_MCTruthFlux_weight_piontotxsec, u);
                        double w_pip = getNuEWeight(nuEScatter_MCTruthFlux_weight_piplus, u);
                        double w_comb = w_horn*w_exp*w_kp*w_km*w_k0*w_nie*w_nqe*w_ntx*w_pim*w_pie*w_piq*w_pit*w_pip;

                        std::vector<double> ws = {w_horn, w_exp, w_kp, w_km, w_k0, w_nie, w_nqe, w_ntx, w_pim, w_pie, w_piq, w_pit, w_pip, w_comb};
                        for(int p = 0; p < nParams; p++) sigCountUniv[p][u] += weights.signalNuE * ws[p];
                    }
                }
            }
        }

        // Determine POT weight for this event
        double potWeight = 0;
        if(signal == 1 && DLCurrent == 5) potWeight = weights.signalNuE;
        else if(signal == 2 && DLCurrent == 5) potWeight = weights.BNBNuE;
        else if(signal == 3 && DLCurrent == 5) potWeight = weights.cosmicsNuE;
        else continue;

        recoilElectron_struct recoilElectron;
        recoilElectron.energy = recoilElectron.angle = recoilElectron.dx = recoilElectron.dy = recoilElectron.dz = -999999;
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

        // Loop through the slices in the event
        int numSlices = (int)reco_sliceID->size();
        for(int slice = 0; slice < numSlices; slice++){
            // Count PFPs in this slice
            double numPFPsSlice = 0, numPrimaryPFPsSlice = 0, numPrimaryPFPs10Slice = 0;
            double numHitsInPFPs = 0, summedEnergy = 0;

            // Loop through particles in slice
            for(size_t pa = 0; pa < reco_particleSliceID->size(); pa++){

                if((int)reco_particleSliceID->at(pa) != (int)reco_sliceID->at(slice)) continue;
                // Skip the particle if it is a clear cosmic
                if(thisCut.clearCosmicCut == 1 && reco_particleClearCosmic->at(pa) == 1) continue;
                
                numPFPsSlice++;
                
                if(reco_particleIsPrimary->at(pa) == 1){
                    numPrimaryPFPsSlice++;
                    if(reco_particleNumHits->at(pa) >= 10) numPrimaryPFPs10Slice++;
                }

                numHitsInPFPs += reco_particleNumHits->at(pa);
                summedEnergy += reco_particleBestPlaneEnergy->at(pa);
            }

            // Count reco neutrinos in this slice
            int numRecoNeutrinos = 0;
            for(size_t nu = 0; nu < reco_neutrinoSliceID->size(); nu++){
                if((int)reco_neutrinoSliceID->at(nu) == (int)reco_sliceID->at(slice)) numRecoNeutrinos++;
            }

            // Find highest energy PFP
            highestEnergyPFP_struct hePFP;
            for(size_t pa = 0; pa < reco_particleSliceID->size(); pa++){
                if((int)reco_particleSliceID->at(pa) != (int)reco_sliceID->at(slice)) continue;
                if(thisCut.clearCosmicCut == 1 && reco_particleClearCosmic->at(pa) == 1) continue;

                double en = reco_particleBestPlaneEnergy->at(pa);
                
                if(en > hePFP.energy){
                    hePFP.energy               = en;
                    hePFP.PFPID                = reco_particleID->at(pa);
                    hePFP.theta                = reco_particleTheta->at(pa);
                    hePFP.dx                   = reco_particleDX->at(pa);
                    hePFP.dy                   = reco_particleDY->at(pa);
                    hePFP.dz                   = reco_particleDZ->at(pa);
                    hePFP.vx                   = reco_particleVX->at(pa);
                    hePFP.vy                   = reco_particleVY->at(pa);
                    hePFP.vz                   = reco_particleVZ->at(pa);
                    hePFP.completeness         = reco_particleCompleteness->at(pa);
                    hePFP.purity               = reco_particlePurity->at(pa);
                    hePFP.trackscore           = reco_particleTrackScore->at(pa);
                    hePFP.primary              = reco_particleIsPrimary->at(pa);
                    hePFP.truePDG              = reco_particleTruePDG->at(pa);
                    hePFP.trueOrigin           = reco_particleTrueOrigin->at(pa);
                    hePFP.trueInt              = reco_particleTrueInteractionType->at(pa);
                    hePFP.bestPlanedEdx        = reco_particleBestPlanedEdx->at(pa);
                    hePFP.razzledPDG11         = reco_particleRazzledPDG11->at(pa);
                    hePFP.razzledPDG13         = reco_particleRazzledPDG13->at(pa);
                    hePFP.razzledPDG22         = reco_particleRazzledPDG22->at(pa);
                    hePFP.razzledPDG211        = reco_particleRazzledPDG211->at(pa);
                    hePFP.razzledPDG2212       = reco_particleRazzledPDG2212->at(pa);
                    hePFP.razzledBestPDG       = reco_particleRazzledBestPDG->at(pa);
                    hePFP.showerLength         = reco_particleShowerLength->at(pa);
                    hePFP.showerOpenAngle      = reco_particleShowerOpenAngle->at(pa);
                    hePFP.showerBestPlaneEnergy= reco_particleShowerBestPlaneEnergy->at(pa);
                    hePFP.trueVX               = reco_particleTrueVX->at(pa);
                    hePFP.trueVY               = reco_particleTrueVY->at(pa);
                    hePFP.trueVZ               = reco_particleTrueVZ->at(pa);
                    hePFP.trueEndX             = reco_particleTrueEndX->at(pa);
                    hePFP.trueEndY             = reco_particleTrueEndY->at(pa);
                    hePFP.trueEndZ             = reco_particleTrueEndZ->at(pa);
                    hePFP.numHits              = reco_particleNumHits->at(pa);
                    hePFP.clearCosmic          = reco_particleClearCosmic->at(pa);
                }
            }

            // Reco vertex
            double recoVX = -999999, recoVY = -999999, recoVZ = -999999;
            for(size_t nu = 0; nu < reco_neutrinoSliceID->size(); nu++){
                if((int)reco_neutrinoSliceID->at(nu) == (int)reco_sliceID->at(slice)){
                    recoVX = reco_neutrinoVX->at(nu);
                    recoVY = reco_neutrinoVY->at(nu);
                    recoVZ = reco_neutrinoVZ->at(nu);
                }
            }

            // PFP 10cm PCA angle
            double pfp10cm_PCAAngle = -999999;
            for(size_t pa = 0; pa < angleRecalculationPCAPFP10cm_pfpID->size(); pa++){
                if(angleRecalculationPCAPFP10cm_pfpID->at(pa) == hePFP.PFPID){
                    pfp10cm_PCAAngle = angleRecalculationPCAPFP10cm_angle->at(pa);
                    break;
                }
            }

            // ----- Apply cuts -----
            const CutConfig& cc = thisCut;
            if(cc.numPFPs0Cut == 1 && numPFPsSlice == 0) continue;

            if(cc.numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0) continue;

            if(cc.CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high)) continue;

            if(cc.FVCut == 1){
                if(!(recoVX < FVCut_xHigh && recoVX > FVCut_xLow && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow  && recoVZ < FVCut_zHigh)) continue;
            }

            if(cc.primaryPFPCut == 1 && numPrimaryPFPs10Slice != primaryPFPCutValue) continue;
            
            if(cc.ETheta2Cut == 1 && (hePFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle > ETheta2High_highestEnergyPFP || hePFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle < ETheta2Low_highestEnergyPFP)) continue;
            
            if(cc.razzledPDG11Cut == 1 && (hePFP.razzledPDG11  > razzled11High_highestEnergyPFP || hePFP.razzledPDG11  < razzled11Low_highestEnergyPFP)) continue;
            
            if(cc.razzledPDG211Cut == 1 && (hePFP.razzledPDG211 > razzled211High_highestEnergyPFP || hePFP.razzledPDG211 < razzled211Low_highestEnergyPFP)) continue;
            
            if(cc.dEdxCut == 1 && (hePFP.bestPlanedEdx > dEdxHigh_highestEnergyPFP || hePFP.bestPlanedEdx < dEdxLow_highestEnergyPFP)) continue;

            // ----- Slice category -----
            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                sliceCategoryPlottingMacro = 0;
            } else if(reco_sliceOrigin->at(slice) == 1){
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    if(cc.FVCut == 0 &&
                       reco_sliceTrueVX->at(slice) < 201.3  && reco_sliceTrueVX->at(slice) > -201.3 &&
                       reco_sliceTrueVY->at(slice) < 203.8  && reco_sliceTrueVY->at(slice) > -203.8 &&
                       reco_sliceTrueVZ->at(slice) > 0      && reco_sliceTrueVZ->at(slice) < 509.5){
                        sliceCategoryPlottingMacro = 1;
                    } else if(cc.FVCut == 1 &&
                               reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow &&
                               std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre &&
                               reco_sliceTrueVY->at(slice) < FVCut_yHigh && reco_sliceTrueVY->at(slice) > FVCut_yLow &&
                               reco_sliceTrueVZ->at(slice) < FVCut_zHigh && reco_sliceTrueVZ->at(slice) > FVCut_zLow){
                        sliceCategoryPlottingMacro = 1;
                    } else {
                        sliceCategoryPlottingMacro = 2;
                    }
                } else {
                    sliceCategoryPlottingMacro = 2;
                }
            } else if(reco_sliceOrigin->at(slice) == 3){
                sliceCategoryPlottingMacro = (reco_sliceCompleteness->at(slice) > 0.5) ? 3.0 : 4.0;
            }
            if(sliceCategoryPlottingMacro == -999999) continue;

            int cat = (int)sliceCategoryPlottingMacro;
            bool isCosmic = (cat == 0);

            // For each of the 25 variables, the value for the slice is computed and it works out which bin it sits in
            for(int v = 0; v < nVars; v++){
                double val = getVariableValue(varConfigs[v].name, reco_sliceCompleteness->at(slice), reco_sliceScore->at(slice), reco_slicePurity->at(slice), numRecoNeutrinos, numPFPsSlice, numPrimaryPFPsSlice, numPrimaryPFPs10Slice, numHitsInPFPs, reco_sliceNumHits->at(slice), summedEnergy, hePFP, pfp10cm_PCAAngle, recoVX, recoVY, recoVZ);

                int nb   = varConfigs[v].nBins;
                double lo = varConfigs[v].min;
                double hi = varConfigs[v].max;
                int b = findBin(val, nb, lo, hi);

                // Nominal (universe index 0) - adds the POT weight
                for(int p = 0; p < nParams; p++)
                    bins[v][cat][p][0][b] += potWeight;

                // Loops through the universes and gets the flux weights for the 13 parameters (or 1 if its a cosmic)
                for(int u = 0; u < NUNIV; u++){
                    double w_horn = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_horncurrent, slice, u, weightsFound);
                    double w_exp = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_expskin, slice, u, weightsFound);
                    double w_kp = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_kplus, slice, u, weightsFound);
                    double w_km = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_kmin, slice, u, weightsFound);
                    double w_k0 = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_kzero, slice, u, weightsFound);
                    double w_nie = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_nucleoninexsec, slice, u, weightsFound);
                    double w_nqe = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_nucleonqexsec, slice, u, weightsFound);
                    double w_ntx = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_nucleontotxsec, slice, u, weightsFound);
                    double w_pim = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_piminus, slice, u, weightsFound);
                    double w_pie = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_pioninexsec, slice, u, weightsFound);
                    double w_piq = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_pionqexsec, slice, u, weightsFound);
                    double w_pit = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_piontotxsec, slice, u, weightsFound);
                    double w_pip = isCosmic ? 1.0 : getSliceWeight(reco_sliceMCTruthFlux_weight_piplus, slice, u, weightsFound);
                    double w_comb = w_horn*w_exp*w_kp*w_km*w_k0*w_nie*w_nqe*w_ntx*w_pim*w_pie*w_piq*w_pit*w_pip;

                    std::vector<double> sw = {w_horn,w_exp,w_kp,w_km,w_k0,w_nie,w_nqe,w_ntx,w_pim,w_pie,w_piq,w_pit,w_pip,w_comb};
                    for(int p = 0; p < nParams; p++)
                        
                        // Adds the POT weight * universe flux parameter weight
                        bins[v][cat][p][u + 1][b] += potWeight * sw[p];
                }
            }
        } 
    } 

    std::cout << "Event loop complete. Writing output TTree..." << std::endl;

    // ============================================================
    // WRITE OUTPUT
    //
    // TTree named after the cut (e.g. "clearCosmic").
    // One entry (one Fill() call).
    //
    // Branch naming: {varName}_{catName}_{paramName}
    // Branch type:   vector<vector<double>>
    //   [0]         = nominal   (nBins entries)
    //   [1..1000]   = universe 0..999 (nBins entries each)
    //
    // A second TTree "sigCount" holds signal counts:
    // Branch naming: sigCount_{paramName}
    // Branch type:   vector<double>  (length NUNIV+1)
    //   [0]         = nominal signal count (actualSignalCount)
    //   [1..1000]   = per-universe signal count for universes 0..999
    // ============================================================

    // Opens the output ROOT file
    TFile *fOut = TFile::Open(outputFile.c_str(), "RECREATE");
    if(!fOut){ std::cerr << "Could not open output file" << std::endl; return; }

    fOut->cd();

    // New TTree is created and named after the cuts being applied
    TTree *outTree = new TTree(thisCut.name.c_str(), ("Histogram bin contents for cut: " + thisCut.name).c_str());

    // Calculates the number of branches that are needed, 1 for each variable, slice category and parameter
    int nBranches = nVars * nCats * nParams;
    std::vector<std::vector<std::vector<double>>*> branchPtrs(nBranches, nullptr);

    // This block of code creates branches. Each branch will hold 1001 histograms for one variable/category/parameter
    // index 0 = nominal, 1-1000 = 1000 universes
    //
    // Loop through the number of variables to plot
    for(int v = 0; v < nVars; v++){
        // Loop through the number of slice categories
        for(int cat = 0; cat < nCats; cat++){
            // Loop through the number of flux parameters
            for(int p = 0; p < nParams; p++){
                int idx = v * nCats * nParams + cat * nParams + p;
                // Put the data into the 2D vector
                branchPtrs[idx] = new std::vector<std::vector<double>>(bins[v][cat][p]);

                // Construct name of the branch
                std::string bname = varConfigs[v].name + "_" + catNames[cat] + "_" + paramNames[p];

                // Attaches the 2D vector as a branch of outTree
                outTree->Branch(bname.c_str(), branchPtrs[idx]);
            }
        }
    }

    // Another treewhich holds the total number of true nu+e elastic scattering events per systematic universe
    TTree *sigTree = new TTree("sigCount", "Per-universe signal counts");
    std::vector<std::vector<double>*> sigBranchPtrs(nParams, nullptr);
    for(int p = 0; p < nParams; p++){
        sigBranchPtrs[p] = new std::vector<double>(NUNIV + 1, 0.0);

        (*sigBranchPtrs[p])[0] = actualSignalCount; // nominal
        
        // Loop through universes
        for(int u = 0; u < NUNIV; u++) (*sigBranchPtrs[p])[u + 1] = sigCountUniv[p][u];
        
        // Construct name of the branch
        std::string bname = "sigCount_" + paramNames[p];

        // Attaches the 2D vector as a branch of sigTree
        sigTree->Branch(bname.c_str(), sigBranchPtrs[p]);
    }

    // Also stores the nominal number of true nu+e elastic scattering events as its own branch
    double nomSigCount = actualSignalCount;
    sigTree->Branch("nominalSignalCount", &nomSigCount);

    // Single Fill for both trees
    outTree->Fill();
    sigTree->Fill();

    // ----- Write systematic summary to text file -----
    std::string summaryPath = "/nashome/c/coackley/systSummary_" + thisCut.name + ".txt";
    std::ofstream systSummary(summaryPath.c_str());
    systSummary << "=== Systematic Uncertainties on nu+e Signal Count ===" << std::endl;
    systSummary << "\nCut: " << thisCut.name << std::endl;
    systSummary << Form("  Nominal: %.2f", actualSignalCount) << std::endl;
    double totalSystSq = 0;
    for(int p = 0; p < nParams - 1; p++){
        // Compute mean and stddev from sigCountUniv[p]
        double sum = 0, sum2 = 0;
        for(int u = 0; u < NUNIV; u++){ sum += sigCountUniv[p][u]; sum2 += sigCountUniv[p][u]*sigCountUniv[p][u]; }
        double mean   = sum / NUNIV;
        double stddev = std::sqrt(sum2/NUNIV - mean*mean);
        double shift  = mean - actualSignalCount;
        systSummary << Form("  %-20s  mean=%.2f  shift=%+.2f (%+.1f%%)  syst=%.2f (%.1f%%)",
            paramNames[p].c_str(), mean, shift,
            actualSignalCount > 0 ? 100.*shift/actualSignalCount : 0,
            stddev,
            actualSignalCount > 0 ? 100.*stddev/actualSignalCount : 0) << std::endl;
        totalSystSq += stddev * stddev;
    }
    double totalSyst = std::sqrt(totalSystSq);
    systSummary << Form("  %-20s  syst=%.2f (%.1f%%)", "TOTAL (quadrature)", totalSyst,
        actualSignalCount > 0 ? 100.*totalSyst/actualSignalCount : 0) << std::endl;
    int pComb = nParams - 1;
    {
        double sum = 0, sum2 = 0;
        for(int u = 0; u < NUNIV; u++){ sum += sigCountUniv[pComb][u]; sum2 += sigCountUniv[pComb][u]*sigCountUniv[pComb][u]; }
        double combMean = sum / NUNIV;
        double combSyst = std::sqrt(sum2/NUNIV - combMean*combMean);
        systSummary << Form("  %-20s  mean=%.2f  syst=%.2f (%.1f%%)", "COMBINED (product)",
            combMean, combSyst,
            actualSignalCount > 0 ? 100.*combSyst/actualSignalCount : 0) << std::endl;
    }
    systSummary.close();

    // Store actual signal count as TParameter
    TParameter<double> *pSigCount = new TParameter<double>("actualSignalCount", actualSignalCount);
    pSigCount->Write();

    fOut->Write();
    fOut->Close();

    // Clean up heap allocations
    for(auto ptr : branchPtrs) delete ptr;
    for(auto ptr : sigBranchPtrs) delete ptr;

    fNuE->Close();
    fNuEWeights->Close();

    std::cout << "\nDone! Output written to: " << outputFile << std::endl;
    std::cout << "Systematic summary written to: " << summaryPath << std::endl;
}
