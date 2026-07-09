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
#include <functional>
#include <array>
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

struct eventCounter_struct{
    double nuE = 0;
    double NCNPi0 = 0;
    double otherNC = 0;
    double CCnumu = 0;
    double CCnue = 0;
    double dirt = 0;
    double nuEDirt = 0;
    double nuEFuzzy = 0;
    double cosmic = 0;
    double other = 0;
};

struct beforeEventCount_struct{
    double signal = 0;
    double background = 0;
    eventCounter_struct splitInt;
};

struct eventCounting_struct{
    double clearCosmicsSig = 0;    double clearCosmicsBack = 0;    eventCounter_struct clearCosmicsIntSplit;
    double numPFPs0Sig = 0;        double numPFPs0Back = 0;        eventCounter_struct numPFPs0IntSplit;
    double numRecoNeut0Sig = 0;    double numRecoNeut0Back = 0;    eventCounter_struct numRecoNeut0IntSplit;
    double FVSig = 0;              double FVBack = 0;              eventCounter_struct FVIntSplit;
    double crumbsSig = 0;          double crumbsBack = 0;          eventCounter_struct crumbsIntSplit;
    double primaryPFPSig = 0;      double primaryPFPBack = 0;      eventCounter_struct primaryPFPIntSplit;
    double razzled11Sig = 0;       double razzled11Back = 0;       eventCounter_struct razzled11IntSplit;
    double razzled211Sig = 0;      double razzled211Back = 0;      eventCounter_struct razzled211IntSplit;
    double ETheta2Sig = 0;         double ETheta2Back = 0;         eventCounter_struct ETheta2IntSplit;
    double dEdxSig = 0;            double dEdxBack = 0;            eventCounter_struct dEdxIntSplit;
};

struct GenieParam_struct{
    std::string mapKey;
    std::string shortName;
    int nUniv;
    bool isMultisim;
};

void nuESelectionNumbersWithXSecSystematics_macro(){

    // Set true to make all plots after each cut has been applied
    bool makePerCutPlots_GENIE = false;

    std::string cutsApplied = "allCuts";
    std::string base_path = "/nashome/c/coackley/systPlotsNumbers9July_GENIE_" + cutsApplied + "/";
    std::string tableFileName = base_path + "table_GENIE.txt";

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

    // AV boundaries
    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0;      double zMax = 509.4;

    if (gSystem->AccessPathName(base_path.c_str())) {
        gSystem->mkdir(base_path.c_str(), kTRUE);
    }

    std::ofstream clearTableFile(tableFileName, std::ios::trunc);
    if (!clearTableFile.is_open()) {
        std::cerr << "Error: could not open or create " << tableFileName << std::endl;
        return;
    }
    clearTableFile.close();

    TFile *fNuE = TFile::Open("/exp/sbnd/data/users/coackley/merged_nu+eIntimeCosmicBNBnue_7Jul_noWeights.root");
    if(!fNuE){ std::cerr << "Error opening the NuE TTree file" << std::endl; return; }

    TDirectory *dirNuE = (TDirectory*)fNuE->Get("ana");
    if(!dirNuE){ std::cerr << "Directory 'ana' not found in NuE file" << std::endl; return; }

    TTree *tree = (TTree*)dirNuE->Get("NuE");
    if(!tree){ std::cerr << "NuE TTree not found" << std::endl; return; }

    TTree *subRunTree = (TTree*)dirNuE->Get("SubRun");
    if(!subRunTree){ std::cerr << "SubRun TTree not found" << std::endl; return; }

    std::vector<std::string> weightFilePaths = {
        "/exp/sbnd/data/users/coackley/merged_nu+eAllJobs_weights_3Jul.root",
        "/exp/sbnd/data/users/coackley/nue_withWeights_5Jul.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_000.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_001.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_002.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_003.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_004.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_005.root",
        "/exp/sbnd/data/users/coackley/hadd_level3/level3_006.root"
    };
    
    const int NWEIGHTFILES = (int)weightFilePaths.size();

    std::vector<TFile*> fNuEWeightsVec(NWEIGHTFILES, nullptr);
    std::vector<TTree*> weightsTreeVec(NWEIGHTFILES, nullptr);

    for(int f = 0; f < NWEIGHTFILES; ++f){
        fNuEWeightsVec[f] = TFile::Open(weightFilePaths[f].c_str());
        if(!fNuEWeightsVec[f] || fNuEWeightsVec[f]->IsZombie()){
            std::cerr << "Error opening weights file: " << weightFilePaths[f] << std::endl;
            return;
        }

        TDirectory *dirNuEWeights = (TDirectory*)fNuEWeightsVec[f]->Get("ana");
        if(!dirNuEWeights){
            std::cerr << "Directory 'ana' not found in " << weightFilePaths[f] << std::endl;
            return;
        }

        weightsTreeVec[f] = (TTree*)dirNuEWeights->Get("NuEWeights");
        if(!weightsTreeVec[f]){
            std::cerr << "NuEWeights TTree not found in " << weightFilePaths[f] << std::endl;
            return;
        }
    }

    TFile *fOut = new TFile("/exp/sbnd/data/users/coackley/selectionNumberSystematicPlots_GENIE_9July.root", "RECREATE");
    if(!fOut || fOut->IsZombie()){ std::cerr << "Error creating output ROOT file" << std::endl; return; }

    beforeEventCount_struct eventsBeforeCuts_DLNuE;
    eventCounting_struct eventsAfterCuts_DLNuE;

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
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsNuENuE;

    double cosmicSpillsSumNuE = 0;
    double BNBSpillsSumNuE = 0;
    double NuESpillsSumNuE = 0;

    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;
    double totalPOTNuENuE = 0;   
 
    double POTSignalNuE_notMissing = 0;
    double POTBNBNuE_notMissing = 0;
    double POTNuENuE_notMissing = 0;

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

        } else if(subRunSignal == 4){
            if(subRunDLCurrent == 5 && seenSubRunsNuENuE.find(key) == seenSubRunsNuENuE.end()){
                totalPOTNuENuE += subRunPOT;
                seenSubRunsNuENuE.insert(key);
            }

            if(subRunDLCurrent == 5) POTNuENuE_notMissing += subRunPOT;
        }
    }

    double targetPOT = 1e21;
    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    double totalPOTNuENuE_notMissing = (POTNuENuE_notMissing + POTBNBNuE_notMissing);

    std::cout << "POT from nue sample = " << POTNuENuE_notMissing << ", POT from BNB sample = " << POTBNBNuE_notMissing << ", total nue POT = " << totalPOTNuENuE_notMissing << std::endl;

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

    std::vector<double> *truth_particleSliceID = nullptr;
    std::vector<double> *truth_particlePDG = nullptr;
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
    std::vector<double> *angleRecalculationPCAPFP10cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dz = nullptr;

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

    tree->SetBranchAddress("truth_particleSliceID", &truth_particleSliceID);
    tree->SetBranchAddress("truth_particlePDG", &truth_particlePDG);
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
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dx", &angleRecalculationPCAPFP10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dy", &angleRecalculationPCAPFP10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dz", &angleRecalculationPCAPFP10cm_dz);

    Long64_t numEntries = tree->GetEntries();

    // NuEWeights Tree branch variable
    UInt_t eventID_weights, runID_weights, subRunID_weights;
    int nuEScatter_weights;
    double nuEScatterTrueVX_weights, nuEScatterTrueVY_weights, nuEScatterTrueVZ_weights;
    double DLCurrent_weights, signal_weights;

    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE = nullptr;
    std::vector<double> *nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE = nullptr;

    std::vector<double> *reco_sliceID_weights = nullptr;  
    std::vector<double> *reco_sliceInteraction_weights = nullptr;  
    std::vector<double> *reco_sliceTrueVX_weights = nullptr;  
    std::vector<double> *reco_sliceTrueVY_weights = nullptr;  
    std::vector<double> *reco_sliceTrueVZ_weights = nullptr;  
    std::vector<double> *reco_sliceOrigin_weights = nullptr;  
    std::vector<double> *reco_sliceTrueCCNC_weights = nullptr;  
    std::vector<double> *reco_sliceTrueNeutrinoType_weights = nullptr;

    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE = nullptr;
    std::vector<std::vector<double>> *reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE = nullptr;

    for(int f = 0; f < NWEIGHTFILES; ++f){
        TTree* wt = weightsTreeVec[f];
        wt->SetBranchAddress("eventID", &eventID_weights);
        wt->SetBranchAddress("runID", &runID_weights);
        wt->SetBranchAddress("subRunID", &subRunID_weights);
        wt->SetBranchAddress("DLCurrent", &DLCurrent_weights);
        wt->SetBranchAddress("signal", &signal_weights);
        
        wt->SetBranchAddress("nuEScatter", &nuEScatter_weights);
        wt->SetBranchAddress("nuEScatterTrueVX", &nuEScatterTrueVX_weights);
        wt->SetBranchAddress("nuEScatterTrueVY", &nuEScatterTrueVY_weights);
        wt->SetBranchAddress("nuEScatterTrueVZ", &nuEScatterTrueVZ_weights);

        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE);
        wt->SetBranchAddress("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE);
        
        wt->SetBranchAddress("reco_sliceID", &reco_sliceID_weights);
        wt->SetBranchAddress("reco_sliceInteraction", &reco_sliceInteraction_weights);
        wt->SetBranchAddress("reco_sliceTrueVX", &reco_sliceTrueVX_weights);
        wt->SetBranchAddress("reco_sliceTrueVY", &reco_sliceTrueVY_weights);
        wt->SetBranchAddress("reco_sliceTrueVZ", &reco_sliceTrueVZ_weights);
        wt->SetBranchAddress("reco_sliceOrigin", &reco_sliceOrigin_weights);
        wt->SetBranchAddress("reco_sliceTrueCCNC", &reco_sliceTrueCCNC_weights);
        wt->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType_weights);

        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE);
        wt->SetBranchAddress("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE);
    }

    struct weightLocation_struct{
        int fileIdx;
        Long64_t entry;
    };

    std::unordered_map<eventKey_struct, weightLocation_struct, eventKeyHash_struct> weightEntryMap;

    for(int f = 0; f < NWEIGHTFILES; ++f){
        TTree* wt = weightsTreeVec[f];
        Long64_t nEntriesThisTree = wt->GetEntries();

        for(Long64_t i = 0; i < nEntriesThisTree; ++i){
            wt->GetEntry(i);
            eventKey_struct key{runID_weights, subRunID_weights, eventID_weights, static_cast<int>(signal_weights), static_cast<int>(DLCurrent_weights)};

            if(weightEntryMap.find(key) != weightEntryMap.end()){
                std::cerr << "Warning: duplicate event key across weight files (run=" << runID_weights
                           << ", subRun=" << subRunID_weights << ", event=" << eventID_weights
                           << ", signal=" << signal_weights << ", DLCurrent=" << DLCurrent_weights
                           << ") - keeping first occurrence" << std::endl;
                continue;
            }

            weightEntryMap[key] = {f, i};
        }
    }

    std::vector<GenieParam_struct> genieParams = {
        // NOvA-style non-resonant pion normalization (23 params)
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

        // Misc interaction systematics (5 params)
        {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu","Misc_C12ToAr40_2p2hScaling_nu",6,false},
        {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar","Misc_C12ToAr40_2p2hScaling_nubar",6,false},
        {"MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression","Misc_SPPLowQ2Suppression",10,false},
        {"MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio","Misc_nuenuebar_xsec_ratio",2,false},
        {"MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio","Misc_nuenumu_xsec_ratio",2,false},

        // MINERvA q0q3 weighting (1 param)
        {"MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement","MINERvAq0q3_Mnv2p2hGaussEnhancement",4,false},

        // MINERvA E2p2h (4 params)
        {"MINERvAE2p2h_SBN_v1_E2p2h_A_nu","MINERvAE2p2h_E2p2h_A_nu",6,false},
        {"MINERvAE2p2h_SBN_v1_E2p2h_A_nubar","MINERvAE2p2h_E2p2h_A_nubar",6,false},
        {"MINERvAE2p2h_SBN_v1_E2p2h_B_nu","MINERvAE2p2h_E2p2h_B_nu",6,false},
        {"MINERvAE2p2h_SBN_v1_E2p2h_B_nubar","MINERvAE2p2h_E2p2h_B_nubar",6,false},

        // GENIEReWeight multisim (30 params, N=100 -> genuine multisim / RMS treatment)
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

        // GENIEReWeight multisigma (52 params -> envelope treatment)
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

    const int NPARAMS_GENIE = (int)genieParams.size(); // should be 115
    std::cout << "Number of GENIE parameters loaded: " << NPARAMS_GENIE << " out of 115" << std::endl;

    // Order must exactly match genieParams above, index-for-index
    std::vector<std::vector<double>*> nuEScatter_GENIEWeightVecs = {
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi,
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi,
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu,
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar,
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression,
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio,
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio,
        nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement,
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu,
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar,
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu,
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE,
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE
    };

    std::vector<std::vector<std::vector<double>>*> reco_slice_GENIEWeightVecs = {
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

    if((int)nuEScatter_GENIEWeightVecs.size() != NPARAMS_GENIE || (int)reco_slice_GENIEWeightVecs.size() != NPARAMS_GENIE){
        std::cerr << "ERROR: GENIE weight-vector array size doesn't match genieParams (" << NPARAMS_GENIE << ")" << std::endl;
        return;
    }

    const int NCUTS = 11;
    std::vector<std::string> cutNames_syst = {"beforeCuts", "clearCosmic", "numPFPs0", "numRecoNeut", "crumbs", "FV", "primaryPFP", "razzled11", "razzled211", "ETheta2", "dEdx"};

    std::vector<double> nomSig_perCut(NCUTS, 0.0);
    std::vector<double> nomBack_perCut(NCUTS, 0.0);

    std::vector<std::vector<double>> count_genie(NPARAMS_GENIE); // count_genie[p][u]: true nu+e count (before cuts) per universe
    std::vector<std::vector<std::vector<double>>> univSig_perCutParam_genie(NPARAMS_GENIE);
    std::vector<std::vector<std::vector<double>>> univBack_perCutParam_genie(NPARAMS_GENIE);

    for(int p = 0; p < NPARAMS_GENIE; p++){
        count_genie[p].assign(genieParams[p].nUniv, 0.0);
        univSig_perCutParam_genie[p].assign(NCUTS, std::vector<double>(genieParams[p].nUniv, 0.0));
        univBack_perCutParam_genie[p].assign(NCUTS, std::vector<double>(genieParams[p].nUniv, 0.0));
    }

    double actualSignalCount = 0.0;

    auto getGenieNuEWeight = [&](std::vector<double>* vec, int u, int expectedN) -> double {
        if(!vec || (int)vec->size() != expectedN) return 1.0;
        if(u < 0 || u >= expectedN) return 1.0;
        return vec->at(u);
    };

    auto getGenieSliceWeight = [&](std::vector<std::vector<double>>* vec, size_t sliceIdx, int u, bool wFound, int expectedN) -> double {
        if(!wFound || !vec || sliceIdx >= vec->size()) return 1.0;
        if((int)vec->at(sliceIdx).size() != expectedN) return 1.0;
        if(u < 0 || u >= expectedN) return 1.0;
        return vec->at(sliceIdx).at(u);
    };

    // Multisim: RMS about nominal, 1/N 
    auto calcSystMultisim = [&](const std::vector<double>& values, double nominal) -> double {
        int N = (int)values.size();
        if(N < 2) return 0.0;
        double sumSq = 0.0;
        for(double x : values) sumSq += (x - nominal)*(x - nominal);
        return std::sqrt(sumSq/N);
    };

    // Multisigma/unisim: envelope estimator (largest absolute deviation among the available points).
    auto calcSystMultisigma = [&](const std::vector<double>& values, double nominal) -> double {
        if(values.empty()) return 0.0;
        double maxDev = 0.0;
        for(double x : values) maxDev = std::max(maxDev, std::fabs(x - nominal));
        return maxDev;
    };

    // Multisigma with the standard 6-point layout [-1sig, +1sig, -2sig, +2sig, -3sig, +3sig].
    // Uses only the +-1sigma points, symmetrized, as is standard practice.
    auto calcSystMultisigma1Sigma = [&](const std::vector<double>& values, double nominal) -> double {
        if(values.size() != 6) return calcSystMultisigma(values, nominal); // fallback for non-standard knobs (e.g. the 1-point unisim knobs)
        double down1sigma = values[0];
        double up1sigma    = values[1];
        return std::fabs(up1sigma - down1sigma) / 2.0;
    };

    auto calcSystGeneric = [&](const std::vector<double>& values, double nominal, bool isMultisim) -> double {
        return isMultisim ? calcSystMultisim(values, nominal) : calcSystMultisigma1Sigma(values, nominal);
    };

    auto calcMeanFromValues = [&](const std::vector<double>& values) -> double {
        if(values.empty()) return 0.0;
        double sum = 0.0;
        for(double x : values) sum += x;
        return sum / values.size();
    };

    // Loop through events
    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);

        bool weightsFound = false;
        Long64_t weightsEntryIdx = -1;

        if(signal == 1 || signal == 2 || signal == 4){
            eventKey_struct key{runID, subRunID, eventID, static_cast<int>(signal), static_cast<int>(DLCurrent)};
            auto it = weightEntryMap.find(key);
            if(it != weightEntryMap.end()){
                weightsTreeVec[it->second.fileIdx]->GetEntry(it->second.entry);
                weightsFound = true;
                weightsEntryIdx = it->second.entry;
            }
        }

        recoilElectron_struct recoilElectron;
        recoilElectron.energy = -999999; recoilElectron.dx = -999999; recoilElectron.dy = -999999; recoilElectron.dz = -999999;
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->at(i) != -999999){
                recoilElectron.energy = truth_recoilElectronEnergy->at(i);
                recoilElectron.dx = truth_recoilElectronDX->at(i);
                recoilElectron.dy = truth_recoilElectronDY->at(i);
                recoilElectron.dz = truth_recoilElectronDZ->at(i);
            }
        }

        // True nu+e elastic scatter: fill count_genie[p][u]
        if(nuEScatter == 1 && signal == 1 && DLCurrent == 5){
            if(recoilElectron.energy > 150){
                bool passesFV = (FVCut == 0 && (((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))) || (FVCut == 1 && (((nuEScatterTrueVX > FVCut_xLow) && (nuEScatterTrueVX < FVCut_xHigh) && (std::abs(nuEScatterTrueVX) > FVCut_xCentre)) && ((nuEScatterTrueVY > FVCut_yLow) && (nuEScatterTrueVY < FVCut_yHigh)) && ((nuEScatterTrueVZ > FVCut_zLow) && (nuEScatterTrueVZ < FVCut_zHigh))));
                if(passesFV){
                    actualSignalCount += weights.signalNuE;

                    if(weightsFound){
                        for(int p = 0; p < NPARAMS_GENIE; p++){
                            int N = genieParams[p].nUniv;
                            for(int u = 0; u < N; u++){
                                double w = getGenieNuEWeight(nuEScatter_GENIEWeightVecs[p], u, N);
                                count_genie[p][u] += weights.signalNuE * w;
                            }
                        }
                    }
                }
            }
        }

        if(reco_sliceID->size() == 0) continue;

        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) == -999999) continue;

            // Assigning a category to the slices
            // 0 = cosmic, 1 = signal, 2 = signal fuzzy, 3 = bnb, 4 = bnb fuzzy, 5 = nue event, 6 = nue fuzzy event
            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                sliceCategoryPlottingMacro = 0;
            } else if(reco_sliceOrigin->at(slice) == 1){
                if(reco_sliceCompleteness->at(slice) > 0.5 && recoilElectron.energy > 150){
                    // Slice must have completeness > 0.5 and nu+e elastic scatter it comes from has true recoil electron energy > 150 MeV
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
                // This is a BNB slice
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    if(reco_sliceTrueCCNC->at(slice) == 0 && reco_sliceTrueNeutrinoType->at(slice) == 12){
                        // This is a slice from a nue event
                        sliceCategoryPlottingMacro = 5;
                    } else{
                        // This is a BNB event (not a nue event)
                        sliceCategoryPlottingMacro = 3;
                        //std::cout << "BNB Slice: sliceCategoryPlottingMacro = 3" << std::endl;
                    }
                } else{
                    if(reco_sliceTrueCCNC->at(slice) == 0 && reco_sliceTrueNeutrinoType->at(slice) == 12){
                        sliceCategoryPlottingMacro = 6;
                    } else{
                        sliceCategoryPlottingMacro = 4;
                        //std::cout << "BNB Fuzzy Slice: sliceCategoryPlottingMacro = 4" << std::endl;
                    }
                }
            }


            // Assigning an interaction category to the slices
            // Event types: Cosmic = 0, nu+e scatter = 1, NC Npi0 = 2, other NC = 3, CC numu = 4, CC nue = 5, Dirt = 6, Dirt nu+e = 7
            // Other = 8, Fuzzy nu+e = 9
            int sliceInteractionType = -999999;
            if(reco_sliceOrigin->at(slice) != 0){
                if(reco_sliceOrigin->at(slice) == 1){
                    if(reco_sliceCompleteness->at(slice) > 0.5 && recoilElectron.energy > 150){
                        // Slice must have completeness > 0.5 and nu+e elastic scatter it comes from has true recoil electron energy > 150 MeV
                        if(FVCut == 0 && (reco_sliceTrueVX->at(slice) > -201.3 && reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVY->at(slice) > -203.8 && reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5)){
                            sliceInteractionType = 1;
                        } else if(FVCut == 1 && ((reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow && std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre) && (reco_sliceTrueVY->at(slice) > FVCut_yLow && reco_sliceTrueVY->at(slice) < FVCut_yHigh) && (reco_sliceTrueVZ->at(slice) > FVCut_zLow && reco_sliceTrueVZ->at(slice) < FVCut_zHigh))){
                            sliceInteractionType = 1;
                        } else{
                            sliceInteractionType = 7;
                        }
                    } else{
                        sliceInteractionType = 9;
                    }
                } else if(reco_sliceOrigin->at(slice) == 3){
                    if((FVCut == 0 && (reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVX->at(slice) > -201.3) && (reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVY->at(slice) > -203.8) && (reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5)) || (FVCut == 1 && (reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow && std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre) && (reco_sliceTrueVY->at(slice) < FVCut_yHigh && reco_sliceTrueVY->at(slice) > FVCut_yLow) && (reco_sliceTrueVZ->at(slice) < FVCut_zHigh && reco_sliceTrueVZ->at(slice) > FVCut_zLow))){
                        if(reco_sliceTrueCCNC->at(slice) == 0){
                            if(reco_sliceTrueNeutrinoType->at(slice) == 12){
                                sliceInteractionType = 5;
                            } else if(reco_sliceTrueNeutrinoType->at(slice) == 14){
                                sliceInteractionType = 4;
                            }
                        } else if(reco_sliceTrueCCNC->at(slice) == 1){
                            int neutralPion = 0;
                            for(size_t trueParticle = 0; trueParticle < truth_particleSliceID->size(); trueParticle++){
                                if(truth_particleSliceID->at(trueParticle) == reco_sliceID->at(slice)){
                                    if(truth_particleStatusCode->at(trueParticle) == 1){
                                        if(truth_particlePDG->at(trueParticle) == 111){
                                            neutralPion++;
                                        }
                                    }
                                }
                            }
                            if(neutralPion > 0){
                                sliceInteractionType = 2;
                            } else{
                                sliceInteractionType = 3;
                            }
                        }
                    } else{
                        sliceInteractionType = 6;
                    }
                }
            } else{
                sliceInteractionType = 0;
            }

            if(sliceInteractionType == -999999){
                sliceInteractionType = 8;
            }
            
            double weight = 0;
            if(signal == 1 && DLCurrent == 5) weight = weights.signalNuE;
            if(signal == 2 && DLCurrent == 5 && sliceCategoryPlottingMacro != 5 && sliceCategoryPlottingMacro != 6) weight = weights.BNBNuE;
            if((signal == 2 || signal == 4) && DLCurrent == 5 && (sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6)) weight = weights.NuENuE;
            if(signal == 3 && DLCurrent == 5) weight = weights.cosmicsNuE;

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
                            highestEnergyPFP.dx = reco_particleDX->at(pfp);
                            highestEnergyPFP.dy = reco_particleDY->at(pfp);
                            highestEnergyPFP.dz = reco_particleDZ->at(pfp);
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
            std::vector<std::vector<double>> sliceUnivWeights_genie(NPARAMS_GENIE);
            for(int p = 0; p < NPARAMS_GENIE; p++) sliceUnivWeights_genie[p].assign(genieParams[p].nUniv, 1.0);

            if(DLCurrent == 5 && weightsFound && signal != 3){
                // Match this slice to its entry in the weights tree by sliceID
                for(size_t ws = 0; ws < reco_sliceID_weights->size(); ++ws){
                    if(reco_sliceID_weights->at(ws) == reco_sliceID->at(slice)){
                        wSliceIdx_cached = ws;
                        sliceWeightValid_cached = true;
                        break;
                    }
                }

                for(int p = 0; p < NPARAMS_GENIE; p++){
                    int N = genieParams[p].nUniv;
                    for(int u = 0; u < N; u++){
                        sliceUnivWeights_genie[p][u] = getGenieSliceWeight(reco_slice_GENIEWeightVecs[p], wSliceIdx_cached, u, sliceWeightValid_cached, N);
                    }
                }
            }

            auto fillSliceSystCounters_genie = [&](int cutIdx){
                // Cat == 0 && signal == 4, cosmic slices from nue sample - skip explicitly (instead of relying on weight = 0)
                if(sliceCategoryPlottingMacro == 0 && signal == 4) return;

                bool isSigSlice = (sliceCategoryPlottingMacro == 1 && signal == 1);
                nomSig_perCut[cutIdx] += isSigSlice ? weight : 0.0;
                nomBack_perCut[cutIdx] += isSigSlice ? 0.0 : weight;

                for(int p = 0; p < NPARAMS_GENIE; p++){
                    int N = genieParams[p].nUniv;
                    for(int u = 0; u < N; u++){
                        double w = weight * sliceUnivWeights_genie[p][u];
                        if(isSigSlice) univSig_perCutParam_genie[p][cutIdx][u] += w;
                        else univBack_perCutParam_genie[p][cutIdx][u] += w;
                    }
                }
            };

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsBeforeCuts_DLNuE.background += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsBeforeCuts_DLNuE.signal += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsBeforeCuts_DLNuE.background += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsBeforeCuts_DLNuE.background += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsBeforeCuts_DLNuE.background += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsBeforeCuts_DLNuE.background += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsBeforeCuts_DLNuE.background += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsBeforeCuts_DLNuE.splitInt.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsBeforeCuts_DLNuE.splitInt.nuE += weight;
                else if(sliceInteractionType == 2) eventsBeforeCuts_DLNuE.splitInt.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsBeforeCuts_DLNuE.splitInt.otherNC += weight;
                else if(sliceInteractionType == 4) eventsBeforeCuts_DLNuE.splitInt.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsBeforeCuts_DLNuE.splitInt.CCnue += weight;
                else if(sliceInteractionType == 6) eventsBeforeCuts_DLNuE.splitInt.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsBeforeCuts_DLNuE.splitInt.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsBeforeCuts_DLNuE.splitInt.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy += weight;

                fillSliceSystCounters_genie(0);
            }

            if(DLCurrent == 5 && clearCosmicCut == 1){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
               
                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEFuzzy += weight; 

                fillSliceSystCounters_genie(1);
            }

            if(numPFPs0Cut == 1 && numPFPsSlice == 0) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.numPFPs0Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numPFPs0IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEFuzzy += weight;                

                fillSliceSystCounters_genie(2);
            }

            if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(3);
            }

            if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high)) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.crumbsSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.crumbsBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.crumbsIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.crumbsIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(4);
            }

            if(FVCut == 1){
                if(!(recoVX < FVCut_xHigh && recoVX > FVCut_xLow && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow && recoVZ < FVCut_zHigh)) continue;
            }

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.FVSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.FVBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.FVIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.FVIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.FVIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.FVIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.FVIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.FVIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(5);
            }

            if(primaryPFPCut == 1 && numPrimaryPFPs10Slice != primaryPFPCutValue) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.primaryPFPBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.primaryPFPIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(6);
            }

            if(razzledPDG11Cut == 1 && ((highestEnergyPFP.razzledPDG11 > razzled11High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG11 < razzled11Low_highestEnergyPFP))) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled11Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled11Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled11IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled11IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(7);
            }

            if(razzledPDG211Cut == 1 && ((highestEnergyPFP.razzledPDG211 > razzled211High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG211 < razzled211Low_highestEnergyPFP))) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled211Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled211Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled211IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled211IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(8);
            }

            if(ETheta2Cut == 1 && ((highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) > ETheta2High_highestEnergyPFP || (highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) < ETheta2Low_highestEnergyPFP)) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.ETheta2Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.ETheta2Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.ETheta2IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(9);
            }

            if(dEdxCut == 1 && (highestEnergyPFP.bestPlanedEdx > dEdxHigh_highestEnergyPFP || highestEnergyPFP.bestPlanedEdx < dEdxLow_highestEnergyPFP)) continue;

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.dEdxSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.dEdxBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.dEdxIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.dEdxIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuEFuzzy += weight;

                fillSliceSystCounters_genie(10);
            }
        }
    }

    std::cout << "\n=== GENIE Systematic Uncertainties on nu+e Signal Count (before cuts) ===" << std::endl;
    std::cout << Form("Nominal: %.2f", actualSignalCount) << std::endl;

    std::vector<double> systValues_beforeCuts(NPARAMS_GENIE, 0.0);

    for(int p = 0; p < NPARAMS_GENIE; p++){
        double stddev = calcSystGeneric(count_genie[p], actualSignalCount, genieParams[p].isMultisim);
        double mean   = calcMeanFromValues(count_genie[p]);
        double shift  = mean - actualSignalCount;
        systValues_beforeCuts[p] = stddev;

        std::cout << Form("%-45s (N=%3d, %-10s)  mean=%.2f  shift=%.2f (%+.1f%%)  syst=%.2f (%.1f%%)", genieParams[p].shortName.c_str(), genieParams[p].nUniv, genieParams[p].isMultisim ? "multisim" : "multisigma", mean, shift, (actualSignalCount != 0 ? 100.*shift/actualSignalCount : 0.), stddev, (actualSignalCount != 0 ? 100.*stddev/actualSignalCount : 0.)) << std::endl;
        TH1D *h = new TH1D(("h_genie_" + genieParams[p].shortName).c_str(), (genieParams[p].shortName + ";Total nu+e count;Universes").c_str(), 60, 0, 600);
        
        for(double v : count_genie[p]) h->Fill(v);

        TCanvas *c = new TCanvas(("c_genie_" + genieParams[p].shortName).c_str(), "", 800, 600);
        c->SetLeftMargin(0.12); c->SetBottomMargin(0.12); c->SetRightMargin(0.05); c->SetTopMargin(0.08);
        h->SetLineColor(kBlue+1); h->SetLineWidth(2); h->SetStats(0);
        h->GetXaxis()->SetTitle("Total nu+e Signal Count");
        h->GetYaxis()->SetTitle("Universes");
        h->Draw("HIST E");

        TLine *nomLine = new TLine(actualSignalCount, 0, actualSignalCount, h->GetMaximum()*1.05);
        nomLine->SetLineColor(kMagenta+1); nomLine->SetLineWidth(2); nomLine->Draw("SAME");

        TLatex paramLabel; paramLabel.SetTextSize(0.04); paramLabel.SetNDC();
        paramLabel.DrawLatex(0.15, 0.85, genieParams[p].shortName.c_str());

        c->Update();
        c->SaveAs((base_path + "nuE_signalCount_GENIE_" + genieParams[p].shortName + ".pdf").c_str());

        fOut->cd();
        TDirectory *dir = fOut->GetDirectory("nuESignalCount_GENIE");
        if(!dir) dir = fOut->mkdir("nuESignalCount_GENIE");
        dir->cd();
        h->Write(("nuESignalCount_GENIE_" + genieParams[p].shortName).c_str());
        fOut->cd();

        delete nomLine;
        delete c;
        delete h;
    }

    double totalSystSq_beforeCuts = 0.0;
    for(double s : systValues_beforeCuts) totalSystSq_beforeCuts += s*s;
    double totalSyst_beforeCuts = std::sqrt(totalSystSq_beforeCuts);

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << Form("%-45s  syst=%.2f (%.1f%%)", "TOTAL GENIE (quadrature, 115 params)", totalSyst_beforeCuts, (actualSignalCount != 0 ? 100.*totalSyst_beforeCuts/actualSignalCount : 0.)) << std::endl;
    std::cout << Form("%-45s  %.2f +/- %.2f (syst)", "Signal count", actualSignalCount, totalSyst_beforeCuts) << std::endl;

    double initialSig = nomSig_perCut[0];

    auto buildUnivVecGenie = [&](int p, int cutIdx, std::function<double(double s, double b, double trueSig, double selSig0)> fn) -> std::vector<double> {
        int N = genieParams[p].nUniv;
        std::vector<double> result(N, 0.0);
        for(int u = 0; u < N; u++){
            double s = univSig_perCutParam_genie[p][cutIdx][u];
            double b = univBack_perCutParam_genie[p][cutIdx][u];
            double trueSig = (u < (int)count_genie[p].size()) ? count_genie[p][u] : 0.0;
            double s_beforeCuts = univSig_perCutParam_genie[p][0][u];
            result[u] = fn(s, b, trueSig, s_beforeCuts);
        }
        return result;
    };

    auto printSystBlock_genie = [&](const std::string& blockName, double nomVal, int cutIdx, std::function<double(double s, double b, double ts, double ss)> fn, bool isPct){
        double scale = isPct ? 100.0 : 1.0;
        std::string unitSuffix = isPct ? "%" : "";

        std::cout << "\n--- " << blockName << " ---" << std::endl;
        std::cout << Form("Nominal: %.4f%s", nomVal * scale, unitSuffix.c_str()) << std::endl;

        std::vector<double> systValues(NPARAMS_GENIE, 0.0);
        for(int p = 0; p < NPARAMS_GENIE; p++){
            std::vector<double> vals = buildUnivVecGenie(p, cutIdx, fn);
            double stddev = calcSystGeneric(vals, nomVal, genieParams[p].isMultisim);
            double mean   = calcMeanFromValues(vals);
            double shift  = mean - nomVal;
            systValues[p] = stddev;
            std::cout << Form("%-45s (N=%3d, %-10s)  mean=%.4f%s  shift=%.4f (%+.1f%%)  syst=%.4f%s (%.1f%%)", genieParams[p].shortName.c_str(), genieParams[p].nUniv, genieParams[p].isMultisim ? "multisim" : "multisigma", mean * scale, unitSuffix.c_str(), shift * scale, (nomVal != 0 ? 100.*shift/nomVal : 0.), stddev * scale, unitSuffix.c_str(), (nomVal != 0 ? 100.*stddev/nomVal : 0.)) << std::endl;
        }
        double totalSystSq = 0.0;
        for(double s : systValues) totalSystSq += s*s;
        double totalSyst = std::sqrt(totalSystSq);

        std::cout << "--------------------------------------------" << std::endl;
        std::cout << Form("%-45s  syst=%.4f%s (%.1f%%)", "TOTAL GENIE (quadrature, 115 params)", totalSyst * scale, unitSuffix.c_str(), (nomVal != 0 ? 100.*totalSyst/nomVal : 0.)) << std::endl;
        std::cout << Form("%-45s  %.4f%s +/- %.4f%s (syst)", blockName.c_str(), nomVal*scale, unitSuffix.c_str(), totalSyst*scale, unitSuffix.c_str()) << std::endl;
    };

    auto getQuadratureSystAllGenie = [&](int cutIdx, std::function<double(double s, double b, double ts, double ss)> fn, double nomVal) -> double {
        double totalSq = 0.0;
        for(int p = 0; p < NPARAMS_GENIE; p++){
            std::vector<double> vals = buildUnivVecGenie(p, cutIdx, fn);
            double s = calcSystGeneric(vals, nomVal, genieParams[p].isMultisim);
            totalSq += s*s;
        }
        return std::sqrt(totalSq);
    };

    auto printCutStage_genie = [&](const std::string& label, int cutIdx){
        double nomS = nomSig_perCut[cutIdx];
        double nomB = nomBack_perCut[cutIdx];
        double nomEff = (actualSignalCount > 0) ? nomS / actualSignalCount : 0.0;
        double nomSel = (initialSig > 0) ? nomS / initialSig : 0.0;
        double nomPur = (nomS + nomB > 0) ? nomS / (nomS + nomB) : 0.0;
        double nomER = nomEff * nomPur;
        double nomSR = nomSel * nomPur;

        int width = (int)label.size() + 8;
        std::string bar(width, '=');
        std::cout << "\n" << bar << std::endl;
        std::cout << "    " << label << " (GENIE)" << std::endl;
        std::cout << bar << std::endl;

        printSystBlock_genie("Signal Count", nomS, cutIdx, [](double s, double b, double ts, double ss){ return s; }, false);
        printSystBlock_genie("Background Count", nomB, cutIdx, [](double s, double b, double ts, double ss){ return b; }, false);
        printSystBlock_genie("Efficiency", nomEff, cutIdx, [](double s, double b, double ts, double ss){ return (ts > 0) ? s/ts : 0.0; }, true);
        printSystBlock_genie("Selection Efficiency", nomSel, cutIdx, [](double s, double b, double ts, double ss){ return (ss > 0) ? s/ss : 0.0; }, true);
        printSystBlock_genie("Purity", nomPur, cutIdx, [](double s, double b, double ts, double ss){ double tot=s+b; return (tot>0)? s/tot : 0.0; }, true);
        printSystBlock_genie("Efficiency x Purity", nomER, cutIdx, [](double s, double b, double ts, double ss){ double tot=s+b; double eff=(ts>0)?s/ts:0.0; double pur=(tot>0)?s/tot:0.0; return eff*pur; }, true);
        printSystBlock_genie("SelEff x Purity", nomSR, cutIdx, [](double s, double b, double ts, double ss){ double tot=s+b; double sel=(ss>0)?s/ss:0.0; double pur=(tot>0)?s/tot:0.0; return sel*pur; }, true);
    };

    printCutStage_genie("Before Cuts", 0);
    if(clearCosmicCut == 1) printCutStage_genie("Cut 1: Clear Cosmic", 1);
    if(numPFPs0Cut == 1) printCutStage_genie("Cut 2: Num PFPs != 0", 2);
    if(numRecoNeutrinosCut == 1) printCutStage_genie("Cut 3: Num Reco Neutrinos", 3);
    if(CRUMBSCut == 1) printCutStage_genie("Cut 4: CRUMBS Score", 4);
    if(FVCut == 1) printCutStage_genie("Cut 5: FV Cut", 5);
    if(primaryPFPCut == 1) printCutStage_genie("Cut 6: Primary PFP", 6);
    if(razzledPDG11Cut == 1) printCutStage_genie("Cut 7: Razzled PDG11", 7);
    if(razzledPDG211Cut == 1) printCutStage_genie("Cut 8: Razzled PDG211", 8);
    if(ETheta2Cut == 1) printCutStage_genie("Cut 9: ETheta2", 9);
    if(dEdxCut == 1) printCutStage_genie("Cut 10: dEdx", 10);

    auto plotPerCutUniverseDist_derived_genie = [&](int p, const std::string& quantityName, const std::string& xAxisTitle,
                                                      std::function<double(double s, double b, double trueSig, double selSig0)> fn, int color, int cutIdx){
        std::vector<double> vals = buildUnivVecGenie(p, cutIdx, fn);
        if(vals.empty()) return;

        double lo = *std::min_element(vals.begin(), vals.end());
        double hi = *std::max_element(vals.begin(), vals.end());
        double range = hi - lo;
        lo = std::max(0.0, lo - 0.1*range);
        hi = hi + 0.1*range;
        if(hi <= lo) hi = lo + 1.0;

        std::string histName = Form("perCut_%s_%s_%s", quantityName.c_str(), genieParams[p].shortName.c_str(), cutNames_syst[cutIdx].c_str());
        TH1D *h = new TH1D(histName.c_str(), "", 50, lo, hi);
        for(double v : vals) h->Fill(v);

        TCanvas *c = new TCanvas(("cPerCut_" + histName).c_str(), "", 800, 600);
        c->SetLeftMargin(0.12); c->SetBottomMargin(0.12); c->SetRightMargin(0.05); c->SetTopMargin(0.08);
        h->SetLineColor(color); h->SetLineWidth(2); h->SetStats(0);
        h->GetXaxis()->SetTitle(xAxisTitle.c_str());
        h->GetYaxis()->SetTitle("Universes");
        h->Draw("HIST E");

        double nomS = nomSig_perCut[cutIdx];
        double nomB = nomBack_perCut[cutIdx];
        double nomVal = fn(nomS, nomB, actualSignalCount, initialSig);
        TLine *ln = new TLine(nomVal, 0, nomVal, h->GetMaximum()*1.05);
        ln->SetLineColor(kMagenta+1); ln->SetLineWidth(2); ln->Draw("SAME");

        TLatex lx; lx.SetTextSize(0.04); lx.SetNDC();
        lx.DrawLatex(0.15, 0.85, (genieParams[p].shortName + " - " + cutNames_syst[cutIdx]).c_str());

        c->Update();
        c->SaveAs((base_path + "perCut_" + quantityName + "_" + genieParams[p].shortName + "_" + cutNames_syst[cutIdx] + ".pdf").c_str());

        fOut->cd();
        std::string dirName = "perCut_GENIE_" + quantityName;
        TDirectory *dir = fOut->GetDirectory(dirName.c_str());
        if(!dir) dir = fOut->mkdir(dirName.c_str());
        dir->cd();
        h->Write(histName.c_str());
        fOut->cd();

        delete ln; delete h; delete c;
    };

    auto plotSevenQuantities_genie = [&](int p, int cutIdx){
        plotPerCutUniverseDist_derived_genie(p, "signalSlice", "Number of Signal Slices", [](double s,double b,double ts,double ss){ return s; }, kBlue+1, cutIdx);
        plotPerCutUniverseDist_derived_genie(p, "backgroundSlice", "Number of Background Slices", [](double s,double b,double ts,double ss){ return b; }, kRed+1, cutIdx);
        plotPerCutUniverseDist_derived_genie(p, "efficiency", "Efficiency", [](double s,double b,double ts,double ss){ return (ts>0)?s/ts:0.0; }, kBlue+1, cutIdx);
        plotPerCutUniverseDist_derived_genie(p, "selEfficiency", "Selection Efficiency", [](double s,double b,double ts,double ss){ return (ss>0)?s/ss:0.0; }, kCyan+1, cutIdx);
        plotPerCutUniverseDist_derived_genie(p, "purity", "Purity", [](double s,double b,double ts,double ss){ double tot=s+b; return (tot>0)?s/tot:0.0; }, kGreen+2, cutIdx);
        plotPerCutUniverseDist_derived_genie(p, "effXpurity", "Efficiency #times Purity", [](double s,double b,double ts,double ss){ double tot=s+b; double eff=(ts>0)?s/ts:0.0; double pur=(tot>0)?s/tot:0.0; return eff*pur; }, kOrange+1, cutIdx);
        plotPerCutUniverseDist_derived_genie(p, "selEffXpurity", "Selection Efficiency #times Purity", [](double s,double b,double ts,double ss){ double tot=s+b; double sel=(ss>0)?s/ss:0.0; double pur=(tot>0)?s/tot:0.0; return sel*pur; }, kViolet+1, cutIdx);
    };

    if(makePerCutPlots_GENIE){
        // ON: all 11 cuts -> 115 * 7 * 11 = 8855 plots
        for(int p = 0; p < NPARAMS_GENIE; p++)
            for(int cut = 0; cut < NCUTS; cut++)
                plotSevenQuantities_genie(p, cut);
    } else {
        // OFF: before cuts + after all cuts -> 115 * 7 * 2 = 1610 plots
        int finalCut = NCUTS - 1; // index 10 = "dEdx"
        for(int p = 0; p < NPARAMS_GENIE; p++){
            plotSevenQuantities_genie(p, 0);        // before cuts
            plotSevenQuantities_genie(p, finalCut); // all cuts applied
        }
    }

    std::ofstream out_tablefile(tableFileName, std::ios::app);
    if(out_tablefile.is_open()){

        auto fmtPMpct = [](double val, double syst, int valPrec, int systPrec) -> std::string {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(valPrec) << 100.*val << " $\\pm$ " << std::fixed << std::setprecision(systPrec) << 100.*syst << "\\%";
            return oss.str();
        };
        auto fmtPM = [](double val, double syst, int valPrec, int systPrec) -> std::string {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(valPrec) << val << " $\\pm$ " << std::fixed << std::setprecision(systPrec) << syst;
            return oss.str();
        };

        out_tablefile << "=========== DL Nu+E Vertexing: GENIE Systematics (quadrature over 115 params) ===========" << std::endl;
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\epsilon$ (\\%)} & \\textbf{Selection $\\epsilon$ (\\%)} & \\textbf{$\\rho$ (\\%)} & \\textbf{$\\epsilon\\rho$} & \\textbf{Selection $\\epsilon\\rho$} & Signal Left & Background Left \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;

        auto writeRow_genie = [&](const std::string& cutLabel, int cutIdx, double nomS, double nomB){
            double eff = (actualSignalCount > 0) ? nomS / actualSignalCount : 0.0;
            double selEff = (eventsBeforeCuts_DLNuE.signal > 0) ? nomS / eventsBeforeCuts_DLNuE.signal : 0.0;
            double pur = (nomS + nomB > 0) ? nomS / (nomS + nomB) : 0.0;
            double epsRho = eff * pur;
            double selEpsRho = selEff * pur;

            double sigSyst   = getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ return s; }, nomS);
            double backSyst  = getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ return b; }, nomB);
            double effSyst   = getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ return (ts>0)?s/ts:0.0; }, eff);
            double selEffSyst= getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ return (ss>0)?s/ss:0.0; }, selEff);
            double purSyst   = getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ double tot=s+b; return (tot>0)?s/tot:0.0; }, pur);
            double epsRhoSyst= getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ double tot=s+b; double e=(ts>0)?s/ts:0.0; double p=(tot>0)?s/tot:0.0; return e*p; }, epsRho);
            double selEpsRhoSyst= getQuadratureSystAllGenie(cutIdx, [](double s,double b,double ts,double ss){ double tot=s+b; double se=(ss>0)?s/ss:0.0; double p=(tot>0)?s/tot:0.0; return se*p; }, selEpsRho);

            out_tablefile
                << cutLabel << " & "
                << fmtPMpct(eff, effSyst, 4, 4) << " & "
                << fmtPMpct(selEff, selEffSyst, 4, 4) << " & "
                << fmtPMpct(pur, purSyst, 4, 4) << " & "
                << fmtPM(epsRho, epsRhoSyst, 6, 6) << " & "
                << fmtPM(selEpsRho, selEpsRhoSyst, 6, 6) << " & "
                << std::fixed << std::setprecision(0) << nomS << " $\\pm$ " << std::fixed << std::setprecision(2) << sigSyst
                << std::defaultfloat << std::setprecision(4) << " (" << 100.*eff << "\\%) & "
                << std::fixed << std::setprecision(0) << nomB << " $\\pm$ " << std::fixed << std::setprecision(2) << backSyst
                << std::defaultfloat << std::setprecision(4) << " (" << 100.*nomB/eventsBeforeCuts_DLNuE.background << "\\%)"
                << " \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        };

        writeRow_genie("No Cut", 0, eventsBeforeCuts_DLNuE.signal, eventsBeforeCuts_DLNuE.background);
        if(clearCosmicCut == 1) writeRow_genie("Remove Clear Cosmic PFPs", 1, eventsAfterCuts_DLNuE.clearCosmicsSig, eventsAfterCuts_DLNuE.clearCosmicsBack);
        if(numPFPs0Cut == 1) writeRow_genie("PFPs in Slice != 0", 2, eventsAfterCuts_DLNuE.numPFPs0Sig, eventsAfterCuts_DLNuE.numPFPs0Back);
        if(numRecoNeutrinosCut == 1) writeRow_genie("1 Reco Neutrino in Slice", 3, eventsAfterCuts_DLNuE.numRecoNeut0Sig, eventsAfterCuts_DLNuE.numRecoNeut0Back);
        if(CRUMBSCut == 1){
            std::ostringstream crumbsLabel;
            crumbsLabel << std::defaultfloat << std::setprecision(7) << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high;
            writeRow_genie(crumbsLabel.str(), 4, eventsAfterCuts_DLNuE.crumbsSig, eventsAfterCuts_DLNuE.crumbsBack);
        }
        if(FVCut == 1) writeRow_genie("FV Cut", 5, eventsAfterCuts_DLNuE.FVSig, eventsAfterCuts_DLNuE.FVBack);
        if(primaryPFPCut == 1){
            std::ostringstream primLabel;
            primLabel << std::defaultfloat << std::setprecision(7) << "Primary PFPs in Slice with $\\geq$ 10 Hits = " << primaryPFPCutValue;
            writeRow_genie(primLabel.str(), 6, eventsAfterCuts_DLNuE.primaryPFPSig, eventsAfterCuts_DLNuE.primaryPFPBack);
        }
        if(razzledPDG11Cut == 1){
            std::ostringstream r11Label;
            r11Label << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Electron Score $\\geq$ " << razzled11Low_highestEnergyPFP;
            writeRow_genie(r11Label.str(), 7, eventsAfterCuts_DLNuE.razzled11Sig, eventsAfterCuts_DLNuE.razzled11Back);
        }
        if(razzledPDG211Cut == 1){
            std::ostringstream r211Label;
            r211Label << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ " << razzled211High_highestEnergyPFP;
            writeRow_genie(r211Label.str(), 8, eventsAfterCuts_DLNuE.razzled211Sig, eventsAfterCuts_DLNuE.razzled211Back);
        }
        if(ETheta2Cut == 1){
            std::ostringstream et2Label;
            et2Label << std::defaultfloat << std::setprecision(7) << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP + PFP Spacepoints 10cm)} $\\leq$ " << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$";
            writeRow_genie(et2Label.str(), 9, eventsAfterCuts_DLNuE.ETheta2Sig, eventsAfterCuts_DLNuE.ETheta2Back);
        }
        if(dEdxCut == 1){
            std::ostringstream dedxLabel;
            dedxLabel << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has " << dEdxLow_highestEnergyPFP << " MeV cm$^{-1}$ $\\leq$ dE/dx $\\leq$ " << dEdxHigh_highestEnergyPFP << " MeV cm$^{-1}$";
            writeRow_genie(dedxLabel.str(), 10, eventsAfterCuts_DLNuE.dEdxSig, eventsAfterCuts_DLNuE.dEdxBack);
        }

        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;
        out_tablefile << "" << std::endl;

        // ----------------------------------------------------------------
        // Table 2: interaction-type breakdown
        // ----------------------------------------------------------------
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c|c| }" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\multicolumn{11}{|c|}{\\textbf{Number of Events Left}} \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\boldsymbol{\\nu+e}$} & \\textbf{NCN$\\boldsymbol{\\pi^0}$} & \\textbf{Other NC} & \\textbf{CC$\\boldsymbol{\\nu_\\mu}$} & \\textbf{CC$\\boldsymbol{\\nu_e}$} & \\textbf{Dirt} & \\textbf{$\\boldsymbol{\\nu+e}$ Dirt} & \\textbf{Cosmic} & \\textbf{Other} & \\textbf{$\\boldsymbol{\\nu+e}$ Fuzzy}\\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "No Cut & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.splitInt.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << std::defaultfloat << std::setprecision(4) << "(" << 100*eventsBeforeCuts_DLNuE.splitInt.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.otherNC << " (" << 100*eventsBeforeCuts_DLNuE.splitInt.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.CCnumu << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.CCnue << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.dirt << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuEDirt << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.cosmic << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.other << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        if(clearCosmicCut == 1){
            out_tablefile << "Remove Clear Cosmic PFPs & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(numPFPs0Cut == 1){
            out_tablefile << "PFPs in Slice != 0 & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(numRecoNeutrinosCut == 1){
            out_tablefile << "1 Reco Neutrino in Slice & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(CRUMBSCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(FVCut == 1){
            out_tablefile << "FV Cut & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.FVIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(primaryPFPCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Primary PFPs in Slice with $\\geq$ 10 Hits = " << primaryPFPCutValue << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(razzledPDG11Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Electron Score $\\geq$ " << razzled11Low_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(razzledPDG211Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ " << razzled211High_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(ETheta2Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP + PFP Spacepoints 10cm)} $\\leq$ " << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$ & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        if(dEdxCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has " << dEdxLow_highestEnergyPFP << " MeV cm^{-1} $\\leq$ dE/dx $\\leq$ " << dEdxHigh_highestEnergyPFP << " MeV cm^{-1} & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;
        out_tablefile << "" << std::endl;
        out_tablefile << "\\newpage" << std::endl;
        out_tablefile << "" << std::endl;
    }

    fOut->Write();
    fOut->Close();
    delete fOut;

    for(int f = 0; f < NWEIGHTFILES; ++f){
        if(fNuEWeightsVec[f]) fNuEWeightsVec[f]->Close();
    }
}
