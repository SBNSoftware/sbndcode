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
    double clearCosmicsSig = 0;
    double clearCosmicsBack = 0;
    eventCounter_struct clearCosmicsIntSplit;
    double numPFPs0Sig = 0;
    double numPFPs0Back = 0;
    eventCounter_struct numPFPs0IntSplit;
    double numRecoNeut0Sig = 0;
    double numRecoNeut0Back = 0;
    eventCounter_struct numRecoNeut0IntSplit;
    double FVSig = 0;
    double FVBack = 0;
    eventCounter_struct FVIntSplit;
    double crumbsSig = 0;
    double crumbsBack = 0;
    eventCounter_struct crumbsIntSplit;
    double primaryPFPSig = 0;
    double primaryPFPBack = 0;
    eventCounter_struct primaryPFPIntSplit;
    double razzled2212Sig = 0;
    double razzled2212Back = 0;
    eventCounter_struct razzled2212IntSplit;
    double razzled13Sig = 0;
    double razzled13Back = 0;
    eventCounter_struct razzled13IntSplit;
    double razzled211Sig = 0;
    double razzled211Back = 0;
    eventCounter_struct razzled211IntSplit;
    double razzled22Sig = 0;
    double razzled22Back = 0;
    eventCounter_struct razzled22IntSplit;
    double razzled11Sig = 0;
    double razzled11Back = 0;
    eventCounter_struct razzled11IntSplit;
    double dEdxSig = 0;
    double dEdxBack = 0;
    eventCounter_struct dEdxIntSplit;
    double fracHitsContainedSig = 0;
    double fracHitsContainedBack = 0;
    eventCounter_struct fracHitsContainedIntSplit;
    double numHitsSig = 0;
    double numHitsBack = 0;
    eventCounter_struct numHitsIntSplit;
    double trackscoreSig = 0;
    double trackscoreBack = 0;
    eventCounter_struct trackscoreIntSplit;
    double ETheta2Sig = 0;    
    double ETheta2Back = 0;    
    eventCounter_struct ETheta2IntSplit;    
    double showerLengthSig = 0;    
    double showerLengthBack = 0;    
    eventCounter_struct showerLengthIntSplit;    
    double showerEnergySig = 0;    
    double showerEnergyBack = 0;    
    eventCounter_struct showerEnergyIntSplit;    
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
    std::sort(fileList.begin(), fileList.end()); // deterministic order
    return fileList;
}

void nuESelectionNumbersWithSystematics_macro(){

    std::string cutsApplied = "allCuts";
    std::string base_path = "/nashome/c/coackley/systPlotsNumbers9July_" + cutsApplied + "/";
    std::string tableFileName = base_path + "table.txt";

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

    // Cut values
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
  
    // Creates the output directory if it doesn't already exist 
    if (gSystem->AccessPathName(base_path.c_str())) {
        gSystem->mkdir(base_path.c_str(), kTRUE);
    }

    // Open and clear the table txt file
    std::ofstream clearTableFile(tableFileName, std::ios::trunc);
    if (!clearTableFile.is_open()) {
        std::cerr << "Error: could not open or create " << tableFileName << std::endl;
        return;
    }
    clearTableFile.close();

    std::string inputDir = "/exp/sbnd/data/users/coackley/testFiles/analysed";
    std::vector<std::string> inputFiles = listRootFiles(inputDir);
    std::cout << "Found " << inputFiles.size() << " input files in " << inputDir << std::endl;
    if(inputFiles.empty()){
        std::cerr << "No input files found in " << inputDir << std::endl;
        return;
    }

    TChain *tree = new TChain("ana/NuE");
    TChain *subRunTree = new TChain("ana/SubRun");
    TChain *weightsTree = new TChain("ana/NuEWeights");

    // Add the SAME file list, in the SAME order, to tree and weightsTree
    // so that global entry i in one corresponds to global entry i in the other.
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

    TFile *fOut = new TFile("/exp/sbnd/data/users/coackley/selectionNumberSystematicPlots_9July.root", "RECREATE");
    if(!fOut || fOut->IsZombie()){
        std::cerr << "Error creating output ROOT file" << std::endl;
        return;
    }
    
    // Initialise the counters used to count number of slices before and after cuts
    beforeEventCount_struct eventsBeforeCuts_DLNuE;
    eventCounting_struct eventsAfterCuts_DLNuE;

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
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsNuENuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBUboone;

    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;
    double totalPOTNuENuE = 0;

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
    double targetSpills = (targetPOT/(5e12));

    double BNBScaledSpills_NuE = ((targetPOT/POTBNBNuE_notMissing) * BNBSpillsSumNuE);
    double SignalScaledSpills_NuE = ((targetPOT/POTSignalNuE_notMissing) * SignalSpillsSumNuE);

    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    double totalPOTNuENuE_notMissing = (POTNuENuE_notMissing + POTBNBNuE_notMissing);

    std::cout << "POT from nue sample = " << POTNuENuE_notMissing << ", POT from BNB sample = " << POTBNBNuE_notMissing << ", total nue POT = " << totalPOTNuENuE_notMissing << std::endl;

    // Weights used to scale everything to 1e21 POT
    weights_struct weights;
    weights.signalNuE = targetPOT / POTSignalNuE_notMissing;
    weights.BNBNuE = targetPOT /POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;
    weights.NuENuE = targetPOT / totalPOTNuENuE_notMissing;

    std::cout << "Weights DLNu+E: BNB = " << weights.BNBNuE << ", Signal = " << weights.signalNuE << ", Intime Cosmics = " << weights.cosmicsNuE << ", nue = " << weights.NuENuE << std::endl;

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
    double nuEScatterTrueVX_weights, nuEScatterTrueVY_weights, nuEScatterTrueVZ_weights;

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

    weightsTree->SetBranchAddress("nuEScatterTrueVX_weightTree", &nuEScatterTrueVX_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVY_weightTree", &nuEScatterTrueVY_weights);
    weightsTree->SetBranchAddress("nuEScatterTrueVZ_weightTree", &nuEScatterTrueVZ_weights);

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

    const int NUNIV = 1000; // number of universes
    long long nuEWeightFallbacks = 0, nuEWeightCalls = 0;
    long long sliceWeightFallbacks = 0, sliceWeightCalls = 0;

    auto getNuEWeight = [&](std::vector<double>* vec, int u) -> double {
        nuEWeightCalls++;
        if(!vec || (int)vec->size() != NUNIV){ nuEWeightFallbacks++; return 1.0; }
        return vec->at(u);
    };

    auto getSliceWeight = [&](std::vector<std::vector<double>>* vec, size_t sliceIdx, int u, bool wFound) -> double {
        sliceWeightCalls++;
        if(!wFound || !vec || sliceIdx >= vec->size() || (int)vec->at(sliceIdx).size() != NUNIV){ sliceWeightFallbacks++; return 1.0; }
        return vec->at(sliceIdx).at(u);
    };

    // Histograms for each flux parameter (13) + combined (1) = 14 total
    // For total number of true nu+e elastic scattering events (before cuts)
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
    TH1D* h_combined = new TH1D("h_combined", "All Parameters Combined", 60, 0, 600); 

    // Covariance matrix setup: reconstructed recoil angle, post-cuts
    const int NANGLEBINS = 18;
    double angleBinLow  = 0.0;
    double angleBinHigh = 90.0;   // <-- adjust to match your actual angle range/units

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

    // Names of the 13 flux parameters + combined
    std::vector<std::string> paramNames = {"horncurrent", "expskin", "kplus", "kmin", "kzero", "nucleoninex", "nucleonqex", "nucleontotx", "piminus", "pioninex", "pionqex", "piontotx", "piplus", "combined_allParams"};
    int nParams = paramNames.size();
    
    // Names of the 5 slice categories
    std::vector<std::string> catNames = {"cosmic", "signal", "signal_fuzzy", "BNB", "BNB_fuzzy"};
    int nCats   = catNames.size();

    // Creates vectors to store nominal histograms (2D array of histogram pointers): nominal[category][parameter]
    std::vector<std::vector<TH1D*>> nominal(nCats, std::vector<TH1D*>(nParams, nullptr));

    // Creates vectors to store universe histograms (3D array of histogram pointers): univ[category][parameter][universe]
    std::vector<std::vector<std::vector<TH1D*>>> univ(nCats, std::vector<std::vector<TH1D*>>(nParams, std::vector<TH1D*>(NUNIV, nullptr)));

    // Running totals across events for number of true nu+e elastic scattering events
    // One entry per universe, i.e. count_horncurrent[u] = total weighted true nu+e signal count in universe u for the horn current parameter variation
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
    std::vector<double> count_combined(NUNIV, 0.0);

    // The cuts to be applied
    const int NCUTS = 11;
    std::vector<std::string> cutNames_syst = {"beforeCuts", "clearCosmic", "numPFPs0", "numRecoNeut", "crumbs", "FV", "primaryPFP", "razzled11", "razzled211", "ETheta2", "dEdx"};

    // Vectors that store running totals of nominal signal and background slice counts at each cut
    std::vector<double> nomSig_perCut(NCUTS, 0.0);
    std::vector<double> nomBack_perCut(NCUTS, 0.0);

    // parameter numbers: horncurrent=0, expskin=1, kplus=2, kmin=3, kzero=4, nucleoninex=5, nucleonqex=6, nucleontotx=7, piminus=8, pioninex=9, 
    // pionqex=10, piontotx=11, piplus=12, combined=13
    const int NPARAMS_SYST = 14;
    std::vector<std::string> paramNames_syst = {
        "horncurrent","expskin","kplus","kmin","kzero",
        "nucleoninex","nucleonqex","nucleontotx","piminus",
        "pioninex","pionqex","piontotx","piplus","combined"
    };
    
    // Creates a 3D array: univSig_perCutParam[parameter][cut][universe]
    // Used to store the (weighted + universe varied) number of signal slices after each cut for each systematic parameter and universe
    std::vector<std::vector<std::vector<double>>> univSig_perCutParam(NPARAMS_SYST, std::vector<std::vector<double>>(NCUTS, std::vector<double>(NUNIV, 0.0)));
    
    // Creates a 3D array: univBack_perCutParam[parameter][cut][universe]
    // Used to store the (weighted + universe varied) number of background slices after each cut for each systematic parameter and universe
    std::vector<std::vector<std::vector<double>>> univBack_perCutParam(NPARAMS_SYST, std::vector<std::vector<double>>(NCUTS, std::vector<double>(NUNIV, 0.0)));

    // Number of true nu+e elastic scattering events
    double actualSignalCount = 0.0;

    // Helper to calculate the systematic given a vector of values. Calculates it as sqrt((1/N * sum from 1 to N(x_j - nominal))
    auto calcSystFromUniverses = [&](const std::vector<double>& values, double nominal) -> double {
        int N = values.size();
        if(N < 2) return 0.0;
        double sumSq = 0.0;

        for(double x : values)
            sumSq += (x - nominal)*(x - nominal);

        return std::sqrt(sumSq/N);
    };

    // Start looping through the events
    for(Long64_t e = 0; e < numEntries; ++e){
        //std::cout << "============= New Event =============" << std::endl;
        tree->GetEntry(e);
        weightsTree->GetEntry(e);   // index-aligned: same physical event as tree entry e

        int trueSignal = 0;       
 
        // Looking at the true recoil electron in the event (if there is one)
        recoilElectron_struct recoilElectron;
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            //if(truth_recoilElectronPDG->size() > 1) std::cout << "More than 1 true recoil electron in event!" << std::endl;
            if(truth_recoilElectronPDG->at(i) != -999999){
                // There is a true recoil electron in the event
                recoilElectron.energy = truth_recoilElectronEnergy->at(i);
                recoilElectron.angle = truth_recoilElectronAngle->at(i);
                recoilElectron.dx = truth_recoilElectronDX->at(i);
                recoilElectron.dy = truth_recoilElectronDY->at(i);
                recoilElectron.dz = truth_recoilElectronDZ->at(i);
            } else if(truth_recoilElectronPDG->size() == 1 && truth_recoilElectronPDG->at(i) == -999999){
                // There is no recoil electron in the event
                recoilElectron.energy = -999999;
                recoilElectron.angle = -999999;
                recoilElectron.dx = -999999;
                recoilElectron.dy = -999999;
                recoilElectron.dz = -999999;
            }
        }


        if(nuEScatter == 1 && signal == 1 && DLCurrent == 5){
            // This is an event with a nu+e elastic scatter in it (from the signal files)
            if(recoilElectron.energy > 150){
                // nu+e elastic scatter must have true recoil electron with energy > 150 MeV
                if((FVCut == 0 && (((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))) || (FVCut == 1 && (((nuEScatterTrueVX > FVCut_xLow) && (nuEScatterTrueVX < FVCut_xHigh) && (std::abs(nuEScatterTrueVX) > FVCut_xCentre)) && ((nuEScatterTrueVY > FVCut_yLow) && (nuEScatterTrueVY < FVCut_yHigh)) && ((nuEScatterTrueVZ > FVCut_zLow) && (nuEScatterTrueVZ < FVCut_zHigh))))){
                    // True nu+e elastic scattering event within the active volume (if FVCut == 0) or FV (if FVCut == 1))
                    actualSignalCount += weights.signalNuE; // nominal value
                    trueSignal = 1;

                    // Check whether there is a weight associated with the true nu+e elastic scatter
                    // if weightsFound == 0 then no event was found to match between the 2 trees
                    // nuEScatter_MCTruthFlux_weight_horncurrent->size() == NUNIV checks whether there are the same number of entries as universes
                    bool nuEWeightsValid = nuEScatter_MCTruthFlux_weight_horncurrent && (nuEScatter_MCTruthFlux_weight_horncurrent->size() == NUNIV);

                    if(nuEWeightsValid){
                        // Loop through the universes
                        for(int u = 0; u < NUNIV; u++){
                            // combinedWeight is the 13 flux parameter weights for that universe multiplied together to get the combined weight in that universe
                            double combinedWeight = getNuEWeight(nuEScatter_MCTruthFlux_weight_horncurrent, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_expskin, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_kplus, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_kmin, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_kzero, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleoninexsec, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleonqexsec, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleontotxsec, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_piminus, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_pioninexsec, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_pionqexsec, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_piontotxsec, u) * getNuEWeight(nuEScatter_MCTruthFlux_weight_piplus, u);

                            // Add to the count (number of true nu+e elastic scattering events) for that universe the POT weight * universe weight for that flux parameter
                            // After looping through all events count_horncurrent[u] will be the number of true nu+e elastic scatters in universe u under a shift of the horn current parameter, etc, etc
                            count_horncurrent[u] += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_horncurrent, u);
                            count_expskin[u]     += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_expskin, u);
                            count_kplus[u]       += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_kplus, u);
                            count_kmin[u]        += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_kmin, u);
                            count_kzero[u]       += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_kzero, u);
                            count_nucleoninex[u] += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleoninexsec, u);
                            count_nucleonqex[u]  += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleonqexsec, u);
                            count_nucleontotx[u] += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_nucleontotxsec, u);
                            count_piminus[u]     += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_piminus, u);
                            count_pioninex[u]    += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_pioninexsec, u);
                            count_pionqex[u]     += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_pionqexsec, u);
                            count_piontotx[u]    += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_piontotxsec, u);
                            count_piplus[u]      += weights.signalNuE * getNuEWeight(nuEScatter_MCTruthFlux_weight_piplus, u);
                            count_combined[u]    += weights.signalNuE * combinedWeight;
                        }
                    }
                }
            }
        }

        // Looking at the reco slices
        if(reco_sliceID->size() == 0) continue;

        //if(weightsFound) std::cout << "Number of slices = " << reco_sliceID_weights->size() << std::endl;
        //std::cout << "--- Slices for event ---" << std::endl;

        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            // Loop through slices in event
            if(reco_sliceID->at(slice) == -999999) continue;

            //std::cout << "Slice " << slice << ": ID = " << reco_sliceID->at(slice) << ", CRUMBS Score = " << reco_sliceScore->at(slice) << std::endl;

            double sliceRecoVX = -999999;
            double sliceRecoVY = -999999;
            double sliceRecoVZ = -999999;

            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    sliceRecoVX = reco_neutrinoVX->at(recoNeut);
                    sliceRecoVY = reco_neutrinoVY->at(recoNeut);
                    sliceRecoVZ = reco_neutrinoVZ->at(recoNeut);
                }
            }

            // Assigning a category to the slices
            // 0 = cosmic, 1 = signal, 2 = signal fuzzy, 3 = bnb, 4 = bnb fuzzy, 5 = nue event, 6 = nue fuzzy event
            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                sliceCategoryPlottingMacro = 0;
            } else if(reco_sliceOrigin->at(slice) == 1){
                if(reco_sliceCompleteness->at(slice) > 0.5 && recoilElectron.energy > 150){
                    // Slice must have completeness > 0.5 and have nu+e elastic scatter it comes from has true recoil electron energy > 150 MeV
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

            //if(sliceCategoryPlottingMacro == 0) std::cout << "Cosmic Slice" << std::endl;
            //if(sliceCategoryPlottingMacro == 1 && signal == 1) std::cout << "Signal Slice" << std::endl;
            //if(sliceCategoryPlottingMacro == 2 && signal == 1) std::cout << "Signal Fuzzy Slice" << std::endl;
            //if(sliceCategoryPlottingMacro == 3) std::cout << "BNB Slice" << std::endl;
            //if(sliceCategoryPlottingMacro == 4) std::cout << "BNB Fuzzy Slice" << std::endl;
            
            for(size_t trueParticle = 0; trueParticle < truth_particleSliceID->size(); trueParticle++){
                if(truth_particleSliceID->at(trueParticle) == reco_sliceID->at(slice)){
                    if(truth_particleStatusCode->at(trueParticle) == 1){
                        //std::cout << "True particle in slice: PDG = " << truth_particlePDG->at(trueParticle) << std::endl;
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
                        // Slice must have completeness > 0.5 and have nu+e elastic scatter it comes from has true recoil electron energy > 150 MeV
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

            double summedEnergy = 0;
            double numPFPsSlice = 0;
            double numPrimaryPFPsSlice = 0;
            double numPrimaryPFPs10Slice = 0;
            double numHitsInPFPs = 0;
            
            highestEnergyPFP_struct highestEnergyPFP; 

            //std::cout << "------ PFPs before cuts ------" << std::endl;
            // Looping through all PFPs in the slice
            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                    // PFP is in the slice

                    // If the clearCosmicCut == 1 then don't look at PFPs that are classed as clear cosmic
                    if((clearCosmicCut == 1 && reco_particleClearCosmic->at(pfp) == 0) || clearCosmicCut == 0){
                        numPFPsSlice++;
                        if(reco_particleIsPrimary->at(pfp) == 1){
                            numPrimaryPFPsSlice++;
                            if(reco_particleNumHits->at(pfp) >= 10) numPrimaryPFPs10Slice++;
                        }
                        
                        numHitsInPFPs += reco_particleNumHits->at(pfp);

                        //std::cout << "PFP " << pfp << ": Energy = " << reco_particleBestPlaneEnergy->at(pfp) << ", Clear Cosmic = " << reco_particleClearCosmic->at(pfp) << ", True PDG = " << reco_particleTruePDG->at(pfp) << ", Vertex = (" << reco_particleVX->at(pfp) << ", " << reco_particleVY->at(pfp) << ", " << reco_particleVZ->at(pfp) << ")" << std::endl;
                        summedEnergy += reco_particleBestPlaneEnergy->at(pfp);

                        if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergyPFP.energy){
                            // This is the highest energy PFP out of the PFPs looked at so far
                            highestEnergyPFP.energy = reco_particleBestPlaneEnergy->at(pfp);
                            highestEnergyPFP.theta = reco_particleTheta->at(pfp);
                            highestEnergyPFP.PFPID = reco_particleID->at(pfp);
                            highestEnergyPFP.dx = reco_particleDX->at(pfp);
                            highestEnergyPFP.dy = reco_particleDY->at(pfp);
                            highestEnergyPFP.dz = reco_particleDZ->at(pfp);
                            highestEnergyPFP.vx = reco_particleVX->at(pfp);
                            highestEnergyPFP.vy = reco_particleVY->at(pfp);
                            highestEnergyPFP.vz = reco_particleVZ->at(pfp);
                            highestEnergyPFP.completeness = reco_particleCompleteness->at(pfp);
                            highestEnergyPFP.purity = reco_particlePurity->at(pfp);
                            highestEnergyPFP.trackscore = reco_particleTrackScore->at(pfp);
                            highestEnergyPFP.primary = reco_particleIsPrimary->at(pfp);
                            highestEnergyPFP.truePDG = reco_particleTruePDG->at(pfp);
                            highestEnergyPFP.trueOrigin = reco_particleTrueOrigin->at(pfp);
                            highestEnergyPFP.trueInt = reco_particleTrueInteractionType->at(pfp);
                            highestEnergyPFP.bestPlanedEdx = reco_particleBestPlanedEdx->at(pfp);
                            highestEnergyPFP.razzledPDG11 = reco_particleRazzledPDG11->at(pfp);
                            highestEnergyPFP.razzledPDG13 = reco_particleRazzledPDG13->at(pfp);
                            highestEnergyPFP.razzledPDG22 = reco_particleRazzledPDG22->at(pfp);
                            highestEnergyPFP.razzledPDG211 = reco_particleRazzledPDG211->at(pfp);
                            highestEnergyPFP.razzledPDG2212 = reco_particleRazzledPDG2212->at(pfp);
                            highestEnergyPFP.razzledBestPDG = reco_particleRazzledBestPDG->at(pfp);
                            highestEnergyPFP.trueVX = reco_particleTrueVX->at(pfp);
                            highestEnergyPFP.trueVY = reco_particleTrueVY->at(pfp);
                            highestEnergyPFP.trueVZ = reco_particleTrueVZ->at(pfp);
                            highestEnergyPFP.trueEndX = reco_particleTrueEndX->at(pfp);
                            highestEnergyPFP.trueEndY = reco_particleTrueEndY->at(pfp);
                            highestEnergyPFP.trueEndZ = reco_particleTrueEndZ->at(pfp);
                            highestEnergyPFP.numHits = reco_particleNumHits->at(pfp);
                            highestEnergyPFP.clearCosmic = reco_particleClearCosmic->at(pfp);

                            if(highestEnergyPFP.trueVX != -999999 && highestEnergyPFP.trueVY != -999999 && highestEnergyPFP.trueVZ != -999999 && highestEnergyPFP.trueEndX != -999999 && highestEnergyPFP.trueEndY != -999999 && highestEnergyPFP.trueEndZ != -999999){
                                double xCoordDiff_length = (highestEnergyPFP.trueVX - highestEnergyPFP.trueEndX);
                                double yCoordDiff_length = (highestEnergyPFP.trueVY - highestEnergyPFP.trueEndY);
                                double zCoordDiff_length = (highestEnergyPFP.trueVZ - highestEnergyPFP.trueEndZ);
                                highestEnergyPFP.trueLength = std::sqrt((xCoordDiff_length * xCoordDiff_length) + (yCoordDiff_length * yCoordDiff_length) + (zCoordDiff_length * zCoordDiff_length));
                            }
                        }
                    }
                }
            }

            double pfp10cm_PCAAngle = -999999;
            double pfp10cm_PCADX = -999999;
            double pfp10cm_PCADY = -999999;
            double pfp10cm_PCADZ = -999999;

            for(size_t pfpAngle = 0; pfpAngle < angleRecalculationPCAPFP10cm_pfpID->size(); ++pfpAngle){
                if(angleRecalculationPCAPFP10cm_pfpID->at(pfpAngle) == highestEnergyPFP.PFPID){
                    pfp10cm_PCAAngle = angleRecalculationPCAPFP10cm_angle->at(pfpAngle);
                    pfp10cm_PCADX = angleRecalculationPCAPFP10cm_dx->at(pfpAngle);
                    pfp10cm_PCADY = angleRecalculationPCAPFP10cm_dy->at(pfpAngle);
                    pfp10cm_PCADZ = angleRecalculationPCAPFP10cm_dz->at(pfpAngle);
                }
            }

            // Looped through all PFPs in the slice and now have the highest energy PFP out
            double angleDifference = -999999;
            double angleDifferencePCAPFP10cm = -999999;

            if((highestEnergyPFP.dx != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((highestEnergyPFP.dx * recoilElectron.dx) + (highestEnergyPFP.dy * recoilElectron.dy) + (highestEnergyPFP.dz * recoilElectron.dz));
                double aMagnitude = std::sqrt((highestEnergyPFP.dx * highestEnergyPFP.dx) + (highestEnergyPFP.dy * highestEnergyPFP.dy) + (highestEnergyPFP.dz * highestEnergyPFP.dz));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifference = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            if((pfp10cm_PCADX != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((pfp10cm_PCADX * recoilElectron.dx) + (pfp10cm_PCADY * recoilElectron.dy) + (pfp10cm_PCADZ * recoilElectron.dz));
                double aMagnitude = std::sqrt((pfp10cm_PCADX * pfp10cm_PCADX) + (pfp10cm_PCADY * pfp10cm_PCADY) + (pfp10cm_PCADZ * pfp10cm_PCADZ));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifferencePCAPFP10cm = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            double recoVX = -999999;
            double recoVY = -999999;
            double recoVZ = -999999;
            int numRecoNeutrinos = 0;

            // Looking for the reco neutrino in the slice
            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    // Reco neutrino is in the slice
                    numRecoNeutrinos++;
                    recoVX = reco_neutrinoVX->at(recoNeut);
                    recoVY = reco_neutrinoVY->at(recoNeut);
                    recoVZ = reco_neutrinoVZ->at(recoNeut);
                }
            }

            size_t wSliceIdx_cached = 999999;
            bool sliceWeightValid_cached = false;
            // Creates a 2D array which will store the systematic weights associated with the slice: sliceUnivWeights[parameter][universe]
            std::vector<std::vector<double>> sliceUnivWeights(NPARAMS_SYST, std::vector<double>(NUNIV, 1.0));
            if(DLCurrent == 5 && signal != 3){
                for(size_t ws = 0; ws < reco_sliceID_weights->size(); ++ws){
                    // Finds the matching slice in the weights tree by comparing sliceID numbers
                    if(reco_sliceID_weights->at(ws) == reco_sliceID->at(slice)){
                        wSliceIdx_cached = ws; // ID of the slice being looked at
                        sliceWeightValid_cached = true;
                        break;
                    }
                }
            
                for(int u = 0; u < NUNIV; u++){
                    // Gets the systematic weights correspinding to the slice for the parameter in the universe
                    double wHorn    = getSliceWeight(reco_sliceMCTruthFlux_weight_horncurrent,    wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wExp     = getSliceWeight(reco_sliceMCTruthFlux_weight_expskin,        wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wKplus   = getSliceWeight(reco_sliceMCTruthFlux_weight_kplus,          wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wKmin    = getSliceWeight(reco_sliceMCTruthFlux_weight_kmin,           wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wKzero   = getSliceWeight(reco_sliceMCTruthFlux_weight_kzero,          wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wNinex   = getSliceWeight(reco_sliceMCTruthFlux_weight_nucleoninexsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wNqex    = getSliceWeight(reco_sliceMCTruthFlux_weight_nucleonqexsec,  wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wNtotx   = getSliceWeight(reco_sliceMCTruthFlux_weight_nucleontotxsec, wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPiminus = getSliceWeight(reco_sliceMCTruthFlux_weight_piminus,        wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPinex   = getSliceWeight(reco_sliceMCTruthFlux_weight_pioninexsec,    wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPiqex   = getSliceWeight(reco_sliceMCTruthFlux_weight_pionqexsec,     wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPitotx  = getSliceWeight(reco_sliceMCTruthFlux_weight_piontotxsec,    wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wPiplus  = getSliceWeight(reco_sliceMCTruthFlux_weight_piplus,         wSliceIdx_cached, u, sliceWeightValid_cached);
                    double wComb    = wHorn*wExp*wKplus*wKmin*wKzero*wNinex*wNqex*wNtotx*wPiminus*wPinex*wPiqex*wPitotx*wPiplus;
                
                    // Puts all the systematic weights corresponding to the slice for the parameter in the universe into the 2D array
                    sliceUnivWeights[0][u]  = wHorn;
                    sliceUnivWeights[1][u]  = wExp;
                    sliceUnivWeights[2][u]  = wKplus;
                    sliceUnivWeights[3][u]  = wKmin;
                    sliceUnivWeights[4][u]  = wKzero;
                    sliceUnivWeights[5][u]  = wNinex;
                    sliceUnivWeights[6][u]  = wNqex;
                    sliceUnivWeights[7][u]  = wNtotx;
                    sliceUnivWeights[8][u]  = wPiminus;
                    sliceUnivWeights[9][u]  = wPinex;
                    sliceUnivWeights[10][u] = wPiqex;
                    sliceUnivWeights[11][u] = wPitotx;
                    sliceUnivWeights[12][u] = wPiplus;
                    sliceUnivWeights[13][u] = wComb;
                }
            } else {
                // Cosmics or no weights: all universe weights stay 1.0 (default value set above)
            }

            // Lambda function which takes 1 input (the cut index)
            // Used to fill the nominal and universe, parameter counts for signal and background slices
            auto fillSliceSystCounters = [&](int cutIdx){
                // Cat == 0 && signal == 4, cosmic slices from nue sample - skip explicitly (instead of relying on weight = 0)
                if(sliceCategoryPlottingMacro == 0 && signal == 4) return;

                // Checks whether the slice is a signal slice
                bool isSigSlice = (sliceCategoryPlottingMacro == 1 && signal == 1);
                // If it is a signal slice then add the POT weight to the nominal signal slice counter, if it isn't signal slice then add 0
                nomSig_perCut[cutIdx] += isSigSlice ? weight : 0.0;

                // If it isn't a signal slice then add the POT weight to the nominal background slice counter, if it is a signal slice then add 0
                nomBack_perCut[cutIdx] += isSigSlice ? 0.0    : weight;

                // Loop through the parameters
                for(int p = 0; p < NPARAMS_SYST; p++){
                    // Loop through the universes
                    for(int u = 0; u < NUNIV; u++){
                        // Weight = POT weight * universe parameter systematic weight
                        double w = weight * sliceUnivWeights[p][u];
                        // If it is a signal slice then add the (POT weight * systematic weight) to the universe signal slice counter, if it isn't signal slice then add 0
                        if(isSigSlice) univSig_perCutParam[p][cutIdx][u] += w;
                        
                        // If it isn't a signal slice then add the (POT weight * systematic weight) to the universe background slice counter, if it is signal slice then add 0
                        else univBack_perCutParam[p][cutIdx][u] += w;
                    }
                }
            };

            // Fill the counters of event categories before any cuts are applied, these are nominal counts
            // eventsBeforeCuts_DLNuE.background, eventsBeforeCuts_DLNuE.signal and the other counters after each cut (eventsAfterCuts_DLNuE.cutSig/Back) should be the same as the nomSig_perCut[cutIdx] and nomBack_perCut[cutIdx] counters 
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 1 && signal == 1){
                    eventsBeforeCuts_DLNuE.signal += weight;
                } else if(sliceCategoryPlottingMacro == 2 && signal == 1){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 3){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 4){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 5){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 6){
                    eventsBeforeCuts_DLNuE.background += weight;
                }

                if(sliceInteractionType == 0 && signal != 4){
                    eventsBeforeCuts_DLNuE.splitInt.cosmic += weight;   
                } else if(sliceInteractionType == 1 && signal == 1){
                    eventsBeforeCuts_DLNuE.splitInt.nuE += weight;
                } else if(sliceInteractionType == 2){
                    eventsBeforeCuts_DLNuE.splitInt.NCNPi0 += weight;
                } else if(sliceInteractionType == 3){
                    eventsBeforeCuts_DLNuE.splitInt.otherNC += weight;
                } else if(sliceInteractionType == 4){
                    eventsBeforeCuts_DLNuE.splitInt.CCnumu += weight;
                } else if(sliceInteractionType == 5){
                    eventsBeforeCuts_DLNuE.splitInt.CCnue += weight;
                } else if(sliceInteractionType == 6){
                    eventsBeforeCuts_DLNuE.splitInt.dirt += weight;
                } else if(sliceInteractionType == 7 && signal == 1){
                    eventsBeforeCuts_DLNuE.splitInt.nuEDirt += weight;
                } else if(sliceInteractionType == 8){
                    eventsBeforeCuts_DLNuE.splitInt.other += weight;
                } else if(sliceInteractionType == 9 && signal == 1){
                    eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy += weight;
                }

                // Fill slice category signal/background counters before cuts (cut index 0 = before cuts)
                fillSliceSystCounters(0);

            }

            // Clear cosmic cut has been applied, add to counters
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
            
                // Fill slice category signal/background counters after clear cosmic cut (cut index 1 = clear cosmic cut)
                fillSliceSystCounters(1);
            }

            // Applying cuts here
            if(numPFPs0Cut == 1 && numPFPsSlice == 0){
                // This is a slice with 0 PFPs in it
                continue;
            }

            // Number of PFPs 0 cut has been applied, add to counters
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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs cuts (cut index 2 = clear cosmic + num PFPs cuts)
                fillSliceSystCounters(2);
            }

            if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0){
                // This is a slice with no reco neutrino
                continue;
            }

            // Number of reco neutrinos cut has been applied, add to counters
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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos cuts (cut index 2 = clear cosmic + num PFPs + num reco neut cuts)
                fillSliceSystCounters(3);
            }

            if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high)){
                // This is a slice with a CRUMBS score outside cut values
                continue;
            }

            // CRUMBS score cut has been applied, add to counters
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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs cuts (cut index 2 = clear cosmic + num PFPs + num reco neut + crumbs cuts)
                fillSliceSystCounters(4);
            }
            
            if(FVCut == 1){
                if(!(recoVX < FVCut_xHigh && recoVX > FVCut_xLow  && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow && recoVZ < FVCut_zHigh)){
                    // Doesn't pass the FV cut values
                    continue;
                }
            }

            // FV cut applied, fill counters
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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs + FV cuts (cut index 2 = clear cosmic + num PFPs + num reco neut + crumbs + FV cuts)
                fillSliceSystCounters(5);
            }

            if(primaryPFPCut == 1 && numPrimaryPFPs10Slice != primaryPFPCutValue){
                // Slice has more than 1 primary PFP in it
                continue;
            }

            // Primary PFP cut has been applied, fill counters
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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs + FV + primary PFP cuts (cut index 2 = clear cosmic + num PFPs + num reco neut + crumbs + FV + primary PFP cuts)
                fillSliceSystCounters(6);
            }

            if(razzledPDG11Cut == 1 && ((highestEnergyPFP.razzledPDG11 > razzled11High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG11 < razzled11Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 11 cut
                continue;
            }

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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs + FV + primary PFP + ETheta2 + razzled 11 cuts (cut index 8 = clear cosmic + num PFPs + num reco neut + crumbs + FV + primary PFP + ETheta2 + razzled11 cuts)
                fillSliceSystCounters(7);
            }

            if(razzledPDG211Cut == 1 && ((highestEnergyPFP.razzledPDG211 > razzled211High_highestEnergyPFP) || (highestEnergyPFP.razzledPDG211 < razzled211Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 211 cut
                continue;
            }

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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs + FV + primary PFP + ETheta2 + razzled 11 + razzled 211 cuts (cut index 9 = clear cosmic + num PFPs + num reco neut + crumbs + FV + primary PFP + ETheta2 + razzled11 + razzled211 cuts)
                fillSliceSystCounters(8);
            }

            if(ETheta2Cut == 1 && ((highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) > ETheta2High_highestEnergyPFP || (highestEnergyPFP.energy * pfp10cm_PCAAngle * pfp10cm_PCAAngle) < ETheta2Low_highestEnergyPFP)){
                // Highest energy PFP in slice doesn't pass the ETheta2 cut
                continue;
            }

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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs + FV + primary PFP + ETheta2 cuts (cut index 7 = clear cosmic + num PFPs + num reco neut + crumbs + FV + primary PFP + ETheta2 cuts)
                fillSliceSystCounters(9);
            }

            if(dEdxCut == 1 && (highestEnergyPFP.bestPlanedEdx > dEdxHigh_highestEnergyPFP || highestEnergyPFP.bestPlanedEdx < dEdxLow_highestEnergyPFP)){
                // Highest energy PFP in slice doesn't pass the dE/dx cut
                continue;
            }

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
            
                // Fill slice category signal/background counters after clear cosmic + num PFPs + num reco neutrinos + crumbs + FV + primary PFP + ETheta2 + razzled 11 + razzled 211 + dE/dx cuts (cut index 10 = clear cosmic + num PFPs + num reco neut + crumbs + FV + primary PFP + ETheta2 + razzled11 + razzled211 + dE/dx cuts)
                fillSliceSystCounters(10);
            }

            // This slice passes all of the cuts applied
            
            // Fill reconstructed-angle histograms (nominal + all 1000 universes)
            if(DLCurrent == 5 && pfp10cm_PCAAngle != -999999){
                h_angle_CV->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);
                if(sliceCategoryPlottingMacro == 1 && signal == 1) h_angle_signal_CV->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);
                else if((sliceCategoryPlottingMacro == 2 && signal == 1) || (sliceCategoryPlottingMacro == 0 && signal != 4) || sliceCategoryPlottingMacro == 3 || sliceCategoryPlottingMacro == 4 || sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6) h_angle_back_CV->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight);
                
                for(int u = 0; u < NUNIV; u++){
                    double wComb = sliceUnivWeights[13][u]; // index 13 = combined flux weight
                    h_angle_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight * wComb);
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) h_angle_signal_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight * wComb);
                    else if((sliceCategoryPlottingMacro == 2 && signal == 1) || (sliceCategoryPlottingMacro == 0 && signal != 4) || sliceCategoryPlottingMacro == 3 || sliceCategoryPlottingMacro == 4 || sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6) h_angle_back_univ[u]->Fill(pfp10cm_PCAAngle*TMath::RadToDeg(), weight * wComb);
                }
            }

        }

        //std::cout << "-------------------------------------------" << std::endl;
    }

    // Fill the universe distribution histograms, there is 1000 entries, 1 for each universe
    // These histograms show the distribution of the total true nu+e elastic scattering events across 1000 universes
    //
    // For example: h_horncurrent is the distribution of the number of true nu+e elastic scattering events in universes under the horncurrent parameter variation
    // count_horncurrent[u] is the number of true nu+e elastic scattering events in universe u under a horncurrent parameter shift 
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
        h_combined->Fill(count_combined[u]);
    }

    // Systematic uncertainty per parameter
    std::cout << "\n=== Systematic Uncertainties on nu+e Signal Count ===" << std::endl;
    std::cout << Form("Nominal: %.2f", actualSignalCount) << std::endl;

    std::vector<double> systValues;

    // Lambda function to calculate the mean, std deviation and shift from a TH1D
    auto computeSyst = [&](const std::string& name, TH1D* h, const std::vector<double>& universes, double nominal){
        // Mean is calculated as 1/N*sum(x)
        double mean = h->GetMean();

        // Std Dev is calculated as sqrt(1/N * sum((x - mean)^2)) - might need to change this
        //double stddev = h->GetStdDev();

        double stddev = calcSystFromUniverses(universes, nominal);

        double shift = mean - nominal;
        std::cout << Form("%-20s  mean=%.2f  shift=%.2f (%+.1f%%)  syst=%.2f (%.1f%%)", name.c_str(), mean, shift, 100.*shift/nominal, stddev, 100.*stddev/nominal) << std::endl;
        systValues.push_back(stddev);
    };

    // Compute the mean, std deviation and shift for each parameter + combined
    computeSyst("horncurrent", h_horncurrent, count_horncurrent, actualSignalCount);
    computeSyst("expskin", h_expskin, count_expskin, actualSignalCount);
    computeSyst("kplus", h_kplus, count_kplus, actualSignalCount);
    computeSyst("kmin", h_kmin, count_kmin, actualSignalCount);
    computeSyst("kzero", h_kzero, count_kzero, actualSignalCount);
    computeSyst("nucleoninex", h_nucleoninex, count_nucleoninex, actualSignalCount);
    computeSyst("nucleonqex", h_nucleonqex, count_nucleonqex, actualSignalCount);
    computeSyst("nucleontotx", h_nucleontotx, count_nucleontotx, actualSignalCount);
    computeSyst("piminus", h_piminus, count_piminus, actualSignalCount);
    computeSyst("pioninex", h_pioninex, count_pioninex, actualSignalCount);
    computeSyst("pionqex", h_pionqex, count_pionqex, actualSignalCount);
    computeSyst("piontotx", h_piontotx, count_piontotx, actualSignalCount);
    computeSyst("piplus", h_piplus, count_piplus, actualSignalCount);

    // Combines the systematics in quadrature (assumes the parameters are uncorrelated)
    double totalSystSq = 0.0;
    for(double s : systValues) totalSystSq += s * s;
    double totalSyst = sqrt(totalSystSq);

    std::cout << "--------------------------------------------" << std::endl;
    std::cout << Form("%-20s  syst=%.2f (%.1f%%)", "TOTAL (quadrature)", totalSyst, 100.*totalSyst/actualSignalCount) << std::endl;
    std::cout << Form("%-20s  %.2f +/- %.2f (syst)", "Signal count", actualSignalCount, totalSyst) << std::endl;

    //double combinedSyst = h_combined->GetStdDev();
    double combinedSyst = calcSystFromUniverses(count_combined, actualSignalCount);
    double combinedMean = h_combined->GetMean();
    std::cout << Form("%-20s  mean=%.2f  shift=%.2f (%+.1f%%)  syst=%.2f (%.1f%%)", "COMBINED (product)", combinedMean, combinedMean - actualSignalCount, 100.*(combinedMean - actualSignalCount)/actualSignalCount, combinedSyst, 100.*combinedSyst/actualSignalCount) << std::endl;

    // For each parameter, the histogram of the universe varied total true nu+e elastic scattering event count is drawn with a vertical line at the nominal value
    // Takes the parameter name, the histogram and the nominal value as inputs
    auto plotUniverseDist = [&](const std::string& paramName, TH1D* h, double nominal){

        TCanvas *c = new TCanvas(("c_" + paramName).c_str(), "", 800, 600);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.08);

        // universe varied histograms being drawn
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
        h->SetStats(0);

        h->Draw("HIST E");

        // Nominal line being drawn
        TLine *nomLine = new TLine(nominal, 0, nominal, h->GetMaximum() * 1.05);
        nomLine->SetLineColor(kMagenta+1);
        nomLine->SetLineWidth(2);
        nomLine->Draw("SAME");

        TLatex latex;
        latex.SetTextColor(kMagenta+1);
        latex.SetTextSize(0.04);
        latex.SetTextFont(62);
        double labelX = nominal - (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) * 0.25;
        double labelY = h->GetMaximum() * 0.65;
        latex.DrawLatex(labelX, labelY, Form("Nominal: %.2f", nominal));

        TLatex potLabel;
        potLabel.SetTextColor(kGray+1);
        potLabel.SetTextSize(0.035);
        potLabel.SetNDC();
        potLabel.DrawLatex(0.70, 0.93, "1#times10^{21} POT");

        TLatex paramLabel;
        paramLabel.SetTextSize(0.04);
        paramLabel.SetNDC();
        paramLabel.DrawLatex(0.15, 0.85, paramName.c_str());

        c->Update();

        std::string outPath = base_path + "nuE_signalCount_" + paramName + ".pdf";
        c->SaveAs(outPath.c_str());

        fOut->cd();
        TDirectory *dir = fOut->GetDirectory("nuESignalCount");
        if(!dir) dir = fOut->mkdir("nuESignalCount");
        dir->cd();
        h->Write(("nuESignalCount_" + paramName).c_str());
        fOut->cd();

        delete nomLine;
        delete c;
    };

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
    plotUniverseDist("combined_allParams", h_combined, actualSignalCount);


    // The number of signal and background slices (nominal) before any cuts have been applied
    double initialSig  = nomSig_perCut[0];
    double initialBack = nomBack_perCut[0];

    // Helper to compute mean and standard deviation of the inputted vector
    auto getMeanStd = [&](const std::vector<double>& v, double nominal) -> std::pair<double,double> {
        double mean = 0; for(double x : v) mean += x; mean /= v.size();
        double var = 0; for(double x : v) var  += (x-nominal)*(x-nominal); var /= (v.size());
        // returns the mean = 1/N * sum(x) and variance = sqrt(1/N * sum((x - nominal)^2)
        return {mean, std::sqrt(var)};
    };

    // Helper to print a block showing the mean, shift, and standard deviation/systematic for all of the 13 parameters + combined
    // Takes the following inputs: blockName = "Signal Count"/"Efficiency"/etc, it is the name of the block
    //                             nomVal = the nominal value
    //                             univVecs = a 2D array, univVec[parameter][universe]
    //                             isPct = flag for whether to print it as a percentage
    auto printSystBlock = [&](const std::string& blockName, double nomVal, const std::vector<std::vector<double>>& univVecs, bool isPct){
        double scale = isPct ? 100.0 : 1.0;
        std::string unitSuffix = isPct ? "%" : "";

        std::cout << "\n--- " << blockName << " ---" << std::endl;
        std::cout << Form("Nominal: %.4f%s", nomVal * scale, unitSuffix.c_str()) << std::endl;

        std::vector<double> systValues;
        std::vector<std::string> pNames = {"horncurrent", "expskin", "kplus", "kmin", "kzero", "nucleoninex", "nucleonqex", "nucleontotx", "piminus", "pioninex", "pionqex", "piontotx", "piplus"};

        // Loop through the 13 flux parameters
        for(int p = 0; p < 13; p++){
            // Use the getMeanStd helper function to get the mean and std deviation
            auto [mean, stddev] = getMeanStd(univVecs[p], nomVal);
            double shift = mean - nomVal;
            std::cout << Form("%-20s  mean=%.4f%s  shift=%.4f (%+.1f%%)  syst=%.4f (%.1f%%)", pNames[p].c_str(), (mean * scale), unitSuffix.c_str(), (shift * scale), (nomVal != 0 ? 100. * shift/nomVal : 0.), (stddev * scale), (nomVal != 0 ? 100. * stddev/nomVal : 0.)) << std::endl;
            systValues.push_back(stddev);
        }

        // Adds the systematic in quadrature to get the total systematic 
        double totalSystSq = 0.0;
        for(double s : systValues) totalSystSq += s * s;
        double totalSyst = std::sqrt(totalSystSq);

        // Gets the mean and standard devaition of the combined parameter shifts
        auto [combMean, combStd] = getMeanStd(univVecs[13], nomVal);
        double combShift = combMean - nomVal;

        std::cout << "--------------------------------------------" << std::endl;
        std::cout << Form("%-20s  syst=%.4f%s (%.1f%%)", "TOTAL (quadrature)", totalSyst * scale, unitSuffix.c_str(), (nomVal != 0 ? 100.*totalSyst/nomVal : 0.)) << std::endl;
        std::cout << Form("%-20s  %.4f%s +/- %.4f%s (syst)", blockName.c_str(), (nomVal * scale), unitSuffix.c_str(), (totalSyst * scale), unitSuffix.c_str()) << std::endl;
        std::cout << Form("%-20s  mean=%.4f%s  shift=%.4f (%+.1f%%)  syst=%.4f (%.1f%%)", "COMBINED (product)", (combMean * scale), unitSuffix.c_str(), (combShift * scale), (nomVal != 0 ? 100. * (combShift/nomVal) : 0.), (combStd * scale), (nomVal != 0 ? 100. * combStd/nomVal   : 0.)) << std::endl;
    };

    // Helper function used to build the per parameter, per universe values of derived quantities such as efficiency and purity
    // Inputs are cutIdx = index of the cut applied, fn = an std::function which takes 4 doubles and returns 1 double
    // Function returns a 2D array with [parameter][universe]
    // Input s = universe varied signal slice count,    b = universe varied background slice count
    //       trueSig = universe varied true nu+e count  selSig0 = nominal 
    auto buildUnivVecs = [&](int cutIdx, std::function<double(double s, double b, double trueSig, double selSig0)> fn) -> std::vector<std::vector<double>>{
        std::vector<std::vector<double>> result(14, std::vector<double>(NUNIV, 0.0));

        // Create true signal vector where each element is a pointer which points to a std::vector<double> containing the number of true nu+e elastic scattering events (1000 entries, 1 for each universe)
        const std::vector<double>* trueSigVecs[14] = {&count_horncurrent, &count_expskin, &count_kplus, &count_kmin, &count_kzero, &count_nucleoninex, &count_nucleonqex, &count_nucleontotx, &count_piminus, &count_pioninex, &count_pionqex, &count_piontotx, &count_piplus, &count_combined};

        // Loop through the 13 + 1 parameters
        for(int p = 0; p < 14; p++){
            // Loop through the 1000 universes
            for(int u = 0; u < NUNIV; u++){
                double s = univSig_perCutParam[p][cutIdx][u]; // number of signal slices passing the cuts for parameter p and universe u
                double b = univBack_perCutParam[p][cutIdx][u]; // number of background slices passing the cuts for parameter p and universe u
                double trueSig = (*trueSigVecs[p])[u]; // number of true nu+e elastic scattering events for parameter p and universe u before cuts
                double s_beforeCuts = univSig_perCutParam[p][0][u]; // number of signal slices for parameter p and universe u before cuts
                result[p][u] = fn(s, b, trueSig, s_beforeCuts); 
            }
        }
        return result;
    };

    // Print one full cut stage
    // Helper function used to compute the nominal values of efficiency, selection efficiency, purity, efficiency x purity and selection efficiency x purity
    // Also generates the universe vectors using buildUnivVecs
    auto printCutStage = [&](const std::string& label, int cutIdx){

        double nomS = nomSig_perCut[cutIdx];
        double nomB = nomBack_perCut[cutIdx];
        double nomEff = (actualSignalCount > 0) ? nomS / actualSignalCount : 0.0;
        double nomSel = (initialSig > 0) ? nomS / initialSig : 0.0;
        double nomPur = (nomS + nomB > 0) ? nomS / (nomS + nomB) : 0.0;
        double nomER = nomEff * nomPur;
        double nomSR = nomSel * nomPur;

        // Prints header
        int width = (int)label.size() + 8;
        std::string bar(width, '=');
        std::cout << "\n" << bar << std::endl;
        std::cout << "    " << label << std::endl;
        std::cout << bar << std::endl;

        // Prints the nominal number of signal slices passing the cut and the mean, shift and syst of the universe varied number of signal slices passing the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the number of signal slices passing the cut with [parameter][universe]
        printSystBlock("Signal Count", nomS, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                    return s;
        }), false);

        // Prints the nominal number of background slices passing the cut and the mean, shift and syst of the universe varied number of background slices passing the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the number of background slices passing the cut with [parameter][universe]
        printSystBlock("Background Count", nomB, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                    return b;
        }), false);

        // Prints the nominal efficiency after applying the cut and the mean, shift and syst of the universe varied efficiency after applying the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the efficiency after applying the cut with [parameter][universe]
        // selection efficiency = s/ts (number of signal slices passing cuts/number of true nu+e elastic scatter events before cuts)
        printSystBlock("Efficiency", nomEff, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                    return (ts > 0) ? s / ts : 0.0;
        }), true);

        // Prints the nominal selection efficiency after applying the cut and the mean, shift and syst of the universe varied selection efficiency after applying the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the selection efficiency after applying the cut with [parameter][universe]
        // selection efficiency = s/ss (number of signal slices passing cuts/number of signal slices before cuts)
        printSystBlock("Selection Efficiency", nomSel, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                    return (ss > 0) ? s / ss : 0.0;
        }), true);

        // Prints the nominal purity after applying the cut and the mean, shift and syst of the universe varied purity after applying the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the purity after applying the cut with [parameter][universe]
        // purity = s/(s+b) (number of signal slices passing cuts/total number of slices passing cuts)
        printSystBlock("Purity", nomPur, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                    double tot = s + b; 
                    return (tot > 0) ? s / tot : 0.0;
        }), true);


        // Prints the nominal eff x purity after applying the cut and the mean, shift and syst of the universe varied eff x purity after applying the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the eff x purity after applying the cut with [parameter][universe]
        // eff x purity = s/ts * s/(s+b) (number of signal slices passing cuts/number of true nu+e elastic scatter events before cuts * number of signal slices passing cuts/total number of slices passing cuts)
        printSystBlock("Efficiency x Purity", nomER, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                double tot = s + b;
                double eff = (ts  > 0) ? s / ts : 0.0;
                double pur = (tot > 0) ? s / tot : 0.0;
                return eff * pur;
        }), true);

        // Prints the nominal selection eff x purity after applying the cut and the mean, shift and syst of the universe varied selection eff x purity after applying the cut (for the 13 + 1 combined parameters)
        // buildUnivVecs returns a 2D array of the selection eff x purity after applying the cut with [parameter][universe]
        // selection eff x purity = s/ss * s/(s+b) (number of signal slices passing cuts/number of signal slices before cuts * number of signal slices passing cuts/total number of slices passing cuts)
        printSystBlock("SelEff x Purity", nomSR, buildUnivVecs(cutIdx, [](double s, double b, double ts, double ss) -> double {
                double tot = s + b;
                double sel = (ss  > 0) ? s / ss : 0.0;
                double pur = (tot > 0) ? s / tot : 0.0;
                return sel * pur;
        }), true);
    };

    // Print all cut stages
    printCutStage("Before Cuts", 0);
    if(clearCosmicCut == 1) printCutStage("Cut 1: Clear Cosmic", 1);
    if(numPFPs0Cut == 1) printCutStage("Cut 2: Num PFPs != 0", 2);
    if(numRecoNeutrinosCut == 1) printCutStage("Cut 3: Num Reco Neutrinos", 3);
    if(CRUMBSCut == 1) printCutStage("Cut 4: CRUMBS Score", 4);
    if(FVCut == 1) printCutStage("Cut 5: FV Cut", 5);
    if(primaryPFPCut == 1) printCutStage("Cut 6: Primary PFP", 6);
    if(razzledPDG11Cut == 1) printCutStage("Cut 7: Razzled PDG11", 7);
    if(razzledPDG211Cut == 1) printCutStage("Cut 8: Razzled PDG211", 8);
    if(ETheta2Cut == 1) printCutStage("Cut 9: ETheta2", 9);
    if(dEdxCut == 1) printCutStage("Cut 10: dEdx", 10);

    auto plotPerCutUniverseDist_derived = [&](int paramIdx, const std::string& paramName, const std::string& quantityName, const std::string& xAxisTitle, std::function<double(double s, double b, double trueSig, double selSig0)> fn, double xLo, double xHi, int color){
        const std::vector<double>* trueSigVecs[14] = {&count_horncurrent, &count_expskin, &count_kplus, &count_kmin, &count_kzero, &count_nucleoninex, &count_nucleonqex, &count_nucleontotx, &count_piminus, &count_pioninex, &count_pionqex, &count_piontotx, &count_piplus, &count_combined};

        for(int cut = 0; cut < NCUTS; cut++){
            TCanvas *c = new TCanvas(Form("cPerCut_%s_%s_c%d", quantityName.c_str(), paramName.c_str(), cut), "", 800, 600);
            c->SetLeftMargin(0.12); 
            c->SetBottomMargin(0.12);
            c->SetRightMargin(0.05); 
            c->SetTopMargin(0.08);

            // Build per-universe values
            std::vector<double> vals(NUNIV, 0.0);
            for(int u = 0; u < NUNIV; u++){
                double s = univSig_perCutParam[paramIdx][cut][u];
                double b = univBack_perCutParam[paramIdx][cut][u];
                double trueSig = (*trueSigVecs[paramIdx])[u];
                double s_beforeCuts = univSig_perCutParam[paramIdx][0][u]; 
                vals[u] = fn(s, b, trueSig, s_beforeCuts);
            }

            double lo = xLo, hi = xHi;
            if(lo == 0.0 && hi == 0.0){
                lo = *std::min_element(vals.begin(), vals.end());
                hi = *std::max_element(vals.begin(), vals.end());
                double range = hi - lo;
                lo = std::max(0.0, lo - 0.1*range);
                hi = hi + 0.1*range;
                if(hi <= lo) hi = lo + 1.0;
            }

            std::string histName = Form("perCut_%s_%s_%s", quantityName.c_str(), paramName.c_str(), cutNames_syst[cut].c_str());

            // Creates and fills the histogram with the 1000 universe values
            TH1D *h = new TH1D(histName.c_str(), "", 50, lo, hi);
            for(double v : vals) h->Fill(v);

            h->SetLineColor(color); h->SetLineWidth(2); h->SetStats(0);
            h->GetXaxis()->SetTitle(xAxisTitle.c_str());
            h->GetYaxis()->SetTitle("Universes");
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetLabelSize(0.04); h->GetYaxis()->SetLabelSize(0.04);
            h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.1);
            h->Draw("HIST E");

            // Nominal value line
            double nomS = nomSig_perCut[cut];
            double nomB = nomBack_perCut[cut];
            double nomVal = fn(nomS, nomB, actualSignalCount, initialSig);
            TLine *ln = new TLine(nomVal, 0, nomVal, h->GetMaximum()*1.05);
            ln->SetLineColor(kMagenta+1); 
            ln->SetLineWidth(2); 
            ln->Draw("SAME");

            TLatex lx;
            lx.SetTextSize(0.04); 
            lx.SetNDC();
            lx.DrawLatex(0.15, 0.85, (paramName + " - " + cutNames_syst[cut]).c_str());

            TLatex nomLabel;
            nomLabel.SetTextColor(kMagenta+1);
            nomLabel.SetTextSize(0.035); nomLabel.SetNDC();
            nomLabel.DrawLatex(0.15, 0.80, Form("Nominal: %.4f", nomVal));

            TLatex potLabel;
            potLabel.SetTextColor(kGray+1); potLabel.SetTextSize(0.035); potLabel.SetNDC();
            potLabel.DrawLatex(0.70, 0.93, "1#times10^{21} POT");

            c->Update();
            c->SaveAs((base_path + "perCut_" + quantityName + "_" + paramName + "_" + cutNames_syst[cut] + ".pdf").c_str());
    
            fOut->cd();
            std::string dirName = "perCut_" + quantityName;
            TDirectory *dir = fOut->GetDirectory(dirName.c_str());
            if(!dir) dir = fOut->mkdir(dirName.c_str());
            dir->cd();
            h->Write(histName.c_str());
            fOut->cd();

            delete ln; 
            delete h; 
            delete c;
        }
    };

    // Call for all 14 parameters and all 5 derived quantities
    for(int p = 0; p < NPARAMS_SYST; p++){
        const std::string& pName = paramNames_syst[p];
        
        plotPerCutUniverseDist_derived(p, pName, "signalSlice", "Number of Signal Slices",[](double s, double b, double ts, double ss) -> double {
                return s;
        }, 0.0, 0.0, kBlue+1);

        plotPerCutUniverseDist_derived(p, pName, "backgroundSlice", "Number of Background Slices",[](double s, double b, double ts, double ss) -> double {
                return b;
        }, 0.0, 0.0, kRed+1);

        // Efficiency = s/ts (number of signal slices after cuts/total number of true nu+e elastic scattering events before cuts)
        plotPerCutUniverseDist_derived(p, pName, "efficiency", "Efficiency",[](double s, double b, double ts, double ss) -> double {
                return (ts > 0) ? s / ts : 0.0;
        }, 0.0, 0.0, kBlue+1);

        // Selection efficiency = s/ss (number of signal slices after cuts/number of signal slices before cuts)
        plotPerCutUniverseDist_derived(p, pName, "selEfficiency", "Selection Efficiency", [](double s, double b, double ts, double ss) -> double {
                return (ss > 0) ? s / ss : 0.0;
        }, 0.0, 0.0, kCyan+1);

        // Purity = s/(s+b) (number of signal slices after cuts/number of slices after cuts)
        plotPerCutUniverseDist_derived(p, pName, "purity", "Purity", [](double s, double b, double ts, double ss) -> double {
                double tot = s + b;
                return (tot > 0) ? s / tot : 0.0;
        }, 0.0, 0.0, kGreen+2);

        // Eff x Pur = s/ts * s/(s+b) (number of signal slices after cuts/total number of true nu+e elastic scattering events before cuts * number of signal slices after cuts/number of slices after cuts)
        plotPerCutUniverseDist_derived(p, pName, "effXpurity", "Efficiency #times Purity", [](double s, double b, double ts, double ss) -> double {
                double tot = s + b;
                double eff = (ts  > 0) ? s / ts  : 0.0;
                double pur = (tot > 0) ? s / tot : 0.0;
                return eff * pur;
        }, 0.0, 0.0, kOrange+1);

        // Selection Eff x Pur = s/ss * s/(s+b) (number of signal slices after cuts/number of signal slices before cuts * number of signal slices after cuts/number of slices after cuts)
        plotPerCutUniverseDist_derived(p, pName, "selEffXpurity", "Selection Efficiency #times Purity", [](double s, double b, double ts, double ss) -> double {
                double tot = s + b;
                double sel = (ss  > 0) ? s / ss  : 0.0;
                double pur = (tot > 0) ? s / tot : 0.0;
                return sel * pur;
        }, 0.0, 0.0, kViolet+1);
    }

    // Helper: given a cut index, returns {sigSyst, backSyst, effSyst, selEffSyst, puritSyst, epsRhoSyst, selEpsRhoSyst}
    // using the combined (paramIdx=13) universe variations
    auto getCombinedSyst = [&](int cutIdx) -> std::array<double,7> {
        const int p = 13; // combined parameter index

        std::vector<double>& svec  = univSig_perCutParam[p][cutIdx];
        std::vector<double>& bvec  = univBack_perCutParam[p][cutIdx];
      
        // Nominal signal and background slices left after cut 
        double s_nom = nomSig_perCut[cutIdx];
        double b_nom = nomBack_perCut[cutIdx];

        std::vector<double> effVec(NUNIV), selEffVec(NUNIV), purVec(NUNIV), epsRhoVec(NUNIV), selEpsRhoVec(NUNIV);
        for(int u = 0; u < NUNIV; u++){
            double s       = svec[u];
            double b       = bvec[u];
            double tot     = s + b;
            double trueSig = count_combined[u];
            double s_beforeCuts = univSig_perCutParam[p][0][u];

            effVec[u]       = (trueSig > 0) ? s / trueSig : 0.0;
            selEffVec[u]    = (s_beforeCuts > 0) ? s / s_beforeCuts : 0.0;
            purVec[u]       = (tot > 0) ? s / tot : 0.0;
            epsRhoVec[u]    = effVec[u] * purVec[u];
            selEpsRhoVec[u] = selEffVec[u] * purVec[u];
        } 
        
        auto getStd = [&](const std::vector<double>& v, double nominal) -> double {
            double mean = 0; for(double x : v) mean += x; mean /= v.size();
            double var  = 0; for(double x : v) var  += (x-nominal)*(x-nominal); var /= (v.size());
            return std::sqrt(var);
        };

        double eff_nom = (s_nom/actualSignalCount);
        double selEff_nom = (s_nom/nomSig_perCut[0]);
        double pur_nom = (s_nom/(s_nom + b_nom));
        double epsRho_nom = (eff_nom * pur_nom);
        double selEpsRho_nom = (selEff_nom * pur_nom);

        return {
            getStd(svec, s_nom),
            getStd(bvec, b_nom),
            getStd(effVec, eff_nom),
            getStd(selEffVec, selEff_nom),
            getStd(purVec, pur_nom),
            getStd(epsRhoVec, epsRho_nom),
            getStd(selEpsRhoVec, selEpsRho_nom)
        };
    };

    std::ofstream out_tablefile(tableFileName, std::ios::app);
    if(out_tablefile.is_open()){

        // Pre-compute systematics for every cut stage
        auto s_before     = getCombinedSyst(0);
        auto s_clearCosmic= getCombinedSyst(1);
        auto s_numPFPs0   = getCombinedSyst(2);
        auto s_numRecoNeut= getCombinedSyst(3);
        auto s_crumbs     = getCombinedSyst(4);
        auto s_FV         = getCombinedSyst(5);
        auto s_primaryPFP = getCombinedSyst(6);
        auto s_razzled11  = getCombinedSyst(7);
        auto s_razzled211 = getCombinedSyst(8);
        auto s_ETheta2    = getCombinedSyst(9);
        auto s_dEdx       = getCombinedSyst(10);

        // Indices into the array returned by getCombinedSyst:
        // 0=sigSyst, 1=backSyst, 2=effSyst, 3=selEffSyst, 4=purSyst, 5=epsRhoSyst, 6=selEpsRhoSyst

        // Convenience macro-lambda: formats "value $\pm$ syst" with chosen precisions
        // valPrec = decimal places for value, systPrec = decimal places for syst
        auto fmtPM = [](double val, double syst, int valPrec, int systPrec) -> std::string {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(valPrec) << val << " $\\pm$ " << std::fixed << std::setprecision(systPrec) << syst;
            return oss.str();
        };

        // Version for percentages: multiplies both by 100
        auto fmtPMpct = [](double val, double syst, int valPrec, int systPrec) -> std::string {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(valPrec) << 100.*val << " $\\pm$ " << std::fixed << std::setprecision(systPrec) << 100.*syst << "\\%";
            return oss.str();
        };

        // ----------------------------------------------------------------
        // Table 1: efficiency / purity / signal / background
        // ----------------------------------------------------------------
        out_tablefile << "=========== DL Nu+E Vertexing ===========" << std::endl;
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\epsilon$ (\\%)} & \\textbf{Selection $\\epsilon$ (\\%)} & \\textbf{$\\rho$ (\\%)} & \\textbf{$\\epsilon\\rho$} & \\textbf{Selection $\\epsilon\\rho$} & Signal Left & Background Left \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;

        // Helper that writes one data row given the counters and the syst array
        // sArr = getCombinedSyst result for this cut
        // nomS, nomB = nominal signal and background for this cut
        // prevNomS   = nominal signal at the previous cut level (for selection efficiency denominator)
        // cutLabel   = latex string for the cut name cell
        auto writeRow = [&](const std::string& cutLabel, double nomS, double nomB, const std::array<double,7>& sArr){
            double eff = (actualSignalCount > 0) ? nomS / actualSignalCount : 0.0;
            double selEff = (eventsBeforeCuts_DLNuE.signal > 0) ? nomS / eventsBeforeCuts_DLNuE.signal : 0.0;
            double pur = (nomS + nomB > 0) ? nomS / (nomS + nomB) : 0.0;
            double epsRho = eff * pur;
            double selEpsRho = selEff * pur;

            out_tablefile
                << cutLabel << " & "
                << fmtPMpct(eff,    sArr[2], 4, 4) << " & "
                << fmtPMpct(selEff, sArr[3], 4, 4) << " & "
                << fmtPMpct(pur,    sArr[4], 4, 4) << " & "
                << fmtPM(epsRho,    sArr[5], 6, 6) << " & "
                << fmtPM(selEpsRho, sArr[6], 6, 6) << " & "
                << std::fixed << std::setprecision(0) << nomS
                << " $\\pm$ " << std::fixed << std::setprecision(2) << sArr[0]
                << std::defaultfloat << std::setprecision(4)
                << " (" << 100.*eff << "\\%) & "
                << std::fixed << std::setprecision(0) << nomB
                << " $\\pm$ " << std::fixed << std::setprecision(2) << sArr[1]
                << std::defaultfloat << std::setprecision(4)
                << " (" << 100.*nomB/eventsBeforeCuts_DLNuE.background << "\\%)"
                << " \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        };

        writeRow("No Cut",
                 eventsBeforeCuts_DLNuE.signal, eventsBeforeCuts_DLNuE.background,
                 s_before);

        if(clearCosmicCut == 1)
            writeRow("Remove Clear Cosmic PFPs",
                     eventsAfterCuts_DLNuE.clearCosmicsSig, eventsAfterCuts_DLNuE.clearCosmicsBack,
                     s_clearCosmic);

        if(numPFPs0Cut == 1)
            writeRow("PFPs in Slice != 0",
                     eventsAfterCuts_DLNuE.numPFPs0Sig, eventsAfterCuts_DLNuE.numPFPs0Back,
                     s_numPFPs0);

        if(numRecoNeutrinosCut == 1)
            writeRow("1 Reco Neutrino in Slice",
                     eventsAfterCuts_DLNuE.numRecoNeut0Sig, eventsAfterCuts_DLNuE.numRecoNeut0Back,
                     s_numRecoNeut);

        if(CRUMBSCut == 1){
            std::ostringstream crumbsLabel;
            crumbsLabel << std::defaultfloat << std::setprecision(7)
                        << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high;
            writeRow(crumbsLabel.str(),
                     eventsAfterCuts_DLNuE.crumbsSig, eventsAfterCuts_DLNuE.crumbsBack,
                     s_crumbs);
        }

        if(FVCut == 1)
            writeRow("FV Cut",
                     eventsAfterCuts_DLNuE.FVSig, eventsAfterCuts_DLNuE.FVBack,
                     s_FV);

        if(primaryPFPCut == 1){
            std::ostringstream primLabel;
            primLabel << std::defaultfloat << std::setprecision(7)
                      << "Primary PFPs in Slice with $\\geq$ 10 Hits = " << primaryPFPCutValue;
            writeRow(primLabel.str(),
                     eventsAfterCuts_DLNuE.primaryPFPSig, eventsAfterCuts_DLNuE.primaryPFPBack,
                     s_primaryPFP);
        }

        if(razzledPDG11Cut == 1){
            std::ostringstream r11Label;
            r11Label << std::defaultfloat << std::setprecision(7)
                     << "Highest Energy PFP in Slice has Electron Score $\\geq$ "
                     << razzled11Low_highestEnergyPFP;
            writeRow(r11Label.str(),
                     eventsAfterCuts_DLNuE.razzled11Sig, eventsAfterCuts_DLNuE.razzled11Back,
                     s_razzled11);
        }

        if(razzledPDG211Cut == 1){
            std::ostringstream r211Label;
            r211Label << std::defaultfloat << std::setprecision(7)
                      << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ "
                      << razzled211High_highestEnergyPFP;
            writeRow(r211Label.str(),
                     eventsAfterCuts_DLNuE.razzled211Sig, eventsAfterCuts_DLNuE.razzled211Back,
                     s_razzled211);
        }

        if(ETheta2Cut == 1){
            std::ostringstream et2Label;
            et2Label << std::defaultfloat << std::setprecision(7)
                     << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP + PFP Spacepoints 10cm)} $\\leq$ "
                     << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$";
            writeRow(et2Label.str(),
                     eventsAfterCuts_DLNuE.ETheta2Sig, eventsAfterCuts_DLNuE.ETheta2Back,
                     s_ETheta2);
        }

        if(dEdxCut == 1){
            std::ostringstream dedxLabel;
            dedxLabel << std::defaultfloat << std::setprecision(7)
                      << "Highest Energy PFP in Slice has "
                      << dEdxLow_highestEnergyPFP
                      << " MeV cm$^{-1}$ $\\leq$ dE/dx $\\leq$ "
                      << dEdxHigh_highestEnergyPFP
                      << " MeV cm$^{-1}$";
            writeRow(dedxLabel.str(),
                     eventsAfterCuts_DLNuE.dEdxSig, eventsAfterCuts_DLNuE.dEdxBack,
                     s_dEdx);
        }

        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;
        out_tablefile << "" << std::endl;
        out_tablefile << "" << std::endl;
        out_tablefile << "" << std::endl;

        // ----------------------------------------------------------------
        // Table 2: interaction-type breakdown  (unchanged from your original)
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

    std::cout << "nomSig_perCut[0] = " << nomSig_perCut[0] << ", eventsBeforeCuts_DLNuE.signal = " << eventsBeforeCuts_DLNuE.signal << std::endl;
    std::cout << "nomSig_perCut[2] = " << nomSig_perCut[2] << ", eventsAfterCuts_DLNuE.numPFPs0Sig = " << eventsAfterCuts_DLNuE.numPFPs0Sig << std::endl;
    std::cout << "nomSig_perCut[4] = " << nomSig_perCut[4] << ", eventsAfterCuts_DLNuE.crumbsSig = " << eventsAfterCuts_DLNuE.crumbsSig << std::endl;
    std::cout << "" << std::endl;
    std::cout << "nomBack_perCut[0] = " << nomBack_perCut[0] << ", eventsBeforeCuts_DLNuE.background = " << eventsBeforeCuts_DLNuE.background << std::endl;
    std::cout << "nomBack_perCut[2] = " << nomBack_perCut[2] << ", eventsAfterCuts_DLNuE.numPFPs0Back = " << eventsAfterCuts_DLNuE.numPFPs0Back << std::endl;
    std::cout << "nomBack_perCut[4] = " << nomBack_perCut[4] << ", eventsAfterCuts_DLNuE.crumbsBack = " << eventsAfterCuts_DLNuE.crumbsBack << std::endl;

    std::cout << "nuE weight fallback: " << nuEWeightFallbacks << "/" << nuEWeightCalls << " (" << (nuEWeightCalls? 100.0*nuEWeightFallbacks/nuEWeightCalls : 0.) << "%)" << std::endl;
    std::cout << "slice weight fallback: " << sliceWeightFallbacks << "/" << sliceWeightCalls << " (" << (sliceWeightCalls? 100.0*sliceWeightFallbacks/sliceWeightCalls : 0.) << "%)" << std::endl;

    // Build flux covariance & correlation matrix: reconstructed angle
    TMatrixDSym covMatrix_angle(NANGLEBINS);
    TMatrixDSym corrMatrix_angle(NANGLEBINS);
    
    TMatrixDSym covMatrix_angle_signal(NANGLEBINS);
    TMatrixDSym corrMatrix_angle_signal(NANGLEBINS);
    
    TMatrixDSym covMatrix_angle_back(NANGLEBINS);
    TMatrixDSym corrMatrix_angle_back(NANGLEBINS);

    std::vector<double> CV_binContent(NANGLEBINS);
    std::vector<double> CV_signal_binContent(NANGLEBINS);
    std::vector<double> CV_back_binContent(NANGLEBINS);
    for(int i = 0; i < NANGLEBINS; i++){
        CV_binContent[i] = h_angle_CV->GetBinContent(i+1);
        CV_signal_binContent[i] = h_angle_signal_CV->GetBinContent(i+1);
        CV_back_binContent[i] = h_angle_back_CV->GetBinContent(i+1);
    }

    for(int i = 0; i < NANGLEBINS; i++){
        for(int j = 0; j < NANGLEBINS; j++){
            double sum = 0.0;
            double sum_signal = 0.0;
            double sum_back = 0.0;
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
    for(int u = 0; u < NUNIV; u++) delete h_angle_univ[u]; // free the 1000 per-universe histograms

    delete h_angle_signal_CV;
    delete h_cov_angle_signal;
    delete h_corr_angle_signal;
    delete c_cov_signal;
    delete c_corr_signal;
    for(int u = 0; u < NUNIV; u++) delete h_angle_signal_univ[u]; // free the 1000 per-universe histograms

    delete h_angle_back_CV;
    delete h_cov_angle_back;
    delete h_corr_angle_back;
    delete c_cov_back;
    delete c_corr_back;
    for(int u = 0; u < NUNIV; u++) delete h_angle_back_univ[u]; // free the 1000 per-universe histograms

    fOut->Write();
    fOut->Close();
    delete fOut;

    for(int f = 0; f < NWEIGHTFILES; ++f){
        if(fNuEWeightsVec[f]) fNuEWeightsVec[f]->Close();
    }

}

