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

struct weights_struct{
    double cosmicsCurrent = 0;
    double cosmicsUboone = 0;
    double cosmicsNuE = 0;
};

void cosmic_macro(){

    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_Intime_DLUbooneNuEBDT_23Feb.root");
    std::string base_path = "/nashome/c/coackley/cosmicPlots/";

    // If the directory already exists, delete everything in it
    // If the directory doesn't exists, create it. 
    if (!gSystem->AccessPathName(base_path.c_str())) {
        gSystem->Exec(Form("rm -rf %s/*", base_path.c_str()));
    }
    gSystem->mkdir(base_path.c_str(), kTRUE);

    if(!file){
        std::cerr << "Error opening the file" << std::endl;
        return;
    }

    TDirectory *dir = (TDirectory*)file->Get("ana");
    if(!dir){
        std::cerr << "Directory 'ana' not found" << std::endl;
        return;
    }

    TTree *subRunTree = (TTree*)dir->Get("SubRun");
    if(!subRunTree){
        std::cerr << "SubRun not found" << std::endl;
        return;
    }
 
    TTree *tree = (TTree*)dir->Get("NuE");
    if(!tree){
        std::cerr << "NuE not found" << std::endl;
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

    double cosmicSpillsSumCurrent = 0;
    double cosmicSpillsSumUboone = 0;
    double cosmicSpillsSumNuE = 0;
    
    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 2) cosmicSpillsSumCurrent += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 0) cosmicSpillsSumUboone += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;
    }

    //double targetPOT = POTSignalBDT_notMissing;
    double targetPOT = 1e21;
    double targetSpills = (targetPOT/(5e12));

    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_BDT = (((1-0.0754) * targetGates)/cosmicSpillsSumCurrent);
    double cosmicsWeights_Uboone = (((1-0.0754) * targetGates)/cosmicSpillsSumUboone);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    weights_struct weights;
    weights.cosmicsCurrent = cosmicsWeights_BDT;
    weights.cosmicsUboone = cosmicsWeights_Uboone;
    weights.cosmicsNuE = cosmicsWeights_NuE;

    printf("Weights:\nCurrent: Cosmics = %f\nUboone: Cosmics = %f\nNu+E: Cosmics = %f\n", weights.cosmicsCurrent, weights.cosmicsUboone, weights.cosmicsNuE);


    // NuETree
    UInt_t eventID, runID, subRunID;
    double DLCurrent, signal;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);
    tree->SetBranchAddress("signal", &signal);

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

    Long64_t numEntries = tree->GetEntries();

    // Make plots here
    TH1D* cosmic_deltaDistance_BDT = new TH1D("cosmic_deltaDistance_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;Distance Between True and Reco Vertex (cm)", 100, 0, 200);
    TH1D* cosmic_deltaDistance_DLUboone = new TH1D("cosmic_deltaDistance_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;Distance Between True and Reco Vertex (cm)", 100, 0, 200);
    TH1D* cosmic_deltaDistance_DLNuE = new TH1D("cosmic_deltaDistance_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;Distance Between True and Reco Vertex (cm)", 100, 0, 200);
    
    TH1D* cosmic_deltaXDistance_BDT = new TH1D("cosmic_deltaXDistance_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;X Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaXDistance_DLUboone = new TH1D("cosmic_deltaXDistance_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;X Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaXDistance_DLNuE = new TH1D("cosmic_deltaXDistance_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;X Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    
    TH1D* cosmic_deltaYDistance_BDT = new TH1D("cosmic_deltaYDistance_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;Y Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaYDistance_DLUboone = new TH1D("cosmic_deltaYDistance_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;Y Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaYDistance_DLNuE = new TH1D("cosmic_deltaYDistance_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;Y Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    
    TH1D* cosmic_deltaZDistance_BDT = new TH1D("cosmic_deltaZDistance_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;Z Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaZDistance_DLUboone = new TH1D("cosmic_deltaZDistance_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;Z Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaZDistance_DLNuE = new TH1D("cosmic_deltaZDistance_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;Z Distance Between True and Reco Vertex (cm)", 200, -200, 200);
    
    TH1D* cosmic_deltaTop_BDT = new TH1D("cosmic_deltaTop_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;Distance Between Plane Intersection and Reco Vertex (cm)", 100, 0, 200);
    TH1D* cosmic_deltaTop_DLUboone = new TH1D("cosmic_deltaTop_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;Distance Between Plane Intersection and Reco Vertex (cm)", 100, 0, 200);
    TH1D* cosmic_deltaTop_DLNuE = new TH1D("cosmic_deltaTop_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;Distance Between Plane Intersection and Reco Vertex (cm)", 100, 0, 200);
    
    TH1D* cosmic_deltaXTop_BDT = new TH1D("cosmic_deltaXTop_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;X Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaXTop_DLUboone = new TH1D("cosmic_deltaXTop_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;X Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaXTop_DLNuE = new TH1D("cosmic_deltaXTop_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;X Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    
    TH1D* cosmic_deltaYTop_BDT = new TH1D("cosmic_deltaYTop_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;Y Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaYTop_DLUboone = new TH1D("cosmic_deltaYTop_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;Y Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaYTop_DLNuE = new TH1D("cosmic_deltaYTop_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;Y Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    
    TH1D* cosmic_deltaZTop_BDT = new TH1D("cosmic_deltaZTop_BDT", "Cosmics Intersecting Top Plane: BDT Vertexing;Z Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaZTop_DLUboone = new TH1D("cosmic_deltaZTop_DLUboone", "Cosmics Intersecting Top Plane: DL Uboone Vertexing;Z Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);
    TH1D* cosmic_deltaZTop_DLNuE = new TH1D("cosmic_deltaZTop_DLNuE", "Cosmics Intersecting Top Plane: DL Nu+E Vertexing;Z Distance Between Plane Intersection and Reco Vertex (cm)", 200, -200, 200);

    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0; double zMax = 509.4;

    double bottomOfTrack = 0; double startOfTrack = 0;
    double bottomOfTrackExit = 0; double bottomOfTrackEnd = 0;

    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);
       
        double weight = 0; 
        if(signal == 3 && DLCurrent == 2) weight = weights.cosmicsCurrent;
        if(signal == 3 && DLCurrent == 0) weight = weights.cosmicsUboone;
        if(signal == 3 && DLCurrent == 5) weight = weights.cosmicsNuE;

        for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
            // Loop through pfps in the event
            //std::cout << "====== NEW PFP ======" << std::endl;
            double pfpTrueOrigin = reco_particleTrueOrigin->at(pfp);
            double pfpVX = reco_particleVX->at(pfp);
            double pfpVY = reco_particleVY->at(pfp);
            double pfpVZ = reco_particleVZ->at(pfp);
            double pfpTrueVX = reco_particleTrueVX->at(pfp);
            double pfpTrueVY = reco_particleTrueVY->at(pfp);
            double pfpTrueVZ = reco_particleTrueVZ->at(pfp);
            double pfpTrueEndX = reco_particleTrueEndX->at(pfp);
            double pfpTrueEndY = reco_particleTrueEndY->at(pfp);
            double pfpTrueEndZ = reco_particleTrueEndZ->at(pfp);
            double pfpIsPrimary = reco_particleIsPrimary->at(pfp);
            double pfpSliceID = reco_particleSliceID->at(pfp);
            double pfpID = reco_particleID->at(pfp);
            double pfpClearCosmic = reco_particleClearCosmic->at(pfp);
            //std::cout << "pfpTrueOrigin = " << pfpTrueOrigin << ", pfpIsPrimary = " << pfpIsPrimary << std::endl; 

            if(pfpIsPrimary == 1){
                // PFP comes from a cosmic in truth and is primary
                double trueDX = pfpTrueEndX - pfpTrueVX; 
                double trueDY = pfpTrueEndY - pfpTrueVY; 
                double trueDZ = pfpTrueEndZ - pfpTrueVZ;
                //std::cout << "Passes pfpTrueOrigin == 2 && pfpIsPrimary == 1" << std::endl;

                if(std::abs(trueDY) > 1e-8){
                    // Track isn't parallel to the top plane
                    //std::cout << "Passes std::abs(trueDY) > 1e-8" << std::endl;
                    
                    double t = (yMax - pfpTrueVY)/trueDY;
                    //std::cout << "t = " << t << std::endl; 

                    if(t > 0.0 && t < 1.0){
                        //std::cout << "Passes t > 0.0 && t < 1.0" << std::endl;
                        double xInt = pfpTrueVX + (t * trueDX);
                        double zInt = pfpTrueVZ + (t * trueDZ);

                        bool insideX = (xInt >= (xMin + 10.0)) && (xInt <= (xMax - 10.0));
                        bool insideZ = (zInt >= (zMin + 10.0)) && (zInt <= (zMax - 10.0));

                        if(insideX && insideZ){
                            // Distance between true and reco vertex
                            //std::cout << "Passes insideX && insideZ" << std::endl;
                            double delta_x = (pfpTrueVX - pfpVX);
                            double delta_y = (pfpTrueVY - pfpVY);
                            double delta_z = (pfpTrueVZ - pfpVZ);
                            double deltaD = std::sqrt((delta_x * delta_x) + (delta_y * delta_y) + (delta_z * delta_z));
                            
                            double deltaTop_x = (xInt - pfpVX);
                            double deltaTop_y = (yMax - pfpVY);
                            double deltaTop_z = (zInt - pfpVZ);
                            double deltaTopD = std::sqrt((deltaTop_x * deltaTop_x) + (deltaTop_y * deltaTop_y) + (deltaTop_z * deltaTop_z));

                            // Add to plot
                            if(DLCurrent == 2){
                                cosmic_deltaDistance_BDT->Fill(deltaD, weight);
                                cosmic_deltaXDistance_BDT->Fill(delta_x, weight);
                                cosmic_deltaYDistance_BDT->Fill(delta_y, weight);
                                cosmic_deltaZDistance_BDT->Fill(delta_z, weight);
                                cosmic_deltaTop_BDT->Fill(deltaTopD, weight);
                                cosmic_deltaXTop_BDT->Fill(deltaTop_x, weight);
                                cosmic_deltaYTop_BDT->Fill(deltaTop_y, weight);
                                cosmic_deltaZTop_BDT->Fill(deltaTop_z, weight);
                            } else if(DLCurrent == 0){
                                cosmic_deltaDistance_DLUboone->Fill(deltaD, weight);
                                cosmic_deltaXDistance_DLUboone->Fill(delta_x, weight);
                                cosmic_deltaYDistance_DLUboone->Fill(delta_y, weight);
                                cosmic_deltaZDistance_DLUboone->Fill(delta_z, weight);
                                cosmic_deltaTop_DLUboone->Fill(deltaTopD, weight);
                                cosmic_deltaXTop_DLUboone->Fill(deltaTop_x, weight);
                                cosmic_deltaYTop_DLUboone->Fill(deltaTop_y, weight);
                                cosmic_deltaZTop_DLUboone->Fill(deltaTop_z, weight);
                            } else if(DLCurrent == 5){
                                cosmic_deltaDistance_DLNuE->Fill(deltaD, weight);
                                cosmic_deltaXDistance_DLNuE->Fill(delta_x, weight);
                                cosmic_deltaYDistance_DLNuE->Fill(delta_y, weight);
                                cosmic_deltaZDistance_DLNuE->Fill(delta_z, weight);
                                cosmic_deltaTop_DLNuE->Fill(deltaTopD, weight);
                                cosmic_deltaXTop_DLNuE->Fill(deltaTop_x, weight);
                                cosmic_deltaYTop_DLNuE->Fill(deltaTop_y, weight);
                                cosmic_deltaZTop_DLNuE->Fill(deltaTop_z, weight);

                                if(deltaTopD < 60 && deltaTopD > 30){
                                    if(pfpClearCosmic == 0) std::cout << "deltaTopD = " << deltaTopD << ", signal = " << signal << ", DLCurrent = " << DLCurrent << ", runID = " << runID << ", subRunID = " << subRunID << ", eventID = " << eventID << ", sliceID = " << pfpSliceID << ", PFPID = " << pfpID << std::endl;
                                }
                            }

                            double endDiffX;
                            double endDiffY;
                            double endDiffZ;

                            double exit = -999999;

                            double tBottom = (yMin - pfpTrueVY) / trueDY;
                            if(tBottom > 0.0 && tBottom < 1.0){
                                //Intersects the bottom plane in truth
                                double xBottom = pfpTrueVX + (tBottom * trueDX);
                                double zBottom = pfpTrueVZ + (tBottom * trueDZ);
                           
                                endDiffX = pfpVX - xBottom;
                                endDiffY = pfpVY - yMin;
                                endDiffZ = pfpVZ - zBottom;

                                exit = 1;

                            } else{
                                // Doesn't intersect the bottom plane in truth - stops in TPC
                                //std::cout << "doesn't interact the bottom plane in truth, true end y coord = " << pfpTrueEndY << ", bottom tpc y = " << yMin << std::endl;
                                endDiffX = pfpVX - pfpTrueEndX;
                                endDiffY = pfpVY - pfpTrueEndY;
                                endDiffZ = pfpVZ - pfpTrueEndZ;

                                exit = 0;
                            }

                            double distanceEndPoint = std::sqrt((endDiffX * endDiffX) + (endDiffY * endDiffY) + (endDiffZ * endDiffZ)); 

                            if(distanceEndPoint > deltaTopD){
                                // Reco vertex is closer to the end point/exit point of cosmic than to the top
                                bottomOfTrack++;

                                if(exit == 1){
                                    bottomOfTrackExit++;
                                } else if(exit == 0){
                                    bottomOfTrackEnd++;
                                }
                            } else{
                                // Reco vertex is close to the top than the end point/exit point of the cosmic
                                startOfTrack++;
                            }
                        }
                    }
                } 
            }
        }

    }

    TCanvas* canvasBDT = new TCanvas("canvasBDT", "Graph Draw Options", 200, 10, 600, 400);
    canvasBDT->SetTickx();
    canvasBDT->SetTicky();
    TH1D* cosmic_deltaDistance_BDT_cf = (TH1D*)cosmic_deltaDistance_BDT->GetCumulative();
    double total_BDT = cosmic_deltaDistance_BDT->Integral();
    cosmic_deltaDistance_BDT_cf->Scale(1.0/total_BDT);
    //cosmic_deltaDistance_BDT->Draw("hist");
    cosmic_deltaDistance_BDT_cf->Draw("hist");
    canvasBDT->SaveAs((base_path + "cosmics_distanceTrueRecoVertex_BDT.pdf").c_str());
    delete canvasBDT;
    
    TCanvas* canvasTopBDT = new TCanvas("canvasTopBDT", "Graph Draw Options", 200, 10, 600, 400);
    canvasTopBDT->SetTickx();
    canvasTopBDT->SetTicky();
    TH1D* cosmic_deltaTop_BDT_cf = (TH1D*)cosmic_deltaTop_BDT->GetCumulative();
    double totalTop_BDT = cosmic_deltaTop_BDT->Integral();
    cosmic_deltaTop_BDT_cf->Scale(1.0/totalTop_BDT);
    //cosmic_deltaDistance_BDT->Draw("hist");
    cosmic_deltaTop_BDT_cf->Draw("hist");
    canvasTopBDT->SaveAs((base_path + "cosmics_distanceTopRecoVertex_BDT.pdf").c_str());
    delete canvasTopBDT;

    TCanvas* canvasDLUboone = new TCanvas("canvasDLUboone", "Graph Draw Options", 200, 10, 600, 400);
    canvasDLUboone->SetTickx();
    canvasDLUboone->SetTicky();
    TH1D* cosmic_deltaDistance_DLUboone_cf = (TH1D*)cosmic_deltaDistance_DLUboone->GetCumulative();
    double total_DLUboone = cosmic_deltaDistance_DLUboone->Integral();
    cosmic_deltaDistance_DLUboone_cf->Scale(1.0/total_DLUboone);
    //cosmic_deltaDistance_DLUboone->Draw("hist");
    cosmic_deltaDistance_DLUboone_cf->Draw("hist");
    canvasDLUboone->SaveAs((base_path + "cosmics_distanceTrueRecoVertex_DLUboone.pdf").c_str());
    delete canvasDLUboone;
    
    TCanvas* canvasTopDLUboone = new TCanvas("canvasTopDLUboone", "Graph Draw Options", 200, 10, 600, 400);
    canvasTopDLUboone->SetTickx();
    canvasTopDLUboone->SetTicky();
    TH1D* cosmic_deltaTop_DLUboone_cf = (TH1D*)cosmic_deltaTop_DLUboone->GetCumulative();
    double totalTop_DLUboone = cosmic_deltaTop_DLUboone->Integral();
    cosmic_deltaTop_DLUboone_cf->Scale(1.0/totalTop_DLUboone);
    //cosmic_deltaDistance_DLUboone->Draw("hist");
    cosmic_deltaTop_DLUboone_cf->Draw("hist");
    canvasTopDLUboone->SaveAs((base_path + "cosmics_distanceTopRecoVertex_DLUboone.pdf").c_str());
    delete canvasTopDLUboone;
    
    TCanvas* canvasDLNuE = new TCanvas("canvasDLNuE", "Graph Draw Options", 200, 10, 600, 400);
    canvasDLNuE->SetTickx();
    canvasDLNuE->SetTicky();
    TH1D* cosmic_deltaDistance_DLNuE_cf = (TH1D*)cosmic_deltaDistance_DLNuE->GetCumulative();
    double total_DLNuE = cosmic_deltaDistance_DLNuE->Integral();
    cosmic_deltaDistance_DLNuE_cf->Scale(1.0/total_DLNuE);
    //cosmic_deltaDistance_DLNuE->Draw("hist");
    cosmic_deltaDistance_DLNuE_cf->Draw("hist");
    canvasDLNuE->SaveAs((base_path + "cosmics_distanceTrueRecoVertex_DLNuE.pdf").c_str());
    delete canvasDLNuE;
    
    TCanvas* canvasTopDLNuE = new TCanvas("canvasTopDLNuE", "Graph Draw Options", 200, 10, 600, 400);
    canvasTopDLNuE->SetTickx();
    canvasTopDLNuE->SetTicky();
    TH1D* cosmic_deltaTop_DLNuE_cf = (TH1D*)cosmic_deltaTop_DLNuE->GetCumulative();
    double totalTop_DLNuE = cosmic_deltaTop_DLNuE->Integral();
    cosmic_deltaTop_DLNuE_cf->Scale(1.0/totalTop_DLNuE);
    //cosmic_deltaDistance_DLNuE->Draw("hist");
    cosmic_deltaTop_DLNuE_cf->Draw("hist");
    canvasTopDLNuE->SaveAs((base_path + "cosmics_distanceTopRecoVertex_DLNuE.pdf").c_str());
    delete canvasTopDLNuE;
    
    TCanvas* canvasAll = new TCanvas("canvasAll", "Graph Draw Options", 200, 10, 600, 400);
    canvasAll->SetTickx();
    canvasAll->SetTicky();

    cosmic_deltaDistance_BDT_cf->SetStats(0);
    cosmic_deltaDistance_DLUboone_cf->SetStats(0);
    cosmic_deltaDistance_DLNuE_cf->SetStats(0);

    cosmic_deltaDistance_BDT_cf->SetLineWidth(2);           cosmic_deltaDistance_BDT_cf->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaDistance_DLUboone_cf->SetLineWidth(2);      cosmic_deltaDistance_DLUboone_cf->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaDistance_DLNuE_cf->SetLineWidth(2);         cosmic_deltaDistance_DLNuE_cf->SetLineColor(TColor::GetColor("#3f90da"));

    cosmic_deltaDistance_BDT->SetLineWidth(2);              cosmic_deltaDistance_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaDistance_DLUboone->SetLineWidth(2);         cosmic_deltaDistance_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaDistance_DLNuE->SetLineWidth(2);            cosmic_deltaDistance_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));

    if(cosmic_deltaDistance_BDT->Integral() != 0)           cosmic_deltaDistance_BDT->Scale(1.0 / cosmic_deltaDistance_BDT->Integral());
    if(cosmic_deltaDistance_DLUboone->Integral() != 0)      cosmic_deltaDistance_DLUboone->Scale(1.0 / cosmic_deltaDistance_DLUboone->Integral());
    if(cosmic_deltaDistance_DLNuE->Integral() != 0)         cosmic_deltaDistance_DLNuE->Scale(1.0 / cosmic_deltaDistance_DLNuE->Integral());

    cosmic_deltaXDistance_BDT->SetLineWidth(2);             cosmic_deltaXDistance_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaXDistance_DLUboone->SetLineWidth(2);        cosmic_deltaXDistance_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaXDistance_DLNuE->SetLineWidth(2);           cosmic_deltaXDistance_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaXDistance_BDT->Integral() != 0)          cosmic_deltaXDistance_BDT->Scale(1.0 / cosmic_deltaXDistance_BDT->Integral());
    if(cosmic_deltaXDistance_DLUboone->Integral() != 0)     cosmic_deltaXDistance_DLUboone->Scale(1.0 / cosmic_deltaXDistance_DLUboone->Integral());
    if(cosmic_deltaXDistance_DLNuE->Integral() != 0)        cosmic_deltaXDistance_DLNuE->Scale(1.0 / cosmic_deltaXDistance_DLNuE->Integral());

    cosmic_deltaYDistance_BDT->SetLineWidth(2);             cosmic_deltaYDistance_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaYDistance_DLUboone->SetLineWidth(2);        cosmic_deltaYDistance_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaYDistance_DLNuE->SetLineWidth(2);           cosmic_deltaYDistance_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaYDistance_BDT->Integral() != 0)          cosmic_deltaYDistance_BDT->Scale(1.0 / cosmic_deltaYDistance_BDT->Integral());
    if(cosmic_deltaYDistance_DLUboone->Integral() != 0)     cosmic_deltaYDistance_DLUboone->Scale(1.0 / cosmic_deltaYDistance_DLUboone->Integral());
    if(cosmic_deltaYDistance_DLNuE->Integral() != 0)        cosmic_deltaYDistance_DLNuE->Scale(1.0 / cosmic_deltaYDistance_DLNuE->Integral());

    cosmic_deltaZDistance_BDT->SetLineWidth(2);             cosmic_deltaZDistance_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaZDistance_DLUboone->SetLineWidth(2);        cosmic_deltaZDistance_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaZDistance_DLNuE->SetLineWidth(2);           cosmic_deltaZDistance_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaZDistance_BDT->Integral() != 0)          cosmic_deltaZDistance_BDT->Scale(1.0 / cosmic_deltaZDistance_BDT->Integral());
    if(cosmic_deltaZDistance_DLUboone->Integral() != 0)     cosmic_deltaZDistance_DLUboone->Scale(1.0 / cosmic_deltaZDistance_DLUboone->Integral());
    if(cosmic_deltaZDistance_DLNuE->Integral() != 0)        cosmic_deltaZDistance_DLNuE->Scale(1.0 / cosmic_deltaZDistance_DLNuE->Integral());

    TLegend* legend = new TLegend(0.69, 0.26, 0.87, 0.137);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(cosmic_deltaDistance_BDT_cf, "BDT", "f");
    legend->AddEntry(cosmic_deltaDistance_DLUboone_cf, "DL Uboone", "f");
    legend->AddEntry(cosmic_deltaDistance_DLNuE_cf, "DL Nu+E", "f");
    legend->SetTextSize(0.01);
    legend->SetMargin(0.12);

    cosmic_deltaDistance_BDT_cf->Draw("hist");
    cosmic_deltaDistance_DLUboone_cf->Draw("histsame");
    cosmic_deltaDistance_DLNuE_cf->Draw("histsame");
    legend->Draw();
    canvasAll->SaveAs((base_path + "cosmics_distanceTrueRecoVertex_cf_all.pdf").c_str());

    canvasAll->Clear();
    canvasAll->Update();
    cosmic_deltaDistance_BDT->Draw("hist");
    cosmic_deltaDistance_DLUboone->Draw("histsame");
    cosmic_deltaDistance_DLNuE->Draw("histsame");
    legend->Draw();
    canvasAll->SaveAs((base_path + "cosmics_distanceTrueRecoVertex_all.pdf").c_str());

    canvasAll->Clear();
    canvasAll->Update();
    cosmic_deltaXDistance_BDT->Draw("hist");
    cosmic_deltaXDistance_DLUboone->Draw("histsame");
    cosmic_deltaXDistance_DLNuE->Draw("histsame");
    legend->Draw();
    canvasAll->SaveAs((base_path + "cosmics_distanceXTrueRecoVertex_all.pdf").c_str());

    canvasAll->Clear();
    canvasAll->Update();
    cosmic_deltaYDistance_BDT->Draw("hist");
    cosmic_deltaYDistance_DLUboone->Draw("histsame");
    cosmic_deltaYDistance_DLNuE->Draw("histsame");
    legend->Draw();
    canvasAll->SaveAs((base_path + "cosmics_distanceYTrueRecoVertex_all.pdf").c_str());

    canvasAll->Clear();
    canvasAll->Update();
    cosmic_deltaZDistance_BDT->Draw("hist");
    cosmic_deltaZDistance_DLUboone->Draw("histsame");
    cosmic_deltaZDistance_DLNuE->Draw("histsame");
    legend->Draw();
    canvasAll->SaveAs((base_path + "cosmics_distanceZTrueRecoVertex_all.pdf").c_str());
    delete canvasAll;

    TCanvas* canvasTopAll = new TCanvas("canvasTopAll", "Graph Draw Options", 200, 10, 600, 400);
    canvasTopAll->SetTickx();
    canvasTopAll->SetTicky();

    cosmic_deltaTop_BDT_cf->SetStats(0);
    cosmic_deltaTop_DLUboone_cf->SetStats(0);
    cosmic_deltaTop_DLNuE_cf->SetStats(0);

    cosmic_deltaTop_BDT_cf->SetLineWidth(2);        cosmic_deltaTop_BDT_cf->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaTop_DLUboone_cf->SetLineWidth(2);   cosmic_deltaTop_DLUboone_cf->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaTop_DLNuE_cf->SetLineWidth(2);      cosmic_deltaTop_DLNuE_cf->SetLineColor(TColor::GetColor("#3f90da"));
    
    cosmic_deltaTop_BDT->SetLineWidth(2);           cosmic_deltaTop_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaTop_DLUboone->SetLineWidth(2);      cosmic_deltaTop_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaTop_DLNuE->SetLineWidth(2);         cosmic_deltaTop_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaTop_BDT->Integral() != 0)        cosmic_deltaTop_BDT->Scale(1.0 / cosmic_deltaTop_BDT->Integral());
    if(cosmic_deltaTop_DLUboone->Integral() != 0)   cosmic_deltaTop_DLUboone->Scale(1.0 / cosmic_deltaTop_DLUboone->Integral());
    if(cosmic_deltaTop_DLNuE->Integral() != 0)      cosmic_deltaTop_DLNuE->Scale(1.0 / cosmic_deltaTop_DLNuE->Integral());
    
    cosmic_deltaXTop_BDT->SetLineWidth(2);          cosmic_deltaXTop_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaXTop_DLUboone->SetLineWidth(2);     cosmic_deltaXTop_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaXTop_DLNuE->SetLineWidth(2);        cosmic_deltaXTop_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaXTop_BDT->Integral() != 0)       cosmic_deltaXTop_BDT->Scale(1.0 / cosmic_deltaXTop_BDT->Integral());
    if(cosmic_deltaXTop_DLUboone->Integral() != 0)  cosmic_deltaXTop_DLUboone->Scale(1.0 / cosmic_deltaXTop_DLUboone->Integral());
    if(cosmic_deltaXTop_DLNuE->Integral() != 0)     cosmic_deltaXTop_DLNuE->Scale(1.0 / cosmic_deltaXTop_DLNuE->Integral());
    
    cosmic_deltaYTop_BDT->SetLineWidth(2);          cosmic_deltaYTop_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaYTop_DLUboone->SetLineWidth(2);     cosmic_deltaYTop_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaYTop_DLNuE->SetLineWidth(2);        cosmic_deltaYTop_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaYTop_BDT->Integral() != 0)       cosmic_deltaYTop_BDT->Scale(1.0 / cosmic_deltaYTop_BDT->Integral());
    if(cosmic_deltaYTop_DLUboone->Integral() != 0)  cosmic_deltaYTop_DLUboone->Scale(1.0 / cosmic_deltaYTop_DLUboone->Integral());
    if(cosmic_deltaYTop_DLNuE->Integral() != 0)     cosmic_deltaYTop_DLNuE->Scale(1.0 / cosmic_deltaYTop_DLNuE->Integral());
    
    cosmic_deltaZTop_BDT->SetLineWidth(2);          cosmic_deltaZTop_BDT->SetLineColor(TColor::GetColor("#bd1f01"));
    cosmic_deltaZTop_DLUboone->SetLineWidth(2);     cosmic_deltaZTop_DLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    cosmic_deltaZTop_DLNuE->SetLineWidth(2);        cosmic_deltaZTop_DLNuE->SetLineColor(TColor::GetColor("#3f90da"));
    
    if(cosmic_deltaZTop_BDT->Integral() != 0)       cosmic_deltaZTop_BDT->Scale(1.0 / cosmic_deltaZTop_BDT->Integral());
    if(cosmic_deltaZTop_DLUboone->Integral() != 0)  cosmic_deltaZTop_DLUboone->Scale(1.0 / cosmic_deltaZTop_DLUboone->Integral());
    if(cosmic_deltaZTop_DLNuE->Integral() != 0)     cosmic_deltaZTop_DLNuE->Scale(1.0 / cosmic_deltaZTop_DLNuE->Integral());

    TLegend* legendTop = new TLegend(0.69, 0.26, 0.87, 0.137);
    legendTop->SetBorderSize(0);
    legendTop->SetFillStyle(0);
    legendTop->AddEntry(cosmic_deltaTop_BDT_cf, "BDT", "f");
    legendTop->AddEntry(cosmic_deltaTop_DLUboone_cf, "DL Uboone", "f");
    legendTop->AddEntry(cosmic_deltaTop_DLNuE_cf, "DL Nu+E", "f");
    legendTop->SetTextSize(0.01);
    legendTop->SetMargin(0.12);

    cosmic_deltaTop_BDT_cf->Draw("hist");
    cosmic_deltaTop_DLUboone_cf->Draw("histsame");
    cosmic_deltaTop_DLNuE_cf->Draw("histsame");
    legendTop->Draw();
    canvasTopAll->SaveAs((base_path + "cosmics_distanceTopRecoVertex_cf_all.pdf").c_str());
    
    canvasTopAll->Clear();
    canvasTopAll->Update();
    cosmic_deltaTop_BDT->Draw("hist");
    cosmic_deltaTop_DLUboone->Draw("histsame");
    cosmic_deltaTop_DLNuE->Draw("histsame");
    legendTop->Draw();
    canvasTopAll->SaveAs((base_path + "cosmics_distanceTopRecoVertex_all.pdf").c_str());
    
    canvasTopAll->Clear();
    canvasTopAll->Update();
    cosmic_deltaXTop_BDT->Draw("hist");
    cosmic_deltaXTop_DLUboone->Draw("histsame");
    cosmic_deltaXTop_DLNuE->Draw("histsame");
    legendTop->Draw();
    canvasTopAll->SaveAs((base_path + "cosmics_distanceXTopRecoVertex_all.pdf").c_str());
    
    canvasTopAll->Clear();
    canvasTopAll->Update();
    cosmic_deltaYTop_BDT->Draw("hist");
    cosmic_deltaYTop_DLUboone->Draw("histsame");
    cosmic_deltaYTop_DLNuE->Draw("histsame");
    legendTop->Draw();
    canvasTopAll->SaveAs((base_path + "cosmics_distanceYTopRecoVertex_all.pdf").c_str());
    
    canvasTopAll->Clear();
    canvasTopAll->Update();
    cosmic_deltaZTop_BDT->Draw("hist");
    cosmic_deltaZTop_DLUboone->Draw("histsame");
    cosmic_deltaZTop_DLNuE->Draw("histsame");
    legendTop->Draw();
    canvasTopAll->SaveAs((base_path + "cosmics_distanceZTopRecoVertex_all.pdf").c_str());
    delete canvasTopAll;

    std::cout << "Number of reco vertices closer to the top y plane intersection = " << startOfTrack << std::endl;
    std::cout << "Number of reco vertices closer to the end/exit point = " << bottomOfTrack << std::endl;
    std::cout << "Number of vertices closer to the end/exit point, EXITED = " << bottomOfTrackExit << std::endl;
    std::cout << "Number of vertices closer to the end/exit point, ENDED = " << bottomOfTrackEnd << std::endl;
}
