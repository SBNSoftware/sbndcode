#include <vector>
#include <map>
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

typedef struct{
    double event;
    double run;
    double subrun;
    double trueNeutrinoVX;
    double trueNeutrinoVY;
    double trueNeutrinoVZ;
    double recoNeutrinoVX;
    double recoNeutrinoVY;
    double recoNeutrinoVZ;
    int trueCCNC;
    int trueNeutrinoType;
    int trueLeptonType;
    int nSlices;
    int nPFPs;
} event_t;

void numberPlots(std::vector<event_t> allEvents_vec){
    size_t numEvents = allEvents_vec.size();
}

void deltaVertices(std::vector<event_t> allEvents_vec){
    size_t numEvents = allEvents_vec.size();
    
    TCanvas *deltaXCanvas = new TCanvas("deltaX_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaX_dist = new TH1F("Delta X", "X Coord Diff", 40, -10, 10);
    deltaX_dist->SetTitle("Neutrino Vertex deltaX;True - Reco X Coordinate;#");

    TCanvas *deltaYCanvas = new TCanvas("deltaY_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaY_dist = new TH1F("Delta Y", "Y Coord Diff", 40, -10, 10);
    deltaY_dist->SetTitle("Neutrino Vertex deltaY;True - Reco Y Coordinate;#");

    TCanvas *deltaZCanvas = new TCanvas("deltaZ_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaZ_dist = new TH1F("Delta Z", "Z Coord Diff", 40, -10, 10);
    deltaZ_dist->SetTitle("Neutrino Vertex deltaZ;True - Reco Z Coordinate;#");

    TCanvas *deltaRCanvas = new TCanvas("deltaR_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaR_dist = new TH1F("Delta R", "R Coord Diff", 40, 0, 20);
    deltaR_dist->SetTitle("Neutrino Vertex deltaR;True - Reco R Coordinate;#");

    for(UInt_t i = 0; i < numEvents; i++){
        event_t event = allEvents_vec.at(i);
        deltaX_dist->Fill(event.trueNeutrinoVX - event.recoNeutrinoVX);
        deltaY_dist->Fill(event.trueNeutrinoVY - event.recoNeutrinoVY);
        deltaZ_dist->Fill(event.trueNeutrinoVZ - event.recoNeutrinoVZ);
        deltaR_dist->Fill(std::sqrt( std::pow(event.trueNeutrinoVX - event.recoNeutrinoVX, 2) + std::pow(event.trueNeutrinoVY - event.recoNeutrinoVY, 2) + std::pow(event.trueNeutrinoVZ - event.recoNeutrinoVZ, 2) ));
    }

    deltaXCanvas->cd();
    deltaX_dist->Draw();
    deltaXCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaX_dist.png");

    deltaYCanvas->cd();
    deltaY_dist->Draw();
    deltaYCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaY_dist.png");

    deltaZCanvas->cd();
    deltaZ_dist->Draw();
    deltaZCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaZ_dist.png");

    deltaRCanvas->cd();
    deltaR_dist->Draw();
    deltaRCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaR_dist.png");
}

void nuE_macro(){
    //TFile *f = TFile::Open("/exp/sbnd/app/users/coackley/nuev09_91_02/srcs/sbndcode/sbndcode/NuE/NuEAnalyserOutput.root", "READ");
    TFile *f = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_03/NuEAnalyserOutput.root", "READ");
    if(!f){
        std::cout << "Failed to read file" << std::endl;
        return;
    }

    TTree *t;
    f->GetObject("NuETree", t);

    std::vector<double> allEvent = std::vector<double>(0);
    std::vector<double> allRun = std::vector<double>(0);
    std::vector<double> allSubrun = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVX = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVY = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVZ = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVX = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVY = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVZ = std::vector<double>(0);
    std::vector<int> allTrueCCNC = std::vector<int>(0);
    std::vector<int> allTrueNeutrinoType = std::vector<int>(0);
    std::vector<int> allTrueLeptonType = std::vector<int>(0);
    std::vector<int> allNumSlices = std::vector<int>(0);
    std::vector<int> allNumPfps = std::vector<int>(0);

    Long64_t numEntries = t->GetEntries();
    int fileNumber = 1;

    for(Long64_t i = 0; i < numEntries; i++){
        std::cout << "File number: " << fileNumber << std::endl;
        std::vector<double> *pEvent = 0;
        std::vector<double> *pRun = 0; 
        std::vector<double> *pSubrun = 0;
        std::vector<double> *pTrueNeutrinoVX = 0;
        std::vector<double> *pTrueNeutrinoVY = 0;
        std::vector<double> *pTrueNeutrinoVZ = 0;
        std::vector<double> *pRecoNeutrinoVX = 0;
        std::vector<double> *pRecoNeutrinoVY = 0;
        std::vector<double> *pRecoNeutrinoVZ = 0;
        std::vector<int> *pTrueCCNC = 0;
        std::vector<int> *pTrueNeutrinoType = 0;
        std::vector<int> *pTrueLeptonType = 0;
        std::vector<int> *pNumSlices = 0;
        std::vector<int> *pNumPfps = 0;

        TBranch *event_branch = 0;
        TBranch *run_branch = 0;
        TBranch *subrun_branch = 0;
        TBranch *trueNeutrinoVX_branch = 0;
        TBranch *trueNeutrinoVY_branch = 0;
        TBranch *trueNeutrinoVZ_branch = 0;
        TBranch *recoNeutrinoVX_branch = 0;
        TBranch *recoNeutrinoVY_branch = 0;
        TBranch *recoNeutrinoVZ_branch = 0;
        TBranch *trueCCNC_branch = 0;
        TBranch *trueNeutrinoType_branch = 0;
        TBranch *trueLeptonType_branch = 0;
        TBranch *numSlices_branch = 0;
        TBranch *numPfps_branch = 0;

        t->SetBranchAddress("event_tree", &pEvent, &event_branch);
        t->SetBranchAddress("run_tree", &pRun, &run_branch);
        t->SetBranchAddress("subrun_tree", &pSubrun, &subrun_branch);
        t->SetBranchAddress("trueNeutrinoVX_tree", &pTrueNeutrinoVX, &trueNeutrinoVX_branch);
        t->SetBranchAddress("trueNeutrinoVY_tree", &pTrueNeutrinoVY, &trueNeutrinoVY_branch);
        t->SetBranchAddress("trueNeutrinoVZ_tree", &pTrueNeutrinoVZ, &trueNeutrinoVZ_branch);
        t->SetBranchAddress("recoNeutrinoVX_tree", &pRecoNeutrinoVX, &recoNeutrinoVX_branch);
        t->SetBranchAddress("recoNeutrinoVY_tree", &pRecoNeutrinoVY, &recoNeutrinoVY_branch);
        t->SetBranchAddress("recoNeutrinoVZ_tree", &pRecoNeutrinoVZ, &recoNeutrinoVZ_branch);
        t->SetBranchAddress("trueCCNC_tree", &pTrueCCNC, &trueCCNC_branch);
        t->SetBranchAddress("trueNeutrinoType_tree", &pTrueNeutrinoType, &trueNeutrinoType_branch);
        t->SetBranchAddress("trueLeptonType_tree", &pTrueLeptonType, &trueLeptonType_branch);
        t->SetBranchAddress("numSlices_tree", &pNumSlices, &numSlices_branch);
        t->SetBranchAddress("numPfps_tree", &pNumPfps, &numPfps_branch);

        event_branch->GetEntry(i);
        run_branch->GetEntry(i);
        subrun_branch->GetEntry(i);
        trueNeutrinoVX_branch->GetEntry(i);
        trueNeutrinoVY_branch->GetEntry(i);
        trueNeutrinoVZ_branch->GetEntry(i);
        recoNeutrinoVX_branch->GetEntry(i);
        recoNeutrinoVY_branch->GetEntry(i);
        recoNeutrinoVZ_branch->GetEntry(i);
        trueCCNC_branch->GetEntry(i);
        trueNeutrinoType_branch->GetEntry(i);
        trueLeptonType_branch->GetEntry(i);
        numSlices_branch->GetEntry(i);
        numPfps_branch->GetEntry(i);

        size_t vecSize = pEvent->size();
        for(UInt_t j = 0; j < vecSize; j++){
            allEvent.push_back(pEvent->at(j));
            allRun.push_back(pRun->at(j));
            allSubrun.push_back(pSubrun->at(j));
            allTrueNeutrinoVX.push_back(pTrueNeutrinoVX->at(j));
            allTrueNeutrinoVY.push_back(pTrueNeutrinoVY->at(j));
            allTrueNeutrinoVZ.push_back(pTrueNeutrinoVZ->at(j));
            allRecoNeutrinoVX.push_back(pRecoNeutrinoVX->at(j));
            allRecoNeutrinoVY.push_back(pRecoNeutrinoVY->at(j));
            allRecoNeutrinoVZ.push_back(pRecoNeutrinoVZ->at(j));
            allTrueCCNC.push_back(pTrueCCNC->at(j));
            allTrueNeutrinoType.push_back(pTrueNeutrinoType->at(j));
            allTrueLeptonType.push_back(pTrueLeptonType->at(j));
            allNumSlices.push_back(pNumSlices->at(j));
            allNumPfps.push_back(pNumPfps->at(j));
        }

        fileNumber++;
    }

    std::vector<event_t> allEvents = std::vector<event_t>();

    size_t newNumEvents = allEvent.size();

    for(UInt_t k = 0; k < newNumEvents; k++){
        event_t event;
        event.event = allEvent.at(k);
        event.run = allRun.at(k);
        event.subrun = allSubrun.at(k);
        event.trueNeutrinoVX = allTrueNeutrinoVX.at(k);
        event.trueNeutrinoVY = allTrueNeutrinoVY.at(k);
        event.trueNeutrinoVZ = allTrueNeutrinoVZ.at(k);
        event.recoNeutrinoVX = allRecoNeutrinoVX.at(k);
        event.recoNeutrinoVY = allRecoNeutrinoVY.at(k);
        event.recoNeutrinoVZ = allRecoNeutrinoVZ.at(k);
        event.trueCCNC = allTrueCCNC.at(k);
        event.trueNeutrinoType = allTrueNeutrinoType.at(k);
        event.trueLeptonType = allTrueLeptonType.at(k);
        event.nSlices = allNumSlices.at(k);
        event.nPFPs = allNumPfps.at(k);
        allEvents.push_back(event);
    }

    deltaVertices(allEvents);
    numberPlots(allEvents);
}
