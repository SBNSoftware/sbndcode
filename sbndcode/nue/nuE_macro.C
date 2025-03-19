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
    int nTracks;
    int nShowers;
    double tpcID;
    int dl_current;
} event_t;

void numberPlots(std::vector<event_t> allEvents_vec){
    size_t numEvents = allEvents_vec.size();

    std::vector<event_t> dlDune = std::vector<event_t>();
    std::vector<event_t> dlUboone = std::vector<event_t>();
    std::vector<event_t> current = std::vector<event_t>();

    for(UInt_t i = 0; i < numEvents; i++){
       event_t event = allEvents_vec.at(i);
       if(event.dl_current == 0){
           dlUboone.push_back(event);
       } else if(event.dl_current == 1){
           dlDune.push_back(event);
       } else if(event.dl_current == 2){
           current.push_back(event);
       }
    }

    TCanvas *numSlicesCanvas = new TCanvas("numSlices_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numSlicesCurrent_dist = new TH1F("Number of Slices", "Num Slices", 5, -0.5, 4.5);
    numSlicesCurrent_dist->SetTitle("Distribution of the Number of Slices in an Event;Number of Slices in Event;# of Events");
    TH1F* numSlicesDLUboone_dist = (TH1F*) numSlicesCurrent_dist->Clone("Number of Slices");
    TH1F* numSlicesDLDune_dist = (TH1F*) numSlicesCurrent_dist->Clone("Number of Slices");

    TCanvas *numPFPsCanvas = new TCanvas("numPFPs_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numPFPsCurrent_dist = new TH1F("Number of PFPs", "Num PFPs", 10, -0.5, 9.5);
    numPFPsCurrent_dist->SetTitle("Distribution of the Number of PFPs in an Event;Number of PFPs in Event;# of Events");
    TH1F* numPFPsDLUboone_dist = (TH1F*) numPFPsCurrent_dist->Clone("Number of PFPs");
    TH1F* numPFPsDLDune_dist = (TH1F*) numPFPsCurrent_dist->Clone("Number of PFPs");

    TCanvas *numTracksCanvas = new TCanvas("numTracks_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numTracksCurrent_dist = new TH1F("Number of Tracks", "Num Tracks", 10, -0.5, 9.5);
    numTracksCurrent_dist->SetTitle("Distribution of the Number of Tracks in an Event;Number of Tracks in Event;# of Events");
    TH1F* numTracksDLUboone_dist = (TH1F*) numTracksCurrent_dist->Clone("Number of Tracks");
    TH1F* numTracksDLDune_dist = (TH1F*) numTracksCurrent_dist->Clone("Number of Tracks");

    TCanvas *numShowersCanvas = new TCanvas("numShowers_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numShowersCurrent_dist = new TH1F("Number of Showers", "Num Showers", 10, -0.5, 9.5);
    numShowersCurrent_dist->SetTitle("Distribution of the Number of Showers in an Event;Number of Showers in Event;# of Events");
    TH1F* numShowersDLUboone_dist = (TH1F*) numShowersCurrent_dist->Clone("Number of Showers");
    TH1F* numShowersDLDune_dist = (TH1F*) numShowersCurrent_dist->Clone("Number of Showers");

    for(UInt_t j = 0; j < current.size(); j++){
        numSlicesCurrent_dist->Fill(current.at(j).nSlices);
        numPFPsCurrent_dist->Fill(current.at(j).nPFPs);
        numTracksCurrent_dist->Fill(current.at(j).nTracks);
        numShowersCurrent_dist->Fill(current.at(j).nShowers);
    }

    for(UInt_t j = 0; j < dlDune.size(); j++){
        numSlicesDLDune_dist->Fill(dlDune.at(j).nSlices);
        numPFPsDLDune_dist->Fill(dlDune.at(j).nPFPs);
        numTracksDLDune_dist->Fill(dlDune.at(j).nTracks);
        numShowersDLDune_dist->Fill(dlDune.at(j).nShowers);
    }

    for(UInt_t j = 0; j < dlUboone.size(); j++){
        numSlicesDLUboone_dist->Fill(dlUboone.at(j).nSlices);
        numPFPsDLUboone_dist->Fill(dlUboone.at(j).nPFPs);
        numTracksDLUboone_dist->Fill(dlUboone.at(j).nTracks);
        numShowersDLUboone_dist->Fill(dlUboone.at(j).nShowers);
    }

    numSlicesCanvas->cd();
    numSlicesCurrent_dist->SetLineWidth(2);
    numSlicesCurrent_dist->SetLineColor(kRed);
    numSlicesDLDune_dist->SetLineWidth(2);
    numSlicesDLDune_dist->SetLineColor(kViolet-5);
    numSlicesDLUboone_dist->SetLineWidth(2);
    numSlicesDLUboone_dist->SetLineColor(kBlue);
    numSlicesCurrent_dist->Draw("hist");
    numSlicesDLDune_dist->Draw("histsame");
    numSlicesDLUboone_dist->Draw("histsame");
    numSlicesCurrent_dist->SetStats(0);
    numSlicesCurrent_dist->GetYaxis()->SetRangeUser(0, 1000);
    
    auto legend = new TLegend(0.56,0.86,0.88,0.70);
    legend->AddEntry(numSlicesDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(numSlicesDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(numSlicesCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();
    numSlicesCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numSlices_dist.pdf");
    
    numPFPsCanvas->cd();
    numPFPsCurrent_dist->SetLineWidth(2);
    numPFPsCurrent_dist->SetLineColor(kRed);
    numPFPsDLDune_dist->SetLineWidth(2);
    numPFPsDLDune_dist->SetLineColor(kViolet-5);
    numPFPsDLUboone_dist->SetLineWidth(2);
    numPFPsDLUboone_dist->SetLineColor(kBlue);
    numPFPsCurrent_dist->Draw("hist");
    numPFPsDLDune_dist->Draw("histsame");
    numPFPsDLUboone_dist->Draw("histsame");
    numPFPsCurrent_dist->SetStats(0);
    numPFPsCurrent_dist->GetYaxis()->SetRangeUser(0, 900);

    auto legend2 = new TLegend(0.56,0.86,0.88,0.70);
    legend2->AddEntry(numPFPsDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend2->AddEntry(numPFPsDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend2->AddEntry(numPFPsCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend2->SetTextSize(0.0225);
    legend2->SetMargin(0.13);
    legend2->Draw();
    numPFPsCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numPFPs_dist.pdf");

    numTracksCanvas->cd();
    numTracksCurrent_dist->SetLineWidth(2);
    numTracksCurrent_dist->SetLineColor(kRed);
    numTracksDLDune_dist->SetLineWidth(2);
    numTracksDLDune_dist->SetLineColor(kViolet-5);
    numTracksDLUboone_dist->SetLineWidth(2);
    numTracksDLUboone_dist->SetLineColor(kBlue);
    numTracksCurrent_dist->Draw("hist");
    numTracksDLDune_dist->Draw("histsame");
    numTracksDLUboone_dist->Draw("histsame");
    numTracksCurrent_dist->SetStats(0);
    numTracksCurrent_dist->GetYaxis()->SetRangeUser(0, 900);
    
    auto legend3 = new TLegend(0.56,0.86,0.88,0.70);
    legend3->AddEntry(numTracksDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3->AddEntry(numTracksDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3->AddEntry(numTracksCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend3->SetTextSize(0.0225);
    legend3->SetMargin(0.13);
    legend3->Draw();
    numTracksCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numTracks_dist.pdf");

    numShowersCanvas->cd();
    numShowersCurrent_dist->SetLineWidth(2);
    numShowersCurrent_dist->SetLineColor(kRed);
    numShowersDLDune_dist->SetLineWidth(2);
    numShowersDLDune_dist->SetLineColor(kViolet-5);
    numShowersDLUboone_dist->SetLineWidth(2);
    numShowersDLUboone_dist->SetLineColor(kBlue);
    numShowersCurrent_dist->Draw("hist");
    numShowersDLDune_dist->Draw("histsame");
    numShowersDLUboone_dist->Draw("histsame");
    numShowersCurrent_dist->SetStats(0);
    numShowersCurrent_dist->GetYaxis()->SetRangeUser(0, 900);

    auto legend4 = new TLegend(0.56,0.86,0.88,0.70);
    legend4->AddEntry(numShowersDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend4->AddEntry(numShowersDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend4->AddEntry(numShowersCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend4->SetTextSize(0.0225);
    legend4->SetMargin(0.13);
    legend4->Draw();
    numShowersCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numShowers_dist.pdf");
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
    deltaX_dist->SetLineWidth(2);
    deltaX_dist->SetStats(0);
    deltaXCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaX_dist_nuE.pdf");

    deltaYCanvas->cd();
    deltaY_dist->Draw();
    deltaY_dist->SetLineWidth(2);
    deltaY_dist->SetStats(0);
    deltaYCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaY_dist_nuE.pdf");

    deltaZCanvas->cd();
    deltaZ_dist->Draw();
    deltaZ_dist->SetLineWidth(2);
    deltaZ_dist->SetStats(0);
    deltaZCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaZ_dist_nuE.pdf");

    deltaRCanvas->cd();
    deltaR_dist->Draw();
    deltaR_dist->SetLineWidth(2);
    deltaR_dist->SetStats(0);
    deltaRCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaR_dist_nuE.pdf");
}

void nuE_macro(){
    //TFile *f = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_03/NuEAnalyserOutput.root", "READ");
    TFile *f = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/merged.root", "READ");
    
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
    std::vector<int> allNumTracks = std::vector<int>(0);
    std::vector<int> allNumShowers = std::vector<int>(0);
    std::vector<double> allTpcID = std::vector<double>(0);
    std::vector<int> allDl_current = std::vector<int>(0);

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
        std::vector<int> *pNumTracks = 0;
        std::vector<int> *pNumShowers = 0;
        std::vector<double> *pTpcID = 0;
        std::vector<int> *pDl_current = 0;

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
        TBranch *numTracks_branch = 0;
        TBranch *numShowers_branch = 0;
        TBranch *tpcID_branch = 0;
        TBranch *dl_current_branch = 0;

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
        t->SetBranchAddress("numTracks_tree", &pNumTracks, &numTracks_branch);
        t->SetBranchAddress("numShowers_tree", &pNumShowers, &numShowers_branch);
        t->SetBranchAddress("tpcID_tree", &pTpcID, &tpcID_branch);
        t->SetBranchAddress("dl_current_tree", &pDl_current, &dl_current_branch);

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
        numTracks_branch->GetEntry(i);
        numShowers_branch->GetEntry(i);
        tpcID_branch->GetEntry(i);
        dl_current_branch->GetEntry(i);

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
            allNumTracks.push_back(pNumTracks->at(j));
            allNumShowers.push_back(pNumShowers->at(j));
            allTpcID.push_back(pTpcID->at(j));
            allDl_current.push_back(pDl_current->at(j));
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
        event.nTracks = allNumTracks.at(k);
        event.nShowers = allNumShowers.at(k);
        event.tpcID = allTpcID.at(k);
        event.dl_current = allDl_current.at(k);
        allEvents.push_back(event);
    }

    //deltaVertices(allEvents);
    numberPlots(allEvents);
}
