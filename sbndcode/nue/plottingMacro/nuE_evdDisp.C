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

typedef struct{
    double uvz;
    double x;
    double plane;
    double slice;
    double pfp;
} hit_t;

typedef struct{
    double x;
    double y;
    double z;
    double TPCValid;
} trueVertex_t;

typedef struct{
    double x;
    double y;
    double z;
} recoVertex_t;

typedef struct{
    std::vector<hit_t> hitVec0;
    std::vector<hit_t> hitVec1;
    std::vector<hit_t> hitVec2;
    std::vector<std::vector<hit_t>> hitVec0PFP;
    std::vector<std::vector<hit_t>> hitVec1PFP;
    std::vector<std::vector<hit_t>> hitVec2PFP;
    std::vector<std::vector<hit_t>> hitVec0Slice;
    std::vector<std::vector<hit_t>> hitVec1Slice;
    std::vector<std::vector<hit_t>> hitVec2Slice;
    trueVertex_t trueVertex;
    recoVertex_t recoVertex;
    int numSlices;
    int numPFPs;
} event_t;

void drawEvtDispSlices(std::vector<event_t>& events, int eventIndex = 0, const std::string& DLCurrent = "Default", int xminU = 0, int xmaxU = 0, int yminU = 0, int ymaxU = 0, int xminV = 0, int xmaxV = 0, int yminV = 0, int ymaxV = 0, int xminZ = 0, int xmaxZ = 0, int yminZ = 0, int ymaxZ = 0){
    if (eventIndex >= events.size()+1) {
        std::cerr << "Invalid event index!" << std::endl;
        return;
    }
    std::cout << "eventIndex: " << eventIndex << std::endl;

    const event_t& evt = events[eventIndex];

    TCanvas* c = new TCanvas("c", "Event Display", 1800, 600);
    c->Divide(3, 1);

    std::vector<int> colors = {kRed, kBlue, kGreen+2, kCyan+2, kOrange, kViolet, kAzure-2, kPink+2, kSpring+5, kMagenta};

    c->cd(1);  // Go to the U plane pad
    TH2F* frameU = new TH2F("frameU", "", 10, xminU, xmaxU, 10, yminU, ymaxU);
    frameU->GetXaxis()->SetTitle("X (cm)");
    frameU->GetYaxis()->SetTitle("U (cm)");
    frameU->Draw();
    gPad->Update();
    gPad->GetFrame()->SetY1(yminU);  // Lower y-axis limit
    gPad->GetFrame()->SetY2(ymaxU);   // Upper y-axis limit
    gPad->Modified();

    c->cd(2);  // Go to the V plane pad
    TH2F* frameV = new TH2F("frameV", "", 10, xminV, xmaxV, 10, yminV, ymaxV);
    frameV->GetXaxis()->SetTitle("X (cm)");
    frameV->GetYaxis()->SetTitle("V (cm)");
    frameV->Draw();
    gPad->Update();
    gPad->GetFrame()->SetY1(yminV);  // Lower y-axis limit
    gPad->GetFrame()->SetY2(ymaxV);   // Upper y-axis limit
    gPad->Modified();

    c->cd(3);  // Go to the Z plane pad
    TH2F* frameZ = new TH2F("frameZ", "", 10, xminZ, xmaxZ, 10, yminZ, ymaxZ);
    frameZ->GetXaxis()->SetTitle("X (cm)");
    frameZ->GetYaxis()->SetTitle("Z (cm)");
    frameZ->Draw();
    gPad->Update();
    gPad->GetFrame()->SetY1(yminZ);  // Lower y-axis limit
    gPad->GetFrame()->SetY2(ymaxZ);   // Upper y-axis limit
    gPad->Modified();

    c->cd(1);
    gPad->SetLeftMargin(0.15);
    TH2F* frame0 = new TH2F("frame0", "Hits in U Plane;x (cm);U (cm)", 100, xminU, xmaxU, 100, yminU, ymaxU);
    frame0->Draw();

    for (size_t sliceID = 0; sliceID < evt.hitVec0Slice.size(); ++sliceID) {
        const auto& hits = evt.hitVec0Slice[sliceID];
        if (hits.empty()) continue;

            TGraph* gr = new TGraph(hits.size());
            for (size_t i = 0; i < hits.size(); ++i) {
                gr->SetPoint(i, hits[i].x, hits[i].uvz);  // x vs uvz
            }
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.5);

            if(sliceID == evt.hitVec0Slice.size() - 1){
            gr->SetMarkerColor(kBlack);
            } else{
                gr->SetMarkerColor(colors[sliceID % colors.size()]);
            }

            gr->SetMinimum(yminU);
            gr->SetMaximum(ymaxU);
            gr->Draw("P SAME");
    }

    c->cd(2);
    gPad->SetLeftMargin(0.15);
    TH2F* frame1 = new TH2F("frame1", "Hits in V Plane;x (cm);V (cm)", 100, xminV, xmaxV, 100, yminV, ymaxV);
    frame1->Draw();

    std::cout << "num slices: " << evt.hitVec1Slice.size()-1 << std::endl;
    for (size_t sliceID = 0; sliceID < evt.hitVec1Slice.size(); ++sliceID) {
        const auto& hits = evt.hitVec1Slice[sliceID];
        if (hits.empty()) continue;

        TGraph* gr = new TGraph(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
            gr->SetPoint(i, hits[i].x, hits[i].uvz);  // x vs uvz
        }
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);

        if(sliceID == evt.hitVec1Slice.size() - 1){
            gr->SetMarkerColor(kBlack);
        } else{
            gr->SetMarkerColor(colors[sliceID % colors.size()]);
        }

        gr->SetMinimum(yminV);
        gr->SetMaximum(ymaxV);
        gr->Draw("P SAME");
    }

    c->cd(3);
    gPad->SetLeftMargin(0.15);
    TH2F* frame2 = new TH2F("frame2", "Hits in Z Plane;x (cm);Z (cm)", 100, xminZ, xmaxZ, 100, yminZ, ymaxZ);
    frame2->Draw();

    for (size_t sliceID = 0; sliceID < evt.hitVec2Slice.size(); ++sliceID) {
        const auto& hits = evt.hitVec2Slice[sliceID];
        if (hits.empty()) continue;

        TGraph* gr = new TGraph(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
        gr->SetPoint(i, hits[i].x, hits[i].uvz);  // x vs uvz
                                        }
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);

        if(sliceID == evt.hitVec2Slice.size() - 1){  
            gr->SetMarkerColor(kBlack);
        } else{
            gr->SetMarkerColor(colors[sliceID % colors.size()]);
        }

        gr->SetMinimum(yminZ);
        gr->SetMaximum(ymaxZ);
        gr->Draw("P SAME");
    }
    gStyle->SetOptStat(0);

    std::ostringstream filename;
    filename << "/nashome/c/coackley/nuEPlots/EvtDisp_Slice_" << DLCurrent << "_" << eventIndex+1 << ".pdf";
    c->SaveAs(filename.str().c_str());

}

void drawEvtDispPFPs(std::vector<event_t>& events, int eventIndex = 0, const std::string& DLCurrent = "Default", int xminU = 0, int xmaxU = 0, int yminU = 0, int ymaxU = 0, int xminV = 0, int xmaxV = 0, int yminV = 0, int ymaxV = 0, int xminZ = 0, int xmaxZ = 0, int yminZ = 0, int ymaxZ = 0){
    if (eventIndex >= events.size()+1) {
        std::cerr << "Invalid event index!" << std::endl;
        return;
    }
    std::cout << "eventIndex: " << eventIndex << std::endl;

    const event_t& evt = events[eventIndex];

    TCanvas* c = new TCanvas("c", "Event Display", 1800, 600);
    c->Divide(3, 1);

    std::vector<int> colors = {kRed, kBlue, kGreen+2, kCyan+2, kOrange, kViolet, kAzure-2, kPink+2, kSpring+5, kMagenta};

    c->cd(1);  // Go to the V plane pad
    TH2F* frameU = new TH2F("frameU", "", 10, xminU, xmaxU, 10, yminU, ymaxU); 
    frameU->GetXaxis()->SetTitle("X (cm)");
    frameU->GetYaxis()->SetTitle("U (cm)");
    frameU->Draw();
    gPad->Update();
    gPad->GetFrame()->SetY1(yminU);  // Lower y-axis limit
    gPad->GetFrame()->SetY2(ymaxU);   // Upper y-axis limit
    gPad->Modified();
    gPad->Update();

    c->cd(2);  // Go to the V plane pad
    TH2F* frameV = new TH2F("frameV", "", 10, xminV, xmaxV, 10, yminV, ymaxV);
    frameV->GetXaxis()->SetTitle("X (cm)");
    frameV->GetYaxis()->SetTitle("V (cm)");
    frameV->Draw();
    gPad->Update();
    gPad->GetFrame()->SetY1(yminV);  // Lower y-axis limit
    gPad->GetFrame()->SetY2(ymaxV);   // Upper y-axis limit
    gPad->Modified();
    gPad->Update();

    c->cd(3);  // Go to the V plane pad
    TH2F* frameZ = new TH2F("frameZ", "", 10, xminZ, xmaxZ, 10, yminZ, ymaxZ);
    frameZ->GetXaxis()->SetTitle("X (cm)");
    frameZ->GetYaxis()->SetTitle("Z (cm)");
    frameZ->Draw();
    gPad->Update();
    gPad->GetFrame()->SetY1(yminZ);  // Lower y-axis limit
    gPad->GetFrame()->SetY2(ymaxZ);   // Upper y-axis limit
    gPad->Modified();
    gPad->Update();
    
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    TH2F* frame0 = new TH2F("frame0", "Hits in U Plane;x (cm);U (cm)", 100, xminU, xmaxU, 100, yminU, ymaxU);
    frame0->Draw();

    for (size_t pfpID = 0; pfpID < evt.hitVec0PFP.size(); ++pfpID) {
        const auto& hits = evt.hitVec0PFP[pfpID];
        if (hits.empty()) continue;

        TGraph* gr = new TGraph(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
            gr->SetPoint(i, hits[i].x, hits[i].uvz);  // x vs uvz
        }
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);
        
        if(pfpID == evt.hitVec0PFP.size() - 1){
            gr->SetMarkerColor(kBlack);
        } else{
            gr->SetMarkerColor(colors[pfpID % colors.size()]);
        }

        gr->SetMinimum(yminU);
        gr->SetMaximum(ymaxU);        
        gr->Draw("P SAME");
    }

    c->cd(2);
    gPad->SetLeftMargin(0.15);
    TH2F* frame1 = new TH2F("frame1", "Hits in V Plane;x (cm);V (cm)", 100, xminV, xmaxV, 100, yminV, ymaxV);
    frame1->Draw();

    std::cout << "num PFPs: " << evt.hitVec1PFP.size()-1 << std::endl;
    for (size_t pfpID = 0; pfpID < evt.hitVec1PFP.size(); ++pfpID) {
        const auto& hits = evt.hitVec1PFP[pfpID];
        if (hits.empty()) continue;

        TGraph* gr = new TGraph(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
            gr->SetPoint(i, hits[i].x, hits[i].uvz);  // x vs uvz
        }
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);

        if(pfpID == evt.hitVec1PFP.size() - 1){
            gr->SetMarkerColor(kBlack);
        } else{
            gr->SetMarkerColor(colors[pfpID % colors.size()]);
        }

        gr->SetMinimum(yminV);
        gr->SetMaximum(ymaxV);        
        gr->Draw("P SAME");
    }

    c->cd(3);
    gPad->SetLeftMargin(0.15);
    TH2F* frame2 = new TH2F("frame2", "Hits in Z Plane;x (cm);Z (cm)", 100, xminZ, xmaxZ, 100, yminZ, ymaxZ);
    frame2->Draw();

    for (size_t pfpID = 0; pfpID < evt.hitVec2PFP.size(); ++pfpID) {
        const auto& hits = evt.hitVec2PFP[pfpID];
        if (hits.empty()) continue;

        TGraph* gr = new TGraph(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
            gr->SetPoint(i, hits[i].x, hits[i].uvz);  // x vs uvz
        }
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);
        
        if(pfpID == evt.hitVec2PFP.size() - 1){
            gr->SetMarkerColor(kBlack);
        } else{
            gr->SetMarkerColor(colors[pfpID % colors.size()]);
        }

        gr->SetMinimum(yminZ);
        gr->SetMaximum(ymaxZ);        
        gr->Draw("P SAME");
    }
    gStyle->SetOptStat(0);
    
    std::ostringstream filename;
    filename << "/nashome/c/coackley/nuEPlots/EvtDisp_PFP_" << DLCurrent << "_" << eventIndex+1 << ".pdf";
    c->SaveAs(filename.str().c_str());
}

void drawEvtDisp(const std::vector<event_t>& events, int eventIndex = 0, const std::string& DLCurrent = "Default"){
    if (eventIndex >= events.size()+1) {
        std::cerr << "Invalid event index!" << std::endl;
        return;
    }
    std::cout << "eventIndex: " << eventIndex << std::endl;
    
    const event_t& evt = events[eventIndex];
    TCanvas* c = new TCanvas("c", "Event Display", 1800, 600);
    c->Divide(3, 1);

    std::vector<hit_t> hitVecs[3] = {evt.hitVec0, evt.hitVec1, evt.hitVec2};
    const char* labels[3] = {"U", "V", "Z"};

    double thetaU = -60 * (TMath::Pi() / 180.0);
    double thetaV = 60 * (TMath::Pi() / 180.0);

    printf("True Vertex: (x, y ,z) = (%f, %f, %f)\n", evt.trueVertex.x, evt.trueVertex.y, evt.trueVertex.z);
    printf("Reco Vertex: (x, y ,z) = (%f, %f, %f)\n", evt.recoVertex.x, evt.recoVertex.y, evt.recoVertex.z);
    printf("(dx, dy ,dz) = (%f, %f, %f)\n", evt.recoVertex.x - evt.trueVertex.x, evt.recoVertex.y - evt.trueVertex.y, evt.recoVertex.z - evt.trueVertex.z);

    double trueVertexX = evt.trueVertex.x;
    double trueVertexU = (evt.trueVertex.z * TMath::Cos(thetaU)) - (evt.trueVertex.y * TMath::Sin(thetaU));
    double trueVertexV = (evt.trueVertex.z * TMath::Cos(thetaV)) - (evt.trueVertex.y * TMath::Sin(thetaV));
    double trueVertexZ = evt.trueVertex.z;

    //double recoVertexX = evt.recoVertex.x;
    //double recoVertexU = (evt.recoVertex.z * TMath::Sin(thetaU)) - (evt.recoVertex.y * TMath::Cos(thetaU));
    //double recoVertexV = (evt.recoVertex.z * TMath::Sin(thetaV)) - (evt.recoVertex.y * TMath::Cos(thetaV));
    //double recoVertexZ = evt.recoVertex.z;
    
    double recoVertexX = -151.137421;
    double recoVertexU = (289.797607 * TMath::Sin(thetaU)) - (-7.071669 * TMath::Cos(thetaU));
    double recoVertexV = (289.797607 * TMath::Sin(thetaV)) - (-7.071669 * TMath::Cos(thetaV));
    double recoVertexZ = 289.797607;

    for (int i = 0; i < 3; ++i) {
        int nHits = hitVecs[i].size();
        if (nHits == 0) continue;
   
        std::cout << "hitVec" << i << " has " << nHits << " hits" << std::endl;

        auto graph = new TGraph(nHits+2);
        auto graphTrue = new TGraph(1);
        auto graphReco = new TGraph(1);

        for (int j = 0; j < nHits; ++j) {
            graph->SetPoint(j, hitVecs[i][j].x, hitVecs[i][j].uvz);
        }
        
        if(i == 0){
            // U Plane
            graphTrue->SetPoint(0, trueVertexX, trueVertexU);
            graph->SetPoint(nHits, trueVertexX, trueVertexU);
            std::cout << "Drew red dot on u Plane at " << trueVertexX << " " << trueVertexU << std::endl;
            

            graphReco->SetPoint(0, recoVertexX, recoVertexU);
            graph->SetPoint(nHits+1, recoVertexX, recoVertexU);
            std::cout << "Drew blue dot on u Plane at " << recoVertexX << " " << recoVertexU << std::endl;
        } else if(i == 1){
            // V Plane
            graphTrue->SetPoint(0, trueVertexX, trueVertexV);
            graph->SetPoint(nHits, trueVertexX, trueVertexV);
            std::cout << "Drew red dot on V Plane at " << trueVertexX << " " << trueVertexV << std::endl;
            
            graphReco->SetPoint(0, recoVertexX, recoVertexV);
            graph->SetPoint(nHits+1, recoVertexX, recoVertexV);
            std::cout << "Drew blue dot on V Plane at " << recoVertexX << " " << recoVertexV << std::endl;
        } else if(i == 2){
            // Z Plane
            graphTrue->SetPoint(0, trueVertexX, trueVertexZ);
            graph->SetPoint(nHits, trueVertexX, trueVertexZ);
            std::cout << "Drew red dot on z Plane at " << trueVertexX << " " << trueVertexZ << std::endl;
            
            graphReco->SetPoint(0, recoVertexX, recoVertexZ);
            graph->SetPoint(nHits+1, recoVertexX, recoVertexZ);
            std::cout << "Drew blue dot on z Plane at " << recoVertexX << " " << recoVertexZ << std::endl;
        }

        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlack);
        graph->SetMarkerSize(0.5);
        
        graphReco->SetMarkerStyle(20);
        graphReco->SetMarkerColor(kBlue);
        graphReco->SetMarkerSize(1.0);
      
        graphTrue->SetMarkerStyle(20);
        graphTrue->SetMarkerColor(kRed);
        graphTrue->SetMarkerSize(1.0);

        c->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        graph->SetTitle(Form("%s-Plane;x (cm);%s (cm)", labels[i], labels[i]));
        graph->Draw("AP");
        graphTrue->Draw("P SAME");
        graphReco->Draw("P SAME");
    }

    c->Update();
    
    std::ostringstream filename;
    filename << "/nashome/c/coackley/nuEPlots/EvtDisp_" << DLCurrent << "_" << eventIndex+1 << ".pdf";
    c->SaveAs(filename.str().c_str());

}

void nuE_evdDisp(){

    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/analysed_Current/NoRefinement/CRUMBS/1.root");
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/analysed_DL_uboone/CRUMBS/1.root");
    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_05/NuEAnalyserOutput.root");
    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_05/srcs/sbndcode/sbndcode/nue/NuEAnalyserOutput.root");
    
    if(!file){
        std::cerr << "Error opening the file" << std::endl;
        return;
    }

    TDirectory *dir = (TDirectory*)file->Get("ana");
    if(!dir){
        std::cerr << "Directory 'ana' not found" << std::endl;
        return;
    }

    TTree *tree = (TTree*)dir->Get("NuE");
    if(!tree){
        std::cerr << "NuE not found" << std::endl;
        return;
    }


    UInt_t eventID, runID, subRunID;
    double DLCurrent;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);

    std::vector<double> *reco_hitPlane = nullptr;
    std::vector<double> *reco_hitX = nullptr;
    std::vector<double> *reco_hitUVZ = nullptr;
    std::vector<double> *reco_hitSlice = nullptr;
    std::vector<double> *reco_hitPFP = nullptr;

    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;

    std::vector<double> *truth_neutrinoVX = nullptr;
    std::vector<double> *truth_neutrinoVY = nullptr;
    std::vector<double> *truth_neutrinoVZ = nullptr;
    std::vector<double> *truth_neutrinoTPCValid = nullptr;

    std::vector<double> *reco_particlePDG = nullptr;
    std::vector<double> *reco_sliceID = nullptr;

    tree->SetBranchAddress("reco_hitPlane", &reco_hitPlane);
    tree->SetBranchAddress("reco_hitX", &reco_hitX);
    tree->SetBranchAddress("reco_hitUVZ", &reco_hitUVZ); 
    tree->SetBranchAddress("reco_hitSlice", &reco_hitSlice);
    tree->SetBranchAddress("reco_hitPFP", &reco_hitPFP);

    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);

    tree->SetBranchAddress("truth_neutrinoVX", &truth_neutrinoVX);
    tree->SetBranchAddress("truth_neutrinoVY", &truth_neutrinoVY);
    tree->SetBranchAddress("truth_neutrinoVZ", &truth_neutrinoVZ);
    tree->SetBranchAddress("truth_neutrinoTPCValid", &truth_neutrinoTPCValid);

    tree->SetBranchAddress("reco_particlePDG", &reco_particlePDG);
    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);

    Long64_t numEntries = tree->GetEntries();
    std::cout << "Number of Events: " << numEntries << std::endl;
    
    std::vector<event_t> events;

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);

        event_t eventInfo;

        std::vector<hit_t> hitsPlane0, hitsPlane1, hitsPlane2;

        std::vector<std::vector<hit_t>> hitsPlane0PFP, hitsPlane1PFP, hitsPlane2PFP;
        std::vector<std::vector<hit_t>> hitsPlane0Slice, hitsPlane1Slice, hitsPlane2Slice;

        int numPFPs = reco_particlePDG->size();
        int numSlices = reco_sliceID->size();

        if(numPFPs == 1 && reco_particlePDG->at(0) == -999999) numPFPs = 0;
        if(numSlices == 1 && reco_sliceID->at(0) == -999999) numSlices = 0;

        hitsPlane0PFP.resize(numPFPs+1);
        hitsPlane1PFP.resize(numPFPs+1);
        hitsPlane2PFP.resize(numPFPs+1);
        hitsPlane0Slice.resize(numSlices+1);
        hitsPlane1Slice.resize(numSlices+1);
        hitsPlane2Slice.resize(numSlices+1);

        std::cout << "reco_hitX size: " << reco_hitX->size() << std::endl;
        for(size_t j = 0; j < reco_hitX->size(); ++j){
            hit_t hit;
            hit.x = reco_hitX->at(j);
            hit.uvz = reco_hitUVZ->at(j);
            hit.plane = reco_hitPlane->at(j);
            hit.slice = reco_hitSlice->at(j);
            hit.pfp = reco_hitPFP->at(j);

            if(reco_hitPlane->at(j) == 0){
                hitsPlane0.push_back(hit);
            } else if(reco_hitPlane->at(j) == 1){
                hitsPlane1.push_back(hit);
            } else if(reco_hitPlane->at(j) == 2){
                hitsPlane2.push_back(hit);
            }

            if(hit.pfp == -999999){
                if(reco_hitPlane->at(j) == 0){
                    hitsPlane0PFP[numPFPs].push_back(hit);
                } else if(reco_hitPlane->at(j) == 1){
                    hitsPlane1PFP[numPFPs].push_back(hit);
                } else if(reco_hitPlane->at(j) == 2){
                    hitsPlane2PFP[numPFPs].push_back(hit);
                }
            } else{
                if(reco_hitPlane->at(j) == 0){
                    hitsPlane0PFP[hit.pfp].push_back(hit);
                } else if(reco_hitPlane->at(j) == 1){
                    hitsPlane1PFP[hit.pfp].push_back(hit);
                } else if(reco_hitPlane->at(j) == 2){
                    hitsPlane2PFP[hit.pfp].push_back(hit);
                }
            }

            if(hit.slice == -999999){
                if(reco_hitPlane->at(j) == 0){
                    hitsPlane0Slice[numSlices].push_back(hit);
                } else if(reco_hitPlane->at(j) == 1){
                    hitsPlane1Slice[numSlices].push_back(hit);
                } else if(reco_hitPlane->at(j) == 2){
                    hitsPlane2Slice[numSlices].push_back(hit);
                }
               
            } else if(hit.slice < 20){
                if(reco_hitPlane->at(j) == 0){
                    hitsPlane0Slice[hit.slice].push_back(hit);
                } else if(reco_hitPlane->at(j) == 1){
                    hitsPlane1Slice[hit.slice].push_back(hit);
                } else if(reco_hitPlane->at(j) == 2){
                    hitsPlane2Slice[hit.slice].push_back(hit);
                }
            }
        }
   
        recoVertex_t recoVertex; 
        for(size_t j = 0; j < reco_neutrinoVX->size(); ++j){
            recoVertex.x = reco_neutrinoVX->at(j);
            recoVertex.y = reco_neutrinoVY->at(j);
            recoVertex.z = reco_neutrinoVZ->at(j);           
            
            if(reco_neutrinoVX->size() > 1){
                recoVertex.x = -999999;
                recoVertex.y = -999999;
                recoVertex.z = -999999;
            }
        }

        trueVertex_t trueVertex;
        for(size_t j = 0; j < truth_neutrinoVX->size(); ++j){
            trueVertex.x = truth_neutrinoVX->at(j);
            trueVertex.y = truth_neutrinoVY->at(j);
            trueVertex.z = truth_neutrinoVZ->at(j);
            trueVertex.TPCValid = truth_neutrinoTPCValid->at(j);
            
            if(truth_neutrinoVX->size() > 1){
                trueVertex.x = -999999;
                trueVertex.y = -999999;
                trueVertex.z = -999999;
                trueVertex.TPCValid = -999999;
            }
        }

        eventInfo.recoVertex = recoVertex;
        eventInfo.trueVertex = trueVertex;
        eventInfo.hitVec0 = hitsPlane0;
        eventInfo.hitVec1 = hitsPlane1;
        eventInfo.hitVec2 = hitsPlane2;
        eventInfo.hitVec0PFP = hitsPlane0PFP;
        eventInfo.hitVec1PFP = hitsPlane1PFP;
        eventInfo.hitVec2PFP = hitsPlane2PFP;
        eventInfo.hitVec0Slice = hitsPlane0Slice;
        eventInfo.hitVec1Slice = hitsPlane1Slice;
        eventInfo.hitVec2Slice = hitsPlane2Slice;
        eventInfo.numPFPs = numPFPs;
        eventInfo.numSlices = numSlices;
        
        events.push_back(eventInfo);
    
        std::cout << "num slices = " << numSlices << ", num PFPs = " << numPFPs << std::endl;
        std::cout << "hitsPlane0 size = " << hitsPlane0.size() << ", hitsPlane1 size = " << hitsPlane1.size() << ", hitsPlane2 size = " << hitsPlane2.size() << std::endl;
    }

    int xminU = 50;
    int xmaxU = 200;
    int yminU = 0;
    int ymaxU = 200;
    
    int xminV = -260;
    int xmaxV = 200;
    int yminV = 70;
    int ymaxV = 280;

    int xminZ = 50;
    int xmaxZ = 200;
    int yminZ = 140;
    int ymaxZ = 450;

    drawEvtDisp(events, 43, "DLUboone");
    drawEvtDispPFPs(events, 43, "DLUboone", xminU, xmaxU, yminU, ymaxU, xminV, xmaxV, yminV, ymaxV, xminZ, xmaxZ, yminZ, ymaxZ);
    drawEvtDispSlices(events, 43, "DLUboone", xminU, xmaxU, yminU, ymaxU, xminV, xmaxV, yminV, ymaxV, xminZ, xmaxZ, yminZ, ymaxZ);
}


