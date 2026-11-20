// ROOT macro to make the nu score plots
// Usage: root -l -b -q 'MacroNuScoreAna.C("input.root")'

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"

// Define cuts
TCut cutFV = "abs(NuvX) > 1.5 && abs(NuvX) < 190 && abs(NuvY) < 190 && NuvZ > 10 && NuvZ < 490";
TCut cutAV = "abs(nuvX) < 200 && abs(nuvY) < 200 && NuvZ > 0 && NuvZ < 500";
TCut cutCC = "IntCCNC == 0";
TCut cutNC = "IntCCNC == 1";
TCut cutQE = "IntMode == 0";
TCut cutRES = "IntMode == 1";
TCut cutDIS = "IntMode == 2";
TCut cutMuon = "abs(IntNuPDG) == 14";
TCut cutLambda = "IntNLambda > 0";

TCut cutLambdaQE = cutFV && cutMuon && cutQE && cutLambda;
TCut cutNuQE =  cutFV && cutMuon && cutQE;

// Colors    
std::vector<int> fColors = {kRed-3,  kAzure-5, kGreen+3, kOrange-3, kMagenta+1, kCyan- 3, kYellow+2, kViolet-1, kTeal-1};
std::vector<int> fLineStyle = {1, 2, 9, 3, 6, 10, 7, 4, 5, 8};
std::vector<int> fMarkerStyle = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

int fNBins = 41;
double fXMin = -10;
double fXMax = 10;

std::map<std::string, TCut> fCuts = {
    {"#nu_{#mu} CC QE", cutNuQE}
    ,{"#bar{#nu}_{#mu} QE #Lambda", cutLambdaQE}
};



void MacroRecoVertexResolution(){
    std::cout << "Loaded MacroNuScoreAna" << std::endl;
    return;
}

THStack * PlotResolution(TTree *tree, std::string var, TPad *pad, std::vector<TH1F*> &hists){

    pad->cd();
    std::string title = "Resolution "+var+";"+var+"_{rec} - "+var+"_{true};AU";
    std::string varName = "Recnuv"+var+" - Nuv"+var;
    THStack *hs = new THStack(title.c_str(), title.c_str());
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    size_t i = 0;
    for (auto const& cut : fCuts){
        std::string histName = var + cut.first;
        TH1F *hRes = new TH1F(histName.c_str(), cut.first.c_str(), fNBins, fXMin, fXMax);
        hRes->SetLineColor(fColors[i]);
        hRes->SetLineStyle(fLineStyle[i]);
        hRes->SetLineWidth(2);    
        // Draw from TTree
        tree->Draw( (varName + ">>"+histName).c_str(), cut.second);
        hRes->Scale(1/hRes->Integral());
        hists.push_back(hRes);

        // Legend
        std::string label = cut.first+", ";
        double std_dev = hRes->GetStdDev();
        // format 1 decimal
        label += Form(" #sigma=%.1f", std_dev);


        leg->AddEntry(hRes, label.c_str(), "l");
        
        hs->Add( hists.back() );
        i++;
    }
    hs->Draw("nostack hist");
    leg->Draw("same");
    return hs;
}

void RunRecoVertexResolution(std::string fFileName, std::string fTreeName="originsAna/LambdaAnaTree"){

    double fNuScoreCut = 0.0;
    // Load tree
    TFile *file = new TFile(fFileName.c_str(), "READ");
    TTree *tree = (TTree*)file->Get(fTreeName.c_str());
    tree->Print();

    // Define histograms signals and backgrounds
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    // Histogram line width
    gStyle->SetHistLineWidth(4);
    // bottom margin
    gStyle->SetPadBottomMargin(0.15);
    // Left margin
    gStyle->SetPadLeftMargin(0.15);
    // Y offset
    gStyle->SetTitleOffset(1.2, "Y");


    // TCanvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    // TPads
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 0.5, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0.5, 1, 1);
    TPad *pad3 = new TPad("pad3", "pad3", 0, 0, 0.5, 0.5);
    TPad *pad4 = new TPad("pad4", "pad4", 0.5, 0, 1, 0.5);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();


    
    //gPad->SetLogy();
   

    // X
    pad1->cd();
    std::vector<TH1F*> histsX;
    THStack *hsResX = PlotResolution(tree, "X", pad1, histsX);

    // Y
    pad2->cd();
    std::vector<TH1F*> histsY;
    THStack *hsResY = PlotResolution(tree, "Y", pad2, histsY);

    // Z
    pad3->cd();
    std::vector<TH1F*> histsZ;
    THStack *hsResZ = PlotResolution(tree, "Z", pad3, histsZ);


   

    c1->cd();
    c1->Update();
    c1->WaitPrimitive();
    c1->Update();
    c1->SaveAs("NuScore.pdf");


}
