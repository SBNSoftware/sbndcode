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
TCut cutFV = "abs(NuScoreMatchedX) > 1.5 && abs(NuScoreMatchedX) < 190 && abs(NuScoreMatchedY) < 190 && NuScoreMatchedZ > 10 && NuScoreMatchedZ < 490";
TCut cutAV = "abs(NuScoreMatchedX) < 200 && abs(NuScoreMatchedY) < 200 && NuScoreMatchedZ > 0 && NuScoreMatchedZ < 500";
TCut cutCC = "NuScoreMatchedCCNC == 0";
TCut cutNC = "NuScoreMatchedCCNC == 1";
TCut cutQE = "NuScoreMatchedIntMode == 0";
TCut cutRES = "NuScoreMatchedIntMode == 1";
TCut cutDIS = "NuScoreMatchedIntMode == 2";
TCut cutLambda = "NuScoreMatchedIntNLambda > 0";

TCut cutLambdaQE = "NuScoreMatchedOrigin == 1" && cutAV && cutQE && cutLambda;
TCut cutNu = "NuScoreMatchedOrigin == 1" && cutAV && !cutLambda;
TCut cutCosmic = "NuScoreMatchedOrigin == 2";
TCut cutDirt = "NuScoreMatchedOrigin == 1" && !cutAV;
TCut cutNuCC = "NuScoreMatchedOrigin == 1" && cutAV && cutCC;
TCut cutNuNC = "NuScoreMatchedOrigin == 1" && cutAV && cutNC;



std::map<std::string, TCut> fCutsGeneral = {
    {"#nu AV", cutNu}
    ,{"Cosmic", cutCosmic}
    //,{"Dirt", cutDirt}
};

std::map<std::string, TCut> fCutsNu = {
    {"#nu CC AV", cutNuCC}
    ,{"#nu NC AV", cutNuNC}
    ,{"#nu QE #Lambda AV", cutLambdaQE}
};

std::map<std::string, TCut> fCuts = fCutsGeneral;



void MacroNuScoreAna(){
    std::cout << "Loaded MacroNuScoreAna" << std::endl;
    return;
}

void RunNuScoreAna(std::string fFileName, std::string fTreeName="originsAna/TreeNuScore"){

    double fNuScoreCut = 0.0;
    // Load tree
    TFile *file = new TFile(fFileName.c_str(), "READ");
    TTree *tree = (TTree*)file->Get(fTreeName.c_str());
    tree->Print();

    // Load branches
    std::vector<float> *NuScoreVect = new std::vector<float>();
    std::vector<int> *NuScoreMatchedOrigin = new std::vector<int>();
    std::vector<int> *NuScoreMatchedCCNC = new std::vector<int>();
    std::vector<int> *NuScoreMatchedNuPDG = new std::vector<int>();
    std::vector<float> *NuScoreMatchedNuEnergy = new std::vector<float>();
    std::vector<int> *NuScoreMatchedNuIntMode = new std::vector<int>();
    std::vector<int> *NuScoreMatchedX = new std::vector<int>();
    std::vector<int> *NuScoreMatchedY = new std::vector<int>();
    std::vector<int> *NuScoreMatchedZ = new std::vector<int>();
    tree->SetBranchAddress("NuScoreVect", &NuScoreVect);
    tree->SetBranchAddress("NuScoreMatchedOrigin", &NuScoreMatchedOrigin);
    tree->SetBranchAddress("NuScoreMatchedCCNC", &NuScoreMatchedCCNC);
    tree->SetBranchAddress("NuScoreMatchedNuPDG", &NuScoreMatchedNuPDG);
    tree->SetBranchAddress("NuScoreMatchedNuEnergy", &NuScoreMatchedNuEnergy);
    tree->SetBranchAddress("NuScoreMatchedIntMode", &NuScoreMatchedNuIntMode);
    tree->SetBranchAddress("NuScoreMatchedX", &NuScoreMatchedX);
    tree->SetBranchAddress("NuScoreMatchedY", &NuScoreMatchedY);
    tree->SetBranchAddress("NuScoreMatchedZ", &NuScoreMatchedZ);

    // Define histograms signals and backgrounds
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    // Histogram line width
    gStyle->SetHistLineWidth(4);
    // bottom margin
    gStyle->SetPadBottomMargin(0.15);

    // Colors    
    std::vector<int> fColors = {kRed-3,  kAzure-5, kGreen+3, kOrange-3, kMagenta+1, kCyan- 3, kYellow+2, kViolet-1, kTeal-1};
    std::vector<int> fLineStyle = {1, 2, 9, 3, 6, 10, 7, 4, 5, 8};
    std::vector<int> fMarkerStyle = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};

    // TCanvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->Divide(1, 1);
    c1->cd(1);
    gPad->SetLogy();
    double fLabelXPos = 0.3;
    double fLabelYPos = 0.3;
    double fLabelYPosStep = 0.05;

    // THStack
    THStack *hsNuScore = new THStack("hsNuScore", ";#nu score;AU");
    // Legend
    TLegend *legend = new TLegend(0.15, 0.9, 0.85, 1.);
    legend->SetNColumns(fCuts.size());
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    // Label
    TPaveText *labelA = new TPaveText(0.15, 0.7, 0.35, 0.85, "NDC");
    labelA->AddText( "Miss" );
    TPaveText *labelB = new TPaveText(0.4, 0.7, 0.6, 0.85, "NDC");
    labelB->AddText( "Unamb. Cosmic" );
    TPaveText *labelC = new TPaveText(0.65, 0.7, 0.85, 0.85, "NDC");
    labelC->AddText( Form("#nu score > %.1f", fNuScoreCut) );
    // Define the histograms
    float binSize = 0.05;
    //std::vector<float> binEdges = {-2, -1, 0};
    //for(float i = binSize; i <= 1.0+binSize; i += binSize) binEdges.push_back(i);

    std::vector<float> binEdges;
    for(float i = -2.2; i <= 1.0+binSize; i += binSize) binEdges.push_back(i);
    
    size_t nBins = binEdges.size()-1;
    int labelYPos = fLabelYPos;
    std::vector<TH1D*> hNuScoreHistVect;
    std::vector<TLatex*> labelVect;
    size_t index = 0;
    for(auto const& cut : fCuts){
        
        // Declare and fill
        std::string hName = "hNuScore_" + cut.first;
        TH1D *hNuScore = new TH1D(Form(hName.c_str(), cut.first.c_str()), Form("NuScore %s", cut.first.c_str()), nBins, binEdges.data());
        tree->Draw( ("NuScoreVect>>" + hName).c_str(), cut.second );
        
        // Set style
        hNuScore->SetLineColor(fColors[index]);
        hNuScore->SetLineStyle(fLineStyle[index]);
        hNuScore->SetMarkerStyle(fMarkerStyle[index]);

        // Normalise
        hNuScore->Scale(1.0/hNuScore->Integral());

        // Add to stack
        hNuScoreHistVect.push_back(hNuScore);
        hsNuScore->Add( hNuScoreHistVect.back() );

        int nEntries = hNuScore->GetEntries();
        
        // get entries in range -2, -1
        double nEventsMiss = hNuScore->Integral(hNuScore->FindBin(-2), hNuScore->FindBin(-2));
        double fracEventsMiss = 100. * nEventsMiss / hNuScore->Integral();

        // get entries in range -1, 0
        double nEventsUnambCosmic = hNuScore->Integral(hNuScore->FindBin(-1), hNuScore->FindBin(-1));
        double fracEventsUnambCosmic = 100. * nEventsUnambCosmic / hNuScore->Integral();

        // get entries in range fNuScoreCut, 1
        double nEvents = hNuScore->Integral(hNuScore->FindBin(fNuScoreCut), hNuScore->FindBin(1));
        double fracEvents = 100. * nEvents / hNuScore->Integral();
        
        std::cout << "Cut: " << cut.first << std::endl;
        std::cout << "Miss: " << nEventsMiss << " (" << fracEventsMiss << "%)" << std::endl;
        std::cout << "Unamb. Cosmic: " << nEventsUnambCosmic << " (" << fracEventsUnambCosmic << "%)" << std::endl;
        std::cout << "NuScore > " << fNuScoreCut << ": " << nEvents << " (" << fracEvents << "%)" << std::endl;

        // Add to legend
        legend->AddEntry(hNuScoreHistVect.back(), cut.first.c_str(), "l");
        if(cut.first != "Cosmic")
            labelA->AddText(Form("%s: %.1f%%", cut.first.c_str(), fracEventsMiss));
        labelB->AddText(Form("%s: %.1f%%", cut.first.c_str(), fracEventsUnambCosmic));
        labelC->AddText(Form("%s: %.1f%%", cut.first.c_str(), fracEvents));
       
        index++;
    }


    hsNuScore->Draw("nostack hist");

    // NDivisions
    hsNuScore->GetXaxis()->SetNdivisions(505);
    /*// NDivisions
    hsNuScore->GetXaxis()->SetBinLabel(1,"Miss");
    hsNuScore->GetXaxis()->SetBinLabel(2,"Unamb. Cosmic");
    hsNuScore->GetXaxis()->SetBinLabel(3,"0");
    //hsNuScore->GetXaxis()->SetBinLabel(3 + (nBins-3)/2,"0.5");
    hsNuScore->GetXaxis()->SetBinLabel(nBins,"1");*/
    
    legend->Draw("same");
    labelA->Draw("same");
    labelB->Draw("same");
    labelC->Draw("same");

    c1->cd();
    c1->Update();
    c1->WaitPrimitive();
    c1->Update();
    c1->SaveAs("NuScore.pdf");


}
