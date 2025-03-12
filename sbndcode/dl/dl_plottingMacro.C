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
    double deltaX;
    double deltaY;
    double deltaZ;
    double deltaR;
    int CCNC;
    int trueNeutrino; // Neutrino type
    int DLCurrent;
    int CosmicSingle;
    int Pass1Pass2;
    unsigned int eventID;
    unsigned int runID;
    unsigned int subRunID;
    int fileNumber; // The file number that this event is in (i.e 1 would be 1.root)
    double tpcID;
    double uncorrectedDeltaX;
    double uncorrectedDeltaR;
    double recoX;
    double recoY;
    double recoZ;
    double trueX;
    double trueY;
    double trueZ;
} vertex_t;

void bigDifferences(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();


    std::vector<vertex_t> dlDuneVerticesPass1 = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVerticesPass1 = std::vector<vertex_t>();
    std::vector<vertex_t> currentVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlDuneVerticesPass2 = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVerticesPass2 = std::vector<vertex_t>();

    double largestValueROC = 0;

    for(UInt_t i = 0; i < numVertices; i++){
        vertex_t vertex = allVertices_vec.at(i);
        if(vertex.Pass1Pass2 == 1){
            if(vertex.DLCurrent == 0){
                dlUbooneVerticesPass1.push_back(vertex);
            } else if(vertex.DLCurrent == 1){
                dlDuneVerticesPass1.push_back(vertex);
            } else if(vertex.DLCurrent == 2){
                currentVertices.push_back(vertex);
            }
        } else if(vertex.Pass1Pass2 == 2){
            if(vertex.DLCurrent == 0){
                dlUbooneVerticesPass2.push_back(vertex);
            } else if(vertex.DLCurrent == 1){
                dlDuneVerticesPass2.push_back(vertex);
            } else if(vertex.DLCurrent == 2){
                currentVertices.push_back(vertex);
            }
        }
    }

}

void pass1pass2(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();

    std::vector<vertex_t> dlDuneVerticesPass1 = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVerticesPass1 = std::vector<vertex_t>();
    std::vector<vertex_t> currentVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlDuneVerticesPass2 = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVerticesPass2 = std::vector<vertex_t>();

    double largestValueROC = 0;

    for(UInt_t i = 0; i < numVertices; i++){
        vertex_t vertex = allVertices_vec.at(i);
        if(vertex.Pass1Pass2 == 1){
            if(vertex.DLCurrent == 0){
                dlUbooneVerticesPass1.push_back(vertex);
            } else if(vertex.DLCurrent == 1){
                dlDuneVerticesPass1.push_back(vertex);
            } else if(vertex.DLCurrent == 2){
                currentVertices.push_back(vertex);
            }
        } else if(vertex.Pass1Pass2 == 2){
            if(vertex.DLCurrent == 0){
                dlUbooneVerticesPass2.push_back(vertex);
            } else if(vertex.DLCurrent == 1){
                dlDuneVerticesPass2.push_back(vertex);
            } else if(vertex.DLCurrent == 2){
                currentVertices.push_back(vertex);
            }
        }

        if(largestValueROC < vertex.deltaR){
            largestValueROC = vertex.deltaR;
        }
    }

    UInt_t numDLDuneVerticesPass1 = dlDuneVerticesPass1.size();
    UInt_t numDLUbooneVerticesPass1 = dlUbooneVerticesPass1.size();
    UInt_t numCurrentVertices = currentVertices.size();
    UInt_t numDLDuneVerticesPass2 = dlDuneVerticesPass2.size();
    UInt_t numDLUbooneVerticesPass2 = dlUbooneVerticesPass2.size();

    double binWidth = 0.47619048;
    double lowerBoundX = -5; // std::floor(smallestValueX); 
    double upperBoundX = 5;  //std::ceil(largestValueX);
    double lowerBoundY = -5; //std::floor(smallestValueY);
    double upperBoundY = 5;//std::ceil(largestValueY);
    double lowerBoundZ = -5;//std::floor(smallestValueZ);
    double upperBoundZ = 5;//std::ceil(largestValueZ);
    double upperBoundR = 5;//std::ceil(largestValueR);
    //double lowerBound = -4;
    //double upperBound = 4;
    double numOfBinsX = std::ceil((upperBoundX - lowerBoundX)/binWidth);
    double numOfBinsY = std::ceil((upperBoundY - lowerBoundY)/binWidth);
    double numOfBinsZ = std::ceil((upperBoundZ - lowerBoundZ)/binWidth);
    double numOfBinsR = std::ceil((upperBoundR)/binWidth);

    // Pass 1
    TH1F *Current_deltaX_freq_pass1 = new TH1F("deltaX Current freq", "Current x Coord Diff freq", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *DL_uboone_deltaX_freq_pass1 = (TH1F*) Current_deltaX_freq_pass1->Clone("deltaX DL uboone freq");
    TH1F *DL_dune_deltaX_freq_pass1 = (TH1F*) Current_deltaX_freq_pass1->Clone("deltaX DL dune freq");

    TH1F *Current_deltaY_freq_pass1 = new TH1F("deltaY DL dune freq", "DL dune Y Coord Diff freq", numOfBinsY, lowerBoundY, upperBoundY);
    TH1F *DL_uboone_deltaY_freq_pass1 = (TH1F*) Current_deltaY_freq_pass1->Clone("deltaY DL uboone freq");
    TH1F *DL_dune_deltaY_freq_pass1 = (TH1F*) Current_deltaY_freq_pass1->Clone("deltaY Current freq");

    TH1F *Current_deltaZ_freq_pass1 = new TH1F("deltaZ Current freq", "Current Z Coord Diff freq", numOfBinsZ, lowerBoundZ, upperBoundZ);
    TH1F *DL_uboone_deltaZ_freq_pass1 = (TH1F*) Current_deltaZ_freq_pass1->Clone("deltaZ DL uboone freq");
    TH1F *DL_dune_deltaZ_freq_pass1 = (TH1F*) Current_deltaZ_freq_pass1->Clone("deltaZ DL dune freq");

    TH1F *Current_deltaR_freq_pass1 = new TH1F("deltaR Current freq", "Current R Coord Diff freq", numOfBinsR, 0, upperBoundR);
    TH1F *DL_uboone_deltaR_freq_pass1 = (TH1F*) Current_deltaR_freq_pass1->Clone("deltaR DL uboone freq");
    TH1F *DL_dune_deltaR_freq_pass1 = (TH1F*) Current_deltaR_freq_pass1->Clone("deltaR DL dune freq");


    // Pass 2
    TH1F *Current_deltaX_freq_pass2 = new TH1F("deltaX Current freq", "Current x Coord Diff freq", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *DL_uboone_deltaX_freq_pass2 = (TH1F*) Current_deltaX_freq_pass2->Clone("deltaX DL uboone freq");
    TH1F *DL_dune_deltaX_freq_pass2 = (TH1F*) Current_deltaX_freq_pass2->Clone("deltaX DL dune freq");

    TH1F *Current_deltaY_freq_pass2 = new TH1F("deltaY DL dune freq", "DL dune Y Coord Diff freq", numOfBinsY, lowerBoundY, upperBoundY);
    TH1F *DL_uboone_deltaY_freq_pass2 = (TH1F*) Current_deltaY_freq_pass2->Clone("deltaY DL uboone freq");
    TH1F *DL_dune_deltaY_freq_pass2 = (TH1F*) Current_deltaY_freq_pass2->Clone("deltaY Current freq");

    TH1F *Current_deltaZ_freq_pass2 = new TH1F("deltaZ Current freq", "Current Z Coord Diff freq", numOfBinsZ, lowerBoundZ, upperBoundZ);
    TH1F *DL_uboone_deltaZ_freq_pass2 = (TH1F*) Current_deltaZ_freq_pass2->Clone("deltaZ DL uboone freq");
    TH1F *DL_dune_deltaZ_freq_pass2 = (TH1F*) Current_deltaZ_freq_pass2->Clone("deltaZ DL dune freq");

    TH1F *Current_deltaR_freq_pass2 = new TH1F("deltaR Current freq", "Current R Coord Diff freq", numOfBinsR, 0, upperBoundR);
    TH1F *DL_uboone_deltaR_freq_pass2 = (TH1F*) Current_deltaR_freq_pass2->Clone("deltaR DL uboone freq");
    TH1F *DL_dune_deltaR_freq_pass2 = (TH1F*) Current_deltaR_freq_pass2->Clone("deltaR DL dune freq");

    std::cout << "DL DUNE ---- Number of pass 2: " << numDLDuneVerticesPass2 << " Number of pass 1: " << numDLDuneVerticesPass1 << std::endl;
    std::cout << "DL UBoone ---- Number of pass 2: " << numDLUbooneVerticesPass2 << " Number of pass 1: " << numDLUbooneVerticesPass1 << std::endl;
    std::cout << "Current ---- Number of vertices: " << numCurrentVertices << std::endl;

    for(UInt_t j = 0; j < numDLDuneVerticesPass1; j++){
        vertex_t dlDuneVertex = dlDuneVerticesPass1.at(j);
        DL_dune_deltaX_freq_pass1->Fill(dlDuneVertex.deltaX);
        DL_dune_deltaY_freq_pass1->Fill(dlDuneVertex.deltaY);
        DL_dune_deltaZ_freq_pass1->Fill(dlDuneVertex.deltaZ);
        DL_dune_deltaR_freq_pass1->Fill(dlDuneVertex.deltaR);
    }

    for(UInt_t j = 0; j < numDLDuneVerticesPass2; j++){
        vertex_t dlDuneVertex = dlDuneVerticesPass2.at(j);
        DL_dune_deltaX_freq_pass2->Fill(dlDuneVertex.deltaX);
        DL_dune_deltaY_freq_pass2->Fill(dlDuneVertex.deltaY);
        DL_dune_deltaZ_freq_pass2->Fill(dlDuneVertex.deltaZ);
        DL_dune_deltaR_freq_pass2->Fill(dlDuneVertex.deltaR);
    }

    for(UInt_t j = 0; j < numDLUbooneVerticesPass1; j++){
        vertex_t dlUbooneVertex = dlUbooneVerticesPass1.at(j);
        DL_uboone_deltaX_freq_pass1->Fill(dlUbooneVertex.deltaX);
        DL_uboone_deltaY_freq_pass1->Fill(dlUbooneVertex.deltaY);
        DL_uboone_deltaZ_freq_pass1->Fill(dlUbooneVertex.deltaZ);
        DL_uboone_deltaR_freq_pass1->Fill(dlUbooneVertex.deltaR);
    }

    for(UInt_t j = 0; j < numDLUbooneVerticesPass2; j++){
        vertex_t dlUbooneVertex = dlUbooneVerticesPass2.at(j);
        DL_uboone_deltaX_freq_pass2->Fill(dlUbooneVertex.deltaX);
        DL_uboone_deltaY_freq_pass2->Fill(dlUbooneVertex.deltaY);
        DL_uboone_deltaZ_freq_pass2->Fill(dlUbooneVertex.deltaZ);
        DL_uboone_deltaR_freq_pass2->Fill(dlUbooneVertex.deltaR);
    }

    for(UInt_t j = 0; j < numCurrentVertices; j++){
        vertex_t currentVertex = currentVertices.at(j);
        Current_deltaX_freq_pass1->Fill(currentVertex.deltaX);
        Current_deltaY_freq_pass1->Fill(currentVertex.deltaY);
        Current_deltaZ_freq_pass1->Fill(currentVertex.deltaZ);
        Current_deltaR_freq_pass1->Fill(currentVertex.deltaR);
        Current_deltaX_freq_pass2->Fill(currentVertex.deltaX);
        Current_deltaY_freq_pass2->Fill(currentVertex.deltaY);
        Current_deltaZ_freq_pass2->Fill(currentVertex.deltaZ);
        Current_deltaR_freq_pass2->Fill(currentVertex.deltaR);
    }

    TCanvas *Canvas100_pass1 = new TCanvas("c100_pass1", "Graph Draw Options", 200, 10, 600, 400);

    // Pass 1
    TH1F *Current_deltaX_perc_pass1 = new TH1F("deltaX DL dune %", "x Coord Diff %", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *DL_uboone_deltaX_perc_pass1 = (TH1F*) Current_deltaX_perc_pass1->Clone("deltaX DL uboone %");
    TH1F *DL_dune_deltaX_perc_pass1 = (TH1F*) Current_deltaX_perc_pass1->Clone("deltaX Current %");

    TH1F *Current_deltaY_perc_pass1 = new TH1F("deltaY DL dune %", "Y Coord Diff %", numOfBinsY, lowerBoundY, upperBoundY);
    TH1F *DL_uboone_deltaY_perc_pass1 = (TH1F*) Current_deltaY_perc_pass1->Clone("deltaY DL uboone %");
    TH1F *DL_dune_deltaY_perc_pass1 = (TH1F*) Current_deltaY_perc_pass1->Clone("deltaY Current %");

    TH1F *Current_deltaZ_perc_pass1 = new TH1F("deltaZ DL dune %", "Z Coord Diff %", numOfBinsZ, lowerBoundZ, upperBoundZ);
    TH1F *DL_uboone_deltaZ_perc_pass1 = (TH1F*) Current_deltaZ_perc_pass1->Clone("deltaZ DL uboone %");
    TH1F *DL_dune_deltaZ_perc_pass1 = (TH1F*) Current_deltaZ_perc_pass1->Clone("deltaZ Current %");

    TH1F *Current_deltaR_perc_pass1 = new TH1F("deltaR DL dune %", "R Coord Diff %", numOfBinsR, 0, upperBoundR);
    TH1F *DL_uboone_deltaR_perc_pass1 = (TH1F*) Current_deltaR_perc_pass1->Clone("deltaR DL uboone %");
    TH1F *DL_dune_deltaR_perc_pass1 = (TH1F*) Current_deltaR_perc_pass1->Clone("deltaR Current %");


    // Pass 2
    TH1F *Current_deltaX_perc_pass2 = new TH1F("deltaX DL dune %", "x Coord Diff %", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *DL_uboone_deltaX_perc_pass2 = (TH1F*) Current_deltaX_perc_pass2->Clone("deltaX DL uboone %");
    TH1F *DL_dune_deltaX_perc_pass2 = (TH1F*) Current_deltaX_perc_pass2->Clone("deltaX Current %");

    TH1F *Current_deltaY_perc_pass2 = new TH1F("deltaY DL dune %", "Y Coord Diff %", numOfBinsY, lowerBoundY, upperBoundY);
    TH1F *DL_uboone_deltaY_perc_pass2 = (TH1F*) Current_deltaY_perc_pass2->Clone("deltaY DL uboone %");
    TH1F *DL_dune_deltaY_perc_pass2 = (TH1F*) Current_deltaY_perc_pass2->Clone("deltaY Current %");

    TH1F *Current_deltaZ_perc_pass2 = new TH1F("deltaZ DL dune %", "Z Coord Diff %", numOfBinsZ, lowerBoundZ, upperBoundZ);
    TH1F *DL_uboone_deltaZ_perc_pass2 = (TH1F*) Current_deltaZ_perc_pass2->Clone("deltaZ DL uboone %");
    TH1F *DL_dune_deltaZ_perc_pass2 = (TH1F*) Current_deltaZ_perc_pass2->Clone("deltaZ Current %");

    TH1F *Current_deltaR_perc_pass2 = new TH1F("deltaR DL dune %", "R Coord Diff %", numOfBinsR, 0, upperBoundR);
    TH1F *DL_uboone_deltaR_perc_pass2 = (TH1F*) Current_deltaR_perc_pass2->Clone("deltaR DL uboone %");
    TH1F *DL_dune_deltaR_perc_pass2 = (TH1F*) Current_deltaR_perc_pass2->Clone("deltaR Current %");

    double binFillingX = lowerBoundX + (binWidth/2);
    double binFillingY = lowerBoundY + (binWidth/2);
    double binFillingZ = lowerBoundZ + (binWidth/2);
    double binFillingR = 0 + (binWidth/2);

    for(int i = 1; i < numOfBinsX+1; i++){
        DL_dune_deltaX_perc_pass1->Fill(binFillingX, (100.0f*DL_dune_deltaX_freq_pass1->GetBinContent(i)/(double)numDLDuneVerticesPass1));
        DL_uboone_deltaX_perc_pass1->Fill(binFillingX, (100.0f*DL_uboone_deltaX_freq_pass1->GetBinContent(i)/(double)numDLUbooneVerticesPass1));
        Current_deltaX_perc_pass1->Fill(binFillingX, (100.0f*Current_deltaX_freq_pass1->GetBinContent(i)/(double)numCurrentVertices));
        DL_dune_deltaX_perc_pass2->Fill(binFillingX, (100.0f*DL_dune_deltaX_freq_pass2->GetBinContent(i)/(double)numDLDuneVerticesPass2));
        DL_uboone_deltaX_perc_pass2->Fill(binFillingX, (100.0f*DL_uboone_deltaX_freq_pass2->GetBinContent(i)/(double)numDLUbooneVerticesPass2));
        Current_deltaX_perc_pass2->Fill(binFillingX, (100.0f*Current_deltaX_freq_pass2->GetBinContent(i)/(double)numCurrentVertices));
        binFillingX += binWidth;
    }

    for(int i = 1; i < numOfBinsY+1; i++){
        DL_dune_deltaY_perc_pass1->Fill(binFillingY, (100.0f*DL_dune_deltaY_freq_pass1->GetBinContent(i)/(double)numDLDuneVerticesPass1));
        DL_uboone_deltaY_perc_pass1->Fill(binFillingY, (100.0f*DL_uboone_deltaY_freq_pass1->GetBinContent(i)/(double)numDLUbooneVerticesPass1));
        Current_deltaY_perc_pass1->Fill(binFillingY, (100.0f*Current_deltaY_freq_pass1->GetBinContent(i)/(double)numCurrentVertices));
        DL_dune_deltaY_perc_pass2->Fill(binFillingY, (100.0f*DL_dune_deltaY_freq_pass2->GetBinContent(i)/(double)numDLDuneVerticesPass2));
        DL_uboone_deltaY_perc_pass2->Fill(binFillingY, (100.0f*DL_uboone_deltaY_freq_pass2->GetBinContent(i)/(double)numDLUbooneVerticesPass2));
        Current_deltaY_perc_pass2->Fill(binFillingY, (100.0f*Current_deltaY_freq_pass2->GetBinContent(i)/(double)numCurrentVertices));
        binFillingY += binWidth;
    }

    for(int i = 1; i < numOfBinsZ+1; i++){
        DL_dune_deltaZ_perc_pass1->Fill(binFillingZ, (100.0f*DL_dune_deltaZ_freq_pass1->GetBinContent(i)/(double)numDLDuneVerticesPass1));
        DL_uboone_deltaZ_perc_pass1->Fill(binFillingZ, (100.0f*DL_uboone_deltaZ_freq_pass1->GetBinContent(i)/(double)numDLUbooneVerticesPass1));
        Current_deltaZ_perc_pass1->Fill(binFillingZ, (100.0f*Current_deltaZ_freq_pass1->GetBinContent(i)/(double)numCurrentVertices));
        DL_dune_deltaZ_perc_pass2->Fill(binFillingZ, (100.0f*DL_dune_deltaZ_freq_pass2->GetBinContent(i)/(double)numDLDuneVerticesPass2));
        DL_uboone_deltaZ_perc_pass2->Fill(binFillingZ, (100.0f*DL_uboone_deltaZ_freq_pass2->GetBinContent(i)/(double)numDLUbooneVerticesPass2));
        Current_deltaZ_perc_pass2->Fill(binFillingZ, (100.0f*Current_deltaZ_freq_pass2->GetBinContent(i)/(double)numCurrentVertices));
        binFillingZ += binWidth;
    }

    for(int i = 1; i < numOfBinsR+1; i++){
        DL_dune_deltaR_perc_pass1->Fill(binFillingR, (100.0f*DL_dune_deltaR_freq_pass1->GetBinContent(i)/(double)numDLDuneVerticesPass1));
        DL_uboone_deltaR_perc_pass1->Fill(binFillingR, (100.0f*DL_uboone_deltaR_freq_pass1->GetBinContent(i)/(double)numDLUbooneVerticesPass1));
        Current_deltaR_perc_pass1->Fill(binFillingR, (100.0f*Current_deltaR_freq_pass1->GetBinContent(i)/(double)numCurrentVertices));
        DL_dune_deltaR_perc_pass2->Fill(binFillingR, (100.0f*DL_dune_deltaR_freq_pass2->GetBinContent(i)/(double)numDLDuneVerticesPass2));
        DL_uboone_deltaR_perc_pass2->Fill(binFillingR, (100.0f*DL_uboone_deltaR_freq_pass2->GetBinContent(i)/(double)numDLUbooneVerticesPass2));
        Current_deltaR_perc_pass2->Fill(binFillingR, (100.0f*Current_deltaR_freq_pass2->GetBinContent(i)/(double)numCurrentVertices));
        binFillingR += binWidth;
    }

    DL_dune_deltaX_perc_pass1->SetLineColor(kViolet-5);
    DL_dune_deltaX_perc_pass1->SetLineWidth(2);
    DL_uboone_deltaX_perc_pass1->SetLineColor(kBlue);
    DL_uboone_deltaX_perc_pass1->SetLineWidth(2);
    Current_deltaX_perc_pass1->SetLineColor(kRed);
    Current_deltaX_perc_pass1->SetLineWidth(2);

    DL_dune_deltaX_perc_pass2->SetLineColor(kViolet-5);
    DL_dune_deltaX_perc_pass2->SetLineWidth(2);
    DL_uboone_deltaX_perc_pass2->SetLineColor(kBlue);
    DL_uboone_deltaX_perc_pass2->SetLineWidth(2);
    Current_deltaX_perc_pass2->SetLineColor(kRed);
    Current_deltaX_perc_pass2->SetLineWidth(2);

    Current_deltaX_perc_pass1->Draw("hist");
    DL_uboone_deltaX_perc_pass1->Draw("histsame");
    DL_dune_deltaX_perc_pass1->Draw("histsame");
    Current_deltaX_perc_pass1->SetStats(0);
    Current_deltaX_perc_pass1->SetTitle("DeltaX Pass 1");
    Current_deltaX_perc_pass1->GetXaxis()->SetTitle("x_{Reco} - x_{True} (cm)");
    Current_deltaX_perc_pass1->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaX_perc_pass1->GetYaxis()->SetRange(22,23);

    auto legend_deltaX_pass1 = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend_deltaX_pass1->AddEntry(DL_dune_deltaX_perc_pass1, "Deep Learning: DUNE/LBNF Tune", "f");
    legend_deltaX_pass1->AddEntry(DL_uboone_deltaX_perc_pass1, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend_deltaX_pass1->AddEntry(Current_deltaX_perc_pass1, "Current SBND Vertexing", "f");
    legend_deltaX_pass1->Draw();

    TCanvas *Canvas100_pass2 = new TCanvas("c100_pass2", "Graph Draw Options", 200, 10, 600, 400);
    Current_deltaX_perc_pass2->Draw("hist");
    DL_uboone_deltaX_perc_pass2->Draw("histsame");
    DL_dune_deltaX_perc_pass2->Draw("histsame");
    Current_deltaX_perc_pass2->SetStats(0);
    Current_deltaX_perc_pass2->SetTitle("DeltaX Pass 2");
    Current_deltaX_perc_pass2->GetXaxis()->SetTitle("x_{Reco} - x_{True} (cm)");
    Current_deltaX_perc_pass2->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaX_perc_pass2->GetYaxis()->SetRange(22,23);

    std::cout << "deltaX Pass 1" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaX_perc_pass1->GetMean() << " RMS: " << Current_deltaX_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaX_perc_pass1->GetMean() << " RMS: " << DL_dune_deltaX_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaX_perc_pass1->GetMean() << "RMS: " << DL_uboone_deltaX_perc_pass1->GetRMS(1) << std::endl;

    auto legend_pass1 = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend_pass1->AddEntry(DL_dune_deltaX_perc_pass1, "Deep Learning: DUNE/LBNF Tune", "f");
    legend_pass1->AddEntry(DL_uboone_deltaX_perc_pass1, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend_pass1->AddEntry(Current_deltaX_perc_pass1, "Current SBND Vertexing", "f");
    legend_pass1->Draw();

    std::cout << "deltaX Pass 2" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaX_perc_pass2->GetMean() << " RMS: " << Current_deltaX_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaX_perc_pass2->GetMean() << " RMS: " << DL_dune_deltaX_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaX_perc_pass2->GetMean() << "RMS: " << DL_uboone_deltaX_perc_pass2->GetRMS(1) << std::endl;

    auto legend_pass2 = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend_pass2->AddEntry(DL_dune_deltaX_perc_pass2, "Deep Learning: DUNE/LBNF Tune", "f");
    legend_pass2->AddEntry(DL_uboone_deltaX_perc_pass2, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend_pass2->AddEntry(Current_deltaX_perc_pass2, "Current SBND Vertexing", "f");
    legend_pass2->Draw();

    TCanvas *Canvas1_pass1 = new TCanvas("c1_pass1", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaY_perc_pass1->SetLineColor(kViolet-5);
    DL_dune_deltaY_perc_pass1->SetLineWidth(2);
    DL_uboone_deltaY_perc_pass1->SetLineColor(kBlue);
    DL_uboone_deltaY_perc_pass1->SetLineWidth(2);
    Current_deltaY_perc_pass1->SetLineColor(kRed);
    Current_deltaY_perc_pass1->SetLineWidth(2);
    Current_deltaY_perc_pass1->Draw("hist");
    DL_uboone_deltaY_perc_pass1->Draw("histsame");
    DL_dune_deltaY_perc_pass1->Draw("histsame");

    Current_deltaY_perc_pass1->SetStats(0);
    Current_deltaY_perc_pass1->SetTitle("DeltaY Pass 1");
    Current_deltaY_perc_pass1->GetXaxis()->SetTitle("#Delta Y (cm)");
    Current_deltaY_perc_pass1->GetYaxis()->SetTitle("Percentage of Events (%)");

    auto legend2_pass1 = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend2_pass1->AddEntry(DL_dune_deltaY_perc_pass1, "Deep Learning: DUNE/LBNF Tune ", "f");
    legend2_pass1->AddEntry(DL_uboone_deltaY_perc_pass1, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend2_pass1->AddEntry(Current_deltaY_perc_pass1, "Current SBND Vertexing", "f");
    legend2_pass1->Draw();

    std::cout << "deltaY Pass 1" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaY_perc_pass1->GetMean() << " RMS: " << Current_deltaY_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaY_perc_pass1->GetMean() << " RMS: " << DL_dune_deltaY_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaY_perc_pass1->GetMean() << "RMS: " << DL_uboone_deltaY_perc_pass1->GetRMS(1) << std::endl;

    TCanvas *Canvas1_pass2 = new TCanvas("c1_pass2", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaY_perc_pass2->SetLineColor(kViolet-5);
    DL_dune_deltaY_perc_pass2->SetLineWidth(2);
    DL_uboone_deltaY_perc_pass2->SetLineColor(kBlue);
    DL_uboone_deltaY_perc_pass2->SetLineWidth(2);
    Current_deltaY_perc_pass2->SetLineColor(kRed);
    Current_deltaY_perc_pass2->SetLineWidth(2);
    Current_deltaY_perc_pass2->Draw("hist");
    DL_uboone_deltaY_perc_pass2->Draw("histsame");
    DL_dune_deltaY_perc_pass2->Draw("histsame");

    Current_deltaY_perc_pass2->SetStats(0);
    Current_deltaY_perc_pass2->SetTitle("DeltaY Pass 2");
    Current_deltaY_perc_pass2->GetXaxis()->SetTitle("#Delta Y (cm)");
    Current_deltaY_perc_pass2->GetYaxis()->SetTitle("Percentage of Events (%)");

    auto legend2_pass2 = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend2_pass2->AddEntry(DL_dune_deltaY_perc_pass2, "Deep Learning: DUNE/LBNF Tune", "f");
    legend2_pass2->AddEntry(DL_uboone_deltaY_perc_pass2, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend2_pass2->AddEntry(Current_deltaY_perc_pass2, "Current SBND Vertexing", "f");
    legend2_pass2->Draw();

    std::cout << "deltaY Pass 2" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaY_perc_pass2->GetMean() << " RMS: " << Current_deltaY_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaY_perc_pass2->GetMean() << " RMS: " << DL_dune_deltaY_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaY_perc_pass2->GetMean() << "RMS: " << DL_uboone_deltaY_perc_pass2->GetRMS(1) << std::endl;

    TCanvas *Canvas3_pass1 = new TCanvas("c3_pass1", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaZ_perc_pass1->SetLineColor(kViolet-5);
    DL_dune_deltaZ_perc_pass1->SetLineWidth(2);
    DL_uboone_deltaZ_perc_pass1->SetLineColor(kBlue);
    DL_uboone_deltaZ_perc_pass1->SetLineWidth(2);
    Current_deltaZ_perc_pass1->SetLineColor(kRed);
    Current_deltaZ_perc_pass1->SetLineWidth(2);
    Current_deltaZ_perc_pass1->Draw("hist");
    DL_uboone_deltaZ_perc_pass1->Draw("histsame");
    DL_dune_deltaZ_perc_pass1->Draw("histsame");

    Current_deltaZ_perc_pass1->SetStats(0);
    Current_deltaZ_perc_pass1->SetTitle("DeltaZ Pass 1");
    Current_deltaZ_perc_pass1->GetXaxis()->SetTitle("#Delta Z (cm)");
    Current_deltaZ_perc_pass1->GetYaxis()->SetTitle("Percentage of Events (%)");

    auto legend3_pass1 = new TLegend(0.1,0.7,0.48,0.9);
    legend3_pass1->AddEntry(DL_dune_deltaZ_perc_pass1, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3_pass1->AddEntry(DL_uboone_deltaZ_perc_pass1, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3_pass1->AddEntry(Current_deltaZ_perc_pass1, "Current SBND Vertexing", "f");
    legend3_pass1->Draw();

    std::cout << "deltaZ Pass 1" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaZ_perc_pass1->GetMean() << " RMS: " << Current_deltaZ_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaZ_perc_pass1->GetMean() << " RMS: " << DL_dune_deltaZ_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaZ_perc_pass1->GetMean() << "RMS: " << DL_uboone_deltaZ_perc_pass1->GetRMS(1) << std::endl;

    TCanvas *Canvas3_pass2 = new TCanvas("c3_pass2", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaZ_perc_pass2->SetLineColor(kViolet-5);
    DL_dune_deltaZ_perc_pass2->SetLineWidth(2);
    DL_uboone_deltaZ_perc_pass2->SetLineColor(kBlue);
    DL_uboone_deltaZ_perc_pass2->SetLineWidth(2);
    Current_deltaZ_perc_pass2->SetLineColor(kRed);
    Current_deltaZ_perc_pass2->SetLineWidth(2);
    Current_deltaZ_perc_pass2->Draw("hist");
    DL_uboone_deltaZ_perc_pass2->Draw("histsame");
    DL_dune_deltaZ_perc_pass2->Draw("histsame");

    Current_deltaZ_perc_pass2->SetStats(0);
    Current_deltaZ_perc_pass2->SetTitle("DeltaZ Pass 2");
    Current_deltaZ_perc_pass2->GetXaxis()->SetTitle("#Delta Z (cm)");
    Current_deltaZ_perc_pass2->GetYaxis()->SetTitle("Percentage of Events (%)");

    auto legend3_pass2 = new TLegend(0.1,0.7,0.48,0.9);
    legend3_pass2->AddEntry(DL_dune_deltaZ_perc_pass2, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3_pass2->AddEntry(DL_uboone_deltaZ_perc_pass2, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3_pass2->AddEntry(Current_deltaZ_perc_pass2, "Current SBND Vertexing", "f");
    legend3_pass2->Draw();

    std::cout << "deltaZ Pass 2" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaZ_perc_pass2->GetMean() << " RMS: " << Current_deltaZ_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaZ_perc_pass2->GetMean() << " RMS: " << DL_dune_deltaZ_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaZ_perc_pass2->GetMean() << "RMS: " << DL_uboone_deltaZ_perc_pass2->GetRMS(1) << std::endl;

    TCanvas *Canvas4_pass1 = new TCanvas("c4_pass1", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaR_perc_pass1->SetLineColor(kViolet-5);
    DL_dune_deltaR_perc_pass1->SetLineWidth(2);
    DL_uboone_deltaR_perc_pass1->SetLineColor(kBlue);
    DL_uboone_deltaR_perc_pass1->SetLineWidth(2);
    Current_deltaR_perc_pass1->SetLineColor(kRed);
    Current_deltaR_perc_pass1->SetLineWidth(2);
    Current_deltaR_perc_pass1->Draw("hist");
    DL_uboone_deltaR_perc_pass1->Draw("histsame");
    DL_dune_deltaR_perc_pass1->Draw("histsame"); 

    Current_deltaR_perc_pass1->SetStats(0);
    Current_deltaR_perc_pass1->SetTitle("DeltaR Pass 1");
    Current_deltaR_perc_pass1->GetXaxis()->SetTitle("|#bar{r}_{Reco} - #bar{r}_{True}| (cm)");
    Current_deltaR_perc_pass1->GetYaxis()->SetTitle("Percentage of Events (%)");

    auto legend4_pass1 = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend4_pass1->AddEntry(DL_dune_deltaR_perc_pass1, "Deep Learning: DUNE/LBNF Tune", "f");
    legend4_pass1->AddEntry(DL_uboone_deltaR_perc_pass1, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend4_pass1->AddEntry(Current_deltaR_perc_pass1, "Current SBND Vertexing", "f");
    legend4_pass1->Draw();

    std::cout << "deltaR Pass 1" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaR_perc_pass1->GetMean() << " RMS: " << Current_deltaR_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaR_perc_pass1->GetMean() << " RMS: " << DL_dune_deltaR_perc_pass1->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaR_perc_pass1->GetMean() << "RMS: " << DL_uboone_deltaR_perc_pass1->GetRMS(1) << std::endl;

    TCanvas *Canvas4_pass2 = new TCanvas("c4_pass2", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaR_perc_pass2->SetLineColor(kViolet-5);
    DL_dune_deltaR_perc_pass2->SetLineWidth(2);
    DL_uboone_deltaR_perc_pass2->SetLineColor(kBlue);
    DL_uboone_deltaR_perc_pass2->SetLineWidth(2);
    Current_deltaR_perc_pass2->SetLineColor(kRed);
    Current_deltaR_perc_pass2->SetLineWidth(2);
    Current_deltaR_perc_pass2->Draw("hist");
    DL_uboone_deltaR_perc_pass2->Draw("histsame");
    DL_dune_deltaR_perc_pass2->Draw("histsame"); 

    Current_deltaR_perc_pass2->SetStats(0);
    Current_deltaR_perc_pass2->SetTitle("DeltaR Pass 2");
    Current_deltaR_perc_pass2->GetXaxis()->SetTitle("|#bar{r}_{Reco} - #bar{r}_{True}| (cm)");
    Current_deltaR_perc_pass2->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaR_perc_pass2->GetYaxis()->SetRangeUser(0, 28);

    auto legend4_pass2 = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend4_pass2->AddEntry(DL_dune_deltaR_perc_pass2, "Deep Learning: DUNE/LBNF Tune", "f");
    legend4_pass2->AddEntry(DL_uboone_deltaR_perc_pass2, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend4_pass2->AddEntry(Current_deltaR_perc_pass2, "Current SBND Vertexing", "f");
    legend4_pass2->Draw();

    std::cout << "deltaR Pass 2" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaR_perc_pass2->GetMean() << " RMS: " << Current_deltaR_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaR_perc_pass2->GetMean() << " RMS: " << DL_dune_deltaR_perc_pass2->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaR_perc_pass2->GetMean() << "RMS: " << DL_uboone_deltaR_perc_pass2->GetRMS(1) << std::endl;


    TH1F *deltaXUboone_pass2 = new TH1F("deltaX uboone % pass 2", "deltaX uboone % pass 2", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *deltaXUboone_pass1 = (TH1F*) deltaXUboone_pass2->Clone("deltaX uboone % pass 1");
    TH1F *deltaXDune_pass2 = (TH1F*) deltaXUboone_pass2->Clone("deltaX dune % pass 2");
    TH1F *deltaXDune_pass1 = (TH1F*) deltaXUboone_pass2->Clone("deltaX dune % pass 1");

    double binFillingX_pass1pass2 = lowerBoundX + (binWidth/2);
    
    for(int i = 1; i < numOfBinsX+1; i++){
        deltaXUboone_pass2->Fill(binFillingX_pass1pass2, DL_uboone_deltaX_perc_pass2->GetBinContent(i));
        deltaXUboone_pass1->Fill(binFillingX_pass1pass2, DL_uboone_deltaX_perc_pass1->GetBinContent(i));
        deltaXDune_pass2->Fill(binFillingX_pass1pass2, DL_dune_deltaX_perc_pass2->GetBinContent(i));
        deltaXDune_pass1->Fill(binFillingX_pass1pass2, DL_dune_deltaX_perc_pass1->GetBinContent(i));
        binFillingX_pass1pass2 += binWidth;
    }

    deltaXUboone_pass2->SetLineColor(kGreen+3);
    deltaXUboone_pass1->SetLineColor(kGreen-3);
    deltaXDune_pass2->SetLineColor(kBlue+1);
    deltaXDune_pass1->SetLineColor(kAzure+2);
    deltaXUboone_pass2->SetLineWidth(2);
    deltaXUboone_pass1->SetLineWidth(2);
    deltaXDune_pass2->SetLineWidth(2);
    deltaXDune_pass1->SetLineWidth(2);

    TCanvas *Canvas639 = new TCanvas("c639", "Graph Draw Options", 200, 10, 600, 400);
    deltaXUboone_pass2->Draw("hist");
    deltaXUboone_pass1->Draw("histsame");
    deltaXDune_pass2->Draw("histsame");
    deltaXDune_pass1->Draw("histsame");

    deltaXUboone_pass2->SetTitle("DeltaX");
    deltaXUboone_pass2->GetXaxis()->SetTitle("x_{Reco} - x_{True} (cm)");
    deltaXUboone_pass2->GetYaxis()->SetTitle("Percentage of Events (%)");
    deltaXUboone_pass2->SetStats(0);
    deltaXUboone_pass2->GetYaxis()->SetRangeUser(0, 38);

    auto legend_deltaXPass12 = new TLegend(0.1,0.7,0.48,0.9);
    legend_deltaXPass12->AddEntry(deltaXDune_pass1, "Deep Learning: DUNE/LBNF Tune Pass 1", "f");
    legend_deltaXPass12->AddEntry(deltaXDune_pass2, "Deep Learning: DUNE/LBNF Tune Pass 2", "f");
    legend_deltaXPass12->AddEntry(deltaXUboone_pass1, "Deep Learning: #muBooNE/BNB Tune Pass 1", "f");
    legend_deltaXPass12->AddEntry(deltaXUboone_pass2, "Deep Learning: #muBooNE/BNB Tune Pass 2", "f");
    legend_deltaXPass12->Draw();
    
    double binWidthROC = 2;
    double numOfBinsROC = std::ceil(largestValueROC)/binWidthROC;

    TH1F *DL_dune_deltaR_4roc_pass1 = new TH1F("deltaR for ROC", "deltaR for ROC", numOfBinsROC, 0, std::ceil(largestValueROC));
    TH1F *DL_dune_deltaR_4roc_pass2 = (TH1F*) DL_dune_deltaR_4roc_pass1->Clone("deltaR for ROC");
    TH1F *DL_uboone_deltaR_4roc_pass1 = (TH1F*) DL_dune_deltaR_4roc_pass1->Clone("deltaR for ROC");
    TH1F *DL_uboone_deltaR_4roc_pass2 = (TH1F*) DL_dune_deltaR_4roc_pass1->Clone("deltaR for ROC");
    TH1F *Current_deltaR_4roc = (TH1F*) DL_dune_deltaR_4roc_pass1->Clone("deltaR for ROC");

    for(UInt_t a = 0; a < numDLDuneVerticesPass1; a++){
        vertex_t vertex = dlDuneVerticesPass1.at(a);
        DL_dune_deltaR_4roc_pass1->Fill(vertex.deltaR);
    }

    for(UInt_t b = 0; b < numDLDuneVerticesPass2; b++){
        vertex_t vertex = dlDuneVerticesPass2.at(b);
        DL_dune_deltaR_4roc_pass2->Fill(vertex.deltaR);
    }

    for(UInt_t c = 0; c < numDLUbooneVerticesPass1; c++){
        vertex_t vertex = dlUbooneVerticesPass1.at(c);
        DL_uboone_deltaR_4roc_pass1->Fill(vertex.deltaR);
    }

    for(UInt_t d = 0; d < numDLUbooneVerticesPass2; d++){
        vertex_t vertex = dlUbooneVerticesPass2.at(d);
        DL_uboone_deltaR_4roc_pass2->Fill(vertex.deltaR);
    }

    for(UInt_t e = 0; e < numCurrentVertices; e++){
        vertex_t vertex = currentVertices.at(e);
        Current_deltaR_4roc->Fill(vertex.deltaR);
    }

    TH1F *DL_dune_deltaR_roc_pass1 = new TH1F("deltaR ROC", "deltaR ROC", numOfBinsROC, 0, std::ceil(largestValueROC));
    TH1F *DL_uboone_deltaR_roc_pass1 = (TH1F*) DL_dune_deltaR_roc_pass1->Clone("deltaR ROC");
    TH1F *Current_deltaR_roc = (TH1F*) DL_dune_deltaR_roc_pass1->Clone("deltaR ROC");
    TH1F *DL_dune_deltaR_roc_pass2 = (TH1F*) DL_dune_deltaR_roc_pass1->Clone("deltaR ROC");
    TH1F *DL_uboone_deltaR_roc_pass2 = (TH1F*) DL_dune_deltaR_roc_pass1->Clone("deltaR ROC");

    double runningTotal_dlDune_pass1 = 0;
    double runningTotal_dlDune_pass2 = 0;
    double runningTotal_dlUboone_pass1 = 0;
    double runningTotal_dlUboone_pass2 = 0;
    double runningTotal_current = 0;

    double binFillingROC = binWidthROC/2;

    for(int i = 1; i < numOfBinsROC+1; i++){
        runningTotal_dlDune_pass1 += DL_dune_deltaR_4roc_pass1->GetBinContent(i);
        runningTotal_dlDune_pass2 += DL_dune_deltaR_4roc_pass2->GetBinContent(i);
        runningTotal_dlUboone_pass1 += DL_uboone_deltaR_4roc_pass1->GetBinContent(i);
        runningTotal_dlUboone_pass2 += DL_uboone_deltaR_4roc_pass2->GetBinContent(i);
        runningTotal_current += Current_deltaR_4roc->GetBinContent(i);

        DL_dune_deltaR_roc_pass1->Fill(binFillingROC, (100.0f*runningTotal_dlDune_pass1/(double)numDLDuneVerticesPass1));
        DL_dune_deltaR_roc_pass2->Fill(binFillingROC, (100.0f*runningTotal_dlDune_pass2/(double)numDLDuneVerticesPass2));
        DL_uboone_deltaR_roc_pass1->Fill(binFillingROC, (100.0f*runningTotal_dlUboone_pass1/(double)numDLUbooneVerticesPass1));
        DL_uboone_deltaR_roc_pass2->Fill(binFillingROC, (100.0f*runningTotal_dlUboone_pass2/(double)numDLUbooneVerticesPass2));
        Current_deltaR_roc->Fill(binFillingROC, (100.0f*runningTotal_current/(double)numCurrentVertices));

        binFillingROC += binWidthROC;
    }

    DL_dune_deltaR_roc_pass1->SetLineWidth(2);
    DL_dune_deltaR_roc_pass2->SetLineWidth(2);
    DL_uboone_deltaR_roc_pass1->SetLineWidth(2);
    DL_uboone_deltaR_roc_pass2->SetLineWidth(2);
    Current_deltaR_roc->SetLineWidth(2);
    DL_dune_deltaR_roc_pass1->SetLineColor(kAzure+2);
    DL_dune_deltaR_roc_pass2->SetLineColor(kBlue+1);
    DL_uboone_deltaR_roc_pass1->SetLineColor(kGreen-3);
    DL_uboone_deltaR_roc_pass2->SetLineColor(kGreen+3);
    Current_deltaR_roc->SetLineColor(kRed);

    TCanvas *Canvas37746 = new TCanvas("c37746", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaR_roc_pass1->Draw("hist");
    DL_dune_deltaR_roc_pass2->Draw("histsame");
    DL_uboone_deltaR_roc_pass1->Draw("histsame");
    DL_uboone_deltaR_roc_pass2->Draw("histsame");
    Current_deltaR_roc->Draw("histsame");
    DL_dune_deltaR_roc_pass1->SetTitle("DeltaR ROC");
    DL_dune_deltaR_roc_pass1->SetStats(0);
    //DL_dune_deltaR_roc_pass1->GetXaxis()->SetRangeUser(0., 150.);
    DL_dune_deltaR_roc_pass1->GetYaxis()->SetRangeUser(0., 105);

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    // legend->SetHeader("legend header", "C");
    legend->AddEntry(DL_dune_deltaR_roc_pass1, "Deep Learning: DUNE/LBNF Tune Pass 1", "f");
    legend->AddEntry(DL_dune_deltaR_roc_pass2, "Deep Learning: DUNE/LBNF Tune Pass 2", "f");
    legend->AddEntry(DL_uboone_deltaR_roc_pass1, "Deep Learning: #muBooNE/BNB Tune Pass 1", "f");
    legend->AddEntry(DL_uboone_deltaR_roc_pass2, "Deep Learning: #muBooNE/BNB Tune Pass 2", "f");
    legend->AddEntry(Current_deltaR_roc, "Current SBND Vertexing", "f");
    legend->Draw();

    gPad->Modified();
    gPad->Update();
    TLine *line = new TLine(gPad->GetUxmin(), 100, gPad->GetUxmax(), 100);
    line->SetLineWidth(2);
    line->SetLineStyle(9);
    line->SetLineColor(kBlack);
    line->Draw();
}

void tpcPlots(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();

    std::vector<vertex_t> dlDuneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> currentVertices = std::vector<vertex_t>();

    for(UInt_t i = 0; i < numVertices; i++){
        vertex_t vertex = allVertices_vec.at(i);

        if(vertex.DLCurrent == 0){
            dlUbooneVertices.push_back(vertex);
        } else if(vertex.DLCurrent == 1){
            dlDuneVertices.push_back(vertex);
        } else if(vertex.DLCurrent == 2){
            currentVertices.push_back(vertex);
        }
    }

    TH1F *TPC0_current = new TH1F("deltaX TPC 0 Current", "deltaX TPC 0 Current", 16, -2, 2);
    TH1F *TPC1_current = (TH1F*) TPC0_current->Clone("deltaX TPC 1 Current");

    TH1F *TPC0_current_corrected = new TH1F("deltaX Corrected TPC 0", "deltaX Corrected TPC 0", 16, -2, 2);
    TH1F *TPC1_current_corrected = (TH1F*) TPC0_current_corrected->Clone("deltaX Corrected TPC 1");

    for(UInt_t j = 0; j < currentVertices.size(); j++){
        if(currentVertices.at(j).tpcID == 0){
            TPC0_current->Fill(currentVertices.at(j).uncorrectedDeltaX);
            TPC0_current_corrected->Fill(currentVertices.at(j).deltaX);
        } else if(currentVertices.at(j).tpcID == 1){
            TPC1_current->Fill(currentVertices.at(j).uncorrectedDeltaX);
            TPC1_current_corrected->Fill(currentVertices.at(j).deltaX);
        }
    }

    TPC0_current->SetLineColor(kBlue);
    TPC0_current->SetLineWidth(2);
    TPC1_current->SetLineColor(kRed);
    TPC1_current->SetLineWidth(2);

    TPC0_current_corrected->SetLineColor(kBlue);
    TPC0_current_corrected->SetLineWidth(2);
    TPC1_current_corrected->SetLineColor(kRed);
    TPC1_current_corrected->SetLineWidth(2);

    TCanvas *Canvas7 = new TCanvas("c7", "Graph Draw Options", 200, 10, 600, 400);
    TPC0_current->Draw("hist");
    TPC1_current->Draw("histsame");

    TCanvas *Canvas8 = new TCanvas("c8", "Graph Draw Options", 200, 10, 600, 400);
    TPC0_current_corrected->Draw("hist");
    TPC1_current_corrected->Draw("histsame");
   
    TPC0_current->SetStats(0);
    TPC0_current->SetTitle("");
    TPC0_current->GetXaxis()->SetTitle("#Delta X (cm)");
    TPC0_current->GetYaxis()->SetTitle("Number of Events");

    auto legend = new TLegend(0.1,0.7,0.48,0.49);
    // legend->SetHeader("legend header", "C");
    legend->AddEntry(TPC0_current, "TPC 0", "f");
    legend->AddEntry(TPC1_current, "TPC 1", "f");
    legend->Draw();
}

void gridComparison(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();
    std::vector<vertex_t> dlDuneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> currentVertices = std::vector<vertex_t>();
    for(UInt_t i = 0; i < numVertices; i++){
        vertex_t vertex = allVertices_vec.at(i);
        if(vertex.DLCurrent == 0){
            dlUbooneVertices.push_back(vertex);
        } else if(vertex.DLCurrent == 1){
            dlDuneVertices.push_back(vertex);
        } else if(vertex.DLCurrent == 2){
            currentVertices.push_back(vertex);
        }
    }

    int numPoints = 0;

    for(UInt_t i = 0; i < dlUbooneVertices.size(); i++){
        vertex_t ubooneVertex = dlUbooneVertices.at(i);
        for(UInt_t j = 0; j < currentVertices.size(); j++){
            vertex_t currentVertex = currentVertices.at(j);
            if((ubooneVertex.eventID == currentVertex.eventID) && (ubooneVertex.runID == currentVertex.runID) && (ubooneVertex.subRunID == currentVertex.subRunID)){
                for(UInt_t k = 0; k < dlDuneVertices.size(); k++){
                    vertex_t duneVertex = dlDuneVertices.at(k);
                    if((duneVertex.eventID == currentVertex.eventID) && (duneVertex.runID == currentVertex.runID) && (duneVertex.subRunID == currentVertex.subRunID)){
                        // If statement for adding points to the plot
                        if((ubooneVertex.deltaR < 10) && (currentVertex.deltaR < 10)){
                            numPoints++;
                        }

                        // If statement for printing out points
                        if(((ubooneVertex.deltaR > 10) && (currentVertex.deltaR < 0.9)) || ((currentVertex.deltaR > 10) && (ubooneVertex.deltaR < 0.9))){
                            std::cout << ubooneVertex.runID << ":" << ubooneVertex.subRunID << ":" << ubooneVertex.eventID << std::endl;
                            std::cout << "uboone deltaR = " << ubooneVertex.deltaR << " current deltaR = " << currentVertex.deltaR << " dune deltaR = " << duneVertex.deltaR << std::endl;
                            std::cout << "file number: " << ubooneVertex.fileNumber << ", " << currentVertex.fileNumber << std::endl;
                            std::cout << "" << std::endl;
                        }
                    }
                }
            }
        }
    }

    double x[numPoints];
    double y[numPoints];
    int counter = 0;
    for(UInt_t i = 0; i < dlUbooneVertices.size(); i++){
        vertex_t ubooneVertex = dlUbooneVertices.at(i);
        for(UInt_t j = 0; j < currentVertices.size(); j++){
            vertex_t currentVertex = currentVertices.at(j);
            if((ubooneVertex.eventID == currentVertex.eventID) && (ubooneVertex.runID == currentVertex.runID) && (ubooneVertex.subRunID == currentVertex.subRunID)){
                for(UInt_t k = 0; k < dlDuneVertices.size(); k++){
                    vertex_t duneVertex = dlDuneVertices.at(k);
                    if((duneVertex.eventID == currentVertex.eventID) && (duneVertex.runID == currentVertex.runID) && (duneVertex.subRunID == currentVertex.subRunID)){
                        // If statement for adding points to plot (must be the same as the above one)
                        if((ubooneVertex.deltaR < 10) && (currentVertex.deltaR < 10)){
                                x[counter] = currentVertex.deltaR;
                                y[counter] = ubooneVertex.deltaR;
                                counter++;
                        }
                    }
                }
            }
            
        }
    }

    TGraph *graph = new TGraph(numPoints, x, y);
    TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);
    graph->Draw("A*"); 
}

void calculations(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();
    std::vector<vertex_t> dlDuneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVertices = std::vector<vertex_t>();

    float below = 10000;
    double numDLDuneVertices = 0;
    double numDLDuneVerticesBelow = 0;
    double numDLUbooneVertices = 0;
    double numDLUbooneVerticesBelow = 0;
    double numCurrentVertices = 0;
    double numCurrentVerticesBelow = 0;

    for(UInt_t i = 0; i < numVertices; i++){
        vertex_t vertex = allVertices_vec.at(i);

        if(vertex.DLCurrent == 0){
            numDLUbooneVertices++;
            if(vertex.deltaR <= below){
                numDLUbooneVerticesBelow++;
            }
        } else if(vertex.DLCurrent == 1){
            numDLDuneVertices++;
            if(vertex.deltaR <= below){
                numDLDuneVerticesBelow++;
            }
        } else if(vertex.DLCurrent == 2){
            numCurrentVertices++;
            if(vertex.deltaR <= below){
                numCurrentVerticesBelow++;
            }
        }
    }
    double dlUboonePerc = ((numDLUbooneVerticesBelow/numDLUbooneVertices)*100);
    double dlDunePerc = ((numDLDuneVerticesBelow/numDLDuneVertices)*100);
    double currentPerc = ((numCurrentVerticesBelow/numCurrentVertices)*100);
    printf("Reco vertex with 10 cm of true vertex: \nDL DUNE: %f out of %f = %f%%\nMicroBooNe: %f out of %f = %f%%\nCurrent: %f out of %f = %f%%\n", numDLDuneVerticesBelow, numDLDuneVertices, dlDunePerc, numDLUbooneVerticesBelow, numDLUbooneVertices, dlUboonePerc, numCurrentVerticesBelow, numCurrentVertices, currentPerc);

}

void rocTogether(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();

    std::vector<vertex_t> dlDuneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> currentVertices = std::vector<vertex_t>();

    double largestValue = 0;

    for(UInt_t i = 0; i < numVertices; i++){
	    vertex_t vertex = allVertices_vec.at(i);
	    if(vertex.DLCurrent == 0){
		    dlUbooneVertices.push_back(vertex);
	    } else if(vertex.DLCurrent == 1){
		    dlDuneVertices.push_back(vertex);
	    } else if(vertex.DLCurrent == 2){
            currentVertices.push_back(vertex);
        }

	    if(largestValue < vertex.deltaR){
		    largestValue = vertex.deltaR;
	    }
    }

    std::cout << "Size of DUNE DL: " << dlDuneVertices.size() << " Size of UBooNE DL: " << dlUbooneVertices.size() << " Size of Current: " << currentVertices.size() << std::endl;

    double binWidth = 3;
    double numberOfBins = std::ceil(largestValue/binWidth);
    //double numberOfBins = std::ceil(150/binWidth);
    
    //TH1F *DL_dune_deltaR_4roc = new TH1F("deltaR for ROC", "R Coord Diff for ROC", numberOfBins, 0, 150);
    TH1F *DL_dune_deltaR_4roc = new TH1F("deltaR for ROC", "R Coord Diff for ROC", numberOfBins, 0, std::ceil(largestValue));
    TH1F *DL_uboone_deltaR_4roc = (TH1F*) DL_dune_deltaR_4roc->Clone("deltaR for ROC");
    TH1F *Current_deltaR_4roc = (TH1F*) DL_dune_deltaR_4roc->Clone("deltaR for ROC");

    double numDLDuneVertices = dlDuneVertices.size();
    double numDLUbooneVertices = dlUbooneVertices.size();
    double numCurrentVertices = currentVertices.size();

    for(UInt_t j = 0; j < numDLDuneVertices; j++){
	    vertex_t vertex = dlDuneVertices.at(j);
	    DL_dune_deltaR_4roc->Fill(vertex.deltaR);
    }

    for(UInt_t j = 0; j < numDLUbooneVertices; j++){
        vertex_t vertex = dlUbooneVertices.at(j);
        DL_uboone_deltaR_4roc->Fill(vertex.deltaR);
    }

    for(UInt_t j = 0; j < numCurrentVertices; j++){
      	vertex_t vertex = currentVertices.at(j);
   	    Current_deltaR_4roc->Fill(vertex.deltaR);
    }

    //TH1F *DL_dune_deltaR_roc = new TH1F("deltaR ROC", "R Coord Diff ROC", numberOfBins, 0, 150);
    
    // Drawing all 3 ROC curves
    //TH1F *DL_dune_deltaR_roc = new TH1F("deltaR ROC", "R Coord Diff ROC", numberOfBins, 0, std::ceil(largestValue));
    //TH1F *DL_uboone_deltaR_roc = (TH1F*) DL_dune_deltaR_roc->Clone("deltaR ROC");
    //TH1F *Current_deltaR_roc = (TH1F*) DL_dune_deltaR_roc->Clone("deltaR ROC");
    
    // Drawing uboone and current ROC curves
    TH1F *Current_deltaR_roc = new TH1F("deltaR ROC", "R Coord Diff ROC", numberOfBins, 0, std::ceil(largestValue));
    TH1F *DL_uboone_deltaR_roc = (TH1F*) Current_deltaR_roc->Clone("deltaR ROC");
    TH1F *DL_dune_deltaR_roc = (TH1F*) Current_deltaR_roc->Clone("deltaR ROC");

    double runningTotal_dlDune = 0;
    double runningTotal_dlUboone = 0;
    double runningTotal_current = 0;

    double binFilling = binWidth/2;

    for(int i = 1; i < numberOfBins+1; i++){
	    runningTotal_dlDune += DL_dune_deltaR_4roc->GetBinContent(i);
        runningTotal_dlUboone += DL_uboone_deltaR_4roc->GetBinContent(i);
	    runningTotal_current += Current_deltaR_4roc->GetBinContent(i);
 	   	
	    DL_dune_deltaR_roc->Fill(binFilling, (100.0f*runningTotal_dlDune/(double)numDLDuneVertices));
        DL_uboone_deltaR_roc->Fill(binFilling, (100.0f*runningTotal_dlUboone/(double)numDLUbooneVertices));
	    Current_deltaR_roc->Fill(binFilling, (100.0f*runningTotal_current/(double)numCurrentVertices));
	
	    binFilling += binWidth;

    }

    double lowestValue = 1000;
    if(DL_dune_deltaR_roc->GetBinContent(1) < lowestValue) lowestValue = DL_dune_deltaR_roc->GetBinContent(1);
    if(DL_uboone_deltaR_roc->GetBinContent(1) < lowestValue) lowestValue = DL_uboone_deltaR_roc->GetBinContent(1);
    if(Current_deltaR_roc->GetBinContent(1) < lowestValue) lowestValue = Current_deltaR_roc->GetBinContent(1);

    TColor *c1 = new TColor(9001,0,8,0);
    TColor *c2 = new TColor(9002,0,0,1);
    TColor *c3 = new TColor(9003,1,0,0);

    TCanvas *Canvas4 = new TCanvas("c4", "Graph Draw Options", 200, 10, 600, 400);

    DL_dune_deltaR_roc->SetLineColor(kViolet-5);
    //Current_deltaR_roc->SetLineColor(kBlack);
    DL_dune_deltaR_roc->SetLineWidth(2);
    DL_uboone_deltaR_roc->SetLineColor(kBlue);
    DL_uboone_deltaR_roc->SetLineWidth(2);
    Current_deltaR_roc->SetLineColor(kRed);
    Current_deltaR_roc->SetLineWidth(2);
   
    Current_deltaR_roc->Draw("hist");
    DL_uboone_deltaR_roc->Draw("histsame"); 
    DL_dune_deltaR_roc->Draw("histsame");

    Current_deltaR_roc->SetStats(0);
    //Current_deltaR_roc->SetTitle("ROC Curves: Neutrino Events with Cosmics");
    Current_deltaR_roc->SetTitle("ROC Curves: Nu + E Elastic Scattering Events");
    Current_deltaR_roc->GetXaxis()->SetTitle("|#bar{r}_{Reco} - #bar{r}_{True}| (cm)");
    Current_deltaR_roc->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaR_roc->GetYaxis()->SetRangeUser(lowestValue - 5, 105.);

    DL_dune_deltaR_4roc->SetLineColor(9001);
    DL_dune_deltaR_4roc->SetLineWidth(2);
    DL_uboone_deltaR_4roc->SetLineColor(kBlue);
    DL_uboone_deltaR_4roc->SetLineWidth(2);
    Current_deltaR_4roc->SetLineColor(kRed);
    Current_deltaR_4roc->SetLineWidth(2);

    auto legend = new TLegend(0.42,0.28,0.85,0.15);
    // legend->SetHeader("legend header", "C");
    legend->AddEntry(DL_dune_deltaR_roc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(DL_uboone_deltaR_roc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(Current_deltaR_roc, "Current SBND Vertexing (without Refinement)", "f");
    legend->Draw();
    
    gPad->Modified();
    gPad->Update();
    TLine *line = new TLine(gPad->GetUxmin(), 100, gPad->GetUxmax(), 100);
    line->SetLineWidth(2);
    line->SetLineStyle(9);
    line->SetLineColor(kBlack);
    line->Draw();

    gPad->SetTickx();
    gPad->SetTicky();

    //std::string fileName = Form("/nashome/c/coackley/DLVertexPlots/ROC_3March25_%.2f_withoutRefinement.pdf", binWidth);
    std::string fileName = Form("/nashome/c/coackley/nuEPlots/ROC_3March25_%.2f_withoutRefinement.pdf", binWidth);
    Canvas4->SaveAs(fileName.c_str());

    line->Delete();

    double lowX = 0.;
    double highX = 5;
    int binNumber = Current_deltaR_roc->FindBin(highX);

    Current_deltaR_roc->GetYaxis()->SetRangeUser(lowestValue - 5, Current_deltaR_roc->GetBinContent(binNumber) + 10);
    Current_deltaR_roc->GetXaxis()->SetRangeUser(lowX, highX);
    //fileName = Form("/nashome/c/coackley/DLVertexPlots/ROC_3March25_%.2f_%.2f-%.2f_withoutRefinement.pdf", binWidth, lowX, highX);
    fileName = Form("/nashome/c/coackley/nuEPlots/ROC_3March25_%.2f_%.2f-%.2f_withoutRefinement.pdf", binWidth, lowX, highX);
    Canvas4->SaveAs(fileName.c_str());

    highX = 10;
    binNumber = Current_deltaR_roc->FindBin(highX);

    Current_deltaR_roc->GetYaxis()->SetRangeUser(lowestValue - 5, Current_deltaR_roc->GetBinContent(binNumber) + 10);
    Current_deltaR_roc->GetXaxis()->SetRangeUser(lowX, highX);
    //fileName = Form("/nashome/c/coackley/DLVertexPlots/ROC_3March25_%.2f_%.2f-%.2f_withoutRefinement.pdf", binWidth, lowX, highX);
    fileName = Form("/nashome/c/coackley/nuEPlots/ROC_3March25_%.2f_%.2f-%.2f_withoutRefinement.pdf", binWidth, lowX, highX);
    Canvas4->SaveAs(fileName.c_str());
}

void plottingTogetherDist(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();

    std::vector<vertex_t> dlDuneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> dlUbooneVertices = std::vector<vertex_t>();
    std::vector<vertex_t> currentVertices = std::vector<vertex_t>();

    double largestValueX = 0;
    double smallestValueX = 0;
    double largestValueY = 0;
    double smallestValueY = 0;
    double largestValueZ = 0;
    double smallestValueZ = 0;
    double largestValueR = 0;
    double smallestValueR = 0;

    for(UInt_t i = 0; i < numVertices; i++){
	    vertex_t vertex = allVertices_vec.at(i);
    	if(vertex.DLCurrent == 0){
		    dlUbooneVertices.push_back(vertex);
    	} else if(vertex.DLCurrent == 1){
            dlDuneVertices.push_back(vertex);
        } else if(vertex.DLCurrent == 2){
		    currentVertices.push_back(vertex);
    	}

	    if(largestValueX < vertex.deltaX){
		    largestValueX = vertex.deltaX;
	    } else if(smallestValueX > vertex.deltaX){
		    smallestValueX = vertex.deltaX;
    	}   
    	if(largestValueY < vertex.deltaY){
    		largestValueY = vertex.deltaY;
	    } else if(smallestValueY > vertex.deltaY){
    		smallestValueY = vertex.deltaY;
    	}
    	if(largestValueZ < vertex.deltaZ){
     		largestValueZ = vertex.deltaZ;
    	} else if(smallestValueZ > vertex.deltaZ){
    		smallestValueZ = vertex.deltaZ;
    	}
    	if(largestValueR < vertex.deltaR){
    		largestValueR = vertex.deltaR;
    	} else if(smallestValueR > vertex.deltaR){
    		smallestValueR = vertex.deltaR;
    	}
    }

    UInt_t numDLDuneVertices = dlDuneVertices.size();
    UInt_t numDLUbooneVertices = dlUbooneVertices.size();
    UInt_t numCurrentVertices = currentVertices.size();

    double binWidth = 0.47619048;
    double lowerBoundX = -5; // std::floor(smallestValueX); 
    double upperBoundX = 5;  //std::ceil(largestValueX);
    double lowerBoundY = -5; //std::floor(smallestValueY);
    double upperBoundY = 5;//std::ceil(largestValueY);
    double lowerBoundZ = -5;//std::floor(smallestValueZ);
    double upperBoundZ = 5;//std::ceil(largestValueZ);
    double upperBoundR = 5;//std::ceil(largestValueR);
    //double lowerBound = -4;
    //double upperBound = 4;
    double numOfBinsX = std::ceil((upperBoundX - lowerBoundX)/binWidth);
    double numOfBinsY = std::ceil((upperBoundY - lowerBoundY)/binWidth);
    double numOfBinsZ = std::ceil((upperBoundZ - lowerBoundZ)/binWidth);
    double numOfBinsR = std::ceil((upperBoundR)/binWidth);

    TH1F *Current_deltaX_freq = new TH1F("deltaX Current freq", "Current x Coord Diff freq", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *DL_uboone_deltaX_freq = (TH1F*) Current_deltaX_freq->Clone("deltaX DL uboone freq");
    TH1F *DL_dune_deltaX_freq = (TH1F*) Current_deltaX_freq->Clone("deltaX DL dune freq");

    TH1F *Current_deltaY_freq = new TH1F("deltaY DL dune freq", "DL dune Y Coord Diff freq", numOfBinsY, lowerBoundY, upperBoundY);
    TH1F *DL_uboone_deltaY_freq = (TH1F*) Current_deltaY_freq->Clone("deltaY DL uboone freq");
    TH1F *DL_dune_deltaY_freq = (TH1F*) Current_deltaY_freq->Clone("deltaY Current freq");

    TH1F *Current_deltaZ_freq = new TH1F("deltaZ Current freq", "Current Z Coord Diff freq", numOfBinsZ, lowerBoundZ, upperBoundZ);
    TH1F *DL_uboone_deltaZ_freq = (TH1F*) Current_deltaZ_freq->Clone("deltaZ DL uboone freq");
    TH1F *DL_dune_deltaZ_freq = (TH1F*) Current_deltaZ_freq->Clone("deltaZ DL dune freq");

    TH1F *Current_deltaR_freq = new TH1F("deltaR Current freq", "Current R Coord Diff freq", numOfBinsR, 0, upperBoundR);
    TH1F *DL_uboone_deltaR_freq = (TH1F*) Current_deltaR_freq->Clone("deltaR DL uboone freq");
    TH1F *DL_dune_deltaR_freq = (TH1F*) Current_deltaR_freq->Clone("deltaR DL dune freq");


    for(UInt_t j = 0; j < numDLDuneVertices; j++){
	    vertex_t dlDuneVertex = dlDuneVertices.at(j);
    	DL_dune_deltaX_freq->Fill(dlDuneVertex.deltaX);
	    DL_dune_deltaY_freq->Fill(dlDuneVertex.deltaY);
	    DL_dune_deltaZ_freq->Fill(dlDuneVertex.deltaZ);
        DL_dune_deltaR_freq->Fill(dlDuneVertex.deltaR);
    }

    for(UInt_t j = 0; j < numDLUbooneVertices; j++){
        vertex_t dlUbooneVertex = dlUbooneVertices.at(j);
        DL_uboone_deltaX_freq->Fill(dlUbooneVertex.deltaX);
        DL_uboone_deltaY_freq->Fill(dlUbooneVertex.deltaY);
        DL_uboone_deltaZ_freq->Fill(dlUbooneVertex.deltaZ);
        DL_uboone_deltaR_freq->Fill(dlUbooneVertex.deltaR);
    }

    for(UInt_t j = 0; j < numCurrentVertices; j++){
    	vertex_t currentVertex = currentVertices.at(j);
 	    Current_deltaX_freq->Fill(currentVertex.deltaX);
	    Current_deltaY_freq->Fill(currentVertex.deltaY);
	    Current_deltaZ_freq->Fill(currentVertex.deltaZ);
	    Current_deltaR_freq->Fill(currentVertex.deltaR);
    }

    TColor *c1 = new TColor(9001,0,1,0);
    TColor *c2 = new TColor(9002,0,0,1);
    TColor *c3 = new TColor(9003,1,0,0);

    TCanvas *Canvas0 = new TCanvas("c0", "Graph Draw Options", 200, 10, 600, 400);

    TH1F *Current_deltaX_perc = new TH1F("deltaX DL dune %", "x Coord Diff %", numOfBinsX, lowerBoundX, upperBoundX);
    TH1F *DL_uboone_deltaX_perc = (TH1F*) Current_deltaX_perc->Clone("deltaX DL uboone %");
    TH1F *DL_dune_deltaX_perc = (TH1F*) Current_deltaX_perc->Clone("deltaX Current %");

    TH1F *Current_deltaY_perc = new TH1F("deltaY DL dune %", "Y Coord Diff %", numOfBinsY, lowerBoundY, upperBoundY);
    TH1F *DL_uboone_deltaY_perc = (TH1F*) Current_deltaY_perc->Clone("deltaY DL uboone %");
    TH1F *DL_dune_deltaY_perc = (TH1F*) Current_deltaY_perc->Clone("deltaY Current %");

    TH1F *Current_deltaZ_perc = new TH1F("deltaZ DL dune %", "Z Coord Diff %", numOfBinsZ, lowerBoundZ, upperBoundZ);
    TH1F *DL_uboone_deltaZ_perc = (TH1F*) Current_deltaZ_perc->Clone("deltaZ DL uboone %");
    TH1F *DL_dune_deltaZ_perc = (TH1F*) Current_deltaZ_perc->Clone("deltaZ Current %");

    TH1F *Current_deltaR_perc = new TH1F("deltaR DL dune %", "R Coord Diff %", numOfBinsR, 0, upperBoundR);
    TH1F *DL_uboone_deltaR_perc = (TH1F*) Current_deltaR_perc->Clone("deltaR DL uboone %");
    TH1F *DL_dune_deltaR_perc = (TH1F*) Current_deltaR_perc->Clone("deltaR Current %");

    double binFillingX = lowerBoundX + (binWidth/2);
    double binFillingY = lowerBoundY + (binWidth/2);
    double binFillingZ = lowerBoundZ + (binWidth/2);
    double binFillingR = 0 + (binWidth/2);

    for(int i = 1; i < numOfBinsX+1; i++){
    	DL_dune_deltaX_perc->Fill(binFillingX, (100.0f*DL_dune_deltaX_freq->GetBinContent(i)/(double)numDLDuneVertices));
        DL_uboone_deltaX_perc->Fill(binFillingX, (100.0f*DL_uboone_deltaX_freq->GetBinContent(i)/(double)numDLUbooneVertices));
	    Current_deltaX_perc->Fill(binFillingX, (100.0f*Current_deltaX_freq->GetBinContent(i)/(double)numCurrentVertices));
	    binFillingX += binWidth;
    }

    for(int i = 1; i < numOfBinsY+1; i++){
        DL_dune_deltaY_perc->Fill(binFillingY, (100.0f*DL_dune_deltaY_freq->GetBinContent(i)/(double)numDLDuneVertices));
        DL_uboone_deltaY_perc->Fill(binFillingY, (100.0f*DL_uboone_deltaY_freq->GetBinContent(i)/(double)numDLUbooneVertices));
        Current_deltaY_perc->Fill(binFillingY, (100.0f*Current_deltaY_freq->GetBinContent(i)/(double)numCurrentVertices));
        binFillingY += binWidth;
    }

    for(int i = 1; i < numOfBinsZ+1; i++){
        DL_dune_deltaZ_perc->Fill(binFillingZ, (100.0f*DL_dune_deltaZ_freq->GetBinContent(i)/(double)numDLDuneVertices));
        DL_uboone_deltaZ_perc->Fill(binFillingZ, (100.0f*DL_uboone_deltaZ_freq->GetBinContent(i)/(double)numDLUbooneVertices));
        Current_deltaZ_perc->Fill(binFillingZ, (100.0f*Current_deltaZ_freq->GetBinContent(i)/(double)numCurrentVertices));
        binFillingZ += binWidth;
    }

    for(int i = 1; i < numOfBinsR+1; i++){
        DL_dune_deltaR_perc->Fill(binFillingR, (100.0f*DL_dune_deltaR_freq->GetBinContent(i)/(double)numDLDuneVertices));
        DL_uboone_deltaR_perc->Fill(binFillingR, (100.0f*DL_uboone_deltaR_freq->GetBinContent(i)/(double)numDLUbooneVertices));
        Current_deltaR_perc->Fill(binFillingR, (100.0f*Current_deltaR_freq->GetBinContent(i)/(double)numCurrentVertices));
        binFillingR += binWidth;
    }

    DL_dune_deltaX_perc->SetLineColor(kViolet-5);
    DL_dune_deltaX_perc->SetLineWidth(2);
    DL_uboone_deltaX_perc->SetLineColor(kBlue);
    DL_uboone_deltaX_perc->SetLineWidth(2);
    Current_deltaX_perc->SetLineColor(kRed);
    Current_deltaX_perc->SetLineWidth(2);

    //TCanvas *Canvas0 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);
    Current_deltaX_perc->Draw("hist");
    DL_uboone_deltaX_perc->Draw("histsame");
    DL_dune_deltaX_perc->Draw("histsame");

    Current_deltaX_perc->SetStats(0);
    Current_deltaX_perc->SetTitle("#Deltax Distribution: Nu + E Elastic Scattering Events");
    Current_deltaX_perc->GetXaxis()->SetTitle("x_{Reco} - x_{True} (cm)");
    Current_deltaX_perc->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaX_perc->GetYaxis()->SetRange(22,23);

    std::cout << "deltaX" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaX_perc->GetMean() << " RMS: " << Current_deltaX_perc->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaX_perc->GetMean() << " RMS: " << DL_dune_deltaX_perc->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaX_perc->GetMean() << "RMS: " << DL_uboone_deltaX_perc->GetRMS(1) << std::endl;

    auto legend = new TLegend(0.56,0.86,0.88,0.70);
    // legend->SetHeader("legend header", "C");
    legend->AddEntry(DL_dune_deltaX_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(DL_uboone_deltaX_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(Current_deltaX_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    TCanvas *Canvas1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaY_perc->SetLineColor(kViolet-5);
    DL_dune_deltaY_perc->SetLineWidth(2);
    DL_uboone_deltaY_perc->SetLineColor(kBlue);
    DL_uboone_deltaY_perc->SetLineWidth(2);
    Current_deltaY_perc->SetLineColor(kRed);
    Current_deltaY_perc->SetLineWidth(2);
    Current_deltaY_perc->Draw("hist");
    DL_uboone_deltaY_perc->Draw("histsame");
    DL_dune_deltaY_perc->Draw("histsame");

    Current_deltaY_perc->SetStats(0);
    Current_deltaY_perc->SetTitle("#Deltay Distribution: Nu + E Elastic Scattering Events");
    Current_deltaY_perc->GetXaxis()->SetTitle("y_{Reco} - y_{True} (cm)");
    Current_deltaY_perc->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaY_perc->GetYaxis()->SetRangeUser(0, 38);

    auto legend2 = new TLegend(0.56,0.86,0.88,0.70);
    // legend->SetHeader("legend header", "C");
    legend2->AddEntry(DL_dune_deltaY_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend2->AddEntry(DL_uboone_deltaY_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend2->AddEntry(Current_deltaY_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend2->SetTextSize(0.0225);
    legend2->SetMargin(0.13);
    legend2->Draw();

    std::cout << "deltaY" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaY_perc->GetMean() << " RMS: " << Current_deltaY_perc->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaY_perc->GetMean() << " RMS: " << DL_dune_deltaY_perc->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaY_perc->GetMean() << "RMS: " << DL_uboone_deltaY_perc->GetRMS(1) << std::endl;

    TCanvas *Canvas3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaZ_perc->SetLineColor(kViolet-5);
    DL_dune_deltaZ_perc->SetLineWidth(2);
    DL_uboone_deltaZ_perc->SetLineColor(kBlue);
    DL_uboone_deltaZ_perc->SetLineWidth(2);
    Current_deltaZ_perc->SetLineColor(kRed);
    Current_deltaZ_perc->SetLineWidth(2);
    Current_deltaZ_perc->Draw("hist");
    DL_uboone_deltaZ_perc->Draw("histsame");
    DL_dune_deltaZ_perc->Draw("histsame");

    Current_deltaZ_perc->SetStats(0);
    Current_deltaZ_perc->SetTitle("#Deltaz Distribution: Nu + E Elastic Scattering Events");
    Current_deltaZ_perc->GetXaxis()->SetTitle("z_{Reco} - z_{True} (cm)");
    Current_deltaZ_perc->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaZ_perc->GetYaxis()->SetRangeUser(0, 35);

    auto legend3 = new TLegend(0.56,0.86,0.88,0.70);
    legend3->AddEntry(DL_dune_deltaZ_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3->AddEntry(DL_uboone_deltaZ_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3->AddEntry(Current_deltaZ_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend3->SetTextSize(0.0225);
    legend3->SetMargin(0.13);
    legend3->Draw();

    std::cout << "deltaZ" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaZ_perc->GetMean() << " RMS: " << Current_deltaZ_perc->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaZ_perc->GetMean() << " RMS: " << DL_dune_deltaZ_perc->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaZ_perc->GetMean() << "RMS: " << DL_uboone_deltaZ_perc->GetRMS(1) << std::endl;

    TCanvas *Canvas4 = new TCanvas("c4", "Graph Draw Options", 200, 10, 600, 400);
    DL_dune_deltaR_perc->SetLineColor(kViolet-5);
    DL_dune_deltaR_perc->SetLineWidth(2);
    DL_uboone_deltaR_perc->SetLineColor(kBlue);
    DL_uboone_deltaR_perc->SetLineWidth(2);
    Current_deltaR_perc->SetLineColor(kRed);
    Current_deltaR_perc->SetLineWidth(2);
    Current_deltaR_perc->Draw("hist");
    DL_uboone_deltaR_perc->Draw("histsame");
    DL_dune_deltaR_perc->Draw("histsame"); 

    Current_deltaR_perc->SetStats(0);
    Current_deltaR_perc->SetTitle("#Delta#bar{r} Distribution: Nu + E Elastic Scattering Events");
    Current_deltaR_perc->GetXaxis()->SetTitle("|#bar{r}_{Reco} - #bar{r}_{True}| (cm)");
    Current_deltaR_perc->GetYaxis()->SetTitle("Percentage of Events (%)");
    Current_deltaR_perc->GetYaxis()->SetRangeUser(0, 27);

    auto legend4 = new TLegend(0.56,0.86,0.88,0.70);
    // legend->SetHeader("legend header", "C");
    legend4->AddEntry(DL_dune_deltaR_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend4->AddEntry(DL_uboone_deltaR_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend4->AddEntry(Current_deltaR_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend4->SetTextSize(0.0225);
    legend4->SetMargin(0.13);
    legend4->Draw();

    std::cout << "deltaR" << std::endl; 
    std::cout << "Current------Mean: " << Current_deltaR_perc->GetMean() << " RMS: " << Current_deltaR_perc->GetRMS(1) << std::endl;
    std::cout << "dune------Mean: " << DL_dune_deltaR_perc->GetMean() << " RMS: " << DL_dune_deltaR_perc->GetRMS(1) << std::endl;
    std::cout << "uboone------Mean: " << DL_uboone_deltaR_perc->GetMean() << " RMS: " << DL_uboone_deltaR_perc->GetRMS(1) << std::endl;

    Canvas0->SaveAs("/nashome/c/coackley/nuEPlots/deltaX_dist.pdf");
    Canvas1->SaveAs("/nashome/c/coackley/nuEPlots/deltaY_dist.pdf");
    Canvas3->SaveAs("/nashome/c/coackley/nuEPlots/deltaZ_dist.pdf");
    Canvas4->SaveAs("/nashome/c/coackley/nuEPlots/deltaR_dist.pdf");

    //TPaveText *t = new TPaveText(0,3,10,13);
    //t->AddText("Current SBND Vertexing"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kRed);
    //t->AddText("Deep Learning: DUNE/LBNF Tune");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kViolet-5);
    //t->AddText("Deep Learning: #muBooNE/BNB Tune"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
    //t->Draw();
    
}

void cosmicVertexHist(std::vector<vertex_t> allVertices_vec){
    size_t numVertices = allVertices_vec.size();

    double largestValue = 0;

    for(UInt_t i = 0; i < numVertices; i++){
	    vertex_t vertex = allVertices_vec.at(i);
    	if(largestValue < vertex.deltaR){
		    largestValue = vertex.deltaR;
	    }
    }

    double binWidth = 5;
    double numberOfBins = std::ceil(largestValue/binWidth);
    std::cout << "largest value: " << largestValue << " number of bins: " << numberOfBins << std::endl; 
    std::cout << "Put in graph, largest value: " <<  std::ceil(largestValue) << std::endl;

    TCanvas *Canvas0 = new TCanvas("c0", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* vertex_deltaX = new TH1F("Delta X", "X Coord Diff", 40, -4, 4);
    TCanvas *Canvas1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* vertex_deltaY = new TH1F("Delta Y", "Y Coord Diff", 40, -4, 4);
    TCanvas *Canvas2 = new TCanvas("c2", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* vertex_deltaZ = new TH1F("Delta Z", "Z Coord Diff", 40, -4, 4);
    TCanvas *Canvas3 = new TCanvas("c3", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* vertex_deltaR = new TH1F("Delta R", "R Coord Diff", 400, 0, 400);

    TH1F* vertex_deltaR_forROC = new TH1F("Delta R, for ROC", "R Coord, for ROC", numberOfBins, 0, std::ceil(largestValue));

    for(UInt_t j = 0; j < numVertices; j++){
        vertex_t vertex = allVertices_vec.at(j);
        vertex_deltaX->Fill(vertex.deltaX);
        vertex_deltaY->Fill(vertex.deltaY);
        vertex_deltaZ->Fill(vertex.deltaZ);
        vertex_deltaR->Fill(vertex.deltaR);
    	vertex_deltaR_forROC->Fill(vertex.deltaR);
     }

    TH1F* roc_deltaR = new TH1F("Delta R ROC", "R Coord ROC", numberOfBins, 0, std::ceil(largestValue));
    
    int runningTotal = 0;

    double binFilling = binWidth/2;

    for(int i = 1; i < numberOfBins+1; i++){
	
	    std::cout << "Bin number: " << i << std::endl;
	    std::cout << "Bin Content: " << vertex_deltaR_forROC->GetBinContent(i) << std::endl;
	    runningTotal += vertex_deltaR_forROC->GetBinContent(i);
	    std::cout << "Running Total: " << (double)runningTotal << " Num Vertices: " << (double)numVertices << std::endl;
	    std::cout << "Percentage: " << ((double)runningTotal/(double)numVertices)*100.0f << std::endl;

	    for(int j = 0; j < abs(((double)runningTotal/(double)numVertices)*100.0f); j++){
		    roc_deltaR->Fill(binFilling);
	    }   
	    binFilling += binWidth; 
    }
    vertex_deltaX->Draw();
    vertex_deltaY->Draw();
    vertex_deltaZ->Draw();
    vertex_deltaR->Draw();
}

void dl_plottingMacro(){
    TFile *f = TFile::Open("/exp/sbnd/data/users/coackley/cosmics/cosmic_3March25_VertexNotRefined_nuE.root", "READ");
    //TFile *f = TFile::Open("/exp/sbnd/app/users/coackley/v09_78_04/srcs/sbndcode/sbndcode/DLVertexAnalyzer/analysisOutput.root");
    if(!f){
        std::cout << "Failed to read file" << std::endl;
        return;
    }

    TTree *t;
    f->GetObject("goodTree", t);
 
    std::vector<double> allDeltaX = std::vector<double>(0);
    std::vector<double> allDeltaY = std::vector<double>(0);
    std::vector<double> allDeltaZ = std::vector<double>(0);
    std::vector<double> allDeltaR = std::vector<double>(0);
    std::vector<int> allCCNC = std::vector<int>(0);
    std::vector<int> allTrueNeutrino = std::vector<int>(0);
    std::vector<int> allDL_Current = std::vector<int>(0);
    std::vector<int> allCosmic_Single = std::vector<int>(0);
    std::vector<int> allPass1Pass2 = std::vector<int>(0);
    std::vector<unsigned int> allEventID = std::vector<unsigned int>(0);
    std::vector<unsigned int> allRunID = std::vector<unsigned int>(0);
    std::vector<unsigned int> allSubRunID = std::vector<unsigned int>(0);
    std::vector<int> allFileNumber = std::vector<int>(0);
    std::vector<double> allTPCID = std::vector<double>(0);
    std::vector<double> allUncorrectedDeltaX = std::vector<double>(0);
    std::vector<double> allUncorrectedDeltaR = std::vector<double>(0);
    std::vector<double> allRecoVertexX = std::vector<double>(0);
    std::vector<double> allRecoVertexY = std::vector<double>(0);
    std::vector<double> allRecoVertexZ = std::vector<double>(0);
    std::vector<double> allMCVertexX = std::vector<double>(0);
    std::vector<double> allMCVertexY = std::vector<double>(0);
    std::vector<double> allMCVertexZ = std::vector<double>(0);

    Long64_t numEntries = t->GetEntries();
    int fileNumber = 1;

    for(Long64_t i = 0; i < numEntries; i++){
        std::cout << fileNumber << std::endl;
	    std::vector<double> *pDeltaX = 0;
        std::vector<double> *pDeltaY = 0;
        std::vector<double> *pDeltaZ = 0;
        std::vector<double> *pDeltaR = 0;
        std::vector<int> *pCCNC = 0;
        std::vector<int> *pTrueNeutrino = 0;
        std::vector<int> *pDL_Current = 0;  
 	    std::vector<int> *pCosmic_Single = 0;
        std::vector<int> *pPass1Pass2 = 0;
        std::vector<unsigned int> *pEventID = 0;
        std::vector<unsigned int> *pRunID = 0;
        std::vector<unsigned int> *pSubRunID = 0;
        std::vector<double> *pTPCID = 0;
        std::vector<double> *pUncorrectedDeltaX = 0;
        std::vector<double> *pUncorrectedDeltaR = 0;
        std::vector<double> *pRecoVertexX = 0;
        std::vector<double> *pRecoVertexY = 0;
        std::vector<double> *pRecoVertexZ = 0;
        std::vector<double> *pMCVertexX = 0;
        std::vector<double> *pMCVertexY = 0;
        std::vector<double> *pMCVertexZ = 0;

    	TBranch *deltaX_branch = 0;
        TBranch *deltaY_branch = 0;
        TBranch *deltaZ_branch = 0; 
        TBranch *deltaR_branch = 0;
        TBranch *CCNC_branch = 0;
        TBranch *trueNeutrino_branch = 0;
        TBranch *dl_current_branch = 0;
	    TBranch *cosmic_single_branch = 0;
        TBranch *pass1pass2_branch = 0;
        TBranch *eventID_branch = 0;
        TBranch *runID_branch = 0;
        TBranch *subRunID_branch = 0;
        TBranch *tpcID_branch = 0;
        TBranch *uncorrectedDeltaX_branch = 0;
        TBranch *uncorrectedDeltaR_branch = 0;
        TBranch *recoVertexX_branch = 0;
        TBranch *recoVertexY_branch = 0;
        TBranch *recoVertexZ_branch = 0;
        TBranch *mcVertexX_branch = 0;
        TBranch *mcVertexY_branch = 0;
        TBranch *mcVertexZ_branch = 0;

	    t->SetBranchAddress("deltaX_goodTree", &pDeltaX, &deltaX_branch);
        t->SetBranchAddress("deltaY_goodTree", &pDeltaY, &deltaY_branch);
        t->SetBranchAddress("deltaZ_goodTree", &pDeltaZ, &deltaZ_branch);
        t->SetBranchAddress("deltaR_goodTree", &pDeltaR, &deltaR_branch);
        t->SetBranchAddress("CCNC_goodTree", &pCCNC, &CCNC_branch);
        t->SetBranchAddress("trueNeutrino_goodTree", &pTrueNeutrino, &trueNeutrino_branch);
        t->SetBranchAddress("dl_current_goodTree", &pDL_Current, &dl_current_branch);
	    t->SetBranchAddress("cosmic_single_goodTree", &pCosmic_Single, &cosmic_single_branch);
        t->SetBranchAddress("pass1_pass2_goodTree", &pPass1Pass2, &pass1pass2_branch);
	    t->SetBranchAddress("eventID_goodTree", &pEventID, &eventID_branch);
        t->SetBranchAddress("runID_goodTree", &pRunID, &runID_branch);
        t->SetBranchAddress("subRunID_goodTree", &pSubRunID, &subRunID_branch);
        t->SetBranchAddress("tpcID_goodTree", &pTPCID, &tpcID_branch);
        t->SetBranchAddress("uncorrectedDeltaX_goodTree", &pUncorrectedDeltaX, &uncorrectedDeltaX_branch);
        t->SetBranchAddress("uncorrectedDeltaR_goodTree", &pUncorrectedDeltaR, &uncorrectedDeltaR_branch);
        t->SetBranchAddress("recoVertexX_goodTree", &pRecoVertexX, &recoVertexX_branch);
        t->SetBranchAddress("recoVertexY_goodTree", &pRecoVertexY, &recoVertexY_branch);
        t->SetBranchAddress("recoVertexZ_goodTree", &pRecoVertexZ, &recoVertexZ_branch);
        t->SetBranchAddress("mcVertexX_goodTree", &pMCVertexX, &mcVertexX_branch);
        t->SetBranchAddress("mcVertexY_goodTree", &pMCVertexY, &mcVertexY_branch);
        t->SetBranchAddress("mcVertexZ_goodTree", &pMCVertexZ, &mcVertexZ_branch);

	    deltaX_branch->GetEntry(i);
        deltaY_branch->GetEntry(i);
        deltaZ_branch->GetEntry(i);
        deltaR_branch->GetEntry(i);    
        CCNC_branch->GetEntry(i);
        trueNeutrino_branch->GetEntry(i);
        dl_current_branch->GetEntry(i);
	    cosmic_single_branch->GetEntry(i);
        pass1pass2_branch->GetEntry(i);
        eventID_branch->GetEntry(i);
        runID_branch->GetEntry(i);
        subRunID_branch->GetEntry(i);
        tpcID_branch->GetEntry(i);
        uncorrectedDeltaX_branch->GetEntry(i);
        uncorrectedDeltaR_branch->GetEntry(i);
        recoVertexX_branch->GetEntry(i);
        recoVertexY_branch->GetEntry(i);
        recoVertexZ_branch->GetEntry(i);
        mcVertexX_branch->GetEntry(i);
        mcVertexY_branch->GetEntry(i);
        mcVertexZ_branch->GetEntry(i);

	    size_t vecSize = pDeltaX->size();
	    for(UInt_t j = 0; j < vecSize; j++){
		    allDeltaX.push_back(pDeltaX->at(j));
		    allDeltaY.push_back(pDeltaY->at(j));
		    allDeltaZ.push_back(pDeltaZ->at(j));
            allDeltaR.push_back(pDeltaR->at(j));
          	allCCNC.push_back(pCCNC->at(j));
          	allTrueNeutrino.push_back(pTrueNeutrino->at(j));
          	allDL_Current.push_back(pDL_Current->at(j));       
	    	allCosmic_Single.push_back(pCosmic_Single->at(j));
            allPass1Pass2.push_back(pPass1Pass2->at(j));
            allEventID.push_back(pEventID->at(j));
            allRunID.push_back(pRunID->at(j));
            allSubRunID.push_back(pSubRunID->at(j));
            allFileNumber.push_back(fileNumber);
            allTPCID.push_back(pTPCID->at(j));
            allUncorrectedDeltaX.push_back(pUncorrectedDeltaX->at(j));
            allUncorrectedDeltaR.push_back(pUncorrectedDeltaR->at(j));
            allRecoVertexX.push_back(pRecoVertexX->at(j));
            allRecoVertexY.push_back(pRecoVertexY->at(j));
            allRecoVertexZ.push_back(pRecoVertexZ->at(j));
            allMCVertexX.push_back(pMCVertexX->at(j));
            allMCVertexY.push_back(pMCVertexY->at(j));
            allMCVertexZ.push_back(pMCVertexZ->at(j));
            //if(fileNumber <= 30){
	           // allFileNumber.push_back(fileNumber);
            //} else if((fileNumber > 30) && (fileNumber <= 60)){
                //allFileNumber.push_back((fileNumber-30));
            //} else if(fileNumber > 60){
                //allFileNumber.push_back((fileNumber-60));
            //}
        }

        fileNumber++;
    }

    std::vector<vertex_t> allVertices = std::vector<vertex_t>();

    size_t newNumEvents = allDeltaX.size();

    for(UInt_t k = 0; k < newNumEvents; k++){
        vertex_t vertex;
	    vertex.deltaX = allDeltaX.at(k);
	    vertex.deltaY = allDeltaY.at(k);
        vertex.deltaZ = allDeltaZ.at(k);
        vertex.deltaR = allDeltaR.at(k);
        vertex.trueNeutrino = allTrueNeutrino.at(k);
        vertex.CCNC = allCCNC.at(k);
        vertex.DLCurrent = allDL_Current.at(k);
        vertex.CosmicSingle = allCosmic_Single.at(k);
        vertex.Pass1Pass2 = allPass1Pass2.at(k);
        vertex.eventID = allEventID.at(k);
        vertex.runID = allRunID.at(k);
        vertex.subRunID = allSubRunID.at(k);
        vertex.fileNumber = allFileNumber.at(k);
        vertex.tpcID = allTPCID.at(k);
        vertex.uncorrectedDeltaX = allUncorrectedDeltaX.at(k);
	    vertex.uncorrectedDeltaR = allUncorrectedDeltaR.at(k);
        vertex.recoX = allRecoVertexX.at(k);
        vertex.recoY = allRecoVertexY.at(k);
        vertex.recoZ = allRecoVertexZ.at(k);
        vertex.trueX = allMCVertexX.at(k);
        vertex.trueY = allMCVertexY.at(k);
        vertex.trueZ = allMCVertexZ.at(k);
        allVertices.push_back(vertex);
    } 

    // for(int h = 0; h < 1; h++){
        // std::cout << allVertices.at(h).trueX << "," << allVertices.at(h).trueY << "," << allVertices.at(h).trueZ << std::endl;
    // }
    // For getting individual delta X, Y, Z, R distributions:
    //cosmicVertexHist(allVertices);
    //
    // For getting distributions (current, DL uboone and dune) plotted together:
    plottingTogetherDist(allVertices);
    //
    // For getting ROC curves (current, DL uboone and dune) plotted together
    //rocTogether(allVertices);
    
    // For calculating percentage of events with reco vertex with a certain distance of true vertex
    //calculations(allVertices);
    //gridComparison(allVertices);
    //tpcPlots(allVertices);
    //pass1pass2(allVertices);

}
