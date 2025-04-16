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
    double pdg;
    double isPrimary;
    double vx;
    double vy;
    double vz;
    double dx;
    double dy;
    double dz;
    double sliceID;
    double bestPlaneEnergy;
    double theta;
    double trackscore;
    double completeness;
    double purity;
} recoParticle;

typedef struct{
    double pdg;
    double isPrimary;
    double vx;
    double vy;
    double vz;
    double sliceID;
} recoNeutrino;

typedef struct{
    double vx;
    double vy;
    double vz;
    double CCNC;
    double pdg;
    double leptonpdg;
    double tpcID;
    double tpcValid;
} trueNeutrino;

typedef struct{
    double pdg;
    double vx;
    double vy;
    double vz;
    double px;
    double py;
    double pz;
    double energy;
    double angle;
    double ETheta2;
    double dx;
    double dy;
    double dz;
} trueParticle;

typedef struct{
    double id;
    double completeness;
    double purity;
    double score;
} recoSlice;

typedef struct{
    std::vector<recoParticle> recoParticleVec;
    std::vector<recoNeutrino> recoNeutrinoVec;
    std::vector<trueNeutrino> trueNeutrinoVec;
    std::vector<trueParticle> trueParticleVec;
    std::vector<recoSlice>    sliceVec;
    double DLCurrent;
} eventData;

typedef struct{
    recoParticle chosenRecoParticle;
    recoNeutrino chosenRecoNeutrino;
    trueNeutrino chosenTrueNeutrino;
    trueParticle chosenTrueParticle;
    recoSlice chosenSlice;
    double totalPFPEnergy;
    double DLCurrent;
} eventAfterCuts;

void styleDraw(TCanvas* canvas, TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, TPaveText* pt = nullptr, int* percentage = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    current->SetLineWidth(2);
    current->SetLineColor(kRed);

    cheated->SetLineWidth(2);
    cheated->SetLineColor(kSpring-5);

    dune->SetLineWidth(2);
    dune->SetLineColor(kViolet-5);

    uboone->SetLineWidth(2);
    uboone->SetLineColor(kBlue);

    current->Draw("hist");
    cheated->Draw("histsame");
    dune->Draw("histsame");
    uboone->Draw("histsame");

    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    current->SetStats(0);
    current->GetXaxis()->SetTickLength(0.04);
    current->GetYaxis()->SetTickLength(0.03);
    current->GetXaxis()->SetTickSize(0.02);
    current->GetYaxis()->SetTickSize(0.02);

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(dune, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(uboone, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(current, "Current SBND Vertexing (without Refinement)", "f");
    legend->AddEntry(cheated, "Cheated SBND Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(pt){
        pt->SetTextSize(legend->GetTextSize());
        pt->SetTextFont(legend->GetTextFont());
        pt->Draw();
    }

    canvas->SaveAs(filename);
} 

void percentage(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax){
    TCanvas *percentageCanvas = new TCanvas("percentage_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentPerc = (TH1F*) current->Clone("perc hist");
    currentPerc->Scale(100.0 * 1.0/sizeCurrent);
    currentPerc->GetYaxis()->SetTitle("Percentage of Events (%)"); 

    TH1F* cheatedPerc = (TH1F*) cheated->Clone("perc hist");
    cheatedPerc->Scale(100.0 * 1.0/sizeCheated);
    cheatedPerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* dunePerc = (TH1F*) dune->Clone("perc hist");
    dunePerc->Scale(100.0 * 1.0/sizeDune);
    dunePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* uboonePerc = (TH1F*) uboone->Clone("perc hist");
    uboonePerc->Scale(100.0 * 1.0/sizeUboone);
    uboonePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 
   
    int funcValue = 1;

    styleDraw(percentageCanvas, currentPerc, cheatedPerc, dunePerc, uboonePerc, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt, &funcValue);
}

recoSlice chooseRecoSlice(std::vector<recoSlice> recoSliceVec){
    int chosenSliceIndex;
    double highestCompleteness = -100;

    for(size_t j = 0; j < recoSliceVec.size(); ++j){
        if(recoSliceVec[j].completeness > highestCompleteness){
            highestCompleteness = recoSliceVec[j].completeness;
            chosenSliceIndex = j;
        }
    }
    
    return recoSliceVec[chosenSliceIndex];
}

recoParticle chooseRecoParticle(std::vector<recoParticle> recoParticleVec, double& totalEnergy){
   totalEnergy = 0; 
   int chosenParticleIndex;
   double highestEnergy = -10000;

    for(size_t j = 0; j < recoParticleVec.size(); ++j){
        totalEnergy += recoParticleVec[j].bestPlaneEnergy;
        if(recoParticleVec[j].bestPlaneEnergy > highestEnergy){
            highestEnergy = recoParticleVec[j].bestPlaneEnergy;
            chosenParticleIndex = j;
        }
    }

    return recoParticleVec[chosenParticleIndex];
}

void numberPlots(std::vector<eventData> allEventData){
    TCanvas *numSlicesCanvas = new TCanvas("numSlices_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numSlicesCurrent_dist = new TH1F("Number of Slices", "Number of Slices", 5, 0, 5);
    numSlicesCurrent_dist->SetTitle("Number of Slices in an Event;Number of Slices;# of Events");
    TH1F* numSlicesDLUboone_dist = (TH1F*) numSlicesCurrent_dist->Clone("numSlices");
    TH1F* numSlicesDLDune_dist = (TH1F*) numSlicesCurrent_dist->Clone("numSlices");
    TH1F* numSlicesCheated_dist = (TH1F*) numSlicesCurrent_dist->Clone("numSlices");  

    TCanvas *numPFPsCanvas = new TCanvas("numPFPs_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numPFPsCurrent_dist = new TH1F("Num PFPs", "Num PFPs", 10, 0, 10);
    numPFPsCurrent_dist->SetTitle("Number of PFPs in an Event;Number of PFPs;# of Events");
    TH1F* numPFPsDLUboone_dist = (TH1F*) numPFPsCurrent_dist->Clone("numPFPs");
    TH1F* numPFPsDLDune_dist = (TH1F*) numPFPsCurrent_dist->Clone("numPFPs");
    TH1F* numPFPsCheated_dist = (TH1F*) numPFPsCurrent_dist->Clone("numPFPs"); 

    TCanvas *numPFPs1SliceCanvas = new TCanvas("numPFPs1Slice_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numPFPs1SliceCurrent_dist = new TH1F("Num PFPs 1 Slice", "Num PFPs 1 Slice", 10, 0, 10);
    numPFPs1SliceCurrent_dist->SetTitle("Number of PFPs in an Event with 1 Slice;Number of PFPs;# of Events");
    TH1F* numPFPs1SliceDLUboone_dist = (TH1F*) numPFPs1SliceCurrent_dist->Clone("numPFPs");
    TH1F* numPFPs1SliceDLDune_dist = (TH1F*) numPFPs1SliceCurrent_dist->Clone("numPFPs");
    TH1F* numPFPs1SliceCheated_dist = (TH1F*) numPFPs1SliceCurrent_dist->Clone("numPFPs"); 
    
    TCanvas *numPFPsSlicesCanvas = new TCanvas("numPFPsSlices_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numPFPsSlicesCurrent_dist = new TH1F("Num PFPs Slices", "Num PFPs Slices", 10, 0, 10);
    numPFPsSlicesCurrent_dist->SetTitle("Number of PFPs in an Event with > 1 Slice;Number of PFPs;# of Events");
    TH1F* numPFPsSlicesDLUboone_dist = (TH1F*) numPFPsSlicesCurrent_dist->Clone("numPFPs");
    TH1F* numPFPsSlicesDLDune_dist = (TH1F*) numPFPsSlicesCurrent_dist->Clone("numPFPs");
    TH1F* numPFPsSlicesCheated_dist = (TH1F*) numPFPsSlicesCurrent_dist->Clone("numPFPs"); 
    
    double sizeCurrent = 0;
    double sizeCheated = 0;
    double sizeDLDune = 0;
    double sizeDLUboone = 0;
    double sizeCurrent1Slice = 0;
    double sizeCheated1Slice = 0;
    double sizeDLDune1Slice = 0;
    double sizeDLUboone1Slice = 0;

    for(const auto& event : allEventData){
        if(!(event.trueNeutrinoVec[0].pdg == -999999 || event.trueNeutrinoVec[0].tpcValid == 0)){
            size_t recoParticleCount = event.recoParticleVec.size();
            size_t recoSliceCount = event.sliceVec.size();

            if(recoParticleCount == 1 && event.recoParticleVec[0].pdg == -999999) {
                recoParticleCount = 0;
            } 
        
            if(recoSliceCount == 1 && event.sliceVec[0].id == -999999){
                recoSliceCount = 0;
            }

            if(event.DLCurrent == 0){
                sizeDLUboone++;
                numPFPsDLUboone_dist->Fill(recoParticleCount);
                numSlicesDLUboone_dist->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1SliceDLUboone_dist->Fill(recoParticleCount);
                    sizeDLUboone1Slice++;
                } else if(recoSliceCount > 1) numPFPsSlicesDLUboone_dist->Fill(recoParticleCount);
               
            } else if(event.DLCurrent == 1){
                sizeDLDune++;
                numPFPsDLDune_dist->Fill(recoParticleCount);
                numSlicesDLDune_dist->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1SliceDLDune_dist->Fill(recoParticleCount);
                    sizeDLDune1Slice++;
                } else if(recoSliceCount > 1) numPFPsSlicesDLDune_dist->Fill(recoParticleCount);
            
            } else if(event.DLCurrent == 2){
                sizeCurrent++;
                numPFPsCurrent_dist->Fill(recoParticleCount);
                numSlicesCurrent_dist->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1SliceCurrent_dist->Fill(recoParticleCount);
                    sizeCurrent1Slice++;
                } else if(recoSliceCount > 1) numPFPsSlicesCurrent_dist->Fill(recoParticleCount); 
            } else if(event.DLCurrent == 3){
                sizeCheated++;
                numPFPsCheated_dist->Fill(recoParticleCount);
                numSlicesCheated_dist->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1SliceCheated_dist->Fill(recoParticleCount);
                    sizeCheated1Slice++;
                } else if(recoSliceCount > 1) numPFPsSlicesCheated_dist->Fill(recoParticleCount);
            }
        }
    }   


    styleDraw(numSlicesCanvas, numSlicesCurrent_dist, numSlicesCheated_dist, numSlicesDLDune_dist, numSlicesDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/numSlices_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(numPFPsCanvas, numPFPsCurrent_dist, numPFPsCheated_dist, numPFPsDLDune_dist, numPFPsDLUboone_dist, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(numPFPs1SliceCanvas, numPFPs1SliceCurrent_dist, numPFPs1SliceCheated_dist, numPFPs1SliceDLDune_dist, numPFPs1SliceDLUboone_dist, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs1Slice_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(numPFPsSlicesCanvas, numPFPsSlicesCurrent_dist, numPFPsSlicesCheated_dist, numPFPsSlicesDLDune_dist, numPFPsSlicesDLUboone_dist, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPsSlices_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    
    percentage(numSlicesCurrent_dist, numSlicesCheated_dist, numSlicesDLDune_dist, numSlicesDLUboone_dist, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/numSlices_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(numPFPsCurrent_dist, numPFPsCheated_dist, numPFPsDLDune_dist, numPFPsDLUboone_dist, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 90, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(numPFPs1SliceCurrent_dist, numPFPs1SliceCheated_dist, numPFPs1SliceDLDune_dist, numPFPs1SliceDLUboone_dist, sizeCurrent1Slice, sizeCheated1Slice, sizeDLDune1Slice, sizeDLUboone1Slice, 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs1Slice_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(numPFPsSlicesCurrent_dist, numPFPsSlicesCheated_dist, numPFPsSlicesDLDune_dist, numPFPsSlicesDLUboone_dist, sizeCurrent - sizeCurrent1Slice, sizeCheated - sizeCheated1Slice, sizeDLDune - sizeDLDune1Slice, sizeDLUboone - sizeDLUboone1Slice, 0, 90, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPsSlices_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    printf("Number of events: Current = %f, Cheated %f, DL Dune = %f, DL Uboone = %f\n", sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone);    
}

void chosenSlicePlots(std::vector<eventAfterCuts> allEvents){


}

void nuE_macroNew(){

    std::vector<eventData> allEventData;
    std::vector<eventAfterCuts> allCutEventData;

    int counterTrueNeut = 0;
    int counterTruePart = 0;
    int counterRecoNeut = 0;
    int counterRecoPart = 0;

    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/merged.root");
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

    std::vector<double> *truth_neutrinoVX = nullptr;
    std::vector<double> *truth_neutrinoVY = nullptr;
    std::vector<double> *truth_neutrinoVZ = nullptr;
    std::vector<double> *truth_CCNC = nullptr;
    std::vector<double> *truth_neutrinoType = nullptr;
    std::vector<double> *truth_chargedLepton = nullptr;
    std::vector<double> *truth_neutrinoTPCID = nullptr;
    std::vector<double> *truth_neutrinoTPCValid = nullptr;

    std::vector<double> *truth_particlePDG = nullptr;
    std::vector<double> *truth_particleVX = nullptr;
    std::vector<double> *truth_particleVY = nullptr;
    std::vector<double> *truth_particleVZ = nullptr;
    std::vector<double> *truth_particlePX = nullptr;
    std::vector<double> *truth_particlePY = nullptr;
    std::vector<double> *truth_particlePZ = nullptr;
    std::vector<double> *truth_particleEnergy = nullptr;
    std::vector<double> *truth_particleAngle = nullptr;
    std::vector<double> *truth_particleETheta2 = nullptr;
    std::vector<double> *truth_particleDirectionX = nullptr;
    std::vector<double> *truth_particleDirectionY = nullptr;
    std::vector<double> *truth_particleDirectionZ = nullptr;

    std::vector<double> *reco_neutrinoPDG = nullptr;
    std::vector<double> *reco_neutrinoIsPrimary = nullptr;
    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;

    std::vector<double> *reco_particlePDG = nullptr;
    std::vector<double> *reco_particleIsPrimary = nullptr;
    std::vector<double> *reco_particleVX = nullptr;
    std::vector<double> *reco_particleVY = nullptr;
    std::vector<double> *reco_particleVZ = nullptr;
    std::vector<double> *reco_particleDirectionX = nullptr;
    std::vector<double> *reco_particleDirectionY = nullptr;
    std::vector<double> *reco_particleDirectionZ = nullptr;
    std::vector<double> *reco_particleSliceID = nullptr;
    std::vector<double> *reco_particleBestPlaneEnergy = nullptr;
    std::vector<double> *reco_particleTheta = nullptr;
    std::vector<double> *reco_particleTrackScore = nullptr;
    std::vector<double> *reco_particleCompleteness = nullptr;
    std::vector<double> *reco_particlePurity = nullptr;

    std::vector<double> *reco_sliceID = nullptr;
    std::vector<double> *reco_sliceCompleteness = nullptr;
    std::vector<double> *reco_slicePurity = nullptr;
    std::vector<double> *reco_sliceScore = nullptr;

    tree->SetBranchAddress("truth_neutrinoVX", &truth_neutrinoVX);
    tree->SetBranchAddress("truth_neutrinoVY", &truth_neutrinoVY);
    tree->SetBranchAddress("truth_neutrinoVZ", &truth_neutrinoVZ);
    tree->SetBranchAddress("truth_CCNC", &truth_CCNC);
    tree->SetBranchAddress("truth_neutrinoType", &truth_neutrinoType);
    tree->SetBranchAddress("truth_chargedLepton", &truth_chargedLepton);
    tree->SetBranchAddress("truth_neutrinoTPCID", &truth_neutrinoTPCID);
    tree->SetBranchAddress("truth_neutrinoTPCValid", &truth_neutrinoTPCValid);

    tree->SetBranchAddress("truth_particlePDG", &truth_particlePDG);
    tree->SetBranchAddress("truth_particleVX", &truth_particleVX);
    tree->SetBranchAddress("truth_particleVY", &truth_particleVY);
    tree->SetBranchAddress("truth_particleVZ", &truth_particleVZ);
    tree->SetBranchAddress("truth_particlePX", &truth_particlePX);
    tree->SetBranchAddress("truth_particlePY", &truth_particlePY);
    tree->SetBranchAddress("truth_particlePZ", &truth_particlePZ);
    tree->SetBranchAddress("truth_particleEnergy", &truth_particleEnergy);
    tree->SetBranchAddress("truth_particleAngle", &truth_particleAngle);
    tree->SetBranchAddress("truth_particleETheta2", &truth_particleETheta2);
    tree->SetBranchAddress("truth_particleDirectionX", &truth_particleDirectionX);
    tree->SetBranchAddress("truth_particleDirectionY", &truth_particleDirectionY);
    tree->SetBranchAddress("truth_particleDirectionZ", &truth_particleDirectionZ);

    tree->SetBranchAddress("reco_neutrinoPDG", &reco_neutrinoPDG);
    tree->SetBranchAddress("reco_neutrinoIsPrimary", &reco_neutrinoIsPrimary);
    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);

    tree->SetBranchAddress("reco_particlePDG", &reco_particlePDG);
    tree->SetBranchAddress("reco_particleIsPrimary", &reco_particleIsPrimary);
    tree->SetBranchAddress("reco_particleVX", &reco_particleVX);
    tree->SetBranchAddress("reco_particleVY", &reco_particleVY);
    tree->SetBranchAddress("reco_particleVZ", &reco_particleVZ);
    tree->SetBranchAddress("reco_particleDirectionX", &reco_particleDirectionX);
    tree->SetBranchAddress("reco_particleDirectionY", &reco_particleDirectionY);
    tree->SetBranchAddress("reco_particleDirectionZ", &reco_particleDirectionZ);
    tree->SetBranchAddress("reco_particleSliceID", &reco_particleSliceID);
    tree->SetBranchAddress("reco_particleBestPlaneEnergy", &reco_particleBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTheta", &reco_particleTheta);
    tree->SetBranchAddress("reco_particleTrackScore", &reco_particleTrackScore);
    tree->SetBranchAddress("reco_particleCompleteness", &reco_particleCompleteness);
    tree->SetBranchAddress("reco_particlePurity", &reco_particlePurity);

    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness", &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_slicePurity", &reco_slicePurity);
    tree->SetBranchAddress("reco_sliceScore", &reco_sliceScore);

    Long64_t numEntries = tree->GetEntries();
    std::cout << "Number of Events: " << numEntries << std::endl;
    //tree->Print();

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);
        //printf("\n_________________________________________________________________________________\n");
        //printf("Run: %d, SubRun: %d, Event: %d, DL/Current: %f\n", runID, subRunID, eventID, DLCurrent);
        
        eventData eData; // A struct containing vectors for each type of particle, i.e. true or reco...
        eventData eDataAfterCuts;

        std::vector<recoParticle> recoParticles;
        std::vector<trueParticle> trueParticles;
        std::vector<recoNeutrino> recoNeutrinos;
        std::vector<trueNeutrino> trueNeutrinos;
        std::vector<recoSlice> slices;
        //printf("Number of true neutrinos = %zu, reco neutrinos = %zu, true particles = %zu, reco particles = %zu, slices = %zu\n\n", truth_neutrinoVX->size(), reco_neutrinoPDG->size(), truth_particlePDG->size(), reco_particlePDG->size(), reco_sliceID->size());
    
        if(truth_neutrinoVX->size() > 1) counterTrueNeut++;
        if(truth_particlePDG->size() > 1) counterTruePart++;
        if(reco_neutrinoPDG->size() < 1) counterRecoNeut++;
        if(reco_particlePDG->size() > 1) counterRecoPart++;

        for(size_t j = 0; j < reco_particlePDG->size(); ++j){
            recoParticle particle;
            particle.pdg = reco_particlePDG->at(j);
            particle.isPrimary = reco_particleIsPrimary->at(j);
            particle.vx = reco_particleVX->at(j);
            particle.vy = reco_particleVY->at(j);
            particle.vz = reco_particleVZ->at(j);
            particle.dx = reco_particleDirectionX->at(j);
            particle.dy = reco_particleDirectionY->at(j);
            particle.dz = reco_particleDirectionZ->at(j);
            particle.sliceID = reco_particleSliceID->at(j);
            particle.bestPlaneEnergy = reco_particleBestPlaneEnergy->at(j);
            particle.theta = reco_particleTheta->at(j);
            particle.trackscore = reco_particleTrackScore->at(j);
            particle.completeness = reco_particleCompleteness->at(j);
            particle.purity = reco_particlePurity->at(j);
            recoParticles.push_back(particle);
            //printf("Reco Particle %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %f, Energy = %f, Theta = %f, Trackscore = %f, Completeness = %f, Purity = %f\n", j, particle.pdg, particle.isPrimary, particle.vx, particle.vy, particle.vz, particle.dx, particle.dy, particle.dz, particle.sliceID, particle.bestPlaneEnergy, particle.theta, particle.trackscore, particle.completeness, particle.purity);
        }

        //printf("\n");
        for(size_t j = 0; j < reco_neutrinoPDG->size(); ++j){
            recoNeutrino neutrino;
            neutrino.pdg = reco_neutrinoPDG->at(j);
            neutrino.isPrimary = reco_neutrinoIsPrimary->at(j);
            neutrino.vx = reco_neutrinoVX->at(j);
            neutrino.vy = reco_neutrinoVY->at(j);
            neutrino.vz = reco_neutrinoVZ->at(j);
            neutrino.sliceID = reco_neutrinoSliceID->at(j);
            recoNeutrinos.push_back(neutrino);
            //printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, neutrino.pdg, neutrino.isPrimary, neutrino.vx, neutrino.vy, neutrino.vz, neutrino.sliceID);
        }

        //printf("\n");
        for(size_t j = 0; j < truth_particlePDG->size(); ++j){
            trueParticle particle;
            particle.pdg = truth_particlePDG->at(j);
            particle.vx = truth_particleVX->at(j);
            particle.vy = truth_particleVY->at(j);
            particle.vz = truth_particleVZ->at(j);
            particle.px = truth_particlePX->at(j);
            particle.py = truth_particlePY->at(j);
            particle.pz = truth_particlePZ->at(j);
            particle.energy = truth_particleEnergy->at(j);
            particle.angle = truth_particleAngle->at(j);
            particle.ETheta2 = truth_particleETheta2->at(j);
            particle.dx = truth_particleDirectionX->at(j);
            particle.dy = truth_particleDirectionY->at(j);
            particle.dz = truth_particleDirectionZ->at(j);
            trueParticles.push_back(particle);
            //printf("Truth Particle %zu: PDG = %f, Vertex = (%f, %f, %f), Momentum = (%f, %f, %f), Energy = %f, Theta = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", j, particle.pdg, particle.vx, particle.vy, particle.vz, particle.px, particle.py, particle.pz, particle.energy, particle.angle, particle.ETheta2, particle.dx, particle.dy, particle.dz);
        }
    
        //printf("\n");
        for(size_t j = 0; j < truth_neutrinoVX->size(); ++j){
            trueNeutrino neutrino;
            neutrino.vx = truth_neutrinoVX->at(j);
            neutrino.vy = truth_neutrinoVY->at(j);
            neutrino.vz = truth_neutrinoVZ->at(j);
            neutrino.CCNC = truth_CCNC->at(j);
            neutrino.pdg = truth_neutrinoType->at(j);
            neutrino.leptonpdg = truth_chargedLepton->at(j);
            neutrino.tpcID = truth_neutrinoTPCID->at(j);
            neutrino.tpcValid = truth_neutrinoTPCValid->at(j);
            trueNeutrinos.push_back(neutrino);
            //printf("Truth Neutrino %zu: PDG = %f, Vertex = (%f, %f, %f), CCNC = %f, Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", j, neutrino.pdg, neutrino.vx, neutrino.vy, neutrino.vz, neutrino.CCNC, neutrino.leptonpdg, neutrino.tpcID, neutrino.tpcValid);
        }
    
        //printf("\n");
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            recoSlice reconstructedSlice;
            reconstructedSlice.id = reco_sliceID->at(j);
            reconstructedSlice.completeness = reco_sliceCompleteness->at(j);
            reconstructedSlice.purity = reco_slicePurity->at(j);
            reconstructedSlice.score = reco_sliceScore->at(j);
            slices.push_back(reconstructedSlice);
            //printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, Score = %f\n", j, reconstructedSlice.id, reconstructedSlice.completeness, reconstructedSlice.purity, reconstructedSlice.score);
        }

        //printf("\nChosen:\n");
        eData.recoParticleVec = recoParticles;
        eData.recoNeutrinoVec = recoNeutrinos;
        eData.trueNeutrinoVec = trueNeutrinos;
        eData.trueParticleVec = trueParticles;
        eData.sliceVec = slices;

        eventAfterCuts eventCut;

        double totalEnergyPFPs = 0;

        if(recoParticles.size() == 1){
            eventCut.chosenRecoParticle = recoParticles[0];     
            totalEnergyPFPs = recoParticles[0].bestPlaneEnergy;
        } else{
            eventCut.chosenRecoParticle = chooseRecoParticle(recoParticles, totalEnergyPFPs);
        }
        eventCut.totalPFPEnergy = totalEnergyPFPs;

        if(recoNeutrinos.size() == 1){
            eventCut.chosenRecoNeutrino = recoNeutrinos[0];
        } else{
        }

        if(trueNeutrinos.size() == 1){
            eventCut.chosenTrueNeutrino = trueNeutrinos[0];
        }
        
        if(trueParticles.size() == 1){
            eventCut.chosenTrueParticle = trueParticles[0];
        }

        if(slices.size() == 1){
            eventCut.chosenSlice = slices[0];
        } else{
            eventCut.chosenSlice = chooseRecoSlice(slices);
        }

        eventCut.DLCurrent = DLCurrent;
        eData.DLCurrent = DLCurrent;

        //printf("Chosen Reco Particle: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %f, Energy = %f, Theta = %f, Trackscore = %f, Completeness = %f, Purity = %f\n", eventCut.chosenRecoParticle.pdg, eventCut.chosenRecoParticle.isPrimary, eventCut.chosenRecoParticle.vx, eventCut.chosenRecoParticle.vy, eventCut.chosenRecoParticle.vz, eventCut.chosenRecoParticle.dx, eventCut.chosenRecoParticle.dy, eventCut.chosenRecoParticle.dz, eventCut.chosenRecoParticle.sliceID, eventCut.chosenRecoParticle.bestPlaneEnergy, eventCut.chosenRecoParticle.theta, eventCut.chosenRecoParticle.trackscore, eventCut.chosenRecoParticle.completeness, eventCut.chosenRecoParticle.purity);
        //printf("Chosen Reco Neutrino: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", eventCut.chosenRecoNeutrino.pdg, eventCut.chosenRecoNeutrino.isPrimary, eventCut.chosenRecoNeutrino.vx, eventCut.chosenRecoNeutrino.vy, eventCut.chosenRecoNeutrino.vz, eventCut.chosenRecoNeutrino.sliceID);
        //printf("Chosen True Neutrino: PDG = %f, Vertex = (%f, %f, %f), CCNC = %f, Charged Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", eventCut.chosenTrueNeutrino.pdg, eventCut.chosenTrueNeutrino.vx, eventCut.chosenTrueNeutrino.vy, eventCut.chosenTrueNeutrino.vz, eventCut.chosenTrueNeutrino.CCNC, eventCut.chosenTrueNeutrino.leptonpdg, eventCut.chosenTrueNeutrino.tpcID, eventCut.chosenTrueNeutrino.tpcValid);
        //printf("Chosen True Particle: PDG = %f, Vertex = (%f, %f, %f), Momentum = (%f, %f, %f), Energy = %f, Theta = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", eventCut.chosenTrueParticle.pdg, eventCut.chosenTrueParticle.vx, eventCut.chosenTrueParticle.vy, eventCut.chosenTrueParticle.vz, eventCut.chosenTrueParticle.px, eventCut.chosenTrueParticle.py, eventCut.chosenTrueParticle.pz, eventCut.chosenTrueParticle.energy, eventCut.chosenTrueParticle.angle, eventCut.chosenTrueParticle.ETheta2, eventCut.chosenTrueParticle.dx, eventCut.chosenTrueParticle.dy, eventCut.chosenTrueParticle.dz);
        //printf("Chosen Slice: ID = %f, Completeness = %f, Purity = %f, Score = %f\n", eventCut.chosenSlice.id, eventCut.chosenSlice.completeness, eventCut.chosenSlice.purity, eventCut.chosenSlice.score);
        //printf("Total Reconstructed Energy = %f\n", totalEnergyPFPs);

        recoParticles.clear();
        recoNeutrinos.clear();
        trueNeutrinos.clear();
        trueParticles.clear();
        slices.clear();

        allEventData.push_back(eData);
        allCutEventData.push_back(eventCut);
    }

    numberPlots(allEventData);
    // allCutEventData
    file->Close();
    printf("\n __________ Number of events with > 1 True Neutrino = %d, True Particle = %d, Reco Neutrino = %d, Reco Particle = %d __________\n", counterTrueNeut, counterTruePart, counterRecoNeut, counterRecoPart);
}
