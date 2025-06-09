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

struct counter{
    int current = 0;
    int cheated = 0;
    int dune = 0;
    int uboone = 0;
    int sbnd = 0;
};

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* current;
    TH1F* cheated;
    TH1F* dune;
    TH1F* uboone;
    TH1F* sbnd;
} histGroup;


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
    double neutrinoParent;
} trueParticle;

typedef struct{
    double id;
    double completeness;
    double purity;
    double score;
} recoSlice;

histGroup createHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);
    
    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_current").c_str()),
        (TH1F*) base->Clone((baseName + "_cheated").c_str()),
        (TH1F*) base->Clone((baseName + "_dldune").c_str()),
        (TH1F*) base->Clone((baseName + "_dluboone").c_str()),
        (TH1F*) base->Clone((baseName + "_dlsbnd").c_str())
    };    
}


void styleDraw(TCanvas* canvas, TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, TH1F* sbnd, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, TPaveText* pt = nullptr, int* percentage = nullptr, int* drawLine = nullptr, int* linePos = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    gStyle->SetPalette(kAvocado);
    gROOT->ForceStyle();
    gPad->Update();

    current->SetLineWidth(2);
    //current->SetLineColor(TColor::GetColorPalette(150));
    current->SetLineColor(TColor::GetColor("#e42536"));

    cheated->SetLineWidth(2);
    //cheated->SetLineColor(TColor::GetColorPalette(200));
    cheated->SetLineColor(TColor::GetColor("#f89c20"));

    dune->SetLineWidth(2);
    //dune->SetLineColor(TColor::GetColorPalette(50));
    dune->SetLineColor(TColor::GetColor("#7a21dd"));

    uboone->SetLineWidth(2);
    //uboone->SetLineColor(TColor::GetColorPalette(100));
    uboone->SetLineColor(TColor::GetColor("#5790fc"));

    sbnd->SetLineWidth(2);
    sbnd->SetLineColor(TColor::GetColor("#9c9ca1"));

    current->Draw("hist");
    cheated->Draw("histsame");
    dune->Draw("histsame");
    uboone->Draw("histsame");
    sbnd->Draw("histsame");

    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    current->SetStats(0);
    current->GetXaxis()->SetTickLength(0.04);
    current->GetYaxis()->SetTickLength(0.03);
    current->GetXaxis()->SetTickSize(0.02);
    current->GetYaxis()->SetTickSize(0.02);

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(dune, "Pandora Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(uboone, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(sbnd, "Pandora Deep Learning: SBND/BNB Tune", "f");
    legend->AddEntry(current, "Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(cheated, "Pandora Cheated SBND Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, current->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, current->GetMaximum() * 0.93, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.1, current->GetMaximum() * 0.93, "2m_{e}");
        }

        latex->SetTextSize(0.035); 
        latex->SetTextAlign(11);
        latex->Draw("same");
    }

    if(pt){
        pt->SetTextSize(legend->GetTextSize());
        pt->SetTextFont(legend->GetTextFont());
        pt->Draw();
    }

    canvas->SaveAs(filename);
}

void percentage(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, TH1F* sbnd, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double sizeSBND, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr, int* linePos = nullptr){
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

    TH1F* sbndPerc = (TH1F*) sbnd->Clone("perc hist");
    sbndPerc->Scale(100.0 * 1.0/sizeSBND);
    sbndPerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.17, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of DL SBND Entries: %d", (int)sizeSBND));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 
   
    int funcValue = 1;

    styleDraw(percentageCanvas, currentPerc, cheatedPerc, dunePerc, uboonePerc, sbndPerc, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt, &funcValue, drawLine, linePos);
}

void efficiency(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, TH1F* sbnd, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double sizeSBND, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr, int* linePos = nullptr, std::string xlabel = ""){
    TCanvas *efficiencyCanvas = new TCanvas("efficiency_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentEff = (TH1F*) current->Clone("eff hist");
    currentEff->Reset();
    currentEff->GetYaxis()->SetTitle("Efficiency"); 
    currentEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* cheatedEff = (TH1F*) cheated->Clone("eff hist");
    cheatedEff->Reset();
    cheatedEff->GetYaxis()->SetTitle("Efficiency");
    cheatedEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* duneEff = (TH1F*) dune->Clone("eff hist");
    duneEff->Reset();
    duneEff->GetYaxis()->SetTitle("Efficiency");
    duneEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneEff = (TH1F*) uboone->Clone("eff hist");
    ubooneEff->Reset();
    ubooneEff->GetYaxis()->SetTitle("Efficiency");
    ubooneEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* sbndEff = (TH1F*) sbnd->Clone("eff hist");
    sbndEff->Reset();
    sbndEff->GetYaxis()->SetTitle("Efficiency");
    sbndEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    int numBins = current->GetNbinsX();
    double currentSum = 0.0;
    double cheatedSum = 0.0;
    double duneSum = 0.0;
    double ubooneSum = 0.0;
    double sbndSum = 0.0;

    for(int i = 1; i <= numBins; ++i){
        currentSum += current->GetBinContent(i);
        cheatedSum += cheated->GetBinContent(i);
        duneSum += dune->GetBinContent(i);
        ubooneSum += uboone->GetBinContent(i);
        sbndSum += sbnd->GetBinContent(i);

        double currentEffValue = currentSum/sizeCurrent;
        double cheatedEffValue = cheatedSum/sizeCheated;
        double duneEffValue = duneSum/sizeDune;
        double ubooneEffValue = ubooneSum/sizeUboone;
        double sbndEffValue = sbndSum/sizeSBND;

        currentEff->SetBinContent(i, currentEffValue);
        cheatedEff->SetBinContent(i, cheatedEffValue);
        duneEff->SetBinContent(i, duneEffValue);
        ubooneEff->SetBinContent(i, ubooneEffValue);
        sbndEff->SetBinContent(i, sbndEffValue);
    }

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of DL SBND Entries: %d", (int)sizeSBND));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 

    int funcValue = 1;

    styleDraw(efficiencyCanvas, currentEff, cheatedEff, duneEff, ubooneEff, sbndEff, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt = nullptr, &funcValue, drawLine, linePos);
}

recoSlice chooseSlice(std::vector<recoSlice> recoSlices, double method){
    int chosenSliceIndex = 0;
    double highestScore = -999999;
    if(method == 0){
        // Choose slice based on completeness
        for(size_t j = 0; j < recoSlices.size(); ++j){
            if(recoSlices[j].completeness > highestScore){
                highestScore = recoSlices[j].completeness;
                chosenSliceIndex = j;
            }
        } 
    } else if(method == 1){
        // Choose sliced based on CRUMBS score
        for(size_t j = 0; j < recoSlices.size(); ++j){
            if(recoSlices[j].score > highestScore){
                highestScore = recoSlices[j].score;
                chosenSliceIndex = j;
            }
        } 
    }
    
    return recoSlices[chosenSliceIndex];
}

recoParticle choosePFP(std::vector<recoParticle> recoParticles, double sliceID, double& totalEnergy, double& numPFPsSlice){
    totalEnergy = 0;
    numPFPsSlice = 0;
    int chosenParticleIndex = 0;
    double highestEnergy = -999999;

    for(size_t j = 0; j < recoParticles.size(); ++j){
        if(recoParticles[j].sliceID == sliceID){
            totalEnergy += recoParticles[j].bestPlaneEnergy;
            numPFPsSlice++;
            if(recoParticles[j].bestPlaneEnergy > highestEnergy){
                highestEnergy = recoParticles[j].bestPlaneEnergy;
                chosenParticleIndex = j;
            }
        }
    }

    return recoParticles[chosenParticleIndex];
}

recoNeutrino chooseRecoNeutrino(std::vector<recoNeutrino> recoNeutrinos, double sliceID){
    int chosenNeutrinoIndex;
    int counter = 0;

    //std::cout << "chosen slice ID: " << sliceID << std::endl;
    //std::cout << "reco Neutrinos size: " << recoNeutrinos.size() << std::endl;
    for(size_t j = 0; j < recoNeutrinos.size(); ++j){
        //std::cout << "reco neutrino " << j << ", slice ID: " << recoNeutrinos[j].sliceID << ", Vertex = (" << recoNeutrinos[j].vx << ", " << recoNeutrinos[j].vy << ", " << recoNeutrinos[j].vz << ")" << std::endl;
        if(recoNeutrinos[j].sliceID == sliceID){
            counter++;
            chosenNeutrinoIndex = j;
        }
    }

    if(counter > 1) printf("%i reco neutrinos in the chosen slice\n", counter);
    return recoNeutrinos[chosenNeutrinoIndex];
}

trueNeutrino chooseTrueNeutrino(std::vector<trueNeutrino> trueNeutrinos){
    // There should just be 1 true neutrino in the event
    //printf("%zu true neutrinos in the event and within the TPC\n", trueNeutrinos.size());
    //printf("True Neutrino Vertex: (%f, %f, %f)\n", trueNeutrinos[0].vx, trueNeutrinos[0].vy, trueNeutrinos[0].vz);
    return trueNeutrinos[0];
}

trueParticle chooseTrueParticle(std::vector<trueParticle> trueParticles, trueNeutrino chosenTrueNeutrino){
    int chosenParticleIndex = 0;
    
    for(size_t j = 0; j < trueParticles.size(); ++j){
        if(trueParticles[j].vx == chosenTrueNeutrino.vx && trueParticles[j].vy == chosenTrueNeutrino.vy && trueParticles[j].vz == chosenTrueNeutrino.vz){
            chosenParticleIndex = j;
        }
    }

    return trueParticles[chosenParticleIndex];
}

void nuE_macro_21May(){
    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_05_00NEW/NuEAnalyserOutput.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/merged.root");
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/merged.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/analysed_Current/NoRefinement/CRUMBS/1.root");
    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_05_00NEW/srcs/sbndcode/sbndcode/nue/NuEAnalyserOutput.root");
    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_05/srcs/sbndcode/sbndcode/nue/NuEAnalyserOutput.root");

    std::string base_path = "/nashome/c/coackley/nuEPlotsWithoutCosmics/";
    //std::string base_path = "/nashome/c/coackley/nuEPlotsWithoutCosmics/";

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
    std::vector<double> *truth_particleNeutrinoParent = nullptr;

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
    tree->SetBranchAddress("truth_particleNeutrinoParent", &truth_particleNeutrinoParent);

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

    // Truth Plots
    auto trueETheta2 = createHistGroup("trueETheta2", "E_{true}#theta_{true}^{2}", "E_{true}#theta_{true}^{2} (MeV rad^{2})", 32, 0, 4.088);

    // Number Plots
    auto numSlices = createHistGroup("numSlices", "Number of Slices in an Event", "Number of Slices", 45, 0, 45);
    // Number of slices with a CRUMBS score
    auto numSlicesCRUMBS = createHistGroup("numSlicesCRUMBS", "Number of Slices with CRUMBS Score in an Event", "Number of Slices with a CRUMBS Score", 10, 0, 10);
    // Number of slices with a completeness > 0
    auto numSlicesCompleteness = createHistGroup("numSlicesCompleteness", "Number of Slices with Completeness > 0 in an Event", "Number of Slices with Completeness > 0", 10, 0, 10);
    auto numRecoNeutrinos = createHistGroup("numRecoNeutrinos", "Number of Reco Neutrinos in an Event", "Number of Reco Neutrinos", 10, 0, 10);
    

    // Highest CRUMBS Score Slice Plots
    auto sliceCompletenessCRUMBS = createHistGroup("sliceCompletenessCRUMBS", "Completeness of the Slice with the Highest CRUMBS Score", "Completeness", 102, 0, 1.02); // Bin width = 0.005
    auto sliceScoreCRUMBS = createHistGroup("sliceScoreCRUMBS", "CRUMBS Score of the Slice with the Highest CRUMBS Score", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCRUMBS = createHistGroup("slicePurityCRUMBS", "Purity of the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto deltaXCRUMBS = createHistGroup("deltaXCRUMBS", "#Deltax Distribution: Slice with Highest CRUMBS Score", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaYCRUMBS = createHistGroup("deltaYCRUMBS", "#Deltay Distribution: Slice with Highest CRUMBS Score", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZCRUMBS = createHistGroup("deltaZCRUMBS", "#Deltaz Distribution: Slice with Highest CRUMBS Score", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaRCRUMBS = createHistGroup("deltaRCRUMBS", "#Delta#bar{r} Distribution: Slice with Highest CRUMBS Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    auto numPFPsCRUMBS = createHistGroup("numPFPsCRUMBS", "Number of PFPs in the Slice with the Highest CRUMBS Score", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCRUMBS = createHistGroup("ratioChosenSummedEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCRUMBS = createHistGroup("ratioChosenTrueEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCRUMBS = createHistGroup("ratioSummedTrueEnergyCRUMBS", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    auto angleDifferenceCRUMBS = createHistGroup("angleDifferenceCRUMBS", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto EtrueThetaRecoCRUMBS = createHistGroup("EtrueThetaRecoCRUMBS", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCRUMBS = createHistGroup("ERecoSumThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBS = createHistGroup("ERecoHighestThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCRUMBS = createHistGroup("ERecoSumThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCRUMBS = createHistGroup("ERecoHighestThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);

    // Highest Completeness Score Slice Plots
    auto sliceCompletenessCompleteness = createHistGroup("sliceCompletenessCompleteness", "Completeness of the Slice with the Highest Completeness", "Completeness Score", 102, 0, 1.02); // Bin width = 0.01
    auto sliceScoreCompleteness = createHistGroup("sliceScoreCompleteness", "CRUMBS Score of the Slice with the Highest Completeness", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCompleteness = createHistGroup("slicePurityCompleteness", "Purity of the Slice with the Highest Completeness", "Purity", 50, 0, 1.02);
    auto deltaXCompleteness = createHistGroup("deltaXCompleteness", "#Deltax Distribution: Slice with Highest Completeness", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaYCompleteness = createHistGroup("deltaYCompleteness", "#Deltay Distribution: Slice with Highest Completeness", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZCompleteness = createHistGroup("deltaZCompleteness", "#Deltaz Distribution: Slice with Highest Completeness", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaRCompleteness = createHistGroup("deltaRCompleteness", "#Delta#bar{r} Distribution: Slice with Highest Completeness Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    auto numPFPsCompleteness = createHistGroup("numPFPsCompleteness", "Number of PFPs in the Slice with the Highest Completeness", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCompleteness = createHistGroup("ratioChosenSummedEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest Completeness", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCompleteness = createHistGroup("ratioChosenTrueEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCompleteness = createHistGroup("ratioSummedTrueEnergyCompleteness", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    auto angleDifferenceCompleteness = createHistGroup("angleDifferenceCompleteness", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest Completeness", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto EtrueThetaRecoCompleteness = createHistGroup("EtrueThetaRecoCompleteness", "E_{true}#theta_{reco}^{2}: Slice with Highest Completeness", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCompleteness = createHistGroup("ERecoSumThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCompleteness = createHistGroup("ERecoHighestThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCompleteness = createHistGroup("ERecoSumThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCompleteness = createHistGroup("ERecoHighestThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    
    counter numEventsTotal;
    counter numEventsTrueNeutrino;
    counter numEventsTrueElectron;
    counter numEventsSlices;
    counter numEventsRecoNeutrino;
    counter numEventsSliceNotEqualNeutrino;
    counter numSlicesCRUMBSCompletenessZero;
    counter numEventsCRUMBSRecoParticle;
    counter sameSliceSelected;

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);

        printf("___________________________________________________________________\n");
        std::vector<recoParticle> recoParticlesInEvent;
        std::vector<recoNeutrino> recoNeutrinosInEvent;
        std::vector<recoSlice> recoSlicesInEvent;

        recoParticle chosenRecoParticleCRUMBS;
        recoParticle chosenRecoParticleCompleteness;
        recoNeutrino chosenRecoNeutrinoCRUMBS;
        recoNeutrino chosenRecoNeutrinoCompleteness;
        recoSlice chosenRecoSliceCompleteness;
        recoSlice chosenRecoSliceCRUMBS;
        trueParticle chosenTrueParticle;
        trueNeutrino chosenTrueNeutrino;
        double energyInSlice = 0;
        double numPFPsInEvent = 0;
        double numPFPsInChosenSlice = 0;

        if(DLCurrent == 0){
            numEventsTotal.uboone++;
        } else if(DLCurrent == 1){
            numEventsTotal.dune++;
        } else if(DLCurrent == 2){
            numEventsTotal.current++;
        } else if(DLCurrent == 3){
            numEventsTotal.cheated++;
        } else if(DLCurrent == 4){
            numEventsTotal.sbnd++;
        }

        int truthNeutrinoInTPC = 0;
        for(size_t j = 0; j < truth_neutrinoVX->size(); ++j){
            if(truth_neutrinoTPCValid->at(j) == 1){
                truthNeutrinoInTPC = 1;
                chosenTrueNeutrino.vx = truth_neutrinoVX->at(j);
                chosenTrueNeutrino.vy = truth_neutrinoVY->at(j);
                chosenTrueNeutrino.vz = truth_neutrinoVZ->at(j); 
                chosenTrueNeutrino.CCNC = truth_CCNC->at(j);
                chosenTrueNeutrino.pdg = truth_neutrinoType->at(j);
                chosenTrueNeutrino.leptonpdg = truth_chargedLepton->at(j);
                chosenTrueNeutrino.tpcID = truth_neutrinoTPCID->at(j);
                chosenTrueNeutrino.tpcValid = truth_neutrinoTPCValid->at(j);
            } else if(truth_neutrinoTPCValid->at(j) != -999999 && truth_neutrinoVX->at(j) < 201.3 && truth_neutrinoVX->at(j) > -201.3 && truth_neutrinoVY->at(j) < 203.732 && truth_neutrinoVY->at(j) > -203.732 && truth_neutrinoVZ->at(j) > 0 && truth_neutrinoVZ->at(j) < 509.4){
                //printf("True Neutrino Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            } else if(truth_neutrinoTPCValid->at(j) != -999999){
                //printf("True Neutrino Not Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            }
        }

        if(truth_neutrinoVX->size() == 1 && truth_neutrinoVX->at(0) == -999999) printf("No True Neutrino in Event\n");

        // If there are no truth neutrinos within the TPC, skip the event
        if(truthNeutrinoInTPC == 0) continue;

        if(DLCurrent == 0){
            numEventsTrueNeutrino.uboone++;
        } else if(DLCurrent == 1){
            numEventsTrueNeutrino.dune++;
        } else if(DLCurrent == 2){
            numEventsTrueNeutrino.current++;
        } else if(DLCurrent == 3){
            numEventsTrueNeutrino.cheated++;
        } else if(DLCurrent == 4){
            numEventsTrueNeutrino.sbnd++;
        }

        printf("True Neutrino Vertex = (%f, %f, %f)\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz);
        // Pick out the true recoil electron
        int truthElectron = 0;
        for(size_t j = 0; j < truth_particlePDG->size(); ++j){
            if(truth_particlePDG->at(j) == 11 && truth_particleNeutrinoParent->at(j) == 1){
                // This is an electron (PDG == 11) whose parent has TrackID of 0 (has neutrino parent == 1)
                printf("True Primary Electron: Vertex = (%f, %f, %f)\n", truth_particleVX->at(j), truth_particleVY->at(j), truth_particleVZ->at(j));
                if(std::abs(truth_particleVX->at(j) - chosenTrueNeutrino.vx) < 0.5 && std::abs(truth_particleVY->at(j) - chosenTrueNeutrino.vy) < 0.5 && std::abs(truth_particleVZ->at(j) - chosenTrueNeutrino.vz) < 0.5){
                    // This is the true recoil electron
                    truthElectron = 1;
                    chosenTrueParticle.pdg = truth_particlePDG->at(j);
                    chosenTrueParticle.vx = truth_particleVX->at(j);
                    chosenTrueParticle.vy = truth_particleVY->at(j);
                    chosenTrueParticle.vz = truth_particleVZ->at(j);
                    chosenTrueParticle.px = truth_particlePX->at(j);
                    chosenTrueParticle.py = truth_particlePY->at(j);
                    chosenTrueParticle.pz = truth_particlePZ->at(j);
                    chosenTrueParticle.energy = truth_particleEnergy->at(j);
                    chosenTrueParticle.angle = truth_particleAngle->at(j);
                    chosenTrueParticle.ETheta2 = truth_particleETheta2->at(j);
                    chosenTrueParticle.dx = truth_particleDirectionX->at(j);
                    chosenTrueParticle.dy = truth_particleDirectionY->at(j);
                    chosenTrueParticle.dz = truth_particleDirectionZ->at(j);
                    chosenTrueParticle.neutrinoParent = truth_particleNeutrinoParent->at(j);
                } else{
                    printf("Doesn't Pass!!\n");
                }      
            }
        }

        // If there are no true electrons from the primary neutrino, skip the event
        if(truthElectron == 0){
            printf("LOOK AT THIS ONE!!!!!!!!!\n");
            continue;
        }

        if(DLCurrent == 0){
            numEventsTrueElectron.uboone++;
            trueETheta2.uboone->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 1){
            numEventsTrueElectron.dune++;
            trueETheta2.dune->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 2){
            numEventsTrueElectron.current++;
            trueETheta2.current->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 3){
            numEventsTrueElectron.cheated++;
            trueETheta2.cheated->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 4){
            numEventsTrueElectron.sbnd++;
            trueETheta2.sbnd->Fill(chosenTrueParticle.ETheta2);
        }

        int slice = 0;
        int crumbsSlice = 0;
        int completenessSlice = 0;
        double numSlicesInEvent = 0;
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            if(reco_sliceID->at(j) != -999999){
                slice = 1;
                numSlicesInEvent++;

                recoSlice recoSlice;
                recoSlice.id = reco_sliceID->at(j);
                recoSlice.completeness = reco_sliceCompleteness->at(j);
                recoSlice.purity = reco_slicePurity->at(j);
                recoSlice.score = reco_sliceScore->at(j);
                recoSlicesInEvent.push_back(recoSlice);
        
                if(recoSlice.score != -999999) crumbsSlice++;  
                if(recoSlice.completeness > 0) completenessSlice++;  
                printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", j, recoSlice.id, recoSlice.completeness, recoSlice.purity, recoSlice.score);
            }
        }    

        // If there are no slices, skip the event
        if(slice == 0) continue;  
     
        // Pick the slice with the highest completeness
        chosenRecoSliceCompleteness = chooseSlice(recoSlicesInEvent, 0);
        printf("\nSlice with Highest Completeness: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", chosenRecoSliceCompleteness.id, chosenRecoSliceCompleteness.completeness, chosenRecoSliceCompleteness.purity, chosenRecoSliceCompleteness.score);

        // Pick the slice with the highest score
        chosenRecoSliceCRUMBS = chooseSlice(recoSlicesInEvent, 1);
        printf("Slice with Highest CRUMBS Score: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n\n", chosenRecoSliceCRUMBS.id, chosenRecoSliceCRUMBS.completeness, chosenRecoSliceCRUMBS.purity, chosenRecoSliceCRUMBS.score);
        
        if(DLCurrent == 0){
            numEventsSlices.uboone++;
            numSlices.uboone->Fill(numSlicesInEvent);
            numSlicesCRUMBS.uboone->Fill(crumbsSlice);
            numSlicesCompleteness.uboone->Fill(completenessSlice);
            
            sliceCompletenessCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.score);
            slicePurityCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.purity);
            
            sliceCompletenessCompleteness.uboone->Fill(chosenRecoSliceCompleteness.completeness);
            sliceScoreCompleteness.uboone->Fill(chosenRecoSliceCompleteness.score);
            slicePurityCompleteness.uboone->Fill(chosenRecoSliceCompleteness.purity);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.uboone++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.uboone++;
        } else if(DLCurrent == 1){
            numEventsSlices.dune++;
            numSlices.dune->Fill(numSlicesInEvent);
            numSlicesCRUMBS.dune->Fill(crumbsSlice);
            numSlicesCompleteness.dune->Fill(completenessSlice);

            sliceCompletenessCRUMBS.dune->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.dune->Fill(chosenRecoSliceCRUMBS.score);
            slicePurityCRUMBS.dune->Fill(chosenRecoSliceCRUMBS.purity);
            
            sliceCompletenessCompleteness.dune->Fill(chosenRecoSliceCompleteness.completeness);
            sliceScoreCompleteness.dune->Fill(chosenRecoSliceCompleteness.score);
            slicePurityCompleteness.dune->Fill(chosenRecoSliceCompleteness.purity);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.dune++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.dune++;
        } else if(DLCurrent == 2){
            numEventsSlices.current++;
            numSlices.current->Fill(numSlicesInEvent);
            numSlicesCRUMBS.current->Fill(crumbsSlice);
            numSlicesCompleteness.current->Fill(completenessSlice);
            
            sliceCompletenessCRUMBS.current->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.current->Fill(chosenRecoSliceCRUMBS.score);
            slicePurityCRUMBS.current->Fill(chosenRecoSliceCRUMBS.purity);
            
            sliceCompletenessCompleteness.current->Fill(chosenRecoSliceCompleteness.completeness);
            sliceScoreCompleteness.current->Fill(chosenRecoSliceCompleteness.score);
            slicePurityCompleteness.current->Fill(chosenRecoSliceCompleteness.purity);
                
            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.current++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.current++;
        } else if(DLCurrent == 3){
            numEventsSlices.cheated++;
            numSlices.cheated->Fill(numSlicesInEvent);
            numSlicesCRUMBS.cheated->Fill(crumbsSlice);
            numSlicesCompleteness.cheated->Fill(completenessSlice);
            
            sliceCompletenessCRUMBS.cheated->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.cheated->Fill(chosenRecoSliceCRUMBS.score);
            slicePurityCRUMBS.cheated->Fill(chosenRecoSliceCRUMBS.purity);
            
            sliceCompletenessCompleteness.cheated->Fill(chosenRecoSliceCompleteness.completeness);
            sliceScoreCompleteness.cheated->Fill(chosenRecoSliceCompleteness.score);
            slicePurityCompleteness.cheated->Fill(chosenRecoSliceCompleteness.purity);
            
            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.cheated++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.cheated++;
        } else if(DLCurrent == 4){
            numEventsSlices.sbnd++;
            numSlices.sbnd->Fill(numSlicesInEvent);
            numSlicesCRUMBS.sbnd->Fill(crumbsSlice);
            numSlicesCompleteness.sbnd->Fill(completenessSlice);

            sliceCompletenessCRUMBS.sbnd->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.sbnd->Fill(chosenRecoSliceCRUMBS.score);
            slicePurityCRUMBS.sbnd->Fill(chosenRecoSliceCRUMBS.purity);
            
            sliceCompletenessCompleteness.sbnd->Fill(chosenRecoSliceCompleteness.completeness);
            sliceScoreCompleteness.sbnd->Fill(chosenRecoSliceCompleteness.score);
            slicePurityCompleteness.sbnd->Fill(chosenRecoSliceCompleteness.purity);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.sbnd++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.sbnd++;
        }
        
        int neutrino = 0;
        int numRecoNeutrinosInEvent = 0;
        std::cout << "num reco neutrinos: " << reco_neutrinoPDG->size() << std::endl;
        for(size_t j = 0; j < reco_neutrinoPDG->size(); ++j){
            if(reco_neutrinoPDG->at(j) != -999999){
                neutrino = 1;
                numRecoNeutrinosInEvent++;
                std::cout << "1" << std::endl;
                
                recoNeutrino recoNeutrino;
                recoNeutrino.pdg = reco_neutrinoPDG->at(j);
                recoNeutrino.isPrimary = reco_neutrinoIsPrimary->at(j);
                recoNeutrino.vx = reco_neutrinoVX->at(j);
                recoNeutrino.vy = reco_neutrinoVY->at(j);
                recoNeutrino.vz = reco_neutrinoVZ->at(j);
                recoNeutrino.sliceID = reco_neutrinoSliceID->at(j);
                recoNeutrinosInEvent.push_back(recoNeutrino);

                printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, recoNeutrino.pdg, recoNeutrino.isPrimary, recoNeutrino.vx, recoNeutrino.vy, recoNeutrino.vz, recoNeutrino.sliceID);
            }
        }
        
        if(DLCurrent == 0) numRecoNeutrinos.uboone->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 1) numRecoNeutrinos.dune->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 2) numRecoNeutrinos.current->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 3) numRecoNeutrinos.cheated->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 4) numRecoNeutrinos.sbnd->Fill(numRecoNeutrinosInEvent);

        // Skip the event if there are no reconstructed neutrinos
        if(neutrino == 0) continue;
    
        // The number of slices with a CRUMBS score != the number of reconstructed neutrinos in the event 
        if(crumbsSlice != numRecoNeutrinosInEvent){
            if(DLCurrent == 0) numEventsSliceNotEqualNeutrino.uboone++; 
            if(DLCurrent == 1) numEventsSliceNotEqualNeutrino.dune++;     
            if(DLCurrent == 2) numEventsSliceNotEqualNeutrino.current++;
            if(DLCurrent == 3) numEventsSliceNotEqualNeutrino.cheated++;
            if(DLCurrent == 4) numEventsSliceNotEqualNeutrino.sbnd++;
        } 

        // Look at reconstruction when we pick the slice with highest CRUMBS score
        chosenRecoNeutrinoCRUMBS = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCRUMBS.id); 

        double deltaXCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vx - chosenTrueNeutrino.vx);
        double deltaYCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vy - chosenTrueNeutrino.vy);
        double deltaZCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vz - chosenTrueNeutrino.vz);
        double deltaRCRUMBSValue = std::sqrt((deltaXCRUMBSValue * deltaXCRUMBSValue) + (deltaYCRUMBSValue * deltaYCRUMBSValue) + (deltaZCRUMBSValue * deltaZCRUMBSValue));

        chosenRecoNeutrinoCompleteness = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCompleteness.id);
        double deltaXCompletenessValue = (chosenRecoNeutrinoCompleteness.vx - chosenTrueNeutrino.vx);
        double deltaYCompletenessValue = (chosenRecoNeutrinoCompleteness.vy - chosenTrueNeutrino.vy);
        double deltaZCompletenessValue = (chosenRecoNeutrinoCompleteness.vz - chosenTrueNeutrino.vz);
        double deltaRCompletenessValue = std::sqrt((deltaXCompletenessValue * deltaXCompletenessValue) + (deltaYCompletenessValue * deltaYCompletenessValue) + (deltaZCompletenessValue * deltaZCompletenessValue));

        if(DLCurrent == 0){
            numEventsRecoNeutrino.uboone++;
            deltaXCRUMBS.uboone->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.uboone->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.uboone->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.uboone->Fill(deltaRCRUMBSValue);

            deltaXCompleteness.uboone->Fill(deltaXCompletenessValue);
            deltaYCompleteness.uboone->Fill(deltaYCompletenessValue);
            deltaZCompleteness.uboone->Fill(deltaZCompletenessValue);
            deltaRCompleteness.uboone->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 1){
            numEventsRecoNeutrino.dune++;
            deltaXCRUMBS.dune->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.dune->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.dune->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.dune->Fill(deltaRCRUMBSValue);
        
            deltaXCompleteness.dune->Fill(deltaXCompletenessValue);
            deltaYCompleteness.dune->Fill(deltaYCompletenessValue);
            deltaZCompleteness.dune->Fill(deltaZCompletenessValue);
            deltaRCompleteness.dune->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 2){
            numEventsRecoNeutrino.current++;
            deltaXCRUMBS.current->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.current->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.current->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.current->Fill(deltaRCRUMBSValue);
        
            deltaXCompleteness.current->Fill(deltaXCompletenessValue);
            deltaYCompleteness.current->Fill(deltaYCompletenessValue);
            deltaZCompleteness.current->Fill(deltaZCompletenessValue);
            deltaRCompleteness.current->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 3){
            numEventsRecoNeutrino.cheated++;
            deltaXCRUMBS.cheated->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.cheated->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.cheated->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.cheated->Fill(deltaRCRUMBSValue);
        
            deltaXCompleteness.cheated->Fill(deltaXCompletenessValue);
            deltaYCompleteness.cheated->Fill(deltaYCompletenessValue);
            deltaZCompleteness.cheated->Fill(deltaZCompletenessValue);
            deltaRCompleteness.cheated->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 4){
            numEventsRecoNeutrino.sbnd++;
            deltaXCRUMBS.sbnd->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.sbnd->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.sbnd->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.sbnd->Fill(deltaRCRUMBSValue);
        
            deltaXCompleteness.sbnd->Fill(deltaXCompletenessValue);
            deltaYCompleteness.sbnd->Fill(deltaYCompletenessValue);
            deltaZCompleteness.sbnd->Fill(deltaZCompletenessValue);
            deltaRCompleteness.sbnd->Fill(deltaRCompletenessValue);
        }

        int recoparticle = 0;
        int recoparticleCRUMBS = 0;
        int recoparticleCompleteness = 0;

        for(size_t j = 0; j < reco_particlePDG->size(); ++j){
            if(reco_particlePDG->at(j) != -999999){
                recoparticle = 1;

                recoParticle recoParticle;
                recoParticle.pdg = reco_particlePDG->at(j);
                recoParticle.isPrimary = reco_particleIsPrimary->at(j);
                recoParticle.vx = reco_particleVX->at(j);
                recoParticle.vy = reco_particleVY->at(j);
                recoParticle.vz = reco_particleVZ->at(j);
                recoParticle.dx = reco_particleDirectionX->at(j);
                recoParticle.dy = reco_particleDirectionY->at(j);
                recoParticle.dz = reco_particleDirectionZ->at(j);
                recoParticle.sliceID = reco_particleSliceID->at(j);
                recoParticle.bestPlaneEnergy = reco_particleBestPlaneEnergy->at(j);
                recoParticle.theta = reco_particleTheta->at(j);
                recoParticle.trackscore = reco_particleTrackScore->at(j);
                recoParticle.completeness = reco_particleCompleteness->at(j);
                recoParticle.purity = reco_particlePurity->at(j);
                recoParticlesInEvent.push_back(recoParticle);
            
                if(recoParticle.sliceID == chosenRecoSliceCRUMBS.id) recoparticleCRUMBS = 1;
                //if(recoParticle.sliceID == chosenRecoSliceCompleteness.id) recoparticleCompleteness = 1;
            }
        }

        // There are no reco particles in the event
        if(recoparticle == 0) continue;
        
        // No reco particles in the slice with the highest CRUMBS score
        if(recoparticleCRUMBS == 0) continue;
        
        // No reco particles in the slice with the highest completeness
        //if(recoparticleCompleteness == 0) continue;

        double totalSliceEnergyCRUMBS = 0;
        double numPFPsSliceCRUMBS = 0;
        chosenRecoParticleCRUMBS = choosePFP(recoParticlesInEvent, chosenRecoSliceCRUMBS.id, totalSliceEnergyCRUMBS, numPFPsSliceCRUMBS);

        double aDOTbCRUMBS = ((chosenRecoParticleCRUMBS.dx * chosenTrueParticle.dx) + (chosenRecoParticleCRUMBS.dy * chosenTrueParticle.dy) + (chosenRecoParticleCRUMBS.dz * chosenTrueParticle.dz));
        double aMagCRUMBS = std::sqrt((chosenRecoParticleCRUMBS.dx * chosenRecoParticleCRUMBS.dx) + (chosenRecoParticleCRUMBS.dy * chosenRecoParticleCRUMBS.dy) + (chosenRecoParticleCRUMBS.dz * chosenRecoParticleCRUMBS.dz));
        double bMagCRUMBS = std::sqrt((chosenTrueParticle.dx * chosenTrueParticle.dx) + (chosenTrueParticle.dy * chosenTrueParticle.dy) + (chosenTrueParticle.dz * chosenTrueParticle.dz));
        double cosAngleCRUMBS = aDOTbCRUMBS/(aMagCRUMBS * bMagCRUMBS);
        if(cosAngleCRUMBS > 1.0) cosAngleCRUMBS = 1.0;
        if(cosAngleCRUMBS < -1.0) cosAngleCRUMBS = -1.0;
        double angleDiffCRUMBS = TMath::ACos(cosAngleCRUMBS) * TMath::RadToDeg();

        double totalSliceEnergyCompleteness = 0;
        double numPFPsSliceCompleteness = 0;
        chosenRecoParticleCompleteness = choosePFP(recoParticlesInEvent, chosenRecoSliceCompleteness.id, totalSliceEnergyCompleteness, numPFPsSliceCompleteness);
        
        double aDOTbCompleteness = ((chosenRecoParticleCompleteness.dx * chosenTrueParticle.dx) + (chosenRecoParticleCompleteness.dy * chosenTrueParticle.dy) + (chosenRecoParticleCompleteness.dz * chosenTrueParticle.dz));
        double aMagCompleteness = std::sqrt((chosenRecoParticleCompleteness.dx * chosenRecoParticleCompleteness.dx) + (chosenRecoParticleCompleteness.dy * chosenRecoParticleCompleteness.dy) + (chosenRecoParticleCompleteness.dz * chosenRecoParticleCompleteness.dz));
        double bMagCompleteness = std::sqrt((chosenTrueParticle.dx * chosenTrueParticle.dx) + (chosenTrueParticle.dy * chosenTrueParticle.dy) + (chosenTrueParticle.dz * chosenTrueParticle.dz));
        double cosAngleCompleteness = aDOTbCompleteness/(aMagCompleteness * bMagCompleteness);
        if(cosAngleCompleteness > 1.0) cosAngleCompleteness = 1.0;
        if(cosAngleCompleteness < -1.0) cosAngleCompleteness = -1.0;
        double angleDiffCompleteness = TMath::ACos(cosAngleCompleteness) * TMath::RadToDeg();
        
        if(DLCurrent == 0){
            numEventsCRUMBSRecoParticle.uboone++;
            numPFPsCRUMBS.uboone->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.uboone->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.uboone->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);

            numPFPsCompleteness.uboone->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.uboone->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.uboone->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.uboone->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.uboone->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.uboone->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 1){
            numEventsCRUMBSRecoParticle.dune++;
            numPFPsCRUMBS.dune->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.dune->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.dune->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.dune->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.dune->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.dune->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        
            numPFPsCompleteness.dune->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.dune->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.dune->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.dune->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.dune->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.dune->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 2){
            numEventsCRUMBSRecoParticle.current++;
            numPFPsCRUMBS.current->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.current->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.current->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.current->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.current->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.current->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
            numPFPsCompleteness.current->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.current->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.current->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.current->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.current->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.current->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 3){
            numEventsCRUMBSRecoParticle.cheated++;
            numPFPsCRUMBS.cheated->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.cheated->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.cheated->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.cheated->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.cheated->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.cheated->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        
            numPFPsCompleteness.cheated->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.cheated->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.cheated->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.cheated->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.cheated->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.cheated->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 4){
            numEventsCRUMBSRecoParticle.sbnd++;
            numPFPsCRUMBS.sbnd->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.sbnd->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.sbnd->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.sbnd->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.sbnd->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.sbnd->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        
            numPFPsCompleteness.sbnd->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.sbnd->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.sbnd->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.sbnd->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.sbnd->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.sbnd->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        }
        
        printf("Chosen True Neutrino: Vertex = (%f, %f, %f), CCNC = %f, PDG = %f, Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz, chosenTrueNeutrino.CCNC, chosenTrueNeutrino.pdg, chosenTrueNeutrino.leptonpdg, chosenTrueNeutrino.tpcID, chosenTrueNeutrino.tpcValid);    
        printf("Chosen True Recoil Electron: Vertex = (%f, %f, %f), PDG = %f, Momentum = (%f, %f, %f), Energy = %f, Angle = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", chosenTrueParticle.vx, chosenTrueParticle.vy, chosenTrueParticle.vz, chosenTrueParticle.pdg, chosenTrueParticle.px, chosenTrueParticle.py, chosenTrueParticle.pz, chosenTrueParticle.energy, chosenTrueParticle.angle, chosenTrueParticle.ETheta2, chosenTrueParticle.dx, chosenTrueParticle.dy, chosenTrueParticle.dz);
        printf("Number of Slices in event: %f\n", numSlicesInEvent);
    }

    printf("Number of events with a true neutrino within the TPC:\nUBoone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsTrueNeutrino.uboone, numEventsTotal.uboone, numEventsTrueNeutrino.dune, numEventsTotal.dune, numEventsTrueNeutrino.current, numEventsTotal.current, numEventsTrueNeutrino.cheated, numEventsTotal.cheated, numEventsTrueNeutrino.sbnd, numEventsTotal.sbnd);
    printf("Number of events with a true recoil electron:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsTrueElectron.uboone, numEventsTotal.uboone, numEventsTrueElectron.dune, numEventsTotal.dune, numEventsTrueElectron.current, numEventsTotal.current, numEventsTrueElectron.cheated, numEventsTotal.cheated, numEventsTrueElectron.sbnd, numEventsTotal.sbnd);
    printf("Number of events with a slice:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsSlices.uboone, numEventsTotal.uboone, numEventsSlices.dune, numEventsTotal.dune, numEventsSlices.current, numEventsTotal.current, numEventsSlices.cheated, numEventsTotal.cheated, numEventsSlices.sbnd, numEventsTotal.sbnd);
    printf("Number of events with a reco neutrino:\nUboone:%i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsRecoNeutrino.uboone, numEventsTotal.uboone, numEventsRecoNeutrino.dune, numEventsTotal.dune, numEventsRecoNeutrino.current, numEventsTotal.current, numEventsRecoNeutrino.cheated, numEventsTotal.cheated, numEventsRecoNeutrino.sbnd, numEventsTotal.sbnd);
    printf("Number of events where the number of slices with a CRUMBS score != number of reco neutrinos:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsSliceNotEqualNeutrino.uboone, numEventsRecoNeutrino.uboone, numEventsSliceNotEqualNeutrino.dune, numEventsRecoNeutrino.dune, numEventsSliceNotEqualNeutrino.current, numEventsRecoNeutrino.current, numEventsSliceNotEqualNeutrino.cheated, numEventsRecoNeutrino.cheated, numEventsSliceNotEqualNeutrino.sbnd, numEventsRecoNeutrino.sbnd);
    printf("Number of events where the chosen Slice has the highest CRUMBS score and completeness:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", sameSliceSelected.uboone, numEventsSlices.uboone, sameSliceSelected.dune, numEventsSlices.dune, sameSliceSelected.current, numEventsSlices.current, sameSliceSelected.cheated, numEventsSlices.cheated, sameSliceSelected.sbnd, numEventsSlices.sbnd);
    printf("Number of events where the slice with the highest CRUMBS score has a completeness of 0:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numSlicesCRUMBSCompletenessZero.uboone, numEventsSlices.uboone, numSlicesCRUMBSCompletenessZero.dune, numEventsSlices.dune, numSlicesCRUMBSCompletenessZero.current, numEventsSlices.current, numSlicesCRUMBSCompletenessZero.cheated, numEventsSlices.cheated, numSlicesCRUMBSCompletenessZero.sbnd, numEventsSlices.sbnd);

    styleDraw(numSlices.canvas, numSlices.current, numSlices.cheated, numSlices.dune, numSlices.uboone, numSlices.sbnd, 0, 900, 999, 999, (base_path + "numSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlices.current, numSlices.cheated, numSlices.dune, numSlices.uboone, numSlices.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 100, 999, 999, (base_path + "numSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
 
    styleDraw(numSlicesCRUMBS.canvas, numSlicesCRUMBS.current, numSlicesCRUMBS.cheated, numSlicesCRUMBS.dune, numSlicesCRUMBS.uboone, numSlicesCRUMBS.sbnd, 0, 900, 999, 999, (base_path + "numCRUMBSSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlicesCRUMBS.current, numSlicesCRUMBS.cheated, numSlicesCRUMBS.dune, numSlicesCRUMBS.uboone, numSlicesCRUMBS.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 100, 999, 999, (base_path + "numCRUMBSSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(numSlicesCompleteness.canvas, numSlicesCompleteness.current, numSlicesCompleteness.cheated, numSlicesCompleteness.dune, numSlicesCompleteness.uboone, numSlicesCompleteness.sbnd, 0, 900, 999, 999, (base_path + "numCompletenessSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlicesCompleteness.current, numSlicesCompleteness.cheated, numSlicesCompleteness.dune, numSlicesCompleteness.uboone, numSlicesCompleteness.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 100, 999, 999, (base_path + "numCompletenessSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    
    styleDraw(numRecoNeutrinos.canvas, numRecoNeutrinos.current, numRecoNeutrinos.cheated, numRecoNeutrinos.dune, numRecoNeutrinos.uboone, numRecoNeutrinos.sbnd, 0, 900, 999, 999, (base_path + "numRecoNeutrinos_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numRecoNeutrinos.current, numRecoNeutrinos.cheated, numRecoNeutrinos.dune, numRecoNeutrinos.uboone, numRecoNeutrinos.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 100, 999, 999, (base_path + "numRecoNeutrinos_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // CRUMBS Slice Score Plots
    styleDraw(sliceCompletenessCRUMBS.canvas, sliceCompletenessCRUMBS.current, sliceCompletenessCRUMBS.cheated, sliceCompletenessCRUMBS.dune, sliceCompletenessCRUMBS.uboone, sliceCompletenessCRUMBS.sbnd, 0, 900, 999, 999, (base_path + "sliceCompletenessCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceCompletenessCRUMBS.current, sliceCompletenessCRUMBS.cheated, sliceCompletenessCRUMBS.dune, sliceCompletenessCRUMBS.uboone, sliceCompletenessCRUMBS.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 100, 999, 999, (base_path + "sliceCompletenessCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceScoreCRUMBS.canvas, sliceScoreCRUMBS.current, sliceScoreCRUMBS.cheated, sliceScoreCRUMBS.dune, sliceScoreCRUMBS.uboone, sliceScoreCRUMBS.sbnd, 0, 160, 999, 999, (base_path + "sliceScoreCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceScoreCRUMBS.current, sliceScoreCRUMBS.cheated, sliceScoreCRUMBS.dune, sliceScoreCRUMBS.uboone, sliceScoreCRUMBS.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 18, 999, 999, (base_path + "sliceScoreCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(slicePurityCRUMBS.canvas, slicePurityCRUMBS.current, slicePurityCRUMBS.cheated, slicePurityCRUMBS.dune, slicePurityCRUMBS.uboone, slicePurityCRUMBS.sbnd, 0, 60, 999, 999, (base_path + "slicePurityCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(slicePurityCRUMBS.current, slicePurityCRUMBS.cheated, slicePurityCRUMBS.dune, slicePurityCRUMBS.uboone, slicePurityCRUMBS.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 10, 999, 999, (base_path + "slicePurityCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    
    // CRUMBS Slice Vertex Plots
    styleDraw(deltaXCRUMBS.canvas, deltaXCRUMBS.current, deltaXCRUMBS.cheated, deltaXCRUMBS.dune, deltaXCRUMBS.uboone, deltaXCRUMBS.sbnd, 0, 460, 999, 999, (base_path + "deltaXCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaYCRUMBS.canvas, deltaYCRUMBS.current, deltaYCRUMBS.cheated, deltaYCRUMBS.dune, deltaYCRUMBS.uboone, deltaYCRUMBS.sbnd, 0, 460, 999, 999, (base_path + "deltaYCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaZCRUMBS.canvas, deltaZCRUMBS.current, deltaZCRUMBS.cheated, deltaZCRUMBS.dune, deltaZCRUMBS.uboone, deltaZCRUMBS.sbnd, 0, 460, 999, 999, (base_path + "deltaZCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaRCRUMBS.canvas, deltaRCRUMBS.current, deltaRCRUMBS.cheated, deltaRCRUMBS.dune, deltaRCRUMBS.uboone, deltaRCRUMBS.sbnd, 0, 860, 999, 999, (base_path + "deltaRCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaXCRUMBS.current, deltaXCRUMBS.cheated, deltaXCRUMBS.dune, deltaXCRUMBS.uboone, deltaXCRUMBS.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 52, 999, 999, (base_path + "deltaXCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaYCRUMBS.current, deltaYCRUMBS.cheated, deltaYCRUMBS.dune, deltaYCRUMBS.uboone, deltaYCRUMBS.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 50, 999, 999, (base_path + "deltaYCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaZCRUMBS.current, deltaZCRUMBS.cheated, deltaZCRUMBS.dune, deltaZCRUMBS.uboone, deltaZCRUMBS.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 50, 999, 999, (base_path + "deltaZCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaRCRUMBS.current, deltaRCRUMBS.cheated, deltaRCRUMBS.dune, deltaRCRUMBS.uboone, deltaRCRUMBS.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 100, 999, 999, (base_path + "deltaRCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // CRUMBS PFP Plots
    styleDraw(numPFPsCRUMBS.canvas, numPFPsCRUMBS.current, numPFPsCRUMBS.cheated, numPFPsCRUMBS.dune, numPFPsCRUMBS.uboone, numPFPsCRUMBS.sbnd, 0, 820, 999, 999, (base_path + "numPFPsCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numPFPsCRUMBS.current, numPFPsCRUMBS.cheated, numPFPsCRUMBS.dune, numPFPsCRUMBS.uboone, numPFPsCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 100, 999, 999, (base_path + "numPFPsCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(ratioChosenSummedEnergyCRUMBS.canvas, ratioChosenSummedEnergyCRUMBS.current, ratioChosenSummedEnergyCRUMBS.cheated, ratioChosenSummedEnergyCRUMBS.dune, ratioChosenSummedEnergyCRUMBS.uboone, ratioChosenSummedEnergyCRUMBS.sbnd, 0, 820, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenSummedEnergyCRUMBS.current, ratioChosenSummedEnergyCRUMBS.cheated, ratioChosenSummedEnergyCRUMBS.dune, ratioChosenSummedEnergyCRUMBS.uboone, ratioChosenSummedEnergyCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 92, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioChosenTrueEnergyCRUMBS.canvas, ratioChosenTrueEnergyCRUMBS.current, ratioChosenTrueEnergyCRUMBS.cheated, ratioChosenTrueEnergyCRUMBS.dune, ratioChosenTrueEnergyCRUMBS.uboone, ratioChosenTrueEnergyCRUMBS.sbnd, 0, 300, 999, 999, (base_path + "ratioChosenTrueEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenTrueEnergyCRUMBS.current, ratioChosenTrueEnergyCRUMBS.cheated, ratioChosenTrueEnergyCRUMBS.dune, ratioChosenTrueEnergyCRUMBS.uboone, ratioChosenTrueEnergyCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 34, 999, 999, (base_path + "ratioChosenTrueEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioSummedTrueEnergyCRUMBS.canvas, ratioSummedTrueEnergyCRUMBS.current, ratioSummedTrueEnergyCRUMBS.cheated, ratioSummedTrueEnergyCRUMBS.dune, ratioSummedTrueEnergyCRUMBS.uboone, ratioSummedTrueEnergyCRUMBS.sbnd, 0, 320, 999, 999, (base_path + "ratioSummedTrueEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioSummedTrueEnergyCRUMBS.current, ratioSummedTrueEnergyCRUMBS.cheated, ratioSummedTrueEnergyCRUMBS.dune, ratioSummedTrueEnergyCRUMBS.uboone, ratioSummedTrueEnergyCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 36, 999, 999, (base_path + "ratioSummedTrueEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    int drawLine = 1;
    int left = 0;
    int right = 1;
    
    styleDraw(angleDifferenceCRUMBS.canvas, angleDifferenceCRUMBS.current, angleDifferenceCRUMBS.cheated, angleDifferenceCRUMBS.dune, angleDifferenceCRUMBS.uboone, angleDifferenceCRUMBS.sbnd, 0, 260, 999, 999, (base_path + "angleDifferenceCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(angleDifferenceCRUMBS.current, angleDifferenceCRUMBS.cheated, angleDifferenceCRUMBS.dune, angleDifferenceCRUMBS.uboone, angleDifferenceCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 38, 999, 999, (base_path + "angleDifferenceCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
   
    styleDraw(EtrueThetaRecoCRUMBS.canvas, EtrueThetaRecoCRUMBS.current, EtrueThetaRecoCRUMBS.cheated, EtrueThetaRecoCRUMBS.dune, EtrueThetaRecoCRUMBS.sbnd, EtrueThetaRecoCRUMBS.uboone, 0, 120, 999, 999, (base_path + "EtrueThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(EtrueThetaRecoCRUMBS.current, EtrueThetaRecoCRUMBS.cheated, EtrueThetaRecoCRUMBS.dune, EtrueThetaRecoCRUMBS.uboone, EtrueThetaRecoCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 14, 999, 999, (base_path + "EtrueThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(EtrueThetaRecoCRUMBS.current, EtrueThetaRecoCRUMBS.cheated, EtrueThetaRecoCRUMBS.dune, EtrueThetaRecoCRUMBS.uboone, EtrueThetaRecoCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "EtrueThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{true}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoSumThetaTrueCRUMBS.canvas, ERecoSumThetaTrueCRUMBS.current, ERecoSumThetaTrueCRUMBS.cheated, ERecoSumThetaTrueCRUMBS.dune, ERecoSumThetaTrueCRUMBS.uboone, ERecoSumThetaTrueCRUMBS.sbnd, 0, 160, 999, 999, (base_path + "ERecoSumThetaTrueCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaTrueCRUMBS.current, ERecoSumThetaTrueCRUMBS.cheated, ERecoSumThetaTrueCRUMBS.dune, ERecoSumThetaTrueCRUMBS.uboone, ERecoSumThetaTrueCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 18, 999, 999, (base_path + "ERecoSumThetaTrueCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaTrueCRUMBS.current, ERecoSumThetaTrueCRUMBS.cheated, ERecoSumThetaTrueCRUMBS.dune, ERecoSumThetaTrueCRUMBS.uboone, ERecoSumThetaTrueCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoSumThetaTrueCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoHighestThetaTrueCRUMBS.canvas, ERecoHighestThetaTrueCRUMBS.current, ERecoHighestThetaTrueCRUMBS.cheated, ERecoHighestThetaTrueCRUMBS.dune, ERecoHighestThetaTrueCRUMBS.uboone, ERecoHighestThetaTrueCRUMBS.sbnd, 0, 180, 999, 999, (base_path + "ERecoHighestThetaTrueCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaTrueCRUMBS.current, ERecoHighestThetaTrueCRUMBS.cheated, ERecoHighestThetaTrueCRUMBS.dune, ERecoHighestThetaTrueCRUMBS.uboone, ERecoHighestThetaTrueCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 20, 999, 999, (base_path + "ERecoHighestThetaTrueCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaTrueCRUMBS.current, ERecoHighestThetaTrueCRUMBS.cheated, ERecoHighestThetaTrueCRUMBS.dune, ERecoHighestThetaTrueCRUMBS.uboone, ERecoHighestThetaTrueCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoHighestThetaTrueCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoSumThetaRecoCRUMBS.canvas, ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.cheated, ERecoSumThetaRecoCRUMBS.dune, ERecoSumThetaRecoCRUMBS.uboone, ERecoSumThetaRecoCRUMBS.sbnd, 0, 180, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.cheated, ERecoSumThetaRecoCRUMBS.dune, ERecoSumThetaRecoCRUMBS.uboone, ERecoSumThetaRecoCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 20, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.cheated, ERecoSumThetaRecoCRUMBS.dune, ERecoSumThetaRecoCRUMBS.uboone, ERecoSumThetaRecoCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoHighestThetaRecoCRUMBS.canvas, ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.cheated, ERecoHighestThetaRecoCRUMBS.dune, ERecoHighestThetaRecoCRUMBS.uboone, ERecoHighestThetaRecoCRUMBS.sbnd, 0, 180, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.cheated, ERecoHighestThetaRecoCRUMBS.dune, ERecoHighestThetaRecoCRUMBS.uboone, ERecoHighestThetaRecoCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 20, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.cheated, ERecoHighestThetaRecoCRUMBS.dune, ERecoHighestThetaRecoCRUMBS.uboone, ERecoHighestThetaRecoCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    // Completeness Slice Score Plots
    styleDraw(sliceCompletenessCompleteness.canvas, sliceCompletenessCompleteness.current, sliceCompletenessCompleteness.cheated, sliceCompletenessCompleteness.dune, sliceCompletenessCompleteness.uboone, sliceCompletenessCompleteness.sbnd, 0, 900, 999, 999, (base_path + "sliceCompletenessCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceCompletenessCompleteness.current, sliceCompletenessCompleteness.cheated, sliceCompletenessCompleteness.dune, sliceCompletenessCompleteness.uboone, sliceCompletenessCompleteness.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 100, 999, 999, (base_path + "sliceCompletenessCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceScoreCompleteness.canvas, sliceScoreCompleteness.current, sliceScoreCompleteness.cheated, sliceScoreCompleteness.dune, sliceScoreCompleteness.uboone, sliceScoreCompleteness.sbnd, 0, 160, 999, 999, (base_path + "sliceScoreCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceScoreCompleteness.current, sliceScoreCompleteness.cheated, sliceScoreCompleteness.dune, sliceScoreCompleteness.uboone, sliceScoreCompleteness.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 16, 999, 999, (base_path + "sliceScoreCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(slicePurityCompleteness.canvas, slicePurityCompleteness.current, slicePurityCompleteness.cheated, slicePurityCompleteness.dune, slicePurityCompleteness.uboone, slicePurityCompleteness.sbnd, 0, 60, 999, 999, (base_path + "slicePurityCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(slicePurityCompleteness.current, slicePurityCompleteness.cheated, slicePurityCompleteness.dune, slicePurityCompleteness.uboone, slicePurityCompleteness.sbnd, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, 0, 6, 999, 999, (base_path + "slicePurityCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    styleDraw(deltaXCompleteness.canvas, deltaXCompleteness.current, deltaXCompleteness.cheated, deltaXCompleteness.dune, deltaXCompleteness.uboone, deltaXCompleteness.sbnd, 0, 460, 999, 999, (base_path + "deltaXCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaYCompleteness.canvas, deltaYCompleteness.current, deltaYCompleteness.cheated, deltaYCompleteness.dune, deltaYCompleteness.uboone, deltaYCompleteness.sbnd, 0, 460, 999, 999, (base_path + "deltaYCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaZCompleteness.canvas, deltaZCompleteness.current, deltaZCompleteness.cheated, deltaZCompleteness.dune, deltaZCompleteness.uboone, deltaZCompleteness.sbnd, 0, 460, 999, 999, (base_path + "deltaZCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaRCompleteness.canvas, deltaRCompleteness.current, deltaRCompleteness.cheated, deltaRCompleteness.dune, deltaRCompleteness.uboone, deltaRCompleteness.sbnd, 0, 860, 999, 999, (base_path + "deltaRCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaXCompleteness.current, deltaXCompleteness.cheated, deltaXCompleteness.dune, deltaXCompleteness.uboone, deltaXCompleteness.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 52, 999, 999, (base_path + "deltaXCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaYCompleteness.current, deltaYCompleteness.cheated, deltaYCompleteness.dune, deltaYCompleteness.uboone, deltaYCompleteness.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 52, 999, 999, (base_path + "deltaYCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaZCompleteness.current, deltaZCompleteness.cheated, deltaZCompleteness.dune, deltaZCompleteness.uboone, deltaZCompleteness.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 52, 999, 999, (base_path + "deltaZCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaRCompleteness.current, deltaRCompleteness.cheated, deltaRCompleteness.dune, deltaRCompleteness.uboone, deltaRCompleteness.sbnd, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, 0, 100, 999, 999, (base_path + "deltaRCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // Completeness PFP Plots
    styleDraw(numPFPsCompleteness.canvas, numPFPsCompleteness.current, numPFPsCompleteness.cheated, numPFPsCompleteness.dune, numPFPsCompleteness.uboone, numPFPsCompleteness.sbnd, 0, 820, 999, 999, (base_path + "numPFPsCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numPFPsCompleteness.current, numPFPsCompleteness.cheated, numPFPsCompleteness.dune, numPFPsCompleteness.uboone, numPFPsCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 92, 999, 999, (base_path + "numPFPsCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(ratioChosenSummedEnergyCompleteness.canvas, ratioChosenSummedEnergyCompleteness.current, ratioChosenSummedEnergyCompleteness.cheated, ratioChosenSummedEnergyCompleteness.dune, ratioChosenSummedEnergyCompleteness.uboone, ratioChosenSummedEnergyCompleteness.sbnd, 0, 820, 999, 999, (base_path + "ratioChosenSummedEnergyCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenSummedEnergyCompleteness.current, ratioChosenSummedEnergyCompleteness.cheated, ratioChosenSummedEnergyCompleteness.dune, ratioChosenSummedEnergyCompleteness.uboone, ratioChosenSummedEnergyCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 92, 999, 999, (base_path + "ratioChosenSummedEnergyCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioChosenTrueEnergyCompleteness.canvas, ratioChosenTrueEnergyCompleteness.current, ratioChosenTrueEnergyCompleteness.cheated, ratioChosenTrueEnergyCompleteness.dune, ratioChosenTrueEnergyCompleteness.uboone, ratioChosenTrueEnergyCompleteness.sbnd, 0, 300, 999, 999, (base_path + "ratioChosenTrueEnergyCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenTrueEnergyCompleteness.current, ratioChosenTrueEnergyCompleteness.cheated, ratioChosenTrueEnergyCompleteness.dune, ratioChosenTrueEnergyCompleteness.uboone, ratioChosenTrueEnergyCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 34, 999, 999, (base_path + "ratioChosenTrueEnergyCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioSummedTrueEnergyCompleteness.canvas, ratioSummedTrueEnergyCompleteness.current, ratioSummedTrueEnergyCompleteness.cheated, ratioSummedTrueEnergyCompleteness.dune, ratioSummedTrueEnergyCompleteness.uboone, ratioSummedTrueEnergyCompleteness.sbnd, 0, 300, 999, 999, (base_path + "ratioSummedTrueEnergyCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioSummedTrueEnergyCompleteness.current, ratioSummedTrueEnergyCompleteness.cheated, ratioSummedTrueEnergyCompleteness.dune, ratioSummedTrueEnergyCompleteness.uboone, ratioSummedTrueEnergyCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 34, 999, 999, (base_path + "ratioSummedTrueEnergyCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    styleDraw(angleDifferenceCompleteness.canvas, angleDifferenceCompleteness.current, angleDifferenceCompleteness.cheated, angleDifferenceCompleteness.dune, angleDifferenceCompleteness.uboone, angleDifferenceCompleteness.sbnd, 0, 260, 999, 999, (base_path + "angleDifferenceCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(angleDifferenceCompleteness.current, angleDifferenceCompleteness.cheated, angleDifferenceCompleteness.dune, angleDifferenceCompleteness.uboone, angleDifferenceCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 28, 999, 999, (base_path + "angleDifferenceCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
   
    styleDraw(EtrueThetaRecoCompleteness.canvas, EtrueThetaRecoCompleteness.current, EtrueThetaRecoCompleteness.cheated, EtrueThetaRecoCompleteness.dune, EtrueThetaRecoCompleteness.uboone, EtrueThetaRecoCompleteness.sbnd, 0, 120, 999, 999, (base_path + "EtrueThetaRecoCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(EtrueThetaRecoCompleteness.current, EtrueThetaRecoCompleteness.cheated, EtrueThetaRecoCompleteness.dune, EtrueThetaRecoCompleteness.uboone, EtrueThetaRecoCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 14, 999, 999, (base_path + "EtrueThetaRecoCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(EtrueThetaRecoCompleteness.current, EtrueThetaRecoCompleteness.cheated, EtrueThetaRecoCompleteness.dune, EtrueThetaRecoCompleteness.uboone, EtrueThetaRecoCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "EtrueThetaRecoCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{true}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoSumThetaTrueCompleteness.canvas, ERecoSumThetaTrueCompleteness.current, ERecoSumThetaTrueCompleteness.cheated, ERecoSumThetaTrueCompleteness.dune, ERecoSumThetaTrueCompleteness.uboone, ERecoSumThetaTrueCompleteness.sbnd, 0, 160, 999, 999, (base_path + "ERecoSumThetaTrueCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaTrueCompleteness.current, ERecoSumThetaTrueCompleteness.cheated, ERecoSumThetaTrueCompleteness.dune, ERecoSumThetaTrueCompleteness.uboone, ERecoSumThetaTrueCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 18, 999, 999, (base_path + "ERecoSumThetaTrueCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaTrueCompleteness.current, ERecoSumThetaTrueCompleteness.cheated, ERecoSumThetaTrueCompleteness.dune, ERecoSumThetaTrueCompleteness.uboone, ERecoSumThetaTrueCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoSumThetaTrueCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoHighestThetaTrueCompleteness.canvas, ERecoHighestThetaTrueCompleteness.current, ERecoHighestThetaTrueCompleteness.cheated, ERecoHighestThetaTrueCompleteness.dune, ERecoHighestThetaTrueCompleteness.uboone, ERecoHighestThetaTrueCompleteness.sbnd, 0, 1000, 999, 999, (base_path + "ERecoHighestThetaTrueCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaTrueCompleteness.current, ERecoHighestThetaTrueCompleteness.cheated, ERecoHighestThetaTrueCompleteness.dune, ERecoHighestThetaTrueCompleteness.uboone, ERecoHighestThetaTrueCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 100, 999, 999, (base_path + "ERecoHighestThetaTrueCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaTrueCompleteness.current, ERecoHighestThetaTrueCompleteness.cheated, ERecoHighestThetaTrueCompleteness.dune, ERecoHighestThetaTrueCompleteness.uboone, ERecoHighestThetaTrueCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoHighestThetaTrueCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoSumThetaRecoCompleteness.canvas, ERecoSumThetaRecoCompleteness.current, ERecoSumThetaRecoCompleteness.cheated, ERecoSumThetaRecoCompleteness.dune, ERecoSumThetaRecoCompleteness.uboone, ERecoSumThetaRecoCompleteness.sbnd, 0, 180, 999, 999, (base_path + "ERecoSumThetaRecoCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaRecoCompleteness.current, ERecoSumThetaRecoCompleteness.cheated, ERecoSumThetaRecoCompleteness.dune, ERecoSumThetaRecoCompleteness.uboone, ERecoSumThetaRecoCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 20, 999, 999, (base_path + "ERecoSumThetaRecoCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaRecoCompleteness.current, ERecoSumThetaRecoCompleteness.cheated, ERecoSumThetaRecoCompleteness.dune, ERecoSumThetaRecoCompleteness.uboone, ERecoSumThetaRecoCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoSumThetaRecoCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoHighestThetaRecoCompleteness.canvas, ERecoHighestThetaRecoCompleteness.current, ERecoHighestThetaRecoCompleteness.cheated, ERecoHighestThetaRecoCompleteness.dune, ERecoHighestThetaRecoCompleteness.uboone, ERecoHighestThetaRecoCompleteness.sbnd, 0, 180, 999, 999, (base_path + "ERecoHighestThetaRecoCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaRecoCompleteness.current, ERecoHighestThetaRecoCompleteness.cheated, ERecoHighestThetaRecoCompleteness.dune, ERecoHighestThetaRecoCompleteness.uboone, ERecoHighestThetaRecoCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 20, 999, 999, (base_path + "ERecoHighestThetaRecoCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaRecoCompleteness.current, ERecoHighestThetaRecoCompleteness.cheated, ERecoHighestThetaRecoCompleteness.dune, ERecoHighestThetaRecoCompleteness.uboone, ERecoHighestThetaRecoCompleteness.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 1, 999, 999, (base_path + "ERecoHighestThetaRecoCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");


    styleDraw(trueETheta2.canvas, trueETheta2.current, trueETheta2.cheated, trueETheta2.dune, trueETheta2.uboone, trueETheta2.sbnd, 0, 120, 999, 999, (base_path + "trueETheta2_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(trueETheta2.current, trueETheta2.cheated, trueETheta2.dune, trueETheta2.uboone, trueETheta2.sbnd, numEventsTrueElectron.current, numEventsTrueElectron.cheated, numEventsTrueElectron.dune, numEventsTrueElectron.uboone, numEventsTrueElectron.sbnd, 0, 14, 999, 999, (base_path + "trueETheta2_perc.pdf").c_str(), 0.56, 0.88, 0.70, 0.86, &drawLine, &right);
    efficiency(trueETheta2.current, trueETheta2.cheated, trueETheta2.dune, trueETheta2.uboone, trueETheta2.sbnd, numEventsTrueElectron.current, numEventsTrueElectron.cheated, numEventsTrueElectron.dune, numEventsTrueElectron.uboone, numEventsTrueElectron.sbnd, 0, 1, 999, 999, (base_path + "trueETheta2_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{true}#theta_{true}^{2} (MeV rad^{2})");
}
