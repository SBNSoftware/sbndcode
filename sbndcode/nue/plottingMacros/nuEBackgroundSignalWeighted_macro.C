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
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* currentCosmic;
    TH1F* ubooneCosmic;
    TH1F* nuECosmic;
    TH1F* currentSignal;
    TH1F* ubooneSignal;
    TH1F* nuESignal;
    TH1F* currentSignalFuzzy;
    TH1F* ubooneSignalFuzzy;
    TH1F* nuESignalFuzzy;
    TH1F* currentBNB;
    TH1F* ubooneBNB;
    TH1F* nuEBNB;
    TH1F* currentBNBFuzzy;
    TH1F* ubooneBNBFuzzy;
    TH1F* nuEBNBFuzzy;
} histGroup_struct;

struct weights_struct{
    double signalCurrent = 0;
    double signalUboone = 0;
    double signalNuE = 0;
    double BNBCurrent = 0;
    double BNBUboone = 0;
    double BNBNuE = 0;
    double cosmicsCurrent = 0;
    double cosmicsUboone = 0;
    double cosmicsNuE = 0;
};

typedef struct{

} PFP_struct;

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
    double etheta2;
    double dx;
    double dy;
    double dz;
} true_recoilElectron_struct;

typedef struct{
    std::vector<PFP_struct> PFPs;
} slice_struct;

histGroup_struct createHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);
    
    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_currentCosmic").c_str()),
        (TH1F*) base->Clone((baseName + "_ubooneCosmic").c_str()),
        (TH1F*) base->Clone((baseName + "_nuECosmic").c_str()),
        (TH1F*) base->Clone((baseName + "_currentSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_ubooneSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignal").c_str()),
        (TH1F*) base->Clone((baseName + "_currentSignalFuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_ubooneSignalFuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignalFuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_currentBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_ubooneBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_nuEBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_currentBNBFuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_ubooneBNBFuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_nuEBNBFuzzy").c_str()),
    };    
}

void styleDrawSignal(histGroup_struct hists, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    hists.currentSignal->SetLineWidth(2);
    hists.currentSignal->SetLineColor(kBlue+1);
    hists.ubooneSignal->SetLineWidth(2);
    hists.ubooneSignal->SetLineColor(kBlue-7);
    hists.nuESignal->SetLineWidth(2);
    hists.nuESignal->SetLineColor(kAzure+5);

    if((ymin != 999) && (ymax != 999)) hists.currentSignal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) hists.currentSignal->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({hists.currentSignal->GetMaximum(), hists.ubooneSignal->GetMaximum(), hists.nuESignal->GetMaximum()});
    double yminVal;
    if((ymin == 999) && (ymax == 999)){
        yminVal = 0;
        double ymaxVal = (maxYValue * 1.1);
        hists.currentSignal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    hists.currentSignal->Draw("hist");
    hists.ubooneSignal->Draw("histsame");
    hists.nuESignal->Draw("histsame");

    hists.currentSignal->SetStats(0);
    hists.currentSignal->GetXaxis()->SetTickLength(0.04);
    hists.currentSignal->GetYaxis()->SetTickLength(0.03);
    hists.currentSignal->GetXaxis()->SetTickSize(0.02);
    hists.currentSignal->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lxmax = 0;
    double Lymin = 0;
    double Lymax = 0;

    if(legendLocation == "topRight"){
        Lxmin = 0.48;
        Lymax = 0.863;
        Lxmax = 0.87;
        Lymin = 0.640;
    } else if(legendLocation == "topLeft"){
        Lxmin = 0.13;
        Lymax = 0.863;
        Lxmax = 0.52;
        Lymin = 0.640;
    } else if(legendLocation == "bottomRight"){
        Lxmin = 0.48;
        Lymax = 0.36;
        Lxmax = 0.87;
        Lymin = 0.137;
    } else if(legendLocation == "bottomLeft"){
        Lxmin = 0.13;
        Lymax = 0.36;
        Lxmax = 0.52;
        Lymin = 0.137;
    }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);

    legend->AddEntry(hists.currentSignal, "Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.ubooneSignal, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuESignal, "Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.12);
    legend->Draw();
    
    if(drawLine){
        TLine* line = new TLine(1.022, yminVal, 1.022, hists.currentSignal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        double ymaxValLine = maxYValue*0.95;
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, ymaxValLine, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.2, ymaxValLine, "2m_{e}");
        }

        latex->SetTextSize(0.035); 
        latex->SetTextAlign(11);
        latex->Draw("same");
    }

    hists.canvas->SaveAs(filename);
}

void styleDrawSignalAll(histGroup_struct hists, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    hists.currentCosmic->SetLineWidth(2);
    hists.currentCosmic->SetLineColor(kPink+9);
    hists.ubooneCosmic->SetLineWidth(2);
    hists.ubooneCosmic->SetLineColor(kPink+1);
    hists.nuECosmic->SetLineWidth(2);
    hists.nuECosmic->SetLineColor(kPink-2);

    hists.currentSignal->SetLineWidth(2);
    hists.currentSignal->SetLineColor(kBlue+1);
    hists.ubooneSignal->SetLineWidth(2);
    hists.ubooneSignal->SetLineColor(kBlue-7);
    hists.nuESignal->SetLineWidth(2);
    hists.nuESignal->SetLineColor(kAzure+5);

    hists.currentSignalFuzzy->SetLineWidth(2);
    hists.currentSignalFuzzy->SetLineColor(kGreen+3);
    hists.ubooneSignalFuzzy->SetLineWidth(2);
    hists.ubooneSignalFuzzy->SetLineColor(kGreen+1);
    hists.nuESignalFuzzy->SetLineWidth(2);
    hists.nuESignalFuzzy->SetLineColor(kGreen-7);
    
    hists.currentBNB->SetLineWidth(2);
    hists.currentBNB->SetLineColor(kOrange+7);
    hists.ubooneBNB->SetLineWidth(2);
    hists.ubooneBNB->SetLineColor(kOrange+6);
    hists.nuEBNB->SetLineWidth(2);
    hists.nuEBNB->SetLineColor(kOrange-5);
    
    hists.currentBNBFuzzy->SetLineWidth(2);
    hists.currentBNBFuzzy->SetLineColor(kViolet+1);
    hists.ubooneBNBFuzzy->SetLineWidth(2);
    hists.ubooneBNBFuzzy->SetLineColor(kViolet-7);
    hists.nuEBNBFuzzy->SetLineWidth(2);
    hists.nuEBNBFuzzy->SetLineColor(kViolet+4);

    if((ymin != 999) && (ymax != 999)) hists.currentSignal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) hists.currentSignal->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({hists.currentSignal->GetMaximum(), hists.ubooneSignal->GetMaximum(), hists.nuESignal->GetMaximum(), hists.currentCosmic->GetMaximum(), hists.ubooneCosmic->GetMaximum(), hists.nuECosmic->GetMaximum(), hists.currentSignalFuzzy->GetMaximum(), hists.ubooneSignalFuzzy->GetMaximum(), hists.nuESignalFuzzy->GetMaximum(), hists.currentBNB->GetMaximum(), hists.ubooneBNB->GetMaximum(), hists.nuEBNB->GetMaximum(), hists.currentBNBFuzzy->GetMaximum(), hists.ubooneBNBFuzzy->GetMaximum(), hists.nuEBNBFuzzy->GetMaximum()});
    double yminVal;
    if((ymin == 999) && (ymax == 999)){
        yminVal = 0;
        double ymaxVal = (maxYValue * 1.1);
        hists.currentSignal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    hists.currentSignal->Draw("hist");
    hists.ubooneSignal->Draw("histsame");
    hists.nuESignal->Draw("histsame");
    hists.currentCosmic->Draw("histsame");
    hists.ubooneCosmic->Draw("histsame");
    hists.nuECosmic->Draw("histsame");
    hists.currentSignalFuzzy->Draw("histsame");
    hists.ubooneSignalFuzzy->Draw("histsame");
    hists.nuESignalFuzzy->Draw("histsame");
    hists.currentBNB->Draw("histsame");
    hists.ubooneBNB->Draw("histsame");
    hists.nuEBNB->Draw("histsame");
    hists.currentBNBFuzzy->Draw("histsame");
    hists.ubooneBNBFuzzy->Draw("histsame");
    hists.nuEBNBFuzzy->Draw("histsame");

    hists.currentSignal->SetStats(0);
    hists.currentSignal->GetXaxis()->SetTickLength(0.04);
    hists.currentSignal->GetYaxis()->SetTickLength(0.03);
    hists.currentSignal->GetXaxis()->SetTickSize(0.02);
    hists.currentSignal->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lxmax = 0;
    double Lymin = 0;
    double Lymax = 0;

    if(legendLocation == "topRight"){
        Lxmin = 0.67;
        Lymax = 0.863;
        Lxmax = 0.87;
        Lymin = 0.640;
    } else if(legendLocation == "topLeft"){
        Lxmin = 0.13;
        Lymax = 0.863;
        Lxmax = 0.33;
        Lymin = 0.640;
    } else if(legendLocation == "bottomRight"){
        Lxmin = 0.67;
        Lymax = 0.36;
        Lxmax = 0.87;
        Lymin = 0.137;
    } else if(legendLocation == "bottomLeft"){
        Lxmin = 0.13;
        Lymax = 0.36;
        Lxmax = 0.33;
        Lymin = 0.137;
    }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(hists.currentSignal, "Signal, Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.ubooneSignal, "Signal, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuESignal, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->AddEntry(hists.currentSignalFuzzy, "Signal Fuzzy, Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.ubooneSignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuESignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->AddEntry(hists.currentBNB, "BNB, Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.ubooneBNB, "BNB, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuEBNB, "BNB, Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->AddEntry(hists.currentBNBFuzzy, "BNB Fuzzy, Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.ubooneBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuEBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->AddEntry(hists.currentCosmic, "Cosmic, Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.ubooneCosmic, "Cosmic, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuECosmic, "Cosmic, Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->SetTextSize(0.01);
    legend->SetMargin(0.12);
    legend->Draw();
    
    if(drawLine){
        TLine* line = new TLine(1.022, yminVal, 1.022, hists.currentSignal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        double ymaxValLine = maxYValue*0.95;
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, ymaxValLine, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.2, ymaxValLine, "2m_{e}");
        }

        latex->SetTextSize(0.035); 
        latex->SetTextAlign(11);
        latex->Draw("same");
    }

    hists.canvas->SaveAs(filename);
}

void nuEBackgroundSignalWeighted_macro(){

    TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nue/srcs/sbndcode/sbndcode/nue/NuEAnalyserOutput.root");
    std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsNew/";

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

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsCosmicsCurrent;
    
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsCosmicsUboone;
    
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsCosmicsNuE;

    double totalPOTSignalCurrent = 0;
    double totalPOTBNBCurrent = 0;
    double totalPOTCosmicsCurrent = 0;
    
    double totalPOTSignalUboone = 0;
    double totalPOTBNBUboone = 0;
    double totalPOTCosmicsUboone = 0;
    
    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;
    double totalPOTCosmicsNuE = 0;

    double cosmicSpillsSumCurrent = 0;
    double cosmicSpillsSumUboone = 0;
    double cosmicSpillsSumNuE = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 2) cosmicSpillsSumCurrent += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 0) cosmicSpillsSumUboone += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;

        std::pair<unsigned int, unsigned int> key = std::make_pair(subRunRun, subRunNumber);

        if(subRunSignal == 1){
            if(subRunDLCurrent == 2 && seenSubRunsSignalCurrent.find(key) == seenSubRunsSignalCurrent.end()){
                totalPOTSignalCurrent += subRunPOT;
                seenSubRunsSignalCurrent.insert(key);
            } else if(subRunDLCurrent == 0 && seenSubRunsSignalUboone.find(key) == seenSubRunsSignalUboone.end()){
                totalPOTSignalUboone += subRunPOT;
                seenSubRunsSignalUboone.insert(key);
            } else if(subRunDLCurrent == 5 && seenSubRunsSignalNuE.find(key) == seenSubRunsSignalNuE.end()){
                totalPOTSignalNuE += subRunPOT;
                seenSubRunsSignalNuE.insert(key);
            }
        } else if(subRunSignal == 2){
            if(subRunDLCurrent == 2 && seenSubRunsBNBCurrent.find(key) == seenSubRunsBNBCurrent.end()){
                totalPOTBNBCurrent += subRunPOT;
                seenSubRunsBNBCurrent.insert(key);
            } else if(subRunDLCurrent == 0 && seenSubRunsBNBUboone.find(key) == seenSubRunsBNBUboone.end()){
                totalPOTBNBUboone += subRunPOT;
                seenSubRunsBNBUboone.insert(key);
            } else if(subRunDLCurrent == 5 && seenSubRunsBNBNuE.find(key) == seenSubRunsBNBNuE.end()){
                totalPOTBNBNuE += subRunPOT;
                seenSubRunsBNBNuE.insert(key);
            }
        } else if(subRunSignal == 3){
            if(subRunDLCurrent == 2 && seenSubRunsCosmicsCurrent.find(key) == seenSubRunsCosmicsCurrent.end()){
                totalPOTCosmicsCurrent += subRunPOT;
                seenSubRunsCosmicsCurrent.insert(key);
            } else if(subRunDLCurrent == 0 && seenSubRunsCosmicsUboone.find(key) == seenSubRunsCosmicsUboone.end()){
                totalPOTCosmicsUboone += subRunPOT;
                seenSubRunsCosmicsUboone.insert(key);
            } else if(subRunDLCurrent == 5 && seenSubRunsCosmicsNuE.find(key) == seenSubRunsCosmicsNuE.end()){
                totalPOTCosmicsNuE += subRunPOT;
                seenSubRunsCosmicsNuE.insert(key);
            }
        }
    }

    double cosmicsPOTCurrent = cosmicSpillsSumCurrent * 5e12;
    double cosmicsPOTUboone = cosmicSpillsSumUboone * 5e12;
    double cosmicsPOTNuE = cosmicSpillsSumNuE * 5e12;

    weights_struct weights;
    weights.signalCurrent = totalPOTSignalCurrent/totalPOTSignalCurrent;
    weights.BNBCurrent = totalPOTSignalCurrent/totalPOTBNBCurrent;
    weights.cosmicsCurrent = totalPOTSignalCurrent/cosmicsPOTCurrent;
    
    weights.signalUboone = totalPOTSignalUboone/totalPOTSignalUboone;
    weights.BNBUboone = totalPOTSignalUboone/totalPOTBNBUboone;
    weights.cosmicsUboone = totalPOTSignalUboone/cosmicsPOTUboone;
    
    weights.signalNuE = totalPOTSignalNuE/totalPOTSignalNuE;
    weights.BNBNuE = totalPOTSignalNuE/totalPOTBNBNuE;
    weights.cosmicsNuE = totalPOTSignalNuE/cosmicsPOTNuE;

    std::cout << "BNB POT Current = " << totalPOTBNBCurrent << ", Uboone = " << totalPOTBNBUboone << ", Nu+E = " << totalPOTBNBNuE << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Signal POT Current = " << totalPOTSignalCurrent << ", Uboone = " << totalPOTSignalUboone << ", Nu+E = " << totalPOTSignalNuE << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Spills Cosmics POT Current = " << cosmicsPOTCurrent << ", Uboone = " << cosmicsPOTUboone << ", Nu+E = " << cosmicsPOTNuE << std::endl;  
    printf("Weights:\nCurrent: Signal = %f, BNB = %f, Cosmics = %f\nUboone: Signal = %f, BNB = %f, Cosmics = %f\nNu+E: Signal = %f, BNB = %f, Cosmics = %f\n", weights.signalCurrent, weights.BNBCurrent, weights.cosmicsCurrent, weights.signalUboone, weights.BNBUboone, weights.cosmicsUboone, weights.signalNuE, weights.BNBNuE, weights.cosmicsNuE);

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

    Long64_t numEntries = tree->GetEntries();

    auto trueETheta2 = createHistGroup("trueETheta2", "E_{true}#theta_{true}^{2}", "E_{true}#theta_{true}^{2} (MeV rad^{2})", 32, 0, 4.088);

    auto sliceCompleteness = createHistGroup("sliceCompleteness", "Slice Completeness", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessDist = createHistGroup("sliceCompletenessDist", "Slice Completeness (Not Weighted)", "Completeness", 102, 0, 1.02);
    auto slicePurity = createHistGroup("slicePurity", "Slice Purity", "Purity", 102, 0, 1.02);
    auto slicePurityDist = createHistGroup("slicePurityDist", "Slice Purity (Not Weighted)", "Purity", 102, 0, 1.02);
    auto sliceCRUMBSScore = createHistGroup("sliceCRUMBSScore", "CRUMBS Score of the Slice", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSScoreDist = createHistGroup("sliceCRUMBSScoreDist", "CRUMBS Score of the Slice (Not Weighted)", "CRUMBS Score", 25, -1, 1);

    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);
        
        // Looking at the true recoil electron
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->size() > 1) std::cout << "More than 1 true recoil electron in event!" << std::endl;
            if(truth_recoilElectronPDG->at(i) != -999999){
                // There is a true recoil electron in the event
                if(DLCurrent == 2) trueETheta2.currentSignal->Fill(truth_recoilElectronETheta2->at(i), weights.signalCurrent); 
                if(DLCurrent == 0) trueETheta2.ubooneSignal->Fill(truth_recoilElectronETheta2->at(i), weights.signalUboone); 
                if(DLCurrent == 5) trueETheta2.nuESignal->Fill(truth_recoilElectronETheta2->at(i), weights.signalNuE); 
            }
        }

        double weight = 0;
        if(signal == 1 && DLCurrent == 2) weight = weights.signalCurrent;
        if(signal == 2 && DLCurrent == 2) weight = weights.BNBCurrent;
        if(signal == 3 && DLCurrent == 2) weight = weights.cosmicsCurrent;
        if(signal == 1 && DLCurrent == 0) weight = weights.signalUboone;
        if(signal == 2 && DLCurrent == 0) weight = weights.BNBUboone;
        if(signal == 3 && DLCurrent == 0) weight = weights.cosmicsUboone;
        if(signal == 1 && DLCurrent == 5) weight = weights.signalNuE;
        if(signal == 2 && DLCurrent == 5) weight = weights.BNBNuE;
        if(signal == 3 && DLCurrent == 5) weight = weights.cosmicsNuE;

        // Looking at the reco slices
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            if(reco_sliceID->at(j) != -999999){
                // There is a reco slice in the event
                if(reco_sliceCategory->at(j) == 0){
                    // This is a cosmic slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentCosmic->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.currentCosmic->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.currentCosmic->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.currentCosmic->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.currentCosmic->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.currentCosmic->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneCosmic->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.ubooneCosmic->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.ubooneCosmic->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.ubooneCosmic->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.ubooneCosmic->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.ubooneCosmic->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuECosmic->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.nuECosmic->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.nuECosmic->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.nuECosmic->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.nuECosmic->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.nuECosmic->Fill(reco_sliceScore->at(j));
                    }
                } else if(reco_sliceCategory->at(j) == 1){
                    // This is a signal slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentSignal->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.currentSignal->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.currentSignal->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.currentSignal->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.currentSignal->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.currentSignal->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneSignal->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.ubooneSignal->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.ubooneSignal->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.ubooneSignal->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.ubooneSignal->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.ubooneSignal->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuESignal->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.nuESignal->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.nuESignal->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.nuESignal->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.nuESignal->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.nuESignal->Fill(reco_sliceScore->at(j));
                    }
                } else if(reco_sliceCategory->at(j) == 2){
                    // This is a fuzzy signal slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentSignalFuzzy->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.currentSignalFuzzy->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.currentSignalFuzzy->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.currentSignalFuzzy->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.currentSignalFuzzy->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.currentSignalFuzzy->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneSignalFuzzy->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.ubooneSignalFuzzy->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.ubooneSignalFuzzy->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.ubooneSignalFuzzy->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.ubooneSignalFuzzy->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.ubooneSignalFuzzy->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuESignalFuzzy->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.nuESignalFuzzy->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.nuESignalFuzzy->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.nuESignalFuzzy->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.nuESignalFuzzy->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.nuESignalFuzzy->Fill(reco_sliceScore->at(j));
                    }
                } else if(reco_sliceCategory->at(j) == 3){
                    // This is a BNB slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentBNB->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.currentBNB->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.currentBNB->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.currentBNB->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.currentBNB->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.currentBNB->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneBNB->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.ubooneBNB->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.ubooneBNB->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.ubooneBNB->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.ubooneBNB->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.ubooneBNB->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuEBNB->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.nuEBNB->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.nuEBNB->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.nuEBNB->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.nuEBNB->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.nuEBNB->Fill(reco_sliceScore->at(j));
                    }
                } else if(reco_sliceCategory->at(j) == 4){
                    // This is a fuzzy BNB slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentBNBFuzzy->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.currentBNBFuzzy->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.currentBNBFuzzy->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.currentBNBFuzzy->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.currentBNBFuzzy->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.currentBNBFuzzy->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneBNBFuzzy->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.ubooneBNBFuzzy->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.ubooneBNBFuzzy->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.ubooneBNBFuzzy->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.ubooneBNBFuzzy->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.ubooneBNBFuzzy->Fill(reco_sliceScore->at(j));
                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuEBNBFuzzy->Fill(reco_sliceCompleteness->at(j), weight);
                        slicePurity.nuEBNBFuzzy->Fill(reco_slicePurity->at(j), weight);
                        sliceCRUMBSScore.nuEBNBFuzzy->Fill(reco_sliceScore->at(j), weight);
                        sliceCompletenessDist.nuEBNBFuzzy->Fill(reco_sliceCompleteness->at(j));
                        slicePurityDist.nuEBNBFuzzy->Fill(reco_slicePurity->at(j));
                        sliceCRUMBSScoreDist.nuEBNBFuzzy->Fill(reco_sliceScore->at(j));
                    }
                }

                // aaaaa
            }
        
        }

        

 
    }

    int drawLine = 1;
    int left = 0;
    int right = 1;

    styleDrawSignal(trueETheta2, 999, 999, 999, 999, (base_path + "trueETheta2_weighted.pdf").c_str(), "bottomRight", &drawLine, &right);
    styleDrawSignalAll(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawSignalAll(sliceCompletenessDist, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawSignalAll(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawSignalAll(slicePurityDist, 999, 999, 999, 999, (base_path + "slicePurity_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawSignalAll(sliceCRUMBSScore, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawSignalAll(sliceCRUMBSScoreDist, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_dist.pdf").c_str(), "topRight", nullptr, &right);
}
