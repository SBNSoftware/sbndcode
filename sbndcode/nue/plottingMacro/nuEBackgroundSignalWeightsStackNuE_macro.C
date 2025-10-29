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

struct trueInfo{
    int numElectrons = 0;
    int numPhotons = 0;
    int numMuons = 0;
    int numPiZero = 0;
    int numChargedPi = 0;
    int numProtons = 0;
};

struct interactionCounter{
    int Unknown = 0;
    int QE = 0;
    int Res = 0;
    int DIS = 0;
    int Coh = 0;
    int CohElastic = 0;
    int ElectronScattering = 0;
    int IMDAnnihilation = 0;
    int InverseBetaDecay = 0;
    int GlashowResonance = 0;
    int AMNuGamma = 0;
    int MEC = 0;
    int Diffractive = 0;
    int EM = 0;
    int WeakMix = 0;
    int NuanceOffset = 0;
    int CCQE = 0;
    int NCQE = 0;
    int NuanceRes = 0;
    int CCDis = 0;
    int NCDis = 0;
    int NuEElastic = 0;
};

struct counter{
    int currentSignal = 0;
    int ubooneSignal = 0;
    int nueSignal = 0;
    int currentBNB = 0;
    int ubooneBNB = 0;
    int nueBNB = 0;
    int currentCosmics = 0;
    int ubooneCosmics = 0;
    int nueCosmics = 0;
};

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* currentSignal;
    TH1F* ubooneSignal;
    TH1F* nueSignal;
    TH1F* currentBNB;
    TH1F* ubooneBNB;
    TH1F* nueBNB;
    TH1F* currentCosmics;
    TH1F* ubooneCosmics;
    TH1F* nueCosmics;
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
    double numHits;
    double numMatchedHits;
    double numTrueHits;
    double truePDG;
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
    double genieMode;
    double genieInteraction;
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
        (TH1F*) base->Clone((baseName + "_currentSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_dlubooneSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_dlnueSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_currentBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_dlubooneBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_dlnueBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_currentCosmics").c_str()),
        (TH1F*) base->Clone((baseName + "_dlubooneCosmics").c_str()),
        (TH1F*) base->Clone((baseName + "_dlnueCosmics").c_str()),
    };    
}


void styleDraw(TCanvas* canvas, TH1F* currentSignal, TH1F* ubooneSignal, TH1F* nueSignal, TH1F* currentBNB, TH1F* ubooneBNB, TH1F* nueBNB, TH1F* currentCosmics, TH1F* ubooneCosmics, TH1F* nueCosmics, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, TPaveText* pt = nullptr, int* weighted = nullptr, int* drawLine = nullptr, int* linePos = nullptr, int* log = nullptr, TH1F* currentPurityHist = nullptr, TH1F* uboonePurityHist = nullptr, TH1F* nuePurityHist = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    gStyle->SetPalette(kAvocado);
    gROOT->ForceStyle();

    if(log && *log){ 
        gPad->SetLogy();
        currentSignal->SetMinimum(0.0000001);
        ubooneSignal->SetMinimum(0.0000001);
        nueSignal->SetMinimum(0.0000001);
        currentBNB->SetMinimum(0.0000001);
        ubooneBNB->SetMinimum(0.0000001);
        nueBNB->SetMinimum(0.0000001);
        currentCosmics->SetMinimum(0.0000001);
        ubooneCosmics->SetMinimum(0.0000001);
        nueCosmics->SetMinimum(0.0000001);
    }

    gPad->Update();

    currentSignal->SetLineWidth(2);
    currentSignal->SetLineColor(kBlue+1);
    ubooneSignal->SetLineWidth(2);
    ubooneSignal->SetLineColor(kBlue-7);
    nueSignal->SetLineWidth(2);
    nueSignal->SetLineColor(kAzure+5);

    currentBNB->SetLineWidth(2);
    currentBNB->SetLineColor(kOrange+7);
    ubooneBNB->SetLineWidth(2);
    ubooneBNB->SetLineColor(kOrange+6);
    nueBNB->SetLineWidth(2);
    nueBNB->SetLineColor(kOrange-5);

    currentCosmics->SetLineWidth(2);
    currentCosmics->SetLineColor(kPink+9);
    ubooneCosmics->SetLineWidth(2);
    ubooneCosmics->SetLineColor(kPink+1);
    nueCosmics->SetLineWidth(2);
    nueCosmics->SetLineColor(kPink-2);


    if(currentPurityHist){
        currentPurityHist->SetLineWidth(2);
        currentPurityHist->SetLineColor(kGreen+3);
        uboonePurityHist->SetLineWidth(2);
        uboonePurityHist->SetLineColor(kGreen+1);
        nuePurityHist->SetLineWidth(2);
        nuePurityHist->SetLineColor(kGreen);
    }

    if((ymin != 999) && (ymax != 999)) currentSignal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) currentSignal->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({currentSignal->GetMaximum(), ubooneSignal->GetMaximum(), nueSignal->GetMaximum(), currentBNB->GetMaximum(), ubooneBNB->GetMaximum(), nueBNB->GetMaximum(), currentCosmics->GetMaximum(), ubooneCosmics->GetMaximum(), nueCosmics->GetMaximum()});

    double yminVal;
    if((ymin == 999) && (ymax == 999)){
        yminVal = (log && *log) ? 0.1 : 0;
        double ymaxVal  = (log && *log) ? maxYValue*10000 : maxYValue*1.1;
        
        if(log && *log && currentPurityHist){
            yminVal = 0.000001;
            ymaxVal = 1;
        }
       
        currentSignal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
        //currentSignal->GetYaxis()->SetRangeUser(0, 1000);
    }

    currentSignal->Draw("hist");
    ubooneSignal->Draw("histsame");
    nueSignal->Draw("histsame");
    currentBNB->Draw("histsame");
    ubooneBNB->Draw("histsame");
    nueBNB->Draw("histsame");
    currentCosmics->Draw("histsame");
    ubooneCosmics->Draw("histsame");
    nueCosmics->Draw("histsame");

    if(currentPurityHist){
        currentPurityHist->Draw("histsame");
        uboonePurityHist->Draw("histsame");
        nuePurityHist->Draw("histsame");
    }

    currentSignal->SetStats(0);
    currentSignal->GetXaxis()->SetTickLength(0.04);
    currentSignal->GetYaxis()->SetTickLength(0.03);
    currentSignal->GetXaxis()->SetTickSize(0.02);
    currentSignal->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

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
    //auto legend = new TLegend(0.48, 0.39, 0.87, 0.167);
    if(!currentPurityHist){
        legend->AddEntry(ubooneSignal, "Signal, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nueSignal, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentSignal, "Signal, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(ubooneBNB, "BNB, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nueBNB, "BNB, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentBNB, "BNB, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(ubooneCosmics, "Cosmics, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nueCosmics, "Cosmics, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentCosmics, "Cosmics, Pandora BDT SBND (without Refinement)", "f");
        legend->SetTextSize(0.0225);
    }

    if(currentPurityHist){
        legend->AddEntry(ubooneSignal, "Signal Efficiency, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nueSignal, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentSignal, "Signal Efficiency, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(ubooneBNB, "BNB Rejection, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nueBNB, "BNB Rejection, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentBNB, "BNB Rejection, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(ubooneCosmics, "Intime Cosmic Rejection, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nueCosmics, "Intime Cosmic Rejection, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentCosmics, "Intime Cosmic Rejection, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(uboonePurityHist, "Selection Purity, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(nuePurityHist, "Selection Purity, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(currentPurityHist, "Selection Purity, Pandora BDT SBND (without Refinement)", "f");
        legend->SetTextSize(0.018);
    }
    
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, yminVal, 1.022, currentSignal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        double ymaxValLine  = (log && *log) ? maxYValue*1000 : maxYValue*0.95;
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, ymaxValLine, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.2, ymaxValLine, "2m_{e}");
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

void styleDrawIndividual(TCanvas* canvas, TH1F* current, TH1F* uboone, TH1F* nue, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, TPaveText* pt = nullptr, int* weighted = nullptr, int* drawLine = nullptr, int* linePos = nullptr, int* log = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    gStyle->SetPalette(kAvocado);
    gROOT->ForceStyle();

    if(log && *log){ 
        gPad->SetLogy();
        current->SetMinimum(0.00001);
        uboone->SetMinimum(0.00001);
        nue->SetMinimum(0.00001);
    }

    gPad->Update();

    current->SetLineWidth(2);
    current->SetLineColor(kBlack+3);
    uboone->SetLineWidth(2);
    uboone->SetLineColor(kBlack+1);
    nue->SetLineWidth(2);
    nue->SetLineColor(kGray);

    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({current->GetMaximum(), uboone->GetMaximum(), nue->GetMaximum()});

    double biggestYVal = 0;
    double smallestYVal = 1000;
    for(int i = 1; i <= current->GetNbinsX(); ++i){
        biggestYVal = std::max({biggestYVal, current->GetBinContent(i), uboone->GetBinContent(i), nue->GetBinContent(i)});
        smallestYVal = std::min({smallestYVal, current->GetBinContent(i), uboone->GetBinContent(i), nue->GetBinContent(i)});

    }

    double yminVal;
    if((ymin == 999) && (ymax == 999)){
        yminVal = (log && *log) ? 0.000001 : smallestYVal*0.9;
        double ymaxVal  = (log && *log) ? biggestYVal*1 : biggestYVal*1.1;
        if(biggestYVal > 0.9) ymaxVal = 1;
        current->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    current->Draw("hist");
    uboone->Draw("histsame");
    nue->Draw("histsame");

    current->SetStats(0);
    current->GetXaxis()->SetTickLength(0.04);
    current->GetYaxis()->SetTickLength(0.03);
    current->GetXaxis()->SetTickSize(0.02);
    current->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

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
    //auto legend = new TLegend(0.48, 0.39, 0.87, 0.167);
    legend->AddEntry(uboone, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(uboone, "Pandora Deep Learning: SBND Nu+E Tune", "f");
    legend->AddEntry(current, "Pandora BDT SBND (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, yminVal, 1.022, current->GetMaximum());
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


void weighted(TH1F* currentSignal, TH1F* ubooneSignal, TH1F* nueSignal, TH1F* currentBNB, TH1F* ubooneBNB, TH1F* nueBNB, TH1F* currentCosmics, TH1F* ubooneCosmics, TH1F* nueCosmics, double signalWeightCurrent, double BNBWeightCurrent, double cosmicsWeightCurrent, double signalWeightUboone, double BNBWeightUboone, double cosmicsWeightUboone, double signalWeightNue, double BNBWeightNue, double cosmicsWeightNue, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    TCanvas *weightedCanvas = new TCanvas("weighted_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentSignalWeighted = (TH1F*) currentSignal->Clone("weighted hist");
    currentSignalWeighted->Scale(signalWeightCurrent);
    currentSignalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 

    TH1F* ubooneSignalWeighted = (TH1F*) ubooneSignal->Clone("weighted hist");
    ubooneSignalWeighted->Scale(signalWeightUboone);
    ubooneSignalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)");
    
    TH1F* nueSignalWeighted = (TH1F*) nueSignal->Clone("weighted hist");
    nueSignalWeighted->Scale(signalWeightNue);
    nueSignalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)");
    
    TH1F* currentBNBWeighted = (TH1F*) currentBNB->Clone("weighted hist");
    currentBNBWeighted->Scale(BNBWeightCurrent);
    currentBNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 

    TH1F* ubooneBNBWeighted = (TH1F*) ubooneBNB->Clone("weighted hist");
    ubooneBNBWeighted->Scale(BNBWeightUboone);
    ubooneBNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)");
    
    TH1F* nueBNBWeighted = (TH1F*) nueBNB->Clone("weighted hist");
    nueBNBWeighted->Scale(BNBWeightNue);
    nueBNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)");

    TH1F* currentCosmicsWeighted = (TH1F*) currentCosmics->Clone("weighted hist");
    currentCosmicsWeighted->Scale(cosmicsWeightCurrent);
    currentCosmicsWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 

    TH1F* ubooneCosmicsWeighted = (TH1F*) ubooneCosmics->Clone("weighted hist");
    ubooneCosmicsWeighted->Scale(cosmicsWeightUboone);
    ubooneCosmicsWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)");
    
    TH1F* nueCosmicsWeighted = (TH1F*) nueCosmics->Clone("weighted hist");
    nueCosmicsWeighted->Scale(cosmicsWeightNue);
    nueCosmicsWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)");
    
    int funcValue = 1;
    int log = 1;

    styleDraw(weightedCanvas, currentSignalWeighted, ubooneSignalWeighted, nueSignalWeighted, currentBNBWeighted, ubooneBNBWeighted, nueBNBWeighted, currentCosmicsWeighted, ubooneCosmicsWeighted, nueCosmicsWeighted, ymin, ymax, xmin, xmax, filename, legendLocation, nullptr, &funcValue, drawLine, linePos, &log);
}

void efficiency(TH1F* currentSignal, TH1F* ubooneSignal, TH1F* nueSignal, TH1F* currentBNB, TH1F* ubooneBNB, TH1F* nueBNB, TH1F* currentCosmics, TH1F* ubooneCosmics, TH1F* nueCosmics, double sizeCurrentSignal, double sizeUbooneSignal, double sizeNueSignal, double sizeCurrentBNB, double sizeUbooneBNB, double sizeNueBNB, double sizeCurrentCosmics, double sizeUbooneCosmics, double sizeNueCosmics, double ymin, double ymax, double xmin, double xmax, const char* filename_eff, const char* filename_rej, const char* filename_effrejpur, const char* filename_pur, const char* filename_effpur, const std::string& legendLocation, double signalWeightCurrent, double BNBWeightCurrent, double cosmicsWeightCurrent, double signalWeightUboone, double BNBWeightUboone, double cosmicsWeightUboone, double signalWeightNue, double BNBWeightNue, double cosmicsWeightNue, int* drawLine = nullptr, int* linePos = nullptr, std::string xlabel = "", double efficiencyWay = 0.0, double* currentMaxEffPurBin = nullptr, double* ubooneMaxEffPurBin = nullptr, double* nueMaxEffPurBin = nullptr){
    TCanvas *efficiencyCanvas = new TCanvas("efficiency_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentSignalEff = (TH1F*) currentSignal->Clone("eff hist");
    currentSignalEff->Reset();
    currentSignalEff->GetYaxis()->SetTitle("Efficiency"); 
    currentSignalEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneSignalEff = (TH1F*) ubooneSignal->Clone("eff hist");
    ubooneSignalEff->Reset();
    ubooneSignalEff->GetYaxis()->SetTitle("Efficiency");
    ubooneSignalEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueSignalEff = (TH1F*) nueSignal->Clone("eff hist");
    nueSignalEff->Reset();
    nueSignalEff->GetYaxis()->SetTitle("Efficiency");
    nueSignalEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* currentBNBEff = (TH1F*) currentBNB->Clone("eff hist");
    currentBNBEff->Reset();
    currentBNBEff->GetYaxis()->SetTitle("Efficiency"); 
    currentBNBEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneBNBEff = (TH1F*) ubooneBNB->Clone("eff hist");
    ubooneBNBEff->Reset();
    ubooneBNBEff->GetYaxis()->SetTitle("Efficiency");
    ubooneBNBEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueBNBEff = (TH1F*) nueBNB->Clone("eff hist");
    nueBNBEff->Reset();
    nueBNBEff->GetYaxis()->SetTitle("Efficiency");
    nueBNBEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* currentCosmicsEff = (TH1F*) currentCosmics->Clone("eff hist");
    currentCosmicsEff->Reset();
    currentCosmicsEff->GetYaxis()->SetTitle("Efficiency"); 
    currentCosmicsEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneCosmicsEff = (TH1F*) ubooneCosmics->Clone("eff hist");
    ubooneCosmicsEff->Reset();
    ubooneCosmicsEff->GetYaxis()->SetTitle("Efficiency");
    ubooneCosmicsEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueCosmicsEff = (TH1F*) nueCosmics->Clone("eff hist");
    nueCosmicsEff->Reset();
    nueCosmicsEff->GetYaxis()->SetTitle("Efficiency");
    nueCosmicsEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TCanvas *rejectionCanvas = new TCanvas("rejection_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* currentSignalRej = (TH1F*) currentSignalEff->Clone("rej hist");
    currentSignalRej->Reset();
    currentSignalRej->GetYaxis()->SetTitle("Rejection");
    currentSignalRej->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneSignalRej = (TH1F*) ubooneSignalEff->Clone("rej hist");
    ubooneSignalRej->Reset();
    ubooneSignalRej->GetYaxis()->SetTitle("Rejection");
    ubooneSignalRej->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueSignalRej = (TH1F*) nueSignalEff->Clone("rej hist");
    nueSignalRej->Reset();
    nueSignalRej->GetYaxis()->SetTitle("Rejection");
    nueSignalRej->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* currentBNBRej = (TH1F*) currentBNBEff->Clone("rej hist");
    currentBNBRej->Reset();
    std::string currentBNBRejTitle = currentBNBRej->GetTitle();
    currentBNBRejTitle += ", BNB";
    currentBNBRej->SetTitle(currentBNBRejTitle.c_str());
    currentBNBRej->GetYaxis()->SetTitle("Rejection");
    currentBNBRej->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneBNBRej = (TH1F*) ubooneBNBEff->Clone("rej hist");
    ubooneBNBRej->Reset();
    std::string ubooneBNBRejTitle = ubooneBNBRej->GetTitle();
    ubooneBNBRejTitle += ", BNB";
    ubooneBNBRej->SetTitle(ubooneBNBRejTitle.c_str());
    ubooneBNBRej->GetYaxis()->SetTitle("Rejection");
    ubooneBNBRej->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueBNBRej = (TH1F*) nueBNBEff->Clone("rej hist");
    nueBNBRej->Reset();
    std::string nueBNBRejTitle = nueBNBRej->GetTitle();
    nueBNBRejTitle += ", BNB";
    nueBNBRej->SetTitle(nueBNBRejTitle.c_str());
    nueBNBRej->GetYaxis()->SetTitle("Rejection");
    nueBNBRej->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* currentCosmicsRej = (TH1F*) currentCosmicsEff->Clone("rej hist");
    currentCosmicsRej->Reset();
    std::string currentCosmicRejTitle = currentCosmicsRej->GetTitle();
    currentCosmicRejTitle += ", Intime Cosmics";
    currentCosmicsRej->SetTitle(currentCosmicRejTitle.c_str());
    currentCosmicsRej->GetYaxis()->SetTitle("Rejection");
    currentCosmicsRej->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneCosmicsRej = (TH1F*) ubooneCosmicsEff->Clone("rej hist");
    ubooneCosmicsRej->Reset();
    std::string ubooneCosmicRejTitle = ubooneCosmicsRej->GetTitle();
    ubooneCosmicRejTitle += ", Intime Cosmics";
    ubooneCosmicsRej->SetTitle(ubooneCosmicRejTitle.c_str());
    ubooneCosmicsRej->GetYaxis()->SetTitle("Rejection");
    ubooneCosmicsRej->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueCosmicsRej = (TH1F*) nueCosmicsEff->Clone("rej hist");
    nueCosmicsRej->Reset();
    std::string nueCosmicRejTitle = nueCosmicsRej->GetTitle();
    nueCosmicRejTitle += ", Intime Cosmics";
    nueCosmicsRej->SetTitle(nueCosmicRejTitle.c_str());
    nueCosmicsRej->GetYaxis()->SetTitle("Rejection");
    nueCosmicsRej->GetXaxis()->SetTitle(xlabel.c_str());
   
    TCanvas *purityCanvas = new TCanvas("purity_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* currentPur = (TH1F*) currentSignalEff->Clone("pur hist");
    currentPur->Reset();
    currentPur->GetYaxis()->SetTitle("Purity");
    currentPur->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* uboonePur = (TH1F*) ubooneSignalEff->Clone("pur hist");
    uboonePur->Reset();
    uboonePur->GetYaxis()->SetTitle("Purity");
    uboonePur->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nuePur = (TH1F*) nueSignalEff->Clone("pur hist");
    nuePur->Reset();
    nuePur->GetYaxis()->SetTitle("Purity");
    nuePur->GetXaxis()->SetTitle(xlabel.c_str());

    TCanvas *effpurCanvas = new TCanvas("effpur_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* currentEffPur = (TH1F*) currentSignalEff->Clone("current eff pur hist");
    //TH1F* currentEffPur = new TH1F("current_eff_pur_hist", "", currentSignalEff->GetNbinsX(), currentSignalEff->GetXaxis()->GetXmin(), currentSignalEff->GetXaxis()->GetXmax()); 
    currentEffPur->Reset();
    currentEffPur->GetYaxis()->SetTitle("Efficiency x Purity");
    currentEffPur->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneEffPur = (TH1F*) ubooneSignalEff->Clone("uboone eff pur hist");
    //TH1F* ubooneEffPur = new TH1F("uboone_eff_pur_hist", "", ubooneSignalEff->GetNbinsX(), ubooneSignalEff->GetXaxis()->GetXmin(), ubooneSignalEff->GetXaxis()->GetXmax());
    ubooneEffPur->Reset();
    ubooneEffPur->GetYaxis()->SetTitle("Efficiency x Purity");
    ubooneEffPur->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueEffPur = (TH1F*) nueSignalEff->Clone("nue eff pur hist");
    //TH1F* ubooneEffPur = new TH1F("uboone_eff_pur_hist", "", ubooneSignalEff->GetNbinsX(), ubooneSignalEff->GetXaxis()->GetXmin(), ubooneSignalEff->GetXaxis()->GetXmax());
    nueEffPur->Reset();
    nueEffPur->GetYaxis()->SetTitle("Efficiency x Purity");
    nueEffPur->GetXaxis()->SetTitle(xlabel.c_str());
    
    TCanvas *EffRejPurCanvas = new TCanvas("efficiency_rejection_canvas", "Graph Draw Options", 200, 10, 600, 400);
    
    TCanvas *numEventsCanvas = new TCanvas("numEvents_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* numEventsCurrentSignal = (TH1F*) currentSignal->Clone("numEventsCurrentSignal hist");
    numEventsCurrentSignal->Reset();
    numEventsCurrentSignal->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsCurrentSignal->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* numEventsUbooneSignal = (TH1F*) currentSignal->Clone("numEventsUbooneSignal hist");
    numEventsUbooneSignal->Reset();
    numEventsUbooneSignal->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsUbooneSignal->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* numEventsNueSignal = (TH1F*) currentSignal->Clone("numEventsNueSignal hist");
    numEventsNueSignal->Reset();
    numEventsNueSignal->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsNueSignal->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* numEventsCurrentBNB = (TH1F*) currentSignal->Clone("numEventsCurrentBNB hist");
    numEventsCurrentBNB->Reset();
    numEventsCurrentBNB->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsCurrentBNB->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* numEventsUbooneBNB = (TH1F*) currentSignal->Clone("numEventsUbooneBNB hist");
    numEventsUbooneBNB->Reset();
    numEventsUbooneBNB->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsUbooneBNB->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* numEventsNueBNB = (TH1F*) currentSignal->Clone("numEventsNueBNB hist");
    numEventsNueBNB->Reset();
    numEventsNueBNB->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsNueBNB->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* numEventsCurrentCosmics = (TH1F*) currentSignal->Clone("numEventsCurrentCosmics hist");
    numEventsCurrentCosmics->Reset();
    numEventsCurrentCosmics->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsCurrentCosmics->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* numEventsUbooneCosmics = (TH1F*) currentSignal->Clone("numEventsUbooneCosmics hist");
    numEventsUbooneCosmics->Reset();
    numEventsUbooneCosmics->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsUbooneCosmics->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* numEventsNueCosmics = (TH1F*) currentSignal->Clone("numEventsNueCosmics hist");
    numEventsNueCosmics->Reset();
    numEventsNueCosmics->GetYaxis()->SetTitle("Number of Events that Pass Cut"); 
    numEventsNueCosmics->GetXaxis()->SetTitle(xlabel.c_str());
    
    int numBins = currentSignal->GetNbinsX();
    double currentSignalSum = 0.0;
    double ubooneSignalSum = 0.0;
    double nueSignalSum = 0.0;
    double currentBNBSum = 0.0;
    double ubooneBNBSum = 0.0;
    double nueBNBSum = 0.0;
    double currentCosmicsSum = 0.0;
    double ubooneCosmicsSum = 0.0;
    double nueCosmicsSum = 0.0;

    double currentSignalTotal = 0.0;
    double ubooneSignalTotal = 0.0;
    double nueSignalTotal = 0.0;
    double currentBNBTotal = 0.0;
    double ubooneBNBTotal = 0.0;
    double nueBNBTotal = 0.0;
    double currentCosmicsTotal = 0.0;
    double ubooneCosmicsTotal = 0.0;
    double nueCosmicsTotal = 0.0;

    // efficiencyWay == -1 includes everything to the right of the cut
    if(efficiencyWay == -1){
        for(int i = 1; i <= numBins; ++i){
            currentSignalTotal += currentSignal->GetBinContent(i);
            ubooneSignalTotal += ubooneSignal->GetBinContent(i);
            nueSignalTotal += nueSignal->GetBinContent(i);
            currentBNBTotal += currentBNB->GetBinContent(i);
            ubooneBNBTotal += ubooneBNB->GetBinContent(i);
            nueBNBTotal += nueBNB->GetBinContent(i);
            currentCosmicsTotal += currentCosmics->GetBinContent(i);
            ubooneCosmicsTotal += ubooneCosmics->GetBinContent(i);
            nueCosmicsTotal += nueCosmics->GetBinContent(i);
        }

        std::cout << "TOTALS:" << std::endl;
        std::cout << "Signal: Current = " << currentSignalTotal << ", Uboone = " << ubooneSignalTotal << ", Nu+E = " << nueSignalTotal << std::endl;
        std::cout << "BNB: Current = " << currentBNBTotal << ", Uboone = " << ubooneBNBTotal << ", Nu+E = " << nueBNBTotal << std::endl;
        std::cout << "Cosmics: Current = " << currentCosmicsTotal << ", Uboone = " << ubooneCosmicsTotal << ", Nu+E = " << nueCosmicsTotal << std::endl;
    
        for(int i = 1; i <= numBins; ++i){
            currentSignalSum += currentSignal->GetBinContent(i);
            ubooneSignalSum += ubooneSignal->GetBinContent(i);
            nueSignalSum += nueSignal->GetBinContent(i);
            currentBNBSum += currentBNB->GetBinContent(i);
            ubooneBNBSum += ubooneBNB->GetBinContent(i);
            nueBNBSum += nueBNB->GetBinContent(i);
            currentCosmicsSum += currentCosmics->GetBinContent(i);
            ubooneCosmicsSum += ubooneCosmics->GetBinContent(i);
            nueCosmicsSum += nueCosmics->GetBinContent(i);

            numEventsCurrentSignal->SetBinContent(i, currentSignalTotal - currentSignalSum);
            numEventsUbooneSignal->SetBinContent(i, ubooneSignalTotal - ubooneSignalSum);
            numEventsNueSignal->SetBinContent(i, nueSignalTotal - nueSignalSum);
            numEventsCurrentBNB->SetBinContent(i, currentBNBTotal - currentBNBSum);
            numEventsUbooneBNB->SetBinContent(i, ubooneBNBTotal - ubooneBNBSum);
            numEventsNueBNB->SetBinContent(i, nueBNBTotal - nueBNBSum);
            numEventsCurrentCosmics->SetBinContent(i, currentCosmicsTotal - currentCosmicsSum);
            numEventsUbooneCosmics->SetBinContent(i, ubooneCosmicsTotal - ubooneCosmicsSum);
            numEventsNueCosmics->SetBinContent(i, nueCosmicsTotal - nueCosmicsSum);

            double currentSignalEffValue = (currentSignalTotal - currentSignalSum)/sizeCurrentSignal;
            double ubooneSignalEffValue = (ubooneSignalTotal - ubooneSignalSum)/sizeUbooneSignal;
            double nueSignalEffValue = (nueSignalTotal - nueSignalSum)/sizeNueSignal;
            double currentBNBEffValue = (currentBNBTotal - currentBNBSum)/sizeCurrentBNB;
            double ubooneBNBEffValue = (ubooneBNBTotal - ubooneBNBSum)/sizeUbooneBNB;
            double nueBNBEffValue = (nueBNBTotal - nueBNBSum)/sizeNueBNB;
            double currentCosmicsEffValue = (currentCosmicsTotal - currentCosmicsSum)/sizeCurrentCosmics;
            double ubooneCosmicsEffValue = (ubooneCosmicsTotal - ubooneCosmicsSum)/sizeUbooneCosmics;
            double nueCosmicsEffValue = (nueCosmicsTotal - nueCosmicsSum)/sizeNueCosmics;

            currentSignalEff->SetBinContent(i, currentSignalEffValue);
            ubooneSignalEff->SetBinContent(i, ubooneSignalEffValue);
            nueSignalEff->SetBinContent(i, nueSignalEffValue);
            currentBNBEff->SetBinContent(i, currentBNBEffValue);
            ubooneBNBEff->SetBinContent(i, ubooneBNBEffValue);
            nueBNBEff->SetBinContent(i, nueBNBEffValue);
            currentCosmicsEff->SetBinContent(i, currentCosmicsEffValue);
            ubooneCosmicsEff->SetBinContent(i, ubooneCosmicsEffValue);
            nueCosmicsEff->SetBinContent(i, nueCosmicsEffValue);
        
            currentSignalRej->SetBinContent(i, 1-currentSignalEffValue);
            ubooneSignalRej->SetBinContent(i, 1-ubooneSignalEffValue);
            nueSignalRej->SetBinContent(i, 1-nueSignalEffValue);
            currentBNBRej->SetBinContent(i, 1-currentBNBEffValue);
            ubooneBNBRej->SetBinContent(i, 1-ubooneBNBEffValue);
            nueBNBRej->SetBinContent(i, 1-nueBNBEffValue);
            currentCosmicsRej->SetBinContent(i, 1-currentCosmicsEffValue);
            ubooneCosmicsRej->SetBinContent(i, 1-ubooneCosmicsEffValue);
            nueCosmicsRej->SetBinContent(i, 1-nueCosmicsEffValue);
   
            double sumCurrentBackground = ((currentBNBTotal - currentBNBSum) * BNBWeightCurrent) + ((currentCosmicsTotal - currentCosmicsSum) * cosmicsWeightCurrent);
            double sumUbooneBackground = ((ubooneBNBTotal - ubooneBNBSum) * BNBWeightUboone) + ((ubooneCosmicsTotal - ubooneCosmicsSum) * cosmicsWeightUboone);
            double sumNueBackground = ((nueBNBTotal - nueBNBSum) * BNBWeightNue) + ((nueCosmicsTotal - nueCosmicsSum) * cosmicsWeightNue);
            double totalSumCurrent = (sumCurrentBackground + (currentSignalTotal - currentSignalSum));
            double totalSumUboone = (sumUbooneBackground + (ubooneSignalTotal - ubooneSignalSum));
            double totalSumNue = (sumNueBackground + (nueSignalTotal - nueSignalSum));

            double currentPurityValue;
            double uboonePurityValue;
            double nuePurityValue;
        
            if(currentSignalSum != 0 && sumCurrentBackground != 0){
                currentPurityValue = ((currentSignalTotal - currentSignalSum) * signalWeightCurrent)/ totalSumCurrent;
            } else{
                currentPurityValue = 0;
            }

            if(ubooneSignalSum != 0 && sumUbooneBackground != 0){
                uboonePurityValue = ((ubooneSignalTotal - ubooneSignalSum) * signalWeightUboone) / totalSumUboone;
            } else{
                uboonePurityValue = 0;
            }
            
            if(nueSignalSum != 0 && sumNueBackground != 0){
                nuePurityValue = ((nueSignalTotal - nueSignalSum) * signalWeightNue) / totalSumNue;
            } else{
                nuePurityValue = 0;
            }
            
            double currentEffPurValue = (currentSignalEffValue * currentPurityValue);
            double ubooneEffPurValue =  (ubooneSignalEffValue * uboonePurityValue);
            double nueEffPurValue =  (nueSignalEffValue * nuePurityValue);
        
            currentPur->SetBinContent(i, currentPurityValue);
            uboonePur->SetBinContent(i, uboonePurityValue);
            nuePur->SetBinContent(i, nuePurityValue);

            currentEffPur->SetBinContent(i, currentEffPurValue);
            ubooneEffPur->SetBinContent(i, ubooneEffPurValue);
            nueEffPur->SetBinContent(i, nueEffPurValue);
        }
    }

    // efficiencyWay == 1 includes everything to the left of the cut
    if(efficiencyWay == 1){
        for(int i = 1; i <= numBins; ++i){
            currentSignalSum += currentSignal->GetBinContent(i);
            ubooneSignalSum += ubooneSignal->GetBinContent(i);
            nueSignalSum += nueSignal->GetBinContent(i);
            currentBNBSum += currentBNB->GetBinContent(i);
            ubooneBNBSum += ubooneBNB->GetBinContent(i);
            nueBNBSum += nueBNB->GetBinContent(i);
            currentCosmicsSum += currentCosmics->GetBinContent(i);
            ubooneCosmicsSum += ubooneCosmics->GetBinContent(i);
            nueCosmicsSum += nueCosmics->GetBinContent(i);

            numEventsCurrentSignal->SetBinContent(i, currentSignalSum);
            numEventsUbooneSignal->SetBinContent(i, ubooneSignalSum);
            numEventsNueSignal->SetBinContent(i, nueSignalSum);
            numEventsCurrentBNB->SetBinContent(i, currentBNBSum);
            numEventsUbooneBNB->SetBinContent(i, ubooneBNBSum);
            numEventsNueBNB->SetBinContent(i, nueBNBSum);
            numEventsCurrentCosmics->SetBinContent(i, currentCosmicsSum);
            numEventsUbooneCosmics->SetBinContent(i, ubooneCosmicsSum);
            numEventsNueCosmics->SetBinContent(i, nueCosmicsSum);

            double currentSignalEffValue = currentSignalSum/sizeCurrentSignal;
            double ubooneSignalEffValue = ubooneSignalSum/sizeUbooneSignal;
            double nueSignalEffValue = nueSignalSum/sizeNueSignal;
            double currentBNBEffValue = currentBNBSum/sizeCurrentBNB;
            double ubooneBNBEffValue = ubooneBNBSum/sizeUbooneBNB;
            double nueBNBEffValue = nueBNBSum/sizeNueBNB;
            double currentCosmicsEffValue = currentCosmicsSum/sizeCurrentCosmics;
            double ubooneCosmicsEffValue = ubooneCosmicsSum/sizeUbooneCosmics;
            double nueCosmicsEffValue = nueCosmicsSum/sizeNueCosmics;

            currentSignalEff->SetBinContent(i, currentSignalEffValue);
            ubooneSignalEff->SetBinContent(i, ubooneSignalEffValue);
            nueSignalEff->SetBinContent(i, nueSignalEffValue);
            currentBNBEff->SetBinContent(i, currentBNBEffValue);
            ubooneBNBEff->SetBinContent(i, ubooneBNBEffValue);
            nueBNBEff->SetBinContent(i, nueBNBEffValue);
            currentCosmicsEff->SetBinContent(i, currentCosmicsEffValue);
            ubooneCosmicsEff->SetBinContent(i, ubooneCosmicsEffValue);
            nueCosmicsEff->SetBinContent(i, nueCosmicsEffValue);
        
            currentSignalRej->SetBinContent(i, 1-currentSignalEffValue);
            ubooneSignalRej->SetBinContent(i, 1-ubooneSignalEffValue);
            nueSignalRej->SetBinContent(i, 1-nueSignalEffValue);
            currentBNBRej->SetBinContent(i, 1-currentBNBEffValue);
            ubooneBNBRej->SetBinContent(i, 1-ubooneBNBEffValue);
            nueBNBRej->SetBinContent(i, 1-nueBNBEffValue);
            currentCosmicsRej->SetBinContent(i, 1-currentCosmicsEffValue);
            ubooneCosmicsRej->SetBinContent(i, 1-ubooneCosmicsEffValue);
            nueCosmicsRej->SetBinContent(i, 1-nueCosmicsEffValue);
   
            double sumCurrentBackground = (currentBNBSum * BNBWeightCurrent) + (currentCosmicsSum * cosmicsWeightCurrent);
            double sumUbooneBackground = (ubooneBNBSum * BNBWeightUboone) + (ubooneCosmicsSum * cosmicsWeightUboone);
            double sumNueBackground = (nueBNBSum * BNBWeightNue) + (ubooneCosmicsSum * cosmicsWeightNue);
            double totalSumCurrent = (sumCurrentBackground + currentSignalSum);
            double totalSumUboone = (sumUbooneBackground + ubooneSignalSum);
            double totalSumNue = (sumNueBackground + nueSignalSum);

            double currentPurityValue;
            double uboonePurityValue;
            double nuePurityValue;
        
            if(currentSignalSum != 0 && sumCurrentBackground != 0){
                currentPurityValue = (currentSignalSum * signalWeightCurrent)/ totalSumCurrent;
            } else{
                currentPurityValue = 0;
            }

            if(ubooneSignalSum != 0 && sumUbooneBackground != 0){
                uboonePurityValue = (ubooneSignalSum * signalWeightUboone) / totalSumUboone;
            } else{
                uboonePurityValue = 0;
            }
            
            if(nueSignalSum != 0 && sumNueBackground != 0){
                nuePurityValue = (nueSignalSum * signalWeightNue) / totalSumNue;
            } else{
                nuePurityValue = 0;
            }
            
            double currentEffPurValue = (currentSignalEffValue * currentPurityValue);
            double ubooneEffPurValue =  (ubooneSignalEffValue * uboonePurityValue);
            double nueEffPurValue =  (nueSignalEffValue * nuePurityValue);
        
            currentPur->SetBinContent(i, currentPurityValue);
            uboonePur->SetBinContent(i, uboonePurityValue);
            nuePur->SetBinContent(i, nuePurityValue);

            currentEffPur->SetBinContent(i, currentEffPurValue);
            ubooneEffPur->SetBinContent(i, ubooneEffPurValue);
            nueEffPur->SetBinContent(i, nueEffPurValue);
        }
    }

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

    if(legendLocation == "topRight"){
        Lxmin = 0.48;
        Lymax = 0.97;
        Lxmax = 0.87;
        Lymin = 0.747;
    } else if(legendLocation == "topLeft"){
        Lxmin = 0.13;
        Lymax = 0.97;
        Lxmax = 0.52;
        Lymin = 0.747;
    } else if(legendLocation == "bottomRight"){
        Lxmin = 0.48;
        Lymax = 0.39;
        Lxmax = 0.87;
        Lymin = 0.167;
    } else if(legendLocation == "bottomLeft"){
        Lxmin = 0.13;
        Lymax = 0.39;
        Lxmax = 0.52;
        Lymin = 0.167;
    }

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of Signal DL Uboone Entries: %d", (int)sizeUbooneSignal));
    pt->AddText(Form("Number of Signal DL Nu+e Entries: %d", (int)sizeNueSignal));
    pt->AddText(Form("Number of Signal Current Entries: %d", (int)sizeCurrentSignal));
    pt->AddText(Form("Number of BNB DL Uboone Entries: %d", (int)sizeUbooneBNB));
    pt->AddText(Form("Number of BNB DL Nu+E Entries: %d", (int)sizeNueBNB));
    pt->AddText(Form("Number of BNB Current Entries: %d", (int)sizeCurrentBNB));
    pt->AddText(Form("Number of Cosmics DL Uboone Entries: %d", (int)sizeUbooneCosmics));
    pt->AddText(Form("Number of Cosmics DL Nu+e Entries: %d", (int)sizeNueCosmics));
    pt->AddText(Form("Number of Cosmics Current Entries: %d", (int)sizeCurrentCosmics));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 

    int funcValue = 1;
    int log = 1;

    currentPur->GetYaxis()->SetNoExponent(false);
    currentPur->GetYaxis()->SetMaxDigits(2);
    uboonePur->GetYaxis()->SetNoExponent(false);
    uboonePur->GetYaxis()->SetMaxDigits(2);
    nuePur->GetYaxis()->SetNoExponent(false);
    nuePur->GetYaxis()->SetMaxDigits(2);
    
    currentEffPur->GetYaxis()->SetNoExponent(false);
    currentEffPur->GetYaxis()->SetMaxDigits(2);
    ubooneEffPur->GetYaxis()->SetNoExponent(false);
    ubooneEffPur->GetYaxis()->SetMaxDigits(2);
    nueEffPur->GetYaxis()->SetNoExponent(false);
    nueEffPur->GetYaxis()->SetMaxDigits(2);

    TCanvas *rejBNBCanvas = new TCanvas("rejBNB_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TCanvas *rejCosmicCanvas = new TCanvas("rejCosmic_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TCanvas *effSignalCanvas = new TCanvas("effSignal_canvas", "Graph Draw Options", 200, 10, 600, 400);

    std::string filename_rej_str(filename_rej);
    std::string filename_rejBNB_str = filename_rej_str;
    std::string filename_rejCosmics_str = filename_rej_str;
    std::string filename_numEvents_str = filename_rej_str;
    size_t posBNB = filename_rejBNB_str.find("rej");
    if (posBNB != std::string::npos) filename_rejBNB_str.replace(posBNB, 3, "rejBNB");
    if (posBNB != std::string::npos) filename_numEvents_str.replace(posBNB, 3, "numEvents");
    size_t posCosmics = filename_rejCosmics_str.find("rej");
    if (posCosmics != std::string::npos) filename_rejCosmics_str.replace(posCosmics, 3, "rejCosmics");
    const char* filename_rejBNB = filename_rejBNB_str.c_str();
    const char* filename_rejCosmics = filename_rejCosmics_str.c_str();
    const char* filename_numEvents = filename_numEvents_str.c_str();
   
    std::string filename_eff_str(filename_eff); 
    std::string filename_effSignal_str = filename_eff_str;
    size_t posSignal = filename_effSignal_str.find("eff");
    if (posSignal != std::string::npos) filename_effSignal_str.replace(posSignal, 3, "effSignal");
    const char* filename_effSignal = filename_effSignal_str.c_str();

    styleDraw(efficiencyCanvas, currentSignalEff, ubooneSignalEff, nueSignalEff, currentBNBEff, ubooneBNBEff, nueBNBEff, currentCosmicsEff, ubooneCosmicsEff, nueCosmicsEff, ymin, ymax, xmin, xmax, filename_eff, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDraw(rejectionCanvas, currentSignalRej, ubooneSignalRej, nueSignalRej, currentBNBRej, ubooneBNBRej, nueBNBRej, currentCosmicsRej, ubooneCosmicsRej, nueCosmicsRej, ymin, ymax, xmin, xmax, filename_rej, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDraw(EffRejPurCanvas, currentSignalEff, ubooneSignalEff, nueSignalEff, currentBNBRej, ubooneBNBRej, nueBNBRej, currentCosmicsRej, ubooneCosmicsRej, nueCosmicsRej, 999, 999, 999, 999, filename_effrejpur, legendLocation, pt = nullptr, &funcValue, drawLine, linePos, &log, currentPur, uboonePur, nuePur);
    styleDrawIndividual(purityCanvas, currentPur, uboonePur, nuePur, 999, 999, 999, 999, filename_pur, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawIndividual(effpurCanvas, currentEffPur, ubooneEffPur, nueEffPur, 999, 999, 999, 999, filename_effpur, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawIndividual(rejBNBCanvas, currentBNBRej, ubooneBNBRej, nueBNBRej, 999, 999, 999, 999, filename_rejBNB, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawIndividual(rejCosmicCanvas, currentCosmicsRej, ubooneCosmicsRej, nueCosmicsRej, 999, 999, 999, 999, filename_rejCosmics, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawIndividual(effSignalCanvas, currentSignalEff, ubooneSignalEff, nueSignalEff, 999, 999, 999, 999, filename_effSignal, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDraw(numEventsCanvas, numEventsCurrentSignal, numEventsUbooneSignal, numEventsNueSignal, numEventsCurrentBNB, numEventsUbooneBNB, numEventsNueBNB, numEventsCurrentCosmics, numEventsUbooneCosmics, numEventsNueCosmics, 999, 999, 999, 999, filename_numEvents, legendLocation, pt = nullptr, &funcValue, drawLine, linePos);

    if(currentMaxEffPurBin) *currentMaxEffPurBin = currentEffPur->GetBinCenter(currentEffPur->GetMaximumBin());
    if(ubooneMaxEffPurBin) *ubooneMaxEffPurBin = ubooneEffPur->GetBinCenter(ubooneEffPur->GetMaximumBin());
    if(nueMaxEffPurBin) *nueMaxEffPurBin = nueEffPur->GetBinCenter(nueEffPur->GetMaximumBin());
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
        // Choose slice based on CRUMBS score
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

double sliceCompletenessCalculator(std::vector<recoParticle> recoParticles, double sliceID){
    int counter = 0;
    double totalNumMatchedHits = 0;
    double totalNumTrueHits = 0;

    for(size_t j = 0; j < recoParticles.size(); ++j){
        if(recoParticles[j].sliceID == sliceID){
            counter++;
            totalNumMatchedHits += recoParticles[j].numMatchedHits;
            totalNumTrueHits += recoParticles[j].numTrueHits;
        }
    }

    double sliceCompleteness;
    if(counter < 1){
        sliceCompleteness = -999999;
    } else {
        sliceCompleteness = (totalNumMatchedHits / totalNumTrueHits);
    }
    //std::cout << "Chosen Slice Completeness = " << sliceCompleteness << std::endl;
    return sliceCompleteness;
}

double slicePurityCalculator(std::vector<recoParticle> recoParticles, double sliceID){
    int counter = 0;
    double totalNumMatchedHits = 0;
    double totalNumHits = 0;

    for(size_t j = 0; j < recoParticles.size(); ++j){
        if(recoParticles[j].sliceID == sliceID){
            counter++;
            totalNumMatchedHits += recoParticles[j].numMatchedHits;
            totalNumHits += recoParticles[j].numHits;
        }
    }

    double slicePurity;
    if(counter < 1){
        slicePurity = -999999;
    } else {
        slicePurity = (totalNumMatchedHits / totalNumHits);
    }
    //std::cout << "Chosen Slice Purity = " << slicePurity << std::endl;
    return slicePurity;
}

void obtainTrueParticles(std::vector<trueParticle> trueParticles, trueNeutrino chosenTrueNeutrino){
    //std::cout << "Number of true particles = " << trueParticles.size() << std::endl;
    trueInfo truth;
    double vx = chosenTrueNeutrino.vx;
    double vy = chosenTrueNeutrino.vy;
    double vz = chosenTrueNeutrino.vz;
    
    for(size_t j = 0; j < trueParticles.size(); ++j){
        
        double vertex = 0;
        if((std::abs(vx - trueParticles[j].vx) < 1e-3) + (std::abs(vy - trueParticles[j].vy) < 1e-3) + (std::abs(vz - trueParticles[j].vz) < 1e-3)) vertex = 1;
        
        if(trueParticles[j].pdg == 11 && trueParticles[j].neutrinoParent == 1 && vertex == 1) truth.numElectrons++;
        if(trueParticles[j].pdg == 22 && trueParticles[j].neutrinoParent == 1 && vertex == 1) truth.numPhotons++;
        if(trueParticles[j].pdg == 13 && trueParticles[j].neutrinoParent == 1 && vertex == 1) truth.numMuons++;
        if(trueParticles[j].pdg == 111 && trueParticles[j].neutrinoParent == 1 && vertex == 1) truth.numPiZero++;
        if(trueParticles[j].pdg == 211 || trueParticles[j].pdg == -211 && trueParticles[j].neutrinoParent == 1 && vertex == 1) truth.numChargedPi++;
        if(trueParticles[j].pdg == 2212&& trueParticles[j].neutrinoParent == 1 && vertex == 1) truth.numProtons++;
    }

    //printf("Number of Electrons = %i, Number of Photons = %i, Number of Muons = %i, Number of PiZero = %i, Number of Charged Pi = %i, Number of Protons = %i\n", truth.numElectrons, truth.numPhotons, truth.numMuons, truth.numPiZero, truth.numChargedPi, truth.numProtons); 
}

void nuEBackgroundSignalWeightsStackNuE_macro(){
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_23July.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_new.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_16Sep.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_22Sep.root");

    // v10_06_00 40k Nu+E elastic scattering with cosmics (SCE ON), BNB+Cosmic, Intime Cosmics
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_BNBIntimeNu+E+Cosmics_29Oct.root");
    
    std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsNuE/";

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
    
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNue;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBNue;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsCosmicsNue;

    double totalPOTSignalCurrent = 0;
    double totalPOTBNBCurrent = 0;
    double totalPOTCosmicsCurrent = 0;
    
    double totalPOTSignalUboone = 0;
    double totalPOTBNBUboone = 0;
    double totalPOTCosmicsUboone = 0;
    
    double totalPOTSignalNue = 0;
    double totalPOTBNBNue = 0;
    double totalPOTCosmicsNue = 0;
    
    double totalPOTSignalCurrentALL = 0;
    double totalPOTBNBCurrentALL = 0;
    double totalPOTCosmicsCurrentALL = 0;
    
    double totalPOTSignalUbooneALL = 0;
    double totalPOTBNBUbooneALL = 0;
    double totalPOTCosmicsUbooneALL = 0;
    
    double totalPOTSignalNueALL = 0;
    double totalPOTBNBNueALL = 0;
    double totalPOTCosmicsNueALL = 0;

    double cosmicSpillsSumCurrent = 0;
    double cosmicSpillsSumUboone = 0;
    double cosmicSpillsSumNue = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 2) cosmicSpillsSumCurrent += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 0) cosmicSpillsSumUboone += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNue += subRunNumGenEvents;

        std::pair<unsigned int, unsigned int> key = std::make_pair(subRunRun, subRunNumber);

        if(subRunSignal == 1){
            if(subRunDLCurrent == 2 && seenSubRunsSignalCurrent.find(key) == seenSubRunsSignalCurrent.end()){
                totalPOTSignalCurrent += subRunPOT;
                seenSubRunsSignalCurrent.insert(key);
            } else if(subRunDLCurrent == 2){
                totalPOTSignalCurrentALL += subRunPOT;
            } else if(subRunDLCurrent == 0 && seenSubRunsSignalUboone.find(key) == seenSubRunsSignalUboone.end()){
                totalPOTSignalUboone += subRunPOT;
                seenSubRunsSignalUboone.insert(key);
            } else if(subRunDLCurrent == 0){
                totalPOTSignalUbooneALL += subRunPOT;
            } else if(subRunDLCurrent == 5 && seenSubRunsSignalNue.find(key) == seenSubRunsSignalNue.end()){
                totalPOTSignalNue += subRunPOT;
                seenSubRunsSignalNue.insert(key);
            } else if(subRunDLCurrent == 5){
                totalPOTSignalNueALL += subRunPOT;
            }
        } else if(subRunSignal == 2){
            if(subRunDLCurrent == 2 && seenSubRunsBNBCurrent.find(key) == seenSubRunsBNBCurrent.end()){
                totalPOTBNBCurrent += subRunPOT;
                seenSubRunsBNBCurrent.insert(key);

            } else if(subRunDLCurrent == 2){
                totalPOTBNBCurrentALL += subRunPOT;
 
           } else if(subRunDLCurrent == 0 && seenSubRunsBNBUboone.find(key) == seenSubRunsBNBUboone.end()){
                totalPOTBNBUboone += subRunPOT;
                seenSubRunsBNBUboone.insert(key);

            } else if(subRunDLCurrent == 0){
                totalPOTBNBUbooneALL += subRunPOT;

            } else if(subRunDLCurrent == 5 && seenSubRunsBNBNue.find(key) == seenSubRunsBNBNue.end()){
                totalPOTBNBNue += subRunPOT;
                seenSubRunsBNBNue.insert(key);

            } else if(subRunDLCurrent == 5){
                totalPOTBNBNueALL += subRunPOT;
            }
        } else if(subRunSignal == 3){
            if(subRunDLCurrent == 2 && seenSubRunsCosmicsCurrent.find(key) == seenSubRunsCosmicsCurrent.end()){
                totalPOTCosmicsCurrent += subRunPOT;
                seenSubRunsCosmicsCurrent.insert(key);
            } else if(subRunDLCurrent == 2){
                totalPOTCosmicsCurrentALL += subRunPOT;
            } else if(subRunDLCurrent == 0 && seenSubRunsCosmicsUboone.find(key) == seenSubRunsCosmicsUboone.end()){
                totalPOTCosmicsUboone += subRunPOT;
                seenSubRunsCosmicsUboone.insert(key);
            } else if(subRunDLCurrent == 0){
                totalPOTCosmicsUbooneALL += subRunPOT;
            } else if(subRunDLCurrent == 5 && seenSubRunsCosmicsNue.find(key) == seenSubRunsCosmicsNue.end()){
                totalPOTCosmicsNue += subRunPOT;
                seenSubRunsCosmicsNue.insert(key);
            } else if(subRunDLCurrent == 5){
                totalPOTCosmicsNueALL += subRunPOT;
            }
        }
    }

    //double numberFiles = 600;
    //double cosmicSpills = numberFiles * 500;
    //std::cout << "cosmicSpillsSum = " << cosmicSpillsSum << std::endl;
    double cosmicsPOTCurrent = cosmicSpillsSumCurrent * 5e12;
    double cosmicsPOTUboone = cosmicSpillsSumUboone * 5e12;
    double cosmicsPOTNue = cosmicSpillsSumNue * 5e12;

    double signalWeightCurrent = totalPOTSignalCurrent/totalPOTSignalCurrent;
    double BNBWeightCurrent = totalPOTSignalCurrent/totalPOTBNBCurrent;
    double cosmicsWeightCurrent = totalPOTSignalCurrent/cosmicsPOTCurrent;
    
    double signalWeightUboone = totalPOTSignalUboone/totalPOTSignalUboone;
    double BNBWeightUboone = totalPOTSignalUboone/totalPOTBNBUboone;
    double cosmicsWeightUboone = totalPOTSignalUboone/cosmicsPOTUboone;
    
    double signalWeightNue = totalPOTSignalNue/totalPOTSignalNue;
    double BNBWeightNue = totalPOTSignalNue/totalPOTBNBNue;
    double cosmicsWeightNue = totalPOTSignalNue/cosmicsPOTNue;

    std::cout << "BNB POT Current = " << totalPOTBNBCurrent << ", Uboone = " << totalPOTBNBUboone << ", Nu+E = " << totalPOTBNBNue << std::endl; 
    std::cout << "ALL BNB POT Current = " << totalPOTBNBCurrentALL << ", Uboone = " << totalPOTBNBUbooneALL << ", Nu+E = " << totalPOTBNBNueALL << std::endl; 
    std::cout << "" << std::endl;
    std::cout << "Signal POT Current = " << totalPOTSignalCurrent << ", Uboone = " << totalPOTSignalUboone << ", Nu+E = " << totalPOTSignalNue << std::endl;
    std::cout << "ALL Signal POT Current = " << totalPOTSignalCurrentALL << ", Uboone = " << totalPOTSignalUbooneALL << ", Nu+E = " << totalPOTSignalNueALL << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Cosmics POT Current = " << totalPOTCosmicsCurrent << ", Uboone = " << totalPOTCosmicsUboone << ", Nu+E = " << totalPOTCosmicsNue << std::endl;    
    std::cout << "ALL Cosmics POT Current = " << totalPOTCosmicsCurrentALL << ", Uboone = " << totalPOTCosmicsUbooneALL << ", Nu+E = " << totalPOTCosmicsNueALL << std::endl;    
    std::cout << "Spills Cosmics POT Current = " << cosmicsPOTCurrent << ", Uboone = " << cosmicsPOTUboone << ", Nu+E = " << cosmicsPOTNue << std::endl;    

    printf("Weights:\nCurrent: Signal = %f, BNB = %f, Cosmics = %f\nUboone: Signal = %f, BNB = %f, Cosmics = %f\nNu+E: Signal = %f, BNB = %f, Cosmics = %f\n", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue);

    UInt_t eventID, runID, subRunID;
    double DLCurrent, signal;
    
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);
    tree->SetBranchAddress("signal", &signal);

    std::vector<double> *truth_neutrinoVX = nullptr;
    std::vector<double> *truth_neutrinoVY = nullptr;
    std::vector<double> *truth_neutrinoVZ = nullptr;
    std::vector<double> *truth_CCNC = nullptr;
    std::vector<double> *truth_neutrinoType = nullptr;
    std::vector<double> *truth_chargedLepton = nullptr;
    std::vector<double> *truth_neutrinoTPCID = nullptr;
    std::vector<double> *truth_neutrinoTPCValid = nullptr;
    std::vector<double> *truth_genieMode = nullptr;
    std::vector<double> *truth_genieInteraction = nullptr;

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
    std::vector<double> *reco_particleNumHits = nullptr;
    std::vector<double> *reco_particleNumMatchedHits = nullptr;
    std::vector<double> *reco_particleNumTrueHits = nullptr;
    std::vector<double> *reco_particleTruePDG = nullptr;

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
    tree->SetBranchAddress("truth_genieMode", &truth_genieMode);
    tree->SetBranchAddress("truth_genieInteraction", &truth_genieInteraction);

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
    tree->SetBranchAddress("reco_particleNumHits", &reco_particleNumHits);
    tree->SetBranchAddress("reco_particleNumMatchedHits", &reco_particleNumMatchedHits);
    tree->SetBranchAddress("reco_particleNumTrueHits", &reco_particleNumTrueHits);
    tree->SetBranchAddress("reco_particleTruePDG", &reco_particleTruePDG);

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
    //auto numSlicesCompleteness = createHistGroup("numSlicesCompleteness", "Number of Slices with Completeness > 0 in an Event", "Number of Slices with Completeness > 0", 10, 0, 10);
    auto numRecoNeutrinos = createHistGroup("numRecoNeutrinos", "Number of Reco Neutrinos in an Event", "Number of Reco Neutrinos", 10, 0, 10);

    // Highest CRUMBS Score Slice Plots
    auto sliceCompletenessCRUMBS = createHistGroup("sliceCompletenessCRUMBS", "Completeness of the Slice with the Highest CRUMBS Score", "Completeness", 102, 0, 1.02); // Bin width = 0.005
    auto sliceScoreCRUMBS = createHistGroup("sliceScoreCRUMBS", "CRUMBS Score of the Slice with the Highest CRUMBS Score", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCRUMBS = createHistGroup("slicePurityCRUMBS", "Purity of the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto highestPFPCompletenessCRUMBS = createHistGroup("highestPFPCompletenessCRUMBS", "Completeness of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "Completeness", 50, 0, 1.02);
    auto highestPFPPurityCRUMBS = createHistGroup("highestPFPPurityCRUMBS", "Purity of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto deltaXCRUMBS = createHistGroup("deltaXCRUMBS", "#Deltax Distribution: Slice with Highest CRUMBS Score", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaYCRUMBS = createHistGroup("deltaYCRUMBS", "#Deltay Distribution: Slice with Highest CRUMBS Score", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZCRUMBS = createHistGroup("deltaZCRUMBS", "#Deltaz Distribution: Slice with Highest CRUMBS Score", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaRCRUMBS = createHistGroup("deltaRCRUMBS", "#Delta#bar{r} Distribution: Slice with Highest CRUMBS Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    auto numPFPsCRUMBS = createHistGroup("numPFPsCRUMBS", "Number of PFPs in the Slice with the Highest CRUMBS Score", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCRUMBS = createHistGroup("ratioChosenSummedEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCRUMBS = createHistGroup("ratioChosenTrueEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCRUMBS = createHistGroup("ratioSummedTrueEnergyCRUMBS", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    //auto angleDifferenceCRUMBS = createHistGroup("angleDifferenceCRUMBS", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto EtrueThetaRecoCRUMBS = createHistGroup("EtrueThetaRecoCRUMBS", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCRUMBS = createHistGroup("ERecoSumThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBS = createHistGroup("ERecoHighestThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCRUMBS = createHistGroup("ERecoSumThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCRUMBS = createHistGroup("ERecoHighestThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);

    auto primaryVertexXCRUMBS = createHistGroup("primaryVertexXCRUMBS", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "X Coordinate (cm)", 82, -201.3, 201.3);
    auto primaryVertexYCRUMBS = createHistGroup("primaryVertexYCRUMBS", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Y Coordinate (cm)", 83, -203.8, 203.8);
    auto primaryVertexZCRUMBS = createHistGroup("primaryVertexZCRUMBS", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Z Coordinate (cm)", 104, 0, 509.4);

    auto primaryVertexXCRUMBSNegative = createHistGroup("primaryVertexXCRUMBSNegative", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "X Coordinate (cm)", 805, -201.3, 201.3);
    auto primaryVertexYCRUMBSNegative = createHistGroup("primaryVertexYCRUMBSNegative", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Y Coordinate (cm)", 815, -203.8, 203.8);
    auto primaryVertexZCRUMBSNegative = createHistGroup("primaryVertexZCRUMBSNegative", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Z Coordinate (cm)", 1019, 0, 509.4);

    auto primaryVertexXCRUMBSPositive = createHistGroup("primaryVertexXCRUMBSPositive", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "X Coordinate (cm)", 805, -201.3, 201.3);
    auto primaryVertexYCRUMBSPositive = createHistGroup("primaryVertexYCRUMBSPositive", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Y Coordinate (cm)", 815, -203.8, 203.8);
    auto primaryVertexZCRUMBSPositive = createHistGroup("primaryVertexZCRUMBSPositive", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Z Coordinate (cm)", 1019, 0, 509.4);
   
    auto allRecoVertexX = createHistGroup("allRecoVertexX", "X Coordinate of the Reco Neutrino Vertex (all reco vertices)", "X Coordinate (cm)", 82, -201.3, 201.3);
    auto allRecoVertexY = createHistGroup("allRecoVertexY", "Y Coordinate of the Reco Neutrino Vertex (all reco vertices)", "Y Coordinate (cm)", 83, -203.8, 203.8);
    auto allRecoVertexZ = createHistGroup("allRecoVertexZ", "Z Coordinate of the Reco Neutrino Vertex (all reco vertices)", "Z Coordinate (cm)", 104, 0, 509.4);

    // Highest Completeness Score Slice Plots
    //auto sliceCompletenessCompleteness = createHistGroup("sliceCompletenessCompleteness", "Completeness of the Slice with the Highest Completeness", "Completeness Score", 102, 0, 1.02); // Bin width = 0.01
    //auto sliceScoreCompleteness = createHistGroup("sliceScoreCompleteness", "CRUMBS Score of the Slice with the Highest Completeness", "CRUMBS Score", 25, -1, 1); 
    //auto slicePurityCompleteness = createHistGroup("slicePurityCompleteness", "Purity of the Slice with the Highest Completeness", "Purity", 50, 0, 1.02);
    //auto deltaXCompleteness = createHistGroup("deltaXCompleteness", "#Deltax Distribution: Slice with Highest Completeness", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    //auto deltaYCompleteness = createHistGroup("deltaYCompleteness", "#Deltay Distribution: Slice with Highest Completeness", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    //auto deltaZCompleteness = createHistGroup("deltaZCompleteness", "#Deltaz Distribution: Slice with Highest Completeness", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    //auto deltaRCompleteness = createHistGroup("deltaRCompleteness", "#Delta#bar{r} Distribution: Slice with Highest Completeness Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    //auto numPFPsCompleteness = createHistGroup("numPFPsCompleteness", "Number of PFPs in the Slice with the Highest Completeness", "Number of PFPs", 10, 0, 10);
    //auto ratioChosenSummedEnergyCompleteness = createHistGroup("ratioChosenSummedEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest Completeness", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    //auto ratioChosenTrueEnergyCompleteness = createHistGroup("ratioChosenTrueEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    //auto ratioSummedTrueEnergyCompleteness = createHistGroup("ratioSummedTrueEnergyCompleteness", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    //auto angleDifferenceCompleteness = createHistGroup("angleDifferenceCompleteness", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest Completeness", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    //auto EtrueThetaRecoCompleteness = createHistGroup("EtrueThetaRecoCompleteness", "E_{true}#theta_{reco}^{2}: Slice with Highest Completeness", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    //auto ERecoSumThetaTrueCompleteness = createHistGroup("ERecoSumThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    //auto ERecoHighestThetaTrueCompleteness = createHistGroup("ERecoHighestThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    //auto ERecoSumThetaRecoCompleteness = createHistGroup("ERecoSumThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    //auto ERecoHighestThetaRecoCompleteness = createHistGroup("ERecoHighestThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    
    counter numEventsTotal;
    counter numEventsTrueNeutrino;
    counter numEventsSlices;
    counter numEventsRecoNeutrino;
    counter numEventsSliceNotEqualNeutrino;
    //counter numSlicesCRUMBSCompletenessZero;
    counter numEventsCRUMBSRecoParticle;
    //counter sameSliceSelected;
    interactionCounter interactionsSignal;
    interactionCounter interactionsBNB;

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);
        //if(signal != 3 || DLCurrent != 0) continue;
        //printf("___________________________________________________________________\n");
        //std::cout << "Signal = " << signal << ", DLCurrent = " << DLCurrent << std::endl;
        std::vector<recoParticle> recoParticlesInEvent;
        std::vector<recoNeutrino> recoNeutrinosInEvent;
        std::vector<recoSlice> recoSlicesInEvent;
        std::vector<trueParticle> trueParticlesInEvent;
        
        recoParticle chosenRecoParticleCRUMBS;
        //recoParticle chosenRecoParticleCompleteness;
        recoNeutrino chosenRecoNeutrinoCRUMBS;
        //recoNeutrino chosenRecoNeutrinoCompleteness;
        //recoSlice chosenRecoSliceCompleteness;
        recoSlice chosenRecoSliceCRUMBS;
        //trueParticle chosenTrueParticle;
        trueNeutrino chosenTrueNeutrino;
        double energyInSlice = 0;
        double numPFPsInEvent = 0;
        double numPFPsInChosenSlice = 0;

        if(DLCurrent == 0 && signal == 1){
            numEventsTotal.ubooneSignal++;
        } else if(DLCurrent == 2 && signal == 1){
            numEventsTotal.currentSignal++;
        } else if(DLCurrent == 5 && signal == 1){
            numEventsTotal.nueSignal++;
        } else if(DLCurrent == 0 && signal == 2){
            numEventsTotal.ubooneBNB++;
        } else if(DLCurrent == 2 && signal == 2){
            numEventsTotal.currentBNB++;
        } else if(DLCurrent == 5 && signal == 2){
            numEventsTotal.nueBNB++;
        } else if(DLCurrent == 0 && signal == 3){
            numEventsTotal.ubooneCosmics++;
        } else if(DLCurrent == 2 && signal == 3){
            numEventsTotal.currentCosmics++;
        } else if(DLCurrent == 5 && signal == 3){
            numEventsTotal.nueCosmics++;
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
                chosenTrueNeutrino.genieMode = truth_genieMode->at(j);
                chosenTrueNeutrino.genieInteraction = truth_genieInteraction->at(j);

                //std::cout << "genieMode: " << chosenTrueNeutrino.genieMode << std::endl;
                //std::cout << "genieInteraction: " << chosenTrueNeutrino.genieInteraction << std::endl;
                if(DLCurrent == 0){
                    if(chosenTrueNeutrino.genieMode == -1 && signal == 1) interactionsSignal.Unknown++;
                    if(chosenTrueNeutrino.genieMode == 0 && signal == 1) interactionsSignal.QE++;
                    if(chosenTrueNeutrino.genieMode == 1 && signal == 1) interactionsSignal.Res++;
                    if(chosenTrueNeutrino.genieMode == 2 && signal == 1) interactionsSignal.DIS++;
                    if(chosenTrueNeutrino.genieMode == 3 && signal == 1) interactionsSignal.Coh++;
                    if(chosenTrueNeutrino.genieMode == 4 && signal == 1) interactionsSignal.CohElastic++;
                    if(chosenTrueNeutrino.genieMode == 5 && signal == 1) interactionsSignal.ElectronScattering++;
                    if(chosenTrueNeutrino.genieMode == 6 && signal == 1) interactionsSignal.IMDAnnihilation++;
                    if(chosenTrueNeutrino.genieMode == 7 && signal == 1) interactionsSignal.InverseBetaDecay++;
                    if(chosenTrueNeutrino.genieMode == 8 && signal == 1) interactionsSignal.GlashowResonance++;
                    if(chosenTrueNeutrino.genieMode == 9 && signal == 1) interactionsSignal.AMNuGamma++;
                    if(chosenTrueNeutrino.genieMode == 10 && signal == 1) interactionsSignal.MEC++;
                    if(chosenTrueNeutrino.genieMode == 11 && signal == 1) interactionsSignal.Diffractive++;
                    if(chosenTrueNeutrino.genieMode == 12 && signal == 1) interactionsSignal.EM++;
                    if(chosenTrueNeutrino.genieMode == 13 && signal == 1) interactionsSignal.WeakMix++;
                    if(chosenTrueNeutrino.genieInteraction == 1000 && signal == 1) interactionsSignal.NuanceOffset++;
                    if(chosenTrueNeutrino.genieInteraction == 1001 && signal == 1) interactionsSignal.CCQE++;
                    if(chosenTrueNeutrino.genieInteraction == 1002 && signal == 1) interactionsSignal.NCQE++;
                    if(chosenTrueNeutrino.genieInteraction > 1002 && chosenTrueNeutrino.genieInteraction < 1091 && signal == 1) interactionsSignal.NuanceRes++;
                    if(chosenTrueNeutrino.genieInteraction == 1091 && signal == 1) interactionsSignal.CCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1092 && signal == 1) interactionsSignal.NCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1098 && signal == 1) interactionsSignal.NuEElastic++;

                    if(chosenTrueNeutrino.genieMode == -1 && signal == 2) interactionsBNB.Unknown++;
                    if(chosenTrueNeutrino.genieMode == 0 && signal == 2) interactionsBNB.QE++;
                    if(chosenTrueNeutrino.genieMode == 1 && signal == 2) interactionsBNB.Res++;
                    if(chosenTrueNeutrino.genieMode == 2 && signal == 2) interactionsBNB.DIS++;
                    if(chosenTrueNeutrino.genieMode == 3 && signal == 2) interactionsBNB.Coh++;
                    if(chosenTrueNeutrino.genieMode == 4 && signal == 2) interactionsBNB.CohElastic++;
                    if(chosenTrueNeutrino.genieMode == 5 && signal == 2) interactionsBNB.ElectronScattering++;
                    if(chosenTrueNeutrino.genieMode == 6 && signal == 2) interactionsBNB.IMDAnnihilation++;
                    if(chosenTrueNeutrino.genieMode == 7 && signal == 2) interactionsBNB.InverseBetaDecay++;
                    if(chosenTrueNeutrino.genieMode == 8 && signal == 2) interactionsBNB.GlashowResonance++;
                    if(chosenTrueNeutrino.genieMode == 9 && signal == 2) interactionsBNB.AMNuGamma++;
                    if(chosenTrueNeutrino.genieMode == 10 && signal == 2) interactionsBNB.MEC++;
                    if(chosenTrueNeutrino.genieMode == 11 && signal == 2) interactionsBNB.Diffractive++;
                    if(chosenTrueNeutrino.genieMode == 12 && signal == 2) interactionsBNB.EM++;
                    if(chosenTrueNeutrino.genieMode == 13 && signal == 2) interactionsBNB.WeakMix++;
                    if(chosenTrueNeutrino.genieInteraction == 1000 && signal == 2) interactionsBNB.NuanceOffset++;
                    if(chosenTrueNeutrino.genieInteraction == 1001 && signal == 2) interactionsBNB.CCQE++;
                    if(chosenTrueNeutrino.genieInteraction == 1002 && signal == 2) interactionsBNB.NCQE++;
                    if(chosenTrueNeutrino.genieInteraction > 1002 && chosenTrueNeutrino.genieInteraction < 1091 && signal == 2) interactionsBNB.NuanceRes++;
                    if(chosenTrueNeutrino.genieInteraction == 1091 && signal == 2) interactionsBNB.CCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1092 && signal == 2) interactionsBNB.NCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1098 && signal == 2) interactionsBNB.NuEElastic++;
                }
            } else if(truth_neutrinoTPCValid->at(j) != -999999 && truth_neutrinoVX->at(j) < 201.3 && truth_neutrinoVX->at(j) > -201.3 && truth_neutrinoVY->at(j) < 203.732 && truth_neutrinoVY->at(j) > -203.732 && truth_neutrinoVZ->at(j) > 0 && truth_neutrinoVZ->at(j) < 509.4){
                //printf("True Neutrino Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            } else if(truth_neutrinoTPCValid->at(j) != -999999){
                //printf("True Neutrino Not Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            }
        }
        
        if(truth_neutrinoVX->size() == 1 && truth_neutrinoVX->at(0) == -999999) printf("No True Neutrino in Event\n");

        // If there are no truth neutrinos within the TPC, skip the event
        if(truthNeutrinoInTPC == 0 && signal != 3) continue;
        
        if(DLCurrent == 0 && signal == 1){
            numEventsTrueNeutrino.ubooneSignal++;
        } else if(DLCurrent == 2 && signal == 1){
            numEventsTrueNeutrino.currentSignal++;
        } else if(DLCurrent == 5 && signal == 1){
            numEventsTrueNeutrino.nueSignal++;
        } else if(DLCurrent == 0 && signal == 2){
            numEventsTrueNeutrino.ubooneBNB++;
        } else if(DLCurrent == 2 && signal == 2){
            numEventsTrueNeutrino.currentBNB++;
        } else if(DLCurrent == 5 && signal == 2){
            numEventsTrueNeutrino.nueBNB++;
        }

        if(signal != 3){
            //printf("True Neutrino Vertex = (%f, %f, %f)\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz);
        }

        //std::cout << "Number of truth particles = " << truth_particlePDG->size() << std::endl;             
    
        for(size_t k = 0; k < truth_particlePDG->size(); ++k){
            trueParticle truthParticle;
            truthParticle.pdg = truth_particlePDG->at(k);
            truthParticle.vx = truth_particleVX->at(k);
            truthParticle.vy = truth_particleVY->at(k);
            truthParticle.vz = truth_particleVZ->at(k);
            truthParticle.px = truth_particlePX->at(k);
            truthParticle.py = truth_particlePY->at(k);
            truthParticle.pz = truth_particlePZ->at(k);
            truthParticle.energy = truth_particleEnergy->at(k);
            truthParticle.angle = truth_particleAngle->at(k);
            truthParticle.ETheta2 = truth_particleETheta2->at(k);
            truthParticle.dx = truth_particleDirectionX->at(k);
            truthParticle.dy = truth_particleDirectionY->at(k);
            truthParticle.dz = truth_particleDirectionZ->at(k);
            truthParticle.neutrinoParent = truth_particleNeutrinoParent->at(k);
            trueParticlesInEvent.push_back(truthParticle);
        }

        obtainTrueParticles(trueParticlesInEvent, chosenTrueNeutrino);

        int slice = 0;
        int crumbsSlice = 0;
        //int completenessSlice = 0;
        double numSlicesInEvent = 0;
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            if(reco_sliceID->at(j) != -999999){
                slice = 1;
                numSlicesInEvent++;

                recoSlice recoSlice;
                recoSlice.id = reco_sliceID->at(j);
                //recoSlice.completeness = reco_sliceCompleteness->at(j);
                //recoSlice.purity = reco_slicePurity->at(j);
                recoSlice.completeness = -999999;
                recoSlice.purity = -999999;
                recoSlice.score = reco_sliceScore->at(j);
                recoSlicesInEvent.push_back(recoSlice);
        
                if(recoSlice.score != -999999) crumbsSlice++;  
                //if(recoSlice.completeness > 0) completenessSlice++;  
                //printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", j, recoSlice.id, recoSlice.completeness, recoSlice.purity, recoSlice.score);
            }
        }    

        // If there are no slices, skip the event
        if(slice == 0) continue;  
     
        // Pick the slice with the highest completeness
        //chosenRecoSliceCompleteness = chooseSlice(recoSlicesInEvent, 0);
        //printf("\nSlice with Highest Completeness: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", chosenRecoSliceCompleteness.id, chosenRecoSliceCompleteness.completeness, chosenRecoSliceCompleteness.purity, chosenRecoSliceCompleteness.score);

        // Pick the slice with the highest score
        chosenRecoSliceCRUMBS = chooseSlice(recoSlicesInEvent, 1);
        //printf("Slice with Highest CRUMBS Score: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n\n", chosenRecoSliceCRUMBS.id, chosenRecoSliceCRUMBS.completeness, chosenRecoSliceCRUMBS.purity, chosenRecoSliceCRUMBS.score);
        
        if(DLCurrent == 0 && signal == 1){
            numEventsSlices.ubooneSignal++;
            numSlices.ubooneSignal->Fill(numSlicesInEvent);
            numSlicesCRUMBS.ubooneSignal->Fill(crumbsSlice);
            sliceScoreCRUMBS.ubooneSignal->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 2 && signal == 1){
            numEventsSlices.currentSignal++;
            numSlices.currentSignal->Fill(numSlicesInEvent);
            numSlicesCRUMBS.currentSignal->Fill(crumbsSlice);
        } else if(DLCurrent == 5 && signal == 1){
            numEventsSlices.nueSignal++;
            numSlices.nueSignal->Fill(numSlicesInEvent);
            numSlicesCRUMBS.nueSignal->Fill(crumbsSlice);
            sliceScoreCRUMBS.nueSignal->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCRUMBS.nueSignal->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 0 && signal == 2){
            numEventsSlices.ubooneBNB++;
            numSlices.ubooneBNB->Fill(numSlicesInEvent);
            numSlicesCRUMBS.ubooneBNB->Fill(crumbsSlice);
            sliceScoreCRUMBS.ubooneBNB->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 2 && signal == 2){
            numEventsSlices.currentBNB++;
            numSlices.currentBNB->Fill(numSlicesInEvent);
            numSlicesCRUMBS.currentBNB->Fill(crumbsSlice);
            sliceScoreCRUMBS.currentBNB->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 5 && signal == 2){
            numEventsSlices.nueBNB++;
            numSlices.nueBNB->Fill(numSlicesInEvent);
            numSlicesCRUMBS.nueBNB->Fill(crumbsSlice);
            sliceScoreCRUMBS.nueBNB->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 0 && signal == 3){
            numEventsSlices.ubooneCosmics++;
            numSlices.ubooneCosmics->Fill(numSlicesInEvent);
            numSlicesCRUMBS.ubooneCosmics->Fill(crumbsSlice);
            sliceScoreCRUMBS.ubooneCosmics->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 2 && signal == 3){
            numEventsSlices.currentCosmics++;
            numSlices.currentCosmics->Fill(numSlicesInEvent);
            numSlicesCRUMBS.currentCosmics->Fill(crumbsSlice);
            sliceScoreCRUMBS.currentCosmics->Fill(chosenRecoSliceCRUMBS.score);
        } else if(DLCurrent == 2 && signal == 3){
            numEventsSlices.nueCosmics++;
            numSlices.nueCosmics->Fill(numSlicesInEvent);
            numSlicesCRUMBS.nueCosmics->Fill(crumbsSlice);
            sliceScoreCRUMBS.nueCosmics->Fill(chosenRecoSliceCRUMBS.score);
        }
        
        int neutrino = 0;
        int numRecoNeutrinosInEvent = 0;
        //std::cout << "num reco neutrinos: " << reco_neutrinoPDG->size() << std::endl;
        for(size_t j = 0; j < reco_neutrinoPDG->size(); ++j){
            if(reco_neutrinoPDG->at(j) != -999999){
                neutrino = 1;
                numRecoNeutrinosInEvent++;
                
                recoNeutrino recoNeutrino;
                recoNeutrino.pdg = reco_neutrinoPDG->at(j);
                recoNeutrino.isPrimary = reco_neutrinoIsPrimary->at(j);
                recoNeutrino.vx = reco_neutrinoVX->at(j);
                recoNeutrino.vy = reco_neutrinoVY->at(j);
                recoNeutrino.vz = reco_neutrinoVZ->at(j);
                recoNeutrino.sliceID = reco_neutrinoSliceID->at(j);
                recoNeutrinosInEvent.push_back(recoNeutrino);

                if(DLCurrent == 0 && signal == 1){
                    allRecoVertexX.ubooneSignal->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.ubooneSignal->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.ubooneSignal->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 2 && signal == 1){
                    allRecoVertexX.currentSignal->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.currentSignal->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.currentSignal->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 5 && signal == 1){
                    allRecoVertexX.nueSignal->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.nueSignal->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.nueSignal->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 0 && signal == 2){
                    allRecoVertexX.ubooneBNB->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.ubooneBNB->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.ubooneBNB->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 2 && signal == 2){
                    allRecoVertexX.currentBNB->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.currentBNB->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.currentBNB->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 5 && signal == 2){
                    allRecoVertexX.nueBNB->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.nueBNB->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.nueBNB->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 0 && signal == 3){
                    allRecoVertexX.ubooneCosmics->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.ubooneCosmics->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.ubooneCosmics->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 2 && signal == 3){
                    allRecoVertexX.currentCosmics->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.currentCosmics->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.currentCosmics->Fill(reco_neutrinoVZ->at(j));
                }
                if(DLCurrent == 5 && signal == 3){
                    allRecoVertexX.nueCosmics->Fill(reco_neutrinoVX->at(j));
                    allRecoVertexY.nueCosmics->Fill(reco_neutrinoVY->at(j));
                    allRecoVertexZ.nueCosmics->Fill(reco_neutrinoVZ->at(j));
                }

                //printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, recoNeutrino.pdg, recoNeutrino.isPrimary, recoNeutrino.vx, recoNeutrino.vy, recoNeutrino.vz, recoNeutrino.sliceID);
            }
        }
        
        if(DLCurrent == 0 && signal == 1) numRecoNeutrinos.ubooneSignal->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 2 && signal == 1) numRecoNeutrinos.currentSignal->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 5 && signal == 1) numRecoNeutrinos.nueSignal->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 0 && signal == 2) numRecoNeutrinos.ubooneBNB->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 2 && signal == 2) numRecoNeutrinos.currentBNB->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 5 && signal == 2) numRecoNeutrinos.nueBNB->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 0 && signal == 3) numRecoNeutrinos.ubooneCosmics->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 2 && signal == 3) numRecoNeutrinos.currentCosmics->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 5 && signal == 3) numRecoNeutrinos.nueCosmics->Fill(numRecoNeutrinosInEvent);

        // Skip the event if there are no reconstructed neutrinos
        if(neutrino == 0) continue;
    
        // The number of slices with a CRUMBS score != the number of reconstructed neutrinos in the event 
        if(crumbsSlice != numRecoNeutrinosInEvent){
            if(DLCurrent == 0 && signal == 1) numEventsSliceNotEqualNeutrino.ubooneSignal++; 
            if(DLCurrent == 2 && signal == 1) numEventsSliceNotEqualNeutrino.currentSignal++;
            if(DLCurrent == 5 && signal == 1) numEventsSliceNotEqualNeutrino.nueSignal++;
            if(DLCurrent == 0 && signal == 2) numEventsSliceNotEqualNeutrino.ubooneBNB++; 
            if(DLCurrent == 2 && signal == 2) numEventsSliceNotEqualNeutrino.currentBNB++;
            if(DLCurrent == 5 && signal == 2) numEventsSliceNotEqualNeutrino.nueBNB++;
            if(DLCurrent == 0 && signal == 3) numEventsSliceNotEqualNeutrino.ubooneCosmics++; 
            if(DLCurrent == 2 && signal == 3) numEventsSliceNotEqualNeutrino.currentCosmics++;
            if(DLCurrent == 5 && signal == 3) numEventsSliceNotEqualNeutrino.nueCosmics++;
        } 

        // Look at reconstruction when we pick the slice with highest CRUMBS score
        chosenRecoNeutrinoCRUMBS = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCRUMBS.id); 

        if(DLCurrent == 0 && signal == 3) numEventsRecoNeutrino.ubooneCosmics++;
        if(DLCurrent == 2 && signal == 3) numEventsRecoNeutrino.currentCosmics++;
        if(DLCurrent == 5 && signal == 3) numEventsRecoNeutrino.nueCosmics++;

        if(signal != 3){
            double deltaXCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vx - chosenTrueNeutrino.vx);
            double deltaYCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vy - chosenTrueNeutrino.vy);
            double deltaZCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vz - chosenTrueNeutrino.vz);
            double deltaRCRUMBSValue = std::sqrt((deltaXCRUMBSValue * deltaXCRUMBSValue) + (deltaYCRUMBSValue * deltaYCRUMBSValue) + (deltaZCRUMBSValue * deltaZCRUMBSValue));

            if(DLCurrent == 0 && signal == 1){
                numEventsRecoNeutrino.ubooneSignal++;
                deltaXCRUMBS.ubooneSignal->Fill(deltaXCRUMBSValue);
                deltaYCRUMBS.ubooneSignal->Fill(deltaYCRUMBSValue);
                deltaZCRUMBS.ubooneSignal->Fill(deltaZCRUMBSValue);
                deltaRCRUMBS.ubooneSignal->Fill(deltaRCRUMBSValue);
            } else if(DLCurrent == 2 && signal == 1){
                numEventsRecoNeutrino.currentSignal++;
                deltaXCRUMBS.currentSignal->Fill(deltaXCRUMBSValue);
                deltaYCRUMBS.currentSignal->Fill(deltaYCRUMBSValue);
                deltaZCRUMBS.currentSignal->Fill(deltaZCRUMBSValue);
                deltaRCRUMBS.currentSignal->Fill(deltaRCRUMBSValue);
            } else if(DLCurrent == 5 && signal == 1){
                numEventsRecoNeutrino.nueSignal++;
                deltaXCRUMBS.nueSignal->Fill(deltaXCRUMBSValue);
                deltaYCRUMBS.nueSignal->Fill(deltaYCRUMBSValue);
                deltaZCRUMBS.nueSignal->Fill(deltaZCRUMBSValue);
                deltaRCRUMBS.nueSignal->Fill(deltaRCRUMBSValue);
            } else if(DLCurrent == 0 && signal == 2){
                numEventsRecoNeutrino.ubooneBNB++;
                deltaXCRUMBS.ubooneBNB->Fill(deltaXCRUMBSValue);
                deltaYCRUMBS.ubooneBNB->Fill(deltaYCRUMBSValue);
                deltaZCRUMBS.ubooneBNB->Fill(deltaZCRUMBSValue);
                deltaRCRUMBS.ubooneBNB->Fill(deltaRCRUMBSValue);
            } else if(DLCurrent == 2 && signal == 2){
                numEventsRecoNeutrino.currentBNB++;
                deltaXCRUMBS.currentBNB->Fill(deltaXCRUMBSValue);
                deltaYCRUMBS.currentBNB->Fill(deltaYCRUMBSValue);
                deltaZCRUMBS.currentBNB->Fill(deltaZCRUMBSValue);
                deltaRCRUMBS.currentBNB->Fill(deltaRCRUMBSValue);
            } else if(DLCurrent == 5 && signal == 2){
                numEventsRecoNeutrino.nueBNB++;
                deltaXCRUMBS.nueBNB->Fill(deltaXCRUMBSValue);
                deltaYCRUMBS.nueBNB->Fill(deltaYCRUMBSValue);
                deltaZCRUMBS.nueBNB->Fill(deltaZCRUMBSValue);
                deltaRCRUMBS.nueBNB->Fill(deltaRCRUMBSValue);
            }
        }

        if(DLCurrent == 0 && signal == 1){
            primaryVertexXCRUMBS.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);

                primaryVertexXCRUMBSNegative.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.ubooneSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);

        } else if(DLCurrent == 2 && signal == 1){
            primaryVertexXCRUMBS.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.currentSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);

        } else if(DLCurrent == 5 && signal == 1){
            primaryVertexXCRUMBS.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.nueSignal->Fill(chosenRecoNeutrinoCRUMBS.vz);

        } else if(DLCurrent == 0 && signal == 2){
            primaryVertexXCRUMBS.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.ubooneBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);

        } else if(DLCurrent == 2 && signal == 2){
            primaryVertexXCRUMBS.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.currentBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);

        } else if(DLCurrent == 5 && signal == 2){
            primaryVertexXCRUMBS.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.nueBNB->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
        } else if(DLCurrent == 0 && signal == 3){
            primaryVertexXCRUMBS.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.ubooneCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);

        } else if(DLCurrent == 2 && signal == 3){
            primaryVertexXCRUMBS.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.currentCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
        } else if(DLCurrent == 5 && signal == 3){
            primaryVertexXCRUMBS.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
            primaryVertexYCRUMBS.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
            primaryVertexZCRUMBS.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
                primaryVertexXCRUMBSNegative.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositive.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexYCRUMBSPositive.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy);
                primaryVertexZCRUMBSNegative.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
                primaryVertexZCRUMBSPositive.nueCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz);
        }

        int recoparticle = 0;
        int recoparticleCRUMBS = 0;

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
                recoParticle.numHits = reco_particleNumHits->at(j);
                recoParticle.numMatchedHits = reco_particleNumMatchedHits->at(j);
                recoParticle.numTrueHits = reco_particleNumTrueHits->at(j);
                recoParticle.truePDG = reco_particleTruePDG->at(j);
                recoParticlesInEvent.push_back(recoParticle);
            
                if(recoParticle.sliceID == chosenRecoSliceCRUMBS.id) recoparticleCRUMBS = 1;
                //if(recoParticle.sliceID == chosenRecoSliceCompleteness.id) recoparticleCompleteness = 1;
            }
        }

        // There are no reco particles in the event
        if(recoparticle == 0) continue;
        
        // No reco particles in the slice with the highest CRUMBS score
        if(recoparticleCRUMBS == 0) continue;
        
        double totalSliceEnergyCRUMBS = 0;
        double numPFPsSliceCRUMBS = 0;
        chosenRecoParticleCRUMBS = choosePFP(recoParticlesInEvent, chosenRecoSliceCRUMBS.id, totalSliceEnergyCRUMBS, numPFPsSliceCRUMBS);
        double chosenSlicePurityCRUMBS = slicePurityCalculator(recoParticlesInEvent, chosenRecoSliceCRUMBS.id);
        double chosenSliceCompletenessCRUMBS = sliceCompletenessCalculator(recoParticlesInEvent, chosenRecoSliceCRUMBS.id);

        //std::cout << "highest energy PFP PDG: " << chosenRecoParticleCRUMBS.truePDG << std::endl;
        //std::cout << chosenRecoParticleCRUMBS.truePDG << std::endl;

        //double aDOTbCRUMBS = ((chosenRecoParticleCRUMBS.dx * chosenTrueParticle.dx) + (chosenRecoParticleCRUMBS.dy * chosenTrueParticle.dy) + (chosenRecoParticleCRUMBS.dz * chosenTrueParticle.dz));
        //double aMagCRUMBS = std::sqrt((chosenRecoParticleCRUMBS.dx * chosenRecoParticleCRUMBS.dx) + (chosenRecoParticleCRUMBS.dy * chosenRecoParticleCRUMBS.dy) + (chosenRecoParticleCRUMBS.dz * chosenRecoParticleCRUMBS.dz));
        //double bMagCRUMBS = std::sqrt((chosenTrueParticle.dx * chosenTrueParticle.dx) + (chosenTrueParticle.dy * chosenTrueParticle.dy) + (chosenTrueParticle.dz * chosenTrueParticle.dz));
        //double cosAngleCRUMBS = aDOTbCRUMBS/(aMagCRUMBS * bMagCRUMBS);
        //if(cosAngleCRUMBS > 1.0) cosAngleCRUMBS = 1.0;
        //if(cosAngleCRUMBS < -1.0) cosAngleCRUMBS = -1.0;
        //double angleDiffCRUMBS = TMath::ACos(cosAngleCRUMBS) * TMath::RadToDeg();
        
        if(DLCurrent == 0 && signal == 1){
            numEventsCRUMBSRecoParticle.ubooneSignal++;
            numPFPsCRUMBS.ubooneSignal->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.ubooneSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.ubooneSignal->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.ubooneSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.ubooneSignal->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.ubooneSignal->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.ubooneSignal->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.ubooneSignal->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 2 && signal == 1){
            numEventsCRUMBSRecoParticle.currentSignal++;
            numPFPsCRUMBS.currentSignal->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.currentSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.currentSignal->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.currentSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.currentSignal->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.currentSignal->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.currentSignal->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.currentSignal->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 5 && signal == 1){
            numEventsCRUMBSRecoParticle.nueSignal++;
            numPFPsCRUMBS.nueSignal->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.nueSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.nueSignal->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.nueSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.nueSignal->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.nueSignal->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.nueSignal->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.nueSignal->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 0 && signal == 2){
            numEventsCRUMBSRecoParticle.ubooneBNB++;
            numPFPsCRUMBS.ubooneBNB->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.ubooneBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.ubooneBNB->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.ubooneBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.ubooneBNB->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.ubooneBNB->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.ubooneBNB->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.ubooneBNB->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 2 && signal == 2){
            numEventsCRUMBSRecoParticle.currentBNB++;
            numPFPsCRUMBS.currentBNB->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.currentBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.currentBNB->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.currentBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.currentBNB->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.currentBNB->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.currentBNB->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.currentBNB->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 5 && signal == 2){
            numEventsCRUMBSRecoParticle.nueBNB++;
            numPFPsCRUMBS.nueBNB->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.nueBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.nueBNB->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.nueBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.nueBNB->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.nueBNB->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.nueBNB->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.nueBNB->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 0 && signal == 3){
            numEventsCRUMBSRecoParticle.ubooneCosmics++;
            numPFPsCRUMBS.ubooneCosmics->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.ubooneCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.ubooneCosmics->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.ubooneCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.ubooneCosmics->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.ubooneCosmics->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.ubooneCosmics->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.ubooneCosmics->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 2 && signal == 3){
            numEventsCRUMBSRecoParticle.currentCosmics++;
            numPFPsCRUMBS.currentCosmics->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.currentCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.currentCosmics->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.currentCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.currentCosmics->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.currentCosmics->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.currentCosmics->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.currentCosmics->Fill(chosenRecoParticleCRUMBS.purity);
        } else if(DLCurrent == 5 && signal == 3){
            numEventsCRUMBSRecoParticle.nueCosmics++;
            numPFPsCRUMBS.nueCosmics->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.nueCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.nueCosmics->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.nueCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            slicePurityCRUMBS.nueCosmics->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.nueCosmics->Fill(chosenSliceCompletenessCRUMBS);
            highestPFPCompletenessCRUMBS.nueCosmics->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.nueCosmics->Fill(chosenRecoParticleCRUMBS.purity);
        }

        //printf("Number of Slices in event: %f\n", numSlicesInEvent);
    }

    //printf("Signal Events:\nNumber of events with a true neutrino within the TPC:\nUBoone: %i out of %i\nCurrent: %i out of %i\n", numEventsTrueNeutrino.ubooneSignal, numEventsTotal.ubooneSignal, numEventsTrueNeutrino.currentSignal, numEventsTotal.currentSignal);
    //printf("Number of events with a slice:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSlices.ubooneSignal, numEventsTotal.ubooneSignal, numEventsSlices.currentSignal, numEventsTotal.currentSignal);
    //printf("Number of events with a reco neutrino:\nUboone:%i out of %i\nCurrent: %i out of %i\n", numEventsRecoNeutrino.ubooneSignal, numEventsTotal.ubooneSignal, numEventsRecoNeutrino.currentSignal, numEventsTotal.currentSignal);
    //printf("Number of events where the number of slices with a CRUMBS score != number of reco neutrinos:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSliceNotEqualNeutrino.ubooneSignal, numEventsRecoNeutrino.ubooneSignal, numEventsSliceNotEqualNeutrino.currentSignal, numEventsRecoNeutrino.currentSignal);
    
    //printf("BNB Events:\nNumber of events with a true neutrino within the TPC:\nUBoone: %i out of %i\nCurrent: %i out of %i\n", numEventsTrueNeutrino.ubooneBNB, numEventsTotal.ubooneBNB, numEventsTrueNeutrino.currentBNB, numEventsTotal.currentBNB);
    //printf("Number of events with a slice:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSlices.ubooneBNB, numEventsTotal.ubooneBNB, numEventsSlices.currentBNB, numEventsTotal.currentBNB);
    //printf("Number of events with a reco neutrino:\nUboone:%i out of %i\nCurrent: %i out of %i\n", numEventsRecoNeutrino.ubooneBNB, numEventsTotal.ubooneBNB, numEventsRecoNeutrino.currentBNB, numEventsTotal.currentBNB);
    //printf("Number of events where the number of slices with a CRUMBS score != number of reco neutrinos:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSliceNotEqualNeutrino.ubooneBNB, numEventsRecoNeutrino.ubooneBNB, numEventsSliceNotEqualNeutrino.currentBNB, numEventsRecoNeutrino.currentBNB);
    
    //printf("Intime Cosmic Events:\nNumber of events with a true neutrino within the TPC:\nUBoone: %i out of %i\nCurrent: %i out of %i\n", numEventsTrueNeutrino.ubooneCosmics, numEventsTotal.ubooneCosmics, numEventsTrueNeutrino.currentCosmics, numEventsTotal.currentCosmics);
    //printf("Number of events with a slice:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSlices.ubooneCosmics, numEventsTotal.ubooneCosmics, numEventsSlices.currentCosmics, numEventsTotal.currentCosmics);
    //printf("Number of events with a reco neutrino:\nUboone:%i out of %i\nCurrent: %i out of %i\n", numEventsRecoNeutrino.ubooneCosmics, numEventsTotal.ubooneCosmics, numEventsRecoNeutrino.currentCosmics, numEventsTotal.currentCosmics);
    //printf("Number of events where the number of slices with a CRUMBS score != number of reco neutrinos:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSliceNotEqualNeutrino.ubooneCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsSliceNotEqualNeutrino.currentCosmics, numEventsRecoNeutrino.currentCosmics);
    
    int drawLine = 1;
    int left = 0;
    int right = 1;

    double currentMaxEffPur = -999999;
    double ubooneMaxEffPur = -999999;
    double nueMaxEffPur = -999999;

    //efficiency(.currentSignal, .ubooneSignal, .currentBNB, .ubooneBNB, .currentCosmics, .ubooneCosmics, counter.currentSignal, .ubooneSignal, .currentBNB, .ubooneBNB, .currentCosmics, .ubooneCosmics, 0, 1, 999, 999, (base_path + "_eff.pdf").c_str(), (base_path + "_rej.pdf").c_str(), (base_path + "_effrejpur.pdf").c_str(), (base_path + "_pur.pdf").c_str(), (base_path + "_effpur.pdf").c_str(), "bottomRight", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "x");

    styleDraw(numSlices.canvas, numSlices.currentSignal, numSlices.ubooneSignal, numSlices.nueSignal, numSlices.currentBNB, numSlices.ubooneBNB, numSlices.nueBNB, numSlices.currentCosmics, numSlices.ubooneCosmics, numSlices.nueCosmics, 999, 999, 999, 999, (base_path + "numSlices_dist.pdf").c_str(), "topRight");
    weighted(numSlices.currentSignal, numSlices.ubooneSignal, numSlices.nueSignal, numSlices.currentBNB, numSlices.ubooneBNB, numSlices.nueBNB, numSlices.currentCosmics, numSlices.ubooneCosmics, numSlices.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "numSlices_weighted.pdf").c_str(), "topRight");
    efficiency(numSlices.currentSignal, numSlices.ubooneSignal, numSlices.nueSignal, numSlices.currentBNB, numSlices.ubooneBNB, numSlices.nueBNB, numSlices.currentCosmics, numSlices.ubooneCosmics, numSlices.nueCosmics, numEventsSlices.currentSignal, numEventsSlices.ubooneSignal, numEventsSlices.nueSignal, numEventsSlices.currentBNB, numEventsSlices.ubooneBNB, numEventsSlices.nueBNB, numEventsSlices.currentCosmics, numEventsSlices.ubooneCosmics, numEventsSlices.nueCosmics, 0, 1, 999, 999, (base_path + "numSlices_eff.pdf").c_str(), (base_path + "numSlices_rej.pdf").c_str(), (base_path + "numSlices_effrejpur.pdf").c_str(), (base_path + "numSlices_pur.pdf").c_str(), (base_path + "numSlices_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Number of Slices", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
 
    styleDraw(numSlicesCRUMBS.canvas, numSlicesCRUMBS.currentSignal, numSlicesCRUMBS.ubooneSignal, numSlicesCRUMBS.nueSignal, numSlicesCRUMBS.currentBNB, numSlicesCRUMBS.ubooneBNB, numSlicesCRUMBS.nueBNB, numSlicesCRUMBS.currentCosmics, numSlicesCRUMBS.ubooneCosmics, numSlicesCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "numCRUMBSSlices_dist.pdf").c_str(), "topRight");
    weighted(numSlicesCRUMBS.currentSignal, numSlicesCRUMBS.ubooneSignal, numSlicesCRUMBS.nueSignal, numSlicesCRUMBS.currentBNB, numSlicesCRUMBS.ubooneBNB, numSlicesCRUMBS.nueBNB, numSlicesCRUMBS.currentCosmics, numSlicesCRUMBS.ubooneCosmics, numSlicesCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "numCRUMBSSlices_weighted.pdf").c_str(), "topRight");
    efficiency(numSlicesCRUMBS.currentSignal, numSlicesCRUMBS.ubooneSignal, numSlicesCRUMBS.nueSignal, numSlicesCRUMBS.currentBNB, numSlicesCRUMBS.ubooneBNB, numSlicesCRUMBS.nueBNB, numSlicesCRUMBS.currentCosmics, numSlicesCRUMBS.ubooneCosmics, numSlicesCRUMBS.nueCosmics, numEventsSlices.currentSignal, numEventsSlices.ubooneSignal, numEventsSlices.nueSignal, numEventsSlices.currentBNB, numEventsSlices.ubooneBNB, numEventsSlices.nueBNB, numEventsSlices.currentCosmics, numEventsSlices.ubooneCosmics, numEventsSlices.nueCosmics, 0, 1, 999, 999, (base_path + "numCRUMBSSlices_eff.pdf").c_str(), (base_path + "numCRUMBSSlices_rej.pdf").c_str(), (base_path + "numCRUMBSSlices_effrejpur.pdf").c_str(), (base_path + "numCRUMBSSlices_pur.pdf").c_str(), (base_path + "numCRUMBSSlices_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Number of Slices with a CRUMBS Score", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
   
    styleDraw(numRecoNeutrinos.canvas, numRecoNeutrinos.currentSignal, numRecoNeutrinos.ubooneSignal, numRecoNeutrinos.nueSignal, numRecoNeutrinos.currentBNB, numRecoNeutrinos.ubooneBNB, numRecoNeutrinos.nueBNB, numRecoNeutrinos.currentCosmics, numRecoNeutrinos.ubooneCosmics, numRecoNeutrinos.nueCosmics, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_dist.pdf").c_str(), "topRight");
    weighted(numRecoNeutrinos.currentSignal, numRecoNeutrinos.ubooneSignal, numRecoNeutrinos.nueSignal, numRecoNeutrinos.currentBNB, numRecoNeutrinos.ubooneBNB, numRecoNeutrinos.nueBNB, numRecoNeutrinos.currentCosmics, numRecoNeutrinos.ubooneCosmics, numRecoNeutrinos.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_weighted.pdf").c_str(), "topRight");
    efficiency(numRecoNeutrinos.currentSignal, numRecoNeutrinos.ubooneSignal, numRecoNeutrinos.nueSignal, numRecoNeutrinos.currentBNB, numRecoNeutrinos.ubooneBNB, numRecoNeutrinos.nueBNB, numRecoNeutrinos.currentCosmics, numRecoNeutrinos.ubooneCosmics, numRecoNeutrinos.nueCosmics, numEventsSlices.currentSignal, numEventsSlices.ubooneSignal, numEventsSlices.nueSignal, numEventsSlices.currentBNB, numEventsSlices.ubooneBNB, numEventsSlices.nueBNB, numEventsSlices.currentCosmics, numEventsSlices.ubooneCosmics, numEventsSlices.nueCosmics, 0, 1, 999, 999, (base_path + "numRecoNeutrinos_eff.pdf").c_str(), (base_path + "numRecoNeutrinos_rej.pdf").c_str(), (base_path + "numRecoNeutrinos_effrejpur.pdf").c_str(), (base_path + "numRecoNeutrinos_pur.pdf").c_str(), (base_path + "numRecoNeutrinos_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Number of Reco Neutrinos", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(primaryVertexXCRUMBS.canvas, primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "primaryVertexX_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "primaryVertexX_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 999, 999, (base_path + "primaryVertexX_eff.pdf").c_str(), (base_path + "primaryVertexX_rej.pdf").c_str(), (base_path + "primaryVertexX_effrejpur.pdf").c_str(), (base_path + "primaryVertexX_pur.pdf").c_str(), (base_path + "primaryVertexX_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex X Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(primaryVertexYCRUMBS.canvas, primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "primaryVertexY_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "primaryVertexY_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 999, 999, (base_path + "primaryVertexY_eff.pdf").c_str(), (base_path + "primaryVertexY_rej.pdf").c_str(), (base_path + "primaryVertexY_effrejpur.pdf").c_str(), (base_path + "primaryVertexY_pur.pdf").c_str(), (base_path + "primaryVertexY_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Y Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(primaryVertexZCRUMBS.canvas, primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "primaryVertexZ_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "primaryVertexZ_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 999, 999, (base_path + "primaryVertexZ_eff.pdf").c_str(), (base_path + "primaryVertexZ_rej.pdf").c_str(), (base_path + "primaryVertexZ_effrejpur.pdf").c_str(), (base_path + "primaryVertexZ_pur.pdf").c_str(), (base_path + "primaryVertexZ_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Z Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    /*
    styleDraw(primaryVertexXCRUMBS.canvas, primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, 999, 999, 186.3, 201.3, (base_path + "primaryVertexXPOS_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 186.3, 201.3, (base_path + "primaryVertexXPOS_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 186.3, 201.3, (base_path + "primaryVertexXPOS_eff.pdf").c_str(), (base_path + "primaryVertexXPOS_rej.pdf").c_str(), (base_path + "primaryVertexXPOS_effrejpur.pdf").c_str(), (base_path + "primaryVertexXPOS_pur.pdf").c_str(), (base_path + "primaryVertexXPOS_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex X Coord", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    styleDraw(primaryVertexYCRUMBS.canvas, primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, 999, 999, 188.8, 203.8, (base_path + "primaryVertexYPOS_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 188.8, 203.8, (base_path + "primaryVertexYPOS_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 188.8, 203.8, (base_path + "primaryVertexYPOS_eff.pdf").c_str(), (base_path + "primaryVertexYPOS_rej.pdf").c_str(), (base_path + "primaryVertexYPOS_effrejpur.pdf").c_str(), (base_path + "primaryVertexYPOS_pur.pdf").c_str(), (base_path + "primaryVertexYPOS_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Y Coord", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    styleDraw(primaryVertexZCRUMBS.canvas, primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, 999, 999, 494.4, 509.4, (base_path + "primaryVertexZPOS_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 494.4, 509.4, (base_path + "primaryVertexZPOS_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 494.4, 509.4, (base_path + "primaryVertexZPOS_eff.pdf").c_str(), (base_path + "primaryVertexZPOS_rej.pdf").c_str(), (base_path + "primaryVertexZPOS_effrejpur.pdf").c_str(), (base_path + "primaryVertexZPOS_pur.pdf").c_str(), (base_path + "primaryVertexZPOS_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Z Coord", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    
    styleDraw(primaryVertexXCRUMBS.canvas, primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, 999, 999, -201.3, -186.3, (base_path + "primaryVertexXNEGNEG_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, -201.3, -186.3, (base_path + "primaryVertexX_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexXCRUMBS.currentSignal, primaryVertexXCRUMBS.ubooneSignal, primaryVertexXCRUMBS.nueSignal, primaryVertexXCRUMBS.currentBNB, primaryVertexXCRUMBS.ubooneBNB, primaryVertexXCRUMBS.nueBNB, primaryVertexXCRUMBS.currentCosmics, primaryVertexXCRUMBS.ubooneCosmics, primaryVertexXCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, -201.3, -186.3, (base_path + "primaryVertexXNEG_eff.pdf").c_str(), (base_path + "primaryVertexXNEG_rej.pdf").c_str(), (base_path + "primaryVertexXNEG_effrejpur.pdf").c_str(), (base_path + "primaryVertexXNEG_pur.pdf").c_str(), (base_path + "primaryVertexXNEG_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex X Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    styleDraw(primaryVertexYCRUMBS.canvas, primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, 999, 999, -203.8, -188.8, (base_path + "primaryVertexYNEGNEG_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, -203.8, -188.8, (base_path + "primaryVertexY_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexYCRUMBS.currentSignal, primaryVertexYCRUMBS.ubooneSignal, primaryVertexYCRUMBS.nueSignal, primaryVertexYCRUMBS.currentBNB, primaryVertexYCRUMBS.ubooneBNB, primaryVertexYCRUMBS.nueBNB, primaryVertexYCRUMBS.currentCosmics, primaryVertexYCRUMBS.ubooneCosmics, primaryVertexYCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, -203.8, -188.8, (base_path + "primaryVertexYNEG_eff.pdf").c_str(), (base_path + "primaryVertexYNEG_rej.pdf").c_str(), (base_path + "primaryVertexYNEG_effrejpur.pdf").c_str(), (base_path + "primaryVertexYNEG_pur.pdf").c_str(), (base_path + "primaryVertexYNEG_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Y Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    styleDraw(primaryVertexZCRUMBS.canvas, primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, 999, 999, 0, 15, (base_path + "primaryVertexZNEGNEG_dist.pdf").c_str(), "topRight");
    weighted(primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 0, 15, (base_path + "primaryVertexZ_weighted.pdf").c_str(), "topRight");
    efficiency(primaryVertexZCRUMBS.currentSignal, primaryVertexZCRUMBS.ubooneSignal, primaryVertexZCRUMBS.nueSignal, primaryVertexZCRUMBS.currentBNB, primaryVertexZCRUMBS.ubooneBNB, primaryVertexZCRUMBS.nueBNB, primaryVertexZCRUMBS.currentCosmics, primaryVertexZCRUMBS.ubooneCosmics, primaryVertexZCRUMBS.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0, 1, 0, 15, (base_path + "primaryVertexZNEG_eff.pdf").c_str(), (base_path + "primaryVertexZNEG_rej.pdf").c_str(), (base_path + "primaryVertexZNEG_effrejpur.pdf").c_str(), (base_path + "primaryVertexZNEG_pur.pdf").c_str(), (base_path + "primaryVertexZNEG_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Z Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    */
    
    styleDraw(primaryVertexXCRUMBSPositive.canvas, primaryVertexXCRUMBSPositive.currentSignal, primaryVertexXCRUMBSPositive.ubooneSignal, primaryVertexXCRUMBSPositive.nueSignal, primaryVertexXCRUMBSPositive.currentBNB, primaryVertexXCRUMBSPositive.ubooneBNB, primaryVertexXCRUMBSPositive.nueBNB, primaryVertexXCRUMBSPositive.currentCosmics, primaryVertexXCRUMBSPositive.ubooneCosmics, primaryVertexXCRUMBSPositive.nueCosmics, 999, 999, 180, 201.3, (base_path + "primaryVertexX_dist_positive.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexXCRUMBSPositive.currentSignal, primaryVertexXCRUMBSPositive.ubooneSignal, primaryVertexXCRUMBSPositive.nueSignal, primaryVertexXCRUMBSPositive.currentBNB, primaryVertexXCRUMBSPositive.ubooneBNB, primaryVertexXCRUMBSPositive.nueBNB, primaryVertexXCRUMBSPositive.currentCosmics, primaryVertexXCRUMBSPositive.ubooneCosmics, primaryVertexXCRUMBSPositive.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 180, 201.3, (base_path + "primaryVertexX_weighted_positive.pdf").c_str(), "bottomLeft");
    efficiency(primaryVertexXCRUMBSPositive.currentSignal, primaryVertexXCRUMBSPositive.ubooneSignal, primaryVertexXCRUMBSPositive.nueSignal, primaryVertexXCRUMBSPositive.currentBNB, primaryVertexXCRUMBSPositive.ubooneBNB, primaryVertexXCRUMBSPositive.nueBNB, primaryVertexXCRUMBSPositive.currentCosmics, primaryVertexXCRUMBSPositive.ubooneCosmics, primaryVertexXCRUMBSPositive.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0.00004, 0.00005, 180, 201.3, (base_path + "primaryVertexX_eff_positive.pdf").c_str(), (base_path + "primaryVertexX_rej_positive.pdf").c_str(), (base_path + "primaryVertexX_effrejpur_positive.pdf").c_str(), (base_path + "primaryVertexX_pur_positive.pdf").c_str(), (base_path + "primaryVertexX_effpur_positive.pdf").c_str(), "bottomLeft", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex X Coord", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    std::cout << "Vertex X Positive: Current = " << currentMaxEffPur << ", Uboone = " << ubooneMaxEffPur << ", Nue = " << nueMaxEffPur << std::endl;
    styleDraw(primaryVertexYCRUMBSPositive.canvas, primaryVertexYCRUMBSPositive.currentSignal, primaryVertexYCRUMBSPositive.ubooneSignal, primaryVertexYCRUMBSPositive.nueSignal, primaryVertexYCRUMBSPositive.currentBNB, primaryVertexYCRUMBSPositive.ubooneBNB, primaryVertexYCRUMBSPositive.nueBNB, primaryVertexYCRUMBSPositive.currentCosmics, primaryVertexYCRUMBSPositive.ubooneCosmics, primaryVertexYCRUMBSPositive.nueCosmics, 999, 999, 180, 203.8, (base_path + "primaryVertexY_dist_positive.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexYCRUMBSPositive.currentSignal, primaryVertexYCRUMBSPositive.ubooneSignal, primaryVertexYCRUMBSPositive.nueSignal, primaryVertexYCRUMBSPositive.currentBNB, primaryVertexYCRUMBSPositive.ubooneBNB, primaryVertexYCRUMBSPositive.nueBNB, primaryVertexYCRUMBSPositive.currentCosmics, primaryVertexYCRUMBSPositive.ubooneCosmics, primaryVertexYCRUMBSPositive.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 180, 203.8, (base_path + "primaryVertexY_weighted_positive.pdf").c_str(), "bottomLeft");
    efficiency(primaryVertexYCRUMBSPositive.currentSignal, primaryVertexYCRUMBSPositive.ubooneSignal, primaryVertexYCRUMBSPositive.nueSignal, primaryVertexYCRUMBSPositive.currentBNB, primaryVertexYCRUMBSPositive.ubooneBNB, primaryVertexYCRUMBSPositive.nueBNB, primaryVertexYCRUMBSPositive.currentCosmics, primaryVertexYCRUMBSPositive.ubooneCosmics, primaryVertexYCRUMBSPositive.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0.000045, 0.000055, 180, 203.8, (base_path + "primaryVertexY_eff_positive.pdf").c_str(), (base_path + "primaryVertexY_rej_positive.pdf").c_str(), (base_path + "primaryVertexY_effrejpur_positive.pdf").c_str(), (base_path + "primaryVertexY_pur_positive.pdf").c_str(), (base_path + "primaryVertexY_effpur_positive.pdf").c_str(), "bottomLeft", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Y Coord", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    std::cout << "Vertex Y Positive: Current = " << currentMaxEffPur << ", Uboone = " << ubooneMaxEffPur << ", Nue = " << nueMaxEffPur << std::endl;
    styleDraw(primaryVertexZCRUMBSPositive.canvas, primaryVertexZCRUMBSPositive.currentSignal, primaryVertexZCRUMBSPositive.ubooneSignal, primaryVertexZCRUMBSPositive.nueSignal, primaryVertexZCRUMBSPositive.currentBNB, primaryVertexZCRUMBSPositive.ubooneBNB, primaryVertexZCRUMBSPositive.nueBNB, primaryVertexZCRUMBSPositive.currentCosmics, primaryVertexZCRUMBSPositive.ubooneCosmics, primaryVertexZCRUMBSPositive.nueCosmics, 999, 999, 490, 509.4, (base_path + "primaryVertexZ_dist_positive.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexZCRUMBSPositive.currentSignal, primaryVertexZCRUMBSPositive.ubooneSignal, primaryVertexZCRUMBSPositive.nueSignal, primaryVertexZCRUMBSPositive.currentBNB, primaryVertexZCRUMBSPositive.ubooneBNB, primaryVertexZCRUMBSPositive.nueBNB, primaryVertexZCRUMBSPositive.currentCosmics, primaryVertexZCRUMBSPositive.ubooneCosmics, primaryVertexZCRUMBSPositive.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 490, 509.4, (base_path + "primaryVertexZ_weighted_positive.pdf").c_str(), "bottomLeft");
    efficiency(primaryVertexZCRUMBSPositive.currentSignal, primaryVertexZCRUMBSPositive.ubooneSignal, primaryVertexZCRUMBSPositive.nueSignal, primaryVertexZCRUMBSPositive.currentBNB, primaryVertexZCRUMBSPositive.ubooneBNB, primaryVertexZCRUMBSPositive.nueBNB, primaryVertexZCRUMBSPositive.currentCosmics, primaryVertexZCRUMBSPositive.ubooneCosmics, primaryVertexZCRUMBSPositive.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0.000045, 0.000055, 490, 509.4, (base_path + "primaryVertexZ_eff_positive.pdf").c_str(), (base_path + "primaryVertexZ_rej_positive.pdf").c_str(), (base_path + "primaryVertexZ_effrejpur_positive.pdf").c_str(), (base_path + "primaryVertexZ_pur_positive.pdf").c_str(), (base_path + "primaryVertexZ_effpur_positive.pdf").c_str(), "bottomLeft", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Z Coord", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    std::cout << "Vertex Z Positive: Current = " << currentMaxEffPur << ", Uboone = " << ubooneMaxEffPur << ", Nue = " << nueMaxEffPur << std::endl;

    styleDraw(primaryVertexXCRUMBSNegative.canvas, primaryVertexXCRUMBSNegative.currentSignal, primaryVertexXCRUMBSNegative.ubooneSignal, primaryVertexXCRUMBSNegative.nueSignal, primaryVertexXCRUMBSNegative.currentBNB, primaryVertexXCRUMBSNegative.ubooneBNB, primaryVertexXCRUMBSNegative.nueBNB, primaryVertexXCRUMBSNegative.currentCosmics, primaryVertexXCRUMBSNegative.ubooneCosmics, primaryVertexXCRUMBSNegative.nueCosmics, 999, 999, -201.3, -180, (base_path + "primaryVertexX_dist_negative.pdf").c_str(), "bottomRight");
    weighted(primaryVertexXCRUMBSNegative.currentSignal, primaryVertexXCRUMBSNegative.ubooneSignal, primaryVertexXCRUMBSNegative.nueSignal, primaryVertexXCRUMBSNegative.currentBNB, primaryVertexXCRUMBSNegative.ubooneBNB, primaryVertexXCRUMBSNegative.nueBNB, primaryVertexXCRUMBSNegative.currentCosmics, primaryVertexXCRUMBSNegative.ubooneCosmics, primaryVertexXCRUMBSNegative.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, -201.3, -180, (base_path + "primaryVertexX_weighted_negative.pdf").c_str(), "bottomRight");
    efficiency(primaryVertexXCRUMBSNegative.currentSignal, primaryVertexXCRUMBSNegative.ubooneSignal, primaryVertexXCRUMBSNegative.nueSignal, primaryVertexXCRUMBSNegative.currentBNB, primaryVertexXCRUMBSNegative.ubooneBNB, primaryVertexXCRUMBSNegative.nueBNB, primaryVertexXCRUMBSNegative.currentCosmics, primaryVertexXCRUMBSNegative.ubooneCosmics, primaryVertexXCRUMBSNegative.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0.000045, 0.000055, -201.3, -180, (base_path + "primaryVertexX_eff_negative.pdf").c_str(), (base_path + "primaryVertexX_rej_negative.pdf").c_str(), (base_path + "primaryVertexX_effrejpur_negative.pdf").c_str(), (base_path + "primaryVertexX_pur_negative.pdf").c_str(), (base_path + "primaryVertexX_effpur_negative.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex X Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    std::cout << "Vertex X Negative: Current = " << currentMaxEffPur << ", Uboone = " << ubooneMaxEffPur << ", Nue = " << nueMaxEffPur << std::endl;
    styleDraw(primaryVertexYCRUMBSNegative.canvas, primaryVertexYCRUMBSNegative.currentSignal, primaryVertexYCRUMBSNegative.ubooneSignal, primaryVertexYCRUMBSNegative.nueSignal, primaryVertexYCRUMBSNegative.currentBNB, primaryVertexYCRUMBSNegative.ubooneBNB, primaryVertexYCRUMBSNegative.nueBNB, primaryVertexYCRUMBSNegative.currentCosmics, primaryVertexYCRUMBSNegative.ubooneCosmics, primaryVertexYCRUMBSNegative.nueCosmics, 999, 999, -203.8, -180, (base_path + "primaryVertexY_dist_negative.pdf").c_str(), "bottomRight");
    weighted(primaryVertexYCRUMBSNegative.currentSignal, primaryVertexYCRUMBSNegative.ubooneSignal, primaryVertexYCRUMBSNegative.nueSignal, primaryVertexYCRUMBSNegative.currentBNB, primaryVertexYCRUMBSNegative.ubooneBNB, primaryVertexYCRUMBSNegative.nueBNB, primaryVertexYCRUMBSNegative.currentCosmics, primaryVertexYCRUMBSNegative.ubooneCosmics, primaryVertexYCRUMBSNegative.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, -203.8, -180, (base_path + "primaryVertexY_weighted_negative.pdf").c_str(), "bottomRight");
    efficiency(primaryVertexYCRUMBSNegative.currentSignal, primaryVertexYCRUMBSNegative.ubooneSignal, primaryVertexYCRUMBSNegative.nueSignal, primaryVertexYCRUMBSNegative.currentBNB, primaryVertexYCRUMBSNegative.ubooneBNB, primaryVertexYCRUMBSNegative.nueBNB, primaryVertexYCRUMBSNegative.currentCosmics, primaryVertexYCRUMBSNegative.ubooneCosmics, primaryVertexYCRUMBSNegative.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0.000045, 0.000055, -203.8, -180, (base_path + "primaryVertexY_eff_negative.pdf").c_str(), (base_path + "primaryVertexY_rej_negative.pdf").c_str(), (base_path + "primaryVertexY_effrejpur_negative.pdf").c_str(), (base_path + "primaryVertexY_pur_negative.pdf").c_str(), (base_path + "primaryVertexY_effpur_negative.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Y Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    std::cout << "Vertex Y Negative: Current = " << currentMaxEffPur << ", Uboone = " << ubooneMaxEffPur << ", Nue = " << nueMaxEffPur << std::endl;
    styleDraw(primaryVertexZCRUMBSNegative.canvas, primaryVertexZCRUMBSNegative.currentSignal, primaryVertexZCRUMBSNegative.ubooneSignal, primaryVertexZCRUMBSNegative.nueSignal, primaryVertexZCRUMBSNegative.currentBNB, primaryVertexZCRUMBSNegative.ubooneBNB, primaryVertexZCRUMBSNegative.nueBNB, primaryVertexZCRUMBSNegative.currentCosmics, primaryVertexZCRUMBSNegative.ubooneCosmics, primaryVertexZCRUMBSNegative.nueCosmics, 999, 999, 0, 20, (base_path + "primaryVertexZ_dist_negative.pdf").c_str(), "bottomRight");
    weighted(primaryVertexZCRUMBSNegative.currentSignal, primaryVertexZCRUMBSNegative.ubooneSignal, primaryVertexZCRUMBSNegative.nueSignal, primaryVertexZCRUMBSNegative.currentBNB, primaryVertexZCRUMBSNegative.ubooneBNB, primaryVertexZCRUMBSNegative.nueBNB, primaryVertexZCRUMBSNegative.currentCosmics, primaryVertexZCRUMBSNegative.ubooneCosmics, primaryVertexZCRUMBSNegative.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 0, 20, (base_path + "primaryVertexZ_weighted_negative.pdf").c_str(), "bottomRight");
    efficiency(primaryVertexZCRUMBSNegative.currentSignal, primaryVertexZCRUMBSNegative.ubooneSignal, primaryVertexZCRUMBSNegative.nueSignal, primaryVertexZCRUMBSNegative.currentBNB, primaryVertexZCRUMBSNegative.ubooneBNB, primaryVertexZCRUMBSNegative.nueBNB, primaryVertexZCRUMBSNegative.currentCosmics, primaryVertexZCRUMBSNegative.ubooneCosmics, primaryVertexZCRUMBSNegative.nueCosmics, numEventsRecoNeutrino.currentSignal, numEventsRecoNeutrino.ubooneSignal, numEventsRecoNeutrino.nueSignal, numEventsRecoNeutrino.currentBNB, numEventsRecoNeutrino.ubooneBNB, numEventsRecoNeutrino.nueBNB, numEventsRecoNeutrino.currentCosmics, numEventsRecoNeutrino.ubooneCosmics, numEventsRecoNeutrino.nueCosmics, 0.000045, 0.000055, 0, 20, (base_path + "primaryVertexZ_eff_negative.pdf").c_str(), (base_path + "primaryVertexZ_rej_negative.pdf").c_str(), (base_path + "primaryVertexZ_effrejpur_negative.pdf").c_str(), (base_path + "primaryVertexZ_pur_negative.pdf").c_str(), (base_path + "primaryVertexZ_effpur_negative.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Primary Vertex Z Coord", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    std::cout << "Vertex Z Negative: Current = " << currentMaxEffPur << ", Uboone = " << ubooneMaxEffPur << ", Nue = " << nueMaxEffPur << std::endl;

    // CRUMBS Slice Score Plots
    styleDraw(sliceScoreCRUMBS.canvas, sliceScoreCRUMBS.currentSignal, sliceScoreCRUMBS.ubooneSignal, sliceScoreCRUMBS.nueSignal, sliceScoreCRUMBS.currentBNB, sliceScoreCRUMBS.ubooneBNB, sliceScoreCRUMBS.nueBNB, sliceScoreCRUMBS.currentCosmics, sliceScoreCRUMBS.ubooneCosmics, sliceScoreCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_dist.pdf").c_str(), "topLeft");
    weighted(sliceScoreCRUMBS.currentSignal, sliceScoreCRUMBS.ubooneSignal, sliceScoreCRUMBS.nueSignal, sliceScoreCRUMBS.currentBNB, sliceScoreCRUMBS.ubooneBNB, sliceScoreCRUMBS.nueBNB, sliceScoreCRUMBS.currentCosmics, sliceScoreCRUMBS.ubooneCosmics, sliceScoreCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_weighted.pdf").c_str(), "topLeft");
    efficiency(sliceScoreCRUMBS.currentSignal, sliceScoreCRUMBS.ubooneSignal, sliceScoreCRUMBS.nueSignal, sliceScoreCRUMBS.currentBNB, sliceScoreCRUMBS.ubooneBNB, sliceScoreCRUMBS.nueBNB, sliceScoreCRUMBS.currentCosmics, sliceScoreCRUMBS.ubooneCosmics, sliceScoreCRUMBS.nueCosmics, numEventsSlices.currentSignal, numEventsSlices.ubooneSignal, numEventsSlices.nueSignal, numEventsSlices.currentBNB, numEventsSlices.ubooneBNB, numEventsSlices.nueBNB, numEventsSlices.currentCosmics, numEventsSlices.ubooneCosmics, numEventsSlices.nueCosmics, 0, 1, 999, 999, (base_path + "sliceScoreCRUMBS_eff.pdf").c_str(), (base_path + "sliceScoreCRUMBS_rej.pdf").c_str(), (base_path + "sliceScoreCRUMBS_effrejpur.pdf").c_str(), (base_path + "sliceScoreCRUMBS_pur.pdf").c_str(), (base_path + "sliceScoreCRUMBS_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "CRUMBS Score", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    
    // CRUMBS PFP Plots
    styleDraw(numPFPsCRUMBS.canvas, numPFPsCRUMBS.currentSignal, numPFPsCRUMBS.ubooneSignal, numPFPsCRUMBS.nueSignal, numPFPsCRUMBS.currentBNB, numPFPsCRUMBS.ubooneBNB, numPFPsCRUMBS.nueBNB, numPFPsCRUMBS.currentCosmics, numPFPsCRUMBS.ubooneCosmics, numPFPsCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_dist.pdf").c_str(), "topRight");
    weighted(numPFPsCRUMBS.currentSignal, numPFPsCRUMBS.ubooneSignal, numPFPsCRUMBS.nueSignal, numPFPsCRUMBS.currentBNB, numPFPsCRUMBS.ubooneBNB, numPFPsCRUMBS.nueBNB, numPFPsCRUMBS.currentCosmics, numPFPsCRUMBS.ubooneCosmics, numPFPsCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_weighted.pdf").c_str(), "topRight");
    efficiency(numPFPsCRUMBS.currentSignal, numPFPsCRUMBS.ubooneSignal, numPFPsCRUMBS.nueSignal, numPFPsCRUMBS.currentBNB, numPFPsCRUMBS.ubooneBNB, numPFPsCRUMBS.nueBNB, numPFPsCRUMBS.currentCosmics, numPFPsCRUMBS.ubooneCosmics, numPFPsCRUMBS.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "numPFPsCRUMBS_eff.pdf").c_str(), (base_path + "numPFPsCRUMBS_rej.pdf").c_str(), (base_path + "numPFPsCRUMBS_effrejpur.pdf").c_str(), (base_path + "numPFPsCRUMBS_pur.pdf").c_str(), (base_path + "numPFPsCRUMBS_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Number of PFPs", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(ratioChosenSummedEnergyCRUMBS.canvas, ratioChosenSummedEnergyCRUMBS.currentSignal, ratioChosenSummedEnergyCRUMBS.ubooneSignal, ratioChosenSummedEnergyCRUMBS.nueSignal, ratioChosenSummedEnergyCRUMBS.currentBNB, ratioChosenSummedEnergyCRUMBS.ubooneBNB, ratioChosenSummedEnergyCRUMBS.nueBNB, ratioChosenSummedEnergyCRUMBS.currentCosmics, ratioChosenSummedEnergyCRUMBS.ubooneCosmics, ratioChosenSummedEnergyCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_dist.pdf").c_str(), "topLeft");
    weighted(ratioChosenSummedEnergyCRUMBS.currentSignal, ratioChosenSummedEnergyCRUMBS.ubooneSignal, ratioChosenSummedEnergyCRUMBS.nueSignal, ratioChosenSummedEnergyCRUMBS.currentBNB, ratioChosenSummedEnergyCRUMBS.ubooneBNB, ratioChosenSummedEnergyCRUMBS.nueBNB, ratioChosenSummedEnergyCRUMBS.currentCosmics, ratioChosenSummedEnergyCRUMBS.ubooneCosmics, ratioChosenSummedEnergyCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_weighted.pdf").c_str(), "topLeft");
    efficiency(ratioChosenSummedEnergyCRUMBS.currentSignal, ratioChosenSummedEnergyCRUMBS.ubooneSignal, ratioChosenSummedEnergyCRUMBS.nueSignal, ratioChosenSummedEnergyCRUMBS.currentBNB, ratioChosenSummedEnergyCRUMBS.ubooneBNB, ratioChosenSummedEnergyCRUMBS.nueBNB, ratioChosenSummedEnergyCRUMBS.currentCosmics, ratioChosenSummedEnergyCRUMBS.ubooneCosmics, ratioChosenSummedEnergyCRUMBS.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_eff.pdf").c_str(), (base_path + "ratioChosenSummedEnergyCRUMBS_rej.pdf").c_str(), (base_path + "ratioChosenSummedEnergyCRUMBS_effrejpur.pdf").c_str(), (base_path + "ratioChosenSummedEnergyCRUMBS_pur.pdf").c_str(), (base_path + "ratioChosenSummedEnergyCRUMBS_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(slicePurityCRUMBS.canvas, slicePurityCRUMBS.currentSignal, slicePurityCRUMBS.ubooneSignal, slicePurityCRUMBS.nueSignal, slicePurityCRUMBS.currentBNB, slicePurityCRUMBS.ubooneBNB, slicePurityCRUMBS.nueBNB, slicePurityCRUMBS.currentCosmics, slicePurityCRUMBS.ubooneCosmics, slicePurityCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "slicePurityCRUMBS_dist.pdf").c_str(), "topLeft");
    weighted(slicePurityCRUMBS.currentSignal, slicePurityCRUMBS.ubooneSignal, slicePurityCRUMBS.nueSignal, slicePurityCRUMBS.currentBNB, slicePurityCRUMBS.ubooneBNB, slicePurityCRUMBS.nueBNB, slicePurityCRUMBS.currentCosmics, slicePurityCRUMBS.ubooneCosmics, slicePurityCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "slicePurityCRUMBS_weighted.pdf").c_str(), "topLeft");
    
    styleDraw(sliceCompletenessCRUMBS.canvas, sliceCompletenessCRUMBS.currentSignal, sliceCompletenessCRUMBS.ubooneSignal, sliceCompletenessCRUMBS.nueSignal, sliceCompletenessCRUMBS.currentBNB, sliceCompletenessCRUMBS.ubooneBNB, sliceCompletenessCRUMBS.nueBNB, sliceCompletenessCRUMBS.currentCosmics, sliceCompletenessCRUMBS.ubooneCosmics, sliceCompletenessCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBS_dist.pdf").c_str(), "topLeft");
    weighted(sliceCompletenessCRUMBS.currentSignal, sliceCompletenessCRUMBS.ubooneSignal, sliceCompletenessCRUMBS.nueSignal, sliceCompletenessCRUMBS.currentBNB, sliceCompletenessCRUMBS.ubooneBNB, sliceCompletenessCRUMBS.nueBNB, sliceCompletenessCRUMBS.currentCosmics, sliceCompletenessCRUMBS.ubooneCosmics, sliceCompletenessCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBS_weighted.pdf").c_str(), "topLeft");

    styleDraw(highestPFPCompletenessCRUMBS.canvas, highestPFPCompletenessCRUMBS.currentSignal, highestPFPCompletenessCRUMBS.ubooneSignal, highestPFPCompletenessCRUMBS.nueSignal, highestPFPCompletenessCRUMBS.currentBNB, highestPFPCompletenessCRUMBS.ubooneBNB, highestPFPCompletenessCRUMBS.nueBNB, highestPFPCompletenessCRUMBS.currentCosmics, highestPFPCompletenessCRUMBS.ubooneCosmics, highestPFPCompletenessCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBS_dist.pdf").c_str(), "topLeft");
    weighted(highestPFPCompletenessCRUMBS.currentSignal, highestPFPCompletenessCRUMBS.ubooneSignal, highestPFPCompletenessCRUMBS.nueSignal, highestPFPCompletenessCRUMBS.currentBNB, highestPFPCompletenessCRUMBS.ubooneBNB, highestPFPCompletenessCRUMBS.nueBNB, highestPFPCompletenessCRUMBS.currentCosmics, highestPFPCompletenessCRUMBS.ubooneCosmics, highestPFPCompletenessCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBS_weighted.pdf").c_str(), "topLeft");
    efficiency(highestPFPCompletenessCRUMBS.currentSignal, highestPFPCompletenessCRUMBS.ubooneSignal, highestPFPCompletenessCRUMBS.nueSignal, highestPFPCompletenessCRUMBS.currentBNB, highestPFPCompletenessCRUMBS.ubooneBNB, highestPFPCompletenessCRUMBS.nueBNB, highestPFPCompletenessCRUMBS.currentCosmics, highestPFPCompletenessCRUMBS.ubooneCosmics, highestPFPCompletenessCRUMBS.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "highestPFPCompletenessCRUMBS_eff.pdf").c_str(), (base_path + "highestPFPCompletenessCRUMBS_rej.pdf").c_str(), (base_path + "highestPFPCompletenessCRUMBS_effrejpur.pdf").c_str(), (base_path + "highestPFPCompletenessCRUMBS_pur.pdf").c_str(), (base_path + "highestPFPCompletenessCRUMBS_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Completeness", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(highestPFPPurityCRUMBS.canvas, highestPFPPurityCRUMBS.currentSignal, highestPFPPurityCRUMBS.ubooneSignal, highestPFPPurityCRUMBS.nueSignal, highestPFPPurityCRUMBS.currentBNB, highestPFPPurityCRUMBS.ubooneBNB, highestPFPPurityCRUMBS.nueBNB, highestPFPPurityCRUMBS.currentCosmics, highestPFPPurityCRUMBS.ubooneCosmics, highestPFPPurityCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBS_dist.pdf").c_str(), "topLeft");
    weighted(highestPFPPurityCRUMBS.currentSignal, highestPFPPurityCRUMBS.ubooneSignal, highestPFPPurityCRUMBS.nueSignal, highestPFPPurityCRUMBS.currentBNB, highestPFPPurityCRUMBS.ubooneBNB, highestPFPPurityCRUMBS.nueBNB, highestPFPPurityCRUMBS.currentCosmics, highestPFPPurityCRUMBS.ubooneCosmics, highestPFPPurityCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBS_weighted.pdf").c_str(), "topLeft");
    efficiency(highestPFPPurityCRUMBS.currentSignal, highestPFPPurityCRUMBS.ubooneSignal, highestPFPPurityCRUMBS.nueSignal, highestPFPPurityCRUMBS.currentBNB, highestPFPPurityCRUMBS.ubooneBNB, highestPFPPurityCRUMBS.nueBNB, highestPFPPurityCRUMBS.currentCosmics, highestPFPPurityCRUMBS.ubooneCosmics, highestPFPPurityCRUMBS.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "highestPFPPurityCRUMBS_eff.pdf").c_str(), (base_path + "highestPFPPurityCRUMBS_rej.pdf").c_str(), (base_path + "highestPFPPurityCRUMBS_effrejpur.pdf").c_str(), (base_path + "highestPFPPurityCRUMBS_pur.pdf").c_str(), (base_path + "highestPFPPurityCRUMBS_effpur.pdf").c_str(), "bottomLeft", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, nullptr, &left, "Purity", -1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(ERecoSumThetaRecoCRUMBS.canvas, ERecoSumThetaRecoCRUMBS.currentSignal, ERecoSumThetaRecoCRUMBS.ubooneSignal, ERecoSumThetaRecoCRUMBS.nueSignal, ERecoSumThetaRecoCRUMBS.currentBNB, ERecoSumThetaRecoCRUMBS.ubooneBNB, ERecoSumThetaRecoCRUMBS.nueBNB, ERecoSumThetaRecoCRUMBS.currentCosmics, ERecoSumThetaRecoCRUMBS.ubooneCosmics, ERecoSumThetaRecoCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_dist.pdf").c_str(), "topRight", nullptr, nullptr, &drawLine, &right);
    weighted(ERecoSumThetaRecoCRUMBS.currentSignal, ERecoSumThetaRecoCRUMBS.ubooneSignal, ERecoSumThetaRecoCRUMBS.nueSignal, ERecoSumThetaRecoCRUMBS.currentBNB, ERecoSumThetaRecoCRUMBS.ubooneBNB, ERecoSumThetaRecoCRUMBS.nueBNB, ERecoSumThetaRecoCRUMBS.currentCosmics, ERecoSumThetaRecoCRUMBS.ubooneCosmics, ERecoSumThetaRecoCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_weighted.pdf").c_str(), "topRight", &drawLine, &right);
    efficiency(ERecoSumThetaRecoCRUMBS.currentSignal, ERecoSumThetaRecoCRUMBS.ubooneSignal, ERecoSumThetaRecoCRUMBS.nueSignal, ERecoSumThetaRecoCRUMBS.currentBNB, ERecoSumThetaRecoCRUMBS.ubooneBNB, ERecoSumThetaRecoCRUMBS.nueBNB, ERecoSumThetaRecoCRUMBS.currentCosmics, ERecoSumThetaRecoCRUMBS.ubooneCosmics, ERecoSumThetaRecoCRUMBS.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_eff.pdf").c_str(), (base_path + "ERecoSumThetaRecoCRUMBS_rej.pdf").c_str(), (base_path + "ERecoSumThetaRecoCRUMBS_effrejpur.pdf").c_str(), (base_path + "ERecoSumThetaRecoCRUMBS_pur.pdf").c_str(), (base_path + "ERecoSumThetaRecoCRUMBS_effpur.pdf").c_str(), "topRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, &drawLine, &right, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    //std::cout << "ERecoHighestThetaReco ____________________________________________________________________" << std::endl;
    styleDraw(ERecoHighestThetaRecoCRUMBS.canvas, ERecoHighestThetaRecoCRUMBS.currentSignal, ERecoHighestThetaRecoCRUMBS.ubooneSignal, ERecoHighestThetaRecoCRUMBS.nueSignal, ERecoHighestThetaRecoCRUMBS.currentBNB, ERecoHighestThetaRecoCRUMBS.ubooneBNB, ERecoHighestThetaRecoCRUMBS.nueBNB, ERecoHighestThetaRecoCRUMBS.currentCosmics, ERecoHighestThetaRecoCRUMBS.ubooneCosmics, ERecoHighestThetaRecoCRUMBS.nueCosmics, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_dist.pdf").c_str(), "topRight", nullptr, nullptr, &drawLine, &right);
    weighted(ERecoHighestThetaRecoCRUMBS.currentSignal, ERecoHighestThetaRecoCRUMBS.ubooneSignal, ERecoHighestThetaRecoCRUMBS.nueSignal, ERecoHighestThetaRecoCRUMBS.currentBNB, ERecoHighestThetaRecoCRUMBS.ubooneBNB, ERecoHighestThetaRecoCRUMBS.nueBNB, ERecoHighestThetaRecoCRUMBS.currentCosmics, ERecoHighestThetaRecoCRUMBS.ubooneCosmics, ERecoHighestThetaRecoCRUMBS.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_weighted.pdf").c_str(), "topRight", &drawLine, &right);
    efficiency(ERecoHighestThetaRecoCRUMBS.currentSignal, ERecoHighestThetaRecoCRUMBS.ubooneSignal, ERecoHighestThetaRecoCRUMBS.nueSignal, ERecoHighestThetaRecoCRUMBS.currentBNB, ERecoHighestThetaRecoCRUMBS.ubooneBNB, ERecoHighestThetaRecoCRUMBS.nueBNB, ERecoHighestThetaRecoCRUMBS.currentCosmics, ERecoHighestThetaRecoCRUMBS.ubooneCosmics, ERecoHighestThetaRecoCRUMBS.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_eff.pdf").c_str(), (base_path + "ERecoHighestThetaRecoCRUMBS_rej.pdf").c_str(), (base_path + "ERecoHighestThetaRecoCRUMBS_effrejpur.pdf").c_str(), (base_path + "ERecoHighestThetaRecoCRUMBS_pur.pdf").c_str(), (base_path + "ERecoHighestThetaRecoCRUMBS_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, &drawLine, &right, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    styleDraw(allRecoVertexX.canvas, allRecoVertexX.currentSignal, allRecoVertexX.ubooneSignal, allRecoVertexX.nueSignal, allRecoVertexX.currentBNB, allRecoVertexX.ubooneBNB, allRecoVertexX.nueBNB, allRecoVertexX.currentCosmics, allRecoVertexX.ubooneCosmics, allRecoVertexX.nueCosmics, 999, 999, 999, 999, (base_path + "allRecoVertexX_dist.pdf").c_str(), "topRight");
    weighted(allRecoVertexX.currentSignal, allRecoVertexX.ubooneSignal, allRecoVertexX.nueSignal, allRecoVertexX.currentBNB, allRecoVertexX.ubooneBNB, allRecoVertexX.nueBNB, allRecoVertexX.currentCosmics, allRecoVertexX.ubooneCosmics, allRecoVertexX.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "allRecoVertexX_weighted.pdf").c_str(), "topRight");
    efficiency(allRecoVertexX.currentSignal, allRecoVertexX.ubooneSignal, allRecoVertexX.nueSignal, allRecoVertexX.currentBNB, allRecoVertexX.ubooneBNB, allRecoVertexX.nueBNB, allRecoVertexX.currentCosmics, allRecoVertexX.ubooneCosmics, allRecoVertexX.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "allRecoVertexX_eff.pdf").c_str(), (base_path + "allRecoVertexX_rej.pdf").c_str(), (base_path + "allRecoVertexX_effrejpur.pdf").c_str(), (base_path + "allRecoVertexX_pur.pdf").c_str(), (base_path + "allRecoVertexX_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, &drawLine, &right, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    
    styleDraw(allRecoVertexY.canvas, allRecoVertexY.currentSignal, allRecoVertexY.ubooneSignal, allRecoVertexY.nueSignal, allRecoVertexY.currentBNB, allRecoVertexY.ubooneBNB, allRecoVertexY.nueBNB, allRecoVertexY.currentCosmics, allRecoVertexY.ubooneCosmics, allRecoVertexY.nueCosmics, 999, 999, 999, 999, (base_path + "allRecoVertexY_dist.pdf").c_str(), "topRight");
    weighted(allRecoVertexY.currentSignal, allRecoVertexY.ubooneSignal, allRecoVertexY.nueSignal, allRecoVertexY.currentBNB, allRecoVertexY.ubooneBNB, allRecoVertexY.nueBNB, allRecoVertexY.currentCosmics, allRecoVertexY.ubooneCosmics, allRecoVertexY.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "allRecoVertexY_weighted.pdf").c_str(), "topRight");
    efficiency(allRecoVertexY.currentSignal, allRecoVertexY.ubooneSignal, allRecoVertexY.nueSignal, allRecoVertexY.currentBNB, allRecoVertexY.ubooneBNB, allRecoVertexY.nueBNB, allRecoVertexY.currentCosmics, allRecoVertexY.ubooneCosmics, allRecoVertexY.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "allRecoVertexY_eff.pdf").c_str(), (base_path + "allRecoVertexY_rej.pdf").c_str(), (base_path + "allRecoVertexY_effrejpur.pdf").c_str(), (base_path + "allRecoVertexY_pur.pdf").c_str(), (base_path + "allRecoVertexY_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, &drawLine, &right, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);
    
    styleDraw(allRecoVertexZ.canvas, allRecoVertexZ.currentSignal, allRecoVertexZ.ubooneSignal, allRecoVertexZ.nueSignal, allRecoVertexZ.currentBNB, allRecoVertexZ.ubooneBNB, allRecoVertexZ.nueBNB, allRecoVertexZ.currentCosmics, allRecoVertexZ.ubooneCosmics, allRecoVertexZ.nueCosmics, 999, 999, 999, 999, (base_path + "allRecoVertexZ_dist.pdf").c_str(), "topRight");
    weighted(allRecoVertexZ.currentSignal, allRecoVertexZ.ubooneSignal, allRecoVertexZ.nueSignal, allRecoVertexZ.currentBNB, allRecoVertexZ.ubooneBNB, allRecoVertexZ.nueBNB, allRecoVertexZ.currentCosmics, allRecoVertexZ.ubooneCosmics, allRecoVertexZ.nueCosmics, signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, 999, 999, 999, 999, (base_path + "allRecoVertexZ_weighted.pdf").c_str(), "topRight");
    efficiency(allRecoVertexZ.currentSignal, allRecoVertexZ.ubooneSignal, allRecoVertexZ.nueSignal, allRecoVertexZ.currentBNB, allRecoVertexZ.ubooneBNB, allRecoVertexZ.nueBNB, allRecoVertexZ.currentCosmics, allRecoVertexZ.ubooneCosmics, allRecoVertexZ.nueCosmics, numEventsCRUMBSRecoParticle.currentSignal, numEventsCRUMBSRecoParticle.ubooneSignal, numEventsCRUMBSRecoParticle.nueSignal, numEventsCRUMBSRecoParticle.currentBNB, numEventsCRUMBSRecoParticle.ubooneBNB, numEventsCRUMBSRecoParticle.nueBNB, numEventsCRUMBSRecoParticle.currentCosmics, numEventsCRUMBSRecoParticle.ubooneCosmics, numEventsCRUMBSRecoParticle.nueCosmics, 0, 1, 999, 999, (base_path + "allRecoVertexZ_eff.pdf").c_str(), (base_path + "allRecoVertexZ_rej.pdf").c_str(), (base_path + "allRecoVertexZ_effrejpur.pdf").c_str(), (base_path + "allRecoVertexZ_pur.pdf").c_str(), (base_path + "allRecoVertexZ_effpur.pdf").c_str(), "bottomRight", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue, &drawLine, &right, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1, &currentMaxEffPur, &ubooneMaxEffPur, &nueMaxEffPur);

    printf("\nInteractions in Signal Sample:\nUnknown = %i\nQE = %i\nRes = %i\nDIS = %i\nCoh = %i\nCoh Elastic = %i\nElastic Scattering = %i\nIMD Annihilation = %i\nInverse Beta Decay = %i\nGlashow Resonance = %i\nAM Nu Gamma = %i\nMEC = %i\nDiffractive = %i\nEM = %i\nWeak Mix = %i\n\n", interactionsSignal.Unknown, interactionsSignal.QE, interactionsSignal.Res, interactionsSignal.DIS, interactionsSignal.Coh, interactionsSignal.CohElastic, interactionsSignal.ElectronScattering, interactionsSignal.IMDAnnihilation, interactionsSignal.InverseBetaDecay, interactionsSignal.GlashowResonance, interactionsSignal.AMNuGamma, interactionsSignal.MEC, interactionsSignal.Diffractive, interactionsSignal.EM, interactionsSignal.WeakMix);
    printf("Nuance Offsets:\nNumber with just Offset = %i\nCCQE = %i\nNCQE = %i\nNuanceRes = %i\nCCDIS = %i\nNCDis = %i\nNuEElastic = %i\n", interactionsSignal.NuanceOffset, interactionsSignal.CCQE, interactionsSignal.NCQE, interactionsSignal.NuanceRes, interactionsSignal.CCDis, interactionsSignal.NCDis, interactionsSignal.NuEElastic);
    printf("\nInteractions in BNB Sample:\nUnknown = %i\nQE = %i\nRes = %i\nDIS = %i\nCoh = %i\nCoh Elastic = %i\nElastic Scattering = %i\nIMD Annihilation = %i\nInverse Beta Decay = %i\nGlashow Resonance = %i\nAM Nu Gamma = %i\nMEC = %i\nDiffractive = %i\nEM = %i\nWeak Mix = %i\n\n", interactionsBNB.Unknown, interactionsBNB.QE, interactionsBNB.Res, interactionsBNB.DIS, interactionsBNB.Coh, interactionsBNB.CohElastic, interactionsBNB.ElectronScattering, interactionsBNB.IMDAnnihilation, interactionsBNB.InverseBetaDecay, interactionsBNB.GlashowResonance, interactionsBNB.AMNuGamma, interactionsBNB.MEC, interactionsBNB.Diffractive, interactionsBNB.EM, interactionsBNB.WeakMix);
    printf("Nuance Offsets:\nNumber with just Offset = %i\nCCQE = %i\nNCQE = %i\nNuanceRes = %i\nCCDIS = %i\nNCDis = %i\nNuEElastic = %i\n\n", interactionsBNB.NuanceOffset, interactionsBNB.CCQE, interactionsBNB.NCQE, interactionsBNB.NuanceRes, interactionsBNB.CCDis, interactionsBNB.NCDis, interactionsBNB.NuEElastic);

    //std::cout << "cosmicSpillsSum = " << cosmicSpillsSum << std::endl;

    //std::cout << "Total POT Signal: " << totalPOTSignal << std::endl;
    //std::cout << "Total POT BNB: " << totalPOTBNB << std::endl;
    //std::cout << "Total POT Cosmics: " << cosmicsPOT << std::endl;
    //printf("\nWeights:\nSignal = %f, BNB = %f, Cosmics = %f\n\n", signalWeightCurrent, BNBWeightCurrent, cosmicsWeightCurrent, signalWeightUboone, BNBWeightUboone, cosmicsWeightUboone, signalWeightNue, BNBWeightNue, cosmicsWeightNue);
   
    //std::cout << "Nu+E before weighting = " << numEventsTotal.ubooneSignal << std::endl;
    //std::cout << "BNB before weighting = " << numEventsTotal.ubooneBNB << std::endl;
    //std::cout << "Intime cosmics before weighting = " << numEventsTotal.ubooneCosmics << std::endl; 
    //std::cout << "Number of Nu+E Elastic Scattering events = " << numEventsTotal.ubooneSignal*signalWeight << std::endl;
    //std::cout << "Number of BNB events = " << numEventsTotal.ubooneBNB*BNBWeight << std::endl;
    //std::cout << "Number of Intime Cosmic events = " << numEventsTotal.ubooneCosmics*cosmicsWeight << std::endl;
    //std::cout << "" << std::endl;
    //double totalNumEvents = (numEventsTotal.ubooneCosmics*cosmicsWeight) + (numEventsTotal.ubooneBNB*BNBWeight) + (numEventsTotal.ubooneSignal*signalWeight);
    //std::cout << "Total number of events: " << totalNumEvents << std::endl;
    //std::cout << "Nu+E Elastic Scattering Events = " << 100*(numEventsTotal.ubooneSignal*signalWeight)/totalNumEvents << "%" << std::endl;
    //std::cout << "BNB Events = " << 100*(numEventsTotal.ubooneBNB*BNBWeight)/totalNumEvents << "%" << std::endl;
    //std::cout << "Intime Cosmics = " << 100*(numEventsTotal.ubooneCosmics*cosmicsWeight)/totalNumEvents << "%" << std::endl;
}
