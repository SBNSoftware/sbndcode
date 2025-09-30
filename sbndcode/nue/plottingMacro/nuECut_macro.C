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
    int keptSignal = 0;
    int cutSignal = 0;
    int keptBNB = 0;
    int cutBNB = 0;
    int keptCosmics = 0;
    int cutCosmics = 0;
};

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* keptSignal;
    TH1F* cutSignal;
    TH1F* keptBNB;
    TH1F* cutBNB;
    TH1F* keptCosmics;
    TH1F* cutCosmics;
} histGroup;

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* Signal;
    TH1F* BNB;
    TH1F* Cosmics;
} cutHistGroup;

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
        (TH1F*) base->Clone((baseName + "_cutSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_keptSignal").c_str()),
        (TH1F*) base->Clone((baseName + "_cutBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_keptBNB").c_str()),
        (TH1F*) base->Clone((baseName + "_cutCosmics").c_str()),
        (TH1F*) base->Clone((baseName + "_keptCosmics").c_str()),
    };    
}

cutHistGroup createCutHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);
    
    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_Signal").c_str()),
        (TH1F*) base->Clone((baseName + "_BNB").c_str()),
        (TH1F*) base->Clone((baseName + "_Cosmics").c_str()),
    };    
}

void styleDraw(histGroup hist, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, TPaveText* pt = nullptr, int* weighted = nullptr, int* drawLine = nullptr, int* linePos = nullptr, int* log = nullptr){
    hist.canvas->cd();
    hist.canvas->SetTickx();
    hist.canvas->SetTicky();

    if(log && *log){
        gPad->SetLogy();
        hist.keptSignal->SetMinimum(0.0000001);
        hist.keptBNB->SetMinimum(0.0000001);
        hist.keptCosmics->SetMinimum(0.0000001);
        hist.cutSignal->SetMinimum(0.0000001);
        hist.cutBNB->SetMinimum(0.0000001);
        hist.cutCosmics->SetMinimum(0.0000001);
    }

    gPad->Update();
    hist.keptSignal->SetLineWidth(2);
    hist.keptSignal->SetLineColor(kBlue+1);
    hist.keptBNB->SetLineWidth(2);
    hist.keptBNB->SetLineColor(kOrange+7);
    hist.keptCosmics->SetLineWidth(2);
    hist.keptCosmics->SetLineColor(kPink+9);
    hist.cutSignal->SetLineWidth(2);
    hist.cutSignal->SetLineColor(kBlue-7);
    hist.cutBNB->SetLineWidth(2);
    hist.cutBNB->SetLineColor(kOrange+6);
    hist.cutCosmics->SetLineWidth(2);
    hist.cutCosmics->SetLineColor(kPink+1);

    if((ymin != 999) && (ymax != 999)) hist.keptSignal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) hist.keptSignal->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({hist.keptSignal->GetMaximum(), hist.cutSignal->GetMaximum(), hist.keptBNB->GetMaximum(), hist.cutBNB->GetMaximum(), hist.keptCosmics->GetMaximum(), hist.cutCosmics->GetMaximum()});

    if((ymin == 999) && (ymax == 999)){
        double yminVal = (log && *log) ? 0.1 : 0;
        double ymaxVal  = (log && *log) ? maxYValue*10000 : maxYValue*1.1;
        
        hist.keptSignal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    hist.keptSignal->Draw("hist");
    hist.keptBNB->Draw("histsame");
    hist.keptCosmics->Draw("histsame");
    hist.cutSignal->Draw("histsame");
    hist.cutBNB->Draw("histsame");
    hist.cutCosmics->Draw("histsame");

    hist.keptSignal->SetStats(0);
    hist.keptSignal->GetXaxis()->SetTickLength(0.04);
    hist.keptSignal->GetYaxis()->SetTickLength(0.03);
    hist.keptSignal->GetXaxis()->SetTickSize(0.02);
    hist.keptSignal->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(hist.keptSignal, "Kept Nu+E + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(hist.keptBNB, "Kept BNB + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(hist.keptCosmics, "Kept Intime Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(hist.cutSignal, "Cut Nu+E + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(hist.cutBNB, "Cut BNB + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(hist.cutCosmics, "Cut Intime Cosmics, Pandora BDT Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, hist.keptSignal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, hist.keptSignal->GetMaximum() * 0.93, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.1, hist.keptSignal->GetMaximum() * 0.93, "2m_{e}");
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

void styleDrawCut(TCanvas* canvas, TH1F* signal, TH1F* BNB, TH1F* cosmics, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, TPaveText* pt = nullptr, int* weighted = nullptr, int* drawLine = nullptr, int* linePos = nullptr, int* log = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    if(log && *log){
        gPad->SetLogy();
        signal->SetMinimum(0.0000001);
        BNB->SetMinimum(0.0000001);
        cosmics->SetMinimum(0.0000001);
    }

    gPad->Update();
    signal->SetLineWidth(2);
    signal->SetLineColor(kBlue+1);
    BNB->SetLineWidth(2);
    BNB->SetLineColor(kOrange+7);
    cosmics->SetLineWidth(2);
    cosmics->SetLineColor(kPink+9);

    if((ymin != 999) && (ymax != 999)) signal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) signal->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({signal->GetMaximum(), BNB->GetMaximum(), cosmics->GetMaximum()});

    if((ymin == 999) && (ymax == 999)){
        double yminVal = (log && *log) ? 0.1 : 0;
        double ymaxVal  = (log && *log) ? maxYValue*10000 : maxYValue*1.1;
        
        signal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    signal->Draw("hist");
    BNB->Draw("histsame");
    cosmics->Draw("histsame");

    signal->SetStats(0);
    signal->GetXaxis()->SetTickLength(0.04);
    signal->GetYaxis()->SetTickLength(0.03);
    signal->GetXaxis()->SetTickSize(0.02);
    signal->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(signal, "Nu+E + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(BNB, "BNB + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(cosmics, "Intime Cosmics, Pandora BDT Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();


    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, signal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, signal->GetMaximum() * 0.93, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.1, signal->GetMaximum() * 0.93, "2m_{e}");
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

void styleDrawCutIndividual(TCanvas* canvas, TH1F* hist, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, TPaveText* pt = nullptr, int* weighted = nullptr, int* drawLine = nullptr, int* linePos = nullptr, int* log = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    if(log && *log){
        gPad->SetLogy();
        hist->SetMinimum(0.0000001);
    }

    gPad->Update();
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue+1);

    if((ymin != 999) && (ymax != 999)) hist->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) hist->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({hist->GetMaximum()});

    if((ymin == 999) && (ymax == 999)){
        double yminVal = (log && *log) ? 0.1 : 0;
        double ymaxVal  = (log && *log) ? maxYValue*10000 : maxYValue*1.1;
        
        hist->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    hist->Draw("hist");

    hist->SetStats(0);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetYaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTickSize(0.02);
    hist->GetYaxis()->SetTickSize(0.02);

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, hist->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, hist->GetMaximum() * 0.93, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.1, hist->GetMaximum() * 0.93, "2m_{e}");
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

void efficiencyCut(cutHistGroup hists, double ymin, double ymax, double xmin, double xmax, const char* filenameBase, const std::string& legendLocation, double signalWeight, double BNBWeight, double cosmicsWeight, int* drawLine = nullptr, int* linePos = nullptr, std::string xlabel = "", double efficiencyWay = 0.0, double* currentMaxEffPurBin = nullptr, double* ubooneMaxEffPurBin = nullptr){

    TCanvas *efficiencyCanvas = new TCanvas("efficiency_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* signalEff = (TH1F*) hists.Signal->Clone("eff hist");
    signalEff->Reset();
    signalEff->GetYaxis()->SetTitle("Efficiency"); 
    signalEff->GetXaxis()->SetTitle(xlabel.c_str());
    TH1F* BNBEff = (TH1F*) hists.BNB->Clone("eff hist");
    BNBEff->Reset();
    BNBEff->GetYaxis()->SetTitle("Efficiency"); 
    BNBEff->GetXaxis()->SetTitle(xlabel.c_str());
    TH1F* cosmicsEff = (TH1F*) hists.Cosmics->Clone("eff hist");
    cosmicsEff->Reset();
    cosmicsEff->GetYaxis()->SetTitle("Efficiency"); 
    cosmicsEff->GetXaxis()->SetTitle(xlabel.c_str());

    TCanvas *rejectionCanvas = new TCanvas("rejection_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* signalRej = (TH1F*) hists.Signal->Clone("rej hist");
    signalRej->Reset();
    signalRej->GetYaxis()->SetTitle("Rejection"); 
    signalRej->GetXaxis()->SetTitle(xlabel.c_str());
    TH1F* BNBRej = (TH1F*) hists.BNB->Clone("rej hist");
    BNBRej->Reset();
    BNBRej->GetYaxis()->SetTitle("Rejection"); 
    BNBRej->GetXaxis()->SetTitle(xlabel.c_str());
    TH1F* cosmicsRej = (TH1F*) hists.Cosmics->Clone("rej hist");
    cosmicsRej->Reset();
    cosmicsRej->GetYaxis()->SetTitle("Rejection"); 
    cosmicsRej->GetXaxis()->SetTitle(xlabel.c_str());

    TCanvas *purityCanvas = new TCanvas("purity_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* pur = (TH1F*) hists.Signal->Clone("pur hist");
    pur->Reset();
    pur->GetYaxis()->SetTitle("Purity"); 
    pur->GetXaxis()->SetTitle(xlabel.c_str());

    TCanvas *effpurCanvas = new TCanvas("effpur_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* effPur = (TH1F*) hists.Signal->Clone("eff pur hist");
    effPur->Reset();
    effPur->GetYaxis()->SetTitle("Efficiency x Purity");
    effPur->GetXaxis()->SetTitle(xlabel.c_str());    

    int numBins = hists.Signal->GetNbinsX();
    double signalSum = 0.0;
    double BNBSum = 0.0;
    double cosmicsSum = 0.0;

    double signalTotal = 0.0;
    double BNBTotal = 0.0;
    double cosmicsTotal = 0.0;

    // efficiencyWay == -1 includes everything to the right of the cut
    if(efficiencyWay == -1){
        for(int i = 1; i <= numBins; ++i){
            signalTotal += hists.Signal->GetBinContent(i);
            BNBTotal += hists.BNB->GetBinContent(i);
            cosmicsTotal += hists.Cosmics->GetBinContent(i);
        }

        std::cout << "TOTALS:" << std::endl;
        std::cout << "Signal = " << signalTotal << std::endl;
        std::cout << "BNB = " << BNBTotal << std::endl;
        std::cout << "Cosmics = " << cosmicsTotal << std::endl;

        double sizeSignal = hists.Signal->GetEntries();
        double sizeBNB = hists.BNB->GetEntries();
        double sizeCosmics = hists.Cosmics->GetEntries();

        for(int i = 1; i <= numBins; ++i){
            signalSum += hists.Signal->GetBinContent(i);
            BNBSum += hists.BNB->GetBinContent(i);
            cosmicsSum += hists.Cosmics->GetBinContent(i);

            double signalEffValue = (signalTotal - signalSum)/sizeSignal;
            double BNBEffValue = (BNBTotal - BNBSum)/sizeBNB;
            double cosmicsEffValue = (cosmicsTotal - cosmicsSum)/sizeCosmics;

            signalEff->SetBinContent(i, signalEffValue);
            BNBEff->SetBinContent(i, BNBEffValue);
            cosmicsEff->SetBinContent(i, cosmicsEffValue);

            signalRej->SetBinContent(i, 1-signalEffValue);
            BNBRej->SetBinContent(i, 1-BNBEffValue);
            cosmicsRej->SetBinContent(i, 1-cosmicsEffValue);

            double sumBackground = ((BNBTotal - BNBSum) * BNBWeight) + ((cosmicsTotal - cosmicsSum) * cosmicsWeight);
            double totalSum = (sumBackground + (signalTotal - signalSum));

            double purityValue;
            if(signalSum != 0 && sumBackground != 0){
                purityValue = ((signalTotal - signalSum) * signalWeight)/ totalSum;
            } else{
                purityValue = 0;
            }

            double effPurValue = (signalEffValue * purityValue);  
            
            pur->SetBinContent(i, purityValue);
            effPur->SetBinContent(i, effPurValue);
        }   
    }


    // efficiencyWay == 1 includes everything to the left of the cut
    if(efficiencyWay == 1){
        for(int i = 1; i <= numBins; ++i){
            signalSum += hists.Signal->GetBinContent(i);
            BNBSum += hists.BNB->GetBinContent(i);
            cosmicsSum += hists.Cosmics->GetBinContent(i);

            double signalEffValue = signalSum/sizeSignal;
            double BNBEffValue = BNBSum/sizeBNB;
            double cosmicsEffValue = cosmicsSum/sizeCosmics;
            BNBEff-SetBinContent(i, currentBNBEffValue);
            cosmicsEff->SetBinContent(i, currentCosmicsEffValue);
        
            signalRej->SetBinContent(i, 1-signalEffValue);
            BNBRej->SetBinContent(i, 1-BNBEffValue);
            cosmicsRej->SetBinContent(i, 1-cosmicsEffValue);

            double sumBackground = (BNBSum * BNBWeight) + (cosmicsSum * cosmicsWeight);
            double totalSum = (sumBackground + signalSum);

            double purityValue;
        
            if(signalSum != 0 && sumBackground != 0){
                purityValue = (signalSum * signalWeight)/ totalSum;
            } else{
                purityValue = 0;
            }

            double effPurValue = (signalEffValue * purityValue);
            pur->SetBinContent(i, purityValue);
            effPur->SetBinContent(i, effPurValue);
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
    pt->AddText(Form("Number of Signal Entries: %d", (int)hists.Signal->GetEntries()));
    pt->AddText(Form("Number of BNB Entries: %d", (int)hists.BNB->GetEntries()));
    pt->AddText(Form("Number of Cosmics Entries: %d", (int)hists.Cosmics->GetEntries()));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0);

    int funcValue = 1;
    int log = 1;

    pur->GetYaxis()->SetNoExponent(false);
    pur->GetYaxis()->SetMaxDigits(2);
    
    effPur->GetYaxis()->SetNoExponent(false);
    effPur->GetYaxis()->SetMaxDigits(2);

    TCanvas *rejBNBCanvas = new TCanvas("rejBNB_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TCanvas *rejCosmicCanvas = new TCanvas("rejCosmic_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TCanvas *effSignalCanvas = new TCanvas("effSignal_canvas", "Graph Draw Options", 200, 10, 600, 400);

    styleDrawCut(efficiencyCanvas, signalEff, BNBEff, cosmicsEff, ymin, ymax, xmin, xmax, (std::string(filenameBase) + "_eff.pdf").c_str(), legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawCut(rejectionCanvas, signalRej, BNBRej, cosmicsRej, ymin, ymax, xmin, xmax, (std::string(filenameBase) + "_rej.pdf").c_str(), legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawCutIndividual(purityCanvas, pur, ymin, ymax, xmin, xmax, (std::string(filenameBase) + "_pur.pdf").c_str(), legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    styleDrawCutIndividual(effpurCanvas, effPur, ymin, ymax, xmin, xmax, (std::string(filenameBase) + "_effpur.pdf").c_str(), legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
}

void weighted(){}

void weightedCut(cutHistGroup hists, double signalWeight, double BNBWeight, double cosmicsWeight, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    TCanvas *weightedCanvas = new TCanvas("weighted_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* signalWeighted = (TH1F*) hists.Signal->Clone("weighted hist");
    signalWeighted->Scale(signalWeight);
    signalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    
    TH1F* BNBWeighted = (TH1F*) hists.BNB->Clone("weighted hist");
    BNBWeighted->Scale(BNBWeight);
    BNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 

    TH1F* cosmicsWeighted = (TH1F*) hists.cosmics->Clone("weighted hist");
    cosmicsWeighted->Scale(cosmicsWeight);
    cosmicsWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    
    int funcValue = 1;
    int log = 1;

    styleDrawCut(weightedCanvas, signalWeighted, BNBWeighted, cosmicsWeighted, ymin, ymax, xmin, xmax, filename, legendLocation, nullptr, &funcValue, drawLine, linePos, &log);
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

void nuECut_macro(){
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_22Sep.root");
    std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeights/";

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

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignal;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNB;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsCosmics;

    double totalPOTSignal = 0;
    double totalPOTBNB = 0;
    double totalPOTCosmics = 0;

    double cosmicSpillsSum = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 0) cosmicSpillsSum += subRunNumGenEvents;

        if(subRunSignal == 1){
            std::pair<unsigned int, unsigned int> keySignal = std::make_pair(subRunRun, subRunNumber);
            if(subRunDLCurrent == 0 && seenSubRunsSignal.find(keySignal) == seenSubRunsSignal.end()){
                totalPOTSignal += subRunPOT;
                seenSubRunsSignal.insert(keySignal);
            }
        } else if(subRunSignal == 2){
            std::pair<unsigned int, unsigned int> keyBNB = std::make_pair(subRunRun, subRunNumber);
            if(subRunDLCurrent == 0 && seenSubRunsBNB.find(keyBNB) == seenSubRunsBNB.end()){
                totalPOTBNB += subRunPOT;
                seenSubRunsBNB.insert(keyBNB);
            } else if(subRunDLCurrent == 0){
                totalPOTBNB += subRunPOT;
            }
        } else if(subRunSignal == 3){
            std::pair<unsigned int, unsigned int> keyCosmics = std::make_pair(subRunRun, subRunNumber);
            if(subRunDLCurrent == 0 && seenSubRunsCosmics.find(keyCosmics) == seenSubRunsCosmics.end()){
                totalPOTCosmics += subRunPOT;
                seenSubRunsCosmics.insert(keyCosmics);
            }
        } 
    }

    //double numberFiles = 600;
    //double cosmicSpills = numberFiles * 500;
    //std::cout << "cosmicSpillsSum = " << cosmicSpillsSum << std::endl;
    double cosmicsPOT = cosmicSpillsSum * 5e12;

    double signalWeight = totalPOTSignal/totalPOTSignal;
    double BNBWeight = totalPOTSignal/totalPOTBNB;
    double cosmicsWeight = totalPOTSignal/cosmicsPOT;
    
    //std::cout << "Total POT Signal: " << totalPOTSignal << std::endl;
    //std::cout << "Total POT BNB: " << totalPOTBNB << std::endl;
    //std::cout << "Total POT Cosmics: " << cosmicsPOT << std::endl;

    //printf("Weights:\nSignal = %f, BNB = %f, Cosmics = %f\n", signalWeight, BNBWeight, cosmicsWeight);

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

    auto numSlices = createHistGroup("numSlices", "Number of Slices in an Event", "Number of Slices", 45, 0, 45);
    auto numSlicesCut = createCutHistGroup("numSlicesCut", "Number of Slices in an Event (After Cuts)", "Number of Slices", 45, 0, 45);

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

        std::vector<recoParticle> recoParticlesInEvent;
        std::vector<recoNeutrino> recoNeutrinosInEvent;
        std::vector<recoSlice> recoSlicesInEvent;
        std::vector<trueParticle> trueParticlesInEvent;
        
        recoParticle chosenRecoParticleCRUMBS;
        recoNeutrino chosenRecoNeutrinoCRUMBS;
        recoSlice chosenRecoSliceCRUMBS;
        trueNeutrino chosenTrueNeutrino;
        
        double energyInSlice = 0;
        double numPFPsInEvent = 0;
        double numPFPsInChosenSlice = 0;

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
                if(DLCurrent == 2){
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

        if(slice == 0) continue; 

        chosenRecoSliceCRUMBS = chooseSlice(recoSlicesInEvent, 1);

        int neutrino = 0;
        int numRecoNeutrinosInEvent = 0;
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

                //printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, recoNeutrino.pdg, recoNeutrino.isPrimary, recoNeutrino.vx, recoNeutrino.vy, recoNeutrino.vz, recoNeutrino.sliceID);
            }
        }

        if(neutrino == 0) continue;

        chosenRecoNeutrinoCRUMBS = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCRUMBS.id); 

        if(signal != 3){
            double deltaXCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vx - chosenTrueNeutrino.vx);
            double deltaYCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vy - chosenTrueNeutrino.vy);
            double deltaZCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vz - chosenTrueNeutrino.vz);
            double deltaRCRUMBSValue = std::sqrt((deltaXCRUMBSValue * deltaXCRUMBSValue) + (deltaYCRUMBSValue * deltaYCRUMBSValue) + (deltaZCRUMBSValue * deltaZCRUMBSValue));
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


        // Cut = 0, Keep = 1
        double cutOrKeep = 0;

        // CUT VALUES
        double FV_xLow = -201.05;
        double FV_xHigh = 200.55;
        double FV_yLow = -200.049;
        double FV_yHigh = 199.549;
        double FV_zLow = 0.249951;
        double FV_zHigh = 490.154;

        if((chosenRecoNeutrinoCRUMBS.vx > FV_xLow && chosenRecoNeutrinoCRUMBS.vx < FV_xHigh) && (chosenRecoNeutrinoCRUMBS.vy > FV_yLow && chosenRecoNeutrinoCRUMBS.vy < FV_yHigh) && (chosenRecoNeutrinoCRUMBS.vz > FV_zLow && chosenRecoNeutrinoCRUMBS.vz < FV_zHigh)) cutOrKeep = 1;

        if(cutOrKeep == 0){
            // Fails cuts
            if(DLCurrent == 2 && signal == 1){
                numSlices.cutSignal->Fill(numSlicesInEvent);
            } else if(DLCurrent == 2 && signal == 2){
                numSlices.cutBNB->Fill(numSlicesInEvent);
            } else if(DLCurrent == 2 && signal == 3){
                numSlices.cutCosmics->Fill(numSlicesInEvent);
            }

        } else if(cutOrKeep == 1){
            // Passes cuts
            if(DLCurrent == 2 && signal == 1){
                numSlices.keptSignal->Fill(numSlicesInEvent);
                numSlicesCut.Signal->Fill(numSlicesInEvent);
            } else if(DLCurrent == 2 && signal == 2){
                numSlices.keptBNB->Fill(numSlicesInEvent);
                numSlicesCut.BNB->Fill(numSlicesInEvent);
            } else if(DLCurrent == 2 && signal == 3){
                numSlices.keptCosmics->Fill(numSlicesInEvent);
                numSlicesCut.Cosmics->Fill(numSlicesInEvent);
            }
        }

    }
   
    int drawLine = 1;
    int left = 0;
    int right = 1;

    styleDrawCut(numSlicesCut.canvas, numSlicesCut.Signal, numSlicesCut.BNB, numSlicesCut.Cosmics, 999, 999, 999, 999, (base_path + "numSlicesCut_dist.pdf").c_str(), legendLocation, pt = nullptr, &funcValue, drawLine, linePos);
    weightedCut();
    efficiencyCut(numSlicesCut, 999, 999, 999, 999, (base_path + "numSlicesCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Number of Slices", 1);

}
