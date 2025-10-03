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

void styleDraw(TCanvas* canvas, TH1F* keptSignal, TH1F* cutSignal, TH1F* keptBNB, TH1F* cutBNB, TH1F* keptCosmics, TH1F* cutCosmics, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, TPaveText* pt = nullptr, int* weighted = nullptr, int* drawLine = nullptr, int* linePos = nullptr, int* log = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    if(log && *log){
        gPad->SetLogy();
        keptSignal->SetMinimum(0.0000001);
        keptBNB->SetMinimum(0.0000001);
        keptCosmics->SetMinimum(0.0000001);
        cutSignal->SetMinimum(0.0000001);
        cutBNB->SetMinimum(0.0000001);
        cutCosmics->SetMinimum(0.0000001);
    }

    gPad->Update();
    keptSignal->SetLineWidth(2);
    keptSignal->SetLineColor(kBlue+1);
    keptBNB->SetLineWidth(2);
    keptBNB->SetLineColor(kOrange+7);
    keptCosmics->SetLineWidth(2);
    keptCosmics->SetLineColor(kPink+9);
    cutSignal->SetLineWidth(2);
    cutSignal->SetLineColor(kBlue-7);
    cutBNB->SetLineWidth(2);
    cutBNB->SetLineColor(kOrange+6);
    cutCosmics->SetLineWidth(2);
    cutCosmics->SetLineColor(kPink+1);

    if((ymin != 999) && (ymax != 999)) keptSignal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) keptSignal->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({keptSignal->GetBinContent(keptSignal->GetMaximumBin()), cutSignal->GetBinContent(cutSignal->GetMaximumBin()), keptBNB->GetBinContent(keptBNB->GetMaximumBin()), cutBNB->GetBinContent(cutBNB->GetMaximumBin()), keptCosmics->GetBinContent(keptCosmics->GetMaximumBin()), cutCosmics->GetBinContent(cutCosmics->GetMaximumBin())});
    std::cout << "maxYValue = " << maxYValue << std::endl;
    std::cout << "Signal Max = " << keptSignal->GetBinContent(keptSignal->GetMaximumBin()) << ", " << cutSignal->GetBinContent(cutSignal->GetMaximumBin()) << ", BNB Max = " << keptBNB->GetBinContent(keptBNB->GetMaximumBin()) << ", " << cutBNB->GetBinContent(cutBNB->GetMaximumBin()) << ", Cosmics Max = " << keptCosmics->GetBinContent(keptCosmics->GetMaximumBin()) << ", " << cutCosmics->GetBinContent(cutCosmics->GetMaximumBin()) << std::endl;

    if((ymin == 999) && (ymax == 999)){
        double yminVal = (log && *log) ? 0.1 : 0;
        double ymaxVal  = (log && *log) ? maxYValue*10000 : maxYValue*1.1;
        
        keptSignal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    keptSignal->Draw("hist");
    keptBNB->Draw("histsame");
    keptCosmics->Draw("histsame");
    cutSignal->Draw("histsame");
    cutBNB->Draw("histsame");
    cutCosmics->Draw("histsame");

    keptSignal->SetStats(0);
    keptSignal->GetXaxis()->SetTickLength(0.04);
    keptSignal->GetYaxis()->SetTickLength(0.03);
    keptSignal->GetXaxis()->SetTickSize(0.02);
    keptSignal->GetYaxis()->SetTickSize(0.02);

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
    legend->AddEntry(keptSignal, "Kept Nu+E + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(cutSignal, "Cut Nu+E + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(keptBNB, "Kept BNB + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(cutBNB, "Cut BNB + Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(keptCosmics, "Kept Intime Cosmics, Pandora BDT Vertexing", "f");
    legend->AddEntry(cutCosmics, "Cut Intime Cosmics, Pandora BDT Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, keptSignal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, keptSignal->GetMaximum() * 0.93, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.1, keptSignal->GetMaximum() * 0.93, "2m_{e}");
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

    double maxYValue = std::max({signal->GetBinContent(signal->GetMaximumBin()), BNB->GetBinContent(BNB->GetMaximumBin()), cosmics->GetBinContent(cosmics->GetMaximumBin())});
    std::cout << "maxYValue = " << maxYValue << std::endl;
    std::cout << "Signal Max = " << signal->GetBinContent(signal->GetMaximumBin()) << ", BNB Max = " << BNB->GetBinContent(BNB->GetMaximumBin()) << ", Cosmics Max = " << cosmics->GetBinContent(cosmics->GetMaximumBin()) << std::endl;

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


    if(legendLocation == "topRight"){
        Lxmin = 0.548;
        Lymax = 0.863;
        Lxmax = 0.87;
        Lymin = 0.7515;
    } else if(legendLocation == "topLeft"){
        Lxmin = 0.13;
        Lymax = 0.863;
        Lxmax = 0.452;
        Lymin = 0.7515;
    } else if(legendLocation == "bottomRight"){
        Lxmin = 0.548;
        Lymax = 0.2485;
        Lxmax = 0.87;
        Lymin = 0.137;
    } else if(legendLocation == "bottomLeft"){
        Lxmin = 0.13;
        Lymax = 0.2485;
        Lxmax = 0.452;
        Lymin = 0.137;
    }

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

    double maxYValue = std::max({hist->GetBinContent(hist->GetMaximumBin())});

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

    //if(pt){
        //pt->SetTextSize(legend->GetTextSize());
        //pt->SetTextFont(legend->GetTextFont());
        //pt->Draw();
    //}

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
            std::cout << "Purity = " << purityValue << ", Efficiency = " << signalEffValue << ", Purity x Efficiency = " << effPurValue << std::endl;
        }   
    }


    // efficiencyWay == 1 includes everything to the left of the cut
    if(efficiencyWay == 1){
        for(int i = 1; i <= numBins; ++i){
            signalSum += hists.Signal->GetBinContent(i);
            BNBSum += hists.BNB->GetBinContent(i);
            cosmicsSum += hists.Cosmics->GetBinContent(i);

            double sizeSignal = hists.Signal->GetEntries();
            double sizeBNB = hists.BNB->GetEntries();
            double sizeCosmics = hists.Cosmics->GetEntries();

            double signalEffValue = signalSum/sizeSignal;
            double BNBEffValue = BNBSum/sizeBNB;
            double cosmicsEffValue = cosmicsSum/sizeCosmics;
            
            signalEff->SetBinContent(i, signalEffValue);
            BNBEff->SetBinContent(i, BNBEffValue);
            cosmicsEff->SetBinContent(i, cosmicsEffValue);
        
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

            std::cout << "Purity = " << purityValue << ", Efficiency = " << signalEffValue << ", Purity x Efficiency = " << effPurValue << std::endl;
        }
    }

    double Lxmin = 0;
    double Lymax = 0;
    double Lxmax = 0;
    double Lymin = 0;

    if(legendLocation == "topRight"){
        Lxmin = 0.548;
        Lymax = 0.863;
        Lxmax = 0.87;
        Lymin = 0.7515;
    } else if(legendLocation == "topLeft"){
        Lxmin = 0.13;
        Lymax = 0.863;
        Lxmax = 0.452;
        Lymin = 0.7515;
    } else if(legendLocation == "bottomRight"){
        Lxmin = 0.548;
        Lymax = 0.2485;
        Lxmax = 0.87;
        Lymin = 0.137;
    } else if(legendLocation == "bottomLeft"){
        Lxmin = 0.13;
        Lymax = 0.2485;
        Lxmax = 0.452;
        Lymin = 0.137;
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

void weighted(histGroup hists, double signalWeight, double BNBWeight, double cosmicsWeight, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    TCanvas *weightedCanvas = new TCanvas("weighted_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* keptSignalWeighted = (TH1F*) hists.keptSignal->Clone("weighted hist");
    keptSignalWeighted->Scale(signalWeight);
    keptSignalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    TH1F* cutSignalWeighted = (TH1F*) hists.cutSignal->Clone("weighted hist");
    cutSignalWeighted->Scale(signalWeight);
    cutSignalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    
    TH1F* keptBNBWeighted = (TH1F*) hists.keptBNB->Clone("weighted hist");
    keptBNBWeighted->Scale(BNBWeight);
    keptBNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    TH1F* cutBNBWeighted = (TH1F*) hists.cutBNB->Clone("weighted hist");
    cutBNBWeighted->Scale(BNBWeight);
    cutBNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 

    TH1F* keptCosmicsWeighted = (TH1F*) hists.keptCosmics->Clone("weighted hist");
    keptCosmicsWeighted->Scale(cosmicsWeight);
    keptCosmicsWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    TH1F* cutCosmicsWeighted = (TH1F*) hists.cutCosmics->Clone("weighted hist");
    cutCosmicsWeighted->Scale(cosmicsWeight);
    cutCosmicsWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    
    int funcValue = 1;
    int log = 1;

    styleDraw(weightedCanvas, keptSignalWeighted, cutSignalWeighted, keptBNBWeighted, cutBNBWeighted, keptCosmicsWeighted, cutCosmicsWeighted, ymin, ymax, xmin, xmax, filename, legendLocation, nullptr, &funcValue, drawLine, linePos, &log);
}

void weightedCut(cutHistGroup hists, double signalWeight, double BNBWeight, double cosmicsWeight, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    TCanvas *weightedCanvas = new TCanvas("weighted_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* signalWeighted = (TH1F*) hists.Signal->Clone("weighted hist");
    signalWeighted->Scale(signalWeight);
    signalWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 
    
    TH1F* BNBWeighted = (TH1F*) hists.BNB->Clone("weighted hist");
    BNBWeighted->Scale(BNBWeight);
    BNBWeighted->GetYaxis()->SetTitle("Number of Events (POT Weighted)"); 

    TH1F* cosmicsWeighted = (TH1F*) hists.Cosmics->Clone("weighted hist");
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
    std::string base_path = "/nashome/c/coackley/nuECutsPlots/";

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

    auto numRecoNeutrinos = createHistGroup("numRecoNeutrinos", "Number of Reco Neutrinos in an Event", "Number of Reco Neutrinos", 10, 0, 10);
    auto numRecoNeutrinosCut = createCutHistGroup("numRecoNeutrinosCut", "Number of Reco Neutrinos in an Event (After Cuts)", "Number of Reco Neutrinos", 10, 0, 10);

    auto sliceCompletenessCRUMBS = createHistGroup("sliceCompletenessCRUMBS", "Completeness of the Slice with the Highest CRUMBS Score", "Completeness", 102, 0, 1.02); // Bin width = 0.005
    auto sliceCompletenessCRUMBSCut = createCutHistGroup("sliceCompletenessCRUMBSCut", "Completeness of the Slice with the Highest CRUMBS Score (After Cuts)", "Completeness", 102, 0, 1.02); // Bin width = 0.005
    auto sliceScoreCRUMBS = createHistGroup("sliceScoreCRUMBS", "CRUMBS Score of the Slice with the Highest CRUMBS Score", "CRUMBS Score", 25, -1, 1); 
    auto sliceScoreCRUMBSCut = createCutHistGroup("sliceScoreCRUMBSCut", "CRUMBS Score of the Slice with the Highest CRUMBS Score (After Cuts)", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCRUMBS = createHistGroup("slicePurityCRUMBS", "Purity of the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto slicePurityCRUMBSCut = createCutHistGroup("slicePurityCRUMBSCut", "Purity of the Slice with the Highest CRUMBS Score (After Cuts)", "Purity", 50, 0, 1.02);
    auto highestPFPCompletenessCRUMBS = createHistGroup("highestPFPCompletenessCRUMBS", "Completeness of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "Completeness", 50, 0, 1.02);
    auto highestPFPCompletenessCRUMBSCut = createCutHistGroup("highestPFPCompletenessCRUMBSCut", "Completeness of the Highest Energy PFP in the Slice with the Highest CRUMBS Score (After Cuts)", "Completeness", 50, 0, 1.02);
    auto highestPFPPurityCRUMBS = createHistGroup("highestPFPPurityCRUMBS", "Purity of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto highestPFPPurityCRUMBSCut = createCutHistGroup("highestPFPPurityCRUMBSCut", "Purity of the Highest Energy PFP in the Slice with the Highest CRUMBS Score (After Cuts)", "Purity", 50, 0, 1.02);
    auto numPFPsCRUMBS = createHistGroup("numPFPsCRUMBS", "Number of PFPs in the Slice with the Highest CRUMBS Score", "Number of PFPs", 10, 0, 10);
    auto numPFPsCRUMBSCut = createCutHistGroup("numPFPsCRUMBSCut", "Number of PFPs in the Slice with the Highest CRUMBS Score (After Cuts)", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCRUMBS = createHistGroup("ratioChosenSummedEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenSummedEnergyCRUMBSCut = createCutHistGroup("ratioChosenSummedEnergyCRUMBSCut", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score (After Cuts)", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCRUMBS = createHistGroup("ratioChosenTrueEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioChosenTrueEnergyCRUMBSCut = createCutHistGroup("ratioChosenTrueEnergyCRUMBSCut", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score and the True Shower Energy (After Cuts)", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCRUMBS = createHistGroup("ratioSummedTrueEnergyCRUMBS", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCRUMBSCut = createCutHistGroup("ratioSummedTrueEnergyCRUMBSCut", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score and the True Shower Energy (After Cuts)", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    auto EtrueThetaRecoCRUMBS = createHistGroup("EtrueThetaRecoCRUMBS", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto EtrueThetaRecoCRUMBSCut = createCutHistGroup("EtrueThetaRecoCRUMBSCut", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score (After Cuts)", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCRUMBS = createHistGroup("ERecoSumThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaTrueCRUMBSCut = createCutHistGroup("ERecoSumThetaTrueCRUMBSCut", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score (After Cuts)", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBS = createHistGroup("ERecoHighestThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBSCut = createCutHistGroup("ERecoHighestThetaTrueCRUMBSCut", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score (After Cuts)", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCRUMBS = createHistGroup("ERecoSumThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoSumThetaRecoCRUMBSCut = createCutHistGroup("ERecoSumThetaRecoCRUMBSCut", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCRUMBS = createHistGroup("ERecoHighestThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCRUMBSCut = createCutHistGroup("ERecoHighestThetaRecoCRUMBSCut", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);

    auto primaryVertexXCRUMBS = createHistGroup("primaryVertexXCRUMBS", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "X Coordinate (cm)", 82, -201.3, 201.3);
    auto primaryVertexYCRUMBS = createHistGroup("primaryVertexYCRUMBS", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Y Coordinate (cm)", 83, -203.8, 203.8);
    auto primaryVertexZCRUMBS = createHistGroup("primaryVertexZCRUMBS", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Z Coordinate (cm)", 104, 0, 509.4);
    auto primaryVertexXCRUMBSCut = createCutHistGroup("primaryVertexXCRUMBSCut", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "X Coordinate (cm)", 82, -201.3, 201.3);
    auto primaryVertexYCRUMBSCut = createCutHistGroup("primaryVertexYCRUMBSCut", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "Y Coordinate (cm)", 83, -203.8, 203.8);
    auto primaryVertexZCRUMBSCut = createCutHistGroup("primaryVertexZCRUMBSCut", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "Z Coordinate (cm)", 104, 0, 509.4);

    auto primaryVertexXCRUMBSNegative = createHistGroup("primaryVertexXCRUMBSNegative", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "X Coordinate (cm)", 805, -201.3, 201.3);
    auto primaryVertexYCRUMBSNegative = createHistGroup("primaryVertexYCRUMBSNegative", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Y Coordinate (cm)", 815, -203.8, 203.8);
    auto primaryVertexZCRUMBSNegative = createHistGroup("primaryVertexZCRUMBSNegative", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Z Coordinate (cm)", 1019, 0, 509.4);
    auto primaryVertexXCRUMBSNegativeCut = createCutHistGroup("primaryVertexXCRUMBSNegativeCut", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "X Coordinate (cm)", 805, -201.3, 201.3);
    auto primaryVertexYCRUMBSNegativeCut = createCutHistGroup("primaryVertexYCRUMBSNegativeCut", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "Y Coordinate (cm)", 815, -203.8, 203.8);
    auto primaryVertexZCRUMBSNegativeCut = createCutHistGroup("primaryVertexZCRUMBSNegativeCut", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "Z Coordinate (cm)", 1019, 0, 509.4);

    auto primaryVertexXCRUMBSPositive = createHistGroup("primaryVertexXCRUMBSPositive", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "X Coordinate (cm)", 805, -201.3, 201.3);
    auto primaryVertexYCRUMBSPositive = createHistGroup("primaryVertexYCRUMBSPositive", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Y Coordinate (cm)", 815, -203.8, 203.8);
    auto primaryVertexZCRUMBSPositive = createHistGroup("primaryVertexZCRUMBSPositive", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score", "Z Coordinate (cm)", 1019, 0, 509.4);
    auto primaryVertexXCRUMBSPositiveCut = createCutHistGroup("primaryVertexXCRUMBSPositiveCut", "X Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "X Coordinate (cm)", 805, -201.3, 201.3);
    auto primaryVertexYCRUMBSPositiveCut = createCutHistGroup("primaryVertexYCRUMBSPositiveCut", "Y Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "Y Coordinate (cm)", 815, -203.8, 203.8);
    auto primaryVertexZCRUMBSPositiveCut = createCutHistGroup("primaryVertexZCRUMBSPositiveCut", "Z Coordinate of the Primary Neutrino Vertex in the Slice with the Highest CRUMBS Score (After Cuts)", "Z Coordinate (cm)", 1019, 0, 509.4);
 
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
    //for(Long64_t i = 0; i < 100; ++i){
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

        double sliceScoreCut = -0.3;
        double numPFPsCut = 1;

        if((chosenRecoNeutrinoCRUMBS.vx > FV_xLow && chosenRecoNeutrinoCRUMBS.vx < FV_xHigh) && (chosenRecoNeutrinoCRUMBS.vy > FV_yLow && chosenRecoNeutrinoCRUMBS.vy < FV_yHigh) && (chosenRecoNeutrinoCRUMBS.vz > FV_zLow && chosenRecoNeutrinoCRUMBS.vz < FV_zHigh)){
            if(chosenRecoSliceCRUMBS.score > sliceScoreCut){
                if(numPFPsCut >= numPFPsSliceCRUMBS){
                    cutOrKeep = 1;
                }
            }
        }

        if(cutOrKeep == 0){
            // Fails cuts
            if(DLCurrent == 2 && signal == 1){
                numSlices.cutSignal->Fill(numSlicesInEvent);
                sliceScoreCRUMBS.cutSignal->Fill(chosenRecoSliceCRUMBS.score);
                numRecoNeutrinos.cutSignal->Fill(numRecoNeutrinosInEvent);
                primaryVertexXCRUMBS.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBS.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBS.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSNegative.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSNegative.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSPositive.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSPositive.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSPositive.cutSignal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                numPFPsCRUMBS.cutSignal->Fill(numPFPsSliceCRUMBS);
                ratioChosenSummedEnergyCRUMBS.cutSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ERecoSumThetaRecoCRUMBS.cutSignal->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBS.cutSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                slicePurityCRUMBS.cutSignal->Fill(chosenSlicePurityCRUMBS);
                sliceCompletenessCRUMBS.cutSignal->Fill(chosenSliceCompletenessCRUMBS);
                highestPFPCompletenessCRUMBS.cutSignal->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPPurityCRUMBS.cutSignal->Fill(chosenRecoParticleCRUMBS.purity);
            } else if(DLCurrent == 2 && signal == 2){
                numSlices.cutBNB->Fill(numSlicesInEvent);
                sliceScoreCRUMBS.cutBNB->Fill(chosenRecoSliceCRUMBS.score);
                numRecoNeutrinos.cutBNB->Fill(numRecoNeutrinosInEvent);
                primaryVertexXCRUMBS.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBS.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBS.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSNegative.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSNegative.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSPositive.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSPositive.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSPositive.cutBNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                numPFPsCRUMBS.cutBNB->Fill(numPFPsSliceCRUMBS);
                ratioChosenSummedEnergyCRUMBS.cutBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ERecoSumThetaRecoCRUMBS.cutBNB->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBS.cutBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                slicePurityCRUMBS.cutBNB->Fill(chosenSlicePurityCRUMBS);
                sliceCompletenessCRUMBS.cutBNB->Fill(chosenSliceCompletenessCRUMBS);
                highestPFPCompletenessCRUMBS.cutBNB->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPPurityCRUMBS.cutBNB->Fill(chosenRecoParticleCRUMBS.purity);
            } else if(DLCurrent == 2 && signal == 3){
                numSlices.cutCosmics->Fill(numSlicesInEvent);
                sliceScoreCRUMBS.cutCosmics->Fill(chosenRecoSliceCRUMBS.score);
                numRecoNeutrinos.cutCosmics->Fill(numRecoNeutrinosInEvent);
                primaryVertexXCRUMBS.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBS.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBS.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSNegative.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSNegative.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSPositive.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSPositive.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSPositive.cutCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                numPFPsCRUMBS.cutCosmics->Fill(numPFPsSliceCRUMBS);
                ratioChosenSummedEnergyCRUMBS.cutCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ERecoSumThetaRecoCRUMBS.cutCosmics->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBS.cutCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                slicePurityCRUMBS.cutCosmics->Fill(chosenSlicePurityCRUMBS);
                sliceCompletenessCRUMBS.cutCosmics->Fill(chosenSliceCompletenessCRUMBS);
                highestPFPCompletenessCRUMBS.cutCosmics->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPPurityCRUMBS.cutCosmics->Fill(chosenRecoParticleCRUMBS.purity);
            }

        } else if(cutOrKeep == 1){
            // Passes cuts
            if(DLCurrent == 2 && signal == 1){
                numSlices.keptSignal->Fill(numSlicesInEvent);
                numSlicesCut.Signal->Fill(numSlicesInEvent);
                sliceScoreCRUMBS.keptSignal->Fill(chosenRecoSliceCRUMBS.score);
                sliceScoreCRUMBSCut.Signal->Fill(chosenRecoSliceCRUMBS.score);
                numRecoNeutrinos.keptSignal->Fill(numRecoNeutrinosInEvent);
                numRecoNeutrinosCut.Signal->Fill(numRecoNeutrinosInEvent);
                primaryVertexXCRUMBS.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBS.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBS.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSNegative.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSNegativeCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSNegativeCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSNegative.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSNegativeCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSPositive.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositiveCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSPositive.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSPositiveCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSPositive.keptSignal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSPositiveCut.Signal->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                numPFPsCRUMBS.keptSignal->Fill(numPFPsSliceCRUMBS);
                numPFPsCRUMBSCut.Signal->Fill(numPFPsSliceCRUMBS);
                ratioChosenSummedEnergyCRUMBS.keptSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ratioChosenSummedEnergyCRUMBSCut.Signal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ERecoSumThetaRecoCRUMBS.keptSignal->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoSumThetaRecoCRUMBSCut.Signal->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBS.keptSignal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBSCut.Signal->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                slicePurityCRUMBS.keptSignal->Fill(chosenSlicePurityCRUMBS);
                slicePurityCRUMBSCut.Signal->Fill(chosenSlicePurityCRUMBS);
                sliceCompletenessCRUMBS.keptSignal->Fill(chosenSliceCompletenessCRUMBS);
                sliceCompletenessCRUMBSCut.Signal->Fill(chosenSliceCompletenessCRUMBS);
                highestPFPCompletenessCRUMBS.keptSignal->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPCompletenessCRUMBSCut.Signal->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPPurityCRUMBS.keptSignal->Fill(chosenRecoParticleCRUMBS.purity);
                highestPFPPurityCRUMBSCut.Signal->Fill(chosenRecoParticleCRUMBS.purity);
            } else if(DLCurrent == 2 && signal == 2){
                numSlices.keptBNB->Fill(numSlicesInEvent);
                numSlicesCut.BNB->Fill(numSlicesInEvent);
                sliceScoreCRUMBS.keptBNB->Fill(chosenRecoSliceCRUMBS.score);
                sliceScoreCRUMBSCut.BNB->Fill(chosenRecoSliceCRUMBS.score);
                numRecoNeutrinos.keptBNB->Fill(numRecoNeutrinosInEvent);
                numRecoNeutrinosCut.BNB->Fill(numRecoNeutrinosInEvent);
                primaryVertexXCRUMBS.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBS.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBS.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSNegative.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSNegativeCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSNegativeCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSNegative.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSNegativeCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSPositive.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositiveCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSPositive.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSPositiveCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSPositive.keptBNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSPositiveCut.BNB->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                numPFPsCRUMBS.keptBNB->Fill(numPFPsSliceCRUMBS);
                numPFPsCRUMBSCut.BNB->Fill(numPFPsSliceCRUMBS);
                ratioChosenSummedEnergyCRUMBS.keptBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ratioChosenSummedEnergyCRUMBSCut.BNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ERecoSumThetaRecoCRUMBS.keptBNB->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoSumThetaRecoCRUMBSCut.BNB->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBS.keptBNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBSCut.BNB->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                slicePurityCRUMBS.keptBNB->Fill(chosenSlicePurityCRUMBS);
                slicePurityCRUMBSCut.BNB->Fill(chosenSlicePurityCRUMBS);
                sliceCompletenessCRUMBS.keptBNB->Fill(chosenSliceCompletenessCRUMBS);
                sliceCompletenessCRUMBSCut.BNB->Fill(chosenSliceCompletenessCRUMBS);
                highestPFPCompletenessCRUMBS.keptBNB->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPCompletenessCRUMBSCut.BNB->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPPurityCRUMBS.keptBNB->Fill(chosenRecoParticleCRUMBS.purity);
                highestPFPPurityCRUMBSCut.BNB->Fill(chosenRecoParticleCRUMBS.purity);
            } else if(DLCurrent == 2 && signal == 3){
                numSlices.keptCosmics->Fill(numSlicesInEvent);
                numSlicesCut.Cosmics->Fill(numSlicesInEvent);
                sliceScoreCRUMBS.keptCosmics->Fill(chosenRecoSliceCRUMBS.score);
                sliceScoreCRUMBSCut.Cosmics->Fill(chosenRecoSliceCRUMBS.score);
                numRecoNeutrinos.keptCosmics->Fill(numRecoNeutrinosInEvent);
                numRecoNeutrinosCut.Cosmics->Fill(numRecoNeutrinosInEvent);
                primaryVertexXCRUMBS.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBS.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBS.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSNegative.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSNegativeCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSNegative.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSNegativeCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSNegative.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSNegativeCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexXCRUMBSPositive.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexXCRUMBSPositiveCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vx);
                primaryVertexYCRUMBSPositive.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexYCRUMBSPositiveCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vy); 
                primaryVertexZCRUMBSPositive.keptCosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                primaryVertexZCRUMBSPositiveCut.Cosmics->Fill(chosenRecoNeutrinoCRUMBS.vz); 
                numPFPsCRUMBS.keptCosmics->Fill(numPFPsSliceCRUMBS);
                numPFPsCRUMBSCut.Cosmics->Fill(numPFPsSliceCRUMBS);
                ratioChosenSummedEnergyCRUMBS.keptCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ratioChosenSummedEnergyCRUMBSCut.Cosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
                ERecoSumThetaRecoCRUMBS.keptCosmics->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoSumThetaRecoCRUMBSCut.Cosmics->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBS.keptCosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                ERecoHighestThetaRecoCRUMBSCut.Cosmics->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
                slicePurityCRUMBS.keptCosmics->Fill(chosenSlicePurityCRUMBS);
                slicePurityCRUMBSCut.Cosmics->Fill(chosenSlicePurityCRUMBS);
                sliceCompletenessCRUMBS.keptCosmics->Fill(chosenSliceCompletenessCRUMBS);
                sliceCompletenessCRUMBSCut.Cosmics->Fill(chosenSliceCompletenessCRUMBS);
                highestPFPCompletenessCRUMBS.keptCosmics->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPCompletenessCRUMBSCut.Cosmics->Fill(chosenRecoParticleCRUMBS.completeness);
                highestPFPPurityCRUMBS.keptCosmics->Fill(chosenRecoParticleCRUMBS.purity);
                highestPFPPurityCRUMBSCut.Cosmics->Fill(chosenRecoParticleCRUMBS.purity);
            }
        }
    }
   
    int drawLine = 1;
    int left = 0;
    int right = 1;

    styleDraw(numSlices.canvas, numSlices.keptSignal, numSlices.cutSignal, numSlices.keptBNB, numSlices.cutBNB, numSlices.keptCosmics, numSlices.cutCosmics, 999, 999, 999, 999, (base_path + "numSlices_dist.pdf").c_str(), "topRight");
    styleDrawCut(numSlicesCut.canvas, numSlicesCut.Signal, numSlicesCut.BNB, numSlicesCut.Cosmics, 999, 999, 999, 999, (base_path + "numSlicesCut_dist.pdf").c_str(), "topRight");
    weighted(numSlices, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "numSlices_weighted.pdf").c_str(), "topRight");
    weightedCut(numSlicesCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "numSlicesCut_weighted.pdf").c_str(), "topRight");
    efficiencyCut(numSlicesCut, 999, 999, 999, 999, (base_path + "numSlicesCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Number of Slices", 1);
   
    styleDraw(sliceScoreCRUMBS.canvas, sliceScoreCRUMBS.keptSignal, sliceScoreCRUMBS.cutSignal, sliceScoreCRUMBS.keptBNB, sliceScoreCRUMBS.cutBNB, sliceScoreCRUMBS.keptCosmics, sliceScoreCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_dist.pdf").c_str(), "topLeft");
    styleDrawCut(sliceScoreCRUMBSCut.canvas, sliceScoreCRUMBSCut.Signal, sliceScoreCRUMBSCut.BNB, sliceScoreCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBSCut_dist.pdf").c_str(), "topLeft");
    weighted(sliceScoreCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(sliceScoreCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBSCut_weighted.pdf").c_str(), "topLeft");
    efficiencyCut(sliceScoreCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "CRUMBS Score", -1);

    styleDraw(numRecoNeutrinos.canvas, numRecoNeutrinos.keptSignal, numRecoNeutrinos.cutSignal, numRecoNeutrinos.keptBNB, numRecoNeutrinos.cutBNB, numRecoNeutrinos.keptCosmics, numRecoNeutrinos.cutCosmics, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_dist.pdf").c_str(), "topRight");
    styleDrawCut(numRecoNeutrinosCut.canvas, numRecoNeutrinosCut.Signal, numRecoNeutrinosCut.BNB, numRecoNeutrinosCut.Cosmics, 999, 999, 999, 999, (base_path + "numRecoNeutrinosCut_dist.pdf").c_str(), "topRight");
    weighted(numRecoNeutrinos, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_weighted.pdf").c_str(), "topRight");
    weightedCut(numRecoNeutrinosCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "numRecoNeutrinosCut_weighted.pdf").c_str(), "topRight");
    efficiencyCut(numRecoNeutrinosCut, 999, 999, 999, 999, (base_path + "numRecoNeutrinosCut").c_str(), "bottomRight", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Number of Reco Neutrinos", 1);
    
    styleDraw(primaryVertexXCRUMBS.canvas, primaryVertexXCRUMBS.keptSignal, primaryVertexXCRUMBS.cutSignal, primaryVertexXCRUMBS.keptBNB, primaryVertexXCRUMBS.cutBNB, primaryVertexXCRUMBS.keptCosmics, primaryVertexXCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBS_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexXCRUMBSCut.canvas, primaryVertexXCRUMBSCut.Signal, primaryVertexXCRUMBSCut.BNB, primaryVertexXCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexXCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(primaryVertexXCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexXCRUMBSCut, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex X Coord", -1);

    styleDraw(primaryVertexYCRUMBS.canvas, primaryVertexYCRUMBS.keptSignal, primaryVertexYCRUMBS.cutSignal, primaryVertexYCRUMBS.keptBNB, primaryVertexYCRUMBS.cutBNB, primaryVertexYCRUMBS.keptCosmics, primaryVertexYCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBS_dist.pdf").c_str(), "topLeft");
    styleDrawCut(primaryVertexYCRUMBSCut.canvas, primaryVertexYCRUMBSCut.Signal, primaryVertexYCRUMBSCut.BNB, primaryVertexYCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSCut_dist.pdf").c_str(), "topLeft");
    weighted(primaryVertexYCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(primaryVertexYCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSCut_weighted.pdf").c_str(), "topLeft");
    efficiencyCut(primaryVertexYCRUMBSCut, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex Y Coord", -1);
    
    styleDraw(primaryVertexZCRUMBS.canvas, primaryVertexZCRUMBS.keptSignal, primaryVertexZCRUMBS.cutSignal, primaryVertexZCRUMBS.keptBNB, primaryVertexZCRUMBS.cutBNB, primaryVertexZCRUMBS.keptCosmics, primaryVertexZCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBS_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexZCRUMBSCut.canvas, primaryVertexZCRUMBSCut.Signal, primaryVertexZCRUMBSCut.BNB, primaryVertexZCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexZCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(primaryVertexZCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexZCRUMBSCut, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex Z Coord", -1);
   
    styleDraw(primaryVertexXCRUMBSNegative.canvas, primaryVertexXCRUMBSNegative.keptSignal, primaryVertexXCRUMBSNegative.cutSignal, primaryVertexXCRUMBSNegative.keptBNB, primaryVertexXCRUMBSNegative.cutBNB, primaryVertexXCRUMBSNegative.keptCosmics, primaryVertexXCRUMBSNegative.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSNegative_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexXCRUMBSNegativeCut.canvas, primaryVertexXCRUMBSNegativeCut.Signal, primaryVertexXCRUMBSNegativeCut.BNB, primaryVertexXCRUMBSNegativeCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSNegativeCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexXCRUMBSNegative, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSNegative_weighted.pdf").c_str(), "bottomLeft");
    weightedCut(primaryVertexXCRUMBSNegativeCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSNegativeCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexXCRUMBSNegativeCut, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSNegativeCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex X Coord", -1);
    
    styleDraw(primaryVertexYCRUMBSNegative.canvas, primaryVertexYCRUMBSNegative.keptSignal, primaryVertexYCRUMBSNegative.cutSignal, primaryVertexYCRUMBSNegative.keptBNB, primaryVertexYCRUMBSNegative.cutBNB, primaryVertexYCRUMBSNegative.keptCosmics, primaryVertexYCRUMBSNegative.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSNegative_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexYCRUMBSNegativeCut.canvas, primaryVertexYCRUMBSNegativeCut.Signal, primaryVertexYCRUMBSNegativeCut.BNB, primaryVertexYCRUMBSNegativeCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSNegativeCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexYCRUMBSNegative, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSNegative_weighted.pdf").c_str(), "bottomLeft");
    weightedCut(primaryVertexYCRUMBSNegativeCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSNegativeCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexYCRUMBSNegativeCut, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSNegativeCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex Y Coord", -1);

    styleDraw(primaryVertexZCRUMBSNegative.canvas, primaryVertexZCRUMBSNegative.keptSignal, primaryVertexZCRUMBSNegative.cutSignal, primaryVertexZCRUMBSNegative.keptBNB, primaryVertexZCRUMBSNegative.cutBNB, primaryVertexZCRUMBSNegative.keptCosmics, primaryVertexZCRUMBSNegative.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSNegative_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexZCRUMBSNegativeCut.canvas, primaryVertexZCRUMBSNegativeCut.Signal, primaryVertexZCRUMBSNegativeCut.BNB, primaryVertexZCRUMBSNegativeCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSNegativeCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexZCRUMBSNegative, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSNegative_weighted.pdf").c_str(), "bottomLeft");
    weightedCut(primaryVertexZCRUMBSNegativeCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSNegativeCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexZCRUMBSNegativeCut, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSNegativeCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex Z Coord", -1);
    
    styleDraw(primaryVertexXCRUMBSPositive.canvas, primaryVertexXCRUMBSPositive.keptSignal, primaryVertexXCRUMBSPositive.cutSignal, primaryVertexXCRUMBSPositive.keptBNB, primaryVertexXCRUMBSPositive.cutBNB, primaryVertexXCRUMBSPositive.keptCosmics, primaryVertexXCRUMBSPositive.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSPositive_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexXCRUMBSPositiveCut.canvas, primaryVertexXCRUMBSPositiveCut.Signal, primaryVertexXCRUMBSPositiveCut.BNB, primaryVertexXCRUMBSPositiveCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSPositiveCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexXCRUMBSPositive, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSPositive_weighted.pdf").c_str(), "bottomLeft");
    weightedCut(primaryVertexXCRUMBSPositiveCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSPositiveCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexXCRUMBSPositiveCut, 999, 999, 999, 999, (base_path + "primaryVertexXCRUMBSPositiveCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex X Coord", 1);
    
    styleDraw(primaryVertexYCRUMBSPositive.canvas, primaryVertexYCRUMBSPositive.keptSignal, primaryVertexYCRUMBSPositive.cutSignal, primaryVertexYCRUMBSPositive.keptBNB, primaryVertexYCRUMBSPositive.cutBNB, primaryVertexYCRUMBSPositive.keptCosmics, primaryVertexYCRUMBSPositive.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSPositive_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexYCRUMBSPositiveCut.canvas, primaryVertexYCRUMBSPositiveCut.Signal, primaryVertexYCRUMBSPositiveCut.BNB, primaryVertexYCRUMBSPositiveCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSPositiveCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexYCRUMBSPositive, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSPositive_weighted.pdf").c_str(), "bottomLeft");
    weightedCut(primaryVertexYCRUMBSPositiveCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSPositiveCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexYCRUMBSPositiveCut, 999, 999, 999, 999, (base_path + "primaryVertexYCRUMBSPositiveCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex Y Coord", 1);
    
    styleDraw(primaryVertexZCRUMBSPositive.canvas, primaryVertexZCRUMBSPositive.keptSignal, primaryVertexZCRUMBSPositive.cutSignal, primaryVertexZCRUMBSPositive.keptBNB, primaryVertexZCRUMBSPositive.cutBNB, primaryVertexZCRUMBSPositive.keptCosmics, primaryVertexZCRUMBSPositive.cutCosmics, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSPositive_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(primaryVertexZCRUMBSPositiveCut.canvas, primaryVertexZCRUMBSPositiveCut.Signal, primaryVertexZCRUMBSPositiveCut.BNB, primaryVertexZCRUMBSPositiveCut.Cosmics, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSPositiveCut_dist.pdf").c_str(), "bottomLeft");
    weighted(primaryVertexZCRUMBSPositive, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSPositive_weighted.pdf").c_str(), "bottomLeft");
    weightedCut(primaryVertexZCRUMBSPositiveCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSPositiveCut_weighted.pdf").c_str(), "bottomLeft");
    efficiencyCut(primaryVertexZCRUMBSPositiveCut, 999, 999, 999, 999, (base_path + "primaryVertexZCRUMBSPositiveCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Primary Vertex Z Coord", 1);
    
    styleDraw(numPFPsCRUMBS.canvas, numPFPsCRUMBS.keptSignal, numPFPsCRUMBS.cutSignal, numPFPsCRUMBS.keptBNB, numPFPsCRUMBS.cutBNB, numPFPsCRUMBS.keptCosmics, numPFPsCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_dist.pdf").c_str(), "topRight");
    styleDrawCut(numPFPsCRUMBSCut.canvas, numPFPsCRUMBSCut.Signal, numPFPsCRUMBSCut.BNB, numPFPsCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "numPFPsCRUMBSCut_dist.pdf").c_str(), "topRight");
    weighted(numPFPsCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_weighted.pdf").c_str(), "topRight");
    weightedCut(numPFPsCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "numPFPsCRUMBSCut_weighted.pdf").c_str(), "topRight");
    efficiencyCut(numPFPsCRUMBSCut, 999, 999, 999, 999, (base_path + "numPFPsCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Number of PFPs", 1);
    
    styleDraw(ratioChosenSummedEnergyCRUMBS.canvas, ratioChosenSummedEnergyCRUMBS.keptSignal, ratioChosenSummedEnergyCRUMBS.cutSignal, ratioChosenSummedEnergyCRUMBS.keptBNB, ratioChosenSummedEnergyCRUMBS.cutBNB, ratioChosenSummedEnergyCRUMBS.keptCosmics, ratioChosenSummedEnergyCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_dist.pdf").c_str(), "topLeft");
    styleDrawCut(ratioChosenSummedEnergyCRUMBSCut.canvas, ratioChosenSummedEnergyCRUMBSCut.Signal, ratioChosenSummedEnergyCRUMBSCut.BNB, ratioChosenSummedEnergyCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBSCut_dist.pdf").c_str(), "topLeft");
    weighted(ratioChosenSummedEnergyCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(ratioChosenSummedEnergyCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBSCut_weighted.pdf").c_str(), "topLeft");
    efficiencyCut(ratioChosenSummedEnergyCRUMBSCut, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", -1);
    
    styleDraw(ERecoSumThetaRecoCRUMBS.canvas, ERecoSumThetaRecoCRUMBS.keptSignal, ERecoSumThetaRecoCRUMBS.cutSignal, ERecoSumThetaRecoCRUMBS.keptBNB, ERecoSumThetaRecoCRUMBS.cutBNB, ERecoSumThetaRecoCRUMBS.keptCosmics, ERecoSumThetaRecoCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_dist.pdf").c_str(), "topRight");
    styleDrawCut(ERecoSumThetaRecoCRUMBSCut.canvas, ERecoSumThetaRecoCRUMBSCut.Signal, ERecoSumThetaRecoCRUMBSCut.BNB, ERecoSumThetaRecoCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBSCut_dist.pdf").c_str(), "topRight");
    weighted(ERecoSumThetaRecoCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_weighted.pdf").c_str(), "topRight");
    weightedCut(ERecoSumThetaRecoCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBSCut_weighted.pdf").c_str(), "topRight");
    efficiencyCut(ERecoSumThetaRecoCRUMBSCut, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1);
    
    styleDraw(ERecoHighestThetaRecoCRUMBS.canvas, ERecoHighestThetaRecoCRUMBS.keptSignal, ERecoHighestThetaRecoCRUMBS.cutSignal, ERecoHighestThetaRecoCRUMBS.keptBNB, ERecoHighestThetaRecoCRUMBS.cutBNB, ERecoHighestThetaRecoCRUMBS.keptCosmics, ERecoHighestThetaRecoCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_dist.pdf").c_str(), "topRight");
    styleDrawCut(ERecoHighestThetaRecoCRUMBSCut.canvas, ERecoHighestThetaRecoCRUMBSCut.Signal, ERecoHighestThetaRecoCRUMBSCut.BNB, ERecoHighestThetaRecoCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBSCut_dist.pdf").c_str(), "topRight");
    weighted(ERecoHighestThetaRecoCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_weighted.pdf").c_str(), "topRight");
    weightedCut(ERecoHighestThetaRecoCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBSCut_weighted.pdf").c_str(), "topRight");
    efficiencyCut(ERecoHighestThetaRecoCRUMBSCut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 1);
    
    styleDraw(slicePurityCRUMBS.canvas, slicePurityCRUMBS.keptSignal, slicePurityCRUMBS.cutSignal, slicePurityCRUMBS.keptBNB, slicePurityCRUMBS.cutBNB, slicePurityCRUMBS.keptCosmics, slicePurityCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "slicePurityCRUMBS_dist.pdf").c_str(), "bottomLeft");
    styleDrawCut(slicePurityCRUMBSCut.canvas, slicePurityCRUMBSCut.Signal, slicePurityCRUMBSCut.BNB, slicePurityCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "slicePurityCRUMBSCut_dist.pdf").c_str(), "bottomLeft");
    weighted(slicePurityCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "slicePurityCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(slicePurityCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "slicePurityCRUMBSCut_weighted.pdf").c_str(), "topLeft");
    efficiencyCut(slicePurityCRUMBSCut, 999, 999, 999, 999, (base_path + "slicePurityCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Purity", -1);
    
    styleDraw(sliceCompletenessCRUMBS.canvas, sliceCompletenessCRUMBS.keptSignal, sliceCompletenessCRUMBS.cutSignal, sliceCompletenessCRUMBS.keptBNB, sliceCompletenessCRUMBS.cutBNB, sliceCompletenessCRUMBS.keptCosmics, sliceCompletenessCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBS_dist.pdf").c_str(), "topRight");
    styleDrawCut(sliceCompletenessCRUMBSCut.canvas, sliceCompletenessCRUMBSCut.Signal, sliceCompletenessCRUMBSCut.BNB, sliceCompletenessCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBSCut_dist.pdf").c_str(), "topRight");
    weighted(sliceCompletenessCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBS_weighted.pdf").c_str(), "topRight");
    weightedCut(sliceCompletenessCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBSCut_weighted.pdf").c_str(), "topRight");
    efficiencyCut(sliceCompletenessCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Completeness", -1);
    
    styleDraw(highestPFPCompletenessCRUMBS.canvas, highestPFPCompletenessCRUMBS.keptSignal, highestPFPCompletenessCRUMBS.cutSignal, highestPFPCompletenessCRUMBS.keptBNB, highestPFPCompletenessCRUMBS.cutBNB, highestPFPCompletenessCRUMBS.keptCosmics, highestPFPCompletenessCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBS_dist.pdf").c_str(), "topLeft");
    styleDrawCut(highestPFPCompletenessCRUMBSCut.canvas, highestPFPCompletenessCRUMBSCut.Signal, highestPFPCompletenessCRUMBSCut.BNB, highestPFPCompletenessCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBSCut_dist.pdf").c_str(), "topLeft");
    weighted(highestPFPCompletenessCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(highestPFPCompletenessCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBSCut_weighted.pdf").c_str(), "topLeft");
    efficiencyCut(highestPFPCompletenessCRUMBSCut, 999, 999, 999, 999, (base_path + "highestPFPCompletenessCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Completeness", -1);
    
    styleDraw(highestPFPPurityCRUMBS.canvas, highestPFPPurityCRUMBS.keptSignal, highestPFPPurityCRUMBS.cutSignal, highestPFPPurityCRUMBS.keptBNB, highestPFPPurityCRUMBS.cutBNB, highestPFPPurityCRUMBS.keptCosmics, highestPFPPurityCRUMBS.cutCosmics, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBS_dist.pdf").c_str(), "topLeft");
    styleDrawCut(highestPFPPurityCRUMBSCut.canvas, highestPFPPurityCRUMBSCut.Signal, highestPFPPurityCRUMBSCut.BNB, highestPFPPurityCRUMBSCut.Cosmics, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBSCut_dist.pdf").c_str(), "topLeft");
    weighted(highestPFPPurityCRUMBS, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBS_weighted.pdf").c_str(), "topLeft");
    weightedCut(highestPFPPurityCRUMBSCut, signalWeight, BNBWeight, cosmicsWeight, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBSCut_weighted.pdf").c_str(), "topLeft");
    efficiencyCut(highestPFPPurityCRUMBSCut, 999, 999, 999, 999, (base_path + "highestPFPPurityCRUMBSCut").c_str(), "bottomLeft", signalWeight, BNBWeight, cosmicsWeight, nullptr, &left, "Purity", -1);
 
}
