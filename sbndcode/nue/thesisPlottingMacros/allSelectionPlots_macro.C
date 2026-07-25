#include <vector>
#include <map>
#include <tuple>
#include <utility>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
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
#include <fstream>
#include <TLegend.h>
#include <THStack.h>
#include <set>
#include <utility>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <iomanip>
#include <TH2D.h>
#include <TProfile.h>
#include <TGaxis.h>
#include "TRandom3.h"
#include <algorithm>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

struct weights_struct{
    double signalNuE = 0;
    double BNBNuE = 0;
    double cosmicsNuE = 0;
    double NuENuE = 0;
};

struct recoilElectron_struct{
    double energy;
    double angle;
    double dx;
    double dy;
    double dz;
};

struct highestEnergyPFP_struct{
    double PFPID = -999999;
    double energy = -999999;
    double theta = -999999;
    double dx = -999999;
    double dy = -999999;
    double dz = -999999;
    double vx = -999999;
    double vy = -999999;
    double vz = -999999;
    double completeness = -999999;
    double purity = -999999;
    double trackscore = -999999;
    double primary = -999999;
    double truePDG = -999999;
    double trueOrigin = -999999;
    double trueInt = -999999;
    double bestPlanedEdx = -999999;
    double razzledPDG11 = -999999;
    double razzledPDG13 = -999999;
    double razzledPDG22 = -999999;
    double razzledPDG211 = -999999;
    double razzledPDG2212 = -999999;
    double razzledBestPDG = -999999;
    double showerLength = -999999;
    double showerOpenAngle = -999999;
    double showerBestPlaneEnergy = -999999;
    double trueVX = -999999;
    double trueVY = -999999;
    double trueVZ = -999999;
    double trueEndX = -999999;
    double trueEndY = -999999;
    double trueEndZ = -999999;
    double trueLength = -999999;
    double numHits = -999999;
    double clearCosmic = -999999;
};

typedef struct{
    TCanvas* canvas;
    TH1D* baseHist;
    TH1D* nuECosmic;
    TH1D* nuESignal;
    TH1D* nuESignalFuzzy;
    TH1D* nuEBNB;
    TH1D* nuEBNBFuzzy;
} histGroup_struct;

typedef struct{
    TCanvas* canvas;
    TH1D* baseHist;
    TH1D* nu_e;
    TH1D* NCNpi0;
    TH1D* otherNC;
    TH1D* CCnumu;
    TH1D* CCnue;
    TH1D* dirt;
    TH1D* nu_eDirt;
    TH1D* cosmic;
    TH1D* other;
    TH1D* nu_eFuzzy;
} splitHistGroup_struct;

struct eventCounter_struct{
    double nuE = 0;
    double NCNPi0 = 0;
    double otherNC = 0;
    double CCnumu = 0;
    double CCnue = 0;
    double dirt = 0;
    double nuEDirt = 0;
    double nuEFuzzy = 0;
    double cosmic = 0;
    double other = 0;
};

struct beforeEventCount_struct{
    double signal = 0;
    double background = 0;
    eventCounter_struct splitInt;
};

struct eventCounting_struct{
    double clearCosmicsSig = 0;
    double clearCosmicsBack = 0;
    eventCounter_struct clearCosmicsIntSplit;
    double numPFPs0Sig = 0;
    double numPFPs0Back = 0;
    eventCounter_struct numPFPs0IntSplit;
    double numRecoNeut0Sig = 0;
    double numRecoNeut0Back = 0;
    eventCounter_struct numRecoNeut0IntSplit;
    double FVSig = 0;
    double FVBack = 0;
    eventCounter_struct FVIntSplit;
    double crumbsSig = 0;
    double crumbsBack = 0;
    eventCounter_struct crumbsIntSplit;
    double primaryPFPSig = 0;
    double primaryPFPBack = 0;
    eventCounter_struct primaryPFPIntSplit;
    double razzled211Sig = 0;
    double razzled211Back = 0;
    eventCounter_struct razzled211IntSplit;
    double razzled11Sig = 0;
    double razzled11Back = 0;
    eventCounter_struct razzled11IntSplit;
    double dEdxSig = 0;
    double dEdxBack = 0;
    eventCounter_struct dEdxIntSplit;
    double ETheta2Sig = 0;    
    double ETheta2Back = 0;    
    eventCounter_struct ETheta2IntSplit;    
};


histGroup_struct createHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);
    
    TH1D* base = new TH1D(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1D*) base->Clone((baseName + "_nuECosmic").c_str()),
        (TH1D*) base->Clone((baseName + "_nuESignal").c_str()),
        (TH1D*) base->Clone((baseName + "_nuESignalFuzzy").c_str()),
        (TH1D*) base->Clone((baseName + "_nuEBNB").c_str()),
        (TH1D*) base->Clone((baseName + "_nuEBNBFuzzy").c_str()),
    };    
}

splitHistGroup_struct createSplitHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);

    TH1D* base = new TH1D(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1D*) base->Clone((baseName + "_nu_e").c_str()),
        (TH1D*) base->Clone((baseName + "_NCNpi0").c_str()),
        (TH1D*) base->Clone((baseName + "_otherNC").c_str()),
        (TH1D*) base->Clone((baseName + "_CCnumu").c_str()),
        (TH1D*) base->Clone((baseName + "_CCnue").c_str()),
        (TH1D*) base->Clone((baseName + "_dirt").c_str()),
        (TH1D*) base->Clone((baseName + "_nu_eDirt").c_str()),
        (TH1D*) base->Clone((baseName + "_cosmic").c_str()),
        (TH1D*) base->Clone((baseName + "_other").c_str()),
        (TH1D*) base->Clone((baseName + "_nu_eFuzzy").c_str())
    };
}

void fillHistogram(histGroup_struct* hist, int DLCurrent, int signal, int type, double value, weights_struct* weight){
    if(!hist || DLCurrent != 5) return;
    
    TH1* target = nullptr;

    switch(type){
        case 0: if(signal != 4) target = hist->nuECosmic; break;
        case 1: if(signal == 1) target = hist->nuESignal; break;
        case 2: if(signal == 1) target = hist->nuESignalFuzzy; break;
        case 3: target = hist->nuEBNB; break;
        case 4: target = hist->nuEBNBFuzzy; break;
        case 5: target = hist->nuEBNB; break;
        case 6: target = hist->nuEBNBFuzzy; break;
    }

    if(!target) return;

    double w = 0.0;

    if(signal == 1) w = weight->signalNuE;
    else if(signal == 2 && (type != 5 && type != 6)) w = weight->BNBNuE;
    else if(signal == 3) w = weight->cosmicsNuE;
    else if((type == 5 || type == 6) && (signal == 2 || signal == 4)) w = weight->NuENuE;

    target->Fill(value, w);
}

void fillSplitIntHistogram(splitHistGroup_struct* hist, int DLCurrent, int signal, int type, double value, weights_struct* weight){
    if(!hist || DLCurrent != 5) return;

    TH1* target = nullptr;

    switch(type){
        case 0: if(signal != 4) target = hist->cosmic; break;
        case 1: if(signal == 1) target = hist->nu_e; break;
        case 2: target = hist->NCNpi0; break;
        case 3: target = hist->otherNC; break;
        case 4: target = hist->CCnumu; break;
        case 5: target = hist->CCnue; break;
        case 6: target = hist->dirt; break;
        case 7: if(signal == 1) target = hist->nu_eDirt; break;
        case 8: target = hist->other; break;
        case 9: if(signal == 1) target = hist->nu_eFuzzy; break;
        case 15: break;
    }

    if(!target) return;

    double w = 0.0;

    if(signal == 1) w = weight->signalNuE;
    else if(signal == 2 && (type != 5)) w = weight->BNBNuE;
    else if(signal == 3) w = weight->cosmicsNuE;
    else if((signal == 2 || signal == 4) && type == 5) w = weight->NuENuE;

    target->Fill(value, w);
}


void styleDrawAll(histGroup_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr,
                  bool includeSignal = true, bool includeSignalFuzzy = true,
                  bool includeBNB = true, bool includeBNBFuzzy = true,
                  bool includeCosmic = true,
                  bool includeDLUboone = true, bool includeDLNuE = true,
                  bool includeBDT = true,
                  bool useLogScale = false, bool bestPDGPlot = false)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    if (useLogScale)
        hists.canvas->SetLogy(1);
    else
        hists.canvas->SetLogy(0);

    std::vector<TH1D*> allHists = {
        hists.currentSignal, hists.ubooneSignal, hists.nuESignal,
        hists.currentSignalFuzzy, hists.ubooneSignalFuzzy, hists.nuESignalFuzzy,
        hists.currentBNB, hists.ubooneBNB, hists.nuEBNB,
        hists.currentBNBFuzzy, hists.ubooneBNBFuzzy, hists.nuEBNBFuzzy,
        hists.currentCosmic, hists.ubooneCosmic, hists.nuECosmic
    };

    if (useLogScale) {
        for (auto* hist : allHists) {
            if (!hist) continue;
            int nBins = hist->GetNbinsX();
            for (int i = 1; i <= nBins; ++i) {
                if (hist->GetBinContent(i) <= 0)
                    hist->SetBinContent(i, 1e-6);
            }
        }
    }

    if(bestPDGPlot){
        for(auto* hist : allHists){
            if(!hist) continue;
            hist->GetXaxis()->SetBinLabel(1, "e^{-}");
            hist->GetXaxis()->SetBinLabel(2, "#mu^{-}");
            hist->GetXaxis()->SetBinLabel(3, "#gamma");
            hist->GetXaxis()->SetBinLabel(4, "#pi^{#pm}");
            hist->GetXaxis()->SetBinLabel(5, "p");
            hist->GetXaxis()->SetBinLabel(6, "Other");
            hist->LabelsOption("h");
        }
    }

    for (auto* hist : allHists)
        if (hist) hist->SetStats(0);

    hists.currentCosmic->SetLineWidth(2);  hists.currentCosmic->SetLineColor(kPink+9);
    hists.ubooneCosmic->SetLineWidth(2);   hists.ubooneCosmic->SetLineColor(kPink+1);
    hists.nuECosmic->SetLineWidth(2);      hists.nuECosmic->SetLineColor(kPink-2);

    hists.currentSignal->SetLineWidth(2);  hists.currentSignal->SetLineColor(kBlue+1);
    hists.ubooneSignal->SetLineWidth(2);   hists.ubooneSignal->SetLineColor(kBlue-7);
    hists.nuESignal->SetLineWidth(2);      hists.nuESignal->SetLineColor(kAzure+5);

    hists.currentSignalFuzzy->SetLineWidth(2);  hists.currentSignalFuzzy->SetLineColor(kGreen+3);
    hists.ubooneSignalFuzzy->SetLineWidth(2);   hists.ubooneSignalFuzzy->SetLineColor(kGreen+1);
    hists.nuESignalFuzzy->SetLineWidth(2);      hists.nuESignalFuzzy->SetLineColor(kGreen-7);

    hists.currentBNB->SetLineWidth(2);  hists.currentBNB->SetLineColor(kOrange+7);
    hists.ubooneBNB->SetLineWidth(2);   hists.ubooneBNB->SetLineColor(kOrange+6);
    hists.nuEBNB->SetLineWidth(2);      hists.nuEBNB->SetLineColor(kOrange-5);

    hists.currentBNBFuzzy->SetLineWidth(2);  hists.currentBNBFuzzy->SetLineColor(kViolet+1);
    hists.ubooneBNBFuzzy->SetLineWidth(2);   hists.ubooneBNBFuzzy->SetLineColor(kViolet-7);
    hists.nuEBNBFuzzy->SetLineWidth(2);      hists.nuEBNBFuzzy->SetLineColor(kViolet+4);

    if((ymin != 999) && (ymax != 999)){
        hists.currentSignal->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.ubooneSignal->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuESignal->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.currentSignalFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.ubooneSignalFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuESignalFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.currentBNB->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.ubooneBNB->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEBNB->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.currentBNBFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.ubooneBNBFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEBNBFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.currentCosmic->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.ubooneCosmic->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuECosmic->GetYaxis()->SetRangeUser(ymin, ymax);
    }
    
    if((xmin != 999) && (xmax != 999)){
        hists.currentSignal->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.ubooneSignal->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuESignal->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.currentSignalFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.ubooneSignalFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuESignalFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.currentBNB->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.ubooneBNB->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEBNB->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.currentBNBFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.ubooneBNBFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEBNBFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.currentCosmic->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.ubooneCosmic->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuECosmic->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    double maxYValue = 0.0;
    for (auto* hist : allHists)
        if (hist && hist->GetMaximum() > maxYValue)
            maxYValue = hist->GetMaximum();

    //std::cout << "maxYValue = " << maxYValue << std::endl;
    double yminVal = useLogScale ? 1e-1 : 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = useLogScale ? (maxYValue * 100.0) : (maxYValue * 1.1);
        //std::cout << "setting yaxis to " << yminVal << ", " << ymaxVal << std::endl;
        for (auto* hist : allHists)
            if (hist) hist->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    bool first = true;
    auto draw = [&](TH1* hist){ if (hist) { hist->Draw(first ? "hist" : "histsame"); first = false; } };

    auto variantAllowed = [&](const std::string& name) {
        bool isBDT = name.find("current") != std::string::npos;
        bool isDLUboone = name.find("uboone") != std::string::npos;
        bool isDLNuE = name.find("nuE") != std::string::npos;

        if (!includeBDT && isBDT) return false;
        if (!includeDLUboone && isDLUboone) return false;
        if (!includeDLNuE && isDLNuE) return false;
        return true;
    };

    if (includeSignal) {
        if (variantAllowed("currentSignal")) draw(hists.currentSignal);
        if (variantAllowed("ubooneSignal")) draw(hists.ubooneSignal);
        if (variantAllowed("nuESignal")) draw(hists.nuESignal);
    }
    if (includeSignalFuzzy) {
        if (variantAllowed("currentSignalFuzzy")) draw(hists.currentSignalFuzzy);
        if (variantAllowed("ubooneSignalFuzzy")) draw(hists.ubooneSignalFuzzy);
        if (variantAllowed("nuESignalFuzzy")) draw(hists.nuESignalFuzzy);
    }
    if (includeBNB) {
        if (variantAllowed("currentBNB")) draw(hists.currentBNB);
        if (variantAllowed("ubooneBNB")) draw(hists.ubooneBNB);
        if (variantAllowed("nuEBNB")) draw(hists.nuEBNB);
    }
    if (includeBNBFuzzy) {
        if (variantAllowed("currentBNBFuzzy")) draw(hists.currentBNBFuzzy);
        if (variantAllowed("ubooneBNBFuzzy")) draw(hists.ubooneBNBFuzzy);
        if (variantAllowed("nuEBNBFuzzy")) draw(hists.nuEBNBFuzzy);
    }
    if (includeCosmic) {
        if (variantAllowed("currentCosmic")) draw(hists.currentCosmic);
        if (variantAllowed("ubooneCosmic")) draw(hists.ubooneCosmic);
        if (variantAllowed("nuECosmic")) draw(hists.nuECosmic);
    }

    hists.currentSignal->SetStats(0);
    hists.currentSignal->GetXaxis()->SetTickLength(0.04);
    hists.currentSignal->GetYaxis()->SetTickLength(0.03);
    hists.currentSignal->GetXaxis()->SetTickSize(0.02);
    hists.currentSignal->GetYaxis()->SetTickSize(0.02);

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    std::vector<std::pair<TH1*, std::string>> legendEntries;

    auto addLegendIf = [&](TH1* hist, const std::string& label, const std::string& name){
        if (hist && variantAllowed(name)) legendEntries.emplace_back(hist, label);
    };

    if (includeSignal) {
        addLegendIf(hists.currentSignal, "Signal, Pandora BDT SBND (without Refinement)", "currentSignal");
        addLegendIf(hists.ubooneSignal, "Signal, Pandora Deep Learning: #muBooNE/BNB Tune", "ubooneSignal");
        addLegendIf(hists.nuESignal, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "nuESignal");
    }
    if (includeSignalFuzzy) {
        addLegendIf(hists.currentSignalFuzzy, "Signal Fuzzy, Pandora BDT SBND (without Refinement)", "currentSignalFuzzy");
        addLegendIf(hists.ubooneSignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: #muBooNE/BNB Tune", "ubooneSignalFuzzy");
        addLegendIf(hists.nuESignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "nuESignalFuzzy");
    }
    if (includeBNB) {
        addLegendIf(hists.currentBNB, "BNB, Pandora BDT SBND (without Refinement)", "currentBNB");
        addLegendIf(hists.ubooneBNB, "BNB, Pandora Deep Learning: #muBooNE/BNB Tune", "ubooneBNB");
        addLegendIf(hists.nuEBNB, "BNB, Pandora Deep Learning: SBND Nu+E Tune", "nuEBNB");
    }
    if (includeBNBFuzzy) {
        addLegendIf(hists.currentBNBFuzzy, "BNB Fuzzy, Pandora BDT SBND (without Refinement)", "currentBNBFuzzy");
        addLegendIf(hists.ubooneBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: #muBooNE/BNB Tune", "ubooneBNBFuzzy");
        addLegendIf(hists.nuEBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "nuEBNBFuzzy");
    }
    if (includeCosmic) {
        addLegendIf(hists.currentCosmic, "Cosmic, Pandora BDT SBND (without Refinement)", "currentCosmic");
        addLegendIf(hists.ubooneCosmic, "Cosmic, Pandora Deep Learning: #muBooNE/BNB Tune", "ubooneCosmic");
        addLegendIf(hists.nuECosmic, "Cosmic, Pandora Deep Learning: SBND Nu+E Tune", "nuECosmic");
    }

    int nEntries = legendEntries.size();
    double height = std::max(0.025 * nEntries, 0.03);

    if(legendLocation == "topRight"){ Lxmin=0.62; Lymax=0.863; Lxmax=0.87; Lymin=Lymax - height; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.38; Lymin=Lymax - height; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.62; Lymin=0.137; Lxmax=0.87; Lymax=Lymin + height; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymin=0.137; Lxmax=0.38; Lymax=Lymin + height; }

    auto legend = new TLegend(Lxmin, Lymin, Lxmax, Lymax);
    for (auto& [hist, label] : legendEntries)
        legend->AddEntry(hist, label.c_str(), "f");

    legend->SetTextSize(0.015);
    legend->SetMargin(0.12);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, yminVal, 1.022, maxYValue);
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        double ymaxValLine = maxYValue * (useLogScale ? 0.5 : 0.95);
        TLatex* latex = new TLatex((*linePos == 0 ? 1.022 - 0.2 : 1.022 + 0.2), ymaxValLine, "2m_{e}");
        latex->SetTextSize(0.035); 
        latex->SetTextAlign(11);
        latex->Draw("same");
    }

    hists.canvas->SaveAs(filename);
}

void styleDrawBackSig(histGroup_struct hists,
                      double ymin, double ymax, double xmin, double xmax,
                      const char* filename, const std::string& legendLocation,
                      bool includeCurrent = true, bool includeUboone = true, bool includeNuE = true,
                      bool useLogScale = false, bool bestPDGPlot = false)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();
    hists.canvas->SetLogy(useLogScale);

    auto combine = [useLogScale](TH1D* a, TH1D* b, TH1D* c, TH1D* d, const char* name) -> TH1D* {
        TH1D* combo = nullptr;
        if (a) combo = (TH1D*)a->Clone(name);
        else if (b) combo = (TH1D*)b->Clone(name);
        else if (c) combo = (TH1D*)c->Clone(name);
        else if (d) combo = (TH1D*)d->Clone(name);
        else return nullptr;

        combo->Reset();
        if (a) combo->Add(a);
        if (b) combo->Add(b);
        if (c) combo->Add(c);
        if (d) combo->Add(d);

        if (useLogScale) {
            for (int i = 1; i <= combo->GetNbinsX(); ++i)
                if (combo->GetBinContent(i) <= 0)
                    combo->SetBinContent(i, 1e-6);
        }

        return combo;
    };

    TH1D* currentSignalCombo     = combine(hists.currentSignal, nullptr, nullptr, nullptr, "currentSignalCombo");
    TH1D* currentBackgroundCombo = combine(hists.currentBNB, hists.currentBNBFuzzy, hists.currentCosmic, hists.currentSignalFuzzy, "currentBackgroundCombo");

    TH1D* ubooneSignalCombo     = combine(hists.ubooneSignal, nullptr, nullptr, nullptr, "ubooneSignalCombo");
    TH1D* ubooneBackgroundCombo = combine(hists.ubooneBNB, hists.ubooneBNBFuzzy, hists.ubooneCosmic, hists.ubooneSignalFuzzy, "ubooneBackgroundCombo");

    TH1D* nuESignalCombo     = combine(hists.nuESignal, nullptr, nullptr, nullptr, "nuESignalCombo");
    TH1D* nuEBackgroundCombo = combine(hists.nuEBNB, hists.nuEBNBFuzzy, hists.nuECosmic, hists.nuESignalFuzzy, "nuEBackgroundCombo");

    std::vector<TH1D*> allHists = {
        currentSignalCombo, currentBackgroundCombo,
        ubooneSignalCombo, ubooneBackgroundCombo,
        nuESignalCombo, nuEBackgroundCombo
    };
    
    if(bestPDGPlot){
        for(auto* hist : allHists){
            if(!hist) continue;
            hist->GetXaxis()->SetBinLabel(1, "e^{-}");
            hist->GetXaxis()->SetBinLabel(2, "#mu^{-}");
            hist->GetXaxis()->SetBinLabel(3, "#gamma");
            hist->GetXaxis()->SetBinLabel(4, "#pi^{#pm}");
            hist->GetXaxis()->SetBinLabel(5, "p");
            hist->GetXaxis()->SetBinLabel(6, "Other");
            hist->LabelsOption("h");
        }
    }

    double maxYValue = 0.0;
    for (auto* hist : allHists)
        if (hist && hist->GetMaximum() > maxYValue)
            maxYValue = hist->GetMaximum();

    double yminVal = useLogScale ? 1e-1 : 0;
    double ymaxVal = (ymin == 999 && ymax == 999)
                     ? (useLogScale ? maxYValue * 100.0 : maxYValue * 1.1)
                     : ymax;

    for (auto* hist : allHists)
        if (hist) hist->GetYaxis()->SetRangeUser(yminVal, ymaxVal);

    if (xmin != 999 && xmax != 999) {
        for (auto* hist : allHists)
            if (hist) hist->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    if (currentSignalCombo)     { currentSignalCombo->SetLineWidth(2); currentSignalCombo->SetLineColor(kBlue+1); }
    if (ubooneSignalCombo)      { ubooneSignalCombo->SetLineWidth(2);  ubooneSignalCombo->SetLineColor(kBlue-7); }
    if (nuESignalCombo)         { nuESignalCombo->SetLineWidth(2);     nuESignalCombo->SetLineColor(kAzure+5); }

    if (currentBackgroundCombo) { currentBackgroundCombo->SetLineWidth(2); currentBackgroundCombo->SetLineColor(kOrange+7); }
    if (ubooneBackgroundCombo)  { ubooneBackgroundCombo->SetLineWidth(2);  ubooneBackgroundCombo->SetLineColor(kOrange+6); }
    if (nuEBackgroundCombo)     { nuEBackgroundCombo->SetLineWidth(2);     nuEBackgroundCombo->SetLineColor(kOrange-5); }

    bool first = true;
    auto draw = [&](TH1* hist){ if (hist) { hist->Draw(first ? "hist" : "histsame"); first = false; } };

    if (includeCurrent) {
        draw(currentBackgroundCombo);
        draw(currentSignalCombo);
    }
    if (includeUboone) {
        draw(ubooneBackgroundCombo);
        draw(ubooneSignalCombo);
    }
    if (includeNuE) {
        draw(nuEBackgroundCombo);
        draw(nuESignalCombo);
    }

    for (auto* hist : allHists) {
        if (!hist) continue;
        hist->GetXaxis()->SetTickLength(0.04);
        hist->GetYaxis()->SetTickLength(0.03);
        hist->GetXaxis()->SetTickSize(0.02);
        hist->GetYaxis()->SetTickSize(0.02);
    }

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    if(legendLocation == "topRight"){ Lxmin=0.62; Lymax=0.863; Lxmax=0.87; Lymin=0.863 - 0.12; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.38; Lymin=0.863 - 0.12; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.62; Lymin=0.137; Lxmax=0.87; Lymax=0.137 + 0.12; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymin=0.137; Lxmax=0.38; Lymax=0.137 + 0.12; }

    auto legend = new TLegend(Lxmin, Lymin, Lxmax, Lymax);

    if (includeCurrent) {
        legend->AddEntry(currentSignalCombo, "Signal, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(currentBackgroundCombo, "Background, Pandora BDT SBND (without Refinement)", "f");
    }
    if (includeUboone) {
        legend->AddEntry(ubooneSignalCombo, "Signal, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(ubooneBackgroundCombo, "Background, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    }
    if (includeNuE) {
        legend->AddEntry(nuESignalCombo, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(nuEBackgroundCombo, "Background, Pandora Deep Learning: SBND Nu+E Tune", "f");
    }

    legend->SetTextSize(0.015);
    legend->SetMargin(0.12);
    legend->Draw();

    hists.canvas->SaveAs(filename);

    for (auto* hist : allHists)
        delete hist;
}


void styleDrawSplit(splitHistGroup_struct hists,
                    double ymin, double ymax, double xmin, double xmax,
                    const char* filename, const std::string& legendLocation,
                    int* drawLine = nullptr, int* linePos = nullptr,
                    bool useLogScale = false, bool bestPDGPlot = false){
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    if (useLogScale)
        hists.canvas->SetLogy(1);
    else
        hists.canvas->SetLogy(0);

    std::vector<TH1D*> allHists = {hists.nu_e, hists.NCNpi0, hists.otherNC, hists.CCnumu, hists.CCnue, hists.dirt, hists.nu_eDirt, hists.cosmic, hists.other, hists.nu_eFuzzy};

    if (useLogScale) {
        for (auto* hist : allHists) {
            if (!hist) continue;
            int nBins = hist->GetNbinsX();
            for (int i = 1; i <= nBins; ++i) {
                if (hist->GetBinContent(i) <= 0)
                    hist->SetBinContent(i, 1e-6);
            }
        }
    }
    
    if(bestPDGPlot){
        for(auto* hist : allHists){
            if(!hist) continue;
            hist->GetXaxis()->SetBinLabel(1, "e^{-}");
            hist->GetXaxis()->SetBinLabel(2, "#mu^{-}");
            hist->GetXaxis()->SetBinLabel(3, "#gamma");
            hist->GetXaxis()->SetBinLabel(4, "#pi^{#pm}");
            hist->GetXaxis()->SetBinLabel(5, "p");
            hist->GetXaxis()->SetBinLabel(6, "Other");
            hist->LabelsOption("h");
        }
    }

    for (auto* hist : allHists){
        if(hist){
            hist->SetStats(0);
            hist->GetXaxis()->SetTickLength(0.04);
            hist->GetYaxis()->SetTickLength(0.03);
            hist->GetXaxis()->SetTickSize(0.02);
            hist->GetYaxis()->SetTickSize(0.02);
        }
    }
    

    hists.nu_e->SetLineWidth(2);        hists.nu_e->SetLineColor(TColor::GetColor("#656364"));
    hists.NCNpi0->SetLineWidth(2);      hists.NCNpi0->SetLineColor(TColor::GetColor("#578dff"));
    hists.otherNC->SetLineWidth(2);     hists.otherNC->SetLineColor(TColor::GetColor("#86c8dd"));
    hists.CCnumu->SetLineWidth(2);      hists.CCnumu->SetLineColor(TColor::GetColor("#adad7d"));
    hists.CCnue->SetLineWidth(2);       hists.CCnue->SetLineColor(TColor::GetColor("#c91f16"));
    hists.dirt->SetLineWidth(2);        hists.dirt->SetLineColor(TColor::GetColor("#ff5e02"));
    hists.nu_eDirt->SetLineWidth(2);    hists.nu_eDirt->SetLineColor(TColor::GetColor("#1845fb"));
    hists.cosmic->SetLineWidth(2);      hists.cosmic->SetLineColor(TColor::GetColor("#c849a9"));
    hists.other->SetLineWidth(2);       hists.other->SetLineColor(TColor::GetColor("#ffa90e"));
    hists.nu_eFuzzy->SetLineWidth(2);   hists.other->SetLineColor(TColor::GetColor("#9c9ca1"));

    if((ymin != 999) && (ymax != 999)){
        for(auto* hist : allHists)
            if (hist) hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }

    if((xmin != 999) && (xmax != 999)){
        for(auto* hist : allHists)
            if (hist) hist->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    double maxYValue = 0.0;
    for (auto* hist : allHists)
        if (hist && hist->GetMaximum() > maxYValue)
            maxYValue = hist->GetMaximum();

    //std::cout << "maxYValue = " << maxYValue << std::endl;
    double yminVal = useLogScale ? 1e-1 : 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = useLogScale ? (maxYValue * 100.0) : (maxYValue * 1.1);
        //std::cout << "setting yaxis to " << yminVal << ", " << ymaxVal << std::endl;
        for (auto* hist : allHists)
            if (hist) hist->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    hists.nu_e->Draw("hist");
    hists.NCNpi0->Draw("histsame");
    hists.otherNC->Draw("histsame");
    hists.CCnumu->Draw("histsame");
    hists.CCnue->Draw("histsame");
    hists.dirt->Draw("histsame");
    hists.nu_eDirt->Draw("histsame");
    hists.cosmic->Draw("histsame");
    hists.other->Draw("histsame");
    hists.nu_eFuzzy->Draw("histsame");

    int nEntries = 9;
    double height = std::max(0.03 * nEntries, 0.03);
    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;

    if(legendLocation == "topRight"){ Lxmin=0.72; Lymax=0.863; Lxmax=0.87; Lymin=Lymax - height; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.28; Lymin=Lymax - height; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.72; Lymin=0.137; Lxmax=0.87; Lymax=Lymin + height; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymin=0.137; Lxmax=0.28; Lymax=Lymin + height; }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    //auto legend = new TLegend(0.48, 0.39, 0.87, 0.167);
    legend->AddEntry(hists.nu_e, "#nu+e Elastic Scatter", "f");
    legend->AddEntry(hists.NCNpi0, "NCN#pi^{0}", "f");
    legend->AddEntry(hists.otherNC, "Other NC", "f");
    legend->AddEntry(hists.CCnumu, "CC#nu_{#mu}", "f");
    legend->AddEntry(hists.CCnue, "CC#nu_{e}", "f");
    legend->AddEntry(hists.dirt, "Dirt", "f");
    legend->AddEntry(hists.nu_eDirt, "#nu+e Dirt", "f");
    legend->AddEntry(hists.cosmic, "Cosmic", "f");
    legend->AddEntry(hists.other, "Other", "f");
    legend->AddEntry(hists.nu_eFuzzy, "#nu+e Fuzzy", "f");
    legend->SetTextSize(0.0225);

    legend->SetMargin(0.13);
    legend->Draw();

    hists.canvas->SaveAs(filename);

    // Drawing the histograms as a stack.
    const char* histsTitle = hists.nu_e->GetTitle();
    const char* xAxisTitle = hists.nu_e->GetXaxis()->GetTitle();
    const char* yAxisTitle = hists.nu_e->GetYaxis()->GetTitle();

    std::string stackTitle = std::string(histsTitle) + ";" + xAxisTitle + ";" + yAxisTitle;

    THStack* stack = new THStack("stack", stackTitle.c_str());
    stack->Add(hists.nu_e);
    stack->Add(hists.NCNpi0);
    stack->Add(hists.otherNC);
    stack->Add(hists.CCnumu);
    stack->Add(hists.CCnue);
    stack->Add(hists.dirt);
    stack->Add(hists.nu_eDirt);
    stack->Add(hists.cosmic);
    stack->Add(hists.other);
    stack->Add(hists.nu_eFuzzy);

    hists.canvas->cd();
    hists.canvas->Clear();

    if(useLogScale) hists.canvas->SetLogy(1);
    else hists.canvas->SetLogy(0);

    // Build a frame with your desired axis limits
    double xminFrame = (xmin != 999 ? xmin : hists.nu_e->GetXaxis()->GetXmin());
    double xmaxFrame = (xmax != 999 ? xmax : hists.nu_e->GetXaxis()->GetXmax());
    double yminFrame = (ymin != 999 ? ymin : 1e-6);
    double ymaxFrame = (ymax != 999 ? ymax : stack->GetMaximum()*1.2);
    
    TH1D* frame = new TH1D("frame", stackTitle.c_str(),
                           1, xminFrame, xmaxFrame);

    frame->SetMinimum(yminFrame);
    frame->SetMaximum(ymaxFrame);
    frame->SetTitle(stackTitle.c_str());
    //frame->GetXaxis()->SetTitle(xAxisTitle);
    //frame->GetYaxis()->SetTitle(yAxisTitle);

    frame->SetLineColor(0);
    frame->SetLineWidth(0);
    frame->SetFillStyle(0);     
    frame->SetStats(0);
    //frame->Draw("axis");
    frame->Draw("HIST");

    hists.canvas->Update();

    stack->Draw("hist same");

    gPad->RedrawAxis();

    auto legendStack = new TLegend(Lxmin, Lymax, Lxmax, Lymin);
    legendStack->AddEntry(hists.nu_e, "#nu+e Elastic Scatter", "f");
    legendStack->AddEntry(hists.NCNpi0, "NCN#pi^{0}", "f");
    legendStack->AddEntry(hists.otherNC, "Other NC", "f");
    legendStack->AddEntry(hists.CCnumu, "CC#nu_{#mu}", "f");
    legendStack->AddEntry(hists.CCnue, "CC#nu_{e}", "f");
    legendStack->AddEntry(hists.dirt, "Dirt", "f");
    legendStack->AddEntry(hists.nu_eDirt, "#nu+e Dirt", "f");
    legendStack->AddEntry(hists.cosmic, "Cosmic", "f");
    legendStack->AddEntry(hists.other, "Other", "f");
    legendStack->AddEntry(hists.nu_eFuzzy, "#nu+e Fuzzy", "f");
    legendStack->SetTextSize(0.0225);
    legendStack->SetMargin(0.13);
    legendStack->Draw();

    std::string fname(filename);
    std::string stackedFname;
    size_t pos = fname.rfind(".pdf");
    if(pos != std::string::npos) stackedFname = fname.substr(0, pos) + "_stacked.pdf";
    else stackedFname = fname + "_stacked.pdf";

    hists.canvas->SaveAs(stackedFname.c_str());

}

void allSelectionPlots_macro(){

    TChain *subRunTree = new TChain("ana/SubRun");
    subRunTree->Add("/exp/sbnd/data/users/coackley/analysisFiles_14Jul/*.root");

    TChain *tree = new TChain("ana/NuE");
    tree->Add("/exp/sbnd/data/users/coackley/analysisFiles_14Jul/*.root");

    std::string base_path = "/nashome/c/coackley/outputDirectory/";

    gROOT->SetBatch(true);

    int clearCosmicCut = 1;
    int numPFPs0Cut = 1;
    int numRecoNeutrinosCut = 1;
    int CRUMBSCut = 1;
    int FVCut = 1;
    int primaryPFPCut = 1;
    int razzledPDG11Cut = 1;
    int razzledPDG211Cut = 1;
    int ETheta2Cut = 1;
    int dEdxCut = 1;

    // Minimum number of hits in primary PFP definition
    int primaryPFPMinHitRequirement = 10;

    // Cut values
    double crumbsScoreCut_low = 0.2;
    double crumbsScoreCut_high = 0.76;

    double FVCut_xHigh = 192; 
    double FVCut_xLow = -192; 
    double FVCut_xCentre = 10; 

    double FVCut_yHigh = 194; 
    double FVCut_yLow = -194; 
    
    double FVCut_zHigh = 450; 
    double FVCut_zLow = 6; 
   
    double primaryPFPCutValue = 1;

    double razzled211High_highestEnergyPFP = 0.0125;
    double razzled211Low_highestEnergyPFP = 0.00;
    double razzled11High_highestEnergyPFP = 1;
    double razzled11Low_highestEnergyPFP = 0.875;
    
    double dEdxHigh_highestEnergyPFP = 3.25;
    double dEdxLow_highestEnergyPFP = 0;

    double ETheta2High_highestEnergyPFP = 3.066;
    double ETheta2Low_highestEnergyPFP = 0;

    double actualSignalCount = 0; 

    beforeEventCount_struct eventsBeforeCuts_DLNuE;
    eventCounting_struct eventsAfterCuts_DLNuE;

    if (!gSystem->AccessPathName(base_path.c_str())) {
        gSystem->Exec(Form("rm -rf %s/*", base_path.c_str()));
    }
    gSystem->mkdir(base_path.c_str(), kTRUE);

    std::string tableFileName = base_path + "table.txt";

    std::ofstream clearTableFile(tableFileName, std::ios::trunc);
    if (!clearTableFile.is_open()) {
        std::cerr << "Error: could not open or create " << tableFileName << std::endl;
        return;
    }
    clearTableFile.close();

    if(subRunTree->GetEntries() == 0){
        std::cerr << "SubRun chain has 0 entries - check the file path/pattern" << std::endl;
        return;
    }

    if(tree->GetEntries() == 0){
        std::cerr << "NuE chain has 0 entries - check the file path/pattern" << std::endl;
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

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsNuENuE;

    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;
    double totalPOTNuENuE = 0;

    double cosmicSpillsSumNuE = 0;
    double BNBSpillsSumNuE = 0;
    double SignalSpillsSumNuE = 0;

    double POTSignalNuE_notMissing = 0;
    double POTBNBNuE_notMissing = 0;
    double POTNuENuE_notMissing = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 2 && subRunDLCurrent == 5) BNBSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 1 && subRunDLCurrent == 5) SignalSpillsSumNuE += subRunNumGenEvents;

        std::pair<unsigned int, unsigned int> key = std::make_pair(subRunRun, subRunNumber);

        if(subRunSignal == 1){
            if(subRunDLCurrent == 5 && seenSubRunsSignalNuE.find(key) == seenSubRunsSignalNuE.end()){
                totalPOTSignalNuE += subRunPOT;
                seenSubRunsSignalNuE.insert(key);
            }
            
            if(subRunDLCurrent == 5) POTSignalNuE_notMissing += subRunPOT;
                
        } else if(subRunSignal == 2){
            if(subRunDLCurrent == 5 && seenSubRunsBNBNuE.find(key) == seenSubRunsBNBNuE.end()){
                totalPOTBNBNuE += subRunPOT;
                seenSubRunsBNBNuE.insert(key);
            }
            
            if(subRunDLCurrent == 5) POTBNBNuE_notMissing += subRunPOT;

        } else if(subRunSignal == 4){
            if(subRunDLCurrent == 5 && seenSubRunsNuENuE.find(key) == seenSubRunsNuENuE.end()){
                totalPOTNuENuE += subRunPOT;
                seenSubRunsNuENuE.insert(key);
            } else{
                // Found a subrun with the same run:subrun:event ID as a previous subrun
                //std::cout << "matching run/subrun with signal = 4, key = (" << key.first << ", " << key.second << ")" << std::endl;
            }

            if(subRunDLCurrent == 5) POTNuENuE_notMissing += subRunPOT;
        }
    }

    double targetPOT = 1e21;

    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    double totalPOTNuENuE_notMissing = (POTNuENuE_notMissing + POTBNBNuE_notMissing);

    weights_struct weights;
    weights.signalNuE = targetPOT / POTSignalNuE_notMissing;
    weights.BNBNuE = targetPOT / POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;
    weights.NuENuE = targetPOT / totalPOTNuENuE_notMissing;


    std::cout << "nu+e POT (not missing) = " << POTSignalNuE_notMissing << ", nu+e POT = " << totalPOTSignalNuE << std::endl;
    std::cout << "BNB POT (not missing) = " << POTBNBNuE_notMissing << ", BNB POT = " << totalPOTBNBNuE << std::endl;
    std::cout << "nue POT (not missing) = " << POTNuENuE_notMissing << ", nue POT = " << totalPOTNuENuE << std::endl;
    std::cout << "Weights DLNu+E: BNB = " << weights.BNBNuE << ", Signal = " << weights.signalNuE << ", Intime Cosmics = " << weights.cosmicsNuE << ", nue = " << weights.NuENuE << std::endl;


    UInt_t eventID, runID, subRunID;
    int nuEScatter;
    double nuEScatterTrueVX, nuEScatterTrueVY, nuEScatterTrueVZ;
    double DLCurrent, signal;

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
    std::vector<double> *reco_sliceTrueVX = nullptr;  
    std::vector<double> *reco_sliceTrueVY = nullptr;  
    std::vector<double> *reco_sliceTrueVZ = nullptr;  
    std::vector<double> *reco_sliceNumHits = nullptr;  
    std::vector<double> *reco_sliceNumHitsTruthMatched = nullptr;  
    std::vector<double> *reco_sliceNumTruthHits = nullptr;  
    std::vector<double> *reco_sliceOrigin = nullptr; 
    std::vector<double> *reco_sliceTrueCCNC = nullptr;  
    std::vector<double> *reco_sliceTrueNeutrinoType = nullptr;

    std::vector<double> *truth_particleSliceID = nullptr;  
    std::vector<double> *truth_particlePrimary = nullptr;  
    std::vector<double> *truth_particleVX = nullptr;  
    std::vector<double> *truth_particleVY = nullptr;  
    std::vector<double> *truth_particleVZ = nullptr;  
    std::vector<double> *truth_particlePDG = nullptr;  
    std::vector<double> *truth_particleTrackID = nullptr;  
    std::vector<double> *truth_particleMother = nullptr;  
    std::vector<double> *truth_particleStatusCode = nullptr;  
    
    std::vector<double> *reco_particlePDG = nullptr;
    std::vector<double> *reco_particleIsPrimary = nullptr;
    std::vector<double> *reco_particleVX = nullptr;
    std::vector<double> *reco_particleVY = nullptr;
    std::vector<double> *reco_particleVZ = nullptr;
    std::vector<double> *reco_particleDX = nullptr;
    std::vector<double> *reco_particleDY = nullptr;
    std::vector<double> *reco_particleDZ = nullptr;
    std::vector<double> *reco_particleSliceID = nullptr;
    std::vector<double> *reco_particleBestPlaneEnergy = nullptr;
    std::vector<double> *reco_particleTheta = nullptr;
    std::vector<double> *reco_particleTrackScore = nullptr;
    std::vector<double> *reco_particleCompleteness = nullptr;
    std::vector<double> *reco_particlePurity = nullptr;
    std::vector<double> *reco_particleID = nullptr;
    std::vector<double> *reco_particleTruePDG = nullptr;
    std::vector<double> *reco_particleTrueOrigin = nullptr;
    std::vector<double> *reco_particleTrueInteractionType = nullptr;
    std::vector<double> *reco_particleNumHits = nullptr;
    std::vector<double> *reco_particleNumHitsTruthMatched = nullptr;
    std::vector<double> *reco_particleNumTruthHits = nullptr;
    std::vector<double> *reco_particleClearCosmic = nullptr;  
    std::vector<double> *reco_particlePlane0dEdx = nullptr;  
    std::vector<double> *reco_particlePlane1dEdx = nullptr;  
    std::vector<double> *reco_particlePlane2dEdx = nullptr;  
    std::vector<double> *reco_particleBestPlanedEdx = nullptr;  
    std::vector<double> *reco_particleRazzledPDG11 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG13 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG22 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG211 = nullptr;  
    std::vector<double> *reco_particleRazzledPDG2212 = nullptr;  
    std::vector<double> *reco_particleRazzledBestPDG = nullptr;  
    std::vector<double> *reco_particleShowerLength = nullptr;  
    std::vector<double> *reco_particleShowerOpenAngle = nullptr;  
    std::vector<double> *reco_particleShowerBestPlaneEnergy = nullptr;  
    std::vector<double> *reco_particleTrueVX = nullptr;  
    std::vector<double> *reco_particleTrueVY = nullptr;  
    std::vector<double> *reco_particleTrueVZ = nullptr;  
    std::vector<double> *reco_particleTrueEndX = nullptr;  
    std::vector<double> *reco_particleTrueEndY = nullptr;  
    std::vector<double> *reco_particleTrueEndZ = nullptr;  
  
    std::vector<double> *reco_neutrinoID = nullptr;
    std::vector<double> *reco_neutrinoPDG = nullptr;
    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;
    
    std::vector<double> *angleRecalculationPCASlice_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice_dz = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_dz = nullptr;
    
    std::vector<double> *angleRecalculationPCAPFP_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_dz = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_dz = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dx = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dy = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_dz = nullptr;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("nuEScatter", &nuEScatter);
    tree->SetBranchAddress("nuEScatterTrueVX", &nuEScatterTrueVX);
    tree->SetBranchAddress("nuEScatterTrueVY", &nuEScatterTrueVY);
    tree->SetBranchAddress("nuEScatterTrueVZ", &nuEScatterTrueVZ);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);
    tree->SetBranchAddress("signal", &signal);

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
    tree->SetBranchAddress("reco_sliceTrueVX", &reco_sliceTrueVX);
    tree->SetBranchAddress("reco_sliceTrueVY", &reco_sliceTrueVY);
    tree->SetBranchAddress("reco_sliceTrueVZ", &reco_sliceTrueVZ);
    tree->SetBranchAddress("reco_sliceNumHits", &reco_sliceNumHits);
    tree->SetBranchAddress("reco_sliceNumHitsTruthMatched", &reco_sliceNumHitsTruthMatched);
    tree->SetBranchAddress("reco_sliceNumTruthHits", &reco_sliceNumTruthHits);
    tree->SetBranchAddress("reco_sliceOrigin", &reco_sliceOrigin);
    tree->SetBranchAddress("reco_sliceTrueCCNC", &reco_sliceTrueCCNC);
    tree->SetBranchAddress("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);

    tree->SetBranchAddress("truth_particleSliceID", &truth_particleSliceID);
    tree->SetBranchAddress("truth_particlePrimary", &truth_particlePrimary);
    tree->SetBranchAddress("truth_particleVX", &truth_particleVX);
    tree->SetBranchAddress("truth_particleVY", &truth_particleVY);
    tree->SetBranchAddress("truth_particleVZ", &truth_particleVZ);
    tree->SetBranchAddress("truth_particlePDG", &truth_particlePDG);
    tree->SetBranchAddress("truth_particleTrackID", &truth_particleTrackID);
    tree->SetBranchAddress("truth_particleMother", &truth_particleMother);
    tree->SetBranchAddress("truth_particleStatusCode", &truth_particleStatusCode);

    tree->SetBranchAddress("reco_particlePDG", &reco_particlePDG);
    tree->SetBranchAddress("reco_particleIsPrimary", &reco_particleIsPrimary);
    tree->SetBranchAddress("reco_particleVX", &reco_particleVX);
    tree->SetBranchAddress("reco_particleVY", &reco_particleVY);
    tree->SetBranchAddress("reco_particleVZ", &reco_particleVZ);
    tree->SetBranchAddress("reco_particleDX", &reco_particleDX);
    tree->SetBranchAddress("reco_particleDY", &reco_particleDY);
    tree->SetBranchAddress("reco_particleDZ", &reco_particleDZ);
    tree->SetBranchAddress("reco_particleSliceID", &reco_particleSliceID);
    tree->SetBranchAddress("reco_particleBestPlaneEnergy", &reco_particleBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTheta", &reco_particleTheta);
    tree->SetBranchAddress("reco_particleTrackScore", &reco_particleTrackScore);
    tree->SetBranchAddress("reco_particleCompleteness", &reco_particleCompleteness);
    tree->SetBranchAddress("reco_particlePurity", &reco_particlePurity);
    tree->SetBranchAddress("reco_particleID", &reco_particleID);
    tree->SetBranchAddress("reco_particleTruePDG", &reco_particleTruePDG);
    tree->SetBranchAddress("reco_particleTrueOrigin", &reco_particleTrueOrigin);
    tree->SetBranchAddress("reco_particleTrueInteractionType", &reco_particleTrueInteractionType);
    tree->SetBranchAddress("reco_particleNumHits", &reco_particleNumHits);
    tree->SetBranchAddress("reco_particleNumHitsTruthMatched", &reco_particleNumHitsTruthMatched);
    tree->SetBranchAddress("reco_particleNumTruthHits", &reco_particleNumTruthHits);
    tree->SetBranchAddress("reco_particleClearCosmic", &reco_particleClearCosmic);
    tree->SetBranchAddress("reco_particlePlane0dEdx", &reco_particlePlane0dEdx);
    tree->SetBranchAddress("reco_particlePlane1dEdx", &reco_particlePlane1dEdx);
    tree->SetBranchAddress("reco_particlePlane2dEdx", &reco_particlePlane2dEdx);
    tree->SetBranchAddress("reco_particleBestPlanedEdx", &reco_particleBestPlanedEdx);
    tree->SetBranchAddress("reco_particleRazzledPDG11", &reco_particleRazzledPDG11);
    tree->SetBranchAddress("reco_particleRazzledPDG13", &reco_particleRazzledPDG13);
    tree->SetBranchAddress("reco_particleRazzledPDG22", &reco_particleRazzledPDG22);
    tree->SetBranchAddress("reco_particleRazzledPDG211", &reco_particleRazzledPDG211);
    tree->SetBranchAddress("reco_particleRazzledPDG2212", &reco_particleRazzledPDG2212);
    tree->SetBranchAddress("reco_particleRazzledBestPDG", &reco_particleRazzledBestPDG);
    tree->SetBranchAddress("reco_particleShowerLength", &reco_particleShowerLength);
    tree->SetBranchAddress("reco_particleShowerOpenAngle", &reco_particleShowerOpenAngle);
    tree->SetBranchAddress("reco_particleShowerBestPlaneEnergy", &reco_particleShowerBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTrueVX", &reco_particleTrueVX);
    tree->SetBranchAddress("reco_particleTrueVY", &reco_particleTrueVY);
    tree->SetBranchAddress("reco_particleTrueVZ", &reco_particleTrueVZ);
    tree->SetBranchAddress("reco_particleTrueEndX", &reco_particleTrueEndX);
    tree->SetBranchAddress("reco_particleTrueEndY", &reco_particleTrueEndY);
    tree->SetBranchAddress("reco_particleTrueEndZ", &reco_particleTrueEndZ);
    
    tree->SetBranchAddress("reco_neutrinoID", &reco_neutrinoID);
    tree->SetBranchAddress("reco_neutrinoPDG", &reco_neutrinoPDG);
    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);
    
    tree->SetBranchAddress("angleRecalculationPCASlice_angle", &angleRecalculationPCASlice_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice_sliceID", &angleRecalculationPCASlice_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice_dx", &angleRecalculationPCASlice_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice_dy", &angleRecalculationPCASlice_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice_dz", &angleRecalculationPCASlice_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_angle", &angleRecalculationPCASlice5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_sliceID", &angleRecalculationPCASlice5cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dx", &angleRecalculationPCASlice5cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dy", &angleRecalculationPCASlice5cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_dz", &angleRecalculationPCASlice5cm_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_angle", &angleRecalculationPCASlice10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_sliceID", &angleRecalculationPCASlice10cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dx", &angleRecalculationPCASlice10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dy", &angleRecalculationPCASlice10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_dz", &angleRecalculationPCASlice10cm_dz);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_angle", &angleRecalculationPCASlice15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_sliceID", &angleRecalculationPCASlice15cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dx", &angleRecalculationPCASlice15cm_dx);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dy", &angleRecalculationPCASlice15cm_dy);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_dz", &angleRecalculationPCASlice15cm_dz);
    
    tree->SetBranchAddress("angleRecalculationPCAPFP_angle", &angleRecalculationPCAPFP_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP_pfpID", &angleRecalculationPCAPFP_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dx", &angleRecalculationPCAPFP_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dy", &angleRecalculationPCAPFP_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP_dz", &angleRecalculationPCAPFP_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_angle", &angleRecalculationPCAPFP5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_pfpID", &angleRecalculationPCAPFP5cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dx", &angleRecalculationPCAPFP5cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dy", &angleRecalculationPCAPFP5cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_dz", &angleRecalculationPCAPFP5cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_angle", &angleRecalculationPCAPFP10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_pfpID", &angleRecalculationPCAPFP10cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dx", &angleRecalculationPCAPFP10cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dy", &angleRecalculationPCAPFP10cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_dz", &angleRecalculationPCAPFP10cm_dz);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_angle", &angleRecalculationPCAPFP15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_pfpID", &angleRecalculationPCAPFP15cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dx", &angleRecalculationPCAPFP15cm_dx);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dy", &angleRecalculationPCAPFP15cm_dy);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_dz", &angleRecalculationPCAPFP15cm_dz);

    Long64_t numEntries = tree->GetEntries();

    // Put Plots Here (Before and After All Cuts)
    auto sliceCompletenessBeforeCuts = createHistGroup("sliceCompletenessBeforeCuts", "Slice Completeness (Before Cuts)", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessBeforeCuts_splitDLNuE = createSplitHistGroup("sliceCompletenessBeforeCuts_splitDLNuE", "Slice Completeness (Before Cuts)", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessAfterCuts = createHistGroup("sliceCompletenessAfterCuts", "Slice Completeness (After Cuts)", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessAfterCuts_splitDLNuE = createSplitHistGroup("sliceCompletenessAfterCuts_splitDLNuE", "Slice Completeness (After Cuts)", "Completeness", 102, 0, 1.02);
    
    auto slicePurityBeforeCuts = createHistGroup("slicePurityBeforeCuts", "Slice Purity (Before Cuts)", "Purity", 102, 0, 1);
    auto slicePurityBeforeCuts_splitDLNuE = createSplitHistGroup("slicePurityBeforeCuts_splitDLNuE", "Slice Purity (Before Cuts)", "Purity", 102, 0, 1);
    auto slicePurityAfterCuts = createHistGroup("slicePurityAfterCuts", "Slice Purity (After Cuts)", "Purity", 102, 0, 1);
    auto slicePurityAfterCuts_splitDLNuE = createSplitHistGroup("slicePurityAfterCuts_splitDLNuE", "Slice Purity (After Cuts)", "CRUMBS Score", 102, 0, 1);
    
    auto sliceCRUMBSBeforeCuts = createHistGroup("sliceCRUMBSBeforeCuts", "Slice CRUMBS Score (Before Cuts)", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSBeforeCuts_splitDLNuE = createSplitHistGroup("sliceCRUMBSBeforeCuts_splitDLNuE", "Slice CRUMBS Score (Before Cuts)", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCuts = createHistGroup("sliceCRUMBSAfterCuts", "Slice CRUMBS Score (After Cuts)", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCuts_splitDLNuE = createSplitHistGroup("sliceCRUMBSAfterCuts_splitDLNuE", "Slice CRUMBS Score (After Cuts)", "CRUMBS Score", 25, -1, 1);

    auto sliceNumRecoNeutBeforeCuts = createHistGroup("sliceNumRecoNeutBeforeCuts", "Number of Reco Neutrinos in Slice (Before Cuts)", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutBeforeCuts_splitDLNuE = createSplitHistGroup("sliceNumRecoNeutBeforeCuts_splitDLNuE", "Number of Reco Neutrinos in Slice (Before Cuts)", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterCuts = createHistGroup("sliceNumRecoNeutAfterCuts", "Number of Reco Neutrinos in Slice (After Cuts)", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumRecoNeutAfterCuts_splitDLNuE", "Number of Reco Neutrinos in Slice (After Cuts)", "Number of Reco Neutrinos", 10, 0, 10);
 
    auto sliceNumPFPsBeforeCuts = createHistGroup("sliceNumPFPsBeforeCuts", "Number of PFPs in Slice (Before Cuts)", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsBeforeCuts_splitDLNuE = createSplitHistGroup("sliceNumPFPsBeforeCuts_splitDLNuE", "Number of PFPs in Slice (Before Cuts)", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts = createHistGroup("sliceNumPFPsAfterCuts", "Number of PFPs in Slice (After Cuts)", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPFPsAfterCuts_splitDLNuE", "Number of PFPs in Slice (After Cuts)", "Number of PFPs", 20, 0, 20);
   
    auto sliceNumPrimaryPFPsBeforeCuts = createHistGroup("sliceNumPrimaryPFPsBeforeCuts", "Number of Primary PFPs in Slice (Before Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsBeforeCuts_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsBeforeCuts_splitDLNuE", "Number of Primary PFPs in Slice (Before Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterCuts = createHistGroup("sliceNumPrimaryPFPsAfterCuts", "Number of Primary PFPs in Slice (After Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsAfterCuts_splitDLNuE", "Number of Primary PFPs in Slice (After Cuts)", "Number of Primary PFPs", 20, 0, 20);
    
    auto sliceNumPrimaryPFPsMinHitBeforeCuts = createHistGroup("sliceNumPrimaryPFPsMinHitBeforeCuts", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement (Before Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsMinHitBeforeCuts_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsMinHitBeforeCuts_splitDLNuE", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement (Before Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsMinHitAfterCuts = createHistGroup("sliceNumPrimaryPFPsMinHitAfterCuts", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement (After Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsMinHitAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsMinHitAfterCuts_splitDLNuE", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement (After Cuts)", "Number of Primary PFPs", 20, 0, 20);
    
    auto ERecoSumThetaRecoBeforeCuts = createHistGroup("ERecoSumThetaRecoBeforeCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoBeforeCuts_splitDLNuE = createSplitHistGroup("ERecoSumThetaRecoBeforeCuts_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoAfterCuts = createHistGroup("ERecoSumThetaRecoAfterCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoAfterCuts_splitDLNuE = createSplitHistGroup("ERecoSumThetaRecoAfterCuts_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    
    auto ERecoHighestThetaRecoBeforeCuts = createHistGroup("ERecoHighestThetaRecoBeforeCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoBeforeCuts_splitDLNuE = createSplitHistGroup("ERecoHighestThetaRecoBeforeCuts_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts = createHistGroup("ERecoHighestThetaRecoAfterCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts_splitDLNuE = createSplitHistGroup("ERecoHighestThetaRecoAfterCuts_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    
    auto ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints = createHistGroup("ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoBeforeCuts_splitDLNuE_pfp10cmPoints = createSplitHistGroup("ERecoHighestThetaRecoBeforeCuts_splitDLNuE_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts_pfp10cmPoints = createHistGroup("ERecoHighestThetaRecoAfterCuts_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts_splitDLNuE_pfp10cmPoints = createSplitHistGroup("ERecoHighestThetaRecoAfterCuts_splitDLNuE_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    
    auto dEdxBeforeCuts = createHistGroup("dEdxBeforeCuts", "dE/dx of the PFP in the Slice with the Highest Energy (Before Cuts)", "dE/dx (MeV/cm)", 40, 0, 10);
    auto dEdxBeforeCuts_splitDLNuE = createSplitHistGroup("dEdxBeforeCuts_splitDLNuE", "dE/dx of the PFP in the Slice with the Highest Energy (Before Cuts)", "dE/dx (MeV/cm)", 40, 0, 10);
    auto dEdxAfterCuts = createHistGroup("dEdxAfterCuts", "dE/dx of the PFP in the Slice with the Highest Energy (After Cuts)", "dE/dx (MeV/cm)", 40, 0, 10);
    auto dEdxAfterCuts_splitDLNuE = createSplitHistGroup("dEdxAfterCuts_splitDLNuE", "dE/dx of the PFP in the Slice with the Highest Energy (After Cuts)", "dE/dx (MeV/cm)", 40, 0, 10);
    
    auto razzledPDG11BeforeCuts = createHistGroup("razzledPDG11BeforeCuts", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG11BeforeCuts_splitDLNuE = createSplitHistGroup("razzledPDG11BeforeCuts_splitDLNuE", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG11AfterCuts = createHistGroup("razzledPDG11AfterCuts", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    auto razzledPDG11AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG11AfterCuts_splitDLNuE", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    
    auto razzledPDG13BeforeCuts = createHistGroup("razzledPDG13BeforeCuts", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG13BeforeCuts_splitDLNuE = createSplitHistGroup("razzledPDG13BeforeCuts_splitDLNuE", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG13AfterCuts = createHistGroup("razzledPDG13AfterCuts", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    auto razzledPDG13AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG13AfterCuts_splitDLNuE", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
   
    auto razzledPDG22BeforeCuts = createHistGroup("razzledPDG22BeforeCuts", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG22BeforeCuts_splitDLNuE = createSplitHistGroup("razzledPDG22BeforeCuts_splitDLNuE", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG22AfterCuts = createHistGroup("razzledPDG22AfterCuts", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    auto razzledPDG22AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG22AfterCuts_splitDLNuE", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    
    auto razzledPDG211BeforeCuts = createHistGroup("razzledPDG211BeforeCuts", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG211BeforeCuts_splitDLNuE = createSplitHistGroup("razzledPDG211BeforeCuts_splitDLNuE", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG211AfterCuts = createHistGroup("razzledPDG211AfterCuts", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    auto razzledPDG211AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG211AfterCuts_splitDLNuE", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    
    auto razzledPDG2212BeforeCuts = createHistGroup("razzledPDG2212BeforeCuts", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG2212BeforeCuts_splitDLNuE = createSplitHistGroup("razzledPDG2212BeforeCuts_splitDLNuE", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 80, 0, 1);
    auto razzledPDG2212AfterCuts = createHistGroup("razzledPDG2212AfterCuts", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    auto razzledPDG2212AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG2212AfterCuts_splitDLNuE", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 80, 0, 1);
    
    auto pfpCompletenessBeforeCuts = createHistGroup("pfpCompletenessBeforeCuts", "Completeness of the PFP in the Slice with the Highest Energy (Before Cuts)", "Completeness", 50, 0, 1);
    auto pfpCompletenessBeforeCuts_splitDLNuE = createSplitHistGroup("pfpCompletenessBeforeCuts_splitDLNuE", "Completeness of the PFP in the Slice with the Highest Energy (Before Cuts)", "Completeness", 50, 0, 1);
    auto pfpCompletenessAfterCuts = createHistGroup("pfpCompletenessAfterCuts", "Completeness of the PFP in the Slice with the Highest Energy (After Cuts)", "Completeness", 50, 0, 1);
    auto pfpCompletenessAfterCuts_splitDLNuE = createSplitHistGroup("pfpCompletenessAfterCuts_splitDLNuE", "Completeness of the PFP in the Slice with the Highest Energy (After Cuts)", "Completeness", 50, 0, 1);
    
    auto pfpPurityBeforeCuts = createHistGroup("pfpPurityBeforeCuts", "Purity of the PFP in the Slice with the Highest Energy (Before Cuts)", "Purity", 50, 0, 1);
    auto pfpPurityBeforeCuts_splitDLNuE = createSplitHistGroup("pfpPurityBeforeCuts_splitDLNuE", "Purity of the PFP in the Slice with the Highest Energy (Before Cuts)", "Purity", 50, 0, 1);
    auto pfpPurityAfterCuts = createHistGroup("pfpPurityAfterCuts", "Purity of the PFP in the Slice with the Highest Energy (After Cuts)", "Purity", 50, 0, 1);
    auto pfpPurityAfterCuts_splitDLNuE = createSplitHistGroup("pfpPurityAfterCuts_splitDLNuE", "Purity of the PFP in the Slice with the Highest Energy (After Cuts)", "Purity", 50, 0, 1);
   
    auto angleDifferenceSignalBeforeCuts = createHistGroup("angleDifferenceSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFPSignalBeforeCuts = createHistGroup("angleDifferencePCAPFPSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFP5cmSignalBeforeCuts = createHistGroup("angleDifferencePCAPFP5cmSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP 5 cm) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFP10cmSignalBeforeCuts = createHistGroup("angleDifferencePCAPFP10cmSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP 10 cm) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFP15cmSignalBeforeCuts = createHistGroup("angleDifferencePCAPFP15cmSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP 15 cm) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASliceSignalBeforeCuts = createHistGroup("angleDifferencePCASliceSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASlice5cmSignalBeforeCuts = createHistGroup("angleDifferencePCASlice5cmSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice 5 cm) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASlice10cmSignalBeforeCuts = createHistGroup("angleDifferencePCASlice10cmSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice 10 cm) (Before Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASlice15cmSignalBeforeCuts = createHistGroup("angleDifferencePCASlice15cmSignalBeforeCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice 15 cm) (Before Cuts)", "Angle (degrees)", 16, 0, 16);    
        
    auto angleDifferenceSignalAfterCuts = createHistGroup("angleDifferenceSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFPSignalAfterCuts = createHistGroup("angleDifferencePCAPFPSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFP5cmSignalAfterCuts = createHistGroup("angleDifferencePCAPFP5cmSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP 5 cm) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFP10cmSignalAfterCuts = createHistGroup("angleDifferencePCAPFP10cmSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP 10 cm) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCAPFP15cmSignalAfterCuts = createHistGroup("angleDifferencePCAPFP15cmSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using PFP 15 cm) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASliceSignalAfterCuts = createHistGroup("angleDifferencePCASliceSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASlice5cmSignalAfterCuts = createHistGroup("angleDifferencePCASlice5cmSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice 5 cm) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASlice10cmSignalAfterCuts = createHistGroup("angleDifferencePCASlice10cmSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice 10 cm) (After Cuts)", "Angle (degrees)", 16, 0, 16);
    auto angleDifferencePCASlice15cmSignalAfterCuts = createHistGroup("angleDifferencePCASlice15cmSignalAfterCuts", "Angle Difference between True and Reconstructed Recoil Electron (PCA using Slice 15 cm) (After Cuts)", "Angle (degrees)", 16, 0, 16);
 
    auto recoVXBeforeCuts = createHistGroup("recoVXBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXBeforeCuts_splitDLNuE = createSplitHistGroup("recoVXBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterCuts = createHistGroup("recoVXAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterCuts_splitDLNuE = createSplitHistGroup("recoVXAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm) ", 202, -202, 202);
   
    auto recoVYBeforeCuts = createHistGroup("recoVYBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYBeforeCuts_splitDLNuE = createSplitHistGroup("recoVYBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterCuts = createHistGroup("recoVYAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterCuts_splitDLNuE = createSplitHistGroup("recoVYAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm) ", 204, -204, 204);
    
    auto recoVZBeforeCuts = createHistGroup("recoVZBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZBeforeCuts_splitDLNuE = createSplitHistGroup("recoVZBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterCuts = createHistGroup("recoVZAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterCuts_splitDLNuE = createSplitHistGroup("recoVZAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 255, 0, 510);
    
    auto recoVXSmallerBinsBeforeCuts = createHistGroup("recoVXSmallerBinsBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsBeforeCuts_splitDLNuE = createSplitHistGroup("recoVXSmallerBinsBeforeCutsBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterCuts = createHistGroup("recoVXSmallerBinsAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("recoVXSmallerBinsAfterCutsAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm) ", 808, -202, 202);
    
    auto recoVYSmallerBinsBeforeCuts = createHistGroup("recoVYSmallerBinsBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsBeforeCuts_splitDLNuE = createSplitHistGroup("recoVYSmallerBinsBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterCuts = createHistGroup("recoVYSmallerBinsAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("recoVYSmallerBinsAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm) ", 816, -204, 204);
    
    auto recoVZSmallerBinsBeforeCuts = createHistGroup("recoVZSmallerBinsBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsBeforeCuts_splitDLNuE = createSplitHistGroup("recoVZSmallerBinsBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterCuts = createHistGroup("recoVZSmallerBinsAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("recoVZSmallerBinsAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 1020, 0, 510);
    
    auto recoVXLowBeforeCuts = createHistGroup("recoVXLowBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, 182);
    auto recoVXLowBeforeCuts_splitDLNuE = createSplitHistGroup("recoVXLowBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterCuts = createHistGroup("recoVXLowAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterCuts_splitDLNuE = createSplitHistGroup("recoVXLowAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm) ", 40, -202, -182);
    
    auto recoVYLowBeforeCuts = createHistGroup("recoVYLowBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowBeforeCuts_splitDLNuE = createSplitHistGroup("recoVYLowBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterCuts = createHistGroup("recoVYLowAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterCuts_splitDLNuE = createSplitHistGroup("recoVYLowAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm) ", 40, -204, -184);
    
    auto recoVZLowBeforeCuts = createHistGroup("recoVZLowBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowBeforeCuts_splitDLNuE = createSplitHistGroup("recoVZLowBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterCuts = createHistGroup("recoVZLowAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterCuts_splitDLNuE = createSplitHistGroup("recoVZLowAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 40, 0, 20);
    
    auto recoVXHighBeforeCuts = createHistGroup("recoVXHighBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighBeforeCuts_splitDLNuE = createSplitHistGroup("recoVXHighBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterCuts = createHistGroup("recoVXHighAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterCuts_splitDLNuE = createSplitHistGroup("recoVXHighAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm)", 40, 182, 202);
    
    auto recoVYHighBeforeCuts = createHistGroup("recoVYHighBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighBeforeCuts_splitDLNuE = createSplitHistGroup("recoVYHighBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterCuts = createHistGroup("recoVYHighAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterCuts_splitDLNuE = createSplitHistGroup("recoVYHighAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm)", 40, 184, 204);
    
    auto recoVZHighBeforeCuts = createHistGroup("recoVZHighBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighBeforeCuts_splitDLNuE = createSplitHistGroup("recoVZHighBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterCuts = createHistGroup("recoVZHighAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterCuts_splitDLNuE = createSplitHistGroup("recoVZHighAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 40, 480, 510);
   
    auto energyAsymmetryBeforeCuts = createHistGroup("energyAsymmetryBeforeCuts", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (Before Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryBeforeCuts_splitDLNuE = createSplitHistGroup("energyAsymmetryBeforeCut_splitDLNuE", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (Before Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts = createHistGroup("energyAsymmetryAfterCuts", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (After Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts_splitDLNuE = createSplitHistGroup("energyAsymmetryAfterCut_splitDLNuE", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (After Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
     

    // Put Plots Here (Dedicated After Previous Cut)
    auto sliceNumPFPsAfterCuts = createHistGroup("sliceNumPFPsAfterCuts", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPFPsAfterCuts_splitDLNuE", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts = createHistGroup("sliceNumPFPsAfterCuts", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPFPsAfterCuts_splitDLNuE", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    
    auto sliceNumRecoNeutAfterNumPFPCut = createHistGroup("sliceNumRecoNeutAfterNumPFPCut", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterNumPFPCut_splitDLNuE = createSplitHistGroup("sliceNumRecoNeutAfterNumPFPCut_splitDLNuE", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterNumNeutrinoCut = createHistGroup("sliceNumRecoNeutAfterNumNeutrinoCut", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterNumNeutrinoCut_splitDLNuE = createSplitHistGroup("sliceNumRecoNeutAfterNumNeutrinoCut_splitDLNuE", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    
    auto sliceCRUMBSAfterNumNeutrinoCut = createHistGroup("sliceCRUMBSAfterNumNeutrinoCut", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterNumNeutrinoCut_splitDLNuE = createSplitHistGroup("sliceCRUMBSAfterNumNeutrinoCut_splitDLNuE", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCRUMBSCut = createHistGroup("sliceCRUMBSAfterCRUMBSCut", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceCRUMBSAfterCRUMBSCut_splitDLNuE", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    
    auto recoVXAfterCRUMBSCut = createHistGroup("recoVXAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVXAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterFVCut = createHistGroup("recoVXAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterFVCut_splitDLNuE = createSplitHistGroup("recoVXAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
   
    auto recoVYAfterCRUMBSCut = createHistGroup("recoVYAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVYAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterFVCut = createHistGroup("recoVYAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterFVCut_splitDLNuE = createSplitHistGroup("recoVYAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    
    auto recoVZAfterCRUMBSCut = createHistGroup("recoVZAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVZAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterFVCut = createHistGroup("recoVZAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterFVCut_splitDLNuE = createSplitHistGroup("recoVZAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    
    auto recoVXSmallerBinsAfterCRUMBSCut = createHistGroup("recoVXSmallerBinsAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVXSmallerBinsAfterCRUMBSCutAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterFVCut = createHistGroup("recoVXSmallerBinsAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterFVCut_splitDLNuE = createSplitHistGroup("recoVXSmallerBinsAfterFVCutAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    
    auto recoVYSmallerBinsAfterCRUMBSCut = createHistGroup("recoVYSmallerBinsAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVYSmallerBinsAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterFVCut = createHistGroup("recoVYSmallerBinsAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterFVCut_splitDLNuE = createSplitHistGroup("recoVYSmallerBinsAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    
    auto recoVZSmallerBinsAfterCRUMBSCut = createHistGroup("recoVZSmallerBinsAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVZSmallerBinsAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterFVCut = createHistGroup("recoVZSmallerBinsAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterFVCut_splitDLNuE = createSplitHistGroup("recoVZSmallerBinsAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    
    auto recoVXLowAfterCRUMBSCut = createHistGroup("recoVXLowAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVXLowAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterFVCut = createHistGroup("recoVXLowAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterFVCut_splitDLNuE = createSplitHistGroup("recoVXLowAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    
    auto recoVYLowAfterCRUMBSCut = createHistGroup("recoVYLowAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVYLowAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterFVCut = createHistGroup("recoVYLowAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterFVCut_splitDLNuE = createSplitHistGroup("recoVYLowAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    
    auto recoVZLowAfterCRUMBSCut = createHistGroup("recoVZLowAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVZLowAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterFVCut = createHistGroup("recoVZLowAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterFVCut_splitDLNuE = createSplitHistGroup("recoVZLowAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    
    auto recoVXHighAfterCRUMBSCut = createHistGroup("recoVXHighAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVXHighAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterFVCut = createHistGroup("recoVXHighAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterFVCut_splitDLNuE = createSplitHistGroup("recoVXHighAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    
    auto recoVYHighAfterCRUMBSCut = createHistGroup("recoVYHighAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVYHighAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterFVCut = createHistGroup("recoVYHighAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterFVCut_splitDLNuE = createSplitHistGroup("recoVYHighAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    
    auto recoVZHighAfterCRUMBSCut = createHistGroup("recoVZHighAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("recoVZHighAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterFVCut = createHistGroup("recoVZHighAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterFVCut_splitDLNuE = createSplitHistGroup("recoVZHighAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
 
    auto sliceNumPrimaryPFPsAfterFVCut = createHistGroup("sliceNumPrimaryPFPsAfterFVCut", "Number of Primary PFPs in Slice", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterFVCut_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsAfterFVCut_splitDLNuE", "Number of Primary PFPs in Slice", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterPrimaryPFPCut = createHistGroup("sliceNumPrimaryPFPsAfterPrimaryPFPCut", "Number of Primary PFPs in Slice", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterPrimaryPFPCut_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsAfterPrimaryPFPCut_splitDLNuE", "Number of Primary PFPs in Slice", "Number of Primary PFPs", 20, 0, 20);
    
    auto sliceNumPrimaryPFPsMinHitAfterFVCut = createHistGroup("sliceNumPrimaryPFPsMinHitAfterFVCut", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsMinHitAfterFVCut_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsMinHitAfterFVCut_splitDLNuE", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut = createHistGroup("sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut_splitDLNuE", "Number of Primary PFPs in Slice with Number of Hits > Hit Requirement", "Number of Primary PFPs", 20, 0, 20);
    
    auto razzledPDG11AfterPrimaryPFPCut = createHistGroup("razzledPDG11AfterPrimaryPFPCut", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    auto razzledPDG11AfterPrimaryPFPCut_splitDLNuE = createSplitHistGroup("razzledPDG11AfterPrimaryPFPCut_splitDLNuE", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    auto razzledPDG11AfterRazzled11Cut = createHistGroup("razzledPDG11AfterRazzled11Cut", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    auto razzledPDG11AfterRazzled11Cut_splitDLNuE = createSplitHistGroup("razzledPDG11AfterRazzled11Cut_splitDLNuE", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    
    auto razzledPDG211AfterRazzled11Cut = createHistGroup("razzledPDG211AfterRazzled11Cut", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    auto razzledPDG211AfterRazzled11Cut_splitDLNuE = createSplitHistGroup("razzledPDG211AfterRazzled11Cut_splitDLNuE", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    auto razzledPDG211AfterRazzled211Cut = createHistGroup("razzledPDG211AfterRazzled211Cut", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    auto razzledPDG211AfterRazzled211Cut_splitDLNuE = createSplitHistGroup("razzledPDG211AfterRazzled211Cut_splitDLNuE", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy", "Score", 80, 0, 1);
    
    auto ERecoHighestThetaRecoAfterRazzled211Cut = createHistGroup("ERecoHighestThetaRecoAfterRazzled211Cut", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE = createSplitHistGroup("ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterETheta2Cut = createHistGroup("ERecoHighestThetaRecoAfterETheta2Cut", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE = createSplitHistGroup("ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    
    auto ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints = createHistGroup("ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE_pfp10cmPoints = createSplitHistGroup("ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterETheta2Cut_pfp10cmPoints = createHistGroup("ERecoHighestThetaRecoAfterETheta2Cut_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE_pfp10cmPoints = createSplitHistGroup("ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE_pfp10cmPoints", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    
    auto dEdxAfterETheta2Cut = createHistGroup("dEdxAfterETheta2Cut", "dE/dx of the PFP in the Slice with the Highest Energy", "dE/dx (MeV/cm)", 40, 0, 10);
    auto dEdxAfterETheta2Cut_splitDLNuE = createSplitHistGroup("dEdxAfterETheta2Cut_splitDLNuE", "dE/dx of the PFP in the Slice with the Highest Energy", "dE/dx (MeV/cm)", 40, 0, 10);
    auto dEdxAfterdEdxCut = createHistGroup("dEdxAfterdEdxCut", "dE/dx of the PFP in the Slice with the Highest Energy", "dE/dx (MeV/cm)", 40, 0, 10);
    auto dEdxAfterdEdxCut_splitDLNuE = createSplitHistGroup("dEdxAfterdEdxCut_splitDLNuE", "dE/dx of the PFP in the Slice with the Highest Energy", "dE/dx (MeV/cm)", 40, 0, 10);

    double numEvents_DLNuECosmic = 0;
    double numEvents_DLNuEBNB = 0;
    double numEvents_DLNuESignal = 0;
    double numEvents_DLNuENuE = 0;

    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);

        if(DLCurrent == 5 && signal == 4) numEvents_DLNuENuE++;
        if(DLCurrent == 5 && signal == 3) numEvents_DLNuECosmic++;
        if(DLCurrent == 5 && signal == 2) numEvents_DLNuEBNB++;
        if(DLCurrent == 5 && signal == 1) numEvents_DLNuESignal++;

        // If there is a true recoil electron in the event, look at it
        recoilElectron_struct recoilElectron; 
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->size() > 1) std::cout << "More than 1 true recoil electron in event!" << std::endl;
            if(truth_recoilElectronPDG->at(i) != -999999){
                // There is a true recoil electron in the event
                recoilElectron.energy = truth_recoilElectronEnergy->at(i);
                recoilElectron.angle = truth_recoilElectronAngle->at(i);
                recoilElectron.dx = truth_recoilElectronDX->at(i);
                recoilElectron.dy = truth_recoilElectronDY->at(i);
                recoilElectron.dz = truth_recoilElectronDZ->at(i);
            } else if(truth_recoilElectronPDG->size() == 1 && truth_recoilElectronPDG->at(i) == -999999){
                // There is no recoil electron in the event
                recoilElectron.energy = -999999;
                recoilElectron.angle = -999999;
                recoilElectron.dx = -999999;
                recoilElectron.dy = -999999;
                recoilElectron.dz = -999999;
            }
        }

        // Check whether there is a true nu+e elastic scatter in the event
        if(nuEScatter == 1 && signal == 1 && DLCurrent == 5){
            // This is an event with a nu+e elastic scatter in it (from the signal files)
            if(recoilElectron.energy == -999999) std::cout << "No True Recoil Electron Energy Found in Nu+E Elastic Scatter Event" << std::endl;
            //std::cout << "energy of true recoil electron in nu+e event = " << recoilElectron.energy << std::endl;

            if(recoilElectron.energy > 150){
                // nu+e elastic scatter must have true recoil electron with energy > 150 MeV to be classed as signal
                if(FVCut == 0 && (((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))){
                    // The true neutrino interaction is within the active volume
                    trueSignal = 1;
                } else if(FVCut == 1 && (((nuEScatterTrueVX > FVCut_xLow) && (nuEScatterTrueVX < FVCut_xHigh) && (std::abs(nuEScatterTrueVX) > FVCut_xCentre)) && ((nuEScatterTrueVY > FVCut_yLow) && (nuEScatterTrueVY < FVCut_yHigh)) && ((nuEScatterTrueVZ > FVCut_zLow) && (nuEScatterTrueVZ < FVCut_zHigh)))){
                    // The true neutrino interaction is within the fiducial volume
                    trueSignal = 1;
                }
            }
        }

        // This event is a true nu+e elastic scattering event from the dedicated nu+e elastic scatter files
        if(trueSignal == 1) actualSignalCount += weights.signalNuE;
        
        // Looking at the reco slices
        if(reco_sliceID->size() == 0) continue;

        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) == -999999) continue;

            double sliceRecoVX = -999999;
            double sliceRecoVY = -999999;
            double sliceRecoVZ = -999999;
            double numRecoNeutrinos = 0;

            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    numRecoNeutrinos++;
                    sliceRecoVX = reco_neutrinoVX->at(recoNeut); 
                    sliceRecoVY = reco_neutrinoVY->at(recoNeut); 
                    sliceRecoVZ = reco_neutrinoVZ->at(recoNeut); 
                }
            }

            // Assigning a category to the slices
            // 0 = cosmic, 1 = signal, 2 = signal fuzzy, 3 = bnb, 4 = bnb fuzzy
            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                // This is a cosmic slice
                sliceCategoryPlottingMacro = 0;
                //std::cout << "Cosmic Slice: sliceCategoryPlottingMacro = 0" << std::endl;
            } else if(reco_sliceOrigin->at(slice) == 1){
                // This is a nu+e elastic scatter slice
                if(reco_sliceCompleteness->at(slice) > 0.5 && recoilElectron.energy > 150){
                    // Slice must have completeness > 0.5 and have nu+e elastic scatter it comes from has true recoil electron energy > 150 MeV
                    if(FVCut == 0 && (reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVX->at(slice) > -201.3) && (reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVY->at(slice) > -203.8) && (reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5)){
                        // True neutrino vertex is within the active volume (this is the signal definition if we aren't using the FV cuts
                        // -201.3 < x < 201.3, -203.8 < y < 203.8, 0 < z < 509.5
                        sliceCategoryPlottingMacro = 1;
                        //std::cout << "Nu+E Slice: sliceCategoryPlottingMacro = 1" << std::endl;
                    } else if(FVCut == 1 && (reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow && std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre) && (reco_sliceTrueVY->at(slice) < FVCut_yHigh && reco_sliceTrueVY->at(slice) > FVCut_yLow) && (reco_sliceTrueVZ->at(slice) < FVCut_zHigh && reco_sliceTrueVZ->at(slice) > FVCut_zLow)){
                        // True neutrino vertex is within the FV (this is signal definition if we are using the FV cut)
                        sliceCategoryPlottingMacro = 1;
                    } else{
                        sliceCategoryPlottingMacro = 2;
                    }
                } else{
                    sliceCategoryPlottingMacro = 2;
                    //std::cout << "Nu+E Fuzzy Slice: sliceCategoryPlottingMacro = 2" << std::endl;
                }
            } else if(reco_sliceOrigin->at(slice) == 3){
                // This is a BNB slice
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    if(reco_sliceTrueCCNC->at(slice) == 0 && reco_sliceTrueNeutrinoType->at(slice) == 12){
                        // This is a slice from a nue event
                        sliceCategoryPlottingMacro = 5;
                    } else{
                        // This is a BNB event (not a nue event)
                        sliceCategoryPlottingMacro = 3;
                        //std::cout << "BNB Slice: sliceCategoryPlottingMacro = 3" << std::endl;
                    }
                } else{
                    if(reco_sliceTrueCCNC->at(slice) == 0 && reco_sliceTrueNeutrinoType->at(slice) == 12){
                        sliceCategoryPlottingMacro = 6;
                    } else{
                        sliceCategoryPlottingMacro = 4;
                        //std::cout << "BNB Fuzzy Slice: sliceCategoryPlottingMacro = 4" << std::endl;
                    }
                }
            }

            double weight = 0;
            if(signal == 1 && DLCurrent == 5) weight = weights.signalNuE;
            if(signal == 2 && DLCurrent == 5 && sliceCategoryPlottingMacro != 5 && sliceCategoryPlottingMacro != 6) weight = weights.BNBNuE;
            if((signal == 2 || signal == 4) && DLCurrent == 5 && (sliceCategoryPlottingMacro == 5 || sliceCategoryPlottingMacro == 6)) weight = weights.NuENuE;
            if(signal == 3 && DLCurrent == 5) weight = weights.cosmicsNuE;


            // Assigning a interaction category to the slices
            // Event types: Cosmic = 0, nu+e scatter = 1, NC Npi0 = 2, other NC = 3, CC numu = 4, CC nue = 5, Dirt = 6, Dirt nu+e = 7
            // Other = 8, Fuzzy nu+e = 9, Dirt nuE = 10
            int sliceInteractionType = -999999;
            if(reco_sliceOrigin->at(slice) != 0){
                // This is a slice that isn't truth-matched to a cosmic
                if(reco_sliceOrigin->at(slice) == 1){
                    // This is a slice that is truth-matched to a nu+e elastic scatter
                    if(reco_sliceCompleteness->at(slice) > 0.5 && recoilElectron.energy > 150){
                        // This is a slice that is truth-matched to a nu+e elastic scatter AND has completeness > 0.5 (signal slice)
                        // The nu+e elastic scatter it is truth matched to must have true recoil electron with energy > 150 MeV
                        //if(FVCut == 0 && (sliceRecoVX < 201.3 && sliceRecoVX > -201.3) && (sliceRecoVY < 203.8 && sliceRecoVY > -203.8) && (sliceRecoVZ > 0 && sliceRecoVZ < 509.5)){
                        if(FVCut == 0 && (reco_sliceTrueVX->at(slice) > -201.3 && reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVY->at(slice) > -203.8 && reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5)){
                            // Interaction happened inside the TPC (in truth)
                            sliceInteractionType = 1;
                        //} else if(FVCut == 1 && (sliceRecoVX < FVCut_xHigh && sliceRecoVX > FVCut_xLow && std::abs(sliceRecoVX) > FVCut_xCentre) && (sliceRecoVY < FVCut_yHigh && sliceRecoVY > FVCut_yLow) && (sliceRecoVZ < FVCut_zHigh && sliceRecoVZ > FVCut_zLow)){
                        } else if(FVCut == 1 && ((reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow && std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre) && (reco_sliceTrueVY->at(slice) > FVCut_yLow && reco_sliceTrueVY->at(slice) < FVCut_yHigh) && (reco_sliceTrueVZ->at(slice) > FVCut_zLow && reco_sliceTrueVZ->at(slice) < FVCut_zHigh))){
                            // Interaction happened inside the FV (in truth)
                            sliceInteractionType = 1;
                        } else{
                            sliceInteractionType = 7;
                        }
                    } else{
                        // This is a slice that is truth-mathced to a nu+e elastic scatter with completeness < 0.5 (not signal)
                        sliceInteractionType = 9;
                    }
                } else if(reco_sliceOrigin->at(slice) == 3){
                    // This is a slice that is truth-matched to a beam neutrino that isn't a nu+e elastic scatter
                    if((FVCut == 0 && (reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVX->at(slice) > -201.3) && (reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVY->at(slice) > -203.8) && (reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5)) || (FVCut == 1 && (reco_sliceTrueVX->at(slice) < FVCut_xHigh && reco_sliceTrueVX->at(slice) > FVCut_xLow && std::abs(reco_sliceTrueVX->at(slice)) > FVCut_xCentre) && (reco_sliceTrueVY->at(slice) < FVCut_yHigh && reco_sliceTrueVY->at(slice) > FVCut_yLow) && (reco_sliceTrueVZ->at(slice) < FVCut_zHigh && reco_sliceTrueVZ->at(slice) > FVCut_zLow))){
                        // Interaction happened inside the TPC/FV (in truth)
                        if(reco_sliceTrueCCNC->at(slice) == 0){
                            // This is a CC process
                            if(reco_sliceTrueNeutrinoType->at(slice) == 12){
                                // This is a CC nue
                                sliceInteractionType = 5;
                            } else if(reco_sliceTrueNeutrinoType->at(slice) == 14){
                                // This is a CC numu
                                sliceInteractionType = 4;
                            }
                        } else if(reco_sliceTrueCCNC->at(slice) == 1){
                            // This is an NC process
                            int neutralPion = 0; // Number of neutral pions with status code = 1
                            for(size_t trueParticle = 0; trueParticle < truth_particleSliceID->size(); trueParticle++){
                                if(truth_particleSliceID->at(trueParticle) == reco_sliceID->at(slice)){
                                    // True particle is in the slice
                                    if(truth_particleStatusCode->at(trueParticle) == 1){
                                        // True particle has status code of 1 -> Tracked by GENIE
                                        if(truth_particlePDG->at(trueParticle) == 111){
                                            // True particle is a neutral pion
                                            neutralPion++;
                                        }
                                    }
                                }
                            }

                            if(neutralPion > 0){
                                // Slice has a true pi0 in it
                                // This is an NC Npi0 process
                                sliceInteractionType = 2;
                            } else{
                                // This is an NC other process
                                sliceInteractionType = 3;
                            }
                        }
                    } else{
                        // Interaction happened outside the TPC/FV - Dirt event
                        sliceInteractionType = 6;
                    }
                }
            } else{
                // This is a cosmic events
                sliceInteractionType = 0;
            }
            
            if(sliceInteractionType == -999999){
                sliceInteractionType = 8;
            }

            double summedEnergy_beforeCuts = 0;
            double numPFPsSlice_beforeCuts = 0;
            double numPrimaryPFPsSlice_beforeCuts = 0;
            double numPrimaryPFPsMinHitSlice_beforeCuts = 0;

            // Looping through all PFPs in slice before clear cosmic cut applied            
            highestEnergyPFP_struct highestEnergyPFP_beforeCuts;
            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                    // PFP is in the slice
                    numPFPsSlice_beforeCuts++;
                    if(reco_particleIsPrimary->at(pfp) == 1){
                        numPrimaryPFPsSlice_beforeCuts++;
                        if(reco_particleNumHits->at(pfp) >= primaryPFPMinHitRequirement) numPrimaryPFPsMinHitSlice_beforeCuts++;
                    }
                    
                    summedEnergy_beforeCuts += reco_particleBestPlaneEnergy->at(pfp);

                    if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergyPFP_beforeCuts.energy){
                        highestEnergyPFP_beforeCuts.energy = reco_particleBestPlaneEnergy->at(pfp);
                        highestEnergyPFP_beforeCuts.theta = reco_particleTheta->at(pfp);
                        highestEnergyPFP_beforeCuts.PFPID = reco_particleID->at(pfp);
                        highestEnergyPFP_beforeCuts.dx = reco_particleDX->at(pfp);
                        highestEnergyPFP_beforeCuts.dy = reco_particleDY->at(pfp);
                        highestEnergyPFP_beforeCuts.dz = reco_particleDZ->at(pfp);
                        highestEnergyPFP_beforeCuts.vx = reco_particleVX->at(pfp);
                        highestEnergyPFP_beforeCuts.vy = reco_particleVY->at(pfp);
                        highestEnergyPFP_beforeCuts.vz = reco_particleVZ->at(pfp);
                        highestEnergyPFP_beforeCuts.completeness = reco_particleCompleteness->at(pfp);
                        highestEnergyPFP_beforeCuts.purity = reco_particlePurity->at(pfp);
                        highestEnergyPFP_beforeCuts.trackscore = reco_particleTrackScore->at(pfp);
                        highestEnergyPFP_beforeCuts.primary = reco_particleIsPrimary->at(pfp);
                        highestEnergyPFP_beforeCuts.truePDG = reco_particleTruePDG->at(pfp);
                        highestEnergyPFP_beforeCuts.trueOrigin = reco_particleTrueOrigin->at(pfp);
                        highestEnergyPFP_beforeCuts.trueInt = reco_particleTrueInteractionType->at(pfp);
                        highestEnergyPFP_beforeCuts.bestPlanedEdx = reco_particleBestPlanedEdx->at(pfp);
                        highestEnergyPFP_beforeCuts.razzledPDG11 = reco_particleRazzledPDG11->at(pfp);
                        highestEnergyPFP_beforeCuts.razzledPDG13 = reco_particleRazzledPDG13->at(pfp);
                        highestEnergyPFP_beforeCuts.razzledPDG22 = reco_particleRazzledPDG22->at(pfp);
                        highestEnergyPFP_beforeCuts.razzledPDG211 = reco_particleRazzledPDG211->at(pfp);
                        highestEnergyPFP_beforeCuts.razzledPDG2212 = reco_particleRazzledPDG2212->at(pfp);
                        highestEnergyPFP_beforeCuts.razzledBestPDG = reco_particleRazzledBestPDG->at(pfp);
                        highestEnergyPFP_beforeCuts.trueVX = reco_particleTrueVX->at(pfp);
                        highestEnergyPFP_beforeCuts.trueVY = reco_particleTrueVY->at(pfp);
                        highestEnergyPFP_beforeCuts.trueVZ = reco_particleTrueVZ->at(pfp);
                        highestEnergyPFP_beforeCuts.trueEndX = reco_particleTrueEndX->at(pfp);
                        highestEnergyPFP_beforeCuts.trueEndY = reco_particleTrueEndY->at(pfp);
                        highestEnergyPFP_beforeCuts.trueEndZ = reco_particleTrueEndZ->at(pfp);
                        highestEnergyPFP_beforeCuts.numHits = reco_particleNumHits->at(pfp);
                        highestEnergyPFP_beforeCuts.clearCosmic = reco_particleClearCosmic->at(pfp);

                        if(highestEnergyPFP_beforeCuts.trueVX != -999999 && highestEnergyPFP_beforeCuts.trueVY != -999999 && highestEnergyPFP_beforeCuts.trueVZ != -999999 && highestEnergyPFP_beforeCuts.trueEndX != -999999 && highestEnergyPFP_beforeCuts.trueEndY != -999999 && highestEnergyPFP_beforeCuts.trueEndZ != -999999){
                            double xCoordDiff_length = (highestEnergyPFP_beforeCuts.trueVX - highestEnergyPFP_beforeCuts.trueEndX);
                            double yCoordDiff_length = (highestEnergyPFP_beforeCuts.trueVY - highestEnergyPFP_beforeCuts.trueEndY);
                            double zCoordDiff_length = (highestEnergyPFP_beforeCuts.trueVZ - highestEnergyPFP_beforeCuts.trueEndZ);
                            highestEnergyPFP_beforeCuts.trueLength = std::sqrt((xCoordDiff_length * xCoordDiff_length) + (yCoordDiff_length * yCoordDiff_length) + (zCoordDiff_length * zCoordDiff_length));
                        }
                    }
                }
            } // End of looping through all PFPs in slice before clear cosmic cut

            // Getting the PCA recalculated angle (using all hits in highest energy PFP with 10 cm of reco vertex)
            double pfp10cm_PCAAngle_beforeCuts = -999999;
            double pfp10cm_PCADX_beforeCuts = -999999;
            double pfp10cm_PCADY_beforeCuts = -999999;
            double pfp10cm_PCADZ_beforeCuts = -999999;
            
            for(size_t pfpAngle = 0; pfpAngle < angleRecalculationPCAPFP10cm_pfpID->size(); ++pfpAngle){
                if(angleRecalculationPCAPFP10cm_pfpID->at(pfpAngle) == highestEnergyPFP_beforeCuts.PFPID){
                    pfp10cm_PCAAngle_beforeCuts = angleRecalculationPCAPFP10cm_angle->at(pfpAngle);
                    pfp10cm_PCADX_beforeCuts = angleRecalculationPCAPFP10cm_dx->at(pfpAngle);
                    pfp10cm_PCADY_beforeCuts = angleRecalculationPCAPFP10cm_dy->at(pfpAngle);
                    pfp10cm_PCADZ_beforeCuts = angleRecalculationPCAPFP10cm_dz->at(pfpAngle);
                    if(pfp10cm_PCAAngle_beforeCuts < smallestAngle_beforeCuts) smallestAngle_beforeCuts = pfp10cm_PCAAngle_beforeCuts;
                }
            }

            // Calculating the angle difference between reconstructed highest energy PFP in slice and true recoil electron
            double angleDifference_beforeCuts = -999999;
            double angleDifferencePCAPFP10cm_beforeCuts = -999999;


            if((highestEnergyPFP_beforeCuts.dx != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((highestEnergyPFP_beforeCuts.dx * recoilElectron.dx) + (highestEnergyPFP_beforeCuts.dy * recoilElectron.dy) + (highestEnergyPFP_beforeCuts.dz * recoilElectron.dz));
                double aMagnitude = std::sqrt((highestEnergyPFP_beforeCuts.dx * highestEnergyPFP_beforeCuts.dx) + (highestEnergyPFP_beforeCuts.dy * highestEnergyPFP_beforeCuts.dy) + (highestEnergyPFP_beforeCuts.dz * highestEnergyPFP_beforeCuts.dz));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifference_beforeCuts = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            if((pfp10cm_PCADX_beforeCuts != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((pfp10cm_PCADX_beforeCuts * recoilElectron.dx) + (pfp10cm_PCADY_beforeCuts * recoilElectron.dy) + (pfp10cm_PCADZ_beforeCuts * recoilElectron.dz));
                double aMagnitude = std::sqrt((pfp10cm_PCADX_beforeCuts * pfp10cm_PCADX_beforeCuts) + (pfp10cm_PCADY_beforeCuts * pfp10cm_PCADY_beforeCuts) + (pfp10cm_PCADZ_beforeCuts * pfp10cm_PCADZ_beforeCuts));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifferencePCAPFP10cm_beforeCuts = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            // Counting number of events before cuts
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 1 && signal == 1){
                    eventsBeforeCuts_DLNuE.signal += weight;
                } else if(sliceCategoryPlottingMacro == 2 && signal == 1){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 3){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 4){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 5){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 6){
                    eventsBeforeCuts_DLNuE.background += weight;
                }

                if(sliceInteractionType == 0 && signal != 4){
                    eventsBeforeCuts_DLNuE.splitInt.cosmic += weight;   
                } else if(sliceInteractionType == 1 && signal == 1){
                    eventsBeforeCuts_DLNuE.splitInt.nuE += weight;
                } else if(sliceInteractionType == 2){
                    eventsBeforeCuts_DLNuE.splitInt.NCNPi0 += weight;
                } else if(sliceInteractionType == 3){
                    eventsBeforeCuts_DLNuE.splitInt.otherNC += weight;
                } else if(sliceInteractionType == 4){
                    eventsBeforeCuts_DLNuE.splitInt.CCnumu += weight;
                } else if(sliceInteractionType == 5){
                    eventsBeforeCuts_DLNuE.splitInt.CCnue += weight;
                } else if(sliceInteractionType == 6){
                    eventsBeforeCuts_DLNuE.splitInt.dirt += weight;
                } else if(sliceInteractionType == 7 && signal == 1){
                    eventsBeforeCuts_DLNuE.splitInt.nuEDirt += weight;
                } else if(sliceInteractionType == 8){
                    eventsBeforeCuts_DLNuE.splitInt.other += weight;
                } else if(sliceInteractionType == 9 && signal == 1){
                    eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy += weight;
                }
            }

            // Add to before cuts plots
            fillHistogram(&sliceCompletenessBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceCompleteness->at(slice), &weights);
            fillHistogram(&slicePurityBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_slicePurity->at(slice), &weights);
            fillHistogram(&sliceCRUMBSBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
            fillHistogram(&sliceNumRecoNeutBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
            fillHistogram(&sliceNumPFPsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_beforeCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_beforeCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsMinHitBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsMinHitSlice_beforeCuts, &weights);
            fillHistogram(&ERecoSumThetaRecoBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_beforeCuts.energy * pfp10cm_PCAAngle_beforeCuts * pfp10cm_PCAAngle_beforeCuts), &weights);
            fillHistogram(&dEdxBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
            fillHistogram(&razzledPDG11BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
            fillHistogram(&razzledPDG13BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
            fillHistogram(&razzledPDG22BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
            fillHistogram(&razzledPDG211BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
            fillHistogram(&razzledPDG2212BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
            fillHistogram(&pfpCompletenessBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.completeness, &weights);
            fillHistogram(&pfpPurityBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.purity, &weights); 
            fillHistogram(&recoVXBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
            fillHistogram(&recoVXSmallerBinsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYSmallerBinsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZSmallerBinsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
            fillHistogram(&recoVXLowBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYLowBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZLowBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
            fillHistogram(&recoVXHighBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYHighBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZHighBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
           
            fillSplitIntHistogram(&sliceCompletenessBeforeCuts, DLCurrent, signal, sliceInteractionType, reco_sliceCompleteness->at(slice), &weights);
            fillSplitIntHistogram(&slicePurityBeforeCuts, DLCurrent, signal, sliceInteractionType, reco_slicePurity->at(slice), &weights);
            fillSplitIntHistogram(&sliceCRUMBSBeforeCuts, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
            fillSplitIntHistogram(&sliceNumRecoNeutBeforeCuts, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);
            fillSplitIntHistogram(&sliceNumPFPsBeforeCuts, DLCurrent, signal, sliceInteractionType, numPFPsSlice_beforeCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsBeforeCuts, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_beforeCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsMinHitBeforeCuts, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsMinHitSlice_beforeCuts, &weights);
            fillSplitIntHistogram(&ERecoSumThetaRecoBeforeCuts, DLCurrent, signal, sliceInteractionType, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoBeforeCuts, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_beforeCuts.energy * pfp10cm_PCAAngle_beforeCuts * pfp10cm_PCAAngle_beforeCuts), &weights);
            fillSplitIntHistogram(&dEdxBeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
            fillSplitIntHistogram(&razzledPDG11BeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
            fillSplitIntHistogram(&razzledPDG13BeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
            fillSplitIntHistogram(&razzledPDG22BeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
            fillSplitIntHistogram(&razzledPDG211BeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
            fillSplitIntHistogram(&razzledPDG2212BeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
            fillSplitIntHistogram(&pfpCompletenessBeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.completeness, &weights);
            fillSplitIntHistogram(&pfpPurityBeforeCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.purity, &weights); 
            fillSplitIntHistogram(&recoVXBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
            fillSplitIntHistogram(&recoVXSmallerBinsBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYSmallerBinsBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZSmallerBinsBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
            fillSplitIntHistogram(&recoVXLowBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYLowBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZLowBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
            fillSplitIntHistogram(&recoVXHighBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYHighBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZHighBeforeCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);

            if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999){
                fillHistogram(&angleDifferenceSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifference_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCAPFPSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP5cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP5cm_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP10cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP10cm_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP15cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP15cm_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCASliceSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCASlice5cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice5cm_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCASlice10cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice10cm_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCASlice15cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice15cm_beforeCuts, &weights);
                fillHistogram(&energyAsymmetryBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy)/recoilElectron.energy), &weights);

                fillSplitIntHistogram(&angleDifferenceSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifference_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFPSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFP5cmSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP5cm_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFP10cmSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP10cm_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFP15cmSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP15cm_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASliceSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASlice5cmSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice5cm_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASlice10cmSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice10cm_beforeCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASlice15cmSignalBeforeCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice15cm_beforeCuts, &weights);
                fillSplitIntHistogram(&energyAsymmetryBeforeCuts, DLCurrent, signal, sliceInteractionType, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy)/recoilElectron.energy), &weights);                
            }           
 
            // Start applying cuts here, this macro has no option to turn off the removal of clear cosmic PFPs
            double summedEnergy_afterCuts = 0;
            double numPFPsSlice_afterCuts = 0;
            double numPrimaryPFPsSlice_afterCuts = 0;
            double numPrimaryPFPsMinHitSlice_afterCuts = 0;

            highestEnergyPFP_struct highestEnergyPFP_afterCuts;
            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                    if(reco_particleClearCosmic->at(pfp) == 0){
                        numPFPsSlice_afterCuts++;
                        if(reco_particleIsPrimary->at(pfp) == 1){
                            numPrimaryPFPsSlice_afterCuts++; // PFP is a primary PFP
                            if(reco_particleNumHits->at(pfp) >= primaryPFPMinHitRequirement) numPrimaryPFPsMinHitSlice_afterCuts++;
                        }

                        summedEnergy_afterCuts += reco_particleBestPlaneEnergy->at(pfp);
                        if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergyPFP_afterCuts.energy){
                            highestEnergyPFP_afterCuts.energy = reco_particleBestPlaneEnergy->at(pfp);
                            highestEnergyPFP_afterCuts.theta = reco_particleTheta->at(pfp);
                            highestEnergyPFP_afterCuts.PFPID = reco_particleID->at(pfp);
                            highestEnergyPFP_afterCuts.dx = reco_particleDX->at(pfp);
                            highestEnergyPFP_afterCuts.dy = reco_particleDY->at(pfp);
                            highestEnergyPFP_afterCuts.dz = reco_particleDZ->at(pfp);
                            highestEnergyPFP_afterCuts.vx = reco_particleVX->at(pfp);
                            highestEnergyPFP_afterCuts.vy = reco_particleVY->at(pfp);
                            highestEnergyPFP_afterCuts.vz = reco_particleVZ->at(pfp);
                            highestEnergyPFP_afterCuts.completeness = reco_particleCompleteness->at(pfp);
                            highestEnergyPFP_afterCuts.purity = reco_particlePurity->at(pfp);
                            highestEnergyPFP_afterCuts.trackscore = reco_particleTrackScore->at(pfp);
                            highestEnergyPFP_afterCuts.primary = reco_particleIsPrimary->at(pfp);
                            highestEnergyPFP_afterCuts.truePDG = reco_particleTruePDG->at(pfp);
                            highestEnergyPFP_afterCuts.trueOrigin = reco_particleTrueOrigin->at(pfp);
                            highestEnergyPFP_afterCuts.trueInt = reco_particleTrueInteractionType->at(pfp);
                            highestEnergyPFP_afterCuts.bestPlanedEdx = reco_particleBestPlanedEdx->at(pfp);
                            highestEnergyPFP_afterCuts.razzledPDG11 = reco_particleRazzledPDG11->at(pfp);
                            highestEnergyPFP_afterCuts.razzledPDG13 = reco_particleRazzledPDG13->at(pfp);
                            highestEnergyPFP_afterCuts.razzledPDG22 = reco_particleRazzledPDG22->at(pfp);
                            highestEnergyPFP_afterCuts.razzledPDG211 = reco_particleRazzledPDG211->at(pfp);
                            highestEnergyPFP_afterCuts.razzledPDG2212 = reco_particleRazzledPDG2212->at(pfp);
                            highestEnergyPFP_afterCuts.razzledBestPDG = reco_particleRazzledBestPDG->at(pfp);
                            highestEnergyPFP_afterCuts.trueVX = reco_particleTrueVX->at(pfp);
                            highestEnergyPFP_afterCuts.trueVY = reco_particleTrueVY->at(pfp);
                            highestEnergyPFP_afterCuts.trueVZ = reco_particleTrueVZ->at(pfp);
                            highestEnergyPFP_afterCuts.trueEndX = reco_particleTrueEndX->at(pfp);
                            highestEnergyPFP_afterCuts.trueEndY = reco_particleTrueEndY->at(pfp);
                            highestEnergyPFP_afterCuts.trueEndZ = reco_particleTrueEndZ->at(pfp);
                            highestEnergyPFP_afterCuts.numHits = reco_particleNumHits->at(pfp);
                            highestEnergyPFP_afterCuts.clearCosmic = reco_particleClearCosmic->at(pfp);

                            if(highestEnergyPFP_afterCuts.trueVX != -999999 && highestEnergyPFP_afterCuts.trueVY != -999999 && highestEnergyPFP_afterCuts.trueVZ != -999999 && highestEnergyPFP_afterCuts.trueEndX != -999999 && highestEnergyPFP_afterCuts.trueEndY != -999999 && highestEnergyPFP_afterCuts.trueEndZ != -999999){
                                double xCoordDiff_length = (highestEnergyPFP_afterCuts.trueVX - highestEnergyPFP_afterCuts.trueEndX);
                                double yCoordDiff_length = (highestEnergyPFP_afterCuts.trueVY - highestEnergyPFP_afterCuts.trueEndY);
                                double zCoordDiff_length = (highestEnergyPFP_afterCuts.trueVZ - highestEnergyPFP_afterCuts.trueEndZ);
                                highestEnergyPFP_afterCuts.trueLength = std::sqrt((xCoordDiff_length * xCoordDiff_length) + (yCoordDiff_length * yCoordDiff_length) + (zCoordDiff_length * zCoordDiff_length));
                            }
                        }
                    }
                }

            } // End of looping through PFPs
                
            double pfp10cm_PCAAngle_afterCuts = -999999;
            double pfp10cm_PCADX_afterCuts = -999999;
            double pfp10cm_PCADY_afterCuts = -999999;
            double pfp10cm_PCADZ_afterCuts = -999999;

            for(size_t pfpAngle = 0; pfpAngle < angleRecalculationPCAPFP10cm_pfpID->size(); ++pfpAngle){
                if(angleRecalculationPCAPFP10cm_pfpID->at(pfpAngle) == highestEnergyPFP_afterCuts.PFPID){
                    pfp10cm_PCAAngle_afterCuts = angleRecalculationPCAPFP10cm_angle->at(pfpAngle);
                    pfp10cm_PCADX_afterCuts = angleRecalculationPCAPFP10cm_dx->at(pfpAngle);
                    pfp10cm_PCADY_afterCuts = angleRecalculationPCAPFP10cm_dy->at(pfpAngle);
                    pfp10cm_PCADZ_afterCuts = angleRecalculationPCAPFP10cm_dz->at(pfpAngle);
                    if(pfp10cm_PCAAngle_afterCuts < smallestAngle_afterCuts) smallestAngle_afterCuts = pfp10cm_PCAAngle_afterCuts;
                }
            }

            double angleDifference_afterCuts = -999999;
            double angleDifferencePCAPFP10cm_afterCuts = -999999;

            if((highestEnergyPFP_afterCuts.dx != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((highestEnergyPFP_afterCuts.dx * recoilElectron.dx) + (highestEnergyPFP_afterCuts.dy * recoilElectron.dy) + (highestEnergyPFP_afterCuts.dz * recoilElectron.dz));
                double aMagnitude = std::sqrt((highestEnergyPFP_afterCuts.dx * highestEnergyPFP_afterCuts.dx) + (highestEnergyPFP_afterCuts.dy * highestEnergyPFP_afterCuts.dy) + (highestEnergyPFP_afterCuts.dz * highestEnergyPFP_afterCuts.dz));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifference_afterCuts = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            if((pfp10cm_PCADX_afterCuts != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((pfp10cm_PCADX_afterCuts * recoilElectron.dx) + (pfp10cm_PCADY_afterCuts * recoilElectron.dy) + (pfp10cm_PCADZ_afterCuts * recoilElectron.dz));
                double aMagnitude = std::sqrt((pfp10cm_PCADX_afterCuts * pfp10cm_PCADX_afterCuts) + (pfp10cm_PCADY_afterCuts * pfp10cm_PCADY_afterCuts) + (pfp10cm_PCADZ_afterCuts * pfp10cm_PCADZ_afterCuts));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifferencePCAPFP10cm_afterCuts = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            // Clear cosmic cut has been applied, add to counters
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEFuzzy += weight;
            }

            if(numPFPs0Cut == 1 && numPFPsSlice_afterCuts == 0){
                // This is a slice with 0 PFPs in it
                continue;
            }
            
            // Number of PFPs 0 cut has been applied, add to counters
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.numPFPs0Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numPFPs0IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEFuzzy += weight;
            }

            if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0){
                // This is a slice with no reco neutrino
                continue;
            }

            // Number of reco neutrinos cut has been applied, add to counters
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEFuzzy += weight;
            }

            fillHistogram(&sliceCRUMBSAfterRecoNeutrinoCut, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);

            if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high)){
                // This is a slice with a CRUMBS score outside cut values
                continue;
            }
            
            // CRUMBS score cut has been applied, add to counters
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.crumbsSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.crumbsBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.crumbsBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.crumbsIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.crumbsIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuEFuzzy += weight;
            }
            
            fillHistogram(&sliceCRUMBSAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);

            if(FVCut == 1){
                if(!(recoVX < FVCut_xHigh && recoVX > FVCut_xLow  && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow && recoVZ < FVCut_zHigh)){
                    // Doesn't pass the FV cut values
                    continue;
                }
            }

            // FV cut applied, fill counters
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.FVSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.FVBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.FVBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.FVIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.FVIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.FVIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.FVIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.FVIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.FVIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuEFuzzy += weight;
            }

            if(primaryPFPCut == 1 && numPrimaryPFPsMinHitSlice_afterCuts != primaryPFPCutValue){
                // Slice has more than 1 primary PFP in it
                continue;
            }

            // Primary PFP cut has been applied, fill counters
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.primaryPFPBack += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.primaryPFPIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEFuzzy += weight;
            }

            if(razzledPDG11Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG11 > razzled11High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG11 < razzled11Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 11 cut
                continue;
            }

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled11Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled11Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled11Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled11IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled11IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuEFuzzy += weight;
            }

            if(razzledPDG22Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG22 > razzled22High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG22 < razzled22Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 22 cut
                continue;
            }

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled22Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled22Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled22Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled22Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled22Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled22Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled22Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled22IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled22IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled22IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled22IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled22IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled22IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled22IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled22IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled22IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled22IntSplit.nuEFuzzy += weight;
            }

            if(razzledPDG13Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG13 > razzled13High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG13 < razzled13Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 13 cut
                continue;
            }

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled13Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled13Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled13Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled13Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled13Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled13Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled13Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled13IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled13IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled13IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled13IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled13IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled13IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled13IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled13IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled13IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled13IntSplit.nuEFuzzy += weight;
            }

            if(razzledPDG2212Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG2212 > razzled2212High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG2212 < razzled2212Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 2212 cut
                continue;
            }
            
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled2212Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled2212Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled2212IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled2212IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled2212IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled2212IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled2212IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled2212IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled2212IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled2212IntSplit.nuEFuzzy += weight;
            }

            if(razzledPDG211Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG211 > razzled211High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG211 < razzled211Low_highestEnergyPFP))){
                // Highest energy PFP in slice doesn't pass the razzled 211 cut
                continue;
            }

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled211Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.razzled211Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.razzled211Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled211IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled211IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuEFuzzy += weight;
            }
            
            if(ETheta2Cut == 1 && ((highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts) > ETheta2High_highestEnergyPFP || (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts) < ETheta2Low_highestEnergyPFP)){
                // Highest energy PFP in slice doesn't pass the ETheta2 cut
                continue;
            }
                
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.ETheta2Sig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.ETheta2Back += weight;

                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.ETheta2IntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEFuzzy += weight;
            }

            if(dEdxCut == 1 && (highestEnergyPFP_afterCuts.bestPlanedEdx > dEdxHigh_highestEnergyPFP || highestEnergyPFP_afterCuts.bestPlanedEdx < dEdxLow_highestEnergyPFP)){
                // Highest energy PFP in slice doesn't pass the dE/dx cut
                continue;
            }

            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0 && signal != 4) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.dEdxSig += weight;
                else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 5) eventsAfterCuts_DLNuE.dEdxBack += weight;
                else if(sliceCategoryPlottingMacro == 6) eventsAfterCuts_DLNuE.dEdxBack += weight;


                if(sliceInteractionType == 0 && signal != 4) eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic += weight;
                else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuE += weight;
                else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0 += weight;
                else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC += weight;
                else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu += weight;
                else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue += weight;
                else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.dEdxIntSplit.dirt += weight;
                else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt += weight;
                else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.dEdxIntSplit.other += weight;
                else if(sliceInteractionType == 9 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuEFuzzy += weight;
            }

            // Fill histograms after cuts here
            fillHistogram(&sliceCompletenessAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceCompleteness->at(slice), &weights);
            fillHistogram(&slicePurityAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_slicePurity->at(slice), &weights);
            fillHistogram(&sliceCRUMBSAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
            fillHistogram(&sliceNumRecoNeutAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
            fillHistogram(&sliceNumPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_afterCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_afterCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsMinHitAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsMinHitSlice_afterCuts, &weights);
            fillHistogram(&ERecoSumThetaRecoAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (summedEnergy_afterCuts * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoAfterCuts_pfp10cmPoints, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);
            fillHistogram(&dEdxAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
            fillHistogram(&razzledPDG11AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
            fillHistogram(&razzledPDG13AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG13, &weights);
            fillHistogram(&razzledPDG22AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG22, &weights);
            fillHistogram(&razzledPDG211AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
            fillHistogram(&razzledPDG2212AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG2212, &weights);
            fillHistogram(&pfpCompletenessAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.completeness, &weights);
            fillHistogram(&pfpPurityAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.purity, &weights); 
            fillHistogram(&recoVXAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
            fillHistogram(&recoVXSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
            fillHistogram(&recoVXLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
            fillHistogram(&recoVXHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
            fillHistogram(&recoVYHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
            fillHistogram(&recoVZHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
           
            fillSplitIntHistogram(&sliceCompletenessAfterCuts, DLCurrent, signal, sliceInteractionType, reco_sliceCompleteness->at(slice), &weights);
            fillSplitIntHistogram(&slicePurityAfterCuts, DLCurrent, signal, sliceInteractionType, reco_slicePurity->at(slice), &weights);
            fillSplitIntHistogram(&sliceCRUMBSAfterCuts, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
            fillSplitIntHistogram(&sliceNumRecoNeutAfterCuts, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);
            fillSplitIntHistogram(&sliceNumPFPsAfterCuts, DLCurrent, signal, sliceInteractionType, numPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsAfterCuts, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsMinHitAfterCuts, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsMinHitSlice_afterCuts, &weights);
            fillSplitIntHistogram(&ERecoSumThetaRecoAfterCuts, DLCurrent, signal, sliceInteractionType, (summedEnergy_afterCuts * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterCuts, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterCuts_pfp10cmPoints, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);
            fillSplitIntHistogram(&dEdxAfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
            fillSplitIntHistogram(&razzledPDG11AfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
            fillSplitIntHistogram(&razzledPDG13AfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG13, &weights);
            fillSplitIntHistogram(&razzledPDG22AfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG22, &weights);
            fillSplitIntHistogram(&razzledPDG211AfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
            fillSplitIntHistogram(&razzledPDG2212AfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG2212, &weights);
            fillSplitIntHistogram(&pfpCompletenessAfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.completeness, &weights);
            fillSplitIntHistogram(&pfpPurityAfterCuts, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.purity, &weights); 
            fillSplitIntHistogram(&recoVXAfterCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYAfterCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZAfterCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
            fillSplitIntHistogram(&recoVXSmallerBinsAfterCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYSmallerBinsAfterCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZSmallerBinsAfterCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
            fillSplitIntHistogram(&recoVXLowAfterCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYLowAfterCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZLowAfterCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
            fillSplitIntHistogram(&recoVXHighAfterCuts, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
            fillSplitIntHistogram(&recoVYHighAfterCuts, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
            fillSplitIntHistogram(&recoVZHighAfterCuts, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);

            if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999){
                fillHistogram(&angleDifferenceSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifference_afterCuts, &weights);
                fillHistogram(&angleDifferencePCAPFPSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP_afterCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP5cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP5cm_afterCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP10cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP10cm_afterCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP15cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP15cm_afterCuts, &weights);
                fillHistogram(&angleDifferencePCASliceSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice_afterCuts, &weights);
                fillHistogram(&angleDifferencePCASlice5cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice5cm_afterCuts, &weights);
                fillHistogram(&angleDifferencePCASlice10cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice10cm_afterCuts, &weights);
                fillHistogram(&angleDifferencePCASlice15cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCASlice15cm_afterCuts, &weights);
                fillHistogram(&energyAsymmetryAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy)/recoilElectron.energy), &weights);

                fillSplitIntHistogram(&angleDifferenceSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifference_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFPSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFP5cmSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP5cm_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFP10cmSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP10cm_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCAPFP15cmSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCAPFP15cm_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASliceSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASlice5cmSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice5cm_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASlice10cmSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice10cm_afterCuts, &weights);
                fillSplitIntHistogram(&angleDifferencePCASlice15cmSignalAfterCuts, DLCurrent, signal, sliceInteractionType, angleDifferencePCASlice15cm_afterCuts, &weights);
                fillSplitIntHistogram(&energyAsymmetryAfterCuts, DLCurrent, signal, sliceInteractionType, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy)/recoilElectron.energy), &weights);                
            }           

        } // End of looping through slices

    } // End of looping through events

    // Draw histograms here
    styleDrawAll(sliceCRUMBSBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCRUMBSAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    // Make table of efficiency, purity, etc for each cut
    std::ofstream out_tablefile(tableFileName, std::ios::app);
    if(out_tablefile.is_open()){
        out_tablefile << "=========== DL Nu+E Vertexing ===========" << std::endl;
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\epsilon$ (\\%)} & \\textbf{Selection $\\epsilon$ (\\%)} & \\textbf{$\\rho$ (\\%)} & \\textbf{$\\epsilon\\rho$} & \\textbf{Selection $\\epsilon\\rho$} & Signal Left & Background Left \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << std::defaultfloat << std::setprecision(7) << "No Cut & " << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.signal/actualSignalCount << " & " << 100*eventsBeforeCuts_DLNuE.signal/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsBeforeCuts_DLNuE.signal/(eventsBeforeCuts_DLNuE.signal+eventsBeforeCuts_DLNuE.background) << " & " << (eventsBeforeCuts_DLNuE.signal/(eventsBeforeCuts_DLNuE.signal+eventsBeforeCuts_DLNuE.background))*(eventsBeforeCuts_DLNuE.signal/actualSignalCount) << " & " << (eventsBeforeCuts_DLNuE.signal/(eventsBeforeCuts_DLNuE.signal+eventsBeforeCuts_DLNuE.background))*(eventsBeforeCuts_DLNuE.signal/eventsBeforeCuts_DLNuE.signal) << " & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.signal << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.signal/actualSignalCount << "\\%)" << std::fixed << std::setprecision(0) << " & " << eventsBeforeCuts_DLNuE.background << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.background/eventsBeforeCuts_DLNuE.background << "\\%)" << " \\\\ " << std::endl;
        out_tablefile << "\\hline" << std::endl;
       
        if(clearCosmicCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Remove Clear Cosmic PFPs & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/(eventsAfterCuts_DLNuE.clearCosmicsSig+eventsAfterCuts_DLNuE.clearCosmicsBack) << " & " << (eventsAfterCuts_DLNuE.clearCosmicsSig/actualSignalCount)*(eventsAfterCuts_DLNuE.clearCosmicsSig/(eventsAfterCuts_DLNuE.clearCosmicsSig+eventsAfterCuts_DLNuE.clearCosmicsBack)) << " & " << (eventsAfterCuts_DLNuE.clearCosmicsSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.clearCosmicsSig/(eventsAfterCuts_DLNuE.clearCosmicsSig+eventsAfterCuts_DLNuE.clearCosmicsBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numPFPs0Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "PFPs in Slice != 0 & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/(eventsAfterCuts_DLNuE.numPFPs0Sig+eventsAfterCuts_DLNuE.numPFPs0Back) << " & " << (eventsAfterCuts_DLNuE.numPFPs0Sig/actualSignalCount)*(eventsAfterCuts_DLNuE.numPFPs0Sig/(eventsAfterCuts_DLNuE.numPFPs0Sig+eventsAfterCuts_DLNuE.numPFPs0Back)) << " & " << (eventsAfterCuts_DLNuE.numPFPs0Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.numPFPs0Sig/(eventsAfterCuts_DLNuE.numPFPs0Sig+eventsAfterCuts_DLNuE.numPFPs0Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numPFPs0Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/actualSignalCount << std::fixed << std::setprecision(0) << "\\%) & " << eventsAfterCuts_DLNuE.numPFPs0Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numRecoNeutrinosCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "1 Reco Neutrino in Slice & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/(eventsAfterCuts_DLNuE.numRecoNeut0Sig+eventsAfterCuts_DLNuE.numRecoNeut0Back) << " & " << (eventsAfterCuts_DLNuE.numRecoNeut0Sig/actualSignalCount)*(eventsAfterCuts_DLNuE.numRecoNeut0Sig/(eventsAfterCuts_DLNuE.numRecoNeut0Sig+eventsAfterCuts_DLNuE.numRecoNeut0Back)) << " & " << (eventsAfterCuts_DLNuE.numRecoNeut0Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.numRecoNeut0Sig/(eventsAfterCuts_DLNuE.numRecoNeut0Sig+eventsAfterCuts_DLNuE.numRecoNeut0Back)) << std::fixed << std::setprecision(0) << " & " << eventsAfterCuts_DLNuE.numRecoNeut0Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numRecoNeut0Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(CRUMBSCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high << " & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.crumbsSig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.crumbsSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.crumbsSig/(eventsAfterCuts_DLNuE.crumbsSig+eventsAfterCuts_DLNuE.crumbsBack) << " & " << (eventsAfterCuts_DLNuE.crumbsSig/actualSignalCount)*(eventsAfterCuts_DLNuE.crumbsSig/(eventsAfterCuts_DLNuE.crumbsSig+eventsAfterCuts_DLNuE.crumbsBack)) << " & " << (eventsAfterCuts_DLNuE.crumbsSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.crumbsSig/(eventsAfterCuts_DLNuE.crumbsSig+eventsAfterCuts_DLNuE.crumbsBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsSig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
       
        if(FVCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "FV Cut & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.FVSig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.FVSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.FVSig/(eventsAfterCuts_DLNuE.FVSig+eventsAfterCuts_DLNuE.FVBack) << " & " << (eventsAfterCuts_DLNuE.FVSig/actualSignalCount)*(eventsAfterCuts_DLNuE.FVSig/(eventsAfterCuts_DLNuE.FVSig+eventsAfterCuts_DLNuE.FVBack)) << " & " << (eventsAfterCuts_DLNuE.FVSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.FVSig/(eventsAfterCuts_DLNuE.FVSig+eventsAfterCuts_DLNuE.FVBack)) << std::fixed << std::setprecision(0) << " & " << eventsAfterCuts_DLNuE.FVSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVSig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.FVBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(primaryPFPCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Primary PFPs in Slice with $\\geq$ " << primaryPFPMinHitRequirement << " Hits = " << primaryPFPCutValue << " & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.primaryPFPSig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.primaryPFPSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.primaryPFPSig/(eventsAfterCuts_DLNuE.primaryPFPSig+eventsAfterCuts_DLNuE.primaryPFPBack) << " & " << (eventsAfterCuts_DLNuE.primaryPFPSig/actualSignalCount)*(eventsAfterCuts_DLNuE.primaryPFPSig/(eventsAfterCuts_DLNuE.primaryPFPSig+eventsAfterCuts_DLNuE.primaryPFPBack)) << " & " << (eventsAfterCuts_DLNuE.primaryPFPSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.primaryPFPSig/(eventsAfterCuts_DLNuE.primaryPFPSig+eventsAfterCuts_DLNuE.primaryPFPBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPSig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG11Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Electron Score $\\geq$ " << razzled11Low_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled11Sig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.razzled11Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled11Sig/(eventsAfterCuts_DLNuE.razzled11Sig+eventsAfterCuts_DLNuE.razzled11Back) << " & " << (eventsAfterCuts_DLNuE.razzled11Sig/actualSignalCount)*(eventsAfterCuts_DLNuE.razzled11Sig/(eventsAfterCuts_DLNuE.razzled11Sig+eventsAfterCuts_DLNuE.razzled11Back)) << " & " << (eventsAfterCuts_DLNuE.razzled11Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled11Sig/(eventsAfterCuts_DLNuE.razzled11Sig+eventsAfterCuts_DLNuE.razzled11Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11Sig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG211Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ " << razzled211High_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled211Sig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.razzled211Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled211Sig/(eventsAfterCuts_DLNuE.razzled211Sig+eventsAfterCuts_DLNuE.razzled211Back) << " & " << (eventsAfterCuts_DLNuE.razzled211Sig/actualSignalCount)*(eventsAfterCuts_DLNuE.razzled211Sig/(eventsAfterCuts_DLNuE.razzled211Sig+eventsAfterCuts_DLNuE.razzled211Back)) << " & " << (eventsAfterCuts_DLNuE.razzled211Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled211Sig/(eventsAfterCuts_DLNuE.razzled211Sig+eventsAfterCuts_DLNuE.razzled211Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211Sig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(ETheta2Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP + PFP Spacepoints 10cm)} $\\leq$ " << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$ & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.ETheta2Sig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.ETheta2Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.ETheta2Sig/(eventsAfterCuts_DLNuE.ETheta2Sig+eventsAfterCuts_DLNuE.ETheta2Back) << " & " << (eventsAfterCuts_DLNuE.ETheta2Sig/actualSignalCount)*(eventsAfterCuts_DLNuE.ETheta2Sig/(eventsAfterCuts_DLNuE.ETheta2Sig+eventsAfterCuts_DLNuE.ETheta2Back)) << " & " << (eventsAfterCuts_DLNuE.ETheta2Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.ETheta2Sig/(eventsAfterCuts_DLNuE.ETheta2Sig+eventsAfterCuts_DLNuE.ETheta2Back)) << std::fixed << std::setprecision(0) << " & " << eventsAfterCuts_DLNuE.ETheta2Sig << " ("  << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.ETheta2Sig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.ETheta2Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(dEdxCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has " << dEdxLow_highestEnergyPFP << " MeV cm^{-1} $\\leq$ dE/dx $\\leq$ " << dEdxHigh_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " MeV cm^{-1} & " << 100*eventsAfterCuts_DLNuE.dEdxSig/actualSignalCount << " & " << 100*eventsAfterCuts_DLNuE.dEdxSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.dEdxSig/(eventsAfterCuts_DLNuE.dEdxSig+eventsAfterCuts_DLNuE.dEdxBack) << " & " << (eventsAfterCuts_DLNuE.dEdxSig/actualSignalCount)*(eventsAfterCuts_DLNuE.dEdxSig/(eventsAfterCuts_DLNuE.dEdxSig+eventsAfterCuts_DLNuE.dEdxBack)) << " & " << (eventsAfterCuts_DLNuE.dEdxSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.dEdxSig/(eventsAfterCuts_DLNuE.dEdxSig+eventsAfterCuts_DLNuE.dEdxBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxSig/actualSignalCount << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;
       
        out_tablefile << "" << std::endl;
        out_tablefile << "" << std::endl;
        out_tablefile << "" << std::endl;
        
        // Put split interaction table here
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c|c| }" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\multicolumn{11}{|c|}{\\textbf{Number of Events Left}} \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\boldsymbol{\\nu+e}$} & \\textbf{NCN$\\boldsymbol{\\pi^0}$} & \\textbf{Other NC} & \\textbf{CC$\\boldsymbol{\\nu_\\mu}$} & \\textbf{CC$\\boldsymbol{\\nu_e}$} & \\textbf{Dirt} & \\textbf{$\\boldsymbol{\\nu+e}$ Dirt} & \\textbf{Cosmic} & \\textbf{Other} & \\textbf{$\\boldsymbol{\\nu+e}$ Fuzzy}\\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "No Cut & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.splitInt.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << std::defaultfloat << std::setprecision(4) << "(" << 100*eventsBeforeCuts_DLNuE.splitInt.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.otherNC << " (" << 100*eventsBeforeCuts_DLNuE.splitInt.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.CCnumu << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.CCnue << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.dirt << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuEDirt << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.cosmic << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.other << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        if(clearCosmicCut == 1){
            out_tablefile << "Remove Clear Cosmic PFPs & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numPFPs0Cut == 1){
            out_tablefile << "PFPs in Slice != 0 & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numRecoNeutrinosCut == 1){
            out_tablefile << "1 Reco Neutrino in Slice & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
       
        if(CRUMBSCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
       
        if(FVCut == 1){
            out_tablefile << "FV Cut & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.FVIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(primaryPFPCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Primary PFPs in Slice with $\\geq$ " << primaryPFPMinHitRequirement << " = " << primaryPFPCutValue << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG11Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Electron Score $\\geq$ " << razzled11Low_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG211Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ " << razzled211High_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(ETheta2Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP + PFP Spacepoints 10cm)} $\\leq$ " << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$ & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(dEdxCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has " << dEdxLow_highestEnergyPFP << " MeV cm^{-1} $\\leq$ dE/dx $\\leq$ " << dEdxHigh_highestEnergyPFP << " MeV cm^{-1} & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.nuEFuzzy << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuEFuzzy/eventsBeforeCuts_DLNuE.splitInt.nuEFuzzy << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;

        out_tablefile << "" << std::endl;
        out_tablefile << "\\newpage" << std::endl;
        out_tablefile << "" << std::endl; 
    }

}
