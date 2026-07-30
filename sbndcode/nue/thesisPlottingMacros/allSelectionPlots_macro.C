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

TFile* outRootFile = nullptr;

void saveCanvasToRootFile(TCanvas* c, const std::string& filename){
    if(!outRootFile || !c) return;

    std::string name = gSystem->BaseName(filename.c_str());
    size_t dotPos = name.find_last_of('.');
    if(dotPos != std::string::npos) name = name.substr(0, dotPos);

    outRootFile->cd();
    c->Write(name.c_str());
}

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
                  const char* filename, const char* rootname, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr,
                  bool includeSignal = true, bool includeSignalFuzzy = true,
                  bool includeBNB = true, bool includeBNBFuzzy = true,
                  bool includeCosmic = true,
                  bool includeDLUboone = true, bool includeDLNuE = true,
                  bool includeBDT = true,
                  bool useLogScale = false, bool bestPDGPlot = false, bool fvPlot = false)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    if (useLogScale)
        hists.canvas->SetLogy(1);
    else
        hists.canvas->SetLogy(0);

    std::vector<TH1D*> allHists = {hists.nuESignal, hists.nuESignalFuzzy, hists.nuEBNB, hists.nuEBNBFuzzy, hists.nuECosmic};

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

    if(fvPlot){
        for(auto* hist : allHists){
            if(!hist) continue;
            hist->GetXaxis()->SetBinLabel(1, "No");
            hist->GetXaxis()->SetBinLabel(2, "Yes");
            hist->LabelsOption("h");
        }
    }

    for (auto* hist : allHists)
        if (hist) hist->SetStats(0);

    hists.nuECosmic->SetLineWidth(2);      hists.nuECosmic->SetLineColor(kPink-2);
    hists.nuESignal->SetLineWidth(2);      hists.nuESignal->SetLineColor(kAzure+5);
    hists.nuESignalFuzzy->SetLineWidth(2);      hists.nuESignalFuzzy->SetLineColor(kGreen-7);
    hists.nuEBNB->SetLineWidth(2);      hists.nuEBNB->SetLineColor(kOrange-5);
    hists.nuEBNBFuzzy->SetLineWidth(2);      hists.nuEBNBFuzzy->SetLineColor(kViolet+4);

    if((ymin != 999) && (ymax != 999)){
        hists.nuESignal->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuESignalFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEBNB->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEBNBFuzzy->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuECosmic->GetYaxis()->SetRangeUser(ymin, ymax);
    }
    
    if((xmin != 999) && (xmax != 999)){
        hists.nuESignal->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuESignalFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEBNB->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEBNBFuzzy->GetXaxis()->SetRangeUser(xmin, xmax);
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
        if (variantAllowed("nuESignal")) draw(hists.nuESignal);
    }
    if (includeSignalFuzzy) {
        if (variantAllowed("nuESignalFuzzy")) draw(hists.nuESignalFuzzy);
    }
    if (includeBNB) {
        if (variantAllowed("nuEBNB")) draw(hists.nuEBNB);
    }
    if (includeBNBFuzzy) {
        if (variantAllowed("nuEBNBFuzzy")) draw(hists.nuEBNBFuzzy);
    }
    if (includeCosmic) {
        if (variantAllowed("nuECosmic")) draw(hists.nuECosmic);
    }

    hists.nuESignal->SetStats(0);
    hists.nuESignal->GetXaxis()->SetTickLength(0.04);
    hists.nuESignal->GetYaxis()->SetTickLength(0.03);
    hists.nuESignal->GetXaxis()->SetTickSize(0.02);
    hists.nuESignal->GetYaxis()->SetTickSize(0.02);

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    std::vector<std::pair<TH1*, std::string>> legendEntries;

    auto addLegendIf = [&](TH1* hist, const std::string& label, const std::string& name){
        if (hist && variantAllowed(name)) legendEntries.emplace_back(hist, label);
    };

    if (includeSignal) {
        addLegendIf(hists.nuESignal, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "nuESignal");
    }
    if (includeSignalFuzzy) {
        addLegendIf(hists.nuESignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "nuESignalFuzzy");
    }
    if (includeBNB) {
        addLegendIf(hists.nuEBNB, "BNB, Pandora Deep Learning: SBND Nu+E Tune", "nuEBNB");
    }
    if (includeBNBFuzzy) {
        addLegendIf(hists.nuEBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "nuEBNBFuzzy");
    }
    if (includeCosmic) {
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
    saveCanvasToRootFile(hists.canvas, rootname);
}

void styleDrawBackSig(histGroup_struct hists,
                      double ymin, double ymax, double xmin, double xmax,
                      const char* filename, const char* rootname, const std::string& legendLocation,
                      bool includeCurrent = true, bool includeUboone = true, bool includeNuE = true,
                      bool useLogScale = false, bool bestPDGPlot = false, bool fvPlot = false)
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

    TH1D* nuESignalCombo     = combine(hists.nuESignal, nullptr, nullptr, nullptr, "nuESignalCombo");
    TH1D* nuEBackgroundCombo = combine(hists.nuEBNB, hists.nuEBNBFuzzy, hists.nuECosmic, hists.nuESignalFuzzy, "nuEBackgroundCombo");

    std::vector<TH1D*> allHists = {nuESignalCombo, nuEBackgroundCombo};
    
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
    
    if(fvPlot){
        for(auto* hist : allHists){
            if(!hist) continue;
            hist->GetXaxis()->SetBinLabel(1, "No");
            hist->GetXaxis()->SetBinLabel(2, "Yes");
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

    if (nuESignalCombo)         { nuESignalCombo->SetLineWidth(2);     nuESignalCombo->SetLineColor(kAzure+5); }
    if (nuEBackgroundCombo)     { nuEBackgroundCombo->SetLineWidth(2);     nuEBackgroundCombo->SetLineColor(kOrange-5); }

    bool first = true;
    auto draw = [&](TH1* hist){ if (hist) { hist->Draw(first ? "hist" : "histsame"); first = false; } };

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

    if (includeNuE) {
        legend->AddEntry(nuESignalCombo, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(nuEBackgroundCombo, "Background, Pandora Deep Learning: SBND Nu+E Tune", "f");
    }

    legend->SetTextSize(0.015);
    legend->SetMargin(0.12);
    legend->Draw();

    hists.canvas->SaveAs(filename);
    saveCanvasToRootFile(hists.canvas, rootname); 

    for (auto* hist : allHists)
        delete hist;
}


void styleDrawSplit(splitHistGroup_struct hists,
                    double ymin, double ymax, double xmin, double xmax,
                    const char* filename, const char* rootname, const std::string& legendLocation,
                    int* drawLine = nullptr, int* linePos = nullptr,
                    bool useLogScale = false, bool bestPDGPlot = false, bool fvPlot = false){
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

    if(fvPlot){
        for(auto* hist : allHists){
            if(!hist) continue;
            hist->GetXaxis()->SetBinLabel(1, "No");
            hist->GetXaxis()->SetBinLabel(2, "Yes");
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
    legend->AddEntry(hists.nu_eDirt, "Dirt #nu+e Elastic Scatter", "f");
    legend->AddEntry(hists.cosmic, "Cosmic", "f");
    legend->AddEntry(hists.other, "Other", "f");
    legend->AddEntry(hists.nu_eFuzzy, "Mis-Reco #nu+e Elastic Scatter", "f");
    legend->SetTextSize(0.0225);

    legend->SetMargin(0.13);
    legend->Draw();

    hists.canvas->SaveAs(filename);
    saveCanvasToRootFile(hists.canvas, rootname);

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

TH1D* makeTotalHist(const TH1D* h){
    TH1D* hc = (TH1D*)h->Clone(Form("%s_totalSum", h->GetName()));
    hc->Reset();

    double totalSum = h->Integral(0, h->GetNbinsX() + 1);
    //std::cout << "total sum = " << totalSum << std::endl; 
    
    for(int i = 0; i <= hc->GetNbinsX() + 1; ++i){
        hc->SetBinContent(i, totalSum);
    }

    return hc;
}

TH1D* makeCumulative(const TH1D* h, bool keepRight){
    TH1D* hc = (TH1D*)h->Clone(Form("%s_cumulative", h->GetName()));
    hc->Reset();

    int n = h->GetNbinsX();

    if (keepRight) {
        double sum = 0.0;
        for (int i = n; i >= 1; --i) {
            sum += h->GetBinContent(i);
            hc->SetBinContent(i, sum);
        }
    } else {
        double sum = 0.0;
        for (int i = 1; i <= n; ++i) {
            sum += h->GetBinContent(i);
            hc->SetBinContent(i, sum);
        }
    }

    return hc;
}

double getMinValueEfficiency(const TEfficiency* eff, double xmin, double xmax, bool includeErrors = false){
    if (!eff) return 0.0;

    const TH1* hTotConst = eff->GetTotalHistogram();
    if (!hTotConst) return 0.0;

    TH1* hTot = const_cast<TH1*>(hTotConst);

    int binMin = hTot->FindBin(xmin);
    int binMax = hTot->FindBin(xmax);

    binMin = std::max(binMin, 1);
    binMax = std::min(binMax, hTot->GetNbinsX());

    double minVal = std::numeric_limits<double>::max();

    for (int i = binMin; i <= binMax; ++i) {
        if (!hTot->GetBinContent(i)) continue;

        double val = eff->GetEfficiency(i);
        if (includeErrors)
            val -= eff->GetEfficiencyErrorLow(i);

        if (val < minVal)
            minVal = val;
    }

    if (minVal == std::numeric_limits<double>::max())
        return 0.0;

    return minVal;
}

double getMaxValueEfficiency(const TEfficiency* eff, bool includeErrors = false){
    double maxVal = 0.0;

    int nBins = eff->GetTotalHistogram()->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        if (!eff->GetTotalHistogram()->GetBinContent(i)) continue;

        double val = eff->GetEfficiency(i);
        if (includeErrors)
            val += eff->GetEfficiencyErrorUp(i);

        if (val > maxVal)
            maxVal = val;
    }
    return maxVal;
}

void drawEfficiencyErrors(TEfficiency* plot, const std::string& filename, const std::string& rootname, double lowY, double highY, const std::string& legendLocation, double xmin, double xmax, double efficiencyWay = 0.0){
    if (!plot) {
        std::cerr << "drawEfficiency: null TEfficiency pointer\n";
        return;
    }

    double maxVal = getMaxValueEfficiency(plot, false);
    double minVal = getMinValueEfficiency(plot, xmin, xmax, false);

    TCanvas* c = new TCanvas("c_eff", "Efficiency comparison", 800, 600);
    c->SetTicks();
    c->SetLeftMargin(0.15);

    plot->SetMarkerColor(kBlack);
    plot->SetMarkerSize(0.7); 
    plot->SetLineWidth(1);
    plot->SetLineColor(kBlack);
    plot->SetMarkerStyle(20);

    const TH1* hTotal = plot->GetTotalHistogram();
    int nBins = hTotal->GetNbinsX();
    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(nBins);    

    double maxEff = 0;
    double maxEffBin = 0;

    for(int i = 1; i <= nBins; ++i){
        double xCenter = hTotal->GetXaxis()->GetBinCenter(i);
        double xErr = (hTotal->GetXaxis()->GetBinUpEdge(i) - hTotal->GetXaxis()->GetBinLowEdge(i)) / 2.0;
        
        double yEff = plot->GetEfficiency(i);
        double yErrLow  = plot->GetEfficiencyErrorLow(i);
        double yErrUp   = plot->GetEfficiencyErrorUp(i);

        if(yEff == maxEff){
            // The eff x pur value is the same as the max value
            if(efficiencyWay == -1) maxEffBin = (xCenter - xErr); 
        } else if(yEff > maxEff){
            maxEff = yEff;
            if(efficiencyWay == 1) maxEffBin = (xCenter + xErr);
            if(efficiencyWay == -1) maxEffBin = (xCenter - xErr);
        }

        gEff->SetPoint(i-1, xCenter, yEff);
        gEff->SetPointError(i-1, xErr, xErr, yErrLow, yErrUp);
    }

    gEff->SetLineColor(kBlack);
    gEff->SetMarkerColor(kBlack);
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(0.7);
    gEff->SetLineWidth(1);
    if(xmin != 999){
        gEff->GetXaxis()->SetLimits(xmin, xmax);
    }
    gEff->GetYaxis()->SetRangeUser(minVal*0.9, maxVal*1.1);

    const TH1* hAxis = plot->GetTotalHistogram();
    gEff->SetTitle(plot->GetTitle());
    gEff->GetXaxis()->SetTitle(hAxis->GetXaxis()->GetTitle());
    gEff->GetYaxis()->SetTitle(hAxis->GetYaxis()->GetTitle());
    gEff->GetYaxis()->SetTitleOffset(1.6);
    gEff->Draw("AP");

    plot->Draw("SAME");
    gPad->Update();

    auto* gBDT = plot->GetPaintedGraph();
    gBDT->SetMarkerSize(0.8);
    gBDT->Draw("PE SAME");

    auto* g = plot->GetPaintedGraph();
    
    if(lowY == -999999 && highY == -999999){
        g->GetYaxis()->SetRangeUser(minVal*0.9, maxVal*1.1);
    } else{
        g->GetYaxis()->SetRangeUser(lowY, highY);
    }

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    if(legendLocation == "topRight"){ Lxmin=0.69; Lymax=0.863; Lxmax=0.87; Lymin=0.74; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.31; Lymin=0.74; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.69; Lymax=0.26; Lxmax=0.87; Lymin=0.137; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymax=0.26; Lxmax=0.31; Lymin=0.137; }

    TLegend* leg = new TLegend(Lxmin, Lymax, Lxmax, Lymin);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(plot, "DL Nu+E", "LEP");

    //leg->Draw();

    c->SaveAs(filename.c_str());
    saveCanvasToRootFile(hists.canvas, rootname.c_str());
    delete c;
}

void drawTEfficiency(TH1D* numerator, TH1D* denominator, const std::string& filename, const std::string& rootname){
    if (!numerator || !denominator) return;

    TEfficiency* eff = new TEfficiency(*numerator, *denominator);
    eff->SetStatisticOption(TEfficiency::kFCP);
    eff->SetUseWeightedEvents(false);

    TGraphAsymmErrors* gEff = eff->CreateGraph();
    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(0.8);
    gEff->SetLineWidth(2);
    gEff->SetLineColor(kBlack);

    {
        TCanvas* c = new TCanvas("c1", "Denominator + Efficiency", 800, 600);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.12);
        gPad->SetTicks(1, 0);

        TH1D* denomCopy = (TH1D*)denominator->Clone("denom_copy1");
        denomCopy->SetFillColor(kAzure-4);
        denomCopy->SetLineColor(kAzure-5);
        denomCopy->SetStats(0);

        denomCopy->SetTitle(Form("%s;%s;# of Events", numerator->GetTitle(), numerator->GetXaxis()->GetTitle()));

        double maxVal = 0;
        for(int i = 1; i <= denomCopy->GetNbinsX(); ++i){
            double val = denomCopy->GetBinContent(i);
            if(val > maxVal) maxVal = val;
        }

        denomCopy->SetMaximum(maxVal * 1.3);

        denomCopy->Draw("HIST");

        gPad->Update();
        double ymin = gPad->GetUymin();
        double ymax = gPad->GetUymax();

        TGraphAsymmErrors* gScaled = (TGraphAsymmErrors*)gEff->Clone("gScaled1");
        for (int i = 0; i < gScaled->GetN(); ++i) {
            double x, y;
            gScaled->GetPoint(i, x, y);

            double yScaled = ymin + y * (ymax - ymin);
            double eyl = gScaled->GetErrorYlow(i) * (ymax - ymin);
            double eyh = gScaled->GetErrorYhigh(i) * (ymax - ymin);

            gScaled->SetPoint(i, x, yScaled);
            gScaled->SetPointError(i, gScaled->GetErrorXlow(i), gScaled->GetErrorXhigh(i), eyl, eyh);
        }

        gScaled->Draw("PE SAME");

        TGaxis* axis = new TGaxis(gPad->GetUxmax(), ymin, gPad->GetUxmax(), ymax, 0.0, 1.0, 510, "+L");

        axis->SetTitle("Efficiency");
        axis->SetTitleOffset(1.2);
        axis->SetLabelFont(42);
        axis->SetTitleFont(42);
        axis->SetLabelSize(0.035);
        axis->SetTitleSize(0.04);
        axis->Draw();

        // Legend
        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(denomCopy, "Before Cuts", "F");
        leg->AddEntry(gScaled, "Efficiency", "PE");
        leg->Draw();
        gPad->RedrawAxis();

        c->SaveAs((filename + "Denominator.pdf").c_str());
        saveCanvasToRootFile(hists.canvas, (rootname + "Denominator.pdf").c_str()

        delete gScaled;
        delete denomCopy;
        delete leg;
        delete c;
    }

    {
        TCanvas* c = new TCanvas("c2", "Efficiency Only", 800, 600);
        c->SetLeftMargin(0.15);

        TH1D* frame = (TH1D*)numerator->Clone("frame");
        frame->Reset();
        frame->SetStats(0);
        frame->SetTitle(Form("%s;%s;Efficiency", numerator->GetTitle(), numerator->GetXaxis()->GetTitle()));
        frame->SetMinimum(0.0);
        frame->SetMaximum(1.0);
        frame->Draw();

        gEff->Draw("PE SAME");

        // Legend
        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gEff, "Efficiency", "PE");
        leg->Draw();
        gPad->RedrawAxis();

        c->SaveAs((filename + ".pdf").c_str());
        saveCanvasToRootFile(hists.canvas, (rootname).c_str()

        delete frame;
        delete leg;
        delete c;
    }

    {
        TCanvas* c = new TCanvas("c3", "Num + Denom + Efficiency", 800, 600);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.12);
        gPad->SetTicks(1, 0);

        TH1D* denomCopy = (TH1D*)denominator->Clone("denom_copy3");
        TH1D* numCopy   = (TH1D*)numerator->Clone("num_copy3");

        denomCopy->SetFillColor(kAzure-4);
        denomCopy->SetLineColor(kAzure-5);
        numCopy->SetFillColor(kOrange+6);
        numCopy->SetLineColor(kOrange+7);

        denomCopy->SetStats(0);
        numCopy->SetStats(0);

        denomCopy->SetTitle(Form("%s;%s;# of Events", numerator->GetTitle(), numerator->GetXaxis()->GetTitle()));

        double maxVal = 0;
        for(int i = 1; i <= denomCopy->GetNbinsX(); ++i){
            double val = denomCopy->GetBinContent(i);
            if(val > maxVal) maxVal = val;
        }

        denomCopy->SetMaximum(maxVal * 1.3);

        denomCopy->Draw("HIST");
        numCopy->Draw("HIST SAME");

        gPad->Update();
        double ymin = gPad->GetUymin();
        double ymax = gPad->GetUymax();

        // Scale efficiency
        TGraphAsymmErrors* gScaled = (TGraphAsymmErrors*)gEff->Clone("gScaled3");
        for (int i = 0; i < gScaled->GetN(); ++i) {
            double x, y;
            gScaled->GetPoint(i, x, y);

            double yScaled = ymin + y * (ymax - ymin);
            double eyl = gScaled->GetErrorYlow(i) * (ymax - ymin);
            double eyh = gScaled->GetErrorYhigh(i) * (ymax - ymin);

            gScaled->SetPoint(i, x, yScaled);
            gScaled->SetPointError(i,
                gScaled->GetErrorXlow(i),
                gScaled->GetErrorXhigh(i),
                eyl, eyh);
        }

        gScaled->Draw("PE SAME");

        TGaxis* axis = new TGaxis(
            gPad->GetUxmax(), ymin,
            gPad->GetUxmax(), ymax,
            0.0, 1.0,
            510, "+L"
        );

        axis->SetTitle("Efficiency");
        axis->SetTitleOffset(1.2);
        axis->SetLabelFont(42);
        axis->SetTitleFont(42);
        axis->SetLabelSize(0.035);
        axis->SetTitleSize(0.04);
        axis->Draw();

        // Legend
        TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(denomCopy, "Before Cuts", "F");
        leg->AddEntry(numCopy, "After Cuts", "F");
        leg->AddEntry(gScaled, "Efficiency", "PE");
        leg->Draw();
        gPad->RedrawAxis();

        c->SaveAs((filename + "DenominatorNumerator.pdf").c_str());
        saveCanvasToRootFile(hists.canvas, (rootname).c_str()));

        delete gScaled;
        delete denomCopy;
        delete numCopy;
        delete leg;
        delete c;
    }

    delete gEff;
    delete eff;
}

void drawEffPurEffPurCombined(TEfficiency* eff, TEfficiency* pur, TEfficiency* effPur,
                               double optimalCut, const std::string& filename, const std::string& rootname,
                               const std::string& legendLocation,
                               double xmin, double xmax){
    if(!eff || !pur || !effPur){
        std::cerr << "drawEffPurEffPurCombined: null TEfficiency pointer\n";
        return;
    }

    TCanvas* c = new TCanvas("c_effPurCombined", "Efficiency, Purity, Efficiency x Purity", 800, 600);
    c->SetTicks();
    c->SetLeftMargin(0.15);

    const TH1* hTotal = eff->GetTotalHistogram();
    int nBins = hTotal->GetNbinsX();

    auto buildScaledGraph = [&](TEfficiency* effObj) -> TGraphAsymmErrors* {
        TGraphAsymmErrors* g = new TGraphAsymmErrors(nBins);
        for(int i = 1; i <= nBins; ++i){
            double xCenter = hTotal->GetXaxis()->GetBinCenter(i);
            double xErr = (hTotal->GetXaxis()->GetBinUpEdge(i) - hTotal->GetXaxis()->GetBinLowEdge(i)) / 2.0;

            double y = effObj->GetEfficiency(i) * 100.0;
            double yErrLow = effObj->GetEfficiencyErrorLow(i) * 100.0;
            double yErrUp = effObj->GetEfficiencyErrorUp(i) * 100.0;

            g->SetPoint(i-1, xCenter, y);
            g->SetPointError(i-1, xErr, xErr, yErrLow, yErrUp);
        }
        return g;
    };

    TGraphAsymmErrors* gEff = buildScaledGraph(eff);
    TGraphAsymmErrors* gPur = buildScaledGraph(pur);
    TGraphAsymmErrors* gEffPur = buildScaledGraph(effPur);

    gEff->SetMarkerColor(kRed+1);      gEff->SetLineColor(kRed+1);
    gEff->SetMarkerStyle(20);          gEff->SetMarkerSize(0.7);      gEff->SetLineWidth(1);

    gPur->SetMarkerColor(kBlue+1);     gPur->SetLineColor(kBlue+1);
    gPur->SetMarkerStyle(20);          gPur->SetMarkerSize(0.7);      gPur->SetLineWidth(1);

    gEffPur->SetMarkerColor(kMagenta+2); gEffPur->SetLineColor(kMagenta+2);
    gEffPur->SetMarkerStyle(20);         gEffPur->SetMarkerSize(0.7);   gEffPur->SetLineWidth(1);

    if(xmin != 999 && xmax != 999){
        gEff->GetXaxis()->SetLimits(xmin, xmax);
    }

    gEff->SetTitle(hTotal->GetTitle());
    gEff->GetXaxis()->SetTitle(hTotal->GetXaxis()->GetTitle());
    gEff->GetYaxis()->SetTitle("Percentage (%)");
    gEff->GetYaxis()->SetRangeUser(0, 130);
    gEff->GetYaxis()->SetTitleOffset(1.4);

    gEff->Draw("AP");

    gPad->Update();
    double yAxisMinRight = gPad->GetUymin();
    double yAxisMaxRight = gPad->GetUymax();

    TGaxis* rightAxis = new TGaxis(gPad->GetUxmax(), yAxisMinRight, gPad->GetUxmax(), yAxisMaxRight, yAxisMinRight, yAxisMaxRight, 510, "+L");
    rightAxis->SetTitle("Percentage (%)");
    rightAxis->SetTitleOffset(1.4);
    rightAxis->SetLabelFont(42);
    rightAxis->SetTitleFont(42);
    rightAxis->SetLabelSize(0.035);
    rightAxis->SetTitleSize(0.04);
    rightAxis->Draw();

    gPur->Draw("P SAME");
    gEffPur->Draw("P SAME");

    // Find the bin for the optimal cut and read off eff/pur/effPur there
    TH1* hTotalNonConst = const_cast<TH1*>(hTotal);
    int optimalBin = hTotalNonConst->FindBin(optimalCut);
    if(optimalBin < 1) optimalBin = 1;
    if(optimalBin > nBins) optimalBin = nBins;

    double effAtOptimal    = eff->GetEfficiency(optimalBin) * 100.0;
    double purAtOptimal    = pur->GetEfficiency(optimalBin) * 100.0;
    double effPurAtOptimal = effPur->GetEfficiency(optimalBin) * 100.0;

    gPad->Update();
    double yAxisMin = gPad->GetUymin();
    double yAxisMax = gPad->GetUymax();

    TLine* line = new TLine(optimalCut, yAxisMin, optimalCut, yAxisMax);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    if(legendLocation == "topRight"){ Lxmin=0.40; Lymax=0.88; Lxmax=0.87; Lymin=0.80; }
    else if(legendLocation == "topLeft"){ Lxmin=0.15; Lymax=0.88; Lxmax=0.62; Lymin=0.80; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.40; Lymax=0.26; Lxmax=0.87; Lymin=0.137; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.15; Lymax=0.26; Lxmax=0.62; Lymin=0.137; }

    TLegend* leg = new TLegend(Lxmin, Lymin, Lxmax, Lymax);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(3);
    leg->AddEntry(gEff, "Efficiency", "PE");
    leg->AddEntry(gPur, "Purity", "PE");
    leg->AddEntry(gEffPur, "Efficiency x Purity", "PE");
    leg->Draw();

    TLatex* text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextAlign(11);

    double xText = 0.16;
    double yText = 0.72;
    double yStep = 0.045;

    text->DrawLatex(xText, yText,           Form("Optimal Cut: %.4f", optimalCut));
    text->DrawLatex(xText, yText - yStep,   Form("Efficiency x Purity: %.2f", effPurAtOptimal));
    text->DrawLatex(xText, yText - 2*yStep, Form("Efficiency: %.2f %%", effAtOptimal));
    text->DrawLatex(xText, yText - 3*yStep, Form("Purity: %.2f %%", purAtOptimal));

    gPad->RedrawAxis();
    c->SaveAs(filename.c_str());
    saveCanvasToRootFile(hists.canvas, rootname.c_str());

    delete gEff;
    delete gPur;
    delete gEffPur;
    delete rightAxis;
    delete line;
    delete leg;
    delete text;
    delete c;
}

void efficiency(double trueSignal, histGroup_struct* histBeforeCuts, histGroup_struct* histAfterCuts, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, double optimalCut, int* drawLine = nullptr, int* linePos = nullptr, double efficiencyWay = 0.0){
    bool keepRight = (efficiencyWay == -1);

    // Total signal before cuts
    //std::cout << "============================ total signal before cuts ============================" << std::endl;
    TH1D* hTotalSignalBeforeCuts = makeTotalHist(histBeforeCuts->nuESignal);
    //std::cout << "======================================================================================" << std::endl;

    TH1D* hTotalSignalNewDefinition = (TH1D*)histBeforeCuts->nuESignal->Clone("hTotalSignalNewDefinition");
    hTotalSignalNewDefinition->Reset();

    for(int i = 1; i <= hTotalSignalNewDefinition->GetNbinsX(); i++){
        hTotalSignalNewDefinition->SetBinContent(i, trueSignal);
    }

    // Total signal kept after cuts (cumulative)
    //std::cout << "============================ total signal kept after cuts ============================" << std::endl;
    TH1D* hPassedSignalAfterCuts = makeCumulative(histAfterCuts->nuESignal, keepRight);
    //std::cout << "======================================================================================" << std::endl;

    // Total background before cuts
    //std::cout << "============================ total background before cuts ============================" << std::endl;
    TH1D* hTotalBackgroundBeforeCutsAdded = (TH1D*) histBeforeCuts->nuECosmic->Clone("hTotalBackgroundBeforeCutsAdded");
    hTotalBackgroundBeforeCutsAdded->Reset();
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuECosmic);
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuEBNB);
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuEBNBFuzzy);
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuESignalFuzzy);
    //std::cout << "Number of entries = " << hTotalBackgroundBeforeCutsAdded->GetEntries() << std::endl;
    TH1D* hTotalBackgroundBeforeCuts = makeTotalHist(hTotalBackgroundBeforeCutsAdded);
    //std::cout << "======================================================================================" << std::endl;

    // Total background kept after cuts (cumulative)
    //std::cout << "============================ total background kept after cuts ============================" << std::endl;
    TH1D* hPassedBackgroundAfterCutsAdded = (TH1D*) histAfterCuts->nuECosmic->Clone("hPassedBackgroundAfterCutsAdded");
    hPassedBackgroundAfterCutsAdded->Reset();
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuECosmic);
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuEBNB);
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuEBNBFuzzy);
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuESignalFuzzy);
    TH1D* hPassedBackgroundAfterCuts = makeCumulative(hPassedBackgroundAfterCutsAdded, keepRight);
    //std::cout << "======================================================================================" << std::endl;

    // Total background rejected after cuts (cumulative)
    //std::cout << "============================ total background rejected after cuts ============================" << std::endl;
    TH1D* hRejectedBackgroundAfterCuts = (TH1D*) hTotalBackgroundBeforeCuts->Clone("hRejectedBackgroundAfterCuts");
    hRejectedBackgroundAfterCuts->Add(hPassedBackgroundAfterCuts, -1.0);
    //std::cout << "======================================================================================" << std::endl;

    // Total background + signal kept after cuts (cumulative)
    //std::cout << "============================ total background + signal after cuts ============================" << std::endl;
    TH1D* hPassedBackgroundSignalAfterCutsAdded = (TH1D*) histAfterCuts->nuECosmic->Clone("hPassedBackgroundSignalAfterCutsAdded");
    hPassedBackgroundSignalAfterCutsAdded->Reset();
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuECosmic);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuEBNB);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuEBNBFuzzy);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuESignalFuzzy);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuESignal);
    TH1D* hPassedBackgroundSignalAfterCuts = makeCumulative(hPassedBackgroundSignalAfterCutsAdded, keepRight);
    //std::cout << "======================================================================================" << std::endl;

    TH1D* hEffPurNumerator = (TH1D*) histAfterCuts->nuESignal->Clone("hEffPurNumerator");
    hEffPurNumerator->Reset();
    hEffPurNumerator->Add(hPassedSignalAfterCuts);
    hEffPurNumerator->Multiply(hPassedSignalAfterCuts);

    TH1D* hEffPurDenominator = (TH1D*) histAfterCuts->nuESignal->Clone("hEffPurDenominator");
    hEffPurDenominator->Reset();
    hEffPurDenominator->Add(hTotalSignalBeforeCuts);
    hEffPurDenominator->Multiply(hPassedBackgroundSignalAfterCuts);

    TH1D* hEffPurDenominator_newDef = (TH1D*) histAfterCuts->nuESignal->Clone("hEffPurDenominator_newDef");
    hEffPurDenominator_newDef->Reset();
    hEffPurDenominator_newDef->Add(hTotalSignalNewDefinition);
    hEffPurDenominator_newDef->Multiply(hPassedBackgroundSignalAfterCuts);

    // Efficiency plot
    TEfficiency* eff = new TEfficiency(*hPassedSignalAfterCuts, *hTotalSignalBeforeCuts);
    TEfficiency* eff_newDef = new TEfficiency(*hPassedSignalAfterCuts, *hTotalSignalNewDefinition);
    TEfficiency* rej = new TEfficiency(*hRejectedBackgroundAfterCuts, *hTotalBackgroundBeforeCuts);
    TEfficiency* pur = new TEfficiency(*hPassedSignalAfterCuts, *hPassedBackgroundSignalAfterCuts);
    TEfficiency* effPur = new TEfficiency(*hEffPurNumerator, *hEffPurDenominator);
    TEfficiency* effPur_newDef = new TEfficiency(*hEffPurNumerator, *hEffPurDenominator_newDef);

    eff->SetTitle(Form("%s;%s;Efficiency", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    eff->SetStatisticOption(TEfficiency::kFNormal);
    
    eff_newDef->SetTitle(Form("%s;%s;Efficiency", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    eff_newDef->SetStatisticOption(TEfficiency::kFNormal);
    
    pur->SetTitle(Form("%s;%s;Purity", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    pur->SetStatisticOption(TEfficiency::kFNormal);
    
    rej->SetTitle(Form("%s;%s;Background Rejection", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    rej->SetStatisticOption(TEfficiency::kFNormal);
    
    effPur->SetTitle(Form("%s;%s;Efficiency x Purity", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    effPur->SetStatisticOption(TEfficiency::kFNormal);
    
    effPur_newDef->SetTitle(Form("%s;%s;Efficiency x Purity", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    effPur_newDef->SetStatisticOption(TEfficiency::kFNormal);

    std::string filenameEff = std::string(filename) + "_eff.pdf";
    std::string filenameEff_newDef = std::string(filename) + "_effNewDef.pdf";
    std::string filenamePur = std::string(filename) + "_pur.pdf";
    std::string filenameRej = std::string(filename) + "_rej.pdf";
    std::string filenameEffPur = std::string(filename) + "_effPur.pdf";
    std::string filenameEffPur_newDef = std::string(filename) + "_effPurNewDef.pdf";

    eff->SetUseWeightedEvents(false);
    eff_newDef->SetUseWeightedEvents(false);
    pur->SetUseWeightedEvents(false);
    rej->SetUseWeightedEvents(false);
    effPur->SetUseWeightedEvents(false);
    effPur_newDef->SetUseWeightedEvents(false);
   
    // Comment out if not wanting individual plots 
    drawEfficiencyErrors(eff, filenameEff, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay);
    drawEfficiencyErrors(eff_newDef, filenameEff_newDef, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay);
    drawEfficiencyErrors(pur, filenamePur, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay);
    drawEfficiencyErrors(rej, filenameRej, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay);
    drawEfficiencyErrors(effPur, filenameEffPur, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay);
    drawEfficiencyErrors(effPur_newDef, filenameEffPur_newDef, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay);

    std::string filenameCombined = std::string(filename) + "_effPurCombined.pdf";
    drawEffPurEffPurCombined(eff_newDef, pur, effPur_newDef, optimalCut, filenameCombined, legendLocation, xmin, xmax);
}

void allSelectionPlots_macro(){

    TChain *subRunTree = new TChain("ana/SubRun");
    subRunTree->Add("/exp/sbnd/data/users/coackley/analysisFiles_14Jul/*.root");

    TChain *tree = new TChain("ana/NuE");
    tree->Add("/exp/sbnd/data/users/coackley/analysisFiles_14Jul/*.root");

    std::string base_path = "/nashome/c/coackley/selectionThesisPlots_26Jul/";

    gROOT->SetBatch(true);

    outRootFile = new TFile((base_path + "allPlots.root").c_str(), "RECREATE");

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

    // Active volume boundaries:
    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0;      double zMax = 509.4;

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
 
    auto sliceRecoVXBeforeCuts = createHistGroup("sliceRecoVXBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto sliceRecoVXBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm) ", 202, -202, 202);
    auto sliceRecoVXAfterCuts = createHistGroup("sliceRecoVXAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto sliceRecoVXAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm) ", 202, -202, 202);
   
    auto sliceRecoVYBeforeCuts = createHistGroup("sliceRecoVYBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto sliceRecoVYBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm) ", 204, -204, 204);
    auto sliceRecoVYAfterCuts = createHistGroup("sliceRecoVYAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto sliceRecoVYAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm) ", 204, -204, 204);
    
    auto sliceRecoVZBeforeCuts = createHistGroup("sliceRecoVZBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto sliceRecoVZBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 255, 0, 510);
    auto sliceRecoVZAfterCuts = createHistGroup("sliceRecoVZAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto sliceRecoVZAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 255, 0, 510);
    
    auto sliceRecoVXSmallerBinsBeforeCuts = createHistGroup("sliceRecoVXSmallerBinsBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto sliceRecoVXSmallerBinsBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXSmallerBinsBeforeCutsBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm) ", 808, -202, 202);
    auto sliceRecoVXSmallerBinsAfterCuts = createHistGroup("sliceRecoVXSmallerBinsAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto sliceRecoVXSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXSmallerBinsAfterCutsAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm) ", 808, -202, 202);
    
    auto sliceRecoVYSmallerBinsBeforeCuts = createHistGroup("sliceRecoVYSmallerBinsBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto sliceRecoVYSmallerBinsBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYSmallerBinsBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm) ", 816, -204, 204);
    auto sliceRecoVYSmallerBinsAfterCuts = createHistGroup("sliceRecoVYSmallerBinsAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto sliceRecoVYSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYSmallerBinsAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm) ", 816, -204, 204);
    
    auto sliceRecoVZSmallerBinsBeforeCuts = createHistGroup("sliceRecoVZSmallerBinsBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto sliceRecoVZSmallerBinsBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZSmallerBinsBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 1020, 0, 510);
    auto sliceRecoVZSmallerBinsAfterCuts = createHistGroup("sliceRecoVZSmallerBinsAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto sliceRecoVZSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZSmallerBinsAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 1020, 0, 510);
    
    auto sliceRecoVXLowBeforeCuts = createHistGroup("sliceRecoVXLowBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, 182);
    auto sliceRecoVXLowBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXLowBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm) ", 40, -202, -182);
    auto sliceRecoVXLowAfterCuts = createHistGroup("sliceRecoVXLowAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto sliceRecoVXLowAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXLowAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm) ", 40, -202, -182);
    
    auto sliceRecoVYLowBeforeCuts = createHistGroup("sliceRecoVYLowBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto sliceRecoVYLowBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYLowBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm) ", 40, -204, -184);
    auto sliceRecoVYLowAfterCuts = createHistGroup("sliceRecoVYLowAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto sliceRecoVYLowAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYLowAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm) ", 40, -204, -184);
    
    auto sliceRecoVZLowBeforeCuts = createHistGroup("sliceRecoVZLowBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto sliceRecoVZLowBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZLowBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 40, 0, 20);
    auto sliceRecoVZLowAfterCuts = createHistGroup("sliceRecoVZLowAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto sliceRecoVZLowAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZLowAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 40, 0, 20);
    
    auto sliceRecoVXHighBeforeCuts = createHistGroup("sliceRecoVXHighBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto sliceRecoVXHighBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXHighBeforeCuts_splitDLNuE", "X Coordinate of Reco Neutrino (Before Cuts)", "x_{Reco} (cm)", 40, 182, 202);
    auto sliceRecoVXHighAfterCuts = createHistGroup("sliceRecoVXHighAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto sliceRecoVXHighAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVXHighAfterCuts_splitDLNuE", "X Coordinate of Reco Neutrino (After Cuts)", "x_{Reco} (cm)", 40, 182, 202);
    
    auto sliceRecoVYHighBeforeCuts = createHistGroup("sliceRecoVYHighBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto sliceRecoVYHighBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYHighBeforeCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (Before Cuts)", "y_{Reco} (cm)", 40, 184, 204);
    auto sliceRecoVYHighAfterCuts = createHistGroup("sliceRecoVYHighAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto sliceRecoVYHighAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVYHighAfterCuts_splitDLNuE", "Y Coordinate of Reco Neutrino (After Cuts)", "y_{Reco} (cm)", 40, 184, 204);
    
    auto sliceRecoVZHighBeforeCuts = createHistGroup("sliceRecoVZHighBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto sliceRecoVZHighBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZHighBeforeCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (Before Cuts)", "z_{Reco} (cm) ", 40, 480, 510);
    auto sliceRecoVZHighAfterCuts = createHistGroup("sliceRecoVZHighAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto sliceRecoVZHighAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoVZHighAfterCuts_splitDLNuE", "Z Coordinate of Reco Neutrino (After Cuts)", "z_{Reco} (cm) ", 40, 480, 510);
    
    auto sliceRecoNeutInFVBeforeCuts = createHistGroup("sliceRecoNeutInFVBeforeCuts", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
    auto sliceRecoNeutInFVBeforeCuts_splitDLNuE = createSplitHistGroup("sliceRecoNeutInFVBeforeCuts_splitDLNuE", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
    auto sliceRecoNeutInFVAfterCuts = createHistGroup("sliceRecoNeutInFVAfterCuts", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
    auto sliceRecoNeutInFVAfterCuts_splitDLNuE = createSplitHistGroup("sliceRecoNeutInFVAfterCuts_splitDLNuE", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
   
    auto energyAsymmetryBeforeCuts = createHistGroup("energyAsymmetryBeforeCuts", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (Before Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryBeforeCuts_splitDLNuE = createSplitHistGroup("energyAsymmetryBeforeCut_splitDLNuE", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (Before Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts = createHistGroup("energyAsymmetryAfterCuts", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (After Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts_splitDLNuE = createSplitHistGroup("energyAsymmetryAfterCut_splitDLNuE", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (After Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
     

    // Put Plots Here (Dedicated After Previous Cut)
    auto sliceNumPFPsAfterClearCosmicCut = createHistGroup("sliceNumPFPsAfterClearCosmicCut", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterClearCosmicCut_splitDLNuE = createSplitHistGroup("sliceNumPFPsAfterClearCosmicCut_splitDLNuE", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterNumPFPCut = createHistGroup("sliceNumPFPsAfterNumPFPCut", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterNumPFPCut_splitDLNuE = createSplitHistGroup("sliceNumPFPsAfterNumPFPCut_splitDLNuE", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    
    auto sliceNumRecoNeutAfterNumPFPCut = createHistGroup("sliceNumRecoNeutAfterNumPFPCut", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterNumPFPCut_splitDLNuE = createSplitHistGroup("sliceNumRecoNeutAfterNumPFPCut_splitDLNuE", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterNumNeutrinoCut = createHistGroup("sliceNumRecoNeutAfterNumNeutrinoCut", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterNumNeutrinoCut_splitDLNuE = createSplitHistGroup("sliceNumRecoNeutAfterNumNeutrinoCut_splitDLNuE", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    
    auto sliceCRUMBSAfterNumNeutrinoCut = createHistGroup("sliceCRUMBSAfterNumNeutrinoCut", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterNumNeutrinoCut_splitDLNuE = createSplitHistGroup("sliceCRUMBSAfterNumNeutrinoCut_splitDLNuE", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCRUMBSCut = createHistGroup("sliceCRUMBSAfterCRUMBSCut", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceCRUMBSAfterCRUMBSCut_splitDLNuE", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    
    auto sliceRecoVXAfterCRUMBSCut = createHistGroup("sliceRecoVXAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto sliceRecoVXAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVXAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto sliceRecoVXAfterFVCut = createHistGroup("sliceRecoVXAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto sliceRecoVXAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVXAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
   
    auto sliceRecoVYAfterCRUMBSCut = createHistGroup("sliceRecoVYAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto sliceRecoVYAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVYAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto sliceRecoVYAfterFVCut = createHistGroup("sliceRecoVYAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto sliceRecoVYAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVYAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    
    auto sliceRecoVZAfterCRUMBSCut = createHistGroup("sliceRecoVZAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto sliceRecoVZAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVZAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto sliceRecoVZAfterFVCut = createHistGroup("sliceRecoVZAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto sliceRecoVZAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVZAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    
    auto sliceRecoVXSmallerBinsAfterCRUMBSCut = createHistGroup("sliceRecoVXSmallerBinsAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto sliceRecoVXSmallerBinsAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVXSmallerBinsAfterCRUMBSCutAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto sliceRecoVXSmallerBinsAfterFVCut = createHistGroup("sliceRecoVXSmallerBinsAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto sliceRecoVXSmallerBinsAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVXSmallerBinsAfterFVCutAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    
    auto sliceRecoVYSmallerBinsAfterCRUMBSCut = createHistGroup("sliceRecoVYSmallerBinsAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto sliceRecoVYSmallerBinsAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVYSmallerBinsAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto sliceRecoVYSmallerBinsAfterFVCut = createHistGroup("sliceRecoVYSmallerBinsAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto sliceRecoVYSmallerBinsAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVYSmallerBinsAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    
    auto sliceRecoVZSmallerBinsAfterCRUMBSCut = createHistGroup("sliceRecoVZSmallerBinsAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto sliceRecoVZSmallerBinsAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVZSmallerBinsAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto sliceRecoVZSmallerBinsAfterFVCut = createHistGroup("sliceRecoVZSmallerBinsAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto sliceRecoVZSmallerBinsAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVZSmallerBinsAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    
    auto sliceRecoVXLowAfterCRUMBSCut = createHistGroup("sliceRecoVXLowAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto sliceRecoVXLowAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVXLowAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto sliceRecoVXLowAfterFVCut = createHistGroup("sliceRecoVXLowAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto sliceRecoVXLowAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVXLowAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    
    auto sliceRecoVYLowAfterCRUMBSCut = createHistGroup("sliceRecoVYLowAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto sliceRecoVYLowAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVYLowAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto sliceRecoVYLowAfterFVCut = createHistGroup("sliceRecoVYLowAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto sliceRecoVYLowAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVYLowAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    
    auto sliceRecoVZLowAfterCRUMBSCut = createHistGroup("sliceRecoVZLowAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto sliceRecoVZLowAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVZLowAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto sliceRecoVZLowAfterFVCut = createHistGroup("sliceRecoVZLowAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto sliceRecoVZLowAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVZLowAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    
    auto sliceRecoVXHighAfterCRUMBSCut = createHistGroup("sliceRecoVXHighAfterCRUMBSCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto sliceRecoVXHighAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVXHighAfterCRUMBSCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto sliceRecoVXHighAfterFVCut = createHistGroup("sliceRecoVXHighAfterFVCut", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto sliceRecoVXHighAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVXHighAfterFVCut_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    
    auto sliceRecoVYHighAfterCRUMBSCut = createHistGroup("sliceRecoVYHighAfterCRUMBSCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto sliceRecoVYHighAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVYHighAfterCRUMBSCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto sliceRecoVYHighAfterFVCut = createHistGroup("sliceRecoVYHighAfterFVCut", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto sliceRecoVYHighAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVYHighAfterFVCut_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    
    auto sliceRecoVZHighAfterCRUMBSCut = createHistGroup("sliceRecoVZHighAfterCRUMBSCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto sliceRecoVZHighAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoVZHighAfterCRUMBSCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto sliceRecoVZHighAfterFVCut = createHistGroup("sliceRecoVZHighAfterFVCut", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto sliceRecoVZHighAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoVZHighAfterFVCut_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    
    auto sliceRecoNeutInFVAfterCRUMBSCut = createHistGroup("sliceRecoNeutInFVAfterCRUMBSCut", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
    auto sliceRecoNeutInFVAfterCRUMBSCut_splitDLNuE = createSplitHistGroup("sliceRecoNeutInFVAfterCRUMBSCut_splitDLNuE", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
    auto sliceRecoNeutInFVAfterFVCut = createHistGroup("sliceRecoNeutInFVAfterFVCut", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
    auto sliceRecoNeutInFVAfterFVCut_splitDLNuE = createSplitHistGroup("sliceRecoNeutInFVAfterFVCut_splitDLNuE", "Reco Neutrino Within FV", "In FV?", 2, 0, 2);
 
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

        int trueSignal = 0;

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
            fillHistogram(&sliceRecoVXBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXSmallerBinsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYSmallerBinsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZSmallerBinsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXLowBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYLowBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZLowBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXHighBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYHighBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZHighBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
           
            fillSplitIntHistogram(&sliceCompletenessBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceCompleteness->at(slice), &weights);
            fillSplitIntHistogram(&slicePurityBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_slicePurity->at(slice), &weights);
            fillSplitIntHistogram(&sliceCRUMBSBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
            fillSplitIntHistogram(&sliceNumRecoNeutBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);
            fillSplitIntHistogram(&sliceNumPFPsBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPFPsSlice_beforeCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_beforeCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsMinHitBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsMinHitSlice_beforeCuts, &weights);
            fillSplitIntHistogram(&ERecoSumThetaRecoBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoBeforeCuts_splitDLNuE_pfp10cmPoints, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_beforeCuts.energy * pfp10cm_PCAAngle_beforeCuts * pfp10cm_PCAAngle_beforeCuts), &weights);
            fillSplitIntHistogram(&dEdxBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
            fillSplitIntHistogram(&razzledPDG11BeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
            fillSplitIntHistogram(&razzledPDG13BeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
            fillSplitIntHistogram(&razzledPDG22BeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
            fillSplitIntHistogram(&razzledPDG211BeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
            fillSplitIntHistogram(&razzledPDG2212BeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
            fillSplitIntHistogram(&pfpCompletenessBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.completeness, &weights);
            fillSplitIntHistogram(&pfpPurityBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.purity, &weights); 
            fillSplitIntHistogram(&sliceRecoVXBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXSmallerBinsBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYSmallerBinsBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZSmallerBinsBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXLowBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYLowBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZLowBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXHighBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYHighBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZHighBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);

            if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999){
                fillHistogram(&angleDifferenceSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifference_beforeCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP10cmSignalBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP10cm_beforeCuts, &weights);
                fillHistogram(&energyAsymmetryBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy)/recoilElectron.energy), &weights);
            }           
 
            // Check whether reco neutrino is in FV
            if(!(sliceRecoVX < FVCut_xHigh && sliceRecoVX > FVCut_xLow  && std::abs(sliceRecoVX) > FVCut_xCentre && sliceRecoVY < FVCut_yHigh && sliceRecoVY > FVCut_yLow && sliceRecoVZ > FVCut_zLow && sliceRecoVZ < FVCut_zHigh)){
                // Reco neutrino isn't in FV
                fillHistogram(&sliceRecoNeutInFVBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, 0.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, 0.5, &weights);
            } else{
                // Reco neutrino is in FV
                fillHistogram(&sliceRecoNeutInFVBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, 1.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVBeforeCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, 1.5, &weights);
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

            fillHistogram(&sliceNumPFPsAfterClearCosmicCut, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPFPsAfterClearCosmicCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPFPsSlice_afterCuts, &weights);

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

            fillHistogram(&sliceNumPFPsAfterNumPFPCut, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPFPsAfterNumPFPCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPFPsSlice_afterCuts, &weights);

            fillHistogram(&sliceNumRecoNeutAfterNumPFPCut, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
            fillSplitIntHistogram(&sliceNumRecoNeutAfterNumPFPCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);

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

            fillHistogram(&sliceNumRecoNeutAfterNumNeutrinoCut, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
            fillSplitIntHistogram(&sliceNumRecoNeutAfterNumNeutrinoCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);

            fillHistogram(&sliceCRUMBSAfterNumNeutrinoCut, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
            fillSplitIntHistogram(&sliceCRUMBSAfterNumNeutrinoCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);

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
            fillSplitIntHistogram(&sliceCRUMBSAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
           
            fillHistogram(&sliceRecoVXAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXSmallerBinsAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYSmallerBinsAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZSmallerBinsAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXLowAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYLowAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZLowAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXHighAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYHighAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZHighAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
           
            fillSplitIntHistogram(&sliceRecoVXAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXSmallerBinsAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYSmallerBinsAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZSmallerBinsAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXLowAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYLowAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZLowAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXHighAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYHighAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZHighAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            
            // Check whether reco neutrino is in FV
            if(!(sliceRecoVX < FVCut_xHigh && sliceRecoVX > FVCut_xLow  && std::abs(sliceRecoVX) > FVCut_xCentre && sliceRecoVY < FVCut_yHigh && sliceRecoVY > FVCut_yLow && sliceRecoVZ > FVCut_zLow && sliceRecoVZ < FVCut_zHigh)){
                // Reco neutrino isn't in FV
                fillHistogram(&sliceRecoNeutInFVAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, 0.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, 0.5, &weights);
            } else{
                // Reco neutrino is in FV
                fillHistogram(&sliceRecoNeutInFVAfterCRUMBSCut, DLCurrent, signal, sliceCategoryPlottingMacro, 1.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVAfterCRUMBSCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, 1.5, &weights);
            }

            if(FVCut == 1){
                if(!(sliceRecoVX < FVCut_xHigh && sliceRecoVX > FVCut_xLow  && std::abs(sliceRecoVX) > FVCut_xCentre && sliceRecoVY < FVCut_yHigh && sliceRecoVY > FVCut_yLow && sliceRecoVZ > FVCut_zLow && sliceRecoVZ < FVCut_zHigh)){
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
          
            fillHistogram(&sliceRecoVXAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXSmallerBinsAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYSmallerBinsAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZSmallerBinsAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXLowAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYLowAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZLowAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXHighAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYHighAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZHighAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
           
            fillSplitIntHistogram(&sliceRecoVXAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXSmallerBinsAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYSmallerBinsAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZSmallerBinsAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXLowAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYLowAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZLowAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXHighAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYHighAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZHighAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            
            // Check whether reco neutrino is in FV
            if(!(sliceRecoVX < FVCut_xHigh && sliceRecoVX > FVCut_xLow  && std::abs(sliceRecoVX) > FVCut_xCentre && sliceRecoVY < FVCut_yHigh && sliceRecoVY > FVCut_yLow && sliceRecoVZ > FVCut_zLow && sliceRecoVZ < FVCut_zHigh)){
                // Reco neutrino isn't in FV
                fillHistogram(&sliceRecoNeutInFVAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, 0.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, 0.5, &weights);
            } else{
                // Reco neutrino is in FV
                fillHistogram(&sliceRecoNeutInFVAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, 1.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, 1.5, &weights);
            }

            fillHistogram(&sliceNumPrimaryPFPsAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_afterCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsMinHitAfterFVCut, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsMinHitSlice_afterCuts, &weights);

            fillSplitIntHistogram(&sliceNumPrimaryPFPsAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsMinHitAfterFVCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsMinHitSlice_afterCuts, &weights);

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

            fillHistogram(&sliceNumPrimaryPFPsAfterPrimaryPFPCut, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_afterCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsMinHitSlice_afterCuts, &weights);

            fillSplitIntHistogram(&sliceNumPrimaryPFPsAfterPrimaryPFPCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsMinHitSlice_afterCuts, &weights);

            fillHistogram(&razzledPDG11AfterPrimaryPFPCut, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
            fillSplitIntHistogram(&razzledPDG11AfterPrimaryPFPCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG11, &weights);

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

            fillHistogram(&razzledPDG11AfterRazzled11Cut, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
            fillSplitIntHistogram(&razzledPDG11AfterRazzled11Cut_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG11, &weights);

            fillHistogram(&razzledPDG211AfterRazzled11Cut, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
            fillSplitIntHistogram(&razzledPDG211AfterRazzled11Cut_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG211, &weights);

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
            
            fillHistogram(&razzledPDG211AfterRazzled211Cut, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
            fillSplitIntHistogram(&razzledPDG211AfterRazzled211Cut_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG211, &weights);

            fillHistogram(&ERecoHighestThetaRecoAfterRazzled211Cut, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE_pfp10cmPoints, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);
 
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

            fillHistogram(&ERecoHighestThetaRecoAfterETheta2Cut, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoAfterETheta2Cut_pfp10cmPoints, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE_pfp10cmPoints, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);

            fillHistogram(&dEdxAfterETheta2Cut, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
            fillSplitIntHistogram(&dEdxAfterETheta2Cut_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);

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
            
            fillHistogram(&dEdxAfterdEdxCut, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
            fillSplitIntHistogram(&dEdxAfterdEdxCut_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);

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
            fillHistogram(&sliceRecoVXAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
            fillHistogram(&sliceRecoVXHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVX, &weights);
            fillHistogram(&sliceRecoVYHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVY, &weights);
            fillHistogram(&sliceRecoVZHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, sliceRecoVZ, &weights);
           
            fillSplitIntHistogram(&sliceCompletenessAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceCompleteness->at(slice), &weights);
            fillSplitIntHistogram(&slicePurityAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_slicePurity->at(slice), &weights);
            fillSplitIntHistogram(&sliceCRUMBSAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
            fillSplitIntHistogram(&sliceNumRecoNeutAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);
            fillSplitIntHistogram(&sliceNumPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_afterCuts, &weights);
            fillSplitIntHistogram(&sliceNumPrimaryPFPsMinHitAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsMinHitSlice_afterCuts, &weights);
            fillSplitIntHistogram(&ERecoSumThetaRecoAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (summedEnergy_afterCuts * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
            fillSplitIntHistogram(&ERecoHighestThetaRecoAfterCuts_splitDLNuE_pfp10cmPoints, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * pfp10cm_PCAAngle_afterCuts * pfp10cm_PCAAngle_afterCuts), &weights);
            fillSplitIntHistogram(&dEdxAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
            fillSplitIntHistogram(&razzledPDG11AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
            fillSplitIntHistogram(&razzledPDG13AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG13, &weights);
            fillSplitIntHistogram(&razzledPDG22AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG22, &weights);
            fillSplitIntHistogram(&razzledPDG211AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
            fillSplitIntHistogram(&razzledPDG2212AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG2212, &weights);
            fillSplitIntHistogram(&pfpCompletenessAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.completeness, &weights);
            fillSplitIntHistogram(&pfpPurityAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.purity, &weights); 
            fillSplitIntHistogram(&sliceRecoVXAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);
            fillSplitIntHistogram(&sliceRecoVXHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVX, &weights);
            fillSplitIntHistogram(&sliceRecoVYHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVY, &weights);
            fillSplitIntHistogram(&sliceRecoVZHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, sliceRecoVZ, &weights);

            if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999){
                fillHistogram(&angleDifferenceSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifference_afterCuts, &weights);
                fillHistogram(&angleDifferencePCAPFP10cmSignalAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifferencePCAPFP10cm_afterCuts, &weights);
                fillHistogram(&energyAsymmetryAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy)/recoilElectron.energy), &weights);
            }           
            
            // Check whether reco neutrino is in FV
            if(!(sliceRecoVX < FVCut_xHigh && sliceRecoVX > FVCut_xLow  && std::abs(sliceRecoVX) > FVCut_xCentre && sliceRecoVY < FVCut_yHigh && sliceRecoVY > FVCut_yLow && sliceRecoVZ > FVCut_zLow && sliceRecoVZ < FVCut_zHigh)){
                // Reco neutrino isn't in FV
                fillHistogram(&sliceRecoNeutInFVAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, 0.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, 0.5, &weights);
            } else{
                // Reco neutrino is in FV
                fillHistogram(&sliceRecoNeutInFVAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, 1.5, &weights);
                fillSplitIntHistogram(&sliceRecoNeutInFVAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, 1.5, &weights);
            }

        } // End of looping through slices

    } // End of looping through events

    int drawLine = 1;
    int left = 0;
    int right = 1;

    // Draw histograms here
    styleDrawAll(sliceCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_beforeCuts.pdf").c_str(), ("sliceCompleteness_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCompletenessBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCompleteness_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts.pdf").c_str(), ("sliceCompleteness_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCompletenessAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(slicePurityBeforeCuts, 999, 999, 999, 999, (base_path + "slicePurity_beforeCuts.pdf").c_str(), ("slicePurity_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(slicePurityBeforeCuts, 999, 999, 999, 999, (base_path + "slicePurity_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(slicePurityBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "slicePurity_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(slicePurityAfterCuts, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts.pdf").c_str(), ("slicePurity_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(slicePurityAfterCuts, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(slicePurityAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceCRUMBSBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_beforeCuts.pdf").c_str(), ("sliceCRUMBS_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCRUMBSBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "beforeCRUMBS_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts.pdf").c_str(), ("sliceCRUMBS_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCRUMBSAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceNumRecoNeutBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_beforeCuts.pdf").c_str(), ("sliceNumRecoNeut_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumRecoNeutBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumRecoNeutBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumRecoNeutAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts.pdf").c_str(), ("sliceNumRecoNeut_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumRecoNeutAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumRecoNeutAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceNumPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_beforeCuts.pdf").c_str(), ("sliceNumPFPs_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPFPsBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts.pdf").c_str(), ("sliceNumPFPs_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPFPsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
 
    styleDrawAll(sliceNumPrimaryPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_beforeCuts.pdf").c_str(), ("sliceNumPrimaryPFPs_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumPrimaryPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts.pdf").c_str(), ("sliceNumPrimaryPFPs_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
 
    styleDrawAll(sliceNumPrimaryPFPsMinHitBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_beforeCuts.pdf").c_str(), ("sliceNumPrimaryPFPsMinHit_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsMinHitBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsMinHitBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumPrimaryPFPsMinHitAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterCuts.pdf").c_str(), ("sliceNumPrimaryPFPsMinHit_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsMinHitAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsMinHitAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(ERecoSumThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_beforeCuts.pdf").c_str(), ("ERecoSumThetaReco_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoSumThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoSumThetaRecoBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(ERecoSumThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts.pdf").c_str(), ("ERecoSumThetaReco_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoSumThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoSumThetaRecoAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(ERecoHighestThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_beforeCuts.pdf").c_str(), ("ERecoHighestThetaReco_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(ERecoHighestThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts.pdf").c_str(), ("ERecoHighestThetaReco_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_pfp10cmPoints_beforeCuts.pdf").c_str(), ("ERecoHighestThetaReco_pfp10cmPoints_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_pfp10cmPoints_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoBeforeCuts_splitDLNuE_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_pfp10cmPoints_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(ERecoHighestThetaRecoAfterCuts_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_pfp10cmPoints_afterCuts.pdf").c_str(), ("ERecoHighestThetaReco_pfp10cmPoints_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterCuts_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_pfp10cmPoints_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterCuts_splitDLNuE_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_pfp10cmPoints_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(dEdxBeforeCuts, 999, 999, 999, 999, (base_path + "dEdx_beforeCuts.pdf").c_str(), ("dEdx_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(dEdxBeforeCuts, 999, 999, 999, 999, (base_path + "dEdx_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(dEdxBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "dEdx_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(dEdxAfterCuts, 999, 999, 999, 999, (base_path + "dEdx_afterCuts.pdf").c_str(), ("dEdx_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(dEdxAfterCuts, 999, 999, 999, 999, (base_path + "dEdx_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(dEdxAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "dEdx_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(razzledPDG11BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_beforeCuts.pdf").c_str(), ("razzledPDG11_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG11BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG11BeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG11_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG11AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts.pdf").c_str(), ("razzledPDG11_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG11AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG11AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(razzledPDG13BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_beforeCuts.pdf").c_str(), ("razzledPDG13_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG13BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG13BeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG13_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG13AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts.pdf").c_str(), ("razzledPDG13_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG13AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG13AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(razzledPDG22BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_beforeCuts.pdf").c_str(), ("razzledPDG22_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG22BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG22BeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG22_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG22AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts.pdf").c_str(), ("razzledPDG22_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG22AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG22AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(razzledPDG211BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_beforeCuts.pdf").c_str(), ("razzledPDG211_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG211BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG211BeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG211_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG211AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts.pdf").c_str(), ("razzledPDG211_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG211AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG211AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(razzledPDG2212BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_beforeCuts.pdf").c_str(), ("razzledPDG2212_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG2212BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG2212BeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG2212_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG2212AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts.pdf").c_str(), ("razzledPDG2212_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG2212AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG2212AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    
    styleDrawAll(pfpCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_beforeCuts.pdf").c_str(), ("pfpCompleteness_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpCompletenessBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpCompleteness_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(pfpCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts.pdf").c_str(), ("pfpCompleteness_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpCompletenessAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(pfpPurityBeforeCuts, 999, 999, 999, 999, (base_path + "pfpPurity_beforeCuts.pdf").c_str(), ("pfpPurity_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpPurityBeforeCuts, 999, 999, 999, 999, (base_path + "pfpPurity_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpPurityBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpPurity_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(pfpPurityAfterCuts, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts.pdf").c_str(), ("pfpPurity_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpPurityAfterCuts, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpPurityAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(angleDifferenceSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferenceSignal_beforeCuts.pdf").c_str(), ("angleDifferenceSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFPSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFPSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCAPFPSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFP5cmSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFP5cmSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCAPFP5cmSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFP10cmSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFP10cmSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCAPFP10cmSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFP15cmSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFP15cmSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCAPFP15cmSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASliceSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASliceSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCASliceSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASlice5cmSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASlice5cmSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCASlice5cmSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASlice10cmSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASlice10cmSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCASlice10cmSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASlice15cmSignalBeforeCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASlice15cmSignal_beforeCuts.pdf").c_str(), ("angleDifferencePCASlice15cmSignal_beforeCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);

    styleDrawAll(angleDifferenceSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferenceSignal_afterCuts.pdf").c_str(), ("angleDifferenceSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFPSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFPSignal_afterCuts.pdf").c_str(), ("angleDifferencePCAPFPSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFP5cmSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFP5cmSignal_afterCuts.pdf").c_str(), ("angleDifferencePCAPFP5cmSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFP10cmSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFP10cmSignal_afterCuts.pdf").c_str(), ("angleDifferencePCAPFP10cmSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCAPFP15cmSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCAPFP15cmSignal_afterCuts.pdf").c_str(), ("angleDifferencePCAPFP15cmSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASliceSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASliceSignal_afterCuts.pdf").c_str(), ("angleDifferencePCASliceSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASlice5cmSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASlice5cmSignal_afterCuts.pdf").c_str(), ("angleDifferencePCASlice5cmSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASlice10cmSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASlice10cmSignal_afterCuts.pdf").c_str(), ("angleDifferencePCASlice10cmSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    styleDrawAll(angleDifferencePCASlice15cmSignalAfterCuts, 999, 999, 999, 999, (base_path + "angleDifferencePCASlice15cmSignal_afterCuts.pdf").c_str(), ("angleDifferencePCASlice15cmSignal_afterCuts").c_str(), "topRight", nullptr, &right, true, false, false, false, false, false, true, false, false);
    
    styleDrawAll(sliceRecoVXBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVX_beforeCuts.pdf").c_str(), ("sliceRecoVX_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVX_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVX_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterCuts.pdf").c_str(), ("sliceRecoVX_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVYBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVY_beforeCuts.pdf").c_str(), ("sliceRecoVY_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVY_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVY_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterCuts.pdf").c_str(), ("sliceRecoVY_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVZBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZ_beforeCuts.pdf").c_str(), ("sliceRecoVZ_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZ_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZ_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterCuts.pdf").c_str(), ("sliceRecoVZ_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVXSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_beforeCuts.pdf").c_str(), ("sliceRecoVXSmallerBins_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXSmallerBinsBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterCuts.pdf").c_str(), ("sliceRecoVXSmallerBins_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXSmallerBinsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVYSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_beforeCuts.pdf").c_str(), ("sliceRecoVYSmallerBins_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYSmallerBinsBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterCuts.pdf").c_str(), ("sliceRecoVYSmallerBins_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYSmallerBinsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVZSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_beforeCuts.pdf").c_str(), ("sliceRecoVZSmallerBins_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZSmallerBinsBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterCuts.pdf").c_str(), ("sliceRecoVZSmallerBins_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZSmallerBinsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVXLowBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_beforeCuts.pdf").c_str(), ("sliceRecoVXLow_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXLowBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXLowBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXLowAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterCuts.pdf").c_str(), ("sliceRecoVXLow_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXLowAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXLowAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    
    styleDrawAll(sliceRecoVYLowBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_beforeCuts.pdf").c_str(), ("sliceRecoVYLow_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYLowBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYLowBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYLowAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterCuts.pdf").c_str(), ("sliceRecoVYLow_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYLowAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYLowAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVZLowBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_beforeCuts.pdf").c_str(), ("sliceRecoVZLow_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZLowBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZLowBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZLowAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterCuts.pdf").c_str(), ("sliceRecoVZLow_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZLowAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZLowAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
 
    styleDrawAll(sliceRecoVXHighBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_beforeCuts.pdf").c_str(), ("sliceRecoVXHigh_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXHighBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_beforeCuts_BackSig.pdf").c_str(), ("sliceRecoVXHigh_beforeCuts_BackSig").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXHighBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXHighAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterCuts.pdf").c_str(), ("sliceRecoVXHigh_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXHighAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXHighAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVYHighBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_beforeCuts.pdf").c_str(), ("sliceRecoVYHigh_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYHighBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYHighBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYHighAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterCuts.pdf").c_str(), ("sliceRecoVYHigh_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYHighAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYHighAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoVZHighBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_beforeCuts.pdf").c_str(), ("sliceRecoVZHigh_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZHighBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_beforeCuts_BackSig.pdf").c_str(), ("sliceRecoVZHigh_beforeCuts_BackSig").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZHighBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZHighAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterCuts.pdf").c_str(), ("sliceRecoVZHigh_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZHighAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZHighAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawAll(sliceRecoNeutInFVBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoNeutInFV_beforeCuts.pdf").c_str(), ("sliceRecoNeutInFV_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true, false, true);
    styleDrawBackSig(sliceRecoNeutInFVBeforeCuts, 999, 999, 999, 999, (base_path + "sliceRecoNeutInFV_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true, false, true);
    styleDrawSplit(sliceRecoNeutInFVBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoNeutInFV_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true, false, true);
    styleDrawAll(sliceRecoNeutInFVAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoNeutInFV_afterCuts.pdf").c_str(), ("sliceRecoNeutInFV_afterCuts").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true, false, true);
    styleDrawBackSig(sliceRecoNeutInFVAfterCuts, 999, 999, 999, 999, (base_path + "sliceRecoNeutInFV_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true, false, true);
    styleDrawSplit(sliceRecoNeutInFVAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoNeutInFV_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true, false, true);

    styleDrawAll(energyAsymmetryBeforeCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_beforeCuts.pdf").c_str(), ("energyAsymmetry_beforeCuts").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false, true);
    styleDrawBackSig(energyAsymmetryBeforeCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(energyAsymmetryBeforeCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "energyAsymmetry_beforeCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(energyAsymmetryAfterCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts.pdf").c_str(), ("energyAsymmetry_afterCuts").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false, true);
    styleDrawBackSig(energyAsymmetryAfterCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(energyAsymmetryAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);

    // Before and after individual cut plots
    styleDrawAll(sliceNumPFPsAfterClearCosmicCut, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterClearCosmicCut.pdf").c_str(), ("sliceNumPFPs_afterClearCosmicCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPFPsAfterClearCosmicCut, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterClearCosmicCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPFPsAfterClearCosmicCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterClearCosmicCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumPFPsAfterNumPFPCut, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterNumPFPCut.pdf").c_str(), ("sliceNumPFPs_afterNumPFPCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPFPsAfterNumPFPCut, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterNumPFPCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPFPsAfterNumPFPCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterNumPFPCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceNumPFPsBeforeCuts, &sliceNumPFPsAfterClearCosmicCut, 999, 999, 999, 999, (base_path + "sliceNumPFPsAfterClearCosmicCut_upperBound").c_str(), "topRight", 1, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceNumPFPsBeforeCuts, &sliceNumPFPsAfterClearCosmicCut, 999, 999, 999, 999, (base_path + "sliceNumPFPsAfterClearCosmicCut_lowerBound").c_str(), "topRight", 1, nullptr, &right, -1);

    styleDrawAll(sliceNumRecoNeutAfterNumPFPCut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterNumPFPCut.pdf").c_str(), ("sliceNumRecoNeut_afterNumPFPCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumRecoNeutAfterNumPFPCut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterNumPFPCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumRecoNeutAfterNumPFPCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterNumPFPCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumRecoNeutAfterNumNeutrinoCut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterNumNeutrinoCut.pdf").c_str(), ("sliceNumRecoNeut_afterNumNeutrinoCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumRecoNeutAfterNumNeutrinoCut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterNumNeutrinoCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumRecoNeutAfterNumNeutrinoCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterNumNeutrinoCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceNumRecoNeutBeforeCuts, &sliceNumRecoNeutAfterNumPFPCut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeutAfterNumPFPCut_upperBound").c_str(), "topRight", 1, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceNumRecoNeutBeforeCuts, &sliceNumRecoNeutAfterNumPFPCut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeutAfterNumPFPCut_lowerBound").c_str(), "topRight", 1, nullptr, &right, -1);

    styleDrawAll(sliceCRUMBSAfterNumNeutrinoCut, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterNumNeutrinoCut.pdf").c_str(), ("sliceCRUMBS_afterNumNeutrinoCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSAfterNumNeutrinoCut, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterNumNeutrinoCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCRUMBSAfterNumNeutrinoCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterNumNeutrinoCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceCRUMBSAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCRUMBSCut.pdf").c_str(), ("sliceCRUMBS_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCRUMBSAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceCRUMBSBeforeCuts, &sliceCRUMBSAfterNumNeutrinoCut, 999, 999, 999, 999, (base_path + "sliceCRUMBSAfterNumNeutrinoCut_upperBound").c_str(), "topRight", 0.76, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceCRUMBSBeforeCuts, &sliceCRUMBSAfterNumNeutrinoCut, 999, 999, 999, 999, (base_path + "sliceCRUMBSAfterNumNeutrinoCut_lowerBound").c_str(), "topRight", 0.2, nullptr, &right, -1);

    styleDrawAll(sliceRecoVXAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVX_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterFVCut.pdf").c_str(), ("sliceRecoVX_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVX_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVXBeforeCuts, &sliceRecoVXAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXAfterCRUMBSCut_upperBound").c_str(), "topRight", 192, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceRecoVXBeforeCuts, &sliceRecoVXAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXAfterCRUMBSCut_lowerBound").c_str(), "topRight", -192, nullptr, &right, -1);

    styleDrawAll(sliceRecoVYAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVY_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterFVCut.pdf").c_str(), ("sliceRecoVY_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVY_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVYBeforeCuts, &sliceRecoVYAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYAfterCRUMBSCut_upperBound").c_str(), "topRight", 194, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceRecoVYBeforeCuts, &sliceRecoVYAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYAfterCRUMBSCut_lowerBound").c_str(), "topRight", -194, nullptr, &right, -1);

    styleDrawAll(sliceRecoVZAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVZ_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterFVCut.pdf").c_str(), ("sliceRecoVZ_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZ_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVZBeforeCuts, &sliceRecoVZAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZAfterCRUMBSCut_upperBound").c_str(), "topRight", 450, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceRecoVZBeforeCuts, &sliceRecoVZAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZAfterCRUMBSCut_lowerBound").c_str(), "topRight", 6, nullptr, &right, -1);

    styleDrawAll(sliceRecoVXSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVXSmallerBins_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXSmallerBinsAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXSmallerBinsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterFVCut.pdf").c_str(), ("sliceRecoVXSmallerBins_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXSmallerBinsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXSmallerBinsAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBins_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVXSmallerBinsBeforeCuts, &sliceRecoVXSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBinsAfterCRUMBSCut_upperBound").c_str(), "topRight", 192, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceRecoVXSmallerBinsBeforeCuts, &sliceRecoVXSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXSmallerBinsAfterCRUMBSCut_lowerBound").c_str(), "topRight", -192, nullptr, &right, -1);

    styleDrawAll(sliceRecoVYSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVYSmallerBins_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYSmallerBinsAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYSmallerBinsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterFVCut.pdf").c_str(), ("sliceRecoVYSmallerBins_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYSmallerBinsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYSmallerBinsAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBins_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVYSmallerBinsBeforeCuts, &sliceRecoVYSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBinsAfterCRUMBSCut_upperBound").c_str(), "topRight", 194, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceRecoVYSmallerBinsBeforeCuts, &sliceRecoVYSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYSmallerBinsAfterCRUMBSCut_lowerBound").c_str(), "topRight", -194, nullptr, &right, -1);

    styleDrawAll(sliceRecoVZSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVZSmallerBins_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZSmallerBinsAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZSmallerBinsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterFVCut.pdf").c_str(), ("sliceRecoVZSmallerBins_afterFVCut").s_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZSmallerBinsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZSmallerBinsAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBins_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVZSmallerBinsBeforeCuts, &sliceRecoVZSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBinsAfterCRUMBSCut_upperBound").c_str(), "topRight", 450, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceRecoVZSmallerBinsBeforeCuts, &sliceRecoVZSmallerBinsAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZSmallerBinsAfterCRUMBSCut_lowerBound").c_str(), "topRight", 6, nullptr, &right, -1);

    styleDrawAll(sliceRecoVXLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVXLow_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXLowAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXLowAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterFVCut.pdf").c_str(), ("sliceRecoVXLow_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXLowAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXLowAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXLow_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVXLowBeforeCuts, &sliceRecoVXLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXLowAfterCRUMBSCut_lowerBound").c_str(), "topRight", -192, nullptr, &right, -1);

    styleDrawAll(sliceRecoVYLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVYLow_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYLowAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYLowAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterFVCut.pdf").c_str(), ("sliceRecoVYLow_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYLowAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYLowAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYLow_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVYLowBeforeCuts, &sliceRecoVYLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYLowAfterCRUMBSCut_lowerBound").c_str(), "topRight", -194, nullptr, &right, -1);

    styleDrawAll(sliceRecoVZLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVZLow_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZLowAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZLowAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterFVCut.pdf").c_str(), ("sliceRecoVZLow_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZLowAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZLowAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZLow_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVZLowBeforeCuts, &sliceRecoVZLowAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZLowAfterCRUMBSCut_lowerBound").c_str(), "topRight", 6, nullptr, &right, -1);

    styleDrawAll(sliceRecoVXHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVXHigh_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXHighAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVXHighAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterFVCut.pdf").c_str(), ("sliceRecoVXHigh_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVXHighAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVXHighAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVXHigh_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVXHighBeforeCuts, &sliceRecoVXHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVXHighAfterCRUMBSCut_upperBound").c_str(), "topRight", 192, nullptr, &right, 1);

    styleDrawAll(sliceRecoVYHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVYHigh_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYHighAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVYHighAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterFVCut.pdf").c_str(), ("sliceRecoVYHigh_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVYHighAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVYHighAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVYHigh_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVYHighBeforeCuts, &sliceRecoVYHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVYHighAfterCRUMBSCut_upperBound").c_str(), "topRight", 194, nullptr, &right, 1);

    styleDrawAll(sliceRecoVZHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterCRUMBSCut.pdf").c_str(), ("sliceRecoVZHigh_afterCRUMBSCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterCRUMBSCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZHighAfterCRUMBSCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterCRUMBSCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceRecoVZHighAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterFVCut.pdf").c_str(), ("sliceRecoVZHigh_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceRecoVZHighAfterFVCut, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceRecoVZHighAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceRecoVZHigh_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceRecoVZHighBeforeCuts, &sliceRecoVZHighAfterCRUMBSCut, 999, 999, 999, 999, (base_path + "sliceRecoVZHighAfterCRUMBSCut_upperBound").c_str(), "topRight", 450, nullptr, &right, 1);

    styleDrawAll(sliceNumPrimaryPFPsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterFVCut.pdf").c_str(), ("sliceNumPrimaryPFPs_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumPrimaryPFPsAfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterPrimaryPFPCut.pdf").c_str(), ("sliceNumPrimaryPFPs_afterPrimaryPFPCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsAfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterPrimaryPFPCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsAfterPrimaryPFPCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterPrimaryPFPCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceNumPrimaryPFPsBeforeCuts, &sliceNumPrimaryPFPsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsAfterFVCut_upperBound").c_str(), "topRight", 1, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceNumPrimaryPFPsBeforeCuts, &sliceNumPrimaryPFPsAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsAfterFVCut_lowerBound").c_str(), "topRight", 1, nullptr, &right, -1);

    styleDrawAll(sliceNumPrimaryPFPsMinHitAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterFVCut.pdf").c_str(), ("sliceNumPrimaryPFPsMinHit_afterFVCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsMinHitAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterFVCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsMinHitAfterFVCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterFVCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterPrimaryPFPCut.pdf").c_str(), ("sliceNumPrimaryPFPsMinHit_afterPrimaryPFPCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterPrimaryPFPCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsMinHitAfterPrimaryPFPCut_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHit_afterPrimaryPFPCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &sliceNumPrimaryPFPsMinHitBeforeCuts, &sliceNumPrimaryPFPsMinHitAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHitAfterFVCut_upperBound").c_str(), "topRight", 1, nullptr, &right, 1);
    efficiency(actualSignalCount, &sliceNumPrimaryPFPsMinHitBeforeCuts, &sliceNumPrimaryPFPsMinHitAfterFVCut, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsMinHitAfterFVCut_lowerBound").c_str(), "topRight", 1, nullptr, &right, -1);

    styleDrawAll(razzledPDG11AfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "razzledPDG11_afterPrimaryPFPCut.pdf").c_str(), ("razzledPDG11_afterPrimaryPFPCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG11AfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "razzledPDG11_afterPrimaryPFPCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG11AfterPrimaryPFPCut_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG11_afterPrimaryPFPCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG11AfterRazzled11Cut, 999, 999, 999, 999, (base_path + "razzledPDG11_afterRazzled11Cut.pdf").c_str(), ("razzledPDG11_afterRazzled11Cut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG11AfterRazzled11Cut, 999, 999, 999, 999, (base_path + "razzledPDG11_afterRazzled11Cut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG11AfterRazzled11Cut_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG11_afterRazzled11Cut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &razzledPDG11BeforeCuts, &razzledPDG11AfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "razzledPDG11AfterPrimaryPFPCut_upperBound").c_str(), "topRight", 1, nullptr, &right, 1);
    efficiency(actualSignalCount, &razzledPDG11BeforeCuts, &razzledPDG11AfterPrimaryPFPCut, 999, 999, 999, 999, (base_path + "razzledPDG11AfterPrimaryPFPCut_lowerBound").c_str(), "topRight", 0.875, nullptr, &right, -1);

    styleDrawAll(razzledPDG211AfterRazzled11Cut, 999, 999, 999, 999, (base_path + "razzledPDG211_afterRazzled11Cut.pdf").c_str(), ("razzledPDG211_afterRazzled11Cut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG211AfterRazzled11Cut, 999, 999, 999, 999, (base_path + "razzledPDG211_afterRazzled11Cut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG211AfterRazzled11Cut_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG211_afterRazzled11Cut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(razzledPDG211AfterRazzled211Cut, 999, 999, 999, 999, (base_path + "razzledPDG211_afterRazzled211Cut.pdf").c_str(), ("razzledPDG211_afterRazzled211Cut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG211AfterRazzled211Cut, 999, 999, 999, 999, (base_path + "razzledPDG211_afterRazzled211Cut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG211AfterRazzled211Cut_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG211_afterRazzled211Cut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &razzledPDG211BeforeCuts, &razzledPDG211AfterRazzled11Cut, 999, 999, 999, 999, (base_path + "razzledPDG211AfterRazzled11Cut_upperBound").c_str(), "topRight", 0.0125, nullptr, &right, 1);
    efficiency(actualSignalCount, &razzledPDG211BeforeCuts, &razzledPDG211AfterRazzled11Cut, 999, 999, 999, 999, (base_path + "razzledPDG211AfterRazzled11Cut_lowerBound").c_str(), "topRight", 0, nullptr, &right, -1);

    styleDrawAll(ERecoHighestThetaRecoAfterRazzled211Cut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterRazzled211Cut.pdf").c_str(), ("ERecoHighestThetaReco_afterRazzled211Cut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterRazzled211Cut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterRazzled211Cut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterRazzled211Cut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(ERecoHighestThetaRecoAfterETheta2Cut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterETheta2Cut.pdf").c_str(), ("ERecoHighestThetaReco_afterETheta2Cut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterETheta2Cut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterETheta2Cut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterETheta2Cut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &ERecoHighestThetaRecoBeforeCuts, &ERecoHighestThetaRecoAfterRazzled211Cut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoAfterRazzled211Cut_upperBound").c_str(), "topRight", 3.066, nullptr, &right, 1);
    efficiency(actualSignalCount, &ERecoHighestThetaRecoBeforeCuts, &ERecoHighestThetaRecoAfterRazzled211Cut, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoAfterRazzled211Cut_lowerBound").c_str(), "topRight", 0, nullptr, &right, -1);

    styleDrawAll(ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterRazzled211Cut_pfp10cmPoints.pdf").c_str(), ("ERecoHighestThetaReco_afterRazzled211Cut_pfp10cmPoints").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterRazzled211Cut_BackSig_pfp10cmPoints.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterRazzled211Cut_splitDLNuE_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterRazzled211Cut_splitInt_pfp10cmPoints.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(ERecoHighestThetaRecoAfterETheta2Cut_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterETheta2Cut_pfp10cmPoints.pdf").c_str(), ("ERecoHighestThetaReco_afterETheta2Cut_pfp10cmPoints").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterETheta2Cut_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterETheta2Cut_BackSig_pfp10cmPoints.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterETheta2Cut_splitDLNuE_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterETheta2Cut_splitInt_pfp10cmPoints.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints, &ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoAfterRazzled211Cut_upperBound_pfp10cmPoints").c_str(), "topRight", 3.066, nullptr, &right, 1);
    efficiency(actualSignalCount, &ERecoHighestThetaRecoBeforeCuts_pfp10cmPoints, &ERecoHighestThetaRecoAfterRazzled211Cut_pfp10cmPoints, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoAfterRazzled211Cut_lowerBound_pfp10cmPoints").c_str(), "topRight", 0, nullptr, &right, -1);

    styleDrawAll(dEdxAfterETheta2Cut, 999, 999, 999, 999, (base_path + "dEdx_afterETheta2Cut.pdf").c_str(), ("dEdx_afterETheta2Cut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(dEdxAfterETheta2Cut, 999, 999, 999, 999, (base_path + "dEdx_afterETheta2Cut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(dEdxAfterETheta2Cut_splitDLNuE, 999, 999, 999, 999, (base_path + "dEdx_afterETheta2Cut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawAll(dEdxAfterdEdxCut, 999, 999, 999, 999, (base_path + "dEdx_afterdEdxCut.pdf").c_str(), ("dEdx_afterdEdxCut").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(dEdxAfterdEdxCut, 999, 999, 999, 999, (base_path + "dEdx_afterdEdxCut_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(dEdxAfterdEdxCut_splitDLNuE, 999, 999, 999, 999, (base_path + "dEdx_afterdEdxCut_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(actualSignalCount, &dEdxBeforeCuts, &dEdxAfterETheta2Cut, 999, 999, 999, 999, (base_path + "dEdxAfterETheta2Cut_upperBound").c_str(), "topRight", 3.25, nullptr, &right, 1);
    efficiency(actualSignalCount, &dEdxBeforeCuts, &dEdxAfterETheta2Cut, 999, 999, 999, 999, (base_path + "dEdxAfterETheta2Cut_lowerBound").c_str(), "topRight", 0, nullptr, &right, -1);

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

    if(outRootFile){
        outRootFile->Write();
        outRootFile->Close();
    }
}
