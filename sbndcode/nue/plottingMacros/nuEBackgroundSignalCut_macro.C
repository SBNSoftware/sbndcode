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
#include <fstream>

typedef struct{
    TCanvas * canvas;
    TH1F* baseHist;
    TH1F* current;
    TH1F* uboone;
    TH1F* nuE;
} purHist_struct;

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

/*
void styleDrawPur(purHist_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr){
    
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    hists.current->SetLineWidth(2);
    hists.current->SetLineColor(TColor::GetColor("#e42536"));
    hists.uboone->SetLineWidth(2);
    hists.uboone->SetLineColor(TColor::GetColor("#5790fc"));
    hists.nuE->SetLineWidth(2);
    hists.nuE->SetLineColor(TColor::GetColor("#f89c20"));

    if((ymin != 999) && (ymax != 999)) hists.current->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) hists.current->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = 1;
    maxYValue = std::max({
        hists.current->GetMaximum(), hists.uboone->GetMaximum(), hists.nuE->GetMaximum()
    });

    std::cout << "Max values: hists.current = " << hists.current->GetMaximum() << ", hists.uboone = " << hists.uboone->GetMaximum() << ", hists.nuE = " << hists.nuE->GetMaximum() << std::endl;

    std::cout << "maxYValue = " << maxYValue << std::endl;

    double yminVal = 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = (maxYValue * 1.1);
        hists.current->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
        std::cout << "Set Y axis range to: " << yminVal << ", " << ymaxVal << std::endl;
    }

    hists.current->Draw("hist");
    hists.uboone->Draw("histsame");
    hists.nuE->Draw("histsame");

    hists.current->SetStats(0);
    hists.current->GetXaxis()->SetTickLength(0.04);
    hists.current->GetYaxis()->SetTickLength(0.03);
    hists.current->GetXaxis()->SetTickSize(0.02);
    hists.current->GetYaxis()->SetTickSize(0.02);

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    if(legendLocation == "topRight"){ Lxmin=0.67; Lymax=0.863; Lxmax=0.87; Lymin=0.640; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.33; Lymin=0.640; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.67; Lymax=0.36; Lxmax=0.87; Lymin=0.137; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymax=0.36; Lxmax=0.33; Lymin=0.137; }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(hists.current, "Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.uboone, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuE, "Pandora Deep Learning: SBND Nu+E Tune", "f");

    legend->SetTextSize(0.01);
    legend->SetMargin(0.12);
    legend->Draw();

    hists.canvas->SaveAs(filename);

}
*/

void styleDrawPur(purHist_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr,
                  bool writeMaxValues = false)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    hists.current->SetLineWidth(2);
    hists.current->SetLineColor(TColor::GetColor("#e42536"));
    hists.uboone->SetLineWidth(2);
    hists.uboone->SetLineColor(TColor::GetColor("#5790fc"));
    hists.nuE->SetLineWidth(2);
    hists.nuE->SetLineColor(TColor::GetColor("#f89c20"));

    if((ymin != 999) && (ymax != 999))
        hists.current->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999))
        hists.current->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({
        hists.current->GetMaximum(),
        hists.uboone->GetMaximum(),
        hists.nuE->GetMaximum()
    });

    std::cout << "Max values: hists.current = " << hists.current->GetMaximum()
              << ", hists.uboone = " << hists.uboone->GetMaximum()
              << ", hists.nuE = " << hists.nuE->GetMaximum() << std::endl;

    std::cout << "maxYValue = " << maxYValue << std::endl;

    double yminVal = 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = (maxYValue * 1.1);
        hists.current->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
        std::cout << "Set Y axis range to: " << yminVal << ", " << ymaxVal << std::endl;
    }

    hists.current->Draw("hist");
    hists.uboone->Draw("histsame");
    hists.nuE->Draw("histsame");

    hists.current->SetStats(0);
    hists.current->GetXaxis()->SetTickLength(0.04);
    hists.current->GetYaxis()->SetTickLength(0.03);
    hists.current->GetXaxis()->SetTickSize(0.02);
    hists.current->GetYaxis()->SetTickSize(0.02);

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    if(legendLocation == "topRight"){ Lxmin=0.69; Lymax=0.863; Lxmax=0.87; Lymin=0.74; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.31; Lymin=0.74; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.69; Lymax=0.26; Lxmax=0.87; Lymin=0.137; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymax=0.26; Lxmax=0.31; Lymin=0.137; }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(hists.current, "Pandora BDT SBND (without Refinement)", "f");
    legend->AddEntry(hists.uboone, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(hists.nuE, "Pandora Deep Learning: SBND Nu+E Tune", "f");

    legend->SetTextSize(0.01);
    legend->SetMargin(0.12);
    legend->Draw();

    hists.canvas->SaveAs(filename);

    if (writeMaxValues) {
        std::ofstream outfile("purity_max_values.txt", std::ios::app);
        if (outfile.is_open()) {
            outfile << "================" << std::endl;
            outfile << filename << std::endl;

            int binCurrMax = hists.current->GetMaximumBin();
            double currMaxVal = hists.current->GetBinContent(binCurrMax);
            double currBinX = hists.current->GetXaxis()->GetBinCenter(binCurrMax);
            outfile << "BDT: EffxPur Max Value = " << currMaxVal
                    << ", Bin Value = " << currBinX << std::endl;

            int binUbooneMax = hists.uboone->GetMaximumBin();
            double ubooneMaxVal = hists.uboone->GetBinContent(binUbooneMax);
            double ubooneBinX = hists.uboone->GetXaxis()->GetBinCenter(binUbooneMax);
            outfile << "DL Uboone Train: EffxPur Max Value = " << ubooneMaxVal
                    << ", Bin Value = " << ubooneBinX << std::endl;

            // nuE
            int binNuEMax = hists.nuE->GetMaximumBin();
            double nuEMaxVal = hists.nuE->GetBinContent(binNuEMax);
            double nuEBinX = hists.nuE->GetXaxis()->GetBinCenter(binNuEMax);
            outfile << "DL Nu+E Train: EffxPur Max Value = " << nuEMaxVal
                    << ", Bin Value = " << nuEBinX << std::endl;

            outfile << "================" << std::endl << std::endl;
            outfile.close();
        } else {
            std::cerr << "Error: could not open purity_max_values.txt for writing." << std::endl;
        }
    }
}

void styleDrawBackSig(histGroup_struct hists,
                      double ymin, double ymax, double xmin, double xmax,
                      const char* filename, const std::string& legendLocation,
                      bool includeCurrent = true, bool includeUboone = true, bool includeNuE = true,
                      bool useLogScale = false)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();
    hists.canvas->SetLogy(useLogScale);

    auto combine = [useLogScale](TH1F* a, TH1F* b, TH1F* c, const char* name) -> TH1F* {
        TH1F* combo = nullptr;
        if (a) combo = (TH1F*)a->Clone(name);
        else if (b) combo = (TH1F*)b->Clone(name);
        else if (c) combo = (TH1F*)c->Clone(name);
        else return nullptr;

        combo->Reset();
        if (a) combo->Add(a);
        if (b) combo->Add(b);
        if (c) combo->Add(c);

        if (useLogScale) {
            for (int i = 1; i <= combo->GetNbinsX(); ++i)
                if (combo->GetBinContent(i) <= 0)
                    combo->SetBinContent(i, 1e-1);
        }

        return combo;
    };

    TH1F* currentSignalCombo     = combine(hists.currentSignal, hists.currentSignalFuzzy, nullptr, "currentSignalCombo");
    TH1F* currentBackgroundCombo = combine(hists.currentBNB, hists.currentBNBFuzzy, hists.currentCosmic, "currentBackgroundCombo");

    TH1F* ubooneSignalCombo     = combine(hists.ubooneSignal, hists.ubooneSignalFuzzy, nullptr, "ubooneSignalCombo");
    TH1F* ubooneBackgroundCombo = combine(hists.ubooneBNB, hists.ubooneBNBFuzzy, hists.ubooneCosmic, "ubooneBackgroundCombo");

    TH1F* nuESignalCombo     = combine(hists.nuESignal, hists.nuESignalFuzzy, nullptr, "nuESignalCombo");
    TH1F* nuEBackgroundCombo = combine(hists.nuEBNB, hists.nuEBNBFuzzy, hists.nuECosmic, "nuEBackgroundCombo");

    std::vector<TH1F*> allHists = {
        currentSignalCombo, currentBackgroundCombo,
        ubooneSignalCombo, ubooneBackgroundCombo,
        nuESignalCombo, nuEBackgroundCombo
    };

    double maxYValue = 0.0;
    for (auto* hist : allHists)
        if (hist && hist->GetMaximum() > maxYValue)
            maxYValue = hist->GetMaximum();

    double yminVal = useLogScale ? 1e-6 : 0;
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


void styleDrawAll(histGroup_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr,
                  bool includeSignal = true, bool includeSignalFuzzy = true,
                  bool includeBNB = true, bool includeBNBFuzzy = true,
                  bool includeCosmic = true,
                  bool includeDLUboone = true, bool includeDLNuE = true,
                  bool includeBDT = true,
                  bool useLogScale = false)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    if (useLogScale)
        hists.canvas->SetLogy(1);
    else
        hists.canvas->SetLogy(0);

    std::vector<TH1F*> allHists = {
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

    std::cout << "maxYValue = " << maxYValue << std::endl;
    double yminVal = useLogScale ? 1e-1 : 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = useLogScale ? (maxYValue * 100.0) : (maxYValue * 1.1);
        std::cout << "setting yaxis to " << yminVal << ", " << ymaxVal << std::endl;
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

void efficiency(histGroup_struct hists, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr, double efficiencyWay = 0.0){
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();
  
    purHist_struct purHists;
    purHists.canvas = hists.canvas;
    purHists.baseHist = hists.baseHist;

    purHist_struct effPurHists;
    effPurHists.canvas = hists.canvas;
    effPurHists.baseHist = hists.baseHist;

    purHists.current = (TH1F*) hists.currentCosmic->Clone("pur_currentCosmic");
    purHists.current->Reset();
    purHists.current->GetYaxis()->SetTitle("Purity");
    purHists.current->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle()); 

    purHists.uboone = (TH1F*) hists.ubooneCosmic->Clone("pur_ubooneCosmic");
    purHists.uboone->Reset();
    purHists.uboone->GetYaxis()->SetTitle("Purity");
    purHists.uboone->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle()); 
    
    purHists.nuE = (TH1F*) hists.nuECosmic->Clone("pur_nuECosmic");
    purHists.nuE->Reset();
    purHists.nuE->GetYaxis()->SetTitle("Purity");
    purHists.nuE->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle()); 
    
    effPurHists.current = (TH1F*) hists.currentCosmic->Clone("effPur_currentCosmic");
    effPurHists.current->Reset();
    effPurHists.current->GetYaxis()->SetTitle("Efficiency x Purity");
    effPurHists.current->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle()); 
    
    effPurHists.uboone = (TH1F*) hists.ubooneCosmic->Clone("effPur_ubooneCosmic");
    effPurHists.uboone->Reset();
    effPurHists.uboone->GetYaxis()->SetTitle("Efficiency x Purity");
    effPurHists.uboone->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle()); 
    
    effPurHists.nuE = (TH1F*) hists.nuECosmic->Clone("effPur_nuECosmic");
    effPurHists.nuE->Reset();
    effPurHists.nuE->GetYaxis()->SetTitle("Efficiency x Purity");
    effPurHists.nuE->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle()); 

    histGroup_struct effHists;
    effHists.canvas = hists.canvas;
    effHists.baseHist = hists.baseHist;

    effHists.currentCosmic = (TH1F*) hists.currentCosmic->Clone("eff_currentCosmic");
    effHists.currentCosmic->Reset();
    effHists.currentCosmic->GetYaxis()->SetTitle("Rejection");
    effHists.currentCosmic->GetXaxis()->SetTitle(hists.currentCosmic->GetXaxis()->GetTitle());
    
    effHists.currentSignal = (TH1F*) hists.currentSignal->Clone("eff_currentSignal");
    effHists.currentSignal->Reset();
    effHists.currentSignal->GetYaxis()->SetTitle("Efficiency");
    effHists.currentSignal->GetXaxis()->SetTitle(hists.currentSignal->GetXaxis()->GetTitle());

    effHists.currentSignalFuzzy = (TH1F*) hists.currentSignalFuzzy->Clone("eff_currentSignalFuzzy");
    effHists.currentSignalFuzzy->Reset();
    effHists.currentSignalFuzzy->GetYaxis()->SetTitle("Efficiency");
    effHists.currentSignalFuzzy->GetXaxis()->SetTitle(hists.currentSignal->GetXaxis()->GetTitle());
    
    effHists.currentBNB = (TH1F*) hists.currentBNB->Clone("eff_currentBNB");
    effHists.currentBNB->Reset();
    effHists.currentBNB->GetYaxis()->SetTitle("Rejection");
    effHists.currentBNB->GetXaxis()->SetTitle(hists.currentBNB->GetXaxis()->GetTitle());

    effHists.currentBNBFuzzy = (TH1F*) hists.currentBNBFuzzy->Clone("eff_currentBNBFuzzy");
    effHists.currentBNBFuzzy->Reset();
    effHists.currentBNBFuzzy->GetYaxis()->SetTitle("Rejection");
    effHists.currentBNBFuzzy->GetXaxis()->SetTitle(hists.currentBNB->GetXaxis()->GetTitle());

    effHists.ubooneCosmic = (TH1F*) hists.ubooneCosmic->Clone("eff_ubooneCosmic");
    effHists.ubooneCosmic->Reset();
    effHists.ubooneCosmic->GetYaxis()->SetTitle("Rejection");
    effHists.ubooneCosmic->GetXaxis()->SetTitle(hists.ubooneCosmic->GetXaxis()->GetTitle());
    
    effHists.ubooneSignal = (TH1F*) hists.ubooneSignal->Clone("eff_ubooneSignal");
    effHists.ubooneSignal->Reset();
    effHists.ubooneSignal->GetYaxis()->SetTitle("Efficiency");
    effHists.ubooneSignal->GetXaxis()->SetTitle(hists.ubooneSignal->GetXaxis()->GetTitle());

    effHists.ubooneSignalFuzzy = (TH1F*) hists.ubooneSignalFuzzy->Clone("eff_ubooneSignalFuzzy");
    effHists.ubooneSignalFuzzy->Reset();
    effHists.ubooneSignalFuzzy->GetYaxis()->SetTitle("Efficiency");
    effHists.ubooneSignalFuzzy->GetXaxis()->SetTitle(hists.ubooneSignal->GetXaxis()->GetTitle());
    
    effHists.ubooneBNB = (TH1F*) hists.ubooneBNB->Clone("eff_ubooneBNB");
    effHists.ubooneBNB->Reset();
    effHists.ubooneBNB->GetYaxis()->SetTitle("Rejection");
    effHists.ubooneBNB->GetXaxis()->SetTitle(hists.ubooneBNB->GetXaxis()->GetTitle());

    effHists.ubooneBNBFuzzy = (TH1F*) hists.ubooneBNBFuzzy->Clone("eff_ubooneBNBFuzzy");
    effHists.ubooneBNBFuzzy->Reset();
    effHists.ubooneBNBFuzzy->GetYaxis()->SetTitle("Rejection");
    effHists.ubooneBNBFuzzy->GetXaxis()->SetTitle(hists.ubooneBNB->GetXaxis()->GetTitle());

    effHists.nuECosmic = (TH1F*) hists.nuECosmic->Clone("eff_nuECosmic");
    effHists.nuECosmic->Reset();
    effHists.nuECosmic->GetYaxis()->SetTitle("Rejection");
    effHists.nuECosmic->GetXaxis()->SetTitle(hists.nuECosmic->GetXaxis()->GetTitle());
    
    effHists.nuESignal = (TH1F*) hists.nuESignal->Clone("eff_nuESignal");
    effHists.nuESignal->Reset();
    effHists.nuESignal->GetYaxis()->SetTitle("Efficiency");
    effHists.nuESignal->GetXaxis()->SetTitle(hists.nuESignal->GetXaxis()->GetTitle());

    effHists.nuESignalFuzzy = (TH1F*) hists.nuESignalFuzzy->Clone("eff_nuESignalFuzzy");
    effHists.nuESignalFuzzy->Reset();
    effHists.nuESignalFuzzy->GetYaxis()->SetTitle("Efficiency");
    effHists.nuESignalFuzzy->GetXaxis()->SetTitle(hists.nuESignal->GetXaxis()->GetTitle());
    
    effHists.nuEBNB = (TH1F*) hists.nuEBNB->Clone("eff_nuEBNB");
    effHists.nuEBNB->Reset();
    effHists.nuEBNB->GetYaxis()->SetTitle("Rejection");
    effHists.nuEBNB->GetXaxis()->SetTitle(hists.nuEBNB->GetXaxis()->GetTitle());

    effHists.nuEBNBFuzzy = (TH1F*) hists.nuEBNBFuzzy->Clone("eff_nuEBNBFuzzy");
    effHists.nuEBNBFuzzy->Reset();
    effHists.nuEBNBFuzzy->GetYaxis()->SetTitle("Rejection");
    effHists.nuEBNBFuzzy->GetXaxis()->SetTitle(hists.nuEBNB->GetXaxis()->GetTitle());

    int numBins = hists.currentSignal->GetNbinsX();

    double currentSignalSum = 0.0;
    double currentSignalTotal = 0.0;
    double currentSignalFuzzySum = 0.0;
    double currentSignalFuzzyTotal = 0.0;
    double ubooneSignalSum = 0.0;
    double ubooneSignalTotal = 0.0;
    double ubooneSignalFuzzySum = 0.0;
    double ubooneSignalFuzzyTotal = 0.0;
    double nuESignalSum = 0.0;
    double nuESignalTotal = 0.0;
    double nuESignalFuzzySum = 0.0;
    double nuESignalFuzzyTotal = 0.0;

    double currentCosmicSum = 0.0;
    double currentCosmicTotal = 0.0;
    double ubooneCosmicSum = 0.0;
    double ubooneCosmicTotal = 0.0;
    double nuECosmicSum = 0.0;
    double nuECosmicTotal = 0.0;
    double currentBNBSum = 0.0;
    double currentBNBTotal = 0.0;
    double ubooneBNBSum = 0.0;
    double ubooneBNBTotal = 0.0;
    double nuEBNBSum = 0.0;
    double nuEBNBTotal = 0.0;
    double currentBNBFuzzySum = 0.0;
    double currentBNBFuzzyTotal = 0.0;
    double ubooneBNBFuzzySum = 0.0;
    double ubooneBNBFuzzyTotal = 0.0;
    double nuEBNBFuzzySum = 0.0;
    double nuEBNBFuzzyTotal = 0.0;

    // efficiencyWay == -1 includes everything to the right of the cut
    for(int i = 0; i <= numBins+1; ++i){
        currentSignalTotal += hists.currentSignal->GetBinContent(i);
        currentSignalFuzzyTotal += hists.currentSignalFuzzy->GetBinContent(i);
        ubooneSignalTotal += hists.ubooneSignal->GetBinContent(i);
        ubooneSignalFuzzyTotal += hists.ubooneSignalFuzzy->GetBinContent(i);
        nuESignalTotal += hists.nuESignal->GetBinContent(i);
        nuESignalFuzzyTotal += hists.nuESignalFuzzy->GetBinContent(i);
    
        currentCosmicTotal += hists.currentCosmic->GetBinContent(i);
        ubooneCosmicTotal += hists.ubooneCosmic->GetBinContent(i);
        nuECosmicTotal += hists.nuECosmic->GetBinContent(i);
        currentBNBTotal += hists.currentBNB->GetBinContent(i);
        ubooneBNBTotal += hists.ubooneBNB->GetBinContent(i);
        nuEBNBTotal += hists.nuEBNB->GetBinContent(i);
        currentBNBFuzzyTotal += hists.currentBNBFuzzy->GetBinContent(i);
        ubooneBNBFuzzyTotal += hists.ubooneBNBFuzzy->GetBinContent(i);
        nuEBNBFuzzyTotal += hists.nuEBNBFuzzy->GetBinContent(i);
    }

    printf("TOTALS:\n Signal: Current = %f, Uboone = %f, NuE = %f\nSignal Fuzzy: Current = %f, Uboone = %f, NuE = %f\n", currentSignalTotal, ubooneSignalTotal, nuESignalTotal, currentSignalFuzzyTotal, ubooneSignalFuzzyTotal, nuESignalFuzzyTotal);
    printf("BNB: Current = %f, Uboone = %f, NuE = %f\nBNB Fuzzy: Current = %f, Uboone = %f, NuE = %f\n", currentBNBTotal, ubooneBNBTotal, nuEBNBTotal, currentBNBFuzzyTotal, ubooneBNBFuzzyTotal, nuEBNBFuzzyTotal);
    printf("Cosmic: Current = %f, Uboone = %f, NuE = %f\n", currentCosmicTotal, ubooneCosmicTotal, nuECosmicTotal);
    
    for(int i = 0; i <= numBins+1; ++i){
        printf("--------------------------\n");
        if(i != 0 || i != numBins+1) std::cout << "Bin " << i << std::endl;
        if(i == 0 || i == numBins+1) std::cout << "Under/Overflow bin" << std::endl;
        currentSignalSum += hists.currentSignal->GetBinContent(i);
        currentSignalFuzzySum += hists.currentSignalFuzzy->GetBinContent(i);
        ubooneSignalSum += hists.ubooneSignal->GetBinContent(i);
        ubooneSignalFuzzySum += hists.ubooneSignalFuzzy->GetBinContent(i);
        nuESignalSum += hists.nuESignal->GetBinContent(i);
        nuESignalFuzzySum += hists.nuESignalFuzzy->GetBinContent(i);

        currentCosmicSum += hists.currentCosmic->GetBinContent(i);
        ubooneCosmicSum += hists.ubooneCosmic->GetBinContent(i);
        nuECosmicSum += hists.nuECosmic->GetBinContent(i);
        currentBNBSum += hists.currentBNB->GetBinContent(i);
        ubooneBNBSum += hists.ubooneBNB->GetBinContent(i);
        nuEBNBSum += hists.nuEBNB->GetBinContent(i);
        currentBNBFuzzySum += hists.currentBNBFuzzy->GetBinContent(i);
        ubooneBNBFuzzySum += hists.ubooneBNBFuzzy->GetBinContent(i);
        nuEBNBFuzzySum += hists.nuEBNBFuzzy->GetBinContent(i);

        printf("Cosmic Sums: Current = %f, Uboone = %f, Nu+E = %f\n", currentCosmicSum, ubooneCosmicSum, nuECosmicSum);
        printf("BNB Sums: Current = %f, Uboone = %f, Nu+E = %f\n", currentBNBSum, ubooneBNBSum, nuEBNBSum);
        printf("BNB Fuzzy Sums: Current = %f, Uboone = %f, Nu+E = %f\n", currentBNBFuzzySum, ubooneBNBFuzzySum, nuEBNBFuzzySum);
        printf("Signal Sums: Current = %f, Uboone = %f, Nu+E = %f\n", currentSignalSum, ubooneSignalSum, nuESignalSum);
        printf("Signal Fuzzy Sums: Current = %f, Uboone = %f, Nu+E = %f\n\n", currentSignalFuzzySum, ubooneSignalFuzzySum, nuESignalFuzzySum);

        double currentSignalEffVal = 0;
        double currentSignalFuzzyEffVal = 0;
        double ubooneSignalEffVal = 0;
        double ubooneSignalFuzzyEffVal = 0;
        double nuESignalEffVal = 0;
        double nuESignalFuzzyEffVal = 0;

        double currentAllSignalEffVal = 0.0;
        double ubooneAllSignalEffVal = 0.0;
        double nuEAllSignalEffVal = 0.0;

        double currentCosmicRejVal = 0;
        double ubooneCosmicRejVal = 0;
        double nuECosmicRejVal = 0;
        double currentBNBRejVal = 0;
        double ubooneBNBRejVal = 0;
        double nuEBNBRejVal = 0;
        double currentBNBFuzzyRejVal = 0;
        double ubooneBNBFuzzyRejVal = 0;
        double nuEBNBFuzzyRejVal = 0;

        double keptSignalCurrent = 0.0;
        double keptBackgroundCurrent = 0.0;
        double keptSignalUboone = 0.0;
        double keptBackgroundUboone = 0.0;
        double keptSignalNuE = 0.0;
        double keptBackgroundNuE = 0.0;

        double currentPurity = 0.0;
        double uboonePurity = 0.0;
        double nuEPurity = 0.0;
       
        double currentEffPur = 0.0;
        double ubooneEffPur = 0.0;
        double nuEEffPur = 0.0;
        
        if(efficiencyWay == -1){
            // efficiencyWay == -1 includes everything to the right of the cut
            currentSignalEffVal = (1 - (currentSignalSum/currentSignalTotal));
            currentSignalFuzzyEffVal = (1 - (currentSignalFuzzySum/currentSignalFuzzyTotal));
            ubooneSignalEffVal = (1 - (ubooneSignalSum/ubooneSignalTotal));
            ubooneSignalFuzzyEffVal = (1 - (ubooneSignalFuzzySum/ubooneSignalFuzzyTotal));
            nuESignalEffVal = (1 - (nuESignalSum/nuESignalTotal));
            nuESignalFuzzyEffVal = (1 - (nuESignalFuzzySum/nuESignalFuzzyTotal));
            
            currentCosmicRejVal = (currentCosmicSum/currentCosmicTotal);
            ubooneCosmicRejVal = (ubooneCosmicSum/ubooneCosmicTotal);
            nuECosmicRejVal = (nuECosmicSum/nuECosmicTotal);
            currentBNBRejVal = (currentBNBSum/currentBNBTotal);
            ubooneBNBRejVal = (ubooneBNBSum/ubooneBNBTotal);
            nuEBNBRejVal = (nuEBNBSum/nuEBNBTotal);
            currentBNBFuzzyRejVal = (currentBNBFuzzySum/currentBNBFuzzyTotal);
            ubooneBNBFuzzyRejVal = (ubooneBNBFuzzySum/ubooneBNBFuzzyTotal);
            nuEBNBFuzzyRejVal = (nuEBNBFuzzySum/nuEBNBFuzzyTotal);

            keptSignalCurrent = (currentSignalTotal - currentSignalSum); 
            keptBackgroundCurrent = ((currentBNBTotal - currentBNBSum) + (currentBNBFuzzyTotal - currentBNBFuzzySum) + (currentCosmicTotal - currentCosmicSum) + (currentSignalFuzzyTotal - currentSignalFuzzySum));
            keptSignalUboone = (ubooneSignalTotal - ubooneSignalSum); 
            keptBackgroundUboone = ((ubooneBNBTotal - ubooneBNBSum) + (ubooneBNBFuzzyTotal - ubooneBNBFuzzySum) + (ubooneCosmicTotal - ubooneCosmicSum) + (ubooneSignalFuzzyTotal - ubooneSignalFuzzySum));
            keptSignalNuE = (nuESignalTotal - nuESignalSum); 
            keptBackgroundNuE = ((nuEBNBTotal - nuEBNBSum) + (nuEBNBFuzzyTotal - nuEBNBFuzzySum) + (nuECosmicTotal - nuECosmicSum) + (nuESignalFuzzyTotal - nuESignalFuzzySum));

            currentAllSignalEffVal = (keptSignalCurrent / (currentSignalTotal));
            ubooneAllSignalEffVal = (keptSignalUboone / (ubooneSignalTotal));
            nuEAllSignalEffVal = (keptSignalNuE / (nuESignalTotal));

        } else if(efficiencyWay == 1){
            // efficiencyWay == 1 includes everything to the left of the cut
            currentSignalEffVal = (currentSignalSum/currentSignalTotal);
            currentSignalFuzzyEffVal = (currentSignalFuzzySum/currentSignalFuzzyTotal);
            ubooneSignalEffVal = (ubooneSignalSum/ubooneSignalTotal);
            ubooneSignalFuzzyEffVal = (ubooneSignalFuzzySum/ubooneSignalFuzzyTotal);
            nuESignalEffVal = (nuESignalSum/nuESignalTotal);
            nuESignalFuzzyEffVal = (nuESignalFuzzySum/nuESignalFuzzyTotal);
            
            currentCosmicRejVal = (1- (currentCosmicSum/currentCosmicTotal));
            ubooneCosmicRejVal = (1 - (ubooneCosmicSum/ubooneCosmicTotal));
            nuECosmicRejVal = (1 - (nuECosmicSum/nuECosmicTotal));
            currentBNBRejVal = (1 - (currentBNBSum/currentBNBTotal));
            ubooneBNBRejVal = (1 - (ubooneBNBSum/ubooneBNBTotal));
            nuEBNBRejVal = (1 - (nuEBNBSum/nuEBNBTotal));
            currentBNBFuzzyRejVal = (1 - (currentBNBFuzzySum/currentBNBFuzzyTotal));
            ubooneBNBFuzzyRejVal = (1 - (ubooneBNBFuzzySum/ubooneBNBFuzzyTotal));
            nuEBNBFuzzyRejVal = (1 - (nuEBNBFuzzySum/nuEBNBFuzzyTotal));
        
            keptSignalCurrent = (currentSignalSum);
            keptBackgroundCurrent = (currentCosmicSum + currentBNBSum + currentBNBFuzzySum + currentSignalFuzzySum);  
            keptSignalUboone = (ubooneSignalSum);
            keptBackgroundUboone = (ubooneCosmicSum + ubooneBNBSum + ubooneBNBFuzzySum + ubooneSignalFuzzySum);  
            keptSignalNuE = (nuESignalSum);
            keptBackgroundNuE = (nuECosmicSum + nuEBNBSum + nuEBNBFuzzySum + nuESignalFuzzySum);  
            
            currentAllSignalEffVal = (keptSignalCurrent / (currentSignalTotal));
            ubooneAllSignalEffVal = (keptSignalUboone / (ubooneSignalTotal));
            nuEAllSignalEffVal = (keptSignalNuE / (nuESignalTotal));
        }

        currentPurity = (keptSignalCurrent / (keptSignalCurrent + keptBackgroundCurrent));
        uboonePurity = (keptSignalUboone / (keptSignalUboone + keptBackgroundUboone));
        nuEPurity = (keptSignalNuE / (keptSignalNuE + keptBackgroundNuE));

        currentEffPur = (currentAllSignalEffVal * currentPurity);
        ubooneEffPur = (ubooneAllSignalEffVal * uboonePurity);
        nuEEffPur = (nuEAllSignalEffVal * nuEPurity);
        printf("All Signal Sums: Current = %f, Uboone = %f, Nu+E = %f\n\n", keptSignalCurrent, keptSignalUboone, keptSignalNuE);

        printf("Efficiency Values:\nSignal: Current = %f, Uboone = %f, Nu+E = %f\nFuzzy Signal: Current = %f, Uboone = %f, Nu+E = %f\n\n", currentSignalEffVal, ubooneSignalEffVal, nuESignalEffVal, currentSignalFuzzyEffVal, ubooneSignalFuzzyEffVal, nuESignalFuzzyEffVal);
        printf("Rejection Values:\nBNB: Current = %f, Uboone = %f, Nu+E = %f\nFuzzy BNB: Current = %f, Uboone = %f, Nu+E = %f\nCosmics: Current = %f, Uboone = %f, Nu+E = %f\n\n", currentBNBRejVal, ubooneBNBRejVal, nuEBNBRejVal, currentBNBFuzzyRejVal, ubooneBNBFuzzyRejVal, nuEBNBFuzzyRejVal, currentCosmicRejVal, ubooneCosmicRejVal, nuECosmicRejVal);
        printf("Purity Values: Current = %f, Uboone = %f, Nu+E = %f\n\n", currentPurity, uboonePurity, nuEPurity);
        printf("Eff x Purity Values: Current = %f, Uboone = %f, Nu+E = %f\n", currentEffPur, ubooneEffPur, nuEEffPur);
        printf("--------------------------\n");
        if(!std::isnan(currentSignalEffVal)) effHists.currentSignal->SetBinContent(i, currentSignalEffVal);
        if(!std::isnan(currentSignalFuzzyEffVal)) effHists.currentSignalFuzzy->SetBinContent(i, currentSignalFuzzyEffVal);
        if(!std::isnan(ubooneSignalEffVal)) effHists.ubooneSignal->SetBinContent(i, ubooneSignalEffVal);
        if(!std::isnan(ubooneSignalFuzzyEffVal)) effHists.ubooneSignalFuzzy->SetBinContent(i, ubooneSignalFuzzyEffVal);
        if(!std::isnan(nuESignalEffVal)) effHists.nuESignal->SetBinContent(i, nuESignalEffVal);
        if(!std::isnan(nuESignalFuzzyEffVal)) effHists.nuESignalFuzzy->SetBinContent(i, nuESignalFuzzyEffVal);
    
        if(!std::isnan(currentCosmicRejVal)) effHists.currentCosmic->SetBinContent(i, currentCosmicRejVal);
        if(!std::isnan(ubooneCosmicRejVal)) effHists.ubooneCosmic->SetBinContent(i, ubooneCosmicRejVal);
        if(!std::isnan(nuECosmicRejVal)) effHists.nuECosmic->SetBinContent(i, nuECosmicRejVal);
        if(!std::isnan(currentBNBRejVal)) effHists.currentBNB->SetBinContent(i, currentBNBRejVal);
        if(!std::isnan(ubooneBNBRejVal)) effHists.ubooneBNB->SetBinContent(i, ubooneBNBRejVal);
        if(!std::isnan(nuEBNBRejVal)) effHists.nuEBNB->SetBinContent(i, nuEBNBRejVal);
        if(!std::isnan(currentBNBFuzzyRejVal)) effHists.currentBNBFuzzy->SetBinContent(i, currentBNBFuzzyRejVal);
        if(!std::isnan(ubooneBNBFuzzyRejVal)) effHists.ubooneBNBFuzzy->SetBinContent(i, ubooneBNBFuzzyRejVal);
        if(!std::isnan(nuEBNBFuzzyRejVal)) effHists.nuEBNBFuzzy->SetBinContent(i, nuEBNBFuzzyRejVal);
    
        if(!std::isnan(currentPurity)) purHists.current->SetBinContent(i, currentPurity);
        if(!std::isnan(uboonePurity)) purHists.uboone->SetBinContent(i, uboonePurity);
        if(!std::isnan(nuEPurity)) purHists.nuE->SetBinContent(i, nuEPurity);
        
        if(!std::isnan(currentEffPur)){ effPurHists.current->SetBinContent(i, currentEffPur); std::cout << "Filled effPurHists.current with " << currentEffPur << std::endl;}
        if(!std::isnan(ubooneEffPur)){ effPurHists.uboone->SetBinContent(i, ubooneEffPur); std::cout << "Filled effPurHists.uboone with " << ubooneEffPur << std::endl;}
        if(!std::isnan(nuEEffPur)){ effPurHists.nuE->SetBinContent(i, nuEEffPur); std::cout << "Filled effPurHists.nuE with " << nuEEffPur << std::endl;}

    }

    purHists.current->SetMaximum(-1111);
    purHists.uboone->SetMaximum(-1111);
    purHists.nuE->SetMaximum(-1111);
    effPurHists.current->SetMaximum(-1111);
    effPurHists.uboone->SetMaximum(-1111);
    effPurHists.nuE->SetMaximum(-1111);

    std::string filenameEff = std::string(filename) + "_eff.pdf";
    std::string filenameRej = std::string(filename) + "_rej.pdf";
    std::string filenamePur = std::string(filename) + "_pur.pdf";
    std::string filenameEffPur = std::string(filename) + "_effPur.pdf";
    
    styleDrawAll(effHists, ymin, ymax, xmin, xmax, filenameEff.c_str(), legendLocation, drawLine, linePos, true, true, false, false, false);
    styleDrawAll(effHists, 0, 1, xmin, xmax, filenameRej.c_str(), legendLocation, drawLine, linePos, false, false, true, true, true);
    styleDrawPur(purHists, 999, 999, 999, 999, filenamePur.c_str(), legendLocation, drawLine, linePos);
    styleDrawPur(effPurHists, 999, 999, 999, 999, filenameEffPur.c_str(), legendLocation, drawLine, linePos, true);
}

void TwoDHistDraw(TH2D* hist, const char* filename, const char* title){
    TCanvas* TwoDHistCanvas = new TCanvas("2dHist_canvas", "2D Histogram", 200, 10, 800, 600);
    TwoDHistCanvas->SetTickx();
    TwoDHistCanvas->SetTicky();

    hist->SetTitle(title);
    hist->Draw("COLZ");
    hist->SetStats(0);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetYaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTickSize(0.02);
    hist->GetYaxis()->SetTickSize(0.02);

    TwoDHistCanvas->SaveAs(filename);

    TProfile* profX = hist->ProfileX("_pfx", 1, -1, "");

    TCanvas* ProfileCanvas = new TCanvas("profile_canvas", "TProfile from TH2D", 300, 50, 800, 600);
    ProfileCanvas->SetTickx();
    ProfileCanvas->SetTicky();

    profX->SetTitle(Form("ProfileX of %s", title));
    profX->SetLineColor(kBlack);
    profX->SetLineWidth(2);
    profX->SetMarkerStyle(20);
    profX->SetMarkerSize(0.8);
    profX->SetMarkerColor(kBlack);
    profX->Draw("E1");

    profX->GetXaxis()->SetTickLength(0.04);
    profX->GetYaxis()->SetTickLength(0.03);
    profX->GetXaxis()->SetTickSize(0.02);
    profX->GetYaxis()->SetTickSize(0.02);
    profX->SetStats(0);

    std::string profileFilename = std::string(filename);
    size_t dotPos = profileFilename.find_last_of(".");
    if (dotPos != std::string::npos)
        profileFilename.insert(dotPos, "_profile");
    else
        profileFilename += "_profile.df"; 

    ProfileCanvas->SaveAs(profileFilename.c_str());

    TwoDHistCanvas->Clear();
    ProfileCanvas->Clear();
}


void nuEBackgroundSignalCut_macro(){
    std::ofstream clearFile("purity_max_values.txt", std::ios::trunc);
    clearFile.close();

    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nue/srcs/sbndcode/sbndcode/nue/mergedAll.root");
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_IntimeBNBNuE_DLUbooneNuEBDT.root");
    std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsFVCuts/";

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

    double POTSignalNuE_notMissing = 0;
    double POTSignalUboone_notMissing = 0;
    double POTSignalBDT_notMissing = 0;
    
    double POTBNBNuE_notMissing = 0;
    double POTBNBUboone_notMissing = 0;
    double POTBNBBDT_notMissing = 0;
        
    double counteraaaaaa = 0;
    
    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 2) cosmicSpillsSumCurrent += subRunNumGenEvents; counteraaaaaa++;
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
            if(subRunDLCurrent == 2) POTSignalBDT_notMissing += subRunPOT;
            if(subRunDLCurrent == 0) POTSignalUboone_notMissing += subRunPOT;
            if(subRunDLCurrent == 5) POTSignalNuE_notMissing += subRunPOT;
                
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
            if(subRunDLCurrent == 2) POTBNBBDT_notMissing += subRunPOT;
            if(subRunDLCurrent == 0) POTBNBUboone_notMissing += subRunPOT;
            if(subRunDLCurrent == 5) POTBNBNuE_notMissing += subRunPOT;

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

    printf("Cosmic Spills: Current = %f, DL Uboone = %f, DL Nu+E = %f\n", cosmicSpillsSumCurrent, cosmicSpillsSumUboone, cosmicSpillsSumNuE);
    std::cout << "cosmicsPOTCurrent = " << cosmicsPOTCurrent << ", " << cosmicsPOTUboone << ", " << cosmicsPOTNuE << std::endl;

    weights_struct weights;
    weights.signalCurrent = POTSignalBDT_notMissing/POTSignalBDT_notMissing;
    weights.BNBCurrent = POTSignalBDT_notMissing/POTBNBBDT_notMissing;
    weights.cosmicsCurrent = POTSignalBDT_notMissing/cosmicsPOTCurrent;
    
    weights.signalUboone = POTSignalBDT_notMissing/POTSignalUboone_notMissing;
    weights.BNBUboone = POTSignalBDT_notMissing/POTBNBUboone_notMissing;
    weights.cosmicsUboone = POTSignalBDT_notMissing/cosmicsPOTUboone;
    
    weights.signalNuE = POTSignalBDT_notMissing/POTSignalNuE_notMissing;
    weights.BNBNuE = POTSignalBDT_notMissing/POTBNBNuE_notMissing;
    weights.cosmicsNuE = POTSignalBDT_notMissing/cosmicsPOTNuE;

    std::cout << "BNB POT Current = " << totalPOTBNBCurrent << ", Uboone = " << totalPOTBNBUboone << ", Nu+E = " << totalPOTBNBNuE << std::endl;
    std::cout << "ALL BNB POT Current = " << POTBNBBDT_notMissing << ", Uboone = " << POTBNBUboone_notMissing << ", Nu+E = " << POTBNBNuE_notMissing << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Signal POT Current = " << totalPOTSignalCurrent << ", Uboone = " << totalPOTSignalUboone << ", Nu+E = " << totalPOTSignalNuE << std::endl;
    std::cout << "ALL Signal POT Current = " << POTSignalBDT_notMissing << ", Uboone = " << POTSignalUboone_notMissing << ", Nu+E = " << POTSignalNuE_notMissing << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Spills Cosmics POT Current = " << cosmicsPOTCurrent << ", Uboone = " << cosmicsPOTUboone << ", Nu+E = " << cosmicsPOTNuE << std::endl;  
    printf("Weights:\nCurrent: Signal = %f, BNB = %f, Cosmics = %f\nUboone: Signal = %f, BNB = %f, Cosmics = %f\nNu+E: Signal = %f, BNB = %f, Cosmics = %f\n", weights.signalCurrent, weights.BNBCurrent, weights.cosmicsCurrent, weights.signalUboone, weights.BNBUboone, weights.cosmicsUboone, weights.signalNuE, weights.BNBNuE, weights.cosmicsNuE);

    printf("\nPOT after weighting:\nSignal: Current = %f, Uboone = %f, Nu+E = %f\n", (totalPOTSignalCurrent*weights.signalCurrent), (totalPOTSignalUboone*weights.signalUboone), (totalPOTSignalNuE*weights.signalNuE));

    if((totalPOTSignalCurrent * weights.signalCurrent == totalPOTSignalUboone * weights.signalUboone) &&
                (totalPOTSignalUboone * weights.signalUboone == totalPOTSignalNuE * weights.signalNuE)) std::cout << "POT IS THE SAME AFTER WEIGHTING" << std::endl;


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
    std::vector<double> *reco_sliceTrueVX = nullptr;  
    std::vector<double> *reco_sliceTrueVY = nullptr;  
    std::vector<double> *reco_sliceTrueVZ = nullptr;  

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
    
    std::vector<double> *reco_neutrinoID = nullptr;
    std::vector<double> *reco_neutrinoPDG = nullptr;
    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;
    
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
    
    tree->SetBranchAddress("reco_neutrinoID", &reco_neutrinoID);
    tree->SetBranchAddress("reco_neutrinoPDG", &reco_neutrinoPDG);
    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);

    Long64_t numEntries = tree->GetEntries();

    auto trueETheta2 = createHistGroup("trueETheta2", "E_{true}#theta_{true}^{2}", "E_{true}#theta_{true}^{2} (MeV rad^{2})", 32, 0, 4.088);

    auto sliceCompleteness = createHistGroup("sliceCompleteness", "Slice Completeness", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessDist = createHistGroup("sliceCompletenessDist", "Slice Completeness (Not Weighted)", "Completeness", 102, 0, 1.02);
    auto slicePurity = createHistGroup("slicePurity", "Slice Purity", "Purity", 102, 0, 1.02);
    auto slicePurityDist = createHistGroup("slicePurityDist", "Slice Purity (Not Weighted)", "Purity", 102, 0, 1.02);
    auto sliceCRUMBSScore = createHistGroup("sliceCRUMBSScore", "CRUMBS Score of the Slice", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSScoreDist = createHistGroup("sliceCRUMBSScoreDist", "CRUMBS Score of the Slice (Not Weighted)", "CRUMBS Score", 25, -1, 1);
    auto sliceNumPFPs = createHistGroup("sliceNumPFPs", "Number of PFPs in the Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsDist = createHistGroup("sliceNumPFPsDist", "Number of PFPs in the Slice (Not Weighted)", "Number of PFPs", 20, 0, 20); 

    auto ERecoSumThetaReco = createHistGroup("ERecoSumThetaReco", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoDist = createHistGroup("ERecoSumThetaRecoDist", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (Not Weighted)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaReco = createHistGroup("ERecoHighestThetaReco", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoDist = createHistGroup("ERecoHighestThetaRecoDist", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Not Weighted)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);

    auto ETrueThetaReco = createHistGroup("ETrueThetaReco", "E_{true}#theta_{reco}^{2}", "E_{true}#theta_{reco}^{2} (MeV rad^{2})", 40, 0, 20.44);
    auto ETrueThetaRecoDist = createHistGroup("ETrueThetaRecoDist", "E_{true}#theta_{reco}^{2} (Not Weighted)", "E_{true}#theta_{reco}^{2} (MeV rad^{2})", 40, 0, 20.44);
    auto ERecoSumThetaTrue = createHistGroup("ERecoSumThetaTrue", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice", "E_{reco}#theta_{true}^{2} (MeV rad^{2})", 24, 0, 3.066);
    auto ERecoSumThetaTrueDist = createHistGroup("ERecoSumThetaTrueDist", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (Not Weighted)", "E_{reco}#theta_{true}^{2} (MeV rad^{2})", 24, 0, 3.066);
    auto ERecoHighestThetaTrue = createHistGroup("ERecoHighestThetaTrue", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{true}^{2} (MeV rad^{2})", 24, 0, 3.066);
    auto ERecoHighestThetaTrueDist = createHistGroup("ERecoHighestThetaTrueDist", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice (Not Weighted)", "E_{reco}#theta_{true}^{2} (MeV rad^{2})", 24, 0, 3.066);
    auto ETrue = createHistGroup("ETrue", "True Energy of the Recoil Electron", "E_{true} (MeV)", 1000, 0, 10000);
    auto ETrueDist = createHistGroup("ETrueDist", "True Energy of the Recoil Electron (Not Weighted)", "E_{true} (MeV)", 1000, 0, 10000);
    auto ThetaTrue = createHistGroup("ThetaTrue", "True Angle of the Recoil Electron", "#theta_{true} (rad)", 20, 0, 0.4);
    auto ThetaTrueDist = createHistGroup("ThetaTrueDist", "True Angle of the Recoil Electron (Not Weighted)", "#theta_{true} (rad)", 20, 0, 0.4);

    auto deltaX = createHistGroup("deltaX", "#Deltax of Neutrino Vertex in Slice", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaXDist = createHistGroup("deltaXDist", "#Deltax of Neutrino Vertex in Slice (Not Weighted)", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaY = createHistGroup("deltaY", "#Deltay of Neutrino Vertex in Slice", "y_{Reco} - y_{True} (cm)", 80, -20, 20);
    auto deltaYDist = createHistGroup("deltaYDist", "#Deltay of Neutrino Vertex in Slice (Not Weighted)", "y_{Reco} - y_{True} (cm)", 80, -20, 20);
    auto deltaZ = createHistGroup("deltaZ", "#Deltaz of Neutrino Vertex in Slice", "z_{Reco} - z_{True} (cm)", 80, -20, 20);
    auto deltaZDist = createHistGroup("deltaZDist", "#Deltaz of Neutrino Vertex in Slice (Not Weighted)", "z_{Reco} - z_{True} (cm)", 80, -20, 20);
    auto deltaR = createHistGroup("deltaR", "#DeltaR of Neutrino Vertex in Slice", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 80, 0, 40);
    auto deltaRDist = createHistGroup("deltaRDist", "#DeltaR of Neutrino Vertex in Slice (Not Weighted)", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 80, 0, 40);
    auto deltaTheta = createHistGroup("deltaTheta", "Angle Between the True Electron and the Highest Energy PFP in the Slice", "#Delta#theta (degrees)", 90, 0, 180);
    auto deltaThetaDist = createHistGroup("deltaThetaDist", "Angle Between the True Electron and the Highest Energy PFP in the Slice (Not Weighted)", "#Delta#theta (degrees)", 90, 0, 180);

    auto pfpCompleteness = createHistGroup("pfpCompleteness", "Completeness of the Highest Energy PFP in the Slice", "Completeness", 50, 0, 1);
    auto pfpCompletenessDist = createHistGroup("pfpCompletenessDist", "Completeness of the Highest Energy PFP in the Slice (Not Weighted)", "Completeness", 50, 0, 1);
    auto pfpPurity = createHistGroup("pfpPurity", "Purity of the Highest Energy PFP in the Slice", "Purity", 50, 0, 1);
    auto pfpPurityDist = createHistGroup("pfpPurityDist", "Purity of the Highest Energy PFP in the Slice (Not Weighted)", "Purity", 50, 0, 1);

    auto recoX = createHistGroup("recoX", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 200, -202, 202);
    auto recoXDist = createHistGroup("recoXDist", "X Coordinate of Reco Neutrino (Not Weighted)", "x_{Reco} (cm)", 200, -202, 202);
    auto recoX_low = createHistGroup("recoX_low", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, -202, -182);
    auto recoXDist_low = createHistGroup("recoXDist_low", "X Coordinate of Reco Neutrino (Not Weighted)", "x_{Reco} (cm)", 40, -202, -182);
    auto recoX_high = createHistGroup("recoX_high", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoXDist_high = createHistGroup("recoXDist_high", "X Coordinate of Reco Neutrino (Not Weighted)", "x_{Reco} (cm)", 40, 182, 202);
    
    auto recoY = createHistGroup("recoY", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 200, -204, 204);
    auto recoYDist = createHistGroup("recoYDist", "Y Coordinate of Reco Neutrino (Not Weighted)", "y_{Reco} (cm)", 200, -204, 204);
    auto recoY_low = createHistGroup("recoY_low", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, -204, -184);
    auto recoYDist_low = createHistGroup("recoYDist_low", "Y Coordinate of Reco Neutrino (Not Weighted)", "y_{Reco} (cm)", 40, -204, -184);
    auto recoY_high = createHistGroup("recoY_high", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoYDist_high = createHistGroup("recoYDist_high", "Y Coordinate of Reco Neutrino (Not Weighted)", "y_{Reco} (cm)", 40, 184, 204);
    
    auto recoZ = createHistGroup("recoZ", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 250, 0, 510);
    auto recoZDist = createHistGroup("recoZDist", "Z Coordinate of Reco Neutrino (Not Weighted)", "z_{Reco} (cm)", 250, 0, 510);
    auto recoZ_low = createHistGroup("recoZ_low", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 40, 0, 20);
    auto recoZDist_low = createHistGroup("recoZDist_low", "Z Coordinate of Reco Neutrino (Not Weighted)", "z_{Reco} (cm)", 40, 0, 20);
    auto recoZ_high = createHistGroup("recoZ_high", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 40, 490, 510);
    auto recoZDist_high = createHistGroup("recoZDist_high", "Z Coordinate of Reco Neutrino (Not Weighted)", "z_{Reco} (cm)", 40, 490, 510);

    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0; double zMax = 509.4;

    TH2D *xCoordAngleDifferenceBDT_low = new TH2D("xCoordAngleDifferenceBDT_low", "", 10, xMin, (xMin + 20), 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDT_low = new TH2D("yCoordAngleDifferenceBDT_low", "", 10, yMin, (yMin + 20), 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDT_low = new TH2D("zCoordAngleDifferenceBDT_low", "", 10, zMin, (zMin + 20), 40, 0, 180);
    TH2D *xCoordAngleDifferenceBDT_high = new TH2D("xCoordAngleDifferenceBDT_high", "", 10, (xMax - 20), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDT_high = new TH2D("yCoordAngleDifferenceBDT_high", "", 10, (yMax - 20), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDT_high = new TH2D("zCoordAngleDifferenceBDT_high", "", 20, (zMax - 40), zMax, 60, 0, 180);
    
    TH2D *xCoordAngleDifferenceDLUboone_low = new TH2D("xCoordAngleDifferenceDLUboone_low", "", 10, xMin, (xMin + 20), 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLUboone_low = new TH2D("yCoordAngleDifferenceDLUboone_low", "", 10, yMin, (yMin + 20), 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLUboone_low = new TH2D("zCoordAngleDifferenceDLUboone_low", "", 10, zMin, (zMin + 20), 40, 0, 180);
    TH2D *xCoordAngleDifferenceDLUboone_high = new TH2D("xCoordAngleDifferenceDLUboone_high", "", 10, (xMax - 20), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLUboone_high = new TH2D("yCoordAngleDifferenceDLUboone_high", "", 10, (yMax - 20), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLUboone_high = new TH2D("zCoordAngleDifferenceDLUboone_high", "", 20, (zMax - 40), zMax, 60, 0, 180);
    
    TH2D *xCoordAngleDifferenceDLNuE_low = new TH2D("xCoordAngleDifferenceDLNuE_low", "", 10, xMin, (xMin + 20), 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLNuE_low = new TH2D("yCoordAngleDifferenceDLNuE_low", "", 10, yMin, (yMin + 20), 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLNuE_low = new TH2D("zCoordAngleDifferenceDLNuE_low", "", 10, zMin, (zMin + 20), 40, 0, 180);
    TH2D *xCoordAngleDifferenceDLNuE_high = new TH2D("xCoordAngleDifferenceDLNuE_high", "", 10, (xMax - 20), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLNuE_high = new TH2D("yCoordAngleDifferenceDLNuE_high", "", 10, (yMax - 20), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLNuE_high = new TH2D("zCoordAngleDifferenceDLNuE_high", "", 20, (zMax - 40), zMax, 60, 0, 180);
    
    TH2D *xCoordAngleDifferenceBDT = new TH2D("xCoordAngleDifferenceBDT", "", (int)round((xMax - xMin)/5), xMin, xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDT = new TH2D("yCoordAngleDifferenceBDT", "", (int)round((yMax - yMin)/5), yMin, yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDT = new TH2D("zCoordAngleDifferenceBDT", "", (int)round((zMax - zMin)/5), zMin, zMax, 40, 0, 180);
    
    TH2D *xCoordAngleDifferenceDLUboone = new TH2D("xCoordAngleDifferenceDLUboone", "", (int)round((xMax - xMin)/5), xMin, xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLUboone = new TH2D("yCoordAngleDifferenceDLUboone", "", (int)round((yMax - yMin)/5), yMin, yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLUboone = new TH2D("zCoordAngleDifferenceDLUboone", "", (int)round((zMax - zMin)/5), zMin, zMax, 40, 0, 180);
    
    TH2D *xCoordAngleDifferenceDLNuE = new TH2D("xCoordAngleDifferenceDLNuE", "", (int)round((xMax - xMin)/5), xMin, xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLNuE = new TH2D("yCoordAngleDifferenceDLNuE", "", (int)round((yMax - yMin)/5), yMin, yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLNuE = new TH2D("zCoordAngleDifferenceDLNuE", "", (int)round((zMax - zMin)/5), zMin, zMax, 40, 0, 180);

    double numEvents_BDTCosmic = 0;
    double numEvents_BDTBNB = 0;
    double numEvents_BDTNuE = 0;
    
    double numEvents_DLNuECosmic = 0;
    double numEvents_DLNuEBNB = 0;
    double numEvents_DLNuENuE = 0;

    for(Long64_t e = 0; e < numEntries; ++e){
        //printf("=============================================================================\n");
        tree->GetEntry(e);

        if(DLCurrent == 2 && signal == 3) numEvents_BDTCosmic++;
        else if(DLCurrent == 2 && signal == 2) numEvents_BDTBNB++;
        else if(DLCurrent == 2 && signal == 1) numEvents_BDTNuE++;
        else if(DLCurrent == 5 && signal == 3) numEvents_DLNuECosmic++;
        else if(DLCurrent == 5 && signal == 2) numEvents_DLNuEBNB++;
        else if(DLCurrent == 5 && signal == 1) numEvents_DLNuENuE++;

        double recoilElectron_energy = -999999;
        double recoilElectron_angle = -999999;
        double recoilElectron_DX = -999999;
        double recoilElectron_DY = -999999;
        double recoilElectron_DZ = -999999;

        // Looking at the true recoil electron
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->size() > 1) std::cout << "More than 1 true recoil electron in event!" << std::endl;
            if(truth_recoilElectronPDG->at(i) != -999999){
                // There is a true recoil electron in the event
                recoilElectron_energy = truth_recoilElectronEnergy->at(i);
                recoilElectron_angle = truth_recoilElectronAngle->at(i);
                recoilElectron_DX = truth_recoilElectronDX->at(i); 
                recoilElectron_DY = truth_recoilElectronDY->at(i); 
                recoilElectron_DZ = truth_recoilElectronDZ->at(i); 
                
                if(signal == 1 && DLCurrent == 2) trueETheta2.currentSignal->Fill(truth_recoilElectronETheta2->at(i), weights.signalCurrent); 
                else if(signal == 1 && DLCurrent == 0) trueETheta2.ubooneSignal->Fill(truth_recoilElectronETheta2->at(i), weights.signalUboone); 
                else if(signal == 1 && DLCurrent == 5) trueETheta2.nuESignal->Fill(truth_recoilElectronETheta2->at(i), weights.signalNuE); 
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
        if(reco_sliceID->size() == 0) continue;
        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) != -999999){
                // There is a reco slice in the event
                //printf("__________________________ NEW SLICE __________________________\n");
                //std::cout << "reco_sliceID->size() = " << reco_sliceID->size() << ", reco_sliceCompleteness->size() = " << reco_sliceCompleteness->size() << ", reco_slicePurity = " << reco_slicePurity->size() << ", reco_sliceScore = " << reco_sliceScore->size() << ", reco_sliceCategory = " << reco_sliceCategory->size() << ", reco_sliceInteraction = " << reco_sliceInteraction->size() << ", reco_sliceTrueVX = " << reco_sliceTrueVX->size() << ", reco_sliceTrueVY = " << reco_sliceTrueVY->size() << ", reco_sliceTrueVZ = " << reco_sliceTrueVZ->size() << std::endl;
                //printf("Slice ID = %f, Category = %f, Interaction = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", reco_sliceID->at(slice), reco_sliceCategory->at(slice), reco_sliceInteraction->at(slice), reco_sliceCompleteness->at(slice), reco_slicePurity->at(slice), reco_sliceScore->at(slice));
                //if(reco_sliceCategory->at(slice) != 0) printf("True Neutrino Vertex = (%f, %f, %f)\n", reco_sliceTrueVX->at(slice), reco_sliceTrueVY->at(slice), reco_sliceTrueVZ->at(slice));
        
                // Loop through all the reco neutrinos in the event          
                int PFPcounter = 0;
                int numPFPsSlice = 0;

                double summedEnergy = 0;
                double highestEnergy_PFPID = -999999;
                double highestEnergy_energy = -999999;
                double highestEnergy_theta = -999999;
                double highestEnergy_DX = -999999;
                double highestEnergy_DY = -999999;
                double highestEnergy_DZ = -999999;
                double highestEnergy_completeness = -999999;
                double highestEnergy_purity = -999999;

                double FVCut_xLow_BDT = -200;
                double FVCut_xHigh_BDT = 199;
                double FVCut_yLow_BDT = -200;
                double FVCut_yHigh_BDT = 182;
                double FVCut_zLow_BDT = 5;
                double FVCut_zHigh_BDT = 491;

                double FVCut_xLow_DLNuE = -199;
                double FVCut_xHigh_DLNuE = 199;
                double FVCut_yLow_DLNuE = -196;
                double FVCut_yHigh_DLNuE = -152;
                double FVCut_zLow_DLNuE = 5.1;
                double FVCut_zHigh_DLNuE = 491;

                for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                    PFPcounter++;
                    if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                        // This PFP is in the slice
                        numPFPsSlice++;
                        //printf("PFP %d: ID = %f, PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Energy = %f, Theta = %f, Track Score = %f, Completeness = %f, Purity = %f\n", PFPcounter, reco_particleID->at(pfp), reco_particlePDG->at(pfp), reco_particleIsPrimary->at(pfp), reco_particleVX->at(pfp), reco_particleVY->at(pfp), reco_particleVZ->at(pfp), reco_particleDX->at(pfp), reco_particleDY->at(pfp), reco_particleDZ->at(pfp), reco_particleBestPlaneEnergy->at(pfp), reco_particleTheta->at(pfp), reco_particleTrackScore->at(pfp), reco_particleCompleteness->at(pfp), reco_particlePurity->at(pfp));
                        
                        summedEnergy += reco_particleBestPlaneEnergy->at(pfp);
                        if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergy_energy){
                            highestEnergy_energy = reco_particleBestPlaneEnergy->at(pfp);
                            highestEnergy_theta = reco_particleTheta->at(pfp);
                            highestEnergy_PFPID = reco_particleID->at(pfp);
                            highestEnergy_DX = reco_particleDX->at(pfp);
                            highestEnergy_DY = reco_particleDY->at(pfp);
                            highestEnergy_DZ = reco_particleDZ->at(pfp);
                            highestEnergy_completeness = reco_particleCompleteness->at(pfp);
                            highestEnergy_purity = reco_particlePurity->at(pfp);
                        }
                    }
                } 

                //printf("PFP with highest energy (%f) has ID %f and Theta = %f\n", highestEnergy_energy, highestEnergy_PFPID, highestEnergy_theta);
                //printf("Summed energy of all PFPs in slice = %f\n", summedEnergy);
                //std::cout << "Number of PFPs in the slice = " << numPFPsSlice << std::endl;

                double angleDifference = -999999;
                if((highestEnergy_DX != -999999) && (recoilElectron_DX != -999999)){
                    double aDOTb = ((highestEnergy_DX * recoilElectron_DX) + (highestEnergy_DY * recoilElectron_DY) + (highestEnergy_DZ * recoilElectron_DZ));
                    double aMagnitude = std::sqrt((highestEnergy_DX * highestEnergy_DX) + (highestEnergy_DY * highestEnergy_DY) + (highestEnergy_DZ * highestEnergy_DZ));
                    double bMagnitude = std::sqrt((recoilElectron_DX * recoilElectron_DX) + (recoilElectron_DY * recoilElectron_DY) + (recoilElectron_DZ * recoilElectron_DZ));
                    double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                    angleDifference = (TMath::ACos(cosAngle) * TMath::RadToDeg());
                }

                double recoVX = -999999;
                double recoVY = -999999;
                double recoVZ = -999999;

                for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                    if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                        // Reco neutrino is in the slice
                        recoVX = reco_neutrinoVX->at(recoNeut);
                        recoVY = reco_neutrinoVY->at(recoNeut);
                        recoVZ = reco_neutrinoVZ->at(recoNeut);
                        //printf("Reco Neutrino in Slice: ID = %f, PDG = %f, Vertex = (%f, %f, %f)\n", reco_neutrinoID->at(recoNeut), reco_neutrinoPDG->at(recoNeut), reco_neutrinoVX->at(recoNeut), reco_neutrinoVY->at(recoNeut), reco_neutrinoVZ->at(recoNeut));
                    }
                }
                
                if(DLCurrent == 2){
                    if (!(recoVX < FVCut_xHigh_BDT && recoVX > FVCut_xLow_BDT && recoVY < FVCut_yHigh_BDT && recoVY > FVCut_yLow_BDT && recoVZ > FVCut_zLow_BDT && recoVZ < FVCut_zHigh_BDT)){ 
                        std::cout << "BDT: DOES NOT PASS CUTS WITH VX = " << recoVX << ", VY = " << recoVY << ", VZ = " << recoVZ << std::endl;
                        continue;
                    }
                } else if(DLCurrent == 5){
                    if (!(recoVX < FVCut_xHigh_DLNuE && recoVX > FVCut_xLow_DLNuE && recoVY < FVCut_yHigh_DLNuE && recoVY > FVCut_yLow_DLNuE && recoVZ > FVCut_zLow_DLNuE && recoVZ < FVCut_zHigh_DLNuE)){ 
                        std::cout << "DLNuE: DOES NOT PASS CUTS WITH VX = " << recoVX << ", VY = " << recoVY << ", VZ = " << recoVZ << std::endl;
                        continue;
                    }
                }

                // Filling Histograms
                if(reco_sliceCategory->at(slice) == 0){
                    // This is a cosmic slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentCosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.currentCosmic->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.currentCosmic->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.currentCosmic->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.currentCosmic->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.currentCosmic->Fill(reco_sliceScore->at(slice));

                        sliceNumPFPs.currentCosmic->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.currentCosmic->Fill(numPFPsSlice);

                        if(recoVX != -999999){
                            recoX.currentCosmic->Fill(recoVX, weight);
                            recoXDist.currentCosmic->Fill(recoVX);
                            recoY.currentCosmic->Fill(recoVY, weight);
                            recoYDist.currentCosmic->Fill(recoVY);
                            recoZ.currentCosmic->Fill(recoVZ, weight);
                            recoZDist.currentCosmic->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.currentCosmic->Fill(recoVX, weight);
                                recoXDist_low.currentCosmic->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.currentCosmic->Fill(recoVX, weight);
                                recoXDist_high.currentCosmic->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.currentCosmic->Fill(recoVY, weight);
                                recoYDist_low.currentCosmic->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.currentCosmic->Fill(recoVY, weight);
                                recoYDist_high.currentCosmic->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.currentCosmic->Fill(recoVZ, weight);
                                recoZDist_low.currentCosmic->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.currentCosmic->Fill(recoVZ, weight);
                                recoZDist_high.currentCosmic->Fill(recoVZ);
                            //}
                        }

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.currentCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.currentCosmic->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.currentCosmic->Fill(highestEnergy_completeness);
                            pfpPurity.currentCosmic->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.currentCosmic->Fill(highestEnergy_purity);
                        }

                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneCosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.ubooneCosmic->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.ubooneCosmic->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.ubooneCosmic->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.ubooneCosmic->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.ubooneCosmic->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.ubooneCosmic->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.ubooneCosmic->Fill(numPFPsSlice);
                        
                        if(recoVX != -999999){
                            recoX.ubooneCosmic->Fill(recoVX, weight);
                            recoXDist.ubooneCosmic->Fill(recoVX);
                            recoY.ubooneCosmic->Fill(recoVY, weight);
                            recoYDist.ubooneCosmic->Fill(recoVY);
                            recoZ.ubooneCosmic->Fill(recoVZ, weight);
                            recoZDist.ubooneCosmic->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.ubooneCosmic->Fill(recoVX, weight);
                                recoXDist_low.ubooneCosmic->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.ubooneCosmic->Fill(recoVX, weight);
                                recoXDist_high.ubooneCosmic->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.ubooneCosmic->Fill(recoVY, weight);
                                recoYDist_low.ubooneCosmic->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.ubooneCosmic->Fill(recoVY, weight);
                                recoYDist_high.ubooneCosmic->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.ubooneCosmic->Fill(recoVZ, weight);
                                recoZDist_low.ubooneCosmic->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.ubooneCosmic->Fill(recoVZ, weight);
                                recoZDist_high.ubooneCosmic->Fill(recoVZ);
                            //}
                        }

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.ubooneCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.ubooneCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.ubooneCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.ubooneCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.ubooneCosmic->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.ubooneCosmic->Fill(highestEnergy_completeness);
                            pfpPurity.ubooneCosmic->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.ubooneCosmic->Fill(highestEnergy_purity);
                        }

                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuECosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.nuECosmic->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.nuECosmic->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.nuECosmic->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.nuECosmic->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.nuECosmic->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.nuECosmic->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.nuECosmic->Fill(numPFPsSlice);
                        
                        if(recoVX != -999999){
                            recoX.nuECosmic->Fill(recoVX, weight);
                            recoXDist.nuECosmic->Fill(recoVX);
                            recoY.nuECosmic->Fill(recoVY, weight);
                            recoYDist.nuECosmic->Fill(recoVY);
                            recoZ.nuECosmic->Fill(recoVZ, weight);
                            recoZDist.nuECosmic->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.nuECosmic->Fill(recoVX, weight);
                                recoXDist_low.nuECosmic->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.nuECosmic->Fill(recoVX, weight);
                                recoXDist_high.nuECosmic->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.nuECosmic->Fill(recoVY, weight);
                                recoYDist_low.nuECosmic->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.nuECosmic->Fill(recoVY, weight);
                                recoYDist_high.nuECosmic->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.nuECosmic->Fill(recoVZ, weight);
                                recoZDist_low.nuECosmic->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.nuECosmic->Fill(recoVZ, weight);
                                recoZDist_high.nuECosmic->Fill(recoVZ);
                            //}
                        }

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.nuECosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuECosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuECosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuECosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.nuECosmic->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.nuECosmic->Fill(highestEnergy_completeness);
                            pfpPurity.nuECosmic->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.nuECosmic->Fill(highestEnergy_purity);
                        }
                    }

                } else if(reco_sliceCategory->at(slice) == 1){
                    // This is a signal slice
                    if(DLCurrent == 2 && signal == 1){
                        sliceCompleteness.currentSignal->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.currentSignal->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.currentSignal->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.currentSignal->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.currentSignal->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.currentSignal->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.currentSignal->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.currentSignal->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            std::cout << "WEIGHT CHECK, should be 1 here - " << weight << std::endl;
                            if(weight != 1){ std::cout << "signal = " << signal << ", DLCurrent = " << DLCurrent << std::endl; std::cout << "Slice Category = " << reco_sliceCategory->at(slice) << ", Slice Interaction = " << reco_sliceInteraction->at(slice) << std::endl;}
                            ERecoSumThetaReco.currentSignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.currentSignal->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.currentSignal->Fill(highestEnergy_completeness);
                            pfpPurity.currentSignal->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.currentSignal->Fill(highestEnergy_purity);
                            deltaTheta.currentSignal->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignal->Fill(angleDifference);
                      
                            ETrueThetaReco.currentSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.currentSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.currentSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.currentSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.currentSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.currentSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));
                            ETrue.currentSignal->Fill(recoilElectron_energy, weight);
                            ETrueDist.currentSignal->Fill(recoilElectron_energy);
                            ThetaTrue.currentSignal->Fill(recoilElectron_angle, weight);
                            ThetaTrueDist.currentSignal->Fill(recoilElectron_angle);
                        
                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceBDT->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceBDT->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceBDT->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceBDT_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceBDT_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceBDT_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceBDT_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceBDT_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceBDT_high->Fill(recoVZ, angleDifference); 
                            }                        
                        }
                        
                        if(recoVX != -999999){ 
                            deltaX.currentSignal->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.currentSignal->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.currentSignal->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.currentSignal->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.currentSignal->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.currentSignal->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.currentSignal->Fill(deltaRVal, weight);
                            deltaRDist.currentSignal->Fill(deltaRVal);
                            
                            recoX.currentSignal->Fill(recoVX, weight);
                            recoXDist.currentSignal->Fill(recoVX);
                            recoY.currentSignal->Fill(recoVY, weight);
                            recoYDist.currentSignal->Fill(recoVY);
                            recoZ.currentSignal->Fill(recoVZ, weight);
                            recoZDist.currentSignal->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.currentSignal->Fill(recoVX, weight);
                                recoXDist_low.currentSignal->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.currentSignal->Fill(recoVX, weight);
                                recoXDist_high.currentSignal->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.currentSignal->Fill(recoVY, weight);
                                recoYDist_low.currentSignal->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.currentSignal->Fill(recoVY, weight);
                                recoYDist_high.currentSignal->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.currentSignal->Fill(recoVZ, weight);
                                recoZDist_low.currentSignal->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.currentSignal->Fill(recoVZ, weight);
                                recoZDist_high.currentSignal->Fill(recoVZ);
                            //}
                        }

                    } else if(DLCurrent == 0 && signal == 1){
                        sliceCompleteness.ubooneSignal->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.ubooneSignal->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.ubooneSignal->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.ubooneSignal->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.ubooneSignal->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.ubooneSignal->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.ubooneSignal->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.ubooneSignal->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.ubooneSignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.ubooneSignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.ubooneSignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.ubooneSignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.ubooneSignal->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.ubooneSignal->Fill(highestEnergy_completeness);
                            pfpPurity.ubooneSignal->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.ubooneSignal->Fill(highestEnergy_purity);
                            deltaTheta.ubooneSignal->Fill(angleDifference, weight);
                            deltaThetaDist.ubooneSignal->Fill(angleDifference);
                      
                            ETrueThetaReco.ubooneSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.ubooneSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.ubooneSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.ubooneSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.ubooneSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.ubooneSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));
                            ETrue.ubooneSignal->Fill(recoilElectron_energy, weight);
                            ETrueDist.ubooneSignal->Fill(recoilElectron_energy);
                            ThetaTrue.ubooneSignal->Fill(recoilElectron_angle, weight);
                            ThetaTrueDist.ubooneSignal->Fill(recoilElectron_angle);
                            
                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceDLUboone->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLUboone->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLUboone->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceDLUboone_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceDLUboone_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceDLUboone_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceDLUboone_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceDLUboone_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceDLUboone_high->Fill(recoVZ, angleDifference); 
                            }                        
                        }
                        
                        if(recoVX != -999999){
                            deltaX.ubooneSignal->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.ubooneSignal->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.ubooneSignal->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.ubooneSignal->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.ubooneSignal->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.ubooneSignal->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.ubooneSignal->Fill(deltaRVal, weight);
                            deltaRDist.ubooneSignal->Fill(deltaRVal);
                            
                            recoX.ubooneSignal->Fill(recoVX, weight);
                            recoXDist.ubooneSignal->Fill(recoVX);
                            recoY.ubooneSignal->Fill(recoVY, weight);
                            recoYDist.ubooneSignal->Fill(recoVY);
                            recoZ.ubooneSignal->Fill(recoVZ, weight);
                            recoZDist.ubooneSignal->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.ubooneSignal->Fill(recoVX, weight);
                                recoXDist_low.ubooneSignal->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.ubooneSignal->Fill(recoVX, weight);
                                recoXDist_high.ubooneSignal->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.ubooneSignal->Fill(recoVY, weight);
                                recoYDist_low.ubooneSignal->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.ubooneSignal->Fill(recoVY, weight);
                                recoYDist_high.ubooneSignal->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.ubooneSignal->Fill(recoVZ, weight);
                                recoZDist_low.ubooneSignal->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.ubooneSignal->Fill(recoVZ, weight);
                                recoZDist_high.ubooneSignal->Fill(recoVZ);
                            //}
                        }

                    } else if(DLCurrent == 5  && signal == 1){
                        sliceCompleteness.nuESignal->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.nuESignal->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.nuESignal->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.nuESignal->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.nuESignal->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.nuESignal->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.nuESignal->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.nuESignal->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.nuESignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.nuESignal->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.nuESignal->Fill(highestEnergy_completeness);
                            pfpPurity.nuESignal->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.nuESignal->Fill(highestEnergy_purity);
                            deltaTheta.nuESignal->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignal->Fill(angleDifference);
                      
                            ETrueThetaReco.nuESignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.nuESignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.nuESignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.nuESignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.nuESignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.nuESignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));
                            ETrue.nuESignal->Fill(recoilElectron_energy, weight);
                            ETrueDist.nuESignal->Fill(recoilElectron_energy);
                            ThetaTrue.nuESignal->Fill(recoilElectron_angle, weight);
                            ThetaTrueDist.nuESignal->Fill(recoilElectron_angle);
                            
                            if(recoVX != -999999){
                
                                xCoordAngleDifferenceDLNuE->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLNuE->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLNuE->Fill(recoVZ, angleDifference);

                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceDLNuE_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceDLNuE_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceDLNuE_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceDLNuE_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceDLNuE_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceDLNuE_high->Fill(recoVZ, angleDifference); 
                            }                        
                        }
                            
                        if(recoVX != -999999){ 
                            deltaX.nuESignal->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.nuESignal->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.nuESignal->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.nuESignal->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.nuESignal->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.nuESignal->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.nuESignal->Fill(deltaRVal, weight);
                            deltaRDist.nuESignal->Fill(deltaRVal);
                            
                            recoX.nuESignal->Fill(recoVX, weight);
                            recoXDist.nuESignal->Fill(recoVX);
                            recoY.nuESignal->Fill(recoVY, weight);
                            recoYDist.nuESignal->Fill(recoVY);
                            recoZ.nuESignal->Fill(recoVZ, weight);
                            recoZDist.nuESignal->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.nuESignal->Fill(recoVX, weight);
                                recoXDist_low.nuESignal->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.nuESignal->Fill(recoVX, weight);
                                recoXDist_high.nuESignal->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.nuESignal->Fill(recoVY, weight);
                                recoYDist_low.nuESignal->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.nuESignal->Fill(recoVY, weight);
                                recoYDist_high.nuESignal->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.nuESignal->Fill(recoVZ, weight);
                                recoZDist_low.nuESignal->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.nuESignal->Fill(recoVZ, weight);
                                recoZDist_high.nuESignal->Fill(recoVZ);
                            //}
                        }
                    }
                } else if(reco_sliceCategory->at(slice) == 2){
                    // This is a fuzzy signal slice
                    if(DLCurrent == 2 && signal == 1){
                        sliceCompleteness.currentSignalFuzzy->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.currentSignalFuzzy->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.currentSignalFuzzy->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.currentSignalFuzzy->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.currentSignalFuzzy->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.currentSignalFuzzy->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.currentSignalFuzzy->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.currentSignalFuzzy->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.currentSignalFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignalFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignalFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignalFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.currentSignalFuzzy->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.currentSignalFuzzy->Fill(highestEnergy_completeness);
                            pfpPurity.currentSignalFuzzy->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.currentSignalFuzzy->Fill(highestEnergy_purity);
                            deltaTheta.currentSignalFuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignalFuzzy->Fill(angleDifference);
                      
                            ETrueThetaReco.currentSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.currentSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.currentSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.currentSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.currentSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.currentSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));
                            ETrue.currentSignalFuzzy->Fill(recoilElectron_energy, weight);
                            ETrueDist.currentSignalFuzzy->Fill(recoilElectron_energy);
                            ThetaTrue.currentSignalFuzzy->Fill(recoilElectron_angle, weight);
                            ThetaTrueDist.currentSignalFuzzy->Fill(recoilElectron_angle);
                            
                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceBDT->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceBDT->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceBDT->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceBDT_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceBDT_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceBDT_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceBDT_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceBDT_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceBDT_high->Fill(recoVZ, angleDifference); 
                            }
                        }
                            
                        if(recoVX != -999999){ 
                            deltaX.currentSignalFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.currentSignalFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.currentSignalFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.currentSignalFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.currentSignalFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.currentSignalFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.currentSignalFuzzy->Fill(deltaRVal, weight);
                            deltaRDist.currentSignalFuzzy->Fill(deltaRVal);
                            
                            recoX.currentSignalFuzzy->Fill(recoVX, weight);
                            recoXDist.currentSignalFuzzy->Fill(recoVX);
                            recoY.currentSignalFuzzy->Fill(recoVY, weight);
                            recoYDist.currentSignalFuzzy->Fill(recoVY);
                            recoZ.currentSignalFuzzy->Fill(recoVZ, weight);
                            recoZDist.currentSignalFuzzy->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.currentSignalFuzzy->Fill(recoVX, weight);
                                recoXDist_low.currentSignalFuzzy->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.currentSignalFuzzy->Fill(recoVX, weight);
                                recoXDist_high.currentSignalFuzzy->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.currentSignalFuzzy->Fill(recoVY, weight);
                                recoYDist_low.currentSignalFuzzy->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.currentSignalFuzzy->Fill(recoVY, weight);
                                recoYDist_high.currentSignalFuzzy->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.currentSignalFuzzy->Fill(recoVZ, weight);
                                recoZDist_low.currentSignalFuzzy->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.currentSignalFuzzy->Fill(recoVZ, weight);
                                recoZDist_high.currentSignalFuzzy->Fill(recoVZ);
                            //}
                        }
                    } else if(DLCurrent == 0 && signal == 1){
                        sliceCompleteness.ubooneSignalFuzzy->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.ubooneSignalFuzzy->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.ubooneSignalFuzzy->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.ubooneSignalFuzzy->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.ubooneSignalFuzzy->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.ubooneSignalFuzzy->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.ubooneSignalFuzzy->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.ubooneSignalFuzzy->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.ubooneSignalFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.ubooneSignalFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.ubooneSignalFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.ubooneSignalFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.ubooneSignalFuzzy->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.ubooneSignalFuzzy->Fill(highestEnergy_completeness);
                            pfpPurity.ubooneSignalFuzzy->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.ubooneSignalFuzzy->Fill(highestEnergy_purity);
                            deltaTheta.ubooneSignalFuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.ubooneSignalFuzzy->Fill(angleDifference);
                      
                            ETrueThetaReco.ubooneSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.ubooneSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.ubooneSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.ubooneSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.ubooneSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.ubooneSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));
                            ETrue.ubooneSignalFuzzy->Fill(recoilElectron_energy, weight);
                            ETrueDist.ubooneSignalFuzzy->Fill(recoilElectron_energy);
                            ThetaTrue.ubooneSignalFuzzy->Fill(recoilElectron_angle, weight);
                            ThetaTrueDist.ubooneSignalFuzzy->Fill(recoilElectron_angle);
                        
                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceDLUboone->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLUboone->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLUboone->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceDLUboone_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceDLUboone_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceDLUboone_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceDLUboone_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceDLUboone_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceDLUboone_high->Fill(recoVZ, angleDifference); 
                            }                        
                        }
                            
                        if(recoVX != -999999){ 
                            deltaX.ubooneSignalFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.ubooneSignalFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.ubooneSignalFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.ubooneSignalFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.ubooneSignalFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.ubooneSignalFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.ubooneSignalFuzzy->Fill(deltaRVal, weight);
                            deltaRDist.ubooneSignalFuzzy->Fill(deltaRVal);
                            
                            recoX.ubooneSignalFuzzy->Fill(recoVX, weight);
                            recoXDist.ubooneSignalFuzzy->Fill(recoVX);
                            recoY.ubooneSignalFuzzy->Fill(recoVY, weight);
                            recoYDist.ubooneSignalFuzzy->Fill(recoVY);
                            recoZ.ubooneSignalFuzzy->Fill(recoVZ, weight);
                            recoZDist.ubooneSignalFuzzy->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.ubooneSignalFuzzy->Fill(recoVX, weight);
                                recoXDist_low.ubooneSignalFuzzy->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.ubooneSignalFuzzy->Fill(recoVX, weight);
                                recoXDist_high.ubooneSignalFuzzy->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.ubooneSignalFuzzy->Fill(recoVY, weight);
                                recoYDist_low.ubooneSignalFuzzy->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.ubooneSignalFuzzy->Fill(recoVY, weight);
                                recoYDist_high.ubooneSignalFuzzy->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.ubooneSignalFuzzy->Fill(recoVZ, weight);
                                recoZDist_low.ubooneSignalFuzzy->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.ubooneSignalFuzzy->Fill(recoVZ, weight);
                                recoZDist_high.ubooneSignalFuzzy->Fill(recoVZ);
                            //}
                        }
                    } else if(DLCurrent == 5 && signal == 1){
                        sliceCompleteness.nuESignalFuzzy->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.nuESignalFuzzy->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.nuESignalFuzzy->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.nuESignalFuzzy->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.nuESignalFuzzy->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.nuESignalFuzzy->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.nuESignalFuzzy->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.nuESignalFuzzy->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.nuESignalFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignalFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignalFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignalFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.nuESignalFuzzy->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.nuESignalFuzzy->Fill(highestEnergy_completeness);
                            pfpPurity.nuESignalFuzzy->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.nuESignalFuzzy->Fill(highestEnergy_purity);
                            deltaTheta.nuESignalFuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignalFuzzy->Fill(angleDifference);
                      
                            ETrueThetaReco.nuESignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.nuESignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.nuESignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.nuESignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.nuESignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.nuESignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));
                            ETrue.nuESignalFuzzy->Fill(recoilElectron_energy, weight);
                            ETrueDist.nuESignalFuzzy->Fill(recoilElectron_energy);
                            ThetaTrue.nuESignalFuzzy->Fill(recoilElectron_angle, weight);
                            ThetaTrueDist.nuESignalFuzzy->Fill(recoilElectron_angle);
                            
                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceDLNuE->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLNuE->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLNuE->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceDLNuE_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceDLNuE_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceDLNuE_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceDLNuE_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceDLNuE_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceDLNuE_high->Fill(recoVZ, angleDifference); 
                            }                        
                        }
                            
                        if(recoVX != -999999){ 
                            deltaX.nuESignalFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.nuESignalFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.nuESignalFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.nuESignalFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.nuESignalFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.nuESignalFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.nuESignalFuzzy->Fill(deltaRVal, weight);
                            deltaRDist.nuESignalFuzzy->Fill(deltaRVal);
                            
                            recoX.nuESignalFuzzy->Fill(recoVX, weight);
                            recoXDist.nuESignalFuzzy->Fill(recoVX);
                            recoY.nuESignalFuzzy->Fill(recoVY, weight);
                            recoYDist.nuESignalFuzzy->Fill(recoVY);
                            recoZ.nuESignalFuzzy->Fill(recoVZ, weight);
                            recoZDist.nuESignalFuzzy->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.nuESignalFuzzy->Fill(recoVX, weight);
                                recoXDist_low.nuESignalFuzzy->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.nuESignalFuzzy->Fill(recoVX, weight);
                                recoXDist_high.nuESignalFuzzy->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.nuESignalFuzzy->Fill(recoVY, weight);
                                recoYDist_low.nuESignalFuzzy->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.nuESignalFuzzy->Fill(recoVY, weight);
                                recoYDist_high.nuESignalFuzzy->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.nuESignalFuzzy->Fill(recoVZ, weight);
                                recoZDist_low.nuESignalFuzzy->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.nuESignalFuzzy->Fill(recoVZ, weight);
                                recoZDist_high.nuESignalFuzzy->Fill(recoVZ);
                            //}
                      
                        }
                    }
                } else if(reco_sliceCategory->at(slice) == 3){
                    // This is a BNB slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentBNB->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.currentBNB->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.currentBNB->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.currentBNB->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.currentBNB->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.currentBNB->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.currentBNB->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.currentBNB->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.currentBNB->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentBNB->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentBNB->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentBNB->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.currentBNB->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.currentBNB->Fill(highestEnergy_completeness);
                            pfpPurity.currentBNB->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.currentBNB->Fill(highestEnergy_purity);
                        }
                            
                        if(recoVX != -999999){ 
                            deltaX.currentBNB->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.currentBNB->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.currentBNB->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.currentBNB->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.currentBNB->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.currentBNB->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.currentBNB->Fill(deltaRVal, weight);
                            deltaRDist.currentBNB->Fill(deltaRVal);
                            
                            recoX.currentBNB->Fill(recoVX, weight);
                            recoXDist.currentBNB->Fill(recoVX);
                            recoY.currentBNB->Fill(recoVY, weight);
                            recoYDist.currentBNB->Fill(recoVY);
                            recoZ.currentBNB->Fill(recoVZ, weight);
                            recoZDist.currentBNB->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.currentBNB->Fill(recoVX, weight);
                                recoXDist_low.currentBNB->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.currentBNB->Fill(recoVX, weight);
                                recoXDist_high.currentBNB->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.currentBNB->Fill(recoVY, weight);
                                recoYDist_low.currentBNB->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.currentBNB->Fill(recoVY, weight);
                                recoYDist_high.currentBNB->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.currentBNB->Fill(recoVZ, weight);
                                recoZDist_low.currentBNB->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.currentBNB->Fill(recoVZ, weight);
                                recoZDist_high.currentBNB->Fill(recoVZ);
                            //}
                        }

                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneBNB->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.ubooneBNB->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.ubooneBNB->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.ubooneBNB->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.ubooneBNB->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.ubooneBNB->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.ubooneBNB->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.ubooneBNB->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.ubooneBNB->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.ubooneBNB->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.ubooneBNB->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.ubooneBNB->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.ubooneBNB->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.ubooneBNB->Fill(highestEnergy_completeness);
                            pfpPurity.ubooneBNB->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.ubooneBNB->Fill(highestEnergy_purity);
                        } 
                           
                        if(recoVX != -999999){ 
                            deltaX.ubooneBNB->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.ubooneBNB->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.ubooneBNB->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.ubooneBNB->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.ubooneBNB->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.ubooneBNB->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.ubooneBNB->Fill(deltaRVal, weight);
                            deltaRDist.ubooneBNB->Fill(deltaRVal);
                            
                            recoX.ubooneBNB->Fill(recoVX, weight);
                            recoXDist.ubooneBNB->Fill(recoVX);
                            recoY.ubooneBNB->Fill(recoVY, weight);
                            recoYDist.ubooneBNB->Fill(recoVY);
                            recoZ.ubooneBNB->Fill(recoVZ, weight);
                            recoZDist.ubooneBNB->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.ubooneBNB->Fill(recoVX, weight);
                                recoXDist_low.ubooneBNB->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.ubooneBNB->Fill(recoVX, weight);
                                recoXDist_high.ubooneBNB->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.ubooneBNB->Fill(recoVY, weight);
                                recoYDist_low.ubooneBNB->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.ubooneBNB->Fill(recoVY, weight);
                                recoYDist_high.ubooneBNB->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.ubooneBNB->Fill(recoVZ, weight);
                                recoZDist_low.ubooneBNB->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.ubooneBNB->Fill(recoVZ, weight);
                                recoZDist_high.ubooneBNB->Fill(recoVZ);
                            //}
                        }

                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuEBNB->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.nuEBNB->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.nuEBNB->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.nuEBNB->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.nuEBNB->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.nuEBNB->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.nuEBNB->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.nuEBNB->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.nuEBNB->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuEBNB->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuEBNB->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuEBNB->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.nuEBNB->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.nuEBNB->Fill(highestEnergy_completeness);
                            pfpPurity.nuEBNB->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.nuEBNB->Fill(highestEnergy_purity);
                        } 
                            
                        if(recoVX != -999999){ 
                            deltaX.nuEBNB->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.nuEBNB->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.nuEBNB->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.nuEBNB->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.nuEBNB->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.nuEBNB->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.nuEBNB->Fill(deltaRVal, weight);
                            deltaRDist.nuEBNB->Fill(deltaRVal);
                            
                            recoX.nuEBNB->Fill(recoVX, weight);
                            recoXDist.nuEBNB->Fill(recoVX);
                            recoY.nuEBNB->Fill(recoVY, weight);
                            recoYDist.nuEBNB->Fill(recoVY);
                            recoZ.nuEBNB->Fill(recoVZ, weight);
                            recoZDist.nuEBNB->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.nuEBNB->Fill(recoVX, weight);
                                recoXDist_low.nuEBNB->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.nuEBNB->Fill(recoVX, weight);
                                recoXDist_high.nuEBNB->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.nuEBNB->Fill(recoVY, weight);
                                recoYDist_low.nuEBNB->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.nuEBNB->Fill(recoVY, weight);
                                recoYDist_high.nuEBNB->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.nuEBNB->Fill(recoVZ, weight);
                                recoZDist_low.nuEBNB->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.nuEBNB->Fill(recoVZ, weight);
                                recoZDist_high.nuEBNB->Fill(recoVZ);
                            //}
                        }
                    }
                } else if(reco_sliceCategory->at(slice) == 4){
                    // This is a fuzzy BNB slice
                    if(DLCurrent == 2){
                        sliceCompleteness.currentBNBFuzzy->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.currentBNBFuzzy->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.currentBNBFuzzy->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.currentBNBFuzzy->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.currentBNBFuzzy->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.currentBNBFuzzy->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.currentBNBFuzzy->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.currentBNBFuzzy->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.currentBNBFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentBNBFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentBNBFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentBNBFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.currentBNBFuzzy->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.currentBNBFuzzy->Fill(highestEnergy_completeness);
                            pfpPurity.currentBNBFuzzy->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.currentBNBFuzzy->Fill(highestEnergy_purity);
                        } 
                            
                        if(recoVX != -999999){ 
                            deltaX.currentBNBFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.currentBNBFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.currentBNBFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.currentBNBFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.currentBNBFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.currentBNBFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.currentBNBFuzzy->Fill(deltaRVal, weight);
                            deltaRDist.currentBNBFuzzy->Fill(deltaRVal);
                            
                            recoX.currentBNBFuzzy->Fill(recoVX, weight);
                            recoXDist.currentBNBFuzzy->Fill(recoVX);
                            recoY.currentBNBFuzzy->Fill(recoVY, weight);
                            recoYDist.currentBNBFuzzy->Fill(recoVY);
                            recoZ.currentBNBFuzzy->Fill(recoVZ, weight);
                            recoZDist.currentBNBFuzzy->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.currentBNBFuzzy->Fill(recoVX, weight);
                                recoXDist_low.currentBNBFuzzy->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.currentBNBFuzzy->Fill(recoVX, weight);
                                recoXDist_high.currentBNBFuzzy->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.currentBNBFuzzy->Fill(recoVY, weight);
                                recoYDist_low.currentBNBFuzzy->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.currentBNBFuzzy->Fill(recoVY, weight);
                                recoYDist_high.currentBNBFuzzy->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.currentBNBFuzzy->Fill(recoVZ, weight);
                                recoZDist_low.currentBNBFuzzy->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.currentBNBFuzzy->Fill(recoVZ, weight);
                                recoZDist_high.currentBNBFuzzy->Fill(recoVZ);
                            //}
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness.ubooneBNBFuzzy->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.ubooneBNBFuzzy->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.ubooneBNBFuzzy->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.ubooneBNBFuzzy->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.ubooneBNBFuzzy->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.ubooneBNBFuzzy->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.ubooneBNBFuzzy->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.ubooneBNBFuzzy->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.ubooneBNBFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.ubooneBNBFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.ubooneBNBFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.ubooneBNBFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.ubooneBNBFuzzy->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.ubooneBNBFuzzy->Fill(highestEnergy_completeness);
                            pfpPurity.ubooneBNBFuzzy->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.ubooneBNBFuzzy->Fill(highestEnergy_purity);
                        } 
                            
                        if(recoVX != -999999){ 
                            deltaX.ubooneBNBFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.ubooneBNBFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.ubooneBNBFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.ubooneBNBFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.ubooneBNBFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.ubooneBNBFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.ubooneBNBFuzzy->Fill(deltaRVal, weight);
                            deltaRDist.ubooneBNBFuzzy->Fill(deltaRVal);
                            
                            recoX.ubooneBNBFuzzy->Fill(recoVX, weight);
                            recoXDist.ubooneBNBFuzzy->Fill(recoVX);
                            recoY.ubooneBNBFuzzy->Fill(recoVY, weight);
                            recoYDist.ubooneBNBFuzzy->Fill(recoVY);
                            recoZ.ubooneBNBFuzzy->Fill(recoVZ, weight);
                            recoZDist.ubooneBNBFuzzy->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.ubooneBNBFuzzy->Fill(recoVX, weight);
                                recoXDist_low.ubooneBNBFuzzy->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.ubooneBNBFuzzy->Fill(recoVX, weight);
                                recoXDist_high.ubooneBNBFuzzy->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.ubooneBNBFuzzy->Fill(recoVY, weight);
                                recoYDist_low.ubooneBNBFuzzy->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.ubooneBNBFuzzy->Fill(recoVY, weight);
                                recoYDist_high.ubooneBNBFuzzy->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.ubooneBNBFuzzy->Fill(recoVZ, weight);
                                recoZDist_low.ubooneBNBFuzzy->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.ubooneBNBFuzzy->Fill(recoVZ, weight);
                                recoZDist_high.ubooneBNBFuzzy->Fill(recoVZ);
                            //}
                        }
                    } else if(DLCurrent == 5){
                        sliceCompleteness.nuEBNBFuzzy->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity.nuEBNBFuzzy->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore.nuEBNBFuzzy->Fill(reco_sliceScore->at(slice), weight);
                        sliceCompletenessDist.nuEBNBFuzzy->Fill(reco_sliceCompleteness->at(slice));
                        slicePurityDist.nuEBNBFuzzy->Fill(reco_slicePurity->at(slice));
                        sliceCRUMBSScoreDist.nuEBNBFuzzy->Fill(reco_sliceScore->at(slice));
                        
                        sliceNumPFPs.nuEBNBFuzzy->Fill(numPFPsSlice, weight);
                        sliceNumPFPsDist.nuEBNBFuzzy->Fill(numPFPsSlice);

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.nuEBNBFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuEBNBFuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuEBNBFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuEBNBFuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                            pfpCompleteness.nuEBNBFuzzy->Fill(highestEnergy_completeness, weight);
                            pfpCompletenessDist.nuEBNBFuzzy->Fill(highestEnergy_completeness);
                            pfpPurity.nuEBNBFuzzy->Fill(highestEnergy_purity, weight);
                            pfpPurityDist.nuEBNBFuzzy->Fill(highestEnergy_purity);
                        } 
                            
                        if(recoVX != -999999){ 
                            deltaX.nuEBNBFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)), weight);
                            deltaXDist.nuEBNBFuzzy->Fill((recoVX - reco_sliceTrueVX->at(slice)));
                            deltaY.nuEBNBFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)), weight);
                            deltaYDist.nuEBNBFuzzy->Fill((recoVY - reco_sliceTrueVY->at(slice)));
                            deltaZ.nuEBNBFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)), weight);
                            deltaZDist.nuEBNBFuzzy->Fill((recoVZ - reco_sliceTrueVZ->at(slice)));
                            double deltaRVal = std::sqrt(((recoVX - reco_sliceTrueVX->at(slice)) * (recoVX - reco_sliceTrueVX->at(slice))) + ((recoVY - reco_sliceTrueVY->at(slice)) * (recoVY - reco_sliceTrueVY->at(slice))) + ((recoVZ - reco_sliceTrueVZ->at(slice))* (recoVZ - reco_sliceTrueVZ->at(slice)))); 
                            deltaR.nuEBNBFuzzy->Fill(deltaRVal, weight);
                            deltaRDist.nuEBNBFuzzy->Fill(deltaRVal);
                            
                            recoX.nuEBNBFuzzy->Fill(recoVX, weight);
                            recoXDist.nuEBNBFuzzy->Fill(recoVX);
                            recoY.nuEBNBFuzzy->Fill(recoVY, weight);
                            recoYDist.nuEBNBFuzzy->Fill(recoVY);
                            recoZ.nuEBNBFuzzy->Fill(recoVZ, weight);
                            recoZDist.nuEBNBFuzzy->Fill(recoVZ);
                            
                            //if(recoVX >= xMin && recoVX <= xMin+20){
                                recoX_low.nuEBNBFuzzy->Fill(recoVX, weight);
                                recoXDist_low.nuEBNBFuzzy->Fill(recoVX);
                            //} else if(recoVX <= xMax && recoVX >= xMax-20){
                                recoX_high.nuEBNBFuzzy->Fill(recoVX, weight);
                                recoXDist_high.nuEBNBFuzzy->Fill(recoVX);
                            //}
                                
                            //if(recoVY >= yMin && recoVY <= yMin+20){
                                recoY_low.nuEBNBFuzzy->Fill(recoVY, weight);
                                recoYDist_low.nuEBNBFuzzy->Fill(recoVY);
                            //} else if(recoVY <= yMax && recoVY >= yMax-20){
                                recoY_high.nuEBNBFuzzy->Fill(recoVY, weight);
                                recoYDist_high.nuEBNBFuzzy->Fill(recoVY);
                            //}
                                
                            //if(recoVZ >= zMin && recoVZ <= zMin+20){
                                recoZ_low.nuEBNBFuzzy->Fill(recoVZ, weight);
                                recoZDist_low.nuEBNBFuzzy->Fill(recoVZ);
                            //} else if(recoVZ <= zMax && recoVZ >= zMax-20){
                                recoZ_high.nuEBNBFuzzy->Fill(recoVZ, weight);
                                recoZDist_high.nuEBNBFuzzy->Fill(recoVZ);
                            //}
                        }
                    }
                }

                //printf("_______________________________________________________________\n");
            }
        
        }

        

 
    }

    int drawLine = 1;
    int left = 0;
    int right = 1;

    styleDrawAll(trueETheta2, 999, 999, 999, 999, (base_path + "trueETheta2_weighted.pdf").c_str(), "bottomRight", &drawLine, &right, true, true, false, false, false, true, false, false);
    
    styleDrawAll(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, true, true);
    styleDrawAll(sliceCompletenessDist, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, true, true);
    styleDrawAll(slicePurityDist, 999, 999, 999, 999, (base_path + "slicePurity_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(sliceCRUMBSScore, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, true, true);
    styleDrawAll(sliceCRUMBSScoreDist, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(sliceCRUMBSScore, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, true, true);
    styleDrawAll(sliceNumPFPsDist, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);

    styleDrawAll(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, true, true);
    styleDrawAll(ERecoSumThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, true, true);
    styleDrawAll(ERecoHighestThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);

    styleDrawAll(ETrueThetaReco, 999, 999, 999, 999, (base_path + "ETrueThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true, false);
    styleDrawAll(ETrueThetaRecoDist, 999, 999, 999, 999, (base_path + "ETrueThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(ERecoSumThetaTrue, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, false, false);
    styleDrawAll(ERecoSumThetaTrueDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, true);
    styleDrawAll(ERecoHighestThetaTrue, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false, false);
    styleDrawAll(ERecoHighestThetaTrueDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, false);
    styleDrawAll(ETrue, 999, 999, 999, 999, (base_path + "ETrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);
    styleDrawAll(ETrueDist, 999, 999, 999, 999, (base_path + "ETrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);
    styleDrawAll(ThetaTrue, 999, 999, 999, 999, (base_path + "ThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);
    styleDrawAll(ThetaTrueDist, 999, 999, 999, 999, (base_path + "ThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);

    styleDrawAll(deltaX, 999, 999, 999, 999, (base_path + "deltaX_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaXDist, 999, 999, 999, 999, (base_path + "deltaX_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false);
    styleDrawAll(deltaXDist, 0, 32000, 999, 999, (base_path + "deltaX_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaXDist, 0, 35000, 999, 999, (base_path + "deltaX_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaXDist, 0, 9000, 999, 999, (base_path + "deltaX_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaXDist, 0, 4000, 999, 999, (base_path + "deltaX_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);
    styleDrawBackSig(deltaX, 999, 999, 999, 999, (base_path + "deltaX_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaY, 999, 999, 999, 999, (base_path + "deltaY_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaYDist, 999, 999, 999, 999, (base_path + "deltaY_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false);
    styleDrawBackSig(deltaY, 999, 999, 999, 999, (base_path + "deltaY_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaYDist, 0, 19000, 999, 999, (base_path + "deltaY_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaYDist, 0, 37000, 999, 999, (base_path + "deltaY_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaYDist, 0, 8000, 999, 999, (base_path + "deltaY_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaYDist, 0, 4000, 999, 999, (base_path + "deltaY_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);
    styleDrawAll(deltaZ, 999, 999, 999, 999, (base_path + "deltaZ_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaZDist, 999, 999, 999, 999, (base_path + "deltaZ_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false);
    styleDrawBackSig(deltaZ, 999, 999, 999, 999, (base_path + "deltaZ_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaZDist, 0, 20000, 999, 999, (base_path + "deltaZ_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaZDist, 0, 36000, 999, 999, (base_path + "deltaZ_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaZDist, 0, 9000, 999, 999, (base_path + "deltaZ_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaZDist, 0, 4000, 999, 999, (base_path + "deltaZ_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);
    styleDrawAll(deltaR, 999, 999, 999, 999, (base_path + "deltaR_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaRDist, 999, 999, 999, 999, (base_path + "deltaR_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false);
    styleDrawBackSig(deltaR, 999, 999, 999, 999, (base_path + "deltaR_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaRDist, 0, 25000, 999, 999, (base_path + "deltaR_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaRDist, 0, 69000, 999, 999, (base_path + "deltaR_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaRDist, 0, 13000, 999, 999, (base_path + "deltaR_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaRDist, 0, 5000, 999, 999, (base_path + "deltaR_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);

    styleDrawAll(recoX, 999, 999, 999, 999, (base_path + "recoX_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoXDist, 999, 999, 999, 999, (base_path + "recoX_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoX, 999, 999, 999, 999, (base_path + "recoX_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    efficiency(recoX, 0, 1, 999, 999, (base_path + "recoX_right").c_str(), "bottomLeft", nullptr, &right, -1);
    efficiency(recoX, 0, 1, 999, 999, (base_path + "recoX_left").c_str(), "topLeft", nullptr, &right, 1);
    styleDrawAll(recoX_low, 999, 999, 999, 999, (base_path + "recoX_low_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoXDist_low, 999, 999, 999, 999, (base_path + "recoX_low_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoX_low, 999, 999, 999, 999, (base_path + "recoX_low_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(recoX_high, 999, 999, 999, 999, (base_path + "recoX_high_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoXDist_high, 999, 999, 999, 999, (base_path + "recoX_high_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoX_high, 999, 999, 999, 999, (base_path + "recoX_high_BackSig_weighted.pdf").c_str(), "bottomLeft", true, true, true, true);
    
    styleDrawAll(recoY, 999, 999, 999, 999, (base_path + "recoY_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoYDist, 999, 999, 999, 999, (base_path + "recoY_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoY, 999, 999, 999, 999, (base_path + "recoY_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    efficiency(recoY, 0, 1, 999, 999, (base_path + "recoY_right").c_str(), "bottomLeft", nullptr, &right, -1);
    efficiency(recoY, 0, 1, 999, 999, (base_path + "recoY_left").c_str(), "topLeft", nullptr, &right, 1);
    styleDrawAll(recoY_low, 999, 999, 999, 999, (base_path + "recoY_low_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoYDist_low, 999, 999, 999, 999, (base_path + "recoY_low_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoY_low, 999, 999, 999, 999, (base_path + "recoY_low_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(recoY_high, 999, 999, 999, 999, (base_path + "recoY_high_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoYDist_high, 999, 999, 999, 999, (base_path + "recoY_high_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoY_high, 999, 999, 999, 999, (base_path + "recoY_high_BackSig_weighted.pdf").c_str(), "bottomLeft", true, true, true, true);
    
    styleDrawAll(recoZ, 999, 999, 999, 999, (base_path + "recoZ_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoZDist, 999, 999, 999, 999, (base_path + "recoZ_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoZ, 999, 999, 999, 999, (base_path + "recoZ_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    efficiency(recoZ, 0, 1, 999, 999, (base_path + "recoZ_right").c_str(), "bottomLeft", nullptr, &right, -1);
    efficiency(recoZ, 0, 1, 999, 999, (base_path + "recoZ_left").c_str(), "topLeft", nullptr, &right, 1);
    styleDrawAll(recoZ_low, 999, 999, 999, 999, (base_path + "recoZ_low_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoZDist_low, 999, 999, 999, 999, (base_path + "recoZ_low_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoZ_low, 999, 999, 999, 999, (base_path + "recoZ_low_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(recoZ_high, 999, 999, 999, 999, (base_path + "recoZ_high_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false, true);
    styleDrawAll(recoZDist_high, 999, 999, 999, 999, (base_path + "recoZ_high_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, false, false);
    styleDrawBackSig(recoZ_high, 999, 999, 999, 999, (base_path + "recoZ_high_BackSig_weighted.pdf").c_str(), "bottomLeft", true, true, true, true);

    styleDrawAll(deltaTheta, 999, 999, 999, 999, (base_path + "deltaTheta_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);
    styleDrawAll(deltaThetaDist, 999, 999, 999, 999, (base_path + "deltaTheta_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false);

    styleDrawAll(pfpCompleteness, 999, 999, 999, 999, (base_path + "pfpCompleteness_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(pfpCompletenessDist, 999, 999, 999, 999, (base_path + "pfpCompleteness_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(pfpCompleteness, 999, 999, 999, 999, (base_path + "pfpCompleteness_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(pfpPurity, 999, 999, 999, 999, (base_path + "pfpPurity_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(pfpPurityDist, 999, 999, 999, 999, (base_path + "pfpPurity_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawBackSig(pfpPurity, 999, 999, 999, 999, (base_path + "pfpPurity_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);

    //Test
    //efficiency(ERecoSumThetaReco, 0, 1, 999, 999, (base_path + "ERecoSumThetaReco").c_str(), "bottomRight", nullptr, &right, 1); 
    
    efficiency(sliceCompleteness, 0, 1, 999, 999, (base_path + "sliceCompleteness").c_str(), "topRight", nullptr, &right, -1);
    efficiency(slicePurity, 0, 1, 999, 999, (base_path + "slicePurity").c_str(), "topRight", nullptr, &right, -1);
    efficiency(sliceCRUMBSScore, 0, 1, 999, 999, (base_path + "sliceCRUMBSScore").c_str(), "topRight", nullptr, &right, -1);
    efficiency(sliceNumPFPs, 0, 1, 999, 999, (base_path + "sliceNumPFPs").c_str(), "bottomRight", nullptr, &right, 1);

    efficiency(ERecoSumThetaReco, 0, 1, 999, 999, (base_path + "ERecoSumThetaReco").c_str(), "bottomRight", nullptr, &right, 1);
    efficiency(ERecoHighestThetaReco, 0, 1, 999, 999, (base_path + "ERecoHighestThetaReco").c_str(), "bottomRight", nullptr, &right, 1);

    efficiency(pfpCompleteness, 0, 1, 999, 999, (base_path + "pfpCompleteness").c_str(), "bottomLeft", nullptr, &right, -1);
    efficiency(pfpPurity, 0, 1, 999, 999, (base_path + "pfpPurity").c_str(), "bottomLeft", nullptr, &right, -1);

    efficiency(recoX_low, 0, 1, 999, 999, (base_path + "recoX_low").c_str(), "bottomLeft", nullptr, &right, -1);
    efficiency(recoY_low, 0, 1, 999, 999, (base_path + "recoY_low").c_str(), "bottomRight", nullptr, &right, -1);
    efficiency(recoZ_low, 0, 1, 999, 999, (base_path + "recoZ_low").c_str(), "bottomRight", nullptr, &right, -1);
    
    efficiency(recoX_high, 0, 1, 999, 999, (base_path + "recoX_high").c_str(), "bottomRight", nullptr, &right, 1);
    efficiency(recoY_high, 0, 1, 999, 999, (base_path + "recoY_high").c_str(), "bottomLeft", nullptr, &right, 1);
    efficiency(recoZ_high, 0, 1, 999, 999, (base_path + "recoZ_high").c_str(), "bottomRight", nullptr, &right, 1);

    TwoDHistDraw(xCoordAngleDifferenceBDT_low, (base_path + "angleDiffPosition_x_BDT_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(xCoordAngleDifferenceBDT_high, (base_path + "angleDiffPosition_x_BDT_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceBDT_low, (base_path + "angleDiffPosition_y_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceBDT_high, (base_path + "angleDiffPosition_y_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceBDT_low, (base_path + "angleDiffPosition_z_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceBDT_high, (base_path + "angleDiffPosition_z_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifferenceDLUboone_low, (base_path + "angleDiffPosition_x_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(xCoordAngleDifferenceDLUboone_high, (base_path + "angleDiffPosition_x_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLUboone_low, (base_path + "angleDiffPosition_y_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLUboone_high, (base_path + "angleDiffPosition_y_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLUboone_low, (base_path + "angleDiffPosition_z_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLUboone_high, (base_path + "angleDiffPosition_z_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifferenceDLNuE_low, (base_path + "angleDiffPosition_x_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(xCoordAngleDifferenceDLNuE_high, (base_path + "angleDiffPosition_x_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLNuE_low, (base_path + "angleDiffPosition_y_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLNuE_high, (base_path + "angleDiffPosition_y_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLNuE_low, (base_path + "angleDiffPosition_z_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLNuE_high, (base_path + "angleDiffPosition_z_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifferenceBDT, (base_path + "angleDiffPosition_x_BDT.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceBDT, (base_path + "angleDiffPosition_y_BDT.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceBDT, (base_path + "angleDiffPosition_z_BDT.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifferenceDLUboone, (base_path + "angleDiffPosition_x_DLUboone.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLUboone, (base_path + "angleDiffPosition_y_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLUboone, (base_path + "angleDiffPosition_z_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifferenceDLNuE, (base_path + "angleDiffPosition_x_DLNuE.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLNuE, (base_path + "angleDiffPosition_y_DLNuE.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLNuE, (base_path + "angleDiffPosition_z_DLNuE.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");

    printf("Number of Events\nUnweighted BDT: Cosmic = %f, BNB = %f, Nu+E = %f\n", numEvents_BDTCosmic, numEvents_BDTBNB, numEvents_BDTNuE);
    printf("Unweighted DL Nu+E: Cosmic = %f, BNB = %f, Nu+E = %f\n", numEvents_DLNuECosmic, numEvents_DLNuEBNB, numEvents_DLNuENuE);
    printf("Weighted BDT: Cosmic = %f, BNB = %f, Nu+E = %f\n", (numEvents_BDTCosmic * weights.cosmicsCurrent), (numEvents_BDTBNB * weights.BNBCurrent), (numEvents_BDTNuE * weights.signalCurrent));
    printf("Weighted DL Nu+E: Cosmic = %f, BNB = %f, Nu+E = %f\n", (numEvents_DLNuECosmic * weights.cosmicsNuE), (numEvents_DLNuEBNB * weights.BNBNuE), (numEvents_DLNuENuE * weights.signalNuE));
    double totalEvent_BDT = ((numEvents_BDTCosmic * weights.cosmicsCurrent) + (numEvents_BDTBNB * weights.BNBCurrent) + (numEvents_BDTNuE * weights.signalCurrent));
    double cosmicPerc_BDT = ((100 * numEvents_BDTCosmic * weights.cosmicsCurrent)/totalEvent_BDT);
    double BNBPerc_BDT = ((100 * numEvents_BDTBNB * weights.BNBCurrent)/totalEvent_BDT);
    double NuEPerc_BDT = ((100 * numEvents_BDTNuE * weights.signalCurrent)/totalEvent_BDT);
    printf("Event Rates:\nBDT: Cosmic = %f, BNB = %f, Nu+E = %f\n", cosmicPerc_BDT, BNBPerc_BDT, NuEPerc_BDT);

    std::cout << "counter = " << counteraaaaaa << std::endl;
}
