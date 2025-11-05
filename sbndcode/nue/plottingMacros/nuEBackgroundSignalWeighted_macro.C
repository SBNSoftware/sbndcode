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

void styleDrawAll(histGroup_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr,
                  bool includeSignal = true, bool includeBNB = true, bool includeCosmic = true)
{
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    //--- Style setup
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

    //--- Axis ranges
    if((ymin != 999) && (ymax != 999)) hists.currentSignal->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) hists.currentSignal->GetXaxis()->SetRangeUser(xmin, xmax);

    //--- Dynamic Y range
    double maxYValue = std::max({
        hists.currentSignal->GetMaximum(), hists.ubooneSignal->GetMaximum(), hists.nuESignal->GetMaximum(),
        hists.currentCosmic->GetMaximum(), hists.ubooneCosmic->GetMaximum(), hists.nuECosmic->GetMaximum(),
        hists.currentSignalFuzzy->GetMaximum(), hists.ubooneSignalFuzzy->GetMaximum(), hists.nuESignalFuzzy->GetMaximum(),
        hists.currentBNB->GetMaximum(), hists.ubooneBNB->GetMaximum(), hists.nuEBNB->GetMaximum(),
        hists.currentBNBFuzzy->GetMaximum(), hists.ubooneBNBFuzzy->GetMaximum(), hists.nuEBNBFuzzy->GetMaximum()
    });
    double yminVal = 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = (maxYValue * 1.1);
        hists.currentSignal->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    //--- Drawing block
    bool first = true;
    auto draw = [&](TH1* hist){
        if (!hist) return;
        hist->Draw(first ? "hist" : "histsame");
        first = false;
    };

    if (includeSignal) {
        draw(hists.currentSignal);
        draw(hists.ubooneSignal);
        draw(hists.nuESignal);
        draw(hists.currentSignalFuzzy);
        draw(hists.ubooneSignalFuzzy);
        draw(hists.nuESignalFuzzy);
    }
    if (includeBNB) {
        draw(hists.currentBNB);
        draw(hists.ubooneBNB);
        draw(hists.nuEBNB);
        draw(hists.currentBNBFuzzy);
        draw(hists.ubooneBNBFuzzy);
        draw(hists.nuEBNBFuzzy);
    }
    if (includeCosmic) {
        draw(hists.currentCosmic);
        draw(hists.ubooneCosmic);
        draw(hists.nuECosmic);
    }

    //--- Axes and stats
    hists.currentSignal->SetStats(0);
    hists.currentSignal->GetXaxis()->SetTickLength(0.04);
    hists.currentSignal->GetYaxis()->SetTickLength(0.03);
    hists.currentSignal->GetXaxis()->SetTickSize(0.02);
    hists.currentSignal->GetYaxis()->SetTickSize(0.02);

    //--- Legend positioning
    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    if(legendLocation == "topRight"){ Lxmin=0.67; Lymax=0.863; Lxmax=0.87; Lymin=0.640; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.33; Lymin=0.640; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.67; Lymax=0.36; Lxmax=0.87; Lymin=0.137; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymax=0.36; Lxmax=0.33; Lymin=0.137; }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    if (includeSignal) {
        legend->AddEntry(hists.currentSignal, "Signal, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(hists.ubooneSignal, "Signal, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(hists.nuESignal, "Signal, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(hists.currentSignalFuzzy, "Signal Fuzzy, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(hists.ubooneSignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(hists.nuESignalFuzzy, "Signal Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "f");
    }
    if (includeBNB) {
        legend->AddEntry(hists.currentBNB, "BNB, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(hists.ubooneBNB, "BNB, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(hists.nuEBNB, "BNB, Pandora Deep Learning: SBND Nu+E Tune", "f");
        legend->AddEntry(hists.currentBNBFuzzy, "BNB Fuzzy, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(hists.ubooneBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(hists.nuEBNBFuzzy, "BNB Fuzzy, Pandora Deep Learning: SBND Nu+E Tune", "f");
    }
    if (includeCosmic) {
        legend->AddEntry(hists.currentCosmic, "Cosmic, Pandora BDT SBND (without Refinement)", "f");
        legend->AddEntry(hists.ubooneCosmic, "Cosmic, Pandora Deep Learning: #muBooNE/BNB Tune", "f");
        legend->AddEntry(hists.nuECosmic, "Cosmic, Pandora Deep Learning: SBND Nu+E Tune", "f");
    }

    legend->SetTextSize(0.01);
    legend->SetMargin(0.12);
    legend->Draw();
    
    //--- Optional vertical line
    if(drawLine){
        TLine* line = new TLine(1.022, yminVal, 1.022, hists.currentSignal->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        double ymaxValLine = maxYValue*0.95;
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
   
    histGroup_struct effHists;
    effHists.canvas = hists.canvas;
    effHists.baseHist = hists.baseHist;

    effHists.currentCosmic = (TH1F*) hists.currentCosmic->Clone("eff_currentCosmic");
    effHists.currentCosmic->Reset();
    effHists.currentCosmic->GetYaxis()->SetTitle("Efficiency");
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
    effHists.currentBNB->GetYaxis()->SetTitle("Efficiency");
    effHists.currentBNB->GetXaxis()->SetTitle(hists.currentBNB->GetXaxis()->GetTitle());

    effHists.currentBNBFuzzy = (TH1F*) hists.currentBNBFuzzy->Clone("eff_currentBNBFuzzy");
    effHists.currentBNBFuzzy->Reset();
    effHists.currentBNBFuzzy->GetYaxis()->SetTitle("Efficiency");
    effHists.currentBNBFuzzy->GetXaxis()->SetTitle(hists.currentBNB->GetXaxis()->GetTitle());

    effHists.ubooneCosmic = (TH1F*) hists.ubooneCosmic->Clone("eff_ubooneCosmic");
    effHists.ubooneCosmic->Reset();
    effHists.ubooneCosmic->GetYaxis()->SetTitle("Efficiency");
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
    effHists.ubooneBNB->GetYaxis()->SetTitle("Efficiency");
    effHists.ubooneBNB->GetXaxis()->SetTitle(hists.ubooneBNB->GetXaxis()->GetTitle());

    effHists.ubooneBNBFuzzy = (TH1F*) hists.ubooneBNBFuzzy->Clone("eff_ubooneBNBFuzzy");
    effHists.ubooneBNBFuzzy->Reset();
    effHists.ubooneBNBFuzzy->GetYaxis()->SetTitle("Efficiency");
    effHists.ubooneBNBFuzzy->GetXaxis()->SetTitle(hists.ubooneBNB->GetXaxis()->GetTitle());

    effHists.nuECosmic = (TH1F*) hists.nuECosmic->Clone("eff_nuECosmic");
    effHists.nuECosmic->Reset();
    effHists.nuECosmic->GetYaxis()->SetTitle("Efficiency");
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
    effHists.nuEBNB->GetYaxis()->SetTitle("Efficiency");
    effHists.nuEBNB->GetXaxis()->SetTitle(hists.nuEBNB->GetXaxis()->GetTitle());

    effHists.nuEBNBFuzzy = (TH1F*) hists.nuEBNBFuzzy->Clone("eff_nuEBNBFuzzy");
    effHists.nuEBNBFuzzy->Reset();
    effHists.nuEBNBFuzzy->GetYaxis()->SetTitle("Efficiency");
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

            keptSignalCurrent = ((currentSignalTotal - currentSignalSum) + (currentSignalFuzzyTotal - currentSignalFuzzySum)); 
            keptBackgroundCurrent = ((currentBNBTotal - currentBNBSum) + (currentBNBFuzzyTotal - currentBNBFuzzySum) + (currentCosmicTotal - currentCosmicSum));
            keptSignalUboone = ((ubooneSignalTotal - ubooneSignalSum) + (ubooneSignalFuzzyTotal - ubooneSignalFuzzySum)); 
            keptBackgroundUboone = ((ubooneBNBTotal - ubooneBNBSum) + (ubooneBNBFuzzyTotal - ubooneBNBFuzzySum) + (ubooneCosmicTotal - ubooneCosmicSum));
            keptSignalNuE = ((nuESignalTotal - nuESignalSum) + (nuESignalFuzzyTotal - nuESignalFuzzySum)); 
            keptBackgroundNuE = ((nuEBNBTotal - nuEBNBSum) + (nuEBNBFuzzyTotal - nuEBNBFuzzySum) + (nuECosmicTotal - nuECosmicSum));

            currentAllSignalEffVal = (keptSignalCurrent / (currentSignalTotal + currentSignalFuzzyTotal));
            ubooneAllSignalEffVal = (keptSignalUboone / (ubooneSignalTotal + ubooneSignalFuzzyTotal));
            nuEAllSignalEffVal = (keptSignalNuE / (nuESignalTotal + nuESignalFuzzyTotal));

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
        
            keptSignalCurrent = (currentSignalSum + currentSignalFuzzySum);
            keptBackgroundCurrent = (currentCosmicSum + currentBNBSum + currentBNBFuzzySum);  
            keptSignalUboone = (ubooneSignalSum + ubooneSignalFuzzySum);
            keptBackgroundUboone = (ubooneCosmicSum + ubooneBNBSum + ubooneBNBFuzzySum);  
            keptSignalNuE = (nuESignalSum + nuESignalFuzzySum);
            keptBackgroundNuE = (nuECosmicSum + nuEBNBSum + nuEBNBFuzzySum);  
            
            currentAllSignalEffVal = (keptSignalCurrent / (currentSignalTotal + currentSignalFuzzyTotal));
            ubooneAllSignalEffVal = (keptSignalUboone / (ubooneSignalTotal + ubooneSignalFuzzyTotal));
            nuEAllSignalEffVal = (keptSignalNuE / (nuESignalTotal + nuESignalFuzzyTotal));
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
    }
   
    std::string filenameEff = std::string(filename) + "_eff.pdf";
    std::string filenameRej = std::string(filename) + "_rej.pdf";
    styleDrawAll(effHists, ymin, ymax, xmin, xmax, filenameEff.c_str(), legendLocation, drawLine, linePos, true, false, false);
    styleDrawAll(effHists, ymin, ymax, xmin, xmax, filenameRej.c_str(), legendLocation, drawLine, linePos, false, true, true);
}

void nuEBackgroundSignalWeighted_macro(){

    TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nue/srcs/sbndcode/sbndcode/nue/NuEAnalyserOutputSignalMerged.root");
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

    auto deltaX = createHistGroup("deltaX", "#Deltax of Neutrino Vertex in Slice", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaXDist = createHistGroup("deltaXDist", "#Deltax of Neutrino Vertex in Slice (Not Weighted)", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaY = createHistGroup("deltaY", "#Deltay of Neutrino Vertex in Slice", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaYDist = createHistGroup("deltaYDist", "#Deltay of Neutrino Vertex in Slice (Not Weighted)", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZ = createHistGroup("deltaZ", "#Deltaz of Neutrino Vertex in Slice", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaZDist = createHistGroup("deltaZDist", "#Deltaz of Neutrino Vertex in Slice (Not Weighted)", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaR = createHistGroup("deltaR", "#DeltaR of Neutrino Vertex in Slice", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    auto deltaRDist = createHistGroup("deltaRDist", "#DeltaR of Neutrino Vertex in Slice (Not Weighted)", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);

    for(Long64_t e = 0; e < numEntries; ++e){
        printf("=============================================================================\n");
        tree->GetEntry(e);
       
        int recoilElectron_energy = -999999;
        int recoilElectron_angle = -999999;

        // Looking at the true recoil electron
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->size() > 1) std::cout << "More than 1 true recoil electron in event!" << std::endl;
            if(truth_recoilElectronPDG->at(i) != -999999){
                recoilElectron_energy = truth_recoilElectronEnergy->at(i);
                recoilElectron_angle = truth_recoilElectronAngle->at(i);
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
        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) != -999999){
                // There is a reco slice in the event
                printf("__________________________ NEW SLICE __________________________\n");
                printf("Slice ID = %f, Category = %f, Interaction = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", reco_sliceID->at(slice), reco_sliceCategory->at(slice), reco_sliceInteraction->at(slice), reco_sliceCompleteness->at(slice), reco_slicePurity->at(slice), reco_sliceScore->at(slice));
                if(reco_sliceCategory->at(slice) != 0) printf("True Neutrino Vertex = (%f, %f, %f)\n", reco_sliceTrueVX->at(slice), reco_sliceTrueVY->at(slice), reco_sliceTrueVZ->at(slice));
        
                // Loop through all the reco neutrinos in the event          
                int PFPcounter = 0;
                int numPFPsSlice = 0;

                double summedEnergy = 0;
                double highestEnergy_PFPID = -999999;
                double highestEnergy_energy = -999999;
                double highestEnergy_theta = -999999;

                for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                    PFPcounter++;
                    if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                        // This PFP is in the slice
                        numPFPsSlice++;
                        printf("PFP %d: ID = %f, PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Energy = %f, Theta = %f, Track Score = %f, Completeness = %f, Purity = %f\n", PFPcounter, reco_particleID->at(pfp), reco_particlePDG->at(pfp), reco_particleIsPrimary->at(pfp), reco_particleVX->at(pfp), reco_particleVY->at(pfp), reco_particleVZ->at(pfp), reco_particleDX->at(pfp), reco_particleDY->at(pfp), reco_particleDZ->at(pfp), reco_particleBestPlaneEnergy->at(pfp), reco_particleTheta->at(pfp), reco_particleTrackScore->at(pfp), reco_particleCompleteness->at(pfp), reco_particlePurity->at(pfp));
                        
                        summedEnergy += reco_particleBestPlaneEnergy->at(pfp);
                        if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergy_energy){
                            highestEnergy_energy = reco_particleBestPlaneEnergy->at(pfp);
                            highestEnergy_theta = reco_particleTheta->at(pfp);
                            highestEnergy_PFPID = reco_particleID->at(pfp);
                        }
                    }
                } 

                printf("PFP with highest energy (%f) has ID %f and Theta = %f\n", highestEnergy_energy, highestEnergy_PFPID, highestEnergy_theta);
                printf("Summed energy of all PFPs in slice = %f\n", summedEnergy);
                std::cout << "Number of PFPs in the slice = " << numPFPsSlice << std::endl;

                double recoVX = -999999;
                double recoVY = -999999;
                double recoVZ = -999999;

                for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                    if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                        // Reco neutrino is in the slice
                        recoVX = reco_neutrinoVX->at(recoNeut);
                        recoVY = reco_neutrinoVX->at(recoNeut);
                        recoVZ = reco_neutrinoVX->at(recoNeut);
                        printf("Reco Neutrino in Slice: ID = %f, PDG = %f, Vertex = (%f, %f, %f)\n", reco_neutrinoID->at(recoNeut), reco_neutrinoPDG->at(recoNeut), reco_neutrinoVX->at(recoNeut), reco_neutrinoVY->at(recoNeut), reco_neutrinoVZ->at(recoNeut));
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

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.currentCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
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

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.ubooneCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.ubooneCosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.ubooneCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.ubooneCosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
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

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            ERecoSumThetaReco.nuECosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuECosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuECosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuECosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    }

                } else if(reco_sliceCategory->at(slice) == 1){
                    // This is a signal slice
                    if(DLCurrent == 2){
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
                            ERecoSumThetaReco.currentSignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                      
                            ETrueThetaReco.currentSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.currentSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.currentSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.currentSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.currentSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.currentSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));

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
                            }
                        }
                    } else if(DLCurrent == 0){
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
                      
                            ETrueThetaReco.ubooneSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.ubooneSignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.ubooneSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.ubooneSignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.ubooneSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.ubooneSignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));

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
                            }
                        }

                    } else if(DLCurrent == 5){
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
                      
                            ETrueThetaReco.nuESignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.nuESignal->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.nuESignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.nuESignal->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.nuESignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.nuESignal->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));

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
                            }
                        }
                    }
                } else if(reco_sliceCategory->at(slice) == 2){
                    // This is a fuzzy signal slice
                    if(DLCurrent == 2){
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
                      
                            ETrueThetaReco.currentSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.currentSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.currentSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.currentSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.currentSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.currentSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));

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
                            }
                        }
                    } else if(DLCurrent == 0){
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
                      
                            ETrueThetaReco.ubooneSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.ubooneSignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.ubooneSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.ubooneSignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.ubooneSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.ubooneSignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));

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
                            }
                        }
                    } else if(DLCurrent == 5){
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
                      
                            ETrueThetaReco.nuESignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ETrueThetaRecoDist.nuESignalFuzzy->Fill((recoilElectron_energy * highestEnergy_theta * highestEnergy_theta));
                            ERecoSumThetaTrue.nuESignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoSumThetaTrueDist.nuESignalFuzzy->Fill((summedEnergy * recoilElectron_angle * recoilElectron_angle));
                            ERecoHighestThetaTrue.nuESignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle), weight);
                            ERecoHighestThetaTrueDist.nuESignalFuzzy->Fill((highestEnergy_energy * recoilElectron_angle * recoilElectron_angle));

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
                            }
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
                            }
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
                            }
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
                            }
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
                            }
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
                            }
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
                            }
                        }
                    }
                }

                // aaaaa
                printf("_______________________________________________________________\n");
            }
        
        }

        

 
    }

    int drawLine = 1;
    int left = 0;
    int right = 1;

    styleDrawAll(trueETheta2, 999, 999, 999, 999, (base_path + "trueETheta2_weighted.pdf").c_str(), "bottomRight", &drawLine, &right, true, false, false);
    
    styleDrawAll(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(sliceCompletenessDist, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(slicePurityDist, 999, 999, 999, 999, (base_path + "slicePurity_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(sliceCRUMBSScore, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(sliceCRUMBSScoreDist, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(sliceNumPFPsDist, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_dist.pdf").c_str(), "topRight", nullptr, &right);

    styleDrawAll(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(ERecoSumThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right);
    styleDrawAll(ERecoHighestThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right);

    styleDrawAll(ETrueThetaReco, 999, 999, 999, 999, (base_path + "ETrueThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, false, false);
    styleDrawAll(ETrueThetaRecoDist, 999, 999, 999, 999, (base_path + "ETrueThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, false, false);
    styleDrawAll(ERecoSumThetaTrue, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, false, false);
    styleDrawAll(ERecoSumThetaTrueDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, false, false);
    styleDrawAll(ERecoHighestThetaTrue, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, false, false);
    styleDrawAll(ERecoHighestThetaTrueDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, false, false);
   
    efficiency(ERecoSumThetaRecoDist, 0, 1, 999, 999, (base_path + "ERecoSumThetaReco").c_str(), "bottomRight", nullptr, &right, -1); 
    
    //efficiency(trueETheta2, 0, 1, 999, 999, (base_path + "trueETheta2").c_str(), "bottomRight", &drawLine, &right, 1);
    
    //efficiency(sliceCompleteness, 0, 1, 999, 999, (base_path + "sliceCompleteness").c_str(), "topRight", nullptr, &right, -1);
    //efficiency(slicePurity, 0, 1, 999, 999, (base_path + "slicePurity").c_str(), "topRight", nullptr, &right, -1);
    //efficiency(sliceCRUMBSScore, 0, 1, 999, 999, (base_path + "sliceCRUMBSScore").c_str(), "topRight", nullptr, &right, -1);
    //efficiency(sliceNumPFPs, 0, 1, 999, 999, (base_path + "sliceNumPFPs").c_str(), "topRight", nullptr, &right, 1);

    //efficiency(ERecoSumThetaReco, 0, 1, 999, 999, (base_path + "ERecoSumThetaReco").c_str(), "topRight", nullptr, &right, 1);
    //efficiency(ERecoHighestThetaReco, 0, 1, 999, 999, (base_path + "ERecoHighestThetaReco").c_str(), "topRight", nullptr, &right, 1);

    //efficiency(ETrueThetaReco, 0, 1, 999, 999, (base_path + "ETrueThetaReco").c_str(), "topRight", nullptr, &right, 1);
    //efficiency(ERecoSumThetaTrue, 0, 1, 999, 999, (base_path + "ERecoSumThetaTrue").c_str(), "topRight", nullptr, &right, 1);
    //efficiency(ERecoHighestThetaTrue, 0, 1, 999, 999, (base_path + "ERecoHighestThetaTrue").c_str(), "topRight", nullptr, &right, 1);
}
