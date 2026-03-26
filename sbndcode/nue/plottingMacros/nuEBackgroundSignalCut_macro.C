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


struct eventCounter_struct{
    double nuE = 0;
    double NCNPi0 = 0;
    double otherNC = 0;
    double CCnumu = 0;
    double CCnue = 0;
    double dirt = 0;
    double nuEDirt = 0;
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
    double razzled2212Sig = 0;
    double razzled2212Back = 0;
    eventCounter_struct razzled2212IntSplit;
    double razzled13Sig = 0;
    double razzled13Back = 0;
    eventCounter_struct razzled13IntSplit;
    double razzled211Sig = 0;
    double razzled211Back = 0;
    eventCounter_struct razzled211IntSplit;
    double razzled22Sig = 0;
    double razzled22Back = 0;
    eventCounter_struct razzled22IntSplit;
    double razzled11Sig = 0;
    double razzled11Back = 0;
    eventCounter_struct razzled11IntSplit;
    double dEdxSig = 0;
    double dEdxBack = 0;
    eventCounter_struct dEdxIntSplit;
    double fracHitsContainedSig = 0;
    double fracHitsContainedBack = 0;
    eventCounter_struct fracHitsContainedIntSplit;
    double numHitsSig = 0;
    double numHitsBack = 0;
    eventCounter_struct numHitsIntSplit;
    double trackscoreSig = 0;
    double trackscoreBack = 0;
    eventCounter_struct trackscoreIntSplit;
    double ETheta2Sig = 0;    
    double ETheta2Back = 0;    
    eventCounter_struct ETheta2IntSplit;    
    double showerLengthSig = 0;    
    double showerLengthBack = 0;    
    eventCounter_struct showerLengthIntSplit;    
    double showerEnergySig = 0;    
    double showerEnergyBack = 0;    
    eventCounter_struct showerEnergyIntSplit;    
};

struct weights_struct{
    double signalNuE = 0;
    double BNBNuE = 0;
    double cosmicsNuE = 0;
    double signalCurrent = 0;
    double BNBCurrent = 0;
    double cosmicsCurrent = 0;
    double signalUboone = 0;
    double BNBUboone = 0;
    double cosmicsUboone = 0;
};

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

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* nu_e;
    TH1F* NCNpi0;
    TH1F* otherNC;
    TH1F* CCnumu;
    TH1F* CCnue;
    TH1F* dirt;
    TH1F* nu_eDirt;
    TH1F* cosmic;
    TH1F* other;
} splitHistGroup_struct;

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* electron;
    TH1F* nuEElectron;
    TH1F* proton;
    TH1F* muon;
    TH1F* cosmicMuon;
    TH1F* cosmicPhoton;
    TH1F* cosmicElectron;
    TH1F* cosmicOther;
    TH1F* pi0;
    TH1F* chargedPi;
    TH1F* photon;
    TH1F* other;
    TH1F* nuEOther;
} splitPFPHistGroup_struct;

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

splitHistGroup_struct createSplitHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);

    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_nu_e").c_str()),
        (TH1F*) base->Clone((baseName + "_NCNpi0").c_str()),
        (TH1F*) base->Clone((baseName + "_otherNC").c_str()),
        (TH1F*) base->Clone((baseName + "_CCnumu").c_str()),
        (TH1F*) base->Clone((baseName + "_CCnue").c_str()),
        (TH1F*) base->Clone((baseName + "_dirt").c_str()),
        (TH1F*) base->Clone((baseName + "_nu_eDirt").c_str()),
        (TH1F*) base->Clone((baseName + "_cosmic").c_str()),
        (TH1F*) base->Clone((baseName + "_other").c_str())
    };
}

splitPFPHistGroup_struct createSplitPFPHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);

    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of PFPs").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_electron").c_str()),
        (TH1F*) base->Clone((baseName + "_nuEElectron").c_str()),
        (TH1F*) base->Clone((baseName + "_proton").c_str()),
        (TH1F*) base->Clone((baseName + "_muon").c_str()),
        (TH1F*) base->Clone((baseName + "_cosmicMuon").c_str()),
        (TH1F*) base->Clone((baseName + "_cosmicPhoton").c_str()),
        (TH1F*) base->Clone((baseName + "_cosmicElectron").c_str()),
        (TH1F*) base->Clone((baseName + "_cosmicOther").c_str()),
        (TH1F*) base->Clone((baseName + "_pi0").c_str()),
        (TH1F*) base->Clone((baseName + "_chargedPi").c_str()),
        (TH1F*) base->Clone((baseName + "_photon").c_str()),
        (TH1F*) base->Clone((baseName + "_other").c_str()),
        (TH1F*) base->Clone((baseName + "_otherNuE").c_str())
    };
}

void fillHistogram(histGroup_struct* hist, int DLCurrent, int signal, int type, double value, weights_struct* weight){
    if(!hist || DLCurrent != 5) return;
    
    TH1* target = nullptr;

    switch(type){
        case 0: target = hist->nuECosmic; break;
        case 1: if(signal == 1) target = hist->nuESignal; break;
        case 2: if(signal == 1) target = hist->nuESignalFuzzy; break;
        case 3: target = hist->nuEBNB; break;
        case 4: target = hist->nuEBNBFuzzy; break;
    }

    if(!target) return;

    double w = 0.0;

    if(signal == 1) w = weight->signalNuE;
    else if(signal == 2) w = weight->BNBNuE;
    else if(signal == 3) w = weight->cosmicsNuE;

    target->Fill(value, w);
}

void fillSplitIntHistogram(splitHistGroup_struct* hist, int DLCurrent, int signal, int type, double value, weights_struct* weight){
    if(!hist || DLCurrent != 5) return;

    TH1* target = nullptr;

    switch(type){
        case 0: target = hist->cosmic; break;
        case 1: if(signal == 1) target = hist->nu_e; break;
        case 2: target = hist->NCNpi0; break;
        case 3: target = hist->otherNC; break;
        case 4: target = hist->CCnumu; break;
        case 5: target = hist->CCnue; break;
        case 6: target = hist->dirt; break;
        case 7: if(signal == 1) target = hist->nu_eDirt; break;
        case 8: target = hist->other; break;
        case 15: break;
    }

    if(!target) return;

    double w = 0.0;

    if(signal == 1) w = weight->signalNuE;
    else if(signal == 2) w = weight->BNBNuE;
    else if(signal == 3) w = weight->cosmicsNuE;

    target->Fill(value, w);
}

void fillSplitPFPHistogram(splitPFPHistGroup_struct* hist, int DLCurrent, int signal, int type, double value, weights_struct* weight){
    if(!hist || DLCurrent != 5) return;

    TH1* target = nullptr;

    switch(type){
        case 0: if(signal == 1) target = hist->nuEElectron; break;
        case 1: if(signal == 1) target = hist->nuEOther; break;
        case 2: target = hist->electron; break;
        case 3: target = hist->proton; break;
        case 4: target = hist->muon; break;
        case 5: target = hist->pi0; break;
        case 6: target = hist->chargedPi; break;
        case 7: target = hist->photon; break;
        case 8: target = hist->other; break;
        case 9: target = hist->cosmicMuon; break;
        case 10: target = hist->cosmicPhoton; break;
        case 11: target = hist->cosmicElectron; break;
        case 12: target = hist->cosmicOther; break;
    }
    
    if(!target) return;

    double w = 0.0;

    if(signal == 1) w = weight->signalNuE;
    else if(signal == 2) w = weight->BNBNuE;
    else if(signal == 3) w = weight->cosmicsNuE;

    target->Fill(value, w);
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

    std::vector<TH1F*> allHists = {hists.nu_e, hists.NCNpi0, hists.otherNC, hists.CCnumu, hists.CCnue, hists.dirt, hists.nu_eDirt, hists.cosmic, hists.other};

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

    hists.canvas->cd();
    hists.canvas->Clear();

    if(useLogScale) hists.canvas->SetLogy(1);
    else hists.canvas->SetLogy(0);

    // Build a frame with your desired axis limits
    double xminFrame = (xmin != 999 ? xmin : hists.nu_e->GetXaxis()->GetXmin());
    double xmaxFrame = (xmax != 999 ? xmax : hists.nu_e->GetXaxis()->GetXmax());
    double yminFrame = (ymin != 999 ? ymin : 1e-6);
    double ymaxFrame = (ymax != 999 ? ymax : stack->GetMaximum()*1.2);
    
    TH1F* frame = new TH1F("frame", stackTitle.c_str(),
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

void styleDrawPFPSplit(splitPFPHistGroup_struct hists,
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

    std::vector<TH1F*> allHists = {hists.electron, hists.nuEElectron, hists.proton, hists.muon, hists.cosmicMuon, hists.cosmicPhoton, hists.cosmicElectron, hists.cosmicOther, hists.pi0, hists.chargedPi, hists.photon, hists.other, hists.nuEOther};

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
    

    hists.nuEElectron->SetLineWidth(2);     hists.nuEElectron->SetLineColor(TColor::GetColor("#656364"));
    hists.electron->SetLineWidth(2);        hists.electron->SetLineColor(TColor::GetColor("#578dff"));
    hists.proton->SetLineWidth(2);          hists.proton->SetLineColor(TColor::GetColor("#86c8dd"));
    hists.muon->SetLineWidth(2);            hists.muon->SetLineColor(TColor::GetColor("#adad7d"));
    hists.cosmicMuon->SetLineWidth(2);      hists.cosmicMuon->SetLineColor(TColor::GetColor("#c91f16"));
    hists.cosmicPhoton->SetLineWidth(2);    hists.cosmicPhoton->SetLineColor(kViolet-4);
    hists.cosmicElectron->SetLineWidth(2);  hists.cosmicElectron->SetLineColor(kGreen+2);
    hists.cosmicOther->SetLineWidth(2);     hists.cosmicOther->SetLineColor(kPink-9);
    hists.pi0->SetLineWidth(2);             hists.pi0->SetLineColor(TColor::GetColor("#1845fb"));
    hists.chargedPi->SetLineWidth(2);       hists.chargedPi->SetLineColor(TColor::GetColor("#c849a9"));
    hists.photon->SetLineWidth(2);          hists.photon->SetLineColor(TColor::GetColor("#7a21dd"));
    hists.other->SetLineWidth(2);           hists.other->SetLineColor(TColor::GetColor("#ffa90e"));
    hists.nuEOther->SetLineWidth(2);        hists.nuEOther->SetLineColor(TColor::GetColor("#a96b59"));

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

    hists.electron->Draw("hist");
    hists.nuEElectron->Draw("histsame");
    hists.proton->Draw("histsame");
    hists.muon->Draw("histsame");
    hists.cosmicMuon->Draw("histsame");
    hists.cosmicPhoton->Draw("histsame");
    hists.cosmicElectron->Draw("histsame");
    hists.cosmicOther->Draw("histsame");
    hists.pi0->Draw("histsame");
    hists.chargedPi->Draw("histsame");
    hists.photon->Draw("histsame");
    hists.other->Draw("histsame");
    hists.nuEOther->Draw("histsame");

    int nEntries = 13;
    double height = std::max(0.03 * nEntries, 0.03);
    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;

    if(legendLocation == "topRight"){ Lxmin=0.72; Lymax=0.863; Lxmax=0.87; Lymin=Lymax - height; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.28; Lymin=Lymax - height; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.72; Lymin=0.137; Lxmax=0.87; Lymax=Lymin + height; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymin=0.137; Lxmax=0.28; Lymax=Lymin + height; }

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    //auto legend = new TLegend(0.48, 0.39, 0.87, 0.167);
    legend->AddEntry(hists.nuEElectron, "#nu+e Electron", "f");
    legend->AddEntry(hists.electron, "Electron", "f");
    legend->AddEntry(hists.proton, "Proton", "f");
    legend->AddEntry(hists.muon, "Muon", "f");
    legend->AddEntry(hists.cosmicMuon, "Cosmic Muon", "f");
    legend->AddEntry(hists.cosmicPhoton, "Cosmic Photon", "f");
    legend->AddEntry(hists.cosmicElectron, "Cosmic Electron", "f");
    legend->AddEntry(hists.cosmicOther, "Cosmic Other", "f");
    legend->AddEntry(hists.pi0, "Neutral Pion", "f");
    legend->AddEntry(hists.chargedPi, "Charged Pion", "f");
    legend->AddEntry(hists.photon, "Photon", "f");
    legend->AddEntry(hists.other, "Other", "f");
    legend->AddEntry(hists.nuEOther, "#nu+e Other", "f");
    legend->SetTextSize(0.0225);

    legend->SetMargin(0.13);
    legend->Draw();

    hists.canvas->SaveAs(filename);

    // Drawing the histograms as a stack.
    const char* histsTitle = hists.electron->GetTitle();
    const char* xAxisTitle = hists.electron->GetXaxis()->GetTitle();
    const char* yAxisTitle = hists.electron->GetYaxis()->GetTitle();

    std::string stackTitle = std::string(histsTitle) + ";" + xAxisTitle + ";" + yAxisTitle;

    THStack* stack = new THStack("stack", stackTitle.c_str());
    stack->Add(hists.electron);
    stack->Add(hists.nuEElectron);
    stack->Add(hists.proton);
    stack->Add(hists.muon);
    stack->Add(hists.cosmicMuon);
    stack->Add(hists.cosmicPhoton);
    stack->Add(hists.cosmicElectron);
    stack->Add(hists.cosmicOther);
    stack->Add(hists.pi0);
    stack->Add(hists.chargedPi);
    stack->Add(hists.photon);
    stack->Add(hists.other);
    stack->Add(hists.nuEOther);

    hists.canvas->cd();
    hists.canvas->Clear();

    if(useLogScale) hists.canvas->SetLogy(1);
    else hists.canvas->SetLogy(0);

    // Build a frame with your desired axis limits
    double xminFrame = (xmin != 999 ? xmin : hists.electron->GetXaxis()->GetXmin());
    double xmaxFrame = (xmax != 999 ? xmax : hists.electron->GetXaxis()->GetXmax());
    double yminFrame = (ymin != 999 ? ymin : 1e-6);
    double ymaxFrame = (ymax != 999 ? ymax : stack->GetMaximum()*1.2);
    
    TH1F* frame = new TH1F("frame", stackTitle.c_str(),
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
    legendStack->AddEntry(hists.nuEElectron, "#nu+e Electron", "f");
    legendStack->AddEntry(hists.electron, "Electron", "f");
    legendStack->AddEntry(hists.proton, "Proton", "f");
    legendStack->AddEntry(hists.muon, "Muon", "f");
    legendStack->AddEntry(hists.cosmicMuon, "Cosmic Muon", "f");
    legendStack->AddEntry(hists.cosmicPhoton, "Cosmic Photon", "f");
    legendStack->AddEntry(hists.cosmicElectron, "Cosmic Electron", "f");
    legendStack->AddEntry(hists.cosmicOther, "Cosmic Other", "f");
    legendStack->AddEntry(hists.pi0, "Neutral Pion", "f");
    legendStack->AddEntry(hists.chargedPi, "Charged Pion", "f");
    legendStack->AddEntry(hists.photon, "Photon", "f");
    legendStack->AddEntry(hists.other, "Other", "f");
    legendStack->AddEntry(hists.nuEOther, "#nu+e Other", "f");
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

    auto combine = [useLogScale](TH1F* a, TH1F* b, TH1F* c, TH1F* d, const char* name) -> TH1F* {
        TH1F* combo = nullptr;
        if (a) combo = (TH1F*)a->Clone(name);
        else if (b) combo = (TH1F*)b->Clone(name);
        else if (c) combo = (TH1F*)c->Clone(name);
        else if (d) combo = (TH1F*)d->Clone(name);
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

    TH1F* currentSignalCombo     = combine(hists.currentSignal, nullptr, nullptr, nullptr, "currentSignalCombo");
    TH1F* currentBackgroundCombo = combine(hists.currentBNB, hists.currentBNBFuzzy, hists.currentCosmic, hists.currentSignalFuzzy, "currentBackgroundCombo");

    TH1F* ubooneSignalCombo     = combine(hists.ubooneSignal, nullptr, nullptr, nullptr, "ubooneSignalCombo");
    TH1F* ubooneBackgroundCombo = combine(hists.ubooneBNB, hists.ubooneBNBFuzzy, hists.ubooneCosmic, hists.ubooneSignalFuzzy, "ubooneBackgroundCombo");

    TH1F* nuESignalCombo     = combine(hists.nuESignal, nullptr, nullptr, nullptr, "nuESignalCombo");
    TH1F* nuEBackgroundCombo = combine(hists.nuEBNB, hists.nuEBNBFuzzy, hists.nuECosmic, hists.nuESignalFuzzy, "nuEBackgroundCombo");

    std::vector<TH1F*> allHists = {
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

TH1F* makeTotalHist(const TH1F* h){
    TH1F* hc = (TH1F*)h->Clone(Form("%s_totalSum", h->GetName()));
    hc->Reset();

    double totalSum = h->Integral(0, h->GetNbinsX() + 1);
    //std::cout << "total sum = " << totalSum << std::endl; 
    
    for(int i = 0; i <= hc->GetNbinsX() + 1; ++i){
        hc->SetBinContent(i, totalSum);
    }

    return hc;
}

TH1F* makeCumulative(const TH1F* h, bool keepRight){
    TH1F* hc = (TH1F*)h->Clone(Form("%s_cumulative", h->GetName()));
    hc->Reset();

    int n = h->GetNbinsX();

    if (keepRight) {
        double sum = 0.0;
        for (int i = n; i >= 1; --i) {
            sum += h->GetBinContent(i);
            hc->SetBinContent(i, sum);
            //std::cout << "Bin " << i << ": bin content = " << h->GetBinContent(i) << ", sum = " << sum << std::endl;
        }
    } else {
        double sum = 0.0;
        for (int i = 1; i <= n; ++i) {
            sum += h->GetBinContent(i);
            hc->SetBinContent(i, sum);
            //std::cout << "Bin " << i << ": bin content = " << h->GetBinContent(i) << ", sum = " << sum << std::endl;
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

void drawEfficiencyErrors(TEfficiency* plot, const std::string& filename, double lowY, double highY, const std::string& legendLocation, double xmin, double xmax, double efficiencyWay = 0.0, bool writeMaxValues = false, const std::string& text_filename = ""){
    if (!plot) {
        std::cerr << "drawEfficiency: null TEfficiency pointer\n";
        return;
    }

    double maxVal = getMaxValueEfficiency(plot, false);
    double minVal = getMinValueEfficiency(plot, xmin, xmax, false);
    //std::cout << "minVal = " << minVal << ", maxVal = " << maxVal << std::endl;

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
 
        if(yEff > maxEff){
            maxEff = yEff;
            if(efficiencyWay == 1) maxEffBin = (xCenter + xErr);
            if(efficiencyWay == -1) maxEffBin = (xCenter - xErr);
        }

        //std::cout << "Bin " << i << ": yEff = " << yEff << ", yErrLow = " << yErrLow << ", yErrUp = " << yErrUp << ", xErr = " << xErr << ", xCenter = " << xCenter << std::endl; 
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

    if(writeMaxValues){
        std::ofstream outfile(text_filename, std::ios::app);
        if(outfile.is_open()){
            outfile << "=======================" << std::endl;
            outfile << filename << std::endl;

            outfile << "DL Nu+E Vertexing: Max Eff x Pur = " << maxEff << ", Bin Value = " << maxEffBin << std::endl;
            outfile.close();
        } else{
            std::cerr << "Error: couldn't open " << text_filename << " for writing." << std::endl;
        }
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
    delete c;
}

void drawTEff(TH1F* numerator, TH1F* denominator, double lowY, double highY, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr){
    TEfficiency* efficiency = new TEfficiency(*numerator, *denominator);
    
    efficiency->SetTitle(Form("%s;%s;Efficiency", numerator->GetTitle(), numerator->GetXaxis()->GetTitle()));
    efficiency->SetStatisticOption(TEfficiency::kFNormal);

    efficiency->SetUseWeightedEvents(false);

    double maxVal = getMaxValueEfficiency(efficiency, false);
    double minVal = getMinValueEfficiency(efficiency, xmin, xmax, false);

    TCanvas* c = new TCanvas("c_eff", "Efficiency comparison", 800, 600);
    c->SetTicks();
    c->SetLeftMargin(0.15);

    efficiency->SetMarkerColor(kBlack);
    efficiency->SetMarkerSize(0.7); 
    efficiency->SetLineWidth(1);
    efficiency->SetLineColor(kBlack);
    efficiency->SetMarkerStyle(20);

    int nBins = numerator->GetNbinsX();
    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(nBins);

    for(int i = 1; i <= nBins; ++i){
        double xCenter = numerator->GetXaxis()->GetBinCenter(i);
        double xErr = (numerator->GetXaxis()->GetBinUpEdge(i) - numerator->GetXaxis()->GetBinLowEdge(i))/ 2.0;
        
        double yEff = efficiency->GetEfficiency(i);
        double yErrLow = efficiency->GetEfficiencyErrorLow(i);
        double yErrUp = efficiency->GetEfficiencyErrorUp(i);

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

    gEff->SetTitle(numerator->GetTitle());
    gEff->GetXaxis()->SetTitle(numerator->GetXaxis()->GetTitle());
    gEff->GetYaxis()->SetTitle(numerator->GetYaxis()->GetTitle());
    gEff->GetYaxis()->SetTitleOffset(1.6);
    gEff->Draw("AP");

    efficiency->Draw("SAME");
    gPad->Update();

    auto gBDT = efficiency->GetPaintedGraph();
    gBDT->SetMarkerSize(0.8);
    gBDT->Draw("PE SAME");

    auto* g = efficiency->GetPaintedGraph();

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
    leg->AddEntry(efficiency, "DL Nu+E", "LEP");

    //leg->Draw();

    c->SaveAs(filename);
    delete c;
}

void efficiency(histGroup_struct* histBeforeCuts, histGroup_struct* histAfterCuts, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr, double efficiencyWay = 0.0, const std::string& text_filename = ""){
    bool keepRight = (efficiencyWay == -1);

    // Total signal before cuts
    //std::cout << "============================ total signal before cuts ============================" << std::endl;
    TH1F* hTotalSignalBeforeCuts = makeTotalHist(histBeforeCuts->nuESignal);
    //std::cout << "======================================================================================" << std::endl;

    // Total signal kept after cuts (cumulative)
    //std::cout << "============================ total signal kept after cuts ============================" << std::endl;
    TH1F* hPassedSignalAfterCuts = makeCumulative(histAfterCuts->nuESignal, keepRight);
    //std::cout << "======================================================================================" << std::endl;

    // Total background before cuts
    //std::cout << "============================ total background before cuts ============================" << std::endl;
    TH1F* hTotalBackgroundBeforeCutsAdded = (TH1F*) histBeforeCuts->nuECosmic->Clone("hTotalBackgroundBeforeCutsAdded");
    hTotalBackgroundBeforeCutsAdded->Reset();
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuECosmic);
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuEBNB);
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuEBNBFuzzy);
    hTotalBackgroundBeforeCutsAdded->Add(histBeforeCuts->nuESignalFuzzy);
    //std::cout << "Number of entries = " << hTotalBackgroundBeforeCutsAdded->GetEntries() << std::endl;
    TH1F* hTotalBackgroundBeforeCuts = makeTotalHist(hTotalBackgroundBeforeCutsAdded);
    //std::cout << "======================================================================================" << std::endl;

    // Total background kept after cuts (cumulative)
    //std::cout << "============================ total background kept after cuts ============================" << std::endl;
    TH1F* hPassedBackgroundAfterCutsAdded = (TH1F*) histAfterCuts->nuECosmic->Clone("hPassedBackgroundAfterCutsAdded");
    hPassedBackgroundAfterCutsAdded->Reset();
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuECosmic);
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuEBNB);
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuEBNBFuzzy);
    hPassedBackgroundAfterCutsAdded->Add(histAfterCuts->nuESignalFuzzy);
    TH1F* hPassedBackgroundAfterCuts = makeCumulative(hPassedBackgroundAfterCutsAdded, keepRight);
    //std::cout << "======================================================================================" << std::endl;

    // Total background rejected after cuts (cumulative)
    //std::cout << "============================ total background rejected after cuts ============================" << std::endl;
    TH1F* hRejectedBackgroundAfterCuts = (TH1F*) hTotalBackgroundBeforeCuts->Clone("hRejectedBackgroundAfterCuts");
    hRejectedBackgroundAfterCuts->Add(hPassedBackgroundAfterCuts, -1.0);
    //std::cout << "======================================================================================" << std::endl;

    // Total background + signal kept after cuts (cumulative)
    //std::cout << "============================ total background + signal after cuts ============================" << std::endl;
    TH1F* hPassedBackgroundSignalAfterCutsAdded = (TH1F*) histAfterCuts->nuECosmic->Clone("hPassedBackgroundSignalAfterCutsAdded");
    hPassedBackgroundSignalAfterCutsAdded->Reset();
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuECosmic);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuEBNB);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuEBNBFuzzy);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuESignalFuzzy);
    hPassedBackgroundSignalAfterCutsAdded->Add(histAfterCuts->nuESignal);
    TH1F* hPassedBackgroundSignalAfterCuts = makeCumulative(hPassedBackgroundSignalAfterCutsAdded, keepRight);
    //std::cout << "======================================================================================" << std::endl;

    TH1F* hEffPurNumerator = (TH1F*) histAfterCuts->nuESignal->Clone("hEffPurNumerator");
    hEffPurNumerator->Reset();
    hEffPurNumerator->Add(hPassedSignalAfterCuts);
    hEffPurNumerator->Multiply(hPassedSignalAfterCuts);

    TH1F* hEffPurDenominator = (TH1F*) histAfterCuts->nuESignal->Clone("hEffPurDenominator");
    hEffPurDenominator->Reset();
    hEffPurDenominator->Add(hTotalSignalBeforeCuts);
    hEffPurDenominator->Multiply(hPassedBackgroundSignalAfterCuts);

    // Efficiency plot
    TEfficiency* eff = new TEfficiency(*hPassedSignalAfterCuts, *hTotalSignalBeforeCuts);
    TEfficiency* rej = new TEfficiency(*hRejectedBackgroundAfterCuts, *hTotalBackgroundBeforeCuts);
    TEfficiency* pur = new TEfficiency(*hPassedSignalAfterCuts, *hPassedBackgroundSignalAfterCuts);
    TEfficiency* effPur = new TEfficiency(*hEffPurNumerator, *hEffPurDenominator);

    eff->SetTitle(Form("%s;%s;Efficiency", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    eff->SetStatisticOption(TEfficiency::kFNormal);
    
    pur->SetTitle(Form("%s;%s;Purity", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    pur->SetStatisticOption(TEfficiency::kFNormal);
    
    rej->SetTitle(Form("%s;%s;Background Rejection", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    rej->SetStatisticOption(TEfficiency::kFNormal);
    
    effPur->SetTitle(Form("%s;%s;Efficiency x Purity", histAfterCuts->nuESignal->GetTitle(), histAfterCuts->nuESignal->GetXaxis()->GetTitle()));
    effPur->SetStatisticOption(TEfficiency::kFNormal);

    std::string filenameEff = std::string(filename) + "_eff.pdf";
    std::string filenamePur = std::string(filename) + "_pur.pdf";
    std::string filenameRej = std::string(filename) + "_rej.pdf";
    std::string filenameEffPur = std::string(filename) + "_effPur.pdf";

    eff->SetUseWeightedEvents(false);
    pur->SetUseWeightedEvents(false);
    rej->SetUseWeightedEvents(false);
    effPur->SetUseWeightedEvents(false);
    
    drawEfficiencyErrors(eff, filenameEff, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay, false, text_filename);
    drawEfficiencyErrors(pur, filenamePur, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay, false, text_filename);
    drawEfficiencyErrors(rej, filenameRej, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay, false, text_filename);
    drawEfficiencyErrors(effPur, filenameEffPur, -999999, -999999, legendLocation, xmin, xmax, efficiencyWay, true, text_filename);
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

    TProfile* profX = hist->ProfileX("_pfx", 1, -1, "s");

    TCanvas* ProfileCanvas = new TCanvas("profile_canvas", "TProfile from TH2D", 300, 50, 800, 600);
    ProfileCanvas->SetTickx();
    ProfileCanvas->SetTicky();

    profX->SetTitle(Form("ProfileX (Standard Deviation) of %s", title));
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
    if(dotPos != std::string::npos){
        profileFilename.insert(dotPos, "_profileSD");
    } else{
        profileFilename += "_profileSD.pdf"; 
    }

    ProfileCanvas->SaveAs(profileFilename.c_str());

    TProfile* profX_errMean = hist->ProfileX("_pfx_errMean", 1, -1, "");
    TCanvas* ProfileErrCanvas = new TCanvas("profileErr_canvas", "TProfile from TH2D (Error on Mean)", 300, 50, 800, 600);
    ProfileErrCanvas->SetTickx();
    ProfileErrCanvas->SetTicky();

    profX_errMean->SetTitle(Form("ProfileX (Error on Mean) of %s", title));
    profX_errMean->SetLineColor(kBlack);
    profX_errMean->SetLineWidth(2);
    profX_errMean->SetMarkerStyle(20);
    profX_errMean->SetMarkerSize(0.8);
    profX_errMean->SetMarkerColor(kBlack);
    profX_errMean->Draw("E1");

    profX_errMean->GetXaxis()->SetTickLength(0.04);
    profX_errMean->GetYaxis()->SetTickLength(0.03);
    profX_errMean->GetXaxis()->SetTickSize(0.02);
    profX_errMean->GetYaxis()->SetTickSize(0.02);
    profX_errMean->SetStats(0);

    std::string profileErrFilename = std::string(filename);
    size_t dotPosErr = profileErrFilename.find_last_of(".");
    if(dotPosErr != std::string::npos){
            profileErrFilename.insert(dotPosErr, "_profileErrorMean");
    } else{
            profileErrFilename += "_profileErrorMean.pdf";
    }

    ProfileErrCanvas->SaveAs(profileErrFilename.c_str());

    TwoDHistCanvas->Clear();
    ProfileCanvas->Clear();
}

void nuEBackgroundSignalCut_macro(){
    //std::string txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_harsherCrumbs_fv_primaryPFP_ETheta2_razzled2212_razzled13_razzled211_razzled22_razzled11_dEdx_containedHits_numHits.txt";
    std::string txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_CRUMBS.txt";
    //std::string txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_harsherCrumbs_fv_primaryPFP_ETheta2.txt";

    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_IntimeBNBNuE_DLNuE_20Feb.root"); 
    //TFile *file = TFile::Open("/exp/sbnd/app/users/coackley/nue/srcs/sbndcode/sbndcode/nue/plottingMacros/merged.root"); 
    //std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsFix_clearCosmic_numPFPs0_recoNeut_harsherCrumbs_fv_primaryPFP_ETheta2_razzled2212_razzled13_razzled211_razzled22_razzled11_dEdx_containedHits_numHits/";
    std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsFix_clearCosmic_numPFPs0_recoNeut_CRUMBS/";
    //std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsFix_clearCosmic_numPFPs0_recoNeut_harsherCrumbs_fv_primaryPFP_ETheta2/";

    gROOT->SetBatch(true);

    int clearCosmicCut = 1;
    int numPFPs0Cut = 1;
    int numRecoNeutrinosCut = 1;
    int CRUMBSCut = 1;
    int FVCut = 0;
    int primaryPFPCut = 0;
    int ETheta2Cut = 0;
    int razzledPDG2212Cut = 0;
    int razzledPDG13Cut = 0;
    int razzledPDG211Cut = 0;
    int razzledPDG22Cut = 0;
    int razzledPDG11Cut = 0;
    int dEdxCut = 0;
    int fracHitsContainedCut = 0;
    int numHitsCut = 0;

    int printLowRazzledMuonMuons = 0;
    int printHighRazzledMuonMuons = 0;

    int printLowRazzledElectronElectrons = 0;
    int printHighRazzledElectronElectrons = 0;

    int printLowRazzledPhotonPhotons = 0;
    int printHighRazzledPhotonPhotons = 0;
    int printHighRazzledElectronPhotons = 0;

    int printLowRazzledChargedPiChargedPis = 0;
    int printHighRazzledChargedPiChargedPis = 0;

    int crumbsPrint = 0;
    int FVCutPrint = 0;
    int primaryPFPPrint = 0;
    int ETheta2Print = 0;
    int razzledProtonPrint = 0;
    int razzledMuonPrint = 0;
    int razzledChargedPiPrint = 0;
    int razzledPhotonPrint = 0;
    int razzledElectronPrint = 0;
    int dEdxPrint = 0;
    int fracHitsPrint = 0;
    int numHitsPrint = 0;

    double numSignalSlicesPFPCountMoreThan1 = 0;
    double numSignalSlicesPFPCountMoreThan1NuE = 0;
    double numSignalSlicesPFPCountMoreThan1NuEPrimaryPlus = 0;
    double numSignalSlicesPFPCount = 0; 

    // Cut values
    double crumbsScoreCut_low = 0.12;
    double crumbsScoreCut_high = 0.6;

    double FVCut_xHigh = 195; 
    double FVCut_xLow = -197; 
    double FVCut_xCentre = 10; 

    double FVCut_yHigh = 194; 
    double FVCut_yLow = -196; 
    
    double FVCut_zHigh = 450; 
    double FVCut_zLow = 6.5; 
   
    double primaryPFPCutValue = 1;

    double razzled2212High_highestEnergyPFP = 0.05;
    double razzled2212Low_highestEnergyPFP = 0;
    double razzled13High_highestEnergyPFP = 0.05;
    double razzled13Low_highestEnergyPFP = 0;
    double razzled211High_highestEnergyPFP = 0.05;
    double razzled211Low_highestEnergyPFP = 0;
    double razzled22High_highestEnergyPFP = 0.15;
    double razzled22Low_highestEnergyPFP = 0;
    double razzled11High_highestEnergyPFP = 1;
    double razzled11Low_highestEnergyPFP = 0.85;
    
    double dEdxHigh_highestEnergyPFP = 4;
    double dEdxLow_highestEnergyPFP = 0.25;

    double numContainedHits_high = 1;
    double numContainedHits_low = 0.7; 
    
    double numHits_high = 2990;
    double numHits_low = 400; 

    double ETheta2High_highestEnergyPFP = 2.044;
    double ETheta2Low_highestEnergyPFP = 0;


    // If the directory already exists, delete everything in it
    // If the directory doesn't exists, create it. 
    if (!gSystem->AccessPathName(base_path.c_str())) {
        gSystem->Exec(Form("rm -rf %s/*", base_path.c_str()));
    }
    gSystem->mkdir(base_path.c_str(), kTRUE);

    std::string tableFileName = base_path + "table.txt";
    
    std::ofstream clearFile(txtFileName, std::ios::trunc);
    if (!clearFile.is_open()) {
        std::cerr << "Error: could not open or create " << txtFileName << std::endl;
        return;
    }
    clearFile.close();
    
    std::ofstream clearTableFile(tableFileName, std::ios::trunc);
    if (!clearTableFile.is_open()) {
        std::cerr << "Error: could not open or create " << tableFileName << std::endl;
        return;
    }
    clearTableFile.close();

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

    beforeEventCount_struct eventsBeforeCuts_DLNuE;
    eventCounting_struct eventsAfterCuts_DLNuE;

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
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBUboone;

    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;

    double cosmicSpillsSumNuE = 0;
    double BNBSpillsSumNuE = 0;
    double NuESpillsSumNuE = 0;

    double POTSignalNuE_notMissing = 0;
    double POTBNBNuE_notMissing = 0;
    
    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 2 && subRunDLCurrent == 5) BNBSpillsSumNuE += subRunNumGenEvents;
        else if(subRunSignal == 1 && subRunDLCurrent == 5) NuESpillsSumNuE += subRunNumGenEvents;

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
        }
    }

    double targetPOT = 1e21;
    double targetSpills = (targetPOT/(5e12));

    double BNBScaledSpills_NuE = ((targetPOT/POTBNBNuE_notMissing) * BNBSpillsSumNuE);
    double SignalScaledSpills_NuE = ((targetPOT/POTSignalNuE_notMissing) * NuESpillsSumNuE);

    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    weights_struct weights;
    weights.signalNuE = targetPOT / POTSignalNuE_notMissing;
    weights.BNBNuE = targetPOT /POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;

    std::cout << "Weights DLNu+E: BNB = " << weights.BNBNuE << ", Signal = " << weights.signalNuE << ", Intime Cosmics = " << weights.cosmicsNuE << std::endl;

    UInt_t eventID, runID, subRunID;
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

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
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

    Long64_t numEntries = tree->GetEntries();

    auto sliceCompletenessBeforeCuts = createHistGroup("sliceCompletenessBeforeCuts", "Slice Completeness (Before Cuts)", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessAfterCuts = createHistGroup("sliceCompletenessAfterCuts", "Slice Completeness (After Cuts)", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessAfterCuts_splitDLNuE = createSplitHistGroup("sliceCompleteness_splitDLNuE", "Slice Completeness: DL Nu+E Vertexing", "Completeness", 102, 0, 1.02);
    auto sliceCompletenessAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceCompleteness_splitPFPDLNuE", "Slice Completeness: DL Nu+E Vertexing", "Completeness", 102, 0, 1.02);    
    
    auto sliceCRUMBSBeforeCuts = createHistGroup("sliceCRUMBSBeforeCuts", "Slice CRUMBS Score (Before Cuts)", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCuts = createHistGroup("sliceCRUMBSAfterCuts", "Slice CRUMBS Score (After Cuts)", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCuts_splitDLNuE = createSplitHistGroup("sliceCRUMBS_splitDLNuE", "Slice CRUMBS Score: DL Nu+E Vertexing", "CRUMBS Score", 25, -1, 1);
    auto sliceCRUMBSAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceCRUMBS_splitPFPDLNuE", "Slice CRUMBS Score: DL Nu+E Vertexing", "CRUMBS Score", 25, -1, 1);    
    
    auto slicePurityBeforeCuts = createHistGroup("slicePurityBeforeCuts", "Slice Purity (Before Cuts)", "Purity", 102, 0, 1);
    auto slicePurityAfterCuts = createHistGroup("slicePurityAfterCuts", "Slice Purity (After Cuts)", "Purity", 102, 0, 1);
    auto slicePurityAfterCuts_splitDLNuE = createSplitHistGroup("slicePurity_splitDLNuE", "Slice Purity: DL Nu+E Vertexing", "CRUMBS Score", 102, 0, 1);
    auto slicePurityAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("slicePurity_splitPFPDLNuE", "Slice Purity: DL Nu+E Vertexing", "CRUMBS Score", 102, 0, 1);    

    auto sliceNumRecoNeutBeforeCuts = createHistGroup("sliceNumRecoNeutBeforeCuts", "Number of Reco Neutrinos in Slice (Before Cuts)", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterCuts = createHistGroup("sliceNumRecoNeutAfterCuts", "Number of Reco Neutrinos in Slice (After Cuts)", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumRecoNeut_splitDLNuE", "Number of Reco Neutrinos in Slice: DL Nu+E Vertexing", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumRecoNeutAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceNumRecoNeut_splitPFPDLNuE", "Number of Reco Neutrinos in Slice: DL Nu+E Vertexing", "Number of Reco Neutrinos", 10, 0, 10);    

    auto sliceNumPFPsBeforeCuts = createHistGroup("sliceNumPFPsBeforeCuts", "Number of PFPs in Slice (Before Cuts)", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts = createHistGroup("sliceNumPFPsAfterCuts", "Number of PFPs in Slice (After Cuts)", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPFPs_splitDLNuE", "Number of PFPs in Slice: DL Nu+E Vertexing", "Number of PFPs", 20, 0, 20);
    auto sliceNumPFPsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceNumPFPs_splitPFPDLNuE", "Number of PFPs in Slice: DL Nu+E Vertexing", "Number of PFPs", 20, 0, 20);    
    
    auto sliceNumPrimaryPFPsBeforeCuts = createHistGroup("sliceNumPrimaryPFPsBeforeCuts", "Number of Primary PFPs in Slice (Before Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterCuts = createHistGroup("sliceNumPrimaryPFPsAfterCuts", "Number of Primary PFPs in Slice (After Cuts)", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPs_splitDLNuE", "Number of Primary PFPs in Slice: DL Nu+E Vertexing", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceNumPrimaryPFPs_splitPFPDLNuE", "Number of Primary PFPs in Slice: DL Nu+E Vertexing", "Number of Primary PFPs", 20, 0, 20);    
    
    auto sliceFracHitsInPFPsBeforeCuts = createHistGroup("sliceFracHitsInPFPsBeforeCuts", "Fraction of Hits in Slice Contained in PFPs (Before Cuts)", "Fraction", 20, 0, 1);
    auto sliceFracHitsInPFPsAfterCuts = createHistGroup("sliceFracHitsInPFPsAfterCuts", "Fraction of Hits in Slice Contained in PFPs (After Cuts)", "Fraction", 20, 0, 1);
    auto sliceFracHitsInPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceFracHitsInPFPs_splitDLNuE", "Fraction of Hits in Slice Contained in PFPs: DL Nu+E Vertexing", "Fraction", 20, 0, 1);
    auto sliceFracHitsInPFPsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceFracHitsInPFPs_splitPFPDLNuE", "Fraction of Hits in Slice Contained in PFPs: DL Nu+E Vertexing", "Fraction", 20, 0, 1);    
    
    auto sliceFracHitsInHighestEnergyPFPsBeforeCuts = createHistGroup("sliceFracHitsInHighestEnergyPFPsBeforeCuts", "Fraction of Hits in Slice Contained in Highest Energy PFP in Slice (Before Cuts)", "Fraction", 20, 0, 1);
    auto sliceFracHitsInHighestEnergyPFPsAfterCuts = createHistGroup("sliceFracHitsInHighestEnergyPFPsAfterCuts", "Fraction of Hits in Slice Contained in Highest Energy PFP in Slice (After Cuts)", "Fraction", 20, 0, 1);
    auto sliceFracHitsInHighestEnergyPFPsAfterCuts_splitDLNuE = createSplitHistGroup("sliceFracHitsInHighestEnergyPFPs_splitDLNuE", "Fraction of Hits in Slice Contained in Highest Energy PFP in Slice: DL Nu+E Vertexing", "Fraction", 20, 0, 1);
    auto sliceFracHitsInHighestEnergyPFPsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceFracHitsInHighestEnergyPFPs_splitPFPDLNuE", "Fraction of Hits in Slice Contained in Highest Energy PFP in Slice: DL Nu+E Vertexing", "Fraction", 20, 0, 1);    
    
    auto ERecoSumThetaRecoBeforeCuts = createHistGroup("ERecoSumThetaRecoBeforeCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoAfterCuts = createHistGroup("ERecoSumThetaRecoAfterCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoAfterCuts_splitDLNuE = createSplitHistGroup("ERecoSumThetaReco_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice: DL Nu+E Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("ERecoSumThetaReco_splitPFPDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice: DL Nu+E Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);    
    
    auto ERecoHighestThetaRecoBeforeCuts = createHistGroup("ERecoHighestThetaRecoBeforeCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Before Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts = createHistGroup("ERecoHighestThetaRecoAfterCuts", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (After Cuts)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts_splitDLNuE = createSplitHistGroup("ERecoHighestThetaReco_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice: DL Nu+E Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("ERecoHighestThetaReco_splitPFPDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice: DL Nu+E Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);    
    
    auto dEdxBeforeCuts = createHistGroup("dEdxBeforeCuts", "dE/dx of the PFP in the Slice with the Highest Energy (Before Cuts)", "dE/dx", 40, 0, 10);
    auto dEdxAfterCuts = createHistGroup("dEdxAfterCuts", "dE/dx of the PFP in the Slice with the Highest Energy (After Cuts)", "dE/dx", 40, 0, 10);
    auto dEdxAfterCuts_splitDLNuE = createSplitHistGroup("dEdx_splitDLNuE", "dE/dx of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "dE/dx", 40, 0, 10);
    auto dEdxAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("dEdx_splitPFPDLNuE", "dE/dx of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "dE/dx", 40, 0, 10);    
    
    auto razzledPDG11BeforeCuts = createHistGroup("razzledPDG11BeforeCuts", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 20, 0, 1);
    auto razzledPDG11AfterCuts = createHistGroup("razzledPDG11AfterCuts", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 20, 0, 1);
    auto razzledPDG11AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG11_splitDLNuE", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);
    auto razzledPDG11AfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("razzledPDG11_splitPFPDLNuE", "Razzled PDG 11 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);    
    
    auto razzledPDG13BeforeCuts = createHistGroup("razzledPDG13BeforeCuts", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 20, 0, 1);
    auto razzledPDG13AfterCuts = createHistGroup("razzledPDG13AfterCuts", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 20, 0, 1);
    auto razzledPDG13AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG13_splitDLNuE", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);
    auto razzledPDG13AfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("razzledPDG13_splitPFPDLNuE", "Razzled PDG 13 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);    
    
    auto razzledPDG22BeforeCuts = createHistGroup("razzledPDG22BeforeCuts", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 20, 0, 1);
    auto razzledPDG22AfterCuts = createHistGroup("razzledPDG22AfterCuts", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 20, 0, 1);
    auto razzledPDG22AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG22_splitDLNuE", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);
    auto razzledPDG22AfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("razzledPDG22_splitPFPDLNuE", "Razzled PDG 22 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);    
    
    auto razzledPDG211BeforeCuts = createHistGroup("razzledPDG211BeforeCuts", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 20, 0, 1);
    auto razzledPDG211AfterCuts = createHistGroup("razzledPDG211AfterCuts", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 20, 0, 1);
    auto razzledPDG211AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG211_splitDLNuE", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);
    auto razzledPDG211AfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("razzledPDG211_splitPFPDLNuE", "Razzled PDG 211 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);    
    
    auto razzledPDG2212BeforeCuts = createHistGroup("razzledPDG2212BeforeCuts", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy (Before Cuts)", "Score", 20, 0, 1);
    auto razzledPDG2212AfterCuts = createHistGroup("razzledPDG2212AfterCuts", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy (After Cuts)", "Score", 20, 0, 1);
    auto razzledPDG2212AfterCuts_splitDLNuE = createSplitHistGroup("razzledPDG2212_splitDLNuE", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);
    auto razzledPDG2212AfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("razzledPDG2212_splitPFPDLNuE", "Razzled PDG 2212 Score of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Score", 20, 0, 1);    

    auto pfpCompletenessBeforeCuts = createHistGroup("pfpCompletenessBeforeCuts", "Completeness of the PFP in the Slice with the Highest Energy (Before Cuts)", "Completeness", 50, 0, 1);
    auto pfpCompletenessAfterCuts = createHistGroup("pfpCompletenessAfterCuts", "Completeness of the PFP in the Slice with the Highest Energy (After Cuts)", "Completeness", 50, 0, 1);
    auto pfpCompletenessAfterCuts_splitDLNuE = createSplitHistGroup("pfpCompleteness_splitDLNuE", "Completeness of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Completeness", 50, 0, 1);
    auto pfpCompletenessAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("pfpCompleteness_splitPFPDLNuE", "Completeness of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Completeness", 50, 0, 1);    

    auto pfpPurityBeforeCuts = createHistGroup("pfpPurityBeforeCuts", "Purity of the PFP in the Slice with the Highest Energy (Before Cuts)", "Purity", 50, 0, 1);
    auto pfpPurityAfterCuts = createHistGroup("pfpPurityAfterCuts", "Purity of the PFP in the Slice with the Highest Energy (After Cuts)", "Purity", 50, 0, 1);
    auto pfpPurityAfterCuts_splitDLNuE = createSplitHistGroup("pfpPurity_splitDLNuE", "Purity of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Purity", 50, 0, 1);
    auto pfpPurityAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("pfpPurity_splitPFPDLNuE", "Purity of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Purity", 50, 0, 1);    
    
    auto pfpNumHitsBeforeCuts = createHistGroup("pfpNumHitsBeforeCuts", "Number of Hits in the PFP in the Slice with the Highest Energy (Before Cuts)", "Number of Hits", 300, 0, 3000);
    auto pfpNumHitsAfterCuts = createHistGroup("pfpNumHitsAfterCuts", "Number of Hits in the PFP in the Slice with the Highest Energy (After Cuts)", "Number of Hits", 300, 0, 3000);
    auto pfpNumHitsAfterCuts_splitDLNuE = createSplitHistGroup("pfpNumHits_splitDLNuE", "Number of Hits in the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Number of Hits", 300, 0, 3000);
    auto pfpNumHitsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("pfpNumHits_splitPFPDLNuE", "Number of Hits in the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Number of Hits", 300, 0, 3000);    
    
    auto sliceNumHitsBeforeCuts = createHistGroup("sliceNumHitsBeforeCuts", "Number of Hits in the Slice (Before Cuts)", "Number of Hits", 300, 0, 3000);
    auto sliceNumHitsAfterCuts = createHistGroup("sliceNumHitsAfterCuts", "Number of Hits in the Slice (After Cuts)", "Number of Hits", 300, 0, 3000);
    auto sliceNumHitsAfterCuts_splitDLNuE = createSplitHistGroup("sliceNumHits_splitDLNuE", "Number of Hits in the Slice: DL Nu+E Vertexing", "Number of Hits", 300, 0, 3000);
    auto sliceNumHitsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("sliceNumHits_splitPFPDLNuE", "Number of Hits in the Slice: DL Nu+E Vertexing", "Number of Hits", 300, 0, 3000);    
    
    auto trueRecoilElectronEnergyBeforeCuts = createHistGroup("trueRecoilElectronEnergyBeforeCuts", "Energy of True Recoil Electron (Before Cuts)", "Energy", 200, 0, 1000);
    auto trueRecoilElectronEnergyAfterCuts = createHistGroup("trueRecoilElectronEnergyAfterCuts", "Energy of True Recoil Electron (After Cuts)", "Energy", 200, 0, 1000);
    
    auto trueRecoilElectronAngleBeforeCuts = createHistGroup("trueRecoilElectronAngleBeforeCuts", "Angle of True Recoil Electron (Before Cuts)", "Angle (degrees)", 20, 0, 20);
    auto trueRecoilElectronAngleAfterCuts = createHistGroup("trueRecoilElectronAngleAfterCuts", "Angle of True Recoil Electron (After Cuts)", "Angle (degrees)", 20, 0, 20);

    auto recoVXBeforeCuts = createHistGroup("recoVXBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterCuts = createHistGroup("recoVXAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterCuts_splitDLNuE = createSplitHistGroup("recoVX_splitDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm) ", 202, -202, 202);
    auto recoVXAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVX_splitPFPDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm) ", 202, -202, 202);    
    
    auto recoVYBeforeCuts = createHistGroup("recoVYBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterCuts = createHistGroup("recoVYAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterCuts_splitDLNuE = createSplitHistGroup("recoVY_splitDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm) ", 204, -204, 204);
    auto recoVYAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVY_splitPFPDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm) ", 204, -204, 204);    
    
    auto recoVZBeforeCuts = createHistGroup("recoVZBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterCuts = createHistGroup("recoVZAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterCuts_splitDLNuE = createSplitHistGroup("recoVZ_splitDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 255, 0, 510);
    auto recoVZAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVZ_splitPFPDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 255, 0, 510);    
    
    auto recoVXSmallerBinsBeforeCuts = createHistGroup("recoVXSmallerBinsBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterCuts = createHistGroup("recoVXSmallerBinsAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("recoVXSmallerBins_splitDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm) ", 808, -202, 202);
    auto recoVXSmallerBinsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVXSmallerBins_splitPFPDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm) ", 808, -202, 202);    
    
    auto recoVYSmallerBinsBeforeCuts = createHistGroup("recoVYSmallerBinsBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterCuts = createHistGroup("recoVYSmallerBinsAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("recoVYSmallerBins_splitDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm) ", 816, -204, 204);
    auto recoVYSmallerBinsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVYSmallerBins_splitPFPDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm) ", 816, -204, 204);    
    
    auto recoVZSmallerBinsBeforeCuts = createHistGroup("recoVZSmallerBinsBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterCuts = createHistGroup("recoVZSmallerBinsAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterCuts_splitDLNuE = createSplitHistGroup("recoVZSmallerBins_splitDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 1020, 0, 510);
    auto recoVZSmallerBinsAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVZSmallerBins_splitPFPDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 1020, 0, 510);    
    
    auto recoVXLowBeforeCuts = createHistGroup("recoVXLowBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, 182);
    auto recoVXLowAfterCuts = createHistGroup("recoVXLowAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterCuts_splitDLNuE = createSplitHistGroup("recoVXLow_splitDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm) ", 40, -202, -182);
    auto recoVXLowAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVXLow_splitPFPDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm) ", 40, -202, -182);    
    
    auto recoVYLowBeforeCuts = createHistGroup("recoVYLowBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterCuts = createHistGroup("recoVYLowAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterCuts_splitDLNuE = createSplitHistGroup("recoVYLow_splitDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm) ", 40, -204, -184);
    auto recoVYLowAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVYLow_splitPFPDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm) ", 40, -204, -184);    
    
    auto recoVZLowBeforeCuts = createHistGroup("recoVZLowBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterCuts = createHistGroup("recoVZLowAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterCuts_splitDLNuE = createSplitHistGroup("recoVZLow_splitDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 40, 0, 20);
    auto recoVZLowAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVZLow_splitPFPDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 40, 0, 20);    
    
    auto recoVXHighBeforeCuts = createHistGroup("recoVXHighBeforeCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterCuts = createHistGroup("recoVXHighAfterCuts", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterCuts_splitDLNuE = createSplitHistGroup("recoVXHigh_splitDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm)", 40, 182, 202);
    auto recoVXHighAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVXHigh_splitPFPDLNuE", "X Coordinate of Reco Neutrino: DL Nu+E Vertexing", "x_{Reco} (cm)", 40, 182, 202);    
    
    auto recoVYHighBeforeCuts = createHistGroup("recoVYHighBeforeCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterCuts = createHistGroup("recoVYHighAfterCuts", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterCuts_splitDLNuE = createSplitHistGroup("recoVYHigh_splitDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm)", 40, 184, 204);
    auto recoVYHighAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVYHigh_splitPFPDLNuE", "Y Coordinate of Reco Neutrino: DL Nu+E Vertexing", "y_{Reco} (cm)", 40, 184, 204);    
    
    auto recoVZHighBeforeCuts = createHistGroup("recoVZHighBeforeCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterCuts = createHistGroup("recoVZHighAfterCuts", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterCuts_splitDLNuE = createSplitHistGroup("recoVZHigh_splitDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 40, 480, 510);
    auto recoVZHighAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("recoVZHigh_splitPFPDLNuE", "Z Coordinate of Reco Neutrino: DL Nu+E Vertexing", "z_{Reco} (cm) ", 40, 480, 510);    
    
    auto energyAsymmetryBeforeCuts = createHistGroup("energyAsymmetryBeforeCuts", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (Before Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts = createHistGroup("energyAsymmetryAfterCuts", "Energy Asymmetry of the PFP in the Slice with the Highest Energy (After Cuts)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts_splitDLNuE = createSplitHistGroup("energyAsymmetry_splitDLNuE", "Energy Asymmetry of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryAfterCuts_splitPFPDLNuE = createSplitPFPHistGroup("energyAsymmetry_splitPFPDLNuE", "Energy Asymmetry of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);    

    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0;      double zMax = 509.4;

    TH2D *xCoordAngleDifference = new TH2D("xCoordAngleDifference", "", (int)round((xMax - xMin)/5), xMin, xMax, 40, 0, 180);
    TH2D *yCoordAngleDifference = new TH2D("yCoordAngleDifference", "", (int)round((yMax - yMin)/5), yMin, yMax, 40, 0, 180);
    TH2D *zCoordAngleDifference = new TH2D("zCoordAngleDifference", "", (int)round((zMax - zMin)/5), zMin, zMax, 40, 0, 180);

    TH2D *xCoordAngleDifference_low = new TH2D("xCoordAngleDifference_low", "", 15, xMin, (xMin + 30), 40, 0, 180);
    TH2D *yCoordAngleDifference_low = new TH2D("yCoordAngleDifference_low", "", 15, yMin, (yMin + 30), 40, 0, 180);
    TH2D *zCoordAngleDifference_low = new TH2D("zCoordAngleDifference_low", "", 15, zMin, (zMin + 30), 40, 0, 180);
    
    TH2D *xCoordAngleDifference_high = new TH2D("xCoordAngleDifference_high", "", 15, (xMax - 30), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifference_high = new TH2D("yCoordAngleDifference_high", "", 15, (yMax - 30), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifference_high = new TH2D("zCoordAngleDifference_high", "", 15, (zMax - 30), zMax, 40, 0, 180);

    TH2D *xCoordEnergyAsymmetry = new TH2D("xCoordEnergyAsymmetry", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetry = new TH2D("yCoordEnergyAsymmetry", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetry = new TH2D("zCoordEnergyAsymmetry", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    
    TH2D *xCoordEnergyAsymmetry_low = new TH2D("xCoordEnergyAsymmetry_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetry_low = new TH2D("yCoordEnergyAsymmetry_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetry_low = new TH2D("zCoordEnergyAsymmetry_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    
    TH2D *xCoordEnergyAsymmetry_high = new TH2D("xCoordEnergyAsymmetry_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetry_high = new TH2D("yCoordEnergyAsymmetry_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetry_high = new TH2D("zCoordEnergyAsymmetry_high", "", 55, (zMax - 110), zMax, 20, -1, 1);

    double numEvents_DLNuECosmic = 0;
    double numEvents_DLNuEBNB = 0;
    double numEvents_DLNuENuE = 0;

    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);

        if(DLCurrent == 5 && signal == 3) numEvents_DLNuECosmic++;
        if(DLCurrent == 5 && signal == 2) numEvents_DLNuEBNB++;
        if(DLCurrent == 5 && signal == 1) numEvents_DLNuENuE++;

        // Looking at the true recoil electron in the event (if there is one)
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

        double weight = 0;
        if(signal == 1 && DLCurrent == 5) weight = weights.signalNuE;
        if(signal == 2 && DLCurrent == 5) weight = weights.BNBNuE;
        if(signal == 3 && DLCurrent == 5) weight = weights.cosmicsNuE;

        // Looking at the reco slices
        if(reco_sliceID->size() == 0) continue;

        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) == -999999) continue;
            // There is a reco slice in the event
            //std::cout << "============================== NEW SLICE ==============================" << std::endl;

            // Assigning a category to the slices
            // 0 = cosmic, 1 = signal, 2 = signal fuzzy, 3 = bnb, 4 = bnb fuzzy
            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                // This is a cosmic slice
                sliceCategoryPlottingMacro = 0;
                //std::cout << "Cosmic Slice: sliceCategoryPlottingMacro = 0" << std::endl;
            } else if(reco_sliceOrigin->at(slice) == 1){
                // This is a nu+e elastic scatter slice
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    sliceCategoryPlottingMacro = 1;
                    //std::cout << "Nu+E Slice: sliceCategoryPlottingMacro = 1" << std::endl;
                } else{
                    sliceCategoryPlottingMacro = 2;
                    //std::cout << "Nu+E Fuzzy Slice: sliceCategoryPlottingMacro = 2" << std::endl;
                }
            } else if(reco_sliceOrigin->at(slice) == 3){
                // This is a BNB slice
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    sliceCategoryPlottingMacro = 3;
                    //std::cout << "BNB Slice: sliceCategoryPlottingMacro = 3" << std::endl;
                } else{
                    sliceCategoryPlottingMacro = 4;
                    //std::cout << "BNB Fuzzy Slice: sliceCategoryPlottingMacro = 4" << std::endl;
                }
            }

            //std::cout << "Slice Origin = " << reco_sliceOrigin->at(slice) << ", CCNC = " << reco_sliceTrueCCNC->at(slice) << ", True Neutrino Type = " << reco_sliceTrueNeutrinoType->at(slice) << ", Vertex = (" << reco_sliceTrueVX->at(slice) << ", " << reco_sliceTrueVY->at(slice) << ", " << reco_sliceTrueVZ->at(slice) << ")" << std::endl;
            for(size_t trueParticle = 0; trueParticle < truth_particleSliceID->size(); trueParticle++){
                if(truth_particleSliceID->at(trueParticle) == reco_sliceID->at(slice)){
                    if(truth_particleStatusCode->at(trueParticle) == 1){
                        //std::cout << "True particle in slice: PDG = " << truth_particlePDG->at(trueParticle) << std::endl;
                    }
                }
            }

            // Assigning a interaction category to the slices
            // Event types: Cosmic = 0, nu+e scatter = 1, NC Npi0 = 2, other NC = 3, CC numu = 4, CC nue = 5, Dirt = 6, Dirt nu+e = 7
            int sliceInteractionType = -999999;
            if(reco_sliceOrigin->at(slice) != 0){
                // This is a slice that isn't truth-matched to a cosmic
                if(reco_sliceOrigin->at(slice) == 1){
                    // This is a slice that is truth-matched to a nu+e elastic scatter
                    if(reco_sliceTrueVX->at(slice) > -201.3 && reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVY->at(slice) > -203.8 && reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5){
                        // Interaction happened inside the TPC
                        sliceInteractionType = 1;
                    } else{
                        sliceInteractionType = 7;
                    }
                } else if(reco_sliceOrigin->at(slice) == 3){
                    // This is a slice that is truth-matched to a beam neutrino that isn't a nu+e elastic scatter
                    if(reco_sliceTrueVX->at(slice) > -201.3 && reco_sliceTrueVX->at(slice) < 201.3 && reco_sliceTrueVY->at(slice) > -203.8 && reco_sliceTrueVY->at(slice) < 203.8 && reco_sliceTrueVZ->at(slice) > 0 && reco_sliceTrueVZ->at(slice) < 509.5){
                        // Interaction happened inside the TPC
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
                        // Interaction happened outside the TPC - Dirt event
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

            //std::cout << "Slice interaction type assigned = " << sliceInteractionType << std::endl;
            if(((sliceInteractionType != 1) && sliceInteractionType != 7) && ((sliceCategoryPlottingMacro == 1) || (sliceCategoryPlottingMacro == 2))){
                //std::cout << "This is a nu+e slice but slice int type is not a nu+e scatter!!" << std::endl;
            }
            
            double summedEnergy_beforeCuts = 0;
            double numPFPsSlice_beforeCuts = 0;
            double numPrimaryPFPsSlice_beforeCuts = 0;
            double numHitsInPFPs_beforeCuts = 0;
            
            highestEnergyPFP_struct highestEnergyPFP_beforeCuts;

            //std::cout << "------ PFPs before cuts ------" << std::endl;
            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                    // PFP is in the slice
                    numPFPsSlice_beforeCuts++;
                    if(reco_particleIsPrimary->at(pfp) == 1) numPrimaryPFPsSlice_beforeCuts++;
                    numHitsInPFPs_beforeCuts += reco_particleNumHits->at(pfp);

                    //std::cout << "PFP " << pfp << ": Energy = " << reco_particleBestPlaneEnergy->at(pfp) << ", Clear Cosmic = " << reco_particleClearCosmic->at(pfp) << ", True PDG = " << reco_particleTruePDG->at(pfp) << ", Vertex = (" << reco_particleVX->at(pfp) << ", " << reco_particleVY->at(pfp) << ", " << reco_particleVZ->at(pfp) << ")" << std::endl;
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
            }
            
            //std::cout << "------------------------------" << std::endl;
            //std::cout << "Before Cuts, PFP with Highest Energy: Energy = " << highestEnergyPFP_beforeCuts.energy << ", Clear Cosmic = " << highestEnergyPFP_beforeCuts.clearCosmic << ", True PDG = " << highestEnergyPFP_beforeCuts.truePDG << ", Vertex = (" << highestEnergyPFP_beforeCuts.vx << ", " << highestEnergyPFP_beforeCuts.vy << ", " << highestEnergyPFP_beforeCuts.vz << ")" << std::endl; 

            // Looped through all PFPs in the slice and now have the highest energy PFP out
            double angleDifference_beforeCuts = -999999;
            if((highestEnergyPFP_beforeCuts.dx != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((highestEnergyPFP_beforeCuts.dx * recoilElectron.dx) + (highestEnergyPFP_beforeCuts.dy * recoilElectron.dy) + (highestEnergyPFP_beforeCuts.dz * recoilElectron.dz));
                double aMagnitude = std::sqrt((highestEnergyPFP_beforeCuts.dx * highestEnergyPFP_beforeCuts.dx) + (highestEnergyPFP_beforeCuts.dy * highestEnergyPFP_beforeCuts.dy) + (highestEnergyPFP_beforeCuts.dz * highestEnergyPFP_beforeCuts.dz));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifference_beforeCuts = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            }

            double recoVX = -999999;
            double recoVY = -999999;
            double recoVZ = -999999;
            int numRecoNeutrinos = 0;

            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    // Reco neutrino is in the slice
                    numRecoNeutrinos++;
                    recoVX = reco_neutrinoVX->at(recoNeut);
                    recoVY = reco_neutrinoVY->at(recoNeut);
                    recoVZ = reco_neutrinoVZ->at(recoNeut);
                }
            }


            // Assigning event type based on true PDG of highest energy PFP in slice
            int slicePFPType_beforeCuts = -999999;
            if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 11 && highestEnergyPFP_beforeCuts.trueInt == 1098 && highestEnergyPFP_beforeCuts.trueOrigin == 1 && signal == 1){
                // This is an electron from a nu+e elastic scatter
                slicePFPType_beforeCuts = 0;
            } else if(highestEnergyPFP_beforeCuts.trueInt == 1098 && highestEnergyPFP_beforeCuts.trueOrigin == 1 && signal == 1){
                // This is something other than an electron from a nu+e elastic scatter
                slicePFPType_beforeCuts = 1;
            } else if(highestEnergyPFP_beforeCuts.trueInt == 1098 && signal != 1){
                // This is a nu+e elastic scatter not from the signal file
                slicePFPType_beforeCuts = 15;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 11 && highestEnergyPFP_beforeCuts.trueOrigin == 1 && highestEnergyPFP_beforeCuts.trueInt != 1098){
                // This is an electron from a beam neutrino
                slicePFPType_beforeCuts = 2;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 2212 && highestEnergyPFP_beforeCuts.trueOrigin == 1){
                // This is a proton from a beam neutrino
                slicePFPType_beforeCuts = 3;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 13 && highestEnergyPFP_beforeCuts.trueOrigin == 1){
                // This is a muon from a beam neutrino
                slicePFPType_beforeCuts = 4;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 111 && highestEnergyPFP_beforeCuts.trueOrigin == 1){
                // This is a pi0 fron a beam neutrino
                slicePFPType_beforeCuts = 5;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 211 && highestEnergyPFP_beforeCuts.trueOrigin == 1){
                // This is a charged pi from a beam neutrino
                slicePFPType_beforeCuts = 6;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 22 && highestEnergyPFP_beforeCuts.trueOrigin == 1){
                // This is a proton from a beam neutrino
                slicePFPType_beforeCuts = 7;
            } else if(highestEnergyPFP_beforeCuts.trueOrigin == 1){
                // This is something else from a beam neutrino
                slicePFPType_beforeCuts = 8;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 13 && highestEnergyPFP_beforeCuts.trueOrigin == 2){
                // This is a muon from a cosmic
                slicePFPType_beforeCuts = 9;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 22 && highestEnergyPFP_beforeCuts.trueOrigin == 2){
                // This is a photon from a cosmic
                slicePFPType_beforeCuts = 10;
            } else if(std::abs(highestEnergyPFP_beforeCuts.truePDG) == 11 && highestEnergyPFP_beforeCuts.trueOrigin == 2){
                // This is an electron from a cosmic
                slicePFPType_beforeCuts = 11;
            } else if(highestEnergyPFP_beforeCuts.trueOrigin == 2){
                // This is something else from a cosmic
                slicePFPType_beforeCuts = 12;
            }

            // Counting number of events before cuts
            if(DLCurrent == 5){
                if(sliceCategoryPlottingMacro == 0){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 1 && signal == 1){
                    eventsBeforeCuts_DLNuE.signal += weight;
                } else if(sliceCategoryPlottingMacro == 2 && signal == 1){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 3){
                    eventsBeforeCuts_DLNuE.background += weight;
                } else if(sliceCategoryPlottingMacro == 4){
                    eventsBeforeCuts_DLNuE.background += weight;
                }

                if(sliceInteractionType == 0){
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
                }
            }

            // Add to plots before cuts
            //std::cout << "Slice Type = " << sliceCategoryPlottingMacro << ", Interaction Type = " << sliceInteractionType << ", True PDG of Highest Energy PFP in Slice = " << slicePFPType_beforeCuts << std::endl;
            //std::cout << "Filling histograms before cuts applied ======================================" << std::endl;
            fillHistogram(&sliceCompletenessBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceCompleteness->at(slice), &weights);
            fillHistogram(&sliceCRUMBSBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
            fillHistogram(&slicePurityBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_slicePurity->at(slice), &weights);
            fillHistogram(&sliceNumRecoNeutBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
            fillHistogram(&sliceNumPFPsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_beforeCuts, &weights);
            fillHistogram(&sliceNumPrimaryPFPsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_beforeCuts, &weights);
            fillHistogram(&sliceFracHitsInPFPsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (numHitsInPFPs_beforeCuts/reco_sliceNumHits->at(slice)), &weights);
            fillHistogram(&sliceFracHitsInHighestEnergyPFPsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_beforeCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
            fillHistogram(&ERecoSumThetaRecoBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillHistogram(&ERecoHighestThetaRecoBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
            fillHistogram(&dEdxBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
            fillHistogram(&razzledPDG11BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
            fillHistogram(&razzledPDG13BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
            fillHistogram(&razzledPDG22BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
            fillHistogram(&razzledPDG211BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
            fillHistogram(&razzledPDG2212BeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
            fillHistogram(&pfpCompletenessBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.completeness, &weights);
            fillHistogram(&pfpPurityBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.purity, &weights);
            fillHistogram(&pfpNumHitsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.numHits, &weights);
            fillHistogram(&sliceNumHitsBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceNumHits->at(slice), &weights);
            if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.energy != -999999) fillHistogram(&trueRecoilElectronEnergyBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoilElectron.energy, &weights);
            if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999) fillHistogram(&trueRecoilElectronAngleBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (recoilElectron.angle*TMath::RadToDeg()), &weights);
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
            //std::cout << "finished =========================================================================================" << std::endl;
            if((sliceCategoryPlottingMacro == 1 || sliceCategoryPlottingMacro == 2) && signal == 1) fillHistogram(&energyAsymmetryBeforeCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy)/recoilElectron.energy), &weights);

            // Add to plots if no clear cosmic score is cut on
            if(clearCosmicCut == 0){
                // Uses all the same variables as the plots before cuts - no PFPs are removed
                fillHistogram(&sliceCompletenessAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceCompleteness->at(slice), &weights);
                fillSplitIntHistogram(&sliceCompletenessAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceCompleteness->at(slice), &weights);
                fillSplitPFPHistogram(&sliceCompletenessAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, reco_sliceCompleteness->at(slice), &weights);
                
                fillHistogram(&sliceCRUMBSAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
                fillSplitIntHistogram(&sliceCRUMBSAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
                fillSplitPFPHistogram(&sliceCRUMBSAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, reco_sliceScore->at(slice), &weights);
                
                fillHistogram(&slicePurityAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_slicePurity->at(slice), &weights);
                fillSplitIntHistogram(&slicePurityAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_slicePurity->at(slice), &weights);
                fillSplitPFPHistogram(&slicePurityAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, reco_slicePurity->at(slice), &weights);
                
                fillHistogram(&sliceNumRecoNeutAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
                fillSplitIntHistogram(&sliceNumRecoNeutAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);
                fillSplitPFPHistogram(&sliceNumRecoNeutAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, numRecoNeutrinos, &weights);
                
                fillHistogram(&sliceNumPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_beforeCuts, &weights);
                fillSplitIntHistogram(&sliceNumPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPFPsSlice_beforeCuts, &weights);
                fillSplitPFPHistogram(&sliceNumPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, numPFPsSlice_beforeCuts, &weights);
                
                fillHistogram(&sliceNumPrimaryPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_beforeCuts, &weights);
                fillSplitIntHistogram(&sliceNumPrimaryPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_beforeCuts, &weights);
                fillSplitPFPHistogram(&sliceNumPrimaryPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, numPrimaryPFPsSlice_beforeCuts, &weights);
                
                fillHistogram(&sliceFracHitsInPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (numHitsInPFPs_beforeCuts/reco_sliceNumHits->at(slice)), &weights);
                fillSplitIntHistogram(&sliceFracHitsInPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (numHitsInPFPs_beforeCuts/reco_sliceNumHits->at(slice)), &weights);
                fillSplitPFPHistogram(&sliceFracHitsInPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, (numHitsInPFPs_beforeCuts/reco_sliceNumHits->at(slice)), &weights);
                
                fillHistogram(&sliceFracHitsInHighestEnergyPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_beforeCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
                fillSplitIntHistogram(&sliceFracHitsInHighestEnergyPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_beforeCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
                fillSplitPFPHistogram(&sliceFracHitsInHighestEnergyPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, (highestEnergyPFP_beforeCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
                
                fillHistogram(&ERecoSumThetaRecoAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
                fillSplitIntHistogram(&ERecoSumThetaRecoAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
                fillSplitPFPHistogram(&ERecoSumThetaRecoAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, (summedEnergy_beforeCuts * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
                
                fillHistogram(&ERecoHighestThetaRecoAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
                fillSplitIntHistogram(&ERecoHighestThetaRecoAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
                fillSplitPFPHistogram(&ERecoHighestThetaRecoAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, (highestEnergyPFP_beforeCuts.energy * highestEnergyPFP_beforeCuts.theta * highestEnergyPFP_beforeCuts.theta), &weights);
                
                fillHistogram(&dEdxAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
                fillSplitIntHistogram(&dEdxAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
                fillSplitPFPHistogram(&dEdxAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.bestPlanedEdx, &weights);
                
                fillHistogram(&razzledPDG11AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
                fillSplitIntHistogram(&razzledPDG11AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
                fillSplitPFPHistogram(&razzledPDG11AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.razzledPDG11, &weights);
                
                fillHistogram(&razzledPDG13AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
                fillSplitIntHistogram(&razzledPDG13AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
                fillSplitPFPHistogram(&razzledPDG13AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.razzledPDG13, &weights);
                
                fillHistogram(&razzledPDG22AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
                fillSplitIntHistogram(&razzledPDG22AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
                fillSplitPFPHistogram(&razzledPDG22AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.razzledPDG22, &weights);
                
                fillHistogram(&razzledPDG211AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
                fillSplitIntHistogram(&razzledPDG211AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
                fillSplitPFPHistogram(&razzledPDG211AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.razzledPDG211, &weights);
                
                fillHistogram(&razzledPDG2212AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
                fillSplitIntHistogram(&razzledPDG2212AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
                fillSplitPFPHistogram(&razzledPDG2212AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.razzledPDG2212, &weights);
                
                fillHistogram(&pfpCompletenessAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.completeness, &weights);
                fillSplitIntHistogram(&pfpCompletenessAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.completeness, &weights);
                fillSplitPFPHistogram(&pfpCompletenessAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.completeness, &weights);
                
                fillHistogram(&pfpPurityAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.purity, &weights);
                fillSplitIntHistogram(&pfpPurityAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.purity, &weights);
                fillSplitPFPHistogram(&pfpPurityAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.purity, &weights);
                
                fillHistogram(&pfpNumHitsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_beforeCuts.numHits, &weights);
                fillSplitIntHistogram(&pfpNumHitsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_beforeCuts.numHits, &weights);
                fillSplitPFPHistogram(&pfpNumHitsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, highestEnergyPFP_beforeCuts.numHits, &weights);
                
                fillHistogram(&sliceNumHitsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceNumHits->at(slice), &weights);
                fillSplitIntHistogram(&sliceNumHitsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceNumHits->at(slice), &weights);
                fillSplitPFPHistogram(&sliceNumHitsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, reco_sliceNumHits->at(slice), &weights);
            
                if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.energy != -999999) fillHistogram(&trueRecoilElectronEnergyAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoilElectron.energy, &weights);
                if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999) fillHistogram(&trueRecoilElectronAngleAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (recoilElectron.angle*TMath::RadToDeg()), &weights);
                
                fillHistogram(&recoVXAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVX, &weights);
                
                fillHistogram(&recoVYAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVY, &weights);
                
                fillHistogram(&recoVZAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVZ, &weights);
                
                fillHistogram(&recoVXSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXSmallerBinsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVX, &weights);
                
                fillHistogram(&recoVYSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYSmallerBinsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVY, &weights);
                
                fillHistogram(&recoVZSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZSmallerBinsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVZ, &weights);
                
                fillHistogram(&recoVXLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXLowAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVX, &weights);
                
                fillHistogram(&recoVYLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYLowAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVY, &weights);
                
                fillHistogram(&recoVZLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZLowAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVZ, &weights);
                
                fillHistogram(&recoVXHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXHighAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVX, &weights);
                
                fillHistogram(&recoVYHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYHighAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVY, &weights);
                
                fillHistogram(&recoVZHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZHighAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, recoVZ, &weights);
                
                if((sliceCategoryPlottingMacro == 1 || sliceCategoryPlottingMacro == 2) && signal == 1){
                    fillHistogram(&energyAsymmetryAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy), &weights);
                    fillSplitIntHistogram(&energyAsymmetryAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy), &weights);
                    fillSplitPFPHistogram(&energyAsymmetryAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy), &weights);
                
                    xCoordEnergyAsymmetry->Fill(recoVX, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                    xCoordAngleDifference->Fill(recoVX, angleDifference_beforeCuts);
                
                    yCoordEnergyAsymmetry->Fill(recoVY, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                    yCoordAngleDifference->Fill(recoVY, angleDifference_beforeCuts);
                
                    zCoordEnergyAsymmetry->Fill(recoVZ, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                    zCoordAngleDifference->Fill(recoVZ, angleDifference_beforeCuts);

                    if(recoVX >= xMin && recoVX <= xMin+30){
                        xCoordEnergyAsymmetry_low->Fill(recoVX, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                        xCoordAngleDifference_low->Fill(recoVX, angleDifference_beforeCuts);
                    } else if(recoVX <= xMax && recoVX >= xMax-30){
                        xCoordEnergyAsymmetry_high->Fill(recoVX, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                        xCoordAngleDifference_high->Fill(recoVX, angleDifference_beforeCuts);
                    }

                    if(recoVY >= yMin && recoVY <= yMin+30){
                        yCoordEnergyAsymmetry_low->Fill(recoVY, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                        yCoordAngleDifference_low->Fill(recoVY, angleDifference_beforeCuts);
                    } else if(recoVY <= yMax && recoVY >= yMax-30){
                        yCoordEnergyAsymmetry_high->Fill(recoVY, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                        yCoordAngleDifference_high->Fill(recoVY, angleDifference_beforeCuts);
                    }

                    if(recoVZ >= zMin && recoVZ <= zMin+30){
                        zCoordEnergyAsymmetry_low->Fill(recoVZ, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                        zCoordAngleDifference_low->Fill(recoVZ, angleDifference_beforeCuts);
                    } else if(recoVZ <= zMax && recoVZ >= zMax-110){
                        zCoordEnergyAsymmetry_high->Fill(recoVZ, ((recoilElectron.energy - highestEnergyPFP_beforeCuts.energy) /recoilElectron.energy));
                        zCoordAngleDifference_high->Fill(recoVZ, angleDifference_beforeCuts);
                    }

                }
            }

            // If applying the clear cosmic score then do this
            // Loop through PFPs again to find the highest energy PFP, also reassign the slicePFPType_afterCuts
            if(clearCosmicCut == 1){
                double summedEnergy_afterCuts = 0;
                double numPFPsSlice_afterCuts = 0;
                double numPrimaryPFPsSlice_afterCuts = 0;
                double numHitsInPFPs_afterCuts = 0;

                highestEnergyPFP_struct highestEnergyPFP_afterCuts;

                //std::cout << "------ PFPs after cuts ------" << std::endl;
                for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                    if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                        // PFP is in the slice
                        if(reco_particleClearCosmic->at(pfp) == 0){
                            // PFP is not a clear cosmic
                            //std::cout << "PFP " << pfp << ": Energy = " << reco_particleBestPlaneEnergy->at(pfp) << ", Clear Cosmic = " << reco_particleClearCosmic->at(pfp) << ", True PDG = " << reco_particleTruePDG->at(pfp) << ", Vertex = (" << reco_particleVX->at(pfp) << ", " << reco_particleVY->at(pfp) << ", " << reco_particleVZ->at(pfp) << ")" << std::endl;
                            numPFPsSlice_afterCuts++;
                            if(reco_particleIsPrimary->at(pfp) == 1) numPrimaryPFPsSlice_afterCuts++; // PFP is a primary PFP
                            numHitsInPFPs_afterCuts += reco_particleNumHits->at(pfp);

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
                }
            
                //std::cout << "------------------------------" << std::endl;
                //std::cout << "After Clear Cosmic Cuts, PFP with Highest Energy: Energy = " << highestEnergyPFP_afterCuts.energy << ", Clear Cosmic = " << highestEnergyPFP_afterCuts.clearCosmic << ", True PDG = " << highestEnergyPFP_afterCuts.truePDG << ", Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << ")" << std::endl; 

                // Looped through all the PFPs that aren't clear cosmics and how have the highest energy PFP out
                double angleDifference_afterCuts = -999999;
                if((highestEnergyPFP_afterCuts.dx != -999999) && (recoilElectron.dx != -999999)){
                    double aDOTb = ((highestEnergyPFP_afterCuts.dx * recoilElectron.dx) + (highestEnergyPFP_afterCuts.dy * recoilElectron.dy) + (highestEnergyPFP_afterCuts.dz * recoilElectron.dz));
                    double aMagnitude = std::sqrt((highestEnergyPFP_afterCuts.dx * highestEnergyPFP_afterCuts.dx) + (highestEnergyPFP_afterCuts.dy * highestEnergyPFP_afterCuts.dy) + (highestEnergyPFP_afterCuts.dz * highestEnergyPFP_afterCuts.dz));
                    double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                    double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                    angleDifference_afterCuts = (TMath::ACos(cosAngle) * TMath::RadToDeg());
                }

                // Assigning event type based on true PDG of highest energy PFP in slice
                int slicePFPType_afterCuts = -999999;
                if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 11 && highestEnergyPFP_afterCuts.trueInt == 1098 && highestEnergyPFP_afterCuts.trueOrigin == 1 && signal == 1){
                    // This is an electron from a nu+e elastic scatter
                    slicePFPType_afterCuts = 0;
                } else if(highestEnergyPFP_afterCuts.trueInt == 1098 && highestEnergyPFP_afterCuts.trueOrigin == 1 && signal == 1){
                    // This is something other than an electron from a nu+e elastic scatter
                    slicePFPType_afterCuts = 1;
                } else if(highestEnergyPFP_afterCuts.trueInt == 1098 && signal != 1){
                    // This is a nu+e elastic scatter that isn't from the signal file
                    slicePFPType_afterCuts = 15;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 11 && highestEnergyPFP_afterCuts.trueOrigin == 1 && highestEnergyPFP_afterCuts.trueInt != 1098){
                    // This is an electron from a beam neutrino
                    slicePFPType_afterCuts = 2;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 2212 && highestEnergyPFP_afterCuts.trueOrigin == 1){
                    // This is a proton from a beam neutrino
                    slicePFPType_afterCuts = 3;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 13 && highestEnergyPFP_afterCuts.trueOrigin == 1){
                    // This is a muon from a beam neutrino
                    slicePFPType_afterCuts = 4;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 111 && highestEnergyPFP_afterCuts.trueOrigin == 1){
                    // This is a pi0 fron a beam neutrino
                    slicePFPType_afterCuts = 5;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 211 && highestEnergyPFP_afterCuts.trueOrigin == 1){
                    // This is a charged pi from a beam neutrino
                    slicePFPType_afterCuts = 6;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 22 && highestEnergyPFP_afterCuts.trueOrigin == 1){
                    // This is a proton from a beam neutrino
                    slicePFPType_afterCuts = 7;
                } else if(highestEnergyPFP_afterCuts.trueOrigin == 1){
                    // This is something else from a beam neutrino
                    slicePFPType_afterCuts = 8;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 13 && highestEnergyPFP_afterCuts.trueOrigin == 2){
                    // This is a muon from a cosmic
                    slicePFPType_afterCuts = 9;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 22 && highestEnergyPFP_afterCuts.trueOrigin == 2){
                    // This is a photon from a cosmic
                    slicePFPType_afterCuts = 10;
                } else if(std::abs(highestEnergyPFP_afterCuts.truePDG) == 11 && highestEnergyPFP_afterCuts.trueOrigin == 2){
                    // This is an electron from a cosmic
                    slicePFPType_afterCuts = 11;
                } else if(highestEnergyPFP_afterCuts.trueOrigin == 2){
                    // This is something else from a cosmic
                    slicePFPType_afterCuts = 12;
                }

                // Clear cosmic cut has been applied, add to counters
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.clearCosmicsBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other += weight;
                }

                // Apply cuts here
                if(numPFPs0Cut == 1 && numPFPsSlice_afterCuts == 0){
                    // This is a slice with 0 PFPs in it
                    continue;
                }
                
                // Number of PFPs 0 cut has been applied, add to counters
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numPFPs0Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numPFPs0IntSplit.other += weight;
                }

                if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0){
                    // This is a slice with no reco neutrino
                    continue;
                }

                // Number of reco neutrinos cut has been applied, add to counters
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numRecoNeut0Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other += weight;
                }

                if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low || reco_sliceScore->at(slice) > crumbsScoreCut_high)){
                    // This is a slice with a CRUMBS score outside cut values
                    if(sliceInteractionType == 1 && signal == 1 && crumbsPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails CRUMBS Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "True Vertex = (" << reco_sliceTrueVX->at(slice) << ", " << reco_sliceTrueVY->at(slice) << ", " << reco_sliceTrueVZ->at(slice) << ")" << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }
                
                // CRUMBS score cut has been applied, add to counters
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.crumbsBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.crumbsSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.crumbsBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.crumbsBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.crumbsBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.crumbsIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.crumbsIntSplit.other += weight;
                }

                if(FVCut == 1){
                    if(!(recoVX < FVCut_xHigh && recoVX > FVCut_xLow  && std::abs(recoVX) > FVCut_xCentre && recoVY < FVCut_yHigh && recoVY > FVCut_yLow && recoVZ > FVCut_zLow && recoVZ < FVCut_zHigh)){
                        // Doesn't pass the FV cut values
                        if(sliceInteractionType == 1 && signal == 1 && FVCutPrint == 1){
                            std::cout << "=============" << std::endl;
                            std::cout << "Signal Fails FV Cut" << std::endl;
                            std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                            std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                            std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                            std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                            std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                            std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                            std::cout << "=============" << std::endl;
                        }
                        continue;
                    }
                }

                // FV cut applied, fill counters
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.FVBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.FVSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.FVBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.FVBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.FVBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.FVIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.FVIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.FVIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.FVIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.FVIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.FVIntSplit.other += weight;
                }

                if(primaryPFPCut == 1 && numPrimaryPFPsSlice_afterCuts != primaryPFPCutValue){
                    // Slice has more than 1 primary PFP in it
                    if(sliceInteractionType == 1 && signal == 1 && primaryPFPPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Primary PFP Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                // Primary PFP cut has been applied, fill counters
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.primaryPFPBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.primaryPFPIntSplit.other += weight;
                }
                
                if(ETheta2Cut == 1 && ((highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) > ETheta2High_highestEnergyPFP || (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) < ETheta2Low_highestEnergyPFP)){
                    // Highest energy PFP in slice doesn't pass the ETheta2 cut
                    if(sliceInteractionType == 1 && signal == 1 && ETheta2Print == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails ETheta2 Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }

                    continue;
                }
                    
                if(sliceInteractionType == 1 && signal == 1 && ETheta2Print == 1){
                    std::cout << "=============" << std::endl;
                    std::cout << "Signal Passes ETheta2 Cut" << std::endl;
                    std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                    std::cout << "=============" << std::endl;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.ETheta2Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.ETheta2Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.ETheta2Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.ETheta2IntSplit.other += weight;
                }


                if(razzledPDG2212Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG2212 > razzled2212High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG2212 < razzled2212Low_highestEnergyPFP))){
                    // Highest energy PFP in slice doesn't pass the razzled 2212 cut
                    if(sliceInteractionType == 1 && signal == 1 && razzledProtonPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Razzled Proton Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }

                    continue;
                }
                
                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled2212Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled2212Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled2212Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.razzled2212IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled2212IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled2212IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled2212IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled2212IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled2212IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled2212IntSplit.other += weight;
                }

                if(razzledPDG13Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG13 > razzled13High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG13 < razzled13Low_highestEnergyPFP))){
                    // Highest energy PFP in slice doesn't pass the razzled 13 cut
                    if(sliceInteractionType == 1 && signal == 1 && razzledMuonPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Razzled Muon Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.razzled13Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled13Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled13Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled13Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled13Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.razzled13IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled13IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled13IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled13IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled13IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled13IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled13IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled13IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled13IntSplit.other += weight;
                }

                if(razzledPDG211Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG211 > razzled211High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG211 < razzled211Low_highestEnergyPFP))){
                    // Highest energy PFP in slice doesn't pass the razzled 211 cut
                    if(sliceInteractionType == 1 && signal == 1 && razzledChargedPiPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Razzled Charged Pi Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.razzled211Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled211Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled211Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled211Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled211Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled211IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled211IntSplit.other += weight;
                }

                if(razzledPDG22Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG22 > razzled22High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG22 < razzled22Low_highestEnergyPFP))){
                    // Highest energy PFP in slice doesn't pass the razzled 22 cut
                    if(sliceInteractionType == 1 && signal == 1 && razzledPhotonPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Razzled Photon Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.razzled22Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled22Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled22Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled22Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled22Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.razzled22IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled22IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled22IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled22IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled22IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled22IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled22IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled22IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled22IntSplit.other += weight;
                }

                if(razzledPDG11Cut == 1 && ((highestEnergyPFP_afterCuts.razzledPDG11 > razzled11High_highestEnergyPFP) || (highestEnergyPFP_afterCuts.razzledPDG11 < razzled11Low_highestEnergyPFP))){
                    // Highest energy PFP in slice doesn't pass the razzled 11 cut
                    if(sliceInteractionType == 1 && signal == 1 && razzledElectronPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Razzled Electron Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.razzled11Back += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled11Sig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.razzled11Back += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.razzled11Back += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.razzled11Back += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.razzled11IntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.razzled11IntSplit.other += weight;
                }

                if(dEdxCut == 1 && (highestEnergyPFP_afterCuts.bestPlanedEdx > dEdxHigh_highestEnergyPFP || highestEnergyPFP_afterCuts.bestPlanedEdx < dEdxLow_highestEnergyPFP)){
                    // Highest energy PFP in slice doesn't pass the dE/dx cut
                    if(sliceInteractionType == 1 && signal == 1 && dEdxPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails dE/dx Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.dEdxBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.dEdxSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.dEdxBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.dEdxBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.dEdxBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.dEdxIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.dEdxIntSplit.other += weight;
                }

                if(fracHitsContainedCut == 1 && (((numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) > numContainedHits_high) || ((numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) < numContainedHits_low))){
                    // Highest energy PFP in slice doesn't pass the number of hits contained in PFPs cut
                    if(sliceInteractionType == 1 && signal == 1 && fracHitsPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Frac Hits Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.fracHitsContainedBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.fracHitsContainedSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.fracHitsContainedBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.fracHitsContainedBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.fracHitsContainedBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.other += weight;
                }

                if(numHitsCut == 1 && ((highestEnergyPFP_afterCuts.numHits > numHits_high) || (highestEnergyPFP_afterCuts.numHits < numHits_low))){
                    // Highest energy PFP in slice doesn't pass the number of hits contained in PFPs cut
                    if(sliceInteractionType == 1 && signal == 1 && numHitsPrint == 1){
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Fails Number of Hits Cut" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        std::cout << "=============" << std::endl;
                    }
                    continue;
                }

                if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) eventsAfterCuts_DLNuE.numHitsBack += weight;
                    else if(sliceCategoryPlottingMacro == 1 && signal == 1) eventsAfterCuts_DLNuE.numHitsSig += weight;
                    else if(sliceCategoryPlottingMacro == 2 && signal == 1) eventsAfterCuts_DLNuE.numHitsBack += weight;
                    else if(sliceCategoryPlottingMacro == 3) eventsAfterCuts_DLNuE.numHitsBack += weight;
                    else if(sliceCategoryPlottingMacro == 4) eventsAfterCuts_DLNuE.numHitsBack += weight;
    
                    if(sliceInteractionType == 0) eventsAfterCuts_DLNuE.numHitsIntSplit.cosmic += weight;
                    else if(sliceInteractionType == 1 && signal == 1) eventsAfterCuts_DLNuE.numHitsIntSplit.nuE += weight;
                    else if(sliceInteractionType == 2) eventsAfterCuts_DLNuE.numHitsIntSplit.NCNPi0 += weight;
                    else if(sliceInteractionType == 3) eventsAfterCuts_DLNuE.numHitsIntSplit.otherNC += weight;
                    else if(sliceInteractionType == 4) eventsAfterCuts_DLNuE.numHitsIntSplit.CCnumu += weight;
                    else if(sliceInteractionType == 5) eventsAfterCuts_DLNuE.numHitsIntSplit.CCnue += weight;
                    else if(sliceInteractionType == 6) eventsAfterCuts_DLNuE.numHitsIntSplit.dirt += weight;
                    else if(sliceInteractionType == 7 && signal == 1) eventsAfterCuts_DLNuE.numHitsIntSplit.nuEDirt += weight;
                    else if(sliceInteractionType == 8) eventsAfterCuts_DLNuE.numHitsIntSplit.other += weight;
                }

                // Fill histograms here
                if((signal == 1 && sliceInteractionType == 1) || (signal == 1 && sliceInteractionType == 7)){
                    if((slicePFPType_afterCuts == 0 || slicePFPType_afterCuts == 1) && (slicePFPType_afterCuts != 1 && slicePFPType_afterCuts != 7)){
                        //std::cout << "DIFFERENT: sliceInteractionType = " << sliceInteractionType << ", slicePFPType_afterCuts = " << slicePFPType_afterCuts << ", slicePFPType_beforeCuts = " << slicePFPType_beforeCuts << std::endl;
                        //std::cout << "highestEnergyPFP_afterCuts.truePDG = " << highestEnergyPFP_afterCuts.truePDG << ", highestEnergyPFP_afterCuts.trueInt = " << highestEnergyPFP_afterCuts.trueInt << ", highestEnergyPFP_afterCuts.trueOrigin = " << highestEnergyPFP_afterCuts.trueOrigin << std::endl;
                        //std::cout << "highestEnergyPFP_beforeCuts.truePDG = " << highestEnergyPFP_beforeCuts.truePDG << ", highestEnergyPFP_beforeCuts.trueInt = " << highestEnergyPFP_beforeCuts.trueInt << ", highestEnergyPFP_beforeCuts.trueOrigin = " << highestEnergyPFP_beforeCuts.trueOrigin << std::endl;
                        //std::cout << "highestEnergyPFP_beforeCuts.clearCosmic = " << highestEnergyPFP_beforeCuts.clearCosmic << ", highestEnergyPFP_afterCuts.clearCosmic = " << highestEnergyPFP_afterCuts.clearCosmic << std::endl;
                    }
                }
               
                if((slicePFPType_afterCuts == 0 || slicePFPType_afterCuts == 1) && (sliceInteractionType != 1 && sliceInteractionType != 7)){
                    std::cout << "HERE!" << std::endl;
                    std::cout << "slicePFPType_afterCuts = " << slicePFPType_afterCuts << ", sliceInteractionType = " << sliceInteractionType << std::endl;
                    std::cout << "reco_sliceOrigin->at(slice) = " << reco_sliceOrigin->at(slice) << std::endl;
                    std::cout << "signal = " << signal << ", sliceCategoryPlottingMacro = " << sliceCategoryPlottingMacro << std::endl;
                    std::cout << "Highest energy PFP before clear cosmic: Energy = " << highestEnergyPFP_beforeCuts.energy << ", true origin = " << highestEnergyPFP_beforeCuts.trueOrigin << ", Vertex = (" << highestEnergyPFP_beforeCuts.vx << ", " << highestEnergyPFP_beforeCuts.vy << ", " << highestEnergyPFP_beforeCuts.vz << "), true pdg = " << highestEnergyPFP_beforeCuts.truePDG << ", true origin = " << highestEnergyPFP_beforeCuts.trueOrigin << ", true int = " << highestEnergyPFP_beforeCuts.trueInt << std::endl;
                    std::cout << "Highest energy PFP after clear cosmic: Enegry = " << highestEnergyPFP_afterCuts.energy << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), true pdg = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << std::endl; 
                }

                if(slicePFPType_afterCuts == 2){
                    std::cout << "Electron!! slicePFPType_afterCuts = " << slicePFPType_afterCuts << ", signal = " << signal << ", highestEnergyPFP_afterCuts.truePDG = " << highestEnergyPFP_afterCuts.truePDG << ", highestEnergyPFP_afterCuts.trueInt = " << highestEnergyPFP_afterCuts.trueInt << ", highestEnergyPFP_afterCuts.trueOrigin = " << highestEnergyPFP_afterCuts.trueOrigin << std::endl;
                }

                
                if(sliceInteractionType == 1 && signal == 1){
                    /*
                    std::cout << "=============" << std::endl;
                    std::cout << "Signal Passes" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                */
                } else if(sliceInteractionType != 1){
                    /*
                    std::cout << "=============" << std::endl;
                    std::cout << "Background Passes" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                    */
                } else{
                    /*
                    std::cout << "=============" << std::endl;
                    std::cout << "signal = " << signal << ", true int = " << highestEnergyPFP_afterCuts.trueInt << std::endl;
                    std::cout << "=============" << std::endl;
                    */
                }

                fillHistogram(&sliceCompletenessAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceCompleteness->at(slice), &weights);
                fillSplitIntHistogram(&sliceCompletenessAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceCompleteness->at(slice), &weights);
                fillSplitPFPHistogram(&sliceCompletenessAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, reco_sliceCompleteness->at(slice), &weights);

                fillHistogram(&sliceCRUMBSAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
                fillSplitIntHistogram(&sliceCRUMBSAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceScore->at(slice), &weights);
                fillSplitPFPHistogram(&sliceCRUMBSAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, reco_sliceScore->at(slice), &weights);
                
                fillHistogram(&slicePurityAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_slicePurity->at(slice), &weights);
                fillSplitIntHistogram(&slicePurityAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_slicePurity->at(slice), &weights);
                fillSplitPFPHistogram(&slicePurityAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, reco_slicePurity->at(slice), &weights);
                
                fillHistogram(&sliceNumRecoNeutAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
                fillSplitIntHistogram(&sliceNumRecoNeutAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numRecoNeutrinos, &weights);
                fillSplitPFPHistogram(&sliceNumRecoNeutAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, numRecoNeutrinos, &weights);
                
                fillHistogram(&sliceNumPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice_afterCuts, &weights);
                fillSplitIntHistogram(&sliceNumPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPFPsSlice_afterCuts, &weights);
                fillSplitPFPHistogram(&sliceNumPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, numPFPsSlice_afterCuts, &weights);

                if(sliceInteractionType == 1 && signal == 1){
                    numSignalSlicesPFPCount++;

                    if(numPFPsSlice_afterCuts > 1){
                        numSignalSlicesPFPCountMoreThan1++;
                        std::cout << "=============" << std::endl;
                        std::cout << "Signal Has More Than 1 PFP in Slice" << std::endl;
                        std::cout << "eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                        std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", PFP completeness = " << highestEnergyPFP_afterCuts.completeness << ", PFP purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true PFP vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true PFP end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                        std::cout << "PFP Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), PFP Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << ", ETheta2 = " << (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta) << std::endl;
                        std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                        std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                        std::cout << "Slice CRUMBS Score = " << reco_sliceScore->at(slice) << ", reco neutrino vertex = (" << recoVX << ", " << recoVY << ", " << recoVZ << "), number of PFPs in slice = " << numPFPsSlice_afterCuts << std::endl;
                        // Loop through the PFPs
                        std::cout << "" << std::endl;
    
                        double numPFPSNuEOrigin = 0;
                        double numPFPsNuEOriginPrimary = 0;
                        int pfpCount = 0;
                        for(size_t pfpLoop = 0; pfpLoop < reco_particlePDG->size(); ++pfpLoop){
                            if(reco_particleSliceID->at(pfpLoop) == reco_sliceID->at(slice)){
                                if(reco_particleClearCosmic->at(pfpLoop) == 0){
                                    pfpCount++;
                                    std::cout << "PFP " << pfpCount << " (ID = " << reco_particleID->at(pfpLoop) << "): True Int = " << reco_particleTrueInteractionType->at(pfpLoop) << ", True Origin = " << reco_particleTrueOrigin->at(pfpLoop) << ", True PDG = " << reco_particleTruePDG->at(pfpLoop) << ", True Vertex = (" << reco_particleTrueVX->at(pfpLoop) << ", " << reco_particleTrueVY->at(pfpLoop) << ", " << reco_particleTrueVZ->at(pfpLoop) << "), True End = (" << reco_particleTrueEndX->at(pfpLoop) << ", " << reco_particleTrueEndZ->at(pfpLoop) << "), Primary = " << reco_particleIsPrimary->at(pfpLoop) << ", PFP completeness = " << reco_particleCompleteness->at(pfpLoop) << ", PFP Purity = " << reco_particlePurity->at(pfpLoop) << ", energy = " << reco_particleBestPlaneEnergy->at(pfpLoop) << ", theta = " << reco_particleTheta->at(pfpLoop) << ", vertex = (" << reco_particleVX->at(pfpLoop) << ", " << reco_particleVY->at(pfpLoop) << ", " << reco_particleVZ->at(pfpLoop) << "), num hits = " << reco_particleNumHits->at(pfpLoop) << ", Direction = (" << reco_particleDX->at(pfpLoop) << ", " << reco_particleDY->at(pfpLoop) << ", " << reco_particleDZ->at(pfpLoop) << ")" << std::endl;
                                    if(reco_particleTrueInteractionType->at(pfpLoop) == 1098){
                                        numPFPSNuEOrigin++;
                                        if(reco_particleIsPrimary->at(pfpLoop) == 1){
                                            numPFPsNuEOriginPrimary++;
                                        }
                                    }
                                }
                            }
                        }

                        if(numPFPSNuEOrigin > 1){
                            numSignalSlicesPFPCountMoreThan1NuE++;
                            if(numPFPsNuEOriginPrimary > 1){
                                numSignalSlicesPFPCountMoreThan1NuEPrimaryPlus++;
                            }
                        }
                        std::cout << "=============" << std::endl;
                    }

                }

                fillHistogram(&sliceNumPrimaryPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice_afterCuts, &weights);
                fillSplitIntHistogram(&sliceNumPrimaryPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, numPrimaryPFPsSlice_afterCuts, &weights);
                fillSplitPFPHistogram(&sliceNumPrimaryPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, numPrimaryPFPsSlice_afterCuts, &weights);
                
                fillHistogram(&sliceFracHitsInPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)), &weights);
                fillSplitIntHistogram(&sliceFracHitsInPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)), &weights);
                fillSplitPFPHistogram(&sliceFracHitsInPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)), &weights);
                
                fillHistogram(&sliceFracHitsInHighestEnergyPFPsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
                fillSplitIntHistogram(&sliceFracHitsInHighestEnergyPFPsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
                fillSplitPFPHistogram(&sliceFracHitsInHighestEnergyPFPsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, (highestEnergyPFP_afterCuts.numHits/reco_sliceNumHits->at(slice)), &weights);
                
                fillHistogram(&ERecoSumThetaRecoAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (summedEnergy_afterCuts * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
                fillSplitIntHistogram(&ERecoSumThetaRecoAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (summedEnergy_afterCuts * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
                fillSplitPFPHistogram(&ERecoSumThetaRecoAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, (summedEnergy_afterCuts * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
                
                fillHistogram(&ERecoHighestThetaRecoAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
                fillSplitIntHistogram(&ERecoHighestThetaRecoAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
                fillSplitPFPHistogram(&ERecoHighestThetaRecoAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, (highestEnergyPFP_afterCuts.energy * highestEnergyPFP_afterCuts.theta * highestEnergyPFP_afterCuts.theta), &weights);
                
                fillHistogram(&dEdxAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
                fillSplitIntHistogram(&dEdxAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
                fillSplitPFPHistogram(&dEdxAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.bestPlanedEdx, &weights);
                
                fillHistogram(&razzledPDG11AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
                fillSplitIntHistogram(&razzledPDG11AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
                fillSplitPFPHistogram(&razzledPDG11AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.razzledPDG11, &weights);
               
                if(printLowRazzledElectronElectrons == 1 && (highestEnergyPFP_afterCuts.razzledPDG11 < 0.06) && (slicePFPType_afterCuts == 2)){
                    // This is a electron with a low razzled electron score
                    std::cout << "=============" << std::endl;
                    std::cout << "Electron with a low razzled electron score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
               
                if(printHighRazzledElectronElectrons == 1 && (highestEnergyPFP_afterCuts.razzledPDG11 > 0.7) && (slicePFPType_afterCuts == 2)){
                    // This is a electron with a high razzled electron score
                    std::cout << "=============" << std::endl;
                    std::cout << "Electron with a high razzled electron score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                fillHistogram(&razzledPDG13AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG13, &weights);
                fillSplitIntHistogram(&razzledPDG13AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG13, &weights);
                fillSplitPFPHistogram(&razzledPDG13AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.razzledPDG13, &weights);
               
                if(printLowRazzledMuonMuons == 1 && (highestEnergyPFP_afterCuts.razzledPDG13 < 0.06) && (slicePFPType_afterCuts == 4)){
                    // This is a muon with a low razzled muon score
                    std::cout << "=============" << std::endl;
                    std::cout << "Muon with a low razzled muon score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                if(printHighRazzledMuonMuons == 1 && (highestEnergyPFP_afterCuts.razzledPDG13 > 0.7) && (slicePFPType_afterCuts == 4)){
                    // This is a muon with a high razzled muon score
                    std::cout << "=============" << std::endl;
                    std::cout << "Muon with a high razzled muon score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                fillHistogram(&razzledPDG22AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG22, &weights);
                fillSplitIntHistogram(&razzledPDG22AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG22, &weights);
                fillSplitPFPHistogram(&razzledPDG22AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.razzledPDG22, &weights);
               
                if(printLowRazzledPhotonPhotons == 1 && (highestEnergyPFP_afterCuts.razzledPDG22 < 0.06) && (slicePFPType_afterCuts == 7)){
                    // This is a photon with a low razzled photon score
                    std::cout << "=============" << std::endl;
                    std::cout << "Photon with a low razzled photon score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                if(printHighRazzledPhotonPhotons == 1 && (highestEnergyPFP_afterCuts.razzledPDG22 < 0.06) && (slicePFPType_afterCuts == 7)){
                    // This is a photon with a low razzled photon score
                    std::cout << "=============" << std::endl;
                    std::cout << "Photon with a high razzled photon score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                if(printHighRazzledElectronPhotons == 1 && (highestEnergyPFP_afterCuts.razzledPDG11 > 0.7) && (slicePFPType_afterCuts == 7)){
                    // This is a photon with a high razzled electron score
                    std::cout << "=============" << std::endl;
                    std::cout << "Photon with a high razzled electron score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                fillHistogram(&razzledPDG211AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
                fillSplitIntHistogram(&razzledPDG211AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
                fillSplitPFPHistogram(&razzledPDG211AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.razzledPDG211, &weights);
                
                if(printLowRazzledChargedPiChargedPis == 1 && (highestEnergyPFP_afterCuts.razzledPDG211 < 0.06) && (slicePFPType_afterCuts == 6)){
                    // This is a charged pion with a low razzled charged pion score
                    std::cout << "=============" << std::endl;
                    std::cout << "Charged pion with a low razzled charged pion score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }
                
                if(printHighRazzledChargedPiChargedPis == 1 && (highestEnergyPFP_afterCuts.razzledPDG211 > 0.7) && (slicePFPType_afterCuts == 6)){
                    // This is a charged pion with a high razzled charged pion score
                    std::cout << "=============" << std::endl;
                    std::cout << "Charged pion with a high razzled charged pion score" << std::endl;
                    std::cout << "Survives Cuts: eventID = " << eventID << ", runID = " << runID << ", subRunID = " << subRunID << ", signal = " << signal << ", slice ID = " << reco_sliceID->at(slice) << std::endl;
                    std::cout << "Highest energy PFP ID = " << highestEnergyPFP_afterCuts.PFPID << ", energy = " << highestEnergyPFP_afterCuts.energy << ", theta = " << highestEnergyPFP_afterCuts.theta << ", completeness = " << highestEnergyPFP_afterCuts.completeness << ", purity = " << highestEnergyPFP_afterCuts.purity << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true origin = " << highestEnergyPFP_afterCuts.trueOrigin << ", true int = " << highestEnergyPFP_afterCuts.trueInt << ", true PDG = " << highestEnergyPFP_afterCuts.truePDG << ", true vertex = (" << highestEnergyPFP_afterCuts.trueVX << ", " << highestEnergyPFP_afterCuts.trueVY << ", " << highestEnergyPFP_afterCuts.trueVZ << "), true end = (" << highestEnergyPFP_afterCuts.trueEndX << ", " << highestEnergyPFP_afterCuts.trueEndY << ", " << highestEnergyPFP_afterCuts.trueEndZ << "), true length = " << highestEnergyPFP_afterCuts.trueLength << std::endl;
                    std::cout << "Direction = (" << highestEnergyPFP_afterCuts.dx << ", " << highestEnergyPFP_afterCuts.dy << ", " << highestEnergyPFP_afterCuts.dz << "), Reco Vertex = (" << highestEnergyPFP_afterCuts.vx << ", " << highestEnergyPFP_afterCuts.vy << ", " << highestEnergyPFP_afterCuts.vz << "), trackscore = " << highestEnergyPFP_afterCuts.trackscore << ", primary = " << highestEnergyPFP_afterCuts.primary << ", shower length = " << highestEnergyPFP_afterCuts.showerLength << ", shower open angle = " << highestEnergyPFP_afterCuts.showerOpenAngle << ", num hits = " << highestEnergyPFP_afterCuts.numHits << ", dE/dx = " << highestEnergyPFP_afterCuts.bestPlanedEdx << std::endl;
                    std::cout << "Razzled scores: electron = " << highestEnergyPFP_afterCuts.razzledPDG11 << ", muon = " << highestEnergyPFP_afterCuts.razzledPDG13 << ", photon = " << highestEnergyPFP_afterCuts.razzledPDG22 << ", charged pi = " << highestEnergyPFP_afterCuts.razzledPDG211 << ", proton = " << highestEnergyPFP_afterCuts.razzledPDG2212 << std::endl;
                    std::cout << "Number of Hits in Slice = " << reco_sliceNumHits->at(slice) << ", Fraction of Hits Contained in PFPs = " << (numHitsInPFPs_afterCuts/reco_sliceNumHits->at(slice)) << std::endl;
                    std::cout << "=============" << std::endl;
                }

                fillHistogram(&razzledPDG2212AfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.razzledPDG2212, &weights);
                fillSplitIntHistogram(&razzledPDG2212AfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.razzledPDG2212, &weights);
                fillSplitPFPHistogram(&razzledPDG2212AfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.razzledPDG2212, &weights);

                fillHistogram(&pfpCompletenessAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.completeness, &weights);
                fillSplitIntHistogram(&pfpCompletenessAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.completeness, &weights);
                fillSplitPFPHistogram(&pfpCompletenessAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.completeness, &weights);
                
                fillHistogram(&pfpPurityAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.purity, &weights);
                fillSplitIntHistogram(&pfpPurityAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.purity, &weights);
                fillSplitPFPHistogram(&pfpPurityAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.purity, &weights);
                
                fillHistogram(&pfpNumHitsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, highestEnergyPFP_afterCuts.numHits, &weights);
                fillSplitIntHistogram(&pfpNumHitsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, highestEnergyPFP_afterCuts.numHits, &weights);
                fillSplitPFPHistogram(&pfpNumHitsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, highestEnergyPFP_afterCuts.numHits, &weights);
                
                fillHistogram(&sliceNumHitsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceNumHits->at(slice), &weights);
                fillSplitIntHistogram(&sliceNumHitsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, reco_sliceNumHits->at(slice), &weights);
                fillSplitPFPHistogram(&sliceNumHitsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, reco_sliceNumHits->at(slice), &weights);
                
                if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.energy != -999999) fillHistogram(&trueRecoilElectronEnergyAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoilElectron.energy, &weights);
                if(signal == 1 && sliceCategoryPlottingMacro == 1 && recoilElectron.angle != -999999) fillHistogram(&trueRecoilElectronAngleAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, (recoilElectron.angle*TMath::RadToDeg()), &weights);
                
                fillHistogram(&recoVXAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVX, &weights);
                
                fillHistogram(&recoVYAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVY, &weights);
                
                fillHistogram(&recoVZAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVZ, &weights);
                
                fillHistogram(&recoVXSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXSmallerBinsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVX, &weights);
                
                fillHistogram(&recoVYSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYSmallerBinsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVY, &weights);
                
                fillHistogram(&recoVZSmallerBinsAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZSmallerBinsAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZSmallerBinsAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVZ, &weights);
                
                fillHistogram(&recoVXLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXLowAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVX, &weights);
                
                fillHistogram(&recoVYLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYLowAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVY, &weights);
                
                fillHistogram(&recoVZLowAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZLowAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZLowAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVZ, &weights);
                
                fillHistogram(&recoVXHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVX, &weights);
                fillSplitIntHistogram(&recoVXHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVX, &weights);
                fillSplitPFPHistogram(&recoVXHighAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVX, &weights);
                
                fillHistogram(&recoVYHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVY, &weights);
                fillSplitIntHistogram(&recoVYHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVY, &weights);
                fillSplitPFPHistogram(&recoVYHighAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVY, &weights);
                
                fillHistogram(&recoVZHighAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, recoVZ, &weights);
                fillSplitIntHistogram(&recoVZHighAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, recoVZ, &weights);
                fillSplitPFPHistogram(&recoVZHighAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_afterCuts, recoVZ, &weights);
                
                if((sliceCategoryPlottingMacro == 1 || sliceCategoryPlottingMacro == 2) && signal == 1){
                    fillHistogram(&energyAsymmetryAfterCuts, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy), &weights);
                    fillSplitIntHistogram(&energyAsymmetryAfterCuts_splitDLNuE, DLCurrent, signal, sliceInteractionType, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy), &weights);
                    fillSplitPFPHistogram(&energyAsymmetryAfterCuts_splitPFPDLNuE, DLCurrent, signal, slicePFPType_beforeCuts, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy), &weights);
                
                    xCoordEnergyAsymmetry->Fill(recoVX, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                    xCoordAngleDifference->Fill(recoVX, angleDifference_afterCuts);
                
                    yCoordEnergyAsymmetry->Fill(recoVY, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                    yCoordAngleDifference->Fill(recoVY, angleDifference_afterCuts);
                
                    zCoordEnergyAsymmetry->Fill(recoVZ, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                    zCoordAngleDifference->Fill(recoVZ, angleDifference_afterCuts);

                    if(recoVX >= xMin && recoVX <= xMin+30){
                        xCoordEnergyAsymmetry_low->Fill(recoVX, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                        xCoordAngleDifference_low->Fill(recoVX, angleDifference_afterCuts);
                    } else if(recoVX <= xMax && recoVX >= xMax-30){
                        xCoordEnergyAsymmetry_high->Fill(recoVX, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                        xCoordAngleDifference_high->Fill(recoVX, angleDifference_afterCuts);
                    }

                    if(recoVY >= yMin && recoVY <= yMin+30){
                        yCoordEnergyAsymmetry_low->Fill(recoVY, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                        yCoordAngleDifference_low->Fill(recoVY, angleDifference_afterCuts);
                    } else if(recoVY <= yMax && recoVY >= yMax-30){
                        yCoordEnergyAsymmetry_high->Fill(recoVY, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                        yCoordAngleDifference_high->Fill(recoVY, angleDifference_afterCuts);
                    }

                    if(recoVZ >= zMin && recoVZ <= zMin+30){
                        zCoordEnergyAsymmetry_low->Fill(recoVZ, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                        zCoordAngleDifference_low->Fill(recoVZ, angleDifference_afterCuts);
                    } else if(recoVZ <= zMax && recoVZ >= zMax-110){
                        zCoordEnergyAsymmetry_high->Fill(recoVZ, ((recoilElectron.energy - highestEnergyPFP_afterCuts.energy) /recoilElectron.energy));
                        zCoordAngleDifference_high->Fill(recoVZ, angleDifference_afterCuts);
                    }

                }
            }

            

        }


    }

    int drawLine = 1;
    int left = 0;
    int right = 1;

    styleDrawAll(sliceCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCompletenessAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceCompletenessAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceCompleteness_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&sliceCompletenessBeforeCuts, &sliceCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "sliceCompletenessHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&sliceCompletenessBeforeCuts, &sliceCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "sliceCompletenessLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);

    styleDrawAll(sliceCRUMBSBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSBeforeCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceCRUMBSAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceCRUMBSAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBS_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    //std::cout << "------------------------------------------------------------------------------ CRUMBS ------------------------------------------------------------------------------" << std::endl;
    efficiency(&sliceCRUMBSBeforeCuts, &sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBSHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    //std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    efficiency(&sliceCRUMBSBeforeCuts, &sliceCRUMBSAfterCuts, 999, 999, 999, 999, (base_path + "sliceCRUMBSLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(slicePurityBeforeCuts, 999, 999, 999, 999, (base_path + "slicePurity_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(slicePurityBeforeCuts, 999, 999, 999, 999, (base_path + "slicePurity_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(slicePurityAfterCuts, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(slicePurityAfterCuts, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(slicePurityAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(slicePurityAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "slicePurity_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&slicePurityBeforeCuts, &slicePurityAfterCuts, 999, 999, 999, 999, (base_path + "slicePurityHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&slicePurityBeforeCuts, &slicePurityAfterCuts, 999, 999, 999, 999, (base_path + "slicePurityLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(sliceNumRecoNeutBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumRecoNeutBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceNumRecoNeutAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumRecoNeutAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumRecoNeutAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceNumRecoNeutAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&sliceNumRecoNeutBeforeCuts, &sliceNumRecoNeutAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeutHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&sliceNumRecoNeutBeforeCuts, &sliceNumRecoNeutAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumRecoNeutLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);

    styleDrawAll(sliceNumPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceNumPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPFPsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceNumPFPsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&sliceNumPFPsBeforeCuts, &sliceNumPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&sliceNumPFPsBeforeCuts, &sliceNumPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPFPsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);

    styleDrawAll(sliceNumPrimaryPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceNumPrimaryPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumPrimaryPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumPrimaryPFPsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceNumPrimaryPFPsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&sliceNumPrimaryPFPsBeforeCuts, &sliceNumPrimaryPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&sliceNumPrimaryPFPsBeforeCuts, &sliceNumPrimaryPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);

    styleDrawAll(sliceFracHitsInPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPs_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceFracHitsInPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPs_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceFracHitsInPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPs_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceFracHitsInPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPs_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceFracHitsInPFPsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPs_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceFracHitsInPFPsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPs_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&sliceFracHitsInPFPsBeforeCuts, &sliceFracHitsInPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&sliceFracHitsInPFPsBeforeCuts, &sliceFracHitsInPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInPFPsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);

    styleDrawAll(sliceFracHitsInHighestEnergyPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPs_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceFracHitsInHighestEnergyPFPsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPs_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceFracHitsInHighestEnergyPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPs_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceFracHitsInHighestEnergyPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPs_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceFracHitsInHighestEnergyPFPsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPs_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceFracHitsInHighestEnergyPFPsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPs_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&sliceFracHitsInHighestEnergyPFPsBeforeCuts, &sliceFracHitsInHighestEnergyPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&sliceFracHitsInHighestEnergyPFPsBeforeCuts, &sliceFracHitsInHighestEnergyPFPsAfterCuts, 999, 999, 999, 999, (base_path + "sliceFracHitsInHighestEnergyPFPsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(ERecoSumThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoSumThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(ERecoSumThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoSumThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoSumThetaRecoAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(ERecoSumThetaRecoAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&ERecoSumThetaRecoBeforeCuts, &ERecoSumThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&ERecoSumThetaRecoBeforeCuts, &ERecoSumThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);

    styleDrawAll(ERecoHighestThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoBeforeCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(ERecoHighestThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(ERecoHighestThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(ERecoHighestThetaRecoAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(ERecoHighestThetaRecoAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&ERecoHighestThetaRecoBeforeCuts, &ERecoHighestThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&ERecoHighestThetaRecoBeforeCuts, &ERecoHighestThetaRecoAfterCuts, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(dEdxBeforeCuts, 999, 999, 999, 999, (base_path + "dEdx_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(dEdxBeforeCuts, 999, 999, 999, 999, (base_path + "dEdx_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(dEdxAfterCuts, 999, 999, 999, 999, (base_path + "dEdx_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(dEdxAfterCuts, 999, 999, 999, 999, (base_path + "dEdx_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(dEdxAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "dEdx_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(dEdxAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "dEdx_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&dEdxBeforeCuts, &dEdxAfterCuts, 999, 999, 999, 999, (base_path + "dEdxHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&dEdxBeforeCuts, &dEdxAfterCuts, 999, 999, 999, 999, (base_path + "dEdxLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(razzledPDG11BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG11BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(razzledPDG11AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG11AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG11AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(razzledPDG11AfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG11_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&razzledPDG11BeforeCuts, &razzledPDG11AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG11High").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&razzledPDG11BeforeCuts, &razzledPDG11AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG11Low").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(razzledPDG13BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG13BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(razzledPDG13AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG13AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG13AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(razzledPDG13AfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG13_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&razzledPDG13BeforeCuts, &razzledPDG13AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG13High").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&razzledPDG13BeforeCuts, &razzledPDG13AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG13Low").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(razzledPDG22BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG22BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(razzledPDG22AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG22AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG22AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(razzledPDG22AfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG22_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&razzledPDG22BeforeCuts, &razzledPDG22AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG22High").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&razzledPDG22BeforeCuts, &razzledPDG22AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG22Low").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(razzledPDG211BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG211BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(razzledPDG211AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG211AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG211AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(razzledPDG211AfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG211_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&razzledPDG211BeforeCuts, &razzledPDG211AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG211High").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&razzledPDG211BeforeCuts, &razzledPDG211AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG211Low").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(razzledPDG2212BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG2212BeforeCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(razzledPDG2212AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(razzledPDG2212AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(razzledPDG2212AfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(razzledPDG2212AfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "razzledPDG2212_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&razzledPDG2212BeforeCuts, &razzledPDG2212AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212High").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&razzledPDG2212BeforeCuts, &razzledPDG2212AfterCuts, 999, 999, 999, 999, (base_path + "razzledPDG2212Low").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(pfpCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpCompletenessBeforeCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(pfpCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpCompletenessAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(pfpCompletenessAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "pfpCompleteness_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&pfpCompletenessBeforeCuts, &pfpCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "pfpCompletenessHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&pfpCompletenessBeforeCuts, &pfpCompletenessAfterCuts, 999, 999, 999, 999, (base_path + "pfpCompletenessLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(pfpPurityBeforeCuts, 999, 999, 999, 999, (base_path + "pfpPurity_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpPurityBeforeCuts, 999, 999, 999, 999, (base_path + "pfpPurity_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(pfpPurityAfterCuts, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpPurityAfterCuts, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpPurityAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(pfpPurityAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "pfpPurity_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&pfpPurityBeforeCuts, &pfpPurityAfterCuts, 999, 999, 999, 999, (base_path + "pfpPurityHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&pfpPurityBeforeCuts, &pfpPurityAfterCuts, 999, 999, 999, 999, (base_path + "pfpPurityLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(pfpNumHitsBeforeCuts, 999, 999, 999, 999, (base_path + "pfpNumHits_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpNumHitsBeforeCuts, 999, 999, 999, 999, (base_path + "pfpNumHits_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(pfpNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "pfpNumHits_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(pfpNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "pfpNumHits_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(pfpNumHitsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "pfpNumHits_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(pfpNumHitsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "pfpNumHits_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&pfpNumHitsBeforeCuts, &pfpNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "pfpNumHitsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&pfpNumHitsBeforeCuts, &pfpNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "pfpNumHitsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(sliceNumHitsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumHits_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumHitsBeforeCuts, 999, 999, 999, 999, (base_path + "sliceNumHits_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumHits_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(sliceNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumHits_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(sliceNumHitsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumHits_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(sliceNumHitsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "sliceNumHits_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&pfpNumHitsBeforeCuts, &sliceNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumHitsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&pfpNumHitsBeforeCuts, &sliceNumHitsAfterCuts, 999, 999, 999, 999, (base_path + "sliceNumHitsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(trueRecoilElectronEnergyBeforeCuts, 999, 999, 999, 999, (base_path + "trueRecoilElectronEnergy_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(trueRecoilElectronEnergyAfterCuts, 999, 999, 999, 999, (base_path + "trueRecoilElectronEnergy_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    
    styleDrawAll(trueRecoilElectronAngleBeforeCuts, 999, 999, 999, 999, (base_path + "trueRecoilElectronAngle_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(trueRecoilElectronAngleAfterCuts, 999, 999, 999, 999, (base_path + "trueRecoilElectronAngle_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    
    styleDrawAll(recoVXBeforeCuts, 999, 999, 999, 999, (base_path + "recoVX_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXBeforeCuts, 999, 999, 999, 999, (base_path + "recoVX_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVXAfterCuts, 999, 999, 999, 999, (base_path + "recoVX_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXAfterCuts, 999, 999, 999, 999, (base_path + "recoVX_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVXAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVX_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVXAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVX_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVXBeforeCuts, &recoVXAfterCuts, 999, 999, 190, 200, (base_path + "recoVXHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVXBeforeCuts, &recoVXAfterCuts, 999, 999, -200, -190, (base_path + "recoVXLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVYBeforeCuts, 999, 999, 999, 999, (base_path + "recoVY_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYBeforeCuts, 999, 999, 999, 999, (base_path + "recoVY_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVYAfterCuts, 999, 999, 999, 999, (base_path + "recoVY_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYAfterCuts, 999, 999, 999, 999, (base_path + "recoVY_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVYAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVY_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVYAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVY_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVYBeforeCuts, &recoVYAfterCuts, 999, 999, 192, 202, (base_path + "recoVYHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVYBeforeCuts, &recoVYAfterCuts, 999, 999, -202, -192, (base_path + "recoVYLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVZBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZ_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZ_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVZAfterCuts, 999, 999, 999, 999, (base_path + "recoVZ_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZAfterCuts, 999, 999, 999, 999, (base_path + "recoVZ_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVZAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVZ_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVZAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVZ_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVZBeforeCuts, &recoVZAfterCuts, 999, 999, 486, 510, (base_path + "recoVZHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVZBeforeCuts, &recoVZAfterCuts, 999, 999, 0, 16, (base_path + "recoVZLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVXSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "recoVXSmallerBins_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "recoVXSmallerBins_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVXSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "recoVXSmallerBins_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "recoVXSmallerBins_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVXSmallerBinsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVXSmallerBins_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVXSmallerBinsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVXSmallerBins_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVXSmallerBinsBeforeCuts, &recoVXSmallerBinsAfterCuts, 999, 999, -202, -190, (base_path + "recoVXSmallerBinsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVXSmallerBinsBeforeCuts, &recoVXSmallerBinsAfterCuts, 999, 999, 190, 202, (base_path + "recoVXSmallerBinsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVYSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "recoVYSmallerBins_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "recoVYSmallerBins_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVYSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "recoVYSmallerBins_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "recoVYSmallerBins_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVYSmallerBinsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVYSmallerBins_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVYSmallerBinsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVYSmallerBins_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVYSmallerBinsBeforeCuts, &recoVYSmallerBinsAfterCuts, 999, 999, -204, -192, (base_path + "recoVYSmallerBinsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVYSmallerBinsBeforeCuts, &recoVYSmallerBinsAfterCuts, 999, 999, 192, 204, (base_path + "recoVYSmallerBinsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVZSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZSmallerBins_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZSmallerBinsBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZSmallerBins_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVZSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "recoVZSmallerBins_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZSmallerBinsAfterCuts, 999, 999, 999, 999, (base_path + "recoVZSmallerBins_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVZSmallerBinsAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVZSmallerBins_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVZSmallerBinsAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVZSmallerBins_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVZSmallerBinsBeforeCuts, &recoVZSmallerBinsAfterCuts, 999, 999, 0, 16, (base_path + "recoVZSmallerBinsHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVZSmallerBinsBeforeCuts, &recoVZSmallerBinsAfterCuts, 999, 999, 486, 510, (base_path + "recoVZSmallerBinsLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVXLowBeforeCuts, 999, 999, 999, 999, (base_path + "recoVXLow_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXLowBeforeCuts, 999, 999, 999, 999, (base_path + "recoVXLow_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVXLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVXLow_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVXLow_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVXLowAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVXLow_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVXLowAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVXLow_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVXLowBeforeCuts, &recoVXLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVXLowHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVXLowBeforeCuts, &recoVXLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVXLowLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVYLowBeforeCuts, 999, 999, 999, 999, (base_path + "recoVYLow_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYLowBeforeCuts, 999, 999, 999, 999, (base_path + "recoVYLow_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVYLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVYLow_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVYLow_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVYLowAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVYLow_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVYLowAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVYLow_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVYLowBeforeCuts, &recoVYLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVYLowHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVYLowBeforeCuts, &recoVYLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVYLowLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVZLowBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZLow_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZLowBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZLow_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVZLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVZLow_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVZLow_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVZLowAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVZLow_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVZLowAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVZLow_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVZLowBeforeCuts, &recoVZLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVZLowHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVZLowBeforeCuts, &recoVZLowAfterCuts, 999, 999, 999, 999, (base_path + "recoVZLowLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVXHighBeforeCuts, 999, 999, 999, 999, (base_path + "recoVXHigh_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXHighBeforeCuts, 999, 999, 999, 999, (base_path + "recoVXHigh_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVXHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVXHigh_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVXHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVXHigh_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVXHighAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVXHigh_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVXHighAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVXHigh_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVXHighBeforeCuts, &recoVXHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVXHighHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVXHighBeforeCuts, &recoVXHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVXHighLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVYHighBeforeCuts, 999, 999, 999, 999, (base_path + "recoVYHigh_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYHighBeforeCuts, 999, 999, 999, 999, (base_path + "recoVYHigh_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVYHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVYHigh_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVYHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVYHigh_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVYHighAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVYHigh_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVYHighAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVYHigh_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVYHighBeforeCuts, &recoVYHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVYHighHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVYHighBeforeCuts, &recoVYHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVYHighLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoVZHighBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZHigh_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZHighBeforeCuts, 999, 999, 999, 999, (base_path + "recoVZHigh_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(recoVZHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVZHigh_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoVZHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVZHigh_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(recoVZHighAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "recoVZHigh_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(recoVZHighAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "recoVZHigh_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);
    efficiency(&recoVZHighBeforeCuts, &recoVZHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVZHighHigh").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(&recoVZHighBeforeCuts, &recoVZHighAfterCuts, 999, 999, 999, 999, (base_path + "recoVZHighLow").c_str(), "topRight", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(energyAsymmetryBeforeCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_beforeCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false, true);
    styleDrawBackSig(energyAsymmetryBeforeCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_beforeCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(energyAsymmetryAfterCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false, true);
    styleDrawBackSig(energyAsymmetryAfterCuts, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts_BackSig.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawSplit(energyAsymmetryAfterCuts_splitDLNuE, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts_splitInt.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawPFPSplit(energyAsymmetryAfterCuts_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "energyAsymmetry_afterCuts_splitPDG.pdf").c_str(), "topRight", nullptr, &right, true);


    TwoDHistDraw(xCoordAngleDifference, (base_path + "xCoordAngleDifference.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifference, (base_path + "yCoordAngleDifference.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifference, (base_path + "zCoordAngleDifference.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifference_low, (base_path + "xCoordAngleDifference_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifference_low, (base_path + "yCoordAngleDifference_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifference_low, (base_path + "zCoordAngleDifference_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    
    TwoDHistDraw(xCoordAngleDifference_high, (base_path + "xCoordAngleDifference_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifference_high, (base_path + "yCoordAngleDifference_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifference_high, (base_path + "zCoordAngleDifference_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");

    TwoDHistDraw(xCoordEnergyAsymmetry, (base_path + "xCoordEnergyAsymmetry.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetry, (base_path + "yCoordEnergyAsymmetry.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetry, (base_path + "zCoordEnergyAsymmetry.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");

    TwoDHistDraw(xCoordEnergyAsymmetry_low, (base_path + "xCoordEnergyAsymmetry_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetry_low, (base_path + "yCoordEnergyAsymmetry_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetry_low, (base_path + "zCoordEnergyAsymmetry_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");

    TwoDHistDraw(xCoordEnergyAsymmetry_high, (base_path + "xCoordEnergyAsymmetry_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetry_high, (base_path + "yCoordEnergyAsymmetry_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetry_high, (base_path + "zCoordEnergyAsymmetry_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP in Slice;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");

    std::cout << "Numbers of events (DL Nu+E): BNB = " << numEvents_DLNuEBNB << ", Intime Cosmics = " << numEvents_DLNuECosmic << ", Nu+E Elastic Scatters = " << numEvents_DLNuENuE << std::endl;

    std::cout << "Cuts applied: clear cosmic = " << clearCosmicCut << ", num PFPs 0 = " << numPFPs0Cut << ", num reco neutrinos 0 = " << numRecoNeutrinosCut << ", CRUMBS = " << CRUMBSCut << ", FV = " << FVCut << std::endl;
    std::cout << "num primary PFPs 1 = " << primaryPFPCut << ", razzled 2212 = " << razzledPDG2212Cut << ", razzled 13 = " << razzledPDG13Cut << ", razzled 211 = " << razzledPDG211Cut << ", razzled 22 = " << razzledPDG22Cut << std::endl;
    std::cout << "razzled 11 = " << razzledPDG11Cut << ", dE/dx = " << dEdxCut << ", ETheta2 = " << ETheta2Cut << ", frac hits contained = " << fracHitsContainedCut << ", num hits = " << numHitsCut << std::endl;

    drawTEff(pfpNumHitsAfterCuts.nuESignal, pfpNumHitsBeforeCuts.nuESignal, 999, 999, 999, 999, (base_path + "pfpNumHitsSignalEfficiency.pdf").c_str(), "topRight", nullptr, &right);
    drawTEff(sliceNumHitsAfterCuts.nuESignal, sliceNumHitsBeforeCuts.nuESignal, 999, 999, 999, 999, (base_path + "sliceNumHitsSignalEfficiency.pdf").c_str(), "topRight", nullptr, &right);
    drawTEff(trueRecoilElectronEnergyAfterCuts.nuESignal, trueRecoilElectronEnergyBeforeCuts.nuESignal, 999, 999, 999, 999, (base_path + "trueRecoilElectronEnergySignalEfficiency.pdf").c_str(), "topRight", nullptr, &right);
    drawTEff(trueRecoilElectronAngleAfterCuts.nuESignal, trueRecoilElectronAngleBeforeCuts.nuESignal, 999, 999, 999, 999, (base_path + "trueRecoilElectronAngleSignalEfficiency.pdf").c_str(), "topRight", nullptr, &right);


    std::ofstream out_file(txtFileName, std::ios::app);
    if(out_file.is_open()){
        out_file << "==================" << std::endl;
        out_file << "Cuts applied: clear cosmic = " << clearCosmicCut << ", num PFPs 0 = " << numPFPs0Cut << ", num reco neutrinos 0 = " << numRecoNeutrinosCut << ", CRUMBS = " << CRUMBSCut << ", FV = " << FVCut << std::endl;
        out_file << "num primary PFPs 1 = " << primaryPFPCut << ", razzled 2212 = " << razzledPDG2212Cut << ", razzled 13 = " << razzledPDG13Cut << ", razzled 211 = " << razzledPDG211Cut << ", razzled 22 = " << razzledPDG22Cut << std::endl;
        out_file << "razzled 11 = " << razzledPDG11Cut << ", dE/dx = " << dEdxCut << ", ETheta2 = " << ETheta2Cut << ", frac hits contained = " << fracHitsContainedCut << ", num hits = " << numHitsCut << std::endl;
        out_file << "==================" << std::endl;
        out_file.close();
    } else{
        std::cerr << "Error: couldn't open txt file" << std::endl;
    }

    std::ofstream out_tablefile(tableFileName, std::ios::app);
    if(out_tablefile.is_open()){
        // DL Nu+E Table 
        out_tablefile << "=========== DL Nu+E Vertexing ===========" << std::endl;
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{|c|c|c|c|c|c|}" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\epsilon$ (\\%)} & \\textbf{$\\rho$ (\\%)} & \\textbf{$\\epsilon\\rho$}& Signal Left & Background Left \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << std::defaultfloat << std::setprecision(7) << "No Cut & " << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.signal/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsBeforeCuts_DLNuE.signal/(eventsBeforeCuts_DLNuE.signal+eventsBeforeCuts_DLNuE.background) << " & " << (eventsBeforeCuts_DLNuE.signal/(eventsBeforeCuts_DLNuE.signal+eventsBeforeCuts_DLNuE.background))*(eventsBeforeCuts_DLNuE.signal/eventsBeforeCuts_DLNuE.signal) << " & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.signal << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.signal/eventsBeforeCuts_DLNuE.signal << "\\%)" << std::fixed << std::setprecision(0) << " & " << eventsBeforeCuts_DLNuE.background << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.background/eventsBeforeCuts_DLNuE.background << "\\%)" << " \\\\ " << std::endl;
        out_tablefile << "\\hline" << std::endl;
       
        if(clearCosmicCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Remove Clear Cosmic PFPs & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/(eventsAfterCuts_DLNuE.clearCosmicsSig+eventsAfterCuts_DLNuE.clearCosmicsBack) << " & " << (eventsAfterCuts_DLNuE.clearCosmicsSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.clearCosmicsSig/(eventsAfterCuts_DLNuE.clearCosmicsSig+eventsAfterCuts_DLNuE.clearCosmicsBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numPFPs0Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "PFPs in Slice != 0 & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/(eventsAfterCuts_DLNuE.numPFPs0Sig+eventsAfterCuts_DLNuE.numPFPs0Back) << " & " << (eventsAfterCuts_DLNuE.numPFPs0Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.numPFPs0Sig/(eventsAfterCuts_DLNuE.numPFPs0Sig+eventsAfterCuts_DLNuE.numPFPs0Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numPFPs0Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0Sig/eventsBeforeCuts_DLNuE.signal << std::fixed << std::setprecision(0) << "\\%) & " << eventsAfterCuts_DLNuE.numPFPs0Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numRecoNeutrinosCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "1 Reco Neutrino in Slice & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/(eventsAfterCuts_DLNuE.numRecoNeut0Sig+eventsAfterCuts_DLNuE.numRecoNeut0Back) << " & " << (eventsAfterCuts_DLNuE.numRecoNeut0Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.numRecoNeut0Sig/(eventsAfterCuts_DLNuE.numRecoNeut0Sig+eventsAfterCuts_DLNuE.numRecoNeut0Back)) << std::fixed << std::setprecision(0) << " & " << eventsAfterCuts_DLNuE.numRecoNeut0Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numRecoNeut0Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(CRUMBSCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high << " & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.crumbsSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.crumbsSig/(eventsAfterCuts_DLNuE.crumbsSig+eventsAfterCuts_DLNuE.crumbsBack) << " & " << (eventsAfterCuts_DLNuE.crumbsSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.crumbsSig/(eventsAfterCuts_DLNuE.crumbsSig+eventsAfterCuts_DLNuE.crumbsBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
       
        if(FVCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "FV Cut & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.FVSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.FVSig/(eventsAfterCuts_DLNuE.FVSig+eventsAfterCuts_DLNuE.FVBack) << " & " << (eventsAfterCuts_DLNuE.FVSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.FVSig/(eventsAfterCuts_DLNuE.FVSig+eventsAfterCuts_DLNuE.FVBack)) << std::fixed << std::setprecision(0) << " & " << eventsAfterCuts_DLNuE.FVSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.FVBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(primaryPFPCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Primary PFPs in Slice = " << primaryPFPCutValue << " & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.primaryPFPSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.primaryPFPSig/(eventsAfterCuts_DLNuE.primaryPFPSig+eventsAfterCuts_DLNuE.primaryPFPBack) << " & " << (eventsAfterCuts_DLNuE.primaryPFPSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.primaryPFPSig/(eventsAfterCuts_DLNuE.primaryPFPSig+eventsAfterCuts_DLNuE.primaryPFPBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(ETheta2Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP)} $\\leq$ " << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$ & " << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.ETheta2Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.ETheta2Sig/(eventsAfterCuts_DLNuE.ETheta2Sig+eventsAfterCuts_DLNuE.ETheta2Back) << " & " << (eventsAfterCuts_DLNuE.ETheta2Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.ETheta2Sig/(eventsAfterCuts_DLNuE.ETheta2Sig+eventsAfterCuts_DLNuE.ETheta2Back)) << std::fixed << std::setprecision(0) << " & " << eventsAfterCuts_DLNuE.ETheta2Sig << " ("  << std::defaultfloat << std::setprecision(4) << 100*eventsAfterCuts_DLNuE.ETheta2Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.ETheta2Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG2212Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Proton Score $\\leq$ " << razzled2212High_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled2212Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled2212Sig/(eventsAfterCuts_DLNuE.razzled2212Sig+eventsAfterCuts_DLNuE.razzled2212Back) << " & " << (eventsAfterCuts_DLNuE.razzled2212Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled2212Sig/(eventsAfterCuts_DLNuE.razzled2212Sig+eventsAfterCuts_DLNuE.razzled2212Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled2212Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled2212Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG13Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Muon Score $\\leq$ " << razzled13High_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled13Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled13Sig/(eventsAfterCuts_DLNuE.razzled13Sig+eventsAfterCuts_DLNuE.razzled13Back) << " & " << (eventsAfterCuts_DLNuE.razzled13Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled13Sig/(eventsAfterCuts_DLNuE.razzled13Sig+eventsAfterCuts_DLNuE.razzled13Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled13Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled13Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG211Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ " << razzled211High_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled211Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled211Sig/(eventsAfterCuts_DLNuE.razzled211Sig+eventsAfterCuts_DLNuE.razzled211Back) << " & " << (eventsAfterCuts_DLNuE.razzled211Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled211Sig/(eventsAfterCuts_DLNuE.razzled211Sig+eventsAfterCuts_DLNuE.razzled211Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG22Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Photon Score $\\leq$ " << razzled22High_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled22Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled22Sig/(eventsAfterCuts_DLNuE.razzled22Sig+eventsAfterCuts_DLNuE.razzled22Back) << " & " << (eventsAfterCuts_DLNuE.razzled22Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled22Sig/(eventsAfterCuts_DLNuE.razzled22Sig+eventsAfterCuts_DLNuE.razzled22Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled22Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled22Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG11Cut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Electron Score $\\geq$ " << razzled11Low_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.razzled11Sig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.razzled11Sig/(eventsAfterCuts_DLNuE.razzled11Sig+eventsAfterCuts_DLNuE.razzled11Back) << " & " << (eventsAfterCuts_DLNuE.razzled11Sig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.razzled11Sig/(eventsAfterCuts_DLNuE.razzled11Sig+eventsAfterCuts_DLNuE.razzled11Back)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11Sig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11Sig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11Back << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11Back/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(dEdxCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has " << dEdxLow_highestEnergyPFP << " MeV cm^{-1} $\\leq$ dE/dx $\\leq$ " << dEdxHigh_highestEnergyPFP << std::defaultfloat << std::setprecision(4) << " MeV cm^{-1} & " << 100*eventsAfterCuts_DLNuE.dEdxSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.dEdxSig/(eventsAfterCuts_DLNuE.dEdxSig+eventsAfterCuts_DLNuE.dEdxBack) << " & " << (eventsAfterCuts_DLNuE.dEdxSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.dEdxSig/(eventsAfterCuts_DLNuE.dEdxSig+eventsAfterCuts_DLNuE.dEdxBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(fracHitsContainedCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Slice has " << numContainedHits_low << " $\\leq$ Fraction of Hits Contained in PFPs $\\leq$ " << numContainedHits_high << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.fracHitsContainedSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.fracHitsContainedSig/(eventsAfterCuts_DLNuE.fracHitsContainedSig+eventsAfterCuts_DLNuE.fracHitsContainedBack) << " & " << (eventsAfterCuts_DLNuE.fracHitsContainedSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.fracHitsContainedSig/(eventsAfterCuts_DLNuE.fracHitsContainedSig+eventsAfterCuts_DLNuE.fracHitsContainedBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.fracHitsContainedSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.fracHitsContainedBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(numHitsCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << numHits_low << " $\\leq$ Number of Hits in Highest Energy PFP in Slice $\\leq$ " << numHits_high << std::defaultfloat << std::setprecision(4) << " & " << 100*eventsAfterCuts_DLNuE.numHitsSig/eventsBeforeCuts_DLNuE.signal << " & " << 100*eventsAfterCuts_DLNuE.numHitsSig/(eventsAfterCuts_DLNuE.numHitsSig+eventsAfterCuts_DLNuE.numHitsBack) << " & " << (eventsAfterCuts_DLNuE.numHitsSig/eventsBeforeCuts_DLNuE.signal)*(eventsAfterCuts_DLNuE.numHitsSig/(eventsAfterCuts_DLNuE.numHitsSig+eventsAfterCuts_DLNuE.numHitsBack)) << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numHitsSig << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsSig/eventsBeforeCuts_DLNuE.signal << "\\%) & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numHitsBack << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsBack/eventsBeforeCuts_DLNuE.background << "\\%) \\\\ " << std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;
       
        out_tablefile << "" << std::endl;
        out_tablefile << "" << std::endl;
        out_tablefile << "" << std::endl;
        // ======================================== 
        // Put split interaction table here
        
        out_tablefile << "\\begin{table}[h!]" << std::endl;
        out_tablefile << "\\centering" << std::endl;
        out_tablefile << "\\resizebox{\\textwidth}{!}{%" << std::endl;
        out_tablefile << "\\begin{tabular}{ |c|c|c|c|c|c|c|c|c|c| }" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\multicolumn{10}{|c|}{\\textbf{Number of Events Left}} \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "\\textbf{Cut Name} & \\textbf{$\\boldsymbol{\\nu+e}$} & \\textbf{NCN$\\boldsymbol{\\pi^0}$} & \\textbf{Other NC} & \\textbf{CC$\\boldsymbol{\\nu_\\mu}$} & \\textbf{CC$\\boldsymbol{\\nu_e}$} & \\textbf{Dirt} & \\textbf{$\\boldsymbol{\\nu+e}$ Dirt} & \\textbf{Cosmic} & \\textbf{Other}\\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        out_tablefile << "No Cut & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsBeforeCuts_DLNuE.splitInt.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << std::defaultfloat << std::setprecision(4) << "(" << 100*eventsBeforeCuts_DLNuE.splitInt.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.otherNC << " (" << 100*eventsBeforeCuts_DLNuE.splitInt.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.CCnumu << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.CCnue << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.dirt << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.nuEDirt << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.cosmic << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) << eventsBeforeCuts_DLNuE.splitInt.other << " (" << std::defaultfloat << std::setprecision(4) << 100*eventsBeforeCuts_DLNuE.splitInt.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\" << std::endl;
        out_tablefile << "\\hline" << std::endl;
        if(clearCosmicCut == 1){
            out_tablefile << "Remove Clear Cosmic PFPs & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.clearCosmicsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numPFPs0Cut == 1){
            out_tablefile << "PFPs in Slice != 0 & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numPFPs0IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numPFPs0IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(numRecoNeutrinosCut == 1){
            out_tablefile << "1 Reco Neutrino in Slice & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numRecoNeut0IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
       
        if(CRUMBSCut == 1){
            out_tablefile << std::defaultfloat << std::setprecision(7) << crumbsScoreCut_low << " $\\leq$ CRUMBS Score $\\leq$ " << crumbsScoreCut_high << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.crumbsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.crumbsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.crumbsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
       
        if(FVCut == 1){
            out_tablefile << "FV Cut & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.FVIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.FVIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.FVIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }

        if(primaryPFPCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Primary PFPs in Slice = " << primaryPFPCutValue << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.primaryPFPIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.primaryPFPIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(ETheta2Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "$\\textrm{E}\\theta^2 \\textrm{ (Highest Energy PFP)} $\\leq$ " << ETheta2High_highestEnergyPFP << "\\textrm{MeV rad}^2$ & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.ETheta2IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.ETheta2IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG2212Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Proton Score $\\leq$ " << razzled2212High_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled2212IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled2212IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled2212IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG13Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Muon Score $\\leq$ " << razzled13High_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled13IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled13IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled13IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG211Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Charged Pion Score $\\leq$ " << razzled211High_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled211IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled211IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled211IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG22Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Photon Score $\\leq$ " << razzled22High_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled22IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled22IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled22IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(razzledPDG11Cut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has Electron Score $\\geq$ " << razzled11Low_highestEnergyPFP << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.razzled11IntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.razzled11IntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.razzled11IntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(dEdxCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Highest Energy PFP in Slice has " << dEdxLow_highestEnergyPFP << " MeV cm^{-1} $\\leq$ dE/dx $\\leq$ " << dEdxHigh_highestEnergyPFP << " MeV cm^{-1} & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.dEdxIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.dEdxIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.dEdxIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(fracHitsContainedCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << "Slice has " << numContainedHits_low << " $\\leq$ Fraction of Hits Contained in PFPs $\\leq$ " << numContainedHits_high << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.fracHitsContainedIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        if(numHitsCut == 1){ 
            out_tablefile << std::defaultfloat << std::setprecision(7) << numHits_low << " $\\leq$ Number of Hits in Highest Energy PFP in Slice $\\leq$ " << numHits_high << " & " << std::fixed << std::setprecision(0) << eventsAfterCuts_DLNuE.numHitsIntSplit.nuE << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.nuE/eventsBeforeCuts_DLNuE.splitInt.nuE << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.NCNPi0 << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.NCNPi0/eventsBeforeCuts_DLNuE.splitInt.NCNPi0 << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.otherNC << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.otherNC/eventsBeforeCuts_DLNuE.splitInt.otherNC << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.CCnumu << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.CCnumu/eventsBeforeCuts_DLNuE.splitInt.CCnumu << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.CCnue << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.CCnue/eventsBeforeCuts_DLNuE.splitInt.CCnue << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.dirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.dirt/eventsBeforeCuts_DLNuE.splitInt.dirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.nuEDirt << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.nuEDirt/eventsBeforeCuts_DLNuE.splitInt.nuEDirt << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.cosmic << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.cosmic/eventsBeforeCuts_DLNuE.splitInt.cosmic << "\\%) & " << std::fixed << std::setprecision(0) <<  eventsAfterCuts_DLNuE.numHitsIntSplit.other << std::defaultfloat << std::setprecision(4) << " (" << 100*eventsAfterCuts_DLNuE.numHitsIntSplit.other/eventsBeforeCuts_DLNuE.splitInt.other << "\\%) \\\\"<< std::endl;
            out_tablefile << "\\hline" << std::endl;
        }
        
        out_tablefile << "\\end{tabular}" << std::endl;
        out_tablefile << "}" << std::endl;
        out_tablefile << "\\end{table}" << std::endl;

        out_tablefile << "" << std::endl;
        out_tablefile << "\\newpage" << std::endl;
        out_tablefile << "" << std::endl; 

    }

    std::cout << "Number of signal slices with more than 1 PFP = " << numSignalSlicesPFPCountMoreThan1 << " out of " << numSignalSlicesPFPCount << " (" << (100*numSignalSlicesPFPCountMoreThan1/numSignalSlicesPFPCount) << "%)" << std::endl;
    std::cout << "Number of signal slices with more than 1 PFP coming from the nu+e scatter = " << numSignalSlicesPFPCountMoreThan1NuE << std::endl;
    std::cout << "Number of signal slices with more than 1 PFP coming from the nu+e scatter + more than 1 Primary PFP = " << numSignalSlicesPFPCountMoreThan1NuEPrimaryPlus << std::endl;
}
