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

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* electron;
    TH1F* nuEElectron;
    TH1F* proton;
    TH1F* muon;
    TH1F* cosmicMuon;
    TH1F* cosmicOther;
    TH1F* pi0;
    TH1F* chargedPi;
    TH1F* other;
    TH1F* nuEOther;
} splitPFPHistGroup_struct;

struct pfpCounter_struct{
    double electron = 0;
    double nuEElectron = 0;
    double proton = 0;
    double muon = 0;
    double cosmicMuon = 0;
    double cosmicOther = 0;
    double pi0 = 0;
    double chargedPi = 0;
    double other = 0;
    double nuEOther = 0;
};

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
        (TH1F*) base->Clone((baseName + "_cosmicOther").c_str()),
        (TH1F*) base->Clone((baseName + "_pi0").c_str()),
        (TH1F*) base->Clone((baseName + "_chargedPi").c_str()),
        (TH1F*) base->Clone((baseName + "_other").c_str()),
        (TH1F*) base->Clone((baseName + "_otherNuE").c_str())
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

void styleDrawPFPSplit(splitPFPHistGroup_struct hists,
                    double ymin, double ymax, double xmin, double xmax,
                    const char* filename, const std::string& legendLocation,
                    int* drawLine = nullptr, int* linePos = nullptr,
                    bool useLogScale = false){
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    if (useLogScale)
        hists.canvas->SetLogy(1);
    else
        hists.canvas->SetLogy(0);

    std::vector<TH1F*> allHists = {hists.electron, hists.nuEElectron, hists.proton, hists.muon, hists.cosmicMuon, hists.cosmicOther, hists.pi0, hists.chargedPi, hists.other, hists.nuEOther};

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

    for (auto* hist : allHists){
        if(hist){
            hist->SetStats(0);
            hist->GetXaxis()->SetTickLength(0.04);
            hist->GetYaxis()->SetTickLength(0.03);
            hist->GetXaxis()->SetTickSize(0.02);
            hist->GetYaxis()->SetTickSize(0.02);
        }
    }
    

    hists.nuEElectron->SetLineWidth(2); hists.nuEElectron->SetLineColor(TColor::GetColor("#656364"));
    hists.electron->SetLineWidth(2);    hists.electron->SetLineColor(TColor::GetColor("#578dff"));
    hists.proton->SetLineWidth(2);      hists.proton->SetLineColor(TColor::GetColor("#86c8dd"));
    hists.muon->SetLineWidth(2);        hists.muon->SetLineColor(TColor::GetColor("#adad7d"));
    hists.cosmicMuon->SetLineWidth(2);  hists.cosmicMuon->SetLineColor(TColor::GetColor("#c91f16"));
    hists.cosmicOther->SetLineWidth(2); hists.cosmicOther->SetLineColor(TColor::GetColor("#ff5e02"));
    hists.pi0->SetLineWidth(2);         hists.pi0->SetLineColor(TColor::GetColor("#1845fb"));
    hists.chargedPi->SetLineWidth(2);   hists.chargedPi->SetLineColor(TColor::GetColor("#c849a9"));
    hists.other->SetLineWidth(2);       hists.other->SetLineColor(TColor::GetColor("#ffa90e"));
    hists.nuEOther->SetLineWidth(2);    hists.nuEOther->SetLineColor(TColor::GetColor("#a96b59"));

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

    std::cout << "maxYValue = " << maxYValue << std::endl;
    double yminVal = useLogScale ? 1e-1 : 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = useLogScale ? (maxYValue * 100.0) : (maxYValue * 1.1);
        std::cout << "setting yaxis to " << yminVal << ", " << ymaxVal << std::endl;
        for (auto* hist : allHists)
            if (hist) hist->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    hists.electron->Draw("hist");
    hists.nuEElectron->Draw("histsame");
    hists.proton->Draw("histsame");
    hists.muon->Draw("histsame");
    hists.cosmicMuon->Draw("histsame");
    hists.cosmicOther->Draw("histsame");
    hists.pi0->Draw("histsame");
    hists.chargedPi->Draw("histsame");
    hists.other->Draw("histsame");
    hists.nuEOther->Draw("histsame");

    int nEntries = 10;
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
    legend->AddEntry(hists.cosmicOther, "Cosmic Other", "f");
    legend->AddEntry(hists.pi0, "Neutral Pion", "f");
    legend->AddEntry(hists.chargedPi, "Charged Pion", "f");
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
    stack->Add(hists.cosmicOther);
    stack->Add(hists.pi0);
    stack->Add(hists.chargedPi);
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
    legendStack->AddEntry(hists.cosmicOther, "Cosmic Other", "f");
    legendStack->AddEntry(hists.pi0, "Neutral Pion", "f");
    legendStack->AddEntry(hists.chargedPi, "Charged Pion", "f");
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

void styleDrawSplit(splitHistGroup_struct hists,
                    double ymin, double ymax, double xmin, double xmax,
                    const char* filename, const std::string& legendLocation,
                    int* drawLine = nullptr, int* linePos = nullptr,
                    bool useLogScale = false){
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

    std::cout << "maxYValue = " << maxYValue << std::endl;
    double yminVal = useLogScale ? 1e-1 : 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = useLogScale ? (maxYValue * 100.0) : (maxYValue * 1.1);
        std::cout << "setting yaxis to " << yminVal << ", " << ymaxVal << std::endl;
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

void styleDrawPur(purHist_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  int* drawLine = nullptr, int* linePos = nullptr,
                  bool writeMaxValues = false, const std::string& textfilename = "")
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
        std::ofstream outfile(textfilename, std::ios::app);
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
            std::cerr << "Error: could not open " << textfilename << " for writing." << std::endl;
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

TH1F* makeCumulative(const TH1F* h, bool keepRight)
{
    TH1F* hc = (TH1F*)h->Clone(Form("%s_cumulative", h->GetName()));
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

void drawEfficiencyErrors(TEfficiency* plot_BDT, TEfficiency* plot_DLUboone, TEfficiency* plot_DLNuE, const std::string& filename, double lowY, double highY, const std::string& legendLocation, bool effPurPlots, double xmin, double xmax, bool writeMaxValues = false, const std::string& textfilename = ""){
    if (!plot_BDT || !plot_DLUboone || !plot_DLNuE) {
        std::cerr << "drawEfficiency: null TEfficiency pointer\n";
        return;
    }

    double maxVal = std::max({getMaxValueEfficiency(plot_BDT, false), getMaxValueEfficiency(plot_DLUboone, false), getMaxValueEfficiency(plot_DLNuE, false)});

    TCanvas* c = new TCanvas("c_eff", "Efficiency comparison", 800, 600);
    c->SetTicks();

    plot_BDT->SetMarkerColor(TColor::GetColor("#e42536"));
    plot_BDT->SetMarkerSize(0.7); 
    plot_BDT->SetLineWidth(1);
    plot_BDT->SetLineColor(kBlack);
    plot_BDT->SetMarkerStyle(20);

    const TH1* hTotal_BDT = plot_BDT->GetTotalHistogram();
    int nBins_BDT = hTotal_BDT->GetNbinsX();
    TGraphAsymmErrors* gEff_BDT = new TGraphAsymmErrors(nBins_BDT);    

    double maxEff_BDT = 0;
    double maxEffBin_BDT = 0;
    for(int i = 1; i <= nBins_BDT; ++i){
        double xCenter = hTotal_BDT->GetXaxis()->GetBinCenter(i);
        double xErr = (hTotal_BDT->GetXaxis()->GetBinUpEdge(i) - hTotal_BDT->GetXaxis()->GetBinLowEdge(i)) / 2.0;
        
        double yEff = plot_BDT->GetEfficiency(i);
        double yErrLow  = plot_BDT->GetEfficiencyErrorLow(i);
        double yErrUp   = plot_BDT->GetEfficiencyErrorUp(i);
   
        if(yEff > maxEff_BDT){
            maxEff_BDT = yEff;
            maxEffBin_BDT = xCenter;
        }

        gEff_BDT->SetPoint(i-1, xCenter, yEff);
        gEff_BDT->SetPointError(i-1, xErr, xErr, yErrLow, yErrUp);
    }

    gEff_BDT->SetLineColor(kBlack);
    gEff_BDT->SetMarkerColor(TColor::GetColor("#e42536"));
    gEff_BDT->SetMarkerStyle(20);
    gEff_BDT->SetMarkerSize(0.7);
    gEff_BDT->SetLineWidth(1);
    if (xmin != 999) {
        gEff_BDT->GetXaxis()->SetLimits(xmin, xmax);
    }

    std::cout << "MAX VAL HERE = " << maxVal << ", *1.1 = " << maxVal*1.1 << std::endl;
    gEff_BDT->GetYaxis()->SetRangeUser(0, maxVal*1.1);

    const TH1* hAxis = plot_BDT->GetTotalHistogram();
    gEff_BDT->SetTitle(plot_BDT->GetTitle());
    gEff_BDT->GetXaxis()->SetTitle(hAxis->GetXaxis()->GetTitle());
    gEff_BDT->GetYaxis()->SetTitle(hAxis->GetYaxis()->GetTitle());
    gEff_BDT->Draw("AP");

    plot_DLUboone->SetMarkerColor(TColor::GetColor("#5790fc"));
    plot_DLUboone->SetMarkerSize(0.7); 
    plot_DLUboone->SetLineWidth(1);
    plot_DLUboone->SetLineColor(kBlack);
    plot_DLUboone->SetMarkerStyle(20);

    const TH1* hTotal_DLUboone = plot_DLUboone->GetTotalHistogram();
    int nBins_DLUboone = hTotal_DLUboone->GetNbinsX();
    TGraphAsymmErrors* gEff_DLUboone = new TGraphAsymmErrors(nBins_DLUboone);    

    double maxEff_DLUboone = 0;
    double maxEffBin_DLUboone = 0;
    for(int i = 1; i <= nBins_DLUboone; ++i){
        double xCenter = hTotal_DLUboone->GetXaxis()->GetBinCenter(i);
        double xErr = (hTotal_DLUboone->GetXaxis()->GetBinUpEdge(i) - hTotal_DLUboone->GetXaxis()->GetBinLowEdge(i)) / 2.0;
        
        double yEff = plot_DLUboone->GetEfficiency(i);
        double yErrLow  = plot_DLUboone->GetEfficiencyErrorLow(i);
        double yErrUp   = plot_DLUboone->GetEfficiencyErrorUp(i);
        
        if(yEff > maxEff_DLUboone){
            maxEff_DLUboone = yEff;
            maxEffBin_DLUboone = xCenter;
        }
   
        gEff_DLUboone->SetPoint(i-1, xCenter, yEff);
        gEff_DLUboone->SetPointError(i-1, xErr, xErr, yErrLow, yErrUp);
    }
    
    gEff_DLUboone->SetLineColor(kBlack);
    gEff_DLUboone->SetMarkerColor(TColor::GetColor("#5790fc"));
    gEff_DLUboone->SetMarkerStyle(20);
    gEff_DLUboone->SetMarkerSize(0.7);
    gEff_DLUboone->SetLineWidth(1);
    gEff_DLUboone->Draw("P SAME");

    plot_DLNuE->SetMarkerColor(TColor::GetColor("#f89c20"));
    plot_DLNuE->SetMarkerSize(0.7); 
    plot_DLNuE->SetLineWidth(1);
    plot_DLNuE->SetLineColor(kBlack);
    plot_DLNuE->SetMarkerStyle(20);
    
    const TH1* hTotal_DLNuE = plot_DLNuE->GetTotalHistogram();
    int nBins_DLNuE = hTotal_DLNuE->GetNbinsX();
    TGraphAsymmErrors* gEff_DLNuE = new TGraphAsymmErrors(nBins_DLNuE);    

    double maxEff_DLNuE = 0;
    double maxEffBin_DLNuE = 0; 
    for(int i = 1; i <= nBins_DLNuE; ++i){
        double xCenter = hTotal_DLNuE->GetXaxis()->GetBinCenter(i);
        double xErr = (hTotal_DLNuE->GetXaxis()->GetBinUpEdge(i) - hTotal_DLNuE->GetXaxis()->GetBinLowEdge(i)) / 2.0;
        
        double yEff = plot_DLNuE->GetEfficiency(i);
        double yErrLow  = plot_DLNuE->GetEfficiencyErrorLow(i);
        double yErrUp   = plot_DLNuE->GetEfficiencyErrorUp(i);
        
        if(yEff > maxEff_DLNuE){
            maxEff_DLNuE = yEff;
            maxEffBin_DLNuE = xCenter;
        }
   
        //if(i == 5) std::cout << "yErrLow = " << yErrLow << ", yErrUp = " << yErrUp << ", xErr = " << xErr << ", xCenter = " << xCenter << std::endl;
        
        gEff_DLNuE->SetPoint(i-1, xCenter, yEff);
        gEff_DLNuE->SetPointError(i-1, xErr, xErr, yErrLow, yErrUp);
    }
    
    gEff_DLNuE->SetLineColor(kBlack);
    gEff_DLNuE->SetMarkerColor(TColor::GetColor("#f89c20"));
    gEff_DLNuE->SetMarkerStyle(20);
    gEff_DLNuE->SetMarkerSize(0.7);
    gEff_DLNuE->SetLineWidth(1);
    gEff_DLNuE->Draw("P SAME");
    
    if(xmin != 999){
        /*  
        plot_BDT->GetTotalHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
        plot_DLUboone->GetTotalHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
        plot_DLNuE->GetTotalHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
        */
    }

    plot_BDT->Draw("SAME");
    plot_DLUboone->Draw("SAME");
    plot_DLNuE->Draw("SAME");
    gPad->Update();

    auto* gBDT      = plot_BDT->GetPaintedGraph();
    auto* gDLUboone = plot_DLUboone->GetPaintedGraph();
    auto* gDLNuE    = plot_DLNuE->GetPaintedGraph();

    gBDT->SetMarkerSize(0.8);
    gDLUboone->SetMarkerSize(0.8);
    gDLNuE->SetMarkerSize(0.8);

    gBDT->Draw("PE SAME");
    gDLUboone->Draw("PE SAME");
    gDLNuE->Draw("PE SAME");

    auto* g = plot_BDT->GetPaintedGraph();
    
    if(lowY == -999999 && highY == -999999){
        g->GetYaxis()->SetRangeUser(0.0, 1.05);
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
    leg->AddEntry(plot_BDT,      "BDT",       "LEP");
    leg->AddEntry(plot_DLUboone, "DL Uboone", "LEP");
    leg->AddEntry(plot_DLNuE,    "DL Nu+E",   "LEP");
    leg->Draw();

    if(effPurPlots){
        std::cout << filename << ":" << std::endl;
        std::cout << "BDT: Max Eff x Pur = " << maxEff_BDT << ", Bin Value = " << maxEffBin_BDT << std::endl;
        std::cout << "DLUboone: Max Eff x Pur = " << maxEff_DLUboone << ", Bin Value = " << maxEffBin_DLUboone << std::endl;
        std::cout << "DLNuE: Max Eff x Pur = " << maxEff_DLNuE << ", Bin Value = " << maxEffBin_DLNuE << std::endl;
        
        if (writeMaxValues){
            std::ofstream outfile(textfilename, std::ios::app);
            if(outfile.is_open()) {
                outfile << "================" << std::endl;
                outfile << filename << std::endl;

                outfile << "BDT: Max Eff x Pur = " << maxEff_BDT << ", Bin Value = " << maxEffBin_BDT << std::endl;
                outfile << "DLUboone: Max Eff x Pur = " << maxEff_DLUboone << ", Bin Value = " << maxEffBin_DLUboone << std::endl;
                outfile << "DLNuE: Max Eff x Pur = " << maxEff_DLNuE << ", Bin Value = " << maxEffBin_DLNuE << std::endl; 
                outfile << "================" << std::endl;
                outfile.close();
            } else{
                std::cerr << "Error: could not open " << textfilename << " for writing." << std::endl;
            }
        } 
    }

    c->SaveAs(filename.c_str());
    delete c;
}

void drawEfficiencyErrorsIndividual(TEfficiency* plot, const std::string& filename, double lowY, double highY, const std::string& legendLocation, const std::string& vertex, double xmin, double xmax){
    if (!plot) {
        std::cerr << "drawEfficiency: null TEfficiency pointer\n";
        return;
    }

    double maxVal = getMaxValueEfficiency(plot, false);
    double minVal = getMinValueEfficiency(plot, xmin, xmax, false);
    std::cout << "minVal = " << minVal << ", maxVal = " << maxVal << std::endl;

    TCanvas* c = new TCanvas("c_eff", "Efficiency comparison", 800, 600);
    c->SetTicks();
    c->SetLeftMargin(0.15);

    if(vertex == "BDT"){
        plot->SetMarkerColor(TColor::GetColor("#e42536"));
        plot->SetMarkerSize(0.7); 
        plot->SetLineWidth(1);
        plot->SetLineColor(kBlack);
        plot->SetMarkerStyle(20);
    } else if(vertex == "DLUboone"){
        plot->SetMarkerColor(TColor::GetColor("#5790fc"));
        plot->SetMarkerSize(0.7); 
        plot->SetLineWidth(1);
        plot->SetLineColor(kBlack);
        plot->SetMarkerStyle(20);
    } else if(vertex == "DLNuE"){
        plot->SetMarkerColor(TColor::GetColor("#f89c20"));
        plot->SetMarkerSize(0.7); 
        plot->SetLineWidth(1);
        plot->SetLineColor(kBlack);
        plot->SetMarkerStyle(20);
    }

    const TH1* hTotal = plot->GetTotalHistogram();
    int nBins = hTotal->GetNbinsX();
    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(nBins);    

    for(int i = 1; i <= nBins; ++i){
        double xCenter = hTotal->GetXaxis()->GetBinCenter(i);
        double xErr = (hTotal->GetXaxis()->GetBinUpEdge(i) - hTotal->GetXaxis()->GetBinLowEdge(i)) / 2.0;
        
        double yEff = plot->GetEfficiency(i);
        double yErrLow  = plot->GetEfficiencyErrorLow(i);
        double yErrUp   = plot->GetEfficiencyErrorUp(i);
  
        std::cout << "Bin " << i << ": yEff = " << yEff << ", yErrLow = " << yErrLow << ", yErrUp = " << yErrUp << ", xErr = " << xErr << ", xCenter = " << xCenter << std::endl; 
        gEff->SetPoint(i-1, xCenter, yEff);
        gEff->SetPointError(i-1, xErr, xErr, yErrLow, yErrUp);
    }

    gEff->SetLineColor(kBlack);
    gEff->SetMarkerColor(TColor::GetColor("#e42536"));
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
    if(vertex == "BDT"){
        leg->AddEntry(plot, "BDT", "LEP");
    } else if(vertex == "DLUboone"){
        leg->AddEntry(plot, "DL Uboone", "LEP");
    } else if(vertex == "DLNuE"){
        leg->AddEntry(plot, "DL Nu+E", "LEP");
    }

    leg->Draw();

    c->SaveAs(filename.c_str());
    delete c;
}

void drawEfficiency(TEfficiency* plot_BDT, TEfficiency* plot_DLUboone, TEfficiency* plot_DLNuE, const std::string& filename, double lowY, double highY, const std::string& legendLocation){
    if (!plot_BDT || !plot_DLUboone || !plot_DLNuE) {
        std::cerr << "drawEfficiency: null TEfficiency pointer\n";
        return;
    }

    TCanvas* c = new TCanvas("c_eff", "Efficiency comparison", 800, 600);
    c->SetTicks();

    plot_BDT->SetMarkerColor(TColor::GetColor("#e42536"));
    plot_BDT->SetMarkerSize(0.7); 
    plot_BDT->SetLineWidth(1);
    plot_BDT->SetLineColor(kBlack);
    plot_BDT->SetMarkerStyle(20);

    plot_DLUboone->SetMarkerColor(TColor::GetColor("#5790fc"));
    plot_DLUboone->SetMarkerSize(0.7); 
    plot_DLUboone->SetLineWidth(1);
    plot_DLUboone->SetLineColor(kBlack);
    plot_DLUboone->SetMarkerStyle(20);

    plot_DLNuE->SetMarkerColor(TColor::GetColor("#f89c20"));
    plot_DLNuE->SetMarkerSize(0.7); 
    plot_DLNuE->SetLineWidth(1);
    plot_DLNuE->SetLineColor(kBlack);
    plot_DLNuE->SetMarkerStyle(20);
    
    plot_BDT->Draw("AP");
    plot_DLUboone->Draw("SAME");
    plot_DLNuE->Draw("SAME");
    gPad->Update();
    
    auto* g = plot_BDT->GetPaintedGraph();
    
    if(lowY == -999999 && highY == -999999){
        g->GetYaxis()->SetRangeUser(0.0, 1.05);
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
    leg->AddEntry(plot_BDT,      "BDT",       "P");
    leg->AddEntry(plot_DLUboone, "DL Uboone", "P");
    leg->AddEntry(plot_DLNuE,    "DL Nu+E",   "P");
    leg->Draw();

    c->SaveAs(filename.c_str());
    delete c;
}

TH1F* makeTotalHist(const TH1F* h){
    TH1F* hc = (TH1F*)h->Clone(Form("%s_totalSum", h->GetName()));
    hc->Reset();

    double totalSum = h->Integral();
    for(int i = 1; i < hc->GetNbinsX() + 1; ++i){
        hc->SetBinContent(i, totalSum);
    }

    return hc;
}

void efficiency(histGroup_struct hists, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation, int* drawLine = nullptr, int* linePos = nullptr, double efficiencyWay = 0.0, const std::string& text_filename = ""){
    bool keepRight = (efficiencyWay == -1);
    // Total signal
    TH1F* hTotalSignal_BDT = (TH1F*)hists.currentSignal->Clone("hTotalSignal_BDT");
    TH1F* hTotalSignal_DLUboone = (TH1F*)hists.ubooneSignal->Clone("hTotalSignal_DLUboone");
    TH1F* hTotalSignal_DLNuE = (TH1F*)hists.nuESignal->Clone("hTotalSignal_DLNuE");

    TH1F* hTotalSummedSignal_BDT = makeTotalHist(hTotalSignal_BDT);
    TH1F* hTotalSummedSignal_DLUboone = makeTotalHist(hTotalSignal_DLUboone);
    TH1F* hTotalSummedSignal_DLNuE = makeTotalHist(hTotalSignal_DLNuE);

    // Total signal (cumulative)
    TH1F* hPassedSignal_BDT = makeCumulative(hists.currentSignal, keepRight);
    TH1F* hPassedSignal_DLUboone = makeCumulative(hists.ubooneSignal, keepRight);
    TH1F* hPassedSignal_DLNuE = makeCumulative(hists.nuESignal, keepRight);
    
    // Total background
    TH1F* hTotalBackground_BDT = (TH1F*) hists.currentCosmic->Clone("hTotalBackground_BDT");
    hTotalBackground_BDT->Reset();
    hTotalBackground_BDT->Add(hists.currentCosmic);
    hTotalBackground_BDT->Add(hists.currentSignalFuzzy);
    hTotalBackground_BDT->Add(hists.currentBNB);
    hTotalBackground_BDT->Add(hists.currentBNBFuzzy);
    TH1F* hTotalSummedBackground_BDT = makeTotalHist(hTotalBackground_BDT);

    TH1F* hTotalBackground_DLUboone = (TH1F*) hists.ubooneCosmic->Clone("hTotalBackground_DLUboone");
    hTotalBackground_DLUboone->Reset();
    hTotalBackground_DLUboone->Add(hists.ubooneCosmic);
    hTotalBackground_DLUboone->Add(hists.ubooneSignalFuzzy);
    hTotalBackground_DLUboone->Add(hists.ubooneBNB);
    hTotalBackground_DLUboone->Add(hists.ubooneBNBFuzzy);
    TH1F* hTotalSummedBackground_DLUboone = makeTotalHist(hTotalBackground_DLUboone);

    TH1F* hTotalBackground_DLNuE = (TH1F*) hists.nuECosmic->Clone("hTotalBackground_DLNuE");
    hTotalBackground_DLNuE->Reset();
    hTotalBackground_DLNuE->Add(hists.nuECosmic);
    hTotalBackground_DLNuE->Add(hists.nuESignalFuzzy);
    hTotalBackground_DLNuE->Add(hists.nuEBNB);
    hTotalBackground_DLNuE->Add(hists.nuEBNBFuzzy);
    TH1F* hTotalSummedBackground_DLNuE = makeTotalHist(hTotalBackground_DLNuE);
    
    // Total background (cumulative)
    TH1F* hPassedBackground_BDT = makeCumulative(hTotalBackground_BDT, keepRight);    
    TH1F* hPassedBackground_DLUboone = makeCumulative(hTotalBackground_DLUboone, keepRight);    
    TH1F* hPassedBackground_DLNuE = makeCumulative(hTotalBackground_DLNuE, keepRight);    
    
    // Total background (cumulative other way)
    TH1F* hRejectedBackground_BDT = makeCumulative(hTotalBackground_BDT, !keepRight);    
    TH1F* hRejectedBackground_DLUboone = makeCumulative(hTotalBackground_DLUboone, !keepRight);    
    TH1F* hRejectedBackground_DLNuE = makeCumulative(hTotalBackground_DLNuE, !keepRight);    
 
    // Total background + signal
    TH1F* hTotalEverything_BDT = (TH1F*) hists.currentCosmic->Clone("hTotalEverything_BDT");
    hTotalEverything_BDT->Reset();
    hTotalEverything_BDT->Add(hists.currentCosmic);
    hTotalEverything_BDT->Add(hists.currentSignal);
    hTotalEverything_BDT->Add(hists.currentSignalFuzzy);
    hTotalEverything_BDT->Add(hists.currentBNB);
    hTotalEverything_BDT->Add(hists.currentBNBFuzzy);
    TH1F* hTotalSummedEverything_BDT = makeTotalHist(hTotalEverything_BDT);

    TH1F* hTotalEverything_DLUboone = (TH1F*) hists.ubooneCosmic->Clone("hTotalEverything_DLUboone");
    hTotalEverything_DLUboone->Reset();
    hTotalEverything_DLUboone->Add(hists.ubooneCosmic);
    hTotalEverything_DLUboone->Add(hists.ubooneSignal);
    hTotalEverything_DLUboone->Add(hists.ubooneSignalFuzzy);
    hTotalEverything_DLUboone->Add(hists.ubooneBNB);
    hTotalEverything_DLUboone->Add(hists.ubooneBNBFuzzy);
    TH1F* hTotalSummedEverything_DLUboone = makeTotalHist(hTotalEverything_DLUboone);

    TH1F* hTotalEverything_DLNuE = (TH1F*) hists.nuECosmic->Clone("hTotalEverything_DLNuE");
    hTotalEverything_DLNuE->Reset();
    hTotalEverything_DLNuE->Add(hists.nuECosmic);
    hTotalEverything_DLNuE->Add(hists.nuESignal);
    hTotalEverything_DLNuE->Add(hists.nuESignalFuzzy);
    hTotalEverything_DLNuE->Add(hists.nuEBNB);
    hTotalEverything_DLNuE->Add(hists.nuEBNBFuzzy);
    TH1F* hTotalSummedEverything_DLNuE = makeTotalHist(hTotalEverything_DLNuE);

    // Total background + signal (cumulative)
    TH1F* hPassedEverything_BDT = makeCumulative(hTotalEverything_BDT, keepRight);
    TH1F* hPassedEverything_DLUboone = makeCumulative(hTotalEverything_DLUboone, keepRight);
    TH1F* hPassedEverything_DLNuE = makeCumulative(hTotalEverything_DLNuE, keepRight);

    /*
    std::cout << "Total num of signal events = " << hTotalSummedSignal_DLNuE->GetBinContent(3) << ", total num of background events = " << hTotalSummedBackground_DLNuE->GetBinContent(3) << ", total num of events = " << hTotalSummedEverything_DLNuE->GetBinContent(3) << std::endl;
    for(int i = 1; i < hTotalSummedSignal_BDT->GetNbinsX()+1; ++i){
        std::cout << "Bin " << i << ": passed signal = " << hPassedSignal_DLNuE->GetBinContent(i) << ", passed background = " << hPassedBackground_DLNuE->GetBinContent(i) << ", passed everything = " << hPassedEverything_DLNuE->GetBinContent(i) << std::endl;
        std::cout << "Efficiency = " << hPassedSignal_DLNuE->GetBinContent(i)/hTotalSummedSignal_DLNuE->GetBinContent(i) << ", Purity = " << hPassedSignal_DLNuE->GetBinContent(i)/(hPassedSignal_DLNuE->GetBinContent(i) + hPassedBackground_DLNuE->GetBinContent(i)) << ", background rejection = " << 1-(hPassedBackground_DLNuE->GetBinContent(i)/hTotalSummedBackground_DLNuE->GetBinContent(i)) << ", Efficiency x Purity = " << (hPassedSignal_DLNuE->GetBinContent(i)/hTotalSummedSignal_DLNuE->GetBinContent(i))*(hPassedSignal_DLNuE->GetBinContent(i)/hPassedEverything_DLNuE->GetBinContent(i)) << std::endl;
    }
    */
   

    TEfficiency* eff_BDT = new TEfficiency(*hPassedSignal_BDT, *hTotalSummedSignal_BDT);
    eff_BDT->SetTitle(Form("%s;%s;Signal Efficiency", hists.currentSignal->GetTitle(), hists.currentSignal->GetXaxis()->GetTitle()));     
    eff_BDT->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* eff_DLUboone = new TEfficiency(*hPassedSignal_DLUboone, *hTotalSummedSignal_DLUboone);
    eff_DLUboone->SetTitle(Form("%s;%s;Signal Efficiency", hists.ubooneSignal->GetTitle(), hists.ubooneSignal->GetXaxis()->GetTitle()));     
    eff_DLUboone->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* eff_DLNuE = new TEfficiency(*hPassedSignal_DLNuE, *hTotalSummedSignal_DLNuE);
    eff_DLNuE->SetTitle(Form("%s;%s;Signal Efficiency", hists.nuESignal->GetTitle(), hists.nuESignal->GetXaxis()->GetTitle()));     
    eff_DLNuE->SetStatisticOption(TEfficiency::kFNormal);

    TEfficiency* rej_BDT = new TEfficiency(*hRejectedBackground_BDT, *hTotalSummedBackground_BDT);
    rej_BDT->SetTitle(Form("%s;%s;Background Rejection", hists.currentCosmic->GetTitle(), hists.currentCosmic->GetXaxis()->GetTitle()));     
    rej_BDT->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* rej_DLUboone = new TEfficiency(*hRejectedBackground_DLUboone, *hTotalSummedBackground_DLUboone);
    rej_DLUboone->SetTitle(Form("%s;%s;Background Rejection", hists.ubooneCosmic->GetTitle(), hists.ubooneCosmic->GetXaxis()->GetTitle()));     
    rej_DLUboone->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* rej_DLNuE = new TEfficiency(*hRejectedBackground_DLNuE, *hTotalSummedBackground_DLNuE);
    rej_DLNuE->SetTitle(Form("%s;%s;Background Rejection", hists.nuECosmic->GetTitle(), hists.nuECosmic->GetXaxis()->GetTitle()));     
    rej_DLNuE->SetStatisticOption(TEfficiency::kFNormal);

    TEfficiency* pur_BDT = new TEfficiency(*hPassedSignal_BDT, *hPassedEverything_BDT); 
    pur_BDT->SetTitle(Form("%s;%s;Signal Purity", hists.currentSignal->GetTitle(), hists.currentSignal->GetXaxis()->GetTitle()));     
    pur_BDT->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* pur_DLUboone = new TEfficiency(*hPassedSignal_DLUboone, *hPassedEverything_DLUboone); 
    pur_DLUboone->SetTitle(Form("%s;%s;Signal Purity", hists.ubooneSignal->GetTitle(), hists.ubooneSignal->GetXaxis()->GetTitle()));     
    pur_DLUboone->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* pur_DLNuE = new TEfficiency(*hPassedSignal_DLNuE, *hPassedEverything_DLNuE); 
    pur_DLNuE->SetTitle(Form("%s;%s;Signal Purity", hists.nuESignal->GetTitle(), hists.nuESignal->GetXaxis()->GetTitle()));     
    pur_DLNuE->SetStatisticOption(TEfficiency::kFNormal);

    TH1F* hEffPurDenominator_BDT = (TH1F*) hTotalSummedSignal_BDT->Clone("hEffPurDenominator_BDT");
    hEffPurDenominator_BDT->Multiply(hPassedEverything_BDT);
    TH1F* hEffPurNumerator_BDT = (TH1F*) hPassedSignal_BDT->Clone("hEffPurNumerator_BDT");
    hEffPurNumerator_BDT->Multiply(hPassedSignal_BDT);

    TH1F* hEffPurDenominator_DLUboone = (TH1F*) hTotalSummedSignal_DLUboone->Clone("hEffPurDenominator_DLUboone");
    hEffPurDenominator_DLUboone->Multiply(hPassedEverything_DLUboone);
    TH1F* hEffPurNumerator_DLUboone = (TH1F*) hPassedSignal_DLUboone->Clone("hEffPurNumerator_DLUboone");
    hEffPurNumerator_DLUboone->Multiply(hPassedSignal_DLUboone);

    TH1F* hEffPurDenominator_DLNuE = (TH1F*) hTotalSummedSignal_DLNuE->Clone("hEffPurDenominator_DLNuE");
    hEffPurDenominator_DLNuE->Multiply(hPassedEverything_DLNuE);
    TH1F* hEffPurNumerator_DLNuE = (TH1F*) hPassedSignal_DLNuE->Clone("hEffPurNumerator_DLNuE");
    hEffPurNumerator_DLNuE->Multiply(hPassedSignal_DLNuE);

    TEfficiency* effPur_BDT = new TEfficiency(*hEffPurNumerator_BDT, *hEffPurDenominator_BDT);
    effPur_BDT->SetTitle(Form("%s;%s;Efficiency x Purity", hists.currentSignal->GetTitle(), hists.currentSignal->GetXaxis()->GetTitle()));     
    effPur_BDT->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* effPur_DLUboone = new TEfficiency(*hEffPurNumerator_DLUboone, *hEffPurDenominator_DLUboone);
    effPur_DLUboone->SetTitle(Form("%s;%s;Efficiency x Purity", hists.ubooneSignal->GetTitle(), hists.ubooneSignal->GetXaxis()->GetTitle()));     
    effPur_DLUboone->SetStatisticOption(TEfficiency::kFNormal);
    TEfficiency* effPur_DLNuE = new TEfficiency(*hEffPurNumerator_DLNuE, *hEffPurDenominator_DLNuE);
    effPur_DLNuE->SetTitle(Form("%s;%s;Efficiency x Purity", hists.nuESignal->GetTitle(), hists.nuESignal->GetXaxis()->GetTitle()));     
    effPur_DLNuE->SetStatisticOption(TEfficiency::kFNormal);

    std::string filenameEff = std::string(filename) + "_eff.pdf";
    std::string filenameRej = std::string(filename) + "_rej.pdf";
    std::string filenamePur = std::string(filename) + "_pur.pdf";
    double maxPurityVal = std::max({getMaxValueEfficiency(pur_BDT, true), getMaxValueEfficiency(pur_DLUboone, true), getMaxValueEfficiency(pur_DLNuE, true)});
    std::string filenameEffPur = std::string(filename) + "_effpur.pdf";
    double maxEffPurityVal = std::max({getMaxValueEfficiency(effPur_BDT, true), getMaxValueEfficiency(effPur_DLUboone, true), getMaxValueEfficiency(effPur_DLNuE, true)});
    std::string filenameEffErrors = std::string(filename) + "_errorseff.pdf";
    std::string filenameRejErrors = std::string(filename) + "_errorsrej.pdf";
    std::string filenamePurErrors = std::string(filename) + "_errorspur.pdf";
    std::string filenameEffPurErrors = std::string(filename) + "_errorseffpur.pdf";

    drawEfficiency(eff_BDT, eff_DLUboone, eff_DLNuE, filenameEff, -999999, -999999, legendLocation);
    drawEfficiency(rej_BDT, rej_DLUboone, rej_DLNuE, filenameRej, -999999, -999999, legendLocation);
    drawEfficiency(pur_BDT, pur_DLUboone, pur_DLNuE, filenamePur, 0, maxPurityVal*1.1, legendLocation);
    drawEfficiency(effPur_BDT, effPur_DLUboone, effPur_DLNuE, filenameEffPur, 0, maxEffPurityVal*1.1, legendLocation);
  
    eff_BDT->SetUseWeightedEvents(false);
    eff_DLUboone->SetUseWeightedEvents(false);
    eff_DLNuE->SetUseWeightedEvents(false);
    rej_BDT->SetUseWeightedEvents(false);
    rej_DLUboone->SetUseWeightedEvents(false);
    rej_DLNuE->SetUseWeightedEvents(false);
    pur_BDT->SetUseWeightedEvents(false);
    pur_DLUboone->SetUseWeightedEvents(false);
    pur_DLNuE->SetUseWeightedEvents(false);
    effPur_BDT->SetUseWeightedEvents(false);
    effPur_DLUboone->SetUseWeightedEvents(false);
    effPur_DLNuE->SetUseWeightedEvents(false);
    
    drawEfficiencyErrors(eff_BDT, eff_DLUboone, eff_DLNuE, filenameEffErrors, -999999, -999999, legendLocation, 0, xmin, xmax, false, text_filename);
    drawEfficiencyErrors(rej_BDT, rej_DLUboone, rej_DLNuE, filenameRejErrors, -999999, -999999, legendLocation, 0, xmin, xmax, false, text_filename);
    drawEfficiencyErrors(pur_BDT, pur_DLUboone, pur_DLNuE, filenamePurErrors, 0, maxPurityVal*1.1, legendLocation, 0, xmin, xmax, false, text_filename);
    drawEfficiencyErrors(effPur_BDT, effPur_DLUboone, effPur_DLNuE, filenameEffPurErrors, 0, maxEffPurityVal*1.1, legendLocation, 1, xmin, xmax, true, text_filename);
    
    std::string filenameEffPurBDT = std::string(filename) + "_errorsBDT_effpur.pdf";
    std::string filenameEffPurDLUboone = std::string(filename) + "_errorsDLUboone_effpur.pdf";
    std::string filenameEffPurDLNuE = std::string(filename) + "_errorsDLNuE_effpur.pdf";
    std::string filenameEffBDT = std::string(filename) + "_errorsBDT_eff.pdf";
    std::string filenameEffDLUboone = std::string(filename) + "_errorsDLUboone_eff.pdf";
    std::string filenameEffDLNuE = std::string(filename) + "_errorsDLNuE_eff.pdf";
    std::string filenameRejBDT = std::string(filename) + "_errorsBDT_rej.pdf";
    std::string filenameRejDLUboone = std::string(filename) + "_errorsDLUboone_rej.pdf";
    std::string filenameRejDLNuE = std::string(filename) + "_errorsDLNuE_rej.pdf";
    std::string filenamePurBDT = std::string(filename) + "_errorsBDT_pur.pdf";
    std::string filenamePurDLUboone = std::string(filename) + "_errorsDLUboone_pur.pdf";
    std::string filenamePurDLNuE = std::string(filename) + "_errorsDLNuE_pur.pdf";

    drawEfficiencyErrorsIndividual(effPur_BDT, filenameEffPurBDT, -999999, -999999, legendLocation, "BDT", xmin, xmax);
    drawEfficiencyErrorsIndividual(effPur_DLUboone, filenameEffPurDLUboone, -999999, -999999, legendLocation, "DLUboone", xmin, xmax);
    drawEfficiencyErrorsIndividual(effPur_DLNuE, filenameEffPurDLNuE, -999999, -999999, legendLocation, "DLNuE", xmin, xmax);
    drawEfficiencyErrorsIndividual(eff_BDT, filenameEffBDT, -999999, -999999, legendLocation, "BDT", xmin, xmax);
    drawEfficiencyErrorsIndividual(eff_DLUboone, filenameEffDLUboone, -999999, -999999, legendLocation, "DLUboone", xmin, xmax);
    drawEfficiencyErrorsIndividual(eff_DLNuE, filenameEffDLNuE, -999999, -999999, legendLocation, "DLNuE", xmin, xmax);
    drawEfficiencyErrorsIndividual(rej_BDT, filenameRejBDT, -999999, -999999, legendLocation, "BDT", xmin, xmax);
    drawEfficiencyErrorsIndividual(rej_DLUboone, filenameRejDLUboone, -999999, -999999, legendLocation, "DLUboone", xmin, xmax);
    drawEfficiencyErrorsIndividual(rej_DLNuE, filenameRejDLNuE, -999999, -999999, legendLocation, "DLNuE", xmin, xmax);
    drawEfficiencyErrorsIndividual(pur_BDT, filenamePurBDT, -999999, -999999, legendLocation, "BDT", xmin, xmax);
    drawEfficiencyErrorsIndividual(pur_DLUboone, filenamePurDLUboone, -999999, -999999, legendLocation, "DLUboone", xmin, xmax);
    drawEfficiencyErrorsIndividual(pur_DLNuE, filenamePurDLNuE, -999999, -999999, legendLocation, "DLNuE", xmin, xmax);
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
    std::string txtFileName = "purity_max_values_withCuts.txt";

    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_IntimeBNBNuE_DLUbooneNuEBDT_1Feb.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_IntimeBNBNuE_DLUbooneNuEBDT_27Jan.root");
    std::string base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT/";

    // If clearCosmicCut == 1 -> cut all PFPs with clearCosmic score of 1. If clearCosmicCut == 0 -> keep all PFPs
    int clearCosmicCut = 1;
    int numPFPs0Cut = 1;
    int numRecoNeutrinosCut = 1;
    int FVCut = 1;
    int CRUMBSCut = 0;
    int primaryPFPCut = 0;
    int trackscorePFPsCut = 0;
    int trackscoreCut = 0;
    int ETheta2Cut = 0;
    int ETheta2SumCut = 0;
   
    int highestEnergyPFPPrimary_DLNuE = 0;
    int totalLeft_DLNuE = 0;
    int highestEnergyPFPPrimaryNuE_DLNuE = 0;
    int totalLeftNuE_DLNuE = 0;

    if(clearCosmicCut == 1 && numPFPs0Cut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic.txt";
    } else if(numPFPs0Cut == 1 && numRecoNeutrinosCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0.txt";
    } else if(numRecoNeutrinosCut == 1 && FVCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut.txt";
    } else if(FVCut == 1 && CRUMBSCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fvCuts.txt";
    } else if(CRUMBSCut == 1 && primaryPFPCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_crumbs_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fv_crumbsCuts.txt";
    } else if(primaryPFPCut == 1 && trackscorePFPsCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFPCuts.txt";
    } else if(trackscorePFPsCut == 1 && trackscoreCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAllCuts.txt";
    } else if(trackscoreCut == 1 && ETheta2Cut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_trackscore_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_trackscoreCuts.txt";
    } else if(ETheta2Cut == 1 && ETheta2SumCut == 0){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_trackscore_etheta2_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_trackscore_etheta2Cuts.txt";
    } else if(ETheta2SumCut == 1){
        base_path = "/nashome/c/coackley/nuEBackgroundSignalPlotsWeightsWithCutsSPLIT_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_trackscore_etheta2_etheta2Sum_cuts/";
        txtFileName = "purity_max_values_withCuts_clearCosmic_numPFPs0_recoNeut_fv_crumbs_primaryPFP_trackscoreAll_trackscore_etheta2_etheta2SumCuts.txt";
    }
    
    
    std::ofstream clearFile(txtFileName, std::ios::trunc);
    if (!clearFile.is_open()) {
        std::cerr << "Error: could not open or create "
                  << txtFileName << std::endl;
        return;
    }
    clearFile.close();

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
    
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBUboone;
    
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNBNuE;

    double totalPOTSignalCurrent = 0;
    double totalPOTBNBCurrent = 0;
    
    double totalPOTSignalUboone = 0;
    double totalPOTBNBUboone = 0;
    
    double totalPOTSignalNuE = 0;
    double totalPOTBNBNuE = 0;

    double cosmicSpillsSumCurrent = 0;
    double cosmicSpillsSumUboone = 0;
    double cosmicSpillsSumNuE = 0;

    double BNBSpillsSumCurrent = 0;
    double BNBSpillsSumUboone = 0;
    double BNBSpillsSumNuE = 0;

    double NuESpillsSumCurrent = 0;
    double NuESpillsSumUboone = 0;
    double NuESpillsSumNuE = 0;

    double POTSignalNuE_notMissing = 0;
    double POTSignalUboone_notMissing = 0;
    double POTSignalBDT_notMissing = 0;
    
    double POTBNBNuE_notMissing = 0;
    double POTBNBUboone_notMissing = 0;
    double POTBNBBDT_notMissing = 0;
    
    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 2) cosmicSpillsSumCurrent += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 0) cosmicSpillsSumUboone += subRunNumGenEvents;
        if(subRunSignal == 3 && subRunDLCurrent == 5) cosmicSpillsSumNuE += subRunNumGenEvents;

        if(subRunSignal == 2 && subRunDLCurrent == 2) BNBSpillsSumCurrent += subRunNumGenEvents;
        if(subRunSignal == 2 && subRunDLCurrent == 0) BNBSpillsSumUboone += subRunNumGenEvents;
        if(subRunSignal == 2 && subRunDLCurrent == 5) BNBSpillsSumNuE += subRunNumGenEvents;
        
        if(subRunSignal == 1 && subRunDLCurrent == 2) NuESpillsSumCurrent += subRunNumGenEvents;
        if(subRunSignal == 1 && subRunDLCurrent == 0) NuESpillsSumUboone += subRunNumGenEvents;
        if(subRunSignal == 1 && subRunDLCurrent == 5) NuESpillsSumNuE += subRunNumGenEvents;

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
        }
    }

    //double targetPOT = POTSignalBDT_notMissing;
    double targetPOT = 1e21;
    double targetSpills = (targetPOT/(5e12));

    double BNBScaledSpills_BDT = ((targetPOT/POTBNBBDT_notMissing) * BNBSpillsSumCurrent);
    double BNBScaledSpills_Uboone = ((targetPOT/POTBNBUboone_notMissing) * BNBSpillsSumUboone);
    double BNBScaledSpills_NuE = ((targetPOT/POTBNBNuE_notMissing) * BNBSpillsSumNuE);

    double SignalScaledSpills_BDT = ((targetPOT/POTSignalBDT_notMissing) * NuESpillsSumCurrent);
    double SignalScaledSpills_Uboone = ((targetPOT/POTSignalUboone_notMissing) * NuESpillsSumUboone);
    double SignalScaledSpills_NuE = ((targetPOT/POTSignalNuE_notMissing) * NuESpillsSumNuE);

    /*
    double cosmicsWeights_BDT = ((targetSpills - BNBScaledSpills_BDT - SignalScaledSpills_BDT) / cosmicSpillsSumCurrent);
    double cosmicsWeights_Uboone = ((targetSpills - BNBScaledSpills_Uboone - SignalScaledSpills_Uboone) / cosmicSpillsSumUboone);
    double cosmicsWeights_NuE = ((targetSpills - BNBScaledSpills_NuE - SignalScaledSpills_NuE) / cosmicSpillsSumCurrent);
    double cosmicsWeights_BDT = ((targetSpills - BNBScaledSpills_BDT) / cosmicSpillsSumCurrent);
    double cosmicsWeights_Uboone = ((targetSpills - BNBScaledSpills_Uboone) / cosmicSpillsSumUboone);
    double cosmicsWeights_NuE = ((targetSpills - BNBScaledSpills_NuE) / cosmicSpillsSumCurrent);
    */

    double targetGates = ((1333568/6.293443e+18)*targetPOT);
    double cosmicsWeights_BDT = (((1-0.0754) * targetGates)/cosmicSpillsSumCurrent);
    double cosmicsWeights_Uboone = (((1-0.0754) * targetGates)/cosmicSpillsSumUboone);
    double cosmicsWeights_NuE = (((1-0.0754) * targetGates)/cosmicSpillsSumNuE);

    weights_struct weights;
    weights.signalCurrent = targetPOT / POTSignalBDT_notMissing;
    weights.BNBCurrent = targetPOT /POTBNBBDT_notMissing;
    weights.cosmicsCurrent = cosmicsWeights_BDT;
    
    weights.signalUboone = targetPOT /POTSignalUboone_notMissing;
    weights.BNBUboone = targetPOT /POTBNBUboone_notMissing;
    weights.cosmicsUboone = cosmicsWeights_Uboone;
    
    weights.signalNuE = targetPOT /POTSignalNuE_notMissing;
    weights.BNBNuE = targetPOT /POTBNBNuE_notMissing;
    weights.cosmicsNuE = cosmicsWeights_NuE;

    std::cout << "BNB POT Current = " << totalPOTBNBCurrent << ", Uboone = " << totalPOTBNBUboone << ", Nu+E = " << totalPOTBNBNuE << std::endl;
    std::cout << "ALL BNB POT Current = " << POTBNBBDT_notMissing << ", Uboone = " << POTBNBUboone_notMissing << ", Nu+E = " << POTBNBNuE_notMissing << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Signal POT Current = " << totalPOTSignalCurrent << ", Uboone = " << totalPOTSignalUboone << ", Nu+E = " << totalPOTSignalNuE << std::endl;
    std::cout << "ALL Signal POT Current = " << POTSignalBDT_notMissing << ", Uboone = " << POTSignalUboone_notMissing << ", Nu+E = " << POTSignalNuE_notMissing << std::endl;
    std::cout << "" << std::endl;
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
    auto sliceNumPrimaryPFPs = createHistGroup("sliceNumPrimaryPFPs", "Number of Primary PFPs in the Slice", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPsDist = createHistGroup("sliceNumPrimaryPFPsDist", "Number of Primary PFPs in the Slice (Not Weighted)", "Number of Primary PFPs", 20, 0, 20); 
    auto sliceNumNeutrinos = createHistGroup("sliceNumNeutrinos", "Number of Reco Neutrinos in the Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumNeutrinosDist = createHistGroup("sliceNumNeutrinosDist", "Number of Reco Neutrinos in the Slice (Not Weighted)", "Number of Reco Neutrinos", 10, 0, 10);

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

    auto trackscoreHighestEnergyPFP = createHistGroup("trackscoreHighestEnergyPFP", "Trackscore of the PFP in the Slice with the Highest Energy", "Trackscore", 20, 0, 1);
    auto trackscoreHighestEnergyPFPDist = createHistGroup("trackscoreHighestEnergyPFPDist", "Trackscore of the PFP in the Slice with the Highest Energy (Not Weighted)", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPs = createHistGroup("trackscoreAllPFPs", "Trackscore of All PFPs in the Slice", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPsDist = createHistGroup("trackscoreAllPFPsDist", "Trackscore of All PFPs in the Slice (Not Weighted)", "Trackscore", 20, 0, 1);
    auto trackscoreHighestScorePFPs = createHistGroup("trackscoreHighestScorePFPs", "Trackscore of the PFP with the Highest Trackscore in the Slice", "Trackscore", 20, 0, 1);
    auto trackscoreHighestScorePFPsDist = createHistGroup("trackscoreHighestScorePFPsDist", "Trackscore of the PFP with the Highest Trackscore in the Slice (Not Weighted)", "Trackscore", 20, 0, 1);

    auto trackscoreAllPFPsPFP = createHistGroup("trackscoreAllPFPsPFP", "Trackscore of All PFPs in the Slice (Split By PFP)", "Trackscore", 20, 0, 1);

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

    auto deltaEnergy = createHistGroup("deltaEnergy", "Energy Asymmetry of the Highest Energy PFP", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto deltaEnergyDist = createHistGroup("deltaEnergyDist", "Energy Asymmetry of the Highest Energy PFP (Not Weighted)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto deltaEnergySum = createHistGroup("deltaEnergySum", "Energy Asymmetry of the Sum of Energies of PFPs in Slice", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto deltaEnergySumDist = createHistGroup("deltaEnergySumDist", "Energy Asymmetry of the Sum of Energies of PFPs in Slice (Not Weighted)", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);

    auto pfpCompleteness = createHistGroup("pfpCompleteness", "Completeness of the Highest Energy PFP in the Slice", "Completeness", 50, 0, 1);
    auto pfpCompletenessDist = createHistGroup("pfpCompletenessDist", "Completeness of the Highest Energy PFP in the Slice (Not Weighted)", "Completeness", 50, 0, 1);
    auto pfpPurity = createHistGroup("pfpPurity", "Purity of the Highest Energy PFP in the Slice", "Purity", 50, 0, 1);
    auto pfpPurityDist = createHistGroup("pfpPurityDist", "Purity of the Highest Energy PFP in the Slice (Not Weighted)", "Purity", 50, 0, 1);

    auto recoX = createHistGroup("recoX", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 200, -202, 202);
    auto recoX_smallerBins = createHistGroup("recoX_smallerBins", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 808, -202, 202);
    auto recoXDist = createHistGroup("recoXDist", "X Coordinate of Reco Neutrino (Not Weighted)", "x_{Reco} (cm)", 200, -202, 202);
    auto recoX_low = createHistGroup("recoX_low", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, -202, -182);
    auto recoXDist_low = createHistGroup("recoXDist_low", "X Coordinate of Reco Neutrino (Not Weighted)", "x_{Reco} (cm)", 40, -202, -182);
    auto recoX_high = createHistGroup("recoX_high", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 40, 182, 202);
    auto recoXDist_high = createHistGroup("recoXDist_high", "X Coordinate of Reco Neutrino (Not Weighted)", "x_{Reco} (cm)", 40, 182, 202);
    
    auto recoY = createHistGroup("recoY", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 200, -204, 204);
    auto recoY_smallerBins = createHistGroup("recoY_smallerBins", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 816, -204, 204);
    auto recoYDist = createHistGroup("recoYDist", "Y Coordinate of Reco Neutrino (Not Weighted)", "y_{Reco} (cm)", 200, -204, 204);
    auto recoY_low = createHistGroup("recoY_low", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, -204, -184);
    auto recoYDist_low = createHistGroup("recoYDist_low", "Y Coordinate of Reco Neutrino (Not Weighted)", "y_{Reco} (cm)", 40, -204, -184);
    auto recoY_high = createHistGroup("recoY_high", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 40, 184, 204);
    auto recoYDist_high = createHistGroup("recoYDist_high", "Y Coordinate of Reco Neutrino (Not Weighted)", "y_{Reco} (cm)", 40, 184, 204);
    
    auto recoZ = createHistGroup("recoZ", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 250, 0, 510);
    auto recoZ_smallerBins = createHistGroup("recoZ_smallerBins", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 1020, 0, 510);
    auto recoZDist = createHistGroup("recoZDist", "Z Coordinate of Reco Neutrino (Not Weighted)", "z_{Reco} (cm)", 250, 0, 510);
    auto recoZ_low = createHistGroup("recoZ_low", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 40, 0, 20);
    auto recoZDist_low = createHistGroup("recoZDist_low", "Z Coordinate of Reco Neutrino (Not Weighted)", "z_{Reco} (cm)", 40, 0, 20);
    auto recoZ_high = createHistGroup("recoZ_high", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 40, 490, 510);
    auto recoZDist_high = createHistGroup("recoZDist_high", "Z Coordinate of Reco Neutrino (Not Weighted)", "z_{Reco} (cm)", 40, 490, 510);

    auto QSquaredHighest = createHistGroup("QSquaredHighest", "Q^{2} Using Highest Energy PFP in Slice", "Q^{2} (GeV^{2})", 100, 0, 0.1);
    auto QSquaredHighestDist = createHistGroup("QSquaredHighestDist", "Q^{2} Using Highest Energy PFP in Slice (Not Weighted)", "Q^{2} (GeV^{2})", 100, 0, 0.1);
    auto QSquaredSum = createHistGroup("QSquaredSum", "Q^{2} Using Sum of PFP Energies in Slice", "Q^{2} (GeV^{2})", 100, 0, 0.1);
    auto QSquaredSumDist = createHistGroup("QSquaredSumDist", "Q^{2} Using Sum of PFP Energies in Slice (Not Weighted)", "Q^{2} (GeV^{2})", 100, 0, 0.1);

    
    // Plots split up into interaction:
    // BDT Vertexing
    auto sliceCompleteness_splitBDT = createSplitHistGroup("sliceCompleteness_splitBDT", "Slice Completeness: BDT Vertexing", "Completeness", 102, 0, 1.02);
    auto slicePurity_splitBDT = createSplitHistGroup("slicePurity_splitBDT", "Slice Purity: BDT Vertexing", "Purity", 102, 0, 1.02);
    auto sliceCRUMBSScore_splitBDT = createSplitHistGroup("sliceCRUMBSScore_splitBDT", "CRUMBS Score of the Slice: BDT Vertexing", "CRUMBS Score", 25, -1, 1);
    auto sliceNumPFPs_splitBDT = createSplitHistGroup("sliceNumPFPs_splitBDT", "Number of PFPs in the Slice: BDT Vertexing", "Number of PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPs_splitBDT = createSplitHistGroup("sliceNumPrimaryPFPs_splitBDT", "Number of Primary PFPs in the Slice: BDT Vertexing", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumNeutrinos_splitBDT = createSplitHistGroup("sliceNumNeutrinos_splitBDT", "Number of Reco Neutrinos in the Slice: BDT Vertexing", "Number of Reco Neutrinos", 10, 0, 10);

    auto ERecoSumThetaReco_splitBDT = createSplitHistGroup("ERecoSumThetaReco_splitBDT", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice: BDT Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaReco_splitBDT = createSplitHistGroup("ERecoHighestThetaReco_splitBDT", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice: BDT Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);

    auto QSquaredHighest_splitBDT = createSplitHistGroup("QSquaredHighest_splitBDT", "Q^{2} Using Highest Energy PFP in Slice: BDT Vertexing", "Q^{2} (GeV^{2})", 100, 0, 0.1);
    auto QSquaredSum_splitBDT = createSplitHistGroup("QSquaredSum_splitBDT", "Q^{2} Using Sum of PFP Energies in Slice: BDT Vertexing", "Q^{2} (GeV^{2})", 100, 0, 0.1);

    auto trackscoreHighestEnergyPFP_splitBDT = createSplitHistGroup("trackscoreHighestEnergyPFP_splitBDT", "Trackscore of the PFP in the Slice with the Highest Energy: BDT Vertexing", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPs_splitBDT = createSplitHistGroup("trackscoreAllPFPs_splitBDT", "Trackscore of All PFPs in the Slice: BDT Vertexing", "Trackscore", 20, 0, 1);
    auto trackscoreHighestScorePFPs_splitBDT = createSplitHistGroup("trackscoreHighestScorePFPs_splitBDT", "Trackscore of the PFP with the Highest Trackscore in the Slice: BDT Vertexing", "Trackscore", 20, 0, 1);

    auto trackscoreAllPFPs_splitPFPBDT = createSplitPFPHistGroup("trackscoreAllPFPs_splitPFPBDT", "Trackscore of All PFPs in the Slice: BDT Vertexing", "Trackscore", 20, 0, 1);

    // DL Uboone
    auto sliceCompleteness_splitDLUboone = createSplitHistGroup("sliceCompleteness_splitDLUboone", "Slice Completeness: DL Uboone Vertexing", "Completeness", 102, 0, 1.02);
    auto slicePurity_splitDLUboone = createSplitHistGroup("slicePurity_splitDLUboone", "Slice Purity: DL Uboone Vertexing", "Purity", 102, 0, 1.02);
    auto sliceCRUMBSScore_splitDLUboone = createSplitHistGroup("sliceCRUMBSScore_splitDLUboone", "CRUMBS Score of the Slice: DL Uboone Vertexing", "CRUMBS Score", 25, -1, 1);
    auto sliceNumPFPs_splitDLUboone = createSplitHistGroup("sliceNumPFPs_splitDLUboone", "Number of PFPs in the Slice: DL Uboone Vertexing", "Number of PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPs_splitDLUboone = createSplitHistGroup("sliceNumPrimaryPFPs_splitDLUboone", "Number of Primary PFPs in the Slice: DL Uboone Vertexing", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumNeutrinos_splitDLUboone = createSplitHistGroup("sliceNumNeutrinos_splitDLUboone", "Number of Reco Neutrinos in the Slice: DL Uboone Vertexing", "Number of Reco Neutrinos", 10, 0, 10);

    auto ERecoSumThetaReco_splitDLUboone = createSplitHistGroup("ERecoSumThetaReco_splitDLUboone", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice: DL Uboone Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaReco_splitDLUboone = createSplitHistGroup("ERecoHighestThetaReco_splitDLUboone", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice: DL Uboone Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);

    auto QSquaredHighest_splitDLUboone = createSplitHistGroup("QSquaredHighest_splitDLUboone", "Q^{2} Using Highest Energy PFP in Slice: DL Uboone Vertexing", "Q^{2} (GeV^{2})", 100, 0, 0.1);
    auto QSquaredSum_splitDLUboone = createSplitHistGroup("QSquaredSum_splitDLUboone", "Q^{2} Using Sum of PFP Energies in Slice: DL Uboone Vertexing", "Q^{2} (GeV^{2})", 100, 0, 0.1);

    auto trackscoreHighestEnergyPFP_splitDLUboone = createSplitHistGroup("trackscoreHighestEnergyPFP_splitDLUboone", "Trackscore of the PFP in the Slice with the Highest Energy: DL Uboone Vertexing", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPs_splitDLUboone = createSplitHistGroup("trackscoreAllPFPs_splitDLUboone", "Trackscore of All PFPs in the Slice: DL Uboone Vertexing", "Trackscore", 20, 0, 1);
    auto trackscoreHighestScorePFPs_splitDLUboone = createSplitHistGroup("trackscoreHighestScorePFPs_splitDLUboone", "Trackscore of the PFP with the Highest Trackscore in the Slice: DL Uboone Vertexing", "Trackscore", 20, 0, 1);
    
    auto trackscoreAllPFPs_splitPFPDLUboone = createSplitPFPHistGroup("trackscoreAllPFPs_splitPFPDLUboone", "Trackscore of All PFPs in the Slice: DL Uboone Vertexing", "Trackscore", 20, 0, 1);

    // DL Nu+E
    auto sliceCompleteness_splitDLNuE = createSplitHistGroup("sliceCompleteness_splitDLNuE", "Slice Completeness: DL Nu+E Vertexing", "Completeness", 102, 0, 1.02);
    auto slicePurity_splitDLNuE = createSplitHistGroup("slicePurity_splitDLNuE", "Slice Purity: DL Nu+E Vertexing", "Purity", 102, 0, 1.02);
    auto sliceCRUMBSScore_splitDLNuE = createSplitHistGroup("sliceCRUMBSScore_splitDLNuE", "CRUMBS Score of the Slice: DL Nu+E Vertexing", "CRUMBS Score", 25, -1, 1);
    auto sliceNumPFPs_splitDLNuE = createSplitHistGroup("sliceNumPFPs_splitDLNuE", "Number of PFPs in the Slice: DL Nu+E Vertexing", "Number of PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPs_splitDLNuE = createSplitHistGroup("sliceNumPrimaryPFPs_splitDLNuE", "Number of Primary PFPs in the Slice: DL Nu+E Vertexing", "Number of Primary PFPs", 20, 0, 20);
    auto sliceNumNeutrinos_splitDLNuE = createSplitHistGroup("sliceNumNeutrinos_splitDLNuE", "Number of Reco Neutrinos in the Slice: DL Nu+E Vertexing", "Number of Reco Neutrinos", 10, 0, 10);

    auto ERecoSumThetaReco_splitDLNuE = createSplitHistGroup("ERecoSumThetaReco_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice: DL Nu+E Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaReco_splitDLNuE = createSplitHistGroup("ERecoHighestThetaReco_splitDLNuE", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice: DL Nu+E Vertexing", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);

    auto QSquaredHighest_splitDLNuE = createSplitHistGroup("QSquaredHighest_splitDLNuE", "Q^{2} Using Highest Energy PFP in Slice: DL Nu+E Vertexing", "Q^{2} (GeV^{2})", 100, 0, 0.1);
    auto QSquaredSum_splitDLNuE = createSplitHistGroup("QSquaredSum_splitDLNuE", "Q^{2} Using Sum of PFP Energies in Slice: DL Nu+E Vertexing", "Q^{2} (GeV^{2})", 100, 0, 0.1);

    auto trackscoreHighestEnergyPFP_splitDLNuE = createSplitHistGroup("trackscoreHighestEnergyPFP_splitDLNuE", "Trackscore of the PFP in the Slice with the Highest Energy: DL Nu+E Vertexing", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPs_splitDLNuE = createSplitHistGroup("trackscoreAllPFPs_splitDLNuE", "Trackscore of All PFPs in the Slice: DL Nu+E Vertexing", "Trackscore", 20, 0, 1);
    auto trackscoreHighestScorePFPs_splitDLNuE = createSplitHistGroup("trackscoreHighestScorePFPs_splitDLNuE", "Trackscore of the PFP with the Highest Trackscore in the Slice: DL Nu+E Vertexing", "Trackscore", 20, 0, 1);
    
    auto trackscoreAllPFPs_splitPFPDLNuE = createSplitPFPHistGroup("trackscoreAllPFPs_splitPFPDLNuE", "Trackscore of All PFPs in the Slice: DL Nu+E Vertexing", "Trackscore", 20, 0, 1);

    auto recoX_low_splitDLNuE = createSplitHistGroup("recoX_low_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 64, -202, -170);
    auto recoX_high_splitDLNuE = createSplitHistGroup("recoX_high_splitDLNuE", "X Coordinate of Reco Neutrino", "x_{Reco} (cm)", 64, 170, 202);
    auto recoY_low_splitDLNuE = createSplitHistGroup("recoY_low_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 68, -204, -170);
    auto recoY_high_splitDLNuE = createSplitHistGroup("recoY_high_splitDLNuE", "Y Coordinate of Reco Neutrino", "y_{Reco} (cm)", 128, 140, 204);
    auto recoZ_low_splitDLNuE = createSplitHistGroup("recoZ_low_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 100, 0, 50);
    auto recoZ_high_splitDLNuE = createSplitHistGroup("recoZ_high_splitDLNuE", "Z Coordinate of Reco Neutrino", "z_{Reco} (cm)", 100, 460, 510);

    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0; double zMax = 509.4;

    TH2D *xCoordAngleDifferenceBDT_low = new TH2D("xCoordAngleDifferenceBDT_low", "", 15, xMin, (xMin + 30), 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDT_low = new TH2D("yCoordAngleDifferenceBDT_low", "", 15, yMin, (yMin + 30), 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDT_low = new TH2D("zCoordAngleDifferenceBDT_low", "", 15, zMin, (zMin + 30), 40, 0, 180);
    TH2D *xCoordAngleDifferenceBDT_high = new TH2D("xCoordAngleDifferenceBDT_high", "", 15, (xMax - 30), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDT_high = new TH2D("yCoordAngleDifferenceBDT_high", "", 15, (yMax - 30), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDT_high = new TH2D("zCoordAngleDifferenceBDT_high", "", 20, (zMax - 40), zMax, 60, 0, 180);
    
    TH2D *xCoordAngleDifferenceDLUboone_low = new TH2D("xCoordAngleDifferenceDLUboone_low", "", 15, xMin, (xMin + 30), 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLUboone_low = new TH2D("yCoordAngleDifferenceDLUboone_low", "", 15, yMin, (yMin + 30), 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLUboone_low = new TH2D("zCoordAngleDifferenceDLUboone_low", "", 15, zMin, (zMin + 30), 40, 0, 180);
    TH2D *xCoordAngleDifferenceDLUboone_high = new TH2D("xCoordAngleDifferenceDLUboone_high", "", 15, (xMax - 30), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLUboone_high = new TH2D("yCoordAngleDifferenceDLUboone_high", "", 15, (yMax - 30), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLUboone_high = new TH2D("zCoordAngleDifferenceDLUboone_high", "", 20, (zMax - 40), zMax, 60, 0, 180);
    
    TH2D *xCoordAngleDifferenceDLNuE_low = new TH2D("xCoordAngleDifferenceDLNuE_low", "", 15, xMin, (xMin + 30), 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLNuE_low = new TH2D("yCoordAngleDifferenceDLNuE_low", "", 15, yMin, (yMin + 30), 40, 0, 180);
    TH2D *zCoordAngleDifferenceDLNuE_low = new TH2D("zCoordAngleDifferenceDLNuE_low", "", 15, zMin, (zMin + 30), 40, 0, 180);
    TH2D *xCoordAngleDifferenceDLNuE_high = new TH2D("xCoordAngleDifferenceDLNuE_high", "", 15, (xMax - 30), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceDLNuE_high = new TH2D("yCoordAngleDifferenceDLNuE_high", "", 15, (yMax - 30), yMax, 40, 0, 180);
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

    TH2D *xCoordEnergyAsymmetryHighestBDT = new TH2D("xCoordEnergyAsymmetryHighestBDT", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestBDT = new TH2D("yCoordEnergyAsymmetryHighestBDT", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestBDT = new TH2D("zCoordEnergyAsymmetryHighestBDT", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    TH2D *xCoordEnergyAsymmetryHighestBDT_low = new TH2D("xCoordEnergyAsymmetryHighestBDT_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestBDT_low = new TH2D("yCoordEnergyAsymmetryHighestBDT_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestBDT_low = new TH2D("zCoordEnergyAsymmetryHighestBDT_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    TH2D *xCoordEnergyAsymmetryHighestBDT_high = new TH2D("xCoordEnergyAsymmetryHighestBDT_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestBDT_high = new TH2D("yCoordEnergyAsymmetryHighestBDT_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestBDT_high = new TH2D("zCoordEnergyAsymmetryHighestBDT_high", "", 55, (zMax - 110), zMax, 20, -1, 1);
    
    TH2D *xCoordEnergyAsymmetryHighestDLUboone = new TH2D("xCoordEnergyAsymmetryHighestDLUboone", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestDLUboone = new TH2D("yCoordEnergyAsymmetryHighestDLUboone", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestDLUboone = new TH2D("zCoordEnergyAsymmetryHighestDLUboone", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    TH2D *xCoordEnergyAsymmetryHighestDLUboone_low = new TH2D("xCoordEnergyAsymmetryHighestDLUboone_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestDLUboone_low = new TH2D("yCoordEnergyAsymmetryHighestDLUboone_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestDLUboone_low = new TH2D("zCoordEnergyAsymmetryHighestDLUboone_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    TH2D *xCoordEnergyAsymmetryHighestDLUboone_high = new TH2D("xCoordEnergyAsymmetryHighestDLUboone_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestDLUboone_high = new TH2D("yCoordEnergyAsymmetryHighestDLUboone_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestDLUboone_high = new TH2D("zCoordEnergyAsymmetryHighestDLUboone_high", "", 55, (zMax - 110), zMax, 20, -1, 1);

    TH2D *xCoordEnergyAsymmetryHighestDLNuE = new TH2D("xCoordEnergyAsymmetryHighestDLNuE", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestDLNuE = new TH2D("yCoordEnergyAsymmetryHighestDLNuE", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestDLNuE = new TH2D("zCoordEnergyAsymmetryHighestDLNuE", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    TH2D *xCoordEnergyAsymmetryHighestDLNuE_low = new TH2D("xCoordEnergyAsymmetryHighestDLNuE_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestDLNuE_low = new TH2D("yCoordEnergyAsymmetryHighestDLNuE_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestDLNuE_low = new TH2D("zCoordEnergyAsymmetryHighestDLNuE_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    TH2D *xCoordEnergyAsymmetryHighestDLNuE_high = new TH2D("xCoordEnergyAsymmetryHighestDLNuE_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetryHighestDLNuE_high = new TH2D("yCoordEnergyAsymmetryHighestDLNuE_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetryHighestDLNuE_high = new TH2D("zCoordEnergyAsymmetryHighestDLNuE_high", "", 55, (zMax - 110), zMax, 20, -1, 1);

    TH2D *xCoordEnergyAsymmetrySummedBDT = new TH2D("xCoordEnergyAsymmetrySummedBDT", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedBDT = new TH2D("yCoordEnergyAsymmetrySummedBDT", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedBDT = new TH2D("zCoordEnergyAsymmetrySummedBDT", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    TH2D *xCoordEnergyAsymmetrySummedBDT_low = new TH2D("xCoordEnergyAsymmetrySummedBDT_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedBDT_low = new TH2D("yCoordEnergyAsymmetrySummedBDT_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedBDT_low = new TH2D("zCoordEnergyAsymmetrySummedBDT_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    TH2D *xCoordEnergyAsymmetrySummedBDT_high = new TH2D("xCoordEnergyAsymmetrySummedBDT_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedBDT_high = new TH2D("yCoordEnergyAsymmetrySummedBDT_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedBDT_high = new TH2D("zCoordEnergyAsymmetrySummedBDT_high", "", 55, (zMax - 110), zMax, 20, -1, 1);
    
    TH2D *xCoordEnergyAsymmetrySummedDLUboone = new TH2D("xCoordEnergyAsymmetrySummedDLUboone", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedDLUboone = new TH2D("yCoordEnergyAsymmetrySummedDLUboone", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedDLUboone = new TH2D("zCoordEnergyAsymmetrySummedDLUboone", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    TH2D *xCoordEnergyAsymmetrySummedDLUboone_low = new TH2D("xCoordEnergyAsymmetrySummedDLUboone_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedDLUboone_low = new TH2D("yCoordEnergyAsymmetrySummedDLUboone_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedDLUboone_low = new TH2D("zCoordEnergyAsymmetrySummedDLUboone_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    TH2D *xCoordEnergyAsymmetrySummedDLUboone_high = new TH2D("xCoordEnergyAsymmetrySummedDLUboone_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedDLUboone_high = new TH2D("yCoordEnergyAsymmetrySummedDLUboone_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedDLUboone_high = new TH2D("zCoordEnergyAsymmetrySummedDLUboone_high", "", 55, (zMax - 110), zMax, 20, -1, 1);

    TH2D *xCoordEnergyAsymmetrySummedDLNuE = new TH2D("xCoordEnergyAsymmetrySummedDLNuE", "", (int)round((xMax - xMin)/5), xMin, xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedDLNuE = new TH2D("yCoordEnergyAsymmetrySummedDLNuE", "", (int)round((yMax - yMin)/5), yMin, yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedDLNuE = new TH2D("zCoordEnergyAsymmetrySummedDLNuE", "", (int)round((zMax - zMin)/5), zMin, zMax, 20, -1, 1);
    TH2D *xCoordEnergyAsymmetrySummedDLNuE_low = new TH2D("xCoordEnergyAsymmetrySummedDLNuE_low", "", 15, xMin, (xMin + 30), 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedDLNuE_low = new TH2D("yCoordEnergyAsymmetrySummedDLNuE_low", "", 15, yMin, (yMin + 30), 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedDLNuE_low = new TH2D("zCoordEnergyAsymmetrySummedDLNuE_low", "", 15, zMin, (zMin + 30), 20, -1, 1);
    TH2D *xCoordEnergyAsymmetrySummedDLNuE_high = new TH2D("xCoordEnergyAsymmetrySummedDLNuE_high", "", 15, (xMax - 30), xMax, 20, -1, 1);
    TH2D *yCoordEnergyAsymmetrySummedDLNuE_high = new TH2D("yCoordEnergyAsymmetrySummedDLNuE_high", "", 15, (yMax - 30), yMax, 20, -1, 1);
    TH2D *zCoordEnergyAsymmetrySummedDLNuE_high = new TH2D("zCoordEnergyAsymmetrySummedDLNuE_high", "", 55, (zMax - 110), zMax, 20, -1, 1);

    double numNuESliceCategory_DLNuE = 0;
    double numNuESliceCategoryPassed_DLNuE = 0;
    double numNuESliceCategoryElse_DLNuE = 0;
    double numNuEIntType_DLNuE = 0;

    double numEvents_BDTCosmic = 0;
    double numEvents_BDTBNB = 0;
    double numEvents_BDTNuE = 0;
    
    double numEvents_DLNuECosmic = 0;
    double numEvents_DLNuEBNB = 0;
    double numEvents_DLNuENuE = 0;
                
    double numSignal_beforeCut_BDT = 0;
    double numSignalFuzzy_beforeCut_BDT = 0;
    double numBNB_beforeCut_BDT = 0;
    double numBNBFuzzy_beforeCut_BDT = 0;
    double numCosmic_beforeCut_BDT = 0;
                
    double numSignal_afterCut_BDT = 0;
    double numSignalFuzzy_afterCut_BDT = 0;
    double numBNB_afterCut_BDT = 0;
    double numBNBFuzzy_afterCut_BDT = 0;
    double numCosmic_afterCut_BDT = 0;
                
    double numSignal_beforeCut_DLUboone = 0;
    double numSignalFuzzy_beforeCut_DLUboone = 0;
    double numBNB_beforeCut_DLUboone = 0;
    double numBNBFuzzy_beforeCut_DLUboone = 0;
    double numCosmic_beforeCut_DLUboone = 0;
                
    double numSignal_afterCut_DLUboone = 0;
    double numSignalFuzzy_afterCut_DLUboone = 0;
    double numBNB_afterCut_DLUboone = 0;
    double numBNBFuzzy_afterCut_DLUboone = 0;
    double numCosmic_afterCut_DLUboone = 0;

    double numSignal_beforeCut_DLNuE = 0;
    double numSignalFuzzy_beforeCut_DLNuE = 0;
    double numBNB_beforeCut_DLNuE = 0;
    double numBNBFuzzy_beforeCut_DLNuE = 0;
    double numCosmic_beforeCut_DLNuE = 0;

    eventCounter_struct numEventCutDLNuE;
    eventCounter_struct numEventBeforeCutDLNuE;
    eventCounter_struct numEventCutWithoutWeightingDLNuE;
    eventCounter_struct numEventBeforeCutWithoutWeightingDLNuE;

    pfpCounter_struct numPFPsBeforeDLNuE;

    double numSignal_afterCut_DLNuE = 0;
    double numSignalFuzzy_afterCut_DLNuE = 0;
    double numBNB_afterCut_DLNuE = 0;
    double numBNBFuzzy_afterCut_DLNuE = 0;
    double numCosmic_afterCut_DLNuE = 0;

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
                //printf("Slice ID = %f, Category = %f, Interaction = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", reco_sliceID->at(slice), sliceCategoryPlottingMacro, reco_sliceInteraction->at(slice), reco_sliceCompleteness->at(slice), reco_slicePurity->at(slice), reco_sliceScore->at(slice));
                //if(sliceCategoryPlottingMacro != 0) printf("True Neutrino Vertex = (%f, %f, %f)\n", reco_sliceTrueVX->at(slice), reco_sliceTrueVY->at(slice), reco_sliceTrueVZ->at(slice));
        
                // Loop through all the reco neutrinos in the event          
                int PFPcounter = 0;
                int numPFPsSlice = 0;
                int numPrimaryPFPsSlice = 0;

                double summedEnergy = 0;
                double highestEnergy_PFPID = -999999;
                double highestEnergy_energy = -999999;
                double highestEnergy_theta = -999999;
                double highestEnergy_DX = -999999;
                double highestEnergy_DY = -999999;
                double highestEnergy_DZ = -999999;
                double highestEnergy_completeness = -999999;
                double highestEnergy_purity = -999999;
                double highestEnergy_trackscore = -999999;
                double highestEnergy_primary = -999999;

                double highestTrackscore = -999999;

                double sliceCategoryPlottingMacro = -999999;

                double FVCut_xLow_BDT = -190;
                double FVCut_xHigh_BDT = 190;
                double FVCut_xCentre_BDT = 10;
                double FVCut_yLow_BDT = -180;
                double FVCut_yHigh_BDT = 180;
                double FVCut_zLow_BDT = 5;
                double FVCut_zHigh_BDT = 400;

                double FVCut_xLow_DLNuE = -185;
                double FVCut_xHigh_DLNuE = 185;
                double FVCut_xCentre_DLNuE = 10;
                double FVCut_yLow_DLNuE = -185;
                double FVCut_yHigh_DLNuE = 155;
                double FVCut_zLow_DLNuE = 7;
                double FVCut_zHigh_DLNuE = 400;

                double FVCut_xLow_DLUboone = -190;
                double FVCut_xHigh_DLUboone = 180;
                double FVCut_xCentre_DLUboone = 10;
                double FVCut_yLow_DLUboone = -185;
                double FVCut_yHigh_DLUboone = 165;
                double FVCut_zLow_DLUboone = 8;
                double FVCut_zHigh_DLUboone = 400;
               
                double crumbsScoreCut_low_BDT = -0.08;
                double crumbsScoreCut_low_DLUboone = 0;
                double crumbsScoreCut_low_DLNuE = 0.24;
                
                double crumbsScoreCut_high_BDT = 0.8;
                double crumbsScoreCut_high_DLUboone = 0.8;
                double crumbsScoreCut_high_DLNuE = 0.8;
    
                double primaryPFPCut_low_BDT = 1;
                double primaryPFPCut_low_DLUboone = 1; 
                double primaryPFPCut_low_DLNuE = 1;
                
                double primaryPFPCut_high_BDT = 0;
                double primaryPFPCut_high_DLUboone = 0;
                double primaryPFPCut_high_DLNuE = 0;

                double trackscorePFPs_upper_BDT = 0.5;
                double trackscorePFPs_lower_BDT = 0.2;

                double EThetaCut_highestPFP_BDT = 2.5;
                double EThetaCut_highestPFP_DLUboone = 3;
                double EThetaCut_highestPFP_DLNuE = 2.5;
                
                double EThetaCut_summedPFP_BDT = 3.577;
                double EThetaCut_summedPFP_DLUboone = 4.6;
                double EThetaCut_summedPFP_DLNuE = 2.5;
                
                double trackscore_highestPFP_high_BDT = 0.4;
                double trackscore_highestPFP_high_DLUboone = 0.4;
                double trackscore_highestPFP_high_DLNuE = 0.4;
                
                double trackscore_highestPFP_low_BDT = 0.25;
                double trackscore_highestPFP_low_DLUboone = 0.25;
                double trackscore_highestPFP_low_DLNuE = 0.25;
                
               
                for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                    if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                        if(clearCosmicCut == 1){
                            if(reco_particleClearCosmic->at(pfp) == 0){
                                if(reco_particleIsPrimary->at(pfp) == 1){
                                    numPrimaryPFPsSlice++;
                                }
                            }
                        } else if(clearCosmicCut == 0){
                            if(reco_particleIsPrimary->at(pfp) == 1){
                                numPrimaryPFPsSlice++;
                            }
                        
                        }
                    }
                } 

                for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                    PFPcounter++;
                    if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                        // This PFP is in the slice
                        if(clearCosmicCut == 1){
                            // Applying a cut that PFPs must not be a clear cosmic 
                            if(reco_particleClearCosmic->at(pfp) == 0){
                                // This is a PFP that isn't a clear cosmic
                                if(trackscorePFPsCut == 1){
                                    if(reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT){
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
                                            highestEnergy_trackscore = reco_particleTrackScore->at(pfp);
                                            highestEnergy_primary = reco_particleIsPrimary->at(pfp);
                                        }

                                        if(reco_particleTrackScore->at(pfp) > highestTrackscore) highestTrackscore = reco_particleTrackScore->at(pfp);
                                    }
                                } else if(trackscorePFPsCut == 0){
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
                                        highestEnergy_trackscore = reco_particleTrackScore->at(pfp);
                                        highestEnergy_primary = reco_particleIsPrimary->at(pfp);
                                    }

                                    if(reco_particleTrackScore->at(pfp) > highestTrackscore) highestTrackscore = reco_particleTrackScore->at(pfp);
                                }
                            }
                        } else if(clearCosmicCut == 0){
                            // Don't apply a cut that PFPs must not be a clear cosmic
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
                                highestEnergy_trackscore = reco_particleTrackScore->at(pfp);
                                highestEnergy_primary = reco_particleIsPrimary->at(pfp);
                            }
                                
                            if(reco_particleTrackScore->at(pfp) > highestTrackscore) highestTrackscore = reco_particleTrackScore->at(pfp);
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

                int numRecoNeutrinos = 0;

                for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                    if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                        // Reco neutrino is in the slice
                        numRecoNeutrinos++;
                        recoVX = reco_neutrinoVX->at(recoNeut);
                        recoVY = reco_neutrinoVY->at(recoNeut);
                        recoVZ = reco_neutrinoVZ->at(recoNeut);
                        //printf("Reco Neutrino in Slice: ID = %f, PDG = %f, Vertex = (%f, %f, %f)\n", reco_neutrinoID->at(recoNeut), reco_neutrinoPDG->at(recoNeut), reco_neutrinoVX->at(recoNeut), reco_neutrinoVY->at(recoNeut), reco_neutrinoVZ->at(recoNeut));
                    }
                }

                double Q2HighestValue = -999999;
                double Q2SumValue = -999999;

                if(highestEnergy_PFPID != -999999){
                    double cosRecoilAngle = TMath::Cos(highestEnergy_theta);
                    double momentumHighest = std::sqrt((highestEnergy_energy * highestEnergy_energy) - (0.511 * 0.511));
                    double momentumSum = std::sqrt((summedEnergy * summedEnergy) - (0.511 * 0.511));

                    //std::cout << "cosRecoilAngle = " << cosRecoilAngle << ", highestEnergy_theta = " << highestEnergy_theta << std::endl;
                    double ENuHighest = (((939.565 * highestEnergy_energy) - ((0.511 * 0.511)/2))/(939.565 - highestEnergy_energy + (momentumHighest * cosRecoilAngle)));
                    Q2HighestValue = (1e-6) * (2 * 939.565 * (ENuHighest - highestEnergy_energy));
                    double ENuSum = (((939.565 * summedEnergy) - ((0.511 * 0.511)/2))/(939.565 - summedEnergy + (momentumSum * cosRecoilAngle)));
                    Q2SumValue = (1e-6) * (2 * 939.565 * (ENuSum - summedEnergy));
                }

                // Assigning new category to the slices
                // 0 = cosmic, 1 = signal, 2 = signal fuzzy, 3 = bnb, 4 = bnb fuzzy
                if(reco_sliceOrigin->at(slice) == 0){
                    // This is a cosmic slice
                    sliceCategoryPlottingMacro = 0;
                } else if(reco_sliceOrigin->at(slice) == 1){
                    // This is a nu+e elastic scatter slice
                    if(DLCurrent == 5) numNuESliceCategory_DLNuE++;
                    if(reco_sliceCompleteness->at(slice) > 0.5){
                        sliceCategoryPlottingMacro = 1;
                        if(DLCurrent == 5) numNuESliceCategoryPassed_DLNuE++;
                    } else{
                        sliceCategoryPlottingMacro = 2;
                        if(DLCurrent == 5) numNuESliceCategoryElse_DLNuE++;
                    }
                } else if(reco_sliceOrigin->at(slice) == 3){
                    // This is a BNB slice
                    if(reco_sliceCompleteness->at(slice) > 0.5){
                        sliceCategoryPlottingMacro = 3;
                    } else{
                        sliceCategoryPlottingMacro = 4;
                    }
                }


                // Assigning event type for split histograms
                int sliceEventType = -999999;
                // Event types: Cosmic = 0, nu+e scatter = 1, NC Npi0 = 2, other NC = 3, CC numu = 4, CC nue = 5, Dirt = 6, Dirt nu+e = 7
                if(reco_sliceOrigin->at(slice) != 0){
                    // This is a slice that isn't truth-matched to a cosmic
                    if(reco_sliceOrigin->at(slice) == 1){
                        // This is a slice that is truth-matched to a nu+e elastic scatter
                        if(DLCurrent == 5) numNuEIntType_DLNuE++;
                        if(reco_sliceTrueVX->at(slice) > -201.3 && reco_sliceTrueVX->at(slice) < 201.3 &&
                           reco_sliceTrueVY->at(slice) > -203.8 && reco_sliceTrueVY->at(slice) < 203.8 &&
                           reco_sliceTrueVZ->at(slice) > 0      && reco_sliceTrueVZ->at(slice) < 509.5){
                            // Interaction happened inside the TPC
                            sliceEventType = 1;
                        } else{
                            // Interaction happened outside the TPC - Dirt nu+e event
                            sliceEventType = 7;
                        }
                    } else if(reco_sliceOrigin->at(slice) == 3){
                        // This is a slice that is truth-matched to a beam neutrino that isn't a nu+e elastic scatter
                        if(reco_sliceTrueVX->at(slice) > -201.3 && reco_sliceTrueVX->at(slice) < 201.3 &&
                           reco_sliceTrueVY->at(slice) > -203.8 && reco_sliceTrueVY->at(slice) < 203.8 &&
                           reco_sliceTrueVZ->at(slice) > 0      && reco_sliceTrueVZ->at(slice) < 509.5){
                            // Interaction happened inside the TPC
                            
                            if(reco_sliceTrueCCNC->at(slice) == 0){
                                // This is a CC process
                                if(reco_sliceTrueNeutrinoType->at(slice) == 12){
                                    // This is a CC nue
                                    sliceEventType = 5;
                                } else if(reco_sliceTrueNeutrinoType->at(slice) == 14){
                                    // This is a CC numu
                                    sliceEventType = 4;
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
                                    sliceEventType = 2;
                                } else{
                                    // This is an NC other process
                                    sliceEventType = 3;
                                }
                            }
                        } else{
                            // Interaction happened outside the TPC - Dirt event
                            sliceEventType = 6;       
                        }
                    }
                } else{
                    // This is a cosmic event
                    sliceEventType = 0;
                    
                }
                
                if(sliceEventType == 0){ //std::cout << "Event type = Cosmic" << std::endl;
                } else if(sliceEventType == 1){ //std::cout << "Event type = nu+e elastic scatter" << std::endl;
                } else if(sliceEventType == 2){ //std::cout << "Event type = NC Npi0" << std::endl;
                } else if(sliceEventType == 3){ //std::cout << "Event type = other NC" << std::endl;
                } else if(sliceEventType == 4){ //std::cout << "Event type = CC numu" << std::endl;
                } else if(sliceEventType == 5){ //std::cout << "Event type = CC nue" << std::endl;
                } else if(sliceEventType == 6){ //std::cout << "Event type = Dirt" << std::endl;
                } else if(sliceEventType == 7){ //std::cout << "Event type = Dirt nu+e" << std::endl;
                } else{
                    //std::cout << "No event type assigned" << std::endl;                
                    sliceEventType = 8;
                }

                // Applying cuts here
                if(DLCurrent == 2){
                    if(sliceCategoryPlottingMacro == 0) numCosmic_beforeCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) numSignal_beforeCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 2 && signal == 1) numSignalFuzzy_beforeCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 3) numBNB_beforeCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 4) numBNBFuzzy_beforeCut_BDT += weight;
             
                    if(numPFPs0Cut == 1 && numPFPsSlice == 0){
                        // This is a slice with 0 PFPs in it
                        continue;
                    }

                    if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0){
                        // This is a slice with no reco neutrino
                        continue;
                    }


                    if(FVCut == 1){
                        if (!(recoVX < FVCut_xHigh_BDT && recoVX > FVCut_xLow_BDT  && std::abs(recoVX) > FVCut_xCentre_BDT && recoVY < FVCut_yHigh_BDT && recoVY > FVCut_yLow_BDT && recoVZ > FVCut_zLow_BDT && recoVZ < FVCut_zHigh_BDT)){ 
                            std::cout << "BDT: DOES NOT PASS CUTS WITH VX = " << recoVX << ", VY = " << recoVY << ", VZ = " << recoVZ << std::endl;
                            if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a BDT signal event" << std::endl;
                            continue;
                        }
                    }

                    if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low_BDT || reco_sliceScore->at(slice) > crumbsScoreCut_high_BDT)){
                        std::cout << "BDT: DOES NOT PASS CUTS WITH CRUMBS SCORE = " << reco_sliceScore->at(slice) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a BDT signal event" << std::endl;
                        continue;
                    }

                    if(primaryPFPCut == 1 && (numPrimaryPFPsSlice > primaryPFPCut_low_BDT || numPrimaryPFPsSlice < primaryPFPCut_high_BDT)){
                        std::cout << "BDT: DOES NOT PASS CUTS WITH NUMBER OF PRIMARY PFPS = " << numPrimaryPFPsSlice << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a BDT signal event" << std::endl;
                        continue;
                    }
                   
                    if(ETheta2Cut == 1 && (highestEnergy_energy * highestEnergy_theta * highestEnergy_theta) > EThetaCut_highestPFP_BDT){
                        std::cout << "BDT: DOES NOT PASS CUTS WITH ETHETA = " << (highestEnergy_energy * highestEnergy_theta * highestEnergy_theta) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a BDT signal event" << std::endl;
                        continue;
                    }
                    
                    if(ETheta2SumCut == 1 && (summedEnergy * highestEnergy_theta * highestEnergy_theta) > EThetaCut_summedPFP_BDT){
                        std::cout << "BDT: DOES NOT PASS CUTS WITH ETHETA (SUMMED) = " << (summedEnergy * highestEnergy_theta * highestEnergy_theta) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a BDT signal event" << std::endl;
                        continue;
                    }
                    
                    if(trackscoreCut == 1 && (highestEnergy_trackscore > trackscore_highestPFP_high_BDT || highestEnergy_trackscore < trackscore_highestPFP_low_BDT)){
                        std::cout << "BDT: DOES NOT PASS CUTS WITH TRACKSCORE = " << highestEnergy_trackscore << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a BDT signal event" << std::endl;
                        continue;
                    }
                    
                    if(sliceCategoryPlottingMacro == 0) numCosmic_afterCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) numSignal_afterCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 2 && signal == 1) numSignalFuzzy_afterCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 3) numBNB_afterCut_BDT += weight;
                    if(sliceCategoryPlottingMacro == 4) numBNBFuzzy_afterCut_BDT += weight;
                
                } else if(DLCurrent == 5){
                    if(sliceCategoryPlottingMacro == 0) numCosmic_beforeCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) numSignal_beforeCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 2 && signal == 1) numSignalFuzzy_beforeCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 3) numBNB_beforeCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 4) numBNBFuzzy_beforeCut_DLNuE += weight;
                   
                    if(sliceEventType == 0){
                        numEventBeforeCutDLNuE.cosmic += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.cosmic++;
                    } else if(sliceEventType == 1 && signal == 1){
                        numEventBeforeCutDLNuE.nuE += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.nuE++;
                    } else if(sliceEventType == 2){
                        numEventBeforeCutDLNuE.NCNPi0 += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.NCNPi0++;
                    } else if(sliceEventType == 3){
                        numEventBeforeCutDLNuE.otherNC += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.otherNC++;
                    } else if(sliceEventType == 4){
                        numEventBeforeCutDLNuE.CCnumu += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.CCnumu++;
                    } else if(sliceEventType == 5){
                        numEventBeforeCutDLNuE.CCnue += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.CCnue++;
                    } else if(sliceEventType == 6){
                        numEventBeforeCutDLNuE.dirt += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.dirt++;
                    } else if(sliceEventType == 7 && signal == 1){
                        numEventBeforeCutDLNuE.nuEDirt += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.nuEDirt++;
                    } else if(sliceEventType == 8){
                        numEventBeforeCutDLNuE.other += weight; 
                        numEventBeforeCutWithoutWeightingDLNuE.other++;
                    }

                    if(numPFPs0Cut == 1 && numPFPsSlice == 0){
                        // This is a slice with 0 PFPs in it
                        continue;
                    }
                    
                    if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0){
                        // This is a slice with no reco neutrino
                        continue;
                    }
               
                    if(FVCut == 1){ 
                        if (!(recoVX < FVCut_xHigh_DLNuE && recoVX > FVCut_xLow_DLNuE && std::abs(recoVX) > FVCut_xCentre_DLNuE && recoVY < FVCut_yHigh_DLNuE && recoVY > FVCut_yLow_DLNuE && recoVZ > FVCut_zLow_DLNuE && recoVZ < FVCut_zHigh_DLNuE)){ 
                            std::cout << "DLNuE: DOES NOT PASS CUTS WITH VX = " << recoVX << ", VY = " << recoVY << ", VZ = " << recoVZ << std::endl;
                            if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLNuE signal event" << std::endl;
                            continue;
                        }
                    }
                    
                    if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low_DLNuE || reco_sliceScore->at(slice) > crumbsScoreCut_high_DLNuE)){
                        std::cout << "DLNuE: DOES NOT PASS CUTS WITH CRUMBS SCORE = " << reco_sliceScore->at(slice) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLNuE signal event" << std::endl;
                        continue;
                    }
                    
                    if(primaryPFPCut == 1 && (numPrimaryPFPsSlice > primaryPFPCut_low_DLNuE || numPrimaryPFPsSlice < primaryPFPCut_high_DLNuE)){
                        std::cout << "DLNuE: DOES NOT PASS CUTS WITH NUMBER OF PRIMARY PFPS = " << numPrimaryPFPsSlice << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLNuE signal event" << std::endl;
                        continue;
                    }
                    
                    if(ETheta2Cut == 1 && (highestEnergy_energy * highestEnergy_theta * highestEnergy_theta) > EThetaCut_highestPFP_DLNuE){
                        std::cout << "DLNuE: DOES NOT PASS CUTS WITH ETHETA = " << (highestEnergy_energy * highestEnergy_theta * highestEnergy_theta) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLNuE signal event" << std::endl;
                        continue;
                    }
                    
                    if(ETheta2SumCut == 1 && (summedEnergy * highestEnergy_theta * highestEnergy_theta) > EThetaCut_summedPFP_DLNuE){
                        std::cout << "DLNuE: DOES NOT PASS CUTS WITH ETHETA (SUMMED) = " << (summedEnergy * highestEnergy_theta * highestEnergy_theta) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLNuE signal event" << std::endl;
                        continue;
                    }
                    
                    if(trackscoreCut == 1 && (highestEnergy_trackscore > trackscore_highestPFP_high_DLNuE || highestEnergy_trackscore < trackscore_highestPFP_low_DLNuE)){
                        std::cout << "DLNuE: DOES NOT PASS CUTS WITH TRACKSCORE = " << highestEnergy_trackscore << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLNuE signal event" << std::endl;
                        continue;
                    }

                    totalLeft_DLNuE++;
                    if(highestEnergy_primary == 1){
                        highestEnergyPFPPrimary_DLNuE++;
                    }
                    
                    if(sliceCategoryPlottingMacro == 0) numCosmic_afterCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 1 && signal == 1){
                        numSignal_afterCut_DLNuE += weight;
                        totalLeftNuE_DLNuE++;
                        if(highestEnergy_primary == 1){
                            highestEnergyPFPPrimaryNuE_DLNuE++;
                        }
                    }
                    if(sliceCategoryPlottingMacro == 2 && signal == 1) numSignalFuzzy_afterCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 3) numBNB_afterCut_DLNuE += weight;
                    if(sliceCategoryPlottingMacro == 4) numBNBFuzzy_afterCut_DLNuE += weight;
                
                } else if(DLCurrent == 0){
                    if(sliceCategoryPlottingMacro == 0) numCosmic_beforeCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) numSignal_beforeCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 2 && signal == 1) numSignalFuzzy_beforeCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 3) numBNB_beforeCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 4) numBNBFuzzy_beforeCut_DLUboone += weight;
                    
                    if(numPFPs0Cut == 1 && numPFPsSlice == 0){
                        // This is a slice with 0 PFPs in it
                        continue;
                    }
                    
                    if(numRecoNeutrinosCut == 1 && numRecoNeutrinos == 0){
                        // This is a slice with no reco neutrino
                        continue;
                    }
                
                    if(FVCut == 1){
                        if (!(recoVX < FVCut_xHigh_DLUboone && recoVX > FVCut_xLow_DLUboone && std::abs(recoVX) > FVCut_xCentre_DLUboone && recoVY < FVCut_yHigh_DLUboone && recoVY > FVCut_yLow_DLUboone && recoVZ > FVCut_zLow_DLUboone && recoVZ < FVCut_zHigh_DLUboone)){ 
                            std::cout << "DLUboone: DOES NOT PASS CUTS WITH VX = " << recoVX << ", VY = " << recoVY << ", VZ = " << recoVZ << std::endl;
                            if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLUboone signal event" << std::endl;
                            continue;
                        }
                    }
                    
                    if(CRUMBSCut == 1 && (reco_sliceScore->at(slice) < crumbsScoreCut_low_DLUboone || reco_sliceScore->at(slice) > crumbsScoreCut_high_DLUboone)){
                        std::cout << "DLUboone: DOES NOT PASS CUTS WITH CRUMBS SCORE = " << reco_sliceScore->at(slice) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLUboone signal event" << std::endl;
                        continue;
                    }
                    
                    if(primaryPFPCut == 1 && (numPrimaryPFPsSlice > primaryPFPCut_low_DLUboone || numPrimaryPFPsSlice < primaryPFPCut_high_DLUboone)){
                        std::cout << "DLUboone: DOES NOT PASS CUTS WITH NUMBER OF PRIMARY PFPS = " << numPrimaryPFPsSlice << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLUboone signal event" << std::endl;
                        continue;
                    }
                    
                    if(ETheta2Cut == 1 && (highestEnergy_energy * highestEnergy_theta * highestEnergy_theta) > EThetaCut_highestPFP_DLUboone){
                        std::cout << "DLUboone: DOES NOT PASS CUTS WITH ETHETA = " << (highestEnergy_energy * highestEnergy_theta * highestEnergy_theta) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLUboone signal event" << std::endl;
                        continue;
                    }

                    if(ETheta2SumCut == 1 && (summedEnergy * highestEnergy_theta * highestEnergy_theta) > EThetaCut_summedPFP_DLUboone){
                        std::cout << "DLUboone: DOES NOT PASS CUTS WITH ETHETA (SUMMED) = " << (summedEnergy * highestEnergy_theta * highestEnergy_theta) << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLUboone signal event" << std::endl;
                        continue;
                    }
                    
                    if(trackscoreCut == 1 && (highestEnergy_trackscore > trackscore_highestPFP_high_DLUboone || highestEnergy_trackscore < trackscore_highestPFP_low_DLUboone)){
                        std::cout << "DLUboone: DOES NOT PASS CUTS WITH TRACKSCORE = " << highestEnergy_trackscore << std::endl;
                        if(sliceCategoryPlottingMacro == 0) std::cout << "Cutting out a DLUboone signal event" << std::endl;
                        continue;
                    }
                    
                    if(sliceCategoryPlottingMacro == 0) numCosmic_afterCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 1 && signal == 1) numSignal_afterCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 2 && signal == 1) numSignalFuzzy_afterCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 3) numBNB_afterCut_DLUboone += weight;
                    if(sliceCategoryPlottingMacro == 4) numBNBFuzzy_afterCut_DLUboone += weight;
                
                }

                // Looping through PFPs
                // Assigning category to PFPs
                for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                    if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                        double pfpCategory = -999999;
                        if(reco_particleTrueOrigin->at(pfp) == 1){
                            // Beam Origin
                            if(reco_particleTrueInteractionType->at(pfp) == 1098){
                                // Comes from a nu+e elastic scatter
                                if(reco_particleTruePDG->at(pfp) == 11){
                                    pfpCategory = 1;
                                    numPFPsBeforeDLNuE.nuEElectron++;
                                } else{
                                    // Something from the nu+e that isn't an electron
                                    pfpCategory = 10;
                                    numPFPsBeforeDLNuE.nuEOther++;
                                }
                            } else{
                                // From the beam but not a nu+e elasic scatter
                                if(std::abs(reco_particleTruePDG->at(pfp)) == 11){
                                    pfpCategory = 2;
                                    numPFPsBeforeDLNuE.electron++;
                                } else if(std::abs(reco_particleTruePDG->at(pfp)) == 2212){
                                    pfpCategory = 3;
                                    numPFPsBeforeDLNuE.proton++;
                                } else if(std::abs(reco_particleTruePDG->at(pfp)) == 13){
                                    pfpCategory = 4;
                                    numPFPsBeforeDLNuE.muon++;
                                } else if(std::abs(reco_particleTruePDG->at(pfp)) == 111){
                                    pfpCategory = 7;
                                    numPFPsBeforeDLNuE.pi0++;
                                } else if(std::abs(reco_particleTruePDG->at(pfp)) == 211){
                                    pfpCategory = 8;
                                    numPFPsBeforeDLNuE.chargedPi++;
                                } else{
                                    pfpCategory = 9;
                                    numPFPsBeforeDLNuE.other++;
                                }
                            }
                        } else if(reco_particleTrueOrigin->at(pfp) == 2){
                            // Cosmic Origin
                            if(std::abs(reco_particleTruePDG->at(pfp)) == 13){
                                pfpCategory = 5;
                                numPFPsBeforeDLNuE.cosmicMuon++;
                            } else{
                                pfpCategory = 6;
                                numPFPsBeforeDLNuE.cosmicOther++;
                            }
                        } else{
                            pfpCategory = 9;
                            numPFPsBeforeDLNuE.other++;
                        }

                        // Categories: nu+e electron = 1, electron = 2, proton = 3, muon = 4, cosmic muon = 5, cosmic other = 6, pi0 = 7
                        // charged pion = 8, other = 9, other from nu+e = 10
                    
                        // PLOTS HERE
                        if(pfpCategory == 1){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.nuEElectron->Fill(reco_particleTrackScore->at(pfp), weight);
                                    if(sliceCategoryPlottingMacro == 1 && signal == 1){
                                        // slice is a nu+e elastic scatter with completeness > 0.5 and PFP is truth-matched to a recoil electron
                                        trackscoreAllPFPsPFP.currentSignal->Fill(reco_particleTrackScore->at(pfp), weight);
                                    } else if(sliceCategoryPlottingMacro == 2 && signal == 1){
                                        trackscoreAllPFPsPFP.currentSignalFuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                    }
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.nuEElectron->Fill(reco_particleTrackScore->at(pfp), weight);
                                    if(sliceCategoryPlottingMacro == 1 && signal == 1){
                                        // slice is a nu+e elastic scatter with completeness > 0.5 and PFP is truth-matched to a recoil electron
                                        trackscoreAllPFPsPFP.ubooneSignal->Fill(reco_particleTrackScore->at(pfp), weight);
                                    } else if(sliceCategoryPlottingMacro == 2 && signal == 1){
                                        trackscoreAllPFPsPFP.ubooneSignalFuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                    }
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.nuEElectron->Fill(reco_particleTrackScore->at(pfp), weight);
                                    if(sliceCategoryPlottingMacro == 1 && signal == 1){
                                        // slice is a nu+e elastic scatter with completeness > 0.5 and PFP is truth-matched to a recoil electron
                                        trackscoreAllPFPsPFP.nuESignal->Fill(reco_particleTrackScore->at(pfp), weight);
                                    } else if(sliceCategoryPlottingMacro == 2 && signal == 1){
                                        trackscoreAllPFPsPFP.nuESignalFuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                    }
                                }
                            }
                        } else if(pfpCategory == 2){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.electron->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.electron->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.electron->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 3){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.proton->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.proton->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.proton->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 4){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.muon->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.muon->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.muon->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 5){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.cosmicMuon->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentCosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.cosmicMuon->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneCosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.cosmicMuon->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuECosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 6){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.cosmicOther->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentCosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.cosmicOther->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneCosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.cosmicOther->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuECosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 7){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.pi0->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.pi0->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.pi0->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 8){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.chargedPi->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.chargedPi->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.chargedPi->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 9){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.other->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.other->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.other->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        } else if(pfpCategory == 10){
                            if(DLCurrent == 2){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPBDT.nuEOther->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.currentBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 0){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLUboone.nuEOther->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.ubooneBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            } else if(DLCurrent == 5){
                                if(reco_particleTrackScore->at(pfp) != -999999 && ((trackscorePFPsCut == 1 && reco_particleTrackScore->at(pfp) < trackscorePFPs_upper_BDT && reco_particleTrackScore->at(pfp) > trackscorePFPs_lower_BDT) || trackscorePFPsCut == 0)){
                                    trackscoreAllPFPs_splitPFPDLNuE.nuEOther->Fill(reco_particleTrackScore->at(pfp), weight);
                                    trackscoreAllPFPsPFP.nuEBNB->Fill(reco_particleTrackScore->at(pfp), weight);
                                }
                            }
                        }
                    }
                } 

                // Filling Split Histograms
                if(sliceEventType == 0){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.cosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.cosmic->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.cosmic->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.cosmic->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.cosmic->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.cosmic->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            trackscoreHighestEnergyPFP_splitBDT.cosmic->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.cosmic->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.cosmic->Fill(highestTrackscore, weight);

                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.cosmic->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.cosmic->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.cosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.cosmic->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.cosmic->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.cosmic->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.cosmic->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.cosmic->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.cosmic->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.cosmic->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }
                            
                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.cosmic->Fill(highestTrackscore, weight);
    
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.cosmic->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.cosmic->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.cosmic += weight;
                        numEventCutWithoutWeightingDLNuE.cosmic++;
                        sliceCompleteness_splitDLNuE.cosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.cosmic->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.cosmic->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.cosmic->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.cosmic->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.cosmic->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLNuE.cosmic->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.cosmic->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.cosmic->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.cosmic->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.cosmic->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.cosmic->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.cosmic->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.cosmic->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.cosmic->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.cosmic->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.cosmic->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 1 && signal == 1){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.nu_e->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.nu_e->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.nu_e->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.nu_e->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.nu_e->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.nu_e->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.nu_e->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.nu_e->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.nu_e->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.nu_e->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.nu_e->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.nu_e->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.nu_e->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.nu_e->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.nu_e->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.nu_e->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.nu_e->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.nu_e->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.nu_e->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.nu_e->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.nu_e->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.nu_e->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.nu_e->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.nu_e->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.nu_e->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.nu_e->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.nuE += weight;
                        numEventCutWithoutWeightingDLNuE.nuE++;
                        sliceCompleteness_splitDLNuE.nu_e->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.nu_e->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.nu_e->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.nu_e->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.nu_e->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.nu_e->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.nu_e->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.nu_e->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                        
                            trackscoreHighestEnergyPFP_splitDLNuE.nu_e->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.nu_e->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.nu_e->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.nu_e->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.nu_e->Fill(Q2SumValue, weight);
                        }
                        
                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.nu_e->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.nu_e->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.nu_e->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.nu_e->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.nu_e->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.nu_e->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 2){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.NCNpi0->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.NCNpi0->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.NCNpi0->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.NCNpi0->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.NCNpi0->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.NCNpi0->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.NCNpi0->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.NCNpi0->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.NCNpi0->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.NCNpi0->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.NCNpi0->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.NCNpi0->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.NCNpi0->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.NCNpi0->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.NCNpi0->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.NCNpi0->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.NCNpi0->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.NCNpi0->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.NCNpi0->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.NCNpi0->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.NCNpi0->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.NCNpi0->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.NCNpi0->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.NCNpi0->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.NCNpi0->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.NCNpi0->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.NCNPi0 += weight;
                        numEventCutWithoutWeightingDLNuE.NCNPi0++;
                        sliceCompleteness_splitDLNuE.NCNpi0->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.NCNpi0->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.NCNpi0->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.NCNpi0->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.NCNpi0->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.NCNpi0->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.NCNpi0->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.NCNpi0->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            
                            trackscoreHighestEnergyPFP_splitDLNuE.NCNpi0->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.NCNpi0->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.NCNpi0->Fill(highestTrackscore, weight);

                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.NCNpi0->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.NCNpi0->Fill(Q2SumValue, weight);
                        }
                        
                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.NCNpi0->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.NCNpi0->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.NCNpi0->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.NCNpi0->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.NCNpi0->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.NCNpi0->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 3){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.otherNC->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.otherNC->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.otherNC->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.otherNC->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.otherNC->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.otherNC->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.otherNC->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.otherNC->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.otherNC->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.otherNC->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.otherNC->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.otherNC->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.otherNC->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.otherNC->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.otherNC->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.otherNC->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.otherNC->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.otherNC->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.otherNC->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.otherNC->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.otherNC->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.otherNC->Fill(highestEnergy_trackscore, weight);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.otherNC->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.otherNC->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.otherNC->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.otherNC->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.otherNC += weight;
                        numEventCutWithoutWeightingDLNuE.otherNC++;
                        sliceCompleteness_splitDLNuE.otherNC->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.otherNC->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.otherNC->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.otherNC->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.otherNC->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.otherNC->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.otherNC->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.otherNC->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLNuE.otherNC->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.otherNC->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.otherNC->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.otherNC->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.otherNC->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.otherNC->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.otherNC->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.otherNC->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.otherNC->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.otherNC->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.otherNC->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 4){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.CCnumu->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.CCnumu->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.CCnumu->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.CCnumu->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.CCnumu->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.CCnumu->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.CCnumu->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.CCnumu->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.CCnumu->Fill(highestEnergy_trackscore, weight);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.CCnumu->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.CCnumu->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.CCnumu->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.CCnumu->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.CCnumu->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.CCnumu->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.CCnumu->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.CCnumu->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.CCnumu->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.CCnumu->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.CCnumu->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.CCnumu->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            
                            trackscoreHighestEnergyPFP_splitDLUboone.CCnumu->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.CCnumu->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.CCnumu->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.CCnumu->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.CCnumu->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.CCnumu += weight;
                        numEventCutWithoutWeightingDLNuE.CCnumu++;
                        sliceCompleteness_splitDLNuE.CCnumu->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.CCnumu->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.CCnumu->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.CCnumu->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.CCnumu->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.CCnumu->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.CCnumu->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.CCnumu->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLNuE.CCnumu->Fill(highestEnergy_trackscore, weight); 

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.CCnumu->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.CCnumu->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.CCnumu->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.CCnumu->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.CCnumu->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.CCnumu->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.CCnumu->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.CCnumu->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.CCnumu->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.CCnumu->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 5){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.CCnue->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.CCnue->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.CCnue->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.CCnue->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.CCnue->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.CCnue->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.CCnue->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.CCnue->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.CCnue->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.CCnue->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.CCnue->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.CCnue->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.CCnue->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.CCnue->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.CCnue->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.CCnue->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.CCnue->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.CCnue->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.CCnue->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.CCnue->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.CCnue->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.CCnue->Fill(highestEnergy_trackscore, weight); 
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.CCnue->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.CCnue->Fill(highestTrackscore, weight);

                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.CCnue->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.CCnue->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.CCnue += weight;
                        numEventCutWithoutWeightingDLNuE.CCnue++;
                        sliceCompleteness_splitDLNuE.CCnue->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.CCnue->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.CCnue->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.CCnue->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.CCnue->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.CCnue->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.CCnue->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.CCnue->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLNuE.CCnue->Fill(highestEnergy_trackscore, weight); 

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.CCnue->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.CCnue->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.CCnue->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.CCnue->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.CCnue->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.CCnue->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.CCnue->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.CCnue->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.CCnue->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.CCnue->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 6){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.dirt->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.dirt->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.dirt->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.dirt->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.dirt->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.dirt->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.dirt->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.dirt->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.dirt->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.dirt->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.dirt->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.dirt->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.dirt->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.dirt->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.dirt->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.dirt->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.dirt->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.dirt->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.dirt->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.dirt->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.dirt->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.dirt->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.dirt->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.dirt->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.dirt->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.dirt->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.dirt += weight;
                        numEventCutWithoutWeightingDLNuE.dirt++;
                        sliceCompleteness_splitDLNuE.dirt->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.dirt->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.dirt->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.dirt->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.dirt->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.dirt->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.dirt->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.dirt->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLNuE.dirt->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.dirt->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.dirt->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.dirt->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.dirt->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.dirt->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.dirt->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.dirt->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.dirt->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.dirt->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.dirt->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 7 && signal == 1){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.nu_eDirt->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.nu_eDirt->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.nu_eDirt->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.nu_eDirt->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.nu_eDirt->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.nu_eDirt->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.nu_eDirt->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.nu_eDirt->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.nu_eDirt->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.nu_eDirt->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.nu_eDirt->Fill(highestTrackscore, weight); 
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.nu_eDirt->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.nu_eDirt->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.nu_eDirt->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.nu_eDirt->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.nu_eDirt->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.nu_eDirt->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.nu_eDirt->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.nu_eDirt->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.nu_eDirt->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.nu_eDirt->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLUboone.nu_eDirt->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.nu_eDirt->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.nu_eDirt->Fill(highestTrackscore, weight); 
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.nu_eDirt->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.nu_eDirt->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.nuEDirt += weight;
                        numEventCutWithoutWeightingDLNuE.nuEDirt++;
                        sliceCompleteness_splitDLNuE.nu_eDirt->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.nu_eDirt->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.nu_eDirt->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.nu_eDirt->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.nu_eDirt->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.nu_eDirt->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.nu_eDirt->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.nu_eDirt->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitDLNuE.nu_eDirt->Fill(highestEnergy_trackscore, weight);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.nu_eDirt->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.nu_eDirt->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.nu_eDirt->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.nu_eDirt->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.nu_eDirt->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.nu_eDirt->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.nu_eDirt->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.nu_eDirt->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.nu_eDirt->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.nu_eDirt->Fill(recoVZ, weight);
                        }
                    }
                } else if(sliceEventType == 8){
                    if(DLCurrent == 2){
                        sliceCompleteness_splitBDT.other->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitBDT.other->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitBDT.other->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitBDT.other->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitBDT.other->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitBDT.other->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitBDT.other->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitBDT.other->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);

                            trackscoreHighestEnergyPFP_splitBDT.other->Fill(highestEnergy_trackscore, weight);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitBDT.other->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitBDT.other->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitBDT.other->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitBDT.other->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 0){
                        sliceCompleteness_splitDLUboone.other->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLUboone.other->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLUboone.other->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLUboone.other->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLUboone.other->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLUboone.other->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLUboone.other->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLUboone.other->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            
                            trackscoreHighestEnergyPFP_splitDLUboone.other->Fill(highestEnergy_trackscore, weight);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLUboone.other->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLUboone.other->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLUboone.other->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLUboone.other->Fill(Q2SumValue, weight);
                        }
                    } else if(DLCurrent == 5){
                        numEventCutDLNuE.other += weight;
                        numEventCutWithoutWeightingDLNuE.other++;
                        sliceCompleteness_splitDLNuE.other->Fill(reco_sliceCompleteness->at(slice), weight);
                        slicePurity_splitDLNuE.other->Fill(reco_slicePurity->at(slice), weight);
                        sliceCRUMBSScore_splitDLNuE.other->Fill(reco_sliceScore->at(slice), weight);
                        sliceNumPFPs_splitDLNuE.other->Fill(numPFPsSlice, weight);
                        sliceNumPrimaryPFPs_splitDLNuE.other->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumNeutrinos_splitDLNuE.other->Fill(numRecoNeutrinos, weight);

                        if(highestEnergy_PFPID != -999999){
                            ERecoSumThetaReco_splitDLNuE.other->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaReco_splitDLNuE.other->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            
                            trackscoreHighestEnergyPFP_splitDLNuE.other->Fill(highestEnergy_trackscore, weight);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs_splitDLNuE.other->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                    }
                                }
                            }

                            if(highestTrackscore != -999999) trackscoreHighestScorePFPs_splitDLNuE.other->Fill(highestTrackscore, weight);
                        }

                        if(Q2HighestValue != -999999){
                            QSquaredHighest_splitDLNuE.other->Fill(Q2HighestValue, weight);
                        }

                        if(Q2SumValue != -999999){    
                            QSquaredSum_splitDLNuE.other->Fill(Q2SumValue, weight);
                        }

                        if(recoVX != -999999){
                            recoX_low_splitDLNuE.other->Fill(recoVX, weight);
                            recoX_high_splitDLNuE.other->Fill(recoVX, weight);
                            recoY_low_splitDLNuE.other->Fill(recoVY, weight);
                            recoY_high_splitDLNuE.other->Fill(recoVY, weight);
                            recoZ_low_splitDLNuE.other->Fill(recoVZ, weight);
                            recoZ_high_splitDLNuE.other->Fill(recoVZ, weight);
                        }
                    }
                }

                // Filling Histograms
                if(sliceCategoryPlottingMacro == 0){
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
                        sliceNumPrimaryPFPs.currentCosmic->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.currentCosmic->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.currentCosmic->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.currentCosmic->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.currentCosmic->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.currentCosmic->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.currentCosmic->Fill(Q2SumValue, weight);
                            QSquaredSumDist.currentCosmic->Fill(Q2SumValue);
                        }

                        if(recoVX != -999999){
                            recoX.currentCosmic->Fill(recoVX, weight);
                            recoXDist.currentCosmic->Fill(recoVX);
                            recoX_smallerBins.currentCosmic->Fill(recoVX, weight);
                            recoY.currentCosmic->Fill(recoVY, weight);
                            recoYDist.currentCosmic->Fill(recoVY);
                            recoY_smallerBins.currentCosmic->Fill(recoVY, weight);
                            recoZ.currentCosmic->Fill(recoVZ, weight);
                            recoZDist.currentCosmic->Fill(recoVZ);
                            recoZ_smallerBins.currentCosmic->Fill(recoVZ, weight);
                            
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

                            trackscoreHighestEnergyPFP.currentCosmic->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentCosmic->Fill(highestEnergy_trackscore);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.currentCosmic->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.currentCosmic->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.currentCosmic->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.currentCosmic->Fill(highestTrackscore);
                            }
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
                        sliceNumPrimaryPFPs.ubooneCosmic->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.ubooneCosmic->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.ubooneCosmic->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.ubooneCosmic->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.ubooneCosmic->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.ubooneCosmic->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.ubooneCosmic->Fill(Q2SumValue, weight);
                            QSquaredSumDist.ubooneCosmic->Fill(Q2SumValue);
                        }                       
 
                        if(recoVX != -999999){
                            recoX.ubooneCosmic->Fill(recoVX, weight);
                            recoXDist.ubooneCosmic->Fill(recoVX);
                            recoX_smallerBins.ubooneCosmic->Fill(recoVX, weight);
                            recoY.ubooneCosmic->Fill(recoVY, weight);
                            recoYDist.ubooneCosmic->Fill(recoVY);
                            recoY_smallerBins.ubooneCosmic->Fill(recoVY, weight);
                            recoZ.ubooneCosmic->Fill(recoVZ, weight);
                            recoZDist.ubooneCosmic->Fill(recoVZ);
                            recoZ_smallerBins.ubooneCosmic->Fill(recoVZ, weight);
                            
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

                            trackscoreHighestEnergyPFP.ubooneCosmic->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.ubooneCosmic->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.ubooneCosmic->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.ubooneCosmic->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.ubooneCosmic->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.ubooneCosmic->Fill(highestTrackscore);
                            }
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
                        sliceNumPrimaryPFPs.nuECosmic->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.nuECosmic->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.nuECosmic->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.nuECosmic->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.nuECosmic->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.nuECosmic->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.nuECosmic->Fill(Q2SumValue, weight);
                            QSquaredSumDist.nuECosmic->Fill(Q2SumValue);
                        }                        

                        if(recoVX != -999999){
                            recoX.nuECosmic->Fill(recoVX, weight);
                            recoXDist.nuECosmic->Fill(recoVX);
                            recoX_smallerBins.nuECosmic->Fill(recoVX, weight);
                            recoY.nuECosmic->Fill(recoVY, weight);
                            recoYDist.nuECosmic->Fill(recoVY);
                            recoY_smallerBins.nuECosmic->Fill(recoVY, weight);
                            recoZ.nuECosmic->Fill(recoVZ, weight);
                            recoZDist.nuECosmic->Fill(recoVZ);
                            recoZ_smallerBins.nuECosmic->Fill(recoVZ, weight);

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
                            
                            trackscoreHighestEnergyPFP.nuECosmic->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuECosmic->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.nuECosmic->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.nuECosmic->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.nuECosmic->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.nuECosmic->Fill(highestTrackscore);
                            }
                        }
                    }

                } else if(sliceCategoryPlottingMacro == 1){
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
                        sliceNumPrimaryPFPs.currentSignal->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.currentSignal->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.currentSignal->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.currentSignal->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.currentSignal->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.currentSignal->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.currentSignal->Fill(Q2SumValue, weight);
                            QSquaredSumDist.currentSignal->Fill(Q2SumValue);
                        }

                        if(highestEnergy_PFPID != -999999){
                            // There is a PFP in the slice, fill the histograms
                            std::cout << "WEIGHT CHECK, should be 1 here - " << weight << std::endl;
                            if(weight != 1){ std::cout << "signal = " << signal << ", DLCurrent = " << DLCurrent << std::endl; std::cout << "Slice Category = " << sliceCategoryPlottingMacro << ", Slice Interaction = " << reco_sliceInteraction->at(slice) << std::endl;}
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
                            deltaEnergy.currentSignal->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy, weight);
                            deltaEnergyDist.currentSignal->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy);
                            deltaEnergySum.currentSignal->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy, weight);
                            deltaEnergySumDist.currentSignal->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy);
                      
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


                            trackscoreHighestEnergyPFP.currentSignal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignal->Fill(highestEnergy_trackscore);
                           
    
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.currentSignal->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.currentSignal->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }                     

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.currentSignal->Fill(highestTrackscore, weight);  
                                trackscoreHighestScorePFPsDist.currentSignal->Fill(highestTrackscore);  
                            }

                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceBDT->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceBDT->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceBDT->Fill(recoVZ, angleDifference);

                                xCoordEnergyAsymmetryHighestBDT->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                xCoordEnergyAsymmetrySummedBDT->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy)); 
                                yCoordEnergyAsymmetryHighestBDT->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                yCoordEnergyAsymmetrySummedBDT->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy)); 
                                zCoordEnergyAsymmetryHighestBDT->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                zCoordEnergyAsymmetrySummedBDT->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));                                
 
                                if(recoVX >= xMin && recoVX <= xMin+30){
                                    xCoordAngleDifferenceBDT_low->Fill(recoVX, angleDifference);
                                    xCoordEnergyAsymmetryHighestBDT_low->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    xCoordEnergyAsymmetrySummedBDT_low->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));                               
 
                                } else if(recoVX <= xMax && recoVX >= xMax-30){
                                    xCoordAngleDifferenceBDT_high->Fill(recoVX, angleDifference); 
                                    xCoordEnergyAsymmetryHighestBDT_high->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    xCoordEnergyAsymmetrySummedBDT_high->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy)); 
                                }

                                if(recoVY >= yMin && recoVY <= yMin+30){
                                    yCoordAngleDifferenceBDT_low->Fill(recoVY, angleDifference);
                                    yCoordEnergyAsymmetryHighestBDT_low->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    yCoordEnergyAsymmetrySummedBDT_low->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                } else if(recoVY <= yMax && recoVY >= yMax-30){
                                    yCoordAngleDifferenceBDT_high->Fill(recoVY, angleDifference);
                                    yCoordEnergyAsymmetryHighestBDT_high->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    yCoordEnergyAsymmetrySummedBDT_high->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));  
                                }

                                if(recoVZ >= zMin && recoVZ <= zMin+30){
                                    zCoordAngleDifferenceBDT_low->Fill(recoVZ, angleDifference);
                                    zCoordEnergyAsymmetryHighestBDT_low->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    zCoordEnergyAsymmetrySummedBDT_low->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                } else if(recoVZ <= zMax && recoVZ >= zMax-110){
                                    zCoordAngleDifferenceBDT_high->Fill(recoVZ, angleDifference);
                                    zCoordEnergyAsymmetryHighestBDT_high->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    zCoordEnergyAsymmetrySummedBDT_high->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy)); 
                                }
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
                            recoX_smallerBins.currentSignal->Fill(recoVX, weight);
                            recoY.currentSignal->Fill(recoVY, weight);
                            recoYDist.currentSignal->Fill(recoVY);
                            recoY_smallerBins.currentSignal->Fill(recoVY, weight);
                            recoZ.currentSignal->Fill(recoVZ, weight);
                            recoZDist.currentSignal->Fill(recoVZ);
                            recoZ_smallerBins.currentSignal->Fill(recoVZ, weight);
                            
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
                        sliceNumPrimaryPFPs.ubooneSignal->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.ubooneSignal->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.ubooneSignal->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.ubooneSignal->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.ubooneSignal->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.ubooneSignal->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.ubooneSignal->Fill(Q2SumValue, weight);
                            QSquaredSumDist.ubooneSignal->Fill(Q2SumValue);
                        }

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
                            deltaEnergy.ubooneSignal->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy, weight);
                            deltaEnergyDist.ubooneSignal->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy);
                            deltaEnergySum.ubooneSignal->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy, weight);
                            deltaEnergySumDist.ubooneSignal->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy);
                      
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

                            trackscoreHighestEnergyPFP.ubooneSignal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.ubooneSignal->Fill(highestEnergy_trackscore);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.ubooneSignal->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.ubooneSignal->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }                           

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.ubooneSignal->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.ubooneSignal->Fill(highestTrackscore);
                            }

                            if(recoVX != -999999){
                                
                                xCoordAngleDifferenceDLUboone->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLUboone->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLUboone->Fill(recoVZ, angleDifference);

                                xCoordEnergyAsymmetryHighestDLUboone->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                yCoordEnergyAsymmetryHighestDLUboone->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                zCoordEnergyAsymmetryHighestDLUboone->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                xCoordEnergyAsymmetrySummedDLUboone->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                yCoordEnergyAsymmetrySummedDLUboone->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                zCoordEnergyAsymmetrySummedDLUboone->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));                               
 
                                if(recoVX >= xMin && recoVX <= xMin+30){
                                    xCoordAngleDifferenceDLUboone_low->Fill(recoVX, angleDifference);
                                    xCoordEnergyAsymmetryHighestDLUboone_low->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    xCoordEnergyAsymmetrySummedDLUboone_low->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));                               
                                } else if(recoVX <= xMax && recoVX >= xMax-30){
                                    xCoordAngleDifferenceDLUboone_high->Fill(recoVX, angleDifference); 
                                    xCoordEnergyAsymmetryHighestDLUboone_high->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    xCoordEnergyAsymmetrySummedDLUboone_high->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                }

                                if(recoVY >= yMin && recoVY <= yMin+30){
                                    yCoordAngleDifferenceDLUboone_low->Fill(recoVY, angleDifference);
                                    yCoordEnergyAsymmetryHighestDLUboone_low->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    yCoordEnergyAsymmetrySummedDLUboone_low->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                } else if(recoVY <= yMax && recoVY >= yMax-30){
                                    yCoordAngleDifferenceDLUboone_high->Fill(recoVY, angleDifference);
                                    yCoordEnergyAsymmetryHighestDLUboone_high->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    yCoordEnergyAsymmetrySummedDLUboone_high->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy)); 
                                }

                                if(recoVZ >= zMin && recoVZ <= zMin+30){
                                    zCoordAngleDifferenceDLUboone_low->Fill(recoVZ, angleDifference);
                                    zCoordEnergyAsymmetryHighestDLUboone_low->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    zCoordEnergyAsymmetrySummedDLUboone_low->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                } else if(recoVZ <= zMax && recoVZ >= zMax-110){
                                    zCoordAngleDifferenceDLUboone_high->Fill(recoVZ, angleDifference);
                                    zCoordEnergyAsymmetryHighestDLUboone_high->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    zCoordEnergyAsymmetrySummedDLUboone_high->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                }
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
                            recoX_smallerBins.ubooneSignal->Fill(recoVX, weight);
                            recoY.ubooneSignal->Fill(recoVY, weight);
                            recoYDist.ubooneSignal->Fill(recoVY);
                            recoY_smallerBins.ubooneSignal->Fill(recoVY, weight);
                            recoZ.ubooneSignal->Fill(recoVZ, weight);
                            recoZDist.ubooneSignal->Fill(recoVZ);
                            recoZ_smallerBins.ubooneSignal->Fill(recoVZ, weight);

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
                        sliceNumPrimaryPFPs.nuESignal->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.nuESignal->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.nuESignal->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.nuESignal->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.nuESignal->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.nuESignal->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.nuESignal->Fill(Q2SumValue, weight);
                            QSquaredSumDist.nuESignal->Fill(Q2SumValue);
                        }

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
                            deltaEnergy.nuESignal->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy, weight);
                            deltaEnergyDist.nuESignal->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy);
                            deltaEnergySum.nuESignal->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy, weight);
                            deltaEnergySumDist.nuESignal->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy);
                      
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

                            trackscoreHighestEnergyPFP.nuESignal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignal->Fill(highestEnergy_trackscore);

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.nuESignal->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.nuESignal->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }                           

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.nuESignal->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.nuESignal->Fill(highestTrackscore);
                            }

                            if(recoVX != -999999){
                
                                xCoordAngleDifferenceDLNuE->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLNuE->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLNuE->Fill(recoVZ, angleDifference);

                                xCoordEnergyAsymmetryHighestDLNuE->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                yCoordEnergyAsymmetryHighestDLNuE->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                zCoordEnergyAsymmetryHighestDLNuE->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                xCoordEnergyAsymmetrySummedDLNuE->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                yCoordEnergyAsymmetrySummedDLNuE->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                zCoordEnergyAsymmetrySummedDLNuE->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));

                                if(recoVX >= xMin && recoVX <= xMin+30){
                                    xCoordAngleDifferenceDLNuE_low->Fill(recoVX, angleDifference);
                                    xCoordEnergyAsymmetryHighestDLNuE_low->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    xCoordEnergyAsymmetrySummedDLNuE_low->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                } else if(recoVX <= xMax && recoVX >= xMax-30){
                                    xCoordAngleDifferenceDLNuE_high->Fill(recoVX, angleDifference); 
                                    xCoordEnergyAsymmetryHighestDLNuE_high->Fill(recoVX, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    xCoordEnergyAsymmetrySummedDLNuE_high->Fill(recoVX, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                }

                                if(recoVY >= yMin && recoVY <= yMin+30){
                                    yCoordAngleDifferenceDLNuE_low->Fill(recoVY, angleDifference);
                                    yCoordEnergyAsymmetryHighestDLNuE_low->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    yCoordEnergyAsymmetrySummedDLNuE_low->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                
                                } else if(recoVY <= yMax && recoVY >= yMax-30){
                                    yCoordAngleDifferenceDLNuE_high->Fill(recoVY, angleDifference); 
                                    yCoordEnergyAsymmetryHighestDLNuE_high->Fill(recoVY, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    yCoordEnergyAsymmetrySummedDLNuE_high->Fill(recoVY, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                }

                                if(recoVZ >= zMin && recoVZ <= zMin+30){
                                    zCoordAngleDifferenceDLNuE_low->Fill(recoVZ, angleDifference);
                                    zCoordEnergyAsymmetryHighestDLNuE_low->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    zCoordEnergyAsymmetrySummedDLNuE_low->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy));
                                } else if(recoVZ <= zMax && recoVZ >= zMax-110){
                                    zCoordAngleDifferenceDLNuE_high->Fill(recoVZ, angleDifference);
                                    zCoordEnergyAsymmetryHighestDLNuE_high->Fill(recoVZ, ((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy));
                                    zCoordEnergyAsymmetrySummedDLNuE_high->Fill(recoVZ, ((recoilElectron_energy - summedEnergy)/recoilElectron_energy)); 
                                }
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
                            recoX_smallerBins.nuESignal->Fill(recoVX, weight);
                            recoY.nuESignal->Fill(recoVY, weight);
                            recoYDist.nuESignal->Fill(recoVY);
                            recoY_smallerBins.nuESignal->Fill(recoVY, weight);
                            recoZ.nuESignal->Fill(recoVZ, weight);
                            recoZDist.nuESignal->Fill(recoVZ);
                            recoZ_smallerBins.nuESignal->Fill(recoVZ, weight);
                            
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
                } else if(sliceCategoryPlottingMacro == 2){
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
                        sliceNumPrimaryPFPs.currentSignalFuzzy->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.currentSignalFuzzy->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.currentSignalFuzzy->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.currentSignalFuzzy->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.currentSignalFuzzy->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.currentSignalFuzzy->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.currentSignalFuzzy->Fill(Q2SumValue, weight);
                            QSquaredSumDist.currentSignalFuzzy->Fill(Q2SumValue);
                        }

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
                            deltaEnergy.currentSignalFuzzy->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy, weight);
                            deltaEnergyDist.currentSignalFuzzy->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy);
                            deltaEnergySum.currentSignalFuzzy->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy, weight);
                            deltaEnergySumDist.currentSignalFuzzy->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy);
                      
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

                            trackscoreHighestEnergyPFP.currentSignalFuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignalFuzzy->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.currentSignalFuzzy->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.currentSignalFuzzy->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }                            

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.currentSignalFuzzy->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.currentSignalFuzzy->Fill(highestTrackscore);
                            }

                            if(recoVX != -999999){
                                /*  
                                xCoordAngleDifferenceBDT->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceBDT->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceBDT->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceBDT_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceBDT_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceBDT_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceBDT_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceBDT_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceBDT_high->Fill(recoVZ, angleDifference); 
                                */ 
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
                            recoX_smallerBins.currentSignalFuzzy->Fill(recoVX, weight);
                            recoY.currentSignalFuzzy->Fill(recoVY, weight);
                            recoYDist.currentSignalFuzzy->Fill(recoVY);
                            recoY_smallerBins.currentSignalFuzzy->Fill(recoVY, weight);
                            recoZ.currentSignalFuzzy->Fill(recoVZ, weight);
                            recoZDist.currentSignalFuzzy->Fill(recoVZ);
                            recoZ_smallerBins.currentSignalFuzzy->Fill(recoVZ, weight);
                            
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
                        sliceNumPrimaryPFPs.ubooneSignalFuzzy->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.ubooneSignalFuzzy->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.ubooneSignalFuzzy->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.ubooneSignalFuzzy->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.ubooneSignalFuzzy->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.ubooneSignalFuzzy->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.ubooneSignalFuzzy->Fill(Q2SumValue, weight);
                            QSquaredSumDist.ubooneSignalFuzzy->Fill(Q2SumValue);
                        }

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
                            deltaEnergy.ubooneSignalFuzzy->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy, weight);
                            deltaEnergyDist.ubooneSignalFuzzy->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy);
                            deltaEnergySum.ubooneSignalFuzzy->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy, weight);
                            deltaEnergySumDist.ubooneSignalFuzzy->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy);
                      
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

                            trackscoreHighestEnergyPFP.ubooneSignalFuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.ubooneSignalFuzzy->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.ubooneSignalFuzzy->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.ubooneSignalFuzzy->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }                       

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.ubooneSignalFuzzy->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.ubooneSignalFuzzy->Fill(highestTrackscore);
                            }

                            if(recoVX != -999999){
                                /* 
                                xCoordAngleDifferenceDLUboone->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLUboone->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLUboone->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceDLUboone_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceDLUboone_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceDLUboone_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceDLUboone_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceDLUboone_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceDLUboone_high->Fill(recoVZ, angleDifference);
                                */ 
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
                            recoX_smallerBins.ubooneSignalFuzzy->Fill(recoVX, weight);
                            recoY.ubooneSignalFuzzy->Fill(recoVY, weight);
                            recoYDist.ubooneSignalFuzzy->Fill(recoVY);
                            recoY_smallerBins.ubooneSignalFuzzy->Fill(recoVY, weight);
                            recoZ.ubooneSignalFuzzy->Fill(recoVZ, weight);
                            recoZDist.ubooneSignalFuzzy->Fill(recoVZ);
                            recoZ_smallerBins.ubooneSignalFuzzy->Fill(recoVZ, weight);
                            
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
                        sliceNumPrimaryPFPs.nuESignalFuzzy->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.nuESignalFuzzy->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.nuESignalFuzzy->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.nuESignalFuzzy->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.nuESignalFuzzy->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.nuESignalFuzzy->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.nuESignalFuzzy->Fill(Q2SumValue, weight);
                            QSquaredSumDist.nuESignalFuzzy->Fill(Q2SumValue);
                        }

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
                            deltaEnergy.nuESignalFuzzy->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy, weight);
                            deltaEnergyDist.nuESignalFuzzy->Fill((recoilElectron_energy - highestEnergy_energy)/recoilElectron_energy);
                            deltaEnergySum.nuESignalFuzzy->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy, weight);
                            deltaEnergySumDist.nuESignalFuzzy->Fill((recoilElectron_energy - summedEnergy)/recoilElectron_energy);
                      
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

                            trackscoreHighestEnergyPFP.nuESignalFuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignalFuzzy->Fill(highestEnergy_trackscore);
                            
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.nuESignalFuzzy->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.nuESignalFuzzy->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }                           

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.nuESignalFuzzy->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.nuESignalFuzzy->Fill(highestTrackscore);
                            }

                            if(recoVX != -999999){
                                /* // Commented out bc it is signal slice with completeness < 0.5 
                                xCoordAngleDifferenceDLNuE->Fill(recoVX, angleDifference);
                                yCoordAngleDifferenceDLNuE->Fill(recoVY, angleDifference);
                                zCoordAngleDifferenceDLNuE->Fill(recoVZ, angleDifference);
                                
                                if(recoVX >= xMin && recoVX <= xMin+20) xCoordAngleDifferenceDLNuE_low->Fill(recoVX, angleDifference);
                                else if(recoVX <= xMax && recoVX >= xMax-20) xCoordAngleDifferenceDLNuE_high->Fill(recoVX, angleDifference); 
                            
                                if(recoVY >= yMin && recoVY <= yMin+20) yCoordAngleDifferenceDLNuE_low->Fill(recoVY, angleDifference);
                                else if(recoVY <= yMax && recoVY >= yMax-20) yCoordAngleDifferenceDLNuE_high->Fill(recoVY, angleDifference); 
                                
                                if(recoVZ >= zMin && recoVZ <= zMin+20) zCoordAngleDifferenceDLNuE_low->Fill(recoVZ, angleDifference);
                                else if(recoVZ <= zMax && recoVZ >= zMax-40) zCoordAngleDifferenceDLNuE_high->Fill(recoVZ, angleDifference); 
                                */                            
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
                            recoX_smallerBins.nuESignalFuzzy->Fill(recoVX, weight);
                            recoY.nuESignalFuzzy->Fill(recoVY, weight);
                            recoYDist.nuESignalFuzzy->Fill(recoVY);
                            recoY_smallerBins.nuESignalFuzzy->Fill(recoVY, weight);
                            recoZ.nuESignalFuzzy->Fill(recoVZ, weight);
                            recoZDist.nuESignalFuzzy->Fill(recoVZ);
                            recoZ_smallerBins.nuESignalFuzzy->Fill(recoVZ, weight);
                            
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
                } else if(sliceCategoryPlottingMacro == 3){
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
                        sliceNumPrimaryPFPs.currentBNB->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.currentBNB->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.currentBNB->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.currentBNB->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.currentBNB->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.currentBNB->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.currentBNB->Fill(Q2SumValue, weight);
                            QSquaredSumDist.currentBNB->Fill(Q2SumValue);
                        }

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

                            trackscoreHighestEnergyPFP.currentBNB->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentBNB->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.currentBNB->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.currentBNB->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.currentBNB->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.currentBNB->Fill(highestTrackscore);
                            }
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
                            recoX_smallerBins.currentBNB->Fill(recoVX, weight);
                            recoY.currentBNB->Fill(recoVY, weight);
                            recoYDist.currentBNB->Fill(recoVY);
                            recoY_smallerBins.currentBNB->Fill(recoVY, weight);
                            recoZ.currentBNB->Fill(recoVZ, weight);
                            recoZDist.currentBNB->Fill(recoVZ);
                            recoZ_smallerBins.currentBNB->Fill(recoVZ, weight);
                            
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
                        sliceNumPrimaryPFPs.ubooneBNB->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.ubooneBNB->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.ubooneBNB->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.ubooneBNB->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.ubooneBNB->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.ubooneBNB->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.ubooneBNB->Fill(Q2SumValue, weight);
                            QSquaredSumDist.ubooneBNB->Fill(Q2SumValue);
                        }

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

                            trackscoreHighestEnergyPFP.ubooneBNB->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.ubooneBNB->Fill(highestEnergy_trackscore);
                           

                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.ubooneBNB->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.ubooneBNB->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.ubooneBNB->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.ubooneBNB->Fill(highestTrackscore);
                            }
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
                            recoX_smallerBins.ubooneBNB->Fill(recoVX, weight);
                            recoY.ubooneBNB->Fill(recoVY, weight);
                            recoYDist.ubooneBNB->Fill(recoVY);
                            recoY_smallerBins.ubooneBNB->Fill(recoVY, weight);
                            recoZ.ubooneBNB->Fill(recoVZ, weight);
                            recoZDist.ubooneBNB->Fill(recoVZ);
                            recoZ_smallerBins.ubooneBNB->Fill(recoVZ, weight);
                            
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
                        sliceNumPrimaryPFPs.nuEBNB->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.nuEBNB->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.nuEBNB->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.nuEBNB->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.nuEBNB->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.nuEBNB->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.nuEBNB->Fill(Q2SumValue, weight);
                            QSquaredSumDist.nuEBNB->Fill(Q2SumValue);
                        }

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

                            trackscoreHighestEnergyPFP.nuEBNB->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuEBNB->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.nuEBNB->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.nuEBNB->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.nuEBNB->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.nuEBNB->Fill(highestTrackscore);
                            }
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
                            recoX_smallerBins.nuEBNB->Fill(recoVX, weight);
                            recoY.nuEBNB->Fill(recoVY, weight);
                            recoYDist.nuEBNB->Fill(recoVY);
                            recoY_smallerBins.nuEBNB->Fill(recoVY, weight);
                            recoZ.nuEBNB->Fill(recoVZ, weight);
                            recoZDist.nuEBNB->Fill(recoVZ);
                            recoZ_smallerBins.nuEBNB->Fill(recoVZ, weight);
                            
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
                } else if(sliceCategoryPlottingMacro == 4){
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
                        sliceNumPrimaryPFPs.currentBNBFuzzy->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.currentBNBFuzzy->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.currentBNBFuzzy->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.currentBNBFuzzy->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.currentBNBFuzzy->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.currentBNBFuzzy->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.currentBNBFuzzy->Fill(Q2SumValue, weight);
                            QSquaredSumDist.currentBNBFuzzy->Fill(Q2SumValue);
                        }

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

                            trackscoreHighestEnergyPFP.currentBNBFuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentBNBFuzzy->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.currentBNBFuzzy->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.currentBNBFuzzy->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.currentBNBFuzzy->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.currentBNBFuzzy->Fill(highestTrackscore);
                            }
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
                            recoX_smallerBins.currentBNBFuzzy->Fill(recoVX, weight);
                            recoY.currentBNBFuzzy->Fill(recoVY, weight);
                            recoYDist.currentBNBFuzzy->Fill(recoVY);
                            recoY_smallerBins.currentBNBFuzzy->Fill(recoVY, weight);
                            recoZ.currentBNBFuzzy->Fill(recoVZ, weight);
                            recoZDist.currentBNBFuzzy->Fill(recoVZ);
                            recoZ_smallerBins.currentBNBFuzzy->Fill(recoVZ, weight);

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
                        sliceNumPrimaryPFPs.ubooneBNBFuzzy->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.ubooneBNBFuzzy->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.ubooneBNBFuzzy->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.ubooneBNBFuzzy->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.ubooneBNBFuzzy->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.ubooneBNBFuzzy->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.ubooneBNBFuzzy->Fill(Q2SumValue, weight);
                            QSquaredSumDist.ubooneBNBFuzzy->Fill(Q2SumValue);
                        }

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

                            trackscoreHighestEnergyPFP.ubooneBNBFuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.ubooneBNBFuzzy->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.ubooneBNBFuzzy->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.ubooneBNBFuzzy->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.ubooneBNBFuzzy->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.ubooneBNBFuzzy->Fill(highestTrackscore);
                            }
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
                            recoX_smallerBins.ubooneBNBFuzzy->Fill(recoVX, weight);
                            recoY.ubooneBNBFuzzy->Fill(recoVY, weight);
                            recoYDist.ubooneBNBFuzzy->Fill(recoVY);
                            recoY_smallerBins.ubooneBNBFuzzy->Fill(recoVY, weight);
                            recoZ.ubooneBNBFuzzy->Fill(recoVZ, weight);
                            recoZDist.ubooneBNBFuzzy->Fill(recoVZ);
                            recoZ_smallerBins.ubooneBNBFuzzy->Fill(recoVZ, weight);
                            
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
                        sliceNumPrimaryPFPs.nuEBNBFuzzy->Fill(numPrimaryPFPsSlice, weight);
                        sliceNumPrimaryPFPsDist.nuEBNBFuzzy->Fill(numPrimaryPFPsSlice);
                        sliceNumNeutrinos.nuEBNBFuzzy->Fill(numRecoNeutrinos, weight);
                        sliceNumNeutrinosDist.nuEBNBFuzzy->Fill(numRecoNeutrinos);

                        if(Q2HighestValue != -999999){
                            QSquaredHighest.nuEBNBFuzzy->Fill(Q2HighestValue, weight);
                            QSquaredHighestDist.nuEBNBFuzzy->Fill(Q2HighestValue);
                        }

                        if(Q2SumValue != -999999){
                            QSquaredSum.nuEBNBFuzzy->Fill(Q2SumValue, weight);
                            QSquaredSumDist.nuEBNBFuzzy->Fill(Q2SumValue);
                        }

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

                            trackscoreHighestEnergyPFP.nuEBNBFuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuEBNBFuzzy->Fill(highestEnergy_trackscore);
                           
                            for(size_t pfpTrack = 0; pfpTrack < reco_particlePDG->size(); ++pfpTrack){
                                if(reco_particleSliceID->at(pfpTrack) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfpTrack) != -999999){
                                        trackscoreAllPFPs.nuEBNBFuzzy->Fill(reco_particleTrackScore->at(pfpTrack), weight);
                                        trackscoreAllPFPsDist.nuEBNBFuzzy->Fill(reco_particleTrackScore->at(pfpTrack));
                                    }
                                }
                            }

                            if(highestTrackscore != -999999){
                                trackscoreHighestScorePFPs.nuEBNBFuzzy->Fill(highestTrackscore, weight);
                                trackscoreHighestScorePFPsDist.nuEBNBFuzzy->Fill(highestTrackscore);
                            }
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
                            recoX_smallerBins.nuEBNBFuzzy->Fill(recoVX, weight);
                            recoY.nuEBNBFuzzy->Fill(recoVY, weight);
                            recoYDist.nuEBNBFuzzy->Fill(recoVY);
                            recoY_smallerBins.nuEBNBFuzzy->Fill(recoVY, weight);
                            recoZ.nuEBNBFuzzy->Fill(recoVZ, weight);
                            recoZDist.nuEBNBFuzzy->Fill(recoVZ);
                            recoZ_smallerBins.nuEBNBFuzzy->Fill(recoVZ, weight);
                            
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

    styleDrawAll(trueETheta2, 999, 999, 999, 999, (base_path + "trueETheta2_weighted.pdf").c_str(), "bottomRight", &drawLine, &right, true, true, false, false, false, false, true, true);
    
    styleDrawAll(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(sliceCompletenessDist, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(slicePurityDist, 999, 999, 999, 999, (base_path + "slicePurity_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(sliceCRUMBSScore, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(sliceCRUMBSScoreDist, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(sliceCRUMBSScore, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_BackSig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(sliceNumPFPsDist, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(sliceNumPrimaryPFPs, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(sliceNumPrimaryPFPsDist, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(sliceNumPrimaryPFPs, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_BackSig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(sliceNumNeutrinos, 999, 999, 999, 999, (base_path + "sliceNumNeutrinos_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(sliceNumNeutrinosDist, 999, 999, 999, 999, (base_path + "sliceNumNeutrinos_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(sliceNumNeutrinos, 999, 999, 999, 999, (base_path + "sliceNumNeutrinos_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);

    styleDrawAll(QSquaredHighest, 999, 999, 999, 999, (base_path + "QSquared_highest_all_lower_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(QSquaredHighestDist, 999, 999, 999, 999, (base_path + "QSquared_highest_all_lower_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawBackSig(QSquaredHighest, 999, 999, 999, 999, (base_path + "QSquared_highest_Backsig_lower_weighted.pdf").c_str(), "topRight", true, true, true, true);
    
    styleDrawAll(QSquaredSum, 999, 999, 999, 999, (base_path + "QSquared_sum_all_lower_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(QSquaredSumDist, 999, 999, 999, 999, (base_path + "QSquared_sum_all_lower_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawBackSig(QSquaredSum, 999, 999, 999, 999, (base_path + "QSquared_sum_Backsig_lower_weighted.pdf").c_str(), "topRight", true, true, true, true);

    styleDrawAll(trackscoreHighestEnergyPFP, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(trackscoreHighestEnergyPFPDist, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawBackSig(trackscoreHighestEnergyPFP, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_Backsig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(trackscoreAllPFPs, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(trackscoreAllPFPsDist, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawBackSig(trackscoreAllPFPs, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_Backsig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    styleDrawAll(trackscoreHighestScorePFPs, 999, 999, 999, 999, (base_path + "trackscoreHighestScorePFPs_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(trackscoreHighestScorePFPsDist, 999, 999, 999, 999, (base_path + "trackscoreHighestScorePFPs_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawBackSig(trackscoreHighestScorePFPs, 999, 999, 999, 999, (base_path + "trackscoreHighestScorePFPs_Backsig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    
    styleDrawAll(trackscoreAllPFPsPFP, 999, 999, 999, 999, (base_path + "trackscoreAllPFPsPFP_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(trackscoreAllPFPsPFP, 999, 999, 999, 999, (base_path + "trackscoreAllPFPsPFP_Backsig_weighted.pdf").c_str(), "topRight", false, false, true, true);

    styleDrawAll(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(ERecoSumThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);
    styleDrawAll(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(ERecoHighestThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);

    styleDrawAll(ETrueThetaReco, 999, 999, 999, 999, (base_path + "ETrueThetaReco_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);
    styleDrawAll(ETrueThetaRecoDist, 999, 999, 999, 999, (base_path + "ETrueThetaReco_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ERecoSumThetaTrue, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ERecoSumThetaTrueDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ERecoHighestThetaTrue, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ERecoHighestThetaTrueDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ETrue, 999, 999, 999, 999, (base_path + "ETrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ETrueDist, 999, 999, 999, 999, (base_path + "ETrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ThetaTrue, 999, 999, 999, 999, (base_path + "ThetaTrue_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(ThetaTrueDist, 999, 999, 999, 999, (base_path + "ThetaTrue_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);

    styleDrawAll(deltaX, 999, 999, 999, 999, (base_path + "deltaX_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaXDist, 999, 999, 999, 999, (base_path + "deltaX_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true);
    styleDrawAll(deltaXDist, 0, 32000, 999, 999, (base_path + "deltaX_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaXDist, 0, 35000, 999, 999, (base_path + "deltaX_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaXDist, 0, 9000, 999, 999, (base_path + "deltaX_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaXDist, 0, 4000, 999, 999, (base_path + "deltaX_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);
    styleDrawBackSig(deltaX, 999, 999, 999, 999, (base_path + "deltaX_BackSig_weighted.pdf").c_str(), "topRight", true, false, true, true);

    styleDrawAll(deltaY, 999, 999, 999, 999, (base_path + "deltaY_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaYDist, 999, 999, 999, 999, (base_path + "deltaY_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true);
    styleDrawBackSig(deltaY, 999, 999, 999, 999, (base_path + "deltaY_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaYDist, 0, 19000, 999, 999, (base_path + "deltaY_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaYDist, 0, 37000, 999, 999, (base_path + "deltaY_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaYDist, 0, 8000, 999, 999, (base_path + "deltaY_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaYDist, 0, 4000, 999, 999, (base_path + "deltaY_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);
    
    styleDrawAll(deltaZ, 999, 999, 999, 999, (base_path + "deltaZ_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaZDist, 999, 999, 999, 999, (base_path + "deltaZ_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true);
    styleDrawBackSig(deltaZ, 999, 999, 999, 999, (base_path + "deltaZ_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaZDist, 0, 20000, 999, 999, (base_path + "deltaZ_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaZDist, 0, 36000, 999, 999, (base_path + "deltaZ_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaZDist, 0, 9000, 999, 999, (base_path + "deltaZ_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaZDist, 0, 4000, 999, 999, (base_path + "deltaZ_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);
    
    styleDrawAll(deltaR, 999, 999, 999, 999, (base_path + "deltaR_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true, true);
    styleDrawAll(deltaRDist, 999, 999, 999, 999, (base_path + "deltaR_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, false, true, true, true);
    styleDrawBackSig(deltaR, 999, 999, 999, 999, (base_path + "deltaR_BackSig_weighted.pdf").c_str(), "topRight", true, true, true, true);
    styleDrawAll(deltaRDist, 0, 25000, 999, 999, (base_path + "deltaR_signalBDT_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, false, true);
    styleDrawAll(deltaRDist, 0, 69000, 999, 999, (base_path + "deltaR_signalDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, false, true, false);
    styleDrawAll(deltaRDist, 0, 13000, 999, 999, (base_path + "deltaR_BNBBDT_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, false, true);
    styleDrawAll(deltaRDist, 0, 5000, 999, 999, (base_path + "deltaR_BNBDLNuE_dist.pdf").c_str(), "topRight", nullptr, &right, false, false, true, true, false, false, true, false);

    styleDrawAll(recoX, 999, 999, 999, 999, (base_path + "recoX_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(recoXDist, 999, 999, 999, 999, (base_path + "recoX_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoX, 999, 999, 999, 999, (base_path + "recoX_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);
    efficiency(recoX, 0, 1, -202, -190, (base_path + "recoX_right").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoX, 0, 1, 190, 202, (base_path + "recoX_left").c_str(), "topLeft", nullptr, &right, 1, txtFileName);
    styleDrawAll(recoX_low, 999, 999, 999, 999, (base_path + "recoX_low_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(recoXDist_low, 999, 999, 999, 999, (base_path + "recoX_low_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoX_low, 999, 999, 999, 999, (base_path + "recoX_low_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);
    styleDrawAll(recoX_high, 999, 999, 999, 999, (base_path + "recoX_high_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(recoXDist_high, 999, 999, 999, 999, (base_path + "recoX_high_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoX_high, 999, 999, 999, 999, (base_path + "recoX_high_BackSig_weighted.pdf").c_str(), "bottomLeft", true, true, true, true);
    styleDrawAll(recoX_smallerBins, 999, 999, 999, 999, (base_path + "recoX_smallerBins_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoX_smallerBins, 999, 999, 999, 999, (base_path + "recoX_smallerBins_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);
    efficiency(recoX_smallerBins, 0, 1, -202, -190, (base_path + "recoX_smallerBins_right").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoX_smallerBins, 0, 1, 190, 202, (base_path + "recoX_smallerBins_left").c_str(), "topLeft", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoY, 999, 999, 999, 999, (base_path + "recoY_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(recoYDist, 999, 999, 999, 999, (base_path + "recoY_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoY, 999, 999, 999, 999, (base_path + "recoY_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);
    efficiency(recoY, 0, 1, -205, 190, (base_path + "recoY_right").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoY, 0, 1, 140, 205, (base_path + "recoY_left").c_str(), "topLeft", nullptr, &right, 1, txtFileName);
    styleDrawAll(recoY_low, 999, 999, 999, 999, (base_path + "recoY_low_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(recoYDist_low, 999, 999, 999, 999, (base_path + "recoY_low_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoY_low, 999, 999, 999, 999, (base_path + "recoY_low_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(recoY_high, 999, 999, 999, 999, (base_path + "recoY_high_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(recoYDist_high, 999, 999, 999, 999, (base_path + "recoY_high_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoY_high, 999, 999, 999, 999, (base_path + "recoY_high_BackSig_weighted.pdf").c_str(), "bottomLeft", true, true, true, true);
    styleDrawAll(recoY_smallerBins, 999, 999, 999, 999, (base_path + "recoY_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoY_smallerBins, 999, 999, 999, 999, (base_path + "recoY_BackSig_weighted.pdf").c_str(), "bottomRight", false, false, true, true);
    efficiency(recoY_smallerBins, 0, 1, -205, -190, (base_path + "recoY_smallerBins_right").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoY_smallerBins, 0, 1, 140, 205, (base_path + "recoY_smallerBins_left").c_str(), "topLeft", nullptr, &right, 1, txtFileName);
    
    styleDrawAll(recoZ, 999, 999, 999, 999, (base_path + "recoZ_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawAll(recoZDist, 999, 999, 999, 999, (base_path + "recoZ_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoZ, 999, 999, 999, 999, (base_path + "recoZ_BackSig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    efficiency(recoZ, 0, 1, 0, 20, (base_path + "recoZ_right").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoZ, 0, 1, 490, 510, (base_path + "recoZ_left").c_str(), "topLeft", nullptr, &right, 1, txtFileName);
    styleDrawAll(recoZ_low, 999, 999, 999, 999, (base_path + "recoZ_low_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(recoZDist_low, 999, 999, 999, 999, (base_path + "recoZ_low_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoZ_low, 999, 999, 999, 999, (base_path + "recoZ_low_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(recoZ_high, 999, 999, 999, 999, (base_path + "recoZ_high_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(recoZDist_high, 999, 999, 999, 999, (base_path + "recoZ_high_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(recoZ_high, 999, 999, 999, 999, (base_path + "recoZ_high_BackSig_weighted.pdf").c_str(), "bottomLeft", true, true, true, true);
    styleDrawAll(recoZ_smallerBins, 999, 999, 999, 999, (base_path + "recoZ_smallerBins_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, false, true, false, true);
    styleDrawBackSig(recoZ_smallerBins, 999, 999, 999, 999, (base_path + "recoZ_smallerBins_BackSig_weighted.pdf").c_str(), "topRight", false, false, true, true);
    efficiency(recoZ_smallerBins, 0, 1, 0, 20, (base_path + "recoZ_smallerBins_right").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoZ_smallerBins, 0, 1, 490, 510, (base_path + "recoZ_smallerBins_left").c_str(), "topLeft", nullptr, &right, 1, txtFileName);

    styleDrawAll(deltaTheta, 999, 999, 999, 999, (base_path + "deltaTheta_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true, false);
    styleDrawAll(deltaThetaDist, 999, 999, 999, 999, (base_path + "deltaTheta_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, false, true, true, true);
    styleDrawAll(deltaEnergy, 999, 999, 999, 999, (base_path + "deltaEnergy_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, false, false, true, true, true);

    styleDrawAll(pfpCompleteness, 999, 999, 999, 999, (base_path + "pfpCompleteness_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(pfpCompletenessDist, 999, 999, 999, 999, (base_path + "pfpCompleteness_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(pfpCompleteness, 999, 999, 999, 999, (base_path + "pfpCompleteness_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);
    styleDrawAll(pfpPurity, 999, 999, 999, 999, (base_path + "pfpPurity_all_weighted.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true, true);
    styleDrawAll(pfpPurityDist, 999, 999, 999, 999, (base_path + "pfpPurity_all_dist.pdf").c_str(), "topRight", nullptr, &right, true, true, true, true, true, true, true, true);
    styleDrawBackSig(pfpPurity, 999, 999, 999, 999, (base_path + "pfpPurity_BackSig_weighted.pdf").c_str(), "bottomRight", true, true, true, true);

    //Test
    //efficiency(ERecoSumThetaReco, 0, 1, 999, 999, (base_path + "ERecoSumThetaReco").c_str(), "bottomRight", nullptr, &right, 1); 
    
    efficiency(sliceCompleteness, 0, 1, 999, 999, (base_path + "sliceCompleteness").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(slicePurity, 0, 1, 999, 999, (base_path + "slicePurity").c_str(), "topRight", nullptr, &right, -1, txtFileName);
    efficiency(sliceCRUMBSScore, 0, 1, -1, 0.8, (base_path + "sliceCRUMBSScoreNegative").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(sliceCRUMBSScore, 0, 1, 0.1, 0.8, (base_path + "sliceCRUMBSScorePositive").c_str(), "bottomLeft", nullptr, &right, 1, txtFileName);
    efficiency(sliceNumPFPs, 0, 1, 999, 999, (base_path + "sliceNumPFPsNegative").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);
    efficiency(sliceNumPFPs, 0, 1, 999, 999, (base_path + "sliceNumPFPsPositive").c_str(), "bottomRight", nullptr, &right, -1, txtFileName);
    efficiency(sliceNumPrimaryPFPs, 0, 1, 999, 999, (base_path + "sliceNumPrimaryPFPsNegative").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);
    efficiency(sliceNumPrimaryPFPs, 0, 1, 999, 999, (base_path + "sliceNumPrimaryPFPsPositive").c_str(), "bottomRight", nullptr, &right, -1, txtFileName);
    efficiency(sliceNumNeutrinos, 0, 1, 0, 2, (base_path + "sliceNumNeutrinos").c_str(), "bottomRight", nullptr, &right, -1, txtFileName);

    efficiency(ERecoSumThetaReco, 0, 1, 999, 999, (base_path + "ERecoSumThetaReco").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);
    efficiency(ERecoHighestThetaReco, 0, 1, 999, 999, (base_path + "ERecoHighestThetaReco").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);

    efficiency(pfpCompleteness, 0, 1, 999, 999, (base_path + "pfpCompleteness").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(pfpPurity, 0, 1, 999, 999, (base_path + "pfpPurity").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);

    efficiency(trackscoreHighestEnergyPFP, 0, 1, 999, 999, (base_path + "trackscoreHighestEnergyPFP").c_str(), "bottomLeft", nullptr, &right, 1, txtFileName);
    efficiency(trackscoreAllPFPs, 0, 1, 999, 999, (base_path + "trackscoreAllPFPs").c_str(), "bottomLeft", nullptr, &right, 1, txtFileName);
    efficiency(trackscoreHighestScorePFPs, 0, 1, 999, 999, (base_path + "trackscoreHighestScorePFPs").c_str(), "bottomLeft", nullptr, &right, 1, txtFileName);
    
    efficiency(trackscoreAllPFPsPFP, 0, 1, 999, 999, (base_path + "trackscoreAllPFPsPFP").c_str(), "bottomLeft", nullptr, &right, 1, txtFileName);

    efficiency(recoX_low, 0, 1, 999, 999, (base_path + "recoX_low").c_str(), "bottomLeft", nullptr, &right, -1, txtFileName);
    efficiency(recoY_low, 0, 1, 999, 999, (base_path + "recoY_low").c_str(), "bottomRight", nullptr, &right, -1, txtFileName);
    efficiency(recoZ_low, 0, 1, 999, 999, (base_path + "recoZ_low").c_str(), "bottomRight", nullptr, &right, -1, txtFileName);
    
    efficiency(recoX_high, 0, 1, 999, 999, (base_path + "recoX_high").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);
    efficiency(recoY_high, 0, 1, 999, 999, (base_path + "recoY_high").c_str(), "bottomLeft", nullptr, &right, 1, txtFileName);
    efficiency(recoZ_high, 0, 1, 999, 999, (base_path + "recoZ_high").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);

    efficiency(QSquaredHighest, 0, 1, 999, 999, (base_path + "QSquared_highest_lower").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);
    efficiency(QSquaredSum, 0, 1, 999, 999, (base_path + "QSquared_sum_lower").c_str(), "bottomRight", nullptr, &right, 1, txtFileName);

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

    TwoDHistDraw(xCoordEnergyAsymmetryHighestBDT, (base_path + "energyAsym_Highest_x_BDT.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestBDT, (base_path + "energyAsym_Highest_y_BDT.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestBDT, (base_path + "energyAsym_Highest_z_BDT.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetryHighestBDT_low, (base_path + "energyAsym_Highest_x_BDT_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestBDT_low, (base_path + "energyAsym_Highest_y_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestBDT_low, (base_path + "energyAsym_Highest_z_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetryHighestBDT_high, (base_path + "energyAsym_Highest_x_BDT_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestBDT_high, (base_path + "energyAsym_Highest_y_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestBDT_high, (base_path + "energyAsym_Highest_z_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedBDT, (base_path + "energyAsym_Summed_x_BDT.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedBDT, (base_path + "energyAsym_Summed_y_BDT.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedBDT, (base_path + "energyAsym_Summed_z_BDT.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedBDT_low, (base_path + "energyAsym_Summed_x_BDT_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedBDT_low, (base_path + "energyAsym_Summed_y_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedBDT_low, (base_path + "energyAsym_Summed_z_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedBDT_high, (base_path + "energyAsym_Summed_x_BDT_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedBDT_high, (base_path + "energyAsym_Summed_y_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedBDT_high, (base_path + "energyAsym_Summed_z_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    
    TwoDHistDraw(xCoordEnergyAsymmetryHighestDLUboone, (base_path + "energyAsym_Highest_x_DLUboone.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestDLUboone, (base_path + "energyAsym_Highest_y_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestDLUboone, (base_path + "energyAsym_Highest_z_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetryHighestDLUboone_low, (base_path + "energyAsym_Highest_x_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestDLUboone_low, (base_path + "energyAsym_Highest_y_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestDLUboone_low, (base_path + "energyAsym_Highest_z_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetryHighestDLUboone_high, (base_path + "energyAsym_Highest_x_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestDLUboone_high, (base_path + "energyAsym_Highest_y_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestDLUboone_high, (base_path + "energyAsym_Highest_z_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedDLUboone, (base_path + "energyAsym_Summed_x_DLUboone.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedDLUboone, (base_path + "energyAsym_Summed_y_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedDLUboone, (base_path + "energyAsym_Summed_z_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedDLUboone_low, (base_path + "energyAsym_Summed_x_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedDLUboone_low, (base_path + "energyAsym_Summed_y_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedDLUboone_low, (base_path + "energyAsym_Summed_z_DLUboone_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedDLUboone_high, (base_path + "energyAsym_Summed_x_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedDLUboone_high, (base_path + "energyAsym_Summed_y_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedDLUboone_high, (base_path + "energyAsym_Summed_z_DLUboone_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    
    TwoDHistDraw(xCoordEnergyAsymmetryHighestDLNuE, (base_path + "energyAsym_Highest_x_DLNuE.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestDLNuE, (base_path + "energyAsym_Highest_y_DLNuE.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestDLNuE, (base_path + "energyAsym_Highest_z_DLNuE.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetryHighestDLNuE_low, (base_path + "energyAsym_Highest_x_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestDLNuE_low, (base_path + "energyAsym_Highest_y_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestDLNuE_low, (base_path + "energyAsym_Highest_z_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetryHighestDLNuE_high, (base_path + "energyAsym_Highest_x_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetryHighestDLNuE_high, (base_path + "energyAsym_Highest_y_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetryHighestDLNuE_high, (base_path + "energyAsym_Highest_z_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Highest Energy PFP: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedDLNuE, (base_path + "energyAsym_Summed_x_DLNuE.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedDLNuE, (base_path + "energyAsym_Summed_y_DLNuE.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedDLNuE, (base_path + "energyAsym_Summed_z_DLNuE.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedDLNuE_low, (base_path + "energyAsym_Summed_x_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedDLNuE_low, (base_path + "energyAsym_Summed_y_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedDLNuE_low, (base_path + "energyAsym_Summed_z_DLNuE_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(xCoordEnergyAsymmetrySummedDLNuE_high, (base_path + "energyAsym_Summed_x_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex X Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(yCoordEnergyAsymmetrySummedDLNuE_high, (base_path + "energyAsym_Summed_y_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Energy Asymmetry");
    TwoDHistDraw(zCoordEnergyAsymmetrySummedDLNuE_high, (base_path + "energyAsym_Summed_z_DLNuE_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Energy Asymmetry of Summed PFP Energy: DL Nu+E Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Energy Asymmetry");

    // Plotting Split Histograms
    // BDT Vertexing
    styleDrawSplit(sliceCompleteness_splitBDT, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(slicePurity_splitBDT, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceCRUMBSScore_splitBDT, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumPFPs_splitBDT, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumPrimaryPFPs_splitBDT, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumNeutrinos_splitBDT, 999, 999, 999, 999, (base_path + "sliceNumNeutrinos_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(QSquaredHighest_splitBDT, 999, 999, 999, 999, (base_path + "QSquared_highest_all_lower_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(QSquaredSum_splitBDT, 999, 999, 999, 999, (base_path + "QSquared_sum_all_lower_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(ERecoSumThetaReco_splitBDT, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(ERecoHighestThetaReco_splitBDT, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreHighestEnergyPFP_splitBDT, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreAllPFPs_splitBDT, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreHighestScorePFPs_splitBDT, 999, 999, 999, 999, (base_path + "trackscoreHighestScorePFPs_all_weighted_splitBDT.pdf").c_str(), "topRight", nullptr, &right, true);

    styleDrawPFPSplit(trackscoreAllPFPs_splitPFPBDT, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted_splitPFPBDT.pdf").c_str(), "topRight", nullptr, &right, true);

    // DL Uboone Vertexing
    styleDrawSplit(sliceCompleteness_splitDLUboone, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(slicePurity_splitDLUboone, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceCRUMBSScore_splitDLUboone, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumPFPs_splitDLUboone, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumPrimaryPFPs_splitDLUboone, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumNeutrinos_splitDLUboone, 999, 999, 999, 999, (base_path + "sliceNumNeutrinos_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(QSquaredHighest_splitDLUboone, 999, 999, 999, 999, (base_path + "QSquared_highest_all_lower_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(QSquaredSum_splitDLUboone, 999, 999, 999, 999, (base_path + "QSquared_sum_all_lower_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(ERecoSumThetaReco_splitDLUboone, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(ERecoHighestThetaReco_splitDLUboone, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreHighestEnergyPFP_splitDLUboone, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreAllPFPs_splitDLUboone, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreHighestScorePFPs_splitDLUboone, 999, 999, 999, 999, (base_path + "trackscoreHighestScorePFPs_all_weighted_splitDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
    
    styleDrawPFPSplit(trackscoreAllPFPs_splitPFPDLUboone, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted_splitPFPDLUboone.pdf").c_str(), "topRight", nullptr, &right, true);
        
    // DL Nu+E Vertexing
    styleDrawSplit(sliceCompleteness_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCompleteness_all_weighted_splitDLNuE.pdf").c_str(), "bottomRight", nullptr, &right, true);
    styleDrawSplit(slicePurity_splitDLNuE, 999, 999, 999, 999, (base_path + "slicePurity_all_weighted_splitDLNuE.pdf").c_str(), "bottomRight", nullptr, &right, true);
    styleDrawSplit(sliceCRUMBSScore_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceCRUMBSScore_all_weighted_splitDLNuE.pdf").c_str(), "topLeft", nullptr, &right, true);
    styleDrawSplit(sliceNumPFPs_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPFPs_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumPrimaryPFPs_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(sliceNumNeutrinos_splitDLNuE, 999, 999, 999, 999, (base_path + "sliceNumNeutrinos_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(QSquaredHighest_splitDLNuE, 999, 999, 999, 999, (base_path + "QSquared_highest_all_lower_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(QSquaredSum_splitDLNuE, 999, 999, 999, 999, (base_path + "QSquared_sum_all_lower_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(ERecoSumThetaReco_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(ERecoHighestThetaReco_splitDLNuE, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreHighestEnergyPFP_splitDLNuE, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreAllPFPs_splitDLNuE, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(trackscoreHighestScorePFPs_splitDLNuE, 999, 999, 999, 999, (base_path + "trackscoreHighestScorePFPs_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(recoX_low_splitDLNuE, 999, 999, 999, 999, (base_path + "recoX_low_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(recoX_high_splitDLNuE, 999, 999, 999, 999, (base_path + "recoX_high_all_weighted_splitDLNuE.pdf").c_str(), "topLeft", nullptr, &right, true);
    styleDrawSplit(recoY_low_splitDLNuE, 999, 999, 999, 999, (base_path + "recoY_low_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(recoY_high_splitDLNuE, 999, 999, 999, 999, (base_path + "recoY_high_all_weighted_splitDLNuE.pdf").c_str(), "topLeft", nullptr, &right, true);
    styleDrawSplit(recoZ_low_splitDLNuE, 999, 999, 999, 999, (base_path + "recoZ_low_all_weighted_splitDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);
    styleDrawSplit(recoZ_high_splitDLNuE, 999, 999, 999, 999, (base_path + "recoZ_high_all_weighted_splitDLNuE.pdf").c_str(), "topLeft", nullptr, &right, true);
    
    styleDrawPFPSplit(trackscoreAllPFPs_splitPFPDLNuE, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_all_weighted_splitPFPDLNuE.pdf").c_str(), "topRight", nullptr, &right, true);

    printf("Number of Events\nUnweighted BDT: Cosmic = %f, BNB = %f, Nu+E = %f\n", numEvents_BDTCosmic, numEvents_BDTBNB, numEvents_BDTNuE);
    printf("Unweighted DL Nu+E: Cosmic = %f, BNB = %f, Nu+E = %f\n", numEvents_DLNuECosmic, numEvents_DLNuEBNB, numEvents_DLNuENuE);
    printf("Weighted BDT: Cosmic = %f, BNB = %f, Nu+E = %f\n", (numEvents_BDTCosmic * weights.cosmicsCurrent), (numEvents_BDTBNB * weights.BNBCurrent), (numEvents_BDTNuE * weights.signalCurrent));
    printf("Weighted DL Nu+E: Cosmic = %f, BNB = %f, Nu+E = %f\n", (numEvents_DLNuECosmic * weights.cosmicsNuE), (numEvents_DLNuEBNB * weights.BNBNuE), (numEvents_DLNuENuE * weights.signalNuE));
    //double totalEvent_BDT = ((numEvents_BDTCosmic * weights.cosmicsCurrent) + (numEvents_BDTBNB * weights.BNBCurrent) + (numEvents_BDTNuE * weights.signalCurrent));
    double totalEvent_BDT = ((numEvents_BDTCosmic * weights.cosmicsCurrent) + (BNBSpillsSumCurrent * weights.BNBCurrent) + (NuESpillsSumCurrent * weights.signalCurrent));
    double cosmicPerc_BDT = ((100 * numEvents_BDTCosmic * weights.cosmicsCurrent)/totalEvent_BDT);
    //double BNBPerc_BDT = ((100 * numEvents_BDTBNB * weights.BNBCurrent)/totalEvent_BDT);
    //double NuEPerc_BDT = ((100 * numEvents_BDTNuE * weights.signalCurrent)/totalEvent_BDT);
    double BNBPerc_BDT = ((100 * BNBSpillsSumCurrent * weights.BNBCurrent)/totalEvent_BDT);
    double NuEPerc_BDT = ((100 * NuESpillsSumCurrent * weights.signalCurrent)/totalEvent_BDT);
    printf("Event Rates:\nBDT: Cosmic = %f, BNB = %f, Nu+E = %f\n", cosmicPerc_BDT, BNBPerc_BDT, NuEPerc_BDT);
    std::cout << "Cosmics = " << numEvents_BDTCosmic << ", BNB = " << BNBSpillsSumCurrent << ", NUE = " << NuESpillsSumCurrent << std::endl;

    //printf("\n______ Number of Events Left After Cuts ______\nBDT\nSignal: Before = %f, After = %f (%f %% left)\nSignal Fuzzy: Before = %f, After = %f (%f %% left)\nBNB: Before = %f, After = %f (%f %% left)\nBNB Fuzzy: Before = %f, After = %f (%f %% left)\nCosmics: Before = %f, After = %f (%f %% left)\n", numSignal_beforeCut_BDT, numSignal_afterCut_BDT, (100*numSignal_beforeCut_BDT/numSignal_afterCut_BDT), numSignalFuzzy_beforeCut_BDT, numSignalFuzzy_afterCut_BDT, (100*numSignalFuzzy_beforeCut_BDT/numSignalFuzzy_afterCut_BDT), numBNB_beforeCut_BDT, numBNB_afterCut_BDT, (100*numBNB_beforeCut_BDT/numBNB_afterCut_BDT), numBNBFuzzy_beforeCut_BDT, numBNBFuzzy_afterCut_BDT, (100*numBNBFuzzy_beforeCut_BDT/numBNBFuzzy_afterCut_BDT), numCosmic_beforeCut_BDT, numCosmic_afterCut_BDT, (100*numCosmic_beforeCut_BDT/numCosmic_afterCut_BDT));
    //printf("\nDL Uboone\nSignal: Before = %f, After = %f (%f %% left)\nSignal Fuzzy: Before = %f, After = %f (%f %% left)\nBNB: Before = %f, After = %f (%f %% left)\nBNB Fuzzy: Before = %f, After = %f (%f %% left)\nCosmics: Before = %f, After = %f (%f %% left)\n", numSignal_beforeCut_DLUboone, numSignal_afterCut_DLUboone, (100*numSignal_beforeCut_DLUboone/numSignal_afterCut_DLUboone), numSignalFuzzy_beforeCut_DLUboone, numSignalFuzzy_afterCut_DLUboone, (100*numSignalFuzzy_beforeCut_DLUboone/numSignalFuzzy_afterCut_DLUboone), numBNB_beforeCut_DLUboone, numBNB_afterCut_DLUboone, (100*numBNB_beforeCut_DLUboone/numBNB_afterCut_DLUboone), numBNBFuzzy_beforeCut_DLUboone, numBNBFuzzy_afterCut_DLUboone, (100*numBNBFuzzy_beforeCut_DLUboone/numBNBFuzzy_afterCut_DLUboone), numCosmic_beforeCut_DLUboone, numCosmic_afterCut_DLUboone, (100*numCosmic_beforeCut_DLUboone/numCosmic_afterCut_DLUboone));
    //printf("\nDL Nu+E\nSignal: Before = %f, After = %f (%f %% left)\nSignal Fuzzy: Before = %f, After = %f (%f %% left)\nBNB: Before = %f, After = %f (%f %% left)\nBNB Fuzzy: Before = %f, After = %f (%f %% left)\nCosmics: Before = %f, After = %f (%f %% left)\n", numSignal_beforeCut_DLNuE, numSignal_afterCut_DLNuE, (100*numSignal_beforeCut_DLNuE/numSignal_afterCut_DLNuE), numSignalFuzzy_beforeCut_DLNuE, numSignalFuzzy_afterCut_DLNuE, (100*numSignalFuzzy_beforeCut_DLNuE/numSignalFuzzy_afterCut_DLNuE), numBNB_beforeCut_DLNuE, numBNB_afterCut_DLNuE, (100*numBNB_beforeCut_DLNuE/numBNB_afterCut_DLNuE), numBNBFuzzy_beforeCut_DLNuE, numBNBFuzzy_afterCut_DLNuE, (100*numBNBFuzzy_beforeCut_DLNuE/numBNBFuzzy_afterCut_DLNuE), numCosmic_beforeCut_DLNuE, numCosmic_afterCut_DLNuE, (100*numCosmic_beforeCut_DLNuE/numCosmic_afterCut_DLNuE));

    printf("\n______ Number of Events Left After Cuts ______\nBDT\nSignal: Before = %f, After = %f (%f %% left)\nSignal Fuzzy: Before = %f, After = %f (%f %% left)\nBNB: Before = %f, After = %f (%f %% left)\nBNB Fuzzy: Before = %f, After = %f (%f %% left)\nCosmics: Before = %f, After = %f (%f %% left)\n", numSignal_beforeCut_BDT, numSignal_afterCut_BDT, (100 * numSignal_afterCut_BDT / numSignal_beforeCut_BDT), numSignalFuzzy_beforeCut_BDT, numSignalFuzzy_afterCut_BDT, (100 * numSignalFuzzy_afterCut_BDT / numSignalFuzzy_beforeCut_BDT), numBNB_beforeCut_BDT, numBNB_afterCut_BDT, (100 * numBNB_afterCut_BDT / numBNB_beforeCut_BDT), numBNBFuzzy_beforeCut_BDT, numBNBFuzzy_afterCut_BDT, (100 * numBNBFuzzy_afterCut_BDT / numBNBFuzzy_beforeCut_BDT), numCosmic_beforeCut_BDT, numCosmic_afterCut_BDT, (100 * numCosmic_afterCut_BDT / numCosmic_beforeCut_BDT));
    printf("\nDL Uboone\nSignal: Before = %f, After = %f (%f %% left)\nSignal Fuzzy: Before = %f, After = %f (%f %% left)\nBNB: Before = %f, After = %f (%f %% left)\nBNB Fuzzy: Before = %f, After = %f (%f %% left)\nCosmics: Before = %f, After = %f (%f %% left)\n", numSignal_beforeCut_DLUboone, numSignal_afterCut_DLUboone, (100 * numSignal_afterCut_DLUboone / numSignal_beforeCut_DLUboone), numSignalFuzzy_beforeCut_DLUboone, numSignalFuzzy_afterCut_DLUboone, (100 * numSignalFuzzy_afterCut_DLUboone / numSignalFuzzy_beforeCut_DLUboone), numBNB_beforeCut_DLUboone, numBNB_afterCut_DLUboone, (100 * numBNB_afterCut_DLUboone / numBNB_beforeCut_DLUboone), numBNBFuzzy_beforeCut_DLUboone, numBNBFuzzy_afterCut_DLUboone, (100 * numBNBFuzzy_afterCut_DLUboone / numBNBFuzzy_beforeCut_DLUboone), numCosmic_beforeCut_DLUboone, numCosmic_afterCut_DLUboone, (100 * numCosmic_afterCut_DLUboone / numCosmic_beforeCut_DLUboone));
    printf("\nDL Nu+E\nSignal: Before = %f, After = %f (%f %% left)\nSignal Fuzzy: Before = %f, After = %f (%f %% left)\nBNB: Before = %f, After = %f (%f %% left)\nBNB Fuzzy: Before = %f, After = %f (%f %% left)\nCosmics: Before = %f, After = %f (%f %% left)\n", numSignal_beforeCut_DLNuE, numSignal_afterCut_DLNuE, (100 * numSignal_afterCut_DLNuE / numSignal_beforeCut_DLNuE), numSignalFuzzy_beforeCut_DLNuE, numSignalFuzzy_afterCut_DLNuE, (100 * numSignalFuzzy_afterCut_DLNuE / numSignalFuzzy_beforeCut_DLNuE), numBNB_beforeCut_DLNuE, numBNB_afterCut_DLNuE, (100 * numBNB_afterCut_DLNuE / numBNB_beforeCut_DLNuE), numBNBFuzzy_beforeCut_DLNuE, numBNBFuzzy_afterCut_DLNuE, (100 * numBNBFuzzy_afterCut_DLNuE / numBNBFuzzy_beforeCut_DLNuE), numCosmic_beforeCut_DLNuE, numCosmic_afterCut_DLNuE, (100 * numCosmic_afterCut_DLNuE / numCosmic_beforeCut_DLNuE));

    double totalBackgroundEvents_before_DLNuE = (numSignalFuzzy_beforeCut_DLNuE + numBNB_beforeCut_DLNuE + numBNBFuzzy_beforeCut_DLNuE + numCosmic_beforeCut_DLNuE);
    double totalBackgroundEvents_after_DLNuE = (numSignalFuzzy_beforeCut_DLNuE + numBNB_afterCut_DLNuE + numBNBFuzzy_afterCut_DLNuE + numCosmic_afterCut_DLNuE);
    double totalBackgroundEvents_before_BDT = (numSignalFuzzy_beforeCut_BDT + numBNB_beforeCut_BDT + numBNBFuzzy_beforeCut_BDT + numCosmic_beforeCut_BDT);
    double totalBackgroundEvents_after_BDT = (numSignalFuzzy_beforeCut_BDT + numBNB_afterCut_BDT + numBNBFuzzy_afterCut_BDT + numCosmic_afterCut_BDT);
    double totalBackgroundEvents_before_DLUboone = (numSignalFuzzy_beforeCut_DLUboone + numBNB_beforeCut_DLUboone + numBNBFuzzy_beforeCut_DLUboone + numCosmic_beforeCut_DLUboone);
    double totalBackgroundEvents_after_DLUboone = (numSignalFuzzy_beforeCut_DLUboone + numBNB_afterCut_DLUboone + numBNBFuzzy_afterCut_DLUboone + numCosmic_afterCut_DLUboone);

    printf("\n\nDL Nu+E:\nNumber of Background Events Left = %f (%f%%)\nNumber of Signal Events Left = %f (%f%%)\n", totalBackgroundEvents_after_DLNuE, (100*totalBackgroundEvents_after_DLNuE/totalBackgroundEvents_before_DLNuE), numSignal_afterCut_DLNuE, (100 * numSignal_afterCut_DLNuE / numSignal_beforeCut_DLNuE));
    printf("\nDL Uboone:\nNumber of Background Events Left = %f (%f%%)\nNumber of Signal Events Left = %f (%f%%)\n", totalBackgroundEvents_after_DLUboone, (100*totalBackgroundEvents_after_DLUboone/totalBackgroundEvents_before_DLUboone), numSignal_afterCut_DLUboone, (100 * numSignal_afterCut_DLUboone / numSignal_beforeCut_DLUboone));
    printf("\nBDT:\nNumber of Background Events Left = %f (%f%%)\nNumber of Signal Events Left = %f (%f%%)\n", totalBackgroundEvents_after_BDT, (100*totalBackgroundEvents_after_BDT/totalBackgroundEvents_before_BDT), numSignal_afterCut_BDT, (100 * numSignal_afterCut_BDT / numSignal_beforeCut_BDT));

    double signalEff_BDT = (numSignal_afterCut_BDT / numSignal_beforeCut_BDT);
    double signalPur_BDT = (numSignal_afterCut_BDT / (numSignal_afterCut_BDT + numSignalFuzzy_afterCut_BDT + numBNB_afterCut_BDT + numBNBFuzzy_afterCut_BDT + numCosmic_afterCut_BDT));

    double signalEff_DLUboone = (numSignal_afterCut_DLUboone / numSignal_beforeCut_DLUboone);
    double signalPur_DLUboone = (numSignal_afterCut_DLUboone / (numSignal_afterCut_DLUboone + numSignalFuzzy_afterCut_DLUboone + numBNB_afterCut_DLUboone + numBNBFuzzy_afterCut_DLUboone + numCosmic_afterCut_DLUboone));

    double signalEff_DLNuE = (numSignal_afterCut_DLNuE / numSignal_beforeCut_DLNuE);
    double signalPur_DLNuE = (numSignal_afterCut_DLNuE / (numSignal_afterCut_DLNuE + numSignalFuzzy_afterCut_DLNuE + numBNB_afterCut_DLNuE + numBNBFuzzy_afterCut_DLNuE + numCosmic_afterCut_DLNuE));

    printf("\nBDT: Signal Eff = %f, Signal Pur = %f, Signal Eff x Pur = %f\n", signalEff_BDT, signalPur_BDT, signalEff_BDT*signalPur_BDT);
    printf("DLUboone: Signal Eff = %f, Signal Pur = %f, Signal Eff x Pur = %f\n", signalEff_DLUboone, signalPur_DLUboone, signalEff_DLUboone*signalPur_DLUboone);
    printf("DLNuE: Signal Eff = %f, Signal Pur = %f, Signal Eff x Pur = %f\n", signalEff_DLNuE, signalPur_DLNuE, signalEff_DLNuE*signalPur_DLNuE);

    printf("\n\n============ With Weighting ============\n");
    printf("Before Cuts:\nnu+e = %f, NCNpi0 = %f, Other NC = %f, CCnumu = %f\nCCnue = %f, Dirt = %f, nu+e Dirt = %f, Cosmic = %f, Other = %f\n\n", numEventBeforeCutDLNuE.nuE, numEventBeforeCutDLNuE.NCNPi0, numEventBeforeCutDLNuE.otherNC, numEventBeforeCutDLNuE.CCnumu, numEventBeforeCutDLNuE.CCnue, numEventBeforeCutDLNuE.dirt, numEventBeforeCutDLNuE.nuEDirt, numEventBeforeCutDLNuE.cosmic, numEventBeforeCutDLNuE.other);
    printf("After Cut:\nnu+e = %f (%f%%), NCNpi0 = %f (%f%%), Other NC = %f (%f%%), CCnumu = %f (%f%%)\n", numEventCutDLNuE.nuE, (100.0*numEventCutDLNuE.nuE/numEventBeforeCutDLNuE.nuE), numEventCutDLNuE.NCNPi0, (100.0*numEventCutDLNuE.NCNPi0/numEventBeforeCutDLNuE.NCNPi0), numEventCutDLNuE.otherNC, (100.0*numEventCutDLNuE.otherNC/numEventBeforeCutDLNuE.otherNC), numEventCutDLNuE.CCnumu, (100.0*numEventCutDLNuE.CCnumu/numEventBeforeCutDLNuE.CCnumu));
    printf("CCnue = %f (%f%%), Dirt = %f (%f%%), nu+e Dirt = %f (%f%%), Cosmic = %f (%f%%), Other = %f (%f%%)\n", numEventCutDLNuE.CCnue, (100.0*numEventCutDLNuE.CCnue/numEventBeforeCutDLNuE.CCnue), numEventCutDLNuE.dirt, (100.0*numEventCutDLNuE.dirt/numEventBeforeCutDLNuE.dirt), numEventCutDLNuE.nuEDirt, (100.0*numEventCutDLNuE.nuEDirt/numEventBeforeCutDLNuE.nuEDirt), numEventCutDLNuE.cosmic, (100.0*numEventCutDLNuE.cosmic/numEventBeforeCutDLNuE.cosmic), numEventCutDLNuE.other, (100.0*numEventCutDLNuE.other/numEventBeforeCutDLNuE.other));

    printf("\n\n============ Without Weighting ============\n");
    printf("Before Cuts:\nnu+e = %f, NCNpi0 = %f, Other NC = %f, CCnumu = %f\nCCnue = %f, Dirt = %f, nu+e Dirt = %f, Cosmic = %f, Other = %f\n\n", numEventBeforeCutWithoutWeightingDLNuE.nuE, numEventBeforeCutWithoutWeightingDLNuE.NCNPi0, numEventBeforeCutWithoutWeightingDLNuE.otherNC, numEventBeforeCutWithoutWeightingDLNuE.CCnumu, numEventBeforeCutWithoutWeightingDLNuE.CCnue, numEventBeforeCutWithoutWeightingDLNuE.dirt, numEventBeforeCutWithoutWeightingDLNuE.nuEDirt, numEventBeforeCutWithoutWeightingDLNuE.cosmic, numEventBeforeCutWithoutWeightingDLNuE.other);
    printf("After Cut:\nnu+e = %f (%f%%), NCNpi0 = %f (%f%%), Other NC = %f (%f%%), CCnumu = %f (%f%%)\n", numEventCutWithoutWeightingDLNuE.nuE, (100.0*numEventCutWithoutWeightingDLNuE.nuE/numEventBeforeCutWithoutWeightingDLNuE.nuE), numEventCutWithoutWeightingDLNuE.NCNPi0, (100.0*numEventCutWithoutWeightingDLNuE.NCNPi0/numEventBeforeCutWithoutWeightingDLNuE.NCNPi0), numEventCutWithoutWeightingDLNuE.otherNC, (100.0*numEventCutWithoutWeightingDLNuE.otherNC/numEventBeforeCutWithoutWeightingDLNuE.otherNC), numEventCutWithoutWeightingDLNuE.CCnumu, (100.0*numEventCutWithoutWeightingDLNuE.CCnumu/numEventBeforeCutWithoutWeightingDLNuE.CCnumu));
    printf("CCnue = %f (%f%%), Dirt = %f (%f%%), nu+e Dirt = %f (%f%%), Cosmic = %f (%f%%), Other = %f (%f%%)\n", numEventCutWithoutWeightingDLNuE.CCnue, (100.0*numEventCutWithoutWeightingDLNuE.CCnue/numEventBeforeCutWithoutWeightingDLNuE.CCnue), numEventCutWithoutWeightingDLNuE.dirt, (100.0*numEventCutWithoutWeightingDLNuE.dirt/numEventBeforeCutWithoutWeightingDLNuE.dirt), numEventCutWithoutWeightingDLNuE.nuEDirt, (100.0*numEventCutWithoutWeightingDLNuE.nuEDirt/numEventBeforeCutWithoutWeightingDLNuE.nuEDirt), numEventCutWithoutWeightingDLNuE.cosmic, (100.0*numEventCutWithoutWeightingDLNuE.cosmic/numEventBeforeCutWithoutWeightingDLNuE.cosmic), numEventCutWithoutWeightingDLNuE.other, (100.0*numEventCutWithoutWeightingDLNuE.other/numEventBeforeCutWithoutWeightingDLNuE.other));

    std::cout << "" << std::endl;
    std::cout << "clear cosmic = " << clearCosmicCut << ", num PFPs 0 = " << numPFPs0Cut << ", num reco neutrinos 0 = " << numRecoNeutrinosCut << ", FV = " << FVCut << std::endl;
    std::cout << "CRUMBS = " << CRUMBSCut << ", Primary PFP = " << primaryPFPCut << ", trackscore of all PFPs = " << trackscorePFPsCut << ", trackscore of highest energy PFP = " << trackscoreCut << std::endl;
    std::cout << "ETheta2 (highest energy PFP) = " << ETheta2Cut << std::endl; 

    std::ofstream out_file(txtFileName, std::ios::app);
    if(out_file.is_open()){
        out_file << "================" << std::endl;
        out_file << "clear cosmic = " << clearCosmicCut << ", num PFPs 0 = " << numPFPs0Cut << ", num reco neutrinos 0 = " << numRecoNeutrinosCut << ", FV = " << FVCut << std::endl;
        out_file << "CRUMBS = " << CRUMBSCut << ", Primary PFP = " << primaryPFPCut << ", trackscore of all PFPs = " << trackscorePFPsCut << ", trackscore of highest energy PFP = " << trackscoreCut << std::endl;
        out_file << "ETheta2 (highest energy PFP) = " << ETheta2Cut << std::endl; 
        out_file << "================" << std::endl;
        out_file.close(); 
    } else{
        std::cerr << "Error: couldnt open txt file" << std::endl;
    }

    printf("\nNumber of: nu+e electrons = %f, nu+e other = %f, electron = %f, proton = %f, muon = %f\npi0 = %f, charged pi = %f, other = %f, cosmic muon = %f, cosmic other = %f\n", numPFPsBeforeDLNuE.nuEElectron, numPFPsBeforeDLNuE.nuEOther, numPFPsBeforeDLNuE.electron, numPFPsBeforeDLNuE.proton, numPFPsBeforeDLNuE.muon, numPFPsBeforeDLNuE.pi0, numPFPsBeforeDLNuE.chargedPi, numPFPsBeforeDLNuE.other, numPFPsBeforeDLNuE.cosmicMuon, numPFPsBeforeDLNuE.cosmicOther);

    std::cout << "" << std::endl;
    std::cout << "Number of slices where the highest energy PFP is primary (after cuts) = " << highestEnergyPFPPrimary_DLNuE << " (" << 100*highestEnergyPFPPrimary_DLNuE/totalLeft_DLNuE << "%)" << std::endl;
    std::cout << "Slices left = " << totalLeft_DLNuE << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Number of nu+e slices where the highest energy PFP is primary (after cuts) = " << highestEnergyPFPPrimaryNuE_DLNuE << " (" << 100*highestEnergyPFPPrimaryNuE_DLNuE/totalLeftNuE_DLNuE << "%)" << std::endl;
    std::cout << "nu+e slices left = " << totalLeftNuE_DLNuE << std::endl;
    
    //printf("\nNum Nu+E from slice category = %f, num Nu+E from int type = %f\n", numNuESliceCategory_DLNuE, numNuEIntType_DLNuE);
    //printf("Num nu+e with completeness > 0.5 = %f, num nu+e in else = %f\n", numNuESliceCategoryPassed_DLNuE, numNuESliceCategoryElse_DLNuE);
}
