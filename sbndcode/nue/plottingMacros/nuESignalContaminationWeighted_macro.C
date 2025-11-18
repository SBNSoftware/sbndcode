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
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* currentSignal_signal;
    TH1F* nuESignal_signal;
    TH1F* currentSignal_signalfuzzy;
    TH1F* nuESignal_signalfuzzy;
    TH1F* currentSignalCosmics_signal;
    TH1F* nuESignalCosmics_signal;
    TH1F* currentSignalCosmics_signalfuzzy;
    TH1F* nuESignalCosmics_signalfuzzy;
    TH1F* currentSignalCosmics_cosmic;
    TH1F* nuESignalCosmics_cosmic;
} histGroup_struct;

struct weights_struct{
    double signalCurrent = 0;
    double signalNuE = 0;
    double signalCosmicsCurrent = 0;
    double signalCosmicsNuE = 0;
};

histGroup_struct createHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);
    
    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_currentSignal_signal").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignal_signal").c_str()),
        (TH1F*) base->Clone((baseName + "_currentSignal_signalfuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignal_signalfuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_currentSignalCosmics_signal").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignalCosmics_signal").c_str()),
        (TH1F*) base->Clone((baseName + "_currentSignalCosmics_signalfuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignalCosmics_signalfuzzy").c_str()),
        (TH1F*) base->Clone((baseName + "_currentSignalCosmics_cosmic").c_str()),
        (TH1F*) base->Clone((baseName + "_nuESignalCosmics_cosmic").c_str()),
    };
}

void styleDrawAll(histGroup_struct hists,
                  double ymin, double ymax, double xmin, double xmax,
                  const char* filename, const std::string& legendLocation,
                  bool includeSignal = true, bool includeSignalFuzzy = true,
                  bool includeCosmic = true, bool includeBDT = true, bool includeDLNuE = true){
    
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    std::vector<TH1F*> allHists;

    if(includeCosmic == true){
        allHists = {
            hists.currentSignal_signal, hists.nuESignal_signal, hists.currentSignal_signalfuzzy, hists.nuESignal_signalfuzzy, 
            hists.currentSignalCosmics_signal, hists.nuESignalCosmics_signal, hists.currentSignalCosmics_signalfuzzy, hists.nuESignalCosmics_signalfuzzy,
            hists.currentSignalCosmics_cosmic, hists.nuESignalCosmics_cosmic
        };
    } else{
        allHists = {
            hists.currentSignal_signal, hists.nuESignal_signal, hists.currentSignal_signalfuzzy, hists.nuESignal_signalfuzzy, 
            hists.currentSignalCosmics_signal, hists.nuESignalCosmics_signal, hists.currentSignalCosmics_signalfuzzy, hists.nuESignalCosmics_signalfuzzy
        };
    }

    for (auto* hist : allHists)
        if (hist) hist->SetStats(0);

    hists.currentSignal_signal->SetLineWidth(2);                hists.currentSignal_signal->SetLineColor(kViolet+1);
    hists.currentSignal_signalfuzzy->SetLineWidth(2);           hists.currentSignal_signalfuzzy->SetLineColor(kOrange+7);
    
    hists.nuESignal_signal->SetLineWidth(2);                    hists.nuESignal_signal->SetLineColor(kBlue+1);
    hists.nuESignal_signalfuzzy->SetLineWidth(2);               hists.nuESignal_signalfuzzy->SetLineColor(kGreen+3);

    hists.currentSignalCosmics_signal->SetLineWidth(2);         hists.currentSignalCosmics_signal->SetLineColor(kViolet+6);
    hists.currentSignalCosmics_signalfuzzy->SetLineWidth(2);    hists.currentSignalCosmics_signalfuzzy->SetLineColor(kOrange+6);

    hists.nuESignalCosmics_signal->SetLineWidth(2);             hists.nuESignalCosmics_signal->SetLineColor(kBlue-7);
    hists.nuESignalCosmics_signalfuzzy->SetLineWidth(2);        hists.nuESignalCosmics_signalfuzzy->SetLineColor(kGreen+1);

    hists.currentSignalCosmics_cosmic->SetLineWidth(2);         hists.currentSignalCosmics_cosmic->SetLineColor(kPink+9);
    hists.nuESignalCosmics_cosmic->SetLineWidth(2);             hists.nuESignalCosmics_cosmic->SetLineColor(kRed);

    if((ymin != 999) && (ymax != 999)){
        for(auto* hist : allHists)
            if(hist) hist->GetYaxis()->SetRangeUser(ymin, ymax);
    }

    if((xmin != 999) && (xmax != 999)){
        for(auto* hist : allHists)
            if(hist) hist->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    double maxYValue = 0.0;
    for(auto* hist : allHists)
        if(hist && hist->GetMaximum() > maxYValue)
            maxYValue = hist->GetMaximum();

    std::cout << "maxYValue = " << maxYValue << std::endl;
    double yminVal = 0;
    if((ymin == 999) && (ymax == 999)){
        double ymaxVal = maxYValue * 1.1;
        std::cout << "setting yaxis to " << yminVal << ", " << ymaxVal << std::endl;
        for(auto* hist : allHists)
            if(hist) hist->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    bool first = true;
    auto draw = [&](TH1* hist){ if (hist) { hist->Draw(first ? "hist" : "histsame"); first = false; } };

    // Decide whether a histogram name is allowed by the boolean flags
    auto variantAllowed = [&](const std::string& name) {
        bool isSignal        = name.size() >= 7  && name.rfind("_signal") == name.size() - 7;
        bool isSignalFuzzy   = name.size() >= 12 && name.rfind("_signalfuzzy") == name.size() - 12;
        bool isCosmic        = name.size() >= 7  && name.rfind("_cosmic") == name.size() - 7;

        bool isBDT  = name.find("current") == 0;   // starts with current
        bool isDLNuE = name.find("nuE") == 0;      // starts with nuE

        if (!includeSignal        && isSignal)        return false;
        if (!includeSignalFuzzy   && isSignalFuzzy)   return false;
        if (!includeCosmic        && isCosmic)        return false;
        if (!includeBDT           && isBDT)           return false;
        if (!includeDLNuE         && isDLNuE)         return false;

        return true;
    };

    if(variantAllowed("currentSignal_signal")) draw(hists.currentSignal_signal);
    if(variantAllowed("nuESignal_signal")) draw(hists.nuESignal_signal);
    if(variantAllowed("currentSignal_signalfuzzy")) draw(hists.currentSignal_signalfuzzy);
    if(variantAllowed("nuESignal_signalfuzzy")) draw(hists.nuESignal_signalfuzzy);
    if(variantAllowed("currentSignalCosmics_signal")) draw(hists.currentSignalCosmics_signal);
    if(variantAllowed("nuESignalCosmics_signal")) draw(hists.nuESignalCosmics_signal);
    if(variantAllowed("currentSignalCosmics_signalfuzzy")) draw(hists.currentSignalCosmics_signalfuzzy);
    if(variantAllowed("nuESignalCosmics_signalfuzzy")) draw(hists.nuESignalCosmics_signalfuzzy);
    if(variantAllowed("currentSignalCosmics_cosmic")) draw(hists.currentSignalCosmics_cosmic);
    if(variantAllowed("nuESignalCosmics_cosmic")) draw(hists.nuESignalCosmics_cosmic);

    for(auto* hist : allHists){
        if(hist){
            hist->SetStats(0);
            hist->GetXaxis()->SetTickLength(0.04);
            hist->GetYaxis()->SetTickLength(0.03);
            hist->GetXaxis()->SetTickSize(0.02);
            hist->GetYaxis()->SetTickSize(0.02);
        }
    }

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0;
    std::vector<std::pair<TH1*, std::string>> legendEntries;

    auto addLegendIf = [&](TH1* hist, const std::string& label, const std::string& name){
        if (hist && variantAllowed(name)) legendEntries.emplace_back(hist, label);
    }; 

    addLegendIf(hists.currentSignal_signal, "Signal, Without Cosmics, BDT Vertexing", "currentSignal_signal");
    addLegendIf(hists.nuESignal_signal, "Signal, Without Cosmics, DL Nu+E Vertexing", "nuESignal_signal");
    addLegendIf(hists.currentSignal_signalfuzzy, "Signal Fuzzy, Without Cosmics, BDT Vertexing", "currentSignal_signalfuzzy");
    addLegendIf(hists.nuESignal_signalfuzzy, "Signal Fuzzy, Without Cosmics, DL Nu+E Vertexing", "nuESignal_signalfuzzy");
    addLegendIf(hists.currentSignalCosmics_signal, "Signal, With Cosmics, BDT Vertexing", "currentSignalCosmics_signal");
    addLegendIf(hists.nuESignalCosmics_signal, "Signal, With Cosmics, DL Nu+E Vertexing", "nuESignalCosmics_signal");
    addLegendIf(hists.currentSignalCosmics_signalfuzzy, "Signal Fuzzy, With Cosmics, BDT Vertexing", "currentSignalCosmics_signalfuzzy");
    addLegendIf(hists.nuESignalCosmics_signalfuzzy, "Signal Fuzzy, With Cosmics, DL Nu+E Vertexing", "nuESignalCosmics_signalfuzzy");
    addLegendIf(hists.currentSignalCosmics_cosmic, "Cosmic, With Cosmics, BDT Vertexing", "currentSignalCosmics_cosmic");
    addLegendIf(hists.nuESignalCosmics_cosmic, "Cosmic, With Cosmics, DL Nu+E Vertexing", "nuESignalCosmics_cosmic");

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

    hists.canvas->SaveAs(filename);
}

void nuESignalContaminationWeighted_macro(){
   
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_Nu+E_WithWithoutCosmics_new.root");
    std::string base_path = "/nashome/c/coackley/nuESignalWithWithoutCosmics/";

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
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalNuE;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCosmicsCurrent;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCosmicsNuE;

    double totalPOTSignalCurrent = 0;
    double totalPOTSignalNuE = 0;
    double POTSignalNuE_notMissing = 0;
    double POTSignalBDT_notMissing = 0;

    double totalPOTSignalCosmicsCurrent = 0;
    double totalPOTSignalCosmicsNuE = 0;
    double POTSignalCosmicsNuE_notMissing = 0;
    double POTSignalCosmicsBDT_notMissing = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        std::pair<unsigned int, unsigned int> key = std::make_pair(subRunRun, subRunNumber);

        if(subRunSignal == 1){
            if(subRunDLCurrent == 2 && seenSubRunsSignalCosmicsCurrent.find(key) == seenSubRunsSignalCosmicsCurrent.end()){
                totalPOTSignalCosmicsCurrent += subRunPOT;
                seenSubRunsSignalCosmicsCurrent.insert(key);             
            } else if(subRunDLCurrent == 5 && seenSubRunsSignalCosmicsNuE.find(key) == seenSubRunsSignalCosmicsNuE.end()){
                totalPOTSignalCosmicsNuE += subRunPOT;
                seenSubRunsSignalCosmicsNuE.insert(key);            
            }

            if(subRunDLCurrent == 2) POTSignalCosmicsBDT_notMissing += subRunPOT;
            if(subRunDLCurrent == 5) POTSignalCosmicsNuE_notMissing += subRunPOT;
        }
        
        if(subRunSignal == 0){
            if(subRunDLCurrent == 2 && seenSubRunsSignalCurrent.find(key) == seenSubRunsSignalCurrent.end()){
                totalPOTSignalCurrent += subRunPOT;
                seenSubRunsSignalCurrent.insert(key);             
            } else if(subRunDLCurrent == 5 && seenSubRunsSignalNuE.find(key) == seenSubRunsSignalNuE.end()){
                totalPOTSignalNuE += subRunPOT;
                seenSubRunsSignalNuE.insert(key);            
            }

            if(subRunDLCurrent == 2) POTSignalBDT_notMissing += subRunPOT;
            if(subRunDLCurrent == 5) POTSignalNuE_notMissing += subRunPOT;
        }
    }

    weights_struct weights;
    weights.signalCurrent = POTSignalBDT_notMissing/POTSignalBDT_notMissing;
    weights.signalNuE = POTSignalBDT_notMissing/POTSignalNuE_notMissing;
    weights.signalCosmicsCurrent = POTSignalBDT_notMissing/POTSignalCosmicsBDT_notMissing;
    weights.signalCosmicsNuE = POTSignalBDT_notMissing/POTSignalCosmicsNuE_notMissing;

    printf("Signal POT: BDT = %f, DL Nu+E = %f\nSignal + Cosmics POT: BDT = %f, DL Nu+E = %f\n", POTSignalBDT_notMissing, POTSignalNuE_notMissing, POTSignalCosmicsBDT_notMissing, POTSignalCosmicsNuE_notMissing);
    printf("Signal Weights: BDT = %f, Signal DL Nu+E = %f\nSignal + Cosmics Weights: BDT = %f, DL Nu+E = %f\n", weights.signalCurrent, weights.signalNuE, weights.signalCosmicsCurrent, weights.signalCosmicsNuE);

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
    auto sliceNumPFPs = createHistGroup("sliceNumPFPs", "Number of PFPs in the Slice", "Number of PFPs", 15, 0, 15);
    auto sliceNumPFPsDist = createHistGroup("sliceNumPFPsDist", "Number of PFPs in the Slice (Not Weighted)", "Number of PFPs", 15, 0, 15); 
    auto trackscoreHighestEnergyPFP = createHistGroup("trackscoreHighestEnergyPFP", "Trackscore of the PFP in the Slice with the Highest Energy", "Trackscore", 20, 0, 1);
    auto trackscoreHighestEnergyPFPDist = createHistGroup("trackscoreHighestEnergyPFPDist", "Trackscore of the PFP in the Slice with the Highest Energy (Not Weighted)", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPs = createHistGroup("trackscoreAllPFPs", "Trackscore of All PFPs in the Slice", "Trackscore", 20, 0, 1);
    auto trackscoreAllPFPsDist = createHistGroup("trackscoreAllPFPsDist", "Trackscore of All PFPs in the Slice (Not Weighted)", "Trackscore", 20, 0, 1);
    
    auto deltaTheta = createHistGroup("deltaTheta", "Angle Between the True Electron and the Highest Energy PFP in the Slice", "#Delta#theta (degrees)", 90, 0, 180);
    auto deltaThetaDist = createHistGroup("deltaThetaDist", "Angle Between the True Electron and the Highest Energy PFP in the Slice (Not Weighted)", "#Delta#theta (degrees)", 90, 0, 180);
    auto ERecoSumThetaReco = createHistGroup("ERecoSumThetaReco", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaRecoDist = createHistGroup("ERecoSumThetaRecoDist", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice (Not Weighted)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaReco = createHistGroup("ERecoHighestThetaReco", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaRecoDist = createHistGroup("ERecoHighestThetaRecoDist", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice (Not Weighted)", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);

    double numEvents_signalBDT = 0;
    double numEvents_signalDLNuE = 0;
    double numEvents_signalCosmicsBDT = 0;
    double numEvents_signalCosmicsDLNuE = 0;
    
    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);
        if(DLCurrent == 2 && signal == 1) numEvents_signalCosmicsBDT++;
        else if(DLCurrent == 5 && signal == 1) numEvents_signalCosmicsDLNuE++;
        else if(DLCurrent == 2 && signal == 0) numEvents_signalBDT++;
        else if(DLCurrent == 5 && signal == 0) numEvents_signalDLNuE++;

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
                // There is a true recoil electron in the event
                if(signal == 0 && DLCurrent == 2) trueETheta2.currentSignal_signal->Fill(truth_recoilElectronETheta2->at(i), weights.signalCurrent); 
                else if(signal == 0 && DLCurrent == 5) trueETheta2.nuESignal_signal->Fill(truth_recoilElectronETheta2->at(i), weights.signalNuE); 
                else if(signal == 1 && DLCurrent == 2) trueETheta2.currentSignalCosmics_signal->Fill(truth_recoilElectronETheta2->at(i), weights.signalCosmicsCurrent); 
                else if(signal == 1 && DLCurrent == 5) trueETheta2.nuESignalCosmics_signal->Fill(truth_recoilElectronETheta2->at(i), weights.signalCosmicsNuE); 
            }
        }

        double weight = 0;
        if(signal == 0 && DLCurrent == 2) weight = weights.signalCurrent;
        else if(signal == 0 && DLCurrent == 5) weight = weights.signalNuE;
        else if(signal == 1 && DLCurrent == 2) weight = weights.signalCosmicsCurrent;
        else if(signal == 1 && DLCurrent == 5) weight = weights.signalCosmicsNuE;

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
                double highestEnergy_trackscore = -999999;

                // Loop through all the PFPs in the event
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
                            highestEnergy_trackscore = reco_particleTrackScore->at(pfp);
                        }
                    }
                }

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

                if(reco_sliceCategory->at(slice) == 0){
                    // This is a cosmic slice
                    if(signal == 1 && DLCurrent == 2){
                        sliceCompleteness.currentSignalCosmics_cosmic->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.currentSignalCosmics_cosmic->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.currentSignalCosmics_cosmic->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.currentSignalCosmics_cosmic->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.currentSignalCosmics_cosmic->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.currentSignalCosmics_cosmic->Fill(numPFPsSlice);   

                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.currentSignalCosmics_cosmic->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignalCosmics_cosmic->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.currentSignalCosmics_cosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.currentSignalCosmics_cosmic->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }

                            deltaTheta.currentSignalCosmics_cosmic->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignalCosmics_cosmic->Fill(angleDifference);
                            ERecoSumThetaReco.currentSignalCosmics_cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignalCosmics_cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignalCosmics_cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignalCosmics_cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        
                        }
                    } else if(signal == 1 && DLCurrent == 5){
                        sliceCompleteness.nuESignalCosmics_cosmic->Fill(reco_sliceCompleteness->at(slice), weight);
                        sliceCompletenessDist.nuESignalCosmics_cosmic->Fill(reco_sliceCompleteness->at(slice));
                        slicePurity.nuESignalCosmics_cosmic->Fill(reco_slicePurity->at(slice), weight);
                        slicePurityDist.nuESignalCosmics_cosmic->Fill(reco_slicePurity->at(slice));
                        sliceNumPFPs.nuESignalCosmics_cosmic->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.nuESignalCosmics_cosmic->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.nuESignalCosmics_cosmic->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignalCosmics_cosmic->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.nuESignalCosmics_cosmic->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.nuESignalCosmics_cosmic->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.nuESignalCosmics_cosmic->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignalCosmics_cosmic->Fill(angleDifference);
                            ERecoSumThetaReco.nuESignalCosmics_cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignalCosmics_cosmic->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignalCosmics_cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignalCosmics_cosmic->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } 
                } else if(reco_sliceCategory->at(slice) == 1){
                    // This is a signal slice
                    if(signal == 1 && DLCurrent == 2){
                        sliceCompleteness.currentSignalCosmics_signal->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.currentSignalCosmics_signal->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.currentSignalCosmics_signal->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.currentSignalCosmics_signal->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.currentSignalCosmics_signal->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.currentSignalCosmics_signal->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.currentSignalCosmics_signal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignalCosmics_signal->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.currentSignalCosmics_signal->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.currentSignalCosmics_signal->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.currentSignalCosmics_signal->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignalCosmics_signal->Fill(angleDifference);
                            ERecoSumThetaReco.currentSignalCosmics_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignalCosmics_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignalCosmics_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignalCosmics_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } else if(signal == 1 && DLCurrent == 5){
                        sliceCompleteness.nuESignalCosmics_signal->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.nuESignalCosmics_signal->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.nuESignalCosmics_signal->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.nuESignalCosmics_signal->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.nuESignalCosmics_signal->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.nuESignalCosmics_signal->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.nuESignalCosmics_signal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignalCosmics_signal->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.nuESignalCosmics_signal->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.nuESignalCosmics_signal->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.nuESignalCosmics_signal->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignalCosmics_signal->Fill(angleDifference);
                            ERecoSumThetaReco.nuESignalCosmics_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignalCosmics_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignalCosmics_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignalCosmics_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } else if(signal == 0 && DLCurrent == 2){
                        sliceCompleteness.currentSignal_signal->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.currentSignal_signal->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.currentSignal_signal->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.currentSignal_signal->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.currentSignal_signal->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.currentSignal_signal->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.currentSignal_signal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignal_signal->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.currentSignal_signal->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.currentSignal_signal->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.currentSignal_signal->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignal_signal->Fill(angleDifference);
                            ERecoSumThetaReco.currentSignal_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignal_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignal_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignal_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } else if(signal == 0 && DLCurrent == 5){
                        sliceCompleteness.nuESignal_signal->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.nuESignal_signal->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.nuESignal_signal->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.nuESignal_signal->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.nuESignal_signal->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.nuESignal_signal->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.nuESignal_signal->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignal_signal->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.nuESignal_signal->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.nuESignal_signal->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.nuESignal_signal->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignal_signal->Fill(angleDifference);
                            ERecoSumThetaReco.nuESignal_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignal_signal->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignal_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignal_signal->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } 
                } else if(reco_sliceCategory->at(slice) == 2){
                    if(signal == 1 && DLCurrent == 2){
                        sliceCompleteness.currentSignalCosmics_signalfuzzy->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.currentSignalCosmics_signalfuzzy->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.currentSignalCosmics_signalfuzzy->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.currentSignalCosmics_signalfuzzy->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.currentSignalCosmics_signalfuzzy->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.currentSignalCosmics_signalfuzzy->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.currentSignalCosmics_signalfuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignalCosmics_signalfuzzy->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.currentSignalCosmics_signalfuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.nuESignalCosmics_signalfuzzy->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.currentSignalCosmics_signalfuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignalCosmics_signalfuzzy->Fill(angleDifference);
                            ERecoSumThetaReco.currentSignalCosmics_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignalCosmics_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignalCosmics_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignalCosmics_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } else if(signal == 1 && DLCurrent == 5){
                        sliceCompleteness.nuESignalCosmics_signalfuzzy->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.nuESignalCosmics_signalfuzzy->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.nuESignalCosmics_signalfuzzy->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.nuESignalCosmics_signalfuzzy->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.nuESignalCosmics_signalfuzzy->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.nuESignalCosmics_signalfuzzy->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.nuESignalCosmics_signalfuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignalCosmics_signalfuzzy->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.nuESignalCosmics_signalfuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.nuESignalCosmics_signalfuzzy->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.nuESignalCosmics_signalfuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignalCosmics_signalfuzzy->Fill(angleDifference);
                            ERecoSumThetaReco.nuESignalCosmics_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignalCosmics_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignalCosmics_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignalCosmics_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } else if(signal == 0 && DLCurrent == 2){
                        sliceCompleteness.currentSignal_signalfuzzy->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.currentSignal_signalfuzzy->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.currentSignal_signalfuzzy->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.currentSignal_signalfuzzy->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.currentSignal_signalfuzzy->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.currentSignal_signalfuzzy->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.currentSignal_signalfuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.currentSignal_signalfuzzy->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.currentSignal_signalfuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.currentSignal_signalfuzzy->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.currentSignal_signalfuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.currentSignal_signalfuzzy->Fill(angleDifference);
                            ERecoSumThetaReco.currentSignal_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.currentSignal_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.currentSignal_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.currentSignal_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } else if(signal == 0 && DLCurrent == 5){
                        sliceCompleteness.nuESignal_signalfuzzy->Fill(reco_sliceCompleteness->at(slice), weight); 
                        sliceCompletenessDist.nuESignal_signalfuzzy->Fill(reco_sliceCompleteness->at(slice)); 
                        slicePurity.nuESignal_signalfuzzy->Fill(reco_slicePurity->at(slice), weight); 
                        slicePurityDist.nuESignal_signalfuzzy->Fill(reco_slicePurity->at(slice)); 
                        sliceNumPFPs.nuESignal_signalfuzzy->Fill(numPFPsSlice, weight);    
                        sliceNumPFPsDist.nuESignal_signalfuzzy->Fill(numPFPsSlice);    
                        
                        if(highestEnergy_PFPID != -999999){
                            trackscoreHighestEnergyPFP.nuESignal_signalfuzzy->Fill(highestEnergy_trackscore, weight);
                            trackscoreHighestEnergyPFPDist.nuESignal_signalfuzzy->Fill(highestEnergy_trackscore);
                        
                            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                                    if(reco_particleTrackScore->at(pfp) == -999999) std::cout << "Trackscore is -999999" << std::endl;
                                    else{
                                        trackscoreAllPFPs.nuESignal_signalfuzzy->Fill(reco_particleTrackScore->at(pfp), weight);
                                        trackscoreAllPFPsDist.nuESignal_signalfuzzy->Fill(reco_particleTrackScore->at(pfp));
                                    }
                                }
                            }
                            
                            deltaTheta.nuESignal_signalfuzzy->Fill(angleDifference, weight);
                            deltaThetaDist.nuESignal_signalfuzzy->Fill(angleDifference);
                            ERecoSumThetaReco.nuESignal_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoSumThetaRecoDist.nuESignal_signalfuzzy->Fill((summedEnergy * highestEnergy_theta * highestEnergy_theta));
                            ERecoHighestThetaReco.nuESignal_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta), weight);
                            ERecoHighestThetaRecoDist.nuESignal_signalfuzzy->Fill((highestEnergy_energy * highestEnergy_theta * highestEnergy_theta));
                        }
                    } 
                } 
            }
        }
    }

    styleDrawAll(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompleteness_weighted.pdf").c_str(), "topLeft", true, true, false, true, true);
    styleDrawAll(sliceCompletenessDist, 999, 999, 999, 999, (base_path + "sliceCompleteness_dist.pdf").c_str(), "topLeft", true, true, false, true, true);
    styleDrawAll(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity_weighted.pdf").c_str(), "topLeft", true, true, false, true, true);
    styleDrawAll(slicePurityDist, 999, 999, 999, 999, (base_path + "slicePurity_dist.pdf").c_str(), "topLeft", true, true, false, true, true);
    styleDrawAll(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs_weighted.pdf").c_str(), "topRight", true, true, false, true, true);
    styleDrawAll(sliceNumPFPsDist, 999, 999, 999, 999, (base_path + "sliceNumPFPs_dist.pdf").c_str(), "topRight", true, true, false, true, true);
    styleDrawAll(trackscoreHighestEnergyPFP, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_weighted.pdf").c_str(), "topRight", true, true, false, true, true);
    styleDrawAll(trackscoreHighestEnergyPFPDist, 999, 999, 999, 999, (base_path + "trackscoreHighestEnergyPFP_dist.pdf").c_str(), "topRight", true, true, false, true, true);
    styleDrawAll(trackscoreAllPFPs, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_weighted.pdf").c_str(), "topRight", true, true, false, true, true);
    styleDrawAll(trackscoreAllPFPsDist, 999, 999, 999, 999, (base_path + "trackscoreAllPFPs_dist.pdf").c_str(), "topRight", true, true, false, true, true);

    styleDrawAll(deltaTheta, 999, 999, 999, 999, (base_path + "deltaTheta_weighted.pdf").c_str(), "topRight", true, false, false, true, true);
    styleDrawAll(deltaThetaDist, 999, 999, 999, 999, (base_path + "deltaTheta_dist.pdf").c_str(), "topRight", true, false, false, true, true);
    styleDrawAll(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_weighted.pdf").c_str(), "topRight", true, false, false, true, true);
    styleDrawAll(ERecoSumThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco_dist.pdf").c_str(), "topRight", true, false, false, true, true);
    styleDrawAll(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_weighted.pdf").c_str(), "topRight", true, false, false, true, true);
    styleDrawAll(ERecoHighestThetaRecoDist, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco_dist.pdf").c_str(), "topRight", true, false, false, true, true);
}
