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
#include <TGaxis.h>

struct weights_struct{
    double signalCheated = 0;
    double signalDLUboone = 0;
    double signalDLDune = 0;
    double signalBDT = 0;
};

typedef struct{
    TCanvas* canvas;
    TH1D* baseHist;
    TH1D* nuECheated;
    TH1D* nuEDLDune;
    TH1D* nuEDLUboone;
    TH1D* nuEBDT;
} histGroup_struct;

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

    TH1D* base = new TH1D(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return{
        canvas,
        base,
        (TH1D*) base->Clone((baseName + "_nuECheated").c_str()),
        (TH1D*) base->Clone((baseName + "_nuEDLDune").c_str()),
        (TH1D*) base->Clone((baseName + "_nuEDLUboone").c_str()),
        (TH1D*) base->Clone((baseName + "_nuEBDT").c_str())
    };
}

void fillHistogram(histGroup_struct* hist, int DLCurrent, int signal, int type, double value, weights_struct* weight){
    if(!hist) return;

    TH1* target = nullptr;
    double w = 0.0;

    if(signal == 1 && DLCurrent == 0 && type == 1){
        w = weight->signalDLUboone;
        target = hist->nuEDLUboone;
    } else if(signal == 1 && DLCurrent == 1 && type == 1){
        w = weight->signalDLDune;
        target = hist->nuEDLDune;
    } else if(signal == 1 && DLCurrent == 2 && type == 1){
        w = weight->signalBDT;
        target = hist->nuEBDT;
    } else if(signal == 1 && DLCurrent == 6 && type == 1){
        w = weight->signalCheated;
        target = hist->nuECheated;
    }
    
    if(!target) return;

    target->Fill(value, w);
}

void styleDrawAll(histGroup_struct hists, double ymin, double ymax, double xmin, double xmax, const char* filename, const std::string& legendLocation){
    hists.canvas->cd();
    hists.canvas->SetTickx();
    hists.canvas->SetTicky();

    std::vector<TH1D*> allHists = {hists.nuECheated, hists.nuEDLDune, hists.nuEDLUboone, hists.nuEBDT};

    for(auto* hist : allHists)
        if(hist) hist->SetStats(0);

    hists.nuECheated->SetLineWidth(2);   hists.nuECheated->SetLineColor(TColor::GetColor("#bd1f01"));
    hists.nuEDLDune->SetLineWidth(2);    hists.nuEDLDune->SetLineColor(TColor::GetColor("#3f90da"));
    hists.nuEDLUboone->SetLineWidth(2);  hists.nuEDLUboone->SetLineColor(TColor::GetColor("#ffa90e"));
    hists.nuEBDT->SetLineWidth(2);       hists.nuEBDT->SetLineColor(TColor::GetColor("#e76300"));

    double minYValue = 0.0;
    double maxYValue = 0.0;
    for (auto* hist : allHists)
        if (hist && hist->GetMaximum() > maxYValue)
            maxYValue = hist->GetMaximum();
    

    if((ymin != 999) && (ymax != 999)){
        hists.nuECheated->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEDLDune->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEDLUboone->GetYaxis()->SetRangeUser(ymin, ymax);
        hists.nuEBDT->GetYaxis()->SetRangeUser(ymin, ymax);
    } else{
        hists.nuECheated->GetYaxis()->SetRangeUser(minYValue, maxYValue*1.1);
        hists.nuEDLDune->GetYaxis()->SetRangeUser(minYValue, maxYValue*1.1);
        hists.nuEDLUboone->GetYaxis()->SetRangeUser(minYValue, maxYValue*1.1);
        hists.nuEBDT->GetYaxis()->SetRangeUser(minYValue, maxYValue*1.1);
    }

    if((xmin != 999) && (xmax != 999)){
        hists.nuECheated->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEDLDune->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEDLUboone->GetXaxis()->SetRangeUser(xmin, xmax);
        hists.nuEBDT->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    hists.nuECheated->Draw("hist");
    hists.nuEDLDune->Draw("histsame");
    hists.nuEDLUboone->Draw("histsame");
    hists.nuEBDT->Draw("histsame");

    hists.nuECheated->SetStats(0);
    hists.nuECheated->GetXaxis()->SetTickLength(0.04);
    hists.nuECheated->GetXaxis()->SetTickSize(0.02);
    hists.nuECheated->GetYaxis()->SetTickLength(0.03);
    hists.nuECheated->GetYaxis()->SetTickSize(0.02);

    double Lxmin=0, Lxmax=0, Lymin=0, Lymax=0; 
    double height = std::max(0.025 * 4, 0.03);

    if(legendLocation == "topRight"){ Lxmin=0.62; Lymax=0.863; Lxmax=0.87; Lymin=Lymax - height; }
    else if(legendLocation == "topLeft"){ Lxmin=0.13; Lymax=0.863; Lxmax=0.38; Lymin=Lymax - height; }
    else if(legendLocation == "bottomRight"){ Lxmin=0.62; Lymin=0.137; Lxmax=0.87; Lymax=Lymin + height; }
    else if(legendLocation == "bottomLeft"){ Lxmin=0.13; Lymin=0.137; Lxmax=0.38; Lymax=Lymin + height; }

    auto legend = new TLegend(Lxmin, Lymin, Lxmax, Lymax);
    legend->AddEntry(hists.nuECheated, "Signal, Pandora Cheated Vertexing", "f");    
    legend->AddEntry(hists.nuEDLDune, "Signal, Pandora DL Vertexing: DUNE Train", "f");    
    legend->AddEntry(hists.nuEDLUboone, "Signal, Pandora DL Vertexing: #muBooNE Train", "f");    
    legend->AddEntry(hists.nuEBDT, "Signal, Pandora BDT Vertexing", "f");

    legend->SetTextSize(0.015);
    legend->SetMargin(0.12);
    legend->Draw();

    hists.canvas->SaveAs(filename);

}

void nuESignalVertexingComparison_macro(){

    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_nue_BDT_DLUboone_DLDune_Cheated_sameEvents.root");
    std::string base_path = "/nashome/c/coackley/nuESignalVertexingComparison/";

    gROOT->SetBatch(true);

    if (!gSystem->AccessPathName(base_path.c_str())) {
        gSystem->Exec(Form("rm -rf %s/*", base_path.c_str()));
    }
    gSystem->mkdir(base_path.c_str(), kTRUE);

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

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalDLDune;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalDLUboone;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalBDT;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignalCheated;

    double totalPOTSignalDLDune = 0;
    double totalPOTSignalDLUboone = 0;
    double totalPOTSignalBDT = 0;
    double totalPOTSignalCheated = 0;
    
    double totalPOTSignalDLDune_notMissing = 0;
    double totalPOTSignalDLUboone_notMissing = 0;
    double totalPOTSignalBDT_notMissing = 0;
    double totalPOTSignalCheated_notMissing = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        std::pair<unsigned int, unsigned int> key = std::make_pair(subRunRun, subRunNumber);

        if(subRunSignal == 1){
            if(subRunDLCurrent == 6 && seenSubRunsSignalCheated.find(key) == seenSubRunsSignalCheated.end()){
                totalPOTSignalCheated += subRunPOT;
                seenSubRunsSignalCheated.insert(key);
            } else if(subRunDLCurrent == 0 && seenSubRunsSignalDLUboone.find(key) == seenSubRunsSignalDLUboone.end()){
                totalPOTSignalDLUboone += subRunPOT;
                seenSubRunsSignalDLUboone.insert(key);
            } else if(subRunDLCurrent == 1 && seenSubRunsSignalDLDune.find(key) == seenSubRunsSignalDLDune.end()){
                totalPOTSignalDLDune += subRunPOT;
                seenSubRunsSignalDLDune.insert(key);
            } else if(subRunDLCurrent == 2 && seenSubRunsSignalBDT.find(key) == seenSubRunsSignalBDT.end()){
                totalPOTSignalBDT += subRunPOT;
                seenSubRunsSignalBDT.insert(key);
            }

            if(subRunDLCurrent == 6){
                totalPOTSignalCheated_notMissing += subRunPOT;
            } else if(subRunDLCurrent == 0){
                totalPOTSignalDLUboone_notMissing += subRunPOT;
            } else if(subRunDLCurrent == 1){
                totalPOTSignalDLDune_notMissing += subRunPOT;
            } else if(subRunDLCurrent == 2){
                totalPOTSignalBDT_notMissing += subRunPOT;
            }
        }

    }

    double targetPOT = 1e21;

    weights_struct weights;
    weights.signalCheated = targetPOT / totalPOTSignalCheated_notMissing;
    weights.signalDLUboone = targetPOT / totalPOTSignalDLUboone_notMissing;
    weights.signalDLDune = targetPOT / totalPOTSignalDLDune_notMissing;
    weights.signalBDT = targetPOT / totalPOTSignalBDT_notMissing;

    if(totalPOTSignalCheated_notMissing != totalPOTSignalCheated) std::cout << "Cheated POT doesn't match: totalPOTSignalCheated = " << totalPOTSignalCheated << ", totalPOTSignalCheated_notMissing = " << totalPOTSignalCheated_notMissing << std::endl;
    if(totalPOTSignalDLUboone_notMissing != totalPOTSignalDLUboone) std::cout << "DL Uboone POT doesn't match: totalPOTSignalDLUboone = " << totalPOTSignalDLUboone << ", totalPOTSignalDLUboone_notMissing = " << totalPOTSignalDLUboone_notMissing << std::endl;
    if(totalPOTSignalDLDune_notMissing != totalPOTSignalDLDune) std::cout << "DL Dune POT doesn't match: totalPOTSignalDLDune = " << totalPOTSignalDLDune << ", totalPOTSignalDLDune_notMissing = " << totalPOTSignalDLDune_notMissing << std::endl;
    if(totalPOTSignalBDT_notMissing != totalPOTSignalBDT) std::cout << "BDT POT doesn't match: totalPOTSignalBDT = " << totalPOTSignalBDT << ", totalPOTSignalBDT_notMissing = " << totalPOTSignalBDT_notMissing << std::endl;

    std::cout << "Weights Signal: Cheated = " << weights.signalCheated << ", DL Uboone = " << weights.signalDLUboone << ", DL Dune = " << weights.signalDLDune << ", BDT = " << weights.signalBDT << std::endl;

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
    std::vector<double> *angleRecalculationPCASlice5cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice5cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice10cm_sliceID = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCASlice15cm_sliceID = nullptr;
    
    std::vector<double> *angleRecalculationPCAPFP_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP5cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP10cm_pfpID = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_angle = nullptr;
    std::vector<double> *angleRecalculationPCAPFP15cm_pfpID = nullptr;

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
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_angle", &angleRecalculationPCASlice5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice5cm_sliceID", &angleRecalculationPCASlice5cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_angle", &angleRecalculationPCASlice10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice10cm_sliceID", &angleRecalculationPCASlice10cm_sliceID);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_angle", &angleRecalculationPCASlice15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCASlice15cm_sliceID", &angleRecalculationPCASlice15cm_sliceID);
    
    tree->SetBranchAddress("angleRecalculationPCAPFP_angle", &angleRecalculationPCAPFP_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP_pfpID", &angleRecalculationPCAPFP_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_angle", &angleRecalculationPCAPFP5cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP5cm_pfpID", &angleRecalculationPCAPFP5cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_angle", &angleRecalculationPCAPFP10cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP10cm_pfpID", &angleRecalculationPCAPFP10cm_pfpID);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_angle", &angleRecalculationPCAPFP15cm_angle);
    tree->SetBranchAddress("angleRecalculationPCAPFP15cm_pfpID", &angleRecalculationPCAPFP15cm_pfpID);

    Long64_t numEntries = tree->GetEntries();

    auto sliceCompleteness = createHistGroup("sliceCompleteness", "Slice Completeness", "Completeness", 102, 0, 1.02);
    auto slicePurity = createHistGroup("slicePurity", "Slice Purity", "Purity", 102, 0, 1);
    auto sliceCRUMBS = createHistGroup("sliceCRUMBS", "Slice CRUMBS Score", "CRUMBS Score", 25, -1, 1);
    auto sliceNumRecoNeut = createHistGroup("sliceNumRecoNeut", "Number of Reco Neutrinos in Slice", "Number of Reco Neutrinos", 10, 0, 10);
    auto sliceNumPFPs = createHistGroup("sliceNumPFPs", "Number of PFPs in Slice", "Number of PFPs", 20, 0, 20);
    auto sliceNumPrimaryPFPs = createHistGroup("sliceNumPrimaryPFPs", "Number of Primary PFPs in Slice", "Number of Primary PFPs", 20, 0, 20);
    auto ERecoSumThetaReco = createHistGroup("ERecoSumThetaReco", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaReco = createHistGroup("ERecoHighestThetaReco", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ETrueThetaReco = createHistGroup("ETrueThetaReco", "E_{true}#theta_{reco}^{2}", "E_{true}#theta_{reco}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoSumThetaTrue = createHistGroup("ERecoSumThetaTrue", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice", "E_{reco}#theta_{true}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto ERecoHighestThetaTrue = createHistGroup("ERecoHighestThetaTrue", "E_{reco}#theta_{true}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice", "E_{reco}#theta_{true}^{2} (MeV rad^{2})", 27, 0, 13.797);
    auto sliceAngleDifference = createHistGroup("angleDifference", "Angle Difference between True and Reconstructed Recoil Electron", "Angle (degrees)", 90, 0, 180);
    auto energyAsymmetrySummed = createHistGroup("energyAsymmetrySummed", "Energy Asymmetry between True and Summed Energies of PFPs in Slice", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto energyAsymmetryHighest = createHistGroup("energyAsymmetryHighest", "Energy Asymmetry between True and Highest Energy PFP in Slice", "(E_{true} - E_{reco})/E_{true}", 20, -1, 1);
    auto deltaVX = createHistGroup("deltaVX", "X Coordinate Difference Between True and Reco Neutrino Vertex", "x_{Reco} - x_{True} (cm)", 50, -5, 5);
    auto deltaVY = createHistGroup("deltaVY", "Y Coordinate Difference Between True and Reco Neutrino Vertex", "y_{Reco} - y_{True} (cm)", 50, -5, 5);
    auto deltaVZ = createHistGroup("deltaVZ", "Z Coordinate Difference Between True and Reco Neutrino Vertex", "z_{Reco} - z_{True} (cm)", 50, -5, 5);
    auto deltaR = createHistGroup("deltaR", "Distance Between True and Reco Neutrino Vertex", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 25, 0, 5);

    double actualSignalCount_cheated = 0;
    double actualSignalCount_DLDune = 0;
    double actualSignalCount_DLUboone = 0;
    double actualSignalCount_BDT = 0;

    double sliceSignalCount_cheated = 0;
    double sliceSignalCount_DLDune = 0;
    double sliceSignalCount_DLUboone = 0;
    double sliceSignalCount_BDT = 0;

    double xMin = -201.3; double xMax = 201.3;
    double yMin = -203.8; double yMax = 203.8;
    double zMin = 0;      double zMax = 509.4;

    for(Long64_t e = 0; e < numEntries; ++e){
        tree->GetEntry(e);

        double weight = 0;
        if(signal == 1 && DLCurrent == 6) weight = weights.signalCheated;
        if(signal == 1 && DLCurrent == 0) weight = weights.signalDLUboone;
        if(signal == 1 && DLCurrent == 1) weight = weights.signalDLDune;
        if(signal == 1 && DLCurrent == 2) weight = weights.signalBDT;

        int trueSignal = 0;

        if(nuEScatter == 1 && signal == 1 && DLCurrent == 6){
            if((((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))){
                actualSignalCount_cheated += weights.signalCheated;
                trueSignal = 1;
            }
        } else if(nuEScatter == 1 && signal == 1 && DLCurrent == 0){
            if((((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))){
                actualSignalCount_DLUboone += weights.signalDLUboone;
            }
        } else if(nuEScatter == 1 && signal == 1 && DLCurrent == 1){
            if((((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))){
                actualSignalCount_DLDune += weights.signalDLDune;
                trueSignal = 1;
            }
        } else if(nuEScatter == 1 && signal == 1 && DLCurrent == 2){
            if((((nuEScatterTrueVX > xMin) && (nuEScatterTrueVX < xMax)) && ((nuEScatterTrueVY > yMin) && (nuEScatterTrueVY < yMax)) && ((nuEScatterTrueVZ > zMin) && (nuEScatterTrueVZ < zMax)))){
                actualSignalCount_BDT += weights.signalBDT;
                trueSignal = 1;
            }
        }

        recoilElectron_struct recoilElectron;
        for(size_t i = 0; i < truth_recoilElectronPDG->size(); ++i){
            if(truth_recoilElectronPDG->size() > 1) std::cout << "More than 1 true recoil electron in event" << std::endl;
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

        // Looking at the reco slices
        if(reco_sliceID->size() == 0) continue;

        for(size_t slice = 0; slice < reco_sliceID->size(); ++slice){
            if(reco_sliceID->at(slice) == -999999) continue;
            // There is a reco slice in the event
            
            double sliceRecoVX = -999999;
            double sliceRecoVY = -999999;
            double sliceRecoVZ = -999999;

            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    sliceRecoVX = reco_neutrinoVX->at(recoNeut); 
                    sliceRecoVY = reco_neutrinoVY->at(recoNeut); 
                    sliceRecoVZ = reco_neutrinoVZ->at(recoNeut); 
                }
            }

            double sliceCategoryPlottingMacro = -999999;
            if(reco_sliceOrigin->at(slice) == 0){
                // This is a cosmic slice
                sliceCategoryPlottingMacro = 0;
            } else if(reco_sliceOrigin->at(slice) == 1){
                // This is a nu+e elastic scatter slice
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    if((sliceRecoVX < 201.3 && sliceRecoVX > -201.3) && (sliceRecoVY < 203.8 && sliceRecoVY > -203.8) && (sliceRecoVZ > 0 && sliceRecoVZ < 509.5)){
                        sliceCategoryPlottingMacro = 1;
                        // Signal slice
                    } else{
                        sliceCategoryPlottingMacro = 2;
                        // Signal fuzzy slice
                    }
                } else{
                    sliceCategoryPlottingMacro = 2;
                    // Signal fuzzy slice
                }
            } else if(reco_sliceOrigin->at(slice) == 3){
                // This is a BNB slice
                if(reco_sliceCompleteness->at(slice) > 0.5){
                    sliceCategoryPlottingMacro = 3;
                    // BNB slice
                } else{
                    sliceCategoryPlottingMacro = 4;
                    // BNB fuzzy slice
                }
            }

            highestEnergyPFP_struct highestEnergyPFP;
            double numPFPsSlice = 0;
            double numPrimaryPFPsSlice = 0;
            double summedPFPEnergy = 0;
            
            for(size_t pfp = 0; pfp < reco_particlePDG->size(); ++pfp){
                if(reco_particleSliceID->at(pfp) == reco_sliceID->at(slice)){
                    // PFP is in the slice
                    numPFPsSlice++;
                    if(reco_particleIsPrimary->at(pfp) == 1) numPrimaryPFPsSlice++;
                
                    summedPFPEnergy += reco_particleBestPlaneEnergy->at(pfp);

                    if(reco_particleBestPlaneEnergy->at(pfp) > highestEnergyPFP.energy){
                        highestEnergyPFP.energy = reco_particleBestPlaneEnergy->at(pfp);
                        highestEnergyPFP.theta = reco_particleTheta->at(pfp);
                        highestEnergyPFP.PFPID = reco_particleID->at(pfp);
                        highestEnergyPFP.dx = reco_particleDX->at(pfp);
                        highestEnergyPFP.dy = reco_particleDY->at(pfp);
                        highestEnergyPFP.dz = reco_particleDZ->at(pfp);
                        highestEnergyPFP.vx = reco_particleVX->at(pfp);
                        highestEnergyPFP.vy = reco_particleVY->at(pfp);
                        highestEnergyPFP.vz = reco_particleVZ->at(pfp);
                        highestEnergyPFP.completeness = reco_particleCompleteness->at(pfp);
                        highestEnergyPFP.purity = reco_particlePurity->at(pfp);
                        highestEnergyPFP.trackscore = reco_particleTrackScore->at(pfp);
                        highestEnergyPFP.primary = reco_particleIsPrimary->at(pfp);
                        highestEnergyPFP.truePDG = reco_particleTruePDG->at(pfp);
                        highestEnergyPFP.trueOrigin = reco_particleTrueOrigin->at(pfp);
                        highestEnergyPFP.trueInt = reco_particleTrueInteractionType->at(pfp);
                        highestEnergyPFP.bestPlanedEdx = reco_particleBestPlanedEdx->at(pfp);
                        highestEnergyPFP.razzledPDG11 = reco_particleRazzledPDG11->at(pfp);
                        highestEnergyPFP.razzledPDG13 = reco_particleRazzledPDG13->at(pfp);
                        highestEnergyPFP.razzledPDG22 = reco_particleRazzledPDG22->at(pfp);
                        highestEnergyPFP.razzledPDG211 = reco_particleRazzledPDG211->at(pfp);
                        highestEnergyPFP.razzledPDG2212 = reco_particleRazzledPDG2212->at(pfp);
                        highestEnergyPFP.razzledBestPDG = reco_particleRazzledBestPDG->at(pfp);
                        highestEnergyPFP.trueVX = reco_particleTrueVX->at(pfp);
                        highestEnergyPFP.trueVY = reco_particleTrueVY->at(pfp);
                        highestEnergyPFP.trueVZ = reco_particleTrueVZ->at(pfp);
                        highestEnergyPFP.trueEndX = reco_particleTrueEndX->at(pfp);
                        highestEnergyPFP.trueEndY = reco_particleTrueEndY->at(pfp);
                        highestEnergyPFP.trueEndZ = reco_particleTrueEndZ->at(pfp);
                        highestEnergyPFP.numHits = reco_particleNumHits->at(pfp);
                        highestEnergyPFP.clearCosmic = reco_particleClearCosmic->at(pfp);
                    }
                }
            }

            double angleDifference = -999999;
            if((highestEnergyPFP.dx != -999999) && (recoilElectron.dx != -999999)){
                double aDOTb = ((highestEnergyPFP.dx * recoilElectron.dx) + (highestEnergyPFP.dy * recoilElectron.dy) + (highestEnergyPFP.dz * recoilElectron.dz));
                double aMagnitude = std::sqrt((highestEnergyPFP.dx * highestEnergyPFP.dx) + (highestEnergyPFP.dy * highestEnergyPFP.dy) + (highestEnergyPFP.dz * highestEnergyPFP.dz));
                double bMagnitude = std::sqrt((recoilElectron.dx * recoilElectron.dx) + (recoilElectron.dy * recoilElectron.dy) + (recoilElectron.dz * recoilElectron.dz));
                double cosAngle = (aDOTb / (aMagnitude * bMagnitude));
                angleDifference = (TMath::ACos(cosAngle) * TMath::RadToDeg());
            } 

            double recoVX = -999999;
            double recoVY = -999999;
            double recoVZ = -999999;
            double numRecoNeutrinos = 0;

            for(size_t recoNeut = 0; recoNeut < reco_neutrinoID->size(); ++recoNeut){
                if(reco_neutrinoSliceID->at(recoNeut) == reco_sliceID->at(slice)){
                    // Reco neutrino is in the slice
                    numRecoNeutrinos++;
                    recoVX = reco_neutrinoVX->at(recoNeut);
                    recoVY = reco_neutrinoVY->at(recoNeut);
                    recoVZ = reco_neutrinoVZ->at(recoNeut);
                }
            }

            // Fill histograms here
            fillHistogram(&sliceCompleteness, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceCompleteness->at(slice), &weights);
            fillHistogram(&slicePurity, DLCurrent, signal, sliceCategoryPlottingMacro, reco_slicePurity->at(slice), &weights);
            fillHistogram(&sliceCRUMBS, DLCurrent, signal, sliceCategoryPlottingMacro, reco_sliceScore->at(slice), &weights);
            fillHistogram(&sliceNumRecoNeut, DLCurrent, signal, sliceCategoryPlottingMacro, numRecoNeutrinos, &weights);
            fillHistogram(&sliceNumPFPs, DLCurrent, signal, sliceCategoryPlottingMacro, numPFPsSlice, &weights);
            fillHistogram(&sliceNumPrimaryPFPs, DLCurrent, signal, sliceCategoryPlottingMacro, numPrimaryPFPsSlice, &weights);
            fillHistogram(&ERecoSumThetaReco, DLCurrent, signal, sliceCategoryPlottingMacro, (summedPFPEnergy*highestEnergyPFP.theta*highestEnergyPFP.theta), &weights);
            fillHistogram(&ERecoHighestThetaReco, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP.energy*highestEnergyPFP.theta*highestEnergyPFP.theta), &weights);
            fillHistogram(&ETrueThetaReco, DLCurrent, signal, sliceCategoryPlottingMacro, (recoilElectron.energy*highestEnergyPFP.theta*highestEnergyPFP.theta), &weights);
            fillHistogram(&ERecoSumThetaTrue, DLCurrent, signal, sliceCategoryPlottingMacro, (summedPFPEnergy*recoilElectron.angle*recoilElectron.angle), &weights);
            fillHistogram(&ERecoHighestThetaTrue, DLCurrent, signal, sliceCategoryPlottingMacro, (highestEnergyPFP.energy*recoilElectron.angle*recoilElectron.angle), &weights);
            fillHistogram(&sliceAngleDifference, DLCurrent, signal, sliceCategoryPlottingMacro, angleDifference, &weights);
            fillHistogram(&energyAsymmetrySummed, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - summedPFPEnergy)/recoilElectron.energy), &weights);
            fillHistogram(&energyAsymmetryHighest, DLCurrent, signal, sliceCategoryPlottingMacro, ((recoilElectron.energy - highestEnergyPFP.energy)/recoilElectron.energy), &weights);

        }

    }

    std::cout << "Numbers of true nu+e events:" << std::endl;
    std::cout << "Cheated = " << actualSignalCount_cheated << std::endl;
    std::cout << "DL Uboone = " << actualSignalCount_DLUboone << std::endl;
    std::cout << "DL Dune = " << actualSignalCount_DLDune << std::endl;
    std::cout << "BDT = " << actualSignalCount_BDT << std::endl; 

    styleDrawAll(sliceCompleteness, 999, 999, 999, 999, (base_path + "sliceCompletness.pdf").c_str(), "topRight");
    styleDrawAll(slicePurity, 999, 999, 999, 999, (base_path + "slicePurity.pdf").c_str(), "topRight");
    styleDrawAll(sliceCRUMBS, 999, 999, 999, 999, (base_path + "sliceCRUMBS.pdf").c_str(), "topRight");
    styleDrawAll(sliceNumRecoNeut, 999, 999, 999, 999, (base_path + "sliceNumRecoNeut.pdf").c_str(), "topRight");
    styleDrawAll(sliceNumPFPs, 999, 999, 999, 999, (base_path + "sliceNumPFPs.pdf").c_str(), "topRight");
    styleDrawAll(sliceNumPrimaryPFPs, 999, 999, 999, 999, (base_path + "sliceNumPrimaryPFPs.pdf").c_str(), "topRight");
    styleDrawAll(ERecoSumThetaReco, 999, 999, 999, 999, (base_path + "ERecoSumThetaReco.pdf").c_str(), "topRight");
    styleDrawAll(ERecoHighestThetaReco, 999, 999, 999, 999, (base_path + "ERecoHighestThetaReco.pdf").c_str(), "topRight");
    styleDrawAll(ETrueThetaReco, 999, 999, 999, 999, (base_path + "ETrueThetaReco.pdf").c_str(), "topRight");
    styleDrawAll(ERecoSumThetaTrue, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrue.pdf").c_str(), "topRight");
    styleDrawAll(ERecoHighestThetaTrue, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrue.pdf").c_str(), "topRight");
    styleDrawAll(sliceAngleDifference, 999, 999, 999, 999, (base_path + "angleDifference.pdf").c_str(), "topRight");
    styleDrawAll(energyAsymmetrySummed, 999, 999, 999, 999, (base_path + "energyAsymmetrySummed.pdf").c_str(), "topRight");
    styleDrawAll(energyAsymmetryHighest, 999, 999, 999, 999, (base_path + "energyAsymmetryHighest.pdf").c_str(), "topRight");

}
