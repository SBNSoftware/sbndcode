#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "LambdaMVAAnalysis.C"
#include "LambdaEvaluateMVAAnalysis.C"
#include "LambdaAnalysis.C"


void LoadAnalysisMacros(){
    std::cout<<" Loaded ana macros\n";
    return;
}

void MacroEvaluateMVAAnalysis(std::string fInputFileName="", bool batchMode=1, std::string method="BDT", std::string configFile = "configMVA.txt", std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{
    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTrees
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, 3.3e20, potScaling, potScalingSignal);
    std::cout<<"POT scaling: "<<potScaling<<" POT scaling signal: "<<potScalingSignal<<std::endl;
    
    RunEvaluateMVAAnalysis(fTree, fTreeHeader, method, "dataset/weights/", configFile, potScaling, potScalingSignal, 3.3e20);
}


void MacroMVAAnalysis(std::string fInputFileName="", double nTrainFrac = -1, std::string configFile="configMVA.txt", std::string method="BDT", bool useBatchMode=false, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree"){
    RunLambdaMVAAnalysis(fInputFileName, nTrainFrac, configFile, method, fTreeDirName, fTreeName);
}


void MacroRunAndEvaluateMVAAnalysis(std::string fInputFileName, std::string fInputFileNameTest, std::string configFile="configMVA.txt", std::string method="BDT", bool batchMode=1, std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTree"){

    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    // First make the training
    RunLambdaMVAAnalysis(fInputFileName, -1, configFile, method, fTreeDirName, fTreeName);

    // Then evaluate the BDT
    //--------- Input TTrees
    TFile *fFile = new TFile(fInputFileNameTest.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, 3.3e20, potScaling, potScalingSignal);
    std::cout<<"POT scaling: "<<potScaling<<" POT scaling signal: "<<potScalingSignal<<std::endl;
    
    RunEvaluateMVAAnalysis(fTree, fTreeHeader, method, "dataset/weights/", configFile, potScaling, potScalingSignal, 3.3e20);
}