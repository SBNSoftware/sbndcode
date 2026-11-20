#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

std::string fTruthInFV = "TruthIsFiducial==1 &&";
std::string fTruthInAV = "abs(NuvX)<200 && abs(NuvY)<200 && NuvZ>0 && NuvZ<500 && ";

//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = {
    {fTruthInAV+"IntOrigin==1 && IntDirt==0", "Nu", false}
    ,{"IntOrigin==1 && IntDirt==1", "Dirt", true}
    ,{"IntOrigin==2", "Cosmics", false}
};

//-------- POT normalization
double fPOTTotalNorm = 3.3e20;

//---------  Phase space cuts
std::vector<PlotDef> phaseSpaceDefs = {};// phaseSpaceVars;

//---------  LATeX output file
std::string fOutputFileName = "CutEfficiencies";
std::string fOutputFileNameNormalized = "CutEfficienciesNormalized";

TCut fCounterCut = "SliceID==0";

std::vector<PlotDef> cutDefsProtons = {
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kNone, 1, {0,2,2}, false, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ slice \\ in \\ FV"}
    ,{"MainMCParticleEnergy", "MainMCParticleEnergy", CutType::kNone, 1, {0,4,20}, false, "Main MCParticle E [GeV]",  "Main \\ MCParticle \\ E [GeV]"}
    ,{"MainMCParticlePDG", "MainMCParticlePDG", CutType::kCenter, 2212, {0,2500,100}, true, "Main MCParticle PDG",  "Main \\ MCParticle PDG \\ = \\ 2212"}
    ,{"!(IntOrigin==1 && IntDirt==0)",  "!(IntOrigin==1 && IntDirt==0)", CutType::kCenter, 1, {0,2,2}, true, "Dirt||Cosmic",  "Dirt||Cosmic"}

    ,{"MainMCParticleEnergy", "MainMCParticleEnergy", CutType::kNone, 1, {0.8,2.,40}, false, "Main MCParticle E [GeV]",  "Main \\ MCParticle \\ E [GeV]"}    
    ,{"MainMCParticleStartX", "MainMCParticleStartX", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle Start X [cm]",  "Main \\ MCParticle \\ Start \\ X [cm]"}
    ,{"MainMCParticleStartY", "MainMCParticleStartY", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle Start Y [cm]",  "Main \\ MCParticle \\ Start \\ Y [cm]"}
    ,{"MainMCParticleStartZ", "MainMCParticleStartZ", CutType::kNone, 1, {0,500,25}, false, "Main MCParticle Start Z [cm]",  "Main \\ MCParticle \\ Start \\ Z [cm]"}

    ,{"MainMCParticleEnergy", "MainMCParticleEnergy", CutType::kLeft, 1.1, {0.8,2.,40}, true, "Main MCParticle E [GeV]",  "Main \\ MCParticle \\ E \\ (GeV)"}    
    ,{"MainMCParticleStartX", "MainMCParticleStartX", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle Start X [cm]",  "Main \\ MCParticle \\ Start \\ X [cm]"}
    ,{"MainMCParticleStartY", "MainMCParticleStartY", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle Start Y [cm]",  "Main \\ MCParticle \\ Start \\ Y [cm]"}
    ,{"MainMCParticleStartZ", "MainMCParticleStartZ", CutType::kNone, 1, {0,500,25}, false, "Main MCParticle Start Z [cm]",  "Main \\ MCParticle \\ Start \\ Z [cm]"}
    
};

std::vector<PlotDef> cutDefsMichel = {
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kNone, 1, {0,2,2}, false, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"abs(RecnuvX)<210 && abs(RecnuvY)<210 && RecnuvZ>-10 && RecnuvZ<510", "(abs(RecnuvX)<210 && abs(RecnuvY)<210 && RecnuvZ>-10 && RecnuvZ<510)", CutType::kCenter, 1, {0,2,2}, true, "Reco in AV",   "Reco \\ slice \\ in \\ AV"} 
    ,{"MainMCParticleEnergy", "MainMCParticleEnergy", CutType::kNone, 1, {0,4,20}, false, "Main MCParticle E [GeV]",  "Main \\ MCParticle \\ E [GeV]"}
    ,{"abs(MainMCParticlePDG)", "abs(MainMCParticlePDG)", CutType::kCenter, 13, {0,2500,100}, true, "Main MCParticle PDG",  "Main \\ MCParticle PDG \\ = \\ 13"}
    //,{"!(IntOrigin==1 && IntDirt==0)",  "!(IntOrigin==1 && IntDirt==0)", CutType::kCenter, 1, {0,2,2}, true, "Dirt||Cosmic",  "Dirt||Cosmic"}

    ,{"MainMCParticleEnergy", "MainMCParticleEnergy", CutType::kNone, 1, {0.8,2.,40}, false, "Main MCParticle E [GeV]",  "Main \\ MCParticle \\ E [GeV]"}    
    ,{"MainMCParticleEndX", "MainMCParticleEndX", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle End X [cm]",  "Main \\ MCParticle \\ End \\ X [cm]"}
    ,{"MainMCParticleEndY", "MainMCParticleEndY", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle End Y [cm]",  "Main \\ MCParticle \\ End \\ Y [cm]"}
    ,{"MainMCParticleEndZ", "MainMCParticleEndZ", CutType::kNone, 1, {0,500,25}, false, "Main MCParticle End Z [cm]",  "Main \\ MCParticle \\ End \\ Z [cm]"}

    ,{"MainMCParticleNDaughters", "MainMCParticleNDaughters", CutType::kNone, 1, {0,10,10}, false, "Main MCParticle N Daughters",  "Main \\ MCParticle \\ N \\ Daughters"}
    ,{"MainMCParticleMichelDecay", "MainMCParticleMichelDecay", CutType::kCenter, 1, {0,2,2}, true, "Main MCParticle Michel Decay",  "Main \\ MCParticle \\ Michel \\ Decay"}
    ,{"MainMCParticleEndX", "MainMCParticleEndX", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle End X [cm]",  "Main \\ MCParticle \\ End \\ X [cm]"}
    ,{"MainMCParticleEndY", "MainMCParticleEndY", CutType::kNone, 1, {-200,200,20}, false, "Main MCParticle End Y [cm]",  "Main \\ MCParticle \\ End \\ Y [cm]"}
    ,{"MainMCParticleEndZ", "MainMCParticleEndZ", CutType::kNone, 1, {0,500,25}, false, "Main MCParticle End Z [cm]",  "Main \\ MCParticle \\ End \\ Z [cm]"}
    
};

std::vector<PlotDef> cutDefs = cutDefsMichel;

//---------  Main function
void DirtAnalysis(){
    std::cout<<"Load Dirt Analysis\n";
}

//---------  Main function
void RunDirtAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName+"/PhaseSpace").c_str());



    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    

    //----------------- POT normalization
    double potScaling = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, fPOTTotalNorm, potScaling, potScalingSignal);


    //--------- Matrix to store the number of events
    std::vector< std::vector<int> > nEventsMatrix;
    std::vector<PlotDef> cutDefsForTable;

    //--------- Loop over the cuts
    std::vector<AnaPlot> anaPlots;
    TCut previousCut("");
    for (size_t i = 0; i < cutDefs.size(); ++i) {

        TCut currentCut = TCut(cutDefs[i].GetCut());
        AnaPlot anaPlot(i, cutDefs[i], sampleDefs, phaseSpaceDefs);

        anaPlot.DrawHistograms(fTree, previousCut);
        anaPlots.push_back(anaPlot);

        if(cutDefs[i].GetAccumulateCut()){
            previousCut = previousCut && currentCut;
            anaPlot.DrawHistograms(fTree, previousCut, 1);

            cutDefsForTable.push_back(cutDefs[i]);
            nEventsMatrix.push_back({});
            size_t cutIndex = nEventsMatrix.size()-1;
            std::map<std::string, int> nEventsMap = anaPlot.GetCountsV();
            for(const auto& sample : sampleDefs){
                nEventsMatrix[cutIndex].push_back(nEventsMap[sample.GetLabelS()]);
            }
        }
    }


    //loop over the samples and set NEvents
    for(auto& sample : sampleDefs){
        int nEvents = fTree->Draw( "", TCut(sample.GetVar())+fCounterCut, "goff");
        sample.SetNEvents(nEvents);
        std::cout<<"Sample: "<<sample.GetLabelS()<<" NEvents: "<<nEvents<<std::endl;
    }
    
    //--------- Create the LaTeX table
    GenerateAndCompileTeXTable(cutDefsForTable, sampleDefs, nEventsMatrix, fOutputFileName, "Cut efficiencies", fOutputDirName, 1, 1, 1, 0);

    //--------- Create the LaTeX table (POT normalized)
    GenerateAndCompileTeXTable(cutDefsForTable, sampleDefs, nEventsMatrix, fOutputFileNameNormalized, "Cut efficiencies POT normalized", fOutputDirName, potScaling, potScalingSignal, fPOTTotalNorm, 0);
   
    //--------- Output hand scans
    CreateHandScanList(fTree, fTreeHeader, previousCut, sampleDefs);
    

    return;
}


