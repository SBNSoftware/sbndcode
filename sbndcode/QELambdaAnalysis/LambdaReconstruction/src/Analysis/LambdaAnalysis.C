#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutDefinitions.C"
#include "CutEfficienciesLATeXInterface.C"

//--------- Settings ---------
//--------- Scale POT
bool fScaleHistogramsToPOT = 1;
//--------- Plot all the variables
bool fPlotAllVars = 0;
//-------- POT normalization
double fPOTTotalNorm = 10e20;
//---------  LATeX output file
std::string fOutputFileName = "CutEfficiencies";
std::string fOutputFileNameNormalized = "CutEfficienciesNormalized";

//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = fSampleDefsStd;
std::vector<SampleDef> sampleDefsExtended = fSampleDefsStd;

//fSampleDefsStd;

//---------  Phase space cuts
std::vector<PlotDef> fPhaseSpaceDefs = {};//psLambdaKinematics;

//---------  Cuts
std::vector<PlotDef> fCutDefs = cutDefsFinalSelection;
//std::vector<PlotDef> fCutDefs = cutDefsTalk2Induction;
//std::vector<PlotDef> fCutDefs = cutDefsOriginsDistributions;

std::vector<PlotDef> fCutDefsCol = cutDefsPIDFull;
std::vector<PlotDef> fCutDefsInd = cutDefsInd;


//---------  Load function
void LambdaAnalysis(){
    std::cout<<"LambdaAnalysis function loaded"<<std::endl;
    return;
}


//---------  Main function
void RunLambdaAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTree")
{

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName+"/PhaseSpace").c_str());
    gSystem->Exec(("mkdir "+fOutputDirName+"/OtherDistributions").c_str());

    //--------- Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //----------------- POT normalization
    double potScalingBg = 1;
    double potScalingSignal = 1;
    double potScalingBgForPlots = 1;
    double potScalingSignalForPlots = 1;

    ReadPOT(fFile, fPOTTotalNorm, potScalingBg, potScalingSignal, (fTreeDirName+"pottree").c_str());
    if(fScaleHistogramsToPOT){
        potScalingBgForPlots = potScalingBg;
        potScalingSignalForPlots = potScalingSignal;
    }
    std::cout<<"POT scaling: "<<potScalingBg<<" POT scaling signal: "<<potScalingSignal<<std::endl;

    //--------- Matrix to store the number of events
    std::vector< std::vector<int> > nEventsMatrix;
    std::vector<PlotDef> cutDefsForTable;

    //--------- Get vector with the cuts to accumulate
    std::vector<PlotDef> cutsToAccumulate;
    if(fPlotAllVars){
        for (size_t i = 0; i < fCutDefs.size(); ++i) {
            if(fCutDefs[i].GetAccumulateCut()){
                cutsToAccumulate.push_back(fCutDefs[i]);
            }
        }
    }

    fTree->Print();
    //--------- Loop over the cuts
    std::vector<AnaPlot> anaPlots;
    TCut previousCut("");
    for (size_t i = 0; i < fCutDefs.size(); ++i) {
        
        // Get the current cut
        TCut currentCut = TCut(fCutDefs[i].GetCut());

        // Create the plot handle
        AnaPlot anaPlot(i, fCutDefs[i], sampleDefs, fPhaseSpaceDefs, cutsToAccumulate);

        // Draw the histograms
        anaPlot.DrawHistograms(fTree, previousCut, 0, potScalingSignalForPlots, potScalingBgForPlots);
        anaPlots.push_back(anaPlot);

        // Id the cut is to be accumulated
        if(fCutDefs[i].GetAccumulateCut()){
            // Store the previous cut
            previousCut = previousCut && currentCut;

            // Draw the histograms again
            anaPlot.DrawHistograms(fTree, previousCut, 1, potScalingSignalForPlots, potScalingBgForPlots);

            // Store the numbers for final efficiency table
            cutDefsForTable.push_back(fCutDefs[i]);
            nEventsMatrix.push_back({});
            size_t cutIndex = nEventsMatrix.size()-1;
            std::map<std::string, int> nEventsMap = anaPlot.GetCountsV();
            for(const auto& sample : sampleDefs){
                nEventsMatrix[cutIndex].push_back(nEventsMap[sample.GetLabelS()]);
            }

        }
    
    }

    // --- Loop over the samples and set NEvents
    for(auto& sample : sampleDefs){
        int nEvents = fTree->Draw( "", TCut(sample.GetVar())+fCounterCut, "goff");
        sample.SetNEvents(nEvents);
        std::cout<<"Sample: "<<sample.GetLabelS()<<" NEvents: "<<nEvents<<std::endl;
    }

    // --- Final plot with all the cuts, all samples
    AnaPlot anaPlotFinal(-1, fCutDefs.back(), sampleDefs, fPhaseSpaceDefs, {});
    anaPlotFinal.DrawHistograms(fTree, previousCut, 0, potScalingSignalForPlots, potScalingBgForPlots);

    //--------- Create the LaTeX table
    MakeCutFlowPlot(cutDefsForTable, sampleDefsExtended, nEventsMatrix, fOutputDirName, potScalingBg, potScalingSignal, fPOTTotalNorm);
    
    //--------- Create the LaTeX table
    GenerateAndCompileTeXTable(cutDefsForTable, sampleDefs, nEventsMatrix, fOutputFileName, "Cut efficiencies", fOutputDirName);

    //--------- Create the LaTeX table (POT normalized)
    GenerateAndCompileTeXTable(cutDefsForTable, sampleDefs, nEventsMatrix, fOutputFileNameNormalized, "Cut efficiencies", fOutputDirName, potScalingBg, potScalingSignal, fPOTTotalNorm);
   
    //--------- Output hand scans
    CreateHandScanList(fTree, fTreeHeader, previousCut, sampleDefs, 0);

    return;
}



void RunCutLoopAnalysis(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTree")
{

    //Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputCutsLoop";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //----------------- POT normalization
    double potScalingBg = 1;
    double potScalingSignal = 1;
    double potScalingBgForPlots = 1;
    double potScalingSignalForPlots = 1;

    ReadPOT(fFile, fPOTTotalNorm, potScalingBg, potScalingSignal, (fTreeDirName+"pottree").c_str());
    if(fScaleHistogramsToPOT){
        potScalingBgForPlots = potScalingBg;
        potScalingSignalForPlots = potScalingSignal;
    }
    std::cout<<"POT scaling: "<<potScalingBg<<" POT scaling signal: "<<potScalingSignal<<std::endl;
    
    
    //---------- Set of cuts to loop over
    std::vector<PlotDef> loopCuts = cutDefsLoop;   
    // Minimal cut
    PlotDef minimalCut = minimalCutForLoop;
    
    // Loop over all possible orders using std::next_permutation
    int permutationId = 0;
    do {
        
        // Start Cut
        TCut currentCut("");

        // Vector of cuts including the minimal cut
        std::vector<PlotDef> loopCutsWithMinimal = loopCuts;
        loopCutsWithMinimal.insert(loopCutsWithMinimal.begin(), minimalCut);
        loopCutsWithMinimal.insert(loopCutsWithMinimal.begin(), startCutForLoop);

        //--------- Matrix to store the number of events
        std::vector< std::vector<int> > nEventsMatrix;
        nEventsMatrix.resize(loopCutsWithMinimal.size());
        for (size_t i = 0; i < loopCutsWithMinimal.size(); ++i) {
            nEventsMatrix[i].resize(sampleDefs.size());
        }

        for (const auto& cut : loopCutsWithMinimal) {
            std::cout << "Adding cut: " << cut.GetCut() << "\n";
            currentCut+=TCut(cut.GetCut());

            for(const auto& sample : sampleDefs){
                TCut currentSampleCut = TCut(sample.GetVar());

                int n = fTree->Draw( "", currentCut+currentSampleCut, "goff");
                nEventsMatrix[&cut - &loopCutsWithMinimal[0]][&sample - &sampleDefs[0]] = n;
            }
            
        }


        //--------- Create the LaTeX table
        GenerateAndCompileTeXTable(loopCutsWithMinimal, sampleDefs, nEventsMatrix, fOutputFileName+std::to_string(permutationId), "Cut efficiencies "+std::to_string(permutationId), fOutputDirName, potScalingBg, potScalingSignal, fPOTTotalNorm);

        permutationId++;


    } while (std::next_permutation(loopCuts.begin(), loopCuts.end()));

    return;
}



//---------  Main function
void RunTreeSelector(std::string fInputFileName="", bool scalePOT=1, bool batchMode=1, std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTreePost")
{


    //--------- Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);
    

    //--------- Input TTree list
    std::vector<std::string> treeNames{
        "originsAnaPost/LambdaAnaTree",
        "originsAnaPostU/LambdaAnaTree",
        "originsAnaPostV/LambdaAnaTree"
    };

    std::string inductionULabel = "originsAnaPostU";
    std::string inductionVLabel = "originsAnaPostV";

    //---------  Maps to store number of events per signal type
    std::map<std::string, int> nEventsMap;
    for(const auto& sample : sampleDefs){
        nEventsMap[sample.GetLabelS()] = 0;
    }

    //--------- Map to store the entries per sample and initialize it
    std::map<std::string, std::set<int>> entriesMapU;
    std::map<std::string, std::set<int>> entriesMapV;
    std::map<std::string, std::set<int>> entriesMapC;

    for(const auto& sample : sampleDefs){
        entriesMapU[sample.GetLabelS()] = {};
        entriesMapV[sample.GetLabelS()] = {};
        entriesMapC[sample.GetLabelS()] = {};
    }

    //--------- Input TFile
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");

    
    //----------------- POT normalization
    double potScalingBg = 1;
    double potScalingSignal = 1;
    ReadPOT(fFile, fPOTTotalNorm, potScalingBg, potScalingSignal, (fTreeDirName+"pottree").c_str());
    std::cout<<"POT scaling: "<<potScalingBg<<" POT scaling signal: "<<potScalingSignal<<std::endl;

    //--------- Loop over the TTrees
    for(const auto& treeName : treeNames){
        TTree *fTree = (TTree *)fFile->Get(treeName.c_str());
        std::cout<<"\n\n Tree: "<<treeName<<std::endl;


        // Fill the nEventsMap, no Cut
        for(const auto& sample : sampleDefs){
            int nEvents = fTree->Draw( "", TCut(sample.GetVar())+TCut("TruthIsFiducial==1"), "goff");
            nEventsMap[sample.GetLabelS()] = nEvents;
            std::cout<<"Sample: "<<sample.GetLabelS()<<" NEvents: "<<nEvents<<std::endl;
        }


        //--------- Get the accumulated cut   
        std::vector<PlotDef> thisCutDefs = fCutDefsCol;
        if(treeName.find(inductionULabel)!=std::string::npos || treeName.find(inductionVLabel)!=std::string::npos){
            std::cout<<"Using induction cuts"<<std::endl;
            thisCutDefs = fCutDefsInd;
        }
        TCut fAllCuts("");
        for (size_t i = 0; i < thisCutDefs.size(); ++i) {
            if(thisCutDefs[i].GetAccumulateCut()){
                fAllCuts+=TCut(thisCutDefs[i].GetCut());
            }
        }

        std::cout<<"Accumulated cuts: "<<fAllCuts<<std::endl;

        //--------- Get the entry list for each sample
        for(const auto& sample : sampleDefs){
            fTree->Draw(">>elist", fAllCuts+TCut(sample.GetVar()), "entrylist");
            TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
            // Cout the entries

            std::cout<<"Entries for "<<sample.GetLabelS()<<": "<<elist->GetN()<<std::endl;
            for (size_t i = 0; i < elist->GetN(); ++i) {
                if(treeName.find(inductionULabel)!=std::string::npos){
                    entriesMapU[sample.GetLabelS()].insert(elist->GetEntry(i));
                }
                else if(treeName.find(inductionVLabel)!=std::string::npos){
                    entriesMapV[sample.GetLabelS()].insert(elist->GetEntry(i));
                }
                else if(treeName.find("originsAna")!=std::string::npos){
                    entriesMapC[sample.GetLabelS()].insert(elist->GetEntry(i));
                    // fill nEventsMap, no Cut


                }
            }
        }
    }

    // cout summary
    std::cout<<"\n\n Summary: \n\n";
    for(const auto& sample : sampleDefs){
        std::cout<<"Sample: "<<sample.GetLabelS()<<" --- NEvents: "<<nEventsMap[sample.GetLabelS()]<<std::endl;

        std::cout<<"Sample: "<<sample.GetLabelS()<<" --- Entries U: "<<entriesMapU[sample.GetLabelS()].size()<<"   eff: "<<100.*entriesMapU[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()]<<std::endl;
        std::cout<<"Sample: "<<sample.GetLabelS()<<" --- Entries V: "<<entriesMapV[sample.GetLabelS()].size()<<"   eff: "<<100.*entriesMapV[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()]<<std::endl;
        std::cout<<"Sample: "<<sample.GetLabelS()<<" --- Entries C: "<<entriesMapC[sample.GetLabelS()].size()<<"   eff: "<<100.*entriesMapC[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()]<<std::endl;
    }


    std::cout<<"\n\n Intersection sets: \n\n";

    // Get intersection of entries for U and V
    std::map<std::string, std::set<int>> intersectionMap;
    for(const auto& sample : sampleDefs){  
        std::set<int> intersection;
        std::set_intersection(entriesMapU[sample.GetLabelS()].begin(), entriesMapU[sample.GetLabelS()].end(), entriesMapV[sample.GetLabelS()].begin(), entriesMapV[sample.GetLabelS()].end(), std::inserter(intersection, intersection.begin()));
        std::cout<<"Sample: "<<sample.GetLabelS()<<" Entries U and V: "<<intersection.size()<<"   eff: "<<100.*intersection.size()/nEventsMap[sample.GetLabelS()]<<std::endl;
        intersectionMap[sample.GetLabelS()] = intersection;
    }


    std::cout<<"\n\n Union sets: \n\n";
    // Get union with C
    std::map<std::string, std::set<int>> unionMap;
    for(const auto& sample : sampleDefs){  
        // union between intersection and C
        std::set<int> unionSet;
        std::set_union(entriesMapC[sample.GetLabelS()].begin(), entriesMapC[sample.GetLabelS()].end(), intersectionMap[sample.GetLabelS()].begin(), intersectionMap[sample.GetLabelS()].end(), std::inserter(unionSet, unionSet.begin()));

        std::cout<<"Sample: "<<sample.GetLabelS()<<" Entries U and V || C: "<<unionSet.size()<<"   eff: "<<100.*unionSet.size()/nEventsMap[sample.GetLabelS()]<<std::endl;
        unionMap[sample.GetLabelS()] = unionSet;
    }

    // Make all the above an output table with efficiences between parenthesis
    int sep = 15;
    std::cout<<"\n\n Summary table: \n\n";
    std::cout << std::setw(2*sep) << "Sample" << std::setw(sep) << "Events" << std::setw(sep) << "U" << std::setw(sep) << "V" << std::setw(sep) << "C" << std::setw(sep) << "U&V" << std::setw(sep) << "(U&V)|C" << std::endl;
    for(const auto& sample : sampleDefs){
       
        double potScale = 1.;
        if(scalePOT==1){
            if(sample.IsSignal()==1){
                potScale = potScalingSignal;
            }
            else{
                potScale = potScalingBg;
            }
        }   

        std::cout << std::setw(2*sep) << sample.GetLabelS() << std::setw(sep) << nEventsMap[sample.GetLabelS()] << std::setw(sep) << potScale*entriesMapU[sample.GetLabelS()].size() << std::setw(sep) << potScale*entriesMapV[sample.GetLabelS()].size() << std::setw(sep) << potScale*entriesMapC[sample.GetLabelS()].size() << std::setw(sep) << potScale*intersectionMap[sample.GetLabelS()].size() << std::setw(sep) << potScale*unionMap[sample.GetLabelS()].size() << std::endl;
    }

    // Make all the above an output table with efficiences between parenthesis, 3 decimal places and with %
     std::cout<<"\n\n Summary efficiency table: \n\n";
    std::cout << std::setw(2*sep) << "Sample" << std::setw(sep) << "Events" << std::setw(sep) << "U" << std::setw(sep) << "V" << std::setw(sep) << "C" << std::setw(sep) << "U&V" << std::setw(sep) << "(U&V)|C" << std::endl;
    for(const auto& sample : sampleDefs){
        std::cout << std::setw(2*sep) << sample.GetLabelS() << std::setw(sep) << nEventsMap[sample.GetLabelS()] << std::setw(sep) << std::setprecision(3) << 100.*entriesMapU[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()] << " %" << std::setw(sep) << std::setprecision(3) << 100.*entriesMapV[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()] << " %" << std::setw(sep) << std::setprecision(3) << 100.*entriesMapC[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()] << " %" << std::setw(sep) << std::setprecision(3) << 100.*intersectionMap[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()] << " %" << std::setw(sep) << std::setprecision(3) << 100.*unionMap[sample.GetLabelS()].size()/nEventsMap[sample.GetLabelS()] << " %" << std::endl;
    }
    



    return;

}
