#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutDefinitions.C"


//--------- Signal and BG definitions
std::vector<SampleDef> fSampleDefs = {
    {fTruthInAV+"IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Signal", true, "Signal"}
    ,{fTruthInAV+"IntOrigin==1 &&  IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "BG  #nu", false, "BG BNB"}
};

//---------  Cuts
std::vector<PlotDef> fCutDefs = cutDefsPIDFull;
PlotDef fMinimalCut = fCutDefs[0];

//---------  Load function
void LambdaPlot2DDistributions(){
    std::cout<<"Lambda 2D function loaded"<<std::endl;
    return;
}

void Make2DPlot(TTree* fTree, std::string outputDirName, int plotIndex, TCut minimalCut, PlotDef cut1, PlotDef cut2, std::vector<SampleDef> sampleDefs){

    // --- Create the TCanvas
    std::string plotLabel = "hAux2D";
    TCanvas *c3 = new TCanvas(("c3"+std::to_string(plotIndex)).c_str(),"2DSelection",0, 0, 1000, 750);
    CutStyler *fStyler = new CutStyler(0);
    c3->cd();

    // --- Cout
    std::cout<<"Making 2D plot "<<plotIndex<<std::endl;
    std::cout<<"  Cut1: "<<cut1.GetVarS()<<std::endl;
    std::cout<<"  Cut2: "<<cut2.GetVarS()<<std::endl;

    // --- Axis labels
    TString axisLabels = ";"+cut1.GetVarLabelS()+";"+cut2.GetVarLabelS();
    double fAlpha = 0.5;

    std::vector<TH2F*> hAux2DVec;
    for (size_t j = 0; j < sampleDefs.size(); ++j) {
    
        // --- Create the histogram
        std::string histLabel = plotLabel+std::to_string(j);
        TH2F *hAux2D = new TH2F( histLabel.c_str(), axisLabels,
                                cut1.GetBins().GetNBins(), cut1.GetBins().GetX1(), cut1.GetBins().GetX2(),
                                cut2.GetBins().GetNBins(), cut2.GetBins().GetX1(), cut2.GetBins().GetX2());
        hAux2D->SetStats(0);
        // --- Sample cut
        TCut sampelCut(sampleDefs[j].GetVar());

        // --- Draw distribution
        fTree->Draw( (cut2.GetVarS()+":"+cut1.GetVarS()+">>"+histLabel).c_str(), minimalCut && sampelCut, "p");

        // --- Draw the histogram
        hAux2DVec.push_back(hAux2D);        
    }

    // --- Draw the histograms
    // Legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Loop
    for (size_t j = 0; j < hAux2DVec.size(); ++j) {
        if (j==0) hAux2DVec[j]->Draw("p");
        else hAux2DVec[j]->Draw("p same");

        // --- Set colors and markers
        hAux2DVec[j]->SetMarkerColorAlpha(fStyler->GetColor(j), fAlpha);
        hAux2DVec[j]->SetMarkerStyle(fStyler->GetMarkerStyle(j));
        hAux2DVec[j]->SetMarkerSize(1.);

        // --- Add to legend
        legend->AddEntry(hAux2DVec[j], sampleDefs[j].GetLabel(), "p");
        
    }
    legend -> Draw("same");
    c3->Update();
    c3->SaveAs( (outputDirName+"/test2D_"+std::to_string(plotIndex)+".pdf").c_str());
    return;

}


//---------  Main function
void RunLambda2D(std::string fInputFileName="", bool batchMode=1, std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTreePost")
{

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "OutputPlots2D";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());

    //--------- Batch mode
    batchMode? gROOT->SetBatch(kTRUE): gROOT->SetBatch(kFALSE);

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    fTree->Print();

    //--------- Minimal cut
    TCut minimalCut(fMinimalCut.GetCut());

    // Make all possible 2D combinations
    int counterCut = 0;
    for (auto cut1 : fCutDefs){
        for (auto cut2 : fCutDefs){
            if (cut1.GetVarS() == cut2.GetVarS()) continue;
            // Make the 2D plot
            Make2DPlot(fTree, fOutputDirName, counterCut, minimalCut, cut1, cut2, fSampleDefs);
            counterCut++;
        }
    }
 

    return 0;
}
