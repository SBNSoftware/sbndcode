#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <vector>

#include "TMVA/Reader.h"

#include "../FRANS/FRANSObj.h"


std::map<std::string, bool> extractLabels(const std::string& xmlDataFile) {
    std::map<std::string, bool> labelVector;

    std::ifstream file(xmlDataFile);

    if (!file.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return labelVector;
    }

    std::string xmlData((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    size_t pos = 0;

    // Search for the beginning of the <Variable> tags
    while ((pos = xmlData.find("<Variable", pos)) != std::string::npos) {
        // Find the Label attribute within the current <Variable> tag
        size_t labelPos = xmlData.find("Label=\"", pos);
        if (labelPos != std::string::npos) {
            labelPos += 7; // Move past the "Label=\"" part
            // Find the closing quote of the Label attribute
            size_t endQuotePos = xmlData.find("\"", labelPos);
            if (endQuotePos != std::string::npos) {
                // Extract the Label value and add it to the vector
                std::string label = xmlData.substr(labelPos, endQuotePos - labelPos);
                labelVector[label] = true;
            }
        }

        // Move to the next position after the current <Variable> tag
        pos += 1;
    }

    return labelVector;
}

int MacroFRANSEvaluate(
    const std::string& fWeightsPath = "FRANSSelectionTMVA_BDT.weights.xml",
    const std::string& fInputFileName = "FRANSAnaOutput.root",
    const std::string& view = "2",
    const std::string fTreeDirName = "FRANSCheatedVx/",
    const std::string fTreeName = "FRANSTree")
{
    // ----- Configuration -----
    const double scoreCut = 0.1;
    bool evaluteSignalOnly = true;

    // ----- Open the input file -----
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    FRANSObj *fFRANSObj = new FRANSObj();
    std::string fFRANSObjName = "FRANSObj"+view;
    fTree->SetBranchAddress(fFRANSObjName.c_str(), &fFRANSObj);
    // Set Branch addresses for MC variables
    int fIsSignal;
    double fGap, fProtonKE, fPionKE;
    fTree->SetBranchAddress("IsSignal", &fIsSignal);
    fTree->SetBranchAddress("Gap", &fGap);
    fTree->SetBranchAddress("ProtonKE", &fProtonKE);
    fTree->SetBranchAddress("PionKE", &fPionKE);


    // ----- Prepare the TMVA reader -----
    TMVA::Reader *fReader = new TMVA::Reader();
    // extract the labels from the xml file
    std::map<std::string, bool> fVarLabels = extractLabels(fWeightsPath);
    // print variables to use
    std::cout << "Variables to use:" << std::endl;
    for (auto const& x : fVarLabels) {
        std::cout << x.first << std::endl;
    }
    float alpha, iota, eta, delta, fitScore;
    float gap, protonKE, pionKE;
    if (fVarLabels["FRANSObj"+view+".fAlpha"]==true) fReader->AddVariable("FRANSObj"+view+".fAlpha", &alpha);
    if (fVarLabels["FRANSObj"+view+".fIota"]==true) fReader->AddVariable("FRANSObj"+view+".fIota", &iota);
    if (fVarLabels["FRANSObj"+view+".fEta"]==true) fReader->AddVariable("FRANSObj"+view+".fEta", &eta);
    if (fVarLabels["FRANSObj"+view+".fDelta"]==true) fReader->AddVariable("FRANSObj"+view+".fDelta", &delta);
    if (fVarLabels["FRANSObj"+view+".fFitScore"]==true) fReader->AddVariable("FRANSObj"+view+".fFitScore", &fitScore);
    fReader->AddSpectator("Gap", &gap);
    fReader->AddSpectator("ProtonKE", &protonKE);
    fReader->AddSpectator("PionKE", &pionKE);

    fReader->BookMVA("FRANS BDT", fWeightsPath.c_str());    


    TEfficiency* efficiency = new TEfficiency("effS", "Selection Efficiency;Gap [cm];#epsilon", 21, -1, 20);
    TH2F hC("hC", ";Gap [cm]; Score", 200, -1, 20, 200, -1., 1.);
    // TH1s for the score signal and background
    TH1F hS("hS", ";Score", 40, -1., 1.);
    TH1F hB("hB", ";Score", 40, -1., 1.);

    for (Long64_t ievt = 0; ievt < fTree->GetEntries(); ievt++) {
        fTree->GetEntry(ievt);

        alpha = fFRANSObj->GetAlpha();
        eta = fFRANSObj->GetEta();
        delta = fFRANSObj->GetDelta();
        fitScore = fFRANSObj->GetFitScore();
        iota = fFRANSObj->GetIota();
        // get the score
        float score = fReader->EvaluateMVA("FRANS BDT");

        if (fIsSignal) {
            hS.Fill(score);
        } else {
            hB.Fill(score);
        }

        if(evaluteSignalOnly && fIsSignal){
            bool passCut = score > scoreCut;
            efficiency->Fill(passCut, fGap);
            hC.Fill(fGap, score);
        }
                
    }


    // --- Draw all the histograms ---    
    TCanvas* c1 = new TCanvas("example", "", 1000, 750);

    // Declare TPad
    TPad *pad1 = new TPad("pad1","pad1",0,0.5,0.5,1);
    TPad *pad2 = new TPad("pad2","pad2",0.5,0.5,1,1);
    TPad *pad3 = new TPad("pad3","pad3",0,0,0.5,0.5);
    TPad *pad4 = new TPad("pad4","pad4",0.5,0,1,0.5);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();

    
    pad1->cd();
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    efficiency->Draw("AP");
    gPad->Update();
    efficiency->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);
    efficiency->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.2);

    pad2->cd();
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    hC.GetXaxis()->SetTitleOffset(1.);
    hC.GetYaxis()->SetTitleOffset(1.);
    hC.Draw();

    pad3->cd();
    hB.Draw();
    hS.Draw("same");
    hS.SetLineColor(kRed);
    hS.SetLineWidth(2);
    hB.SetLineColor(kBlue);
    hB.SetLineWidth(2);
    hB.SetStats(0);
    //Legend
    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(&hS,"Signal","l");
    legend->AddEntry(&hB,"Background","l");
    legend->Draw("same");


    c1->cd();
    c1->Update();
    c1->WaitPrimitive();

    delete efficiency;
    return 0;
}
