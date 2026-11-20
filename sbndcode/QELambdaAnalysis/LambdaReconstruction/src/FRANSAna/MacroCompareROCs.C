#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <vector>

std::vector<int> colors = {kAzure-3, kOrange+7, kRed, kGreen+3, kMagenta, kBlack, kBlue, kOrange, kRed+2, kGreen+2, kMagenta+2, kBlack+2, kBlue+2, kOrange-3, kRed-3, kGreen-3, kMagenta-3, kBlack-3, kBlue-3, kOrange-6, kRed-6, kGreen-6, kMagenta-6, kBlack-6, kBlue-6, kOrange-9, kRed-9, kGreen-9, kMagenta-9, kBlack-9, kBlue-9};

std::vector<int> lines = {1, 2, 3, 4, 5, 6, 7, 8};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMTMVAResultsNew2.root", "TMVAResultsFinal/TMTMVAResultsNew4.root", "TMVAResultsFinal/TMTMVAResultsOld.root", "TMVAResultsFinal/TMTMVAResultsNew4NoAlpha.root"}; std::vector<std::string> myLabels = {"New2", "New4", "Old", "New4 no alpha"};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMTMVAResultsResNew4Alpha.root", "TMVAResultsFinal/TMTMVAResultsResNew4NoAlpha.root", "TMVAResultsFinal/TMTMVAResultsResOld.root"};
//std::vector<std::string> myLabels = {"Alpha", "NoAlpha", "Old"};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMTMVAResultsFullNew2.root", "TMVAResultsFinal/TMTMVAResultsFullNew4.root", "TMVAResultsFinal/TMTMVAResultsFullOld.root"}; std::vector<std::string> myLabels = {"New2", "New4", "Old"};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMVAResultsStd/TMVAResults_reco.root"};


/*std::vector<std::string> myFiles = {"TMTMVAResults_C.root", "TMTMVAResults_C_Reco.root", "TMTMVAResults_C_WithWidth.root", "TMTMVAResults_C_Iota.root", "TMTMVAResults_C_NoAlpha.root", "TMTMVAResults_C_WithWidthIota.root", "TMTMVAResults_C_WithWidthIota_Inclusive.root"};
std::vector<std::string> myLabels = {"C", "C Reco", "C wWidth", "C Iota", "C NoAlpha", "C wWidthIota", "C wWidthIota Inclusive"};*/


TString GetLabel(std::string label){
    
    // replace eta, delta, alpha, iota, fitscore
    if (label.find("Eta") != std::string::npos) label.replace(label.find("Eta"), 3, ": #eta");
    if (label.find("Delta") != std::string::npos) label.replace(label.find("Delta"), 5, ": #Delta");
    if (label.find("Alpha") != std::string::npos) label.replace(label.find("Alpha"), 5, ": #alpha");
    if (label.find("Iota") != std::string::npos) label.replace(label.find("Iota"), 4, ": #iota");
    if (label.find("FitScore") != std::string::npos) label.replace(label.find("FitScore"), 8, ": FitScore");

    return label.c_str();
}


void CompareROCCurves(std::string path="./", std::string keyLabel="") {

    
    // map with label and file name
    std::map<std::string, std::string> fileMap;

    TSystemDirectory dir(path.c_str(), path.c_str());
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            std::cout<<file->GetName()<<std::endl;
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith("TMVAResults")) {
                // check fname contains keyLabel
                if (keyLabel != "" && !fname.Contains(keyLabel.c_str())) continue;
                std:string label = fname.Data();
                label.erase(0, 11);
                label.erase(label.size()-5, label.size());
                fileMap[label] = fname.Data();
            }
        }
    }


    TCanvas* c = new TCanvas("c", "ROC Comparison", 800, 800);
    c->Divide(1,1);

    c->cd(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);

    std::map<std::string, TH1F*> rocs;
    // loop over the file map
    for (auto it = fileMap.begin(); it != fileMap.end(); ++it) {

        std::cout<<it->first<<" "<<it->second<<std::endl;

        TFile* file1 = TFile::Open((path + it->second).c_str());
        std::string fPlotName;
        TH1F* roc;
        
        // check if the name contains Fisher
        if (it->first.find("Fisher") != std::string::npos) {
            std::cout<<"Fisher"<<std::endl;
            fPlotName = "dataset/Method_Fisher/Fisher/MVA_Fisher_rejBvsS";
            roc = (TH1F*)file1->Get(fPlotName.c_str());
            rocs[it->first] = roc;
        }
        else if(it->first.find("Cuts") != std::string::npos) {
            std::cout<<"Cuts"<<std::endl;
            fPlotName = "dataset/Method_Cuts/Cuts/MVA_Cuts_rejBvsS";
            roc = (TH1F*)file1->Get(fPlotName.c_str());
            rocs[it->first] = roc;
        }
        else if(it->first.find("BDT") != std::string::npos) {
            std::cout<<"BDT"<<std::endl;
            fPlotName = "dataset/Method_BDT/BDT/MVA_BDT_rejBvsS";
            roc = (TH1F*)file1->Get(fPlotName.c_str());
            rocs[it->first] = roc;
        }
    }

    TH2F hFrame("hFrame", ";Signal efficiency; Background rejection", 500, 0, 1, 500, 0, 1);
    hFrame.SetStats(0);
    hFrame.Draw();
    hFrame.GetXaxis()->SetTitleOffset(1.2);
    hFrame.GetYaxis()->SetTitleOffset(1.2);

    TLegend* legend = new TLegend(0.2, 0.2, 0.5, 0.4);

    int rocCounter = 0;
    for (auto it = rocs.begin(); it != rocs.end(); ++it) {

        std::cout<<it->first<<" "<<it->second->Integral()<<std::endl;

        rocs[it->first]->SetLineColor(colors.at(rocCounter));
        rocs[it->first]->SetLineStyle(lines.at(rocCounter));
        rocs[it->first]->SetLineWidth(3);
        rocs[it->first]->Draw("L SAME");        

        std::cout<<it->first<<" "<<it->second->Integral()<<std::endl;
        std::cout<<GetLabel(it->first)<<std::endl;
        legend->AddEntry(it->second, GetLabel(it->first), "l");
        legend->Draw();
        rocCounter++;
    }

    c->cd();
    c->Update();
    c->WaitPrimitive();
    c->SaveAs("ROC_Comparison.pdf");
}