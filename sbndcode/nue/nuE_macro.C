#include <vector>
#include <map>
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

typedef struct{
    double event;
    double run;
    double subrun;
    double trueNeutrinoVX;
    double trueNeutrinoVY;
    double trueNeutrinoVZ;
    double recoNeutrinoVX;
    double recoNeutrinoVY;
    double recoNeutrinoVZ;
    int trueCCNC;
    int trueNeutrinoType;
    int trueLeptonType;
    int nSlices;
    int nPFPs;
    int nTracks;
    int nShowers;
    double tpcID;
    int dl_current;
    double sliceCompleteness;
    double trackScore;
    double showerEnergy;
    double showerTheta;
    double showerSmallestTheta;
    double showerETheta2;
} event_t;

void shower(std::vector<event_t> allEvents_vec){
    std::vector<event_t> dlDune = std::vector<event_t>();
    std::vector<event_t> dlUboone = std::vector<event_t>();
    std::vector<event_t> current = std::vector<event_t>();

    for(UInt_t i = 0; i < allEvents_vec.size(); i++){
        event_t event = allEvents_vec.at(i);
        if(event.dl_current == 0){
            dlUboone.push_back(event);
        } else if(event.dl_current == 1){
            dlDune.push_back(event);
        } else if(event.dl_current == 2){
            current.push_back(event);
        }
    }

    TCanvas *showerECanvas = new TCanvas("showerE_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* showerECurrent_dist = new TH1F("Shower Energy", "Energy of Shower", 60, 0, 6000);
    showerECurrent_dist->SetTitle("Energy of the Recoil Electron Shower in an Event;Energy of Shower (MeV);# of Events");
    TH1F* showerEDLUboone_dist = (TH1F*) showerECurrent_dist->Clone("Shower Energy");
    TH1F* showerEDLDune_dist = (TH1F*) showerECurrent_dist->Clone("Shower Energy");

    TCanvas *showerThetaCanvas = new TCanvas("showerTheta_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* showerThetaCurrent_dist = new TH1F("Shower Theta", "Theta of Shower", 90, 0, 360);
    showerThetaCurrent_dist->SetTitle("Angle of the Recoil Electron Shower Relative to Beam Neutrino Direction in an Event;Theta of Shower (degrees);# of Events");
    TH1F* showerThetaDLUboone_dist = (TH1F*) showerThetaCurrent_dist->Clone("Shower Theta");
    TH1F* showerThetaDLDune_dist = (TH1F*) showerThetaCurrent_dist->Clone("Shower Theta");

    TCanvas *showerSmallThetaCanvas = new TCanvas("showerSmallTheta_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* showerSmallThetaCurrent_dist = new TH1F("Shower Small Theta", "Small Theta of Shower", 90, 0, 360);
    showerSmallThetaCurrent_dist->SetTitle("Recoil Electron Shower with the Smallest Theta Relative to Beam Neutrino Direction in an Event;Smallest Shower Theta(degrees);# of Events");
    TH1F* showerSmallThetaDLUboone_dist = (TH1F*) showerSmallThetaCurrent_dist->Clone("Small Shower Theta");
    TH1F* showerSmallThetaDLDune_dist = (TH1F*) showerSmallThetaCurrent_dist->Clone("Small Shower Theta");

    TCanvas *showerETheta2Canvas = new TCanvas("showerETheta2_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* showerETheta2Current_dist = new TH1F("Shower ETheta2", "ETheta2 of Shower", 500, 0, 1000);
    showerETheta2Current_dist->SetTitle("E#theta^{2} of the Recoil Electron Shower in an Event;E#theta^{2} of Shower (MeV);# of Events");
    TH1F* showerETheta2DLUboone_dist = (TH1F*) showerETheta2Current_dist->Clone("Shower ETheta2");
    TH1F* showerETheta2DLDune_dist = (TH1F*) showerETheta2Current_dist->Clone("Shower ETheta2");

    for(UInt_t j = 0; j < current.size(); j++){
        showerECurrent_dist->Fill(current.at(j).showerEnergy);
        showerThetaCurrent_dist->Fill(TMath::RadToDeg() * current.at(j).showerTheta);
        showerSmallThetaCurrent_dist->Fill(TMath::RadToDeg() * current.at(j).showerSmallestTheta);
        showerETheta2Current_dist->Fill(current.at(j).showerETheta2);
    }

    for(UInt_t j = 0; j < dlUboone.size(); j++){
        showerEDLUboone_dist->Fill(dlUboone.at(j).showerEnergy);
        showerThetaDLUboone_dist->Fill(TMath::RadToDeg() * dlUboone.at(j).showerTheta);
        showerSmallThetaDLUboone_dist->Fill(TMath::RadToDeg() * dlUboone.at(j).showerSmallestTheta);
        showerETheta2DLUboone_dist->Fill(dlUboone.at(j).showerETheta2);
    }

    for(UInt_t j = 0; j < dlDune.size(); j++){
        showerEDLDune_dist->Fill(dlDune.at(j).showerEnergy);
        showerThetaDLDune_dist->Fill(TMath::RadToDeg() * dlDune.at(j).showerTheta);
        showerSmallThetaDLDune_dist->Fill(TMath::RadToDeg() * dlDune.at(j).showerSmallestTheta);
        showerETheta2DLDune_dist->Fill(dlDune.at(j).showerETheta2);
    }

    showerECanvas->cd();
    showerECurrent_dist->SetLineWidth(2);
    showerECurrent_dist->SetLineColor(kRed);
    showerEDLDune_dist->SetLineWidth(2);
    showerEDLDune_dist->SetLineColor(kViolet-5);
    showerEDLUboone_dist->SetLineWidth(2);
    showerEDLUboone_dist->SetLineColor(kBlue);
    showerECurrent_dist->Draw("hist");
    showerEDLDune_dist->Draw("histsame");
    showerEDLUboone_dist->Draw("histsame");
    showerECurrent_dist->SetStats(0);
    
    auto legend = new TLegend(1-0.86,0.86,1-0.54,0.70);
    legend->AddEntry(showerEDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(showerEDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(showerECurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();
    showerECanvas->SaveAs("/nashome/c/coackley/nuEPlots/showerE_dist.pdf");

    showerThetaCanvas->cd();
    showerThetaCurrent_dist->SetLineWidth(2);
    showerThetaCurrent_dist->SetLineColor(kRed);
    showerThetaDLDune_dist->SetLineWidth(2);
    showerThetaDLDune_dist->SetLineColor(kViolet-5);
    showerThetaDLUboone_dist->SetLineWidth(2);
    showerThetaDLUboone_dist->SetLineColor(kBlue);
    showerThetaCurrent_dist->Draw("hist");
    showerThetaDLDune_dist->Draw("histsame");
    showerThetaDLUboone_dist->Draw("histsame");
    showerThetaCurrent_dist->SetStats(0);

    auto legend1 = new TLegend(1-0.86,0.86,1-0.54,0.70);
    legend1->AddEntry(showerThetaDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend1->AddEntry(showerThetaDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend1->AddEntry(showerThetaCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend1->SetTextSize(0.0225);
    legend1->SetMargin(0.13);
    legend1->Draw();
    showerThetaCanvas->SaveAs("/nashome/c/coackley/nuEPlots/showerTheta_dist.pdf");

    showerSmallThetaCanvas->cd();
    showerSmallThetaCurrent_dist->SetLineWidth(2);
    showerSmallThetaCurrent_dist->SetLineColor(kRed);
    showerSmallThetaDLDune_dist->SetLineWidth(2);
    showerSmallThetaDLDune_dist->SetLineColor(kViolet-5);
    showerSmallThetaDLUboone_dist->SetLineWidth(2);
    showerSmallThetaDLUboone_dist->SetLineColor(kBlue);
    showerSmallThetaCurrent_dist->Draw("hist");
    showerSmallThetaDLDune_dist->Draw("histsame");
    showerSmallThetaDLUboone_dist->Draw("histsame");
    showerSmallThetaCurrent_dist->SetStats(0);

    auto legend2 = new TLegend(1-0.86,0.86,1-0.54,0.70);
    legend2->AddEntry(showerSmallThetaDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend2->AddEntry(showerSmallThetaDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend2->AddEntry(showerSmallThetaCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend2->SetTextSize(0.0225);
    legend2->SetMargin(0.13);
    legend2->Draw();
    showerSmallThetaCanvas->SaveAs("/nashome/c/coackley/nuEPlots/showerSmallTheta_dist.pdf");

    showerETheta2Canvas->cd();
    showerETheta2Current_dist->SetLineWidth(2);
    showerETheta2Current_dist->SetLineColor(kRed);
    showerETheta2DLDune_dist->SetLineWidth(2);
    showerETheta2DLDune_dist->SetLineColor(kViolet-5);
    showerETheta2DLUboone_dist->SetLineWidth(2);
    showerETheta2DLUboone_dist->SetLineColor(kBlue);
    showerETheta2Current_dist->Draw("hist");
    showerETheta2DLDune_dist->Draw("histsame");
    showerETheta2DLUboone_dist->Draw("histsame");
    showerETheta2Current_dist->SetStats(0);

    auto legend3 = new TLegend(1-0.86,0.86,1-0.54,0.70);
    legend3->AddEntry(showerETheta2DLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3->AddEntry(showerETheta2DLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3->AddEntry(showerETheta2Current_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend3->SetTextSize(0.0225);
    legend3->SetMargin(0.13);
    legend3->Draw();
    showerETheta2Canvas->SaveAs("/nashome/c/coackley/nuEPlots/showerETheta2_dist.pdf");
}


void scores(std::vector<event_t> allEvents_vec){
    size_t numEvents = allEvents_vec.size();
    std::vector<event_t> dlDune = std::vector<event_t>();
    std::vector<event_t> dlUboone = std::vector<event_t>();
    std::vector<event_t> current = std::vector<event_t>();

    for(UInt_t i = 0; i < numEvents; i++){
        event_t event = allEvents_vec.at(i);
        if(event.dl_current == 0){
            dlUboone.push_back(event);
        } else if(event.dl_current == 1){
            dlDune.push_back(event);
        } else if(event.dl_current == 2){
            current.push_back(event);
        }
    }

    TCanvas *sliceCompletenessCanvas = new TCanvas("sliceCompleteness_canvas", "Graph Draw Options", 200, 10, 600, 400);
    //TH1F* sliceCompletenessCurrent_dist = new TH1F("Slice Completeness", "Completeness of Slice", 35, 0.655, 1.005);
    TH1F* sliceCompletenessCurrent_dist = new TH1F("Slice Completeness", "Completeness of Slice", 35, 0.935, 1.005);
    sliceCompletenessCurrent_dist->SetTitle("Completeness of Slice in an Event;Completeness Score;# of Events");
    TH1F* sliceCompletenessDLUboone_dist = (TH1F*) sliceCompletenessCurrent_dist->Clone("Slice Completeness");
    TH1F* sliceCompletenessDLDune_dist = (TH1F*) sliceCompletenessCurrent_dist->Clone("Slice Completeness");

    TCanvas *trackScoreCanvas = new TCanvas("trackscore_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* trackScoreCurrent_dist = new TH1F("Track Score", "Score of Track", 21, -0.025, 1.025);
    trackScoreCurrent_dist->SetTitle("Recoil Electron Track Score in an Event (with Number of PFPs = 2);Track Score;# of Events");
    TH1F* trackScoreDLUboone_dist = (TH1F*) trackScoreCurrent_dist->Clone("Track Score");
    TH1F* trackScoreDLDune_dist = (TH1F*) trackScoreCurrent_dist->Clone("Track Score");
   
    double sizeCurrentSlice = 0;
    double sizeCurrentTrack = 0;
    double sizeDLDuneSlice = 0;
    double sizeDLDuneTrack = 0;
    double sizeDLUbooneSlice = 0;
    double sizeDLUbooneTrack = 0;

    std::cout << "current size: " << current.size() << std::endl; 
    for(UInt_t j = 0; j < current.size(); j++){
        //if(current.at(j).nSlices == 1){
            sliceCompletenessCurrent_dist->Fill(current.at(j).sliceCompleteness);
            sizeCurrentSlice++;
        //}

        if(current.at(j).nPFPs == 2){
            trackScoreCurrent_dist->Fill(current.at(j).trackScore);
            sizeCurrentTrack++;
        }
    }

    std::cout << "dl uboone: " << dlUboone.size() << std::endl;
    for(UInt_t j = 0; j < dlUboone.size(); j++){
        //if(dlUboone.at(j).nSlices == 1){ 
            sliceCompletenessDLUboone_dist->Fill(dlUboone.at(j).sliceCompleteness);
            sizeDLUbooneSlice++;
        //}

        if(dlUboone.at(j).nPFPs == 2){
            trackScoreDLUboone_dist->Fill(dlUboone.at(j).trackScore);
            sizeDLUbooneTrack++;
        }
    }

    std::cout << "dl dune size: " << dlDune.size() << std::endl;
    for(UInt_t j = 0; j < dlDune.size(); j++){
        //if(dlDune.at(j).nSlices == 1){
            sliceCompletenessDLDune_dist->Fill(dlDune.at(j).sliceCompleteness);
            sizeDLDuneSlice++;
        //}

        if(dlDune.at(j).nPFPs == 2){
            trackScoreDLDune_dist->Fill(dlDune.at(j).trackScore);
            sizeDLDuneTrack++;
        }
    }

    sliceCompletenessCanvas->cd();
    sliceCompletenessCurrent_dist->SetLineWidth(2);
    sliceCompletenessCurrent_dist->SetLineColor(kRed);
    sliceCompletenessDLDune_dist->SetLineWidth(2);
    sliceCompletenessDLDune_dist->SetLineColor(kViolet-5);
    sliceCompletenessDLUboone_dist->SetLineWidth(2);
    sliceCompletenessDLUboone_dist->SetLineColor(kBlue);
    sliceCompletenessCurrent_dist->Draw("hist");
    sliceCompletenessDLDune_dist->Draw("histsame");
    sliceCompletenessDLUboone_dist->Draw("histsame");
    sliceCompletenessCurrent_dist->SetStats(0);
    sliceCompletenessCurrent_dist->GetYaxis()->SetRangeUser(0, 800);
    //sliceCompletenessCurrent_dist->GetXaxis()->SetRangeUser(0.8, 1.0052);
    sliceCompletenessCurrent_dist->GetXaxis()->SetRangeUser(0.935, 1.001);

    auto legend = new TLegend(1-0.86,0.86,1-0.54,0.70);
    legend->AddEntry(sliceCompletenessDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(sliceCompletenessDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(sliceCompletenessCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();
    //sliceCompletenessCanvas->SaveAs("/nashome/c/coackley/nuEPlots/sliceCompleteness_dist.pdf");
    sliceCompletenessCanvas->SaveAs("/nashome/c/coackley/nuEPlots/sliceCompleteness_dist_zoomed.pdf");

    trackScoreCanvas->cd();
    trackScoreCurrent_dist->SetLineWidth(2);
    trackScoreCurrent_dist->SetLineColor(kRed);
    trackScoreDLDune_dist->SetLineWidth(2);
    trackScoreDLDune_dist->SetLineColor(kViolet-5);
    trackScoreDLUboone_dist->SetLineWidth(2);
    trackScoreDLUboone_dist->SetLineColor(kBlue);
    trackScoreCurrent_dist->Draw("hist");
    trackScoreDLDune_dist->Draw("histsame");
    trackScoreDLUboone_dist->Draw("histsame");
    trackScoreCurrent_dist->SetStats(0);
    trackScoreCurrent_dist->GetYaxis()->SetRangeUser(0, 200);
    
    auto legend1 = new TLegend(0.56,0.86,0.88,0.70);
    legend1->AddEntry(trackScoreDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend1->AddEntry(trackScoreDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend1->AddEntry(trackScoreCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend1->SetTextSize(0.0225);
    legend1->SetMargin(0.13);
    legend1->Draw();
    trackScoreCanvas->SaveAs("/nashome/c/coackley/nuEPlots/trackScore_dist.pdf");

    TCanvas *sliceCompletenessPercCanvas = new TCanvas("sliceCompletenessPerc_canvas", "Graph Draw Options", 200, 10, 600, 400);
    //TH1F* sliceCompletenessCurrent_dist_perc = new TH1F("Slice Completeness Perc", "Completeness of Slice Perc", 35, 0.655, 1.005);
    TH1F* sliceCompletenessCurrent_dist_perc = new TH1F("Slice Completeness Perc", "Completeness of Slice Perc", 35, 0.935, 1.005);
    sliceCompletenessCurrent_dist_perc->SetTitle("Completeness of Slice in an Event;Completeness Score;Percentage of Events (%)");
    TH1F* sliceCompletenessDLUboone_dist_perc = (TH1F*) sliceCompletenessCurrent_dist_perc->Clone("Slice Completeness Perc");
    TH1F* sliceCompletenessDLDune_dist_perc = (TH1F*) sliceCompletenessCurrent_dist_perc->Clone("Slice Completeness Perc");

    TCanvas *trackScorePercCanvas = new TCanvas("trackscorePerc_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* trackScoreCurrent_dist_perc = new TH1F("Track Score Perc", "Score of Track Perc", 21, -0.025, 1.025);
    trackScoreCurrent_dist_perc->SetTitle("Recoil Electron Track Score in an Event (with Number of PFPs = 2);Track Score;Percentage of Events (%)");
    TH1F* trackScoreDLUboone_dist_perc = (TH1F*) trackScoreCurrent_dist_perc->Clone("Track Score Perc");
    TH1F* trackScoreDLDune_dist_perc = (TH1F*) trackScoreCurrent_dist_perc->Clone("Track Score Perc");

    //double binFilling = 0.66; // lowest value + bin width/2 (0.655 + 0.005)
    double binFilling = 0.936; // lowest value + bin width/2 (0.935 + 0.001)
    for(int i = 1; i < 36; i++){
        //std::cout << "bin filling: " << binFilling << std::endl;
        sliceCompletenessCurrent_dist_perc->Fill(binFilling, (100.0f*sliceCompletenessCurrent_dist->GetBinContent(i)/sizeCurrentSlice));
        sliceCompletenessDLUboone_dist_perc->Fill(binFilling, (100.0f*sliceCompletenessDLUboone_dist->GetBinContent(i)/sizeDLUbooneSlice));
        sliceCompletenessDLDune_dist_perc->Fill(binFilling, (100.0f*sliceCompletenessDLDune_dist->GetBinContent(i)/sizeDLDuneSlice));
        //binFilling += 0.01;
        binFilling += 0.002; // bin width
        //std::cout << sliceCompletenessCurrent_dist->GetBinContent(i) << ", " << sliceCompletenessDLUboone_dist->GetBinContent(i) << ", " << sliceCompletenessDLDune_dist->GetBinContent(i) << std::endl;
    }

    binFilling = 0; // lowest value + bin width/2 (-0.025 + 0.025)
    for(int i = 1; i < 22; i++){
        trackScoreCurrent_dist_perc->Fill(binFilling, (100.0f*trackScoreCurrent_dist->GetBinContent(i)/sizeCurrentTrack));
        trackScoreDLUboone_dist_perc->Fill(binFilling, (100.0f*trackScoreDLUboone_dist->GetBinContent(i)/sizeDLUbooneTrack));
        trackScoreDLDune_dist_perc->Fill(binFilling, (100.0f*trackScoreDLDune_dist->GetBinContent(i)/sizeDLDuneTrack));
        binFilling += 0.05; // bin width
    }

    sliceCompletenessPercCanvas->cd();
    sliceCompletenessCurrent_dist_perc->SetLineWidth(2);
    sliceCompletenessCurrent_dist_perc->SetLineColor(kRed);
    sliceCompletenessDLDune_dist_perc->SetLineWidth(2);
    sliceCompletenessDLDune_dist_perc->SetLineColor(kViolet-5);
    sliceCompletenessDLUboone_dist_perc->SetLineWidth(2);
    sliceCompletenessDLUboone_dist_perc->SetLineColor(kBlue);
    sliceCompletenessCurrent_dist_perc->Draw("hist");
    sliceCompletenessDLDune_dist_perc->Draw("histsame");
    sliceCompletenessDLUboone_dist_perc->Draw("histsame");
    sliceCompletenessCurrent_dist_perc->SetStats(0);
    sliceCompletenessCurrent_dist_perc->GetYaxis()->SetRangeUser(0, 80);
    //sliceCompletenessCurrent_dist_perc->GetXaxis()->SetRangeUser(0.8, 1.0052);
    sliceCompletenessCurrent_dist_perc->GetXaxis()->SetRangeUser(0.935, 1.001);

    auto legend3 = new TLegend(1-0.86,0.86,1-0.54,0.70);
    legend3->AddEntry(sliceCompletenessDLDune_dist_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3->AddEntry(sliceCompletenessDLUboone_dist_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3->AddEntry(sliceCompletenessCurrent_dist_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend3->SetTextSize(0.0225);
    legend3->SetMargin(0.13);
    legend3->Draw();
    //sliceCompletenessPercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/sliceCompleteness_dist_perc.pdf");
    sliceCompletenessPercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/sliceCompleteness_dist_perc_zoomed.pdf");

    trackScorePercCanvas->cd();
    trackScoreCurrent_dist_perc->SetLineWidth(2);
    trackScoreCurrent_dist_perc->SetLineColor(kRed);
    trackScoreDLDune_dist_perc->SetLineWidth(2);
    trackScoreDLDune_dist_perc->SetLineColor(kViolet-5);
    trackScoreDLUboone_dist_perc->SetLineWidth(2);
    trackScoreDLUboone_dist_perc->SetLineColor(kBlue);
    trackScoreCurrent_dist_perc->Draw("hist");
    trackScoreDLDune_dist_perc->Draw("histsame");
    trackScoreDLUboone_dist_perc->Draw("histsame");
    trackScoreCurrent_dist_perc->SetStats(0);
    trackScoreCurrent_dist_perc->GetYaxis()->SetRangeUser(0, 30);
    
    auto legend4 = new TLegend(0.56,0.86,0.88,0.70);
    legend4->AddEntry(trackScoreDLDune_dist_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend4->AddEntry(trackScoreDLUboone_dist_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend4->AddEntry(trackScoreCurrent_dist_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend4->SetTextSize(0.0225);
    legend4->SetMargin(0.13);
    legend4->Draw();
    trackScorePercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/trackScore_dist_perc.pdf");
}

void numberPlots(std::vector<event_t> allEvents_vec){
    size_t numEvents = allEvents_vec.size();

    std::vector<event_t> dlDune = std::vector<event_t>();
    std::vector<event_t> dlUboone = std::vector<event_t>();
    std::vector<event_t> current = std::vector<event_t>();

    for(UInt_t i = 0; i < numEvents; i++){
       event_t event = allEvents_vec.at(i);
       if(event.dl_current == 0){
           dlUboone.push_back(event);
       } else if(event.dl_current == 1){
           dlDune.push_back(event);
       } else if(event.dl_current == 2){
           current.push_back(event);
       }
    }

    TCanvas *numSlicesCanvas = new TCanvas("numSlices_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numSlicesCurrent_dist = new TH1F("Number of Slices", "Num Slices", 5, -0.5, 4.5);
    numSlicesCurrent_dist->SetTitle("Distribution of the Number of Slices in an Event;Number of Slices in Event;# of Events");
    TH1F* numSlicesDLUboone_dist = (TH1F*) numSlicesCurrent_dist->Clone("Number of Slices");
    TH1F* numSlicesDLDune_dist = (TH1F*) numSlicesCurrent_dist->Clone("Number of Slices");

    TCanvas *numPFPsCanvas = new TCanvas("numPFPs_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numPFPsCurrent_dist = new TH1F("Number of PFPs", "Num PFPs", 10, -0.5, 9.5);
    numPFPsCurrent_dist->SetTitle("Distribution of the Number of PFPs in an Event;Number of PFPs in Event;# of Events");
    TH1F* numPFPsDLUboone_dist = (TH1F*) numPFPsCurrent_dist->Clone("Number of PFPs");
    TH1F* numPFPsDLDune_dist = (TH1F*) numPFPsCurrent_dist->Clone("Number of PFPs");

    TCanvas *numTracksCanvas = new TCanvas("numTracks_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numTracksCurrent_dist = new TH1F("Number of Tracks", "Num Tracks", 10, -0.5, 9.5);
    numTracksCurrent_dist->SetTitle("Distribution of the Number of Tracks in an Event;Number of Tracks in Event;# of Events");
    TH1F* numTracksDLUboone_dist = (TH1F*) numTracksCurrent_dist->Clone("Number of Tracks");
    TH1F* numTracksDLDune_dist = (TH1F*) numTracksCurrent_dist->Clone("Number of Tracks");

    TCanvas *numShowersCanvas = new TCanvas("numShowers_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numShowersCurrent_dist = new TH1F("Number of Showers", "Num Showers", 10, -0.5, 9.5);
    numShowersCurrent_dist->SetTitle("Distribution of the Number of Showers in an Event;Number of Showers in Event;# of Events");
    TH1F* numShowersDLUboone_dist = (TH1F*) numShowersCurrent_dist->Clone("Number of Showers");
    TH1F* numShowersDLDune_dist = (TH1F*) numShowersCurrent_dist->Clone("Number of Showers");

    for(UInt_t j = 0; j < current.size(); j++){
        numSlicesCurrent_dist->Fill(current.at(j).nSlices);
        numPFPsCurrent_dist->Fill(current.at(j).nPFPs);
        numTracksCurrent_dist->Fill(current.at(j).nTracks);
        numShowersCurrent_dist->Fill(current.at(j).nShowers);
    }

    for(UInt_t j = 0; j < dlDune.size(); j++){
        numSlicesDLDune_dist->Fill(dlDune.at(j).nSlices);
        numPFPsDLDune_dist->Fill(dlDune.at(j).nPFPs);
        numTracksDLDune_dist->Fill(dlDune.at(j).nTracks);
        numShowersDLDune_dist->Fill(dlDune.at(j).nShowers);
    }

    for(UInt_t j = 0; j < dlUboone.size(); j++){
        numSlicesDLUboone_dist->Fill(dlUboone.at(j).nSlices);
        numPFPsDLUboone_dist->Fill(dlUboone.at(j).nPFPs);
        numTracksDLUboone_dist->Fill(dlUboone.at(j).nTracks);
        numShowersDLUboone_dist->Fill(dlUboone.at(j).nShowers);
    }

    numSlicesCanvas->cd();
    numSlicesCurrent_dist->SetLineWidth(2);
    numSlicesCurrent_dist->SetLineColor(kRed);
    numSlicesDLDune_dist->SetLineWidth(2);
    numSlicesDLDune_dist->SetLineColor(kViolet-5);
    numSlicesDLUboone_dist->SetLineWidth(2);
    numSlicesDLUboone_dist->SetLineColor(kBlue);
    numSlicesCurrent_dist->Draw("hist");
    numSlicesDLDune_dist->Draw("histsame");
    numSlicesDLUboone_dist->Draw("histsame");
    numSlicesCurrent_dist->SetStats(0);
    numSlicesCurrent_dist->GetYaxis()->SetRangeUser(0, 1000);
    
    auto legend = new TLegend(0.56,0.86,0.88,0.70);
    legend->AddEntry(numSlicesDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(numSlicesDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(numSlicesCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();
    numSlicesCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numSlices_dist.pdf");
    
    numPFPsCanvas->cd();
    numPFPsCurrent_dist->SetLineWidth(2);
    numPFPsCurrent_dist->SetLineColor(kRed);
    numPFPsDLDune_dist->SetLineWidth(2);
    numPFPsDLDune_dist->SetLineColor(kViolet-5);
    numPFPsDLUboone_dist->SetLineWidth(2);
    numPFPsDLUboone_dist->SetLineColor(kBlue);
    numPFPsCurrent_dist->Draw("hist");
    numPFPsDLDune_dist->Draw("histsame");
    numPFPsDLUboone_dist->Draw("histsame");
    numPFPsCurrent_dist->SetStats(0);
    numPFPsCurrent_dist->GetYaxis()->SetRangeUser(0, 900);

    auto legend2 = new TLegend(0.56,0.86,0.88,0.70);
    legend2->AddEntry(numPFPsDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend2->AddEntry(numPFPsDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend2->AddEntry(numPFPsCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend2->SetTextSize(0.0225);
    legend2->SetMargin(0.13);
    legend2->Draw();
    numPFPsCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numPFPs_dist.pdf");

    numTracksCanvas->cd();
    numTracksCurrent_dist->SetLineWidth(2);
    numTracksCurrent_dist->SetLineColor(kRed);
    numTracksDLDune_dist->SetLineWidth(2);
    numTracksDLDune_dist->SetLineColor(kViolet-5);
    numTracksDLUboone_dist->SetLineWidth(2);
    numTracksDLUboone_dist->SetLineColor(kBlue);
    numTracksCurrent_dist->Draw("hist");
    numTracksDLDune_dist->Draw("histsame");
    numTracksDLUboone_dist->Draw("histsame");
    numTracksCurrent_dist->SetStats(0);
    numTracksCurrent_dist->GetYaxis()->SetRangeUser(0, 900);
    
    auto legend3 = new TLegend(0.56,0.86,0.88,0.70);
    legend3->AddEntry(numTracksDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend3->AddEntry(numTracksDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend3->AddEntry(numTracksCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend3->SetTextSize(0.0225);
    legend3->SetMargin(0.13);
    legend3->Draw();
    numTracksCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numTracks_dist.pdf");

    numShowersCanvas->cd();
    numShowersCurrent_dist->SetLineWidth(2);
    numShowersCurrent_dist->SetLineColor(kRed);
    numShowersDLDune_dist->SetLineWidth(2);
    numShowersDLDune_dist->SetLineColor(kViolet-5);
    numShowersDLUboone_dist->SetLineWidth(2);
    numShowersDLUboone_dist->SetLineColor(kBlue);
    numShowersCurrent_dist->Draw("hist");
    numShowersDLDune_dist->Draw("histsame");
    numShowersDLUboone_dist->Draw("histsame");
    numShowersCurrent_dist->SetStats(0);
    numShowersCurrent_dist->GetYaxis()->SetRangeUser(0, 900);

    auto legend4 = new TLegend(0.56,0.86,0.88,0.70);
    legend4->AddEntry(numShowersDLDune_dist, "Deep Learning: DUNE/LBNF Tune", "f");
    legend4->AddEntry(numShowersDLUboone_dist, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend4->AddEntry(numShowersCurrent_dist, "Current SBND Vertexing (without Refinement)", "f");
    legend4->SetTextSize(0.0225);
    legend4->SetMargin(0.13);
    legend4->Draw();
    numShowersCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numShowers_dist.pdf");

    TCanvas *numSlicesPercCanvas = new TCanvas("numSlicesPerc_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numSlicesCurrent_dist_perc = new TH1F("Number of Slices Perc", "Num Slices Perc", 5, -0.5, 4.5);
    numSlicesCurrent_dist_perc->SetTitle("Distribution of the Number of Slices in an Event;Number of Slices in Event;Percentage of Events (%)");
    TH1F* numSlicesDLUboone_dist_perc = (TH1F*) numSlicesCurrent_dist_perc->Clone("Number of Slices Perc");
    TH1F* numSlicesDLDune_dist_perc = (TH1F*) numSlicesCurrent_dist_perc->Clone("Number of Slices Perc");

    TCanvas *numPFPsPercCanvas = new TCanvas("numPFPsPerc_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numPFPsCurrent_dist_perc = new TH1F("Number of PFPs Perc", "Num PFPs Perc", 10, -0.5, 9.5);
    numPFPsCurrent_dist_perc->SetTitle("Distribution of the Number of PFPs in an Event;Number of PFPs in Event;Percentage of Events (%)");
    TH1F* numPFPsDLUboone_dist_perc = (TH1F*) numPFPsCurrent_dist_perc->Clone("Number of PFPs Perc");
    TH1F* numPFPsDLDune_dist_perc = (TH1F*) numPFPsCurrent_dist_perc->Clone("Number of PFPs Perc");    

    TCanvas *numTracksPercCanvas = new TCanvas("numTracksPerc_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numTracksCurrent_dist_perc = new TH1F("Number of Tracks Perc", "Num Tracks Perc", 10, -0.5, 9.5);
    numTracksCurrent_dist_perc->SetTitle("Distribution of the Number of Tracks in an Event;Number of Tracks in Event;Percentage of Events (%)");
    TH1F* numTracksDLUboone_dist_perc = (TH1F*) numTracksCurrent_dist_perc->Clone("Number of Tracks Perc");
    TH1F* numTracksDLDune_dist_perc = (TH1F*) numTracksCurrent_dist_perc->Clone("Number of Tracks Perc");

    TCanvas *numShowersPercCanvas = new TCanvas("numShowersPerc_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* numShowersCurrent_dist_perc = new TH1F("Number of Showers Perc", "Num Showers Perc", 10, -0.5, 9.5);
    numShowersCurrent_dist_perc->SetTitle("Distribution of the Number of Showers in an Event;Number of Showers in Event;Percentage of Events (%)");
    TH1F* numShowersDLUboone_dist_perc = (TH1F*) numShowersCurrent_dist_perc->Clone("Number of Showers Perc");
    TH1F* numShowersDLDune_dist_perc = (TH1F*) numShowersCurrent_dist_perc->Clone("Number of Showers Perc");

    double binFilling = 0; // lowest value + bin width/2 (-0.5 + 0.5)
    for(int i = 1; i < 11; i++){
        numPFPsCurrent_dist_perc->Fill(binFilling, (100.0f*numPFPsCurrent_dist->GetBinContent(i)/(double)current.size()));
        numPFPsDLUboone_dist_perc->Fill(binFilling, (100.0f*numPFPsDLUboone_dist->GetBinContent(i)/(double)dlUboone.size()));
        numPFPsDLDune_dist_perc->Fill(binFilling, (100.0f*numPFPsDLDune_dist->GetBinContent(i)/(double)dlDune.size()));
        numTracksCurrent_dist_perc->Fill(binFilling, (100.0f*numTracksCurrent_dist->GetBinContent(i)/(double)current.size()));
        numTracksDLUboone_dist_perc->Fill(binFilling, (100.0f*numTracksDLUboone_dist->GetBinContent(i)/(double)dlUboone.size()));
        numTracksDLDune_dist_perc->Fill(binFilling, (100.0f*numTracksDLDune_dist->GetBinContent(i)/(double)dlDune.size()));
        numShowersCurrent_dist_perc->Fill(binFilling, (100.0f*numShowersCurrent_dist->GetBinContent(i)/(double)current.size()));
        numShowersDLUboone_dist_perc->Fill(binFilling, (100.0f*numShowersDLUboone_dist->GetBinContent(i)/(double)dlUboone.size()));
        numShowersDLDune_dist_perc->Fill(binFilling, (100.0f*numShowersDLDune_dist->GetBinContent(i)/(double)dlDune.size()));
        binFilling += 1;
    }

    binFilling = 0; // lowest value + bin width/2 (-0.5 + 0.5)
    for(int i = 1; i < 6; i++){
        numSlicesCurrent_dist_perc->Fill(binFilling, (100.0f*numSlicesCurrent_dist->GetBinContent(i)/(double)current.size()));
        numSlicesDLUboone_dist_perc->Fill(binFilling, (100.0f*numSlicesDLUboone_dist->GetBinContent(i)/(double)dlUboone.size()));
        numSlicesDLDune_dist_perc->Fill(binFilling, (100.0f*numSlicesDLDune_dist->GetBinContent(i)/(double)dlDune.size()));
        binFilling += 1;
    }

    numSlicesPercCanvas->cd();
    numSlicesCurrent_dist_perc->SetLineWidth(2);
    numSlicesCurrent_dist_perc->SetLineColor(kRed);
    numSlicesDLDune_dist_perc->SetLineWidth(2);
    numSlicesDLDune_dist_perc->SetLineColor(kViolet-5);
    numSlicesDLUboone_dist_perc->SetLineWidth(2);
    numSlicesDLUboone_dist_perc->SetLineColor(kBlue);
    numSlicesCurrent_dist_perc->Draw("hist");
    numSlicesDLDune_dist_perc->Draw("histsame");
    numSlicesDLUboone_dist_perc->Draw("histsame");
    numSlicesCurrent_dist_perc->SetStats(0);
    numSlicesCurrent_dist_perc->GetYaxis()->SetRangeUser(0, 110);

    auto legend5 = new TLegend(0.56,0.86,0.88,0.70);
    legend5->AddEntry(numSlicesDLDune_dist_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend5->AddEntry(numSlicesDLUboone_dist_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend5->AddEntry(numSlicesCurrent_dist_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend5->SetTextSize(0.0225);
    legend5->SetMargin(0.13);
    legend5->Draw();
    numSlicesPercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numSlices_dist_percentage.pdf");

    numPFPsPercCanvas->cd();
    numPFPsCurrent_dist_perc->SetLineWidth(2);
    numPFPsCurrent_dist_perc->SetLineColor(kRed);
    numPFPsDLDune_dist_perc->SetLineWidth(2);
    numPFPsDLDune_dist_perc->SetLineColor(kViolet-5);
    numPFPsDLUboone_dist_perc->SetLineWidth(2);
    numPFPsDLUboone_dist_perc->SetLineColor(kBlue);
    numPFPsCurrent_dist_perc->Draw("hist");
    numPFPsDLDune_dist_perc->Draw("histsame");
    numPFPsDLUboone_dist_perc->Draw("histsame");
    numPFPsCurrent_dist_perc->SetStats(0);
    numPFPsCurrent_dist_perc->GetYaxis()->SetRangeUser(0, 80);

    auto legend6 = new TLegend(0.56,0.86,0.88,0.70);
    legend6->AddEntry(numPFPsDLDune_dist_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend6->AddEntry(numPFPsDLUboone_dist_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend6->AddEntry(numPFPsCurrent_dist_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend6->SetTextSize(0.0225);
    legend6->SetMargin(0.13);
    legend6->Draw();
    numPFPsPercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numPFPs_dist_percentage.pdf");

    numTracksPercCanvas->cd();
    numTracksCurrent_dist_perc->SetLineWidth(2);
    numTracksCurrent_dist_perc->SetLineColor(kRed);
    numTracksDLDune_dist_perc->SetLineWidth(2);
    numTracksDLDune_dist_perc->SetLineColor(kViolet-5);
    numTracksDLUboone_dist_perc->SetLineWidth(2);
    numTracksDLUboone_dist_perc->SetLineColor(kBlue);
    numTracksCurrent_dist_perc->Draw("hist");
    numTracksDLDune_dist_perc->Draw("histsame");
    numTracksDLUboone_dist_perc->Draw("histsame");
    numTracksCurrent_dist_perc->SetStats(0);
    numTracksCurrent_dist_perc->GetYaxis()->SetRangeUser(0, 80);
    
    auto legend7 = new TLegend(0.56,0.86,0.88,0.70);
    legend7->AddEntry(numTracksDLDune_dist_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend7->AddEntry(numTracksDLUboone_dist_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend7->AddEntry(numTracksCurrent_dist_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend7->SetTextSize(0.0225);
    legend7->SetMargin(0.13);
    legend7->Draw();
    numTracksPercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numTracks_dist_percentage.pdf");

    numShowersPercCanvas->cd();
    numShowersCurrent_dist_perc->SetLineWidth(2);
    numShowersCurrent_dist_perc->SetLineColor(kRed);
    numShowersDLDune_dist_perc->SetLineWidth(2);
    numShowersDLDune_dist_perc->SetLineColor(kViolet-5);
    numShowersDLUboone_dist_perc->SetLineWidth(2);
    numShowersDLUboone_dist_perc->SetLineColor(kBlue);
    numShowersCurrent_dist_perc->Draw("hist");
    numShowersDLDune_dist_perc->Draw("histsame");
    numShowersDLUboone_dist_perc->Draw("histsame");
    numShowersCurrent_dist_perc->SetStats(0);
    numShowersCurrent_dist_perc->GetYaxis()->SetRangeUser(0, 80);

    auto legend8 = new TLegend(0.56,0.86,0.88,0.70);
    legend8->AddEntry(numShowersDLDune_dist_perc, "Deep Learning: DUNE/LBNF Tune", "f");
    legend8->AddEntry(numShowersDLUboone_dist_perc, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend8->AddEntry(numShowersCurrent_dist_perc, "Current SBND Vertexing (without Refinement)", "f");
    legend8->SetTextSize(0.0225);
    legend8->SetMargin(0.13);
    legend8->Draw();
    numShowersPercCanvas->SaveAs("/nashome/c/coackley/nuEPlots/numShowers_dist_percentage.pdf");
}

void deltaVertices(std::vector<event_t> allEvents_vec){
    size_t numEvents = allEvents_vec.size();
    
    TCanvas *deltaXCanvas = new TCanvas("deltaX_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaX_dist = new TH1F("Delta X", "X Coord Diff", 40, -10, 10);
    deltaX_dist->SetTitle("Neutrino Vertex deltaX;True - Reco X Coordinate;#");

    TCanvas *deltaYCanvas = new TCanvas("deltaY_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaY_dist = new TH1F("Delta Y", "Y Coord Diff", 40, -10, 10);
    deltaY_dist->SetTitle("Neutrino Vertex deltaY;True - Reco Y Coordinate;#");

    TCanvas *deltaZCanvas = new TCanvas("deltaZ_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaZ_dist = new TH1F("Delta Z", "Z Coord Diff", 40, -10, 10);
    deltaZ_dist->SetTitle("Neutrino Vertex deltaZ;True - Reco Z Coordinate;#");

    TCanvas *deltaRCanvas = new TCanvas("deltaR_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaR_dist = new TH1F("Delta R", "R Coord Diff", 40, 0, 20);
    deltaR_dist->SetTitle("Neutrino Vertex deltaR;True - Reco R Coordinate;#");

    for(UInt_t i = 0; i < numEvents; i++){
        event_t event = allEvents_vec.at(i);
        deltaX_dist->Fill(event.trueNeutrinoVX - event.recoNeutrinoVX);
        deltaY_dist->Fill(event.trueNeutrinoVY - event.recoNeutrinoVY);
        deltaZ_dist->Fill(event.trueNeutrinoVZ - event.recoNeutrinoVZ);
        deltaR_dist->Fill(std::sqrt( std::pow(event.trueNeutrinoVX - event.recoNeutrinoVX, 2) + std::pow(event.trueNeutrinoVY - event.recoNeutrinoVY, 2) + std::pow(event.trueNeutrinoVZ - event.recoNeutrinoVZ, 2) ));
    }

    deltaXCanvas->cd();
    deltaX_dist->Draw();
    deltaX_dist->SetLineWidth(2);
    deltaX_dist->SetStats(0);
    deltaXCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaX_dist_nuE.pdf");

    deltaYCanvas->cd();
    deltaY_dist->Draw();
    deltaY_dist->SetLineWidth(2);
    deltaY_dist->SetStats(0);
    deltaYCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaY_dist_nuE.pdf");

    deltaZCanvas->cd();
    deltaZ_dist->Draw();
    deltaZ_dist->SetLineWidth(2);
    deltaZ_dist->SetStats(0);
    deltaZCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaZ_dist_nuE.pdf");

    deltaRCanvas->cd();
    deltaR_dist->Draw();
    deltaR_dist->SetLineWidth(2);
    deltaR_dist->SetStats(0);
    deltaRCanvas->SaveAs("/nashome/c/coackley/nuEPlots/deltaR_dist_nuE.pdf");
}

void nuE_macro(){
    TFile *f = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_03/NuEAnalyserOutput.root", "READ");
    //TFile *f = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/merged.root", "READ");
    
    if(!f){
        std::cout << "Failed to read file" << std::endl;
        return;
    }

    TTree *t;
    f->GetObject("NuETree", t);

    std::vector<double> allEvent = std::vector<double>(0);
    std::vector<double> allRun = std::vector<double>(0);
    std::vector<double> allSubrun = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVX = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVY = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVZ = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVX = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVY = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVZ = std::vector<double>(0);
    std::vector<int> allTrueCCNC = std::vector<int>(0);
    std::vector<int> allTrueNeutrinoType = std::vector<int>(0);
    std::vector<int> allTrueLeptonType = std::vector<int>(0);
    std::vector<int> allNumSlices = std::vector<int>(0);
    std::vector<int> allNumPfps = std::vector<int>(0);
    std::vector<int> allNumTracks = std::vector<int>(0);
    std::vector<int> allNumShowers = std::vector<int>(0);
    std::vector<double> allTpcID = std::vector<double>(0);
    std::vector<int> allDl_current = std::vector<int>(0);
    std::vector<double> allSliceCompleteness = std::vector<double>(0);
    std::vector<double> allTrackScore = std::vector<double>(0);
    std::vector<double> allShowerEnergy = std::vector<double>(0);
    std::vector<double> allShowerTheta = std::vector<double>(0);
    std::vector<double> allShowerSmallestTheta = std::vector<double>(0);
    std::vector<double> allShowerETheta2 = std::vector<double>(0);

    Long64_t numEntries = t->GetEntries();
    int fileNumber = 1;

    for(Long64_t i = 0; i < numEntries; i++){
        std::cout << "File number: " << fileNumber << std::endl;
        std::vector<double> *pEvent = 0;
        std::vector<double> *pRun = 0; 
        std::vector<double> *pSubrun = 0;
        std::vector<double> *pTrueNeutrinoVX = 0;
        std::vector<double> *pTrueNeutrinoVY = 0;
        std::vector<double> *pTrueNeutrinoVZ = 0;
        std::vector<double> *pRecoNeutrinoVX = 0;
        std::vector<double> *pRecoNeutrinoVY = 0;
        std::vector<double> *pRecoNeutrinoVZ = 0;
        std::vector<int> *pTrueCCNC = 0;
        std::vector<int> *pTrueNeutrinoType = 0;
        std::vector<int> *pTrueLeptonType = 0;
        std::vector<int> *pNumSlices = 0;
        std::vector<int> *pNumPfps = 0;
        std::vector<int> *pNumTracks = 0;
        std::vector<int> *pNumShowers = 0;
        std::vector<double> *pTpcID = 0;
        std::vector<int> *pDl_current = 0;
        std::vector<double> *pSliceCompleteness = 0;
        std::vector<double> *pTrackScore = 0;
        std::vector<double> *pShowerEnergy = 0;
        std::vector<double> *pShowerTheta = 0;
        std::vector<double> *pShowerSmallestTheta = 0;
        std::vector<double> *pShowerETheta2 = 0;

        TBranch *event_branch = 0;
        TBranch *run_branch = 0;
        TBranch *subrun_branch = 0;
        TBranch *trueNeutrinoVX_branch = 0;
        TBranch *trueNeutrinoVY_branch = 0;
        TBranch *trueNeutrinoVZ_branch = 0;
        TBranch *recoNeutrinoVX_branch = 0;
        TBranch *recoNeutrinoVY_branch = 0;
        TBranch *recoNeutrinoVZ_branch = 0;
        TBranch *trueCCNC_branch = 0;
        TBranch *trueNeutrinoType_branch = 0;
        TBranch *trueLeptonType_branch = 0;
        TBranch *numSlices_branch = 0;
        TBranch *numPfps_branch = 0;
        TBranch *numTracks_branch = 0;
        TBranch *numShowers_branch = 0;
        TBranch *tpcID_branch = 0;
        TBranch *dl_current_branch = 0;
        TBranch *sliceCompleteness_branch = 0;
        TBranch *trackScore_branch = 0;
        TBranch *showerEnergy_branch = 0;
        TBranch *showerTheta_branch = 0;
        TBranch *showerSmallestTheta_branch = 0;
        TBranch *showerETheta2_branch = 0;

        t->SetBranchAddress("event_tree", &pEvent, &event_branch);
        t->SetBranchAddress("run_tree", &pRun, &run_branch);
        t->SetBranchAddress("subrun_tree", &pSubrun, &subrun_branch);
        t->SetBranchAddress("trueNeutrinoVX_tree", &pTrueNeutrinoVX, &trueNeutrinoVX_branch);
        t->SetBranchAddress("trueNeutrinoVY_tree", &pTrueNeutrinoVY, &trueNeutrinoVY_branch);
        t->SetBranchAddress("trueNeutrinoVZ_tree", &pTrueNeutrinoVZ, &trueNeutrinoVZ_branch);
        t->SetBranchAddress("recoNeutrinoVX_tree", &pRecoNeutrinoVX, &recoNeutrinoVX_branch);
        t->SetBranchAddress("recoNeutrinoVY_tree", &pRecoNeutrinoVY, &recoNeutrinoVY_branch);
        t->SetBranchAddress("recoNeutrinoVZ_tree", &pRecoNeutrinoVZ, &recoNeutrinoVZ_branch);
        t->SetBranchAddress("trueCCNC_tree", &pTrueCCNC, &trueCCNC_branch);
        t->SetBranchAddress("trueNeutrinoType_tree", &pTrueNeutrinoType, &trueNeutrinoType_branch);
        t->SetBranchAddress("trueLeptonType_tree", &pTrueLeptonType, &trueLeptonType_branch);
        t->SetBranchAddress("numSlices_tree", &pNumSlices, &numSlices_branch);
        t->SetBranchAddress("numPfps_tree", &pNumPfps, &numPfps_branch);
        t->SetBranchAddress("numTracks_tree", &pNumTracks, &numTracks_branch);
        t->SetBranchAddress("numShowers_tree", &pNumShowers, &numShowers_branch);
        t->SetBranchAddress("tpcID_tree", &pTpcID, &tpcID_branch);
        t->SetBranchAddress("dl_current_tree", &pDl_current, &dl_current_branch);
        t->SetBranchAddress("sliceCompleteness_tree", &pSliceCompleteness, &sliceCompleteness_branch);
        t->SetBranchAddress("trackScore_tree", &pTrackScore, &trackScore_branch);
        t->SetBranchAddress("showerEnergy_tree", &pShowerEnergy, &showerEnergy_branch);
        t->SetBranchAddress("showerTheta_tree", &pShowerTheta, &showerTheta_branch);
        t->SetBranchAddress("showerSmallestTheta_tree", &pShowerSmallestTheta, &showerSmallestTheta_branch);
        t->SetBranchAddress("showerETheta2_tree", &pShowerETheta2, &showerETheta2_branch);

        event_branch->GetEntry(i);
        run_branch->GetEntry(i);
        subrun_branch->GetEntry(i);
        trueNeutrinoVX_branch->GetEntry(i);
        trueNeutrinoVY_branch->GetEntry(i);
        trueNeutrinoVZ_branch->GetEntry(i);
        recoNeutrinoVX_branch->GetEntry(i);
        recoNeutrinoVY_branch->GetEntry(i);
        recoNeutrinoVZ_branch->GetEntry(i);
        trueCCNC_branch->GetEntry(i);
        trueNeutrinoType_branch->GetEntry(i);
        trueLeptonType_branch->GetEntry(i);
        numSlices_branch->GetEntry(i);
        numPfps_branch->GetEntry(i);
        numTracks_branch->GetEntry(i);
        numShowers_branch->GetEntry(i);
        tpcID_branch->GetEntry(i);
        dl_current_branch->GetEntry(i);
        sliceCompleteness_branch->GetEntry(i);
        trackScore_branch->GetEntry(i);
        showerEnergy_branch->GetEntry(i);
        showerTheta_branch->GetEntry(i);
        showerSmallestTheta_branch->GetEntry(i);
        showerETheta2_branch->GetEntry(i);

        size_t vecSize = pEvent->size();
        for(UInt_t j = 0; j < vecSize; j++){
            allEvent.push_back(pEvent->at(j));
            allRun.push_back(pRun->at(j));
            allSubrun.push_back(pSubrun->at(j));
            allTrueNeutrinoVX.push_back(pTrueNeutrinoVX->at(j));
            allTrueNeutrinoVY.push_back(pTrueNeutrinoVY->at(j));
            allTrueNeutrinoVZ.push_back(pTrueNeutrinoVZ->at(j));
            allRecoNeutrinoVX.push_back(pRecoNeutrinoVX->at(j));
            allRecoNeutrinoVY.push_back(pRecoNeutrinoVY->at(j));
            allRecoNeutrinoVZ.push_back(pRecoNeutrinoVZ->at(j));
            allTrueCCNC.push_back(pTrueCCNC->at(j));
            allTrueNeutrinoType.push_back(pTrueNeutrinoType->at(j));
            allTrueLeptonType.push_back(pTrueLeptonType->at(j));
            allNumSlices.push_back(pNumSlices->at(j));
            allNumPfps.push_back(pNumPfps->at(j));
            allNumTracks.push_back(pNumTracks->at(j));
            allNumShowers.push_back(pNumShowers->at(j));
            allTpcID.push_back(pTpcID->at(j));
            allDl_current.push_back(pDl_current->at(j));
            allSliceCompleteness.push_back(pSliceCompleteness->at(j));
            allTrackScore.push_back(pTrackScore->at(j));
            allShowerEnergy.push_back(pShowerEnergy->at(j));
            allShowerTheta.push_back(pShowerTheta->at(j));
            allShowerSmallestTheta.push_back(pShowerSmallestTheta->at(j));
            allShowerETheta2.push_back(pShowerETheta2->at(j));
        }

        fileNumber++;
    }

    std::vector<event_t> allEvents = std::vector<event_t>();

    size_t newNumEvents = allEvent.size();

    for(UInt_t k = 0; k < newNumEvents; k++){
        event_t event;
        event.event = allEvent.at(k);
        event.run = allRun.at(k);
        event.subrun = allSubrun.at(k);
        event.trueNeutrinoVX = allTrueNeutrinoVX.at(k);
        event.trueNeutrinoVY = allTrueNeutrinoVY.at(k);
        event.trueNeutrinoVZ = allTrueNeutrinoVZ.at(k);
        event.recoNeutrinoVX = allRecoNeutrinoVX.at(k);
        event.recoNeutrinoVY = allRecoNeutrinoVY.at(k);
        event.recoNeutrinoVZ = allRecoNeutrinoVZ.at(k);
        event.trueCCNC = allTrueCCNC.at(k);
        event.trueNeutrinoType = allTrueNeutrinoType.at(k);
        event.trueLeptonType = allTrueLeptonType.at(k);
        event.nSlices = allNumSlices.at(k);
        event.nPFPs = allNumPfps.at(k);
        event.nTracks = allNumTracks.at(k);
        event.nShowers = allNumShowers.at(k);
        event.tpcID = allTpcID.at(k);
        event.dl_current = allDl_current.at(k);
        event.sliceCompleteness = allSliceCompleteness.at(k);
        event.trackScore = allTrackScore.at(k);
        event.showerEnergy = allShowerEnergy.at(k);
        event.showerTheta = allShowerTheta.at(k);
        event.showerSmallestTheta = allShowerSmallestTheta.at(k);
        event.showerETheta2 = allShowerETheta2.at(k);
        allEvents.push_back(event);
    }

    //deltaVertices(allEvents);
    //numberPlots(allEvents);
    //scores(allEvents);
    shower(allEvents);
}
