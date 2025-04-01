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
    double fullyReco;
	int DLCurrent;
	double event;
    double run;
    double subrun;
    double trueNeutrinoVX;
    double trueNeutrinoVY;
    double trueNeutrinoVZ;
    double trueCCNC;
    double trueNeutrinoType;
    double trueLeptonType;
    double recoNeutrinoVX;
    double recoNeutrinoVY;
    double recoNeutrinoVZ;
    double nSlices;
    double sliceCompleteness;
	double sliceNumPFPs;
	double nPFPs;
    double showerEnergy;
    double showerTheta;
	double showerTrackScore;
	double nShowers;
    double showerETheta2;
} event_t;

void styleDraw(TCanvas* canvas, TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax){
    canvas->cd();
    
    current->SetLineWidth(2);
    current->SetLineColor(kRed);

    cheated->SetLineWidth(2);
    cheated->SetLineColor(kSpring-5);

    dune->SetLineWidth(2);
    dune->SetLineColor(kViolet-5);

    uboone->SetLineWidth(2);
    uboone->SetLineColor(kBlue);

    current->Draw("hist");
    cheated->Draw("histsame");
    dune->Draw("histsame");
    uboone->Draw("histsame");

    current->SetStats(0);
    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(dune, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(uboone, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(current, "Current SBND Vertexing (without Refinement)", "f");
    legend->AddEntry(cheated, "Cheated SBND Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    canvas->SaveAs(filename);
} 

void percentage(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax){
    TCanvas *percentageCanvas = new TCanvas("percentage_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentPerc = (TH1F*) current->Clone("perc hist");
    currentPerc->Scale(100.0 * 1.0/sizeCurrent);
    currentPerc->GetYaxis()->SetTitle("Percentage of Events (%)"); 

    TH1F* cheatedPerc = (TH1F*) cheated->Clone("perc hist");
    cheatedPerc->Scale(100.0 * 1.0/sizeCheated);
    cheatedPerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* dunePerc = (TH1F*) dune->Clone("perc hist");
    dunePerc->Scale(100.0 * 1.0/sizeDune);
    dunePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* uboonePerc = (TH1F*) uboone->Clone("perc hist");
    uboonePerc->Scale(100.0 * 1.0/sizeUboone);
    uboonePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    styleDraw(percentageCanvas, currentPerc, cheatedPerc, dunePerc, uboonePerc, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax);
}

void shower(std::vector<event_t> dlUboone, std::vector<event_t> dlDune, std::vector<event_t> current, std::vector<event_t> cheated){
    
    TCanvas *showerECanvas = new TCanvas("showerE_canvas", "Graph Draw Options", 200, 10, 600, 400);
    //TH1F* showerECurrent_dist = new TH1F("Shower Energy", "Energy of Shower", 80, 0, 4000);
    TH1F* showerECurrent_dist = new TH1F("Shower Energy", "Energy of Shower", 60, 0, 3000);
    showerECurrent_dist->SetTitle("Energy of the Recoil Electron Shower in an Event;Energy of Shower (MeV);# of Events");
    TH1F* showerEDLUboone_dist = (TH1F*) showerECurrent_dist->Clone("Shower Energy");
    TH1F* showerEDLDune_dist = (TH1F*) showerECurrent_dist->Clone("Shower Energy");
    TH1F* showerECheated_dist = (TH1F*) showerECurrent_dist->Clone("Shower Energy");

    TCanvas *showerThetaCanvas = new TCanvas("showerTheta_canvas", "Graph Draw Options", 200, 10, 600, 400);
    //TH1F* showerThetaCurrent_dist = new TH1F("Shower Theta", "Theta of Shower", 90, 0, 360);
    TH1F* showerThetaCurrent_dist = new TH1F("Shower Theta", "Theta of Shower", 45, 0, 180);
    showerThetaCurrent_dist->SetTitle("Angle of the Recoil Electron Shower Relative to Beam Neutrino Direction in an Event;Theta of Shower (degrees);# of Events");
    TH1F* showerThetaDLUboone_dist = (TH1F*) showerThetaCurrent_dist->Clone("Shower Theta");
    TH1F* showerThetaDLDune_dist = (TH1F*) showerThetaCurrent_dist->Clone("Shower Theta");
    TH1F* showerThetaCheated_dist = (TH1F*) showerThetaCurrent_dist->Clone("Shower Theta");


    TCanvas *showerETheta2Canvas = new TCanvas("showerETheta2_canvas", "Graph Draw Options", 200, 10, 600, 400);
    //TH1F* showerETheta2Current_dist = new TH1F("Shower ETheta2", "ETheta2 of Shower", 200, 0, 204.4);
    TH1F* showerETheta2Current_dist = new TH1F("Shower ETheta2", "ETheta2 of Shower", 20, 0, 20.44);
    showerETheta2Current_dist->SetTitle("E#theta^{2} of the Recoil Electron Shower in an Event;E#theta^{2} of Shower (MeV);# of Events");
    TH1F* showerETheta2DLUboone_dist = (TH1F*) showerETheta2Current_dist->Clone("Shower ETheta2");
    TH1F* showerETheta2DLDune_dist = (TH1F*) showerETheta2Current_dist->Clone("Shower ETheta2");
    TH1F* showerETheta2Cheated_dist = (TH1F*) showerETheta2Current_dist->Clone("Shower ETheta2");

    TCanvas *showerTrackScoreCanvas = new TCanvas("showerTrackscore_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* showerTrackScoreCurrent_dist = new TH1F("Shower Trackscore", "Trackscore of Shower", 51, -0.01, 1.01);
    showerTrackScoreCurrent_dist->SetTitle("Trackscore of the Recoil Electron Shower in an Event;Trackscore;# of Events");
    TH1F* showerTrackScoreDLUboone_dist = (TH1F*) showerTrackScoreCurrent_dist->Clone("Shower Trackscore");
    TH1F* showerTrackScoreDLDune_dist = (TH1F*) showerTrackScoreCurrent_dist->Clone("Shower Trackscore");
    TH1F* showerTrackScoreCheated_dist = (TH1F*) showerTrackScoreCurrent_dist->Clone("Shower Trackscore");


    TCanvas *showerNumCanvas = new TCanvas("showerNum_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* showerNumCurrent_dist = new TH1F("Shower Num", "Number of Showers", 11, -0.5, 10.5);
    showerNumCurrent_dist->SetTitle("Number of Showers in an Event;Number of Showers;# of Events");
    TH1F* showerNumDLUboone_dist = (TH1F*) showerNumCurrent_dist->Clone("Shower Num");
    TH1F* showerNumDLDune_dist = (TH1F*) showerNumCurrent_dist->Clone("Shower Num");
    TH1F* showerNumCheated_dist = (TH1F*) showerNumCurrent_dist->Clone("Shower Num");

    for(UInt_t j = 0; j < current.size(); j++){
        if(current.at(j).nShowers != 0){
            showerNumCurrent_dist->Fill(current.at(j).nShowers);
            if(!isnan(current.at(j).showerEnergy)) showerECurrent_dist->Fill(current.at(j).showerEnergy);
            if(!isnan(current.at(j).showerTheta)) showerThetaCurrent_dist->Fill(TMath::RadToDeg() * current.at(j).showerTheta);
            if(!isnan(current.at(j).showerETheta2)) showerETheta2Current_dist->Fill(current.at(j).showerETheta2);
            if(!isnan(current.at(j).showerTrackScore)) showerTrackScoreCurrent_dist->Fill(current.at(j).showerTrackScore);
        }
    }
 
    for(UInt_t j = 0; j < dlDune.size(); j++){
        if(dlDune.at(j).nShowers != 0){
            showerNumDLDune_dist->Fill(dlDune.at(j).nShowers);
            if(!isnan(dlDune.at(j).showerEnergy)) showerEDLDune_dist->Fill(dlDune.at(j).showerEnergy);
            if(!isnan(dlDune.at(j).showerTheta)) showerThetaDLDune_dist->Fill(TMath::RadToDeg() * dlDune.at(j).showerTheta);
            if(!isnan(dlDune.at(j).showerETheta2)) showerETheta2DLDune_dist->Fill(dlDune.at(j).showerETheta2);
            if(!isnan(dlDune.at(j).showerTrackScore)) showerTrackScoreDLDune_dist->Fill(dlDune.at(j).showerTrackScore);
        }
    }

    for(UInt_t j = 0; j < dlUboone.size(); j++){
        if(dlUboone.at(j).nShowers != 0){
            showerNumDLUboone_dist->Fill(dlUboone.at(j).nShowers);
            if(!isnan(dlUboone.at(j).showerEnergy)) showerEDLUboone_dist->Fill(dlUboone.at(j).showerEnergy);
            if(!isnan(dlUboone.at(j).showerTheta)) showerThetaDLUboone_dist->Fill(TMath::RadToDeg() * dlUboone.at(j).showerTheta);
            if(!isnan(dlUboone.at(j).showerETheta2)) showerETheta2DLUboone_dist->Fill(dlUboone.at(j).showerETheta2);
            if(!isnan(dlUboone.at(j).showerTrackScore)) showerTrackScoreDLUboone_dist->Fill(dlUboone.at(j).showerTrackScore);
        }
    }

    for(UInt_t j = 0; j < cheated.size(); j++){
        if(cheated.at(j).nShowers != 0){
            showerNumCheated_dist->Fill(cheated.at(j).nShowers);
            if(!isnan(cheated.at(j).showerEnergy)) showerECheated_dist->Fill(cheated.at(j).showerEnergy);
            if(!isnan(cheated.at(j).showerTheta)) showerThetaCheated_dist->Fill(TMath::RadToDeg() * cheated.at(j).showerTheta);
            if(!isnan(cheated.at(j).showerETheta2)) showerETheta2Cheated_dist->Fill(cheated.at(j).showerETheta2);
            if(!isnan(cheated.at(j).showerTrackScore)) showerTrackScoreCheated_dist->Fill(cheated.at(j).showerTrackScore);
        }
    }

    styleDraw(showerECanvas, showerECurrent_dist, showerECheated_dist, showerEDLDune_dist, showerEDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/showerE_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(showerThetaCanvas, showerThetaCurrent_dist, showerThetaCheated_dist, showerThetaDLDune_dist, showerThetaDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/showerTheta_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(showerETheta2Canvas, showerETheta2Current_dist, showerETheta2Cheated_dist, showerETheta2DLDune_dist, showerETheta2DLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/showerETheta2_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(showerNumCanvas, showerNumCurrent_dist, showerNumCheated_dist, showerNumDLDune_dist, showerNumDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/showerNum_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(showerTrackScoreCanvas, showerTrackScoreCurrent_dist, showerTrackScoreCheated_dist, showerTrackScoreDLDune_dist, showerTrackScoreDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/showerTrackScore_dist.pdf", 0.56, 0.88, 0.70, 0.86);

    percentage(showerECurrent_dist, showerECheated_dist, showerEDLDune_dist, showerEDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/showerE_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(showerThetaCurrent_dist, showerThetaCheated_dist, showerThetaDLDune_dist, showerThetaDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/showerTheta_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(showerETheta2Current_dist, showerETheta2Cheated_dist, showerETheta2DLDune_dist, showerETheta2DLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/showerETheta2_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(showerNumCurrent_dist, showerNumCheated_dist, showerNumDLDune_dist, showerNumDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/showerNumShowers_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(showerTrackScoreCurrent_dist, showerTrackScoreCheated_dist, showerTrackScoreDLDune_dist, showerTrackScoreDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/showerTrackScore_perc.pdf", 0.56, 0.88, 0.70, 0.86);   
}

void vertices(std::vector<event_t> dlUboone, std::vector<event_t> dlDune, std::vector<event_t> current, std::vector<event_t> cheated){

    TCanvas *deltaXCanvas = new TCanvas("deltaX_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaXCurrent_dist = new TH1F("Delta X", "X Coord Diff", 40, -5, 5);
    deltaXCurrent_dist->SetTitle("#Deltax Distribution: Nu + E Elastic Scattering Events;x_{Reco} - x_{True} (cm);# of Events");
    TH1F* deltaXDLUboone_dist = (TH1F*) deltaXCurrent_dist->Clone("deltaX");
    TH1F* deltaXDLDune_dist = (TH1F*) deltaXCurrent_dist->Clone("deltaX");
    TH1F* deltaXCheated_dist = (TH1F*) deltaXCurrent_dist->Clone("deltaX");

    TCanvas *deltaYCanvas = new TCanvas("deltaY_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaYCurrent_dist = new TH1F("DeltaYX", "Y Coord Diff", 40, -5, 5);
    deltaYCurrent_dist->SetTitle("#Deltay Distribution: Nu + E Elastic Scattering Events;y_{Reco} - y_{True} (cm);# of Events");
    TH1F* deltaYDLUboone_dist = (TH1F*) deltaYCurrent_dist->Clone("deltaY");
    TH1F* deltaYDLDune_dist = (TH1F*) deltaYCurrent_dist->Clone("deltaY");
    TH1F* deltaYCheated_dist = (TH1F*) deltaYCurrent_dist->Clone("deltaY");

    TCanvas *deltaZCanvas = new TCanvas("deltaZ_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaZCurrent_dist = new TH1F("Delta Z", "Z Coord Diff", 40, -5, 5);
    deltaZCurrent_dist->SetTitle("#Deltaz Distribution: Nu + E Elastic Scattering Events;z_{Reco} - z_{True} (cm);# of Events");
    TH1F* deltaZDLUboone_dist = (TH1F*) deltaZCurrent_dist->Clone("deltaZ");
    TH1F* deltaZDLDune_dist = (TH1F*) deltaZCurrent_dist->Clone("deltaZ");
    TH1F* deltaZCheated_dist = (TH1F*) deltaZCurrent_dist->Clone("deltaZ");
    
    TCanvas *deltaRCanvas = new TCanvas("deltaR_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TH1F* deltaRCurrent_dist = new TH1F("Delta R", "R Coord Diff", 20, 0, 5);
    deltaRCurrent_dist->SetTitle("#Delta#bar{r} Distribution: Nu + E Elastic Scattering Events;|#bar{r}_{Reco} - #bar{r}_{True}| (cm);# of Events");
    TH1F* deltaRDLUboone_dist = (TH1F*) deltaRCurrent_dist->Clone("deltaR");
    TH1F* deltaRDLDune_dist = (TH1F*) deltaRCurrent_dist->Clone("deltaR");
    TH1F* deltaRCheated_dist = (TH1F*) deltaRCurrent_dist->Clone("deltaR");

    for(UInt_t j = 0; j < current.size(); j++){
        if(!isnan(current.at(j).trueNeutrinoVX) && !isnan(current.at(j).recoNeutrinoVX)) deltaXCurrent_dist->Fill(current.at(j).recoNeutrinoVX - current.at(j).trueNeutrinoVX);
        if(!isnan(current.at(j).trueNeutrinoVY) && !isnan(current.at(j).recoNeutrinoVY)) deltaYCurrent_dist->Fill(current.at(j).recoNeutrinoVY - current.at(j).trueNeutrinoVY);
        if(!isnan(current.at(j).trueNeutrinoVZ) && !isnan(current.at(j).recoNeutrinoVZ)) deltaZCurrent_dist->Fill(current.at(j).recoNeutrinoVZ - current.at(j).trueNeutrinoVZ);
        
        if(!isnan(current.at(j).trueNeutrinoVX) && !isnan(current.at(j).recoNeutrinoVX) &&
           !isnan(current.at(j).trueNeutrinoVY) && !isnan(current.at(j).recoNeutrinoVY) &&
           !isnan(current.at(j).trueNeutrinoVZ) && !isnan(current.at(j).recoNeutrinoVZ)){
            double deltaR = std::sqrt( std::pow(current.at(j).trueNeutrinoVX - current.at(j).recoNeutrinoVX, 2) + std::pow(current.at(j).trueNeutrinoVY - current.at(j).recoNeutrinoVY, 2) + std::pow(current.at(j).trueNeutrinoVZ - current.at(j).recoNeutrinoVZ, 2));
            deltaRCurrent_dist->Fill(deltaR);
        }
    }

    for(UInt_t j = 0; j < dlDune.size(); j++){
        if(!isnan(dlDune.at(j).trueNeutrinoVX) && !isnan(dlDune.at(j).recoNeutrinoVX)) deltaXDLDune_dist->Fill(dlDune.at(j).recoNeutrinoVX - dlDune.at(j).trueNeutrinoVX);
        if(!isnan(dlDune.at(j).trueNeutrinoVY) && !isnan(dlDune.at(j).recoNeutrinoVY)) deltaYDLDune_dist->Fill(dlDune.at(j).recoNeutrinoVY - dlDune.at(j).trueNeutrinoVY);
        if(!isnan(dlDune.at(j).trueNeutrinoVZ) && !isnan(dlDune.at(j).recoNeutrinoVZ)) deltaZDLDune_dist->Fill(dlDune.at(j).recoNeutrinoVZ - dlDune.at(j).trueNeutrinoVZ);
        
        if(!isnan(dlDune.at(j).trueNeutrinoVX) && !isnan(dlDune.at(j).recoNeutrinoVX) &&
           !isnan(dlDune.at(j).trueNeutrinoVY) && !isnan(dlDune.at(j).recoNeutrinoVY) &&
           !isnan(dlDune.at(j).trueNeutrinoVZ) && !isnan(dlDune.at(j).recoNeutrinoVZ)){
            double deltaR = std::sqrt( std::pow(dlDune.at(j).trueNeutrinoVX - dlDune.at(j).recoNeutrinoVX, 2) + std::pow(dlDune.at(j).trueNeutrinoVY - dlDune.at(j).recoNeutrinoVY, 2) + std::pow(dlDune.at(j).trueNeutrinoVZ - dlDune.at(j).recoNeutrinoVZ, 2));
            deltaRDLDune_dist->Fill(deltaR);
        }
    }
    
    for(UInt_t j = 0; j < dlUboone.size(); j++){
        if(!isnan(dlUboone.at(j).trueNeutrinoVX) && !isnan(dlUboone.at(j).recoNeutrinoVX)) deltaXDLUboone_dist->Fill(dlUboone.at(j).recoNeutrinoVX - dlUboone.at(j).trueNeutrinoVX);
        if(!isnan(dlUboone.at(j).trueNeutrinoVY) && !isnan(dlUboone.at(j).recoNeutrinoVY)) deltaYDLUboone_dist->Fill(dlUboone.at(j).recoNeutrinoVY - dlUboone.at(j).trueNeutrinoVY);
        if(!isnan(dlUboone.at(j).trueNeutrinoVZ) && !isnan(dlUboone.at(j).recoNeutrinoVZ)) deltaZDLUboone_dist->Fill(dlUboone.at(j).recoNeutrinoVZ - dlUboone.at(j).trueNeutrinoVZ);
        
        if(!isnan(dlUboone.at(j).trueNeutrinoVX) && !isnan(dlUboone.at(j).recoNeutrinoVX) &&
           !isnan(dlUboone.at(j).trueNeutrinoVY) && !isnan(dlUboone.at(j).recoNeutrinoVY) &&
           !isnan(dlUboone.at(j).trueNeutrinoVZ) && !isnan(dlUboone.at(j).recoNeutrinoVZ)){
            double deltaR = std::sqrt( std::pow(dlUboone.at(j).trueNeutrinoVX - dlUboone.at(j).recoNeutrinoVX, 2) + std::pow(dlUboone.at(j).trueNeutrinoVY - dlUboone.at(j).recoNeutrinoVY, 2) + std::pow(dlUboone.at(j).trueNeutrinoVZ - dlUboone.at(j).recoNeutrinoVZ, 2));
            deltaRDLUboone_dist->Fill(deltaR);
        }
    }

    for(UInt_t j = 0; j < cheated.size(); j++){
        if(!isnan(cheated.at(j).trueNeutrinoVX) && !isnan(cheated.at(j).recoNeutrinoVX)) deltaXCheated_dist->Fill(cheated.at(j).recoNeutrinoVX - cheated.at(j).trueNeutrinoVX);
        if(!isnan(cheated.at(j).trueNeutrinoVY) && !isnan(cheated.at(j).recoNeutrinoVY)) deltaYCheated_dist->Fill(cheated.at(j).recoNeutrinoVY - cheated.at(j).trueNeutrinoVY);
        if(!isnan(cheated.at(j).trueNeutrinoVZ) && !isnan(cheated.at(j).recoNeutrinoVZ)) deltaZCheated_dist->Fill(cheated.at(j).recoNeutrinoVZ - cheated.at(j).trueNeutrinoVZ);
        
        if(!isnan(cheated.at(j).trueNeutrinoVX) && !isnan(cheated.at(j).recoNeutrinoVX) &&
           !isnan(cheated.at(j).trueNeutrinoVY) && !isnan(cheated.at(j).recoNeutrinoVY) &&
           !isnan(cheated.at(j).trueNeutrinoVZ) && !isnan(cheated.at(j).recoNeutrinoVZ)){
            double deltaR = std::sqrt( std::pow(cheated.at(j).trueNeutrinoVX - cheated.at(j).recoNeutrinoVX, 2) + std::pow(cheated.at(j).trueNeutrinoVY - cheated.at(j).recoNeutrinoVY, 2) + std::pow(cheated.at(j).trueNeutrinoVZ - cheated.at(j).recoNeutrinoVZ, 2));
            deltaRCheated_dist->Fill(deltaR);
        }
    }

    styleDraw(deltaXCanvas, deltaXCurrent_dist, deltaXCheated_dist, deltaXDLDune_dist, deltaXDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/deltaX_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(deltaYCanvas, deltaYCurrent_dist, deltaYCheated_dist, deltaYDLDune_dist, deltaYDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/deltaY_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(deltaZCanvas, deltaZCurrent_dist, deltaZCheated_dist, deltaZDLDune_dist, deltaZDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/deltaZ_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(deltaRCanvas, deltaRCurrent_dist, deltaRCheated_dist, deltaRDLDune_dist, deltaRDLUboone_dist, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/deltaR_dist.pdf", 0.56, 0.88, 0.70, 0.86);

    percentage(deltaXCurrent_dist, deltaXCheated_dist, deltaXDLDune_dist, deltaXDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/deltaX_perc.pdf", 0.56, 0.88, 0.70, 0.86); 
    percentage(deltaYCurrent_dist, deltaYCheated_dist, deltaYDLDune_dist, deltaYDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/deltaY_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(deltaZCurrent_dist, deltaZCheated_dist, deltaZDLDune_dist, deltaZDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/deltaZ_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(deltaRCurrent_dist, deltaRCheated_dist, deltaRDLDune_dist, deltaRDLUboone_dist, current.size(), cheated.size(), dlDune.size(), dlUboone.size(), 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/deltaR_perc.pdf", 0.56, 0.88, 0.70, 0.86);
}

void slices(std::vector<event_t> dlUboone, std::vector<event_t> dlDune, std::vector<event_t> current, std::vector<event_t> cheated){
}

void pfps(std::vector<event_t> dlUboone, std::vector<event_t> dlDune, std::vector<event_t> current, std::vector<event_t> cheated){
}

void nuE_macro(){
    //TFile *f = TFile::Open("/exp/sbnd/app/users/coackley/nuev10_04_05/NuEAnalyserOutput.root", "READ");
    TFile *f = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/merged.root", "READ");
    
    if(!f){
        std::cout << "Failed to read file" << std::endl;
        return;
    }

    TTree *t;
    f->GetObject("NuETree", t);

    std::vector<double> allFullyReco = std::vector<double>(0);
    std::vector<int> allDLCurrent = std::vector<int>(0);    
    std::vector<double> allEvent = std::vector<double>(0);
    std::vector<double> allRun = std::vector<double>(0);
    std::vector<double> allSubrun = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVX = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVY = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoVZ = std::vector<double>(0);
    std::vector<double> allTrueCCNC = std::vector<double>(0);
    std::vector<double> allTrueNeutrinoType = std::vector<double>(0); 
    std::vector<double> allTrueLeptonType = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVX = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVY = std::vector<double>(0);
    std::vector<double> allRecoNeutrinoVZ = std::vector<double>(0);
    std::vector<double> allNumSlices = std::vector<double>(0);
    std::vector<double> allSliceCompleteness = std::vector<double>(0);
    std::vector<double> allSliceNumPFPs = std::vector<double>(0);
    std::vector<double> allNumPFPs = std::vector<double>(0);
    std::vector<double> allShowerEnergy = std::vector<double>(0);
    std::vector<double> allShowerTheta = std::vector<double>(0);
    std::vector<double> allShowerTrackScore = std::vector<double>(0);
    std::vector<double> allNumShowers = std::vector<double>(0);
    std::vector<double> allShowerETheta2 = std::vector<double>(0);

    Long64_t numEntries = t->GetEntries();
    int fileNumber = 1;

    for(Long64_t i = 0; i < numEntries; i++){
        std::cout << "File number: " << fileNumber << std::endl;
        std::vector<double> *pFullyReco = 0;
        std::vector<int> *pDLCurrent = 0;
        std::vector<double> *pEvent = 0;
        std::vector<double> *pRun = 0; 
        std::vector<double> *pSubrun = 0;
        std::vector<double> *pTrueNeutrinoVX = 0;
        std::vector<double> *pTrueNeutrinoVY = 0;
        std::vector<double> *pTrueNeutrinoVZ = 0;
        std::vector<double> *pTrueCCNC = 0;
        std::vector<double> *pTrueNeutrinoType = 0;
        std::vector<double> *pTrueLeptonType = 0;
        std::vector<double> *pRecoNeutrinoVX = 0;
        std::vector<double> *pRecoNeutrinoVY = 0;
        std::vector<double> *pRecoNeutrinoVZ = 0;
        std::vector<double> *pNumSlices = 0;
        std::vector<double> *pSliceCompleteness = 0;
        std::vector<double> *pSliceNumPFPs = 0;
        std::vector<double> *pNumPFPs = 0;
        std::vector<double> *pShowerEnergy = 0;
        std::vector<double> *pShowerTheta = 0;
        std::vector<double> *pShowerTrackScore = 0;
        std::vector<double> *pNumShowers = 0;
        std::vector<double> *pShowerETheta2 = 0;

        TBranch *fullyReco_branch = 0;
        TBranch *DLCurrent_branch = 0;
        TBranch *event_branch = 0;
        TBranch *run_branch = 0;
        TBranch *subrun_branch = 0;
        TBranch *trueNeutrinoVX_branch = 0;
        TBranch *trueNeutrinoVY_branch = 0;
        TBranch *trueNeutrinoVZ_branch = 0;
        TBranch *trueCCNC_branch = 0;
        TBranch *trueNeutrinoType_branch = 0;
        TBranch *trueLeptonType_branch = 0;
        TBranch *recoNeutrinoVX_branch = 0;
        TBranch *recoNeutrinoVY_branch = 0;
        TBranch *recoNeutrinoVZ_branch = 0;
        TBranch *numSlices_branch = 0;
        TBranch *sliceCompleteness_branch = 0;
        TBranch *sliceNumPFPs_branch = 0;
        TBranch *numPFPs_branch = 0;
        TBranch *showerEnergy_branch = 0;
        TBranch *showerTheta_branch = 0;
        TBranch *showerTrackScore_branch = 0;
        TBranch *numShowers_branch = 0;
        TBranch *showerETheta2_branch = 0;

        t->SetBranchAddress("fullyReco_tree", &pFullyReco, &fullyReco_branch);
        t->SetBranchAddress("DLCurrent_tree", &pDLCurrent, &DLCurrent_branch);
        t->SetBranchAddress("event_tree", &pEvent, &event_branch);
        t->SetBranchAddress("run_tree", &pRun, &run_branch);
        t->SetBranchAddress("subrun_tree", &pSubrun, &subrun_branch);
        t->SetBranchAddress("trueNeutrinoVX_tree", &pTrueNeutrinoVX, &trueNeutrinoVX_branch);
        t->SetBranchAddress("trueNeutrinoVY_tree", &pTrueNeutrinoVY, &trueNeutrinoVY_branch);
        t->SetBranchAddress("trueNeutrinoVZ_tree", &pTrueNeutrinoVZ, &trueNeutrinoVZ_branch);
        t->SetBranchAddress("trueCCNC_tree", &pTrueCCNC, &trueCCNC_branch);
        t->SetBranchAddress("trueNeutrinoType_tree", &pTrueNeutrinoType, &trueNeutrinoType_branch);
        t->SetBranchAddress("trueLeptonType_tree", &pTrueLeptonType, &trueLeptonType_branch);
        t->SetBranchAddress("recoNeutrinoVX_tree", &pRecoNeutrinoVX, &recoNeutrinoVX_branch);
        t->SetBranchAddress("recoNeutrinoVY_tree", &pRecoNeutrinoVY, &recoNeutrinoVY_branch);
        t->SetBranchAddress("recoNeutrinoVZ_tree", &pRecoNeutrinoVZ, &recoNeutrinoVZ_branch);
        t->SetBranchAddress("numSlices_tree", &pNumSlices, &numSlices_branch);
        t->SetBranchAddress("sliceCompleteness_tree", &pSliceCompleteness, &sliceCompleteness_branch);
        t->SetBranchAddress("sliceNumPFPs_tree", &pSliceNumPFPs, &sliceNumPFPs_branch);
        t->SetBranchAddress("numPFPs_tree", &pNumPFPs, &numPFPs_branch);
        t->SetBranchAddress("showerEnergy_tree", &pShowerEnergy, &showerEnergy_branch);
        t->SetBranchAddress("showerTheta_tree", &pShowerTheta, &showerTheta_branch);
        t->SetBranchAddress("showerTrackScore_tree", &pShowerTrackScore, &showerTrackScore_branch);
        t->SetBranchAddress("numShowers_tree", &pNumShowers, &numShowers_branch);
        t->SetBranchAddress("showerETheta2_tree", &pShowerETheta2, &showerETheta2_branch);

        fullyReco_branch->GetEntry(i);
        DLCurrent_branch->GetEntry(i);
        event_branch->GetEntry(i);
        run_branch->GetEntry(i);
        subrun_branch->GetEntry(i);
        trueNeutrinoVX_branch->GetEntry(i);
        trueNeutrinoVY_branch->GetEntry(i);
        trueNeutrinoVZ_branch->GetEntry(i);
        trueCCNC_branch->GetEntry(i);
        trueNeutrinoType_branch->GetEntry(i);
        trueLeptonType_branch->GetEntry(i);
        recoNeutrinoVX_branch->GetEntry(i);
        recoNeutrinoVY_branch->GetEntry(i);
        recoNeutrinoVZ_branch->GetEntry(i);
        numSlices_branch->GetEntry(i);
        sliceCompleteness_branch->GetEntry(i);
        sliceNumPFPs_branch->GetEntry(i);
        numPFPs_branch->GetEntry(i);
        showerEnergy_branch->GetEntry(i);
        showerTheta_branch->GetEntry(i);
        showerTrackScore_branch->GetEntry(i);
        numShowers_branch->GetEntry(i);
        showerETheta2_branch->GetEntry(i);

        size_t vecSize = pEvent->size();
        for(UInt_t j = 0; j < vecSize; j++){
            allFullyReco.push_back(pFullyReco->at(j));
            allDLCurrent.push_back(pDLCurrent->at(j));
            allEvent.push_back(pEvent->at(j));
            allRun.push_back(pRun->at(j));
            allSubrun.push_back(pSubrun->at(j));
            allTrueNeutrinoVX.push_back(pTrueNeutrinoVX->at(j));
            allTrueNeutrinoVY.push_back(pTrueNeutrinoVY->at(j));
            allTrueNeutrinoVZ.push_back(pTrueNeutrinoVZ->at(j));
            allTrueCCNC.push_back(pTrueCCNC->at(j));
            allTrueNeutrinoType.push_back(pTrueNeutrinoType->at(j));
            allTrueLeptonType.push_back(pTrueLeptonType->at(j));
            allRecoNeutrinoVX.push_back(pRecoNeutrinoVX->at(j));
            allRecoNeutrinoVY.push_back(pRecoNeutrinoVY->at(j));
            allRecoNeutrinoVZ.push_back(pRecoNeutrinoVZ->at(j));
            allNumSlices.push_back(pNumSlices->at(j));
            allSliceCompleteness.push_back(pSliceCompleteness->at(j));
            allSliceNumPFPs.push_back(pSliceNumPFPs->at(j));
            allNumPFPs.push_back(pNumPFPs->at(j));
            allShowerEnergy.push_back(pShowerEnergy->at(j));
            allShowerTheta.push_back(pShowerTheta->at(j));
            allShowerTrackScore.push_back(pShowerTrackScore->at(j));
            allNumShowers.push_back(pNumShowers->at(j));
            allShowerETheta2.push_back(pShowerETheta2->at(j));
        }

        fileNumber++;
    }

    std::vector<event_t> allEvents = std::vector<event_t>();
    std::vector<event_t> dlDuneEvents = std::vector<event_t>();
    std::vector<event_t> dlUbooneEvents = std::vector<event_t>();
    std::vector<event_t> currentEvents = std::vector<event_t>();
    std::vector<event_t> cheatedEvents = std::vector<event_t>();

    size_t newNumEvents = allEvent.size();

    for(UInt_t k = 0; k < newNumEvents; k++){
        event_t event;
        event.fullyReco = allFullyReco.at(k);
        event.DLCurrent = allDLCurrent.at(k);
        event.event = allEvent.at(k);
        event.run = allRun.at(k);
        event.subrun = allSubrun.at(k);
        event.trueNeutrinoVX = allTrueNeutrinoVX.at(k);
        event.trueNeutrinoVY = allTrueNeutrinoVY.at(k);
        event.trueNeutrinoVZ = allTrueNeutrinoVZ.at(k);
        event.trueCCNC = allTrueCCNC.at(k);
        event.trueNeutrinoType = allTrueNeutrinoType.at(k);
        event.trueLeptonType = allTrueLeptonType.at(k);
        event.recoNeutrinoVX = allRecoNeutrinoVX.at(k);
        event.recoNeutrinoVY = allRecoNeutrinoVY.at(k);
        event.recoNeutrinoVZ = allRecoNeutrinoVZ.at(k);
        event.nSlices = allNumSlices.at(k);
        event.sliceCompleteness = allSliceCompleteness.at(k);
        event.sliceNumPFPs = allSliceNumPFPs.at(k);
        event.nPFPs = allSliceNumPFPs.at(k);
        event.showerEnergy = allShowerEnergy.at(k);
        event.showerTheta = allShowerTheta.at(k);
        event.showerTrackScore = allShowerTrackScore.at(k);
        event.nShowers = allNumShowers.at(k);
        event.showerETheta2 = allShowerETheta2.at(k);
        allEvents.push_back(event);

        if(event.DLCurrent == 0) dlUbooneEvents.push_back(event);
        if(event.DLCurrent == 1) dlDuneEvents.push_back(event);
        if(event.DLCurrent == 2) currentEvents.push_back(event);
        if(event.DLCurrent == 3) cheatedEvents.push_back(event); 
    }

    std::cout << "uboone: " << dlUbooneEvents.size() << ", dune: " << dlDuneEvents.size() << ", current: " << currentEvents.size() << ", cheated: " << cheatedEvents.size() << std::endl;

    shower(dlUbooneEvents, dlDuneEvents, currentEvents, cheatedEvents);
    vertices(dlUbooneEvents, dlDuneEvents, currentEvents, cheatedEvents);
}
