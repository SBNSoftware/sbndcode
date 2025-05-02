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

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* current;
    TH1F* cheated;
    TH1F* dune;
    TH1F* uboone;
} histGroup;

typedef struct{
    double pdg;
    double isPrimary;
    double vx;
    double vy;
    double vz;
    double dx;
    double dy;
    double dz;
    double sliceID;
    double bestPlaneEnergy;
    double theta;
    double trackscore;
    double completeness;
    double purity;
} recoParticle;

typedef struct{
    double pdg;
    double isPrimary;
    double vx;
    double vy;
    double vz;
    double sliceID;
} recoNeutrino;

typedef struct{
    double vx;
    double vy;
    double vz;
    double CCNC;
    double pdg;
    double leptonpdg;
    double tpcID;
    double tpcValid;
} trueNeutrino;

typedef struct{
    double pdg;
    double vx;
    double vy;
    double vz;
    double px;
    double py;
    double pz;
    double energy;
    double angle;
    double ETheta2;
    double dx;
    double dy;
    double dz;
} trueParticle;

typedef struct{
    double id;
    double completeness;
    double purity;
    double score;
} recoSlice;

typedef struct{
    std::vector<recoParticle> recoParticleVec;
    std::vector<recoNeutrino> recoNeutrinoVec;
    std::vector<trueNeutrino> trueNeutrinoVec;
    std::vector<trueParticle> trueParticleVec;
    std::vector<recoSlice>    sliceVec;
    double DLCurrent;
} eventData;

typedef struct{
    recoParticle chosenRecoParticle;
    recoNeutrino chosenRecoNeutrino;
    trueNeutrino chosenTrueNeutrino;
    trueParticle chosenTrueParticle;
    recoSlice chosenSlice;
    double totalPFPEnergy;
    double DLCurrent;
} eventAfterCuts;

typedef struct{
    eventData eventBeforeCut;
    eventAfterCuts eventAfterCut;
    double event;
    double run;
    double subrun;
} allEventData;

histGroup createHistGroup(const std::string& baseName, const std::string& title, const std::string& xAxisTitle, int bins, float xlow, float xup){
    TCanvas* canvas = new TCanvas((baseName + "_canvas").c_str(), "Graph Draw Options", 200, 10, 600, 400);
    
    TH1F* base = new TH1F(baseName.c_str(), title.c_str(), bins, xlow, xup);
    base->SetTitle((title + ";" + xAxisTitle + ";# of Events").c_str());

    return {
        canvas,
        base,
        (TH1F*) base->Clone((baseName + "_current").c_str()),
        (TH1F*) base->Clone((baseName + "_cheated").c_str()),
        (TH1F*) base->Clone((baseName + "_dldune").c_str()),
        (TH1F*) base->Clone((baseName + "_dluboone").c_str())
    };    
}

void styleDraw(TCanvas* canvas, TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, TPaveText* pt = nullptr, int* percentage = nullptr, int* drawLine = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    gStyle->SetPalette(kAvocado);
    gROOT->ForceStyle();
    gPad->Update();

    std::cout << "Palette size: " << gStyle->GetNumberOfColors() << std::endl;
    current->SetLineWidth(2);
    current->SetLineColor(TColor::GetColorPalette(150));

    cheated->SetLineWidth(2);
    cheated->SetLineColor(TColor::GetColorPalette(200));

    dune->SetLineWidth(2);
    dune->SetLineColor(TColor::GetColorPalette(50));

    uboone->SetLineWidth(2);
    uboone->SetLineColor(TColor::GetColorPalette(100));

    current->Draw("hist");
    cheated->Draw("histsame");
    dune->Draw("histsame");
    uboone->Draw("histsame");

    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    current->SetStats(0);
    current->GetXaxis()->SetTickLength(0.04);
    current->GetYaxis()->SetTickLength(0.03);
    current->GetXaxis()->SetTickSize(0.02);
    current->GetYaxis()->SetTickSize(0.02);

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(dune, "Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(uboone, "Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(current, "Current SBND Vertexing (without Refinement)", "f");
    legend->AddEntry(cheated, "Cheated SBND Vertexing", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, current->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    if(pt){
        pt->SetTextSize(legend->GetTextSize());
        pt->SetTextFont(legend->GetTextFont());
        pt->Draw();
    }

    canvas->SaveAs(filename);
} 

void percentage(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr){
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

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 
   
    int funcValue = 1;

    styleDraw(percentageCanvas, currentPerc, cheatedPerc, dunePerc, uboonePerc, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt, &funcValue, drawLine);
}

void efficiency(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr){
    TCanvas *efficiencyCanvas = new TCanvas("efficiency_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentEff = (TH1F*) current->Clone("eff hist");
    currentEff->Reset();
    currentEff->GetYaxis()->SetTitle("Efficiency"); 

    TH1F* cheatedEff = (TH1F*) cheated->Clone("eff hist");
    cheatedEff->Reset();
    cheatedEff->GetYaxis()->SetTitle("Efficiency");

    TH1F* duneEff = (TH1F*) dune->Clone("eff hist");
    duneEff->Reset();
    duneEff->GetYaxis()->SetTitle("Efficiency");

    TH1F* ubooneEff = (TH1F*) uboone->Clone("eff hist");
    ubooneEff->Reset();
    ubooneEff->GetYaxis()->SetTitle("Efficiency");

    int numBins = current->GetNbinsX();
    double currentSum = 0.0;
    double cheatedSum = 0.0;
    double duneSum = 0.0;
    double ubooneSum = 0.0;

    for(int i = 1; i <= numBins; ++i){
        currentSum += current->GetBinContent(i);
        cheatedSum += cheated->GetBinContent(i);
        duneSum += dune->GetBinContent(i);
        ubooneSum += uboone->GetBinContent(i);

        double currentEffValue = currentSum/sizeCurrent;
        double cheatedEffValue = cheatedSum/sizeCheated;
        double duneEffValue = duneSum/sizeDune;
        double ubooneEffValue = ubooneSum/sizeUboone;

        currentEff->SetBinContent(i, currentEffValue);
        cheatedEff->SetBinContent(i, cheatedEffValue);
        duneEff->SetBinContent(i, duneEffValue);
        ubooneEff->SetBinContent(i, ubooneEffValue);
    }

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 

    int funcValue = 1;

    styleDraw(efficiencyCanvas, currentEff, cheatedEff, duneEff, ubooneEff, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt = nullptr, &funcValue, drawLine);
}

recoSlice chooseRecoSlice(std::vector<recoSlice> recoSliceVec){
    int chosenSliceIndex;
    double highestCompleteness = -100;

    for(size_t j = 0; j < recoSliceVec.size(); ++j){
        if(recoSliceVec[j].completeness > highestCompleteness){
            highestCompleteness = recoSliceVec[j].completeness;
            chosenSliceIndex = j;
        }
    }
    
    return recoSliceVec[chosenSliceIndex];
}

recoParticle chooseRecoParticle(std::vector<recoParticle> recoParticleVec, double& totalEnergy){
   totalEnergy = 0; 
   int chosenParticleIndex;
   double highestEnergy = -10000;

    for(size_t j = 0; j < recoParticleVec.size(); ++j){
        totalEnergy += recoParticleVec[j].bestPlaneEnergy;
        if(recoParticleVec[j].bestPlaneEnergy > highestEnergy){
            highestEnergy = recoParticleVec[j].bestPlaneEnergy;
            chosenParticleIndex = j;
        }
    }

    return recoParticleVec[chosenParticleIndex];
}

void numberPlots(std::vector<allEventData> allEventsData){
 
    auto numSlices = createHistGroup("numSlices", "Number of Slices in an Event", "Number of Slices", 5, 0, 5);
    auto numPFPs = createHistGroup("numPFPs", "Number of PFPs in an Event", "Number of PFPs", 10, 0, 10);
    auto numPFPs1Slice = createHistGroup("numPFPs1Slice", "Number of PFPs in 1-Slice Events", "Number of PFPs", 10, 0, 10);
    auto numPFPsSlices = createHistGroup("numPFPsSlices", "Number of PFPs in Multi-Slice Events", "Number of PFPs", 10, 0, 10);
   
    double sizeCurrent = 0;
    double sizeCheated = 0;
    double sizeDLDune = 0;
    double sizeDLUboone = 0;
    double sizeCurrent1Slice = 0;
    double sizeCheated1Slice = 0;
    double sizeDLDune1Slice = 0;
    double sizeDLUboone1Slice = 0;
    double sizeCurrentSlices = 0;
    double sizeCheatedSlices = 0;
    double sizeDLDuneSlices = 0;
    double sizeDLUbooneSlices = 0;
    

    for(const auto& event : allEventsData){
        if(!(event.eventBeforeCut.trueNeutrinoVec[0].pdg == -999999 || event.eventBeforeCut.trueNeutrinoVec[0].tpcValid == 0)){
            size_t recoParticleCount = event.eventBeforeCut.recoParticleVec.size();
            size_t recoSliceCount = event.eventBeforeCut.sliceVec.size();

            if(recoParticleCount == 1 && event.eventBeforeCut.recoParticleVec[0].pdg == -999999) {
                recoParticleCount = 0;
            } 
        
            if(recoSliceCount == 1 && event.eventBeforeCut.sliceVec[0].id == -999999){
                recoSliceCount = 0;
            }

            //if(recoParticleCount > 1 && recoSliceCount == 1){
            //if(event.event == 44){
                //printf("Event: %f, Run %f, SubRun: %f\n", event.event, event.run, event.subrun);
                //printf("Number of PFPs: %zu\n", recoParticleCount);
            //}

            //if(recoSliceCount > 1){
            //if(event.event == 44){
                //printf("Event: %f, Run %f, SubRun: %f\n", event.event, event.run, event.subrun);
                //printf("Number of Slices: %zu, Number of PFPs: %zu\n", recoSliceCount, recoParticleCount);
            //}

            if(event.eventBeforeCut.DLCurrent == 0){
                sizeDLUboone++;
                numPFPs.uboone->Fill(recoParticleCount);
                numSlices.uboone->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1Slice.uboone->Fill(recoParticleCount);
                    sizeDLUboone1Slice++;
                } else if(recoSliceCount > 1){
                    numPFPsSlices.uboone->Fill(recoParticleCount);
                    sizeDLUbooneSlices++;
                }
            } else if(event.eventBeforeCut.DLCurrent == 1){
                sizeDLDune++;
                numPFPs.dune->Fill(recoParticleCount);
                numSlices.dune->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1Slice.dune->Fill(recoParticleCount);
                    sizeDLDune1Slice++;
                } else if(recoSliceCount > 1){
                    numPFPsSlices.dune->Fill(recoParticleCount);
                    sizeDLDuneSlices++;
                }
            } else if(event.eventBeforeCut.DLCurrent == 2){
                sizeCurrent++;
                numPFPs.current->Fill(recoParticleCount);
                numSlices.current->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1Slice.current->Fill(recoParticleCount);
                    sizeCurrent1Slice++;
                } else if(recoSliceCount > 1){
                    numPFPsSlices.current->Fill(recoParticleCount); 
                    sizeCurrentSlices++;
                }
            } else if(event.eventBeforeCut.DLCurrent == 3){
                sizeCheated++;
                numPFPs.cheated->Fill(recoParticleCount);
                numSlices.cheated->Fill(recoSliceCount);
                if(recoSliceCount == 1){
                    numPFPs1Slice.cheated->Fill(recoParticleCount);
                    sizeCheated1Slice++;
                } else if(recoSliceCount > 1){
                    numPFPsSlices.cheated->Fill(recoParticleCount);
                    sizeCheatedSlices++;
                }
            }
        }
    }   


    styleDraw(numSlices.canvas, numSlices.current, numSlices.cheated, numSlices.dune, numSlices.uboone, 0, 1000, 999, 999, "/nashome/c/coackley/nuEPlots/numSlices_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(numPFPs.canvas, numPFPs.current, numPFPs.cheated, numPFPs.dune, numPFPs.uboone, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(numPFPs1Slice.canvas, numPFPs1Slice.current, numPFPs1Slice.cheated, numPFPs1Slice.dune, numPFPs1Slice.uboone, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs1Slice_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(numPFPsSlices.canvas, numPFPsSlices.current, numPFPsSlices.cheated, numPFPsSlices.dune, numPFPsSlices.uboone, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPsSlices_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    
    percentage(numSlices.current, numSlices.cheated, numSlices.dune, numSlices.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/numSlices_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(numPFPs.current, numPFPs.cheated, numPFPs.dune, numPFPs.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 95, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(numPFPs1Slice.current, numPFPs1Slice.cheated, numPFPs1Slice.dune, numPFPs1Slice.uboone, sizeCurrent1Slice, sizeCheated1Slice, sizeDLDune1Slice, sizeDLUboone1Slice, 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPs1Slice_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(numPFPsSlices.current, numPFPsSlices.cheated, numPFPsSlices.dune, numPFPsSlices.uboone, sizeCurrentSlices, sizeCheatedSlices, sizeDLDuneSlices, sizeDLUbooneSlices, 0, 95, 999, 999, "/nashome/c/coackley/nuEPlots/numPFPsSlices_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    //printf("Number of events: Current = %f, Cheated %f, DL Dune = %f, DL Uboone = %f\n", sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone);    

}

void chosenSlicePlots(std::vector<allEventData> allEvents){
    auto sliceCompleteness = createHistGroup("sliceCompleteness", "Completeness of the Chosen Slice in an Event", "Slice Completeness", 44, 0.8, 1.02);
    auto sliceCompleteness1Slice = createHistGroup("sliceCompleteness1Slice", "Completeness of Slice in an Event with 1 Slice", "Slice Completeness", 44, 0.8, 1.02);
    auto sliceCompletenessSlices = createHistGroup("sliceCompletenessSlices", "Completeness of the Chosen Slice in an Event with > 1 Slice", "Slice Completeness", 51, 0, 1.02);
    auto sliceScore = createHistGroup("sliceScore", "Score of the Chosen Slice in an Event", "Slice Score", 25, 0, 1);

    double sizeCurrent = 0;
    double sizeCheated = 0;
    double sizeDLDune = 0;
    double sizeDLUboone = 0;
    double sizeCurrent1Slice = 0;
    double sizeCheated1Slice = 0;
    double sizeDLDune1Slice = 0;
    double sizeDLUboone1Slice = 0;

    for(const auto& event : allEvents){
        if(!(event.eventBeforeCut.trueNeutrinoVec[0].pdg == -999999 || event.eventBeforeCut.trueNeutrinoVec[0].tpcValid == 0 || (event.eventBeforeCut.sliceVec.size() == 1 && event.eventBeforeCut.sliceVec[0].id == -999999))){
            if(event.eventBeforeCut.DLCurrent == 0){
                sizeDLUboone++;
                sliceCompleteness.uboone->Fill(event.eventAfterCut.chosenSlice.completeness);
                sliceScore.uboone->Fill(event.eventAfterCut.chosenSlice.score);
                if(event.eventBeforeCut.sliceVec.size() == 1){
                    sizeDLUboone1Slice++;
                    sliceCompleteness1Slice.uboone->Fill(event.eventAfterCut.chosenSlice.completeness);
                } else if(event.eventBeforeCut.sliceVec.size() > 1){
                    sliceCompletenessSlices.uboone->Fill(event.eventAfterCut.chosenSlice.completeness);
                }
            } else if(event.eventBeforeCut.DLCurrent == 1){
                sizeDLDune++;
                sliceCompleteness.dune->Fill(event.eventAfterCut.chosenSlice.completeness);
                sliceScore.dune->Fill(event.eventAfterCut.chosenSlice.score);
                if(event.eventBeforeCut.sliceVec.size() == 1){
                    sizeDLDune1Slice++;
                    sliceCompleteness1Slice.dune->Fill(event.eventAfterCut.chosenSlice.completeness);
                } else if(event.eventBeforeCut.sliceVec.size() > 1){
                    sliceCompletenessSlices.dune->Fill(event.eventAfterCut.chosenSlice.completeness);
                }
            } else if(event.eventBeforeCut.DLCurrent == 2){
                sizeCurrent++;
                sliceCompleteness.current->Fill(event.eventAfterCut.chosenSlice.completeness);
                sliceScore.current->Fill(event.eventAfterCut.chosenSlice.score);
                if(event.eventBeforeCut.sliceVec.size() == 1){
                    sizeCurrent1Slice++;
                    sliceCompleteness1Slice.current->Fill(event.eventAfterCut.chosenSlice.completeness);
                } else if(event.eventBeforeCut.sliceVec.size() > 1){
                    sliceCompletenessSlices.current->Fill(event.eventAfterCut.chosenSlice.completeness);
                }
            } else if(event.eventBeforeCut.DLCurrent == 3){
                sizeCheated++;
                sliceCompleteness.cheated->Fill(event.eventAfterCut.chosenSlice.completeness);
                sliceScore.cheated->Fill(event.eventAfterCut.chosenSlice.score);
                if(event.eventBeforeCut.sliceVec.size() == 1){
                    sizeCheated1Slice++;
                    sliceCompleteness1Slice.cheated->Fill(event.eventAfterCut.chosenSlice.completeness);
                } else if(event.eventBeforeCut.sliceVec.size() > 1){
                    sliceCompletenessSlices.cheated->Fill(event.eventAfterCut.chosenSlice.completeness);
                }
            }
        }
    }

    styleDraw(sliceCompleteness.canvas, sliceCompleteness.current, sliceCompleteness.cheated, sliceCompleteness.dune, sliceCompleteness.uboone, 0, 1000, 0.8, 1.005, "/nashome/c/coackley/nuEPlots/chosenSliceCompleteness_dist_zoomed.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceCompleteness1Slice.canvas, sliceCompleteness1Slice.current, sliceCompleteness1Slice.cheated, sliceCompleteness1Slice.dune, sliceCompleteness1Slice.uboone, 0, 1000, 0.8, 1.005, "/nashome/c/coackley/nuEPlots/chosenSliceCompleteness1Slice_dist_zoomed.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceCompletenessSlices.canvas, sliceCompletenessSlices.current, sliceCompletenessSlices.cheated, sliceCompletenessSlices.dune, sliceCompletenessSlices.uboone, 0, 1000, 0, 1.05, "/nashome/c/coackley/nuEPlots/chosenSliceCompletenessSlices_dist.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceScore.canvas, sliceScore.current, sliceScore.cheated, sliceScore.dune, sliceScore.uboone, 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/chosenSliceScore_dist_bigBins.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    
    percentage(sliceCompleteness.current, sliceCompleteness.cheated, sliceCompleteness.dune, sliceCompleteness.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 100, 0.8, 1.005, "/nashome/c/coackley/nuEPlots/chosenSliceCompleteness_perc_zoomed.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceCompleteness1Slice.current, sliceCompleteness1Slice.cheated, sliceCompleteness1Slice.dune, sliceCompleteness1Slice.uboone, sizeCurrent1Slice, sizeCheated1Slice, sizeDLDune1Slice, sizeDLUboone1Slice, 0, 100, 0.8, 1.005, "/nashome/c/coackley/nuEPlots/chosenSliceCompleteness1Slice_perc_zoomed.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceCompletenessSlices.current, sliceCompletenessSlices.cheated, sliceCompletenessSlices.dune, sliceCompletenessSlices.uboone, sizeCurrent - sizeCurrent1Slice, sizeCheated - sizeCheated1Slice, sizeDLDune - sizeDLDune1Slice, sizeDLUboone - sizeDLUboone1Slice, 0, 100, 0, 1.05, "/nashome/c/coackley/nuEPlots/chosenSliceCompletenessSlices_perc.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceScore.current, sliceScore.cheated, sliceScore.dune, sliceScore.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 20, 999, 999, "/nashome/c/coackley/nuEPlots/chosenSliceScore_perc_bigBins.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
}

void deltaVertex(std::vector<allEventData> allEventsData){
    auto deltaX = createHistGroup("deltaX", "#Deltax Distribution: Nu + E Elastic Scattering Events", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaY = createHistGroup("deltaY", "#Deltay Distribution: Nu + E Elastic Scattering Events", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZ = createHistGroup("deltaZ", "#Deltaz Distribution: Nu + E Elastic Scattering Events", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaR = createHistGroup("deltaR", "#Delta#bar{r} Distribution: Nu + E Elastic Scattering Events", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    
    int sizeCurrent = 0;
    int sizeCheated = 0;
    int sizeDune = 0;
    int sizeUboone = 0;

    for(const auto& event : allEventsData){
        if(!(event.eventBeforeCut.trueNeutrinoVec[0].pdg == -999999 || event.eventBeforeCut.trueNeutrinoVec[0].tpcValid == 0 || (event.eventBeforeCut.recoNeutrinoVec.size() == 1 &&  event.eventBeforeCut.recoNeutrinoVec[0].pdg == -999999))){
             double deltaXValue = event.eventAfterCut.chosenRecoNeutrino.vx - event.eventAfterCut.chosenTrueNeutrino.vx;
             double deltaYValue = event.eventAfterCut.chosenRecoNeutrino.vy - event.eventAfterCut.chosenTrueNeutrino.vy;
             double deltaZValue = event.eventAfterCut.chosenRecoNeutrino.vz - event.eventAfterCut.chosenTrueNeutrino.vz;
             double deltaRValue = std::sqrt((deltaXValue * deltaXValue) + (deltaYValue * deltaYValue) + (deltaZValue * deltaZValue));

             //if(deltaZValue > 10){
             //if(event.event == 44){
                 //printf("Event: %f, Run %f, SubRun: %f\n", event.event, event.run, event.subrun);
                 //printf("True Vertex: (%f, %f, %f), Reco Vertex: (%f, %f, %f), Differences: (%f, %f, %f)\n", event.eventAfterCut.chosenTrueNeutrino.vx, event.eventAfterCut.chosenTrueNeutrino.vy, event.eventAfterCut.chosenTrueNeutrino.vz, event.eventAfterCut.chosenRecoNeutrino.vx, event.eventAfterCut.chosenRecoNeutrino.vy, event.eventAfterCut.chosenRecoNeutrino.vz, deltaXValue, deltaYValue, deltaZValue);
             //}

             if(event.eventBeforeCut.DLCurrent == 0){
                 deltaX.uboone->Fill(deltaXValue);
                 deltaY.uboone->Fill(deltaYValue);
                 deltaZ.uboone->Fill(deltaZValue);
                 deltaR.uboone->Fill(deltaRValue); 
                 sizeUboone++;
             } else if(event.eventBeforeCut.DLCurrent == 1){
                 deltaX.dune->Fill(deltaXValue);
                 deltaY.dune->Fill(deltaYValue);
                 deltaZ.dune->Fill(deltaZValue);
                 deltaR.dune->Fill(deltaRValue);
                 sizeDune++;
             } else if(event.eventBeforeCut.DLCurrent == 2){
                 deltaX.current->Fill(deltaXValue);
                 deltaY.current->Fill(deltaYValue);
                 deltaZ.current->Fill(deltaZValue);
                 deltaR.current->Fill(deltaRValue);
                 sizeCurrent++;
             } else if(event.eventBeforeCut.DLCurrent == 3){
                 deltaX.cheated->Fill(deltaXValue);
                 deltaY.cheated->Fill(deltaYValue);
                 deltaZ.cheated->Fill(deltaZValue);
                 deltaR.cheated->Fill(deltaRValue);
                 sizeCheated++;
             }

        } else{
            //std::cout << "true neutrino: pdg = " << event.eventBeforeCut.trueNeutrinoVec[0].pdg << ", valid in tpc = " << event.eventBeforeCut.trueNeutrinoVec[0].tpcValid << std::endl;
            //std::cout << "reco neutrino: vec size = " << event.eventBeforeCut.recoNeutrinoVec.size() << ", pdg num = " << event.eventBeforeCut.recoNeutrinoVec[0].pdg << std::endl;
        }
    }

    styleDraw(deltaX.canvas, deltaX.current, deltaX.cheated, deltaX.dune, deltaX.uboone, 0, 500, 999, 999, "/nashome/c/coackley/nuEPlots/deltaX_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(deltaY.canvas, deltaY.current, deltaY.cheated, deltaY.dune, deltaY.uboone, 0, 500, 999, 999, "/nashome/c/coackley/nuEPlots/deltaY_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(deltaZ.canvas, deltaZ.current, deltaZ.cheated, deltaZ.dune, deltaZ.uboone, 0, 500, 999, 999, "/nashome/c/coackley/nuEPlots/deltaZ_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(deltaR.canvas, deltaR.current, deltaR.cheated, deltaR.dune, deltaR.uboone, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/deltaR_dist.pdf", 0.56, 0.88, 0.70, 0.86);
              
    percentage(deltaX.current, deltaX.cheated, deltaX.dune, deltaX.uboone, sizeCurrent, sizeCheated, sizeDune, sizeUboone, 0, 50, 999, 999, "/nashome/c/coackley/nuEPlots/deltaX_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(deltaY.current, deltaY.cheated, deltaY.dune, deltaY.uboone, sizeCurrent, sizeCheated, sizeDune, sizeUboone, 0, 50, 999, 999, "/nashome/c/coackley/nuEPlots/deltaY_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(deltaZ.current, deltaZ.cheated, deltaZ.dune, deltaZ.uboone, sizeCurrent, sizeCheated, sizeDune, sizeUboone, 0, 50, 999, 999, "/nashome/c/coackley/nuEPlots/deltaZ_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(deltaR.current, deltaR.cheated, deltaR.dune, deltaR.uboone, sizeCurrent, sizeCheated, sizeDune, sizeUboone, 0, 90, 999, 999, "/nashome/c/coackley/nuEPlots/deltaR_perc.pdf", 0.56, 0.88, 0.70, 0.86);
}

void energyAngle(std::vector<allEventData> allEventsData){
    auto deltaEAllEnergy = createHistGroup("deltaE All Energy", "#DeltaE Distribution for E_{reco} Being Sum of PFPs Energies", "E_{reco} - E_{true} (MeV)", 35, -600, 100);
    auto EDividedAllEnergy = createHistGroup("ERecoETrue All Energy", "E_{reco}/E_{true} for E_{reco} Being Sum of PFPs Energies", "E_{reco}/E_{true}", 24, 0, 1.2);
    auto deltaEChosenShower = createHistGroup("deltaE Chosen Shower", "#DeltaE Distribution for E_{reco} Being Chosen Shower Energy", "E_{reco} - E_{true} (MeV)", 31, -600, 20);
    auto EDividedChosenShower = createHistGroup("ERecoETrue Chosen Shoert", "E_{reco}/E_{true} for E_{reco} Being Chosen Shower Energy", "E_{reco}/E_{true}", 24, 0, 1.2);
    auto EChosenShowerEAllEnergyDivided = createHistGroup("ERecoEReco All Energy and Chosen Shower", "Ratio of the Energy of the Chosen Shower and the Summed PFP Energies", "E_{reco, chosen shower}/E_{reco, PFP energies}", 24, 0, 1.2);
    auto ETrueThetaTrue2 = createHistGroup("ETrueThetaTrue2", "E_{true}#theta_{true}^{2}", "E_{true}#theta_{true}^{2} (MeV)", 32, 0, 4.088);
    auto ETrueThetaReco2 = createHistGroup("ETrueThetaReco2", "E_{true}#theta_{reco}^{2}", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoThetaTrue2AllEnergy = createHistGroup("ERecoThetaTrue All", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of PFPs Energies", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoThetaTrue2ChosenShower = createHistGroup("ERecoThetaTrue Chosen Shower", "E_{reco}#theta_{true}^{2} for E_{reco} Being Chosen Shower Energy", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoThetaReco2AllEnergy = createHistGroup("ERecoThetaReco All", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of PFPs Energies", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoThetaReco2ChosenShower = createHistGroup("ERecoThetaReco Chosen", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Chosen Shower Energy", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto angleBetweenDirections = createHistGroup("AngleDifference", "Angle Between the Direction of the Chosen PFP and Direction of the True Particle", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto ETrueThetaTrue2BeforeCuts = createHistGroup("ETrueThetaTrue2 Before Cuts", "E_{true}#theta_{true}^{2} (Before Cuts)", "E_{true}#theta_{true}^{2} (MeV)", 32, 0, 4.088);

    double sizeCurrentBefore = 0;
    double sizeCheatedBefore = 0;
    double sizeDLDuneBefore = 0;
    double sizeDLUbooneBefore = 0;
    double sizeCurrent = 0;
    double sizeCheated = 0;
    double sizeDLDune = 0;
    double sizeDLUboone = 0;

    for(const auto& event : allEventsData){
        double deltaEAllEnergyValue = event.eventAfterCut.totalPFPEnergy - event.eventAfterCut.chosenTrueParticle.energy;
        double deltaEChosenShowerValue = event.eventAfterCut.chosenRecoParticle.bestPlaneEnergy - event.eventAfterCut.chosenTrueParticle.energy;
        double EDividedAllEnergyValue = event.eventAfterCut.totalPFPEnergy/event.eventAfterCut.chosenTrueParticle.energy;
        double EDividedChosenShowerValue = event.eventAfterCut.chosenRecoParticle.bestPlaneEnergy/event.eventAfterCut.chosenTrueParticle.energy;
        double EChosenShowerEAllEnergyDividedValue = event.eventAfterCut.chosenRecoParticle.bestPlaneEnergy/event.eventAfterCut.totalPFPEnergy;
        double ETrueThetaTrue2Value = event.eventAfterCut.chosenTrueParticle.energy * event.eventAfterCut.chosenTrueParticle.angle * event.eventAfterCut.chosenTrueParticle.angle;
        double ETrueThetaReco2Value = event.eventAfterCut.chosenTrueParticle.energy * event.eventAfterCut.chosenRecoParticle.theta * event.eventAfterCut.chosenRecoParticle.theta;
        double ERecoThetaTrue2AllEnergyValue = event.eventAfterCut.totalPFPEnergy * event.eventAfterCut.chosenTrueParticle.angle * event.eventAfterCut.chosenTrueParticle.angle;
        double ERecoThetaTrue2ChosenShowerValue = event.eventAfterCut.chosenRecoParticle.bestPlaneEnergy * event.eventAfterCut.chosenTrueParticle.angle * event.eventAfterCut.chosenTrueParticle.angle;
        double ERecoThetaReco2AllEnergyValue = event.eventAfterCut.totalPFPEnergy * event.eventAfterCut.chosenRecoParticle.theta * event.eventAfterCut.chosenRecoParticle.theta;
        double ERecoThetaReco2ChosenShowerValue = event.eventAfterCut.chosenRecoParticle.bestPlaneEnergy * event.eventAfterCut.chosenRecoParticle.theta * event.eventAfterCut.chosenRecoParticle.theta;
        
        if(!(event.eventBeforeCut.trueNeutrinoVec[0].pdg == -999999 || event.eventBeforeCut.trueNeutrinoVec[0].tpcValid == 0 || event.eventBeforeCut.trueParticleVec[0].pdg == -999999)){
            // Events with a true neutrino within the TPC volume and a true particle
            if(event.eventBeforeCut.DLCurrent == 0){
                ETrueThetaTrue2BeforeCuts.uboone->Fill(ETrueThetaTrue2Value);
                sizeDLUbooneBefore++;
            } else if(event.eventBeforeCut.DLCurrent == 1){
                ETrueThetaTrue2BeforeCuts.dune->Fill(ETrueThetaTrue2Value);
                sizeDLDuneBefore++;
            } else if(event.eventBeforeCut.DLCurrent == 2){
                ETrueThetaTrue2BeforeCuts.current->Fill(ETrueThetaTrue2Value);
                sizeCurrentBefore++;
            } else if(event.eventBeforeCut.DLCurrent == 3){
                ETrueThetaTrue2BeforeCuts.cheated->Fill(ETrueThetaTrue2Value);
                sizeCheatedBefore++;
            }
        }
        
        // Cut: Must be a true neutrino within the TPC volume, must be at least one true particle, must be at least one PFP, must be at least one reconstructed neutrino, the summed PFP energies > 0
        if(!(event.eventBeforeCut.trueNeutrinoVec[0].pdg == -999999 || event.eventBeforeCut.trueNeutrinoVec[0].tpcValid == 0 || event.eventBeforeCut.trueParticleVec[0].pdg == -999999 || (event.eventBeforeCut.recoParticleVec.size() == 1 && event.eventBeforeCut.recoParticleVec[0].pdg == -999999) || (event.eventBeforeCut.recoNeutrinoVec.size() == 1 && event.eventBeforeCut.recoNeutrinoVec[0].pdg == -999999) || event.eventAfterCut.totalPFPEnergy == 0)){

            double aDOTb = ((event.eventAfterCut.chosenRecoParticle.dx * event.eventAfterCut.chosenTrueParticle.dx) + (event.eventAfterCut.chosenRecoParticle.dy * event.eventAfterCut.chosenTrueParticle.dy) + (event.eventAfterCut.chosenRecoParticle.dz * event.eventAfterCut.chosenTrueParticle.dz));
            double aMag = std::sqrt((event.eventAfterCut.chosenRecoParticle.dx * event.eventAfterCut.chosenRecoParticle.dx) + (event.eventAfterCut.chosenRecoParticle.dy * event.eventAfterCut.chosenRecoParticle.dy) + (event.eventAfterCut.chosenRecoParticle.dz * event.eventAfterCut.chosenRecoParticle.dz));
            double bMag = std::sqrt((event.eventAfterCut.chosenTrueParticle.dx * event.eventAfterCut.chosenTrueParticle.dx) + (event.eventAfterCut.chosenTrueParticle.dy * event.eventAfterCut.chosenTrueParticle.dy) + (event.eventAfterCut.chosenTrueParticle.dz * event.eventAfterCut.chosenTrueParticle.dz));
            double cosAngle = aDOTb/(aMag * bMag);
            if(cosAngle > 1.0) cosAngle = 1.0;
            if(cosAngle < -1.0) cosAngle = -1.0;
            double angleDiff = TMath::ACos(cosAngle) * TMath::RadToDeg();

            //if(angleDiff > 150){
            //if(event.event == 44){
                //std::cout << "Look" << std::endl;
                //printf("Event: %f, Run %f, SubRun: %f\n", event.event, event.run, event.subrun);
                //printf("(%f, %f, %f), (%f, %f, %f) -> angle = %f degrees\n", event.eventAfterCut.chosenRecoParticle.dx, event.eventAfterCut.chosenRecoParticle.dy, event.eventAfterCut.chosenRecoParticle.dz, event.eventAfterCut.chosenTrueParticle.dx, event.eventAfterCut.chosenTrueParticle.dy, event.eventAfterCut.chosenTrueParticle.dz, angleDiff);
            //}

            if(event.eventBeforeCut.DLCurrent == 0){
                sizeDLUboone++;
                deltaEAllEnergy.uboone->Fill(deltaEAllEnergyValue);
                EDividedAllEnergy.uboone->Fill(EDividedAllEnergyValue);
                deltaEChosenShower.uboone->Fill(deltaEChosenShowerValue);
                EDividedChosenShower.uboone->Fill(EDividedChosenShowerValue);
                EChosenShowerEAllEnergyDivided.uboone->Fill(EChosenShowerEAllEnergyDividedValue);
                ETrueThetaTrue2.uboone->Fill(ETrueThetaTrue2Value);
                ETrueThetaReco2.uboone->Fill(ETrueThetaReco2Value);
                ERecoThetaTrue2AllEnergy.uboone->Fill(ERecoThetaTrue2AllEnergyValue);
                ERecoThetaTrue2ChosenShower.uboone->Fill(ERecoThetaTrue2ChosenShowerValue);
                ERecoThetaReco2AllEnergy.uboone->Fill(ERecoThetaReco2AllEnergyValue);
                ERecoThetaReco2ChosenShower.uboone->Fill(ERecoThetaReco2ChosenShowerValue);
                angleBetweenDirections.uboone->Fill(angleDiff);

            } else if(event.eventBeforeCut.DLCurrent == 1){
                sizeDLDune++;
                deltaEAllEnergy.dune->Fill(deltaEAllEnergyValue);
                EDividedAllEnergy.dune->Fill(EDividedAllEnergyValue);
                deltaEChosenShower.dune->Fill(deltaEChosenShowerValue);
                EDividedChosenShower.dune->Fill(EDividedChosenShowerValue);
                EChosenShowerEAllEnergyDivided.dune->Fill(EChosenShowerEAllEnergyDividedValue);
                ETrueThetaTrue2.dune->Fill(ETrueThetaTrue2Value);
                ETrueThetaReco2.dune->Fill(ETrueThetaReco2Value);
                ERecoThetaTrue2AllEnergy.dune->Fill(ERecoThetaTrue2AllEnergyValue);
                ERecoThetaTrue2ChosenShower.dune->Fill(ERecoThetaTrue2ChosenShowerValue);
                ERecoThetaReco2AllEnergy.dune->Fill(ERecoThetaReco2AllEnergyValue);
                ERecoThetaReco2ChosenShower.dune->Fill(ERecoThetaReco2ChosenShowerValue);
                angleBetweenDirections.dune->Fill(angleDiff);

            } else if(event.eventBeforeCut.DLCurrent == 2){
                sizeCurrent++;
                deltaEAllEnergy.current->Fill(deltaEAllEnergyValue);
                EDividedAllEnergy.current->Fill(EDividedAllEnergyValue);
                deltaEChosenShower.current->Fill(deltaEChosenShowerValue);
                EDividedChosenShower.current->Fill(EDividedChosenShowerValue);
                EChosenShowerEAllEnergyDivided.current->Fill(EChosenShowerEAllEnergyDividedValue);
                ETrueThetaTrue2.current->Fill(ETrueThetaTrue2Value);
                ETrueThetaReco2.current->Fill(ETrueThetaReco2Value);
                ERecoThetaTrue2AllEnergy.current->Fill(ERecoThetaTrue2AllEnergyValue);
                ERecoThetaTrue2ChosenShower.current->Fill(ERecoThetaTrue2ChosenShowerValue);
                ERecoThetaReco2AllEnergy.current->Fill(ERecoThetaReco2AllEnergyValue);
                ERecoThetaReco2ChosenShower.current->Fill(ERecoThetaReco2ChosenShowerValue);
                angleBetweenDirections.current->Fill(angleDiff);

            } else if(event.eventBeforeCut.DLCurrent == 3){
                sizeCheated++;
                deltaEAllEnergy.cheated->Fill(deltaEAllEnergyValue);
                EDividedAllEnergy.cheated->Fill(EDividedAllEnergyValue);
                deltaEChosenShower.cheated->Fill(deltaEChosenShowerValue);
                EDividedChosenShower.cheated->Fill(EDividedChosenShowerValue);
                EChosenShowerEAllEnergyDivided.cheated->Fill(EChosenShowerEAllEnergyDividedValue);
                ETrueThetaTrue2.cheated->Fill(ETrueThetaTrue2Value);
                ETrueThetaReco2.cheated->Fill(ETrueThetaReco2Value);
                ERecoThetaTrue2AllEnergy.cheated->Fill(ERecoThetaTrue2AllEnergyValue);
                ERecoThetaTrue2ChosenShower.cheated->Fill(ERecoThetaTrue2ChosenShowerValue);
                ERecoThetaReco2AllEnergy.cheated->Fill(ERecoThetaReco2AllEnergyValue);
                ERecoThetaReco2ChosenShower.cheated->Fill(ERecoThetaReco2ChosenShowerValue);
                angleBetweenDirections.cheated->Fill(angleDiff);
            }
        }
    }

    int drawLine = 1;

    styleDraw(deltaEAllEnergy.canvas, deltaEAllEnergy.current, deltaEAllEnergy.cheated, deltaEAllEnergy.dune, deltaEAllEnergy.uboone, 0, 200, 999, 999, "/nashome/c/coackley/nuEPlots/deltaEAllEnergy_dist.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(EDividedAllEnergy.canvas, EDividedAllEnergy.current, EDividedAllEnergy.cheated, EDividedAllEnergy.dune, EDividedAllEnergy.uboone, 0, 350, 999, 999, "/nashome/c/coackley/nuEPlots/EDividedAllEnergy_dist.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(deltaEChosenShower.canvas, deltaEChosenShower.current, deltaEChosenShower.cheated, deltaEChosenShower.dune, deltaEChosenShower.uboone, 0, 200, 999, 999, "/nashome/c/coackley/nuEPlots/deltaEChosenShower_dist.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(EDividedChosenShower.canvas, EDividedChosenShower.current, EDividedChosenShower.cheated, EDividedChosenShower.dune, EDividedChosenShower.uboone, 0, 350, 999, 999, "/nashome/c/coackley/nuEPlots/EDividedChosenShower_dist.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(EChosenShowerEAllEnergyDivided.canvas, EChosenShowerEAllEnergyDivided.current, EChosenShowerEAllEnergyDivided.cheated, EChosenShowerEAllEnergyDivided.dune, EChosenShowerEAllEnergyDivided.uboone, 0, 900, 999, 999, "/nashome/c/coackley/nuEPlots/EChosenShowerEAllEnergyDivided_dist.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ETrueThetaTrue2.canvas, ETrueThetaTrue2.current, ETrueThetaTrue2.cheated, ETrueThetaTrue2.dune, ETrueThetaTrue2.uboone, 0, 150, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaTrue2_dist.pdf", 0.56, 0.88, 0.70, 0.86, nullptr, nullptr, &drawLine);
    styleDraw(ETrueThetaReco2.canvas, ETrueThetaReco2.current, ETrueThetaReco2.cheated, ETrueThetaReco2.dune, ETrueThetaReco2.uboone, 0, 150, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaReco2_dist.pdf", 0.56, 0.88, 0.70, 0.86, nullptr, nullptr, &drawLine);
    styleDraw(ERecoThetaTrue2AllEnergy.canvas, ERecoThetaTrue2AllEnergy.current, ERecoThetaTrue2AllEnergy.cheated, ERecoThetaTrue2AllEnergy.dune, ERecoThetaTrue2AllEnergy.uboone, 0, 200, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaTrue2AllEnergy_dist.pdf", 0.56, 0.88, 0.70, 0.86, nullptr, nullptr, &drawLine);
    styleDraw(ERecoThetaTrue2ChosenShower.canvas, ERecoThetaTrue2ChosenShower.current, ERecoThetaTrue2ChosenShower.cheated, ERecoThetaTrue2ChosenShower.dune, ERecoThetaTrue2ChosenShower.uboone, 0, 200, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaTrue2ChosenShower_dist.pdf", 0.56, 0.88, 0.70, 0.86, nullptr, nullptr, &drawLine);
    styleDraw(ERecoThetaReco2AllEnergy.canvas, ERecoThetaReco2AllEnergy.current, ERecoThetaReco2AllEnergy.cheated, ERecoThetaReco2AllEnergy.dune, ERecoThetaReco2AllEnergy.uboone, 0, 200, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaReco2AllEnergy_dist.pdf", 0.56, 0.88, 0.70, 0.86, nullptr, nullptr, &drawLine);
    styleDraw(ERecoThetaReco2ChosenShower.canvas, ERecoThetaReco2ChosenShower.current, ERecoThetaReco2ChosenShower.cheated, ERecoThetaReco2ChosenShower.dune, ERecoThetaReco2ChosenShower.uboone, 0, 200, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaReco2ChosenShower_dist.pdf", 0.56, 0.88, 0.70, 0.86, nullptr, nullptr, &drawLine);
    styleDraw(angleBetweenDirections.canvas, angleBetweenDirections.current, angleBetweenDirections.cheated, angleBetweenDirections.dune, angleBetweenDirections.uboone, 0, 250, 999, 999, "/nashome/c/coackley/nuEPlots/angleDiffDirections_dist.pdf", 0.56, 0.88, 0.70, 0.86);
    styleDraw(ETrueThetaTrue2BeforeCuts.canvas, ETrueThetaTrue2BeforeCuts.current, ETrueThetaTrue2BeforeCuts.cheated, ETrueThetaTrue2BeforeCuts.dune, ETrueThetaTrue2BeforeCuts.uboone, 0, 150, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaTrue2BeforeCuts_dist.pdf", 0.56, 0.88, 0.70, 0.86);

    percentage(deltaEAllEnergy.current, deltaEAllEnergy.cheated, deltaEAllEnergy.dune, deltaEAllEnergy.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 20, 999, 999, "/nashome/c/coackley/nuEPlots/deltaEAllEnergy_perc.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(EDividedAllEnergy.current, EDividedAllEnergy.cheated, EDividedAllEnergy.dune, EDividedAllEnergy.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 40, 999, 999, "/nashome/c/coackley/nuEPlots/EDividedAllEnergy_perc.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(deltaEChosenShower.current, deltaEChosenShower.cheated, deltaEChosenShower.dune, deltaEChosenShower.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 20, 999, 999, "/nashome/c/coackley/nuEPlots/deltaEChosenShower_perc.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(EDividedChosenShower.current, EDividedChosenShower.cheated, EDividedChosenShower.dune, EDividedChosenShower.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 40, 999, 999, "/nashome/c/coackley/nuEPlots/EDividedChosenShower_perc.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(EChosenShowerEAllEnergyDivided.current, EChosenShowerEAllEnergyDivided.cheated, EChosenShowerEAllEnergyDivided.dune, EChosenShowerEAllEnergyDivided.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 100, 999, 999, "/nashome/c/coackley/nuEPlots/EChosenShowerEAllEnergyDivided_perc.pdf", 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ETrueThetaTrue2.current, ETrueThetaTrue2.cheated, ETrueThetaTrue2.dune, ETrueThetaTrue2.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 15, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaTrue2_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);
    percentage(ETrueThetaReco2.current, ETrueThetaReco2.cheated, ETrueThetaReco2.dune, ETrueThetaReco2.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 15, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaReco2_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);
    percentage(ERecoThetaTrue2AllEnergy.current, ERecoThetaTrue2AllEnergy.cheated, ERecoThetaTrue2AllEnergy.dune, ERecoThetaTrue2AllEnergy.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 20, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaTrue2AllEnergy_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);
    percentage(ERecoThetaTrue2ChosenShower.current, ERecoThetaTrue2ChosenShower.cheated, ERecoThetaTrue2ChosenShower.dune, ERecoThetaTrue2ChosenShower.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 20, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaTrue2ChosenShower_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);
    percentage(ERecoThetaReco2AllEnergy.current, ERecoThetaReco2AllEnergy.cheated, ERecoThetaReco2AllEnergy.dune, ERecoThetaReco2AllEnergy.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 25, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaReco2AllEnergy_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);
    percentage(ERecoThetaReco2ChosenShower.current, ERecoThetaReco2ChosenShower.cheated, ERecoThetaReco2ChosenShower.dune, ERecoThetaReco2ChosenShower.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 25, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaReco2ChosenShower_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);
    percentage(angleBetweenDirections.current, angleBetweenDirections.cheated, angleBetweenDirections.dune, angleBetweenDirections.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 30, 999, 999, "/nashome/c/coackley/nuEPlots/angleDiffDirections_perc.pdf", 0.56, 0.88, 0.70, 0.86);
    percentage(ETrueThetaTrue2BeforeCuts.current, ETrueThetaTrue2BeforeCuts.cheated, ETrueThetaTrue2BeforeCuts.dune, ETrueThetaTrue2BeforeCuts.uboone, sizeCurrentBefore, sizeCheatedBefore, sizeDLDuneBefore, sizeDLUbooneBefore, 0, 15, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaTrue2BeforeCuts_perc.pdf", 0.56, 0.88, 0.70, 0.86, &drawLine);

    efficiency(ETrueThetaTrue2.current, ETrueThetaTrue2.cheated, ETrueThetaTrue2.dune, ETrueThetaTrue2.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaTrue2_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine);
    efficiency(ETrueThetaReco2.current, ETrueThetaReco2.cheated, ETrueThetaReco2.dune, ETrueThetaReco2.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaReco2_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine);
    efficiency(ERecoThetaTrue2AllEnergy.current, ERecoThetaTrue2AllEnergy.cheated, ERecoThetaTrue2AllEnergy.dune, ERecoThetaTrue2AllEnergy.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaTrue2AllEnergy_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine); 
    efficiency(ERecoThetaTrue2ChosenShower.current, ERecoThetaTrue2ChosenShower.cheated, ERecoThetaTrue2ChosenShower.dune, ERecoThetaTrue2ChosenShower.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaTrue2ChosenShower_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine);
    efficiency(ERecoThetaReco2AllEnergy.current, ERecoThetaReco2AllEnergy.cheated, ERecoThetaReco2AllEnergy.dune, ERecoThetaReco2AllEnergy.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaReco2AllEnergy_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine);
    efficiency(ERecoThetaReco2ChosenShower.current, ERecoThetaReco2ChosenShower.cheated, ERecoThetaReco2ChosenShower.dune, ERecoThetaReco2ChosenShower.uboone, sizeCurrent, sizeCheated, sizeDLDune, sizeDLUboone, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ERecoThetaReco2ChosenShower_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine);
    efficiency(ETrueThetaTrue2BeforeCuts.current, ETrueThetaTrue2BeforeCuts.cheated, ETrueThetaTrue2BeforeCuts.dune, ETrueThetaTrue2BeforeCuts.uboone, sizeCurrentBefore, sizeCheatedBefore, sizeDLDuneBefore, sizeDLUbooneBefore, 0, 1, 999, 999, "/nashome/c/coackley/nuEPlots/ETrueThetaTrue2BeforeCuts_eff.pdf", 0.56, 0.88, 0.14, 0.3, &drawLine);

}

void nuE_macro(){

    std::vector<allEventData> eventsBeforeAfterCuts;

    int counterTrueNeut = 0;
    int counterTruePart = 0;
    int counterRecoNeut = 0;
    int counterRecoPart = 0;

    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/merged.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/analysed_DL_uboone/CRUMBS/1.root");
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/analysed_Current/NoRefinement/CRUMBS/1.root");
    if(!file){
        std::cerr << "Error opening the file" << std::endl;
        return;
    }

    TDirectory *dir = (TDirectory*)file->Get("ana");
    if(!dir){
        std::cerr << "Directory 'ana' not found" << std::endl;
        return;
    }

    TTree *tree = (TTree*)dir->Get("NuE");
    if(!tree){
        std::cerr << "NuE not found" << std::endl;
        return;
    }

    UInt_t eventID, runID, subRunID;
    double DLCurrent;
    
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);

    std::vector<double> *truth_neutrinoVX = nullptr;
    std::vector<double> *truth_neutrinoVY = nullptr;
    std::vector<double> *truth_neutrinoVZ = nullptr;
    std::vector<double> *truth_CCNC = nullptr;
    std::vector<double> *truth_neutrinoType = nullptr;
    std::vector<double> *truth_chargedLepton = nullptr;
    std::vector<double> *truth_neutrinoTPCID = nullptr;
    std::vector<double> *truth_neutrinoTPCValid = nullptr;

    std::vector<double> *truth_particlePDG = nullptr;
    std::vector<double> *truth_particleVX = nullptr;
    std::vector<double> *truth_particleVY = nullptr;
    std::vector<double> *truth_particleVZ = nullptr;
    std::vector<double> *truth_particlePX = nullptr;
    std::vector<double> *truth_particlePY = nullptr;
    std::vector<double> *truth_particlePZ = nullptr;
    std::vector<double> *truth_particleEnergy = nullptr;
    std::vector<double> *truth_particleAngle = nullptr;
    std::vector<double> *truth_particleETheta2 = nullptr;
    std::vector<double> *truth_particleDirectionX = nullptr;
    std::vector<double> *truth_particleDirectionY = nullptr;
    std::vector<double> *truth_particleDirectionZ = nullptr;

    std::vector<double> *reco_neutrinoPDG = nullptr;
    std::vector<double> *reco_neutrinoIsPrimary = nullptr;
    std::vector<double> *reco_neutrinoVX = nullptr;
    std::vector<double> *reco_neutrinoVY = nullptr;
    std::vector<double> *reco_neutrinoVZ = nullptr;
    std::vector<double> *reco_neutrinoSliceID = nullptr;

    std::vector<double> *reco_particlePDG = nullptr;
    std::vector<double> *reco_particleIsPrimary = nullptr;
    std::vector<double> *reco_particleVX = nullptr;
    std::vector<double> *reco_particleVY = nullptr;
    std::vector<double> *reco_particleVZ = nullptr;
    std::vector<double> *reco_particleDirectionX = nullptr;
    std::vector<double> *reco_particleDirectionY = nullptr;
    std::vector<double> *reco_particleDirectionZ = nullptr;
    std::vector<double> *reco_particleSliceID = nullptr;
    std::vector<double> *reco_particleBestPlaneEnergy = nullptr;
    std::vector<double> *reco_particleTheta = nullptr;
    std::vector<double> *reco_particleTrackScore = nullptr;
    std::vector<double> *reco_particleCompleteness = nullptr;
    std::vector<double> *reco_particlePurity = nullptr;

    std::vector<double> *reco_sliceID = nullptr;
    std::vector<double> *reco_sliceCompleteness = nullptr;
    std::vector<double> *reco_slicePurity = nullptr;
    std::vector<double> *reco_sliceScore = nullptr;

    tree->SetBranchAddress("truth_neutrinoVX", &truth_neutrinoVX);
    tree->SetBranchAddress("truth_neutrinoVY", &truth_neutrinoVY);
    tree->SetBranchAddress("truth_neutrinoVZ", &truth_neutrinoVZ);
    tree->SetBranchAddress("truth_CCNC", &truth_CCNC);
    tree->SetBranchAddress("truth_neutrinoType", &truth_neutrinoType);
    tree->SetBranchAddress("truth_chargedLepton", &truth_chargedLepton);
    tree->SetBranchAddress("truth_neutrinoTPCID", &truth_neutrinoTPCID);
    tree->SetBranchAddress("truth_neutrinoTPCValid", &truth_neutrinoTPCValid);

    tree->SetBranchAddress("truth_particlePDG", &truth_particlePDG);
    tree->SetBranchAddress("truth_particleVX", &truth_particleVX);
    tree->SetBranchAddress("truth_particleVY", &truth_particleVY);
    tree->SetBranchAddress("truth_particleVZ", &truth_particleVZ);
    tree->SetBranchAddress("truth_particlePX", &truth_particlePX);
    tree->SetBranchAddress("truth_particlePY", &truth_particlePY);
    tree->SetBranchAddress("truth_particlePZ", &truth_particlePZ);
    tree->SetBranchAddress("truth_particleEnergy", &truth_particleEnergy);
    tree->SetBranchAddress("truth_particleAngle", &truth_particleAngle);
    tree->SetBranchAddress("truth_particleETheta2", &truth_particleETheta2);
    tree->SetBranchAddress("truth_particleDirectionX", &truth_particleDirectionX);
    tree->SetBranchAddress("truth_particleDirectionY", &truth_particleDirectionY);
    tree->SetBranchAddress("truth_particleDirectionZ", &truth_particleDirectionZ);

    tree->SetBranchAddress("reco_neutrinoPDG", &reco_neutrinoPDG);
    tree->SetBranchAddress("reco_neutrinoIsPrimary", &reco_neutrinoIsPrimary);
    tree->SetBranchAddress("reco_neutrinoVX", &reco_neutrinoVX);
    tree->SetBranchAddress("reco_neutrinoVY", &reco_neutrinoVY);
    tree->SetBranchAddress("reco_neutrinoVZ", &reco_neutrinoVZ);
    tree->SetBranchAddress("reco_neutrinoSliceID", &reco_neutrinoSliceID);

    tree->SetBranchAddress("reco_particlePDG", &reco_particlePDG);
    tree->SetBranchAddress("reco_particleIsPrimary", &reco_particleIsPrimary);
    tree->SetBranchAddress("reco_particleVX", &reco_particleVX);
    tree->SetBranchAddress("reco_particleVY", &reco_particleVY);
    tree->SetBranchAddress("reco_particleVZ", &reco_particleVZ);
    tree->SetBranchAddress("reco_particleDirectionX", &reco_particleDirectionX);
    tree->SetBranchAddress("reco_particleDirectionY", &reco_particleDirectionY);
    tree->SetBranchAddress("reco_particleDirectionZ", &reco_particleDirectionZ);
    tree->SetBranchAddress("reco_particleSliceID", &reco_particleSliceID);
    tree->SetBranchAddress("reco_particleBestPlaneEnergy", &reco_particleBestPlaneEnergy);
    tree->SetBranchAddress("reco_particleTheta", &reco_particleTheta);
    tree->SetBranchAddress("reco_particleTrackScore", &reco_particleTrackScore);
    tree->SetBranchAddress("reco_particleCompleteness", &reco_particleCompleteness);
    tree->SetBranchAddress("reco_particlePurity", &reco_particlePurity);

    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness", &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_slicePurity", &reco_slicePurity);
    tree->SetBranchAddress("reco_sliceScore", &reco_sliceScore);

    Long64_t numEntries = tree->GetEntries();
    //std::cout << "Number of Events: " << numEntries << std::endl;
    //tree->Print();

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);
        //printf("\n_________________________________________________________________________________\n");
        //printf("Run: %d, SubRun: %d, Event: %d, DL/Current: %f\n", runID, subRunID, eventID, DLCurrent);
        
        eventData eData; // A struct containing vectors for each type of particle, i.e. true or reco...
        eventData eDataAfterCuts;

        std::vector<recoParticle> recoParticles;
        std::vector<trueParticle> trueParticles;
        std::vector<recoNeutrino> recoNeutrinos;
        std::vector<trueNeutrino> trueNeutrinos;
        std::vector<recoSlice> slices;
        //printf("Number of true neutrinos = %zu, reco neutrinos = %zu, true particles = %zu, reco particles = %zu, slices = %zu\n\n", truth_neutrinoVX->size(), reco_neutrinoPDG->size(), truth_particlePDG->size(), reco_particlePDG->size(), reco_sliceID->size());
    
        if(truth_neutrinoVX->size() > 1) counterTrueNeut++;
        if(truth_particlePDG->size() > 1) counterTruePart++;
        if(reco_neutrinoPDG->size() < 1) counterRecoNeut++;
        if(reco_particlePDG->size() > 1) counterRecoPart++;

        for(size_t j = 0; j < reco_particlePDG->size(); ++j){
            recoParticle particle;
            particle.pdg = reco_particlePDG->at(j);
            particle.isPrimary = reco_particleIsPrimary->at(j);
            particle.vx = reco_particleVX->at(j);
            particle.vy = reco_particleVY->at(j);
            particle.vz = reco_particleVZ->at(j);
            particle.dx = reco_particleDirectionX->at(j);
            particle.dy = reco_particleDirectionY->at(j);
            particle.dz = reco_particleDirectionZ->at(j);
            particle.sliceID = reco_particleSliceID->at(j);
            particle.bestPlaneEnergy = reco_particleBestPlaneEnergy->at(j);
            particle.theta = reco_particleTheta->at(j);
            particle.trackscore = reco_particleTrackScore->at(j);
            particle.completeness = reco_particleCompleteness->at(j);
            particle.purity = reco_particlePurity->at(j);
            recoParticles.push_back(particle);
            //printf("Reco Particle %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %f, Energy = %f, Theta = %f, Trackscore = %f, Completeness = %f, Purity = %f\n", j, particle.pdg, particle.isPrimary, particle.vx, particle.vy, particle.vz, particle.dx, particle.dy, particle.dz, particle.sliceID, particle.bestPlaneEnergy, particle.theta, particle.trackscore, particle.completeness, particle.purity);
        }

        //printf("\n");
        for(size_t j = 0; j < reco_neutrinoPDG->size(); ++j){
            recoNeutrino neutrino;
            neutrino.pdg = reco_neutrinoPDG->at(j);
            neutrino.isPrimary = reco_neutrinoIsPrimary->at(j);
            neutrino.vx = reco_neutrinoVX->at(j);
            neutrino.vy = reco_neutrinoVY->at(j);
            neutrino.vz = reco_neutrinoVZ->at(j);
            neutrino.sliceID = reco_neutrinoSliceID->at(j);
            recoNeutrinos.push_back(neutrino);
            //printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, neutrino.pdg, neutrino.isPrimary, neutrino.vx, neutrino.vy, neutrino.vz, neutrino.sliceID);
        }

        //printf("\n");
        for(size_t j = 0; j < truth_particlePDG->size(); ++j){
            trueParticle particle;
            particle.pdg = truth_particlePDG->at(j);
            particle.vx = truth_particleVX->at(j);
            particle.vy = truth_particleVY->at(j);
            particle.vz = truth_particleVZ->at(j);
            particle.px = truth_particlePX->at(j);
            particle.py = truth_particlePY->at(j);
            particle.pz = truth_particlePZ->at(j);
            particle.energy = truth_particleEnergy->at(j);
            particle.angle = truth_particleAngle->at(j);
            particle.ETheta2 = truth_particleETheta2->at(j);
            particle.dx = truth_particleDirectionX->at(j);
            particle.dy = truth_particleDirectionY->at(j);
            particle.dz = truth_particleDirectionZ->at(j);
            trueParticles.push_back(particle);
            //printf("Truth Particle %zu: PDG = %f, Vertex = (%f, %f, %f), Momentum = (%f, %f, %f), Energy = %f, Theta = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", j, particle.pdg, particle.vx, particle.vy, particle.vz, particle.px, particle.py, particle.pz, particle.energy, particle.angle, particle.ETheta2, particle.dx, particle.dy, particle.dz);
        }
    
        //printf("\n");
        for(size_t j = 0; j < truth_neutrinoVX->size(); ++j){
            trueNeutrino neutrino;
            neutrino.vx = truth_neutrinoVX->at(j);
            neutrino.vy = truth_neutrinoVY->at(j);
            neutrino.vz = truth_neutrinoVZ->at(j);
            neutrino.CCNC = truth_CCNC->at(j);
            neutrino.pdg = truth_neutrinoType->at(j);
            neutrino.leptonpdg = truth_chargedLepton->at(j);
            neutrino.tpcID = truth_neutrinoTPCID->at(j);
            neutrino.tpcValid = truth_neutrinoTPCValid->at(j);
            trueNeutrinos.push_back(neutrino);
            //printf("Truth Neutrino %zu: PDG = %f, Vertex = (%f, %f, %f), CCNC = %f, Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", j, neutrino.pdg, neutrino.vx, neutrino.vy, neutrino.vz, neutrino.CCNC, neutrino.leptonpdg, neutrino.tpcID, neutrino.tpcValid);
        }
    
        //printf("\n");
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            recoSlice reconstructedSlice;
            reconstructedSlice.id = reco_sliceID->at(j);
            reconstructedSlice.completeness = reco_sliceCompleteness->at(j);
            reconstructedSlice.purity = reco_slicePurity->at(j);
            reconstructedSlice.score = reco_sliceScore->at(j);
            slices.push_back(reconstructedSlice);
            //printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, Score = %f\n", j, reconstructedSlice.id, reconstructedSlice.completeness, reconstructedSlice.purity, reconstructedSlice.score);
        }

        //printf("\nChosen:\n");
        eData.recoParticleVec = recoParticles;
        eData.recoNeutrinoVec = recoNeutrinos;
        eData.trueNeutrinoVec = trueNeutrinos;
        eData.trueParticleVec = trueParticles;
        eData.sliceVec = slices;

        eventAfterCuts eventCut;

        double totalEnergyPFPs = 0;

        if(recoParticles.size() == 1){
            eventCut.chosenRecoParticle = recoParticles[0];     
            totalEnergyPFPs = recoParticles[0].bestPlaneEnergy;
        } else{
            eventCut.chosenRecoParticle = chooseRecoParticle(recoParticles, totalEnergyPFPs);
        }
        eventCut.totalPFPEnergy = totalEnergyPFPs;

        if(recoNeutrinos.size() == 1){
            eventCut.chosenRecoNeutrino = recoNeutrinos[0];
        } else{
        }

        if(trueNeutrinos.size() == 1){
            eventCut.chosenTrueNeutrino = trueNeutrinos[0];
        }
        
        if(trueParticles.size() == 1){
            eventCut.chosenTrueParticle = trueParticles[0];
        }

        if(slices.size() == 1){
            eventCut.chosenSlice = slices[0];
        } else{
            eventCut.chosenSlice = chooseRecoSlice(slices);
        }

        eventCut.DLCurrent = DLCurrent;
        eData.DLCurrent = DLCurrent;

        //printf("Chosen Reco Particle: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Direction = (%f, %f, %f), Slice ID = %f, Energy = %f, Theta = %f, Trackscore = %f, Completeness = %f, Purity = %f\n", eventCut.chosenRecoParticle.pdg, eventCut.chosenRecoParticle.isPrimary, eventCut.chosenRecoParticle.vx, eventCut.chosenRecoParticle.vy, eventCut.chosenRecoParticle.vz, eventCut.chosenRecoParticle.dx, eventCut.chosenRecoParticle.dy, eventCut.chosenRecoParticle.dz, eventCut.chosenRecoParticle.sliceID, eventCut.chosenRecoParticle.bestPlaneEnergy, eventCut.chosenRecoParticle.theta, eventCut.chosenRecoParticle.trackscore, eventCut.chosenRecoParticle.completeness, eventCut.chosenRecoParticle.purity);
        //printf("Chosen Reco Neutrino: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", eventCut.chosenRecoNeutrino.pdg, eventCut.chosenRecoNeutrino.isPrimary, eventCut.chosenRecoNeutrino.vx, eventCut.chosenRecoNeutrino.vy, eventCut.chosenRecoNeutrino.vz, eventCut.chosenRecoNeutrino.sliceID);
        //printf("Chosen True Neutrino: PDG = %f, Vertex = (%f, %f, %f), CCNC = %f, Charged Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", eventCut.chosenTrueNeutrino.pdg, eventCut.chosenTrueNeutrino.vx, eventCut.chosenTrueNeutrino.vy, eventCut.chosenTrueNeutrino.vz, eventCut.chosenTrueNeutrino.CCNC, eventCut.chosenTrueNeutrino.leptonpdg, eventCut.chosenTrueNeutrino.tpcID, eventCut.chosenTrueNeutrino.tpcValid);
        //printf("Chosen True Particle: PDG = %f, Vertex = (%f, %f, %f), Momentum = (%f, %f, %f), Energy = %f, Theta = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", eventCut.chosenTrueParticle.pdg, eventCut.chosenTrueParticle.vx, eventCut.chosenTrueParticle.vy, eventCut.chosenTrueParticle.vz, eventCut.chosenTrueParticle.px, eventCut.chosenTrueParticle.py, eventCut.chosenTrueParticle.pz, eventCut.chosenTrueParticle.energy, eventCut.chosenTrueParticle.angle, eventCut.chosenTrueParticle.ETheta2, eventCut.chosenTrueParticle.dx, eventCut.chosenTrueParticle.dy, eventCut.chosenTrueParticle.dz);
        //printf("Chosen Slice: ID = %f, Completeness = %f, Purity = %f, Score = %f\n", eventCut.chosenSlice.id, eventCut.chosenSlice.completeness, eventCut.chosenSlice.purity, eventCut.chosenSlice.score);
        //printf("Total Reconstructed Energy = %f\n", totalEnergyPFPs);

        recoParticles.clear();
        recoNeutrinos.clear();
        trueNeutrinos.clear();
        trueParticles.clear();
        slices.clear();
        
        allEventData eventAllDataWithin;
        eventAllDataWithin.eventBeforeCut = eData;
        eventAllDataWithin.eventAfterCut = eventCut;
        eventAllDataWithin.event = eventID;
        eventAllDataWithin.run = runID;
        eventAllDataWithin.subrun = subRunID;
        eventsBeforeAfterCuts.push_back(eventAllDataWithin);
    }

    numberPlots(eventsBeforeAfterCuts);
    chosenSlicePlots(eventsBeforeAfterCuts);
    deltaVertex(eventsBeforeAfterCuts);
    energyAngle(eventsBeforeAfterCuts);
    file->Close();
    //printf("\n __________ Number of events with > 1 True Neutrino = %d, True Particle = %d, Reco Neutrino = %d, Reco Particle = %d __________\n", counterTrueNeut, counterTruePart, counterRecoNeut, counterRecoPart);

}
