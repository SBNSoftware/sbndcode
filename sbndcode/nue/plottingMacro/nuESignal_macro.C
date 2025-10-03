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

struct counter{
    int current = 0;
    int cheated = 0;
    int dune = 0;
    int uboone = 0;
    int sbnd = 0;
    int nue = 0;
    int nue_100k = 0;
};

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* current;
    TH1F* cheated;
    TH1F* dune;
    TH1F* uboone;
    TH1F* sbnd;
    TH1F* nue;
    TH1F* nue_100k;
} histGroup;

struct edgeSlice {
    histGroup hg;
    double low;
    double high;
};

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
    double numHits;
    double numMatchedHits;
    double numTrueHits;
    double truePDG;
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
    double neutrinoParent;
} trueParticle;

typedef struct{
    double id;
    double completeness;
    double purity;
    double score;
} recoSlice;

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
        (TH1F*) base->Clone((baseName + "_dluboone").c_str()),
        (TH1F*) base->Clone((baseName + "_dlsbnd").c_str()),
        (TH1F*) base->Clone((baseName + "_dlnue").c_str()),
        (TH1F*) base->Clone((baseName + "_dlnue_100k").c_str())
    };    
}

void TwoDHistDraw(TH2D* hist, const char* filename, const char* title){
    TCanvas* TwoDHistCanvas = new TCanvas("2dHist_canvas", "Graph Draw Options", 200, 10, 600, 400);
    TwoDHistCanvas->SetTickx();
    TwoDHistCanvas->SetTicky();

    hist->SetTitle(title);
    hist->Draw("COLZ");
    hist->SetStats(0);
    hist->GetXaxis()->SetTickLength(0.04);
    hist->GetYaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTickSize(0.02);
    hist->GetYaxis()->SetTickSize(0.02);
       
    TwoDHistCanvas->SaveAs(filename);
    TwoDHistCanvas->Clear();
}

void ProfileDraw(TProfile* profile, const char* filename, const char* title){
    TCanvas* ProfileCanvas = new TCanvas("Profile_canvas", "Graph Draw Options", 200, 10, 600, 400);
    ProfileCanvas->SetTickx();
    ProfileCanvas->SetTicky();

    //TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    //pt->AddText(Form("Number of Entries: %f", (double)profile->GetEntries()));
    //pt->AddText(Form("Angle Difference Mean: %f", (double)profile->GetMean(2)));
    //pt->AddText(Form("Angle Difference Std Dev: %f", (double)profile->GetStdDev(2)));
    //pt->SetFillColor(kWhite);
    //pt->SetFillStyle(1001);
    //pt->SetBorderSize(0); 
	//pt->Draw();

    profile->SetTitle(title);
    profile->SetErrorOption("i");
    profile->SetLineWidth(2);
    profile->SetMarkerStyle(20);
    profile->SetMarkerSize(0.8);
    profile->SetMarkerColor(kBlack);
    profile->SetLineColor(kBlack);
    profile->Draw();
    profile->GetXaxis()->SetTickLength(0.04);
    profile->GetYaxis()->SetTickLength(0.03);
    profile->GetXaxis()->SetTickSize(0.02);
    profile->GetYaxis()->SetTickSize(0.02);
       
    ProfileCanvas->SaveAs(filename);
    ProfileCanvas->Clear();
}

void makeEqualStatProfile(std::vector<double>& xVals, std::vector<double>& yVals, int nbins, double xMin, double xMax, const char* name, const char* title, const char* filename){
    if (xVals.size() != yVals.size()) {
        std::cerr << "xVals and yVals must be the same size!" << std::endl;
        return;
    }

    // Step 1: Fill a temporary histogram with only values inside [xMin,xMax]
    TH1F *hx = new TH1F("hx","x distribution", 1000, xMin, xMax);
    for (size_t i = 0; i < xVals.size(); i++) {
        if (xVals[i] >= xMin && xVals[i] <= xMax)
            hx->Fill(xVals[i]);
    }

    // Step 2: Compute quantiles in this restricted range
    std::vector<double> q(nbins+1);
    std::vector<double> p(nbins+1);
    for (int i = 0; i <= nbins; i++) {
        p[i] = double(i)/nbins;  // cumulative probabilities: 0, 0.1, ..., 1
    }
    hx->GetQuantiles(nbins+1, &q[0], &p[0]);

    // Step 3: Build profile with variable-width bins
    TProfile *prof = new TProfile(name, title, nbins, &q[0], 0, 180);

    // Step 4: Fill only entries in the range
    for (size_t i = 0; i < xVals.size(); i++) {
        if (xVals[i] >= xMin && xVals[i] <= xMax)
            prof->Fill(xVals[i], yVals[i]);
    }

    // Step 5: Draw and save (your style)
    TCanvas* c = new TCanvas("c","Profile",600,400);
    c->SetTickx();
    c->SetTicky();

    prof->SetErrorOption("i");
    prof->SetLineWidth(2);
    prof->SetMarkerStyle(20);
    prof->SetMarkerSize(0.8);
    prof->SetMarkerColor(kBlack);
    prof->SetLineColor(kBlack);
    gStyle->SetOptStat(0);

    prof->Draw();
    prof->GetXaxis()->SetTickLength(0.04);
    prof->GetYaxis()->SetTickLength(0.03);
    prof->GetXaxis()->SetTickSize(0.02);
    prof->GetYaxis()->SetTickSize(0.02);

    c->SaveAs(filename);

    // Debug: print bin edges
    std::cout << "Bin edges for " << name << ":\n";
    for (size_t i = 0; i < q.size(); i++) std::cout << q[i] << " ";
    std::cout << std::endl;

    delete hx;
}

void styleDraw(TCanvas* canvas, TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, TH1F* sbnd, TH1F* nue, TH1F* nue_100k, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, TPaveText* pt = nullptr, int* percentage = nullptr, int* drawLine = nullptr, int* linePos = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    //gStyle->SetPalette(kAvocado);
    //gROOT->ForceStyle();
    //gPad->Update();

    current->SetLineWidth(2);
    //current->SetLineColor(TColor::GetColorPalette(150));
    current->SetLineColor(TColor::GetColor("#e42536"));

    cheated->SetLineWidth(2);
    //cheated->SetLineColor(TColor::GetColorPalette(200));
    cheated->SetLineColor(TColor::GetColor("#f89c20"));
    
    nue->SetLineWidth(2);
    nue->SetLineColor(TColor::GetColor("#7a21dd"));

    nue_100k->SetLineWidth(2);
    nue_100k->SetLineColor(TColor::GetColor("#f89c20"));

    dune->SetLineWidth(2);
    //dune->SetLineColor(TColor::GetColorPalette(50));
    dune->SetLineColor(TColor::GetColor("#7a21dd"));

    uboone->SetLineWidth(2);
    //uboone->SetLineColor(TColor::GetColorPalette(100));
    uboone->SetLineColor(TColor::GetColor("#5790fc"));

    sbnd->SetLineWidth(2);
    sbnd->SetLineColor(TColor::GetColor("#9c9ca1"));

    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({current->GetMaximum(), cheated->GetMaximum(), dune->GetMaximum(), uboone->GetMaximum(), sbnd->GetMaximum(), nue->GetMaximum(), nue_100k->GetMaximum()});

    if((ymin == 999) && (ymax == 999)){
        double yminVal = 0;
        double ymaxVal  = maxYValue * 1.1;
        current->GetYaxis()->SetRangeUser(yminVal, ymaxVal);
    }

    current->Draw("hist");
    //cheated->Draw("histsame");
    //dune->Draw("histsame");
    uboone->Draw("histsame");
    //sbnd->Draw("histsame");
    //nue->Draw("histsame");
    //nue_100k->Draw("histsame");
    
    current->SetStats(0);
    current->GetXaxis()->SetTickLength(0.04);
    current->GetYaxis()->SetTickLength(0.03);
    current->GetXaxis()->SetTickSize(0.02);
    current->GetYaxis()->SetTickSize(0.02);

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(dune, "Pandora Deep Learning: DUNE/LBNF Tune", "f");
    legend->AddEntry(uboone, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(sbnd, "Pandora Deep Learning: SBND/BNB Tune", "f");
    legend->AddEntry(current, "Pandora BDT SBND (without Refinement)", "f");
    //legend->AddEntry(cheated, "Pandora Cheated SBND Vertexing", "f");
    legend->AddEntry(nue, "Pandora Deep Learning: SBND Nu+E Elastic (50k Events)", "f");
    legend->AddEntry(nue_100k, "Pandora Deep Learning: SBND Nu+E Elastic (100k Events)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    //legend->Draw();

    if(drawLine){
        TLine* line = new TLine(1.022, 0, 1.022, current->GetMaximum());
        line->SetLineColor(kGray+2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            latex = new TLatex(1.022 - 0.2, current->GetMaximum() * 0.93, "2m_{e}");
        } else{
            latex = new TLatex(1.022 + 0.1, current->GetMaximum() * 0.93, "2m_{e}");
        }

        latex->SetTextSize(0.035); 
        latex->SetTextAlign(11);
        latex->Draw("same");
    }

    if(pt){
        pt->SetTextSize(legend->GetTextSize());
        pt->SetTextFont(legend->GetTextFont());
        pt->Draw();
    }

    canvas->SaveAs(filename);
}

void percentage(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, TH1F* sbnd, TH1F* nue, TH1F* nue_100k, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double sizeSBND, double sizeNuE, double sizeNuE_100k, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr, int* linePos = nullptr){
    TCanvas *percentageCanvas = new TCanvas("percentage_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentPerc = (TH1F*) current->Clone("perc hist");
    currentPerc->Scale(100.0 * 1.0/sizeCurrent);
    currentPerc->GetYaxis()->SetTitle("Percentage of Events (%)"); 

    TH1F* cheatedPerc = (TH1F*) cheated->Clone("perc hist");
    //cheatedPerc->Scale(100.0 * 1.0/sizeCheated);
    cheatedPerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* dunePerc = (TH1F*) dune->Clone("perc hist");
    //dunePerc->Scale(100.0 * 1.0/sizeDune);
    dunePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* uboonePerc = (TH1F*) uboone->Clone("perc hist");
    uboonePerc->Scale(100.0 * 1.0/sizeUboone);
    uboonePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* sbndPerc = (TH1F*) sbnd->Clone("perc hist");
    sbndPerc->Scale(100.0 * 1.0/sizeSBND);
    sbndPerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TH1F* nuePerc = (TH1F*) nue->Clone("perc hist");
    nuePerc->Scale(100.0 * 1.0/sizeNuE);
    nuePerc->GetYaxis()->SetTitle("Percentage of Events (%)");
    
    TH1F* nue_100kPerc = (TH1F*) nue_100k->Clone("perc hist");
    nue_100kPerc->Scale(100.0 * 1.0/sizeNuE_100k);
    nue_100kPerc->GetYaxis()->SetTitle("Percentage of Events (%)");
    
    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.17, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of DL SBND Entries: %d", (int)sizeSBND));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of DL Nu+E Entries (50k): %d", (int)sizeNuE));
    pt->AddText(Form("Number of DL Nu+E Entries (100k): %d", (int)sizeNuE_100k));
    //pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 
   
    int funcValue = 1;

    styleDraw(percentageCanvas, currentPerc, cheatedPerc, dunePerc, uboonePerc, sbndPerc, nuePerc, nue_100kPerc, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt, &funcValue, drawLine, linePos);
}

void efficiency(TH1F* current, TH1F* cheated, TH1F* dune, TH1F* uboone, TH1F* sbnd, TH1F* nue, TH1F* nue_100k, double sizeCurrent, double sizeCheated, double sizeDune, double sizeUboone, double sizeSBND, double sizeNuE, double sizeNuE_100k, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr, int* linePos = nullptr, std::string xlabel = ""){
    TCanvas *efficiencyCanvas = new TCanvas("efficiency_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentEff = (TH1F*) current->Clone("eff hist");
    currentEff->Reset();
    currentEff->GetYaxis()->SetTitle("Efficiency"); 
    currentEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* cheatedEff = (TH1F*) cheated->Clone("eff hist");
    //cheatedEff->Reset();
    cheatedEff->GetYaxis()->SetTitle("Efficiency");
    cheatedEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* duneEff = (TH1F*) dune->Clone("eff hist");
    //duneEff->Reset();
    duneEff->GetYaxis()->SetTitle("Efficiency");
    duneEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneEff = (TH1F*) uboone->Clone("eff hist");
    ubooneEff->Reset();
    ubooneEff->GetYaxis()->SetTitle("Efficiency");
    ubooneEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* sbndEff = (TH1F*) sbnd->Clone("eff hist");
    sbndEff->Reset();
    sbndEff->GetYaxis()->SetTitle("Efficiency");
    sbndEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    TH1F* nueEff = (TH1F*) nue->Clone("eff hist");
    nueEff->Reset();
    nueEff->GetYaxis()->SetTitle("Efficiency");
    nueEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* nue_100kEff = (TH1F*) nue_100k->Clone("eff hist");
    nue_100kEff->Reset();
    nue_100kEff->GetYaxis()->SetTitle("Efficiency");
    nue_100kEff->GetXaxis()->SetTitle(xlabel.c_str());

    int numBins = current->GetNbinsX();
    double currentSum = 0.0;
    double cheatedSum = 0.0;
    double duneSum = 0.0;
    double ubooneSum = 0.0;
    double sbndSum = 0.0;
    double nueSum = 0.0;
    double nue_100kSum = 0.0;

    for(int i = 1; i <= numBins; ++i){
        currentSum += current->GetBinContent(i);
        //cheatedSum += cheated->GetBinContent(i);
        //duneSum += dune->GetBinContent(i);
        ubooneSum += uboone->GetBinContent(i);
        sbndSum += sbnd->GetBinContent(i);
        nueSum += nue->GetBinContent(i);
        nue_100kSum += nue_100k->GetBinContent(i);

        double currentEffValue = currentSum/sizeCurrent;
        //double cheatedEffValue = cheatedSum/sizeCheated;
        //double duneEffValue = duneSum/sizeDune;
        double ubooneEffValue = ubooneSum/sizeUboone;
        double sbndEffValue = sbndSum/sizeSBND;
        double nueEffValue = nueSum/sizeNuE;
        double nue_100kEffValue = nue_100kSum/sizeNuE_100k;

        currentEff->SetBinContent(i, currentEffValue);
        //cheatedEff->SetBinContent(i, cheatedEffValue);
        //duneEff->SetBinContent(i, duneEffValue);
        ubooneEff->SetBinContent(i, ubooneEffValue);
        sbndEff->SetBinContent(i, sbndEffValue);
        nueEff->SetBinContent(i, nueEffValue);
        nue_100kEff->SetBinContent(i, nue_100kEffValue);
    }

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    //pt->AddText(Form("Number of DL Dune Entries: %d", (int)sizeDune));
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of DL SBND Entries: %d", (int)sizeSBND));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->AddText(Form("Number of Nu+E Elastic Entries (50k): %d", (int)sizeNuE));
    pt->AddText(Form("Number of Nu+E Elastic Entries (100k): %d", (int)sizeNuE_100k));
    
    //pt->AddText(Form("Number of Cheated Entries: %d", (int)sizeCheated));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 

    int funcValue = 1;

    styleDraw(efficiencyCanvas, currentEff, cheatedEff, duneEff, ubooneEff, sbndEff, nueEff, nue_100kEff, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt = nullptr, &funcValue, drawLine, linePos);
}

recoSlice chooseSlice(std::vector<recoSlice> recoSlices, double method){
    int chosenSliceIndex = 0;
    double highestScore = -999999;
    if(method == 0){
        // Choose slice based on completeness
        for(size_t j = 0; j < recoSlices.size(); ++j){
            if(recoSlices[j].completeness > highestScore){
                highestScore = recoSlices[j].completeness;
                chosenSliceIndex = j;
            }
        } 
    } else if(method == 1){
        // Choose sliced based on CRUMBS score
        for(size_t j = 0; j < recoSlices.size(); ++j){
            if(recoSlices[j].score > highestScore){
                highestScore = recoSlices[j].score;
                chosenSliceIndex = j;
            }
        } 
    }
    
    return recoSlices[chosenSliceIndex];
}

recoParticle choosePFP(std::vector<recoParticle> recoParticles, double sliceID, double& totalEnergy, double& numPFPsSlice){
    totalEnergy = 0;
    numPFPsSlice = 0;
    int chosenParticleIndex = 0;
    double highestEnergy = -999999;

    for(size_t j = 0; j < recoParticles.size(); ++j){
        if(recoParticles[j].sliceID == sliceID){
            totalEnergy += recoParticles[j].bestPlaneEnergy;
            numPFPsSlice++;
            if(recoParticles[j].bestPlaneEnergy > highestEnergy){
                highestEnergy = recoParticles[j].bestPlaneEnergy;
                chosenParticleIndex = j;
            }
        }
    }

    return recoParticles[chosenParticleIndex];
}

recoNeutrino chooseRecoNeutrino(std::vector<recoNeutrino> recoNeutrinos, double sliceID){
    int chosenNeutrinoIndex;
    int counter = 0;

    //std::cout << "chosen slice ID: " << sliceID << std::endl;
    //std::cout << "reco Neutrinos size: " << recoNeutrinos.size() << std::endl;
    for(size_t j = 0; j < recoNeutrinos.size(); ++j){
        //std::cout << "reco neutrino " << j << ", slice ID: " << recoNeutrinos[j].sliceID << ", Vertex = (" << recoNeutrinos[j].vx << ", " << recoNeutrinos[j].vy << ", " << recoNeutrinos[j].vz << ")" << std::endl;
        if(recoNeutrinos[j].sliceID == sliceID){
            counter++;
            chosenNeutrinoIndex = j;
        }
    }

    if(counter > 1) printf("%i reco neutrinos in the chosen slice\n", counter);
    return recoNeutrinos[chosenNeutrinoIndex];
}

trueNeutrino chooseTrueNeutrino(std::vector<trueNeutrino> trueNeutrinos){
    // There should just be 1 true neutrino in the event
    //printf("%zu true neutrinos in the event and within the TPC\n", trueNeutrinos.size());
    //printf("True Neutrino Vertex: (%f, %f, %f)\n", trueNeutrinos[0].vx, trueNeutrinos[0].vy, trueNeutrinos[0].vz);
    return trueNeutrinos[0];
}

trueParticle chooseTrueParticle(std::vector<trueParticle> trueParticles, trueNeutrino chosenTrueNeutrino){
    int chosenParticleIndex = 0;
    
    for(size_t j = 0; j < trueParticles.size(); ++j){
        if(trueParticles[j].vx == chosenTrueNeutrino.vx && trueParticles[j].vy == chosenTrueNeutrino.vy && trueParticles[j].vz == chosenTrueNeutrino.vz){
            chosenParticleIndex = j;
        }
    }

    return trueParticles[chosenParticleIndex];
}

double sliceCompletenessCalculator(std::vector<recoParticle> recoParticles, double sliceID){
    int counter = 0;
    double totalNumMatchedHits = 0;
    double totalNumTrueHits = 0;

    for(size_t j = 0; j < recoParticles.size(); ++j){
        if(recoParticles[j].sliceID == sliceID){
            counter++;
            totalNumMatchedHits += recoParticles[j].numMatchedHits;
            totalNumTrueHits += recoParticles[j].numTrueHits;
        }
    }

    double sliceCompleteness;
    if(counter < 1){
        sliceCompleteness = -999999;
    } else {
        sliceCompleteness = (totalNumMatchedHits / totalNumTrueHits);
    }
    //std::cout << "Chosen Slice Completeness = " << sliceCompleteness << std::endl;
    return sliceCompleteness;
}

double slicePurityCalculator(std::vector<recoParticle> recoParticles, double sliceID){
    int counter = 0;
    double totalNumMatchedHits = 0;
    double totalNumHits = 0;

    for(size_t j = 0; j < recoParticles.size(); ++j){
        if(recoParticles[j].sliceID == sliceID){
            counter++;
            totalNumMatchedHits += recoParticles[j].numMatchedHits;
            totalNumHits += recoParticles[j].numHits;
        }
    }

    double slicePurity;
    if(counter < 1){
        slicePurity = -999999;
    } else {
        slicePurity = (totalNumMatchedHits / totalNumHits);
    }
    //std::cout << "Chosen Slice Purity = " << slicePurity << std::endl;
    return slicePurity;
}

void nuESignal_macro(){
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/merged_SCEON_24Sep.root");
    //std::string base_path = "/nashome/c/coackley/nuEPlotsWithCosmicsSCEON_24Sep/";
    
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/22Sep_SCEON/merged_22Sep.root");
    //std::string base_path = "/nashome/c/coackley/nuEPlotsWithCosmics/";
    
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E/analysed_noSCE/merged_noSCE_8Sep.root");
    //std::string base_path = "/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/";
    
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/22Sep_SCEON/merged_22Sep.root");
    std::string base_path = "/nashome/c/coackley/nuEPlotsWithCosmicsSCEON_40kEvents/";
    
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics_v10_09_00/analysed/enuelastic_v10_09_00_gen_g4_detsim_reco1_reco2_analysed_80221889_50_Analysed_BDT_output-9c4b127e-5141-4f5b-82cf-3153df9693ae.root");
    //std::string base_path = "/nashome/c/coackley/v10_09_00_through/";
    
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/analysed_v10_09_00/BDT/merged.root");
    //std::string base_path = "/nashome/c/coackley/nuEPlotsWithCosmicsSCEON_v10_09_00/";
    
    //TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/Nu+E_Cosmics/analysed_v10_06_00/BDT/merged.root");
    //std::string base_path = "/nashome/c/coackley/nuEPlotsWithCosmicsSCEON_v10_06_00/";
    
    std::string base_path_perSlice = base_path + "perSlicePlots/"; 
    bool makePerSlicePlots = false;   // set to false to skip making per-slice plots


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
    std::vector<double> *truth_particleNeutrinoParent = nullptr;

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
    std::vector<double> *reco_particleNumHits = nullptr;
    std::vector<double> *reco_particleNumMatchedHits = nullptr;
    std::vector<double> *reco_particleNumTrueHits = nullptr;
    std::vector<double> *reco_particleTruePDG = nullptr;

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
    tree->SetBranchAddress("truth_particleNeutrinoParent", &truth_particleNeutrinoParent);

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
    tree->SetBranchAddress("reco_particleNumHits", &reco_particleNumHits);
    tree->SetBranchAddress("reco_particleNumMatchedHits", &reco_particleNumMatchedHits);
    tree->SetBranchAddress("reco_particleNumTrueHits", &reco_particleNumTrueHits);
    tree->SetBranchAddress("reco_particleTruePDG", &reco_particleTruePDG);

    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness", &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_slicePurity", &reco_slicePurity);
    tree->SetBranchAddress("reco_sliceScore", &reco_sliceScore);

    Long64_t numEntries = tree->GetEntries();

    // Truth Plots
    auto trueETheta2 = createHistGroup("trueETheta2", "E_{true}#theta_{true}^{2}", "E_{true}#theta_{true}^{2} (MeV rad^{2})", 32, 0, 4.088);

    // Number Plots
    auto numSlices = createHistGroup("numSlices", "Number of Slices in an Event", "Number of Slices", 45, 0, 45);
    // Number of slices with a CRUMBS score
    auto numSlicesCRUMBS = createHistGroup("numSlicesCRUMBS", "Number of Slices with CRUMBS Score in an Event", "Number of Slices with a CRUMBS Score", 10, 0, 10);
    // Number of slices with a completeness > 0
    auto numSlicesCompleteness = createHistGroup("numSlicesCompleteness", "Number of Slices with Completeness > 0 in an Event", "Number of Slices with Completeness > 0", 10, 0, 10);
    auto numRecoNeutrinos = createHistGroup("numRecoNeutrinos", "Number of Reco Neutrinos in an Event", "Number of Reco Neutrinos", 10, 0, 10);
    

    // Highest CRUMBS Score Slice Plots
    auto sliceCompletenessCRUMBS = createHistGroup("sliceCompletenessCRUMBS", "Completeness of the Slice with the Highest CRUMBS Score", "Completeness", 102, 0, 1.02); // Bin width = 0.005
    auto sliceScoreCRUMBS = createHistGroup("sliceScoreCRUMBS", "CRUMBS Score of the Slice with the Highest CRUMBS Score", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCRUMBS = createHistGroup("slicePurityCRUMBS", "Purity of the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto deltaXCRUMBS = createHistGroup("deltaXCRUMBS", "#Deltax Distribution: Slice with Highest CRUMBS Score", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaYCRUMBS = createHistGroup("deltaYCRUMBS", "#Deltay Distribution: Slice with Highest CRUMBS Score", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZCRUMBS = createHistGroup("deltaZCRUMBS", "#Deltaz Distribution: Slice with Highest CRUMBS Score", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaRCRUMBS = createHistGroup("deltaRCRUMBS", "#Delta#bar{r} Distribution: Slice with Highest CRUMBS Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    auto numPFPsCRUMBS = createHistGroup("numPFPsCRUMBS", "Number of PFPs in the Slice with the Highest CRUMBS Score", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCRUMBS = createHistGroup("ratioChosenSummedEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCRUMBS = createHistGroup("ratioChosenTrueEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCRUMBS = createHistGroup("ratioSummedTrueEnergyCRUMBS", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    auto angleDifferenceCRUMBS = createHistGroup("angleDifferenceCRUMBS", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto EtrueThetaRecoCRUMBS = createHistGroup("EtrueThetaRecoCRUMBS", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCRUMBS = createHistGroup("ERecoSumThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBS = createHistGroup("ERecoHighestThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCRUMBS = createHistGroup("ERecoSumThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCRUMBS = createHistGroup("ERecoHighestThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);

    auto trueZCRUMBS = createHistGroup("trueZCRUMBS", "Truth z Distribution: Slice with Highest CRUMBS Score", "z_{True} (cm)", 40, 0, 500);
    auto recoZCRUMBS = createHistGroup("recoZCRUMBS", "Reco z Distribution: Slice with Highest CRUMBS Score", "z_{Reco} (cm)", 40, 0, 500);
    
    // Highest Completeness Score Slice Plots
    auto sliceCompletenessCompleteness = createHistGroup("sliceCompletenessCompleteness", "Completeness of the Slice with the Highest Completeness", "Completeness Score", 102, 0, 1.02); // Bin width = 0.01
    auto sliceScoreCompleteness = createHistGroup("sliceScoreCompleteness", "CRUMBS Score of the Slice with the Highest Completeness", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCompleteness = createHistGroup("slicePurityCompleteness", "Purity of the Slice with the Highest Completeness", "Purity", 50, 0, 1.02);
    auto deltaXCompleteness = createHistGroup("deltaXCompleteness", "#Deltax Distribution: Slice with Highest Completeness", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    auto deltaYCompleteness = createHistGroup("deltaYCompleteness", "#Deltay Distribution: Slice with Highest Completeness", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    auto deltaZCompleteness = createHistGroup("deltaZCompleteness", "#Deltaz Distribution: Slice with Highest Completeness", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    auto deltaRCompleteness = createHistGroup("deltaRCompleteness", "#Delta#bar{r} Distribution: Slice with Highest Completeness Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    auto numPFPsCompleteness = createHistGroup("numPFPsCompleteness", "Number of PFPs in the Slice with the Highest Completeness", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCompleteness = createHistGroup("ratioChosenSummedEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest Completeness", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCompleteness = createHistGroup("ratioChosenTrueEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCompleteness = createHistGroup("ratioSummedTrueEnergyCompleteness", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    auto angleDifferenceCompleteness = createHistGroup("angleDifferenceCompleteness", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest Completeness", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto EtrueThetaRecoCompleteness = createHistGroup("EtrueThetaRecoCompleteness", "E_{true}#theta_{reco}^{2}: Slice with Highest Completeness", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCompleteness = createHistGroup("ERecoSumThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCompleteness = createHistGroup("ERecoHighestThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCompleteness = createHistGroup("ERecoSumThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCompleteness = createHistGroup("ERecoHighestThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
  
    const double xMin = -201.3, xMax = 201.3;
    const double yMin = -203.8, yMax = 203.8;
    const double zMin = 0.0,    zMax = 509.4;
    const double step = 10.0;  // slice width in cm
    
    TH2D *xCoordAngleDifferenceBDTCRUMBS = new TH2D("xCoordAngleDifferenceBDTCRUMBS", "", 50, xMin, xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDTCRUMBS = new TH2D("yCoordAngleDifferenceBDTCRUMBS", "", 50, yMin, yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDTCRUMBS = new TH2D("zCoordAngleDifferenceBDTCRUMBS", "", 50, zMin, zMax, 60, 0, 180);
    TH2D *xCoordAngleDifferenceDLUbooneCRUMBS = new TH2D("xCoordAngleDifferenceDLUbooneCRUMBS", "", 50, xMin, xMax, 60, 0, 180);
    TH2D *yCoordAngleDifferenceDLUbooneCRUMBS = new TH2D("yCoordAngleDifferenceDLUbooneCRUMBS", "", 50, yMin, yMax, 60, 0, 180);
    TH2D *zCoordAngleDifferenceDLUbooneCRUMBS = new TH2D("zCoordAngleDifferencDLUbooneCRUMBS", "", 50, zMin, zMax, 60, 0, 180);

    TH2D *xCoordAngleDifferenceBDTCRUMBS_low = new TH2D("xCoordAngleDifferenceBDTCRUMBS_low", "", 10, xMin, (xMin + 20), 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDTCRUMBS_low = new TH2D("yCoordAngleDifferenceBDTCRUMBS_low", "", 10, yMin, (yMin + 20), 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDTCRUMBS_low = new TH2D("zCoordAngleDifferenceBDTCRUMBS_low", "", 10, zMin, (zMin + 20), 40, 0, 180);
    TH2D *xCoordAngleDifferenceBDTCRUMBS_high = new TH2D("xCoordAngleDifferenceBDTCRUMBS_high", "", 10, (xMax - 20), xMax, 40, 0, 180);
    TH2D *yCoordAngleDifferenceBDTCRUMBS_high = new TH2D("yCoordAngleDifferenceBDTCRUMBS_high", "", 10, (yMax - 20), yMax, 40, 0, 180);
    TH2D *zCoordAngleDifferenceBDTCRUMBS_high = new TH2D("zCoordAngleDifferenceBDTCRUMBS_high", "", 20, (zMax - 40), zMax, 60, 0, 180);

    TProfile *xCoordAngleDifferenceBDTCRUMBSProfile = new TProfile("xCoordAngleDifferenceBDTCRUMBSProfile", "", 50, xMin, xMax, 0, 180);
    TProfile *yCoordAngleDifferenceBDTCRUMBSProfile = new TProfile("yCoordAngleDifferenceBDTCRUMBSProfile", "", 50, yMin, yMax, 0, 180);
    TProfile *zCoordAngleDifferenceBDTCRUMBSProfile = new TProfile("zCoordAngleDifferenceBDTCRUMBSProfile", "", 50, zMin, zMax, 0, 180);
    TProfile *xCoordAngleDifferenceDLUbooneCRUMBSProfile = new TProfile("xCoordAngleDifferenceDLUbooneCRUMBSProfile", "", 50, xMin, xMax, 0, 180);
    TProfile *yCoordAngleDifferenceDLUbooneCRUMBSProfile = new TProfile("yCoordAngleDifferenceDLUbooneCRUMBSProfile", "", 50, yMin, yMax, 0, 180);
    TProfile *zCoordAngleDifferenceDLUbooneCRUMBSProfile = new TProfile("zCoordAngleDifferencDLUbooneCRUMBSProfile", "", 50, zMin, zMax, 0, 180);

    TProfile *xCoordAngleDifferenceBDTCRUMBSProfile_low = new TProfile("xCoordAngleDifferenceBDTCRUMBSProfile_low", "", 10, xMin, (xMin + 20), 0, 180);
    TProfile *yCoordAngleDifferenceBDTCRUMBSProfile_low = new TProfile("yCoordAngleDifferenceBDTCRUMBSProfile_low", "", 10, yMin, (yMin + 20), 0, 180);
    TProfile *zCoordAngleDifferenceBDTCRUMBSProfile_low = new TProfile("zCoordAngleDifferenceBDTCRUMBSProfile_low", "", 10, zMin, (zMin + 20), 0, 180);
    TProfile *xCoordAngleDifferenceBDTCRUMBSProfile_high = new TProfile("xCoordAngleDifferenceBDTCRUMBSProfile_high", "", 10, (xMax - 20), xMax, 0, 180);
    TProfile *yCoordAngleDifferenceBDTCRUMBSProfile_high = new TProfile("yCoordAngleDifferenceBDTCRUMBSProfile_high", "", 10, (yMax - 20), yMax, 0, 180);
    TProfile *zCoordAngleDifferenceBDTCRUMBSProfile_high = new TProfile("zCoordAngleDifferenceBDTCRUMBSProfile_high", "", 20, (zMax - 40), zMax, 0, 180);

    std::vector<double> xVals_xCoord;
    std::vector<double> xVals_yCoord;
    std::vector<double> xVals_zCoord;
    std::vector<double> yVals;

    std::vector<edgeSlice> xEdgeSlices, yEdgeSlices, zEdgeSlices;

    for(double low = xMin; low < xMax; low += step){
        double high = std::min(low + step, xMax);
        std::string name = "x_slice_" + std::to_string(xEdgeSlices.size());
        std::string title = "Angle Difference CRUMBS (X in [" + std::to_string(low) + "," + std::to_string(high) + "])";
        xEdgeSlices.push_back({ createHistGroup(name, title, "Angle Difference (deg)", 45, 0, 180), low, high });
    }

    for(double low = yMin; low < yMax; low += step){
        double high = std::min(low + step, yMax);
        std::string name = "y_slice_" + std::to_string(yEdgeSlices.size());
        std::string title = "Angle Difference CRUMBS (Y in [" + std::to_string(low) + "," + std::to_string(high) + "])";
        yEdgeSlices.push_back({ createHistGroup(name, title, "Angle Difference (deg)", 45, 0, 180), low, high });
    }

    for(double low = zMin; low < zMax; low += step){
        double high = std::min(low + step, zMax);
        std::string name = "z_slice_" + std::to_string(zEdgeSlices.size());
        std::string title = "Angle Difference CRUMBS (Z in [" + std::to_string(low) + "," + std::to_string(high) + "])";
        zEdgeSlices.push_back({ createHistGroup(name, title, "Angle Difference (deg)", 45, 0, 180), low, high });
    }
        
    // TO HERE!

    counter numEventsTotal;
    counter numEventsTrueNeutrino;
    counter numEventsTrueElectron;
    counter numEventsSlices;
    counter numEventsRecoNeutrino;
    counter numEventsSliceNotEqualNeutrino;
    counter numSlicesCRUMBSCompletenessZero;
    counter numEventsCRUMBSRecoParticle;
    counter sameSliceSelected;

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);

        printf("___________________________________________________________________\n");
        std::vector<recoParticle> recoParticlesInEvent;
        std::vector<recoNeutrino> recoNeutrinosInEvent;
        std::vector<recoSlice> recoSlicesInEvent;

        recoParticle chosenRecoParticleCRUMBS;
        recoParticle chosenRecoParticleCompleteness;
        recoNeutrino chosenRecoNeutrinoCRUMBS;
        recoNeutrino chosenRecoNeutrinoCompleteness;
        recoSlice chosenRecoSliceCompleteness;
        recoSlice chosenRecoSliceCRUMBS;
        trueParticle chosenTrueParticle;
        trueNeutrino chosenTrueNeutrino;
        double energyInSlice = 0;
        double numPFPsInEvent = 0;
        double numPFPsInChosenSlice = 0;

        if(DLCurrent == 0){
            numEventsTotal.uboone++;
        } else if(DLCurrent == 1){
            numEventsTotal.dune++;
        } else if(DLCurrent == 2){
            numEventsTotal.current++;
        } else if(DLCurrent == 3){
            numEventsTotal.cheated++;
        } else if(DLCurrent == 4){
            numEventsTotal.sbnd++;
        } else if(DLCurrent == 5){
            numEventsTotal.nue++;
        } else if(DLCurrent == 6){
            numEventsTotal.nue_100k++; 
        }

        int truthNeutrinoInTPC = 0;
        for(size_t j = 0; j < truth_neutrinoVX->size(); ++j){
            if(truth_neutrinoTPCValid->at(j) == 1){
                truthNeutrinoInTPC = 1;
                chosenTrueNeutrino.vx = truth_neutrinoVX->at(j);
                chosenTrueNeutrino.vy = truth_neutrinoVY->at(j);
                chosenTrueNeutrino.vz = truth_neutrinoVZ->at(j); 
                chosenTrueNeutrino.CCNC = truth_CCNC->at(j);
                chosenTrueNeutrino.pdg = truth_neutrinoType->at(j);
                chosenTrueNeutrino.leptonpdg = truth_chargedLepton->at(j);
                chosenTrueNeutrino.tpcID = truth_neutrinoTPCID->at(j);
                chosenTrueNeutrino.tpcValid = truth_neutrinoTPCValid->at(j);
            } else if(truth_neutrinoTPCValid->at(j) != -999999 && truth_neutrinoVX->at(j) < 201.3 && truth_neutrinoVX->at(j) > -201.3 && truth_neutrinoVY->at(j) < 203.732 && truth_neutrinoVY->at(j) > -203.732 && truth_neutrinoVZ->at(j) > 0 && truth_neutrinoVZ->at(j) < 509.4){
                //printf("True Neutrino Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            } else if(truth_neutrinoTPCValid->at(j) != -999999){
                //printf("True Neutrino Not Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            }
        }

        if(truth_neutrinoVX->size() == 1 && truth_neutrinoVX->at(0) == -999999) printf("No True Neutrino in Event\n");

        // If there are no truth neutrinos within the TPC, skip the event
        if(truthNeutrinoInTPC == 0) continue;

        if(DLCurrent == 0){
            numEventsTrueNeutrino.uboone++;
        } else if(DLCurrent == 1){
            numEventsTrueNeutrino.dune++;
        } else if(DLCurrent == 2){
            numEventsTrueNeutrino.current++;
        } else if(DLCurrent == 3){
            numEventsTrueNeutrino.cheated++;
        } else if(DLCurrent == 4){
            numEventsTrueNeutrino.sbnd++;
        } else if(DLCurrent == 5){
            numEventsTrueNeutrino.nue++;
        } else if(DLCurrent == 6){
            numEventsTrueNeutrino.nue_100k++;
        }

        printf("True Neutrino Vertex = (%f, %f, %f)\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz);
        // Pick out the true recoil electron
        int truthElectron = 0;
        for(size_t j = 0; j < truth_particlePDG->size(); ++j){
            if(truth_particlePDG->at(j) == 11 && truth_particleNeutrinoParent->at(j) == 1){
                // This is an electron (PDG == 11) whose parent has TrackID of 0 (has neutrino parent == 1)
                printf("True Primary Electron: Vertex = (%f, %f, %f)\n", truth_particleVX->at(j), truth_particleVY->at(j), truth_particleVZ->at(j));
                if(std::abs(truth_particleVX->at(j) - chosenTrueNeutrino.vx) < 0.5 && std::abs(truth_particleVY->at(j) - chosenTrueNeutrino.vy) < 0.5 && std::abs(truth_particleVZ->at(j) - chosenTrueNeutrino.vz) < 0.5){
                    // This is the true recoil electron
                    truthElectron = 1;
                    chosenTrueParticle.pdg = truth_particlePDG->at(j);
                    chosenTrueParticle.vx = truth_particleVX->at(j);
                    chosenTrueParticle.vy = truth_particleVY->at(j);
                    chosenTrueParticle.vz = truth_particleVZ->at(j);
                    chosenTrueParticle.px = truth_particlePX->at(j);
                    chosenTrueParticle.py = truth_particlePY->at(j);
                    chosenTrueParticle.pz = truth_particlePZ->at(j);
                    chosenTrueParticle.energy = truth_particleEnergy->at(j);
                    chosenTrueParticle.angle = truth_particleAngle->at(j);
                    chosenTrueParticle.ETheta2 = truth_particleETheta2->at(j);
                    chosenTrueParticle.dx = truth_particleDirectionX->at(j);
                    chosenTrueParticle.dy = truth_particleDirectionY->at(j);
                    chosenTrueParticle.dz = truth_particleDirectionZ->at(j);
                    chosenTrueParticle.neutrinoParent = truth_particleNeutrinoParent->at(j);
                } else{
                    printf("Doesn't Pass!!\n");
                }      
            }
        }

        // If there are no true electrons from the primary neutrino, skip the event
        if(truthElectron == 0){
            printf("LOOK AT THIS ONE!!!!!!!!!\n");
            continue;
        }

        if(DLCurrent == 0){
            numEventsTrueElectron.uboone++;
            trueETheta2.uboone->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 1){
            numEventsTrueElectron.dune++;
            trueETheta2.dune->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 2){
            numEventsTrueElectron.current++;
            trueETheta2.current->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 3){
            numEventsTrueElectron.cheated++;
            trueETheta2.cheated->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 4){
            numEventsTrueElectron.sbnd++;
            trueETheta2.sbnd->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 5){
            numEventsTrueElectron.nue++;
            trueETheta2.nue->Fill(chosenTrueParticle.ETheta2);
        } else if(DLCurrent == 6){
            numEventsTrueElectron.nue_100k++;
            trueETheta2.nue_100k->Fill(chosenTrueParticle.ETheta2);
        }

        int slice = 0;
        int crumbsSlice = 0;
        int completenessSlice = 0;
        double numSlicesInEvent = 0;
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            if(reco_sliceID->at(j) != -999999){
                slice = 1;
                numSlicesInEvent++;

                recoSlice recoSlice;
                recoSlice.id = reco_sliceID->at(j);
                recoSlice.completeness = reco_sliceCompleteness->at(j);
                recoSlice.purity = reco_slicePurity->at(j);
                recoSlice.score = reco_sliceScore->at(j);
                recoSlicesInEvent.push_back(recoSlice);
        
                if(recoSlice.score != -999999) crumbsSlice++;  
                if(recoSlice.completeness > 0) completenessSlice++;  
                printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", j, recoSlice.id, recoSlice.completeness, recoSlice.purity, recoSlice.score);
            }
        }    

        // If there are no slices, skip the event
        if(slice == 0) continue;  
     
        // Pick the slice with the highest completeness
        chosenRecoSliceCompleteness = chooseSlice(recoSlicesInEvent, 0);
        printf("\nSlice with Highest Completeness: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", chosenRecoSliceCompleteness.id, chosenRecoSliceCompleteness.completeness, chosenRecoSliceCompleteness.purity, chosenRecoSliceCompleteness.score);

        // Pick the slice with the highest score
        chosenRecoSliceCRUMBS = chooseSlice(recoSlicesInEvent, 1);
        printf("Slice with Highest CRUMBS Score: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n\n", chosenRecoSliceCRUMBS.id, chosenRecoSliceCRUMBS.completeness, chosenRecoSliceCRUMBS.purity, chosenRecoSliceCRUMBS.score);
       
 
        if(DLCurrent == 0){
            numEventsSlices.uboone++;
            numSlices.uboone->Fill(numSlicesInEvent);
            numSlicesCRUMBS.uboone->Fill(crumbsSlice);
            numSlicesCompleteness.uboone->Fill(completenessSlice);
            
            sliceScoreCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.uboone->Fill(chosenRecoSliceCompleteness.score);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.uboone++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.uboone++;
        } else if(DLCurrent == 1){
            numEventsSlices.dune++;
            numSlices.dune->Fill(numSlicesInEvent);
            numSlicesCRUMBS.dune->Fill(crumbsSlice);
            numSlicesCompleteness.dune->Fill(completenessSlice);

            sliceScoreCRUMBS.dune->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.dune->Fill(chosenRecoSliceCompleteness.score);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.dune++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.dune++;
        } else if(DLCurrent == 2){
            numEventsSlices.current++;
            numSlices.current->Fill(numSlicesInEvent);
            numSlicesCRUMBS.current->Fill(crumbsSlice);
            numSlicesCompleteness.current->Fill(completenessSlice);
            
            sliceScoreCRUMBS.current->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.current->Fill(chosenRecoSliceCompleteness.score);
                
            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.current++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.current++;
        } else if(DLCurrent == 3){
            numEventsSlices.cheated++;
            numSlices.cheated->Fill(numSlicesInEvent);
            numSlicesCRUMBS.cheated->Fill(crumbsSlice);
            numSlicesCompleteness.cheated->Fill(completenessSlice);
            
            sliceScoreCRUMBS.cheated->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.cheated->Fill(chosenRecoSliceCompleteness.score);
            
            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.cheated++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.cheated++;
        } else if(DLCurrent == 4){
            numEventsSlices.sbnd++;
            numSlices.sbnd->Fill(numSlicesInEvent);
            numSlicesCRUMBS.sbnd->Fill(crumbsSlice);
            numSlicesCompleteness.sbnd->Fill(completenessSlice);

            sliceScoreCRUMBS.sbnd->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.sbnd->Fill(chosenRecoSliceCompleteness.score);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.sbnd++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.sbnd++;
        } else if(DLCurrent == 5){
            numEventsSlices.nue++;
            numSlices.nue->Fill(numSlicesInEvent);
            numSlicesCRUMBS.nue->Fill(crumbsSlice);
            numSlicesCompleteness.nue->Fill(completenessSlice);

            sliceScoreCRUMBS.nue->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.nue->Fill(chosenRecoSliceCompleteness.score);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.nue++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.nue++;
        } else if(DLCurrent == 6){
            numEventsSlices.nue_100k++;
            numSlices.nue_100k->Fill(numSlicesInEvent);
            numSlicesCRUMBS.nue_100k->Fill(crumbsSlice);
            numSlicesCompleteness.nue_100k->Fill(completenessSlice);

            sliceScoreCRUMBS.nue_100k->Fill(chosenRecoSliceCRUMBS.score);
            sliceScoreCompleteness.nue_100k->Fill(chosenRecoSliceCompleteness.score);

            if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.nue_100k++;
            if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.nue_100k++;
        }

        int neutrino = 0;
        int numRecoNeutrinosInEvent = 0;
        std::cout << "num reco neutrinos: " << reco_neutrinoPDG->size() << std::endl;
        for(size_t j = 0; j < reco_neutrinoPDG->size(); ++j){
            if(reco_neutrinoPDG->at(j) != -999999){
                neutrino = 1;
                numRecoNeutrinosInEvent++;
                std::cout << "1" << std::endl;
                
                recoNeutrino recoNeutrino;
                recoNeutrino.pdg = reco_neutrinoPDG->at(j);
                recoNeutrino.isPrimary = reco_neutrinoIsPrimary->at(j);
                recoNeutrino.vx = reco_neutrinoVX->at(j);
                recoNeutrino.vy = reco_neutrinoVY->at(j);
                recoNeutrino.vz = reco_neutrinoVZ->at(j);
                recoNeutrino.sliceID = reco_neutrinoSliceID->at(j);
                recoNeutrinosInEvent.push_back(recoNeutrino);

                printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, recoNeutrino.pdg, recoNeutrino.isPrimary, recoNeutrino.vx, recoNeutrino.vy, recoNeutrino.vz, recoNeutrino.sliceID);
            }
        }
        
        if(DLCurrent == 0) numRecoNeutrinos.uboone->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 1) numRecoNeutrinos.dune->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 2) numRecoNeutrinos.current->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 3) numRecoNeutrinos.cheated->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 4) numRecoNeutrinos.sbnd->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 5) numRecoNeutrinos.nue->Fill(numRecoNeutrinosInEvent);
        if(DLCurrent == 6) numRecoNeutrinos.nue_100k->Fill(numRecoNeutrinosInEvent);

        // Skip the event if there are no reconstructed neutrinos
        if(neutrino == 0) continue;
    
        // The number of slices with a CRUMBS score != the number of reconstructed neutrinos in the event 
        if(crumbsSlice != numRecoNeutrinosInEvent){
            if(DLCurrent == 0) numEventsSliceNotEqualNeutrino.uboone++; 
            if(DLCurrent == 1) numEventsSliceNotEqualNeutrino.dune++;     
            if(DLCurrent == 2) numEventsSliceNotEqualNeutrino.current++;
            if(DLCurrent == 3) numEventsSliceNotEqualNeutrino.cheated++;
            if(DLCurrent == 4) numEventsSliceNotEqualNeutrino.sbnd++;
            if(DLCurrent == 5) numEventsSliceNotEqualNeutrino.nue++;
            if(DLCurrent == 6) numEventsSliceNotEqualNeutrino.nue_100k++;
        } 

        // Look at reconstruction when we pick the slice with highest CRUMBS score
        chosenRecoNeutrinoCRUMBS = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCRUMBS.id); 

        double deltaXCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vx - chosenTrueNeutrino.vx);
        double deltaYCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vy - chosenTrueNeutrino.vy);
        double deltaZCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vz - chosenTrueNeutrino.vz);
        double deltaRCRUMBSValue = std::sqrt((deltaXCRUMBSValue * deltaXCRUMBSValue) + (deltaYCRUMBSValue * deltaYCRUMBSValue) + (deltaZCRUMBSValue * deltaZCRUMBSValue));

        chosenRecoNeutrinoCompleteness = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCompleteness.id);
        double deltaXCompletenessValue = (chosenRecoNeutrinoCompleteness.vx - chosenTrueNeutrino.vx);
        double deltaYCompletenessValue = (chosenRecoNeutrinoCompleteness.vy - chosenTrueNeutrino.vy);
        double deltaZCompletenessValue = (chosenRecoNeutrinoCompleteness.vz - chosenTrueNeutrino.vz);
        double deltaRCompletenessValue = std::sqrt((deltaXCompletenessValue * deltaXCompletenessValue) + (deltaYCompletenessValue * deltaYCompletenessValue) + (deltaZCompletenessValue * deltaZCompletenessValue));

        if(DLCurrent == 0){
            numEventsRecoNeutrino.uboone++;
            deltaXCRUMBS.uboone->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.uboone->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.uboone->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.uboone->Fill(deltaRCRUMBSValue);

            trueZCRUMBS.uboone->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.uboone->Fill(chosenRecoNeutrinoCRUMBS.vz);

            deltaXCompleteness.uboone->Fill(deltaXCompletenessValue);
            deltaYCompleteness.uboone->Fill(deltaYCompletenessValue);
            deltaZCompleteness.uboone->Fill(deltaZCompletenessValue);
            deltaRCompleteness.uboone->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 1){
            numEventsRecoNeutrino.dune++;
            deltaXCRUMBS.dune->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.dune->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.dune->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.dune->Fill(deltaRCRUMBSValue);
        
            trueZCRUMBS.dune->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.dune->Fill(chosenRecoNeutrinoCRUMBS.vz);
            
            deltaXCompleteness.dune->Fill(deltaXCompletenessValue);
            deltaYCompleteness.dune->Fill(deltaYCompletenessValue);
            deltaZCompleteness.dune->Fill(deltaZCompletenessValue);
            deltaRCompleteness.dune->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 2){
            numEventsRecoNeutrino.current++;
            deltaXCRUMBS.current->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.current->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.current->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.current->Fill(deltaRCRUMBSValue);
            
            trueZCRUMBS.current->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.current->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
            deltaXCompleteness.current->Fill(deltaXCompletenessValue);
            deltaYCompleteness.current->Fill(deltaYCompletenessValue);
            deltaZCompleteness.current->Fill(deltaZCompletenessValue);
            deltaRCompleteness.current->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 3){
            numEventsRecoNeutrino.cheated++;
            deltaXCRUMBS.cheated->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.cheated->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.cheated->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.cheated->Fill(deltaRCRUMBSValue);
            
            trueZCRUMBS.cheated->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.cheated->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
            deltaXCompleteness.cheated->Fill(deltaXCompletenessValue);
            deltaYCompleteness.cheated->Fill(deltaYCompletenessValue);
            deltaZCompleteness.cheated->Fill(deltaZCompletenessValue);
            deltaRCompleteness.cheated->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 4){
            numEventsRecoNeutrino.sbnd++;
            deltaXCRUMBS.sbnd->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.sbnd->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.sbnd->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.sbnd->Fill(deltaRCRUMBSValue);
            
            trueZCRUMBS.sbnd->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.sbnd->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
            deltaXCompleteness.sbnd->Fill(deltaXCompletenessValue);
            deltaYCompleteness.sbnd->Fill(deltaYCompletenessValue);
            deltaZCompleteness.sbnd->Fill(deltaZCompletenessValue);
            deltaRCompleteness.sbnd->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 5){
            numEventsRecoNeutrino.nue++;
            deltaXCRUMBS.nue->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.nue->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.nue->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.nue->Fill(deltaRCRUMBSValue);
            
            trueZCRUMBS.nue->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.nue->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
            deltaXCompleteness.nue->Fill(deltaXCompletenessValue);
            deltaYCompleteness.nue->Fill(deltaYCompletenessValue);
            deltaZCompleteness.nue->Fill(deltaZCompletenessValue);
            deltaRCompleteness.nue->Fill(deltaRCompletenessValue);
        } else if(DLCurrent == 6){
            numEventsRecoNeutrino.nue_100k++;
            deltaXCRUMBS.nue_100k->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.nue_100k->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.nue_100k->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.nue_100k->Fill(deltaRCRUMBSValue);
            
            trueZCRUMBS.nue_100k->Fill(chosenTrueNeutrino.vz);
            recoZCRUMBS.nue_100k->Fill(chosenRecoNeutrinoCRUMBS.vz);
        
            deltaXCompleteness.nue_100k->Fill(deltaXCompletenessValue);
            deltaYCompleteness.nue_100k->Fill(deltaYCompletenessValue);
            deltaZCompleteness.nue_100k->Fill(deltaZCompletenessValue);
            deltaRCompleteness.nue_100k->Fill(deltaRCompletenessValue);
        }

        int recoparticle = 0;
        int recoparticleCRUMBS = 0;
        int recoparticleCompleteness = 0;

        for(size_t j = 0; j < reco_particlePDG->size(); ++j){
            if(reco_particlePDG->at(j) != -999999){
                recoparticle = 1;

                recoParticle recoParticle;
                recoParticle.pdg = reco_particlePDG->at(j);
                recoParticle.isPrimary = reco_particleIsPrimary->at(j);
                recoParticle.vx = reco_particleVX->at(j);
                recoParticle.vy = reco_particleVY->at(j);
                recoParticle.vz = reco_particleVZ->at(j);
                recoParticle.dx = reco_particleDirectionX->at(j);
                recoParticle.dy = reco_particleDirectionY->at(j);
                recoParticle.dz = reco_particleDirectionZ->at(j);
                recoParticle.sliceID = reco_particleSliceID->at(j);
                recoParticle.bestPlaneEnergy = reco_particleBestPlaneEnergy->at(j);
                recoParticle.theta = reco_particleTheta->at(j);
                recoParticle.trackscore = reco_particleTrackScore->at(j);
                recoParticle.completeness = reco_particleCompleteness->at(j);
                recoParticle.purity = reco_particlePurity->at(j);
                recoParticle.numHits = reco_particleNumHits->at(j);
                recoParticle.numMatchedHits = reco_particleNumMatchedHits->at(j);
                recoParticle.numTrueHits = reco_particleNumTrueHits->at(j);
                recoParticle.truePDG = reco_particleTruePDG->at(j);
                recoParticlesInEvent.push_back(recoParticle);
            
                if(recoParticle.sliceID == chosenRecoSliceCRUMBS.id) recoparticleCRUMBS = 1;
                //if(recoParticle.sliceID == chosenRecoSliceCompleteness.id) recoparticleCompleteness = 1;
            }
        }

        // There are no reco particles in the event
        if(recoparticle == 0) continue;
        
        // No reco particles in the slice with the highest CRUMBS score
        if(recoparticleCRUMBS == 0) continue;
        
        // No reco particles in the slice with the highest completeness
        //if(recoparticleCompleteness == 0) continue;

        double totalSliceEnergyCRUMBS = 0;
        double numPFPsSliceCRUMBS = 0;
        chosenRecoParticleCRUMBS = choosePFP(recoParticlesInEvent, chosenRecoSliceCRUMBS.id, totalSliceEnergyCRUMBS, numPFPsSliceCRUMBS);
        double chosenSlicePurityCRUMBS = slicePurityCalculator(recoParticlesInEvent, chosenRecoSliceCRUMBS.id);
        double chosenSliceCompletenessCRUMBS = sliceCompletenessCalculator(recoParticlesInEvent, chosenRecoSliceCRUMBS.id);

        double aDOTbCRUMBS = ((chosenRecoParticleCRUMBS.dx * chosenTrueParticle.dx) + (chosenRecoParticleCRUMBS.dy * chosenTrueParticle.dy) + (chosenRecoParticleCRUMBS.dz * chosenTrueParticle.dz));
        double aMagCRUMBS = std::sqrt((chosenRecoParticleCRUMBS.dx * chosenRecoParticleCRUMBS.dx) + (chosenRecoParticleCRUMBS.dy * chosenRecoParticleCRUMBS.dy) + (chosenRecoParticleCRUMBS.dz * chosenRecoParticleCRUMBS.dz));
        double bMagCRUMBS = std::sqrt((chosenTrueParticle.dx * chosenTrueParticle.dx) + (chosenTrueParticle.dy * chosenTrueParticle.dy) + (chosenTrueParticle.dz * chosenTrueParticle.dz));
        double cosAngleCRUMBS = aDOTbCRUMBS/(aMagCRUMBS * bMagCRUMBS);
        if(cosAngleCRUMBS > 1.0) cosAngleCRUMBS = 1.0;
        if(cosAngleCRUMBS < -1.0) cosAngleCRUMBS = -1.0;
        double angleDiffCRUMBS = TMath::ACos(cosAngleCRUMBS) * TMath::RadToDeg();

        double totalSliceEnergyCompleteness = 0;
        double numPFPsSliceCompleteness = 0;
        chosenRecoParticleCompleteness = choosePFP(recoParticlesInEvent, chosenRecoSliceCompleteness.id, totalSliceEnergyCompleteness, numPFPsSliceCompleteness);
        double chosenSlicePurityCompleteness = slicePurityCalculator(recoParticlesInEvent, chosenRecoSliceCompleteness.id);
        double chosenSliceCompletenessCompleteness = sliceCompletenessCalculator(recoParticlesInEvent, chosenRecoSliceCompleteness.id);
        
        double aDOTbCompleteness = ((chosenRecoParticleCompleteness.dx * chosenTrueParticle.dx) + (chosenRecoParticleCompleteness.dy * chosenTrueParticle.dy) + (chosenRecoParticleCompleteness.dz * chosenTrueParticle.dz));
        double aMagCompleteness = std::sqrt((chosenRecoParticleCompleteness.dx * chosenRecoParticleCompleteness.dx) + (chosenRecoParticleCompleteness.dy * chosenRecoParticleCompleteness.dy) + (chosenRecoParticleCompleteness.dz * chosenRecoParticleCompleteness.dz));
        double bMagCompleteness = std::sqrt((chosenTrueParticle.dx * chosenTrueParticle.dx) + (chosenTrueParticle.dy * chosenTrueParticle.dy) + (chosenTrueParticle.dz * chosenTrueParticle.dz));
        double cosAngleCompleteness = aDOTbCompleteness/(aMagCompleteness * bMagCompleteness);
        if(cosAngleCompleteness > 1.0) cosAngleCompleteness = 1.0;
        if(cosAngleCompleteness < -1.0) cosAngleCompleteness = -1.0;
        double angleDiffCompleteness = TMath::ACos(cosAngleCompleteness) * TMath::RadToDeg();
        
        if(DLCurrent == 0){
            numEventsCRUMBSRecoParticle.uboone++;
            numPFPsCRUMBS.uboone->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.uboone->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.uboone->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            xCoordAngleDifferenceDLUbooneCRUMBS->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceDLUbooneCRUMBS->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceDLUbooneCRUMBS->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 

            slicePurityCRUMBS.uboone->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.uboone->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.uboone->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.uboone->Fill(chosenSliceCompletenessCompleteness);

            numPFPsCompleteness.uboone->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.uboone->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.uboone->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.uboone->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.uboone->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.uboone->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.uboone->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 1){
            numEventsCRUMBSRecoParticle.dune++;
            numPFPsCRUMBS.dune->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.dune->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.dune->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.dune->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.dune->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.dune->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.dune->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);

            slicePurityCRUMBS.dune->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.dune->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.dune->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.dune->Fill(chosenSliceCompletenessCompleteness);

            numPFPsCompleteness.dune->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.dune->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.dune->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.dune->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.dune->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.dune->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.dune->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 2){
            numEventsCRUMBSRecoParticle.current++;
            numPFPsCRUMBS.current->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.current->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.current->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.current->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.current->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.current->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
            xCoordAngleDifferenceBDTCRUMBS->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceBDTCRUMBS->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceBDTCRUMBS->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 
            xCoordAngleDifferenceBDTCRUMBS_low->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceBDTCRUMBS_low->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceBDTCRUMBS_low->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 
            xCoordAngleDifferenceBDTCRUMBS_high->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceBDTCRUMBS_high->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceBDTCRUMBS_high->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 

            xVals_xCoord.push_back(chosenRecoNeutrinoCRUMBS.vx);
            xVals_yCoord.push_back(chosenRecoNeutrinoCRUMBS.vy);
            xVals_zCoord.push_back(chosenRecoNeutrinoCRUMBS.vz);
            yVals.push_back(angleDiffCRUMBS);

            xCoordAngleDifferenceBDTCRUMBSProfile->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceBDTCRUMBSProfile->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceBDTCRUMBSProfile->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 
            xCoordAngleDifferenceBDTCRUMBSProfile_low->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceBDTCRUMBSProfile_low->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceBDTCRUMBSProfile_low->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 
            xCoordAngleDifferenceBDTCRUMBSProfile_high->Fill(chosenRecoNeutrinoCRUMBS.vx, angleDiffCRUMBS); 
            yCoordAngleDifferenceBDTCRUMBSProfile_high->Fill(chosenRecoNeutrinoCRUMBS.vy, angleDiffCRUMBS); 
            zCoordAngleDifferenceBDTCRUMBSProfile_high->Fill(chosenRecoNeutrinoCRUMBS.vz, angleDiffCRUMBS); 

            slicePurityCRUMBS.current->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.current->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.current->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.current->Fill(chosenSliceCompletenessCompleteness);
            
            numPFPsCompleteness.current->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.current->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.current->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.current->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.current->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.current->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.current->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 3){
            numEventsCRUMBSRecoParticle.cheated++;
            numPFPsCRUMBS.cheated->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.cheated->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.cheated->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.cheated->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.cheated->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.cheated->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.cheated->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
            slicePurityCRUMBS.cheated->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.cheated->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.cheated->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.cheated->Fill(chosenSliceCompletenessCompleteness);
        
            numPFPsCompleteness.cheated->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.cheated->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.cheated->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.cheated->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.cheated->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.cheated->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.cheated->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 4){
            numEventsCRUMBSRecoParticle.sbnd++;
            numPFPsCRUMBS.sbnd->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.sbnd->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.sbnd->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.sbnd->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.sbnd->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.sbnd->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.sbnd->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
            slicePurityCRUMBS.sbnd->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.sbnd->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.sbnd->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.sbnd->Fill(chosenSliceCompletenessCompleteness);
        
            numPFPsCompleteness.sbnd->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.sbnd->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.sbnd->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.sbnd->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.sbnd->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.sbnd->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.sbnd->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 5){
            numEventsCRUMBSRecoParticle.nue++;
            numPFPsCRUMBS.nue->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.nue->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.nue->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.nue->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.nue->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.nue->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.nue->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.nue->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.nue->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.nue->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
            slicePurityCRUMBS.nue->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.nue->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.nue->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.nue->Fill(chosenSliceCompletenessCompleteness);
        
            numPFPsCompleteness.nue->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.nue->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.nue->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.nue->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.nue->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.nue->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.nue->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.nue->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.nue->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.nue->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        } else if(DLCurrent == 6){
            numEventsCRUMBSRecoParticle.nue_100k++;
            numPFPsCRUMBS.nue_100k->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.nue_100k->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ratioChosenTrueEnergyCRUMBS.nue_100k->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCRUMBS.nue_100k->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            angleDifferenceCRUMBS.nue_100k->Fill(angleDiffCRUMBS);
            EtrueThetaRecoCRUMBS.nue_100k->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoSumThetaTrueCRUMBS.nue_100k->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCRUMBS.nue_100k->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.nue_100k->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.nue_100k->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
            slicePurityCRUMBS.nue_100k->Fill(chosenSlicePurityCRUMBS);
            sliceCompletenessCRUMBS.nue_100k->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCompleteness.nue_100k->Fill(chosenSlicePurityCompleteness);
            sliceCompletenessCompleteness.nue_100k->Fill(chosenSliceCompletenessCompleteness);
        
            numPFPsCompleteness.nue_100k->Fill(numPFPsSliceCompleteness);
            ratioChosenSummedEnergyCompleteness.nue_100k->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / totalSliceEnergyCompleteness);
            ratioChosenTrueEnergyCompleteness.nue_100k->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy / chosenTrueParticle.energy);
            ratioSummedTrueEnergyCompleteness.nue_100k->Fill(totalSliceEnergyCompleteness / chosenTrueParticle.energy);
            angleDifferenceCompleteness.nue_100k->Fill(angleDiffCompleteness);
            EtrueThetaRecoCompleteness.nue_100k->Fill(chosenTrueParticle.energy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoSumThetaTrueCompleteness.nue_100k->Fill(totalSliceEnergyCompleteness * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoHighestThetaTrueCompleteness.nue_100k->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCompleteness.nue_100k->Fill(totalSliceEnergyCompleteness * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
            ERecoHighestThetaRecoCompleteness.nue_100k->Fill(chosenRecoParticleCompleteness.bestPlaneEnergy * chosenRecoParticleCompleteness.theta * chosenRecoParticleCompleteness.theta);
        }

        // HERE!!

        if(makePerSlicePlots){
            double vx = chosenRecoParticleCRUMBS.vx;
            double vy = chosenRecoParticleCRUMBS.vy;
            double vz = chosenRecoParticleCRUMBS.vz;

            auto fillHistGroup = [&](histGroup& hg){
                if(DLCurrent == 0){      
                    hg.uboone->Fill(angleDiffCRUMBS);
                    std::cout << "FILLED DLUBoone_________________________________________________" << std::endl;
                }
                else if(DLCurrent == 1) hg.dune->Fill(angleDiffCRUMBS);
                else if(DLCurrent == 2){
                    hg.current->Fill(angleDiffCRUMBS);
                    std::cout << "FILLED Current_________________________________________________" << std::endl;
                }
                else if(DLCurrent == 3) hg.cheated->Fill(angleDiffCRUMBS);
                else if(DLCurrent == 4) hg.sbnd->Fill(angleDiffCRUMBS);
                else if(DLCurrent == 5) hg.nue->Fill(angleDiffCRUMBS);
                else if(DLCurrent == 6) hg.nue_100k->Fill(angleDiffCRUMBS);
            };

            // Check X edges
            for(auto& slice : xEdgeSlices){
                if(vx >= slice.low && vx < slice.high){
                    fillHistGroup(slice.hg);
                    break;
                }
            }
            // Check Y edges
            for(auto& slice : yEdgeSlices){
                if(vy >= slice.low && vy < slice.high){
                    fillHistGroup(slice.hg);
                    break;
                }
            }
            // Check Z edges
            for(auto& slice : zEdgeSlices){
                if(vz >= slice.low && vz < slice.high){
                    fillHistGroup(slice.hg);
                    break;
                }
            }
        }

        // TO HERE!!!

        printf("Chosen True Neutrino: Vertex = (%f, %f, %f), CCNC = %f, PDG = %f, Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz, chosenTrueNeutrino.CCNC, chosenTrueNeutrino.pdg, chosenTrueNeutrino.leptonpdg, chosenTrueNeutrino.tpcID, chosenTrueNeutrino.tpcValid);    
        printf("Chosen True Recoil Electron: Vertex = (%f, %f, %f), PDG = %f, Momentum = (%f, %f, %f), Energy = %f, Angle = %f, ETheta2 = %f, Direction = (%f, %f, %f)\n", chosenTrueParticle.vx, chosenTrueParticle.vy, chosenTrueParticle.vz, chosenTrueParticle.pdg, chosenTrueParticle.px, chosenTrueParticle.py, chosenTrueParticle.pz, chosenTrueParticle.energy, chosenTrueParticle.angle, chosenTrueParticle.ETheta2, chosenTrueParticle.dx, chosenTrueParticle.dy, chosenTrueParticle.dz);
        printf("Number of Slices in event: %f\n", numSlicesInEvent);
    }

    printf("Number of events with a true neutrino within the TPC:\nUBoone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsTrueNeutrino.uboone, numEventsTotal.uboone, numEventsTrueNeutrino.dune, numEventsTotal.dune, numEventsTrueNeutrino.current, numEventsTotal.current, numEventsTrueNeutrino.cheated, numEventsTotal.cheated, numEventsTrueNeutrino.sbnd, numEventsTotal.sbnd);
    printf("Number of events with a true recoil electron:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsTrueElectron.uboone, numEventsTotal.uboone, numEventsTrueElectron.dune, numEventsTotal.dune, numEventsTrueElectron.current, numEventsTotal.current, numEventsTrueElectron.cheated, numEventsTotal.cheated, numEventsTrueElectron.sbnd, numEventsTotal.sbnd);
    printf("Number of events with a slice:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsSlices.uboone, numEventsTotal.uboone, numEventsSlices.dune, numEventsTotal.dune, numEventsSlices.current, numEventsTotal.current, numEventsSlices.cheated, numEventsTotal.cheated, numEventsSlices.sbnd, numEventsTotal.sbnd);
    printf("Number of events with a reco neutrino:\nUboone:%i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsRecoNeutrino.uboone, numEventsTotal.uboone, numEventsRecoNeutrino.dune, numEventsTotal.dune, numEventsRecoNeutrino.current, numEventsTotal.current, numEventsRecoNeutrino.cheated, numEventsTotal.cheated, numEventsRecoNeutrino.sbnd, numEventsTotal.sbnd);
    printf("Number of events where the number of slices with a CRUMBS score != number of reco neutrinos:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numEventsSliceNotEqualNeutrino.uboone, numEventsRecoNeutrino.uboone, numEventsSliceNotEqualNeutrino.dune, numEventsRecoNeutrino.dune, numEventsSliceNotEqualNeutrino.current, numEventsRecoNeutrino.current, numEventsSliceNotEqualNeutrino.cheated, numEventsRecoNeutrino.cheated, numEventsSliceNotEqualNeutrino.sbnd, numEventsRecoNeutrino.sbnd);
    printf("Number of events where the chosen Slice has the highest CRUMBS score and completeness:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", sameSliceSelected.uboone, numEventsSlices.uboone, sameSliceSelected.dune, numEventsSlices.dune, sameSliceSelected.current, numEventsSlices.current, sameSliceSelected.cheated, numEventsSlices.cheated, sameSliceSelected.sbnd, numEventsSlices.sbnd);
    printf("Number of events where the slice with the highest CRUMBS score has a completeness of 0:\nUboone: %i out of %i\nDune: %i out of %i\nCurrent: %i out of %i\nCheated: %i out of %i\nDL SBND: %i out of %i\n", numSlicesCRUMBSCompletenessZero.uboone, numEventsSlices.uboone, numSlicesCRUMBSCompletenessZero.dune, numEventsSlices.dune, numSlicesCRUMBSCompletenessZero.current, numEventsSlices.current, numSlicesCRUMBSCompletenessZero.cheated, numEventsSlices.cheated, numSlicesCRUMBSCompletenessZero.sbnd, numEventsSlices.sbnd);

    styleDraw(numSlices.canvas, numSlices.current, numSlices.cheated, numSlices.dune, numSlices.uboone, numSlices.sbnd, numSlices.nue, numSlices.nue_100k, 999, 999, 999, 999, (base_path + "numSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlices.current, numSlices.cheated, numSlices.dune, numSlices.uboone, numSlices.sbnd, numSlices.nue, numSlices.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "numSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
 
    styleDraw(numSlicesCRUMBS.canvas, numSlicesCRUMBS.current, numSlicesCRUMBS.cheated, numSlicesCRUMBS.dune, numSlicesCRUMBS.uboone, numSlicesCRUMBS.sbnd, numSlicesCRUMBS.nue, numSlicesCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "numCRUMBSSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlicesCRUMBS.current, numSlicesCRUMBS.cheated, numSlicesCRUMBS.dune, numSlicesCRUMBS.uboone, numSlicesCRUMBS.sbnd, numSlicesCRUMBS.nue, numSlicesCRUMBS.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "numCRUMBSSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(numSlicesCompleteness.canvas, numSlicesCompleteness.current, numSlicesCompleteness.cheated, numSlicesCompleteness.dune, numSlicesCompleteness.uboone, numSlicesCompleteness.sbnd, numSlicesCompleteness.nue, numSlicesCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "numCompletenessSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlicesCompleteness.current, numSlicesCompleteness.cheated, numSlicesCompleteness.dune, numSlicesCompleteness.uboone, numSlicesCompleteness.sbnd, numSlicesCompleteness.nue, numSlicesCompleteness.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "numCompletenessSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    
    styleDraw(numRecoNeutrinos.canvas, numRecoNeutrinos.current, numRecoNeutrinos.cheated, numRecoNeutrinos.dune, numRecoNeutrinos.uboone, numRecoNeutrinos.sbnd, numRecoNeutrinos.nue, numRecoNeutrinos.nue_100k, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numRecoNeutrinos.current, numRecoNeutrinos.cheated, numRecoNeutrinos.dune, numRecoNeutrinos.uboone, numRecoNeutrinos.sbnd, numRecoNeutrinos.nue, numRecoNeutrinos.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // CRUMBS Slice Score Plots
    styleDraw(sliceCompletenessCRUMBS.canvas, sliceCompletenessCRUMBS.current, sliceCompletenessCRUMBS.cheated, sliceCompletenessCRUMBS.dune, sliceCompletenessCRUMBS.uboone, sliceCompletenessCRUMBS.sbnd, sliceCompletenessCRUMBS.nue, sliceCompletenessCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceCompletenessCRUMBS.current, sliceCompletenessCRUMBS.cheated, sliceCompletenessCRUMBS.dune, sliceCompletenessCRUMBS.uboone, sliceCompletenessCRUMBS.sbnd, sliceCompletenessCRUMBS.nue, sliceCompletenessCRUMBS.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "sliceCompletenessCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceScoreCRUMBS.canvas, sliceScoreCRUMBS.current, sliceScoreCRUMBS.cheated, sliceScoreCRUMBS.dune, sliceScoreCRUMBS.uboone, sliceScoreCRUMBS.sbnd, sliceScoreCRUMBS.nue, sliceScoreCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceScoreCRUMBS.current, sliceScoreCRUMBS.cheated, sliceScoreCRUMBS.dune, sliceScoreCRUMBS.uboone, sliceScoreCRUMBS.sbnd, sliceScoreCRUMBS.nue, sliceScoreCRUMBS.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(slicePurityCRUMBS.canvas, slicePurityCRUMBS.current, slicePurityCRUMBS.cheated, slicePurityCRUMBS.dune, slicePurityCRUMBS.uboone, slicePurityCRUMBS.sbnd, slicePurityCRUMBS.nue, slicePurityCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "slicePurityCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(slicePurityCRUMBS.current, slicePurityCRUMBS.cheated, slicePurityCRUMBS.dune, slicePurityCRUMBS.uboone, slicePurityCRUMBS.sbnd, slicePurityCRUMBS.nue, slicePurityCRUMBS.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "slicePurityCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    
    // CRUMBS Slice Vertex Plots
    styleDraw(deltaXCRUMBS.canvas, deltaXCRUMBS.current, deltaXCRUMBS.cheated, deltaXCRUMBS.dune, deltaXCRUMBS.uboone, deltaXCRUMBS.sbnd, deltaXCRUMBS.nue, deltaXCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "deltaXCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaYCRUMBS.canvas, deltaYCRUMBS.current, deltaYCRUMBS.cheated, deltaYCRUMBS.dune, deltaYCRUMBS.uboone, deltaYCRUMBS.sbnd, deltaYCRUMBS.nue, deltaYCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "deltaYCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaZCRUMBS.canvas, deltaZCRUMBS.current, deltaZCRUMBS.cheated, deltaZCRUMBS.dune, deltaZCRUMBS.uboone, deltaZCRUMBS.sbnd, deltaZCRUMBS.nue, deltaZCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "deltaZCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaRCRUMBS.canvas, deltaRCRUMBS.current, deltaRCRUMBS.cheated, deltaRCRUMBS.dune, deltaRCRUMBS.uboone, deltaRCRUMBS.sbnd, deltaRCRUMBS.nue, deltaRCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "deltaRCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaXCRUMBS.current, deltaXCRUMBS.cheated, deltaXCRUMBS.dune, deltaXCRUMBS.uboone, deltaXCRUMBS.sbnd, deltaXCRUMBS.nue, deltaXCRUMBS.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaXCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaYCRUMBS.current, deltaYCRUMBS.cheated, deltaYCRUMBS.dune, deltaYCRUMBS.uboone, deltaYCRUMBS.sbnd, deltaYCRUMBS.nue, deltaYCRUMBS.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaYCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaZCRUMBS.current, deltaZCRUMBS.cheated, deltaZCRUMBS.dune, deltaZCRUMBS.uboone, deltaZCRUMBS.sbnd, deltaZCRUMBS.nue, deltaZCRUMBS.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaZCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaRCRUMBS.current, deltaRCRUMBS.cheated, deltaRCRUMBS.dune, deltaRCRUMBS.uboone, deltaRCRUMBS.sbnd, deltaRCRUMBS.nue, deltaRCRUMBS.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaRCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(trueZCRUMBS.canvas, trueZCRUMBS.current, trueZCRUMBS.cheated, trueZCRUMBS.dune, trueZCRUMBS.uboone, trueZCRUMBS.sbnd, trueZCRUMBS.nue, trueZCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "trueZCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(recoZCRUMBS.canvas, recoZCRUMBS.current, recoZCRUMBS.cheated, recoZCRUMBS.dune, recoZCRUMBS.uboone, recoZCRUMBS.sbnd, recoZCRUMBS.nue, recoZCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "recoZCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // CRUMBS PFP Plots
    styleDraw(numPFPsCRUMBS.canvas, numPFPsCRUMBS.current, numPFPsCRUMBS.cheated, numPFPsCRUMBS.dune, numPFPsCRUMBS.uboone, numPFPsCRUMBS.sbnd, numPFPsCRUMBS.nue, numPFPsCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numPFPsCRUMBS.current, numPFPsCRUMBS.cheated, numPFPsCRUMBS.dune, numPFPsCRUMBS.uboone, numPFPsCRUMBS.sbnd, numPFPsCRUMBS.nue, numPFPsCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(ratioChosenSummedEnergyCRUMBS.canvas, ratioChosenSummedEnergyCRUMBS.current, ratioChosenSummedEnergyCRUMBS.cheated, ratioChosenSummedEnergyCRUMBS.dune, ratioChosenSummedEnergyCRUMBS.uboone, ratioChosenSummedEnergyCRUMBS.sbnd, ratioChosenSummedEnergyCRUMBS.nue, ratioChosenSummedEnergyCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenSummedEnergyCRUMBS.current, ratioChosenSummedEnergyCRUMBS.cheated, ratioChosenSummedEnergyCRUMBS.dune, ratioChosenSummedEnergyCRUMBS.uboone, ratioChosenSummedEnergyCRUMBS.sbnd, ratioChosenSummedEnergyCRUMBS.nue, ratioChosenSummedEnergyCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioChosenTrueEnergyCRUMBS.canvas, ratioChosenTrueEnergyCRUMBS.current, ratioChosenTrueEnergyCRUMBS.cheated, ratioChosenTrueEnergyCRUMBS.dune, ratioChosenTrueEnergyCRUMBS.uboone, ratioChosenTrueEnergyCRUMBS.sbnd, ratioChosenTrueEnergyCRUMBS.nue, ratioChosenTrueEnergyCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenTrueEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenTrueEnergyCRUMBS.current, ratioChosenTrueEnergyCRUMBS.cheated, ratioChosenTrueEnergyCRUMBS.dune, ratioChosenTrueEnergyCRUMBS.uboone, ratioChosenTrueEnergyCRUMBS.sbnd, ratioChosenTrueEnergyCRUMBS.nue, ratioChosenTrueEnergyCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenTrueEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioSummedTrueEnergyCRUMBS.canvas, ratioSummedTrueEnergyCRUMBS.current, ratioSummedTrueEnergyCRUMBS.cheated, ratioSummedTrueEnergyCRUMBS.dune, ratioSummedTrueEnergyCRUMBS.uboone, ratioSummedTrueEnergyCRUMBS.sbnd, ratioSummedTrueEnergyCRUMBS.nue, ratioSummedTrueEnergyCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ratioSummedTrueEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioSummedTrueEnergyCRUMBS.current, ratioSummedTrueEnergyCRUMBS.cheated, ratioSummedTrueEnergyCRUMBS.dune, ratioSummedTrueEnergyCRUMBS.uboone, ratioSummedTrueEnergyCRUMBS.sbnd, ratioSummedTrueEnergyCRUMBS.nue, ratioSummedTrueEnergyCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ratioSummedTrueEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    int drawLine = 1;
    int left = 0;
    int right = 1;
    
    styleDraw(angleDifferenceCRUMBS.canvas, angleDifferenceCRUMBS.current, angleDifferenceCRUMBS.cheated, angleDifferenceCRUMBS.dune, angleDifferenceCRUMBS.uboone, angleDifferenceCRUMBS.sbnd, angleDifferenceCRUMBS.nue, angleDifferenceCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "angleDifferenceCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(angleDifferenceCRUMBS.current, angleDifferenceCRUMBS.cheated, angleDifferenceCRUMBS.dune, angleDifferenceCRUMBS.uboone, angleDifferenceCRUMBS.sbnd, angleDifferenceCRUMBS.nue, angleDifferenceCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "angleDifferenceCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
   
    styleDraw(EtrueThetaRecoCRUMBS.canvas, EtrueThetaRecoCRUMBS.current, EtrueThetaRecoCRUMBS.cheated, EtrueThetaRecoCRUMBS.dune, EtrueThetaRecoCRUMBS.sbnd, EtrueThetaRecoCRUMBS.uboone, EtrueThetaRecoCRUMBS.nue, EtrueThetaRecoCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "EtrueThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(EtrueThetaRecoCRUMBS.current, EtrueThetaRecoCRUMBS.cheated, EtrueThetaRecoCRUMBS.dune, EtrueThetaRecoCRUMBS.uboone, EtrueThetaRecoCRUMBS.sbnd, EtrueThetaRecoCRUMBS.nue, EtrueThetaRecoCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "EtrueThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(EtrueThetaRecoCRUMBS.current, EtrueThetaRecoCRUMBS.cheated, EtrueThetaRecoCRUMBS.dune, EtrueThetaRecoCRUMBS.uboone, EtrueThetaRecoCRUMBS.sbnd, EtrueThetaRecoCRUMBS.nue, EtrueThetaRecoCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "EtrueThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{true}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoSumThetaTrueCRUMBS.canvas, ERecoSumThetaTrueCRUMBS.current, ERecoSumThetaTrueCRUMBS.cheated, ERecoSumThetaTrueCRUMBS.dune, ERecoSumThetaTrueCRUMBS.uboone, ERecoSumThetaTrueCRUMBS.sbnd, ERecoSumThetaTrueCRUMBS.nue, ERecoSumThetaTrueCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrueCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaTrueCRUMBS.current, ERecoSumThetaTrueCRUMBS.cheated, ERecoSumThetaTrueCRUMBS.dune, ERecoSumThetaTrueCRUMBS.uboone, ERecoSumThetaTrueCRUMBS.sbnd, ERecoSumThetaTrueCRUMBS.nue, ERecoSumThetaTrueCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrueCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaTrueCRUMBS.current, ERecoSumThetaTrueCRUMBS.cheated, ERecoSumThetaTrueCRUMBS.dune, ERecoSumThetaTrueCRUMBS.uboone, ERecoSumThetaTrueCRUMBS.sbnd, ERecoSumThetaTrueCRUMBS.nue, ERecoSumThetaTrueCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoSumThetaTrueCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoHighestThetaTrueCRUMBS.canvas, ERecoHighestThetaTrueCRUMBS.current, ERecoHighestThetaTrueCRUMBS.cheated, ERecoHighestThetaTrueCRUMBS.dune, ERecoHighestThetaTrueCRUMBS.uboone, ERecoHighestThetaTrueCRUMBS.sbnd, ERecoHighestThetaTrueCRUMBS.nue, ERecoHighestThetaTrueCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrueCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaTrueCRUMBS.current, ERecoHighestThetaTrueCRUMBS.cheated, ERecoHighestThetaTrueCRUMBS.dune, ERecoHighestThetaTrueCRUMBS.uboone, ERecoHighestThetaTrueCRUMBS.sbnd, ERecoHighestThetaTrueCRUMBS.nue, ERecoHighestThetaTrueCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrueCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaTrueCRUMBS.current, ERecoHighestThetaTrueCRUMBS.cheated, ERecoHighestThetaTrueCRUMBS.dune, ERecoHighestThetaTrueCRUMBS.uboone, ERecoHighestThetaTrueCRUMBS.sbnd, ERecoHighestThetaTrueCRUMBS.nue, ERecoHighestThetaTrueCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoHighestThetaTrueCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoSumThetaRecoCRUMBS.canvas, ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.cheated, ERecoSumThetaRecoCRUMBS.dune, ERecoSumThetaRecoCRUMBS.uboone, ERecoSumThetaRecoCRUMBS.sbnd, ERecoSumThetaRecoCRUMBS.nue, ERecoSumThetaRecoCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.cheated, ERecoSumThetaRecoCRUMBS.dune, ERecoSumThetaRecoCRUMBS.uboone, ERecoSumThetaRecoCRUMBS.sbnd, ERecoSumThetaRecoCRUMBS.nue, ERecoSumThetaRecoCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.cheated, ERecoSumThetaRecoCRUMBS.dune, ERecoSumThetaRecoCRUMBS.uboone, ERecoSumThetaRecoCRUMBS.sbnd, ERecoSumThetaRecoCRUMBS.nue, ERecoSumThetaRecoCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoHighestThetaRecoCRUMBS.canvas, ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.cheated, ERecoHighestThetaRecoCRUMBS.dune, ERecoHighestThetaRecoCRUMBS.uboone, ERecoHighestThetaRecoCRUMBS.sbnd, ERecoHighestThetaRecoCRUMBS.nue, ERecoHighestThetaRecoCRUMBS.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.cheated, ERecoHighestThetaRecoCRUMBS.dune, ERecoHighestThetaRecoCRUMBS.uboone, ERecoHighestThetaRecoCRUMBS.sbnd, ERecoHighestThetaRecoCRUMBS.nue, ERecoHighestThetaRecoCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.cheated, ERecoHighestThetaRecoCRUMBS.dune, ERecoHighestThetaRecoCRUMBS.uboone, ERecoHighestThetaRecoCRUMBS.sbnd, ERecoHighestThetaRecoCRUMBS.nue, ERecoHighestThetaRecoCRUMBS.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    // Completeness Slice Score Plots
    styleDraw(sliceCompletenessCompleteness.canvas, sliceCompletenessCompleteness.current, sliceCompletenessCompleteness.cheated, sliceCompletenessCompleteness.dune, sliceCompletenessCompleteness.uboone, sliceCompletenessCompleteness.sbnd, sliceCompletenessCompleteness.nue, sliceCompletenessCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "sliceCompletenessCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceCompletenessCompleteness.current, sliceCompletenessCompleteness.cheated, sliceCompletenessCompleteness.dune, sliceCompletenessCompleteness.uboone, sliceCompletenessCompleteness.sbnd, sliceCompletenessCompleteness.nue, sliceCompletenessCompleteness.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "sliceCompletenessCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(sliceScoreCompleteness.canvas, sliceScoreCompleteness.current, sliceScoreCompleteness.cheated, sliceScoreCompleteness.dune, sliceScoreCompleteness.uboone, sliceScoreCompleteness.sbnd, sliceScoreCompleteness.nue, sliceScoreCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "sliceScoreCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceScoreCompleteness.current, sliceScoreCompleteness.cheated, sliceScoreCompleteness.dune, sliceScoreCompleteness.uboone, sliceScoreCompleteness.sbnd, sliceScoreCompleteness.nue, sliceScoreCompleteness.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "sliceScoreCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(slicePurityCompleteness.canvas, slicePurityCompleteness.current, slicePurityCompleteness.cheated, slicePurityCompleteness.dune, slicePurityCompleteness.uboone, slicePurityCompleteness.sbnd, slicePurityCompleteness.nue, slicePurityCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "slicePurityCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(slicePurityCompleteness.current, slicePurityCompleteness.cheated, slicePurityCompleteness.dune, slicePurityCompleteness.uboone, slicePurityCompleteness.sbnd, slicePurityCompleteness.nue, slicePurityCompleteness.nue_100k, numEventsSlices.current, numEventsSlices.cheated, numEventsSlices.dune, numEventsSlices.uboone, numEventsSlices.sbnd, numEventsSlices.nue, numEventsSlices.nue_100k, 999, 999, 999, 999, (base_path + "slicePurityCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    styleDraw(deltaXCompleteness.canvas, deltaXCompleteness.current, deltaXCompleteness.cheated, deltaXCompleteness.dune, deltaXCompleteness.uboone, deltaXCompleteness.sbnd, deltaXCompleteness.nue, deltaXCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "deltaXCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaYCompleteness.canvas, deltaYCompleteness.current, deltaYCompleteness.cheated, deltaYCompleteness.dune, deltaYCompleteness.uboone, deltaYCompleteness.sbnd, deltaYCompleteness.nue, deltaYCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "deltaYCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaZCompleteness.canvas, deltaZCompleteness.current, deltaZCompleteness.cheated, deltaZCompleteness.dune, deltaZCompleteness.uboone, deltaZCompleteness.sbnd, deltaZCompleteness.nue, deltaZCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "deltaZCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaRCompleteness.canvas, deltaRCompleteness.current, deltaRCompleteness.cheated, deltaRCompleteness.dune, deltaRCompleteness.uboone, deltaRCompleteness.sbnd, deltaRCompleteness.nue, deltaRCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "deltaRCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaXCompleteness.current, deltaXCompleteness.cheated, deltaXCompleteness.dune, deltaXCompleteness.uboone, deltaXCompleteness.sbnd, deltaXCompleteness.nue, deltaXCompleteness.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaXCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaYCompleteness.current, deltaYCompleteness.cheated, deltaYCompleteness.dune, deltaYCompleteness.uboone, deltaYCompleteness.sbnd, deltaYCompleteness.nue, deltaYCompleteness.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaYCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaZCompleteness.current, deltaZCompleteness.cheated, deltaZCompleteness.dune, deltaZCompleteness.uboone, deltaZCompleteness.sbnd, deltaZCompleteness.nue, deltaZCompleteness.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaZCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaRCompleteness.current, deltaRCompleteness.cheated, deltaRCompleteness.dune, deltaRCompleteness.uboone, deltaRCompleteness.sbnd, deltaRCompleteness.nue, deltaRCompleteness.nue_100k, numEventsRecoNeutrino.current, numEventsRecoNeutrino.cheated, numEventsRecoNeutrino.dune, numEventsRecoNeutrino.uboone, numEventsRecoNeutrino.sbnd, numEventsRecoNeutrino.nue, numEventsRecoNeutrino.nue_100k, 999, 999, 999, 999, (base_path + "deltaRCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // Completeness PFP Plots
    styleDraw(numPFPsCompleteness.canvas, numPFPsCompleteness.current, numPFPsCompleteness.cheated, numPFPsCompleteness.dune, numPFPsCompleteness.uboone, numPFPsCompleteness.sbnd, numPFPsCompleteness.nue, numPFPsCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "numPFPsCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numPFPsCompleteness.current, numPFPsCompleteness.cheated, numPFPsCompleteness.dune, numPFPsCompleteness.uboone, numPFPsCompleteness.sbnd, numPFPsCompleteness.nue, numPFPsCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "numPFPsCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(ratioChosenSummedEnergyCompleteness.canvas, ratioChosenSummedEnergyCompleteness.current, ratioChosenSummedEnergyCompleteness.cheated, ratioChosenSummedEnergyCompleteness.dune, ratioChosenSummedEnergyCompleteness.uboone, ratioChosenSummedEnergyCompleteness.sbnd, ratioChosenSummedEnergyCompleteness.nue, ratioChosenSummedEnergyCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenSummedEnergyCompleteness.current, ratioChosenSummedEnergyCompleteness.cheated, ratioChosenSummedEnergyCompleteness.dune, ratioChosenSummedEnergyCompleteness.uboone, ratioChosenSummedEnergyCompleteness.sbnd, ratioChosenSummedEnergyCompleteness.nue, ratioChosenSummedEnergyCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioChosenTrueEnergyCompleteness.canvas, ratioChosenTrueEnergyCompleteness.current, ratioChosenTrueEnergyCompleteness.cheated, ratioChosenTrueEnergyCompleteness.dune, ratioChosenTrueEnergyCompleteness.uboone, ratioChosenTrueEnergyCompleteness.sbnd, ratioChosenTrueEnergyCompleteness.nue, ratioChosenTrueEnergyCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenTrueEnergyCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenTrueEnergyCompleteness.current, ratioChosenTrueEnergyCompleteness.cheated, ratioChosenTrueEnergyCompleteness.dune, ratioChosenTrueEnergyCompleteness.uboone, ratioChosenTrueEnergyCompleteness.sbnd, ratioChosenTrueEnergyCompleteness.nue, ratioChosenTrueEnergyCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ratioChosenTrueEnergyCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    styleDraw(ratioSummedTrueEnergyCompleteness.canvas, ratioSummedTrueEnergyCompleteness.current, ratioSummedTrueEnergyCompleteness.cheated, ratioSummedTrueEnergyCompleteness.dune, ratioSummedTrueEnergyCompleteness.uboone, ratioSummedTrueEnergyCompleteness.sbnd, ratioSummedTrueEnergyCompleteness.nue, ratioSummedTrueEnergyCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ratioSummedTrueEnergyCompleteness_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioSummedTrueEnergyCompleteness.current, ratioSummedTrueEnergyCompleteness.cheated, ratioSummedTrueEnergyCompleteness.dune, ratioSummedTrueEnergyCompleteness.uboone, ratioSummedTrueEnergyCompleteness.sbnd, ratioSummedTrueEnergyCompleteness.nue, ratioSummedTrueEnergyCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ratioSummedTrueEnergyCompleteness_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    styleDraw(angleDifferenceCompleteness.canvas, angleDifferenceCompleteness.current, angleDifferenceCompleteness.cheated, angleDifferenceCompleteness.dune, angleDifferenceCompleteness.uboone, angleDifferenceCompleteness.sbnd, angleDifferenceCompleteness.nue, angleDifferenceCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "angleDifferenceCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(angleDifferenceCompleteness.current, angleDifferenceCompleteness.cheated, angleDifferenceCompleteness.dune, angleDifferenceCompleteness.uboone, angleDifferenceCompleteness.sbnd, angleDifferenceCompleteness.nue, angleDifferenceCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "angleDifferenceCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
   
    styleDraw(EtrueThetaRecoCompleteness.canvas, EtrueThetaRecoCompleteness.current, EtrueThetaRecoCompleteness.cheated, EtrueThetaRecoCompleteness.dune, EtrueThetaRecoCompleteness.uboone, EtrueThetaRecoCompleteness.sbnd, EtrueThetaRecoCompleteness.nue, EtrueThetaRecoCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "EtrueThetaRecoCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(EtrueThetaRecoCompleteness.current, EtrueThetaRecoCompleteness.cheated, EtrueThetaRecoCompleteness.dune, EtrueThetaRecoCompleteness.uboone, EtrueThetaRecoCompleteness.sbnd, EtrueThetaRecoCompleteness.nue, EtrueThetaRecoCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "EtrueThetaRecoCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(EtrueThetaRecoCompleteness.current, EtrueThetaRecoCompleteness.cheated, EtrueThetaRecoCompleteness.dune, EtrueThetaRecoCompleteness.uboone, EtrueThetaRecoCompleteness.sbnd, EtrueThetaRecoCompleteness.nue, EtrueThetaRecoCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "EtrueThetaRecoCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{true}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoSumThetaTrueCompleteness.canvas, ERecoSumThetaTrueCompleteness.current, ERecoSumThetaTrueCompleteness.cheated, ERecoSumThetaTrueCompleteness.dune, ERecoSumThetaTrueCompleteness.uboone, ERecoSumThetaTrueCompleteness.sbnd, ERecoSumThetaTrueCompleteness.nue, ERecoSumThetaTrueCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrueCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaTrueCompleteness.current, ERecoSumThetaTrueCompleteness.cheated, ERecoSumThetaTrueCompleteness.dune, ERecoSumThetaTrueCompleteness.uboone, ERecoSumThetaTrueCompleteness.sbnd, ERecoSumThetaTrueCompleteness.nue, ERecoSumThetaTrueCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaTrueCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaTrueCompleteness.current, ERecoSumThetaTrueCompleteness.cheated, ERecoSumThetaTrueCompleteness.dune, ERecoSumThetaTrueCompleteness.uboone, ERecoSumThetaTrueCompleteness.sbnd, ERecoSumThetaTrueCompleteness.nue, ERecoSumThetaTrueCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoSumThetaTrueCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoHighestThetaTrueCompleteness.canvas, ERecoHighestThetaTrueCompleteness.current, ERecoHighestThetaTrueCompleteness.cheated, ERecoHighestThetaTrueCompleteness.dune, ERecoHighestThetaTrueCompleteness.uboone, ERecoHighestThetaTrueCompleteness.sbnd, ERecoHighestThetaTrueCompleteness.nue, ERecoHighestThetaTrueCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrueCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaTrueCompleteness.current, ERecoHighestThetaTrueCompleteness.cheated, ERecoHighestThetaTrueCompleteness.dune, ERecoHighestThetaTrueCompleteness.uboone, ERecoHighestThetaTrueCompleteness.sbnd, ERecoHighestThetaTrueCompleteness.nue, ERecoHighestThetaTrueCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaTrueCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaTrueCompleteness.current, ERecoHighestThetaTrueCompleteness.cheated, ERecoHighestThetaTrueCompleteness.dune, ERecoHighestThetaTrueCompleteness.uboone, ERecoHighestThetaTrueCompleteness.sbnd, ERecoHighestThetaTrueCompleteness.nue, ERecoHighestThetaTrueCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoHighestThetaTrueCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{true}^{2} (MeV)");

    styleDraw(ERecoSumThetaRecoCompleteness.canvas, ERecoSumThetaRecoCompleteness.current, ERecoSumThetaRecoCompleteness.cheated, ERecoSumThetaRecoCompleteness.dune, ERecoSumThetaRecoCompleteness.uboone, ERecoSumThetaRecoCompleteness.sbnd, ERecoSumThetaRecoCompleteness.nue, ERecoSumThetaRecoCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaRecoCompleteness.current, ERecoSumThetaRecoCompleteness.cheated, ERecoSumThetaRecoCompleteness.dune, ERecoSumThetaRecoCompleteness.uboone, ERecoSumThetaRecoCompleteness.sbnd, ERecoSumThetaRecoCompleteness.nue, ERecoSumThetaRecoCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaRecoCompleteness.current, ERecoSumThetaRecoCompleteness.cheated, ERecoSumThetaRecoCompleteness.dune, ERecoSumThetaRecoCompleteness.uboone, ERecoSumThetaRecoCompleteness.sbnd, ERecoSumThetaRecoCompleteness.nue, ERecoSumThetaRecoCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoSumThetaRecoCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoHighestThetaRecoCompleteness.canvas, ERecoHighestThetaRecoCompleteness.current, ERecoHighestThetaRecoCompleteness.cheated, ERecoHighestThetaRecoCompleteness.dune, ERecoHighestThetaRecoCompleteness.uboone, ERecoHighestThetaRecoCompleteness.sbnd, ERecoHighestThetaRecoCompleteness.nue, ERecoHighestThetaRecoCompleteness.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCompleteness_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaRecoCompleteness.current, ERecoHighestThetaRecoCompleteness.cheated, ERecoHighestThetaRecoCompleteness.dune, ERecoHighestThetaRecoCompleteness.uboone, ERecoHighestThetaRecoCompleteness.sbnd, ERecoHighestThetaRecoCompleteness.nue, ERecoHighestThetaRecoCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCompleteness_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaRecoCompleteness.current, ERecoHighestThetaRecoCompleteness.cheated, ERecoHighestThetaRecoCompleteness.dune, ERecoHighestThetaRecoCompleteness.uboone, ERecoHighestThetaRecoCompleteness.sbnd, ERecoHighestThetaRecoCompleteness.nue, ERecoHighestThetaRecoCompleteness.nue_100k, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, numEventsCRUMBSRecoParticle.nue, numEventsCRUMBSRecoParticle.nue_100k, 0, 1, 999, 999, (base_path + "ERecoHighestThetaRecoCompleteness_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");


    styleDraw(trueETheta2.canvas, trueETheta2.current, trueETheta2.cheated, trueETheta2.dune, trueETheta2.uboone, trueETheta2.sbnd, trueETheta2.nue, trueETheta2.nue_100k, 999, 999, 999, 999, (base_path + "trueETheta2_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(trueETheta2.current, trueETheta2.cheated, trueETheta2.dune, trueETheta2.uboone, trueETheta2.sbnd, trueETheta2.nue, trueETheta2.nue_100k, numEventsTrueElectron.current, numEventsTrueElectron.cheated, numEventsTrueElectron.dune, numEventsTrueElectron.uboone, numEventsTrueElectron.sbnd, numEventsTrueElectron.nue, numEventsTrueElectron.nue_100k, 999, 999, 999, 999, (base_path + "trueETheta2_perc.pdf").c_str(), 0.56, 0.88, 0.70, 0.86, &drawLine, &right);
    efficiency(trueETheta2.current, trueETheta2.cheated, trueETheta2.dune, trueETheta2.uboone, trueETheta2.sbnd, trueETheta2.nue, trueETheta2.nue_100k, numEventsTrueElectron.current, numEventsTrueElectron.cheated, numEventsTrueElectron.dune, numEventsTrueElectron.uboone, numEventsTrueElectron.sbnd, numEventsTrueElectron.nue, numEventsTrueElectron.nue_100k, 0, 1, 999, 999, (base_path + "trueETheta2_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{true}#theta_{true}^{2} (MeV rad^{2})");
   
    if(makePerSlicePlots){
        for(auto& slice : xEdgeSlices){
            styleDraw(slice.hg.canvas, slice.hg.current, slice.hg.cheated, slice.hg.dune, slice.hg.uboone, slice.hg.sbnd,
                      slice.hg.nue, slice.hg.nue_100k, 999, 999, 0, 180,
                      (base_path_perSlice + std::string(slice.hg.baseHist->GetName()) + ".pdf").c_str(),
                       0.6, 0.9, 0.7, 0.9);
            std::cout << "slice.hg.current entries: " << slice.hg.current->GetBinContent(3) << std::endl;
            std::cout << "Integral (all bins): " << slice.hg.current->Integral(0, slice.hg.current->GetNbinsX()+1) << std::endl;
        }
        for(auto& slice : yEdgeSlices){
            styleDraw(slice.hg.canvas, slice.hg.current, slice.hg.cheated, slice.hg.dune, slice.hg.uboone, slice.hg.sbnd,
                      slice.hg.nue, slice.hg.nue_100k, 999, 999, 0, 180,
                      (base_path_perSlice + std::string(slice.hg.baseHist->GetName()) + ".pdf").c_str(),
                       0.6, 0.9, 0.7, 0.9);
        }
        for(auto& slice : zEdgeSlices){
            styleDraw(slice.hg.canvas, slice.hg.current, slice.hg.cheated, slice.hg.dune, slice.hg.uboone, slice.hg.sbnd,
                      slice.hg.nue, slice.hg.nue_100k, 999, 999, 0, 180,
                      (base_path_perSlice + std::string(slice.hg.baseHist->GetName()) + ".pdf").c_str(),
                       0.6, 0.9, 0.7, 0.9);
        }
    }

    if(makePerSlicePlots){
        
        // Summary Plots for X Coordinate: BDT Vertexing
        std::vector<double> xvals, yPeaks, yWidths;
        for(const auto& slice : xEdgeSlices){
            TH1F* hist = slice.hg.current;
            double peakCenter = -1.0;
            double width = -1.0;
            std::cout << "slice.hg.current entries later: " << slice.hg.current->GetEntries() << std::endl;
            if(slice.hg.current->GetEntries() > 0){
                int peakBin = slice.hg.current->GetMaximumBin();
                peakCenter = hist->GetBinCenter(peakBin);
                width = hist->GetRMS();
            }
            double center = 0.5 * (slice.low + slice.high);
            xvals.push_back(center);
            yPeaks.push_back(peakCenter);
            yWidths.push_back(width);
        }

        
        TGraph* gXpeak = new TGraph(xvals.size(), xvals.data(), yPeaks.data());
        gXpeak->SetTitle("X Peak Positions: BDT Vertexing;X Coordinate [cm];Angle Difference Peak [deg]");
        gXpeak->SetMarkerStyle(20);
        TCanvas* cXpeak = new TCanvas("cXpeak", "X slice peaks", 800, 600);
        gXpeak->Draw("AP");
        //cXpeak->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/X_BDT_peaks.pdf");
        cXpeak->SaveAs((base_path + "X_BDT_peaks.pdf").c_str());
        
        TGraph* gXwidth = new TGraph(xvals.size(), xvals.data(), yWidths.data());
        gXwidth->SetTitle("X Peak Widths: BDT Vertexing;X Coordinate [cm];Angle Difference RMS [deg]");
        gXwidth->SetMarkerStyle(21);
        TCanvas* cXwidth = new TCanvas("cXwidth", "X slice widths", 800, 600);
        gXwidth->Draw("AP");
        cXwidth->SaveAs((base_path + "X_BDT_widths.pdf").c_str());
        //cXwidth->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/X_BDT_widths.pdf");

        // Summary Plots for Y Coordinate: BDT Vertexing
        xvals.clear(); yPeaks.clear(); yWidths.clear();
        for(const auto& slice : yEdgeSlices){
            TH1F* hist = slice.hg.current;
            double peakCenter = -1.0;
            double width = -1.0;
            if(hist->GetEntries() > 0){
                int peakBin = hist->GetMaximumBin();
                peakCenter = hist->GetBinCenter(peakBin);
                width = hist->GetRMS();
            }
            double center = 0.5 * (slice.low + slice.high);
            xvals.push_back(center);
            yPeaks.push_back(peakCenter);
            yWidths.push_back(width);
        }
        
        TGraph* gYpeak = new TGraph(xvals.size(), xvals.data(), yPeaks.data());
        gYpeak->SetTitle("Y Peak Positions: BDT Vertexing;Y Coordinate [cm];Angle Difference Peak [deg]");
        gYpeak->SetMarkerStyle(20);
        TCanvas* cYpeak = new TCanvas("cYpeak", "Y slice peaks", 800, 600);
        gYpeak->Draw("AP");
        //cYpeak->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Y_BDT_peaks.pdf");
        cYpeak->SaveAs((base_path + "Y_BDT_peaks.pdf").c_str());
        
        TGraph* gYwidth = new TGraph(xvals.size(), xvals.data(), yWidths.data());
        gYwidth->SetTitle("Y Peak Widths: BDT Vertexing;Y Coordinate [cm];Angle Difference RMS [deg]");
        gYwidth->SetMarkerStyle(21);
        TCanvas* cYwidth = new TCanvas("cYwidth", "Y slice widths", 800, 600);
        gYwidth->Draw("AP");
        //cYwidth->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Y_BDT_widths.pdf");
        cYwidth->SaveAs((base_path + "Y_BDT_width.pdf").c_str());

        // Summary Plots for Z Coordinate
        xvals.clear(); yPeaks.clear(); yWidths.clear();
        for(const auto& slice : zEdgeSlices){
            TH1F* hist = slice.hg.current;
            double peakCenter = -1.0;
            double width = -1.0;
            if(hist->GetEntries() > 0){
                int peakBin = hist->GetMaximumBin();
                peakCenter = hist->GetBinCenter(peakBin);
                width = hist->GetRMS();
            }
            double center = 0.5 * (slice.low + slice.high);
            xvals.push_back(center);
            yPeaks.push_back(peakCenter);
            yWidths.push_back(width);
        }
        
        TGraph* gZpeak = new TGraph(xvals.size(), xvals.data(), yPeaks.data());
        gZpeak->SetTitle("Z Peak Positions: BDT Vertexing;Z coordinate [cm];Angle Difference Peak [deg]");
        gZpeak->SetMarkerStyle(20);
        TCanvas* cZpeak = new TCanvas("cZpeak", "Z slice peaks", 800, 600);
        gZpeak->Draw("AP");
        //cZpeak->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Z_BDT_peaks.pdf");
        cZpeak->SaveAs((base_path + "Z_BDT_peaks.pdf").c_str());
        
        TGraph* gZwidth = new TGraph(xvals.size(), xvals.data(), yWidths.data());
        gZwidth->SetTitle("Z Peak Widths (RMS): BDT Vertexing;Z coordinate [cm];Angle Difference RMS [deg]");
        gZwidth->SetMarkerStyle(21);
        TCanvas* cZwidth = new TCanvas("cZwidth", "Z slice widths", 800, 600);
        gZwidth->Draw("AP");
        //cZwidth->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Z_BDT_widths.pdf");
        cZwidth->SaveAs((base_path + "Z_BDT_widths.pdf").c_str());
    
    
        // Summary Plots for X Coordinate: DL Uboone Vertexing
        std::vector<double> xvalsUboone, yPeaksUboone, yWidthsUboone;
        for(const auto& slice : xEdgeSlices){
            TH1F* histUboone = slice.hg.uboone;
            double peakCenter = -1.0;
            double width = -1.0;
            if(slice.hg.uboone->GetEntries() > 0){
                int peakBin = slice.hg.uboone->GetMaximumBin();
                peakCenter = histUboone->GetBinCenter(peakBin);
                width = histUboone->GetRMS();
            }
            double center = 0.5 * (slice.low + slice.high);
            xvalsUboone.push_back(center);
            yPeaksUboone.push_back(peakCenter);
            yWidthsUboone.push_back(width);
        }
        
        TGraph* gXpeakUboone = new TGraph(xvalsUboone.size(), xvalsUboone.data(), yPeaksUboone.data());
        gXpeakUboone->SetTitle("X Peak Positions: DL Uboone Vertexing;X Coordinate [cm];Angle Difference Peak [deg]");
        gXpeakUboone->SetMarkerStyle(20);
        TCanvas* cXpeakUboone = new TCanvas("cXpeakUboone", "X slice peaks", 800, 600);
        gXpeakUboone->Draw("AP");
        //cXpeakUboone->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/X_DLUboone_peaks.pdf");
        cXpeakUboone->SaveAs((base_path + "X_DLUboone_peaks.pdf").c_str());
        
        TGraph* gXwidthUboone = new TGraph(xvalsUboone.size(), xvalsUboone.data(), yWidthsUboone.data());
        gXwidthUboone->SetTitle("X Peak Widths: DL Uboone Vertexing;X Coordinate [cm];Angle Difference RMS [deg]");
        gXwidthUboone->SetMarkerStyle(21);
        TCanvas* cXwidthUboone = new TCanvas("cXwidthUboone", "X slice widths", 800, 600);
        gXwidthUboone->Draw("AP");
        //cXwidthUboone->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/X_DLUboone_widths.pdf");
        cXwidthUboone->SaveAs((base_path + "X_DLUboone_widths.pdf").c_str());

        // Summary Plots for Y Coordinate: DL Uboone Vertexing
        xvalsUboone.clear(); yPeaksUboone.clear(); yWidthsUboone.clear();
        for(const auto& slice : yEdgeSlices){
            TH1F* histUboone = slice.hg.uboone;
            double peakCenter = -1.0;
            double width = -1.0;
            if(histUboone->GetEntries() > 0){
                int peakBin = histUboone->GetMaximumBin();
                peakCenter = histUboone->GetBinCenter(peakBin);
                width = histUboone->GetRMS();
            }
            double center = 0.5 * (slice.low + slice.high);
            xvalsUboone.push_back(center);
            yPeaksUboone.push_back(peakCenter);
            yWidthsUboone.push_back(width);
        }
        
        TGraph* gYpeakUboone = new TGraph(xvalsUboone.size(), xvalsUboone.data(), yPeaksUboone.data());
        gYpeakUboone->SetTitle("Y Peak Positions: DL Uboone Vertexing;Y Coordinate [cm];Angle Difference Peak [deg]");
        gYpeakUboone->SetMarkerStyle(20);
        TCanvas* cYpeakUboone = new TCanvas("cYpeakUboone", "Y slice peaks", 800, 600);
        gYpeakUboone->Draw("AP");
        //cYpeakUboone->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Y_DLUboone_peaks.pdf");
        cYpeakUboone->SaveAs((base_path + "Y_DLUboone_peaks.pdf").c_str());
        
        TGraph* gYwidthUboone = new TGraph(xvalsUboone.size(), xvalsUboone.data(), yWidthsUboone.data());
        gYwidthUboone->SetTitle("Y Peak Widths: DL Uboone Vertexing;Y Coordinate [cm];Angle Difference RMS [deg]");
        gYwidthUboone->SetMarkerStyle(21);
        TCanvas* cYwidthUboone = new TCanvas("cYwidthUboone", "Y slice widths", 800, 600);
        gYwidthUboone->Draw("AP");
        //cYwidthUboone->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Y_DLUboone_widths.pdf");
        cYwidthUboone->SaveAs((base_path + "Y_DLUboone_widths.pdf").c_str());

        // Summary Plots for Z Coordinate
        xvalsUboone.clear(); yPeaksUboone.clear(); yWidthsUboone.clear();
        for(const auto& slice : zEdgeSlices){
            TH1F* histUboone = slice.hg.uboone;
            double peakCenter = -1.0;
            double width = -1.0;
            if(histUboone->GetEntries() > 0){
                int peakBin = histUboone->GetMaximumBin();
                peakCenter = histUboone->GetBinCenter(peakBin);
                width = histUboone->GetRMS();
            }
            double center = 0.5 * (slice.low + slice.high);
            xvalsUboone.push_back(center);
            yPeaksUboone.push_back(peakCenter);
            yWidthsUboone.push_back(width);
        }
        
        TGraph* gZpeakUboone = new TGraph(xvalsUboone.size(), xvalsUboone.data(), yPeaksUboone.data());
        gZpeakUboone->SetTitle("Z Peak Positions: DL Uboone Vertexing;Z coordinate [cm];Angle Difference Peak [deg]");
        gZpeakUboone->SetMarkerStyle(20);
        TCanvas* cZpeakUboone = new TCanvas("cZpeakUboone", "Z slice peaks", 800, 600);
        gZpeakUboone->Draw("AP");
        //cZpeakUboone->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Z_DLUboone_peaks.pdf");
        cZpeakUboone->SaveAs((base_path + "Z_DLUboone_peaks.pdf").c_str());
        
        TGraph* gZwidthUboone = new TGraph(xvalsUboone.size(), xvalsUboone.data(), yWidthsUboone.data());
        gZwidthUboone->SetTitle("Z Peak Widths (RMS): DL Uboone Vertexing;Z coordinate [cm];Angle Difference RMS [deg]");
        gZwidthUboone->SetMarkerStyle(21);
        TCanvas* cZwidthUboone = new TCanvas("cZwidthUboone", "Z slice widths", 800, 600);
        gZwidthUboone->Draw("AP");
        //cZwidthUboone->SaveAs("/nashome/c/coackley/nuEPlotsWithoutCosmicsSCEOFF_50k+100k/Z_DLUboone_widths.pdf");
        cZwidthUboone->SaveAs((base_path + "Z_DLUboone_widths.pdf").c_str());
    }

    TwoDHistDraw(xCoordAngleDifferenceBDTCRUMBS, (base_path + "angleDiffPosition_x_BDT.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceBDTCRUMBS, (base_path + "angleDiffPosition_y_BDT.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceBDTCRUMBS, (base_path + "angleDiffPosition_z_BDT.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(xCoordAngleDifferenceBDTCRUMBS_low, (base_path + "angleDiffPosition_x_BDT_low.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceBDTCRUMBS_low, (base_path + "angleDiffPosition_y_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceBDTCRUMBS_low, (base_path + "angleDiffPosition_z_BDT_low.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(xCoordAngleDifferenceBDTCRUMBS_high, (base_path + "angleDiffPosition_x_BDT_high.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceBDTCRUMBS_high, (base_path + "angleDiffPosition_y_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceBDTCRUMBS_high, (base_path + "angleDiffPosition_z_BDT_high.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");

    ProfileDraw(xCoordAngleDifferenceBDTCRUMBSProfile, (base_path + "angleDiffPositionProfile_x_BDT.pdf").c_str(), "Profile of Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(yCoordAngleDifferenceBDTCRUMBSProfile, (base_path + "angleDiffPositionProfile_y_BDT.pdf").c_str(), "Profile of Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(zCoordAngleDifferenceBDTCRUMBSProfile, (base_path + "angleDiffPositionProfile_z_BDT.pdf").c_str(), "Profile of Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(xCoordAngleDifferenceBDTCRUMBSProfile_low, (base_path + "angleDiffPositionProfile_x_BDT_low.pdf").c_str(), "Profile of Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(yCoordAngleDifferenceBDTCRUMBSProfile_low, (base_path + "angleDiffPositionProfile_y_BDT_low.pdf").c_str(), "Profile of Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(zCoordAngleDifferenceBDTCRUMBSProfile_low, (base_path + "angleDiffPositionProfile_z_BDT_low.pdf").c_str(), "Profile of Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(xCoordAngleDifferenceBDTCRUMBSProfile_high, (base_path + "angleDiffPositionProfile_x_BDT_high.pdf").c_str(), "Profile of Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(yCoordAngleDifferenceBDTCRUMBSProfile_high, (base_path + "angleDiffPositionProfile_y_BDT_high.pdf").c_str(), "Profile of Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    ProfileDraw(zCoordAngleDifferenceBDTCRUMBSProfile_high, (base_path + "angleDiffPositionProfile_z_BDT_high.pdf").c_str(), "Profile of Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: BDT Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");

    TwoDHistDraw(xCoordAngleDifferenceDLUbooneCRUMBS, (base_path + "angleDiffPosition_x_DLUboone.pdf").c_str(), "Reco Neutrino Vertex X Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex X Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(yCoordAngleDifferenceDLUbooneCRUMBS, (base_path + "angleDiffPosition_y_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Y Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Y Coordinate (cm);Angle Difference (degrees)");
    TwoDHistDraw(zCoordAngleDifferenceDLUbooneCRUMBS, (base_path + "angleDiffPosition_z_DLUboone.pdf").c_str(), "Reco Neutrino Vertex Z Coordinate vs Angle Between True and Reco Track: DL Uboone Vertexing;Reco Neutrino Vertex Z Coordinate (cm);Angle Difference (degrees)");

    makeEqualStatProfile(xVals_xCoord, yVals, 25, -201.3, 201.3, "xCoordAngleDifferenceBDTCRUMBSProfile_equalStats", "Profile of Reco Neutrino Vertex X Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino X Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_x_BDT_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_yCoord, yVals, 25, -203.8, 203.8, "yCoordAngleDifferenceBDTCRUMBSProfile_equalStats", "Profile of Reco Neutrino Vertex Y Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino Y Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_y_BDT_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_zCoord, yVals, 25, 0, 509.5, "zCoordAngleDifferenceBDTCRUMBSProfile_equalStats", "Profile of Reco Neutrino Vertex Z Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino Z Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_z_BDT_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_xCoord, yVals, 25, -201.3, -171.3, "xCoordAngleDifferenceBDTCRUMBSProfile_low_equalStats", "Profile of Reco Neutrino Vertex X Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino X Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_x_BDT_low_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_yCoord, yVals, 25, -203.8, -173.8, "yCoordAngleDifferenceBDTCRUMBSProfile_low_equalStats", "Profile of Reco Neutrino Vertex Y Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino Y Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_y_BDT_low_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_zCoord, yVals, 25, 0, 40, "zCoordAngleDifferenceBDTCRUMBSProfile_low_equalStats", "Profile of Reco Neutrino Vertex Z Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino Z Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_z_BDT_low_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_xCoord, yVals, 27, 171.3, 201.3, "xCoordAngleDifferenceBDTCRUMBSProfile_high_equalStats", "Profile of Reco Neutrino Vertex X Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino X Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_x_BDT_high_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_yCoord, yVals, 30, 173.8, 203.8, "yCoordAngleDifferenceBDTCRUMBSProfile_high_equalStats", "Profile of Reco Neutrino Vertex Y Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino Y Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_y_BDT_high_equalStats.pdf").c_str());
    makeEqualStatProfile(xVals_zCoord, yVals, 30, 469.5, 509.5, "zCoordAngleDifferenceBDTCRUMBSProfile_high_equalStats", "Profile of Reco Neutrino Vertex Z Coordinate vs Angle Betwee True and Reco Track: BDT Vertexing;Reco Neutrino Z Coordinate (cm);Angle Difference (degrees)", (base_path + "angleDiffPositionProfile_z_BDT_high_equalStats.pdf").c_str());
}
