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

struct trueInfo{
    int numElectrons = 0;
    int numPhotons = 0;
    int numMuons = 0;
    int numPiZero = 0;
    int numChargedPi = 0;
    int numProtons = 0;
};

struct interactionCounter{
    int Unknown = 0;
    int QE = 0;
    int Res = 0;
    int DIS = 0;
    int Coh = 0;
    int CohElastic = 0;
    int ElectronScattering = 0;
    int IMDAnnihilation = 0;
    int InverseBetaDecay = 0;
    int GlashowResonance = 0;
    int AMNuGamma = 0;
    int MEC = 0;
    int Diffractive = 0;
    int EM = 0;
    int WeakMix = 0;
    int NuanceOffset = 0;
    int CCQE = 0;
    int NCQE = 0;
    int NuanceRes = 0;
    int CCDis = 0;
    int NCDis = 0;
    int NuEElastic = 0;
};

struct counter{
    int current = 0;
    int uboone = 0;
};

typedef struct{
    TCanvas* canvas;
    TH1F* baseHist;
    TH1F* current;
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
    double genieMode;
    double genieInteraction;
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
        (TH1F*) base->Clone((baseName + "_dluboone").c_str()),
    };    
}


void styleDraw(TCanvas* canvas, TH1F* current, TH1F* uboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, TPaveText* pt = nullptr, int* percentage = nullptr, int* drawLine = nullptr, int* linePos = nullptr){
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();

    gStyle->SetPalette(kAvocado);
    gROOT->ForceStyle();
    gPad->Update();

    current->SetLineWidth(2);
    current->SetLineColor(TColor::GetColor("#e42536"));

    uboone->SetLineWidth(2);
    uboone->SetLineColor(TColor::GetColor("#5790fc"));

    current->Draw("hist");
    uboone->Draw("histsame");

    if((ymin != 999) && (ymax != 999)) current->GetYaxis()->SetRangeUser(ymin, ymax);
    
    if((xmin != 999) && (xmax != 999)) current->GetXaxis()->SetRangeUser(xmin, xmax);

    double maxYValue = std::max({current->GetMaximum(), uboone->GetMaximum()});
    double minYValue = std::min({current->GetMinimum(), uboone->GetMinimum()});

    std::cout << minYValue << std::endl;
    if(ymax == 999 && ymin == 999) current->GetYaxis()->SetRangeUser(minYValue*0.95, maxYValue*1.05);
    if(ymax == 999 && ymin != 999) current->GetYaxis()->SetRangeUser(ymin, maxYValue*1.05);
    if(ymax != 999 && ymin == 999) current->GetYaxis()->SetRangeUser(minYValue*0.95, ymax);

    current->SetStats(0);
    current->GetXaxis()->SetTickLength(0.04);
    current->GetYaxis()->SetTickLength(0.03);
    current->GetXaxis()->SetTickSize(0.02);
    current->GetYaxis()->SetTickSize(0.02);

    auto legend = new TLegend(Lxmin,Lymax,Lxmax,Lymin);
    legend->AddEntry(uboone, "Pandora Deep Learning: #muBooNE/BNB Tune", "f");
    legend->AddEntry(current, "Pandora BDT SBND (without Refinement)", "f");
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.13);
    legend->Draw();

    if(drawLine){
        TLine* line = nullptr;
        
        if(ymax != 999 && ymin != 999) line = new TLine(1.022, ymin, 1.022, ymax);
        if(ymax == 999 && ymin == 999) line = new TLine(1.022, minYValue*0.95, 1.022, maxYValue*1.05);
        if(ymax == 999 && ymin != 999) line = new TLine(1.022, ymin, 1.022, maxYValue*1.05);
        if(ymax != 999 && ymin == 999) line = new TLine(1.022, minYValue*0.95, 1.022, ymax);
        
        if(line){
            line->SetLineColor(kGray+2);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");
        }

        TLatex* latex = nullptr;    
        // Labels line on the left
        if(*linePos == 0){
            if(ymax != 999) latex = new TLatex(1.022 - 0.7, ymax * 0.98, "2m_{e}");
            if(ymax == 999) latex = new TLatex(1.022 - 0.7, maxYValue * 0.98, "2m_{e}");
        } else{
            if(ymax != 999) latex = new TLatex(1.022 + 0.7, ymax * 0.98, "2m_{e}");
            if(ymax == 999) latex = new TLatex(1.022 + 0.7, maxYValue * 0.98, "2m_{e}");
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

void percentage(TH1F* current, TH1F* uboone, double sizeCurrent, double sizeUboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr, int* linePos = nullptr){
    TCanvas *percentageCanvas = new TCanvas("percentage_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentPerc = (TH1F*) current->Clone("perc hist");
    currentPerc->Scale(100.0 * 1.0/sizeCurrent);
    currentPerc->GetYaxis()->SetTitle("Percentage of Events (%)"); 

    TH1F* uboonePerc = (TH1F*) uboone->Clone("perc hist");
    uboonePerc->Scale(100.0 * 1.0/sizeUboone);
    uboonePerc->GetYaxis()->SetTitle("Percentage of Events (%)");

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.17, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 
   
    int funcValue = 1;

    styleDraw(percentageCanvas, currentPerc, uboonePerc, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt, &funcValue, drawLine, linePos);
}

void efficiency(TH1F* current, TH1F* uboone, double sizeCurrent, double sizeUboone, double ymin, double ymax, double xmin, double xmax, const char* filename, double Lxmin, double Lxmax, double Lymin, double Lymax, int* drawLine = nullptr, int* linePos = nullptr, std::string xlabel = ""){
    TCanvas *efficiencyCanvas = new TCanvas("efficiency_canvas", "Graph Draw Options", 200, 10, 600, 400); 
    TH1F* currentEff = (TH1F*) current->Clone("eff hist");
    currentEff->Reset();
    currentEff->GetYaxis()->SetTitle("Background Rejection"); 
    currentEff->GetXaxis()->SetTitle(xlabel.c_str());

    TH1F* ubooneEff = (TH1F*) uboone->Clone("eff hist");
    ubooneEff->Reset();
    ubooneEff->GetYaxis()->SetTitle("Background Rejection");
    ubooneEff->GetXaxis()->SetTitle(xlabel.c_str());
    
    int numBins = current->GetNbinsX();
    double currentSum = 0.0;
    double ubooneSum = 0.0;

    for(int i = 1; i <= numBins; ++i){
        currentSum += current->GetBinContent(i);
        ubooneSum += uboone->GetBinContent(i);

        double currentEffValue = currentSum/sizeCurrent;
        double ubooneEffValue = ubooneSum/sizeUboone;

        currentEff->SetBinContent(i, 1-currentEffValue);
        ubooneEff->SetBinContent(i, 1-ubooneEffValue);
    }

    TPaveText* pt = new TPaveText(Lxmin, Lymin - 0.02 - 0.15, Lxmax, Lymin - 0.02, "NDC");
    pt->AddText(Form("Number of DL Uboone Entries: %d", (int)sizeUboone));
    pt->AddText(Form("Number of Current Entries: %d", (int)sizeCurrent));
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(0); 

    int funcValue = 1;

    styleDraw(efficiencyCanvas, currentEff, ubooneEff, ymin, ymax, xmin, xmax, filename, Lxmin, Lxmax, Lymin, Lymax, pt = nullptr, &funcValue, drawLine, linePos);
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
        // Choose slice based on CRUMBS score
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

void obtainTrueParticles(std::vector<trueParticle> trueParticles){
    std::cout << "Number of true particles = " << trueParticles.size() << std::endl;
    trueInfo truth;
    for(size_t j = 0; j < trueParticles.size(); ++j){
        if(trueParticles[j].pdg == 11) truth.numElectrons++;
        if(trueParticles[j].pdg == 22) truth.numPhotons++;
        if(trueParticles[j].pdg == 13) truth.numMuons++;
        if(trueParticles[j].pdg == 111) truth.numPiZero++;
        if(trueParticles[j].pdg == 211 || trueParticles[j].pdg == -211) truth.numChargedPi++;
        if(trueParticles[j].pdg == 2212) truth.numProtons++;
    }

    printf("Number of Electrons = %i, Number of Photons = %i, Number of Muons = %i, Number of PiZero = %i, Number of Charged Pi = %i, Number of Protons = %i\n", truth.numElectrons, truth.numPhotons, truth.numMuons, truth.numPiZero, truth.numChargedPi, truth.numProtons); 
}

void nuEBackgroundCosmics_macro(){
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/IntimeCosmics/merged.root");

    std::string base_path = "/nashome/c/coackley/nuEBackgroundPlotsCosmics/";

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
    double DLCurrent, signal;
    
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("runID", &runID);
    tree->SetBranchAddress("subRunID", &subRunID);
    tree->SetBranchAddress("DLCurrent", &DLCurrent);
    tree->SetBranchAddress("signal", &signal);

    std::vector<double> *truth_neutrinoVX = nullptr;
    std::vector<double> *truth_neutrinoVY = nullptr;
    std::vector<double> *truth_neutrinoVZ = nullptr;
    std::vector<double> *truth_CCNC = nullptr;
    std::vector<double> *truth_neutrinoType = nullptr;
    std::vector<double> *truth_chargedLepton = nullptr;
    std::vector<double> *truth_neutrinoTPCID = nullptr;
    std::vector<double> *truth_neutrinoTPCValid = nullptr;
    std::vector<double> *truth_genieMode = nullptr;
    std::vector<double> *truth_genieInteraction = nullptr;

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
    tree->SetBranchAddress("truth_genieMode", &truth_genieMode);
    tree->SetBranchAddress("truth_genieInteraction", &truth_genieInteraction);

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
    //auto numSlicesCompleteness = createHistGroup("numSlicesCompleteness", "Number of Slices with Completeness > 0 in an Event", "Number of Slices with Completeness > 0", 10, 0, 10);
    auto numRecoNeutrinos = createHistGroup("numRecoNeutrinos", "Number of Reco Neutrinos in an Event", "Number of Reco Neutrinos", 10, 0, 10);

    // Highest CRUMBS Score Slice Plots
    //auto sliceCompletenessCRUMBS = createHistGroup("sliceCompletenessCRUMBS", "Completeness of the Slice with the Highest CRUMBS Score", "Completeness", 102, 0, 1.02); // Bin width = 0.005
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
    //auto angleDifferenceCRUMBS = createHistGroup("angleDifferenceCRUMBS", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    auto EtrueThetaRecoCRUMBS = createHistGroup("EtrueThetaRecoCRUMBS", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    auto ERecoSumThetaTrueCRUMBS = createHistGroup("ERecoSumThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBS = createHistGroup("ERecoHighestThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCRUMBS = createHistGroup("ERecoSumThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    auto ERecoHighestThetaRecoCRUMBS = createHistGroup("ERecoHighestThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);

    // Highest Completeness Score Slice Plots
    //auto sliceCompletenessCompleteness = createHistGroup("sliceCompletenessCompleteness", "Completeness of the Slice with the Highest Completeness", "Completeness Score", 102, 0, 1.02); // Bin width = 0.01
    //auto sliceScoreCompleteness = createHistGroup("sliceScoreCompleteness", "CRUMBS Score of the Slice with the Highest Completeness", "CRUMBS Score", 25, -1, 1); 
    //auto slicePurityCompleteness = createHistGroup("slicePurityCompleteness", "Purity of the Slice with the Highest Completeness", "Purity", 50, 0, 1.02);
    //auto deltaXCompleteness = createHistGroup("deltaXCompleteness", "#Deltax Distribution: Slice with Highest Completeness", "x_{Reco} - x_{True} (cm)", 40, -5, 5);
    //auto deltaYCompleteness = createHistGroup("deltaYCompleteness", "#Deltay Distribution: Slice with Highest Completeness", "y_{Reco} - y_{True} (cm)", 40, -5, 5);
    //auto deltaZCompleteness = createHistGroup("deltaZCompleteness", "#Deltaz Distribution: Slice with Highest Completeness", "z_{Reco} - z_{True} (cm)", 40, -5, 5);
    //auto deltaRCompleteness = createHistGroup("deltaRCompleteness", "#Delta#bar{r} Distribution: Slice with Highest Completeness Score", "|#bar{r}_{Reco} - #bar{r}_{True}| (cm)", 20, 0, 5);
    //auto numPFPsCompleteness = createHistGroup("numPFPsCompleteness", "Number of PFPs in the Slice with the Highest Completeness", "Number of PFPs", 10, 0, 10);
    //auto ratioChosenSummedEnergyCompleteness = createHistGroup("ratioChosenSummedEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest Completeness", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    //auto ratioChosenTrueEnergyCompleteness = createHistGroup("ratioChosenTrueEnergyCompleteness", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    //auto ratioSummedTrueEnergyCompleteness = createHistGroup("ratioSummedTrueEnergyCompleteness", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest Completeness and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
    //auto angleDifferenceCompleteness = createHistGroup("angleDifferenceCompleteness", "Angle Between the True Electron and the Highest Energy PFP in the Slice with the Highest Completeness", "arccos#left(|#vec{A} #upoint #vec{B}|/|#vec{A}||#vec{B}|#right) (degrees)", 93, 0, 186);
    //auto EtrueThetaRecoCompleteness = createHistGroup("EtrueThetaRecoCompleteness", "E_{true}#theta_{reco}^{2}: Slice with Highest Completeness", "E_{true}#theta_{reco}^{2} (MeV)", 40, 0, 20.44);
    //auto ERecoSumThetaTrueCompleteness = createHistGroup("ERecoSumThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    //auto ERecoHighestThetaTrueCompleteness = createHistGroup("ERecoHighestThetaTrueCompleteness", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    //auto ERecoSumThetaRecoCompleteness = createHistGroup("ERecoSumThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    //auto ERecoHighestThetaRecoCompleteness = createHistGroup("ERecoHighestThetaRecoCompleteness", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest Completeness", "E_{reco}#theta_{reco}^{2} (MeV)", 27, 0, 13.797);
    
    counter numEventsTotal;
    counter numEventsTrueNeutrino;
    counter numEventsSlices;
    counter numEventsRecoNeutrino;
    counter numEventsSliceNotEqualNeutrino;
    //counter numSlicesCRUMBSCompletenessZero;
    counter numEventsCRUMBSRecoParticle;
    //counter sameSliceSelected;
    interactionCounter interactions;

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);

        printf("___________________________________________________________________\n");
        std::vector<recoParticle> recoParticlesInEvent;
        std::vector<recoNeutrino> recoNeutrinosInEvent;
        std::vector<recoSlice> recoSlicesInEvent;
        std::vector<trueParticle> trueParticlesInEvent;

        recoParticle chosenRecoParticleCRUMBS;
        //recoParticle chosenRecoParticleCompleteness;
        recoNeutrino chosenRecoNeutrinoCRUMBS;
        //recoNeutrino chosenRecoNeutrinoCompleteness;
        //recoSlice chosenRecoSliceCompleteness;
        recoSlice chosenRecoSliceCRUMBS;
        //trueParticle chosenTrueParticle;
        trueNeutrino chosenTrueNeutrino;
        double energyInSlice = 0;
        double numPFPsInEvent = 0;
        double numPFPsInChosenSlice = 0;

        if(DLCurrent == 0){
            numEventsTotal.uboone++;
        } else if(DLCurrent == 2){
            numEventsTotal.current++;
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
                chosenTrueNeutrino.genieMode = truth_genieMode->at(j);
                chosenTrueNeutrino.genieInteraction = truth_genieInteraction->at(j);

                if(DLCurrent == 0){
                    if(chosenTrueNeutrino.genieMode == -1) interactions.Unknown++;
                    if(chosenTrueNeutrino.genieMode == 0) interactions.QE++;
                    if(chosenTrueNeutrino.genieMode == 1) interactions.Res++;
                    if(chosenTrueNeutrino.genieMode == 2) interactions.DIS++;
                    if(chosenTrueNeutrino.genieMode == 3) interactions.Coh++;
                    if(chosenTrueNeutrino.genieMode == 4) interactions.CohElastic++;
                    if(chosenTrueNeutrino.genieMode == 5) interactions.ElectronScattering++;
                    if(chosenTrueNeutrino.genieMode == 6) interactions.IMDAnnihilation++;
                    if(chosenTrueNeutrino.genieMode == 7) interactions.InverseBetaDecay++;
                    if(chosenTrueNeutrino.genieMode == 8) interactions.GlashowResonance++;
                    if(chosenTrueNeutrino.genieMode == 9) interactions.AMNuGamma++;
                    if(chosenTrueNeutrino.genieMode == 10) interactions.MEC++;
                    if(chosenTrueNeutrino.genieMode == 11) interactions.Diffractive++;
                    if(chosenTrueNeutrino.genieMode == 12) interactions.EM++;
                    if(chosenTrueNeutrino.genieMode == 13) interactions.WeakMix++;
                    if(chosenTrueNeutrino.genieInteraction == 1000) interactions.NuanceOffset++;
                    if(chosenTrueNeutrino.genieInteraction == 1001) interactions.CCQE++;
                    if(chosenTrueNeutrino.genieInteraction == 1002) interactions.NCQE++;
                    if(chosenTrueNeutrino.genieInteraction > 1002 && chosenTrueNeutrino.genieInteraction < 1091) interactions.NuanceRes++;
                    if(chosenTrueNeutrino.genieInteraction == 1091) interactions.CCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1092) interactions.NCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1098) interactions.NuEElastic++;
                } 
            } else if(truth_neutrinoTPCValid->at(j) != -999999 && truth_neutrinoVX->at(j) < 201.3 && truth_neutrinoVX->at(j) > -201.3 && truth_neutrinoVY->at(j) < 203.732 && truth_neutrinoVY->at(j) > -203.732 && truth_neutrinoVZ->at(j) > 0 && truth_neutrinoVZ->at(j) < 509.4){
                //printf("True Neutrino Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            } else if(truth_neutrinoTPCValid->at(j) != -999999){
                //printf("True Neutrino Not Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            }
        }
        
        if(truth_neutrinoVX->size() == 1 && truth_neutrinoVX->at(0) == -999999) printf("No True Neutrino in Event\n");

        // If there are no truth neutrinos within the TPC, skip the event
        //if(truthNeutrinoInTPC == 0) continue;

        if(DLCurrent == 0){
            numEventsTrueNeutrino.uboone++;
        } else if(DLCurrent == 2){
            numEventsTrueNeutrino.current++;
        }

        printf("True Neutrino Vertex = (%f, %f, %f)\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz);
        std::cout << "Number of truth particles = " << truth_particlePDG->size() << std::endl;             
    
        for(size_t k = 0; k < truth_particlePDG->size(); ++k){
            trueParticle truthParticle;
            truthParticle.pdg = truth_particlePDG->at(k);
            truthParticle.vx = truth_particleVX->at(k);
            truthParticle.vy = truth_particleVY->at(k);
            truthParticle.vz = truth_particleVZ->at(k);
            truthParticle.px = truth_particlePX->at(k);
            truthParticle.py = truth_particlePY->at(k);
            truthParticle.pz = truth_particlePZ->at(k);
            truthParticle.energy = truth_particleEnergy->at(k);
            truthParticle.angle = truth_particleAngle->at(k);
            truthParticle.ETheta2 = truth_particleETheta2->at(k);
            truthParticle.dx = truth_particleDirectionX->at(k);
            truthParticle.dy = truth_particleDirectionY->at(k);
            truthParticle.dz = truth_particleDirectionZ->at(k);
            truthParticle.neutrinoParent = truth_particleNeutrinoParent->at(k);
            trueParticlesInEvent.push_back(truthParticle);
        }

        obtainTrueParticles(trueParticlesInEvent);

        int slice = 0;
        int crumbsSlice = 0;
        //int completenessSlice = 0;
        double numSlicesInEvent = 0;
        for(size_t j = 0; j < reco_sliceID->size(); ++j){
            if(reco_sliceID->at(j) != -999999){
                slice = 1;
                numSlicesInEvent++;

                recoSlice recoSlice;
                recoSlice.id = reco_sliceID->at(j);
                //recoSlice.completeness = reco_sliceCompleteness->at(j);
                //recoSlice.purity = reco_slicePurity->at(j);
                recoSlice.completeness = -999999;
                recoSlice.purity = -999999;
                recoSlice.score = reco_sliceScore->at(j);
                recoSlicesInEvent.push_back(recoSlice);
        
                if(recoSlice.score != -999999) crumbsSlice++;  
                //if(recoSlice.completeness > 0) completenessSlice++;  
                printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", j, recoSlice.id, recoSlice.completeness, recoSlice.purity, recoSlice.score);
            }
        }    

        // If there are no slices, skip the event
        if(slice == 0) continue;  
     
        // Pick the slice with the highest completeness
        //chosenRecoSliceCompleteness = chooseSlice(recoSlicesInEvent, 0);
        //printf("\nSlice with Highest Completeness: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", chosenRecoSliceCompleteness.id, chosenRecoSliceCompleteness.completeness, chosenRecoSliceCompleteness.purity, chosenRecoSliceCompleteness.score);

        // Pick the slice with the highest score
        chosenRecoSliceCRUMBS = chooseSlice(recoSlicesInEvent, 1);
        printf("Slice with Highest CRUMBS Score: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n\n", chosenRecoSliceCRUMBS.id, chosenRecoSliceCRUMBS.completeness, chosenRecoSliceCRUMBS.purity, chosenRecoSliceCRUMBS.score);
        
        if(DLCurrent == 0){
            numEventsSlices.uboone++;
            numSlices.uboone->Fill(numSlicesInEvent);
            numSlicesCRUMBS.uboone->Fill(crumbsSlice);
            //numSlicesCompleteness.uboone->Fill(completenessSlice);
            
            //sliceCompletenessCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.score);
            //slicePurityCRUMBS.uboone->Fill(chosenRecoSliceCRUMBS.purity);
            
            //sliceCompletenessCompleteness.uboone->Fill(chosenRecoSliceCompleteness.completeness);
            //sliceScoreCompleteness.uboone->Fill(chosenRecoSliceCompleteness.score);
            //slicePurityCompleteness.uboone->Fill(chosenRecoSliceCompleteness.purity);

            //if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.uboone++;
            //if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.uboone++;
        } else if(DLCurrent == 2){
            numEventsSlices.current++;
            numSlices.current->Fill(numSlicesInEvent);
            numSlicesCRUMBS.current->Fill(crumbsSlice);
            //numSlicesCompleteness.current->Fill(completenessSlice);
            
            //sliceCompletenessCRUMBS.current->Fill(chosenRecoSliceCRUMBS.completeness);
            sliceScoreCRUMBS.current->Fill(chosenRecoSliceCRUMBS.score);
            //slicePurityCRUMBS.current->Fill(chosenRecoSliceCRUMBS.purity);
            
            //sliceCompletenessCompleteness.current->Fill(chosenRecoSliceCompleteness.completeness);
            //sliceScoreCompleteness.current->Fill(chosenRecoSliceCompleteness.score);
            //slicePurityCompleteness.current->Fill(chosenRecoSliceCompleteness.purity);
                
            //if(chosenRecoSliceCRUMBS.completeness == 0) numSlicesCRUMBSCompletenessZero.current++;
            //if(chosenRecoSliceCRUMBS.id == chosenRecoSliceCompleteness.id) sameSliceSelected.current++;
        }
        
        int neutrino = 0;
        int numRecoNeutrinosInEvent = 0;
        std::cout << "num reco neutrinos: " << reco_neutrinoPDG->size() << std::endl;
        for(size_t j = 0; j < reco_neutrinoPDG->size(); ++j){
            if(reco_neutrinoPDG->at(j) != -999999){
                neutrino = 1;
                numRecoNeutrinosInEvent++;
                
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
        if(DLCurrent == 2) numRecoNeutrinos.current->Fill(numRecoNeutrinosInEvent);

        // Skip the event if there are no reconstructed neutrinos
        if(neutrino == 0) continue;
    
        // The number of slices with a CRUMBS score != the number of reconstructed neutrinos in the event 
        if(crumbsSlice != numRecoNeutrinosInEvent){
            if(DLCurrent == 0) numEventsSliceNotEqualNeutrino.uboone++; 
            if(DLCurrent == 2) numEventsSliceNotEqualNeutrino.current++;
        } 

        // Look at reconstruction when we pick the slice with highest CRUMBS score
        chosenRecoNeutrinoCRUMBS = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCRUMBS.id); 

        double deltaXCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vx - chosenTrueNeutrino.vx);
        double deltaYCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vy - chosenTrueNeutrino.vy);
        double deltaZCRUMBSValue = (chosenRecoNeutrinoCRUMBS.vz - chosenTrueNeutrino.vz);
        double deltaRCRUMBSValue = std::sqrt((deltaXCRUMBSValue * deltaXCRUMBSValue) + (deltaYCRUMBSValue * deltaYCRUMBSValue) + (deltaZCRUMBSValue * deltaZCRUMBSValue));

        if(DLCurrent == 0){
            numEventsRecoNeutrino.uboone++;
            deltaXCRUMBS.uboone->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.uboone->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.uboone->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.uboone->Fill(deltaRCRUMBSValue);
        } else if(DLCurrent == 2){
            numEventsRecoNeutrino.current++;
            deltaXCRUMBS.current->Fill(deltaXCRUMBSValue);
            deltaYCRUMBS.current->Fill(deltaYCRUMBSValue);
            deltaZCRUMBS.current->Fill(deltaZCRUMBSValue);
            deltaRCRUMBS.current->Fill(deltaRCRUMBSValue);
        }

        int recoparticle = 0;
        int recoparticleCRUMBS = 0;

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
                //recoParticle.completeness = reco_particleCompleteness->at(j);
                //recoParticle.purity = reco_particlePurity->at(j);
                recoParticle.completeness = -999999;
                recoParticle.purity = -999999;
                recoParticlesInEvent.push_back(recoParticle);
            
                if(recoParticle.sliceID == chosenRecoSliceCRUMBS.id) recoparticleCRUMBS = 1;
                //if(recoParticle.sliceID == chosenRecoSliceCompleteness.id) recoparticleCompleteness = 1;
            }
        }

        // There are no reco particles in the event
        if(recoparticle == 0) continue;
        
        // No reco particles in the slice with the highest CRUMBS score
        if(recoparticleCRUMBS == 0) continue;
        
        double totalSliceEnergyCRUMBS = 0;
        double numPFPsSliceCRUMBS = 0;
        chosenRecoParticleCRUMBS = choosePFP(recoParticlesInEvent, chosenRecoSliceCRUMBS.id, totalSliceEnergyCRUMBS, numPFPsSliceCRUMBS);

        //double aDOTbCRUMBS = ((chosenRecoParticleCRUMBS.dx * chosenTrueParticle.dx) + (chosenRecoParticleCRUMBS.dy * chosenTrueParticle.dy) + (chosenRecoParticleCRUMBS.dz * chosenTrueParticle.dz));
        //double aMagCRUMBS = std::sqrt((chosenRecoParticleCRUMBS.dx * chosenRecoParticleCRUMBS.dx) + (chosenRecoParticleCRUMBS.dy * chosenRecoParticleCRUMBS.dy) + (chosenRecoParticleCRUMBS.dz * chosenRecoParticleCRUMBS.dz));
        //double bMagCRUMBS = std::sqrt((chosenTrueParticle.dx * chosenTrueParticle.dx) + (chosenTrueParticle.dy * chosenTrueParticle.dy) + (chosenTrueParticle.dz * chosenTrueParticle.dz));
        //double cosAngleCRUMBS = aDOTbCRUMBS/(aMagCRUMBS * bMagCRUMBS);
        //if(cosAngleCRUMBS > 1.0) cosAngleCRUMBS = 1.0;
        //if(cosAngleCRUMBS < -1.0) cosAngleCRUMBS = -1.0;
        //double angleDiffCRUMBS = TMath::ACos(cosAngleCRUMBS) * TMath::RadToDeg();
        
        if(DLCurrent == 0){
            numEventsCRUMBSRecoParticle.uboone++;
            numPFPsCRUMBS.uboone->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            //ratioChosenTrueEnergyCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            //ratioSummedTrueEnergyCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            //angleDifferenceCRUMBS.uboone->Fill(angleDiffCRUMBS);
            //EtrueThetaRecoCRUMBS.uboone->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            //ERecoSumThetaTrueCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            //ERecoHighestThetaTrueCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.uboone->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.uboone->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);

        } else if(DLCurrent == 2){
            numEventsCRUMBSRecoParticle.current++;
            numPFPsCRUMBS.current->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            //ratioChosenTrueEnergyCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / chosenTrueParticle.energy);
            //ratioSummedTrueEnergyCRUMBS.current->Fill(totalSliceEnergyCRUMBS / chosenTrueParticle.energy);
            //angleDifferenceCRUMBS.current->Fill(angleDiffCRUMBS);
            //EtrueThetaRecoCRUMBS.current->Fill(chosenTrueParticle.energy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            //ERecoSumThetaTrueCRUMBS.current->Fill(totalSliceEnergyCRUMBS * chosenTrueParticle.angle * chosenTrueParticle.angle);
            //ERecoHighestThetaTrueCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenTrueParticle.angle * chosenTrueParticle.angle);
            ERecoSumThetaRecoCRUMBS.current->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.current->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            
        }

        printf("Chosen True Neutrino: Vertex = (%f, %f, %f), CCNC = %f, PDG = %f, Lepton PDG = %f, TPC ID = %f, TPC Valid = %f\n", chosenTrueNeutrino.vx, chosenTrueNeutrino.vy, chosenTrueNeutrino.vz, chosenTrueNeutrino.CCNC, chosenTrueNeutrino.pdg, chosenTrueNeutrino.leptonpdg, chosenTrueNeutrino.tpcID, chosenTrueNeutrino.tpcValid);    
        printf("Number of Slices in event: %f\n", numSlicesInEvent);
    }

    printf("Number of events with a true neutrino within the TPC:\nUBoone: %i out of %i\nCurrent: %i out of %i\n", numEventsTrueNeutrino.uboone, numEventsTotal.uboone, numEventsTrueNeutrino.current, numEventsTotal.current);
    printf("Number of events with a slice:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSlices.uboone, numEventsTotal.uboone, numEventsSlices.current, numEventsTotal.current);
    printf("Number of events with a reco neutrino:\nUboone:%i out of %i\nCurrent: %i out of %i\n", numEventsRecoNeutrino.uboone, numEventsTotal.uboone, numEventsRecoNeutrino.current, numEventsTotal.current);
    printf("Number of events where the number of slices with a CRUMBS score != number of reco neutrinos:\nUboone: %i out of %i\nCurrent: %i out of %i\n", numEventsSliceNotEqualNeutrino.uboone, numEventsRecoNeutrino.uboone, numEventsSliceNotEqualNeutrino.current, numEventsRecoNeutrino.current);

    styleDraw(numSlices.canvas, numSlices.current, numSlices.uboone, 999, 999, 999, 999, (base_path + "numSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlices.current, numSlices.uboone, numEventsSlices.current, numEventsSlices.uboone, 999, 999, 999, 999, (base_path + "numSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
 
    styleDraw(numSlicesCRUMBS.canvas, numSlicesCRUMBS.current, numSlicesCRUMBS.uboone, 999, 999, 999, 999, (base_path + "numCRUMBSSlices_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numSlicesCRUMBS.current, numSlicesCRUMBS.uboone, numEventsSlices.current, numEventsSlices.uboone, 999, 999, 999, 999, (base_path + "numCRUMBSSlices_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    
    styleDraw(numRecoNeutrinos.canvas, numRecoNeutrinos.current, numRecoNeutrinos.uboone, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numRecoNeutrinos.current, numRecoNeutrinos.uboone, numEventsSlices.current, numEventsSlices.uboone, 999, 999, 999, 999, (base_path + "numRecoNeutrinos_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // CRUMBS Slice Score Plots
    styleDraw(sliceScoreCRUMBS.canvas, sliceScoreCRUMBS.current, sliceScoreCRUMBS.uboone, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(sliceScoreCRUMBS.current, sliceScoreCRUMBS.uboone, numEventsSlices.current, numEventsSlices.uboone, 999, 999, 999, 999, (base_path + "sliceScoreCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    
    // CRUMBS Slice Vertex Plots
    styleDraw(deltaXCRUMBS.canvas, deltaXCRUMBS.current, deltaXCRUMBS.uboone, 999, 999, 999, 999, (base_path + "deltaXCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaYCRUMBS.canvas, deltaYCRUMBS.current, deltaYCRUMBS.uboone, 999, 999, 999, 999, (base_path + "deltaYCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaZCRUMBS.canvas, deltaZCRUMBS.current, deltaZCRUMBS.uboone, 999, 999, 999, 999, (base_path + "deltaZCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    styleDraw(deltaRCRUMBS.canvas, deltaRCRUMBS.current, deltaRCRUMBS.uboone, 999, 999, 999, 999, (base_path + "deltaRCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaXCRUMBS.current, deltaXCRUMBS.uboone, numEventsRecoNeutrino.current, numEventsRecoNeutrino.uboone, 999, 999, 999, 999, (base_path + "deltaXCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaYCRUMBS.current, deltaYCRUMBS.uboone, numEventsRecoNeutrino.current, numEventsRecoNeutrino.uboone, 999, 999, 999, 999, (base_path + "deltaYCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaZCRUMBS.current, deltaZCRUMBS.uboone, numEventsRecoNeutrino.current, numEventsRecoNeutrino.uboone, 999, 999, 999, 999, (base_path + "deltaZCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(deltaRCRUMBS.current, deltaRCRUMBS.uboone, numEventsRecoNeutrino.current, numEventsRecoNeutrino.uboone, 999, 999, 999, 999, (base_path + "deltaRCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    // CRUMBS PFP Plots
    styleDraw(numPFPsCRUMBS.canvas, numPFPsCRUMBS.current, numPFPsCRUMBS.uboone, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    percentage(numPFPsCRUMBS.current, numPFPsCRUMBS.uboone, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.uboone, 999, 999, 999, 999, (base_path + "numPFPsCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(ratioChosenSummedEnergyCRUMBS.canvas, ratioChosenSummedEnergyCRUMBS.current, ratioChosenSummedEnergyCRUMBS.uboone, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_dist.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);
    percentage(ratioChosenSummedEnergyCRUMBS.current, ratioChosenSummedEnergyCRUMBS.uboone, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.uboone, 999, 999, 999, 999, (base_path + "ratioChosenSummedEnergyCRUMBS_perc.pdf").c_str(), 1-0.86, 1-0.54, 0.70, 0.86);

    int drawLine = 1;
    int left = 0;
    int right = 1;
    
    //styleDraw(angleDifferenceCRUMBS.canvas, angleDifferenceCRUMBS.current, angleDifferenceCRUMBS.cheated, angleDifferenceCRUMBS.dune, angleDifferenceCRUMBS.uboone, angleDifferenceCRUMBS.sbnd, 0, 260, 999, 999, (base_path + "angleDifferenceCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);
    //percentage(angleDifferenceCRUMBS.current, angleDifferenceCRUMBS.cheated, angleDifferenceCRUMBS.dune, angleDifferenceCRUMBS.uboone, angleDifferenceCRUMBS.sbnd, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.cheated, numEventsCRUMBSRecoParticle.dune, numEventsCRUMBSRecoParticle.uboone, numEventsCRUMBSRecoParticle.sbnd, 0, 25, 999, 999, (base_path + "angleDifferenceCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86);

    styleDraw(ERecoSumThetaRecoCRUMBS.canvas, ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.uboone, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.uboone, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.uboone, 999, 999, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoSumThetaRecoCRUMBS.current, ERecoSumThetaRecoCRUMBS.uboone, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.uboone, 0.9, 1, 999, 999, (base_path + "ERecoSumThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    styleDraw(ERecoHighestThetaRecoCRUMBS.canvas, ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.uboone, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_dist.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, nullptr, nullptr, &drawLine, &right);
    percentage(ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.uboone, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.uboone, 999, 999, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_perc.pdf").c_str(), 0.56, 0.88, 0.7, 0.86, &drawLine, &right);
    efficiency(ERecoHighestThetaRecoCRUMBS.current, ERecoHighestThetaRecoCRUMBS.uboone, numEventsCRUMBSRecoParticle.current, numEventsCRUMBSRecoParticle.uboone, 0.9, 1, 999, 999, (base_path + "ERecoHighestThetaRecoCRUMBS_eff.pdf").c_str(), 0.56, 0.88, 0.14, 0.3, &drawLine, &left, "E_{reco}#theta_{reco}^{2} (MeV)");

    printf("Interactions in Sample:\nUnknown = %i\nQE = %i\nRes = %i\nDIS = %i\nCoh = %i\nCoh Elastic = %i\nElastic Scattering = %i\nIMD Annihilation = %i\nInverse Beta Decay = %i\nGlashow Resonance = %i\nAM Nu Gamma = %i\nMEC = %i\nDiffractive = %i\nEM = %i\nWeak Mix = %i\n\n", interactions.Unknown, interactions.QE, interactions.Res, interactions.DIS, interactions.Coh, interactions.CohElastic, interactions.ElectronScattering, interactions.IMDAnnihilation, interactions.InverseBetaDecay, interactions.GlashowResonance, interactions.AMNuGamma, interactions.MEC, interactions.Diffractive, interactions.EM, interactions.WeakMix);
    printf("Nuance Offsets:\nNumber with just Offset = %i\nCCQE = %i\nNCQE = %i\nNuanceRes = %i\nCCDIS = %i\nNCDis = %i\nNuEElastic = %i\n", interactions.NuanceOffset, interactions.CCQE, interactions.NCQE, interactions.NuanceRes, interactions.CCDis, interactions.NCDis, interactions.NuEElastic);
}
