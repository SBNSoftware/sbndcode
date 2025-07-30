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
    int electron = 0;
    int positron = 0;
    int muon = 0;
    int antimuon = 0;
    int chargedpion = 0;
    int negativepion = 0;
    int photon = 0;
    int proton = 0; 
};

typedef struct{
    TCanvas* canvas;
    TH1F* Unknown;
    TH1F* CCQE;
    TH1F* NCQE;
    TH1F* Res;
    TH1F* CCDIS;
    TH1F* NCDIS;
    TH1F* Coh;
    TH1F* NuEElastic;
    TH1F* MEC;
    TH1F* other;
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
        (TH1F*) base->Clone((baseName + "_Unknown").c_str()),
        (TH1F*) base->Clone((baseName + "_CCQE").c_str()),
        (TH1F*) base->Clone((baseName + "_NCQE").c_str()),
        (TH1F*) base->Clone((baseName + "_Res").c_str()),
        (TH1F*) base->Clone((baseName + "_CCDIS").c_str()),
        (TH1F*) base->Clone((baseName + "_NCDIS").c_str()),
        (TH1F*) base->Clone((baseName + "_Coh").c_str()),
        (TH1F*) base->Clone((baseName + "_NuEElastic").c_str()),
        (TH1F*) base->Clone((baseName + "_MEC").c_str()),
        (TH1F*) base->Clone((baseName + "_other").c_str())
    };    
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

void styleDraw(histGroup hist, const char* filename){
    hist.canvas->cd();
    hist.canvas->SetTickx();
    hist.canvas->SetTicky();

    gPad->Update();
    hist.Unknown->SetLineWidth(2);
    hist.Unknown->SetLineColor(TColor::GetColor("#92dadd"));
    hist.CCQE->SetLineWidth(2);
    hist.CCQE->SetLineColor(TColor::GetColor("#717581"));
    hist.NCQE->SetLineWidth(2);
    hist.NCQE->SetLineColor(TColor::GetColor("#b9ac70"));
    hist.Res->SetLineWidth(2);
    hist.Res->SetLineColor(TColor::GetColor("#e76300"));
    hist.CCDIS->SetLineWidth(2);
    hist.CCDIS->SetLineColor(TColor::GetColor("#832db6"));
    hist.NCDIS->SetLineWidth(2);
    hist.NCDIS->SetLineColor(TColor::GetColor("#bd1f01"));
    hist.Coh->SetLineWidth(2);
    hist.Coh->SetLineColor(TColor::GetColor("#ffa90e"));
    hist.NuEElastic->SetLineWidth(2);
    hist.NuEElastic->SetLineColor(TColor::GetColor("#3f90da"));
    hist.MEC->SetLineWidth(2);
    hist.MEC->SetLineColor(TColor::GetColor("#a96b59"));
    hist.other->SetLineWidth(2);
    hist.other->SetLineColor(TColor::GetColor("#94a4a2"));

    const char* histTitle = hist.Unknown->GetTitle();
    const char* xAxisTitle = hist.Unknown->GetXaxis()->GetTitle();
    const char* yAxisTitle = hist.Unknown->GetYaxis()->GetTitle();

    std::string stackTitle = std::string(histTitle) + ";" + xAxisTitle + ";" + yAxisTitle;

    THStack* stack = new THStack("stack", stackTitle.c_str());
    stack->Add(hist.Unknown);
    stack->Add(hist.CCQE);
    stack->Add(hist.NCQE);
    stack->Add(hist.Res);
    stack->Add(hist.CCDIS);
    stack->Add(hist.NCDIS);
    stack->Add(hist.Coh);
    stack->Add(hist.NuEElastic);
    stack->Add(hist.MEC);
    stack->Add(hist.other);

    stack->Draw("hist");

    TLegend* legend = new TLegend(0.4, 0.8, 0.9, 0.9);
    legend->SetHeader("GENIE Interaction", "C");
    legend->SetNColumns(5);
    legend->AddEntry(hist.Unknown, "Unknown", "f");
    legend->AddEntry(hist.CCQE, "CCQE", "f");
    legend->AddEntry(hist.NCQE, "NCQE", "f");
    legend->AddEntry(hist.Res, "Res", "f");
    legend->AddEntry(hist.CCDIS, "CCDIS", "f");
    legend->AddEntry(hist.NCDIS, "NCDIS", "f");
    legend->AddEntry(hist.Coh, "Coh", "f");
    legend->AddEntry(hist.NuEElastic, "NuEElastic", "f");
    legend->AddEntry(hist.MEC, "MEC", "f");
    legend->AddEntry(hist.other, "Other", "f");
    legend->SetEntrySeparation(0.025);
    legend->SetTextSize(0.0225);
    legend->SetMargin(0.2);
    legend->Draw();

    hist.canvas->SaveAs(filename);

    delete stack;
}

void nuEStackGenieInt_macro(){
    TFile *file = TFile::Open("/exp/sbnd/data/users/coackley/merged_23July.root");
    //std::string base_path = "/nashome/c/coackley/nuEPlotsGenieIntStack/Signal_Uboone_";
    //std::string base_path = "/nashome/c/coackley/nuEPlotsGenieIntStack/BNB_Uboone_";
    std::string base_path = "/nashome/c/coackley/nuEPlotsGenieIntStack/IntimeCosmic_Uboone_";
    //std::string base_path = "/nashome/c/coackley/nuEPlotsGenieIntStack/Signal_Current_";
    //std::string base_path = "/nashome/c/coackley/nuEPlotsGenieIntStack/BNB_Current_";
    //std::string base_path = "/nashome/c/coackley/nuEPlotsGenieIntStack/IntimeCosmic_Current_";

    // Change based on what events you want to see
    // DL: 0 = DL Vertexing UBoone, 2 = BDT Vertexing
    // Signal: 1 = Nu+E Elastic, 2 = BNB, 3 = Intime Cosmics
    double DLCurrentReq = 0;
    double SignalReq = 3;
    
    if(!file){
        std::cerr << "Error opening the file" << std::endl;
        return;
    }

    TDirectory *dir = (TDirectory*)file->Get("ana");
    if(!dir){
        std::cerr << "Directory 'ana' not found" << std::endl;
        return;
    }

    TTree *subRunTree = (TTree*)dir->Get("SubRun");
    if(!subRunTree){
        std::cerr << "SubRun not found" << std::endl;
        return;
    }

    TTree *tree = (TTree*)dir->Get("NuE");
    if(!tree){
        std::cerr << "NuE not found" << std::endl;
        return;
    }

    double subRunSignal, subRunDLCurrent, subRunPOT;
    int subRunSpills, subRunNumGenEvents;
    unsigned int subRunNumber, subRunRun;

    subRunTree->SetBranchAddress("signal", &subRunSignal);
    subRunTree->SetBranchAddress("DLCurrent", &subRunDLCurrent);
    subRunTree->SetBranchAddress("pot", &subRunPOT);
    subRunTree->SetBranchAddress("spills", &subRunSpills);
    subRunTree->SetBranchAddress("numGenEvents", &subRunNumGenEvents);
    subRunTree->SetBranchAddress("subRun", &subRunNumber);
    subRunTree->SetBranchAddress("run", &subRunRun);

    Long64_t numEntriesSubRun = subRunTree->GetEntries();

    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsSignal;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsBNB;
    std::set<std::pair<unsigned int, unsigned int>> seenSubRunsCosmics;

    double totalPOTSignal = 0;
    double totalPOTBNB = 0;
    double totalPOTCosmics = 0;

    double cosmicSpillsSum = 0;

    for(Long64_t i = 0; i < numEntriesSubRun; ++i){
        subRunTree->GetEntry(i);

        if(subRunSignal == 3 && subRunDLCurrent == 0) cosmicSpillsSum += subRunNumGenEvents;

        if(subRunSignal == 1){
            std::pair<unsigned int, unsigned int> keySignal = std::make_pair(subRunRun, subRunNumber);
            if(subRunDLCurrent == 0 && seenSubRunsSignal.find(keySignal) == seenSubRunsSignal.end()){
                totalPOTSignal += subRunPOT;
                seenSubRunsSignal.insert(keySignal);
            }
        } else if(subRunSignal == 2){
            std::pair<unsigned int, unsigned int> keyBNB = std::make_pair(subRunRun, subRunNumber);
            if(subRunDLCurrent == 0 && seenSubRunsBNB.find(keyBNB) == seenSubRunsBNB.end()){
                totalPOTBNB += subRunPOT;
                seenSubRunsBNB.insert(keyBNB);
            }
        } else if(subRunSignal == 3){
            std::pair<unsigned int, unsigned int> keyCosmics = std::make_pair(subRunRun, subRunNumber);
            if(subRunDLCurrent == 0 && seenSubRunsCosmics.find(keyCosmics) == seenSubRunsCosmics.end()){
                totalPOTCosmics += subRunPOT;
                seenSubRunsCosmics.insert(keyCosmics);
            }
        } 
    }

    double cosmicsPOT = cosmicSpillsSum * 5e12;

    double signalWeight = totalPOTSignal/totalPOTSignal;
    double BNBWeight = totalPOTSignal/totalPOTBNB;
    double cosmicsWeight = totalPOTSignal/cosmicsPOT;

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
    tree->SetBranchAddress("reco_particleNumHits", &reco_particleNumHits);
    tree->SetBranchAddress("reco_particleNumMatchedHits", &reco_particleNumMatchedHits);
    tree->SetBranchAddress("reco_particleNumTrueHits", &reco_particleNumTrueHits);
    tree->SetBranchAddress("reco_particleTruePDG", &reco_particleTruePDG);

    tree->SetBranchAddress("reco_sliceID", &reco_sliceID);
    tree->SetBranchAddress("reco_sliceCompleteness", &reco_sliceCompleteness);
    tree->SetBranchAddress("reco_slicePurity", &reco_slicePurity);
    tree->SetBranchAddress("reco_sliceScore", &reco_sliceScore);

    Long64_t numEntries = tree->GetEntries();

    interactionCounter interactionsSignal;
    interactionCounter interactionsBNB;

    auto numSlices = createHistGroup("numSlices", "Number of Slices in an Event", "Number of Slices", 45, 0, 45);
    auto numSlicesCRUMBS = createHistGroup("numSlicesCRUMBS", "Number of Slices with CRUMBS Score in an Event", "Number of Slices with a CRUMBS Score", 10, 0, 10);
    auto numRecoNeutrinos = createHistGroup("numRecoNeutrinos", "Number of Reco Neutrinos in an Event", "Number of Reco Neutrinos", 10, 0, 10);

    auto sliceCompletenessCRUMBS = createHistGroup("sliceCompletenessCRUMBS", "Completeness of the Slice with the Highest CRUMBS Score", "Completeness", 102, 0, 1.02); // Bin width = 0.005
    auto sliceScoreCRUMBS = createHistGroup("sliceScoreCRUMBS", "CRUMBS Score of the Slice with the Highest CRUMBS Score", "CRUMBS Score", 25, -1, 1); 
    auto slicePurityCRUMBS = createHistGroup("slicePurityCRUMBS", "Purity of the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);
    auto highestPFPCompletenessCRUMBS = createHistGroup("highestPFPCompletenessCRUMBS", "Completeness of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "Completeness", 50, 0, 1.02);
    auto highestPFPPurityCRUMBS = createHistGroup("highestPFPPurityCRUMBS", "Purity of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "Purity", 50, 0, 1.02);

    auto numPFPsCRUMBS = createHistGroup("numPFPsCRUMBS", "Number of PFPs in the Slice with the Highest CRUMBS Score", "Number of PFPs", 10, 0, 10);
    auto ratioChosenSummedEnergyCRUMBS = createHistGroup("ratioChosenSummedEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP and the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score", "E_{reco, highest energy PFP}/E_{reco, summed PFP energies}", 21, 0, 1.05);
    auto ratioChosenTrueEnergyCRUMBS = createHistGroup("ratioChosenTrueEnergyCRUMBS", "Ratio of the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, highest energy PFP}/E_{true}", 24, 0, 1.2);
    auto ratioSummedTrueEnergyCRUMBS = createHistGroup("ratioSummedTrueEnergyCRUMBS", "Ratio of the Summed Energy of the PFPs in the Slice with the Highest CRUMBS Score and the True Shower Energy", "E_{reco, summed PFP energies}/E_{true}", 24, 0, 1.2);
   
    auto EtrueThetaRecoCRUMBS = createHistGroup("EtrueThetaRecoCRUMBS", "E_{true}#theta_{reco}^{2}: Slice with Highest CRUMBS Score", "E_{true}#theta_{reco}^{2} (MeV)", 20, 0, 20.44);
    auto ERecoSumThetaTrueCRUMBS = createHistGroup("ERecoSumThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoHighestThetaTrueCRUMBS = createHistGroup("ERecoHighestThetaTrueCRUMBS", "E_{reco}#theta_{true}^{2} for E_{reco} Being the Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{true}^{2} (MeV)", 24, 0, 3.066);
    auto ERecoSumThetaRecoCRUMBS = createHistGroup("ERecoSumThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Sum of Energies of PFPs in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 14, 0, 14.308);
    auto ERecoHighestThetaRecoCRUMBS = createHistGroup("ERecoHighestThetaRecoCRUMBS", "E_{reco}#theta_{reco}^{2} for E_{reco} Being Energy of the Highest Energy PFP in the Slice with the Highest CRUMBS Score", "E_{reco}#theta_{reco}^{2} (MeV)", 14, 0, 14.308);
    

    for(Long64_t i = 0; i < numEntries; ++i){
        tree->GetEntry(i);

        // Change based on what events you want to see
        // DL: 0 = DL Vertexing UBoone, 2 = BDT Vertexing
        // Signal: 1 = Nu+E Elastic, 2 = BNB, 3 = Intime Cosmics
        if(DLCurrent != DLCurrentReq || signal != SignalReq) continue;

        printf("___________________________________________________________________\n");

        std::cout << "DL Current = " << DLCurrent << ", Signal = " << signal << std::endl;
        std::vector<recoParticle> recoParticlesInEvent;
        std::vector<recoNeutrino> recoNeutrinosInEvent;
        std::vector<recoSlice> recoSlicesInEvent;
        std::vector<trueParticle> trueParticlesInEvent;
        
        recoParticle chosenRecoParticleCRUMBS;
        recoNeutrino chosenRecoNeutrinoCRUMBS;
        recoSlice chosenRecoSliceCRUMBS;
        trueNeutrino chosenTrueNeutrino;
        double energyInSlice = 0;
        double numPFPsInEvent = 0;
        double numPFPsInChosenSlice = 0;


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

                //std::cout << "genieMode: " << chosenTrueNeutrino.genieMode << std::endl;
                //std::cout << "genieInteraction: " << chosenTrueNeutrino.genieInteraction << std::endl;
                if(DLCurrent == 0){
                    if(chosenTrueNeutrino.genieMode == -1 && signal == 1) interactionsSignal.Unknown++;
                    if(chosenTrueNeutrino.genieMode == 0 && signal == 1) interactionsSignal.QE++;
                    if(chosenTrueNeutrino.genieMode == 1 && signal == 1) interactionsSignal.Res++;
                    if(chosenTrueNeutrino.genieMode == 2 && signal == 1) interactionsSignal.DIS++;
                    if(chosenTrueNeutrino.genieMode == 3 && signal == 1) interactionsSignal.Coh++;
                    if(chosenTrueNeutrino.genieMode == 4 && signal == 1) interactionsSignal.CohElastic++;
                    if(chosenTrueNeutrino.genieMode == 5 && signal == 1) interactionsSignal.ElectronScattering++;
                    if(chosenTrueNeutrino.genieMode == 6 && signal == 1) interactionsSignal.IMDAnnihilation++;
                    if(chosenTrueNeutrino.genieMode == 7 && signal == 1) interactionsSignal.InverseBetaDecay++;
                    if(chosenTrueNeutrino.genieMode == 8 && signal == 1) interactionsSignal.GlashowResonance++;
                    if(chosenTrueNeutrino.genieMode == 9 && signal == 1) interactionsSignal.AMNuGamma++;
                    if(chosenTrueNeutrino.genieMode == 10 && signal == 1) interactionsSignal.MEC++;
                    if(chosenTrueNeutrino.genieMode == 11 && signal == 1) interactionsSignal.Diffractive++;
                    if(chosenTrueNeutrino.genieMode == 12 && signal == 1) interactionsSignal.EM++;
                    if(chosenTrueNeutrino.genieMode == 13 && signal == 1) interactionsSignal.WeakMix++;
                    if(chosenTrueNeutrino.genieInteraction == 1000 && signal == 1) interactionsSignal.NuanceOffset++;
                    if(chosenTrueNeutrino.genieInteraction == 1001 && signal == 1) interactionsSignal.CCQE++;
                    if(chosenTrueNeutrino.genieInteraction == 1002 && signal == 1) interactionsSignal.NCQE++;
                    if(chosenTrueNeutrino.genieInteraction > 1002 && chosenTrueNeutrino.genieInteraction < 1091 && signal == 1) interactionsSignal.NuanceRes++;
                    if(chosenTrueNeutrino.genieInteraction == 1091 && signal == 1) interactionsSignal.CCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1092 && signal == 1) interactionsSignal.NCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1098 && signal == 1) interactionsSignal.NuEElastic++;

                    if(chosenTrueNeutrino.genieMode == -1 && signal == 2) interactionsBNB.Unknown++;
                    if(chosenTrueNeutrino.genieMode == 0 && signal == 2) interactionsBNB.QE++;
                    if(chosenTrueNeutrino.genieMode == 1 && signal == 2) interactionsBNB.Res++;
                    if(chosenTrueNeutrino.genieMode == 2 && signal == 2) interactionsBNB.DIS++;
                    if(chosenTrueNeutrino.genieMode == 3 && signal == 2) interactionsBNB.Coh++;
                    if(chosenTrueNeutrino.genieMode == 4 && signal == 2) interactionsBNB.CohElastic++;
                    if(chosenTrueNeutrino.genieMode == 5 && signal == 2) interactionsBNB.ElectronScattering++;
                    if(chosenTrueNeutrino.genieMode == 6 && signal == 2) interactionsBNB.IMDAnnihilation++;
                    if(chosenTrueNeutrino.genieMode == 7 && signal == 2) interactionsBNB.InverseBetaDecay++;
                    if(chosenTrueNeutrino.genieMode == 8 && signal == 2) interactionsBNB.GlashowResonance++;
                    if(chosenTrueNeutrino.genieMode == 9 && signal == 2) interactionsBNB.AMNuGamma++;
                    if(chosenTrueNeutrino.genieMode == 10 && signal == 2) interactionsBNB.MEC++;
                    if(chosenTrueNeutrino.genieMode == 11 && signal == 2) interactionsBNB.Diffractive++;
                    if(chosenTrueNeutrino.genieMode == 12 && signal == 2) interactionsBNB.EM++;
                    if(chosenTrueNeutrino.genieMode == 13 && signal == 2) interactionsBNB.WeakMix++;
                    if(chosenTrueNeutrino.genieInteraction == 1000 && signal == 2) interactionsBNB.NuanceOffset++;
                    if(chosenTrueNeutrino.genieInteraction == 1001 && signal == 2) interactionsBNB.CCQE++;
                    if(chosenTrueNeutrino.genieInteraction == 1002 && signal == 2) interactionsBNB.NCQE++;
                    if(chosenTrueNeutrino.genieInteraction > 1002 && chosenTrueNeutrino.genieInteraction < 1091 && signal == 2) interactionsBNB.NuanceRes++;
                    if(chosenTrueNeutrino.genieInteraction == 1091 && signal == 2) interactionsBNB.CCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1092 && signal == 2) interactionsBNB.NCDis++;
                    if(chosenTrueNeutrino.genieInteraction == 1098 && signal == 2) interactionsBNB.NuEElastic++;
                }
            } else if(truth_neutrinoTPCValid->at(j) != -999999 && truth_neutrinoVX->at(j) < 201.3 && truth_neutrinoVX->at(j) > -201.3 && truth_neutrinoVY->at(j) < 203.732 && truth_neutrinoVY->at(j) > -203.732 && truth_neutrinoVZ->at(j) > 0 && truth_neutrinoVZ->at(j) < 509.4){
                //printf("True Neutrino Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            } else if(truth_neutrinoTPCValid->at(j) != -999999){
                //printf("True Neutrino Not Within TPC: TPC Valid = %f, Vertex = (%f, %f, %f), PDG = %f, Lepton PDG = %f\n", truth_neutrinoTPCValid->at(j), truth_neutrinoVX->at(j), truth_neutrinoVY->at(j), truth_neutrinoVZ->at(j), truth_neutrinoType->at(j), truth_chargedLepton->at(j));
            }
        }

        if(truth_neutrinoVX->size() == 1 && truth_neutrinoVX->at(0) == -999999) printf("No True Neutrino in Event\n");
        if(truthNeutrinoInTPC == 0 && signal != 3) continue;

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

        int slice = 0;
        int crumbsSlice = 0;
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
                //printf("Slice %zu: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n", j, recoSlice.id, recoSlice.completeness, recoSlice.purity, recoSlice.score);
            }
        }    

        // If there are no slices, skip the event
        if(slice == 0) continue;  

        // Pick the slice with the highest score
        chosenRecoSliceCRUMBS = chooseSlice(recoSlicesInEvent, 1);
        //printf("Slice with Highest CRUMBS Score: ID = %f, Completeness = %f, Purity = %f, CRUMBS Score = %f\n\n", chosenRecoSliceCRUMBS.id, chosenRecoSliceCRUMBS.completeness, chosenRecoSliceCRUMBS.purity, chosenRecoSliceCRUMBS.score);
        
        int neutrino = 0;
        int numRecoNeutrinosInEvent = 0;

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

                //printf("Reco Neutrino %zu: PDG = %f, Is Primary = %f, Vertex = (%f, %f, %f), Slice ID = %f\n", j, recoNeutrino.pdg, recoNeutrino.isPrimary, recoNeutrino.vx, recoNeutrino.vy, recoNeutrino.vz, recoNeutrino.sliceID);
            }
        }

        // Skip the event if there are no reconstructed neutrinos
        if(neutrino == 0) continue;

        // Look at reconstruction when we pick the slice with highest CRUMBS score
        chosenRecoNeutrinoCRUMBS = chooseRecoNeutrino(recoNeutrinosInEvent, chosenRecoSliceCRUMBS.id); 

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
                recoParticle.completeness = reco_particleCompleteness->at(j);
                recoParticle.purity = reco_particlePurity->at(j);
                recoParticle.numHits = reco_particleNumHits->at(j);
                recoParticle.numMatchedHits = reco_particleNumMatchedHits->at(j);
                recoParticle.numTrueHits = reco_particleNumTrueHits->at(j);
                recoParticle.truePDG = reco_particleTruePDG->at(j);
                recoParticlesInEvent.push_back(recoParticle);
            
                if(recoParticle.sliceID == chosenRecoSliceCRUMBS.id) recoparticleCRUMBS = 1;
            }
        }

        // There are no reco particles in the event
        if(recoparticle == 0) continue;
        
        // No reco particles in the slice with the highest CRUMBS score
        if(recoparticleCRUMBS == 0) continue;

        double totalSliceEnergyCRUMBS = 0;
        double numPFPsSliceCRUMBS = 0;
        chosenRecoParticleCRUMBS = choosePFP(recoParticlesInEvent, chosenRecoSliceCRUMBS.id, totalSliceEnergyCRUMBS, numPFPsSliceCRUMBS);
        double chosenSlicePurityCRUMBS = slicePurityCalculator(recoParticlesInEvent, chosenRecoSliceCRUMBS.id);
        double chosenSliceCompletenessCRUMBS = sliceCompletenessCalculator(recoParticlesInEvent, chosenRecoSliceCRUMBS.id);

        //std::cout << "highest energy PFP PDG: " << chosenRecoParticleCRUMBS.truePDG << std::endl;
        //std::cout << chosenRecoParticleCRUMBS.truePDG << std::endl;

        if(chosenTrueNeutrino.genieMode == -1){
            numSlices.Unknown->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.Unknown->Fill(crumbsSlice);
            sliceScoreCRUMBS.Unknown->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.Unknown->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.Unknown->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.Unknown->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.Unknown->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.Unknown->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.Unknown->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.Unknown->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.Unknown->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.Unknown->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 0 && chosenTrueNeutrino.genieInteraction == 1001){
            numSlices.CCQE->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.CCQE->Fill(crumbsSlice);
            sliceScoreCRUMBS.CCQE->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.CCQE->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.CCQE->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.CCQE->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.CCQE->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.CCQE->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.CCQE->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.CCQE->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.CCQE->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.CCQE->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 0 && chosenTrueNeutrino.genieInteraction == 1002){
            numSlices.NCQE->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.NCQE->Fill(crumbsSlice);
            sliceScoreCRUMBS.NCQE->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.NCQE->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.NCQE->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.NCQE->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.NCQE->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.NCQE->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.NCQE->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.NCQE->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.NCQE->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.NCQE->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 1){
            numSlices.Res->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.Res->Fill(crumbsSlice);
            sliceScoreCRUMBS.Res->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.Res->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.Res->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.Res->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.Res->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.Res->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.Res->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.Res->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.Res->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.Res->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 2 && chosenTrueNeutrino.genieInteraction == 1091){
            numSlices.CCDIS->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.CCDIS->Fill(crumbsSlice);
            sliceScoreCRUMBS.CCDIS->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.CCDIS->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.CCDIS->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.CCDIS->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.CCDIS->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.CCDIS->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.CCDIS->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.CCDIS->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.CCDIS->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.CCDIS->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 2 && chosenTrueNeutrino.genieInteraction == 1092){
            numSlices.NCDIS->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.NCDIS->Fill(crumbsSlice);
            sliceScoreCRUMBS.NCDIS->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.NCDIS->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.NCDIS->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.NCDIS->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.NCDIS->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.NCDIS->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.NCDIS->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.NCDIS->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.NCDIS->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.NCDIS->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 3){
            numSlices.Coh->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.Coh->Fill(crumbsSlice);
            sliceScoreCRUMBS.Coh->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.Coh->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.Coh->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.Coh->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.Coh->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.Coh->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.Coh->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.Coh->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.Coh->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.Coh->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 5 && chosenTrueNeutrino.genieInteraction == 1098){
            numSlices.NuEElastic->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.NuEElastic->Fill(crumbsSlice);
            sliceScoreCRUMBS.NuEElastic->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.NuEElastic->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.NuEElastic->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.NuEElastic->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.NuEElastic->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.NuEElastic->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.NuEElastic->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.NuEElastic->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.NuEElastic->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.NuEElastic->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else if(chosenTrueNeutrino.genieMode == 10){
            numSlices.MEC->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.MEC->Fill(crumbsSlice);
            sliceScoreCRUMBS.MEC->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.MEC->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.MEC->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.MEC->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.MEC->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.MEC->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.MEC->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.MEC->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.MEC->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.MEC->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        } else{
            numSlices.other->Fill(numSlicesInEvent); 
            numSlicesCRUMBS.other->Fill(crumbsSlice);
            sliceScoreCRUMBS.other->Fill(chosenRecoSliceCRUMBS.score);
            numRecoNeutrinos.other->Fill(numRecoNeutrinosInEvent);
            sliceCompletenessCRUMBS.other->Fill(chosenSliceCompletenessCRUMBS);
            slicePurityCRUMBS.other->Fill(chosenSlicePurityCRUMBS);
            highestPFPCompletenessCRUMBS.other->Fill(chosenRecoParticleCRUMBS.completeness);
            highestPFPPurityCRUMBS.other->Fill(chosenRecoParticleCRUMBS.purity);
            numPFPsCRUMBS.other->Fill(numPFPsSliceCRUMBS);
            ratioChosenSummedEnergyCRUMBS.other->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy / totalSliceEnergyCRUMBS);
            ERecoSumThetaRecoCRUMBS.other->Fill(totalSliceEnergyCRUMBS * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
            ERecoHighestThetaRecoCRUMBS.other->Fill(chosenRecoParticleCRUMBS.bestPlaneEnergy * chosenRecoParticleCRUMBS.theta * chosenRecoParticleCRUMBS.theta);
        }
    }

    styleDraw(numSlices, (base_path + "numSlices_dist.pdf").c_str());
    styleDraw(numSlicesCRUMBS, (base_path + "numSlicesCRUMBS_dist.pdf").c_str());
    styleDraw(sliceScoreCRUMBS, (base_path + "sliceScoreCRUMBS_dist.pdf").c_str());
    styleDraw(numRecoNeutrinos, (base_path + "numRecoNeutrinos_dist.pdf").c_str());
    styleDraw(sliceCompletenessCRUMBS, (base_path + "sliceCompletenessCRUMBS_dist.pdf").c_str());
    styleDraw(slicePurityCRUMBS, (base_path + "slicePurityCRUMBS_dist.pdf").c_str());
    styleDraw(highestPFPCompletenessCRUMBS, (base_path + "highestPFPCompleteness_dist.pdf").c_str());
    styleDraw(highestPFPPurityCRUMBS, (base_path + "highestPFPPurityCRUMBS_dist.pdf").c_str());
    styleDraw(numPFPsCRUMBS, (base_path + "numPFPsCRUMBS_dist.pdf").c_str());
    styleDraw(ratioChosenSummedEnergyCRUMBS, (base_path + "ratioChosenSummedEnergyCRUMBS_dist.pdf").c_str());
    styleDraw(ERecoSumThetaRecoCRUMBS, (base_path + "ERecoSumThetaRecoCRUMBS_dist.pdf").c_str());
    styleDraw(ERecoHighestThetaRecoCRUMBS, (base_path + "ERecoHighestThetaRecoCRUMBS_dist.pdf").c_str());
}
