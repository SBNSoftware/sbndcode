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

typedef struct{
    uint16_t mac5;
    uint16_t flags;
    uint32_t ts0;
    uint32_t ts1;
    uint32_t unixs;
    std::vector<uint16_t> adcVal;
    uint32_t coinc;
} FEBData_t;

void crtMacro(){

    TFile *f = TFile::Open("/exp/sbnd/data/users/coackley/output_crtanalyser.root");
    if(!f){
        std::cout << "Failed to read file" << std::endl;
        return;
    }

    TTree *t;
    f->GetObject("fTree", t);

    std::vector<uint16_t> allMac5 = std::vector<uint16_t>(0);
    std::vector<uint16_t> allFlags = std::vector<uint16_t>(0);
    std::vector<uint32_t> allTs0 = std::vector<uint32_t>(0);
    std::vector<uint32_t> allTs1 = std::vector<uint32_t>(0);
    std::vector<uint32_t> allUnixs = std::vector<uint32_t>(0);
    std::vector<std::vector<uint16_t>> allADC = std::vector<std::vector<uint16_t>>(0);
    std::vector<uint32_t> allCoinc = std::vector<uint32_t>(0);

    Long64_t numEntries = t->GetEntries();
    int fileNumber = 1;

    for(Long64_t i = 0; i < numEntries; i++){
        std::cout << "file num: " << fileNumber << std::endl;
        
        std::vector<uint16_t> *pMac5 = 0;
        TBranch *mac5_branch = 0;
        t->SetBranchAddress("feb_mac5", &pMac5, &mac5_branch);
        mac5_branch->GetEntry(i);

        std::vector<uint16_t> *pFlags = 0;
        TBranch *flags_branch = 0;
        t->SetBranchAddress("feb_flags", &pFlags, &flags_branch);
        flags_branch->GetEntry(i);

        std::vector<uint32_t> *pTs0 = 0;
        TBranch *ts0_branch = 0;
        t->SetBranchAddress("feb_ts0", &pTs0, &ts0_branch);
        ts0_branch->GetEntry(i);

        std::vector<uint32_t> *pTs1 = 0;
        TBranch *ts1_branch = 0;
        t->SetBranchAddress("feb_ts1", &pTs1, &ts1_branch);
        ts1_branch->GetEntry(i);

        std::vector<uint32_t> *pUnixs = 0;
        TBranch *unixs_branch = 0;
        t->SetBranchAddress("feb_unixs", &pUnixs, &unixs_branch);
        unixs_branch->GetEntry(i);

        std::vector<std::vector<uint16_t>> *pADC = 0;
        TBranch *ADC_branch = 0;
        t->SetBranchAddress("feb_adc", &pADC, &ADC_branch);
        ADC_branch->GetEntry(i);

        std::vector<uint32_t> *pCoinc = 0;
        TBranch *coinc_branch = 0;
        t->SetBranchAddress("feb_coinc", &pCoinc, &coinc_branch);
        coinc_branch->GetEntry(i);

        size_t vecSize = pCoinc->size();
        for(UInt_t j = 0; j < vecSize; j++){
            allMac5.push_back(pMac5->at(j));
            allFlags.push_back(pFlags->at(j));
            allTs0.push_back(pTs0->at(j));
            allTs1.push_back(pTs1->at(j));
            allUnixs.push_back(pUnixs->at(j));
            allADC.push_back(pADC->at(j));
            allCoinc.push_back(pCoinc->at(j));
        }
        fileNumber++;
    }

    std::vector<FEBData_t> allFEB = std::vector<FEBData_t>();

    size_t newNumEvents = allMac5.size();

    for(UInt_t k = 0; k < newNumEvents; k++){
        FEBData_t FEB;
        FEB.mac5 = allMac5.at(k);
        FEB.flags = allFlags.at(k);
        FEB.ts0 = allTs0.at(k);
        FEB.ts1 = allTs1.at(k);
        FEB.unixs = allUnixs.at(k);
        FEB.adcVal = allADC.at(k);
        FEB.coinc = allCoinc.at(k);
        allFEB.push_back(FEB);

        //std::cout << "ADC Vec Size: " << FEB.adcVal.size() << std::endl;
        //std::cout << "coinc: " << FEB.coinc << std::endl;
    }


}
