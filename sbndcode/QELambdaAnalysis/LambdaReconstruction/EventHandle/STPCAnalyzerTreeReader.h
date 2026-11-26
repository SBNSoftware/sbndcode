////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_STPCANALYZERTREEREADER_H
#define TPC_STPCANALYZERTREEREADER_H

#include <vector>

// ROOT includess
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleHits.h"
#else
#include "TPCSimpleHits.h"
#endif


// Class to read the TPCAnalyzer TTree
class MyTPCTreeReader {
private:
    TFile *file;    // Pointer to the ROOT file containing the TTree
    TTree *tree;    // Pointer to the TTree object

    int fNEntries;
    
public:
    
    // --- Data members ---
    // Event information
    int eventID;
    int subrunID;
    int runID;
    
    // True Vertex information
    double nuvE, nuvT, nuvX, nuvY, nuvZ;
    int nuvU, nuvV, nuvC, nuvTimeTick;

    // True interaction mode
    int intMode;
    int intNLambda;

    // MC PDGs
    std::vector<int> *truePrimeriesPDG = new std::vector<int>;
    
    // Reco Vertex information
    double recnuvX;
    double recnuvY;
    double recnuvZ;
    int recnuvU;
    int recnuvV;
    int recnuvC;
    int recnuvTimeTick;

    // Hits information
    std::vector<int> *hitsView = new std::vector<int>;
    std::vector<double> *hitsIntegral = new std::vector<double>;
    std::vector<double> *hitsSummedADC = new std::vector<double>;
    std::vector<double> *hitsPeakTime = new std::vector<double>;
    std::vector<int> *hitsChannel = new std::vector<int>;
    std::vector<double> *hitsRMS = new std::vector<double>;
    std::vector<double> *hitsStartT = new std::vector<double>;
    std::vector<double> *hitsEndT = new std::vector<double>;
    std::vector<double> *hitsChi2 = new std::vector<double>;
    std::vector<double> *hitsNDF = new std::vector<double>;
    std::vector<int> *hitsClusterID = new std::vector<int>;
    std::vector<double> *hitsX = new std::vector<double>;
    std::vector<double> *hitsY = new std::vector<double>;
    std::vector<double> *hitsZ = new std::vector<double>;
    std::vector<double> *lambdaProtonPDir = new std::vector<double>;
    std::vector<double> *lambdaPionPDir = new std::vector<double>;

    // --- Methods ---
    MyTPCTreeReader(TString fileName, std::string treeName);
    ~MyTPCTreeReader();
    bool GetEntry(int entry);
    int NEntries();
    std::vector<SHit> GetHitsInView(int view);


};

#endif