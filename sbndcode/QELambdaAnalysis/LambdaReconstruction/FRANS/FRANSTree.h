///////////////////////////////////////////////////////////////////////
/// File: ChargeDensity.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////
#ifndef SBND_FRANSTREEHANDLE_H
#define SBND_FRANSTREEHANDLE_H

#include <iostream>
#include <TTree.h>
#include "FRANSObj.h"

class FRANSTree {
public:
    // constructors
    FRANSTree();
    
    FRANSTree(TTree* tree, bool readMode);

    void SetTree(TTree* tree, bool readMode);

    // Fill functions
    void FillDataMC(int runID, int subRunID, int eventID, int nnuints, int inttype,  int intmode,
                    double trueVx, double trueVy, double trueVz,
                    double gap, int issignal, double keproton, double kepion, double kelambda);
    
    void FillData(int view, FRANSObj fransObj);
    
    void FillTree();

private:

    // Set branch function
    template <typename T>
    void SetBranch(const char* name, T* variable, bool setAddress);

    // Initialize function
    void InitializeTree(bool readMode);

    TTree* fTree;
    std::vector<int> fViews;

    // Event label variables
    int fRunID;
    int fSubRunID;
    int fEventID;

    // Lambda MC variables
    int fNNuInts;
    int fIntMode;
    int fIntType;
    double fGap;
    double fLambdaKE;
    double fProtonKE;
    double fPionKE;
    int fIsSignal;
    double fTrueVx;
    double fTrueVy;
    double fTrueVz;

    // FRANS objects
    FRANSObj fFRANSObj0;
    FRANSObj fFRANSObj1;
    FRANSObj fFRANSObj2;
};

#endif
