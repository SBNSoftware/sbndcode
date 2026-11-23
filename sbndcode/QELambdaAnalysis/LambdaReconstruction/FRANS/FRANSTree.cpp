///////////////////////////////////////////////////////////////////////
/// File: FRANSTreeHandle.cpp
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////

#include "FRANSTree.h"

template <typename T>
void FRANSTree::SetBranch(const char* name, T* variable, bool setAddress) {
    if (setAddress) {
        fTree->SetBranchAddress(name, variable);
    } else {
        fTree->Branch(name, variable);
    }
}

// instance types
template void FRANSTree::SetBranch(const char* name, int* variable, bool setAddress);
template void FRANSTree::SetBranch(const char* name, float* variable, bool setAddress);
template void FRANSTree::SetBranch(const char* name, double* variable, bool setAddress);
template void FRANSTree::SetBranch(const char* name, bool* variable, bool setAddress);
template void FRANSTree::SetBranch(const char* name, std::string* variable, bool setAddress);
template void FRANSTree::SetBranch(const char* name, FRANSObj* variable, bool setAddress);


FRANSTree::FRANSTree()
    : fTree(nullptr)
{}


FRANSTree::FRANSTree(TTree* tree, bool readMode)
    : fTree(tree)
{
    InitializeTree(readMode);
}

void FRANSTree::SetTree(TTree* tree, bool readMode)
{
    fTree = tree;
    InitializeTree(readMode);
}


void FRANSTree::FillDataMC( int runID, int subRunID, int eventID, int nnuints, int inttype, int intmode,
                            double trueVx, double trueVy, double trueVz,
                            double gap, int issignal, double keproton, double kepion, double kelambda)
{
    fRunID = runID;
    fSubRunID = subRunID;
    fEventID = eventID;
    
    fNNuInts = nnuints;
    fIntType = inttype;
    fIntMode = intmode;

    fTrueVx = trueVx;
    fTrueVy = trueVy;
    fTrueVz = trueVz;
    
    fGap = gap;
    fIsSignal = issignal;
    fProtonKE=keproton;
    fPionKE=kepion;
    fLambdaKE=kelambda;
}

void FRANSTree::FillData(int view, FRANSObj fransObj){
    if (view == 0) fFRANSObj0 = fransObj;
    else if (view == 1) fFRANSObj1 = fransObj;
    else if (view == 2) fFRANSObj2 = fransObj;
    else std::cout << "ERROR: FRANSTree::FillData() - view not valid" << std::endl;
}

void FRANSTree::FillTree(){ fTree->Fill(); }


void FRANSTree::InitializeTree(bool readMode)
{
    // Event label variables
    SetBranch("RunID", &fRunID, readMode);
    SetBranch("SubRunID", &fSubRunID, readMode);
    SetBranch("EventID", &fEventID, readMode);

    // Lambda MC variables
    SetBranch("NNuInts", &fNNuInts, readMode);
    SetBranch("IntType", &fIntType, readMode);
    SetBranch("IntMode", &fIntMode, readMode);
    SetBranch("TrueVx", &fTrueVx, readMode);
    SetBranch("TrueVy", &fTrueVy, readMode);
    SetBranch("TrueVz", &fTrueVz, readMode);
    SetBranch("Gap", &fGap, readMode);
    SetBranch("IsSignal", &fIsSignal, readMode);
    SetBranch("ProtonKE", &fProtonKE, readMode);
    SetBranch("PionKE", &fPionKE, readMode);
    SetBranch("LambdaKE", &fLambdaKE, readMode);

    // TTree for BDT training
    SetBranch("FRANSObj0", &fFRANSObj0, readMode);
    SetBranch("FRANSObj1", &fFRANSObj1, readMode);
    SetBranch("FRANSObj2", &fFRANSObj2, readMode);
}