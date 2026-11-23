////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "STPCAnalyzerTreeReader.h"

// --- Destructor ---
MyTPCTreeReader::~MyTPCTreeReader(){
    // Clean up when the class is destroyed
    delete file;
}

// --- Constructor ---
MyTPCTreeReader::MyTPCTreeReader(TString fileName, std::string treeName) {
    // Open the ROOT file
    file = new TFile(fileName);

    // Get the TTree from the file (assuming the TTree is named "myTree")
    tree = (TTree *)file->Get(treeName.c_str());

    fNEntries = tree->GetEntries();

    // Set branch addresses for event information
    tree->SetBranchAddress("EventID", &eventID);
    tree->SetBranchAddress("SubRunID", &subrunID);
    tree->SetBranchAddress("RunID", &runID);

    // Set branch addresses for true vertex information
    tree->SetBranchAddress("TrueVEnergy", &nuvE);
    tree->SetBranchAddress("TrueVt", &nuvT);
    tree->SetBranchAddress("TrueVx", &nuvX);
    tree->SetBranchAddress("TrueVy", &nuvY);
    tree->SetBranchAddress("TrueVz", &nuvZ);
    tree->SetBranchAddress("TrueVU", &nuvU);
    tree->SetBranchAddress("TrueVV", &nuvV);
    tree->SetBranchAddress("TrueVC", &nuvC);
    tree->SetBranchAddress("TrueVTimeTick", &nuvTimeTick);

    tree->SetBranchAddress("IntMode", &intMode);
    tree->SetBranchAddress("IntNLambda", &intNLambda);

    // Set branch addresses for reco vertex information
    tree->SetBranchAddress("RecoVx", &recnuvX);
    tree->SetBranchAddress("RecoVy", &recnuvY);
    tree->SetBranchAddress("RecoVz", &recnuvZ);
    tree->SetBranchAddress("RecoVU", &recnuvU);
    tree->SetBranchAddress("RecoVV", &recnuvV);
    tree->SetBranchAddress("RecoVC", &recnuvC);
    tree->SetBranchAddress("RecoVTimeTick", &recnuvTimeTick);

    // True PDGs
    tree->SetBranchAddress("TruePrimariesPDG", &truePrimeriesPDG);
    
    // Set branch addresses for hits information
    tree->SetBranchAddress("HitsView", &hitsView);
    tree->SetBranchAddress("HitsIntegral", &hitsIntegral);
    tree->SetBranchAddress("HitsSummedADC", &hitsSummedADC);
    tree->SetBranchAddress("HitsPeakTime", &hitsPeakTime);
    tree->SetBranchAddress("HitsChannel", &hitsChannel);
    tree->SetBranchAddress("HitsRMS", &hitsRMS);
    tree->SetBranchAddress("HitsStartT", &hitsStartT);
    tree->SetBranchAddress("HitsEndT", &hitsEndT);
    tree->SetBranchAddress("HitsChi2", &hitsChi2);
    tree->SetBranchAddress("HitsNDF", &hitsNDF);
    tree->SetBranchAddress("HitsClusterID", &hitsClusterID);
    tree->SetBranchAddress("HitsX", &hitsX);
    tree->SetBranchAddress("HitsY", &hitsY);
    tree->SetBranchAddress("HitsZ", &hitsZ);

    tree->SetBranchAddress("LambdaProtonPDir", &lambdaProtonPDir);
    tree->SetBranchAddress("LambdaPionPDir", &lambdaPionPDir);
}


// --- Get the number of entries ---
int MyTPCTreeReader::NEntries(){return fNEntries;}


// --- Get the entry ---
bool MyTPCTreeReader::GetEntry(int entry) {
    if (entry < 0 || entry >= tree->GetEntries())
        return false;
    
    tree->GetEntry(entry);
    return true;
}


// --- Get the hits in a given view ---
std::vector<SHit> MyTPCTreeReader::GetHitsInView(int view)
{
    // set variables
    std::vector<SHit> hitList;
    int nTotalHits = hitsChannel->size();

    // loop over the hits
    for (int i = 0; i < nTotalHits; i++) {

        // filter channels for the view        
        if ( hitsView->at(i)==view ) {
            //SHit hit(-1, _X->at(i), _Y->at(i), _Wi->at(i), _Int->at(i), _ST->at(i), _ET->at(i), _Chi2->at(i), _ClusterId->at(i));
            SHit hit(i, hitsChannel->at(i), hitsPeakTime->at(i), hitsRMS->at(i), hitsIntegral->at(i),
                    hitsStartT->at(i), hitsEndT->at(i), hitsChi2->at(i), hitsClusterID->at(i),
                    hitsX->at(i), hitsY->at(i), hitsZ->at(i) );
            hitList.push_back(hit);
        }
    }

    return hitList;

}
