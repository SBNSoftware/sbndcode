////////////////////////////////////////////////////////////////////////
// Class:       MichelClassifier
// File:        MichelClassifier.cc
//
// Created at Thu Jul 2 2026 by Sergio Yebes Baena 
////////////////////////////////////////////////////////////////////////

#pragma once

#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelObj.hh"

// ROOT and C++ includes
#include <TTree.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include "TMVA/RBDT.hxx"
#include "TMVA/Reader.h"

#include <string.h>
#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <numeric>
#include <optional>


class MichelClassifier {
    private:
        // Inputs
        TMVA::Reader fReader;

        bool fLoaded = false;

        float fNMuTracks;
        float fLowerChi2Ratio;
        float fAngleQ80;
        float fDropRadius;
        float fNBlips;

    public:

        // --------- Constructor --------------
        MichelClassifier() : fReader("!Color:!Silent") {}

        void Load(std::string const& weightsFile);

        // Getters
        bool Loaded() const { return fLoaded; }

        // --------- Main function --------------
        float Evaluate (const MichelFeatures& michelFeatures);
};