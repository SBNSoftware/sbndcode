////////////////////////////////////////////////////////////////////////
// Class:       MichelGeo
// File:        MichelGeo.hh
//
// Created at Thu Jul 2 2026 by Sergio Yebes Baena 
////////////////////////////////////////////////////////////////////////

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

/**
 * @brief Class to perform Bragg identification and muon PCA using calo points of a given plane
 */
class MichelGeo {
    private:
        // Inputs
        GeoParams fParams;
        std::vector<MuCalo> fmuCalos;
        int fplane;

        Bragg fbragg;
        PCAResult fpcaResult;

        // --------- Bragg finder --------------
        void findBraggPeak();
        void getMuonDirection();

    public:

        // --------- Constructor --------------
        MichelGeo(
            const GeoParams& params, 
            const std::vector<MuCalo>& muCalos)
            : fParams(params)
            , fmuCalos(muCalos)
        {}

        // --------- Getters --------------
        Bragg getBragg() const { return fbragg; }
        PCAResult getPCAResult() const { return fpcaResult; }

        // --------- Main function --------------
        void Compute();
};