////////////////////////////////////////////////////////////////////////
// Class:       MichelRadial
// File:        MichelRadial.hh
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
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
#include <Math/QuantFuncMathCore.h>

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
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>


/**
 * @brief Class to obtain the variable Angle q80 from hits
 */
class MichelRadial {
    private:
        // Inputs
        RadialParams fParams;
        std::vector<MuHit> fmuHits;
        Bragg fbragg;
        PCAResult fpca;

        double fdropRadius;
        std::vector<double> fradii;
        std::vector<double> fcharge;
        std::vector<double> fdQdr;
        std::vector<double> fsmoothdQdr;
        int fdropIdx;
        double fmindQdr;
        int fdropThr;
        bool fdropFound;
        
    public:
        // --------- Constructor --------------
        MichelRadial(
            const RadialParams& params,
            const std::vector<MuHit>& muHits, 
            const Bragg& bragg,
            const PCAResult& pca
            )
            : fParams(params)
            , fmuHits(muHits)
            , fbragg(bragg)
            , fpca(pca)
        {}

        // --------- Getters --------------
        double getDropRadius() { return fdropRadius; }
        std::vector<double> getRadii() { return fradii; }
        std::vector<double> getCharge() { return fcharge; }
        std::vector<double> getdQdr() { return fdQdr; }
        std::vector<double> getSmoothdQdr() { return fsmoothdQdr; }
        int getDropIdx() { return fdropIdx; }
        double getMindQdr() { return fmindQdr; }

        
        // --------- Main functions --------------
        std::pair<std::vector<double>, std::vector<double>> getChargeProfile(const std::vector<MuHit>& hits);
        void Compute();
};
