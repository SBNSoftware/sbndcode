////////////////////////////////////////////////////////////////////////
// Class:       MichelAngles
// File:        MichelAngles.hh
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
#include <TF1.h>

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


/**
 * @brief Class to obtain the variable Angle q80 from hits
 */
class MichelAngles {
    private:
        // Inputs
        double fmaxRadius;
        std::vector<MuHit> fmuHits;
        Bragg fbragg;
        PCAResult fpca;

        std::vector<double> fangles;
        
    public:
        // --------- Constructor --------------
        MichelAngles(
            const double maxRadius,
            const std::vector<MuHit>& muHits, 
            const Bragg& bragg,
            const PCAResult& pca
            )
            : fmaxRadius(maxRadius)
            , fmuHits(muHits)
            , fbragg(bragg)
            , fpca(pca)
        {}

        std::vector<double> getAngles() { return fangles; }

        double AnglesQ80();
};
