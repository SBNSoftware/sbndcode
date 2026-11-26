////////////////////////////////////////////////////////////////////////
// Class:       MichelBlips
// File:        MichelBlips.hh
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
#include <algorithm>



/**
 * @brief Class to obtain the variable Angle q80 from hits
 */
class MichelBlips {
    private:
        // Inputs
        int fmaxWireSpan = 3;
        float fmaxBlipDist;
        std::vector<MuBlip> fmuBlips;
        std::vector<MuTrack> fmuTrack;
        Bragg fbragg;
        double frecoVtxX;
        double frecoVtxY;
        double frecoVtxZ;
        
    public:
        // --------- Constructor --------------
        MichelBlips(
            const int maxWireSpan,
            const float maxBlipDist,
            const std::vector<MuBlip>& muBlips,
            const std::vector<MuTrack>& muTrack,
            const Bragg& bragg,
            const double recoVtxX,
            const double recoVtxY,
            const double recoVtxZ)
            : fmaxWireSpan(maxWireSpan)
            , fmaxBlipDist(maxBlipDist)
            , fmuBlips(muBlips)
            , fmuTrack(muTrack)
            , fbragg(bragg)
            , frecoVtxX(recoVtxX)
            , frecoVtxY(recoVtxY)
            , frecoVtxZ(recoVtxZ)
        {}


        int NBlips();
};
