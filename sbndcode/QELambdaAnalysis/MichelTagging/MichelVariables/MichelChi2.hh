////////////////////////////////////////////////////////////////////////
// Class:       MichelChi2
// File:        MichelChi2.hh
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


/**
 * @brief Class to obtain the variable Lower Chi2 ratio from the electron calo points
 */
class MichelChi2 {
    private:
        // Inputs
        GeoParams fParams;
        std::vector<MuCalo> fmuCalos;
        Bragg fbragg;

        
    public:
        // --------- Constructor --------------
        MichelChi2(
            const GeoParams& params,
            const std::vector<MuCalo>& muCalos, 
            const Bragg& bragg
            )
            : fParams(params)
            , fmuCalos(muCalos)
            , fbragg(bragg)
        {}

        // --------- Main function --------------
        double LowerChi2Ratio();
        double AllChi2Ratio();


};