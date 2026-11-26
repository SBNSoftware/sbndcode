////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesHoughs
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_HOUGH_H
#define TPC_LINES_HOUGH_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <numeric>


#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleHits.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleLines.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/TPCLinesAlgo/TPCLinesParameters.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/TPCLinesAlgo/TPCLinesDisplay.h"

#else
#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"
#include "TPCLinesParameters.h"
#include "TPCLinesDisplay.h"
#endif





class TPCLinesHough {
    private:
        HoughAlgorithmPsetType fTPCLinesHoughPset;

        std::vector<SHit> GetHitsInBall(std::vector<SHit> hitList, SVertex vertex, double d);
        LineEquation ComputeLineEquationFromHough(double th, SPoint p);

        // weight function
        double HoughWeightDistanceStep(LineEquation line, std::vector<SHit> hitList);
        double HoughWeightDistance(LineEquation line, std::vector<SHit> hitList);
        double HoughWeightDistanceIntegral(LineEquation line, std::vector<SHit> hitList);

        // alternative weight functions
        double HoughWeightAverageDistance(LineEquation line, std::vector<SHit> hitList);
        double GetLinearR2(std::vector<double> hypoV, std::vector<double> dataV);
        double HoughWeightLinearR2(LineEquation line, std::vector<SHit> hitList);

        // display app directory
        TPCLinesDisplay fDisplay;
        std::string fDisplayAppPath;


    public:
        // constructor
        TPCLinesHough(HoughAlgorithmPsetType tpcLinesHoughPset);

        // main function
        HoughLine GetBestHoughLine(std::vector<SHit> hitList, SVertex vertex);
        
};

#endif // TPC_SIMPLE_LINES_H