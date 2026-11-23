////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////


#ifndef TPC_SIMPLE_CALO_H
#define TPC_SIMPLE_CALO_H


#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include <TVectorD.h>
#include <TMatrixD.h>
#include "TVirtualFitter.h"
#include <TArrow.h>

#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleHits.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleTriangles.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/TPCLinesAlgo/TPCLinesParameters.h"
#else
#include "TPCSimpleHits.h"
#include "TPCSimpleTriangles.h"
#include "TPCLinesParameters.h"
#endif

// --- Drawing options ---
const double fMaxMarkerSize = 5;
const double fAlpha = 0.6;
const double fColor1 = kAzure+7;
const double fColor2 = kOrange+8;

class SCalo {
public:
    // Constructor to initialize the collection
    SCalo(const std::vector<SHit>& points = {}, double angle=0, double yAngle=0);

    // Method to display the calculated path lengths
    void Display();

    std::vector<double> GetResidualRange() const {return fResidualRange;};
    std::vector<double> GetDepositedEnergy() const {return fDepositedEnergy;};
    double GetTrackLength() const {return fTrackLength;};
    double GetResidualRangeRMS() const {return fResidualRangeRMS;};

private:

    // Conversion factors
    double fHitIntegralToEnergy;
    double fStepXToLength;
    double fStepYToLength;
    
    // Track legnth
    double fTrackLength;
    // TrackOrientationAngle
    double fCosTrackAngle;
    double fCosY;
    double fTrackAngle;
    // Define the hit points and path lengths as private members
    std::vector<SHit> fHitList;
    std::vector<double> fPathLengths;
    std::vector<double> fDepositedEnergy;
    std::vector<double> fResidualRange;
    // Residual range RMS
    double fResidualRangeRMS;

    // Function to calculate the distance between two hit points
    double CalculateDistance(const SHit& p1, const SHit& p2);

    // Method to calculate and store the path lengths
    void CalculatePathLengths();

    // Function to calculate residual range and deposited energy
    void CalculateResidualRange();
    
};



class STriangleCalo {
    
    private:
        // Algorithm parameters
        CaloAlgorithmPsetType fCaloPset;

        STriangle fTriangle;
        SCalo fCalo1;
        SCalo fCalo2;

        double fRatioUpperLimit;

        std::vector<SHit> fAllHits, fHitsTrack1, fHitsTrack2, fVertexHits, fResidualHits, fFitHits;
        // Fit parameters
        bool fPassFit;
        double fM1;
        double fN1;
        double fM2;
        double fXB;
        double fN2;

        // Range
        double fXVertex;
        double fYVertex;
        double fXFitVtx;
        double fYFitVtx;
        double fMinX;
        double fMaxX;
        double fMinY;
        double fMaxY;
        double fMaxIntegral;

        double fChargeRatioAverage;
        double fChargeDifferenceAverage;
        double fChargeRelativeDifferenceAverage;

        double fChargeRatioIntegral;
        double fChargeDifferenceIntegral;
        double fChargeRelativeDifferenceIntegral;

        bool fPassChargeFit;
        double fBandOverlap;
        double fBandCrossHits;
        double fChargeRatioFit;
        double fChargeDifferenceFit;
        double fChargeRelativeDifferenceFit;

        int fNVertexHits;
        int fNBulkHits;
        double fVertexHitIntegralRatio;
        double fVertexHitIntegralDifference;
        double fVertexHitIntegralRelativeDifference;

        double fTrackLength1;
        double fTrackLength2;
        double fTrackLengthRatio;
        double fResidualRange1RMS;
        double fResidualRange2RMS;
        double fResidualRangeMinRMS;
        double fResidualRangeMaxRMS;
        double fResidualRange1AngleRMS;
        double fResidualRange2AngleRMS;
        double fResidualRangeMinAngleRMS;
        double fResidualRangeMaxAngleRMS;

        double fTwoLinesChi2;

        void MakeEnergyLossVsResidualRangePlot(SCalo calo1, SCalo calo2, TPad *pad);

        // Function to print calorimetry results
        void PrintResults();
    public:
        STriangleCalo(CaloAlgorithmPsetType caloPset);
        
        void SetTriangle(STriangle triangle);

        // Function to display dEdx
        void CreateEnergyLossVsResidualRangePlot();
        // Triangle joint fit analysis
        void JointFitAnalysis(std::vector<double> *pProton=nullptr, std::vector<double> *pPion=nullptr);
        void JointFitAnalysisFisher();
        void Display(TCanvas *c1);

        // Return functions
        double ChargeRatioAverage() {return fChargeRatioAverage;};
        double ChargeDifferenceAverage() {return fChargeDifferenceAverage;};
        double ChargeRelativeDifferenceAverage() {return fChargeRelativeDifferenceAverage;};
        double ChargeRatioIntegral() {return fChargeRatioIntegral;};
        double ChargeDifferenceIntegral() {return fChargeDifferenceIntegral;};
        double ChargeRelativeDifferenceIntegral() {return fChargeRelativeDifferenceIntegral;};
        bool PassFit() {return fPassFit;};
        double ChargeRatioFit() {return fChargeRatioFit;};
        double ChargeDifferenceFit() {return fChargeDifferenceFit;};
        double ChargeRelativeDifferenceFit() {return fChargeRelativeDifferenceFit;};
        int NVertexHits() {return fNVertexHits;};
        int NBulkHits() {return fNBulkHits;};
        double VertexHitIntegralRatio() {return fVertexHitIntegralRatio;};
        double VertexHitIntegralDifference() {return fVertexHitIntegralDifference;};
        double VertexHitIntegralRelativeDifference() {return fVertexHitIntegralRelativeDifference;};
        bool PassChargeFit() {return fPassChargeFit;};
        double BandOverlap() {return fBandOverlap;};
        double BandCrossHits() {return fBandCrossHits;};
        double TwoLinesChi2() {return fTwoLinesChi2;};
        double TrackLengthRatio() {return fTrackLengthRatio;};
        double TrackLength1() {return fTrackLength1;};
        double TrackLength2() {return fTrackLength2;};
        double ResidualRange1RMS() {return fResidualRange1RMS;};
        double ResidualRange2RMS() {return fResidualRange2RMS;};
        double ResidualRangeMinRMS() {return fResidualRangeMinRMS;};
        double ResidualRangeMaxRMS() {return fResidualRangeMaxRMS;};
        double ResidualRange1AngleRMS() {return fResidualRange1AngleRMS;};
        double ResidualRange2AngleRMS() {return fResidualRange2AngleRMS;};
        double ResidualRangeMinAngleRMS() {return fResidualRangeMinAngleRMS;};
        double ResidualRangeMaxAngleRMS() {return fResidualRangeMaxAngleRMS;};
};



#endif // TPC_SIMPLE_CALO_H