///////////////////////////////////////////////////////////////////////
/// File: ChargeDensity.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////
#ifndef SBND_CHARGEDENSITY_H
#define SBND_CHARGEDENSITY_H

#ifndef LAMBDAANA_LARSOFT
#define LAMBDAANA_LARSOFT 1
#endif

#include <vector>
#include <limits>
#include <filesystem>
#include <numeric>
#include <map>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "TGraphErrors.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TSystem.h"
#include <TROOT.h>
#include "TFile.h"
#include "TMVA/Reader.h"

#include "LeastSquares.h"
#include "VertexView.h"
#include "ChargeDensityPset.h"

#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/TPCSimpleHits.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/FRANS/FRANSObj.h"
#else
#include "TPCSimpleHits.h"
#include "FRANSObj.h"
#endif

#define DefaultMaxZSize 2000

class ChargeDensity{

  public:
    ChargeDensity(ChargeDensity const&) = delete;
    ChargeDensity(ChargeDensity&&) = delete;
    ChargeDensity & operator=(ChargeDensity const&) = delete;
    ChargeDensity & operator=(ChargeDensity &&) = delete;

    ChargeDensity(FRANSPsetType const& config, int view=-1);

    //void Fill(std::vector<art::Ptr<recob::Hit>> hitsVect, VertexView vertex);
    void Fill(std::vector<SHit> hitsVect, SVertex vertex);
    FRANSObj GetFRANSResult();
    void Display(TCanvas *c);

    // Access class members
    double Delta(){return fDelta;}
    double Eta(){return fEta;}
    double FitScore(){return fFitScore;}
    double Alpha(){return fAlpha;}
    double Omega(){return fOmega;}
    double Tau(){return fTau;}
    double Iota(){return fIota;}
    double Score(){return fScore;}
    double AverageHitChi2(){return fMeanChi2;}
    double HitDensity(){return fHitDensity;}
    double NHits(){return fNHits;}

  private:

    // configuration parameters
    FRANSPsetType const fFRANSPset;

    // view information
    unsigned int fView;
    unsigned int fTPC;

    // input vertex
    SVertex fVertex;

    // vectors storing the charge profile
    std::vector<double> fZ;
    std::vector<double> fZCounter;
    std::vector<double> fRho;
    std::vector<double> fZCum;
    std::vector<double> fZCumStart;
    std::vector<double> fZCumDer;
    std::vector<double> fZCumDerErr;

    double fNormUnAvSmooth;
    double fZetaNorm;
    double fMaxCumulative;

    // analysis parameters
    float fDelta;
    float fEta;
    float fFitScore;
    float fAlpha;
    float fOmega;
    float fTau;
    float fIota;
    float fScore;

    int fNHits;
    float fHitDensity;
    float fMeanChi2;

    float fGap;
    float fProtonKE;
    float fPionKE;

    TMVA::Reader fTMVAReader;

    void Reset();
    double GetDistance(int x, double y, int x0, double y0);
    int GetDistanceWires(int x, int x0);
    void UpdateMetrics();
    void ApplyExpoAvSmoothing(std::vector<double>& wf);
    void ApplyUnAvSmoothing(std::vector<double>& wf);
    void SlidingWindow(std::vector<double>& wf);
    void FillCumulative();
    void Rescale(std::vector<double>& wf, double kappa);

};


#endif
