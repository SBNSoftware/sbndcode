///////////////////////////////////////////////////////////////////////
/// File: FRANSObj.cpp
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////

#include "FRANSObj.h"



FRANSObj::FRANSObj(int view, double score, double hitDensity, double meanChi2,
                   double delta, double eta, double fitScore,
                   double alpha, double omega, double tau, double iota,
                   int nOrigins, int nOriginsM1, int nOriginsM2, int nOriginsM3)
:   fView(view),
    fScore(score),
    fHitDensity(hitDensity),
    fMeanChi2(meanChi2),
    fDelta(delta),
    fEta(eta),
    fFitScore(fitScore),
    fAlpha(alpha),
    fOmega(omega),
    fTau(tau),
    fIota(iota),
    fNOrigins(nOrigins),
    fNOriginsM1(nOriginsM1),
    fNOriginsM2(nOriginsM2),
    fNOriginsM3(nOriginsM3)
{}

// Getters
int FRANSObj::GetView() const { return fView; }
double FRANSObj::GetHitDensity() const { return fHitDensity; }
double FRANSObj::GetMeanChi2() const { return fMeanChi2; }
double FRANSObj::GetScore() const { return fScore; }
double FRANSObj::GetDelta() const { return fDelta; }
double FRANSObj::GetEta() const { return fEta; }
double FRANSObj::GetFitScore() const { return fFitScore; }
double FRANSObj::GetAlpha() const { return fAlpha; }
double FRANSObj::GetOmega() const { return fOmega; }
double FRANSObj::GetTau() const { return fTau; }
double FRANSObj::GetIota() const { return fIota; }
int FRANSObj::GetNOrigins() const { return fNOrigins; }
int FRANSObj::GetNOriginsM1() const { return fNOriginsM1; }
int FRANSObj::GetNOriginsM2() const { return fNOriginsM2; }
int FRANSObj::GetNOriginsM3() const { return fNOriginsM3; }


// Setters
void FRANSObj::SetView(int view) { fView = view; }
void FRANSObj::SetHitDensity(double hitDensity) { fHitDensity = hitDensity; }
void FRANSObj::SetMeanChi2(double meanChi2) { fMeanChi2 = meanChi2; }
void FRANSObj::SetScore(double score) { fScore = score; }
void FRANSObj::SetDelta(double delta) { fDelta = delta; }
void FRANSObj::SetEta(double eta) { fEta = eta; }
void FRANSObj::SetFitScore(double fitScore) { fFitScore = fitScore; }
void FRANSObj::SetAlpha(double alpha) { fAlpha = alpha; }
void FRANSObj::SetOmega(double omega) { fOmega = omega; }
void FRANSObj::SetTau(double tau) { fTau = tau; }
void FRANSObj::SetIota(double iota) { fIota = iota; }
void FRANSObj::SetNOrigins(int nOrigins) { fNOrigins = nOrigins; }
void FRANSObj::SetNOriginsM1(int nOriginsM1) { fNOriginsM1 = nOriginsM1; }
void FRANSObj::SetNOriginsM2(int nOriginsM2) { fNOriginsM2 = nOriginsM2; }
void FRANSObj::SetNOriginsM3(int nOriginsM3) { fNOriginsM3 = nOriginsM3; }




