///////////////////////////////////////////////////////////////////////
/// File: ChargeDensity.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////
#ifndef SBND_FRANSObj_H
#define SBND_FRANSObj_H

#include <TTree.h>

#define defValue -9999

class FRANSObj {

private:

    // FRANS Variables
    int fView;
    double fScore;
    double fHitDensity;
    double fMeanChi2;
    double fDelta;
    double fEta;
    double fFitScore;
    double fAlpha;
    double fOmega;
    double fTau;
    double fIota;
    int fNOrigins;
    int fNOriginsM1;
    int fNOriginsM2;
    int fNOriginsM3;
    
public:

    FRANSObj(int view=defValue, double score=defValue, double hitDensity=defValue, double meanChi2=defValue,
             double delta=defValue, double eta=defValue, double fitScore=defValue,
             double alpha=defValue, double omega=defValue, double tau=defValue, double iota=defValue,
             int nOrigins=defValue, int nOriginsM1=defValue, int nOriginsM2=defValue, int nOriginsM3=defValue);

    // Getters
    int GetView() const;
    double GetHitDensity() const;
    double GetMeanChi2() const;
    double GetScore() const;
    double GetDelta() const;
    double GetEta() const;
    double GetFitScore() const;
    double GetAlpha() const;
    double GetOmega() const;
    double GetTau() const;
    double GetIota() const;
    int GetNOrigins() const;
    int GetNOriginsM1() const;
    int GetNOriginsM2() const;
    int GetNOriginsM3() const;
    

    // Setters
    void SetView(int view);
    void SetHitDensity(double hitDensity);
    void SetMeanChi2(double meanChi2);
    void SetScore(double score);
    void SetDelta(double delta);
    void SetEta(double eta);
    void SetFitScore(double fitScore);
    void SetAlpha(double alpha);
    void SetOmega(double omega);
    void SetTau(double tau);
    void SetIota(double iota);
    void SetNOrigins(int nOrigins);
    void SetNOriginsM1(int nOriginsM1);
    void SetNOriginsM2(int nOriginsM2);
    void SetNOriginsM3(int nOriginsM3);


};

#endif
