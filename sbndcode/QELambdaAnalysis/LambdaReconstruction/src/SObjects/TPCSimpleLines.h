////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineObjects.h
//
// \brief Definition of TPCSimpleLines
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_LINES_H
#define TPC_SIMPLE_LINES_H

#include <cmath>
#include "TPCSimpleHits.h" 

class LineEquation {
    private:
        float fM;
        float fN;
        float fGoodness;

    public:
        LineEquation(float slope=0, float intercept=0, float goodness=0);

        void SetGoodness(float goodness);
        
        float Slope() {return fM;}
        float Intercept() {return fN;}
        float Goodness() {return fGoodness;}

        SPoint GetLineClosestPoint(double a, double b, double c, SPoint p);
        SPoint GetLineClosestPoint(SPoint P);
        SPoint GetLineClosestPoint(SHit P);
        float GetDistance(SPoint P);
        float Evaluate(SPoint p);
        float EvaluateX(double x);
};


class HoughLine {
    private:
        LineEquation fEquation;
        float fScore;
        int fNHits;
    public:
        HoughLine(LineEquation line = LineEquation(0, 0), float score = -1, int nHits = 0);
        
        void SetLineEquation(LineEquation eq){ fEquation = eq;}
        void SetScore(float sc){ fScore = sc;}
        void SetNHits(float nhits){ fNHits = nhits;}

        float Score(){return fScore;}
        float NHits(){return fNHits;}
        LineEquation GetLineEquation(){return fEquation;}
};


#endif // TPC_SIMPLE_LINES_H