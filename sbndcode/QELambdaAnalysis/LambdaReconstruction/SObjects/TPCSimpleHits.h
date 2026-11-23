////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////


#ifndef TPC_SIMPLE_HITS_H
#define TPC_SIMPLE_HITS_H


#include <iostream>
#include <string>
#include <iosfwd>
#include <vector>
#include <cmath>

class SPoint {
    private:
        float fX;
        float fY;

    public:
        SPoint(int x=0, int y=0);
        SPoint(double x, double y);
        SPoint(float x, float y);
        float X() const {return fX;}
        float Y() const {return fY;}
        void AssignX(float x) { fX = x;}
        void AssignY(float y) { fY = y;}

        friend std::ostream& operator<<(std::ostream& out, SPoint const& p);
};

class SVertex {
    private:
        SPoint fP;
        std::string fView;
        bool fActive;

    public:
        SVertex();
        SVertex(SPoint p, std::string view="");
        
        SPoint Point(){return fP;}
        std::string View() const {return fView;}
        bool IsActive() const {return fActive;}

        double X() const {return fP.X();}
        double Y() const 
        {return fP.Y();}

        friend std::ostream& operator<<(std::ostream& out, SVertex const& v);
};

class SHit {
    private:
        int fId;
        SPoint fP;
        float fWidth;
        float fStartT;
        float fEndT;
        float fIntegral;
        float fChi2;
        int fClusterId;
        float fSPX;
        float fSPY;
        float fSPZ;
        float fXProj;
        float fYProj;
        float fCompactness;
        float fConnectedness;
        float fConnectednes1D;

    public:
        // Constructor
        SHit(int id=-1, float x=0, float y=0, float w=0, float integral=0, float st=0, float et=0, float chi2=0, int clusterID=-1, float spx=0, float spy=0, float spz=0);
        SHit(float x, float y);

        // Getters
        int Id(){return fId;}
        float X() const {return fP.X();}
        float Y() const {return fP.Y();}
        float Width() const {return fWidth;}
        float StartT() const {return fStartT;}
        float EndT() const {return fEndT;}
        float Integral() const {return fIntegral;}
        float Chi2() const {return fChi2;}
        int ClusterId() const {return fClusterId;}
        float SPX() const {return fSPX;}
        float SPY() const {return fSPY;}
        float SPZ() const {return fSPZ;}
        float XProj() const {return fXProj;}
        float YProj() const {return fYProj;}
        float Compactness() const {return fCompactness;}
        float Connectednes() const {return fConnectedness;}
        float Connectednes1D() const {return fConnectednes1D;}
        SPoint GetPoint() const {return fP;}

        // Setters
        void SetX(float x) { fP.AssignX(x); };
        void SetY(float y) { fP.AssignY(y); };
        void SetWidth(float w) { fWidth = w; };


        // Set the hit connectivity
        void SetHitConnectivity(float comp, float conn, float conn1d);

        // Overload cout
        friend std::ostream& operator<<(std::ostream& out, SHit const& h);
};


#endif // TPC_SIMPLE_HITS_H