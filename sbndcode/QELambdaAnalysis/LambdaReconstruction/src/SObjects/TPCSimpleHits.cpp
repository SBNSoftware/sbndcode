////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleHits.cpp
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleHits.h"

// SPoint class functions
SPoint::SPoint(int x, int y) 
    : fX((float)x),
    fY((float)y)
{}

SPoint::SPoint(float x, float y) 
    : fX(x),
    fY(y)
{}

SPoint::SPoint(double x, double y) 
    : fX((float)x),
    fY((float)y)
{}

std::ostream& operator<<(std::ostream& out, SPoint const& p)
{
    out << " SPoint -- x= " << p.fX << " y=" <<p.fY <<std::endl;
    return out;
}

// SVertex class functions
SVertex::SVertex() 
  : fP( SPoint(-1, -1)),
  fView(""),
  fActive(false)
{
}

SVertex::SVertex(SPoint p, std::string view) 
  : fP(p),
  fView(view),
  fActive(true)
{
}

std::ostream& operator<<(std::ostream& out, SVertex const& v)
{
    out << " SVertex -- " << v.fP << std::endl;
    return out;
}

// SHit class functions
SHit::SHit(int id, float x, float y, float w, float integral, float st, float et, float chi2, int clusterId, float spx, float spy, float spz)
    : fId(id),
    fP( SPoint(x, y) ),
    fWidth(w),
    fStartT(st),
    fEndT(et),
    fIntegral(integral),
    fChi2(chi2),
    fClusterId(clusterId),
    fSPX(spx),
    fSPY(spy),
    fSPZ(spz),
    fXProj(-1000),
    fYProj(-1000),
    fCompactness(-1000),
    fConnectedness(-1000),
    fConnectednes1D(-1000)
{
}

SHit::SHit(float x, float y)
    : fId(-1),
    fP( SPoint(x, y) ),
    fWidth(-1),
    fStartT(-1),
    fEndT(-1),
    fIntegral(-1),
    fChi2(-1),
    fXProj(-1000),
    fYProj(-1000),
    fCompactness(-1000),
    fConnectedness(-1000),
    fConnectednes1D(-1000)
{
}


void SHit::SetHitConnectivity(float comp, float conn, float conn1d) {
    fCompactness = comp;
    fConnectedness = conn;
    fConnectednes1D = conn1d;
}

std::ostream& operator<<(std::ostream& out, SHit const& h)
{
    out << " SHit -- x=" << h.X() << " y=" << h.Y() << std::endl;
    return out;
}