#ifndef HITANALYSISALG_H
#define HITANALYSISALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       HitAnalysisAlg
// Module Type: producer
// File:        HitAnalysisAlg.h
//
//              The intent of this module is to provide methods for
//              "analyzing" hits on waveforms
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//
// Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"

#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/Hit.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace HitAnalysis
{
    
// The following typedefs will, obviously, be useful
using HitVec           = std::vector<recob::Hit>;
using PlaneHitMap      = std::map<size_t,HitVec>;
using TrackPlaneHitMap = std::map<int,PlaneHitMap>;
    
class HitAnalysisAlg
{
public:

    // Copnstructors, destructor.
    HitAnalysisAlg(fhicl::ParameterSet const & pset);
    ~HitAnalysisAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void setup(const geo::GeometryCore&, TDirectory*);
    void endJob(int numEvents);
    
    void fillHistograms(const TrackPlaneHitMap&) const;
    void fillHistograms(const HitVec&)           const;
    
private:

    // Fcl parameters.
    std::string fLocalDirName;     ///< Fraction for truncated mean
    TDirectory* fRootDirectory;
    
    // Pointers to the histograms we'll create.
    std::unique_ptr<TH1D>     fHitsByWire[3];
    std::unique_ptr<TH1D>     fDriftTimes[3];
    std::unique_ptr<TH1D>     fHitsByTime[3];
    std::unique_ptr<TH1D>     fPulseHeight[3];
    std::unique_ptr<TH1D>     fPulseHeightSingle[3];
    std::unique_ptr<TH1D>     fPulseHeightMulti[3];
    std::unique_ptr<TH1D>     fChi2DOF[3];
    std::unique_ptr<TH1D>     fNumDegFree[3];
    std::unique_ptr<TH1D>     fChi2DOFSingle[3];
    std::unique_ptr<TH1D>     fHitMult[3];
    std::unique_ptr<TH1D>     fHitCharge[3];
    std::unique_ptr<TH1D>     fFitWidth[3];
    std::unique_ptr<TH1D>     fHitSumADC[3];
    std::unique_ptr<TH2D>     fNDFVsChi2[3];
    std::unique_ptr<TH2D>     fPulseHVsWidth[3];
    std::unique_ptr<TH2D>     fPulseHVsCharge[3];
    std::unique_ptr<TH1D>     fBadWPulseHeight;
    std::unique_ptr<TH2D>     fBadWPulseHVsWidth;
    std::unique_ptr<TH1D>     fBadWHitsByWire;
    std::unique_ptr<TProfile> fPulseHVsHitNo[3];
    std::unique_ptr<TProfile> fChargeVsHitNo[3];
    std::unique_ptr<TProfile> fChargeVsHitNoS[3];
    
    std::unique_ptr<TH2D>     fSPHvsIdx[3];
    std::unique_ptr<TH2D>     fSWidVsIdx[3];
    std::unique_ptr<TH2D>     f1PPHvsWid[3];
    std::unique_ptr<TH2D>     fSPPHvsWid[3];
    std::unique_ptr<TH2D>     fSOPHvsWid[3];
    std::unique_ptr<TH2D>     fPHRatVsIdx[3];
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
};

} // end of namespace caldata

#endif
