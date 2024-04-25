#ifndef CRTHITRECOALG_H_SEEN
#define CRTHITRECOALG_H_SEEN

#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Persistency/Common/Ptr.h" 

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"

namespace sbnd{

  struct CRTStripHit {
     CRTStripHit(uint32_t _channel, uint32_t _t0, uint32_t _t1, uint32_t _s, double _x,
                 double _ex, uint16_t _adc1, uint16_t _adc2, uint16_t _pedestal1, uint16_t _pedestal2, uint8_t _febdataindex 
                 /*,std::array<float, 2> _channels_below_threshold*/) 
     : channel      (_channel)
     , t0           (_t0)
     , t1           (_t1)
     , s            (_s)
     , x            (_x)
     , ex           (_ex)
     , adc1         (_adc1)
     , adc2         (_adc2)
     , pedestal1    (_pedestal1)
     , pedestal2    (_pedestal2)
     , febdataindex (_febdataindex)
     //, channels_below_threshold (_channels_below_threshold)
    {};

    uint32_t channel;      // Channel ID for 1st SiPM
    uint32_t t0;           // T0 counter [ns]
    uint32_t t1;           // T1 counter [ns]
    uint32_t s;            // Unixtime of event [s]
    double   x;            // Lateral position within strip [cm]
    double   ex;           // Error on lateral position [cm]
    uint16_t adc1;         // ADC 1st SiPM
    uint16_t adc2;         // ADC 2nd SiPM
    uint16_t pedestal1;    // pedestal 1st SiPM
    uint16_t pedestal2;    // pedestal 2nd SiPM
    uint8_t  febdataindex; // Index of FEBData used to make strip hit
    //std::array<float, 2> channels_below_threshold; /*format [total # strips, total # of ADC]*/
  };


  class CRTHitRecoAlg {
  public:

    CRTHitRecoAlg(const fhicl::ParameterSet &p);

    CRTHitRecoAlg();

    ~CRTHitRecoAlg();

    std::map<std::string, std::vector<std::vector<CRTStripHit>>>
      ProduceStripHits(const std::vector<art::Ptr<sbnd::crt::FEBData>> &FEBdataVec);
    
    std::vector<sbn::crt::CRTHit> ProduceCRTHits(const std::map<std::string, std::vector<std::vector<CRTStripHit>>> &taggerStripHits);

    std::vector<std::pair<std::pair<unsigned, unsigned>, sbn::crt::CRTHit>> ProduceCRTHitCandidates(const std::string &tagger, 
                                                                                                    const std::vector<CRTStripHit> &hitsOrien0,
                                                                                                    const std::vector<CRTStripHit> &hitsOrien1);

    std::vector<double> FindOverlap(const CRTStripHit &hit0,   const CRTStripHit &hit1,
                                    const CRTStripGeo &strip0, const CRTStripGeo &strip1);

    void CentralPosition(const std::vector<double> overlap, TVector3 &pos, TVector3 &err);

    //void ReconstructPE(const TVector3 &pos, const CRTStripHit &hit0,
                       //const CRTStripHit &hit1, double &pe0, double &pe1);
    void ReconstructPE(const TVector3 &pos, const CRTStripHit &hit0, const CRTStripHit &hit1, 
                       double &pe0, double &pe1, 
                       double &uncorrected_pe, double &uncorrected_pe1, 
                       double &dist0, double &dist1);

    //double ReconstructPE(const double &dist, const CRTStripHit &hit);
    double ReconstructPE(const double &dist, const CRTStripHit &hit, double &uncorrected_pe);

    double ADCtoPE(const uint16_t &adc);

    void CorrectADC(const TVector3 &pos, const CRTStripHit &hit0,
                    const CRTStripHit &hit1, std::array<uint16_t, 4> &corr_adc);

    void CorrectTimings(const TVector3 &pos, const CRTStripHit &hit0,
                        const CRTStripHit &hit1, const double &pe0, const double &pe1,
                        uint32_t &t0, uint32_t &t1, double &diff);

    uint32_t TimingCorrectionOffset(const double &dist, const double &pe);

  private:
    
    CRTGeoAlg fCRTGeoAlg;
    uint16_t  fADCThreshold;
    uint16_t  fADCSaturation;
    double    fQPedestal;
    double    fQSlope;
    double    fPEAttenuation;
    double    fPropDelay;
    double    fTimeWalkNorm;
    double    fTimeWalkShift;
    double    fTimeWalkSigma;
    double    fTimeWalkOffset;
    uint32_t  fHitCoincidenceRequirement;
    uint32_t  fT1Offset;
  };

}

#endif