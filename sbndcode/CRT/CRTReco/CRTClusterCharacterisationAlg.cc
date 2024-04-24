#include "CRTClusterCharacterisationAlg.h"

namespace sbnd::crt {
  
  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(const fhicl::ParameterSet& pset)
    : fCRTGeoAlg(pset.get<fhicl::ParameterSet>("GeoAlg", fhicl::ParameterSet()))
    , fUseT1(pset.get<bool>("UseT1"))
    , fTimeOffset(pset.get<double>("TimeOffset"))
    , fOverlapBuffer(pset.get<double>("OverlapBuffer"))
    , fPEAttenuation(pset.get<double>("PEAttenuation"))
    , fPropDelay(pset.get<double>("PropDelay"))
    , fTimeWalkNorm(pset.get<double>("TimeWalkNorm"))
    , fTimeWalkScale(pset.get<double>("TimeWalkScale"))
  {
  }

  CRTClusterCharacterisationAlg::CRTClusterCharacterisationAlg(){}
  
  CRTClusterCharacterisationAlg::~CRTClusterCharacterisationAlg(){}

  CRTSpacePoint CRTClusterCharacterisationAlg::CharacteriseSingleHitCluster(const art::Ptr<CRTCluster> &cluster, const art::Ptr<CRTStripHit> &stripHit)
  {
    const std::array<double, 6> hitPos = fCRTGeoAlg.StripHit3DPos(stripHit->Channel(), stripHit->Pos(), stripHit->Error());

    const double pe = ADCToPE(stripHit->Channel(), stripHit->ADC1(), stripHit->ADC2());
    
    geo::Point_t pos, err;
    CentralPosition(hitPos, pos, err);

    return CRTSpacePoint(pos, err, pe, stripHit->Ts1() + fTimeOffset, 0., false);
  }

  bool CRTClusterCharacterisationAlg::CharacteriseDoubleHitCluster(const art::Ptr<CRTCluster> &cluster, const std::vector<art::Ptr<CRTStripHit>> &stripHits, CRTSpacePoint &spacepoint)
  {
    const art::Ptr<CRTStripHit> &hit0 = stripHits[0];
    const art::Ptr<CRTStripHit> &hit1 = stripHits[1];

    return TwoHitSpacePoint(hit0, hit1, spacepoint);
  }

  bool CRTClusterCharacterisationAlg::TwoHitSpacePoint(const art::Ptr<CRTStripHit> hit0, const art::Ptr<CRTStripHit> hit1, CRTSpacePoint &spacepoint)
  {
    const bool threeD = fCRTGeoAlg.ChannelToOrientation(hit0->Channel()) != fCRTGeoAlg.ChannelToOrientation(hit1->Channel());

    if(threeD)
      {
        if(fCRTGeoAlg.CheckOverlap(hit0->Channel(), hit1->Channel(), fOverlapBuffer))
          {
            std::array<double, 6> overlap = FindOverlap(hit0, hit1);
            
            geo::Point_t pos, err;
            CentralPosition(overlap, pos, err);

            const double pe   = ReconstructPE(hit0, hit1, pos);
            double time, etime;
            CorrectTime(hit0, hit1, pos, time, etime);

            spacepoint = CRTSpacePoint(pos, err, pe, time + fTimeOffset, etime, true);
            return true;
          }
        return false;
      }
    else
      {
        if(fCRTGeoAlg.AdjacentStrips(hit0->Channel(), hit1->Channel(), fOverlapBuffer))
          {
            std::array<double, 6> overlap = FindAdjacentPosition(hit0, hit1);

            geo::Point_t pos, err;
            CentralPosition(overlap, pos, err);

            const double pe   = ADCToPE(hit0->Channel(), hit0->ADC1(), hit0->ADC2()) + ADCToPE(hit1->Channel(), hit1->ADC1(), hit1->ADC2());
            const double time  = (hit0->Ts1() + hit1->Ts1()) / 2.;
            const double etime = std::abs((double)hit0->Ts1() - (double)hit1->Ts1()) / 2.;

            spacepoint = CRTSpacePoint(pos, err, pe, time + fTimeOffset, etime, false);
            return true;
          }
        return false;
      }
  }

  bool CRTClusterCharacterisationAlg::CharacteriseMultiHitCluster(const art::Ptr<CRTCluster> &cluster, const std::vector<art::Ptr<CRTStripHit>> &stripHits, CRTSpacePoint &spacepoint)
  {
    std::vector<CRTSpacePoint> spacepoints, complete_spacepoints;

    for(unsigned i = 0; i < stripHits.size(); ++i)
      {
        const art::Ptr<CRTStripHit> hit0 = stripHits[i];

        for(unsigned ii = i + 1; ii < stripHits.size(); ++ii)
          {
            const art::Ptr<CRTStripHit> hit1 = stripHits[ii];

            CRTSpacePoint sp;
            if(TwoHitSpacePoint(hit0, hit1, sp))
              {
                spacepoints.push_back(sp);
        
                if(sp.Complete()) 
                  complete_spacepoints.push_back(sp);
              }
          }
      }

    if(complete_spacepoints.size() == 0)
      return false;
    
    double pe = 0.;
    std::vector<double> times;

    for(auto const &sp : complete_spacepoints)
      {
        pe += sp.PE();
        times.push_back(sp.Time());
      }

    geo::Point_t pos, err;
    AggregatePositions(complete_spacepoints, pos, err);

    double time, etime;
    TimeErrorCalculator(times, time, etime);
    
    spacepoint = CRTSpacePoint(pos, err, pe, time, etime, true);
    return true;
  }

  double CRTClusterCharacterisationAlg::ADCToPE(const uint16_t channel, const uint16_t adc1, const uint16_t adc2)
  {
    return ADCToPE(channel, adc1) + ADCToPE(channel+1, adc2);
  }

  double CRTClusterCharacterisationAlg::ADCToPE(const uint16_t channel, const uint16_t adc)
  {
    return fCRTGeoAlg.GetSiPM(channel).gain * adc;
  }

  std::array<double, 6> CRTClusterCharacterisationAlg::FindOverlap(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1)
  {
    const std::array<double, 6> hit0pos = fCRTGeoAlg.StripHit3DPos(hit0->Channel(), hit0->Pos(), hit0->Error());
    const std::array<double, 6> hit1pos = fCRTGeoAlg.StripHit3DPos(hit1->Channel(), hit1->Pos(), hit1->Error());

    std::array<double, 6> overlap({std::max(hit0pos[0], hit1pos[0]),
                                   std::min(hit0pos[1], hit1pos[1]),
                                   std::max(hit0pos[2], hit1pos[2]),
                                   std::min(hit0pos[3], hit1pos[3]),
                                   std::max(hit0pos[4], hit1pos[4]),
                                   std::min(hit0pos[5], hit1pos[5])});

    return overlap;
  }

  std::array<double, 6> CRTClusterCharacterisationAlg::FindAdjacentPosition(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1)
  {
    const std::array<double, 6> hit0pos = fCRTGeoAlg.StripHit3DPos(hit0->Channel(), hit0->Pos(), hit0->Error());
    const std::array<double, 6> hit1pos = fCRTGeoAlg.StripHit3DPos(hit1->Channel(), hit1->Pos(), hit1->Error());

    std::array<double, 6> overlap({std::min(hit0pos[0], hit1pos[0]),
                                   std::max(hit0pos[1], hit1pos[1]),
                                   std::min(hit0pos[2], hit1pos[2]),
                                   std::max(hit0pos[3], hit1pos[3]),
                                   std::min(hit0pos[4], hit1pos[4]),
                                   std::max(hit0pos[5], hit1pos[5])});

    return overlap;
  }

  void CRTClusterCharacterisationAlg::CentralPosition(const std::array<double, 6> overlap, 
                                                      geo::Point_t &pos, geo::Point_t &err)
  {
    pos = geo::Point_t((overlap[0] + overlap[1])/2.,
                       (overlap[2] + overlap[3])/2.,
                       (overlap[4] + overlap[5])/2.);
    
    err = geo::Point_t(std::abs((overlap[0] - overlap[1])/2.),
                       std::abs((overlap[2] - overlap[3])/2.),
                       std::abs((overlap[4] - overlap[5])/2.));
  }

  double CRTClusterCharacterisationAlg::ReconstructPE(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1, const geo::Point_t &pos)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0->Channel());
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1->Channel());

    return ReconstructPE(hit0, dist0) + ReconstructPE(hit1, dist1);
  }

  double CRTClusterCharacterisationAlg::ReconstructPE(const art::Ptr<CRTStripHit> &hit, const double dist)
  {
    const double pe         = ADCToPE(hit->Channel(), hit->ADC1(), hit->ADC2());
    const double correction = std::pow(dist - fPEAttenuation, 2) / std::pow(fPEAttenuation, 2);

    return pe * correction;
  }

  void CRTClusterCharacterisationAlg::CorrectTime(const art::Ptr<CRTStripHit> &hit0, const art::Ptr<CRTStripHit> &hit1, const geo::Point_t &pos,
                                                  double &time, double &etime)
  {
    const double dist0 = fCRTGeoAlg.DistanceDownStrip(pos, hit0->Channel());
    const double dist1 = fCRTGeoAlg.DistanceDownStrip(pos, hit1->Channel());

    const double pe0 = ReconstructPE(hit0, dist0);
    const double pe1 = ReconstructPE(hit1, dist1);

    const double corr0 = TimingCorrectionOffset(dist0, pe0);
    const double corr1 = TimingCorrectionOffset(dist1, pe1);

    if(fUseT1)
      {
        time  = ((double)hit0->Ts1() - corr0 + (double)hit1->Ts1() - corr1) / 2.;
        etime = std::abs(((double)hit0->Ts1() - corr0) - ((double)hit1->Ts1() - corr1)) / 2.;
      }
    else
      {
        time  = ((double)hit0->Ts0() - corr0 + (double)hit1->Ts0() - corr1) / 2.;
        etime = std::abs(((double)hit0->Ts0() - corr0) - ((double)hit1->Ts0() - corr1)) / 2.;
      }
  }

  double CRTClusterCharacterisationAlg::TimingCorrectionOffset(const double &dist, const double &pe)
  {

    double t_TimeWalk  = fTimeWalkNorm * std::exp(- fTimeWalkScale * pe);
    double t_PropDelay = fPropDelay * dist;
    return t_PropDelay + t_TimeWalk;
  }

  void CRTClusterCharacterisationAlg::AggregatePositions(const std::vector<CRTSpacePoint> &complete_spacepoints, geo::Point_t &pos, geo::Point_t &err)
  {
    double sum_x = 0., sum_y = 0., sum_z = 0.;
    const unsigned n_sp = complete_spacepoints.size();

    for(auto const &sp : complete_spacepoints)
      {
        const geo::Point_t sp_pos = sp.Pos();
        sum_x += sp_pos.X();
        sum_y += sp_pos.Y();
        sum_z += sp_pos.Z();
      }

    pos = geo::Point_t(sum_x / n_sp, sum_y / n_sp, sum_z / n_sp);

    double sum_var_x = 0., sum_var_y = 0., sum_var_z = 0.;

    for(auto const &sp : complete_spacepoints)
      {
        const geo::Point_t sp_pos = sp.Pos();
        sum_var_x += std::pow(sp_pos.X() - pos.X(), 2);
        sum_var_y += std::pow(sp_pos.Y() - pos.Y(), 2);
        sum_var_z += std::pow(sp_pos.Z() - pos.Z(), 2);
      }

    err = geo::Point_t(std::sqrt(sum_var_x / n_sp), std::sqrt(sum_var_y / n_sp), std::sqrt(sum_var_z / n_sp));

    if(err.X() < std::numeric_limits<double>::epsilon())
      err.SetX(complete_spacepoints[0].XErr());
    if(err.Y() < std::numeric_limits<double>::epsilon())
      err.SetY(complete_spacepoints[0].YErr());
    if(err.Z() < std::numeric_limits<double>::epsilon())
      err.SetZ(complete_spacepoints[0].ZErr());
  }

  void CRTClusterCharacterisationAlg::TimeErrorCalculator(const std::vector<double> &times, double &mean, double &err)
  {
    double sum = 0.;
    for(auto const &time : times)
      sum += time;

    mean = sum / times.size();

    double summed_var = 0.;
    for(auto const &time : times)
      summed_var += std::pow((time - mean), 2);

    err = std::sqrt(summed_var / times.size());
  }
}
