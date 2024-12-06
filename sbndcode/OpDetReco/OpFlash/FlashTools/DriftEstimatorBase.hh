///////////////////////////////////////////////////////////////////////
/// File: DriftEstimatorBase.hh
///
/// Interface class for drift coordinate estimation from
/// reco objects
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#ifndef SBND_DRIFTESTIMATORBASE_H
#define SBND_DRIFTESTIMATORBASE_H

namespace lightana
{
  class DriftEstimatorBase{

  public:
    // Default destructor
    virtual ~DriftEstimatorBase() noexcept = default;

    // Method giving the estimated drift coordinate
    virtual double GetDriftPosition(std::vector<double> PE_v) = 0;

    // Method giving the photon propagation
    virtual double GetPropagationTime(double drift) = 0;


  };
}

#endif
