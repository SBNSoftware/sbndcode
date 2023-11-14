#include "SecondShowerFinderAlg.h"

SecondShowerFinderAlg::SecondShowerFinderAlg()
{

}

SecondShowerFinderAlg::SecondShowerFinderAlg(fhicl::ParameterSet const& p)
{

}

bool FindSecondShower(const art::Ptr<recob::Shower> &shower, const std::vector<art::Ptr<recob::Hit>> &hits)
{
  return false;
}
