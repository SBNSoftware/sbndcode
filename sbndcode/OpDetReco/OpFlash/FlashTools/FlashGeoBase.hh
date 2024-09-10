///////////////////////////////////////////////////////////////////////
/// File: FlashGeoBase.h
///
/// Interfacce class for a tool to calculate the recob::OpFlash
/// Y and Z coordinates (PDS plane)
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#ifndef SBND_FLASHGEOBASE_H
#define SBND_FLASHGEOBASE_H

#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderFMWKInterface.h"

namespace lightana
{
  class FlashGeoBase{

  public:
    // Default destructor
    virtual ~FlashGeoBase() noexcept = default;

    // Method to calculate flash geometric properties
    virtual void GetFlashLocation(std::vector<double> pePerOpChannel,
                                  double& Ycenter, double& Zcenter,
                                  double& Ywidth, double& Zwidth) = 0 ;

  private:

  };
}

#endif
