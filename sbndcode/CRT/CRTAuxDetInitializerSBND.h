#ifndef SBND_CRT_CRTAUXDETINITIALIZERSBND_H
#define SBND_CRT_CRTAUXDETINITIALIZERSBND_H

#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"

namespace sbnd::crt {

  class CRTAuxDetInitializerSBND : public geo::AuxDetInitializer {
  public:
    explicit CRTAuxDetInitializerSBND(fhicl::ParameterSet const&);

  private:
    geo::AuxDetReadoutInitializers
    initialize(std::vector<geo::AuxDetGeo> const& adgeo) const override;
  };

} // namespace sbnd::crt

#endif // SBND_CRT_CRTAUXDETINITIALIZERSBND_H
