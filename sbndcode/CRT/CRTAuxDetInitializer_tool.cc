#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"

#include "art/Utilities/ToolMacros.h"

namespace {
  class CRTAuxDetInitializer : public geo::AuxDetInitializer {
  public:
    explicit CRTAuxDetInitializer(fhicl::ParameterSet const&) {}

  private:
    geo::AuxDetReadoutInitializers
    initialize(std::vector<geo::AuxDetGeo> const& adgeo) const override
    {
      geo::AuxDetReadoutInitializers result;
      // Map the AuxDetGeo names to their position in the sorted vector
      //
      // Each tagger is composed of scintillator modules, composed of 16 strips.
      // In the geometry, CRTStripArrays are AuxDets and CRTStrips are the
      // AuxDetSensitives. Each strip has two SiPM channels, one per optical
      // fiber (not in the geometry).
      //
      // 2x TaggerTop:  (5x2  +5x2)x16 = 320 strips,  640 channels
      // 2x TaggerSide: (4x2  +5x2)x16 = 288 strips,  576 channels
      // 1x TaggerFace: (4x2  +4x2)x16 = 256 strips,  512 channels
      // 1x TaggerFace: (4x2-1+4x2)x16 = 240 strips,  480 channels
      // 1x TaggerBot:           34x16 = 544 strips, 1088 channels

      for (size_t a=0; a<adgeo.size(); a++){
        std::string volName(adgeo[a].TotalVolume()->GetName());

        long unsigned int number_scintillating_strips = 0;

        if (strncmp(((adgeo[a].TotalVolume())->GetShape())->GetName(), "CRTstripMINOSArray", 18) == 0) {
          number_scintillating_strips = 20;    //To account for the MINOS modules.
        }
        else {number_scintillating_strips = 16;}

        size_t nsv = adgeo[a].NSensitiveVolume();
        if (nsv != number_scintillating_strips) {
          throw cet::exception("CRTChannelMap")
            << "Wrong number of sensitive volumes for CRT volume "
            << volName << " (got " << nsv << ", expected 16)" << std::endl;
        }

        result.ADGeoToName[a] = volName;
        result.NameToADGeo[volName] = a;

        if (volName.find("CRTStripArray") != std::string::npos) {
          for (size_t svID=0; svID<number_scintillating_strips; svID++) {
            for (size_t ich=0; ich<2; ich++) {
              size_t chID = 2 * svID + ich;
              result.ADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
            }
          }
        }
      }
      return result;
    }
  };
}

DEFINE_ART_CLASS_TOOL(CRTAuxDetInitializer)
