// SBND optical path tool

#ifndef SBNDOpticalPath_H
#define SBNDOpticalPath_H

#include "art/Utilities/ToolMacros.h" 
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include <iostream>

namespace phot {
    class SBNDOpticalPath : public phot::OpticalPath {
    public:
        explicit SBNDOpticalPath(fhicl::ParameterSet const& ps) {};
        ~SBNDOpticalPath() noexcept override = default;

        const bool isOpDetVisible(geo::Point_t const& ScintPoint, geo::Point_t const& OpDetPoint) override {           
            // special case for SBND 
            // check x coordinate has same sign or is close to zero
            if ((ScintPoint.X() < 0.) != (OpDetPoint.X() < 0.)) return false;            
            else return true; 
        }
    };
}

#endif