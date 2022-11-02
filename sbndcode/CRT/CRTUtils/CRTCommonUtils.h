#ifndef CRTCOMMONUTILS_H_SEEN
#define CRTCOMMONUTILS_H_SEEN

///////////////////////////////////////////////
// CRTCommonUtils.h
//
// Common functions for CRT reconstruction
///////////////////////////////////////////////

#include <string>

namespace sbnd {
  enum CRTPlane {
    kCRTNotDefined = -1,   ///< Not defined
    kCRTBot = 0,           ///< Bottom
    kCRTFaceSouth = 1,     ///< Face South (Front)
    kCRTFaceNorth,         ///< Face North (Back)
    kCRTPosMax
  };
}

namespace sbnd{
namespace CRTCommonUtils{

  // Returns the CRT plane index given the tagger name
  enum ::sbnd::CRTPlane GetPlaneIndex(std::string tagger);
}
}

#endif
