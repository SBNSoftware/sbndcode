#ifndef CRTCOMMONUTILS_H_SEEN
#define CRTCOMMONUTILS_H_SEEN

///////////////////////////////////////////////
// CRTCommonUtils.h
//
// Common functions for CRT reconstruction
///////////////////////////////////////////////

#include <string>
#include "sbnobj/SBND/CRT/CRTEnums.hh"

namespace sbnd::crt {

  namespace CRTCommonUtils{

    // Returns the CRT plane index given the tagger name
    enum CRTTagger GetTaggerEnum(std::string tagger);

    // Returns the coordinate constrained by virtue of the tagger's position
    enum CoordSet GetTaggerDefinedCoordinate(const CRTTagger tagger);

    // Returns the coordinate direction of the strip's width given the tagger & local orientation
    enum CoordSet GetStripWidthGlobalCoordinate(const CRTTagger tagger, const uint16_t orientation);

    // Returns whether the tagger is horizontal (top or bottom), useful for track building
    bool IsHorizontalTagger(const CRTTagger tagger);
  }
}

#endif
