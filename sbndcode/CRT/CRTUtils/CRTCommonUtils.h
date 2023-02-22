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

    // Returns the CRT tagger index given the tagger name
    enum CRTTagger GetTaggerEnum(const std::string tagger);

    // Returns the CRT tagger name given the tagger index
    std::string GetTaggerName(const CRTTagger tagger);

    // Returns the coordinate constrained by virtue of the tagger's position
    enum CoordSet GetTaggerDefinedCoordinate(const CRTTagger tagger);

    // Returns the coordinate direction of the strip's width given the tagger & local orientation
    enum CoordSet GetStripWidthGlobalCoordinate(const CRTTagger tagger, const uint16_t orientation);

    // Returns whether the tagger is one of the top taggers, useful for track building
    bool IsTopTagger(const CRTTagger tagger);

    // Returns whether a set of three taggers includes both of the top taggers, useful for track building
    bool CoverTopTaggers(const CRTTagger tagger1, const CRTTagger tagger2, const CRTTagger tagger3);
  }
}

#endif
