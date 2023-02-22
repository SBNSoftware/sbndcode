#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>

namespace sbnd::crt {

  enum CRTTagger CRTCommonUtils::GetTaggerEnum(const std::string tagger)
  {
    if      (tagger == "volTaggerBot_0"     ) return kBottomTagger;
    else if (tagger == "volTaggerSouth_0"   ) return kSouthTagger;
    else if (tagger == "volTaggerNorth_0"   ) return kNorthTagger;
    else if (tagger == "volTaggerWest_0"    ) return kWestTagger;
    else if (tagger == "volTaggerEast_0"    ) return kEastTagger;
    else if (tagger == "volTaggerTopLow_0"  ) return kTopLowTagger;
    else if (tagger == "volTaggerTopHigh_0" ) return kTopHighTagger;
    else {
      mf::LogWarning("CRTCommonUtils") << "CRT tagger unknown: " << tagger << std::endl;
      return kUndefinedTagger;
    }
  }

  std::string CRTCommonUtils::GetTaggerName(const CRTTagger tagger)
  {
    if      (tagger == kBottomTagger)  return "volTaggerBot_0";
    else if (tagger == kSouthTagger)   return "volTaggerSouth_0";
    else if (tagger == kNorthTagger)   return "volTaggerNorth_0";
    else if (tagger == kWestTagger)    return "volTaggerWest_0";
    else if (tagger == kEastTagger)    return "volTaggerEast_0";
    else if (tagger == kTopLowTagger)  return "volTaggerTopLow_0";
    else if (tagger == kTopHighTagger) return "volTaggerTopHigh_0";
    else {
      mf::LogWarning("CRTCommonUtils") << "CRT tagger unknown: " << tagger << std::endl;
      return "";
    }
  }

  enum CoordSet CRTCommonUtils::GetTaggerDefinedCoordinate(const CRTTagger tagger)
  {
    switch(tagger)
    {
      case kWestTagger:
      case kEastTagger:
        return kX;
      case kBottomTagger:
      case kTopLowTagger:
      case kTopHighTagger:
        return kY;
      case kSouthTagger:
      case kNorthTagger:
        return kZ;
      case kUndefinedTagger:
        return kUndefinedSet;
    }

    return kUndefinedSet;
  }
  
  enum CoordSet CRTCommonUtils::GetStripWidthGlobalCoordinate(const CRTTagger tagger, const uint16_t orientation)
  {
    if((tagger == kBottomTagger && orientation == 1) ||
       ((tagger == kTopLowTagger || tagger == kTopHighTagger) && orientation == 0) ||
       ((tagger == kSouthTagger || tagger == kNorthTagger) && orientation == 0))
      return kX;
    else if(((tagger == kSouthTagger || tagger == kNorthTagger) && orientation == 1) ||
            ((tagger == kWestTagger || tagger == kEastTagger) && orientation == 1))
      return kY;
    else if((tagger == kBottomTagger && orientation == 0) ||
            ((tagger == kTopLowTagger || tagger == kTopHighTagger) && orientation == 1) ||
            ((tagger == kWestTagger || tagger == kEastTagger) && orientation == 0))
      return kZ;

    return kUndefinedSet;
  }

  bool CRTCommonUtils::IsTopTagger(const CRTTagger tagger)
  {
    return tagger == kTopLowTagger || tagger == kTopHighTagger;
  }

  bool CRTCommonUtils::CoverTopTaggers(const CRTTagger tagger1, const CRTTagger tagger2, const CRTTagger tagger3)
  {
    return ( tagger1 == kTopLowTagger || tagger2 == kTopLowTagger || tagger3 == kTopLowTagger)
      && ( tagger1 == kTopHighTagger || tagger2 == kTopHighTagger || tagger3 == kTopHighTagger);
  }
}
