#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>

namespace sbnd::crt {

  enum CRTTagger CRTCommonUtils::GetTaggerEnum(std::string tagger) {

    if      (tagger == "volTaggerBot_0"     ) return kBottomTagger;
    else if (tagger == "volTaggerSouth_0"   ) return kSouthTagger;
    else if (tagger == "volTaggerNorth_0"   ) return kNorthTagger;
    else if (tagger == "volTaggerWest_0"    ) return kWestTagger;
    else if (tagger == "volTaggerEast_0"    ) return kEastTagger;
    else if (tagger == "volTaggerTopLow_0"  ) return kTopLowTagger;
    else if (tagger == "volTaggerTopHigh_0" ) return kTopHighTagger;
    else {
      mf::LogWarning("CRTCommonUtils") << "CRT tagger unkown: " << tagger << std::endl;
      return kUndefinedTagger;
    }
  }
}
