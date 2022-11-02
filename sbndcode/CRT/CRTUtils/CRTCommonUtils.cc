#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>

namespace sbnd{

enum CRTPlane CRTCommonUtils::GetPlaneIndex(std::string tagger) {

  if      (tagger == "volTaggerBot_0"      ) return kCRTBot;
  else if (tagger == "volTESTTaggerSouth_0") return kCRTFaceSouth;
  else if (tagger == "volTESTTaggerNorth_0") return kCRTFaceNorth;
  else {
    mf::LogWarning("CRTCommonUtils") << "CRT tagger unkown: " << tagger << std::endl;
    return kCRTNotDefined;
  }
}
}

