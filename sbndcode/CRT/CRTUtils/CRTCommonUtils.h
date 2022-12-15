#ifndef CRTCOMMONUTILS_H_SEEN
#define CRTCOMMONUTILS_H_SEEN

///////////////////////////////////////////////
// CRTCommonUtils.h
//
// Common functions for CRT reconstruction
///////////////////////////////////////////////

#include <string>
#include "sbnobj/SBND/CRT/CRTCluster.hh"

namespace sbnd::crt {

  namespace CRTCommonUtils{

    // Returns the CRT plane index given the tagger name
    enum CRTTagger GetTaggerEnum(std::string tagger);
  }
}

#endif
