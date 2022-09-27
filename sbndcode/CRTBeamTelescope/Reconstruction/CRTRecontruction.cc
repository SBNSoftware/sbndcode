#ifndef SBND_CRTRECONTRUCTION_CXX
#define SBND_CRTRECONTRUCTION_CXX

#include "CRTRecontruction.hh"

#include "sbndcode/CRTBeamTelescope/Reconstruction/LibCRTImport/CRTEvent.h"
#include "cetlib/search_path.h"
// #include "cetlib_except/exception.h"

namespace sbnd{
namespace crt{

  CRTRecontruction::CRTRecontruction():
    fMac5(0)
  {
    cet::search_path sp("FW_SEARCH_PATH");

    std::string fname_cabledelay;
    sp.find_file("SBND_CableDelay-V3.txt", fname_cabledelay);

    std::string fname_positions;
    sp.find_file("SBND_CRTpositionsSiPM-nogaps.txt", fname_positions);

    std::string fname_calibs;
    sp.find_file("SBND_ADCCalib-V4.txt", fname_calibs);

    _calibs = CRTCalibs(fname_cabledelay.c_str(), fname_positions.c_str(), fname_calibs.c_str());
  }


  CRTRecontruction::~CRTRecontruction() {}

  uint16_t CRTRecontruction::Mac5() const
  {
    return fMac5;
  }

} // namespace crt
} // namespace sbnd

#endif
