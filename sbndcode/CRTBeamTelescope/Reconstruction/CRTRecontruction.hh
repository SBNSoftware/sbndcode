/**
 * \brief Class for CRT reconstruction
 *
 * \author Marco Del Tutto
 */

#ifndef SBND_CRTRECONTRUCTION_HH
#define SBND_CRTRECONTRUCTION_HH

#include <stdint.h>
#include <vector>
#include <utility>
#include <array>

#include "LibCRTImport/CRTEvent.h"


namespace sbnd::crt {

  class CRTRecontruction {

    uint16_t fMac5; ///< ID of the FEB
    CRTCalibs _calibs;

   public:

    CRTRecontruction();


    virtual ~CRTRecontruction();

    uint16_t Mac5() const;

  };

} // namespace sbnd::crt

#endif