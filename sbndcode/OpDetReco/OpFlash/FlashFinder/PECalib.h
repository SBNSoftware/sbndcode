/**
 * \file PECalib.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class PECalib
 *
 * @author drinkingkazu
 */
/** \addtogroup UBFlashFinder 
 *
 *     @{*/

#ifndef PECALIB_H
#define PECALIB_H

#include "FlashFinderTypes.h"
#include "FlashFinderFMWKInterface.h"
#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>

namespace lightana{
/**
 \class PECalib
 User defined class PECalib ... these comments are used to generate
 doxygen documentation!
 */

  class PECalib {
    
  public:
    
    /// Default constructor
    PECalib();
        
    /// Default destructor
    ~PECalib(){}

    void Configure(const Config_t &pset);

    double Calibrate(const size_t opdet, const double area) const;

    protected:

    std::vector<double> _spe_area_gain_v;
    std::vector<double> _relative_qe_v;

  };

} 
#endif
/** @} */ // end of doxygen group
