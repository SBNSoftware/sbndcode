#ifndef LCVN_SBNDPIXELMAP_H
#define LCVN_SBNDPIXELMAP_H

#include  <iostream>

#include "larrecodnn/CVN/func/PixelMap.h"

namespace lcvn
{
 class SBNDPixelMap : public PixelMap
 {
  public:
    SBNDPixelMap(unsigned int nWire, unsigned int nTdc, const Boundary& bound): PixelMap::PixelMap(nWire, nTdc, bound) {std::cout << "============ Calling the function SBNDPixelMap::SBNDPixelMap() ==============\n";}
    SBNDPixelMap(): PixelMap::PixelMap(){std::cout << "============ Calling the function SBNDPixelMap::SBNDPixelMap() ==============\n";}
    int fSliceID = -9999;
 };	
}

#endif  // LCVN_SBNDPIXELMAP_H
