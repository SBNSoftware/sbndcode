#ifndef __FLASHFINDERFMWKINTERFACE_H__
#define __FLASHFINDERFMWKINTERFACE_H__

//#include "FhiclLite/ConfigManager.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include <stdlib.h>

namespace lightana {

  //typedef ::fcllite::PSet Config_t;
  typedef fhicl::ParameterSet Config_t;

  std::vector<size_t> ListOpChannels(int cryostat=-1);
  
  std::vector<size_t> ListOpDets(int cryostat=-1);

  size_t NOpDets(int cryostat=-1);

  size_t OpDetFromOpChannel(size_t opch);

  void OpDetCenterFromOpChannel(size_t opch, double *xyz);

}
#endif

