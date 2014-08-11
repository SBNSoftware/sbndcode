
#include "NuAnaAlg.h"

namespace lar1nd{

  NuAnaAlg::NuAnaAlg(){
    
  }


  void NuAnaAlg::configureGeometry(art::ServiceHandle<geo::Geometry> ){
  
    xlow  = - geom -> DetHalfWidth();
    xhigh = geom -> DetHalfWidth();
    ylow  = - geom -> DetHalfHeight();
    yhigh = geom -> DetHalfHeight();
    zlow  = 0.0
    zhigh = geom -> DetLength();

    return;
  }


}