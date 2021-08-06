// Framework includes 
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

// SBN/SBND includes
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"

// ROOT includes 
#include "TRandom3.h"

// C++ includes
#include <map>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <bitset>

class testSimpleFilter : public art::EDFilter {
public:
   explicit testSimpleFilter(fhicl::ParameterSet const& p);
   virtual bool filter(art::Event& e) override;
   void    reconfigure(fhicl::ParameterSet const& p);
   virtual ~testSimpleFilter() { }
   
private:
   int ncts;
   // services 
   art::ServiceHandle<art::TFileService> tfs;
};

testSimpleFilter::testSimpleFilter(fhicl::ParameterSet const& p): EDFilter{p} 
{}

void testSimpleFilter::reconfigure(fhicl::ParameterSet const& p)
{}

bool testSimpleFilter::filter(art::Event& evt)
{ 
   bool pass = false;
   ncts = 0;
   art::Handle<std::vector<sbn::crt::CRTTrack> > crtTrackListHandle;
   std::vector<art::Ptr<sbn::crt::CRTTrack> > ctrklist;
   if (evt.getByLabel("crttrack", crtTrackListHandle))  {
      art::fill_ptr_vector(ctrklist, crtTrackListHandle);
      ncts = ctrklist.size();
   }
   if (ncts > 0){
      pass = true;
   }
   return pass; 
}
// A macro required for a JobControl module.
DEFINE_ART_MODULE(testSimpleFilter)