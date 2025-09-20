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

// SBN/SBND includes
#include "sbnobj/SBND/Commissioning/MuonTrack.hh"

// ROOT includes 

// C++ includes
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <iostream>

class MuonTrackFilter : public art::EDFilter {
public:
   explicit MuonTrackFilter(fhicl::ParameterSet const& p);
   virtual bool filter(art::Event& e) override;
   void    reconfigure(fhicl::ParameterSet const& p);
   virtual ~MuonTrackFilter() { }
   
private:
   // services 
   art::ServiceHandle<art::TFileService> tfs;
};

MuonTrackFilter::MuonTrackFilter(fhicl::ParameterSet const& p): EDFilter{p} 
{}

void MuonTrackFilter::reconfigure(fhicl::ParameterSet const& p)
{}

bool MuonTrackFilter::filter(art::Event& evt)
{ 
   bool pass = false;
   int ncts = 0;
   art::Handle<std::vector<sbnd::comm::MuonTrack> > muonTrackListHandle;
   std::vector<art::Ptr<sbnd::comm::MuonTrack> > muontrklist;
   if (evt.getByLabel("MuonTrackProducer", muonTrackListHandle))  {
      art::fill_ptr_vector(muontrklist, muonTrackListHandle);
      ncts = muontrklist.size();
   }
   if (ncts > 0){
      pass = true;
   }
   return pass; 
}
// A macro required for a JobControl module.
DEFINE_ART_MODULE(MuonTrackFilter)