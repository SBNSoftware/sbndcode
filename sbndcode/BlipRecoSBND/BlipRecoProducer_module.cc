////////////////////////////////////////////////////////////////////////
// Class:       BlipRecoProducer
// Plugin Type: producer (art v2_11_03)
// File:        BlipRecoProducer_module.cc
//
// W. Foreman
// May 2022
////////////////////////////////////////////////////////////////////////

// Framework includes 
#include "art/Framework/Core/EDProducer.h"
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
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/Utilities/AssociationUtil.h"

// C++ includes
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <bitset>
#include <memory>

// SBND-specific includes
#include "sbndcode/BlipRecoSBND/Alg/BlipRecoAlg.h"


class BlipRecoProducer;

class BlipRecoProducer : public art::EDProducer 
{
  public:
  explicit BlipRecoProducer(fhicl::ParameterSet const & p);
  BlipRecoProducer(BlipRecoProducer const &) = delete;
  BlipRecoProducer(BlipRecoProducer &&) = delete;
  BlipRecoProducer & operator = (BlipRecoProducer const &) = delete;
  BlipRecoProducer & operator = (BlipRecoProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  private:
  blip::BlipRecoAlg*      fBlipAlg;
  std::string             fHitProducer;

};



//###################################################
//  BlipRecoProducer constructor and destructor
//###################################################
BlipRecoProducer::BlipRecoProducer(fhicl::ParameterSet const & pset)
   : EDProducer(pset)
{
  // Read in fcl parameters for blip reco alg
  fhicl::ParameterSet pset_blipalg = pset.get<fhicl::ParameterSet>("BlipAlg");
  fBlipAlg        = new blip::BlipRecoAlg( pset_blipalg );
  
  fHitProducer    = pset_blipalg.get<std::string>   ("HitProducer");
 
  // declare what we're going to produce
  produces< std::vector<  recob::SpacePoint > >();
  produces< art::Assns <  recob::Hit, recob::SpacePoint> >();
  
  //produces< std::vector<  blip::Blip > >();
  //produces< art::Assns <  recob::SpacePoint,  blip::Blip> >();
  //produces< art::Assns <  recob::Hit,         blip::Blip> >();

  //produces< std::vector<  recob::Cluster    > >();
  //produces< art::Assns <  recob::Cluster,   recob::SpacePoint> >();
  //produces< art::Assns <  recob::Hit,       recob::Cluster> >();
  
}



//###################################################
//  Main produce function
//###################################################
void BlipRecoProducer::produce(art::Event & evt)
{
  
  std::cout<<"\n"
  <<"******************** BlipReco ****************************\n"
  <<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"\n";

  //============================================
  // Make unique pointers to the vectors of objects 
  // and associations we will create
  //============================================
  std::unique_ptr< std::vector< recob::SpacePoint> > SpacePoint_v(new std::vector<recob::SpacePoint>);
  std::unique_ptr< art::Assns <recob::Hit, recob::SpacePoint> >  assn_hit_sps_v(new art::Assns<recob::Hit,recob::SpacePoint> );

  //std::unique_ptr< std::vector< blip::Blip > > Blip_v(new std::vector<blip::Blip>);
  
  //============================================
  // Get hits from input module
  //============================================
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitProducer,hitHandle))
    art::fill_ptr_vector(hitlist, hitHandle);
  
  
  //============================================
  // Run blip reconstruction: 
  //============================================
  fBlipAlg->RunBlipReco(evt);
  
  
  //===========================================
  // Make blip::Blips
  //===========================================
  /*
  for(size_t i=0; i<fBlipAlg->blips.size(); i++){
    auto& b = fBlipAlg->blips[i];
    blip::Blip nb = b;
    Blip_v->emplace_back(nb);
  }
  */
  
  
  //===========================================
  // Make recob::SpacePoints out of the blip::Blips
  //===========================================
  for(size_t i=0; i<fBlipAlg->blips.size(); i++){
    auto& b = fBlipAlg->blips[i];
    
    Double32_t xyz[3];
    Double32_t xyz_err[6];
    Double32_t chiSquare = 0;
    Double32_t err = 0.; //b.MaxIntersectDiff?
    xyz[0]      = b.Position.X();
    xyz[1]      = b.Position.Y();
    xyz[2]      = b.Position.Z();
    xyz_err[0]  = err;
    xyz_err[1]  = err;
    xyz_err[2]  = err;
    xyz_err[3]  = err;
    xyz_err[4]  = err;
    xyz_err[5]  = err;
    
    recob::SpacePoint newpt(xyz,xyz_err,chiSquare);
    SpacePoint_v->emplace_back(newpt);
    
    
    // Hit associations 
    for(auto& hc : b.clusters ) {
      for(auto& ihit : hc.HitIDs ) {
        auto& hitptr = hitlist[ihit];
        util::CreateAssn(*this, evt, *SpacePoint_v, hitptr, *assn_hit_sps_v);
      }
    }
  
  }
    
  //===========================================
  // Put them on the event
  //===========================================
  evt.put(std::move(SpacePoint_v));
  evt.put(std::move(assn_hit_sps_v));
  
  //evt.put(std::move(Blip_v));

}//END EVENT LOOP

DEFINE_ART_MODULE(BlipRecoProducer)
