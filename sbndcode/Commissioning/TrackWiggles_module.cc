#ifndef TRACKWIGGLES_H
#define TRACKWIGGLES_H

// Framework includes
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDFilter.h" 
#include "canvas/Persistency/Common/FindMany.h" 
#include "art/Framework/Services/Registry/ServiceRegistry.h"


// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" 
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "cetlib/search_path.h"
#include "canvas/Persistency/Common/Assns.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/RecoBase/Wire.h"


//sbnd includes
#include "sbndcode/Utilities/DigitalNoiseChannelStatus.h"

// C++ includes
#include <cstring>
#include <utility>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cmath>
#include <unistd.h>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"



class TrackWiggles;

class TrackWiggles : public art::EDFilter{
    public:
    explicit TrackWiggles(fhicl::ParameterSet const & p);
    virtual bool filter(art::Event& e) override;
    virtual ~TrackWiggles() { };
    bool CheckForHitsInChannels(std::vector<art::Ptr<recob::Hit>, std::allocator<art::Ptr<recob::Hit> > > hitList, int ChannelOfInterest);
    std::vector<int> GetIntersectingChannels(int ChID);
    bool CheckForNoiseBlock(art::Event& e);
    //Store Fcl Params somewhere
    private:
    //Data members for analysis/output
    //Start with fcl params
    // Declare member data here.
    std::string fHitProducer; //For tracks
    std::string fTrackProducer; //For tracks
    std::string fTrackHitAssociator; //For tracks
    std::vector<int> fChannelsToInclude;
    geo::GeometryCore const&                fGeom;
    std::string fRecobWireLabel;
  };

TrackWiggles::TrackWiggles(fhicl::ParameterSet const& p): 
    EDFilter{p}, fGeom{*lar::providerFrom<geo::Geometry>()}
    // More initializers here.
  {
    //Read in assorted fcl parameters
    fHitProducer = p.get< std::string >("HitProducer", "gaushit" );
    fTrackProducer = p.get< std::string >("TrackProducer", "pandoraTrack" );
    fTrackHitAssociator = p.get< std::string>("TrackHitAssosciator", "pandoraTrack");
    fChannelsToInclude = p.get< std::vector<int> >("Channels" );
    fRecobWireLabel = p.get<std::string>("RecobWireLabel","caldata");
    std::cout <<"Including Channels";
    for(int i=0; i<int(fChannelsToInclude.size()); i++)
    {
        std::cout << " " << fChannelsToInclude[i] << " ";
    }
    std::cout<<std::endl;
    //auto TempTest = GetIntersectingChannels( 0);
    //auto TempTest2 = GetIntersectingChannels( 3000);
    std::cout << "all done setting up" << std::endl;
    //made list of interesection channels
    //art::ServiceHandle<geo::Geometry> geometry;
    double y =-195;//-213 is the best fit location from PDS, but thats out of TPC;
    double z = 20;
    unsigned int Plane_0=0;
    unsigned int Plane_1=1;
    unsigned int Plane_2=2;
    unsigned int TPC =1;
    unsigned int Cryostat=0;
    geo::PlaneID TempPlane_0(Cryostat, TPC, Plane_0);
    geo::PlaneID TempPlane_1(Cryostat, TPC, Plane_1);
    geo::PlaneID TempPlane_2(Cryostat, TPC, Plane_2);
    geo::Point_t TempPoint(0.0, y, z);
    std::cout << "Closest wire location in TPC 1 plane 0 "<< fGeom.NearestChannel(TempPoint, TempPlane_0) << std::endl;
    std::cout << "Closest wire location in TPC 1 plane 1 "<< fGeom.NearestChannel(TempPoint, TempPlane_1) << std::endl;
    std::cout << "Closest wire location in TPC 1 plane 2 "<< fGeom.NearestChannel(TempPoint, TempPlane_2) << std::endl;

    TPC=0;
    y=195;
    z=375;
    geo::PlaneID TPC0_TempPlane_0(Cryostat, TPC, Plane_0);
    geo::PlaneID TPC0_TempPlane_1(Cryostat, TPC, Plane_1);
    geo::PlaneID TPC0_TempPlane_2(Cryostat, TPC, Plane_2);
    geo::Point_t TempPoint_2(0.0, y, z);
    std::cout << "What about the point in the opposite most corner" << std::endl;
    std::cout << "Closest wire location in TPC 0 plane 0 "<< fGeom.NearestChannel(TempPoint_2, TPC0_TempPlane_0) << std::endl;
    std::cout << "Closest wire location in TPC 0 plane 1 "<< fGeom.NearestChannel(TempPoint_2, TPC0_TempPlane_1) << std::endl;
    std::cout << "Closest wire location in TPC 0 plane 2 "<< fGeom.NearestChannel(TempPoint_2, TPC0_TempPlane_2) << std::endl;
    sleep(20);
  }


std::vector<int> TrackWiggles::GetIntersectingChannels(int ChID)
{
  //Build vectors of intersecting channels
    double y=0, z=0;
    art::ServiceHandle<geo::Geometry> geometry;
    std::vector<int> TempIntersections;
    int Start=0;
    int TotalWires=11264; //Only check wires from the same TPC
    if(ChID>=5632) Start=5632;
    else TotalWires=5632; 
    for(int j=Start; j<TotalWires; j++) 
    {
      if( (geometry->View(j) == geometry->View(ChID)) ) 
      {
        continue;
      }
      if( geometry->ChannelsIntersect(ChID,j,y,z) )
      {
        TempIntersections.push_back(j);
      } 
    }
    return TempIntersections;
}

bool TrackWiggles::filter(art::Event& evt)
{ 
   art::ServiceHandle<geo::Geometry> geometry;
   sbnd::DigitalNoiseChannelStatus * DigitalNoiseSVC =  &*art::ServiceHandle<sbnd::DigitalNoiseChannelStatus>();
   bool pass = false;
   //Load in hits and tracks and associations between those objects
   art::Handle< std::vector<recob::Hit> > hitHandle; //User handle for vector of OpDetWaveforms
   std::vector<art::Ptr<recob::Hit>> hitList;
   evt.getByLabel(fHitProducer, hitHandle);
   //art::Handle< std::vector<recob::Track> > trackHandle;
   //std::vector<art::Ptr<recob::Track>> trackList;
  // evt.getByLabel(fTrackProducer,trackHandle);
   //art::Handle< art::Assns<recob::Track,recob::Hit,void> > handle_TrackHitAssociation;
   //art::Assns<recob::Track,recob::Hit,void> TrackHitAssociation;
   //evt.getByLabel(fTrackHitAssociator,handle_TrackHitAssociation);
   //Replace artHandles with Normal vectors 
   art::fill_ptr_vector(hitList, hitHandle);
   //Stuff using actual tracks
   //art::fill_ptr_vector(trackList, trackHandle);
   //art::fill_ptr_vector(TrackHitAssociation, handle_TrackHitAssociation); // Probably doesn't work?
   //art::FindMany<recob::Hit> FoundTrackHit(trackHandle, evt, fTrackHitAssociator); //Maybe dont need to use association directly
   //Extract involved in making each track and compare to the requested list of channels
   //for(int i=0; i<int(trackList.size()); i++)
   //{
   // auto hits_for_track  = FoundTrackHit.at(i); 
    //See if this track has any hits in the relevant sections
   // for(int j=0; j<int(hits_for_track.size()); j++)
   // {
   //     unsigned int ChID = (hits_for_track)[j]->Channel(); 
        //std::cout << ChID << std::endl;
        //Check if ChID is in our included channels vector
   //     pass = std::any_of(fChannelsToInclude.begin(), fChannelsToInclude.end(), [ChID](int Ch) { return Ch == int(ChID); });
   //     if(pass) break; //stop both loops
   // }
   // if(pass) break; //exiting this loop too 
   // }

   //Stuff Just considering hits
   //Require that whole list of channels to include is matched
   std::vector<bool> PassChannels;
   for(int i=0; i<int(fChannelsToInclude.size()); i++)
   {
    PassChannels.push_back(false);
   }

    for(int j=0; j<int(fChannelsToInclude.size()); j++)
    {
       bool BadChannel = DigitalNoiseSVC->IsBad( raw::ChannelID_t(fChannelsToInclude[j]) );
       if(BadChannel)
       {
        std::cout << "One of the requested channels is bad, cutting event" << std::endl;
        return false;
       }
       bool HitsChannelOfInterest =CheckForHitsInChannels(hitList, fChannelsToInclude[j]);
       bool HitOtherPlaneChannel = false;
       if(HitsChannelOfInterest)
       {
        auto TempOverlapWires = GetIntersectingChannels( fChannelsToInclude[j]);
        int NumOtherHits=0;
        for(int k=0; k<int(TempOverlapWires.size()); k++)
        {
          HitOtherPlaneChannel = CheckForHitsInChannels(hitList, TempOverlapWires[k]);
          if(HitOtherPlaneChannel) NumOtherHits++;
          if(NumOtherHits>10) break;
        }
       }
       PassChannels[j] = (HitOtherPlaneChannel && HitsChannelOfInterest);
    }
   int MinChannels=12;
   pass = (std::all_of(PassChannels.begin(), PassChannels.end(), [](bool FoundChannel) { return FoundChannel; }) || 
            (std::count(PassChannels.begin(), PassChannels.end(), true)>=MinChannels) );
    if(pass)
    {
      bool NoiseBlock = CheckForNoiseBlock(evt);
      if(NoiseBlock) return false; // else return pass (true)
    }
    std::cout << "About to return pass " << pass << std::endl;
   return pass;
}   


bool TrackWiggles::CheckForNoiseBlock(art::Event& evt)
{
  auto const& rbwires = evt.getProduct<std::vector<recob::Wire>>(fRecobWireLabel);
  bool FoundNoiseBlock = false;
  for(int j=0; j<int(fChannelsToInclude.size()); j++)
  {
    const auto& rw = rbwires[fChannelsToInclude[j]];
    rw.NSignal(); 
    //do waveform analysis on some neighboring channels to find a large rectangular block
  }
  return FoundNoiseBlock;

}


bool TrackWiggles::CheckForHitsInChannels(std::vector<art::Ptr<recob::Hit>, std::allocator<art::Ptr<recob::Hit> > > hitList, int ChannelOfInterest)
{
  bool FoundGoodHit=false;
  for(int i=0; i<int((hitList).size()); i++)
   {
    unsigned int ChID = (hitList)[i]->Channel();
    if(int(ChID)==ChannelOfInterest) 
    {
      FoundGoodHit = true;
    }
   }
   return FoundGoodHit;
}
// Magic ART line that has to exist
DEFINE_ART_MODULE(TrackWiggles)
#endif
