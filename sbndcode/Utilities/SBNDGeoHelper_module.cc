////////////////////////////////////////////////////////////////////////
// Class:       SBNDGeoHelper
// Plugin Type: analyzer (Unknown Unknown)
// File:        SBNDGeoHelper_module.cc
//
// Generated at Wed Jan 17 13:10:15 2024 by Tingjun Yang using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"

#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"

#include <iostream>

using namespace std;

namespace util {
  class SBNDGeoHelper;
}


class util::SBNDGeoHelper : public art::EDAnalyzer {
public:
  explicit SBNDGeoHelper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDGeoHelper(SBNDGeoHelper const&) = delete;
  SBNDGeoHelper(SBNDGeoHelper&&) = delete;
  SBNDGeoHelper& operator=(SBNDGeoHelper const&) = delete;
  SBNDGeoHelper& operator=(SBNDGeoHelper&&) = delete;

private:
  void analyze(art::Event const& e) override {}
};


util::SBNDGeoHelper::SBNDGeoHelper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  art::ServiceHandle<SBND::TPCChannelMapService> channelMap;

  int input = -1;
  while (input!=0){
    cout<<"0: quit"<<endl;
    cout<<"1: convert tpc/plane/wire to offline channel number"<<endl;
    cout<<"2: find intersection of two wires"<<endl;
    cout<<"Type a number: ";
    cin >> input;
    if (input == 1){
      cout<<"Convert tpc/plane/wire to offline channel number"<<endl;
      int tpc, plane, wire;
      cout<<"TPC number: ";
      cin>>tpc;
      cout<<"Plane number: ";
      cin>>plane;
      cout<<"Wire number: ";
      cin>>wire;
      geo::WireID wireid(0,tpc,plane,wire);
      auto const & channel = wireReadout.PlaneWireToChannel(wireid);
      auto const & chaninfo = channelMap->GetChanInfoFromOfflChan(channel);
      cout<<"##############################"<<endl;
      cout<<"TPC             = "<<tpc<<endl;
      cout<<"Plane           = "<<plane<<endl;
      cout<<"Wire            = "<<wire<<endl;
      cout<<"Offline channel = "<<channel<<endl;
      cout<<"WIBCrate (1-4)  = "<<chaninfo.WIBCrate<<endl;
      cout<<"WIB (1-6)       = "<<chaninfo.WIB<<endl;
      cout<<"FEMBOnWib (0-3) = "<<chaninfo.FEMBOnWIB<<endl;
      cout<<"FEMBCh (0-127)  = "<<chaninfo.FEMBCh<<endl;
      cout<<"FEMCrate (1-11) = "<<chaninfo.FEMCrate<<endl;
      cout<<"FEM      (1-16) = "<<chaninfo.FEM<<endl;
      cout<<"FEMCh (0-63)    = "<<chaninfo.FEMCh<<endl;
    }
    else if (input == 2){
      cout<<"Find wire intersection of two wires"<<endl;
      int tpc1, plane1, wire1, tpc2, plane2, wire2;
      cout<<"TPC1 number: ";
      cin>>tpc1;
      cout<<"Plane1 number: ";
      cin>>plane1;
      cout<<"Wire1 number: ";
      cin>>wire1;
      cout<<"TPC2 number: ";
      cin>>tpc2;
      cout<<"Plane2 number: ";
      cin>>plane2;
      cout<<"Wire2 number: ";
      cin>>wire2;
      geo::WireID wid1(0,tpc1,plane1,wire1);
      geo::WireID wid2(0,tpc2,plane2,wire2);
      geo::Point_t intersection;
      auto intersect = wireReadout.WireIDsIntersect(wid1,wid2,intersection);
      cout<<"WireID 1:"<<wid1.toString()<<endl;
      cout<<"WireID 2:"<<wid2.toString()<<endl;
      if (intersect){
        cout<<"The intersection point is ("<<intersection.Y()<<","<<intersection.Z()<<")."<<endl;
        cout<<"It is inside the detector."<<endl;
      }
      else{
        cout<<"It is outside the detector."<<endl;
      }
    }
  }

}

DEFINE_ART_MODULE(util::SBNDGeoHelper)
