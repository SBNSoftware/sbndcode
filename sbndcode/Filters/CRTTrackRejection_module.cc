////////////////////////////////////////////////////////////////////////
// Class:       CRTTrackFilter
// Plugin Type: filter
// File:        CRTTrackFilter_module.cc
//
// Filter module to find events with north-south crossing muon tracks that are close and nearly parallel to the anode planes
// Useful for wire plane transparency analysis being led by Beth
//
// Author:      Lauren Yates (yatesla@fnal.gov) and Bethany McCusker (b.mccusker@lancaster.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include <cmath>


class CRTTrackRejection : public art::EDFilter {
public:
  explicit CRTTrackRejection(fhicl::ParameterSet const& p);
  //virtual ~CRTTrackFilter() { }
  virtual bool filter(art::Event& e) override;
  //void    reconfigure(fhicl::ParameterSet const& p);
  
private:
  // fcl parameters
  std::string fCRTTrackModuleLabel;

  double fXPosMax;     // maximum absolute value of x position of CRT track endpoints, in cm
  double fXPosMin;
  bool fVerbose;


  // useful constants
  const std::set<sbnd::crt::CRTTagger> fNSTaggers = { sbnd::crt::kSouthTagger, sbnd::crt::kNorthTagger };
  const double PI = 3.14159265358979323;


};

CRTTrackRejection::CRTTrackRejection(fhicl::ParameterSet const& p) : EDFilter{p} {
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel","crttracks");
  fXPosMax             = p.get<double>("XPosMax",300.);
  fXPosMin             = p.get<double>("XPosMin",300.);
  fVerbose             = p.get<bool>("Verbose",false);
}

bool CRTTrackRejection::filter(art::Event& e) {
  
  if(fVerbose) std::cout << "Event " << e.id().event() << "..." << std::endl;
  
  // Get the CRT track data product
  art::Handle<std::vector<sbnd::crt::CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  auto const& CRTTrackVec(*CRTTrackHandle);
 
  bool good_crt=false;   
  // Loop over CRT tracks in the event...
  for(size_t i_ct=0; i_ct<CRTTrackVec.size(); ++i_ct) {
    
   

    // Get this CRT track
    auto const& crt_trk = CRTTrackVec[i_ct];
    if(crt_trk.Ts0()<-2500 || crt_trk.Ts0()>-500) continue;
    if(fVerbose) std::cout << "    ts0 " << crt_trk.Ts0() << "ns" << std::endl;
    //got the trigger track now                                                     
    // Check whether it is a north-south crossing muon track

    if(fVerbose) std::cout << "  Track " << i_ct << "..." << std::endl;
    if(fVerbose) {
      std::cout << "    Taggers: ";
      for(auto tagger : crt_trk.Taggers()) std::cout << tagger << " ";
      std::cout << std::endl;
    }
    if(crt_trk.Taggers() != fNSTaggers) continue;
    if(fVerbose) std::cout << "    Found a north-south crossing muon track" << std::endl;
    
    // Get (x, y, z) coordinates of the south and north CRT track points
    double x_S, y_S, z_S, x_N, y_N, z_N;
    for(auto point : crt_trk.Points()) {
      if(fVerbose) std::cout << "      Point: " << "(" << point.X() << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
      if(point.Z() < 0.) {
	x_S = point.X();
	y_S = point.Y();
	z_S = point.Z();
      }
      if(point.Z() > 400.) {
	x_N = point.X();
	y_N = point.Y();
	z_N = point.Z();
      }
    }
    if(fVerbose) {
      std::cout << "    South point: (" << x_S << ", " << y_S << ", " << z_S << ")" << std::endl;
      std::cout << "    North point: (" << x_N << ", " << y_N << ", " << z_N << ")" << std::endl;
    }
    

   

    if( not(abs(x_S) > fXPosMin and abs(x_S) < fXPosMax) ) continue;

    if( not(abs(x_N) > fXPosMin and abs(x_N) < fXPosMax) ) continue;
    if(fVerbose){
      std::cout <<" good track" <<std::endl;
      }


    good_crt=true;

  } // end loop over CRT tracks

 
  if(good_crt) return true;
 else  return false;

}

// A macro required for a JobControl module.
DEFINE_ART_MODULE(CRTTrackRejection)

