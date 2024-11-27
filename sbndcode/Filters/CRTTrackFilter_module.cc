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


class CRTTrackFilter : public art::EDFilter {
public:
  explicit CRTTrackFilter(fhicl::ParameterSet const& p);
  //virtual ~CRTTrackFilter() { }
  virtual bool filter(art::Event& e) override;
  //void    reconfigure(fhicl::ParameterSet const& p);
  
private:
  // fcl parameters
  std::string fCRTTrackModuleLabel;
  double fXPosMin;     // minimum absolute value of x position of CRT track endpoints, in cm
  double fXPosMax;     // maximum absolute value of x position of CRT track endpoints, in cm
  double fThetaXZMax;  // maximum angle theta_XZ of CRT track, in degrees
  bool fUseLengthMin;  // flag to use the fancy calculation of projected CRT track length the active volume, and to cut on that
  double fLengthMin;   // if above is true, then this is minimum length in the active volume, in cm
  bool fVerbose;
  std::string fHitsModuleLabel;     ///< Label for Hit dataproduct (to be set via fcl)
  // useful constants
  const std::set<sbnd::crt::CRTTagger> fNSTaggers = { sbnd::crt::kSouthTagger, sbnd::crt::kNorthTagger };
  const double PI = 3.14159265358979323;
  const double z_AV_S   =    0.;
  const double z_AV_N   =  500.;
  const double y_AV_top =  200.;
  const double y_AV_bot = -200.;
  int _nhits;
};

CRTTrackFilter::CRTTrackFilter(fhicl::ParameterSet const& p) : EDFilter{p} {
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel","crttracks");
  fXPosMin             = p.get<double>("XPosMin",150.);
  fXPosMax             = p.get<double>("XPosMax",200.);
  fThetaXZMax          = p.get<double>("ThetaXZMax",10.);
  fUseLengthMin        = p.get<bool>("UseLengthMin",true);
  fLengthMin           = p.get<double>("LengthMin",300.);
  fVerbose             = p.get<bool>("Verbose",false);
  fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel");
}

bool CRTTrackFilter::filter(art::Event& e) {
  
  if(fVerbose) std::cout << "Event " << e.id().event() << "..." << std::endl;
  
  // Get the CRT track data product
  art::Handle<std::vector<sbnd::crt::CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  auto const& CRTTrackVec(*CRTTrackHandle);

  //Get Number of hits
  art::Handle<std::vector<recob::Hit>> hitListHandle;
  std::vector<art::Ptr<recob::Hit>> hitlist;
  if (e.getByLabel(fHitsModuleLabel,hitListHandle)) {
    art::fill_ptr_vector(hitlist, hitListHandle);
  }
  if(fVerbose) std::cout << "    Number of hits " << hitlist.size() << std::endl;
  if(hitlist.size()>21000) return false;
  
  // Loop over CRT tracks in the event...
  for(size_t i_ct=0; i_ct<CRTTrackVec.size(); ++i_ct) {
    
    if(fVerbose) std::cout << "  Track " << i_ct << "..." << std::endl;

    // Get this CRT track
    auto const& crt_trk = CRTTrackVec[i_ct];
    
    // Check whether it is a north-south crossing muon track
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
    
    // Calculate (signed) difference between north and south endpoints of the CRT track
    double Delta_x = x_N-x_S;
    double Delta_y = y_N-y_S;
    double Delta_z = z_N-z_S;

    // Check whether both CRT track points satisfy the x position requirements
    if( not(abs(x_S) > fXPosMin and abs(x_S) < fXPosMax) ) continue;
    if(fVerbose) std::cout << "    Found track that meets x position requirements on south point" << std::endl;
    if( not(abs(x_N) > fXPosMin and abs(x_N) < fXPosMax) ) continue;
    if(fVerbose) std::cout << "    Found track that meets x position requirements on both points" << std::endl;
    // Check whether Delta x between the CRT track points satisfies the requirement: Delta x < tan(10deg) * Delta z
    if( not(abs(Delta_x) < tan(fThetaXZMax*PI/180)*abs(Delta_z)) ) continue;
    if(fVerbose) std::cout << "    Found track that meets Delta x requirement" << std::endl;
    // If not using the fancy cut on a minimum length in the active volume, then at least...
    // Check whether at least one CRT track point is below the top of the TPC active volume
    if( not(fUseLengthMin) ) {
      if( not(y_N < y_AV_top or y_S < y_AV_top) ) continue;
      if(fVerbose) std::cout << "    Found track that meets simple y position requirement" << std::endl;
    }
    // ...Else, do the fancy thing...
    // Check whether track length inside the active volume satisfies the minimum length requirement...
    //   Note: All of the below assumes that the track does not exit through the anode, which *should* be true given fXPosMax = 200cm
    if(fUseLengthMin) {
      if(fVerbose) std::cout << "      Evaluating points where CRT track is projected to cross active volume boundaries" << std::endl;
      // Initialize vector for points where the CRT track crosses the boundary of the active volume in (y,z) space
      std::vector<geo::Point_t> crossing_points;
      // Calculate if/where the track crosses the south face of the TPC
      double y_AV_S = (Delta_y/Delta_z)*(z_AV_S-z_S)+y_S;
      if( y_AV_S > y_AV_bot and y_AV_S < y_AV_top ) {
	double x_AV_S = (Delta_x/Delta_z)*(z_AV_S-z_S)+x_S;
	crossing_points.push_back( {x_AV_S, y_AV_S, z_AV_S} );
	if(fVerbose) std::cout << "      Found crossing point on south face of active volume at (" << x_AV_S << ", " << y_AV_S << ", " << z_AV_S << ")" << std::endl;
      }
      // Calculate if/where the track crosses the north face of the TPC
      double y_AV_N = (Delta_y/Delta_z)*(z_AV_N-z_S)+y_S;
      if( y_AV_N > y_AV_bot and y_AV_N < y_AV_top ) {
	double x_AV_N = (Delta_x/Delta_z)*(z_AV_N-z_S)+x_S;
	crossing_points.push_back( {x_AV_N, y_AV_N, z_AV_N} );
	if(fVerbose) std::cout << "      Found crossing point on north face of active volume at (" << x_AV_N << ", " << y_AV_N << ", " << z_AV_N << ")" << std::endl;
      }
      // Calculate if/where the track crosses the top of the TPC
      double z_AV_top = (Delta_z/Delta_y)*(y_AV_top-y_S)+z_S;
      if( z_AV_top > z_AV_S and z_AV_top < z_AV_N ) {
	double x_AV_top = (Delta_x/Delta_y)*(y_AV_top-y_S)+x_S;
	crossing_points.push_back( {x_AV_top, y_AV_top, z_AV_top} );
	if(fVerbose) std::cout << "      Found crossing point on top face of active volume at (" << x_AV_top << ", " << y_AV_top << ", " << z_AV_top << ")" << std::endl;
      }
      // Calculate if/where the track crosses the bottom of the TPC
      double z_AV_bot = (Delta_z/Delta_y)*(y_AV_bot-y_S)+z_S;
      if( z_AV_bot > z_AV_S and z_AV_bot < z_AV_N ) {
	double x_AV_bot = (Delta_x/Delta_y)*(y_AV_bot-y_S)+x_S;
	crossing_points.push_back( {x_AV_bot, y_AV_bot, z_AV_bot} );
	if(fVerbose) std::cout << "      Found crossing point on bottom face of active volume at (" << x_AV_bot << ", " << y_AV_bot << ", " << z_AV_bot << ")" << std::endl;
      }
      // Print a warning if there aren't exactly zero or two crossing points, which would imply a bug in the above
      if( not(crossing_points.size()==0 or crossing_points.size()==2)) std::cout << "    WARNING: found " << crossing_points.size() << " active volume crossings (expected exactly 0 or exactly 2)" << std::endl;
      // Calculate the length of the track inside the active volume
      double trk_len_AV = 0.;
      if(crossing_points.size()==2) trk_len_AV = (crossing_points[1] - crossing_points[0]).R();
      if(fVerbose) std::cout << "      Estimated track length in active volume is " << trk_len_AV << "cm" << std::endl;
      // Check whether track length inside the active volume satisfies the minimum length requirement
      if( trk_len_AV < fLengthMin ) continue;
      if(fVerbose) std::cout << "    Found track that meets track length in active volume requirement" << std::endl;
    }
    if(crt_trk.Ts0()<-2500 || crt_trk.Ts0()>-500) continue;    
    if(abs(crt_trk.ToF()-(crt_trk.Length()/29.9792))>15) continue;

    // If we've made it this far, then this track is good and this event passes the filter
    if(fVerbose) std::cout << "    Found a good track!" << std::endl;
    if(fVerbose) std::cout << "    ToF " << crt_trk.ToF() << "ns" << std::endl;
    if(fVerbose) std::cout << "    ts0 " << crt_trk.Ts0() << "ns" << std::endl;
    if(fVerbose) std::cout << "    abs(tof - (length/29.9792...)) " << abs(crt_trk.ToF()-(crt_trk.Length()/29.9792))<<std::endl;
    if(fVerbose) std::cout << "    Number of hits " << hitlist.size() << std::endl;  
    return true;

  } // end loop over CRT tracks

  // If none of the CRT tracks in this event pass the criteria above, then this event does not pass the filter
  return false;

}

// A macro required for a JobControl module.
DEFINE_ART_MODULE(CRTTrackFilter)
