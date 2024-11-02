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
#include <cmath>


class CRTTrackFilter : public art::EDFilter {
public:
  explicit CRTTrackFilter(fhicl::ParameterSet const& p);
  //virtual ~CRTTrackFilter() { }
  virtual bool filter(art::Event& e) override;
  //void    reconfigure(fhicl::ParameterSet const& p);
  
private:
  std::string fCRTTrackModuleLabel;
  double fXPosMin;
  double fXPosMax;
  bool fVerbose;
  const std::set<sbnd::crt::CRTTagger> fNSTaggers = { sbnd::crt::kSouthTagger, sbnd::crt::kNorthTagger };
};

CRTTrackFilter::CRTTrackFilter(fhicl::ParameterSet const& p) : EDFilter{p} {
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel","crttracks");
  fXPosMin = p.get<double>("XPosMin",150.);
  fXPosMax = p.get<double>("XPosMax",200.);
  fVerbose = p.get<bool>("Verbose",true);
}

bool CRTTrackFilter::filter(art::Event& e) {

  // Get the CRT track data product
  art::Handle<std::vector<sbnd::crt::CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  auto const& CRTTrackVec(*CRTTrackHandle);

  if(fVerbose) std::cout << "Event " << e.id().event() << "..." << std::endl;

  // Loop over CRT tracks in the event...
  for(size_t i_ct=0; i_ct<CRTTrackVec.size(); ++i_ct) {
    
    std::cout << "  Track " << i_ct << "..." << std::endl;

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

    // Get the (x, y, z) coordinates of the south and north CRT track points
    double x_S, y_S, z_S, x_N, y_N, z_N;
    for(auto point : crt_trk.Points()) {
      if(fVerbose) std::cout << "      Point N: " << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
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
      std::cout << "    South point: " << x_S << " " << y_S << " " << z_S << std::endl;
      std::cout << "    North point: " << x_N << " " << y_N << " " << z_N << std::endl;
    }

    // Check whether both CRT track points satisfy the x position requirements
    if( not(abs(x_N) > fXPosMin and abs(x_N) < fXPosMax) ) continue;
    if(fVerbose) std::cout << "    Found track that meets x position requirements on north point" << std::endl;
    if( not(abs(x_S) > fXPosMin and abs(x_S) < fXPosMax) ) continue;
    if(fVerbose) std::cout << "    Found track that meets x position requirements on both points" << std::endl;
    // Check whether Delta x between the CRT track points satisfies the requirement: Delta x < tan(10deg) * Delta z
    if( not(abs(x_N-x_S) < 0.176327*abs(z_N-z_S)) ) continue;
    if(fVerbose) std::cout << "    Found track that meets Delta x requirement" << std::endl;
    // Check whether at least one CRT track point is below the top of the TPC active volume
    // ... part of trying to make sure that there's a decent length of track in the active volume
    if( not(y_N < 200. or y_S < 200.) ) continue;
    if(fVerbose) std::cout << "    Found track that meets y position requirement" << std::endl;

    // If we've made it this far, then this track is good
    if(fVerbose) std::cout << "    Found a good track!" << std::endl;
    if(fVerbose) std::cout << "    ToF " << crt_trk.ToF() << "ns" << std::endl;
    return true;

  } // end loop over CRT tracks

  // If none of the CRT tracks in this event pass the criteria above, then this event does not pass the filter
  return false;

}

// A macro required for a JobControl module.
DEFINE_ART_MODULE(CRTTrackFilter)
