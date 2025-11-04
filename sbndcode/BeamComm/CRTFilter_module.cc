////////////////////////////////////////////////////////////////////////
// Class:       CRTFilter
// Plugin Type: filter (Unknown Unknown)
// File:        CRTFilter_module.cc
//
// Generated at Fri Dec 13 06:19:04 2024 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "art_root_io/TFileService.h"

#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

#include <memory>

namespace sbnd {
  class CRTFilter;
}


class sbnd::CRTFilter : public art::EDFilter {
public:
  explicit CRTFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTFilter(CRTFilter const&) = delete;
  CRTFilter(CRTFilter&&) = delete;
  CRTFilter& operator=(CRTFilter const&) = delete;
  CRTFilter& operator=(CRTFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:
  art::InputTag fCrtSpacePointLabel;
  art::InputTag fFrameLabel;

  float ftime_min; // us
  float ftime_max; // us
  
  // Declare member data here.

};


sbnd::CRTFilter::CRTFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fCrtSpacePointLabel = p.get<art::InputTag>("CrtSpacePointLabel", "crtspacepoints");
  fFrameLabel = p.get<art::InputTag>("FrameLabel", "frameshift");
  ftime_min = p.get<float>("time_min",-800);
  ftime_max = p.get<float>("time_max", 2800);
}

bool sbnd::CRTFilter::filter(art::Event& e)
{
  bool pass = false;

  //Get frame
  double frame_shift = 0;

  art::Handle<sbnd::timing::FrameShiftInfo> frameHandle;
  e.getByLabel(fFrameLabel, frameHandle);

  if (!frameHandle.isValid()){
      return pass; //no frame = fail filter
  }
  else{
      sbnd::timing::FrameShiftInfo const& frame(*frameHandle);
      frame_shift = frame.FrameApplyAtCaf();
  }

  art::Handle<std::vector<sbnd::crt::CRTSpacePoint>> crtSpacePointHandle;
  std::vector<art::Ptr<sbnd::crt::CRTSpacePoint>> crt_sp_v;
  e.getByLabel(fCrtSpacePointLabel, crtSpacePointHandle);

  //Get CRT SpacePoint
  if (!crtSpacePointHandle.isValid() || crtSpacePointHandle->size() == 0){
      return pass; //no crt spacepoint = fail filter
  }
  else{
      art::fill_ptr_vector(crt_sp_v, crtSpacePointHandle);
      
      for (auto const& crt_sp: crt_sp_v){

          if (!crt_sp->Complete()) continue; //skip bad spacepoint

          double ts = crt_sp->Ts0();
          ts += frame_shift;

          //require at least 1 spacepoint within beam spill
          if ((ts > ftime_min) && (ts < ftime_max)) pass = true;
      }
  }

  return pass;
}

DEFINE_ART_MODULE(sbnd::CRTFilter)
