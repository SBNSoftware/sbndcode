////////////////////////////////////////////////////////////////////////
// Class:       MetricFilter
// Plugin Type: filter (Unknown Unknown)
// File:        MetricFilter_module.cc
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

#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"

#include <memory>

namespace sbndaq {
  class MetricFilter;
}


class sbndaq::MetricFilter : public art::EDFilter {
public:
  explicit MetricFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MetricFilter(MetricFilter const&) = delete;
  MetricFilter(MetricFilter&&) = delete;
  MetricFilter& operator=(MetricFilter const&) = delete;
  MetricFilter& operator=(MetricFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:
  std::string fmetric_instance_name;
  std::vector<std::string> fmetric_module_labels;              

  float ftime_min; // us
  float ftime_max; // us
  float fpe_thresh;
  
  // Declare member data here.

};


sbndaq::MetricFilter::MetricFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
  // More initializers here.
{
  fmetric_module_labels = p.get<std::vector<std::string>>("metric_module_labels",{"pmtmetricproducer",
                                                                                 "pmtmetricbnbzero",
                                                                                 "pmtmetricbnblight",
                                                                                 "pmtmetricoffbeamzero",
                                                                                 "pmtmetricoffbeamlight",
                                                                                 "pmtmetriccrossingmuon"});
  fmetric_instance_name = p.get<std::string>("metric_instance_name","");
  ftime_min = p.get<float>("time_min",0);
  ftime_max = p.get<float>("time_max",3);
  fpe_thresh = p.get<float>("pe_thresh",50);
}

bool sbndaq::MetricFilter::filter(art::Event& e)
{
  bool pass = false;

  bool found_metrics = false;
  size_t nlabels = fmetric_module_labels.size();
  for (size_t i=0; i<nlabels; i++){
    auto label = fmetric_module_labels[i];

    art::Handle<std::vector<sbnd::trigger::pmtSoftwareTrigger>> metricHandle;
    if (fmetric_instance_name.empty()){
      e.getByLabel(label,metricHandle);
    }
    else{
      e.getByLabel(label,fmetric_instance_name,metricHandle);
    }

    if (!metricHandle.isValid() || metricHandle->size() == 0)
      continue;
    found_metrics = true;
    const std::vector<sbnd::trigger::pmtSoftwareTrigger> metric_v(*metricHandle);
    auto metric = metric_v[0];
    auto flash_peakpe   = metric.peakPE;
    auto flash_peaktime = metric.peaktime; // us

    bool pass_timing = (flash_peaktime > ftime_min) && (flash_peaktime < ftime_max);
    bool pass_pe = flash_peakpe > fpe_thresh;
    pass = (pass_timing && pass_pe);
  }
  if (!found_metrics){
    std::cout << "No Metrics found" << std::endl;
    pass=false;
  }
  return pass;
}

DEFINE_ART_MODULE(sbndaq::MetricFilter)
