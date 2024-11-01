#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "sbnobj/SBND/CRT/CRTTrack.hh"


class CRTTrigFilter : public art::EDFilter {
public:
  explicit CRTTrigFilter(fhicl::ParameterSet const& p);
  //virtual ~CRTTrigFilter() { }
  virtual bool filter(art::Event& e) override;
  //void    reconfigure(fhicl::ParameterSet const& p);
  
private:
  std::string fCRTTrackModuleLabel;
  
};

CRTTrigFilter::CRTTrigFilter(fhicl::ParameterSet const& p) : EDFilter{p} {
  fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel","crt");
}

bool CRTTrigFilter::filter(art::Event& e) {
  
  bool KeepMe = false;

  

  return KeepMe;

}

// A macro required for a JobControl module.
DEFINE_ART_MODULE(CRTTrigFilter)
