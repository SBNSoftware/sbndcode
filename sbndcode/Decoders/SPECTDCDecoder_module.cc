////////////////////////////////////////////////////////////////////////
// Class:       SPECTDCDecoder
// Plugin Type: producer
// File:        SPECTDCDecoder_module.cc
// Author:      Henry Lay         (h.lay@lancaster.ac.uk)
// Author:      Vu Chi Lan Nguyen (vclnguyen1@sheffield.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "artdaq-core/Data/ContainerFragment.hh"
#include "sbndaq-artdaq-core/Overlays/FragmentType.hh"
//#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

class SPECTDCDecoder;


class SPECTDCDecoder : public art::EDProducer {
public:
  explicit SPECTDCDecoder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SPECTDCDecoder(SPECTDCDecoder const&) = delete;
  SPECTDCDecoder(SPECTDCDecoder&&) = delete;
  SPECTDCDecoder& operator=(SPECTDCDecoder const&) = delete;
  SPECTDCDecoder& operator=(SPECTDCDecoder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  std::vector<sbnd::timing::DAQTimestamp> FragToDAQTimestamp(const artdaq::Fragment &frag);

private:

  std::vector<std::string> fSPECTDCModuleLabels;
};


SPECTDCDecoder::SPECTDCDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fSPECTDCModuleLabels(p.get<std::vector<std::string>>("SPECTDCModuleLabels"))
{
  produces<std::vector<sbnd::timing::DAQTimestamp>>();
  // produces<art::Assns<artdaq::Fragment,sbnd::timing::DAQTimestamp>>();
}

void SPECTDCDecoder::produce(art::Event& e)
{
  auto daqTimestampVec      = std::make_unique<std::vector<sbnd::timing::DAQTimestamp>>();
  // auto daqTimestampFragAssn = std::make_unique<art::Assns<artdaq::Fragment,sbnd::timing::DAQTimestamp>>();
  
  for(const std::string &SPECTDCModuleLabel : fSPECTDCModuleLabels)
    {
      art::Handle<std::vector<artdaq::Fragment>> fragmentHandle;
      e.getByLabel(SPECTDCModuleLabel, fragmentHandle);

      if(!fragmentHandle.isValid() || fragmentHandle->size() == 0)
        continue;

      if(fragmentHandle->front().type() == artdaq::Fragment::ContainerFragmentType)
        {
          for(auto cont : *fragmentHandle)
            {
              artdaq::ContainerFragment contf(cont);
              if(contf.fragment_type() == sbndaq::detail::FragmentType::TDCTIMESTAMP)
                {
                  for(unsigned i = 0; i < contf.block_count(); ++i)
                    {
                      std::vector<sbnd::timing::DAQTimestamp> newDAQTimestamps = FragToDAQTimestamp(*contf[i].get());
                      daqTimestampVec->insert(daqTimestampVec->end(), newDAQTimestamps.begin(), newDAQTimestamps.end());
                    }
                }
            }
        }
      else
	{
	  if(fragmentHandle->front().type() == sbndaq::detail::FragmentType::TDCTIMESTAMP)
	    {
	      for(auto frag : *fragmentHandle)
		{
		  std::vector<sbnd::timing::DAQTimestamp> newDAQTimestamps = FragToDAQTimestamp(frag);
		  daqTimestampVec->insert(daqTimestampVec->end(), newDAQTimestamps.begin(), newDAQTimestamps.end());
		}
	    }
	}
    }
  
  e.put(std::move(daqTimestampVec));
}

std::vector<sbnd::timing::DAQTimestamp> SPECTDCDecoder::FragToDAQTimestamp(const artdaq::Fragment &frag)
{
  std::vector<sbnd::timing::DAQTimestamp> daqTimestamps;

  return daqTimestamps;
}

DEFINE_ART_MODULE(SPECTDCDecoder)
