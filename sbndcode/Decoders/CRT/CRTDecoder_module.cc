////////////////////////////////////////////////////////////////////////
// Class:       CRTDecoder
// Plugin Type: producer
// File:        CRTDecoder_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
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
#include "sbndaq-artdaq-core/Overlays/Common/BernCRTFragmentV2.hh"

#include "sbnobj/SBND/CRT/FEBData.hh"

class CRTDecoder;


class CRTDecoder : public art::EDProducer {
public:
  explicit CRTDecoder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTDecoder(CRTDecoder const&) = delete;
  CRTDecoder(CRTDecoder&&) = delete;
  CRTDecoder& operator=(CRTDecoder const&) = delete;
  CRTDecoder& operator=(CRTDecoder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  std::vector<sbnd::crt::FEBData> FragToFEB(const artdaq::Fragment &frag);

private:

  std::vector<std::pair<unsigned, unsigned>> fMac5ToGeoIDVec;
  std::map<unsigned, unsigned>               fMac5ToGeoID;

  std::string              fCRTModuleLabel;
  std::vector<std::string> fCRTInstanceLabels;
};


CRTDecoder::CRTDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTModuleLabel(p.get<std::string>("CRTModuleLabel", "daq"))
  , fCRTInstanceLabels(p.get<std::vector<std::string>>("CRTInstanceLabels"))
  {
    produces<std::vector<sbnd::crt::FEBData>>();
    //  produces<art::Assns<artdaq::Fragment,sbnd::crt::FEBData>>();

    fMac5ToGeoIDVec = p.get<std::vector<std::pair<unsigned, unsigned>>>("FEBMac5ToGeometryIDMap");
    fMac5ToGeoID    = std::map<unsigned, unsigned>(fMac5ToGeoIDVec.begin(), fMac5ToGeoIDVec.end());
  }

void CRTDecoder::produce(art::Event& e)
{
  auto febDataVec      = std::make_unique<std::vector<sbnd::crt::FEBData>>();
  //  auto febDataFragAssn = std::make_unique<art::Assns<artdaq::Fragment,sbnd::crt::FEBData>>();

  for(const std::string &CRTInstanceLabel : fCRTInstanceLabels)
    {
      art::Handle<std::vector<artdaq::Fragment>> fragmentHandle;
      e.getByLabel(fCRTModuleLabel, CRTInstanceLabel, fragmentHandle);

      if(!fragmentHandle.isValid() || fragmentHandle->size() == 0)
        continue;

      if(fragmentHandle->front().type() == artdaq::Fragment::ContainerFragmentType)
        {
          for(auto cont : *fragmentHandle)
            {
              artdaq::ContainerFragment contf(cont);
              if(contf.fragment_type() == sbndaq::detail::FragmentType::BERNCRTV2)
                {
                  for(unsigned i = 0; i < contf.block_count(); ++i)
                    {
                      std::vector<sbnd::crt::FEBData> newFebDatas = FragToFEB(*contf[i].get());
                      febDataVec->insert(febDataVec->end(), newFebDatas.begin(), newFebDatas.end());
                    }
                }
            }
        }
      else if(fragmentHandle->front().type() == sbndaq::detail::FragmentType::BERNCRTV2)
        {
          for(auto frag : *fragmentHandle)
            {
              std::vector<sbnd::crt::FEBData> newFebDatas = FragToFEB(frag);
              febDataVec->insert(febDataVec->end(), newFebDatas.begin(), newFebDatas.end());
            }
        }
    }
  
  e.put(std::move(febDataVec));
}

std::vector<sbnd::crt::FEBData> CRTDecoder::FragToFEB(const artdaq::Fragment &frag)
{
  std::vector<sbnd::crt::FEBData> feb_datas;

  const sbndaq::BernCRTFragmentV2 bern_frag(frag);
  const sbndaq::BernCRTFragmentMetadataV2* bern_frag_meta = bern_frag.metadata();
  
  for(unsigned i = 0; i < bern_frag_meta->hits_in_fragment(); ++i)
    {
      const sbndaq::BernCRTHitV2 *bern_hit = bern_frag.eventdata(i);

      std::array<uint16_t, 32> adc_array;
      unsigned ii = 0;
      for(auto const &adc : bern_hit->adc)
        {
          adc_array[ii] = adc;
          ++ii;
        }

      feb_datas.emplace_back(fMac5ToGeoID[bern_frag_meta->MAC5()],
                             bern_hit->flags,
                             bern_hit->ts0,
                             bern_hit->ts1,
                             bern_hit->timestamp / 1e9,
                             adc_array,
                             bern_hit->coinc);

      std::string adc_string = "";

      for(auto const &adc : feb_datas.back().ADC())
        {
          adc_string += std::to_string(adc);
          adc_string += ", ";
        }
      adc_string.resize(adc_string.size() - 2);

      mf::LogInfo("CRTDecoder")
        << "Creating FEBData object from BernCRT Fragment\n"
        << "Mac5: " << feb_datas.back().Mac5() << '\n'
        << "Flags: " << feb_datas.back().Flags() << '\n'
        << "Ts0: " << feb_datas.back().Ts0() << '\n'
        << "Ts1: " << feb_datas.back().Ts1() << '\n'
        << "UnixS: " << feb_datas.back().UnixS() << '\n'
        << "ADC: [" << adc_string << "]\n"
        << "Coinc: " << feb_datas.back().Coinc() << '\n' << std::endl;
    }

  return feb_datas;
}

DEFINE_ART_MODULE(CRTDecoder)
