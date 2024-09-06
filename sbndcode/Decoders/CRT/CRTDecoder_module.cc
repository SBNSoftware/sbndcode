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

#include "sbndcode/ChannelMaps/CRT/CRTChannelMapService.h"

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

  void endJob() override;

private:

  std::string              fCRTModuleLabel;
  std::vector<std::string> fCRTInstanceLabels;
  bool                     fDebug;

  art::ServiceHandle<SBND::CRTChannelMapService> fCRTChannelMapService;

  std::set<unsigned int> fUnfoundMAC5s;
};


CRTDecoder::CRTDecoder(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTModuleLabel(p.get<std::string>("CRTModuleLabel", "daq"))
  , fCRTInstanceLabels(p.get<std::vector<std::string>>("CRTInstanceLabels"))
  , fDebug(p.get<bool>("Debug", false))
{
  produces<std::vector<sbnd::crt::FEBData>>();
}

void CRTDecoder::produce(art::Event& e)
{
  auto febDataVec = std::make_unique<std::vector<sbnd::crt::FEBData>>();

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
      SBND::CRTChannelMapService::ModuleInfo_t module = fCRTChannelMapService->GetModuleInfoFromFEBMAC5(bern_frag_meta->MAC5());

      if(!module.valid)
        {
          mf::LogInfo("CRTDecoder") << "===========================================================\n"
                                    << "ERROR: Cannot find simulation module for MAC5: "
                                    << unsigned(bern_frag_meta->MAC5()) << '\n'
                                    << "===========================================================\n"
                                    << std::endl;

          fUnfoundMAC5s.insert(bern_frag_meta->MAC5());

          continue;
        }

      const sbndaq::BernCRTHitV2 *bern_hit = bern_frag.eventdata(i);
      // Fill ADC Array. If channel order is swapped in the GDML
      // compared to reality then we fill the array in reverse.
      std::array<uint16_t, 32> adc_array;
      unsigned ii = module.channel_order_swapped ? 31 : 0;
      for(auto const &adc : bern_hit->adc)
        {
          adc_array[ii] = adc;

          if(module.channel_order_swapped)
            --ii;
          else
            ++ii;
        }

      // Timestamp field stores the T0 time corrected by the cable length
      // stored in the original fcl. We use it to get the cable delay for
      // the module and store it in the unused coinc slot.
      const int64_t whole_second_ns = static_cast<int64_t>(1e9);
      int cable_length = bern_hit->timestamp % whole_second_ns - bern_hit->ts0;

      if(cable_length < 0)
        cable_length += whole_second_ns;

      if(cable_length > 1000 || cable_length < 0)
        throw std::runtime_error("Why is the cable length: " + std::to_string(cable_length) + "?");

      feb_datas.emplace_back(module.offline_module_id,
                             bern_hit->flags,
                             bern_hit->ts0,
                             bern_hit->ts1,
                             bern_hit->timestamp / 1e9,
                             adc_array,
                             cable_length);

      std::string adc_string = "";

      for(auto const &adc : feb_datas.back().ADC())
        {
          adc_string += std::to_string(adc);
          adc_string += ", ";
        }
      adc_string.resize(adc_string.size() - 2);

      if(fDebug)
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

void CRTDecoder::endJob()
{
  if(fUnfoundMAC5s.size())
    {
      std::cout << "\n=================================================\n"
                << "CRT Decoder finished\n"
                << "Following MAC5s not found in channel map:\n";

      for(const unsigned int mac5 : fUnfoundMAC5s)
        std::cout << '\t' << mac5 << '\n';

      std::cout << "================================================="
                << std::endl;
    }
}

DEFINE_ART_MODULE(CRTDecoder)
