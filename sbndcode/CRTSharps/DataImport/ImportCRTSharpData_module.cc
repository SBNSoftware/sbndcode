////////////////////////////////////////////////////////////////////////
// Class:       ImportCRTSharpData
// Plugin Type: producer
// File:        ImportCRTSharpData_module.cc
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

class ImportCRTSharpData;


class ImportCRTSharpData : public art::EDProducer {
public:
  explicit ImportCRTSharpData(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ImportCRTSharpData(ImportCRTSharpData const&) = delete;
  ImportCRTSharpData(ImportCRTSharpData&&) = delete;
  ImportCRTSharpData& operator=(ImportCRTSharpData const&) = delete;
  ImportCRTSharpData& operator=(ImportCRTSharpData&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  std::vector<sbnd::crt::FEBData> FragToFEB(const artdaq::Fragment &frag);

private:

  std::vector<std::pair<unsigned, unsigned>> fMac5ToGeoIDVec;
  std::map<unsigned, unsigned>               fMac5ToGeoID;
};


ImportCRTSharpData::ImportCRTSharpData(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  produces<std::vector<sbnd::crt::FEBData>>();
  //  produces<art::Assns<artdaq::Fragment,sbnd::crt::FEBData>>();

  fMac5ToGeoIDVec = p.get<std::vector<std::pair<unsigned, unsigned>>>("FEBMac5ToGeometryIDMap");
  fMac5ToGeoID    = std::map<unsigned, unsigned>(fMac5ToGeoIDVec.begin(), fMac5ToGeoIDVec.end());
}

void ImportCRTSharpData::produce(art::Event& e)
{
  auto febDataVec      = std::make_unique<std::vector<sbnd::crt::FEBData>>();
  //  auto febDataFragAssn = std::make_unique<art::Assns<artdaq::Fragment,sbnd::crt::FEBData>>();
  
  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles = e.getMany<std::vector<artdaq::Fragment>>();

  for(auto handle : fragmentHandles)
    {
      if(!handle.isValid() || handle->size() == 0)
        continue;

      if(handle->front().type() == artdaq::Fragment::ContainerFragmentType)
        {
          for(auto cont : *handle)
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
    }
  
  e.put(std::move(febDataVec));
}

std::vector<sbnd::crt::FEBData> ImportCRTSharpData::FragToFEB(const artdaq::Fragment &frag)
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

      mf::LogInfo("ImportCRTSharpData")
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

DEFINE_ART_MODULE(ImportCRTSharpData)
