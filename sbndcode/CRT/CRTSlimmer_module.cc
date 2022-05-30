/**
 * \brief Makes sbnd::crt::CRTData from sbnd::crt::FEBData
 *
 * \details This class loops over sbnd::crt::FEBData objects and only selects
 * strips which have both SiPM ADC values above a certain (configurable) threshold.
 * The two SiPMs in the selected strips are then saved in sbnd::crt::CRTData objects.
 *
 * \author Marco Del Tutto
 */

////////////////////////////////////////////////////////////////////////
// Class:       CRTSlimmer
// Plugin Type: producer (Unknown Unknown)
// File:        CRTSlimmer_module.cc
//
// Generated at Thu Mar 24 11:27:06 2022 by Marco Del Tutto using cetskelgen
// from  version .
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"


#include <memory>
#include <math.h>

#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"

namespace sbnd {
  namespace crt {
    class CRTSlimmer;
  }
}


class sbnd::crt::CRTSlimmer : public art::EDProducer {
public:
  explicit CRTSlimmer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTSlimmer(CRTSlimmer const&) = delete;
  CRTSlimmer(CRTSlimmer&&) = delete;
  CRTSlimmer& operator=(CRTSlimmer const&) = delete;
  CRTSlimmer& operator=(CRTSlimmer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string _feb_data_producer;
  uint16_t _adc_threshold;

};


sbnd::crt::CRTSlimmer::CRTSlimmer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , _feb_data_producer(p.get<std::string>("FEBDataProducer"))
  , _adc_threshold(p.get<uint16_t>("ADCThreshold"))
{
  produces<std::vector<sbnd::crt::CRTData>>();
  produces<art::Assns<sbnd::crt::CRTData, sbnd::crt::FEBData>>();
  produces<art::Assns<sbnd::crt::CRTData, sim::AuxDetIDE>>();

  consumes<std::vector<sbnd::crt::FEBData>>(_feb_data_producer);
  consumes<art::Assns<sbnd::crt::FEBData, sim::AuxDetIDE, sbnd::crt::FEBTruthInfo>>(_feb_data_producer);
}

void sbnd::crt::CRTSlimmer::produce(art::Event& e)
{
  // Implementation of required member function here.
  std::unique_ptr<std::vector<sbnd::crt::CRTData>> crt_data_v(new std::vector<sbnd::crt::CRTData>);
  std::unique_ptr<art::Assns<sbnd::crt::CRTData, sbnd::crt::FEBData>> crtdata_to_febdata_assns
      (new art::Assns<sbnd::crt::CRTData, sbnd::crt::FEBData>);
  std::unique_ptr<art::Assns<sbnd::crt::CRTData, sim::AuxDetIDE>> crtdata_to_ide_assns
      (new art::Assns<sbnd::crt::CRTData, sim::AuxDetIDE>);

  art::Handle<std::vector<sbnd::crt::FEBData>> feb_data_h;
  e.getByLabel(_feb_data_producer, feb_data_h);

  // make sure hits look good
  if (!feb_data_h.isValid()) {
    throw art::Exception(art::errors::Configuration) << "could not locate FEBData." << std::endl;;
  }

  std::vector<art::Ptr<sbnd::crt::FEBData>> feb_data_v;
  art::fill_ptr_vector(feb_data_v, feb_data_h);

  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> febdata_to_ides (feb_data_h, e, _feb_data_producer);

  art::PtrMaker<sbnd::crt::CRTData> makeDataPtr(e);

  for (size_t feb_i = 0; feb_i < feb_data_v.size(); feb_i++) {

    auto const feb_data = feb_data_v[feb_i];
    mf::LogDebug("CRTSlimmer") << "FEB " << feb_i << " with mac " << feb_data->Mac5() << std::endl;

    auto ides = febdata_to_ides.at(feb_data.key());
    mf::LogDebug("CRTSlimmer") << "We have " << ides.size() << " IDEs." << std::endl;

    // Construct a map to go from SiPM ID to the index of the IDE
    // FEBTruthInfo stores, for each AuxDetIDE, a vector containing
    // the SiPM IDs contributing to that energy deposit
    std::map<int, std::vector<int>> sipm_to_ideids;
    for (size_t j = 0; j < 32; j++) sipm_to_ideids[j] = std::vector<int>();
    for (size_t ide_i = 0; ide_i < febdata_to_ides.data(feb_i).size(); ide_i++) {
      const sbnd::crt::FEBTruthInfo *fti = febdata_to_ides.data(feb_i)[ide_i];
      sipm_to_ideids[fti->GetChannel()].push_back(ide_i);
      mf::LogDebug("CRTSlimmer") << "ide_i " << ide_i << " ene " << ides[ide_i]->energyDeposited << " fti->GetChannel() " << fti->GetChannel() << std::endl;
    }

    auto adcs = feb_data->ADC();

    for (size_t i = 0; i < adcs.size(); i+=2) {

      // Both sipms need to be above threshold to save them
      if (not (adcs[i] >= _adc_threshold and adcs[i+1] >= _adc_threshold)) {
        continue;
      }

      for (size_t sipm = 0; sipm < 2; sipm++) {
        sbnd::crt::CRTData crt_data = sbnd::crt::CRTData(feb_data->Mac5() * 32 + i + sipm,
                                                         feb_data->Ts0(),
                                                         feb_data->Ts1(),
                                                         adcs[i+sipm]);

        mf::LogDebug("CRTSlimmer") << "Adding SiPM " << i + sipm << " with mac " << feb_data->Mac5()
                                   << " mapped to channel " << crt_data.Channel() << std::endl;

        crt_data_v->emplace_back(std::move(crt_data));

        art::Ptr<sbnd::crt::CRTData> crt_data_p = makeDataPtr(crt_data_v->size() - 1);

        // Create the association between CRTData and FEBData
        crtdata_to_febdata_assns->addSingle(crt_data_p, feb_data);

        // Create the association between CRTData and AuxDetIDEs
        auto & ide_ids = sipm_to_ideids[i];
        for (auto ide_id: ide_ids) {
          auto ide = ides[ide_id];
          mf::LogDebug("CRTSlimmer") << "Adding IDE with ID " << ide_id << std::endl;
          crtdata_to_ide_assns->addSingle(crt_data_p, ide);
        }
      }
    }
  }

  mf::LogInfo("CRTSlimmer") << "Creating " << crt_data_v->size()
                            << " CRTData products." << std::endl;

  e.put(std::move(crt_data_v));
  e.put(std::move(crtdata_to_febdata_assns));
  e.put(std::move(crtdata_to_ide_assns));
}



DEFINE_ART_MODULE(sbnd::crt::CRTSlimmer)



