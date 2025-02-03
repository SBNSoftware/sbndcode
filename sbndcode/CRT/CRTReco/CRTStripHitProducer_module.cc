////////////////////////////////////////////////////////////////////////
// Class:       CRTStripHitProducer
// Plugin Type: producer
// File:        CRTStripHitProducer_module.cc
//
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

#include <memory>

namespace sbnd::crt {
  class CRTStripHitProducer;
}


class sbnd::crt::CRTStripHitProducer : public art::EDProducer {
public:
  explicit CRTStripHitProducer(fhicl::ParameterSet const& p);

  CRTStripHitProducer(CRTStripHitProducer const&) = delete;
  CRTStripHitProducer(CRTStripHitProducer&&) = delete;
  CRTStripHitProducer& operator=(CRTStripHitProducer const&) = delete;
  CRTStripHitProducer& operator=(CRTStripHitProducer&&) = delete;

  void produce(art::Event& e) override;

  std::vector<CRTStripHit> CreateStripHits(art::Ptr<FEBData> &data, const uint32_t ref_unix,
                                           const uint64_t etrig_time);
  std::set<uint32_t> UnixSet(const std::vector<art::Ptr<FEBData>> &datas);

private:

  CRTGeoAlg           fCRTGeoAlg;
  std::string         fFEBDataModuleLabel;
  uint16_t            fADCThreshold;
  std::vector<double> fErrorCoeff;
  bool                fAllowFlag1;
  bool                fApplyTs1Window;
  int64_t             fTs1Min;
  int64_t             fTs1Max;
  bool                fCorrectForDifferentSecond;
  bool                fReferenceTs0ToETrig;
  std::string         fSPECTDCModuleLabel;
  uint32_t            fSPECTDCETrigChannel;
};


sbnd::crt::CRTStripHitProducer::CRTStripHitProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
  , fFEBDataModuleLabel(p.get<std::string>("FEBDataModuleLabel"))
  , fADCThreshold(p.get<uint16_t>("ADCThreshold"))
  , fErrorCoeff(p.get<std::vector<double>>("ErrorCoeff"))
  , fAllowFlag1(p.get<bool>("AllowFlag1"))
  , fApplyTs1Window(p.get<bool>("ApplyTs1Window"))
  , fTs1Min(p.get<int64_t>("Ts1Min", 0))
  , fTs1Max(p.get<int64_t>("Ts1Max", std::numeric_limits<int64_t>::max()))
  , fCorrectForDifferentSecond(p.get<bool>("CorrectForDifferentSecond"))
  , fReferenceTs0ToETrig(p.get<bool>("ReferenceTs0ToETrig"))
  , fSPECTDCModuleLabel(p.get<std::string>("SPECTDCModuleLabel", ""))
  , fSPECTDCETrigChannel(p.get<uint32_t>("SPECTDCETrigChannel", 4))
{
  produces<std::vector<CRTStripHit>>();
  produces<art::Assns<FEBData, CRTStripHit>>();
}

void sbnd::crt::CRTStripHitProducer::produce(art::Event& e)
{
  auto stripHitVec      = std::make_unique<std::vector<CRTStripHit>>();
  auto stripHitDataAssn = std::make_unique<art::Assns<FEBData, CRTStripHit>>();
  
  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  
  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // We get a set of the unix times of all the FEBDatas. In data they should only ever take up to 2 values.
  // This comes if the PPS reset arrives mid-event. In the strip hit creation method we correct for this.
  // These lines define our prescription for doing so. We take the later unix time as the reference point.
  // The actual correction should always be turned off for simulation (by setting CorrectForDifferentSecond to false).
  std::set<uint32_t> unix_set = UnixSet(FEBDataVec);
  const uint32_t ref_unix = unix_set.size() ? *unix_set.rbegin() : 0;

  uint64_t etrig_time = 0;

  if(fReferenceTs0ToETrig)
    {
      art::Handle<std::vector<sbnd::timing::DAQTimestamp>> TDCHandle;
      e.getByLabel(fSPECTDCModuleLabel, TDCHandle);

      if(TDCHandle.isValid() && TDCHandle->size() != 0)
        {
          std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> TDCVec;
          art::fill_ptr_vector(TDCVec, TDCHandle);

          for(auto ts : TDCVec)
            {
              if(ts->Channel() == fSPECTDCETrigChannel)
                {
                  etrig_time = ts->Timestamp() % static_cast<uint64_t>(1e9);

                  if(fCorrectForDifferentSecond)
                    {
                      const uint64_t etrig_unix = ts->Timestamp() / static_cast<uint64_t>(1e9);
                      const int64_t unix_diff = ref_unix - etrig_unix;
                      if(unix_diff != 0 && unix_diff != 1)
                        throw std::runtime_error(Form("Unix timestamps differ by more than 1 (%li)", unix_diff));

                      const bool previous_second = unix_diff == 1;

                      if(previous_second)
                        etrig_time -= static_cast<int>(1e9);
                    }
                }
            }
        }
    }

  for(auto data : FEBDataVec)
    {
      std::vector<CRTStripHit> newStripHits = CreateStripHits(data, ref_unix, etrig_time);

      for(auto hit : newStripHits)
        {
          stripHitVec->push_back(hit);
          util::CreateAssn(*this, e, *stripHitVec, data, *stripHitDataAssn);
        }
    }

  e.put(std::move(stripHitVec));
  e.put(std::move(stripHitDataAssn));
}

std::vector<sbnd::crt::CRTStripHit> sbnd::crt::CRTStripHitProducer::CreateStripHits(art::Ptr<FEBData> &data, const uint32_t ref_unix,
                                                                                    const uint64_t etrig_time)
{
  std::vector<CRTStripHit> stripHits;

  const uint32_t mac5  = data->Mac5();
  uint32_t unixs       = data->UnixS();

  // Only consider "real data" readouts, not clock resets etc
  if(!(data->Flags() == 3 || (fAllowFlag1 && data->Flags() == 1)))
    return stripHits;
  
  const CRTModuleGeo module = fCRTGeoAlg.GetModule(mac5 * 32);

  // Correct for FEB readout cable length
  // (time is FEB-by-FEB not channel-by-channel)
  int64_t t0 = (int)data->Ts0() + (int)module.t0CableDelayCorrection;
  int64_t t1 = (int)data->Ts1() + (int)module.t1CableDelayCorrection;

  if(fCorrectForDifferentSecond)
    {
      // For data events we correct the t0 time used for clustering if the
      // events fell either side of the PPS reset.

      const int64_t unix_diff = ref_unix - unixs;
      if(unix_diff != 0 && unix_diff != 1)
        throw std::runtime_error(Form("Unix timestamps differ by more than 1 (%li)", unix_diff));

      const bool previous_second = unix_diff == 1;

      if(previous_second)
        {
          t0 -= static_cast<int>(1e9);
          unixs += 1;
        }
    }

  if(fReferenceTs0ToETrig)
    t0 -= etrig_time;

  if(fApplyTs1Window && (t1 < fTs1Min || t1 > fTs1Max))
    return stripHits;

  // Iterate via strip (2 SiPMs per strip)
  const auto &sipm_adcs = data->ADC();
  for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
    {
      // Calculate SiPM channel number
      const uint16_t channel = mac5 * 32 + adc_i;

      const CRTStripGeo strip = fCRTGeoAlg.GetStrip(channel);
      const CRTSiPMGeo sipm1  = fCRTGeoAlg.GetSiPM(channel);
      const CRTSiPMGeo sipm2  = fCRTGeoAlg.GetSiPM(channel+1);

      // Subtract channel pedestals
      const uint16_t adc1 = sipm1.pedestal < sipm_adcs[adc_i]   ? sipm_adcs[adc_i] - sipm1.pedestal   : 0;
      const uint16_t adc2 = sipm2.pedestal < sipm_adcs[adc_i+1] ? sipm_adcs[adc_i+1] - sipm2.pedestal : 0;

      // Keep hit if both SiPMs above threshold
      if(adc1 > fADCThreshold && adc2 > fADCThreshold)
        {
          // Access width of strip from the geometry algorithm
          const double width = strip.width;

          // Use light ratio to infer lateral position
          const double pos = width / 2. * tanh(log(1. * adc2/adc1)) + width / 2.;
          double err       = fErrorCoeff[0] * width + fErrorCoeff[1] * pos + fErrorCoeff[2] * pos * pos;

          // Ensure error does not allow positions outside the strip geometry
          if(pos + err > width)
            err = width - pos;

          if(pos - err < 0)
            err = pos;

          stripHits.emplace_back(channel, t0, t1, ref_unix, pos, err, adc1, adc2);
        }
    }

  return stripHits;
}

std::set<uint32_t> sbnd::crt::CRTStripHitProducer::UnixSet(const std::vector<art::Ptr<FEBData>> &datas)
{
  std::set<uint32_t> set;

  for(auto const& data : datas)
    {
      set.insert(data->UnixS());
    }

  return set;
}

DEFINE_ART_MODULE(sbnd::crt::CRTStripHitProducer)
