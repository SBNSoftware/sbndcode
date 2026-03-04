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

#include "artdaq-core/Data/RawEvent.hh"

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoService.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbndcode/ChannelMaps/CRT/CRTChannelMapService.h"

#include <memory>
#include <bitset>

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

  std::vector<CRTStripHit> CreateStripHits(art::Ptr<FEBData> &data, const uint32_t ref_time_s,
					   const uint32_t ref_time_ns, const uint32_t etrig_time_s,
					   const uint32_t etrig_time_ns);
  std::set<uint32_t> UnixSet(const std::vector<art::Ptr<FEBData>> &datas);

private:

  art::ServiceHandle<CRTGeoService>              fCRTGeoService;
  art::ServiceHandle<SBND::CRTChannelMapService> fCRTChannelMapService;

  std::string           fFEBDataModuleLabel;
  std::string           fFrameShiftModuleLabel;
  uint16_t              fADCThreshold;
  std::vector<double>   fErrorCoeff;
  bool                  fAllowFlag1;
  bool                  fApplyETrigWindow;
  double                fETrigMin;
  double                fETrigMax;
  bool                  fApplyTs0Window;
  double                fTs0Min;
  double                fTs0Max;
  bool                  fApplyTs1Window;
  double                fTs1Min;
  double                fTs1Max;
  bool                  fCorrectForDifferentSecond;
  bool                  fReferenceTs0;
};


sbnd::crt::CRTStripHitProducer::CRTStripHitProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fFEBDataModuleLabel(p.get<std::string>("FEBDataModuleLabel"))
  , fFrameShiftModuleLabel(p.get<std::string>("FrameShiftModuleLabel"))
  , fADCThreshold(p.get<uint16_t>("ADCThreshold"))
  , fErrorCoeff(p.get<std::vector<double>>("ErrorCoeff"))
  , fAllowFlag1(p.get<bool>("AllowFlag1"))
  , fApplyETrigWindow(p.get<bool>("ApplyETrigWindow"))
  , fETrigMin(p.get<double>("ETrigMin", 0))
  , fETrigMax(p.get<double>("ETrigMax", std::numeric_limits<double>::max()))
  , fApplyTs0Window(p.get<bool>("ApplyTs0Window"))
  , fTs0Min(p.get<double>("Ts0Min", 0))
  , fTs0Max(p.get<double>("Ts0Max", std::numeric_limits<double>::max()))
  , fApplyTs1Window(p.get<bool>("ApplyTs1Window"))
  , fTs1Min(p.get<double>("Ts1Min", 0))
  , fTs1Max(p.get<double>("Ts1Max", std::numeric_limits<double>::max()))
  , fCorrectForDifferentSecond(p.get<bool>("CorrectForDifferentSecond"))
  , fReferenceTs0(p.get<bool>("ReferenceTs0"))
{
  produces<std::vector<CRTStripHit>>();
  produces<art::Assns<FEBData, CRTStripHit>>();

  if(fReferenceTs0)
    produces<raw::TimingReferenceInfo>();
}

void sbnd::crt::CRTStripHitProducer::produce(art::Event& e)
{
  auto stripHitVec         = std::make_unique<std::vector<CRTStripHit>>();
  auto stripHitDataAssn    = std::make_unique<art::Assns<FEBData, CRTStripHit>>();
  auto timingReferenceInfo = std::make_unique<raw::TimingReferenceInfo>();

  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  
  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  uint32_t ref_time_s = 0, ref_time_ns = 0, etrig_time_s = 0, etrig_time_ns = 0;

  if(fReferenceTs0)
    {
      art::Handle<sbnd::timing::FrameShiftInfo> frameShiftHandle;
      e.getByLabel(fFrameShiftModuleLabel, frameShiftHandle);

      if(!frameShiftHandle.isValid())
        throw std::runtime_error("Frame Shift Info object is invalid, check data quality");

      uint64_t ref_time = frameShiftHandle->FrameDefault();

      ref_time_s  = ref_time / sbnd::timing::kSecondInNanoseconds;
      ref_time_ns = ref_time % sbnd::timing::kSecondInNanoseconds;

      timingReferenceInfo->timingType    = frameShiftHandle->TimingTypeDefault();
      timingReferenceInfo->timingChannel = frameShiftHandle->TimingChannelDefault();

      if(timingReferenceInfo->timingType == sbnd::timing::kNoShiftType)
	{
	  std::set<uint32_t> unix_set = UnixSet(FEBDataVec);
          ref_time_s = unix_set.size() ? *unix_set.rbegin(): 0;
	}

      uint64_t etrig_time = frameShiftHandle->FrameEtrig();

      etrig_time_s  = etrig_time / sbnd::timing::kSecondInNanoseconds;
      etrig_time_ns = etrig_time % sbnd::timing::kSecondInNanoseconds;
    }

  for(auto data : FEBDataVec)
    {
      std::vector<CRTStripHit> newStripHits = CreateStripHits(data, ref_time_s, ref_time_ns, etrig_time_s, etrig_time_ns);

      for(auto hit : newStripHits)
        {
          stripHitVec->push_back(hit);
          util::CreateAssn(*this, e, *stripHitVec, data, *stripHitDataAssn);
        }
    }

  e.put(std::move(stripHitVec));
  e.put(std::move(stripHitDataAssn));

  if(fReferenceTs0)
    e.put(std::move(timingReferenceInfo));
}

std::vector<sbnd::crt::CRTStripHit> sbnd::crt::CRTStripHitProducer::CreateStripHits(art::Ptr<FEBData> &data, const uint32_t ref_time_s,
                                                                                    const uint32_t ref_time_ns, const uint32_t etrig_time_s,
                                                                                    const uint32_t etrig_time_ns)
{
  std::vector<CRTStripHit> stripHits;

  const uint32_t offline_module_id = data->Mac5();
  uint32_t unixs                   = data->UnixS();

  // Only consider "real data" readouts, not clock resets etc
  if(!(data->Flags() == 3 || (fAllowFlag1 && data->Flags() == 1)))
    return stripHits;

  const uint32_t offline_channel_id = fCRTChannelMapService->ConstructOfflineChannelIDFromOfflineModuleIDAndOfflineLocalChannel(offline_module_id, 0);
  const CRTModuleGeo module         = fCRTGeoService->GetModule(offline_channel_id);

  // Correct for FEB readout cable length
  // (time is FEB-by-FEB not channel-by-channel)
  double t0 = (int)data->Ts0() + module.t0DelayCorrection;
  double t1 = (int)data->Ts1() + module.t1DelayCorrection;

  if(fCorrectForDifferentSecond)
    {
      // For data events we correct the t0 time used for clustering if the
      // events fell either side of the PPS reset.

      const int64_t unix_diff = static_cast<int64_t>(ref_time_s) - static_cast<int64_t>(unixs);

      if(unix_diff < -1 || unix_diff > 1)
        {
          throw std::runtime_error(Form("Unix timestamps differ by more than 1 (%li)", unix_diff));
        }

      if(unix_diff == 1)
        {
          t0 -= 1e9;
          unixs += 1;
        }
      else if(unix_diff == -1)
        {
          t0 += 1e9;
          unixs -= 1;
        }
    }

  if(fReferenceTs0)
    t0 -= ref_time_ns;

  if(fApplyETrigWindow)
    {
      double t0_etrig         = (int)data->Ts0() + module.t0DelayCorrection;
      const int64_t unix_diff = static_cast<int64_t>(ref_time_s) - static_cast<int64_t>(etrig_time_s);

      if(unix_diff < -1 || unix_diff > 1)
        {
          throw std::runtime_error(Form("Unix timestamps differ by more than 1 (%li)", unix_diff));
        }

      t0_etrig -= etrig_time_ns;

      if(unix_diff == 1)
	t0_etrig -= 1e9;
      else if(unix_diff == -1)
	t0_etrig += 1e9;

      if(t0_etrig < fETrigMin || t0_etrig > fETrigMax)
	return stripHits;
    }

  if(fApplyTs0Window && (t0 < fTs0Min || t0 > fTs0Max))
    return stripHits;

  if(fApplyTs1Window && (t1 < fTs1Min || t1 > fTs1Max))
    return stripHits;

  // Iterate via strip (2 SiPMs per strip)
  const auto &sipm_adcs = data->ADC();
  for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
    {
      // Calculate SiPM channel number
      const uint16_t offline_channel_id = fCRTChannelMapService->ConstructOfflineChannelIDFromOfflineModuleIDAndOfflineLocalChannel(offline_module_id, adc_i);

      const CRTStripGeo strip = fCRTGeoService->GetStrip(offline_channel_id);
      const CRTSiPMGeo sipm1  = fCRTGeoService->GetSiPM(offline_channel_id);
      const CRTSiPMGeo sipm2  = fCRTGeoService->GetSiPM(offline_channel_id+1);

      if(sipm1.status == CRTChannelStatus::kDeadChannel || sipm2.status == CRTChannelStatus::kDeadChannel)
        continue;

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

          stripHits.emplace_back(offline_channel_id, t0, t1, ref_time_s, pos, err, adc1, adc2);
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
