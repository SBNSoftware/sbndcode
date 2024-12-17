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
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"

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
                                           const uint32_t ref_time_ns);
  std::set<uint32_t> UnixSet(const std::vector<art::Ptr<FEBData>> &datas);
  bool SPECTDCReference(art::Event& e, const uint64_t &raw_ts, uint64_t &ref_time);
  bool PTBHLTReference(art::Event& e, const uint64_t &raw_ts, uint64_t &ref_time, uint32_t &hlt_code);
  std::bitset<32> TriggerWordBitset(uint32_t trig_word);

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
  bool                fReferenceTs0;
  int                 fTimingType;
  std::string         fDAQHeaderModuleLabel;
  std::string         fDAQHeaderInstanceLabel;
  uint32_t            fRawTSCorrection;
  uint32_t            fMaxAllowedRefTimeDiff;
  std::string         fSPECTDCModuleLabel;
  uint32_t            fSPECTDCETrigChannel;
  std::string         fPTBModuleLabel;
  std::set<uint32_t>  fAllowedPTBHLTs;
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
  , fReferenceTs0(p.get<bool>("ReferenceTs0"))
  , fTimingType(p.get<int>("TimingType", 0))
  , fDAQHeaderModuleLabel(p.get<std::string>("DAQHeaderModuleLabel", ""))
  , fDAQHeaderInstanceLabel(p.get<std::string>("DAQHeaderInstanceLabel", ""))
  , fRawTSCorrection(p.get<uint32_t>("RawTSCorrection", 0))
  , fMaxAllowedRefTimeDiff(p.get<uint32_t>("MaxAllowedRefTimeDiff", 0))
  , fSPECTDCModuleLabel(p.get<std::string>("SPECTDCModuleLabel", ""))
  , fSPECTDCETrigChannel(p.get<uint32_t>("SPECTDCETrigChannel", 4))
  , fPTBModuleLabel(p.get<std::string>("PTBModuleLabel", ""))
  , fAllowedPTBHLTs(p.get<std::vector<uint32_t>>("AllowedPTBHLTs", {}).begin(), p.get<std::vector<uint32_t>>("AllowedPTBHLTs", {}).end())
{
  produces<std::vector<CRTStripHit>>();
  produces<art::Assns<FEBData, CRTStripHit>>();
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

  uint64_t raw_ts = 0, ref_time = 0;
  uint32_t ref_time_s = 0, ref_time_ns = 0;

  if(fReferenceTs0)
    {
      art::Handle<artdaq::detail::RawEventHeader> DAQHeaderHandle;
      e.getByLabel(fDAQHeaderModuleLabel, fDAQHeaderInstanceLabel, DAQHeaderHandle);

      if(DAQHeaderHandle.isValid())
        {
          artdaq::RawEvent rawHeaderEvent = artdaq::RawEvent(*DAQHeaderHandle);
          raw_ts = rawHeaderEvent.timestamp() - fRawTSCorrection;
        }

      int timingType = fTimingType;
      int timingCh   = 0;

      if(timingType == 0)
        {
          uint64_t spec_tdc_ref_time = 0;

          if(SPECTDCReference(e, raw_ts, spec_tdc_ref_time))
            {
              ref_time = spec_tdc_ref_time;
              timingCh = fSPECTDCETrigChannel;
            }
          else
            ++timingType;
        }

      if(timingType == 1)
        {
          uint64_t ptb_hlt_ref_time = 0;
          uint32_t hlt_code         = 0;

          if(PTBHLTReference(e, raw_ts, ptb_hlt_ref_time, hlt_code))
            {
              ref_time = ptb_hlt_ref_time;
              timingCh = hlt_code;
            }
          else
            ++timingType;
        }

      if(timingType == 2)
        {
          std::set<uint32_t> unix_set = UnixSet(FEBDataVec);
          ref_time = unix_set.size() ? *unix_set.rbegin() * static_cast<uint64_t>(1e9) : 0;
        }

      ref_time_s  = ref_time / static_cast<uint64_t>(1e9);
      ref_time_ns = ref_time % static_cast<uint64_t>(1e9);

      timingReferenceInfo->timingType    = timingType;
      timingReferenceInfo->timingChannel = timingCh;
    }

  for(auto data : FEBDataVec)
    {
      std::vector<CRTStripHit> newStripHits = CreateStripHits(data, ref_time_s, ref_time_ns);

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
                                                                                    const uint32_t ref_time_ns)
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

      const int64_t unix_diff = static_cast<int64_t>(ref_time_s) - static_cast<int64_t>(unixs);

      if(unix_diff < -1 || unix_diff > 1)
        {
          throw std::runtime_error("Unix timestamps differ by more than 1" + unix_diff);
        }

      if(unix_diff == 1)
        {
          t0 -= static_cast<int>(1e9);
          unixs += 1;
        }
      else if(unix_diff == -1)
        {
          t0 += static_cast<int>(1e9);
          unixs -= 1;
        }
    }

  if(fReferenceTs0)
    t0 -= ref_time_ns;

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

          stripHits.emplace_back(channel, t0, t1, ref_time_s, pos, err, adc1, adc2);
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

bool sbnd::crt::CRTStripHitProducer::SPECTDCReference(art::Event& e, const uint64_t &raw_ts, uint64_t &ref_time)
{
  bool found = false;

  art::Handle<std::vector<sbnd::timing::DAQTimestamp>> TDCHandle;
  e.getByLabel(fSPECTDCModuleLabel, TDCHandle);

  if(!TDCHandle.isValid() || TDCHandle->size() == 0)
    return found;

  std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> TDCVec;
  art::fill_ptr_vector(TDCVec, TDCHandle);

  int64_t min_diff     = std::numeric_limits<int64_t>::max();
  uint64_t min_diff_ts = 0;

  for(auto ts : TDCVec)
    {
      if(ts->Channel() == fSPECTDCETrigChannel)
        {
          int64_t diff = raw_ts > ts->Timestamp() ? raw_ts - ts->Timestamp() : ts->Timestamp() - raw_ts;

          if(diff < min_diff)
            {
              min_diff    = diff;
              min_diff_ts = ts->Timestamp();
              found       = true;
            }
        }
    }

  if(min_diff > fMaxAllowedRefTimeDiff)
    return false;

  ref_time = min_diff_ts;

  return found;
}

bool sbnd::crt::CRTStripHitProducer::PTBHLTReference(art::Event& e, const uint64_t &raw_ts, uint64_t &ref_time, uint32_t &hlt_code)
{
  bool found = false;

  art::Handle<std::vector<raw::ptb::sbndptb>> PTBHandle;
  e.getByLabel(fPTBModuleLabel, PTBHandle);

  if(!PTBHandle.isValid() || PTBHandle->size() == 0)
    return found;

  std::vector<art::Ptr<raw::ptb::sbndptb>> PTBVec;
  art::fill_ptr_vector(PTBVec, PTBHandle);

  int64_t min_diff     = std::numeric_limits<int64_t>::max();
  uint64_t min_diff_ts = 0;

  for(auto ptb : PTBVec)
    {
      for(auto hlt : ptb->GetHLTriggers())
        {
          uint64_t hlt_timestamp          = (hlt.timestamp * 20);
          std::bitset<32> hlt_word_bitset = TriggerWordBitset(hlt.trigger_word);

          for(uint32_t allowed_hlt : fAllowedPTBHLTs)
            {
              if(hlt_word_bitset[allowed_hlt])
                {
                  int64_t diff = raw_ts > hlt_timestamp ? raw_ts - hlt_timestamp : hlt_timestamp - raw_ts;

                  if(diff < min_diff)
                    {
                      min_diff    = diff;
                      min_diff_ts = hlt_timestamp;
                      hlt_code    = allowed_hlt;
                      found       = true;
                    }
                }
            }
        }
    }

  if(min_diff > fMaxAllowedRefTimeDiff)
    return false;

  ref_time = min_diff_ts;

  return found;
}

std::bitset<32> sbnd::crt::CRTStripHitProducer::TriggerWordBitset(uint32_t trig_word)
{
  uint32_t trig_word_dec;
  std::stringstream ss;

  ss << std::hex << trig_word << std::dec;
  ss >> trig_word_dec;

  return std::bitset<32>(trig_word_dec);
}

DEFINE_ART_MODULE(sbnd::crt::CRTStripHitProducer)
