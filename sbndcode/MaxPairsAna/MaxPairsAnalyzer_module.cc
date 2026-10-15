// Art analyzer: software maxPairs from PMT waveforms + PTB summary,
// plus MTCA analog mon-pulse for comparison (CASBAna-style).

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr size_t kHoldLength = 48;

const std::vector<int> kPair1 = {
    7,   9,   11,  13,  15,  17,  37,  39,  41,  61,  63,  67,  69,  71,
    85,  87,  89,  91,  93,  95,  115, 117, 119, 139, 141, 145, 147, 149,
    163, 165, 169, 171, 173, 193, 195, 197, 217, 219, 221, 223, 225, 227,
    241, 243, 247, 249, 251, 271, 273, 275, 295, 297, 299, 301, 303, 305};
const std::vector<int> kPair2 = {
    6,   8,   10,  12,  14,  16,  36,  38,  40,  60,  62,  66,  68,  70,
    84,  86,  88,  90,  92,  94,  114, 116, 118, 138, 140, 144, 146, 148,
    162, 164, 168, 170, 172, 192, 194, 196, 216, 218, 220, 222, 224, 226,
    240, 242, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304};
const std::vector<int> kUnpaired = {65, 64, 143, 142, 167, 166, 245, 244};

std::vector<bool> ConstructBinaryResponse(const raw::OpDetWaveform& wvf,
                                          int monThreshold) {
  const auto& wfData = wvf.Waveform();
  std::vector<bool> binaryResponse(wfData.size(), false);
  size_t holdCounter = 0;

  for (size_t t = 0; t < wfData.size(); ++t) {
    if (wfData[t] < monThreshold) {
      binaryResponse[t] = true;
      holdCounter = kHoldLength;
    }
    if (holdCounter > 0) {
      binaryResponse[t] = true;
      --holdCounter;
    }
  }
  return binaryResponse;
}

// Modal non-empty waveform length in the event (0 if none).
size_t ModalWaveformLength(const std::vector<raw::OpDetWaveform>& waveforms) {
  std::map<size_t, int> counts;
  for (const auto& wf : waveforms) {
    if (!wf.empty()) {
      ++counts[wf.size()];
    }
  }
  size_t modal = 0;
  int best = -1;
  for (const auto& kv : counts) {
    if (kv.second > best) {
      best = kv.second;
      modal = kv.first;
    }
  }
  return modal;
}

// Pick one fragment per channel:
//  1) closest TimeStamp to refTs (if hasRefTs)
//  2) else prefer modal length, then deepest min ADC
// Optional TimeStamp window (us) restricts which fragments are considered.
// Keeps CASBAna pair OR / hold / channel lists; drops last-wins + rigid size skip.
std::vector<int> ConstructMonPulse(const std::vector<raw::OpDetWaveform>& waveforms,
                                   int monThreshold,
                                   bool hasRefTs,
                                   double refTs,
                                   bool restrictTs,
                                   double tsMin,
                                   double tsMax) {
  if (waveforms.empty()) {
    return {};
  }

  const size_t modalLen = ModalWaveformLength(waveforms);

  std::map<int, size_t> channelToIndex;
  for (size_t idx = 0; idx < waveforms.size(); ++idx) {
    const auto& wf = waveforms[idx];
    if (wf.empty()) {
      continue;
    }
    if (restrictTs) {
      const double ts = wf.TimeStamp();
      if (ts < tsMin || ts > tsMax) {
        continue;
      }
    }
    const int chan = static_cast<int>(wf.ChannelNumber());
    const int minADC =
        static_cast<int>(*std::min_element(wf.Waveform().begin(), wf.Waveform().end()));
    const bool modalMatch = (modalLen > 0 && wf.size() == modalLen);

    auto it = channelToIndex.find(chan);
    if (it == channelToIndex.end()) {
      channelToIndex[chan] = idx;
      continue;
    }

    const auto& prev = waveforms[it->second];
    const int prevMin =
        static_cast<int>(*std::min_element(prev.Waveform().begin(), prev.Waveform().end()));
    const bool prevModal = (modalLen > 0 && prev.size() == modalLen);

    bool takeNew = false;
    if (hasRefTs) {
      const double dNew = std::abs(wf.TimeStamp() - refTs);
      const double dPrev = std::abs(prev.TimeStamp() - refTs);
      if (dNew < dPrev) {
        takeNew = true;
      } else if (dNew == dPrev) {
        // Same time: prefer modal length, then deeper pulse.
        if (modalMatch && !prevModal) {
          takeNew = true;
        } else if (modalMatch == prevModal && minADC < prevMin) {
          takeNew = true;
        }
      }
    } else {
      if (modalMatch && !prevModal) {
        takeNew = true;
      } else if (modalMatch == prevModal && minADC < prevMin) {
        takeNew = true;
      } else if (modalMatch == prevModal && minADC == prevMin &&
                 wf.size() > prev.size()) {
        takeNew = true;
      }
    }

    if (takeNew) {
      channelToIndex[chan] = idx;
    }
  }

  if (channelToIndex.empty()) {
    return {};
  }

  // Output length from selected fragments (no rigid first-waveform size skip).
  size_t wfSize = 0;
  for (const auto& kv : channelToIndex) {
    wfSize = std::max(wfSize, waveforms[kv.second].size());
  }
  if (wfSize == 0) {
    return {};
  }

  std::vector<int> pairsAbove(wfSize, 0);
  std::vector<bool> binTemp(wfSize, false);

  auto orBinaryInto = [&](const raw::OpDetWaveform& wf) {
    if (wf.empty()) {
      return;
    }
    auto bin = ConstructBinaryResponse(wf, monThreshold);
    const size_t n = std::min(bin.size(), wfSize);
    for (size_t t = 0; t < n; ++t) {
      binTemp[t] = binTemp[t] || bin[t];
    }
  };

  for (size_t i = 0; i < kPair1.size(); ++i) {
    auto it1 = channelToIndex.find(kPair1[i]);
    auto it2 = channelToIndex.find(kPair2[i]);
    if (it1 == channelToIndex.end() && it2 == channelToIndex.end()) {
      continue;
    }

    std::fill(binTemp.begin(), binTemp.end(), false);

    if (it1 != channelToIndex.end()) {
      orBinaryInto(waveforms[it1->second]);
    }
    if (it2 != channelToIndex.end()) {
      orBinaryInto(waveforms[it2->second]);
    }

    for (size_t t = 0; t < wfSize; ++t) {
      if (binTemp[t]) {
        ++pairsAbove[t];
      }
    }
  }

  for (int chan : kUnpaired) {
    auto it = channelToIndex.find(chan);
    if (it == channelToIndex.end()) {
      continue;
    }
    std::fill(binTemp.begin(), binTemp.end(), false);
    orBinaryInto(waveforms[it->second]);
    for (size_t t = 0; t < wfSize; ++t) {
      if (binTemp[t]) {
        ++pairsAbove[t];
      }
    }
  }

  return pairsAbove;
}

// Mode of ADC samples (CASBAna baseline estimator).
int CalculateMode(const raw::OpDetWaveform& waveform) {
  const auto& wf = waveform.Waveform();
  if (wf.empty()) {
    return 0;
  }

  std::map<uint16_t, int> counts;
  for (const auto& val : wf) {
    counts[val]++;
  }

  int maxCount = 0;
  uint16_t mode = wf[0];
  for (const auto& pair : counts) {
    if (pair.second > maxCount) {
      maxCount = pair.second;
      mode = pair.first;
    }
  }
  return static_cast<int>(mode);
}

constexpr double kSampleUs = 0.002;  // 2 ns PMT samples

// Local maxima in digital mon-pulse with height >= minHeight.
std::vector<size_t> FindLocalMaxima(const std::vector<int>& mon, int minHeight) {
  std::vector<size_t> peaks;
  if (mon.size() < 3) {
    return peaks;
  }
  for (size_t i = 1; i + 1 < mon.size(); ++i) {
    if (mon[i] < minHeight) {
      continue;
    }
    if (mon[i] > mon[i - 1] && mon[i] >= mon[i + 1]) {
      peaks.push_back(i);
    }
  }
  return peaks;
}

}  // namespace

class MaxPairsAnalyzer : public art::EDAnalyzer {
public:
  explicit MaxPairsAnalyzer(fhicl::ParameterSet const& p);
  void beginJob() override;
  void analyze(art::Event const& e) override;
  void endJob() override;

private:
  art::InputTag fPMTTag;
  art::InputTag fPTBTag;
  art::InputTag fTimingTag;
  int fPairThreshold;
  int fMTCAChannelIndex;
  double fMTCA1NHITAmplitude;
  int fDigMonPulseNHitCut;  // fill ADC hist only where digital mon-pulse > this
  double fADCPlotXMin;      // x-axis zoom for MTCA ADC canvas/PNG
  double fADCPlotXMax;
  double fRawPMTPlotYMin;   // raw PMT debug y limits; ignored if min >= max (auto)
  double fRawPMTPlotYMax;
  int fRawPMTDebugIndex; // 1-based index among mon-pulse debug events (0 = off)
  int fDebugEvent;     // -1 = any event; >=0 = only that event number
  int fNDebugEvents;   // save first N HLT2 digital/MTCA overlays as PNG
  int fNDebugIneffEvents;      // maxPairs > cut, no HLT2 overlays
  int fDebugIneffMaxPairsMin;  // min maxPairs for inefficiency debug sample
  int fHLTEffBit;      // HLT bit for efficiency (2 = BNB light+beam)
  double fHLT2SinceHLT1XMin;  // x-axis for hHLT2SinceHLT1 (PTB ticks)
  double fHLT2SinceHLT1XMax;
  int fHLT2SinceHLT1NBins;
  double fMaxPairsSampleXMin;  // x-axis for hMaxPairsSample
  double fMaxPairsSampleXMax;
  int fMaxPairsSampleNBins;
  double fFlashTsXMin;  // per-event debug flash-TS axis (us)
  double fFlashTsXMax;
  int fFlashTsNBins;
  double fFlashTsJobXMin;  // job-level peak-time hist axis (us)
  double fFlashTsJobXMax;
  int fFlashTsJobNBins;
  double fEffPeakTsMin;  // efficiency: peak time window (us vs ETRIG/ZB)
  double fEffPeakTsMax;
  int fPeakMinHeight;  // min digital pairs to count as a local-max peak
  int fNFlashTsDebugEvents;       // per-event flash-TS PNGs (HLT2)
  int fNFlashTsDebugIneffEvents;  // flash-TS PNGs (inefficiency sample)
  bool fSavePlotPNGs;
  std::string fPlotOutputPrefix;

  TTree* fTree = nullptr;
  TH1I* fHMaxPairs = nullptr;
  TH2D* fHMaxPairsVsNHitMTCA = nullptr;
  TH1D* fHMaxPairsMinusNHitMTCA = nullptr;
  TH1D* fHBaselineSubMTCAADC = nullptr;
  TH1D* fHEffDenomMaxPairs = nullptr;
  TH1D* fHEffNumerMaxPairs = nullptr;
  TH1D* fHEffDenomNHitMTCA = nullptr;
  TH1D* fHEffNumerNHitMTCA = nullptr;
  TH1I* fHHLTBits = nullptr;
  TH1D* fHHLT2SinceHLT1 = nullptr;  // HLT2 ts - HLT1 ts (events with HLT2)
  TH1D* fHMaxPairsSample = nullptr;  // sample index of max digital mon-pulse
  TH1D* fHPeakTsRelZB = nullptr;         // all local-max peak times vs ZB
  TH1D* fHNPeaks = nullptr;        // max local-max peak height per event (any time)
  TH1D* fHNPeaksWindow = nullptr;  // brightest maxPairs in EffPeakTs window
  int fNDebugPlotsSaved = 0;
  int fNDebugIneffPlotsSaved = 0;
  int fNFlashTsDebugSaved = 0;
  int fNFlashTsIneffDebugSaved = 0;
  bool fRawPMTDebugSaved = false;

  int fRun = 0;
  int fSubRun = 0;
  int fEvent = 0;
  int fMaxPairs = 0;
  int fMaxPairsSample = -1;  // sample index of max digital mon-pulse
  int fNWaveforms = 0;
  int fNUniquePmtFlashes = 0;
  int fNUniqueMsumFlashes = 0;
  int fHasPTB = 0;
  int fHasMTCA = 0;
  int fHasInWindowFlash = 0;  // 1 if maxPairs>0 peak found in EffPeakTs window
  int fHasHLTBit = 0;  // 1 if fHLTEffBit set in any HLT word
  int fNHLT = 0;
  int fNLLT = 0;
  int fNPeaks = 0;       // local-max peaks (digital mon-pulse) in the event
  int fMaxPeakPairs = 0; // brightest local-max peak height in the event
  double fAmpMTCA = 0.;
  double fNHitMTCA = 0.;
  double fDiffMaxPairsNHitMTCA = 0.;
  double fSelectedFlashTs = 0.;  // OpDetWaveform TimeStamp of fragment with peak
  double fSelectedPeakTs = 0.;   // absolute peak time (us vs ETRIG/ZB)
  ULong64_t fHLTWordOr = 0;
  ULong64_t fLLTWordOr = 0;
  ULong64_t fHLTTsMin = 0;
  ULong64_t fLLTTsMin = 0;
  std::vector<ULong64_t> fHLTWords;
  std::vector<ULong64_t> fHLTTimestamps;
  std::vector<ULong64_t> fLLTWords;
  std::vector<ULong64_t> fLLTTimestamps;
  std::vector<int> fDigitalMonPulse;
  std::vector<double> fAnalogMTCANhit;
  std::vector<double> fPmtFlashTs;   // unique PMT flash TimeStamps (us)
  std::vector<double> fMsumFlashTs;  // unique MSUM flash TimeStamps (us)
};

MaxPairsAnalyzer::MaxPairsAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer(p)
  , fPMTTag(p.get<std::string>("PMTInputTag", "pmtdecoder:PMTChannels:PMTDecoder"))
  , fPTBTag(p.get<std::string>("PTBModuleLabel", "ptbdecoder"))
  , fTimingTag(p.get<std::string>("TimingInputTag",
                                  "pmtdecoder:TimingChannels:TriggerEffDecoder"))
  , fPairThreshold(p.get<int>("PairThreshold", 14230))
  , fMTCAChannelIndex(p.get<int>("MTCAChannelIndex", 3))
  , fMTCA1NHITAmplitude(p.get<double>("MTCA1NHITAmplitude", 14.0))
  , fDigMonPulseNHitCut(p.get<int>("DigMonPulseNHitCut", 3))
  , fADCPlotXMin(p.get<double>("ADCPlotXMin", -10.))
  , fADCPlotXMax(p.get<double>("ADCPlotXMax", 200.))
  , fRawPMTPlotYMin(p.get<double>("RawPMTPlotYMin", -1.))
  , fRawPMTPlotYMax(p.get<double>("RawPMTPlotYMax", -1.))
  , fRawPMTDebugIndex(p.get<int>("RawPMTDebugIndex", 1))
  , fDebugEvent(p.get<int>("DebugEvent", -1))
  , fNDebugEvents(p.get<int>("NDebugEvents", 10))
  , fNDebugIneffEvents(p.get<int>("NDebugIneffEvents", 10))
  , fDebugIneffMaxPairsMin(p.get<int>("DebugIneffMaxPairsMin", 10))
  , fHLTEffBit(p.get<int>("HLTEffBit", 2))
  , fHLT2SinceHLT1XMin(p.get<double>("HLT2SinceHLT1XMin", -50.))
  , fHLT2SinceHLT1XMax(p.get<double>("HLT2SinceHLT1XMax", 500.))
  , fHLT2SinceHLT1NBins(p.get<int>("HLT2SinceHLT1NBins", 550))
  , fMaxPairsSampleXMin(p.get<double>("MaxPairsSampleXMin", 0.))
  , fMaxPairsSampleXMax(p.get<double>("MaxPairsSampleXMax", 7000.))
  , fMaxPairsSampleNBins(p.get<int>("MaxPairsSampleNBins", 700))
  , fFlashTsXMin(p.get<double>("FlashTsXMin", -10.))
  , fFlashTsXMax(p.get<double>("FlashTsXMax", 10.))
  , fFlashTsNBins(p.get<int>("FlashTsNBins", 40))
  , fFlashTsJobXMin(p.get<double>("FlashTsJobXMin", -2000.))
  , fFlashTsJobXMax(p.get<double>("FlashTsJobXMax", 2000.))
  , fFlashTsJobNBins(p.get<int>("FlashTsJobNBins", 400))
  , fEffPeakTsMin(p.get<double>("EffPeakTsMin", 0.))
  , fEffPeakTsMax(p.get<double>("EffPeakTsMax", 3.))
  , fPeakMinHeight(p.get<int>("PeakMinHeight", 1))
  , fNFlashTsDebugEvents(p.get<int>("NFlashTsDebugEvents", 10))
  , fNFlashTsDebugIneffEvents(p.get<int>("NFlashTsDebugIneffEvents", 10))
  , fSavePlotPNGs(p.get<bool>("SavePlotPNGs", true))
  , fPlotOutputPrefix(p.get<std::string>("PlotOutputPrefix", "maxpairs"))
{}

void MaxPairsAnalyzer::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("maxpairs", "maxPairs + PTB + MTCA analog mon-pulse");
  fTree->Branch("run", &fRun, "run/I");
  fTree->Branch("subrun", &fSubRun, "subrun/I");
  fTree->Branch("event", &fEvent, "event/I");
  fTree->Branch("maxPairs", &fMaxPairs, "maxPairs/I");
  fTree->Branch("maxPairsSample", &fMaxPairsSample, "maxPairsSample/I");
  fTree->Branch("nWaveforms", &fNWaveforms, "nWaveforms/I");
  fTree->Branch("nUniquePmtFlashes", &fNUniquePmtFlashes, "nUniquePmtFlashes/I");
  fTree->Branch("nUniqueMsumFlashes", &fNUniqueMsumFlashes,
                 "nUniqueMsumFlashes/I");
  fTree->Branch("hasPTB", &fHasPTB, "hasPTB/I");
  fTree->Branch("hasMTCA", &fHasMTCA, "hasMTCA/I");
  fTree->Branch("hasInWindowFlash", &fHasInWindowFlash, "hasInWindowFlash/I");
  fTree->Branch("hasHLTBit", &fHasHLTBit, "hasHLTBit/I");
  fTree->Branch("nPeaks", &fNPeaks, "nPeaks/I");
  fTree->Branch("maxPeakPairs", &fMaxPeakPairs, "maxPeakPairs/I");
  fTree->Branch("amp_MTCA", &fAmpMTCA, "amp_MTCA/D");
  fTree->Branch("nHIT_MTCA", &fNHitMTCA, "nHIT_MTCA/D");
  fTree->Branch("diff_maxPairs_nHIT_MTCA", &fDiffMaxPairsNHitMTCA,
                 "diff_maxPairs_nHIT_MTCA/D");
  fTree->Branch("selectedFlashTs", &fSelectedFlashTs, "selectedFlashTs/D");
  fTree->Branch("selectedPeakTs", &fSelectedPeakTs, "selectedPeakTs/D");
  fTree->Branch("nHLT", &fNHLT, "nHLT/I");
  fTree->Branch("nLLT", &fNLLT, "nLLT/I");
  fTree->Branch("hlt_word_or", &fHLTWordOr, "hlt_word_or/l");
  fTree->Branch("llt_word_or", &fLLTWordOr, "llt_word_or/l");
  fTree->Branch("hlt_ts_min", &fHLTTsMin, "hlt_ts_min/l");
  fTree->Branch("llt_ts_min", &fLLTTsMin, "llt_ts_min/l");
  fTree->Branch("hlt_words", &fHLTWords);
  fTree->Branch("hlt_timestamps", &fHLTTimestamps);
  fTree->Branch("llt_words", &fLLTWords);
  fTree->Branch("llt_timestamps", &fLLTTimestamps);
  fTree->Branch("digital_monpulse", &fDigitalMonPulse);
  fTree->Branch("analog_mtca_nhit", &fAnalogMTCANhit);
  fTree->Branch("pmt_flash_ts", &fPmtFlashTs);
  fTree->Branch("msum_flash_ts", &fMsumFlashTs);

  fHMaxPairs = tfs->make<TH1I>("hMaxPairs", ";maxPairs;Events", 70, 0, 70);
  fHMaxPairsVsNHitMTCA = tfs->make<TH2D>(
      "hMaxPairsVsNHitMTCA",
      "Software maxPairs vs MTCA analog NHIT;nHIT_{MTCA};maxPairs",
      80, 0., 40., 70, 0., 70.);
  fHMaxPairsMinusNHitMTCA = tfs->make<TH1D>(
      "hMaxPairsMinusNHitMTCA",
      "maxPairs - nHIT_{MTCA};maxPairs - nHIT_{MTCA};Events", 80, -20., 20.);

  // CASBAna-style ADC range; fill only samples with digital mon-pulse > cut.
  constexpr int adcRangeMin = -10000;
  constexpr int adcRangeMax = 10000;
  constexpr int nBinsADC = adcRangeMax - adcRangeMin;
  {
    const std::string title =
        "Baseline-Subtracted MTCA ADC (digital mon-pulse > " +
        std::to_string(fDigMonPulseNHitCut) + ");Amplitude (ADC);Count";
    fHBaselineSubMTCAADC = tfs->make<TH1D>(
        "hBaselineSubtractedMTCA_ADC_digCut", title.c_str(), nBinsADC,
        adcRangeMin, adcRangeMax);
  }
  fHBaselineSubMTCAADC->SetLineColor(kGray + 1);
  fHBaselineSubMTCAADC->SetLineWidth(2);

  // Trigger efficiency numer/denom (HLT bit set / in-window flashes in bin).
  const std::string effBitLabel = "HLT bit " + std::to_string(fHLTEffBit);
  fHEffDenomMaxPairs = tfs->make<TH1D>(
      "hEffDenomMaxPairs", ";maxPairs;Events", 70, 0., 70.);
  {
    const std::string numerTitle =
        ";maxPairs;Events with " + effBitLabel;
    fHEffNumerMaxPairs = tfs->make<TH1D>(
        "hEffNumerMaxPairs", numerTitle.c_str(), 70, 0., 70.);
  }
  fHEffDenomNHitMTCA = tfs->make<TH1D>(
      "hEffDenomNHitMTCA", ";nHIT_{MTCA};Events", 80, 0., 40.);
  {
    const std::string numerTitle =
        ";nHIT_{MTCA};Events with " + effBitLabel;
    fHEffNumerNHitMTCA = tfs->make<TH1D>(
        "hEffNumerNHitMTCA", numerTitle.c_str(), 80, 0., 40.);
  }
  fHEffDenomMaxPairs->Sumw2();
  fHEffNumerMaxPairs->Sumw2();
  fHEffDenomNHitMTCA->Sumw2();
  fHEffNumerNHitMTCA->Sumw2();

  fHHLTBits = tfs->make<TH1I>(
      "hHLTBits", "HLT bits fired (all events);HLT bit number;Count", 32, 0,
      32);
  fHHLTBits->SetLineColor(kBlack);
  fHHLTBits->SetLineWidth(2);

  // Beam-accept timing: for events with >=1 HLT2, HLT2 time relative to HLT1 (ZB).
  {
    int nBins = fHLT2SinceHLT1NBins;
    if (nBins < 1) {
      nBins = 1;
    }
    double xMin = fHLT2SinceHLT1XMin;
    double xMax = fHLT2SinceHLT1XMax;
    if (!(xMax > xMin)) {
      xMin = -50.;
      xMax = 500.;
    }
    fHHLT2SinceHLT1 = tfs->make<TH1D>(
        "hHLT2SinceHLT1",
        "HLT2 time since HLT1 (events with HLT2);"
        "t_{HLT2} - t_{HLT1} (PTB ticks);HLT2 triggers",
        nBins, xMin, xMax);
  }
  fHHLT2SinceHLT1->SetLineColor(kBlue + 1);
  fHHLT2SinceHLT1->SetLineWidth(2);

  {
    int nBins = fMaxPairsSampleNBins;
    if (nBins < 1) {
      nBins = 1;
    }
    double xMin = fMaxPairsSampleXMin;
    double xMax = fMaxPairsSampleXMax;
    if (!(xMax > xMin)) {
      xMin = 0.;
      xMax = 7000.;
    }
    fHMaxPairsSample = tfs->make<TH1D>(
        "hMaxPairsSample",
        "Sample of max digital mon-pulse;sample index;Events", nBins, xMin,
        xMax);
  }
  fHMaxPairsSample->SetLineColor(kBlue + 1);
  fHMaxPairsSample->SetLineWidth(2);

  // Job-level: all local-max peak times vs event ref (ETRIG ~ ZB at 0).
  {
    int nBins = fFlashTsJobNBins;
    if (nBins < 1) {
      nBins = 1;
    }
    double xMin = fFlashTsJobXMin;
    double xMax = fFlashTsJobXMax;
    if (!(xMax > xMin)) {
      xMin = -2000.;
      xMax = 2000.;
    }
    fHPeakTsRelZB = tfs->make<TH1D>(
        "hPeakTsRelZB",
        "All digital mon-pulse local-max peak times vs ZB;"
        "t_{peak} - t_{ETRIG/ZB} (#mus);Peaks",
        nBins, xMin, xMax);
  }
  fHPeakTsRelZB->SetLineColor(kMagenta + 1);
  fHPeakTsRelZB->SetLineWidth(2);

  fHNPeaks = tfs->make<TH1D>(
      "hNPeaks",
      "Max local-max peak height per event (any time);"
      "maxPairs of brightest peak;Events",
      70, 0., 70.);
  fHNPeaksWindow = tfs->make<TH1D>(
      "hNPeaksWindow",
      "Brightest maxPairs in efficiency peak window;"
      "maxPairs in EffPeakTs window;Events",
      70, 0., 70.);
}

void MaxPairsAnalyzer::analyze(art::Event const& e) {
  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.event();

  fHLTWords.clear();
  fHLTTimestamps.clear();
  fLLTWords.clear();
  fLLTTimestamps.clear();
  fDigitalMonPulse.clear();
  fAnalogMTCANhit.clear();
  fPmtFlashTs.clear();
  fMsumFlashTs.clear();
  fHLTWordOr = 0;
  fLLTWordOr = 0;
  fHLTTsMin = 0;
  fLLTTsMin = 0;
  fNHLT = 0;
  fNLLT = 0;
  fHasPTB = 0;
  fHasMTCA = 0;
  fHasInWindowFlash = 0;
  fHasHLTBit = 0;
  fAmpMTCA = 0.;
  fNHitMTCA = 0.;
  fDiffMaxPairsNHitMTCA = 0.;
  fSelectedFlashTs = 0.;
  fSelectedPeakTs = 0.;
  fMaxPairs = -1;
  fMaxPairsSample = -1;
  fNPeaks = 0;
  fMaxPeakPairs = 0;
  fNWaveforms = 0;
  fNUniquePmtFlashes = 0;
  fNUniqueMsumFlashes = 0;

  art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle;
  if (e.getByLabel(fPTBTag, ptbHandle) && ptbHandle.isValid()) {
    fHasPTB = 1;
    ULong64_t hltMin = std::numeric_limits<ULong64_t>::max();
    ULong64_t lltMin = std::numeric_limits<ULong64_t>::max();

    for (auto const& ptb : *ptbHandle) {
      for (size_t j = 0; j < ptb.GetNHLTriggers(); ++j) {
        auto const& trig = ptb.GetHLTrigger(j);
        ULong64_t word = trig.trigger_word & 0xFFFFFFFFull;
        fHLTWords.push_back(word);
        fHLTTimestamps.push_back(trig.timestamp);
        fHLTWordOr |= word;
        hltMin = std::min(hltMin, static_cast<ULong64_t>(trig.timestamp));
      }
      for (size_t j = 0; j < ptb.GetNLLTriggers(); ++j) {
        auto const& trig = ptb.GetLLTrigger(j);
        ULong64_t word = trig.trigger_word & 0xFFFFFFFFull;
        fLLTWords.push_back(word);
        fLLTTimestamps.push_back(trig.timestamp);
        fLLTWordOr |= word;
        lltMin = std::min(lltMin, static_cast<ULong64_t>(trig.timestamp));
      }
    }
    fNHLT = static_cast<int>(fHLTWords.size());
    fNLLT = static_cast<int>(fLLTWords.size());
    if (fNHLT > 0) {
      fHLTTsMin = hltMin;
    }
    if (fNLLT > 0) {
      fLLTTsMin = lltMin;
    }
  }
  if (fHLTEffBit >= 0 && fHLTEffBit < 64 &&
      (fHLTWordOr & (1ULL << fHLTEffBit))) {
    fHasHLTBit = 1;
  }
  // Debug: count which HLT bits appear (once per event via OR of all HLT words).
  if (fHHLTBits && fHLTWordOr != 0) {
    for (int b = 0; b < 32; ++b) {
      if (fHLTWordOr & (1ULL << b)) {
        fHHLTBits->Fill(b);
      }
    }
  }

  // HLT2 times relative to HLT1 (ZB), for beam-accept window studies.
  // Require at least one HLT2 in the event; for each HLT2, use the nearest
  // HLT1 at or before that HLT2 (else nearest HLT1 overall).
  constexpr int kHLT1Bit = 1;  // Zero Bias
  constexpr int kHLT2Bit = 2;  // Light + Beam
  if (fHHLT2SinceHLT1 && (fHLTWordOr & (1ULL << kHLT2Bit))) {
    std::vector<ULong64_t> hlt1Ts;
    hlt1Ts.reserve(fHLTTimestamps.size());
    for (size_t i = 0; i < fHLTWords.size(); ++i) {
      if (fHLTWords[i] & (1ULL << kHLT1Bit)) {
        hlt1Ts.push_back(fHLTTimestamps[i]);
      }
    }
    if (!hlt1Ts.empty()) {
      std::sort(hlt1Ts.begin(), hlt1Ts.end());
      for (size_t i = 0; i < fHLTWords.size(); ++i) {
        if (!(fHLTWords[i] & (1ULL << kHLT2Bit))) {
          continue;
        }
        const ULong64_t t2 = fHLTTimestamps[i];
        // Nearest HLT1 with t1 <= t2; else earliest HLT1.
        auto it = std::upper_bound(hlt1Ts.begin(), hlt1Ts.end(), t2);
        ULong64_t t1 = hlt1Ts.front();
        if (it != hlt1Ts.begin()) {
          t1 = *(--it);
        }
        const double dt = static_cast<double>(t2) - static_cast<double>(t1);
        fHHLT2SinceHLT1->Fill(dt);
      }
    }
  }

  // Timing / MSUM: list flash starts; MSUM fragment chosen later to match the
  // PMT digitizer window that owns the brightest peak in EffPeakTs [min,max].
  art::Handle<std::vector<raw::OpDetWaveform>> timingHandle;
  bool hasTiming = e.getByLabel(fTimingTag, timingHandle) && timingHandle.isValid();
  int msumChannel = -1;
  const raw::OpDetWaveform* selectedMsum = nullptr;

  if (hasTiming &&
      timingHandle->size() > static_cast<size_t>(fMTCAChannelIndex) &&
      !(*timingHandle)[fMTCAChannelIndex].empty()) {
    msumChannel =
        static_cast<int>((*timingHandle)[fMTCAChannelIndex].ChannelNumber());
  }

  if (hasTiming) {
    std::set<double> msumTs;
    for (auto const& wf : *timingHandle) {
      if (wf.empty()) {
        continue;
      }
      if (msumChannel >= 0 &&
          static_cast<int>(wf.ChannelNumber()) != msumChannel) {
        continue;
      }
      if (msumChannel < 0) {
        msumChannel = static_cast<int>(wf.ChannelNumber());
      }
      msumTs.insert(wf.TimeStamp());
    }
    fMsumFlashTs.assign(msumTs.begin(), msumTs.end());
    fNUniqueMsumFlashes = static_cast<int>(fMsumFlashTs.size());
  }

  // Digital mon-pulse: scan every unique PMT flash start. Efficiency uses the
  // brightest maxPairs sample whose absolute peak time is in EffPeakTs window.
  // No cut on digitizer window start time.
  art::Handle<std::vector<raw::OpDetWaveform>> pmtHandle;
  bool hasPMT = e.getByLabel(fPMTTag, pmtHandle) && pmtHandle.isValid() &&
                !pmtHandle->empty();
  if (hasPMT) {
    fNWaveforms = static_cast<int>(pmtHandle->size());
    std::set<double> pmtTs;
    for (auto const& wf : *pmtHandle) {
      if (wf.empty()) {
        continue;
      }
      pmtTs.insert(wf.TimeStamp());
    }
    fPmtFlashTs.assign(pmtTs.begin(), pmtTs.end());
    fNUniquePmtFlashes = static_cast<int>(fPmtFlashTs.size());

    int bestInWindowPairs = -1;
    int bestInWindowSample = -1;
    double bestInWindowPeakTs = 0.;
    double bestInWindowFlashTs = 0.;
    std::vector<int> bestInWindowMon;
    std::vector<std::pair<double, int>> allPeaks;  // (peakTs, maxPairs)

    for (double ts : fPmtFlashTs) {
      std::vector<int> mon =
          ConstructMonPulse(*pmtHandle, fPairThreshold, false, 0., true, ts, ts);
      if (mon.empty()) {
        continue;
      }

      // Efficiency: brightest sample whose absolute time is in EffPeakTs window.
      // Typical waveform ~14 us; only search fragments that can overlap.
      const bool canOverlapPeakWindow =
          !(ts > fEffPeakTsMax || ts < fEffPeakTsMin - 20.);
      if (canOverlapPeakWindow) {
        for (size_t i = 0; i < mon.size(); ++i) {
          const double peakTs = ts + kSampleUs * static_cast<double>(i);
          if (peakTs < fEffPeakTsMin || peakTs > fEffPeakTsMax) {
            continue;
          }
          if (mon[i] > bestInWindowPairs) {
            bestInWindowPairs = mon[i];
            bestInWindowSample = static_cast<int>(i);
            bestInWindowPeakTs = peakTs;
            bestInWindowFlashTs = ts;
            bestInWindowMon = mon;
          }
        }
      }

      // Bias diagnostic: local maxima in every digitizer window.
      for (size_t ip : FindLocalMaxima(mon, fPeakMinHeight)) {
        const double peakTs = ts + kSampleUs * static_cast<double>(ip);
        allPeaks.emplace_back(peakTs, mon[ip]);
      }
    }

    fNPeaks = static_cast<int>(allPeaks.size());
    fMaxPeakPairs = 0;
    for (const auto& pk : allPeaks) {
      if (pk.second > fMaxPeakPairs) {
        fMaxPeakPairs = pk.second;
      }
      if (fHPeakTsRelZB) {
        fHPeakTsRelZB->Fill(pk.first);
      }
    }
    if (fHNPeaks) {
      fHNPeaks->Fill(fMaxPeakPairs);
    }

    // Require maxPairs > 0 in the peak window for efficiency sample.
    if (bestInWindowPairs > 0) {
      fHasInWindowFlash = 1;
      fMaxPairs = bestInWindowPairs;
      fMaxPairsSample = bestInWindowSample;
      fSelectedFlashTs = bestInWindowFlashTs;
      fSelectedPeakTs = bestInWindowPeakTs;
      fDigitalMonPulse = std::move(bestInWindowMon);
      fHMaxPairs->Fill(fMaxPairs);
      if (fHMaxPairsSample && fMaxPairsSample >= 0) {
        fHMaxPairsSample->Fill(fMaxPairsSample);
      }
      if (fHNPeaksWindow) {
        fHNPeaksWindow->Fill(fMaxPairs);
      }
    } else {
      fDigitalMonPulse.clear();
    }
  }

  // Analog MTCA / MSUM from the digitizer window that owns the selected peak.
  if (fHasInWindowFlash && hasTiming) {
    double bestDts = 1e300;
    for (auto const& wf : *timingHandle) {
      if (wf.empty()) {
        continue;
      }
      if (msumChannel >= 0 &&
          static_cast<int>(wf.ChannelNumber()) != msumChannel) {
        continue;
      }
      const double dts = std::abs(wf.TimeStamp() - fSelectedFlashTs);
      if (dts < bestDts) {
        bestDts = dts;
        selectedMsum = &wf;
      }
    }
  }

  if (selectedMsum && !selectedMsum->empty() && fHasInWindowFlash) {
    fHasMTCA = 1;
    const auto& wfMTCA = *selectedMsum;
    const auto& wfData = wfMTCA.Waveform();
    const int baseline = CalculateMode(wfMTCA);
    fAnalogMTCANhit.resize(wfData.size());
    for (size_t t = 0; t < wfData.size(); ++t) {
      fAnalogMTCANhit[t] =
          static_cast<double>(baseline - static_cast<int>(wfData[t])) /
          std::max(fMTCA1NHITAmplitude, 1e-9);
    }
    // Analog peak lags the digital mon-pulse: take max NHIT among samples
    // whose absolute time is in the efficiency peak window (same cut as maxPairs).
    const double msumTs = wfMTCA.TimeStamp();
    const double invNh = 1.0 / std::max(fMTCA1NHITAmplitude, 1e-9);
    bool foundInWindow = false;
    double bestAmp = 0.;
    for (size_t t = 0; t < wfData.size(); ++t) {
      const double peakTs = msumTs + kSampleUs * static_cast<double>(t);
      if (peakTs < fEffPeakTsMin || peakTs > fEffPeakTsMax) {
        continue;
      }
      const double amp =
          static_cast<double>(baseline - static_cast<int>(wfData[t]));
      if (!foundInWindow || amp > bestAmp) {
        bestAmp = amp;
        foundInWindow = true;
      }
    }
    if (foundInWindow) {
      fAmpMTCA = bestAmp;
      fNHitMTCA = bestAmp * invNh;
    } else {
      const int minADC = *std::min_element(wfData.begin(), wfData.end());
      fAmpMTCA = static_cast<double>(baseline - minADC);
      fNHitMTCA = fAmpMTCA * invNh;
    }
    if (!fDigitalMonPulse.empty()) {
      const size_t nCommon =
          std::min(fDigitalMonPulse.size(), wfData.size());
      for (size_t t = 0; t < nCommon; ++t) {
        if (fDigitalMonPulse[t] > fDigMonPulseNHitCut) {
          fHBaselineSubMTCAADC->Fill(
              static_cast<double>(baseline - static_cast<int>(wfData[t])));
        }
      }
    }
    if (fMaxPairs >= 0) {
      fDiffMaxPairsNHitMTCA = static_cast<double>(fMaxPairs) - fNHitMTCA;
      fHMaxPairsVsNHitMTCA->Fill(fNHitMTCA, fMaxPairs);
      fHMaxPairsMinusNHitMTCA->Fill(fDiffMaxPairsNHitMTCA);
    }
  }

  // Per-event flash TimeStamp vs ZB/ETRIG (0 at center).
  // Two samples: HLT2, and inefficiency (peak in [EffPeakTs], maxPairs > cut, no HLT2).
  const bool isHLT2 = (fHLTWordOr & (1ULL << kHLT2Bit)) != 0;
  const bool isIneffCandidate =
      fHasInWindowFlash && fMaxPairs > fDebugIneffMaxPairsMin && !isHLT2;
  const bool debugEventOk =
      (fDebugEvent < 0) || (fEvent == fDebugEvent);
  const bool wantFlashTsHLT2 =
      hasPMT && !fPmtFlashTs.empty() && isHLT2 && debugEventOk &&
      fNFlashTsDebugSaved < fNFlashTsDebugEvents;
  const bool wantFlashTsIneff =
      hasPMT && !fPmtFlashTs.empty() && isIneffCandidate && debugEventOk &&
      fNFlashTsIneffDebugSaved < fNFlashTsDebugIneffEvents;
  const bool wantFlashTsDebug = wantFlashTsHLT2 || wantFlashTsIneff;
  if (wantFlashTsDebug) {
    art::ServiceHandle<art::TFileService> tfs;
    const bool isIneffDebug = wantFlashTsIneff && !wantFlashTsHLT2;
    const std::string sampleTag = isIneffDebug ? "ineff_" : "";
    const std::string tag =
        sampleTag + "r" + std::to_string(fRun) + "e" + std::to_string(fEvent) +
        "s" + std::to_string(fSubRun);
    const std::string titleBase =
        std::string(isIneffDebug ? "INEFF debug" : "Debug") + " event " +
        std::to_string(fRun) + ":" + std::to_string(fSubRun) + ":" +
        std::to_string(fEvent) +
        (isIneffDebug ? (" maxPairs=" + std::to_string(fMaxPairs) + " no HLT2")
                      : "");

    int nBins = fFlashTsNBins;
    if (nBins < 1) {
      nBins = 1;
    }
    double xMin = fFlashTsXMin;
    double xMax = fFlashTsXMax;
    if (!(xMax > xMin)) {
      xMin = -100.;
      xMax = 100.;
    }

    // HLT2 times relative to HLT1, in us (PTB tick = 20 ns).
    constexpr double kPtbTickUs = 0.02;
    std::vector<double> hlt2RelUs;
    {
      std::vector<ULong64_t> hlt1Ts;
      hlt1Ts.reserve(fHLTTimestamps.size());
      for (size_t i = 0; i < fHLTWords.size(); ++i) {
        if (fHLTWords[i] & (1ULL << kHLT1Bit)) {
          hlt1Ts.push_back(fHLTTimestamps[i]);
        }
      }
      std::sort(hlt1Ts.begin(), hlt1Ts.end());
      if (!hlt1Ts.empty()) {
        for (size_t i = 0; i < fHLTWords.size(); ++i) {
          if (!(fHLTWords[i] & (1ULL << kHLT2Bit))) {
            continue;
          }
          const ULong64_t t2 = fHLTTimestamps[i];
          auto it = std::upper_bound(hlt1Ts.begin(), hlt1Ts.end(), t2);
          ULong64_t t1 = hlt1Ts.front();
          if (it != hlt1Ts.begin()) {
            t1 = *(--it);
          }
          hlt2RelUs.push_back(kPtbTickUs * (static_cast<double>(t2) -
                                             static_cast<double>(t1)));
        }
      }
    }

    TH1D* hPmt = tfs->make<TH1D>(
        ("hDebugPmtFlashTsRelZB_" + tag).c_str(),
        (titleBase +
         ";t_{flash start} - t_{ETRIG/ZB} (#mus);Flashes")
            .c_str(),
        nBins, xMin, xMax);
    TH1D* hMsum = tfs->make<TH1D>(
        ("hDebugMsumFlashTsRelZB_" + tag).c_str(),
        (titleBase +
         ";t_{flash start} - t_{ETRIG/ZB} (#mus);Flashes")
            .c_str(),
        nBins, xMin, xMax);
    hPmt->SetLineColor(kBlue + 1);
    hPmt->SetLineWidth(2);
    hPmt->SetMarkerColor(kBlue + 1);
    hPmt->SetMarkerStyle(20);
    hPmt->SetMarkerSize(1.1);
    hMsum->SetLineColor(kRed + 1);
    hMsum->SetLineWidth(2);
    hMsum->SetFillStyle(0);

    for (double t : fPmtFlashTs) {
      hPmt->Fill(t);
    }
    for (double t : fMsumFlashTs) {
      hMsum->Fill(t);
    }

    TCanvas* cTs = tfs->make<TCanvas>(
        ("cDebugFlashTsRelZB_" + tag).c_str(),
        (titleBase + " flash TimeStamp vs ZB/ETRIG").c_str(), 900, 600);
    cTs->cd();
    const double ymax = std::max(hPmt->GetMaximum(), hMsum->GetMaximum());
    hMsum->SetMaximum(std::max(1.0, ymax) * 1.45);
    hMsum->SetMinimum(0.);
    hMsum->SetTitle((titleBase + " flash starts vs ZB/ETRIG (0);"
                                 "t - t_{ETRIG/ZB} (#mus);Flashes / markers")
                        .c_str());
    // MSUM first; PMT with markers on top (same bin otherwise hides blue).
    hMsum->Draw("HIST");
    hPmt->Draw("HIST P SAME");

    const double yLine = hMsum->GetMaximum();
    TLine* lineZB = new TLine(0., 0., 0., yLine);
    lineZB->SetLineColor(kBlack);
    lineZB->SetLineStyle(1);
    lineZB->SetLineWidth(2);
    lineZB->Draw();

    // Efficiency flash-peak window [EffPeakTsMin, EffPeakTsMax].
    TLine* lineEffLo =
        new TLine(fEffPeakTsMin, 0., fEffPeakTsMin, yLine);
    lineEffLo->SetLineColor(kOrange + 7);
    lineEffLo->SetLineStyle(2);
    lineEffLo->SetLineWidth(2);
    lineEffLo->Draw();
    TLine* lineEffHi =
        new TLine(fEffPeakTsMax, 0., fEffPeakTsMax, yLine);
    lineEffHi->SetLineColor(kOrange + 7);
    lineEffHi->SetLineStyle(2);
    lineEffHi->SetLineWidth(2);
    lineEffHi->Draw();

    TLine* lineHlt2Leg = nullptr;
    for (double tHlt2 : hlt2RelUs) {
      TLine* lineHlt2 = new TLine(tHlt2, 0., tHlt2, yLine);
      lineHlt2->SetLineColor(kCyan + 2);
      lineHlt2->SetLineStyle(1);
      lineHlt2->SetLineWidth(3);
      lineHlt2->Draw();
      if (!lineHlt2Leg) {
        lineHlt2Leg = lineHlt2;
      }
    }

    auto inRange = [&](double t) { return t >= xMin && t <= xMax; };
    TLine* lineSel = nullptr;
    TLine* linePeak = nullptr;
    if (fHasInWindowFlash && inRange(fSelectedFlashTs)) {
      lineSel = new TLine(fSelectedFlashTs, 0., fSelectedFlashTs, yLine);
      lineSel->SetLineColor(kGreen + 2);
      lineSel->SetLineStyle(3);
      lineSel->SetLineWidth(2);
      lineSel->Draw();
    }
    if (fHasInWindowFlash && inRange(fSelectedPeakTs)) {
      linePeak = new TLine(fSelectedPeakTs, 0., fSelectedPeakTs, yLine);
      linePeak->SetLineColor(kMagenta + 1);
      linePeak->SetLineStyle(3);
      linePeak->SetLineWidth(2);
      linePeak->Draw();
    }

    auto* legTs = new TLegend(0.55, 0.62, 0.88, 0.88);
    legTs->AddEntry(hPmt, "PMT flash starts", "lp");
    legTs->AddEntry(hMsum, "MSUM flash starts", "l");
    legTs->AddEntry(lineZB, "ZB/ETRIG (0)", "l");
    legTs->AddEntry(lineEffLo, "Eff peak window", "l");
    if (lineHlt2Leg) {
      legTs->AddEntry(lineHlt2Leg, "HLT2 (vs HLT1)", "l");
    }
    if (lineSel) {
      legTs->AddEntry(lineSel, "Selected flash start", "l");
    }
    if (linePeak) {
      legTs->AddEntry(linePeak, "Selected flash peak", "l");
    }
    legTs->Draw();

    cTs->Write();
    if (fSavePlotPNGs) {
      const std::string png =
          fPlotOutputPrefix + "_debug_flash_ts_" + tag + ".png";
      cTs->SaveAs(png.c_str());
      if (isIneffDebug) {
        std::cout << "MaxPairsAnalyzer: wrote " << png << " (ineff "
                  << (fNFlashTsIneffDebugSaved + 1) << "/"
                  << fNFlashTsDebugIneffEvents << ") maxPairs=" << fMaxPairs
                  << " nPmtFlashes=" << fNUniquePmtFlashes
                  << " nMsumFlashes=" << fNUniqueMsumFlashes
                  << " selectedTs=" << fSelectedFlashTs
                  << " peakTs=" << fSelectedPeakTs << std::endl;
      } else {
        std::cout << "MaxPairsAnalyzer: wrote " << png << " ("
                  << (fNFlashTsDebugSaved + 1) << "/" << fNFlashTsDebugEvents
                  << ") nPmtFlashes=" << fNUniquePmtFlashes
                  << " nMsumFlashes=" << fNUniqueMsumFlashes
                  << " selectedTs=" << fSelectedFlashTs
                  << " peakTs=" << fSelectedPeakTs
                  << " nHLT2=" << hlt2RelUs.size();
        for (double tHlt2 : hlt2RelUs) {
          std::cout << " hlt2_us=" << tHlt2;
        }
        std::cout << std::endl;
      }
    }
    if (isIneffDebug) {
      ++fNFlashTsIneffDebugSaved;
    } else {
      ++fNFlashTsDebugSaved;
    }
  }

  // Debug overlays: HLT2 sample and inefficiency sample (maxPairs > cut, no HLT2).
  const bool monPulseReady =
      hasPMT && fHasMTCA && !fDigitalMonPulse.empty() &&
      !fAnalogMTCANhit.empty() && debugEventOk;
  const bool wantDebugHLT2 =
      monPulseReady && isHLT2 && fNDebugPlotsSaved < fNDebugEvents;
  const bool wantDebugIneff =
      monPulseReady && isIneffCandidate &&
      fNDebugIneffPlotsSaved < fNDebugIneffEvents;
  const bool wantDebug = wantDebugHLT2 || wantDebugIneff;
  if (wantDebug) {
    art::ServiceHandle<art::TFileService> tfs;
    const bool isIneffDebug = wantDebugIneff && !wantDebugHLT2;
    const size_t nDig = fDigitalMonPulse.size();
    const size_t nAna = fAnalogMTCANhit.size();
    const size_t nBins = std::max(nDig, nAna);
    const std::string sampleTag = isIneffDebug ? "ineff_" : "";
    const std::string tag =
        sampleTag + "r" + std::to_string(fRun) + "e" + std::to_string(fEvent) +
        "s" + std::to_string(fSubRun);
    const std::string titleBase =
        std::string(isIneffDebug ? "INEFF debug" : "Debug") + " event " +
        std::to_string(fRun) + ":" + std::to_string(fSubRun) + ":" +
        std::to_string(fEvent) +
        (isIneffDebug ? (" maxPairs=" + std::to_string(fMaxPairs) + " no HLT2")
                      : "");

    TH1D* hDig = tfs->make<TH1D>(
        ("hDebugDigitalMonPulse_" + tag).c_str(),
        (titleBase + " digital mon-pulse;sample;pairs").c_str(),
        static_cast<int>(nBins), 0., static_cast<double>(nBins));
    TH1D* hAna = tfs->make<TH1D>(
        ("hDebugAnalogMTCA_" + tag).c_str(),
        (titleBase + " MTCA analog (NHIT);sample;NHIT").c_str(),
        static_cast<int>(nBins), 0., static_cast<double>(nBins));
    hDig->SetLineColor(kBlue);
    hAna->SetLineColor(kRed);

    for (size_t t = 0; t < nDig; ++t) {
      hDig->SetBinContent(static_cast<int>(t) + 1, fDigitalMonPulse[t]);
    }
    for (size_t t = 0; t < nAna; ++t) {
      hAna->SetBinContent(static_cast<int>(t) + 1, fAnalogMTCANhit[t]);
    }

    TCanvas* cDebug = tfs->make<TCanvas>(
        ("cDebugMonPulseOverlay_" + tag).c_str(),
        (titleBase + " digital vs MTCA analog").c_str(), 900, 600);
    cDebug->cd();
    const double ymax = std::max(hDig->GetMaximum(), hAna->GetMaximum());
    hDig->SetMaximum(ymax * 1.2);
    hDig->SetTitle((titleBase + ";sample;pairs / NHIT").c_str());
    hDig->Draw("HIST");
    hAna->Draw("HIST SAME");
    auto* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(hDig, "digital mon-pulse", "l");
    leg->AddEntry(hAna, "MTCA analog (NHIT)", "l");
    leg->Draw();
    cDebug->Write();
    if (fSavePlotPNGs) {
      const std::string png =
          fPlotOutputPrefix + "_debug_monpulse_" + tag + ".png";
      cDebug->SaveAs(png.c_str());
      if (isIneffDebug) {
        std::cout << "MaxPairsAnalyzer: wrote " << png << " (ineff "
                  << (fNDebugIneffPlotsSaved + 1) << "/" << fNDebugIneffEvents
                  << ") maxPairs=" << fMaxPairs << std::endl;
      } else {
        std::cout << "MaxPairsAnalyzer: wrote " << png << " ("
                  << (fNDebugPlotsSaved + 1) << "/" << fNDebugEvents << ")"
                  << std::endl;
      }
    }

    if (isIneffDebug) {
      ++fNDebugIneffPlotsSaved;
    } else {
      ++fNDebugPlotsSaved;
    }
  }

  // Raw PMT overlay for the Nth HLT2 mon-pulse debug event (1-based RawPMTDebugIndex).
  const bool wantRawPMT =
      !fRawPMTDebugSaved && hasPMT && pmtHandle.isValid() &&
      fRawPMTDebugIndex > 0 && wantDebugHLT2 &&
      fNDebugPlotsSaved == fRawPMTDebugIndex;
  if (wantRawPMT) {
    art::ServiceHandle<art::TFileService> tfs;
    const std::string tag =
        "r" + std::to_string(fRun) + "e" + std::to_string(fEvent) + "s" +
        std::to_string(fSubRun);
    const std::string titleBase =
        "Debug event " + std::to_string(fRun) + ":" + std::to_string(fSubRun) +
        ":" + std::to_string(fEvent);

    size_t nBinsRaw = 0;
    double yMin = 1e300;
    double yMax = -1e300;
    for (auto const& wf : *pmtHandle) {
      if (wf.empty()) {
        continue;
      }
      nBinsRaw = std::max(nBinsRaw, wf.size());
      for (auto adc : wf.Waveform()) {
        yMin = std::min(yMin, static_cast<double>(adc));
        yMax = std::max(yMax, static_cast<double>(adc));
      }
    }

    if (nBinsRaw > 0) {
      TCanvas* cRaw = tfs->make<TCanvas>(
          ("cDebugRawPMTWaveforms_" + tag).c_str(),
          (titleBase + " raw PMT waveforms").c_str(), 1200, 700);
      cRaw->cd();

      TH1D* hFrame = tfs->make<TH1D>(
          ("hDebugRawPMTFrame_" + tag).c_str(),
          (titleBase + " raw PMT waveforms;sample;ADC").c_str(),
          static_cast<int>(nBinsRaw), 0., static_cast<double>(nBinsRaw));
      hFrame->SetStats(0);
      if (fRawPMTPlotYMin < fRawPMTPlotYMax) {
        hFrame->SetMinimum(fRawPMTPlotYMin);
        hFrame->SetMaximum(fRawPMTPlotYMax);
      } else {
        yMin = std::min(yMin, static_cast<double>(fPairThreshold));
        yMax = std::max(yMax, static_cast<double>(fPairThreshold));
        const double pad = 0.05 * (yMax - yMin + 1.);
        hFrame->SetMinimum(yMin - pad);
        hFrame->SetMaximum(yMax + pad);
      }
      hFrame->Draw("AXIS");

      int colorIdx = 0;
      size_t wfIdx = 0;
      for (auto const& wf : *pmtHandle) {
        const auto& data = wf.Waveform();
        if (data.empty()) {
          ++wfIdx;
          continue;
        }
        TH1D* hCh = tfs->make<TH1D>(
            ("hDebugRawPMT_i" + std::to_string(wfIdx) + "_ch" +
             std::to_string(wf.ChannelNumber()) + "_" + tag)
                .c_str(),
            ";sample;ADC", static_cast<int>(data.size()), 0.,
            static_cast<double>(data.size()));
        hCh->SetStats(0);
        hCh->SetLineColor(2 + (colorIdx % 8));
        hCh->SetLineWidth(1);
        for (size_t t = 0; t < data.size(); ++t) {
          hCh->SetBinContent(static_cast<int>(t) + 1, data[t]);
        }
        hCh->Draw("HIST SAME");
        ++colorIdx;
        ++wfIdx;
      }

      TLine* threshLine = new TLine(0., fPairThreshold,
                                    static_cast<double>(nBinsRaw),
                                    fPairThreshold);
      threshLine->SetLineColor(kBlack);
      threshLine->SetLineStyle(2);
      threshLine->SetLineWidth(2);
      threshLine->Draw();

      auto* legRaw = new TLegend(0.55, 0.82, 0.88, 0.88);
      legRaw->AddEntry(threshLine,
                       ("PairThreshold = " + std::to_string(fPairThreshold))
                           .c_str(),
                       "l");
      legRaw->Draw();

      cRaw->Write();
      if (fSavePlotPNGs) {
        const std::string png =
            fPlotOutputPrefix + "_debug_raw_pmt_" + tag + ".png";
        cRaw->SaveAs(png.c_str());
        std::cout << "MaxPairsAnalyzer: wrote " << png << " (" << colorIdx
                  << " PMT waveforms)" << std::endl;
      }
      fRawPMTDebugSaved = true;
    }
  }

  // Trigger efficiency: only events with maxPairs > 0 peak in EffPeakTs window.
  if (fHasInWindowFlash && fMaxPairs > 0 && fHEffDenomMaxPairs &&
      fHEffNumerMaxPairs) {
    fHEffDenomMaxPairs->Fill(fMaxPairs);
    if (fHasHLTBit) {
      fHEffNumerMaxPairs->Fill(fMaxPairs);
    }
  }
  if (fHasInWindowFlash && fHasMTCA && fHEffDenomNHitMTCA &&
      fHEffNumerNHitMTCA) {
    fHEffDenomNHitMTCA->Fill(fNHitMTCA);
    if (fHasHLTBit) {
      fHEffNumerNHitMTCA->Fill(fNHitMTCA);
    }
  }

  fTree->Fill();
}

void MaxPairsAnalyzer::endJob() {
  art::ServiceHandle<art::TFileService> tfs;

  if (fHMaxPairsVsNHitMTCA) {
    TCanvas* cScatter = tfs->make<TCanvas>(
        "cMaxPairsVsNHitMTCA", "Software maxPairs vs MTCA analog NHIT", 800, 600);
    cScatter->cd();
    fHMaxPairsVsNHitMTCA->Draw("COLZ");
    cScatter->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_maxpairs_vs_mtca.png";
      cScatter->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png << std::endl;
    }
  }

  if (fHMaxPairsMinusNHitMTCA) {
    TCanvas* cDiff = tfs->make<TCanvas>(
        "cMaxPairsMinusNHitMTCA", "maxPairs - nHIT_MTCA per event", 800, 600);
    cDiff->cd();
    fHMaxPairsMinusNHitMTCA->SetLineColor(kBlack);
    fHMaxPairsMinusNHitMTCA->SetLineWidth(2);
    fHMaxPairsMinusNHitMTCA->Draw("HIST");
    cDiff->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_maxpairs_minus_nhit_mtca.png";
      cDiff->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png << std::endl;
    }
  }

  if (fHBaselineSubMTCAADC) {
    const std::string ctitle =
        "Baseline-sub MTCA ADC (digital mon-pulse > " +
        std::to_string(fDigMonPulseNHitCut) + ")";
    TCanvas* cADC = tfs->make<TCanvas>(
        "cBaselineSubtractedMTCA_ADC_digCut", ctitle.c_str(), 800, 600);
    cADC->cd();
    cADC->SetLogy();
    fHBaselineSubMTCAADC->GetXaxis()->SetRangeUser(fADCPlotXMin, fADCPlotXMax);
    fHBaselineSubMTCAADC->Draw("HIST");
    cADC->Write();
    if (fSavePlotPNGs) {
      const std::string png =
          fPlotOutputPrefix + "_mtca_baseline_sub_adc_digcut.png";
      cADC->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png
                << " (entries=" << fHBaselineSubMTCAADC->GetEntries() << ")"
                << std::endl;
    }
  }

  if (fHHLTBits) {
    TCanvas* cHLT = tfs->make<TCanvas>("cHLTBits", "HLT bits fired", 800, 600);
    cHLT->cd();
    cHLT->SetLogy();
    fHHLTBits->SetMinimum(0.5);
    fHHLTBits->Draw("HIST");
    cHLT->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_hlt_bits.png";
      cHLT->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png
                << " (entries=" << fHHLTBits->GetEntries() << ")" << std::endl;
    }
  }

  if (fHHLT2SinceHLT1) {
    TCanvas* cDt = tfs->make<TCanvas>(
        "cHLT2SinceHLT1", "HLT2 time since HLT1", 800, 600);
    cDt->cd();
    fHHLT2SinceHLT1->Draw("HIST");
    cDt->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_hlt2_since_hlt1.png";
      cDt->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png
                << " (entries=" << fHHLT2SinceHLT1->GetEntries() << ")"
                << std::endl;
    }
  }

  // Keep sample-index hist in the ROOT file (no PNG).
  if (fHMaxPairsSample) {
    TCanvas* cSamp = tfs->make<TCanvas>(
        "cMaxPairsSample", "Sample of max digital mon-pulse", 800, 600);
    cSamp->cd();
    fHMaxPairsSample->Draw("HIST");
    cSamp->Write();
  }

  auto drawFlashTsCanvas = [&](TH1D* h, const std::string& cname,
                               const std::string& pngSuffix) {
    if (!h) {
      return;
    }
    TCanvas* c = tfs->make<TCanvas>(cname.c_str(), h->GetTitle(), 800, 600);
    c->cd();
    h->Draw("HIST");
    const double yLine = std::max(1.0, h->GetMaximum());
    TLine* lineZB = new TLine(0., 0., 0., yLine);
    lineZB->SetLineColor(kBlack);
    lineZB->SetLineWidth(2);
    lineZB->Draw();
    TLine* lineEffLo =
        new TLine(fEffPeakTsMin, 0., fEffPeakTsMin, yLine);
    lineEffLo->SetLineColor(kOrange + 7);
    lineEffLo->SetLineStyle(2);
    lineEffLo->SetLineWidth(2);
    lineEffLo->Draw();
    TLine* lineEffHi =
        new TLine(fEffPeakTsMax, 0., fEffPeakTsMax, yLine);
    lineEffHi->SetLineColor(kOrange + 7);
    lineEffHi->SetLineStyle(2);
    lineEffHi->SetLineWidth(2);
    lineEffHi->Draw();
    c->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + pngSuffix;
      c->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png
                << " (entries=" << h->GetEntries() << ")" << std::endl;
    }
  };
  // All local-max peak times across the job.
  drawFlashTsCanvas(fHPeakTsRelZB, "cPeakTsRelZB", "_peak_times.png");

  if (fHNPeaks) {
    TCanvas* cN = tfs->make<TCanvas>(
        "cNPeaks", "Max local-max peak height per event (any time)", 800, 600);
    cN->cd();
    fHNPeaks->Draw("HIST");
    cN->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_npeaks.png";
      cN->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png << std::endl;
    }
  }
  if (fHNPeaksWindow) {
    TCanvas* cNW = tfs->make<TCanvas>(
        "cNPeaksWindow",
        "Brightest maxPairs in efficiency peak window", 800, 600);
    cNW->cd();
    fHNPeaksWindow->Draw("HIST");
    cNW->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_npeaks_window.png";
      cNW->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png << std::endl;
    }
  }

  auto makeEffPercent = [&](TH1D* numer, TH1D* denom, const std::string& name,
                            const std::string& title) -> TH1D* {
    if (!numer || !denom) {
      return nullptr;
    }
    TH1D* heff = tfs->make<TH1D>(name.c_str(), title.c_str(), denom->GetNbinsX(),
                                 denom->GetXaxis()->GetXmin(),
                                 denom->GetXaxis()->GetXmax());
    heff->Divide(numer, denom, 1.0, 1.0, "B");
    heff->Scale(100.);
    heff->SetMinimum(0.);
    heff->SetMaximum(105.);
    heff->SetLineColor(kBlue + 1);
    heff->SetMarkerColor(kBlue + 1);
    heff->SetMarkerStyle(20);
    heff->SetMarkerSize(0.8);
    const std::string ytitle =
        "HLT bit " + std::to_string(fHLTEffBit) + " efficiency (%)";
    heff->GetYaxis()->SetTitle(ytitle.c_str());
    return heff;
  };

  if (TH1D* heffMP = makeEffPercent(
          fHEffNumerMaxPairs, fHEffDenomMaxPairs, "hEffHLT_vs_maxPairs",
          ("HLT bit " + std::to_string(fHLTEffBit) +
           " (BNB light+beam) efficiency vs maxPairs;maxPairs;Efficiency (%)")
              .c_str())) {
    TCanvas* cEffMP = tfs->make<TCanvas>("cEffHLT_vs_maxPairs",
                                         "HLT efficiency vs maxPairs", 800, 600);
    cEffMP->cd();
    heffMP->Draw("E1");
    cEffMP->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_eff_hlt_vs_maxpairs.png";
      cEffMP->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png << std::endl;
    }
  }

  if (TH1D* heffNH = makeEffPercent(
          fHEffNumerNHitMTCA, fHEffDenomNHitMTCA, "hEffHLT_vs_nHIT_MTCA",
          ("HLT bit " + std::to_string(fHLTEffBit) +
           " (BNB light+beam) efficiency vs nHIT_MTCA;nHIT_{MTCA};Efficiency (%)")
              .c_str())) {
    TCanvas* cEffNH = tfs->make<TCanvas>("cEffHLT_vs_nHIT_MTCA",
                                         "HLT efficiency vs nHIT_MTCA", 800, 600);
    cEffNH->cd();
    heffNH->Draw("E1");
    cEffNH->Write();
    if (fSavePlotPNGs) {
      const std::string png = fPlotOutputPrefix + "_eff_hlt_vs_nhit_mtca.png";
      cEffNH->SaveAs(png.c_str());
      std::cout << "MaxPairsAnalyzer: wrote " << png << std::endl;
    }
  }
}

DEFINE_ART_MODULE(MaxPairsAnalyzer)
