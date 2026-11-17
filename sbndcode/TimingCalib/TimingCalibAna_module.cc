////////////////////////////////////////////////////////////////////
// Class:       TimingCalibAna
// Plugin Type: analyzer
// File:        TimingCalibAna_module.cc
//
// Fills a ROOT TTree once per slice with slice-level OpFlash timing and
// event-level frame-shift information for use in timing calibration studies;
// in MC (detected via art::Event::isRealData) also the slice truth match,
// with the dk2nu v01_11_00 vertex-time fix of sbncode PR #675 (CAFMaker
// CorrectMCTiming) applied to the true and matched-slice flash times.
// A second tree (subrunTree) holds the exposure per subrun as gen2 CAFMaker
// fills hdr.pot / hdr.nbnbinfo / hdr.noffbeambnb (POT and gates for the
// cafpybara off-beam scaling), plus the number of events processed.
////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "dk2nu/tree/dk2nu.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"
#include "sbnobj/Common/Reco/CorrectedOpFlashTiming.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"
#include "sbnobj/SBND/Timing/TimingInfo.hh"

#include "TTree.h"

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace calib {
  class TimingCalibAna;
}

class calib::TimingCalibAna : public art::EDAnalyzer {
public:
  struct Config {
    fhicl::Atom<art::InputTag> SliceLabel{
      fhicl::Name("SliceLabel"),
      fhicl::Comment("Producer label for the recob::Slice collection"),
      art::InputTag("pandora")};

    fhicl::Atom<art::InputTag> PFParticleLabel{
      fhicl::Name("PFParticleLabel"),
      fhicl::Comment("Producer label for the recob::PFParticle collection, "
                      "and for the Slice->PFParticle, "
                      "PFParticle->PFParticleMetadata and "
                      "PFParticle->Vertex associations"),
      art::InputTag("pandora")};

    fhicl::Atom<art::InputTag> CorrectedOpFlashLabel{
      fhicl::Name("CorrectedOpFlashLabel"),
      fhicl::Comment("Producer label for the Slice->CorrectedOpFlashTiming "
                      "association (e.g. sbndcode's LightPropagationCorrection "
                      "module, instance \"lightpropagationcorrectionSCE\")"),
      art::InputTag("lightpropagationcorrectionSCE")};

    fhicl::Atom<art::InputTag> BarycenterFMLabel{
      fhicl::Name("BarycenterFMLabel"),
      fhicl::Comment("Producer label for the sbn::TPCPMTBarycenterMatch "
                      "collection and its Slice->TPCPMTBarycenterMatch and "
                      "TPCPMTBarycenterMatch->OpFlash associations (sbndcode's "
                      "TPCPMTBarycenterMatching module; must be the one "
                      "LightPropagationCorrection used)"),
      art::InputTag("tpcpmtbarycentermatchingSCE")};

    fhicl::Sequence<art::InputTag> OpFlashLabels{
      fhicl::Name("OpFlashLabels"),
      fhicl::Comment("Producer labels for the recob::OpFlash collections "
                      "(one per TPC), used for the event-level flash count "
                      "and the nearest other flash to the matched one"),
      std::vector<art::InputTag>{"opflashtpc0", "opflashtpc1"}};

    fhicl::Atom<art::InputTag> GeneratorLabel{
      fhicl::Name("GeneratorLabel"),
      fhicl::Comment("MC only: producer label for the neutrino simb::MCTruth "
                      "collection (beam interactions, incl. dirt) and the subrun "
                      "sumdata::POTSummary (subrunTree pot, as CAFMaker's GenLabel)"),
      art::InputTag("generator")};

    fhicl::Atom<art::InputTag> SimChannelLabel{
      fhicl::Name("SimChannelLabel"),
      fhicl::Comment("MC only: label of the sim::SimChannel collection, used for "
                      "the slice truth-match completeness denominator (as "
                      "CAFMaker's SimChannelLabel: true visible energy per neutrino)"),
      art::InputTag("simtpc2d:simpleSC")};

    fhicl::Atom<art::InputTag> SPECTDCLabel{
      fhicl::Name("SPECTDCLabel"),
      fhicl::Comment("Producer label for the sbnd::timing::DAQTimestamp "
                      "(SPEC-TDC) collection, used to require a channel-2 "
                      "(RWM) timestamp before filling the tree"),
      art::InputTag("tdcdecoder")};

    fhicl::Atom<art::InputTag> FrameShiftLabel{
      fhicl::Name("FrameShiftLabel"),
      fhicl::Comment("Producer label for the sbnd::timing::FrameShiftInfo "
                      "product (one per event, from sbndcode's FrameShift "
                      "module)"),
      art::InputTag("frameshift")};

    fhicl::Atom<art::InputTag> TimingInfoLabel{
      fhicl::Name("TimingInfoLabel"),
      fhicl::Comment("Producer label for the sbnd::timing::TimingInfo "
                      "product (one per event, from sbndcode's FrameShift "
                      "module)"),
      art::InputTag("frameshift")};

    fhicl::Atom<double> NuScoreThreshold{
      fhicl::Name("NuScoreThreshold"),
      fhicl::Comment(
        "Minimum Pandora NuScore, read from the PFParticleMetadata of the "
        "slice's primary neutrino-hypothesis PFParticle (|PdgCode| == 12 "
        "or 14), required to keep a slice. Mirrors sbndcode's "
        "LightPropagationCorrection module's slice selection (must match "
        "its NuScoreThreshold fhicl value for this dataset), so that this "
        "analyzer only considers slices for which CorrectedOpFlashTiming "
        "was actually produced"),
      0.1};

    fhicl::Atom<double> FMScoreThreshold{
      fhicl::Name("FMScoreThreshold"),
      fhicl::Comment(
        "Minimum BFM score LightPropagationCorrection requires before "
        "building a CorrectedOpFlashTiming (must match its FMScoreThreshold "
        "fhicl value for this dataset); only used for corr_flash_status"),
      0.001};

    fhicl::Atom<art::InputTag> BNBPOTDataLabel{
      fhicl::Name("BNBPOTDataLabel"),
      fhicl::Comment("Data only: producer label for the subrun "
        "std::vector<sbn::BNBSpillInfo> (on-beam POT + gates), as CAFMaker's "
        "BNBPOTDataLabel"),
      art::InputTag("bnbinfo")};

    fhicl::Atom<art::InputTag> OffbeamBNBCountDataLabel{
      fhicl::Name("OffbeamBNBCountDataLabel"),
      fhicl::Comment("Data only: producer label for the subrun "
        "std::vector<sbn::EXTCountInfo> (off-beam BNB gates), as CAFMaker's "
        "OffbeamBNBCountDataLabel"),
      art::InputTag("bnbextinfo")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit TimingCalibAna(Parameters const& config);

  TimingCalibAna(TimingCalibAna const&) = delete;
  TimingCalibAna(TimingCalibAna&&) = delete;
  TimingCalibAna& operator=(TimingCalibAna const&) = delete;
  TimingCalibAna& operator=(TimingCalibAna&&) = delete;

  void analyze(art::Event const& e) override;
  void beginSubRun(art::SubRun const& sr) override;
  void endSubRun(art::SubRun const& sr) override;

private:
  art::InputTag fSliceLabel;
  art::InputTag fPFParticleLabel;
  art::InputTag fCorrectedOpFlashLabel;
  art::InputTag fBarycenterFMLabel;
  std::vector<art::InputTag> fOpFlashLabels;
  art::InputTag fGeneratorLabel;
  art::InputTag fSimChannelLabel;
  art::InputTag fSPECTDCLabel;
  art::InputTag fFrameShiftLabel;
  art::InputTag fTimingInfoLabel;
  double fNuScoreThreshold;
  double fFMScoreThreshold;
  art::InputTag fBNBPOTDataLabel;
  art::InputTag fOffbeamBNBCountDataLabel;

  TTree* fTree{nullptr};
  TTree* fSubRunTree{nullptr};

  // subrunTree: exposure per subrun, filled as gen2 CAFMaker's beginSubRun
  // (first found of: BNBSpillInfo -> on-beam, EXTCountInfo -> off-beam,
  // MC POTSummary under GeneratorLabel). pot = sum of BNBSpillInfo::POT()
  // (TOR875, units 1e12) in on-beam data = hdr.pot, POTSummary::totgoodpot
  // in MC; n_bnb_spills = number of BNBSpillInfo = hdr.nbnbinfo (on-beam
  // gates); n_offbeam_gates = sum of EXTCountInfo::gates_since_last_trigger
  // = hdr.noffbeambnb; exposure_type 0 = none found, 1 = on-beam BNB,
  // 2 = off-beam BNB, 3 = MC POTSummary; n_events = events this module
  // processed in the subrun (all events, not only those with slices)
  unsigned int fSRRun{0};
  unsigned int fSRSubRun{0};
  double fSRPOT{0.};
  int fSRNBNBSpills{0};
  double fSRNOffbeamGates{0.};
  int fSRExposureType{0};
  int fSRNEvents{0};

  // Event identifiers
  unsigned int fRun{0};
  unsigned int fSubRun{0};
  unsigned int fEvent{0};

  // Slice identifier: index into the recob::Slice collection for this entry
  unsigned int fSliceId{0};

  // Slice-level reco: NuScore and vertex of the primary neutrino-hypothesis
  // PFParticle (the one with the highest NuScore), vertex in cm, -9999 if
  // the PFParticle has no vertex
  double fNuScore{-9999.};
  double fNuVtxX{-9999.};
  double fNuVtxY{-9999.};
  double fNuVtxZ{-9999.};

  // One tree entry per slice with a neutrino-hypothesis PFParticle
  // (clear-cosmic slices are not filled), with or without a
  // CorrectedOpFlashTiming. corr_flash_status mirrors
  // LightPropagationCorrection's skips: 0 = corrected, 1 = NuScore below
  // NuScoreThreshold, 2 = data event without a SPEC-TDC RWM (channel 2)
  // timestamp, 3 = no BFM match, 4 = best BFM score below FMScoreThreshold,
  // 5 = other (should not happen)
  int fCorrFlashStatus{-1};
  int fHasCorrFlash{0};

  // sbn::CorrectedOpFlashTiming, converted us -> ns (-9999 if none)
  double fOpFlashT0Ns{-9999.};
  double fNuToFLightNs{-9999.};
  double fNuToFChargeNs{-9999.};
  double fOpFlashT0CorrectedNs{-9999.};
  double fParticlePropTimeNs{-9999.};
  double fPhotonPropTimeNs{-9999.};

  // slice flash = the CorrectedOpFlashTiming OpFlash if any, else the flash
  // of the slice's highest-score BFM match (the one LightPropagationCorrection
  // would have used). raw_flash_t0_ns = its OpFlash::Time shifted as gen2
  // CAFMaker shifts opflash.time (SBNDShiftPMTReference: + frame_apply_at_caf_ns
  // in data, no shift in MC), -9999 if no flash. n_slices_same_flash =
  // number of other slices (incl. clear cosmic) whose highest-score BFM flash
  // is this slice flash
  double fRawFlashT0Ns{-9999.};
  int fNSlicesSameFlash{0};

  // sbn::TPCPMTBarycenterMatch of this slice pointing at the slice flash
  // (above): matched-flash block, times converted us -> ns, positions in cm,
  // -9999 if no such match
  double fBFMFlashTimeNs{-9999.};
  double fBFMFlashPEs{-9999.};
  double fBFMFlashCenterY{-9999.};
  double fBFMFlashCenterZ{-9999.};
  double fBFMFlashWidthY{-9999.};
  double fBFMFlashWidthZ{-9999.};

  // same BFM match: charge block and match quality (cm, deg, us -> ns), plus
  // the number of BFM matches of the slice and the best score among the
  // other ones (-9999 if none)
  int fBFMNMatch{0};
  double fBFMScore{-9999.};
  double fBFMSecondScore{-9999.};
  double fBFMChi2{-9999.};
  double fBFMDeltaY{-9999.};
  double fBFMDeltaZ{-9999.};
  double fBFMRadius{-9999.};
  double fBFMAngle{-9999.};
  double fBFMDeltaTNs{-9999.};
  double fBFMRadiusTrigger{-9999.};
  double fBFMChargeTotal{-9999.};
  double fBFMChargeCenterX{-9999.};
  double fBFMChargeCenterY{-9999.};
  double fBFMChargeCenterZ{-9999.};
  double fBFMChargeWidthX{-9999.};
  double fBFMChargeWidthY{-9999.};
  double fBFMChargeWidthZ{-9999.};

  // event-level pile-up: slices passing the NuScore cut, OpFlashes in all
  // OpFlashLabels collections, index (into OpFlashLabels) of the collection
  // holding the slice flash (-1 if not found), and the nearest other flash
  // (OpFlash::Time, signed other - slice flash, us -> ns, and its total PE) in
  // the same collection and in the other collections (-9999 if none)
  int fNNuSlices{0};
  int fNFlashes{0};
  int fFlashTPC{-1};
  double fNearFlashSameDtNs{-9999.};
  double fNearFlashSamePE{-9999.};
  double fNearFlashOtherDtNs{-9999.};
  double fNearFlashOtherPE{-9999.};

  // MC truth, filled only when !e.isRealData() (is_mc = 0 -> all -9999 / -1).
  // Slice <-> neutrino match mirrors sbncode's CAFMaker MatchSlice2Truth:
  // hit energy backtracked to GeneratorLabel MCTruths, index = the MCTruth
  // with the largest share if > 0.5 (else -1); pur = that share (filled even
  // when unmatched), eff = its slice energy / its true visible energy (IDE
  // energy of its G4 particles in all SimChannels, all planes: CAFMaker's
  // sum of SRTrueInteraction plane[c][p].visE),
  // cosmic_frac = slice energy share not from any GeneratorLabel MCTruth;
  // tmatch_vise_in_slc_gev / tmatch_vise_cosmic_gev = CAFMaker SRTruthMatch
  // visEinslc / visEcosmic: backtracked hit energy of the slice (all planes,
  // CAFRecoUtils::TotalHitEnergy) and its part not from any GeneratorLabel
  // MCTruth, GeV; filled as CAFMaker: matched slices, or events without
  // neutrino MCTruths (both = the slice total), else -9999.
  // True neutrino: GeV, cm, time = G4 ns + mc_time_offset_ns;
  // near_true_nu_dt_ns = signed time to the nearest other GeneratorLabel
  // neutrino in the event (both offset-corrected).
  // dk2nu v01_11_00 SetT(0) fix, as sbncode PR #675 CorrectMCTiming: per GENIE
  // neutrino (MCTruth with a GTruth), offset = prod_time (meson ancestor
  // Dk2Nu::ancestor.back().startt, ns) + dk2gen (MCFlux::fdk2gen, decay -> flux
  // window, m -> cm) / c; added to true_nu_time_ns and to the matched slice's
  // opflash_t0_ns, opflash_t0_corrected_ns, raw_flash_t0_ns and bfm_flash_time_ns.
  // mc_time_offset_ns = offset applied to this slice (0 if none);
  // true_nu_g4_time_ns = generator vertex time before the offset (MCNeutrino T;
  // = bunch time + z_vtx / c, the flight from detector z = 0 included);
  // true_nu_gen2vtx_cm = MCFlux::fgen2vtx (flux window -> vertex), saved only:
  // gen2vtx / c - z_vtx / c = the flux window -> z = 0 leg, not in the offset;
  // true_nu_prod_time_ns / true_nu_dk2gen_cm / true_nu_gen2vtx_cm = -9999 if not available
  int fIsMC{0};
  int fNTrueNu{-1};
  int fTMatchIndex{-1};
  double fTMatchPur{-9999.};
  double fTMatchEff{-9999.};
  double fTMatchCosmicFrac{-9999.};
  double fTMatchVisEInSlcGeV{-9999.};
  double fTMatchVisECosmicGeV{-9999.};
  int fTrueNuPdg{0};
  int fTrueNuCCNC{-1};
  int fTrueNuMode{-1};
  double fTrueNuE{-9999.};
  double fTrueNuVtxX{-9999.};
  double fTrueNuVtxY{-9999.};
  double fTrueNuVtxZ{-9999.};
  double fTrueNuTimeNs{-9999.};
  double fNearTrueNuDtNs{-9999.};
  double fTrueNuProdTimeNs{-9999.};
  double fTrueNuDk2GenCm{-9999.};
  double fTrueNuGen2VtxCm{-9999.};
  double fTrueNuG4TimeNs{-9999.};
  double fMCTimeOffsetNs{0.};

  // sbnd::timing::FrameShiftInfo: one per event, already stored in ns
  // (data only, 0 in MC)
  double fFrameTdcRwmNs{0.};
  double fFrameHltBeamGateNs{0.};

  // sbnd::timing::TimingInfo: one per event, already stored in ns
  // (data only, 0 in MC)
  uint64_t fTdcRwmNs{0};
  uint64_t fTdcBesNs{0};
  uint64_t fHltBeamGateNs{0};

  // DAQTimestamp::TimestampPs() of the SPEC-TDC RWM (channel 2) / BES
  // (channel 1) entry whose Timestamp() equals TimingInfo::TdcRwm() /
  // TdcBes(): ps within the UTC second (coarse * 8000 + frac * 8000 / 4096);
  // TimingInfo has no ps field (data only, 0 in MC or when no entry matches)
  uint64_t fTdcRwmPs{0};
  uint64_t fTdcBesPs{0};

  // FrameShiftInfo::FrameApplyAtCaf, ns: shift gen2 CAFMaker applies to PMT
  // times (opflash) in data; 0 in MC, or when NaN (CAFMaker then skips
  // the shift)
  double fFrameApplyAtCafNs{0.};
};

calib::TimingCalibAna::TimingCalibAna(Parameters const& config)
  : EDAnalyzer(config)
  , fSliceLabel(config().SliceLabel())
  , fPFParticleLabel(config().PFParticleLabel())
  , fCorrectedOpFlashLabel(config().CorrectedOpFlashLabel())
  , fBarycenterFMLabel(config().BarycenterFMLabel())
  , fOpFlashLabels(config().OpFlashLabels())
  , fGeneratorLabel(config().GeneratorLabel())
  , fSimChannelLabel(config().SimChannelLabel())
  , fSPECTDCLabel(config().SPECTDCLabel())
  , fFrameShiftLabel(config().FrameShiftLabel())
  , fTimingInfoLabel(config().TimingInfoLabel())
  , fNuScoreThreshold(config().NuScoreThreshold())
  , fFMScoreThreshold(config().FMScoreThreshold())
  , fBNBPOTDataLabel(config().BNBPOTDataLabel())
  , fOffbeamBNBCountDataLabel(config().OffbeamBNBCountDataLabel())
{
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("timingCalibTree", "Timing calibration tree");
  fTree->Branch("run", &fRun);
  fTree->Branch("subrun", &fSubRun);
  fTree->Branch("event", &fEvent);
  fTree->Branch("slice_id", &fSliceId);
  fTree->Branch("nu_score", &fNuScore);
  fTree->Branch("nu_vtx_x", &fNuVtxX);
  fTree->Branch("nu_vtx_y", &fNuVtxY);
  fTree->Branch("nu_vtx_z", &fNuVtxZ);
  fTree->Branch("corr_flash_status", &fCorrFlashStatus);
  fTree->Branch("has_corr_flash", &fHasCorrFlash);
  fTree->Branch("opflash_t0_ns", &fOpFlashT0Ns);
  fTree->Branch("nu_tof_light_ns", &fNuToFLightNs);
  fTree->Branch("nu_tof_charge_ns", &fNuToFChargeNs);
  fTree->Branch("opflash_t0_corrected_ns", &fOpFlashT0CorrectedNs);
  fTree->Branch("particle_prop_time_ns", &fParticlePropTimeNs);
  fTree->Branch("photon_prop_time_ns", &fPhotonPropTimeNs);
  fTree->Branch("raw_flash_t0_ns", &fRawFlashT0Ns);
  fTree->Branch("n_slices_same_flash", &fNSlicesSameFlash);
  fTree->Branch("bfm_flash_time_ns", &fBFMFlashTimeNs);
  fTree->Branch("bfm_flash_pes", &fBFMFlashPEs);
  fTree->Branch("bfm_flash_center_y", &fBFMFlashCenterY);
  fTree->Branch("bfm_flash_center_z", &fBFMFlashCenterZ);
  fTree->Branch("bfm_flash_width_y", &fBFMFlashWidthY);
  fTree->Branch("bfm_flash_width_z", &fBFMFlashWidthZ);
  fTree->Branch("bfm_n_match", &fBFMNMatch);
  fTree->Branch("bfm_score", &fBFMScore);
  fTree->Branch("bfm_second_score", &fBFMSecondScore);
  fTree->Branch("bfm_chi2", &fBFMChi2);
  fTree->Branch("bfm_delta_y", &fBFMDeltaY);
  fTree->Branch("bfm_delta_z", &fBFMDeltaZ);
  fTree->Branch("bfm_radius", &fBFMRadius);
  fTree->Branch("bfm_angle", &fBFMAngle);
  fTree->Branch("bfm_delta_t_ns", &fBFMDeltaTNs);
  fTree->Branch("bfm_radius_trigger", &fBFMRadiusTrigger);
  fTree->Branch("bfm_charge_total", &fBFMChargeTotal);
  fTree->Branch("bfm_charge_center_x", &fBFMChargeCenterX);
  fTree->Branch("bfm_charge_center_y", &fBFMChargeCenterY);
  fTree->Branch("bfm_charge_center_z", &fBFMChargeCenterZ);
  fTree->Branch("bfm_charge_width_x", &fBFMChargeWidthX);
  fTree->Branch("bfm_charge_width_y", &fBFMChargeWidthY);
  fTree->Branch("bfm_charge_width_z", &fBFMChargeWidthZ);
  fTree->Branch("n_nu_slices", &fNNuSlices);
  fTree->Branch("n_flashes", &fNFlashes);
  fTree->Branch("flash_tpc", &fFlashTPC);
  fTree->Branch("near_flash_same_dt_ns", &fNearFlashSameDtNs);
  fTree->Branch("near_flash_same_pe", &fNearFlashSamePE);
  fTree->Branch("near_flash_other_dt_ns", &fNearFlashOtherDtNs);
  fTree->Branch("near_flash_other_pe", &fNearFlashOtherPE);
  fTree->Branch("is_mc", &fIsMC);
  fTree->Branch("n_true_nu", &fNTrueNu);
  fTree->Branch("tmatch_index", &fTMatchIndex);
  fTree->Branch("tmatch_pur", &fTMatchPur);
  fTree->Branch("tmatch_eff", &fTMatchEff);
  fTree->Branch("tmatch_cosmic_frac", &fTMatchCosmicFrac);
  fTree->Branch("tmatch_vise_in_slc_gev", &fTMatchVisEInSlcGeV);
  fTree->Branch("tmatch_vise_cosmic_gev", &fTMatchVisECosmicGeV);
  fTree->Branch("true_nu_pdg", &fTrueNuPdg);
  fTree->Branch("true_nu_ccnc", &fTrueNuCCNC);
  fTree->Branch("true_nu_mode", &fTrueNuMode);
  fTree->Branch("true_nu_e", &fTrueNuE);
  fTree->Branch("true_nu_vtx_x", &fTrueNuVtxX);
  fTree->Branch("true_nu_vtx_y", &fTrueNuVtxY);
  fTree->Branch("true_nu_vtx_z", &fTrueNuVtxZ);
  fTree->Branch("true_nu_time_ns", &fTrueNuTimeNs);
  fTree->Branch("near_true_nu_dt_ns", &fNearTrueNuDtNs);
  fTree->Branch("true_nu_prod_time_ns", &fTrueNuProdTimeNs);
  fTree->Branch("true_nu_dk2gen_cm", &fTrueNuDk2GenCm);
  fTree->Branch("true_nu_gen2vtx_cm", &fTrueNuGen2VtxCm);
  fTree->Branch("true_nu_g4_time_ns", &fTrueNuG4TimeNs);
  fTree->Branch("mc_time_offset_ns", &fMCTimeOffsetNs);
  fTree->Branch("frame_tdc_rwm_ns", &fFrameTdcRwmNs);
  fTree->Branch("frame_hlt_beam_gate_ns", &fFrameHltBeamGateNs);
  fTree->Branch("tdc_rwm_ns", &fTdcRwmNs);
  fTree->Branch("tdc_rwm_ps", &fTdcRwmPs);
  fTree->Branch("tdc_bes_ns", &fTdcBesNs);
  fTree->Branch("tdc_bes_ps", &fTdcBesPs);
  fTree->Branch("hlt_beam_gate_ns", &fHltBeamGateNs);
  fTree->Branch("frame_apply_at_caf_ns", &fFrameApplyAtCafNs);

  fSubRunTree = tfs->make<TTree>("subrunTree", "Exposure per subrun (POT, gates)");
  fSubRunTree->Branch("run", &fSRRun);
  fSubRunTree->Branch("subrun", &fSRSubRun);
  fSubRunTree->Branch("pot", &fSRPOT);
  fSubRunTree->Branch("n_bnb_spills", &fSRNBNBSpills);
  fSubRunTree->Branch("n_offbeam_gates", &fSRNOffbeamGates);
  fSubRunTree->Branch("exposure_type", &fSRExposureType);
  fSubRunTree->Branch("n_events", &fSRNEvents);
}

void
calib::TimingCalibAna::beginSubRun(art::SubRun const& sr)
{
  fSRRun = sr.run();
  fSRSubRun = sr.subRun();
  fSRPOT = 0.;
  fSRNBNBSpills = 0;
  fSRNOffbeamGates = 0.;
  fSRExposureType = 0;
  fSRNEvents = 0;

  auto const bnb_spill = sr.getHandle<std::vector<sbn::BNBSpillInfo>>(fBNBPOTDataLabel);
  auto const bnb_offbeam_spill = sr.getHandle<std::vector<sbn::EXTCountInfo>>(fOffbeamBNBCountDataLabel);

  // as CAFMaker: at most one exposure source per subrun
  if (bnb_spill && bnb_offbeam_spill)
    throw cet::exception("TimingCalibAna")
      << "Found both " << fBNBPOTDataLabel << " and " << fOffbeamBNBCountDataLabel
      << " in run " << fSRRun << " subrun " << fSRSubRun << "\n";

  if (bnb_spill) {
    fSRExposureType = 1;
    fSRNBNBSpills = bnb_spill->size();
    for (sbn::BNBSpillInfo const& info : *bnb_spill)
      fSRPOT += info.POT();
  }
  else if (bnb_offbeam_spill) {
    fSRExposureType = 2;
    for (sbn::EXTCountInfo const& spill : *bnb_offbeam_spill)
      fSRNOffbeamGates += spill.gates_since_last_trigger;
  }
  else if (auto const pot_handle = sr.getHandle<sumdata::POTSummary>(fGeneratorLabel)) {
    fSRExposureType = 3;
    fSRPOT = pot_handle->totgoodpot;
  }
  else {
    mf::LogWarning("TimingCalibAna")
      << "No exposure found in run " << fSRRun << " subrun " << fSRSubRun << " (" << fBNBPOTDataLabel
      << ", " << fOffbeamBNBCountDataLabel << ", " << fGeneratorLabel << " POTSummary)";
  }
}

void
calib::TimingCalibAna::endSubRun(art::SubRun const&)
{
  fSubRunTree->Fill();
}

void
calib::TimingCalibAna::analyze(art::Event const& e)
{
  fIsMC = e.isRealData() ? 0 : 1;
  ++fSRNEvents;

  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.event();

  // SPEC-TDC and FrameShift products exist in data only (as in
  // LightPropagationCorrection's IsMC mode); in MC the frame branches stay 0
  fFrameTdcRwmNs = fFrameHltBeamGateNs = 0.;
  fTdcRwmNs = fTdcBesNs = fHltBeamGateNs = fTdcRwmPs = fTdcBesPs = 0;
  fFrameApplyAtCafNs = 0.;
  bool hasRWM = fIsMC;
  if (!fIsMC) {
    // events without a SPEC-TDC RWM timestamp are kept, with
    // corr_flash_status 2 for every slice, since LightPropagationCorrection
    // can't build OpFlashT0 / OpFlashT0Corrected there
    auto const tdcHandle = e.getHandle<std::vector<sbnd::timing::DAQTimestamp>>(fSPECTDCLabel);
    if (tdcHandle.isValid()) {
      for (auto const& tdc : *tdcHandle) {
        if (tdc.Channel() == 2) hasRWM = true;
      }
    }

    auto const& frameShiftHandle =
      e.getValidHandle<sbnd::timing::FrameShiftInfo>(fFrameShiftLabel);
    fFrameTdcRwmNs = frameShiftHandle->FrameTdcRwm();
    fFrameHltBeamGateNs = frameShiftHandle->FrameHltBeamGate();
    if (!std::isnan(frameShiftHandle->FrameApplyAtCaf()))
      fFrameApplyAtCafNs = frameShiftHandle->FrameApplyAtCaf();

    auto const& timingInfoHandle =
      e.getValidHandle<sbnd::timing::TimingInfo>(fTimingInfoLabel);
    fTdcRwmNs = timingInfoHandle->TdcRwm();
    fTdcBesNs = timingInfoHandle->TdcBes();
    fHltBeamGateNs = timingInfoHandle->HltBeamGate();

    // FrameShift picks one RWM / BES by its ns Timestamp(); take the ps of
    // that same entry
    if (tdcHandle.isValid()) {
      for (auto const& tdc : *tdcHandle) {
        if (fTdcRwmNs != 0 && tdc.Channel() == 2 && tdc.Timestamp() == fTdcRwmNs)
          fTdcRwmPs = tdc.TimestampPs();
        if (fTdcBesNs != 0 && tdc.Channel() == 1 && tdc.Timestamp() == fTdcBesNs)
          fTdcBesPs = tdc.TimestampPs();
      }
    }
  }

  auto const& sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  auto const& pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

  art::FindManyP<recob::PFParticle> slicePFPAssns(sliceHandle, e, fPFParticleLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssns(
    pfpHandle, e, fPFParticleLabel);
  art::FindManyP<recob::Vertex> pfpVertexAssns(pfpHandle, e, fPFParticleLabel);
  art::FindManyP<sbn::CorrectedOpFlashTiming> sliceOpFlashAssns(
    sliceHandle, e, fCorrectedOpFlashLabel);

  auto const& corrOpFlashHandle =
    e.getValidHandle<std::vector<sbn::CorrectedOpFlashTiming>>(fCorrectedOpFlashLabel);
  art::FindManyP<recob::OpFlash> corrOpFlashFlashAssns(
    corrOpFlashHandle, e, fCorrectedOpFlashLabel);

  auto const& bfmHandle =
    e.getValidHandle<std::vector<sbn::TPCPMTBarycenterMatch>>(fBarycenterFMLabel);
  art::FindManyP<sbn::TPCPMTBarycenterMatch> sliceBFMAssns(
    sliceHandle, e, fBarycenterFMLabel);
  art::FindManyP<recob::OpFlash> bfmFlashAssns(bfmHandle, e, fBarycenterFMLabel);

  std::vector<art::Handle<std::vector<recob::OpFlash>>> flashHandles;
  fNFlashes = 0;
  for (auto const& label : fOpFlashLabels) {
    flashHandles.push_back(e.getHandle<std::vector<recob::OpFlash>>(label));
    if (flashHandles.back().isValid()) fNFlashes += flashHandles.back()->size();
  }

  constexpr double kUsToNs = 1000.;
  // us -> ns, keeping the producers' -9999 "not set" default as -9999
  auto usToNs = [](double us) { return us == -9999. ? -9999. : us * kUsToNs; };

  // MC truth setup: neutrino MCTruths, and true visible energy (SimChannel
  // IDEs) per MCTruth (completeness denominator)
  fNTrueNu = -1;
  std::vector<art::Ptr<simb::MCTruth>> nuTruths;
  std::optional<art::FindManyP<recob::Hit>> sliceHitAssns;
  std::map<art::Ptr<simb::MCTruth>, double> eventTruthEnergy;
  // per nuTruths entry: dk2nu prod_time (ns), dk2gen and gen2vtx (cm), and the
  // timing offset (ns, NaN if not a GENIE neutrino with dk2nu info)
  std::vector<double> nuProdTimeNs, nuDk2GenCm, nuGen2VtxCm, nuTimeOffsetNs;
  detinfo::DetectorClocksData const clockData =
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  // backtracked hit energy (MeV) per GeneratorLabel MCTruth; total = all
  // backtracked energy in the hits
  auto const truthEnergies = [&](std::vector<art::Ptr<recob::Hit>> const& hits, double& total) {
    art::ServiceHandle<cheat::BackTrackerService const> bt;
    art::ServiceHandle<cheat::ParticleInventoryService const> pi;
    std::map<art::Ptr<simb::MCTruth>, double> energies;
    total = 0.;
    for (auto const& hit : hits) {
      for (auto const& ide : bt->HitToTrackIDEs(clockData, hit)) {
        total += ide.energy;
        art::Ptr<simb::MCTruth> truth;
        try {
          truth = pi->TrackIdToMCTruth_P(std::abs(ide.trackID));
        }
        // ignore track IDs that cannot be looked up (as CAFMaker does)
        catch (...) {
          continue;
        }
        if (std::find(nuTruths.begin(), nuTruths.end(), truth) != nuTruths.end())
          energies[truth] += ide.energy;
      }
    }
    return energies;
  };
  if (fIsMC) {
    auto const truthHandle = e.getHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);
    if (truthHandle.isValid()) art::fill_ptr_vector(nuTruths, truthHandle);
    fNTrueNu = nuTruths.size();
    sliceHitAssns.emplace(sliceHandle, e, fSliceLabel);

    // dk2nu timing offset per neutrino (MCFlux / Dk2Nu indexed as the MCTruths)
    nuProdTimeNs.assign(nuTruths.size(), -9999.);
    nuDk2GenCm.assign(nuTruths.size(), -9999.);
    nuGen2VtxCm.assign(nuTruths.size(), -9999.);
    nuTimeOffsetNs.assign(nuTruths.size(), std::numeric_limits<double>::quiet_NaN());
    auto const fluxHandle = e.getHandle<std::vector<simb::MCFlux>>(fGeneratorLabel);
    auto const dk2nuHandle = e.getHandle<std::vector<bsim::Dk2Nu>>(fGeneratorLabel);
    std::optional<art::FindManyP<simb::GTruth>> truthGTruthAssns;
    if (e.getHandle<art::Assns<simb::MCTruth, simb::GTruth>>(fGeneratorLabel).isValid())
      truthGTruthAssns.emplace(truthHandle, e, fGeneratorLabel);
    for (unsigned int iNu = 0; iNu < nuTruths.size(); ++iNu) {
      if (dk2nuHandle.isValid() && iNu < dk2nuHandle->size() &&
          !dk2nuHandle->at(iNu).ancestor.empty())
        nuProdTimeNs[iNu] = dk2nuHandle->at(iNu).ancestor.back().startt;
      if (fluxHandle.isValid() && iNu < fluxHandle->size() && fluxHandle->at(iNu).fdk2gen > -998.)
        nuDk2GenCm[iNu] = fluxHandle->at(iNu).fdk2gen * util::kMeterToCentimeter;
      if (fluxHandle.isValid() && iNu < fluxHandle->size() && fluxHandle->at(iNu).fgen2vtx > -998.)
        nuGen2VtxCm[iNu] = fluxHandle->at(iNu).fgen2vtx * util::kMeterToCentimeter;
      bool const isGENIE = truthGTruthAssns && !truthGTruthAssns->at(iNu).empty();
      if (isGENIE && nuProdTimeNs[iNu] != -9999. && nuDk2GenCm[iNu] != -9999.)
        nuTimeOffsetNs[iNu] = nuProdTimeNs[iNu] + nuDk2GenCm[iNu] / util::kc;
    }
    // true visible energy per neutrino MCTruth, as CAFMaker (PrepSimChannels +
    // visE per particle, summed over the interaction's particles): IDE energy
    // (MeV) of every SimChannel, all planes, track ID |id| -> MCTruth
    if (!nuTruths.empty()) {
      art::ServiceHandle<cheat::ParticleInventoryService const> pi;
      auto const simChannelHandle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);
      std::map<int, art::Ptr<simb::MCTruth>> idToTruth;     // lookup cache; null = not found
      for (auto const& sc : *simChannelHandle) {
        for (auto const& item : sc.TDCIDEMap()) {
          for (auto const& ide : item.second) {
            int const id = std::abs(ide.trackID);
            auto it = idToTruth.find(id);
            if (it == idToTruth.end()) {
              art::Ptr<simb::MCTruth> truth;
              try {
                truth = pi->TrackIdToMCTruth_P(id);
              }
              // ignore track IDs that cannot be looked up (as CAFMaker does)
              catch (...) {}
              it = idToTruth.emplace(id, truth).first;
            }
            if (it->second.isNonnull() &&
                std::find(nuTruths.begin(), nuTruths.end(), it->second) != nuTruths.end())
              eventTruthEnergy[it->second] += ide.energy;
          }
        }
      }
    }
  }

  // Mirror sbndcode's LightPropagationCorrection_module.cc slice
  // selection: find the primary neutrino-hypothesis PFParticle
  // (|PdgCode| == 12 or 14) and require its NuScore to clear a
  // threshold. Clear-cosmic slices have a track-like (pdg 13) primary
  // instead, so they never contribute a NuScore and are not filled,
  // consistent with them never getting a CorrectedOpFlashTiming.
  // Done for all slices first so the event-level n_nu_slices is known, along
  // with each slice's highest-score BFM match (LightPropagationCorrection's
  // HighestBFMScoreIdx: first one wins ties) and its flash.
  std::vector<double> sliceNuScores(sliceHandle->size(), -9999.);
  std::vector<art::Ptr<recob::PFParticle>> sliceNuPFPs(sliceHandle->size());
  std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>> sliceBestBFMs(sliceHandle->size());
  std::vector<art::Ptr<recob::OpFlash>> sliceBestBFMFlashes(sliceHandle->size());
  fNNuSlices = 0;
  for (unsigned int iSlice = 0; iSlice < sliceHandle->size(); ++iSlice) {
    for (auto const& bfm : sliceBFMAssns.at(iSlice)) {
      if (sliceBestBFMs[iSlice].isNull() || bfm->score > sliceBestBFMs[iSlice]->score)
        sliceBestBFMs[iSlice] = bfm;
    }
    if (sliceBestBFMs[iSlice].isNonnull()) {
      auto const& bestFlash = bfmFlashAssns.at(sliceBestBFMs[iSlice].key());
      if (!bestFlash.empty()) sliceBestBFMFlashes[iSlice] = bestFlash.front();
    }

    double& sliceMaxNuScore = sliceNuScores[iSlice];
    art::Ptr<recob::PFParticle>& nuPFP = sliceNuPFPs[iSlice];
    for (auto const& pfp : slicePFPAssns.at(iSlice)) {
      if (!pfp->IsPrimary()) continue;
      if (std::abs(pfp->PdgCode()) != 12 && std::abs(pfp->PdgCode()) != 14) continue;
      for (auto const& meta : pfpMetadataAssns.at(pfp.key())) {
        auto const& props = meta->GetPropertiesMap();
        auto const it = props.find("NuScore");
        if (it != props.end() && it->second > sliceMaxNuScore) {
          sliceMaxNuScore = it->second;
          nuPFP = pfp;
        }
      }
    }
    if (sliceMaxNuScore >= fNuScoreThreshold) ++fNNuSlices;
  }

  for (unsigned int iSlice = 0; iSlice < sliceHandle->size(); ++iSlice) {
    double const sliceMaxNuScore = sliceNuScores[iSlice];
    art::Ptr<recob::PFParticle> const& nuPFP = sliceNuPFPs[iSlice];
    if (nuPFP.isNull()) continue; // clear cosmic

    fSliceId = iSlice;

    fNuScore = sliceMaxNuScore;
    fNuVtxX = fNuVtxY = fNuVtxZ = -9999.;
    std::vector<art::Ptr<recob::Vertex>> const& nuVertex = pfpVertexAssns.at(nuPFP.key());
    if (!nuVertex.empty()) {
      auto const& pos = nuVertex.front()->position();
      fNuVtxX = pos.X();
      fNuVtxY = pos.Y();
      fNuVtxZ = pos.Z();
    }

    // corrected OpFlash timing associated with this slice (0 or 1 per
    // slice), ignored if never filled by the producer (NaN sentinel)
    art::Ptr<sbn::CorrectedOpFlashTiming> corrPtr;
    std::vector<art::Ptr<sbn::CorrectedOpFlashTiming>> const& sliceOpFlash =
      sliceOpFlashAssns.at(iSlice);
    if (!sliceOpFlash.empty() && !std::isnan(sliceOpFlash.front()->OpFlashT0))
      corrPtr = sliceOpFlash.front();
    fHasCorrFlash = corrPtr.isNonnull() ? 1 : 0;

    fOpFlashT0Ns = fNuToFLightNs = fNuToFChargeNs = fOpFlashT0CorrectedNs = -9999.;
    fParticlePropTimeNs = fPhotonPropTimeNs = -9999.;
    art::Ptr<recob::OpFlash> corrFlash;
    if (corrPtr.isNonnull()) {
      fOpFlashT0Ns = corrPtr->OpFlashT0 * kUsToNs;
      fNuToFLightNs = corrPtr->NuToFLight * kUsToNs;
      fNuToFChargeNs = corrPtr->NuToFCharge * kUsToNs;
      fOpFlashT0CorrectedNs = corrPtr->OpFlashT0Corrected * kUsToNs;
      fParticlePropTimeNs = corrPtr->ParticlePropagationTime * kUsToNs;
      fPhotonPropTimeNs = corrPtr->PhotonPropagationTime * kUsToNs;
      std::vector<art::Ptr<recob::OpFlash>> const& corrFlashes =
        corrOpFlashFlashAssns.at(corrPtr.key());
      if (!corrFlashes.empty()) corrFlash = corrFlashes.front();
    }

    std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>> const& sliceBFMs =
      sliceBFMAssns.at(iSlice);
    art::Ptr<sbn::TPCPMTBarycenterMatch> const& bestBFM = sliceBestBFMs[iSlice];

    // same order as LightPropagationCorrection's skips
    if (sliceMaxNuScore < fNuScoreThreshold) fCorrFlashStatus = 1;
    else if (!hasRWM) fCorrFlashStatus = 2;
    else if (sliceBFMs.empty()) fCorrFlashStatus = 3;
    else if (bestBFM->score < fFMScoreThreshold) fCorrFlashStatus = 4;
    else if (corrPtr.isNull()) fCorrFlashStatus = 5;
    else fCorrFlashStatus = 0;

    // slice flash: the corrected one, else the highest-score BFM flash
    art::Ptr<recob::OpFlash> const sliceFlash =
      corrFlash.isNonnull() ? corrFlash : sliceBestBFMFlashes[iSlice];

    fRawFlashT0Ns = -9999.;
    if (sliceFlash.isNonnull()) fRawFlashT0Ns = sliceFlash->Time() * kUsToNs + fFrameApplyAtCafNs;

    fNSlicesSameFlash = 0;
    if (sliceFlash.isNonnull()) {
      for (unsigned int jSlice = 0; jSlice < sliceHandle->size(); ++jSlice) {
        if (jSlice != iSlice && sliceBestBFMFlashes[jSlice] == sliceFlash) ++fNSlicesSameFlash;
      }
    }

    // BFM match of this slice pointing at the slice flash
    // (LightPropagationCorrection's BarycenterFM mode builds the
    // CorrectedOpFlashTiming from the slice's highest-score BFM flash)
    fBFMFlashTimeNs = fBFMFlashPEs = -9999.;
    fBFMFlashCenterY = fBFMFlashCenterZ = fBFMFlashWidthY = fBFMFlashWidthZ = -9999.;
    fBFMScore = fBFMSecondScore = fBFMChi2 = -9999.;
    fBFMDeltaY = fBFMDeltaZ = fBFMRadius = fBFMAngle = fBFMDeltaTNs = -9999.;
    fBFMRadiusTrigger = fBFMChargeTotal = -9999.;
    fBFMChargeCenterX = fBFMChargeCenterY = fBFMChargeCenterZ = -9999.;
    fBFMChargeWidthX = fBFMChargeWidthY = fBFMChargeWidthZ = -9999.;
    fBFMNMatch = sliceBFMs.size();
    if (sliceFlash.isNonnull()) {
      for (auto const& bfm : sliceBFMs) {
        std::vector<art::Ptr<recob::OpFlash>> const& bfmFlash = bfmFlashAssns.at(bfm.key());
        if (bfmFlash.empty() || bfmFlash.front() != sliceFlash) {
          if (bfm->score > fBFMSecondScore) fBFMSecondScore = bfm->score;
          continue;
        }
        if (fBFMScore != -9999.) continue; // already filled
        fBFMFlashTimeNs = bfm->flashTime * kUsToNs;
        fBFMFlashPEs = bfm->flashPEs;
        fBFMFlashCenterY = bfm->flashCenter.Y();
        fBFMFlashCenterZ = bfm->flashCenter.Z();
        fBFMFlashWidthY = bfm->flashWidth.Y();
        fBFMFlashWidthZ = bfm->flashWidth.Z();
        fBFMScore = bfm->score;
        fBFMChi2 = bfm->chi2;
        fBFMDeltaY = bfm->deltaY;
        fBFMDeltaZ = bfm->deltaZ;
        fBFMRadius = bfm->radius;
        fBFMAngle = bfm->angle;
        fBFMDeltaTNs = usToNs(bfm->deltaT);
        fBFMRadiusTrigger = bfm->radius_Trigger;
        fBFMChargeTotal = bfm->chargeTotal;
        fBFMChargeCenterX = bfm->chargeCenter.X();
        fBFMChargeCenterY = bfm->chargeCenter.Y();
        fBFMChargeCenterZ = bfm->chargeCenter.Z();
        fBFMChargeWidthX = bfm->chargeWidth.X();
        fBFMChargeWidthY = bfm->chargeWidth.Y();
        fBFMChargeWidthZ = bfm->chargeWidth.Z();
      }
    }

    // nearest other OpFlash to the slice flash, in the same collection
    // (same TPC) and in the other collections
    fFlashTPC = -1;
    fNearFlashSameDtNs = fNearFlashSamePE = -9999.;
    fNearFlashOtherDtNs = fNearFlashOtherPE = -9999.;
    if (sliceFlash.isNonnull()) {
      art::Ptr<recob::OpFlash> const& matchedFlash = sliceFlash;
      for (unsigned int iLabel = 0; iLabel < flashHandles.size(); ++iLabel) {
        if (flashHandles[iLabel].isValid() && flashHandles[iLabel].id() == matchedFlash.id())
          fFlashTPC = iLabel;
      }
      double const matchedTime = matchedFlash->Time();
      for (unsigned int iLabel = 0; iLabel < flashHandles.size(); ++iLabel) {
        auto const& flashHandle = flashHandles[iLabel];
        if (!flashHandle.isValid()) continue;
        bool const same = (static_cast<int>(iLabel) == fFlashTPC);
        double& nearDt = same ? fNearFlashSameDtNs : fNearFlashOtherDtNs;
        double& nearPE = same ? fNearFlashSamePE : fNearFlashOtherPE;
        for (unsigned int iFlash = 0; iFlash < flashHandle->size(); ++iFlash) {
          if (same && iFlash == matchedFlash.key()) continue;
          recob::OpFlash const& flash = flashHandle->at(iFlash);
          double const dt = (flash.Time() - matchedTime) * kUsToNs;
          if (nearDt == -9999. || std::abs(dt) < std::abs(nearDt)) {
            nearDt = dt;
            nearPE = flash.TotalPE();
          }
        }
      }
    }

    // MC truth match of the slice
    fTMatchIndex = -1;
    fTMatchPur = fTMatchEff = fTMatchCosmicFrac = fTMatchVisEInSlcGeV = fTMatchVisECosmicGeV = -9999.;
    fTrueNuPdg = 0;
    fTrueNuCCNC = fTrueNuMode = -1;
    fTrueNuE = fTrueNuVtxX = fTrueNuVtxY = fTrueNuVtxZ = fTrueNuTimeNs = -9999.;
    fNearTrueNuDtNs = -9999.;
    fTrueNuProdTimeNs = fTrueNuDk2GenCm = fTrueNuGen2VtxCm = fTrueNuG4TimeNs = -9999.;
    fMCTimeOffsetNs = 0.;
    if (fIsMC) {
      double sliceTotal = 0.;
      auto const sliceEnergies = truthEnergies(sliceHitAssns->at(iSlice), sliceTotal);
      if (sliceTotal > 0.) {
        double nuEnergy = 0.;
        art::Ptr<simb::MCTruth> best;
        double bestEnergy = 0.;
        for (auto const& [truth, energy] : sliceEnergies) {
          nuEnergy += energy;
          if (energy > bestEnergy) {
            bestEnergy = energy;
            best = truth;
          }
        }
        fTMatchCosmicFrac = (sliceTotal - nuEnergy) / sliceTotal;
        fTMatchPur = bestEnergy / sliceTotal;
        if (nuTruths.empty())     // CAFMaker: no truths -> all cosmic
          fTMatchVisEInSlcGeV = fTMatchVisECosmicGeV = sliceTotal / 1000.;
        if (best.isNonnull() && fTMatchPur > 0.5) {
          fTMatchVisEInSlcGeV = sliceTotal / 1000.;     // MeV -> GeV
          fTMatchVisECosmicGeV = (sliceTotal - nuEnergy) / 1000.;
          fTMatchIndex =
            std::distance(nuTruths.begin(), std::find(nuTruths.begin(), nuTruths.end(), best));
          if (eventTruthEnergy[best] > 0.) fTMatchEff = bestEnergy / eventTruthEnergy[best];
          fTrueNuProdTimeNs = nuProdTimeNs[fTMatchIndex];
          fTrueNuDk2GenCm = nuDk2GenCm[fTMatchIndex];
          fTrueNuGen2VtxCm = nuGen2VtxCm[fTMatchIndex];
          if (!std::isnan(nuTimeOffsetNs[fTMatchIndex])) {
            fMCTimeOffsetNs = nuTimeOffsetNs[fTMatchIndex];
            if (fOpFlashT0Ns != -9999.) {
              fOpFlashT0Ns += fMCTimeOffsetNs;
              fOpFlashT0CorrectedNs += fMCTimeOffsetNs;
            }
            if (fRawFlashT0Ns != -9999.) fRawFlashT0Ns += fMCTimeOffsetNs;
            if (fBFMFlashTimeNs != -9999.) fBFMFlashTimeNs += fMCTimeOffsetNs;
          }
          if (best->NeutrinoSet()) {
            simb::MCNeutrino const& nu = best->GetNeutrino();
            fTrueNuPdg = nu.Nu().PdgCode();
            fTrueNuCCNC = nu.CCNC();
            fTrueNuMode = nu.Mode();
            fTrueNuE = nu.Nu().E();
            fTrueNuVtxX = nu.Nu().Vx();
            fTrueNuVtxY = nu.Nu().Vy();
            fTrueNuVtxZ = nu.Nu().Vz();
            fTrueNuG4TimeNs = nu.Nu().T();
            fTrueNuTimeNs = nu.Nu().T() + fMCTimeOffsetNs;
            for (unsigned int iNu = 0; iNu < nuTruths.size(); ++iNu) {
              auto const& other = nuTruths[iNu];
              if (other == best || !other->NeutrinoSet()) continue;
              double const otherOffsetNs = std::isnan(nuTimeOffsetNs[iNu]) ? 0. : nuTimeOffsetNs[iNu];
              double const dt = other->GetNeutrino().Nu().T() + otherOffsetNs - fTrueNuTimeNs;
              if (fNearTrueNuDtNs == -9999. || std::abs(dt) < std::abs(fNearTrueNuDtNs))
                fNearTrueNuDtNs = dt;
            }
          }
        }
      }
    }

    fTree->Fill();
  }
}

DEFINE_ART_MODULE(calib::TimingCalibAna)
