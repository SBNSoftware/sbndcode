#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include <memory>
#include <vector>
#include <string>

namespace evgen {

  class GENIEOscPileup : public art::EDProducer {
  public:

    struct Config {
      fhicl::Atom<art::InputTag> FullOscTag {
        fhicl::Name("FullOscTag"),
        fhicl::Comment("InputTag for the full-osc GENIEGen module")
      };
      fhicl::Atom<art::InputTag> UnOscTag {
        fhicl::Name("UnOscTag"),
        fhicl::Comment("InputTag for the un-osc GENIEGen module")
      };
    };

    using Parameters = art::EDProducer::Table<Config>;

    explicit GENIEOscPileup(Parameters const& p);

    void produce(art::Event& evt) override;
    void beginRun(art::Run& run) override;
    void endSubRun(art::SubRun& subrun) override;

  private:
    art::InputTag fFullOscTag;
    art::InputTag fUnOscTag;

    void pileuposc(art::Event& evt,
                   std::vector<simb::MCTruth>&          pileupTruth,
                   std::vector<simb::MCFlux>&           pileupFlux,
                   std::vector<simb::GTruth>&           pileupGTruth,
                   art::Assns<simb::MCTruth, simb::MCFlux>& pileupTFassn,
                   art::Assns<simb::MCTruth, simb::GTruth>& pileupTGTassn,
                   std::vector<sim::BeamGateInfo>&      pileupGate) const;
  };


  GENIEOscPileup::GENIEOscPileup(Parameters const& p)
    : EDProducer{p}
    , fFullOscTag{p().FullOscTag()}
    , fUnOscTag{p().UnOscTag()}
  {
    produces<std::vector<simb::MCTruth>>();
    produces<std::vector<simb::MCFlux>>();
    produces<std::vector<simb::GTruth>>();
    produces< sumdata::RunData, art::InRun >();
    produces< sumdata::POTSummary, art::InSubRun >();
    produces<art::Assns<simb::MCTruth, simb::MCFlux>>();
    produces<art::Assns<simb::MCTruth, simb::GTruth>>();
    produces<std::vector<sim::BeamGateInfo>>();
  }

  void GENIEOscPileup::beginRun(art::Run& run)
  {
    art::Handle<sumdata::RunData> hRunData;
    run.getByLabel(fFullOscTag, hRunData);
    auto fullOscRunData = std::make_unique<sumdata::RunData>(*hRunData);
    run.put(std::move(fullOscRunData), art::fullRun());
  }

  void GENIEOscPileup::endSubRun(art::SubRun& subrun)
  {
    art::Handle<sumdata::POTSummary> hPOT;
    subrun.getByLabel(fFullOscTag, hPOT);
    auto fullOscPOT = std::make_unique<sumdata::POTSummary>(*hPOT);
    subrun.put(std::move(fullOscPOT), art::subRunFragment());
  }

  void GENIEOscPileup::produce(art::Event& evt)
  {
    auto pileupTruth   = std::make_unique<std::vector<simb::MCTruth>>();
    auto pileupFlux    = std::make_unique<std::vector<simb::MCFlux>>();
    auto pileupGTruth  = std::make_unique<std::vector<simb::GTruth>>();
    auto pileupTFassn  = std::make_unique<art::Assns<simb::MCTruth, simb::MCFlux>>();
    auto pileupTGTassn = std::make_unique<art::Assns<simb::MCTruth, simb::GTruth>>();
    auto pileupGate    = std::make_unique<std::vector<sim::BeamGateInfo>>();

    pileuposc(evt, *pileupTruth, *pileupFlux, *pileupGTruth,
              *pileupTFassn, *pileupTGTassn, *pileupGate);

    evt.put(std::move(pileupTruth));
    evt.put(std::move(pileupFlux));
    evt.put(std::move(pileupGTruth));
    evt.put(std::move(pileupTFassn));
    evt.put(std::move(pileupTGTassn));
    evt.put(std::move(pileupGate));
  }


  void GENIEOscPileup::pileuposc(art::Event& evt,
                                 std::vector<simb::MCTruth>&          pileupTruth,
                                 std::vector<simb::MCFlux>&           pileupFlux,
                                 std::vector<simb::GTruth>&           pileupGTruth,
                                 art::Assns<simb::MCTruth, simb::MCFlux>& pileupTFassn,
                                 art::Assns<simb::MCTruth, simb::GTruth>& pileupTGTassn,
                                 std::vector<sim::BeamGateInfo>&      pileupGate) const
  {
    art::PtrMaker<simb::MCTruth> makeTruthPtr{evt};
    art::PtrMaker<simb::MCFlux>  makeFluxPtr{evt};
    art::PtrMaker<simb::GTruth>  makeGTruthPtr{evt};

    // ===== Full-osc sample =====
    art::Handle<std::vector<simb::MCTruth>>       hFullOscTruth;
    art::Handle<std::vector<simb::MCFlux>>        hFullOscFlux;
    art::Handle<std::vector<simb::GTruth>>        hFullOscGTruth;
    art::Handle<std::vector<sim::BeamGateInfo>>   hFullOscGate;

    evt.getByLabel(fFullOscTag, hFullOscTruth);
    evt.getByLabel(fFullOscTag, hFullOscFlux);
    evt.getByLabel(fFullOscTag, hFullOscGTruth);
    evt.getByLabel(fFullOscTag, hFullOscGate);

    if (hFullOscTruth->size() != hFullOscFlux->size() ||
        hFullOscTruth->size() != hFullOscGTruth->size()) {
      throw art::Exception(art::errors::DataCorruption)
        << "GENIEOscPileup: full-osc collections have mismatched sizes in event "
        << evt.id() << " (MCTruth="  << hFullOscTruth->size()
        << ", MCFlux="               << hFullOscFlux->size()
        << ", GTruth="               << hFullOscGTruth->size() << ")\n";
    }

    if (hFullOscTruth->size() > 1) {
      mf::LogWarning("GENIEOscPileup")
        << "Full-osc collections have " << hFullOscTruth->size()
        << " entries in event " << evt.id() << "; taking only the first one.";
    }

    // Take index 0 from fullosc collection
    std::size_t idx0 = pileupTruth.size();

    pileupTruth.emplace_back((*hFullOscTruth)[0]);
    pileupFlux.emplace_back((*hFullOscFlux)[0]);
    pileupGTruth.emplace_back((*hFullOscGTruth)[0]);
    pileupGate.emplace_back((*hFullOscGate)[0]);

    // Associations
    art::Ptr<simb::MCTruth> t0 = makeTruthPtr(idx0);
    art::Ptr<simb::MCFlux>  f0 = makeFluxPtr(idx0);
    art::Ptr<simb::GTruth>  g0 = makeGTruthPtr(idx0);

    pileupTFassn.addSingle(t0, f0);
    pileupTGTassn.addSingle(t0, g0);

    // ===== Un-osc sample =====
    art::Handle<std::vector<simb::MCTruth>>       hUnOscTruth;
    art::Handle<std::vector<simb::MCFlux>>        hUnOscFlux;
    art::Handle<std::vector<simb::GTruth>>        hUnOscGTruth;

    evt.getByLabel(fUnOscTag, hUnOscTruth);
    evt.getByLabel(fUnOscTag, hUnOscFlux);
    evt.getByLabel(fUnOscTag, hUnOscGTruth);

    if (hUnOscTruth->size() != hUnOscFlux->size() ||
        hUnOscTruth->size() != hUnOscGTruth->size()) {
      throw art::Exception(art::errors::DataCorruption)
        << "GENIEOscPileup: un-osc collections have mismatched sizes in event "
        << evt.id() << " (MCTruth="  << hUnOscTruth->size()
        << ", MCFlux="               << hUnOscFlux->size()
        << ", GTruth="               << hUnOscGTruth->size() << ")\n";
    }

    // Append all unoscillated neutrinos after first neutrino (1..N)
    for (std::size_t i = 1; i < hUnOscTruth->size(); ++i) {
      std::size_t idx = pileupTruth.size();

      pileupTruth.emplace_back((*hUnOscTruth)[i]);
      pileupFlux.emplace_back((*hUnOscFlux)[i]);
      pileupGTruth.emplace_back((*hUnOscGTruth)[i]);
      art::Ptr<simb::MCTruth> t = makeTruthPtr(idx);
      art::Ptr<simb::MCFlux>  f = makeFluxPtr(idx);
      art::Ptr<simb::GTruth>  g = makeGTruthPtr(idx);

      pileupTFassn.addSingle(t, f);
      pileupTGTassn.addSingle(t, g);
    }

  }

}

DEFINE_ART_MODULE(evgen::GENIEOscPileup)
