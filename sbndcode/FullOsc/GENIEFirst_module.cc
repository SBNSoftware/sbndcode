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

  class GENIEFirst : public art::EDProducer {
  public:

    struct Config {
      fhicl::Atom<art::InputTag> InputTag {
        fhicl::Name("InputTag"),
        fhicl::Comment("InputTag for the GENIEGen module")
      };
    };

    using Parameters = art::EDProducer::Table<Config>;

    explicit GENIEFirst(Parameters const& p);

    void produce(art::Event& evt) override;
    void beginRun(art::Run& run) override;
    void endSubRun(art::SubRun& subrun) override;

  private:
    art::InputTag fInputTag;

    void firstosc(art::Event& evt,
                   std::vector<simb::MCTruth>&          firstTruth,
                   std::vector<simb::MCFlux>&           firstFlux,
                   std::vector<simb::GTruth>&           firstGTruth,
                   art::Assns<simb::MCTruth, simb::MCFlux>& firstTFassn,
                   art::Assns<simb::MCTruth, simb::GTruth>& firstTGTassn,
                   std::vector<sim::BeamGateInfo>&      firstGate) const;
  };


  GENIEFirst::GENIEFirst(Parameters const& p)
    : EDProducer{p}
    , fInputTag{p().InputTag()}
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

  void GENIEFirst::beginRun(art::Run& run)
  {
    art::Handle<sumdata::RunData> hRunData;
    run.getByLabel(fInputTag, hRunData);
    auto fullOscRunData = std::make_unique<sumdata::RunData>(*hRunData);
    run.put(std::move(fullOscRunData), art::fullRun());
  }

  void GENIEFirst::endSubRun(art::SubRun& subrun)
  {
    art::Handle<sumdata::POTSummary> hPOT;
    subrun.getByLabel(fInputTag, hPOT);
    auto fullOscPOT = std::make_unique<sumdata::POTSummary>(*hPOT);
    subrun.put(std::move(fullOscPOT), art::subRunFragment());
  }

  void GENIEFirst::produce(art::Event& evt)
  {
    auto firstTruth   = std::make_unique<std::vector<simb::MCTruth>>();
    auto firstFlux    = std::make_unique<std::vector<simb::MCFlux>>();
    auto firstGTruth  = std::make_unique<std::vector<simb::GTruth>>();
    auto firstTFassn  = std::make_unique<art::Assns<simb::MCTruth, simb::MCFlux>>();
    auto firstTGTassn = std::make_unique<art::Assns<simb::MCTruth, simb::GTruth>>();
    auto firstGate    = std::make_unique<std::vector<sim::BeamGateInfo>>();

    firstosc(evt, *firstTruth, *firstFlux, *firstGTruth,
              *firstTFassn, *firstTGTassn, *firstGate);

    evt.put(std::move(firstTruth));
    evt.put(std::move(firstFlux));
    evt.put(std::move(firstGTruth));
    evt.put(std::move(firstTFassn));
    evt.put(std::move(firstTGTassn));
    evt.put(std::move(firstGate));
  }


  void GENIEFirst::firstosc(art::Event& evt,
                                 std::vector<simb::MCTruth>&          firstTruth,
                                 std::vector<simb::MCFlux>&           firstFlux,
                                 std::vector<simb::GTruth>&           firstGTruth,
                                 art::Assns<simb::MCTruth, simb::MCFlux>& firstTFassn,
                                 art::Assns<simb::MCTruth, simb::GTruth>& firstTGTassn,
                                 std::vector<sim::BeamGateInfo>&      firstGate) const
  {
    art::PtrMaker<simb::MCTruth> makeTruthPtr{evt};
    art::PtrMaker<simb::MCFlux>  makeFluxPtr{evt};
    art::PtrMaker<simb::GTruth>  makeGTruthPtr{evt};

    // ===== Full-osc sample =====
    art::Handle<std::vector<simb::MCTruth>>       hFullOscTruth;
    art::Handle<std::vector<simb::MCFlux>>        hFullOscFlux;
    art::Handle<std::vector<simb::GTruth>>        hFullOscGTruth;
    art::Handle<std::vector<sim::BeamGateInfo>>   hFullOscGate;

    evt.getByLabel(fInputTag, hFullOscTruth);
    evt.getByLabel(fInputTag, hFullOscFlux);
    evt.getByLabel(fInputTag, hFullOscGTruth);
    evt.getByLabel(fInputTag, hFullOscGate);

    if (hFullOscTruth->size() != hFullOscFlux->size() ||
        hFullOscTruth->size() != hFullOscGTruth->size()) {
      throw art::Exception(art::errors::DataCorruption)
        << "GENIEFirst: full-osc collections have mismatched sizes in event "
        << evt.id() << " (MCTruth="  << hFullOscTruth->size()
        << ", MCFlux="               << hFullOscFlux->size()
        << ", GTruth="               << hFullOscGTruth->size() << ")\n";
    }

    if (hFullOscTruth->size() > 1) {
      mf::LogWarning("GENIEFirst")
        << "Full-osc collections have " << hFullOscTruth->size()
        << " entries in event " << evt.id() << "; taking only the first one.";
    }

    // Take index 0 from fullosc collection
    std::size_t idx0 = firstTruth.size();

    firstTruth.emplace_back((*hFullOscTruth)[0]);
    firstFlux.emplace_back((*hFullOscFlux)[0]);
    firstGTruth.emplace_back((*hFullOscGTruth)[0]);
    firstGate.emplace_back((*hFullOscGate)[0]);

    // Associations
    art::Ptr<simb::MCTruth> t0 = makeTruthPtr(idx0);
    art::Ptr<simb::MCFlux>  f0 = makeFluxPtr(idx0);
    art::Ptr<simb::GTruth>  g0 = makeGTruthPtr(idx0);

    firstTFassn.addSingle(t0, f0);
    firstTGTassn.addSingle(t0, g0);


  }

}

DEFINE_ART_MODULE(evgen::GENIEFirst)
