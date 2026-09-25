////////////////////////////////////////////////////////////////////////
// Class:       POTDump
// Plugin Type: analyzer
// File:        POTDump_module.cc
//
// Dumps the exposure products of the SubRun to CSV, whichever are present:
//   sumdata::POTSummary            (MC generator)       -> <prefix>_pot.csv
//   std::vector<sbn::BNBSpillInfo> (SBNDBNBRetriever)   -> <prefix>_spills.csv
//   std::vector<sbn::EXTCountInfo> (SBNDBNBEXTRetriever)-> <prefix>_gates.csv
// and every event it sees                               -> <prefix>_events.csv
//
// Runs in the same job as the retriever producers (end path), so the reco1
// file never has to be rewritten, or alone over MC reco1 files.
// EXTCountInfo carries no event id: the retriever appends one entry per
// event in order, so it is joined to the events of the subrun by position.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace sbnd {
  class POTDump;
}

class sbnd::POTDump : public art::EDAnalyzer {
public:
  explicit POTDump(fhicl::ParameterSet const& p);

  void beginSubRun(art::SubRun const& sr) override;
  void analyze(art::Event const& e) override;
  void endSubRun(art::SubRun const& sr) override;

private:
  art::InputTag fPOTSummaryTag;
  art::InputTag fBNBSpillTag;
  art::InputTag fEXTCountTag;
  std::string fPrefix;

  std::ofstream fEvents, fPOT, fSpills, fGates;
  std::vector<art::EventID> fEventIDs;   // the events of the current subrun, in order

  // opened on first use, so a job leaves only the files it has rows for
  std::ofstream& out(std::ofstream& f, std::string const& name, std::string const& header);
};

sbnd::POTDump::POTDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fPOTSummaryTag{p.get<std::string>("POTSummaryLabel")}
  , fBNBSpillTag{p.get<std::string>("BNBSpillInfoLabel")}
  , fEXTCountTag{p.get<std::string>("EXTCountInfoLabel")}
  , fPrefix{p.get<std::string>("OutputPrefix")}
{}

std::ofstream& sbnd::POTDump::out(std::ofstream& f, std::string const& name, std::string const& header)
{
  if (!f.is_open()) {
    f.open(fPrefix + "_" + name + ".csv");
    f << header << '\n';
  }
  return f;
}

void sbnd::POTDump::beginSubRun(art::SubRun const&)
{
  fEventIDs.clear();
}

void sbnd::POTDump::analyze(art::Event const& e)
{
  fEventIDs.push_back(e.id());
  out(fEvents, "events", "run,subrun,event") << e.run() << ',' << e.subRun() << ',' << e.event() << '\n';
}

void sbnd::POTDump::endSubRun(art::SubRun const& sr)
{
  auto const run = sr.run();
  auto const subrun = sr.subRun();

  if (auto h = sr.getHandle<sumdata::POTSummary>(fPOTSummaryTag); h.isValid()) {
    out(fPOT, "pot", "run,subrun,events,totpot,totgoodpot,totspills,goodspills")
      << run << ',' << subrun << ',' << fEventIDs.size() << ','
      << std::setprecision(12) << h->totpot << ',' << h->totgoodpot << ','
      << h->totspills << ',' << h->goodspills << '\n';
  }

  if (auto h = sr.getHandle<std::vector<sbn::BNBSpillInfo>>(fBNBSpillTag); h.isValid()) {
    auto& f = out(fSpills, "spills", "run,subrun,event,spill_time_s,spill_time_ns,TOR860,TOR875,FOM");
    for (auto const& s : *h) {
      f << run << ',' << subrun << ',' << s.event << ',' << s.spill_time_s << ',' << s.spill_time_ns << ','
        << std::setprecision(9) << s.TOR860 << ',' << s.TOR875 << ',' << s.FOM << '\n';
    }
  }

  if (auto h = sr.getHandle<std::vector<sbn::EXTCountInfo>>(fEXTCountTag); h.isValid()) {
    if (h->size() != fEventIDs.size()) {
      mf::LogWarning("POTDump") << "subrun " << run << "/" << subrun << ": " << h->size()
                                << " EXTCountInfo entries for " << fEventIDs.size() << " events";
    }
    auto& f = out(fGates, "gates", "run,subrun,event,gates");
    for (std::size_t i = 0; i < h->size(); ++i) {
      long const event = i < fEventIDs.size() ? static_cast<long>(fEventIDs[i].event()) : -1;
      f << run << ',' << subrun << ',' << event << ',' << std::setprecision(9)
        << (*h)[i].gates_since_last_trigger << '\n';
    }
  }
}

DEFINE_ART_MODULE(sbnd::POTDump)
