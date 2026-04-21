////////////////////////////////////////////////////////////////////////////////
// Class:       CathodeSimEnergyDepositFilter
// Plugin Type: producer (art v3_xx / LArSoft v09+)
//
// Reads one or more sim::SimEnergyDeposit collections, drops every deposit
// whose MidPoint (and, optionally, Start or End) lies inside the data-driven
// SBND cathode volume (sbnd::SBNDCathode), and writes the surviving deposits
// back out under this module's label with user-specified instance names.
//
// Configuring multiple (input, output-instance) pairs in a single job mirrors
// ionandscint's behaviour of emitting both "" (default instance) and
// "priorSCE" collections.  With a single-element list the module behaves
// like a simple 1:1 filter.
//
// FHiCL parameters
// -----------------
//   SimEnergyDepositLabels (vector<InputTag>)
//       default { "ionandscint" }
//       input collections, one per produced output.
//   OutputInstances        (vector<string>)
//       default { "" }
//       instance name used for each produced collection.  Must have the
//       same length as SimEnergyDepositLabels.
//   CathodeVolumeFile      (string)   required
//   UseStepEndpoints       (bool)     default false
//   NanPolicy              (string)   default "outside" ("outside"|"inside"|"throw")
//   Interp                 (string)   default "bilinear" ("bilinear"|"nearest")
//   Verbose                (bool)     default false
//
// Produces
// --------
//   std::vector<sim::SimEnergyDeposit>   (one per OutputInstances entry)
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "sbndcode/CathodeFilter/SBNDCathode.h"

#include <memory>
#include <string>
#include <vector>

namespace sbnd {

class CathodeSimEnergyDepositFilter : public art::EDProducer {
public:
  explicit CathodeSimEnergyDepositFilter(fhicl::ParameterSet const& p);

  CathodeSimEnergyDepositFilter(CathodeSimEnergyDepositFilter const&)            = delete;
  CathodeSimEnergyDepositFilter(CathodeSimEnergyDepositFilter&&)                 = delete;
  CathodeSimEnergyDepositFilter& operator=(CathodeSimEnergyDepositFilter const&) = delete;
  CathodeSimEnergyDepositFilter& operator=(CathodeSimEnergyDepositFilter&&)      = delete;

  void produce(art::Event& e) override;

private:
  std::vector<art::InputTag> fInputLabels;
  std::vector<std::string>   fOutputInstances;
  std::string   fCathodeVolumeFile;
  bool          fUseStepEndpoints;
  std::string   fNanPolicyStr;
  std::string   fInterpStr;
  bool          fVerbose;

  std::unique_ptr<SBNDCathode> fCathode;

  bool depositInsideCathode(sim::SimEnergyDeposit const& d) const;
};

CathodeSimEnergyDepositFilter::CathodeSimEnergyDepositFilter(fhicl::ParameterSet const& p)
  : EDProducer(p)
  , fInputLabels      (p.get<std::vector<art::InputTag>>("SimEnergyDepositLabels",
                                                        { art::InputTag("ionandscint") }))
  , fOutputInstances  (p.get<std::vector<std::string>>("OutputInstances",
                                                      { std::string{} }))
  , fCathodeVolumeFile(p.get<std::string>("CathodeVolumeFile"))
  , fUseStepEndpoints (p.get<bool>       ("UseStepEndpoints", false))
  , fNanPolicyStr     (p.get<std::string>("NanPolicy",        "outside"))
  , fInterpStr        (p.get<std::string>("Interp",           "bilinear"))
  , fVerbose          (p.get<bool>       ("Verbose",          false))
{
  if (fInputLabels.empty())
    throw cet::exception("CathodeSimEnergyDepositFilter")
      << "SimEnergyDepositLabels must contain at least one entry";
  if (fInputLabels.size() != fOutputInstances.size())
    throw cet::exception("CathodeSimEnergyDepositFilter")
      << "SimEnergyDepositLabels and OutputInstances must have the same length "
      << "(" << fInputLabels.size() << " vs " << fOutputInstances.size() << ")";

  fCathode = std::make_unique<SBNDCathode>(fCathodeVolumeFile);

  if      (fNanPolicyStr == "outside") fCathode->setNanPolicy(SBNDCathode::NanPolicy::Outside);
  else if (fNanPolicyStr == "inside" ) fCathode->setNanPolicy(SBNDCathode::NanPolicy::Inside);
  else if (fNanPolicyStr == "throw"  ) fCathode->setNanPolicy(SBNDCathode::NanPolicy::Throw);
  else throw cet::exception("CathodeSimEnergyDepositFilter")
         << "Unknown NanPolicy '" << fNanPolicyStr << "'";

  if      (fInterpStr == "bilinear") fCathode->setInterp(SBNDCathode::Interp::Bilinear);
  else if (fInterpStr == "nearest" ) fCathode->setInterp(SBNDCathode::Interp::Nearest);
  else throw cet::exception("CathodeSimEnergyDepositFilter")
         << "Unknown Interp '" << fInterpStr << "'";

  mf::LogInfo("CathodeSimEnergyDepositFilter") << fCathode->describe();
  for (size_t i = 0; i < fInputLabels.size(); ++i) {
    mf::LogInfo("CathodeSimEnergyDepositFilter")
      << "filter [" << i << "]: "
      << fInputLabels[i].encode()
      << "  ->  <this>:\"" << fOutputInstances[i] << "\"";
    produces<std::vector<sim::SimEnergyDeposit>>(fOutputInstances[i]);
  }
}

bool CathodeSimEnergyDepositFilter::depositInsideCathode(
    sim::SimEnergyDeposit const& d) const {
  const auto mid = d.MidPoint();
  if (fCathode->contains(mid.X(), mid.Y(), mid.Z())) return true;

  if (fUseStepEndpoints) {
    if (fCathode->contains(d.StartX(), d.StartY(), d.StartZ())) return true;
    if (fCathode->contains(d.EndX(),   d.EndY(),   d.EndZ()  )) return true;
  }
  return false;
}

void CathodeSimEnergyDepositFilter::produce(art::Event& e) {
  for (size_t i = 0; i < fInputLabels.size(); ++i) {
    auto const& in =
      e.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fInputLabels[i]);

    auto out = std::make_unique<std::vector<sim::SimEnergyDeposit>>();
    out->reserve(in->size());

    size_t nKept = 0, nDropped = 0;
    double eKept = 0.0, eDropped = 0.0;
    for (sim::SimEnergyDeposit const& dep : *in) {
      if (depositInsideCathode(dep)) {
        ++nDropped;
        eDropped += dep.Energy();
        continue;
      }
      out->push_back(dep);
      ++nKept;
      eKept += dep.Energy();
    }

    if (fVerbose) {
      mf::LogInfo("CathodeSimEnergyDepositFilter")
        << "[" << fInputLabels[i].encode() << "  ->  :"
        << fOutputInstances[i] << "] "
        << "kept " << nKept << " (" << eKept << " MeV) / "
        << "dropped " << nDropped << " (" << eDropped << " MeV) / "
        << "total " << in->size();
    }

    e.put(std::move(out), fOutputInstances[i]);
  }
}

} // namespace sbnd

DEFINE_ART_MODULE(sbnd::CathodeSimEnergyDepositFilter)
