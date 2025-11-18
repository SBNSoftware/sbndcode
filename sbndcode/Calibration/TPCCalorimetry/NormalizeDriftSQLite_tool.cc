// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// Tool include
#include "larreco/Calorimetry/INormalizeCharge.h"

// Services
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Lab helpers
#include "wda.h"

// C++
#include <string>
#include <optional>
#include <cassert>

namespace sbnd {
  namespace calo {

class NormalizeDriftSQLite : public INormalizeCharge
{
public:
  NormalizeDriftSQLite(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  void setup(const art::Event& e) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  std::string fDBFileName;
  std::string fDBTag;
  bool fVerbose;

  lariov::DBFolder fDB;

  std::optional<detinfo::DetectorClocksData> fClockData; // need delayed construction

  // Class to hold data from DB
  class RunInfo {
  public:
    double tau_E;
    double tau_W;
  };

  // Helpers
  RunInfo GetRunInfo(uint64_t run);

  // Cache run requests
  std::map<uint32_t, RunInfo> fRunInfos;
};

DEFINE_ART_CLASS_TOOL(NormalizeDriftSQLite)

  } // end namespace calo
} // end namespace sbnd


sbnd::calo::NormalizeDriftSQLite::NormalizeDriftSQLite(fhicl::ParameterSet const &pset):
  fDBFileName(pset.get<std::string>("DBFileName")),
  fDBTag(pset.get<std::string>("DBTag")),
  fVerbose(pset.get<bool>("Verbose", false)),
  fDB(fDBFileName, "", "", fDBTag, true, false)
{}

void sbnd::calo::NormalizeDriftSQLite::configure(const fhicl::ParameterSet& pset) {}

void sbnd::calo::NormalizeDriftSQLite::setup(const art::Event& e) {
  fClockData.emplace(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));
}

sbnd::calo::NormalizeDriftSQLite::RunInfo sbnd::calo::NormalizeDriftSQLite::GetRunInfo(uint64_t run) {
  // check the cache
  if (fRunInfos.count(run)) {
    return fRunInfos.at(run);
  }

  // Look up the run
  //
  // Translate the run into a fake "timestamp"
  fDB.UpdateData((run+1000000000)*1000000000);

  RunInfo thisrun;

  double this_tau_E, this_tau_W;
  fDB.GetNamedChannelData(0, "etau_sce_spatial_east", this_tau_E);
  fDB.GetNamedChannelData(0, "etau_sce_spatial_west", this_tau_W);
  thisrun.tau_E = this_tau_E;
  thisrun.tau_W = this_tau_W;

  if (fVerbose) std::cout << "NormalizeDriftSQLite Tool -- Lifetime Data:" << "\nTPC East: " << thisrun.tau_E << "\nTPC West: " << thisrun.tau_W << std::endl;

  // Set the cache
  fRunInfos[run] = thisrun;

  return thisrun;
}

double sbnd::calo::NormalizeDriftSQLite::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {

  if (!fClockData) {
    std::cout << "Error: fClockData is not valid" << std::endl;
    throw cet::exception("fClockData is not valid");
  }

  // Get the info
  RunInfo runelifetime = GetRunInfo(e.id().runID().run());

  // lookup the TPC
  double thiselifetime = -1;
  unsigned tpc = hit.WireID().TPC;
  unsigned cryo = hit.WireID().Cryostat;

  // East
  if (cryo == 0 && tpc == 0) thiselifetime = runelifetime.tau_E;

  // West
  if (cryo == 0 && tpc == 1) thiselifetime = runelifetime.tau_W;

  // Get the hit time
  double thit = fClockData->TPCTick2TrigTime(hit.PeakTime()) - t0;
  thit = thit * 1.e-3;

  if (fVerbose) std::cout << "NormalizeDriftSQLite Tool -- Norm factor: " << exp(thit / thiselifetime) << " at TPC: " << tpc << " Cryo: " << cryo << " Time: " << thit << " Track T0: " << t0 << ", x: " << location.X() << std::endl;

  // Scale
  if (thiselifetime > 0) {
    dQdx = dQdx*exp(thit / thiselifetime);
  }
  // Throw exception if thiselifetime is not updated to non-zero value
  else {
    std::cout << "sbnd::calo::NormalizeDriftSQLite::Normalize electron lifetime is not found for run " << e.id().runID().run() << std::endl;
    throw cet::exception("Electron lifetime is not found");
  }

  return dQdx;
}

