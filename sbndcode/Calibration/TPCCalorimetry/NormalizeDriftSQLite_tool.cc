#include "sbndcode/Calibration/TPCCalorimetry/NormalizeDriftSQLite_tool.h"

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

