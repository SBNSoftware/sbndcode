// Framework Includes
//#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"


#include "fhiclcpp/ParameterSet.h"

#include "larevt/CalibrationDBI/Providers/DBFolder.h"

#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larreco/Calorimetry/INormalizeCharge.h"
// Services
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Lab helpers
//#include "wda.h"

// C++
#include <string>
#include <optional>
#include <cassert>
#include <map>
#include <cstdint>

namespace sbnd {
  namespace calo {

    class NormalizeDriftSQLite : public INormalizeCharge
    {
    public:
      NormalizeDriftSQLite(fhicl::ParameterSet const &pset);

      void configure(const fhicl::ParameterSet& pset) override;
      void setup(const art::Event& e) override;
      double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;
      // Class to hold data from DB
      class RunInfo {
      public:
    double tau_E;
    double tau_W;
      };

      // Helpers
      RunInfo GetRunInfo(uint64_t run);

    private:
      // Configuration
      std::string fDBFileName;
      std::string fDBTag;
      bool fVerbose;

      lariov::DBFolder fDB;

      std::optional<detinfo::DetectorClocksData> fClockData; // need delayed construction

      // Cache run requests
      std::map<uint32_t, RunInfo> fRunInfos;
    };

    //DEFINE_ART_CLASS_TOOL(NormalizeDriftSQLite)

      } // end namespace calo
} // end namespace sbnd