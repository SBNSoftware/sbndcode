/**
 * \brief Basic unit tests for the CRT detector simulation
 *
 * \author Marco Del Tutto
 */

#include <iostream> // for output before message facility is set up
#include <string>
#include <memory> // std::unique_ptr<>
#include <utility> // std::move(), std::forward()
#include <map>
#include <type_traits> // std::add_rvalue_reference()
#include <stdexcept> // std::logic_error

#include "larcorealg/TestUtils/unit_test_base.h"
#include "CLHEP/Random/JamesRandom.h"
#include "sbndcode/CRT/CRTSimulation/CRTDetSimAlg.h"
#include "sbndcode/CRT/CRTSimulation/CRTDetSimParams.h"

constexpr long SEED = 0;

using Parameters = fhicl::Table<sbnd::crt::CRTDetSimParams>;

Parameters get_parameters(int argc, char const** argv) {

    int iParam = 0;

    // first argument: configuration file (mandatory)
    std::string config_path;
    if (++iParam < argc) config_path = argv[iParam];


    char const* fhicl_env = getenv("FHICL_FILE_PATH");
    std::string search_path = fhicl_env? std::string(fhicl_env) + ":": ".:";
    testing::details::FirstAbsoluteOrLookupWithDotPolicy policy(search_path);

    fhicl::intermediate_table table;
    table = fhicl::parse_document(argv[iParam], policy);

    // translate into a parameter set
    fhicl::ParameterSet params;
    params = fhicl::ParameterSet::make(table);

    fhicl::ParameterSet p = params.get<fhicl::ParameterSet>("testcrtsim");

    // std::cout << "p " << p.to_string() << std::endl;

    Parameters detsim_params(p.template get<fhicl::ParameterSet>("DetSimParams"));

    return detsim_params;
}

int main(int argc, char const** argv) {

    int errors = 0;

    //
    // Read parameters
    //
    Parameters detsim_params = get_parameters(argc, argv);


    //
    // Create DetSimAlg class
    //
    CLHEP::HepJamesRandom engine(SEED);

    sbnd::crt::CRTDetSimAlg detsim_alg(detsim_params,
                                       engine,
                                       0.);

    //
    // Time response
    //
    uint32_t ts = detsim_alg.getChannelTriggerTicks(10000, // true time, ns
                                                    30, // PEs
                                                    50); // distance to readout

    std::cout << "ts " << ts << std::endl;

    // if (ts) {
    //     std::cout << "Signal saturated with a 1.5 MeV deposited energy?" << std::endl;
    //     errors++;
    // }


    //
    // Charge response
    //
    long npe0, npe1;
    double q0, q1;
    detsim_alg.ChargeResponse(1.5*1e-3, // 1.5 MeV deposited energy
                              5, // 5 cm distance to fiber 0
                              5, // 5 cm distance to fiber 1
                              50, // 50 cm distance to readout
                              npe0, npe1, q0, q1);

    std::cout << "npe0 " << npe0
              << "\nnpe1 " << npe1
              << "\nq0 " << q0
              << "\nq1 " << q1 << std::endl;

    if (q0 == detsim_alg.Params().AdcSaturation() or q1 == detsim_alg.Params().AdcSaturation()) {
        std::cout << "Signal saturated with a 1.5 MeV deposited energy?" << std::endl;
        errors++;
    }

    if (q0 < 0 or q1 < 0) {
        std::cout << "Predicted ADCs are negative?" << std::endl;
        errors++;
    }


    //
    // Waveform emulation
    //

    double original_adc = 2000;
    uint16_t adc = detsim_alg.WaveformEmulation(10, // time delay
                                                original_adc); // adc

    std::cout << "original_adc " << original_adc
              << "\nadc " << adc << std::endl;

    if (adc > static_cast<uint16_t>(original_adc)) {
        std::cout << "ADC value after waveform emulation is larger?" << std::endl;
        errors++;
    }




    return errors;

}


