

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

    Parameters detsim_params = get_parameters(argc, argv);


    CLHEP::HepJamesRandom engine(SEED);

    sbnd::crt::CRTDetSimAlg detsim_alg(detsim_params,
                                       engine,
                                       0.);

    detsim_alg.getChannelTriggerTicks(10000, 100, 50);



    return 0;

}


