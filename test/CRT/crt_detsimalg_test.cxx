

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


int main(int argc, char const** argv) {

    int iParam = 0;

    // first argument: configuration file (mandatory)
    std::cout << "simpleCRT_test SetConfigurationPath" << std::endl;
    std::string config_path;
    if (++iParam < argc) config_path = argv[iParam];


    char const* fhicl_env = getenv("FHICL_FILE_PATH");
    std::string search_path = fhicl_env? std::string(fhicl_env) + ":": ".:";
    testing::details::FirstAbsoluteOrLookupWithDotPolicy policy(search_path);

    // parse a configuration file; obtain intermediate form
    fhicl::intermediate_table table;
    table = fhicl::parse_document(argv[iParam], policy);

    // translate into a parameter set
    fhicl::ParameterSet params;
    params = fhicl::ParameterSet::make(table);

    std::cout << "simpleCRT_test " << params.to_string() << std::endl;


    long seed = 0;
    CLHEP::HepJamesRandom engine(seed);

    fhicl::ParameterSet p = params.get<fhicl::ParameterSet>("testcrtsim");

    std::cout << "Just before CRTDetSimAlg:" << std::endl;
    std::cout << "p " << p.to_string() << std::endl;
    std::cout << "p.get<fhicl::ParameterSet>(DetSimParams) " << p.get<fhicl::ParameterSet>("DetSimParams").to_string() << std::endl;

    using Parameters = fhicl::Table<sbnd::crt::CRTDetSimParams>;
    Parameters detsim_params(p.template get<fhicl::ParameterSet>("DetSimParams"));



    sbnd::crt::CRTDetSimAlg detsim_alg(detsim_params,
                                       engine,
                                       0.);



    return 0;

}


