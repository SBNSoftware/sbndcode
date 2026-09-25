cet_build_plugin(POTDump art::module
  LIBRARIES PRIVATE
    sbnobj::Common_POTAccounting
    larcoreobj::SummaryData
    art::Framework_Core
    art::Framework_Principal
    canvas::canvas
    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
)

install_fhicl()
install_source()
