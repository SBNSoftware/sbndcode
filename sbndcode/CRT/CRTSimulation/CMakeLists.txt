art_make_library( SOURCE CRTDetSimAlg.cxx
                LIBRARIES  larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           lardataobj::Simulation
                           lardataobj::RawData
                           lardataobj::RecoBase
                           lardata::RecoObjects
                           lardataalg::DetectorInfo
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           ROOT::Core
                           ROOT::MathMore
                           sbnobj::SBND_CRT
        )

install_headers()
install_fhicl()
install_source()
