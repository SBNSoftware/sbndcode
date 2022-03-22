art_make(   LIB_LIBRARIES  larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation
                           lardataobj_Simulation
                           lardataobj_RawData
                           lardataobj_RecoBase
                           lardata_RecoObjects
                           lardataalg_DetectorInfo
                           nusimdata::SimulationBase
                           cetlib::cetlib cetlib_except
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art_root_io::tfile_support ROOT::Core
                           art_root_io::TFileService_service
                           art::Persistency_Common canvas
                           art::Persistency_Provenance canvas
                           art::Utilities canvas
                           messagefacility::MF_MessageLogger

                           fhiclcpp::fhiclcpp
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ${ROOT_BASIC_LIB_LIST}

                           sbndcode_CRT
                           sbnobj_Common_CRT
                           sbnobj_SBND_CRT
                           sbndcode_GeoWrappers
                           sbndcode_RecoUtils
        )

install_headers()
install_fhicl()
install_source()
