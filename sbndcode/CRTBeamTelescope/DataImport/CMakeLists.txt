cet_enable_asserts()

simple_plugin( ImportCRTData
               "source"
               sbndcode_CRTBeamTelescope_DataImport
               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
)

art_make( EXCLUDE ImportCRTData_source.cc
          LIB_LIBRARIES nusimdata::SimulationBase
                        sbnobj_SBND_CRT
                        larcoreobj::SummaryData
                        lardata::Utilities
                        art::Framework_Core
                        art::Framework_IO_Sources
                        art::Framework_Principal

                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        art::Persistency_Common canvas::canvas
                        art::Persistency_Provenance canvas::canvas
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp

                        ${ROOT_BASIC_LIB_LIST}
                        cetlib::cetlib
                        cetlib_except::cetlib_except

          MODULE_LIBRARIES sbnobj_SBND_CRT
                        nusimdata::SimulationBase
                        larcoreobj::SummaryData
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}

                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        fhiclcpp::fhiclcpp
                        # GENIE STUFF #####################################
                        xml2
                        log4cpp::log4cpp
                        GENIE::GFwEG
                        GENIE::GFwAlg
                        GENIE::GFwMsg
                        GENIE::GFwGHEP
                        GENIE::GTlFlx
                        GENIE::GFwParDat
                        ${GBASE}
                        ${GEVGCORE}
                        ${GEVGDRIVERS}
                        ${GFLUXDRIVERS}
                        ${GINTERACTION}
                        ${DK2NUDATA}
                        ${DK2NUGENIE}
                        ####################################################
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_source()
install_fhicl()
