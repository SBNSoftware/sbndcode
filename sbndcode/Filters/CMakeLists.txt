
art_make(
        MODULE_LIBRARIES
                larcorealg_Geometry
                larcore_Geometry_Geometry_service
                larsim_Simulation lardataobj_Simulation
                larsim_MCCheater_BackTrackerService_service
                larsim_MCCheater_ParticleInventoryService_service
                lardata_Utilities
                larevt_Filters
                lardataobj_RawData
                lardataobj_RecoBase
                lardata_RecoObjects
                larpandora_LArPandoraInterface
                nusimdata::SimulationBase
#                sbndcode_CRT
#                sbndcode_CRTData
                sbnobj_Common_CRT
                sbnobj_SBND_CRT
                sbndcode_CRTUtils
                art::Framework_Core
                ${ART_FRAMEWORK_MODULES}
                art::Framework_Principal
                art::Framework_Services_Registry
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                
                art::Persistency_Common canvas
                art::Persistency_Provenance canvas
                art::Utilities canvas
                messagefacility::MF_MessageLogger
                
                fhiclcpp::fhiclcpp
                ROOT::Geom
                ROOT::XMLIO
                ROOT::Gdml
                ${ROOT_BASIC_LIB_LIST}
                sbndcode_RecoUtils
        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(fcls)
