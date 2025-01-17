
set(
        MODULE_LIBRARIES
                larcorealg::Geometry
                larcore::Geometry_Geometry_service
                larsim::Simulation
                lardataobj::Simulation
                larsim::MCCheater_BackTrackerService_service
                larsim::MCCheater_ParticleInventoryService_service
                lardata::Utilities
                larevt::Filters
                lardataobj::RawData
                lardataobj::RecoBase
                lardata::RecoObjects
                larpandora::LArPandoraInterface
                nusimdata::SimulationBase
                sbnobj::Common_CRT
                sbnobj::SBND_CRT
                sbndcode_CRTUtils
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art::Persistency_Common
                art::Persistency_Provenance
                art::Utilities
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                ROOT::Geom
                ROOT::XMLIO
                ROOT::Gdml
                ROOT::Core
                sbndcode_RecoUtils
        )
cet_build_plugin(CRTTrigFilter art::module SOURCE CRTTrigFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(GenCRTFilter art::module SOURCE GenCRTFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(GenNuFilter art::module SOURCE GenNuFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(LArG4CRTFilter art::module SOURCE LArG4CRTFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(LArG4FakeTriggerFilter art::module SOURCE LArG4FakeTriggerFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})
cet_build_plugin(SimEnergyDepFakeTriggerFilter art::module SOURCE SimEnergyDepFakeTriggerFilter_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()

add_subdirectory(fcls)
