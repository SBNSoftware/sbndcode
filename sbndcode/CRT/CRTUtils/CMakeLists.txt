art_make(
         LIBRARY_NAME sbndcode_CRTUtils

         LIB_LIBRARIES
         messagefacility::MF_MessageLogger
         art::Persistency_Common
         lardataobj::Simulation
         larcore::Geometry_Geometry_service
         nusimdata::SimulationBase
         sbnobj::SBND_CRT
)

install_headers()
install_fhicl()
install_source()
