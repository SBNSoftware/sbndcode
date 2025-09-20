install_headers()
install_fhicl()
install_source()

file(GLOB calib_values_crt *.txt)
install_fw( LIST ${calib_values_crt} )

art_make(SERVICE_LIBRARIES 
                           art::Framework_Services_Registry
                           art::Framework_Principal
                           art::Framework_Core
                           art::Persistency_Provenance
                           messagefacility::MF_MessageLogger
                           ROOT::Core
                           sbnobj::SBND_CRT
)
