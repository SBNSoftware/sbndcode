cet_build_plugin(PMTGainFluctuations1Dynode art::tool 
      SOURCE
            PMTGainFluctuations1Dynode_tool.cc 
      LIBRARIES 
            nurandom::RandomUtils_NuRandomService_service
)


cet_build_plugin(PMTNonLinearityTF1 art::tool 
      SOURCE
            PMTNonLinearityTF1_tool.cc 
      LIBRARIES
            ROOT::Hist
)


install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
