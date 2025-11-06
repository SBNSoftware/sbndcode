cet_build_plugin(PMTGainFluctuations1Dynode art::tool 
      SOURCE
            PMTGainFluctuations1Dynode_tool.cc 
      LIBRARIES 
            nurandom::RandomUtils_NuRandomService_service
)

cet_build_plugin(PMTGaussianGainFluctuation art::tool 
      SOURCE
            PMTGaussianGainFluctuation_tool.cc 
      LIBRARIES 
            nurandom::RandomUtils_NuRandomService_service
            larcore::Geometry_Geometry_service
)

cet_build_plugin(PMTNonLinearityTF1 art::tool 
      SOURCE
            PMTNonLinearityTF1_tool.cc 
      LIBRARIES
            ROOT::Hist
)

#cet_build_plugin(PMTNonLinearityTF1ChannelByChannel art::tool 
#      SOURCE
#            PMTNonLinearityTF1ChannelByChannel_tool.cc 
#      LIBRARIES
#            art::Framework_Services_Registry
#            ROOT::Hist
#            larcore::Geometry_Geometry_service
#)

install_headers()
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)