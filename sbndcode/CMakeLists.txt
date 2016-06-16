
#art_make( )
 
#install_headers()
#install_source()
cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )


# source
add_subdirectory(LArSoftConfigurations)
add_subdirectory(JobConfigurations)
add_subdirectory(Geo)
add_subdirectory(Utilities)
add_subdirectory(Calibration)
add_subdirectory(DetectorSim)
add_subdirectory(AnalysisTree)
add_subdirectory(MCTruthExtractor)
add_subdirectory(SBNDPandora)

