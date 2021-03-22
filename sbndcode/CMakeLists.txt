cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

add_subdirectory(JobConfigurations)
add_subdirectory(LArSoftConfigurations)
add_subdirectory(Utilities)
add_subdirectory(Geometry)
add_subdirectory(Calibration)
add_subdirectory(CRT)
add_subdirectory(DetectorSim)
add_subdirectory(MCTruthExtractor)
add_subdirectory(SBNDPandora)
add_subdirectory(SimulationFilters)
add_subdirectory(AnalysisTree)
add_subdirectory(RecoUtils)
add_subdirectory(OpDetSim)
add_subdirectory(OpDetReco)
add_subdirectory(OpT0Finder)
add_subdirectory(CosmicId)
add_subdirectory(ShowerAna)
add_subdirectory(Commissioning)

# subdirectory for gallery stuff
add_subdirectory(gallery)

# for Space Charge
add_subdirectory(SpaceChargeServices)
add_subdirectory(SpaceCharge)
add_subdirectory(SpaceChargeTest)

add_subdirectory(FlashMatch)

add_subdirectory(Filters)
