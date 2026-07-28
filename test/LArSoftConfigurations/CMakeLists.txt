# the executables are tests from lardata
cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} )
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})

cet_test( DetectorClocksSBND_test
  HANDBUILT
  DATAFILES dettest_sbnd.fcl
  TEST_EXEC DetectorClocksStandard_test
  TEST_ARGS ./dettest_sbnd.fcl
  OPTIONAL_GROUPS BROKEN
)

cet_test( LArPropertiesSBND_test
  HANDBUILT
  DATAFILES dettest_sbnd.fcl
  TEST_EXEC LArPropertiesStandard_test
  TEST_ARGS ./dettest_sbnd.fcl
  OPTIONAL_GROUPS BROKEN
)

cet_test( DetectorPropertiesSBND_test
  HANDBUILT
  DATAFILES dettest_sbnd.fcl
  TEST_EXEC DetectorPropertiesStandard_test
  TEST_ARGS ./dettest_sbnd.fcl
  OPTIONAL_GROUPS BROKEN
)

