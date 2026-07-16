include(CetTest)
cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} ${CMAKE_CURRENT_LIST_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR} ${sbnobj_BINARY_DIR}/lib)
cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR} ${sbndcode_BINARY_DIR}/scripts)
cet_test_env("SBNDCODE_DIR=${PROJECT_SOURCE_DIR}")

# for PREBUILT tests, the test name is the executable name, and no argument is allowed
cet_test(fast_simulation_reconstruction_test_electrons.sh PREBUILT)
cet_test(fast_simulation_reconstruction_test_muons.sh PREBUILT)

cet_test(prodcosmics_cry_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config prodcosmics_cry_sbnd.fcl -n 5
  )
