
cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} ${CMAKE_CURRENT_LIST_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})

cet_test(CheckGENIEGeomScan HANDBUILT
     TEST_EXEC python
     TEST_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/check_genie_geom_scan.py
     DATAFILES check_genie_geom_scan.py
)

