# the executables are tests from lardata
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

