
# Test fcl files.

cet_test( testFHiCLfiles_sbnd.sh PREBUILT
  TEST_ARGS "--quiet" "--no-color" "--exclude-from=SkipTesting.txt"
  DATAFILES testFHiCLfiles.sh SkipTesting.txt
  )

install_scripts()
