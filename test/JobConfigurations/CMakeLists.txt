# for PREBUILT tests, the test name is the executable name, and no argument is allowed
cet_test(fast_simulation_reconstruction_test_electrons.sh PREBUILT)
cet_test(fast_simulation_reconstruction_test_muons.sh PREBUILT)

cet_test(prodcosmics_cry_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config prodcosmics_cry_sbnd.fcl -n 5
  )
