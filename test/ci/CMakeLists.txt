#
# The intention is to
# 
# 1. keep the integration test scripts in their own directory
# 2. have the tests run in both development (`mrbsetenv`) and run environments
# 
# In order to achieve that:
# 
# 1. there is no way (I know of) with cetbuildtools of installing the scripts
#    into `${SBNDCODE_DIR}/test/ci` directory: `install_scripts(AS_TEST)`
#    installs all directly into `${SBNDCODE_DIR}/test`
# 2. the integration test configuration (`ci_tests.cfg`) stores the command for
#    executing the test scripts, which will be passed to a subshell for
#    execution; there is no way to specify multiple commands, therefore the
#    single command must describe a suitable command for both environments
# 3. `install_source()` installs under `${SBNDCODE_DIR}/source` directory,
#    and the directory `source` has no equivalent in the development area;
#    also, the files installed by `install_source()` are stripped of their
#    executable flag, which the test runner script needs
# 4. given that there seems to be no way to obtain the same formal path for the
#    script in the two environments, we are left with the option of not
#    specifying any path, and let the subshell find the location of the script
#    following its `PATH` settings
# 5. the test scripts can be installed, as anything else, in a location that is
#    guaranteed to be in the execution path (`PATH`) by using
#    `install_scripts()`: this location is `${SBNDCODE_FQ_DIR}/bin` under
#    running environment, and `${MRB_BUILDDIR}/sbndcode/bin` under development
#    environment; to fill the latter directory, `mrb build` will suffice, the
#    installation phase is not even required. The presence of a script in a
#    directory under `PATH` allows the subshell from the test runner script to
#    find it.
# 6. the integration test configuration still needs to be found, and the test
#    runner script expects it to be `test/ci_tests.cfg` under some directory
#    specified by a UPS-like variable `XXX_DIR`
# 7. the integration test configuration needs to be installed with
#    `install_scripts(AS_TEST)`, so that it's saved in `${SBNDCODE_DIR}/test`
#    and available under execution environment
# 8. under development environment, the integration test configuration is not
#    in the magic location, and therefore the configuration file path must be
#    specified explicitly on the command line:
#    `test_runner -c "${SBNDCODE_DIR}/test/ci/ci_tests.cfg"` or equivalent
# 9. the test job configuration files (FHiCL) must also be stored in a location
#    known by the FHiCL parsers, which is achieved by `install_fhicl()`
#    
# Points 1-5 are relevant only when a custom script is needed.
# Point 8 is relevant only until issue #18270 is resolved.
# 

# install_scripts() # see point 5 above
install_scripts(AS_TEST)

install_fhicl()
