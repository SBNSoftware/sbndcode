file(GLOB pythonfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
message(STATUS "Installing python test files: ${pythonfiles}")
install_scripts(LIST ${pythonfiles})
install_scripts()


cet_test(test_opt0finder.py PREBUILT)

#find_package(PythonInterp REQUIRED)

#cet_test(opt0finder_fhicl_check HANDBUILT
#         TEST_EXEC ${PYTHON_EXECUTABLE}
#         TEST_ARGS test_opt0finder.py
#        )



