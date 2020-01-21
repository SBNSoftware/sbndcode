# We follow the pattern of `icarusutil`:
# we in principle distinguish executable and non-executable (module) python files.
# The former are installed with execute permission.
# The latter are installed without execute permission.
# In practice, everything that is not explicitly listed as executable is
# considered module.

# Executable python files (install in a PATH-accessible directory CET-knows-where)

set (exes) # none so far
if(exes)
  message(STATUS "Executable python scripts: ${exes}")
  install_scripts(LIST ${exes})
endif(exes)

# Non-executable python files (install in our custom PYTHONPATH-accessible directory)
file(GLOB nonexes RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
if(nonexes AND exes)
  list(REMOVE_ITEM nonexes ${exes})
endif(nonexes AND exes)
if(nonexes)
  message(STATUS "Non-executable python modules: ${nonexes}")
  install(FILES ${nonexes} DESTINATION "${product}/${version}/python")
endif(nonexes)
