#!/usr/bin/env python2
#
# This script changes C++ code and CMake files to use sbndcode/Geomery
# in place of the old sbndcode/Geo .
# 
# Change log:
# 20170405 (petrillo@fnal.gov)
#   original version
#

import sys, re

import SerialSubstitution
from SerialSubstitution import AddProcessor, RunSubstitutor


################################################################################
if __name__ == "__main__":

   #############################################################################
   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # CMakeLists.txt
   #
   Subst = AddProcessor(SerialSubstitution.ProcessorClass("cmake"))

   Subst.AddFileNamePattern("CMakeLists.txt")

   Subst.AddWord         ("sbndcode_Geo", "sbndcode_Geometry")

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # C++ source code (including modules and services)
   #
   Subst = AddProcessor(SerialSubstitution.ProcessorClass("code"))

   Subst.AddFileType("h", "cc", "cpp", "cxx")

   Subst.AddWord         ("sbndcode/Geo/", "sbndcode/Geometry/")

   # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   #############################################################################

   sys.exit(RunSubstitutor())
# 
