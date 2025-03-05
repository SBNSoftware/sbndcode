#!/usr/bin/env python


__doc__ = """
Collection of utilities to interface SBND with python, gallery and LArSoft.

This module requires ROOT.
"""

__all__ = [
  'loadSBNDgeometry',
  'loadWireReadout',
  'loadSBNDauxDetgeometry',
]

import LArSoftUtils
from ROOTutils import ROOT


################################################################################
### Geometry
###
def loadSBNDgeometry(config = None, registry = None):
  """Loads and returns SBND geometry with the standard SBND channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  SourceCode = LArSoftUtils.SourceCode # alias
  
  SourceCode.loadHeaderFromUPS('sbndcode/Geometry/GeoObjectSorterSBND.h')
  SourceCode.loadLibrary('sbndcode_Geometry')
  return LArSoftUtils.loadGeometry \
    (config=config, registry=registry, sorterClass=ROOT.geo.GeoObjectSorterSBND)
# loadSBNDgeometry()


def loadSBNDwireReadout(config = None, registry = None):
  """Loads and returns SBND wire readout with the standard SBND channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  assert registry, "Registry is required" # because we'll load Geometry from it
  SourceCode = LArSoftUtils.SourceCode # alias
  
  SourceCode.loadHeaderFromUPS('sbndcode/Geometry/WireReadoutSorterSBND.h')
  SourceCode.loadLibrary('sbndcode_Geometry') # should be already loaded by now
  return LArSoftUtils.loadWireReadout \
    (config=config, registry=registry, sorterClass=ROOT.geo.WireReadoutSorterSBND)
# loadSBNDwireReadout()


def loadSBNDauxDetgeometry(config = None, registry = None):
  """Loads and returns SBND geometry with the standard SBND channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  SourceCode = LArSoftUtils.SourceCode # alias
  
  # SourceCode.loadHeaderFromUPS('sbndcode/CRT/CRTGeoObjectSorter.h')
  # SourceCode.loadLibrary('sbndcode_CRTData')
  # SourceCode.loadHeaderFromUPS('sbndcode/CRT/CRTAuxDetInitializerSBND.h')
  # SourceCode.loadLibrary('libsbndcode_CRT_CRTAuxDetInitializerSBND_tool')

  return LArSoftUtils.loadAuxDetGeometry(config=config, registry=registry)
#    , auxDetReadoutInitClass=ROOT.sbnd.crt.CRTAuxDetInitializerSBND)


################################################################################
