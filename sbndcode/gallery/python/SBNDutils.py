#!/usr/bin/env python


__doc__ = """
Collection of utilities to interface SBND with python, gallery and LArSoft.

This module requires ROOT.
"""

__all__ = [
  'loadSBNDgeometry',
  'justLoadSBNDgeometry',
  'loadWireReadout',
  'justLoadSBNDwireReadout',
  'loadSBNDauxDetgeometry',
  'justLoadSBNDauxDetgeometry',
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
    (config=config, registry=registry, sorter=ROOT.geo.GeoObjectSorterSBND)
# loadSBNDgeometry()


def justLoadSBNDgeometry(configFile, mapping = None):
  """Loads and returns SBND geometry from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  if mapping is not None:
    raise NotImplementedError("Support for non-standard mapping not implemented yet.")
  return loadSBNDgeometry(config=LArSoftUtils.ConfigurationClass(configFile))
# justLoadSBNDgeometry()


def loadSBNDwireReadout(config = None, registry = None):
  """Loads and returns SBND wire readout with the standard SBND channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  SourceCode = LArSoftUtils.SourceCode # alias
  
  SourceCode.loadHeaderFromUPS('sbndcode/Geometry/GeoObjectSorterSBND.h')
  SourceCode.loadHeaderFromUPS('sbndcode/Geometry/WireReadoutSorterSBND.h')
  SourceCode.loadLibrary('sbndcode_Geometry')
  return LArSoftUtils.loadWireReadout \
    (config=config, registry=registry, mapping=ROOT.geo.WireReadoutSorterSBND, sorter=ROOT.geo.GeoObjectSorterSBND)
# loadSBNDwireReadout()


def justLoadSBNDwireReadout(configFile, mapping = None):
  """Loads and returns SBND wire readout from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  if mapping is not None:
    raise NotImplementedError("Support for non-standard mapping not implemented yet.")
  return loadSBNDwireReadout(config=LArSoftUtils.ConfigurationClass(configFile))
# justLoadSBNDwireReadout()


def loadSBNDauxDetgeometry(config = None, registry = None):
  """Loads and returns SBND geometry with the standard SBND channel mapping.
  
  See `loadGeometry()` for the meaning of the arguments.
  """
  SourceCode = LArSoftUtils.SourceCode # alias
  
  SourceCode.loadHeaderFromUPS('sbndcode/CRT/CRTGeoObjectSorter.h')
  SourceCode.loadHeaderFromUPS('sbndcode/CRT/CRTAuxDetInitializer.h')
  SourceCode.loadLibrary('sbndcode_CRTData')

  return LArSoftUtils.loadAuxDetGeometry \
    (config=config, registry=registry, sorter=None, auxdetinit=ROOT.sbnd.crt.CRTAuxDetInitializer) #auxdetinit=ROOT.geo.AuxDetReadoutGeom) 
    # (config=config, registry=registry, sorter=ROOT.geo.AuxDetGeoObjectSorterStandard, auxdetinit=ROOT.sbnd.crt.CRTAuxDetInitializer) #auxdetinit=ROOT.geo.AuxDetReadoutGeom) 
# loadSBNDauxDetgeometry()


def justLoadSBNDauxDetgeometry(configFile, mapping = None):
  """Loads and returns SBND geometry from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  if mapping is not None:
    raise NotImplementedError("Support for non-standard mapping not implemented yet.")
  return loadSBNDauxDetgeometry(config=LArSoftUtils.ConfigurationClass(configFile))
# justLoadSBNDauxDetgeometry()

################################################################################
