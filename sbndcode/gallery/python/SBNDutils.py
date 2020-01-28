#!/usr/bin/env python


__doc__ = """
Collection of utilities to interface SBND with python, gallery and LArSoft.

This module requires ROOT.
"""

__all__ = [
  'loadSBNDgeometry',
  'justLoadSBNDgeometry',
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
  
  SourceCode.loadHeaderFromUPS('sbndcode/Geometry/ChannelMapSBNDAlg.h')
  SourceCode.loadLibrary('sbnd_Geometry')
  return LArSoftUtils.loadGeometry \
    (config=config, registry=registry, mapping=ROOT.geo.ChannelMapSBNDAlg)
# loadSBNDgeometry()


def justLoadSBNDgeometry(configFile, mapping = None):
  """Loads and returns SBND geometry from the specified configuration file.
  
  This is a one-stop procedure recommended only when running interactively.
  """
  if mapping is not None:
    raise NotImplementedError("Support for non-standard mapping not implemented yet.")
  return loadSBNDgeometry(config=LArSoftUtils.ConfigurationClass(configFile))
# justLoadSBNDgeometry()


################################################################################
