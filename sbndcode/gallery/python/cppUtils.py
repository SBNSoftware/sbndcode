#!/usr/bin/env python

from __future__ import print_function

__doc__ = """
Collection of utilities to interface C++ code with Python via PyROOT.

This module requires ROOT.
"""

__all__ = [
  'readHeader',
  'SourceCode',
  ]

import sys, os
from ROOTutils import ROOT

# Make sure <span> is not included in the range v3 library
ROOT.gROOT.ProcessLine('#define RANGES_WORKAROUND_MSVC_UNUSABLE_SPAN 1')

################################################################################
def readHeader(headerPath):
    """Make the ROOT C++ jit compiler read the specified header."""
    ROOT.gROOT.ProcessLine('#include "%s"' % headerPath)
# readHeader()


################################################################################
class SourceCentral:
  """
  A class keeping track of the sources and where to look for them.
  """
  AllPlatformInfo = {
    'Linux': {
      'Name':       'Linux',
      'LibSuffix':  '.so',
      'LibEnvPath': 'LD_LIBRARY_PATH',
    },
    'Darwin': {
      'Name':       'Darwin',
      'LibSuffix':  '.dylib',
      'LibEnvPath': 'DYLD_LIBRARY_PATH', # might be not honoured
    },
  } # AllPlatformInfo
  PlatformInfo = AllPlatformInfo[os.uname()[0]]
  
  def __init__(self, *includePaths):
    self.headers = {}
    self.libraries = {}
    self.includePaths = []
    self.addIncPaths(*includePaths)
    # for
  # __init__()
  
  def addIncPath(self, path, force=False):
    expPath = os.path.expandvars(path)
    if not os.path.isdir(expPath):
      print(
        "Warning: include path '%s'" % path,
        (" ( => '%s')" % expPath if path != expPath else ""),
        " does not exist.",
        sep='',
        file=sys.stderr
        )
    if force or expPath not in self.includePaths:
      self.includePaths.append(expPath)
  # addIncPath()
  
  def addIncPathEnv(self, varName, force = False):
    self.addIncPath(os.environ[varName], force=force)
  
  def addIncPaths(self, *paths):
    for path in paths: self.addIncPath(path)

  def addIncPathEnvs(self, *varNames):
    self.addIncPaths(*map((lambda varName: os.environ[varName]), varNames))

  def find(self, relPath, extraPaths = []):
    return self.findLibrary(relPath, extraPaths=extraPaths) if self.isLibrary(relPath) else self.findHeader(relPath, extraPaths=extraPaths)
  # find()
  
  def findLibrary(self, libName, extraPaths = []):
    expLibName = SourceCentral.expandLibraryName(libName)
    for path in reversed(
     SourceCentral.LibraryPaths() + list(map(os.path.expandvars, extraPaths))
     ):
      candidate = os.path.join(path, expLibName)
      if os.path.exists(candidate): return candidate
    else: return None
  # findLibrary()
  
  def findHeader(self, relPath, extraPaths = []):
    for path in reversed(self.includePaths + list(map(os.path.expandvars, extraPaths))):
      candidate = os.path.join(path, relPath)
      if os.path.exists(candidate): return candidate
    else: return None
  # findHeader()
  
  def loadLibrary(self, relPath, extraPaths = [], force = False):
    expandedName = self.expandLibraryName(relPath)
    res = ROOT.gSystem.Load(expandedName)
    if res == 0: self.libraries[relPath] = expandedName
    return res
  # loadLibrary()
  
  def loadHeader(self, headerRelPath, extraPaths = [], force = False):
    try: return self.headers[headerRelPath]
    except KeyError: pass
    headerPath = self.findHeader(headerRelPath, extraPaths=extraPaths)
    if not headerPath: raise RuntimeError("Can't locate header file '%s'" % headerRelPath)
    readHeader(headerPath)
    self.headers[headerRelPath] = headerPath
    return headerPath
  # loadHeader()
  
  
  def loadHeaderFromUPS(self, headerRelPath, extraPaths = [], force = False):
    """
    Loads a C++ header from a UPS product.
    
    Assumptions:
    * the specified relative path of the header is under the include directory
      of its UPS product
    * the include directory path is set in a environment variable named with
      the standard UPS pattern (`PRODUCTNAME_INC`)
    * the header relative path starts with a directory that reflects the name
      of the UPS product, `productname/relative/package/path/header.h`
    
    For example, for a `headerRelPath` of `larcorealg/Geometry/GeometryCore.h`,
    the full path must be represented by
    `${LARCOREALG_INC}/larcorealg/Geometry/GeometryCore.h`, with the content
    of `LARCOREALG_INC` variable being an absolute path.
    """
    # make sure that if there is a INC variable for the package, that one is included
    return self.loadHeader(
      headerRelPath,
      extraPaths
        =([ '$' + self.packageVarNameFromHeaderPath('INC', headerRelPath) ] + extraPaths),
      force=force
      )
  # loadHeaderFromUPS()

  def load(self, relPath, extraPaths = [], force = False):
    return (self.loadLibrary if self.isLibrary(relPath) else self.loadHeaderFromUPS)(relPath, extraPaths=extraPaths, force=force)
  # load()
  
  def isLibrary(self, path):
    return os.path.splitext(path)[-1] in [ self.PlatformInfo['LibSuffix'], '' ]
  
  def expandLibraryName(self, name):
    if not name.startswith('lib'): name = 'lib' + name
    LibSuffix = self.PlatformInfo['LibSuffix']
    if not name.endswith(LibSuffix): name += LibSuffix
    return name
  # expandLibraryName()
  
  @staticmethod
  def packageNameFromHeaderPath(headerPath):
    return os.path.split(os.path.dirname(headerPath))[0]
  
  @staticmethod
  def packageVarNameFromHeaderPath(varSuffix, headerPath):
    return SourceCentral.packageNameFromHeaderPath(headerPath).upper() + '_' + varSuffix
  
  @staticmethod
  def LibraryPaths():
    return os.getenv(SourceCentral.PlatformInfo['LibEnvPath']) \
     .split(SourceCentral.PlatformInfo.get('LibEnvPathSep', ':'))
  # LibraryPaths()
  
# class SourceCentral

################################################################################

# global instance of source tracking class
SourceCode = SourceCentral()

################################################################################
