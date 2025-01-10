#!/usr/bin/env python


__doc__ = """
Collection of utilities to ease interaction with ROOT.

Unsurprisingly, this module requires ROOT.
"""

__all__ = [
  "splitROOTpath", "createROOTpath", "getROOTclass",
  "activateDirectory",
  "ROOT",
  ]

################################################################################
###
### Try to save the command line arguments from unconsiderate ROOT behaviour
###
def ROOTloader():
  """
  The ROOT interpreter likes to peek at the command line arguments and interpret
  them its own way.
  For example, an option `-b` will be interpreted to set ROOT in batch mode.
  Likewise, if a `--help` argument is present on the command line, the
  interpreter will print the ROOT help message and, even worse, halt the script.
  This is clearly not very friendly to the script.
  
  The initialization of the interpreter happens lazily the first time anything
  is read out of `ROOT` module: `import ROOT` is not enough, but pretty much
  everything after that involving ROOT is (exception: `ROOT.gROOT` will not
  trigger the interpreter, but trying to get anything out of `ROOT.gROOT` will).
  That makes it complicate to control when that happens.
  
  This function triggers the interpreter initialization after having removed all
  command line options, and finally restores them (we use a context manager to
  show that we know Python). The loaded module is returned.
  
  This function is called as soon as this module is imported. It is important
  that this happens as early as possible, possibly as a replacement of ROOT
  import, as:
      
      from ROOTutils import ROOT
      
      from ROOTutils import *
      
      import ROOTutils
      import ROOT
      
  or equivalent.
  """
  import sys, logging
  try:              alreadyLoaded = 'gInterpreter' in dir(ROOT)
  except NameError: alreadyLoaded = False
  if alreadyLoaded:
    logging.warning(
      "ROOT module was loaded before ROOTutils.py: command line arguments may be garbled"
      )
    return sys.modules['ROOT']
  # if already loaded 
  
  class EmptyArgs:
    def __enter__(self):
      self.args = sys.argv
      logging.debug("Saving command line: %s", self.args)
      sys.argv = sys.argv[0:1]
      logging.debug(
       "Replaced command line %s with %s before loading ROOT module",
       self.args, sys.argv)
    def __exit__(self, exc_type, exc_value, traceback):
      sys.argv = self.args
      logging.debug("Restored command line %s", sys.argv)
  # class EmptyArgs
  
  with EmptyArgs():
    import ROOT
    ROOT.gInterpreter # make sure that an interpreter is initialized
    return ROOT
  # with
  
# ROOTloader()

ROOT = ROOTloader() # import ROOT...
del ROOTloader

################################################################################
### Print vectors easily
def TVector2ToString(v):
  return "( %g, %g )" % (v.X(), v.Y())
def TVector3ToString(v):
  return "( %g, %g, %g )" % (v.X(), v.Y(), v.Z())
def TLorentzVectorToString(v):
  return "( %g, %g, %g; %g )" % (v.X(), v.Y(), v.Z(), v.T())

ROOT.TVector2.__str__ = TVector2ToString
ROOT.TVector3.__str__ = TVector3ToString
ROOT.TLorentzVector.__str__ = TLorentzVectorToString


################################################################################
###  File management
###  
def splitROOTpath(path):
  """
  Returns the specified path split into file path and ROOT directory path.
  
  The `path` is in the form:
  "/UNIX/path/to/file.root:internal/ROOT/directory/and/object".
  The returned value is a pair `(filePath, dirPath)`: in the example, that
  would be `("/UNIX/path/to/file.root", "internal/ROOT/directory/and/object")`.
  
  Note: for compatibility with some ROOT tradition, the separator ':' can be
  replaced by '/'
  """
  
  # this implementation is not robust, but I am in a hurry :-P
  try:
    filePath, ROOTpath = path.rsplit('.root')
  except ValueError:
    raise RuntimeError("Path '{}' does not include a ROOT file.".format(path))
  filePath += '.root'
  ROOTpath.lstrip(':')
  ROOTpath.lstrip('/')
  return filePath, ROOTpath
  
# splitROOTpath()


def createROOTpath(path, fileMode = "UPDATE"):
  """
  Creates a complete ROOT directory path.
  
  The `path` is in the form:
  "/UNIX/path/to/file.root:internal/ROOT/directory/structure".
  The ROOT file `/UNIX/path/to/file.root` will be created with the specified
  `fileMode` (path may be relative), then the `TDirectoryFile` hierarchy
  `internal/ROOT/directory/structure` will be created under it.
  The return value is a pair `(file, dir)`, where `file` is a open `TFile`
  for `/UNIX/path/to/file.root` and `dir` is the `TDirectory` object of
  `structure`.
  
  Remember to keep track of `file`, or else python may close it compromising
  `dir` as well.
  """
  
  filePath, ROOTpath = splitROOTpath(path)
  
  ROOTfile = ROOT.TFile(filePath, fileMode)
  if not ROOTfile.IsOpen():
    raise RuntimeError \
      ("Can't open ROOT file '{}' in '{}' mode".format(filePath, fileMode))
  
  # instead of using `TDirectory.mkdir()`, we do that manually
  ROOTpathElements = ROOTpath.split('/')
  ROOTdir = ROOTfile
  for ROOTdirName in ROOTpathElements:
    if not ROOTdirName: continue # empty name does nothing
    daughterDir = ROOTdir.GetDirectory(ROOTdirName)
    if not daughterDir:
      daughterDir = ROOTdir.CreateDirectory(ROOTdirName)
    if not daughterDir:
      raise RuntimeError("Can't access directory '{}' under '{}'".format
       (ROOTdirName, ROOTdir.GetPath()))
    ROOTdir = daughterDir
  # for
  
  return ROOTfile, ROOTdir
  
# createROOTpath()


class DirectoryChanger:
  """
  Object changing ROOT directory while on scope.
  
  The purpose is to make a ROOT directory current only as long as it is needed.
  The most typical uses of this objects include the automatic restoration of
  the previous directory as the object falls out of scope.
  Two methods are supported:
  1. function scope:
        
        def writeEverythingInto(dir, everything):
          dirChanger = ROOTutils.DirectoryChanger(dir)
          for item in everything: item.Write()
        # writeEverythingInto()
        
  2. local scope (equivalent to using `activateDirectory()`):
        
        with DirectoryChanger(dir):
          for item in everything: item.Write()
        # with
        
  
  """
  def __init__(self, newDir = None, saveDir = None):
    if saveDir: self.saveDir(saveDir)
    else:       self.saveCurrentDir()
    self.newDir = newDir
    self.changeDir()
  # __init__()
  
  def saveCurrentDir(self): self.saveDir(ROOT.gDirectory)
  
  def saveDir(self, ROOTdir): self.oldDir = ROOTdir
  
  def changeDir(self):
    if self.newDir: self.newDir.cd()
  
  def restoreDir(self):
    if self.oldDir: self.oldDir.cd()
  
  def forget(self): self.oldDir = None
  
  def __del__(self):
    self.restoreDir()
  
  def __enter__(self):
    self.changeDir()
  
  def __exit__(self, exc_type, exc_value, traceback):
    self.restoreDir()
    self.forget()
  
# DirectoryChanger()


def activateDirectory(ROOTdir):
  """
  Sets a directory with `DirectoryChanger`.
  
  Example:
        
        dir = outputFile.GetDirectory("plots")
        with activateDirectory(dir):
          for plot in plots: item.Write()
        # with
        
  """
  return DirectoryChanger(ROOTdir)


################################################################################
def getROOTclass(classPath):
  """Returns the object specified by `classPath` within ROOT module.
  
  Throws `AttributeError` if any object in the path is not available.
  
  Example: `getROOTclass('geo::GeometryCore')` returns `ROOT.geo.GeometryCore`.
  """
  classPath = classPath.replace('::', '.').lstrip('.')
  base = ROOT
  for objName in classPath.split('.'):
    base = getattr(base, objName)
  return base
# getROOTclass()


################################################################################
