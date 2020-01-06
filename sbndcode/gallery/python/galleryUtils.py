#!/usr/bin/env python


__doc__ = """
Collection of utilities to interface gallery with python.

This module requires ROOT.
"""

__all__ = [
  'readHeader', # imported from `cppUtils`
  'SourceCode', # imported from `cppUtils`
  'make_getValidHandle',
  'makeFileList',
  'forEach',
  'eventLoop',
  'findFHiCL',
  'loadConfiguration',
  'ConfigurationClass',
  'startMessageFacility',
  'ServiceRegistryClass',
  'ConfigurationHelper',
  ]

import sys, os
from ROOTutils import ROOT
import cppUtils
import warnings


################################################################################
### Pass-through
### 
SourceCode = cppUtils.SourceCode
readHeader = cppUtils.readHeader


################################################################################
###  gallery
################################################################################

class HandleMaker:
  """Make the ROOT C++ jit compiler instantiate the
  Event::getValidHandle member template for template parameter klass.
  
  Needs to keep track of what was done already, because Cling will protest if
  the same class is asked twice.
  """
  AlreadyMade = set()
  
  def __call__(self, klass):
    if klass in HandleMaker.AlreadyMade: return
    res = HandleMaker.make(klass)
    if res != ROOT.TInterpreter.kNoError:
      raise RuntimeError(
       "Could not create `ROOT.gallery.Event.getValidHandle` for '%s' (code: %d)"
       % (klass, res)
       )
    # if
    HandleMaker.AlreadyMade.add(klass)
  # __call__()
  
  @staticmethod
  def make(klass):
    ROOT.gROOT.ProcessLine('template gallery::ValidHandle<%(name)s> gallery::Event::getValidHandle<%(name)s>(art::InputTag const&) const;' % {'name' : klass})
  
# class HandleMaker
make_getValidHandle = HandleMaker()



def makeFileList(*filePaths):
  """Creates a file list suitable for `gallery::Event`."""
  files = ROOT.vector(ROOT.string)()
  for path in filePaths: files.push_back(path)
  return files
# makeFileList()


def forEach(event):
  """
  Simplifies sequential event looping.
  
  This function handles the state of `event`, moving the current event around.
  Therefore, it has side effects.
  The function returns `event` itself, after it has set to the next available
  event.
  
  This function is actually a generator.
  
  Example of loop:
      
      for iEvent, event in enumerate(forEach(event)):
        ...
      
  
  """
  while (not event.atEnd()):
    yield event
    event.next()
  # while
# forEach()


class EventIterator:
  """
  Iterator for an event sequence.
  
  It can be used directly as:
      
      for event in EventIterator(event): ...
      
  or it can be plugged in directly in the gallery.Event class, making the latter
  be iterable:
      
      for event in event: ...
      
  (or `for evt in event: ...`).
  """
  def __init__(self, event):
    self._event = event
    self._index = None
  def __iter__(self):
    self._index = None
    return self
  def next(self):
    if self._index is not None:
      self._event.next()
      self._index += 1
    if self._event.atEnd(): raise StopIteration
    return self._event
  # __next__
  
# class EventIterator


def eventLoop(inputFiles,
 process,
 options = {},
 ):
  """
  Applies the `process` function to each and every event from the specified
  input files, in sequence.
  
  The `inputFiles` list may be a single file, or a list (or any iterable object)
  of files, or a `std::vector<std::string>` object.
  
  The `process` callable is executed with two arguments: the number of argument
  in the loop, and the event itself. No information on which file the event is
  taken from is provided. If a call returns exactly `False`, it is considered
  to have failed and an error counter is incremented. Exceptions raised in
  `process` are not handled.
  
  The error counter is returned at the end of the execution.
  
  Options:
  - 'nEvents': number of events to be processed (does not include skipped ones)
  - 'nSkip': number of events from the beginning of the sample to be skipped
  """
  
  # option reading
  nSkip = options.get('nSkip', 0)
  nEvents = options.get('nEvents', None)
  
  # make sure the input file list is in the right format
  if not isinstance(inputFiles, ROOT.vector(ROOT.string)):
    if isinstance(inputFiles, str): inputFiles = [ inputFiles, ]
    inputFiles = makeFileList(*inputFiles)
  # if
  
  event = ROOT.gallery.Event(inputFiles)
  
  # ROOT.gStyle.SetOptStat(0)
  
  iEvent = 0
  nProcessedEvents = 0
  nErrors = 0
  
  iFile = None
  for iEvent, event in enumerate(forEach(event)):
    
    if iFile != event.fileEntry():
      iFile = event.fileEntry()
      print "Opening: '%s'" % inputFiles[iFile]
    # if new file
    
    # event flow control
    if iEvent < nSkip: continue
    if (nEvents is not None) and (nProcessedEvents >= nEvents): break
    nProcessedEvents += 1
    
    ###
    ###
    ###
    res = process(event, iEvent)
    if isinstance(res, bool) and not res: nErrors += 1
    
    ###
    ### all done
    ###
    
  # for
  if nErrors > 0:
    print >>sys.stderr, "Encountered %d/%d errors." % (nErrors, nProcessedEvents)
  return nErrors
# eventLoop()



# this does not really work...
# ROOT.gallery.Event.__iter__ = lambda self: EventIterator(self)


################################################################################
### Infrastructure
################################################################################

def findFHiCL(configRelPath, extraDirs = []):
  
  if os.path.isfile(configRelPath):
    return os.path.join(os.getcwd(), configRelPath)
  for path in extraDirs + os.environ.get('FHICL_FILE_PATH', "").split(':'):
    candidate = os.path.join(path, configRelPath)
    if os.path.isfile(candidate): return candidate
  else: return None

# findFHiCL()


################################################################################
class ConfigurationHelper:
  """
  Provides more user-friendly interface for the configuration access.
  
  While FHiCL C++ interface is friendly enough, some of that is lost in Python
  bindings.
  
  This helper class keeps a reference to the FHiCL configuration object, and
  provides some `get()` overloads to make the interface easier.
  """
  def __init__(self, config):
    self.config = config
  def get(self, key, default=None, klass=None):
    """
    Returns the value associated with the specified key.
    
    Parameters
    -----------
    
    key _(string)_
      the key of the configuration parameter to be read
    default
      the value to provide if no key is present; if default is `None`,
      no default is present and if key is missing a `KeyError` exception will
      be raised
    klass _(class type)_
      type of object returned (if `None`, it's taken from `default` if
      specified)
    
    Returns
    --------
    
    An object of type `klass`, initialised with the value associated to `key`
    or with `default`.
    
    Raises
    -------
    
    `KeyError` if the `key` is not present in the configuration and `default` is
    specified as `None`.
    
    """
    if klass is None:
      assert default is not None
      klass = default.__class__
    # first try to directly return the key
    try: return self.config.get(klass)(key)
    except Exception: pass
    # if that failed, act depending on the default value
    if default is None: raise KeyError(key)
    return klass(default)
  # get()
  
  def has(self, key):
    """
    Returns whether there is a value associated with the specified key.
    
    Parameters
    -----------
    
    key _(string)_
      the key of the configuration parameter to be read
    
    """
    return self.config.has(klass)
  # has()
  
  __call__ = get # alias
  
  def pset(self): return self.config
  
# class ConfigurationHelper


class TemporaryFile:
  def __init__(self, data = None):
    with warnings.catch_warnings():
      # os.tempnam() is a potential security risk: ACK
      warnings.filterwarnings("ignore", ".*tempnam .*", RuntimeWarning)
      self._file = open(os.tempnam(), "w+")
    self.name = self._file.name
    if data is not None:
      self._file.write(str(data))
      self._file.flush() # we are not going to close this file...
    # 
  # __init__()
  def __del__(self):
    if not self._file: return
    del self._file
    os.remove(self.name)
  # __del__()
  def file_(self): return self._file
  def __str__(self): return self.name
# class TemporaryFile


def loadConfiguration(configSpec):
  # this utility actually relies on generic utilities that while not LArSoft
  # specific, are nevertheless distributed with LArSoft (`larcorealg`).
  SourceCode.loadHeaderFromUPS("larcorealg/Geometry/StandaloneBasicSetup.h")
  
  if isinstance(configSpec, ConfigurationString):
    configFile = TemporaryFile(configSpec)
    configPath = configFile.name
  else:
    configFile = None
    configPath = configSpec
  # if
  
  fullPath = findFHiCL(configPath)
  if not fullPath:
    raise RuntimeError("Couldn't find configuration file '%s'" % configPath)
  return ROOT.lar.standalone.ParseConfiguration(fullPath)
# loadConfiguration()


class ConfigurationString:
  """Wrapper to a string that should be considered configuration text."""
  def __init__(self, config): self.config = config
  def __str__(self): return self.config
# class ConfigurationString


class ConfigurationClass:
  def __init__(self, configPath):
    self.config = loadConfiguration(configPath)
  def paramsFor(self, FHiCLpath):
    return self.config.get(ROOT.fhicl.ParameterSet)(FHiCLpath)
  def service(self, serviceName):
    return self.paramsFor("services." + serviceName)
  def producer(self, moduleName):
    return self.paramsFor("physics.producers." + moduleName)
# class ConfigurationClass


################################################################################
class ServiceRegistryClass:
  def __init__(self, config):
    self.fullConfig = config if isinstance(config, ConfigurationClass) else ConfigurationClass(config)
    self.services = {}
  # __init__(self)
  
  def config(self, serviceName): return self.fullConfig.service(serviceName)
  
  def has(self, serviceName): return serviceName in self.services
  
  def register(self, serviceName, service):
    self.services[serviceName] = service
    return service
  
  def registeredServiceNames(self): return self.services.keys()
  
  def registry(self): return self # behave like a manager (LArSoftUtils concept)
  
  def get(self, serviceName): return self.services[serviceName]
  def __call__(self, serviceName): return self.get(serviceName)
  
  def create(self, serviceName, serviceClass, *otherServiceConstructorArgs):
    serviceConfig = self.config(serviceName)
    if not serviceConfig:
      raise RuntimeError("Couldn't find the configuration for service '%s'" % serviceName)
    service = serviceClass(serviceConfig, *otherServiceConstructorArgs)
    if not service:
      raise RuntimeError("Failed to create service '%s' (type: %s)" % (serviceName, serviceClass.__class__.__name__))
    return self.register(serviceName, service)
  # create()
  
# class ServiceRegistryClass


################################################################################
class startMessageFacility:
  """Use it as a function: `startMessageFacility(config, applName)`.
  
  `config` either is `ConfigurationClass` object, or it is a parameter set
  with the configuration for the service.
  """
  Init = False
  def __init__(self, config, applName = None):
    if not startMessageFacility.Init: self.init(config, applName)
  def init(self, config, applName):
    if not applName: applName = os.path.basename(sys.argv[0])
    if isinstance(config, ConfigurationClass): config = config.service("message")
    print "Starting message facility for %s..." % applName
    ROOT.mf.StartMessageFacility(config, applName)
    startMessageFacility.Init = True
  # init()
# class startMessageFacility


################################################################################
