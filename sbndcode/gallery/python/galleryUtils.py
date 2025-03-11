#!/usr/bin/env python

from __future__ import print_function

__doc__ = """
Collection of utilities to interface gallery with python.

This module requires ROOT.

An example of a interactive session counting the number of muons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import galleryUtils
import ROOT
sampleEvents = galleryUtils.makeEvent("data.root")
LArG4tag = ROOT.art.InputTag("largeant")

getParticleHandle \\
  = galleryUtils.make_getValidHandle("std::vector<simb::MCParticle>", sampleEvents)
for event in galleryUtils.forEach(sampleEvents):
  particles = getParticleHandle(LArG4tag).product()
  
  nMuons = sum(1 for part in particles if abs(part.PdgCode()) == 13)
  print("%s: %d muons" % (event.eventAuxiliary().id(), nMuons))
  
# for all events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
from ROOTutils import ROOT, expandFileList
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

# override conversion of art::EventID into a string (used by `print()`)
ROOT.art.EventID.__str__ \
  = lambda self: "R:%d S:%d E:%d" % (self.run(), self.subRun(), self.event())


class HandleMaker:
  """Make the ROOT C++ jit compiler instantiate the
  Event::getValidHandle member template for template parameter klass.
  
  Needs to keep track of what was done already, because Cling will protest if
  the same class is asked twice.
  
  See the documentation of `__call__` for examples of usage of instances of
  this class.
  
  """
  class ManyByTypeProc:
    def __init__(self, klass):
      self.handleVectorClass = ROOT.std.vector[f'gallery::Handle<{klass}>']
      self.getHandles = ROOT.gallery.Event.getManyByType[klass]
    def __call__(self, event):
      handles = self.handleVectorClass()
      self.getHandles(event, handles)
      return handles
  # ManyByTypeProc
  
  
  def validHandle(self,
    klass: "data product class to prepare for (as a C++ class name string)",
    event: "(optional) the event object to retrieve the data products from " = None,
    ) -> "if `event` is specified, bound `getValidHandle<klass>`, unbound otherwise":
    # this has been tested with ROOT 6.22;
    # big improvements in cppyy make this **way** simpler than it used to be
    getHandle = ROOT.gallery.Event.getValidHandle[klass]
    return (lambda tag: getHandle(event, tag)) if event else getHandle
  # validHandle()
  
  def manyByType(self,
    klass: "data product class to prepare for (as a C++ class name string)",
    event: "(optional) the event object to retrieve the data products from " = None,
    ) -> """if `event` is specified, a list (std::vector) of handles,
            otherwise a callable that when called on an event returns such list""":
    getManyByType = HandleMaker.ManyByTypeProc(klass)
    return getManyByType(event) if event else (lambda event: getManyByType(event))
  # manyByType()
  
  def __call__(self,
    klass: "data product class to prepare for (as a C++ class name string)",
    event: "(optional) the event object to retrieve the data products from " = None,
    ) -> "if `event` is specified, bound `getValidHandle<klass>`, unbound otherwise":
    """Prepares for reading data products of the specified class.
    
    This function causes the instantiation of the `gallery::Event::getValidHandle()`
    method specialized to read the specified `klass`.
    It also returns a bound or unbound function to actually retrieve data products.
    
    If an `event` instance is specified, the returned value is a bound function
    that can be used directly.
        
        from galleryUtils import makeEvent, make_getValidHandle
        import ROOT
        
        event = makeEvent("test.root")
        inputTag = ROOT.art.InputTag("largeant")
        
        getParticlesFrom = make_getValidHandle("std::vector<simb::MCParticle>")
        # note: the following is currently mostly broken (see below)
        particles1 = getParticlesFrom(event, inputTag).product()
        
        getParticles = make_getValidHandle("std::vector<simb::MCParticle>", event)
        particles2 = getParticles(inputTag).product()
        
    Both `particles1` and `particles2` point to the same data product.
    
    Exception (`RuntimeError`) is raised on error.
    """
    return self.validHandle(klass, event)
  # __call__()
  
# class HandleMaker
make_getValidHandle = HandleMaker()


def makeFileList(
 *filePaths: "a list of input files"
 ) -> "a list of files suitable to construct a `gallery.Event object":
  """Creates a file list suitable for `gallery::Event`.
  
  If a file ends with `.root`, it is added directly to the list.
  Otherwise, it is interpreted as a file list and treated as such
  (see `ROOTutils.expandFileList()`).
  File list recursion is disabled.
  """
  files = ROOT.vector(ROOT.string)()
  for path in filePaths:
    entries = [ path ] if path.endswith('.root') else expandFileList(path)
    for entry in entries: files.push_back(entry)
  # for
  return files
# makeFileList()


def makeEvent(
 files: "files or file lists, as supported by `ROOTutils.makeFileList()`",
 **kwargs: "additional arguments to `gallery::Event` constructor"
 ) -> "a gallery.Event object reading the specified input files":
  return ROOT.gallery.Event(makeFileList(files), **kwargs)
# makeEvent()


def forEach(
  event,
  fromStart: "always restarts from the beginning of the dataset" = True,
  maxEvents: "stop after processing this many events (None = all)" = None,
  skipEvents: "skip these many events at the beginning of the dataset" = 0,
  ):
  """
  Simplifies sequential event looping.
  
  This function handles the state of `event`, moving the current event around.
  Therefore, it has side effects.
  The function returns `event` itself, after it has set to the next available
  event.
  Normally the function restarts from the beginning of the event list.
  If `fromStart` is `False`, though, it will start after the last event
  processed by `event`.
  A number of entries to skip (either from the start or from the current event,
  depending on the value of `fromStart`) and a maximum number of events to
  process can be specified. However, there is no feedback on whether processing
  stopped because of a reached limit or because of no more events are available.
  
  This function is actually a generator.
  
  Example of loop:
      
      for iEvent, event in enumerate(forEach(event, maxEvents=5)):
        ...
      
  will process up to 5 events from the start of the sample;
      
      for iEvent, event in enumerate(forEach(event, maxEvents=5, skipEvents=8), 8):
        ...
      
  will process events from 8 to 12 (included) at most; `iEvent` starts from `8`
  because we explicitly requested that in the `enumerate()` call, or else it
  would have started from `0` (i.e. `enumerate()` does not know about the events
  being skipped).
  """
  if fromStart:
    if skipEvents > 0: event.goToEntry(skipEvents)
    else:              event.toBegin()
  elif skipEvents > 0:
    for _ in range(skipEvents): event.next() # probably slower than goToEntry()
  nEvents = 0
  while (not event.atEnd() and ((maxEvents is None) or (nEvents < maxEvents))):
    yield event
    nEvents += 1
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
      print("Opening: '%s'" % inputFiles[iFile])
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
    print("Encountered %d/%d errors." % (nErrors, nProcessedEvents),file=sys.stderr)
  return nErrors
# eventLoop()



# this does not really work...
# ROOT.gallery.Event.__iter__ = lambda self: EventIterator(self)


################################################################################
### Infrastructure
################################################################################

def findFHiCL(configRelPath, extraDirs = []):
  
  if os.path.isfile(configRelPath):
    return os.path.join(os.getcwd(), str(configRelPath))
  for path in extraDirs + os.environ.get('FHICL_FILE_PATH', "").split(':'):
    candidate = os.path.join(path, str(configRelPath))
    if os.path.isfile(candidate): return candidate
  else: return None

# findFHiCL()


def getTableIfPresent(
 pset: "fhicl.ParameterSet object containing the desired configuration table",
 key: "name of the table to be extracted",
 defValue: "value returned if the key is not present (None by default)" = None,
 type_: "type to look for (instead of ROOT.fhicl.ParameterSet)" = None,
 ) -> "ParameterSet with `key` in `pset`, or `defValue` if not present":
  
  try: return pset.get(type_ if type_ else ROOT.fhicl.ParameterSet)(key)
  except Exception: return defValue
  
# getTableIfPresent()


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
    try: return self.config.get[klass](key)
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


def loadConfiguration(configSpec):
  # this utility actually relies on generic utilities that while not LArSoft
  # specific, are nevertheless distributed with LArSoft (`larcorealg`).
  SourceCode.loadHeaderFromUPS("larcorealg/Geometry/StandaloneBasicSetup.h")
 
  if isinstance(configSpec, ConfigurationString):
    import tempfile
    configFile = tempfile.NamedTemporaryFile("w+")
    configFile.write(str(configSpec))
    configFile.flush()
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
    return self.config.get[ROOT.fhicl.ParameterSet](FHiCLpath)
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
    print("Starting message facility for %s..." % applName)
    ROOT.mf.StartMessageFacility(config, applName)
    startMessageFacility.Init = True
  # init()
# class startMessageFacility


################################################################################
