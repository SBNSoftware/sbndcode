#!/usr/bin/env python

from __future__ import print_function

__doc__ = """
Collection of utilities to interface LArSoft with python and gallery.

This module requires ROOT.
"""

__all__ = [
  'readHeader',
  'SourceCode',
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
  'loadGeometry',
  'justLoadGeometry',
  'loadSimpleService',
  ]

import sys, os
import ROOTutils
from ROOTutils import ROOT
import galleryUtils
from galleryUtils import (
  make_getValidHandle, makeFileList, forEach, eventLoop,
  findFHiCL, loadConfiguration, ConfigurationHelper, ConfigurationClass,
  startMessageFacility, ServiceRegistryClass,
  )


################################################################################
### Pass-through
###
SourceCode = galleryUtils.SourceCode
readHeader = galleryUtils.readHeader


################################################################################
### LArSoft
################################################################################
def loadGeometry(config=None, registry=None, mapping=None):
  """The argument `config` is an instance of `ConfigurationClass`.

  If a config object is provided, configurations will be read from there.
  Otherwise, they will be read from the registry.
  If a registry is provided, the services will be registered in there.
  """
  assert(config or registry)
  serviceName = 'Geometry'

  # deal with messagefacility first
  if not (registry and registry.has("message")):
    messageConfig = config.service("message") if config else registry.config("message")
    startMessageFacility(messageConfig) # use default name
    if registry: registry.register("message", None) # there is no direct access, sorry
  # if need to load message facility

  geometryConfig = config.service(serviceName) if config else registry.config(serviceName)
  if geometryConfig is None:
    raise RuntimeError("Failed to retrieve the configuration for %s service" % serviceName)

  if not mapping:
    SourceCode.loadHeaderFromUPS('larcorealg/Geometry/ChannelMapStandardAlg.h')
    mapping = ROOT.geo.ChannelMapStandardAlg
  SourceCode.loadHeaderFromUPS("larcorealg/Geometry/StandaloneGeometrySetup.h")
  SourceCode.loadLibrary("larcorealg_Geometry")
  service = ROOT.lar.standalone.SetupGeometry(mapping)(geometryConfig)
  if registry: registry.register(serviceName, service)

  # make it easy to print points and vectors in python
  for varName in ( 'Point_t', 'Vector_t', ):
    try: klass = getattr(ROOT.geo, varName)
    except AttributeError: continue
    klass.__str__ = ROOTutils.TVector3ToString
  # ... and IDs...
  for varName in ( 'CryostatID', 'TPCID', 'PlaneID', 'WireID', ):
    try: klass = getattr(ROOT.geo, varName)
    except AttributeError: continue
    klass.__str__ = klass.toString
  # for ID
  for varName in ( 'CryostatID', 'TPCsetID', 'ROPID', ):
    try: klass = getattr(ROOT.readout, varName)
    except AttributeError: continue
    klass.__str__ = klass.toString
  # for ID
  # ... and geometry objects
  for varName in ( 'CryostatGeo', 'TPCGeo', 'PlaneGeo', 'WireGeo', 'OpDetGeo', 'AuxDetGeo', ):
    try: klass = getattr(ROOT.geo, varName)
    except AttributeError: continue
    klass.__str__ = getattr(klass, varName[:-3] + "Info")
  # for geo object

  return service
# loadGeometry()


################################################################################
def justLoadGeometry(configFile, mapping=None):
  """Loads and returns the geometry from the specified configuration file.

  This is a one-stop procedure recommended only when running interactively.
  """
  return loadGeometry(config=ConfigurationClass(configFile), mapping=mapping)
# justLoadGeometry()


################################################################################
def loadSimpleService \
  (serviceClass, config=None, registry=None, interfaceClass=None, args = []):
  """Loads a service assuming some simple requirements:
   * no dependency from other services
   * constructor accepts a FHiCL parameter set

  If a config object is provided, configurations will be read from there.
  Otherwise, they will be read from the registry.
  If a registry is provided, the services will be registered in there.

  The service configuration is read from an item called as `interfaceClass`,
  or `serviceClass` itself if `interfaceClass` is None, with "Service" appended.

  As a first attempt, the test helpers are attempted to load the service, and
  the arguments are used to construct a `providers_type` object.
  If there is no test helper for `serviceClass`, then direct construction is
  attempted using all the specified arguments.
  """

  serviceName = (interfaceClass if interfaceClass else serviceClass).__name__
  configKey = serviceName + "Service"

  if not isinstance(config, ROOT.fhicl.ParameterSet):
    config = config.service(configKey) if config else registry.config(configKey)
  if config is None:
    raise RuntimeError("Failed to retrieve the configuration for %s service" % serviceName)

  TestSetupClass = ROOT.testing.ProviderSetupClass(serviceClass)
  if TestSetupClass:
    serviceSetup = TestSetupClass.setup
    try: providers = [ serviceClass.providers_type(*args), ]
    except:
      # we received arguments, called expected them to be used, we did not:
      if args: raise # there must have been some problem
      providers = []
    service = serviceSetup(config, *providers)
  else:
    service = serviceClass(config, *args)

  if registry: registry.register(serviceName, service)
  return service
# loadSimpleService()


################################################################################
### special known services
###
###
###

def makeStringList(l): return [ l, ] if isinstance(l, str) else l


class SimpleServiceLoader:
  """Class storing the parameters needed to load a "simple" service.

  So far a "simple" service is one that can be loaded with `loadSimpleService()`
  allowing service provider dependencies and simple configuration adaptions.
  """
  def __init__(self,
    serviceClass,
    interfaceClass = None,
    headers = [],
    libraries = [],
    dependencies = [],
    purgeConfig = [],
    addConfig = {}
    ):
    assert serviceClass is not None
    self.serviceClass = serviceClass
    self.interfaceClass = interfaceClass
    self.headers = makeStringList(headers)
    self.libraries = makeStringList(libraries)
    self.serviceDependencies = makeStringList(dependencies)
    self.purgeConfig = makeStringList(purgeConfig)
    self.addConfig = addConfig
  # __init__()

  def _needsSpecialConfig(self): return self.addConfig or self.purgeConfig

  def _makeConfig(self, registry):
    for configKey in ( self.serviceKey(), self.serviceKey() + "Service", ):
      try:
        config = registry.config(configKey)
        break
      except Exception: pass
    else: config = None
    if not config:
      raise RuntimeError \
        ("No configuration for service '{}'".format(self.serviceKey()))
    # if
    for key in self.purgeConfig: config.erase(key)
    for key, value in self.addConfig.items(): config.put(key, str(value))
    return config
  # _makeConfig()

  def serviceKey(self):
    source = (self.interfaceClass if self.interfaceClass else self.serviceClass)
    if not isinstance(source, str): source = source.__name__
    return source.replace('::', '.').split('.')[-1]
  # serviceKey()

  def load(self, manager):
    # load the required service dependencies
    dependencies = [
      manager(dependency) for dependency in self.serviceDependencies
      ]

    # load the required headers and libraries
    for header in self.headers:
      galleryUtils.SourceCode.loadHeaderFromUPS(header)
    for library in self.libraries:
      galleryUtils.SourceCode.loadLibrary(library)

    # loads the actual classes from ROOT
    self.expandClass('serviceClass')
    if self.interfaceClass is not None: self.expandClass('interfaceClass')

    # if we don't need a special configuration,
    # we let loadSimpleService() find it
    registry = manager.registry()
    config = self._makeConfig(registry) if self._needsSpecialConfig() else None

    return loadSimpleService(
      self.serviceClass, config=config, registry=registry,
      interfaceClass=self.interfaceClass,
      args=dependencies,
      )
  # get()

  def __call__(self, manager): return self.load(manager)

  def expandClass(self, attrName):
    classOrName = getattr(self, attrName)
    if isinstance(classOrName, str):
      classOrName = ROOTutils.getROOTclass(classOrName)
      setattr(self, attrName, classOrName)
    return classOrName
  # expandClass()


# class SimpleServiceLoader


class GeometryServiceGetter:

  def serviceKey(self): return 'Geometry'

  def load(self, manager): return loadGeometry(registry=manager.registry())

  def get(self, manager):
    try: return manager.registry().get(self.serviceKey())
    except KeyError:
      return self.load(manager)
  # get()

  def __call__(self, manager): return self.get(manager)

# class GeometryServiceGetter



class ServiceManagerInterface:
  """The interface of a service manager implementation."""

  def registry(self):
    """Returns the service registry."""
    return NotImplementedError
  # registry()

  def loaders(self):
    """Returns the service registry."""
    return NotImplementedError
  # loaders()


  def loaded(self):
    """Returns a list names of service already loaded."""
    return self.registry().registeredServiceNames()
  # loaded()


  def supported(self):
    """Returns a list names of services known to be supported."""
    return self.loaders().keys()
  # supported()


  def registerLoader(self, serviceKey, loader):
    """Registers a service provider loader, that can then be invoked to create
    the service.
    """
    self.loaders()[serviceKey] = loader
    return self
  # registerLoader()


  def get(self, serviceName, interfaceClass = None):
    """Return (and load first when needed) the specified service.

    The service can be specified by name or by class.
    In the former case, if the service is not already configured, an exception
    will be raised.
    """
    return NotImplementedError
  # get()

  def __call__(self, serviceName, interfaceClass = None):
    return self.get(serviceName, interfaceClass=interfaceClass)


# class ServiceManagerInterface


class ServiceManagerClass(ServiceManagerInterface):

  def __init__(self, config, loadingTable = {}, preload = []):

    #
    # prepare the service registry
    #
    self.registry_ = galleryUtils.ServiceRegistryClass(config)

    #
    # register the standard services
    #
    self.loadingTable = loadingTable.copy()

    #
    # message facility
    #
    galleryUtils.startMessageFacility(self.registry().config("message"))
    self.registry().register("message", None) # there is no direct access, sorry

    #
    # preload services
    #
    for serviceKey in preload: self.get(serviceKey)

  # __init__()

  def registry(self): return self.registry_

  def loaders(self): return self.loadingTable


  def get(self, serviceKey, interfaceClass = None):
    """Return the specified service.

    The service can be specified by name or by class.
    In the former case, if the service is not already configured, an exception
    will be raised.
    If the service is specified by class instead, an attempt is made to load it,
    in which case `interfaceClass` is also passed to
    `LArSoftUtils.loadSimpleService()`.

    """

    # if it is already cached, let's use it
    try: return self.registry().get(serviceKey)
    except KeyError: pass

    # first try to see if there is a loader available for this service
    try:
      loader = self.loadingTable[serviceKey]
    except KeyError:
      loader = SimpleServiceLoader(serviceKey, interfaceClass=interfaceClass)

    print("Loading service provider: '{}'".format(serviceKey))
    return loader(self)

  # get()


# class ServiceManagerClass



################################################################################
class ServiceManagerInstance(ServiceManagerInterface):
  """A service manager with a lazy setup.

  """

  StandardLoadingTable = {

    'Geometry': GeometryServiceGetter(),

    'LArProperties': SimpleServiceLoader(
      "detinfo::LArPropertiesStandard",
      interfaceClass="detinfo::LArProperties",
      headers = 'lardataalg/DetectorInfo/LArPropertiesStandardTestHelpers.h',
      libraries = 'lardataalg_DetectorInfo',
      ),

    'DetectorClocks': SimpleServiceLoader(
      "detinfo::DetectorClocksStandard",
      interfaceClass="detinfo::DetectorClocks",
      headers = 'lardataalg/DetectorInfo/DetectorClocksStandardTestHelpers.h',
      libraries = 'lardataalg_DetectorInfo',
      ),

    'DetectorProperties': SimpleServiceLoader(
      "detinfo::DetectorPropertiesStandard",
      interfaceClass="detinfo::DetectorProperties",
      headers = 'lardataalg/DetectorInfo/DetectorPropertiesStandardTestHelpers.h',
      libraries = 'lardataalg_DetectorInfo',
      dependencies = [ 'Geometry', 'LArProperties' ],
      ),

  } # StandardLoadingTable


  class ConfigurationInfo:
    __slots__ = [ 'configPath', 'serviceTable' ]
    def __init__(self, configPath = None, serviceTable = None):
      self.configPath = configPath
      self.serviceTable = serviceTable
    def fullConfig(self): return self.serviceTable is None
    def isValid(self): return self.configPath is not None

    def __str__(self):
      if self.fullConfig(): return self.configPath
      else: return "{{services={}}}@{}".format(self.serviceTable, self.configPath)
    # __str__()

  # class ConfigurationInfo


  def __init__(self):
    self.manager = None
    self.configuration = None

  def registry(self):
    """Returns the service registry."""
    return self.manager.registry()
  # registry()

  def registerLoader(self, serviceKey, loader):
    """Registers a service provider loader, that can then be invoked to create
    the service.

    It returns the service manager.
    """
    return self.manager.registerLoader(serviceKey, loader)
  # registerLoader()


  def get(self, serviceKey, interfaceClass = None):
    """Return (and load first when needed) the specified service.

    The service can be specified by name or by class.
    In the former case, if the service is not already configured, an exception
    will be raised.
    """
    if not self.manager: self.setup()
    return self.manager(serviceKey)
  # get()


  def defaultConfiguration(self): return None

  def setConfiguration(self, configFile, serviceTable = None):
    """Sets which configuration to use for setup.

    If `serviceTable` is not `None`, a new configuration is created with the
    service table as `serviceTable`, and `configPath` is included in that
    configuration (presumably to define `serviceTable`).

    If `serviceTable` is `None` instead, the configuration file in
    `configPath` is included directly, and it is assumed that it already
    properly defines a `services` table.
    """
    assert configFile is not None
    self.configuration \
      = ServiceManagerInstance.ConfigurationInfo(configFile, serviceTable)
  # setConfiguration()

  def setup(self):
    """Prepares for service provider access in python/Gallery."""

    if self.manager is not None: return self.manager

    #
    # configuration "file"
    #
    configurationInfo = self.defaultConfiguration() \
      if self.configuration is None else self.configuration

    # if assertion fails, then `setConfiguration()` was not correctly called.
    assert self.configuration.isValid()

    if configurationInfo.fullConfig():
      config = configurationInfo.configPath
    else:
      config = galleryUtils.ConfigurationString(
          '#include "{configPath}"'
        '\n'
        '\nservices: @local::{serviceTable}'
        '\n'
        .format(
          configPath=configurationInfo.configPath,
          serviceTable=configurationInfo.serviceTable
          )
        )
    # if ... else

    #
    # prepare the service registry and manager
    #
    self.manager = ServiceManagerClass \
      (config, loadingTable=ServiceManagerInstance.StandardLoadingTable)

    #
    # register the services we know about;
    # some are already known
    # (`LArSoftUtils.ServiceManagerClass.StandardLoadingTable`), including
    # 'Geometry', 'LArProperties', 'DetectorClocks' and 'DetectorProperties',
    # but se override the former with our flavor of it
    #

    return self.manager

  # setup()

# class ServiceManagerInstance


################################################################################