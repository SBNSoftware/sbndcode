#!/usr/bin/env python

__doc__ = """
Provides a service manager preconfigured with SBND service providers.

A `ServiceManager` object bridges the access to the service providers.
In the most straightforward cases, getting a service provider is as simple as
calling `ServiceManager(providerName)` with the name of the service as string
(e.g. `larProp = ServiceManager('LArProperties')`).



If this module is executed as a script, the service manager is loaded and a
python interactive session is started. The main code also show how to override
the service manager setup by choosing a different service configuration.
"""

__all__ = [ 'ServiceManager', 'geometry', ]


import SBNDutils  # loadSBNDgeometry()
import LArSoftUtils


################################################################################
### special known services
###

class SBNDGeometryServiceGetter(LArSoftUtils.SimpleServiceLoader):
  
  def __init__(self):
    LArSoftUtils.SimpleServiceLoader.__init__(self, 'Geometry')
  
  def load(self, manager):
    return SBNDutils.loadSBNDgeometry(registry=manager.registry())
  
# class SBNDGeometryServiceGetter


################################################################################
###  module setup - part I
###  
###  Where global service manager is set up.
###  

class SBNDserviceManagerClass(LArSoftUtils.ServiceManagerInstance):
  
  DefaultConfigPath = "simulationservices_sbnd.fcl"
  DefaultServiceTable = "sbnd_simulation_services"
  
  def defaultConfiguration(self):
    """
    
    Configuration:
    
    If `serviceTable` is not `None`, a new configuration is created with the
    service table as `serviceTable`, and `configPath` is included in that
    configuration (presumably to define `serviceTable`). In this case, if
    `configPath` is `None`, "simulationservices_sbnd.fcl" is included.
    
    If both `configPath` and `serviceTable` are `None`, the configuration is
    created as above, using "sbnd_simulation_services" as `serviceTable`.
    
    Finally, if only `serviceTable` is `None`, the configuration file in
    `configPath` is included directly, and it is assumed that it already
    properly defines a `services` table.
    """
    return DefaultConfigPath, DefaultServiceTable
  # defaultConfiguration()
  
  def __init__(self):
    LArSoftUtils.ServiceManagerInstance.__init__(self)
    self.setConfiguration(
      configFile=SBNDserviceManagerClass.DefaultConfigPath,
      serviceTable=SBNDserviceManagerClass.DefaultServiceTable,
      )
  # __init__()
  
  def setup(self):
    """Prepares for SBND service provider access in python/Gallery."""
    
    LArSoftUtils.ServiceManagerInstance.setup(self)

    #
    # register the services we know about;
    # some are already known
    # (`LArSoftUtils.ServiceManagerClass.StandardLoadingTable`), including
    # 'Geometry', 'LArProperties', 'DetectorClocks' and 'DetectorProperties',
    # but se override the former with our flavor of it
    #
    
    self.manager.registerLoader('Geometry', SBNDGeometryServiceGetter())
    
    return self.manager
    
  # setup()

# class SBNDserviceManagerClass


ServiceManager = SBNDserviceManagerClass()


################################################################################

def geometry(): return ServiceManager.get('Geometry')


################################################################################

if __name__ == "__main__":
  
  #
  # unfortunately, ROOT module interferes with command line arguments.
  #
  import argparse
  
  Parser = argparse.ArgumentParser(
    description=
      "Starts a python interactive session with `ServiceManager` available."
    )
  
  Parser.add_argument("--config", "-c", dest="configPath",
    help="configuration file path (must define `services` or host serviceTable below)"
    )
  Parser.add_argument("--servicetable", "-T", dest="serviceTable",
    help="name of the FHiCL table where all services are configured")
  
  args = Parser.parse_args()
  
  if args.configPath is not None:
    ServiceManager.setConfiguration(args.configPath, args.serviceTable)
  
  # we want ROOT module known in the interactive session;
  # and we keep the good habit of loading it via ROOTutils
  from ROOTutils import ROOT
  
  try:
    import IPython
    IPython.embed()
  except ImportError:
    import code
    code.interact(local=dict(globals(), **locals()))
  # try ... except
  
# main
