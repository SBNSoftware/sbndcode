% SBND and LArSoft utilities for interface with Python
% Gianluca Petrillo (petrillo@slac.stanford.edu)
% June 23, 2019

This document describes the utilities provided in `sbndcode` to facilitate the
use of LArSoft and SBND-specific code developed for _art_/_gallery_ into
Python.

This kind of documents tend to become outdated very quickly, so the reader is
invited, when in doubt, to refer to the code itself. Most of the utilities that
are designed to be part of the interface to users have some terse documentation
in the form of python documentation strings.


Format of this document
------------------------

This file (`README.md`) is written in a dialect of [Markdown] format that is
compatible with [`pandoc`][pandoc]. The format is designed to br readable as
plain text: in that case, the reader will find all the links to external
resources at the end of the document.
The program `pandoc` allows to render this file in other formats. For example,
to render it into Portable Document Format:
    
    pandoc --toc -o README.pdf README.md
    
and to render it into HyperText Transfer Protocol format:
    
    pandoc --toc -o README.html README.md
    


Introduction to the interface to Python
========================================

LArSoft and SBND code, as well as _gallery_ code, rely on ROOT to interface
with Python. Despite some relevant limitations, and despite what we are all used
to think about ROOT, its interface with Python is able to expose into Python a
great deal of the C++ code, including some more abstruse constructs that have
no equivalent in Python.

ROOT can import C++ objects, functions, namespaces and variables from C++
compiled code, and it places it under `ROOT` python module. Therefore, the
object `geo::Point_t` of `larcoreobj` is importer as `ROOT.geo.Point_t` and
the whole namespace `geo` is available as `ROOT.geo`. Note that for template
classes ROOT invents a special pattern, where the template class name, e.g.
`ROOT.lar.sparse_vector`, corresponds to `lar::sparse_vector` (which is almost
never used in C++ code), while the actual class name, e.g.
`lar::sparse_vector<int>`, is accessed as `ROOT.lar.sparse_vector("int")`.

ROOT needs information from two sources: the interface (e.g. the definition of
`geo::GeometryCore`) is learned from a C++ header file, while the actual code
(e.g. what to do when `geo::GeometryCore::IteratePlanes` is called) is learned
from a compiled library. The fact that the building system of _art_/_gallery_
(including `cet_build_tools`, UPS, MRB) follows some conventions helps to
automate the discovery of all the needed components, but this is not fool-proof.


> Currently this code is written in Python 2 (2.7.14 or newer expected).
> An update to Python 3 will be attempted in due time, when LArSoft switch to
> it. In the examples, it is assumed that the command `python` accesses the
> correct version. If this is not the case on a certain platform, the use of
> executables like `python2` or `python2.7` might be necessary.



Components of the Python support libraries
===========================================

The support libraries provided here are organized according to the software they
interact with:

* `ROOTutils.py` simplifies some common tasks in handling ROOT objects, like
  analyzing a ROOT object path and printing standard ROOT vectors
* `cppUtils.py` facilitates the management of C++ source code interacting with
  ROOT
* `galleryUtils.py` simplifies some common tasks like creating a data product
  handle, loading a FHiCL configuration file, keep track of service providers,
  and more
* `LArSoftUtils.py` simplifies the creation of LArSoft service providers, and
  some special support for the Geometry service provider (`geo::GeometryCore`)
* `SBNDutils.py` includes functions to load SBND-specific service providers
  (like SBND geometry implementation)
* `SBNDservices.py` facilitates the management of service providers in the
  context of an SBND execution environment, with fitting defaults.
  It can also used as an executable to start an interactive section where some
  initialization has already automatically taken place.

In the following sections, details are provided on the utilities from these
libraries that are expected to be most commonly used, starting from the highest
level ones.



`SBNDservices.py`: service providers (and more) in SBND context
--------------------------------------------------------------------

`SBNDservices.py` may be used as library in a Python script, or directly to
launch an interactive session.
It provides easier access to LArSoft service providers in Python.


### Interactive Python session configured for SBND

The interactive section can be started by just executing
    
    python SBNDservices.py
    
Curious to see something about the geometry of the wire planes? or detector
timing settings? the standard methods of service providers can be easily
accessed:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
geom = ServiceManager('Geometry')
for plane in geom.IteratePlanes(): print plane.PlaneInfo()

detClocks = ServiceManager('DetectorClocks')
detClocks.debugReport()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supported service providers include: `Geometry`, `DetectorClocks`,
`LArProperties`, `DetectorProperties`.


### Library utilities

`ServiceManager`
  
:   object to load and keep track of services.
    It is configured by FHiCL and its configuration must include all the
    services that are going to be requested, in a way similar to a regular
    LArSoft job configuration (e.g. in a `services` table).
  
    `ServiceManager.setConfiguration(configPath, serviceTable)`
    :   selects the configuration file (`configPath`) and the FHiCL table
        (`"services"` by default) where the configuration of the services is
        stored; all the following requests for new services will use this
        configuration
    
    An example of using `ServiceManager`:
    
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
    from SBNDservices import ServiceManager
    ServiceManager.setConfiguration \
      ('services_sbnd_simulation.fcl', 'icarus_simulation_services')
    detClocks = ServiceManager('DetectorClocks')
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    uses one of the most complete standard service configuration tables of
    SBND, `icarus_simulation_services`, from `services_sbnd_simulation.fcl`
    configuration file. That file is found in the usual way via
    `FHICL_FILE_PATH`, which should have been set already during the environment 
    setup.

`geometry()`

:   shortcut for `ServiceManager('Geometry')`: loads and returns an instance of
    LArSoft geometry service provider (`geo::GeometryCore`) with SBND
    configuration as found in the `Geometry` entry of the service configuration
    table described above.
  


`SBNDutils.py`: lower level SBND-specific functions
--------------------------------------------------------

`SBNDutils.py` provides SBND-specific functions:

* `Geometry` service provider for SBND:

    `loadSBNDgeometry`

    :   loads and returns SBND geometry with the standard SBND channel
        mapping. It is usually easier to go through `ServiceManager` in `SBNDservices.py`, which uses this facility anyway.

    `justLoadSBNDgeometry`

    :   Loads and returns SBND geometry from the specified configuration file.
        This is a one-stop procedure recommended only when running interactively
        and needing just the geometry service provider. Otherwise, it is better
        to have a centralized system to manage the configuration and the
        services: see `SBNDservices.ServiceManager`.
  


`LArSoftUtils.py`: one-stop shop for LArSoft/_gallery_ scripts
---------------------------------------------------------------

The module `LArSoftUtils.py` provides most of the facilities that are typically
needed for a standalone LArSoft script or a _gallery_ script.
Most of the utilities are actually exposed directly from `galleryUtils.py`.

* source code management: exposes `galleryUtils.SourceCode` and
  `galleryUtils.readHeader`
* configuration helpers: exposes `galleryUtils.findFHiCL`,
  `galleryUtils.ConfigurationClass`, `galleryUtils.ConfigurationHelper` and
  `galleryUtils.loadConfiguration`
* data product support: exposes `galleryUtils.make_getValidHandle`
* _art_ event processing support
  exposes `galleryUtils.makeFileList`, `galleryUtils.forEach` and
  `galleryUtils.eventLoop`
* service providers:
  exposes `galleryUtils.ServiceRegistryClass` and
  `galleryUtils.startMessageFacility`. It also provides:
    
  `loadSimpleService`
  
  :   a simple service is a service whose provider that can be initialized
      in a standardised way. LArSoft already uses some protocol to initialize
      services designed for simple unit tests not pulling _art_ in.
      If a service supports this type of setup (by specialization of
      `testing::ProviderSetupClass` template), then 'loadSimpleService' should
      be able to initialize it.
  
  `loadGeometry`, `justLoadGeometry`
  
  :   `Geometry` service is not simple, and its initialization involves
      another "helper" service. Because of this, a special initialization
      function is provided for the `Geometry` service provider. The two
      versions are designed one in the context of an environment where
      multiple services are used, and the other in an environment where just
      the geometry is needed. It is recommended that the former is preferred.

  

`galleryUtils.py`: helpers for data product and event access
-------------------------------------------------------------

The module `galleryUtils.py` providers helper functions to access _art_ data
products, manage FHiCL configuration and service providers. It also exposes
utilities for source code management:

* source code management: exposes `cppUtils.readHeader` and
  `cppUtils.SourceCode`
* data product access:
  
  `make_getValidHandle`
  
  : getting a data product handle via Python requires some obscure preliminary
    steps, involving template member functions. This functions takes the name
    of the data type in the data product and performs those steps. After this,
    it will be possible to call `gallery.Event.getValidHandle` to read the data
    products from the events. Example:
    `galleryUtils.make_getValidHandle('std::vector<raw::RawDigit>')`.
  
* _gallery_ event loop:

  `makeFileList`
  
  :   `gallery.Event` constructor requires a file list as argument, in the form
      of `std::vector<std::string>`; this function creates such an object
      starting from a sequence of Python strings.
    
  
  `forEach`
  
  :   an event loop helper for _gallery_, allowing the use of a `gallery.Event`
      in a `for` loop. The loop may look like:
      
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
      for iEvent, event in enumerate(galleryUtils.forEach(event)):
        ...
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  `eventLoop`
  
  :   a complete event loop, initialized with input file list, a function to
      call on each event, and some options (number of events to skip and to
      process). Useful if you are in a hurry.
      
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
      import galleryUtils
      from ROOTutils import ROOT
      
      def countChannels(iEvent, event):
        rawDigits = event.getValidHandle(ROOT.vector(ROOT.raw.RawDigit))('daq')
        print "Event %d has %d channels" % (iEvent, len(rawDigits))
      # countChannels()
      
      galleryUtils.make_getValidHandle('std::vector<raw::RawDigit>')
      galleryUtils.eventLoop('detsim.root', countChannels)
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
* configuration management:
  
  `findFHiCL`
  
  :   returns the full path of a single FHiCL file. The algorithms used is not
      (necessarily) the same as in _art_: the file is searched in the current
      directory and then in all directories specified in `FHICL_FILE_PATH`
      environment variable.
    
  `loadConfiguration`
  
  :   parses a FHiCL configuration file and returns a `fhicl::ParameterSet`
      object.

  `ConfigurationHelper`
  
  :   a simple object wrapping a `fhicl::ParameterSet` object and mimicking its
      `get` and `has` methods. The `get()` function (as well as a direct call to
      the object itself) allows to specify either a class type or a default
      value (in which case the class type is deduced from it), or both.
      Example:
      
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
      pset = galleryUtils.loadConfiguration('parameters.fcl')
      config = galleryUtils.ConfigurationHelper(pset)
      threshold = config('threshold', 0.0)
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  `ConfigurationClass`
  
  :   an object that reads a job configuration FHiCL file and allows for quick
      access to service and module configuration. Example:
      
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
      config = galleryUtils.ConfigurationClass('standard_g4_icarus.fcl')
      LArG4Parameters = config.service('LArG4Parameters')
      LArG4 = config.producer('largeant')
      Output = config.paramsFor('outputs.out1')
      ProcessName = config.config.get('std::string')('process_name')
      
      print "Process name: '%s'" % ProcessName
      print "LArG4 configuration:"
      print LArG4.to_indented_string()
      print "LArG4 parameters:"
      print LArG4Parameters.to_indented_string()
      print "Output:"
      print Output.to_indented_string()
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
* service providers:

  `ServiceRegistryClass`
  
  :   an object keeping track of both configuration and service providers.
      Services are referred to by the name of their configuration
      (e.g. `Geometry`, `LArPropertiesService`).
      It is initialized with a FHiCL configuration file, and provides:
      
      `config`
      
      :   configuration of the specified service
      
      `register`
      
      :   registers (or replaces) an object in the service registry
      
      `has`
      
      :   whether the specified service is registered
      
      `registeredServiceNames`
      
      :   list of names of all registered services
      
      `get`, or direct call
      
      :   returns the object registered with the specified name
      
      `create`
      
      :   constructs and registers an object. The configuration is the one
          returned by `config()`.
          Example:
          
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
          from galleryUtils import SourceCode, ServiceRegistryClass
          from ROOTutils import ROOT
          
          SourceCode.loadHeaderFromUPS('lardataalg/DetectorInfo/LArPropertiesStandard.h')
          SourceCode.loadLibrary('lardataalg_DetectorInfo')
          
          ServiceRegistry = ServiceRegistryClass('standard_g4_icarus.fcl')
          larProp = ServiceRegistry.create \
            ('LArPropertiesService', ROOT.detinfo.LArPropertiesStandard)
          print larProp.RadiationLength()
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          Note that this is equivalent to using `SBNDservices.py`:
          
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
          from SBNDservices import ServiceManager
          ServiceManager.setConfiguration('standard_g4_icarus.fcl')
          larProp = ServiceManager('LArProperties')
          print larProp.RadiationLength()
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          which also ensures that services required by `LArProperties` are
          loaded.
      
  
  `startMessageFacility`
  
  :   configures the `messages` service, typically used by LArSoft algorithms
      and service providers, as well as by _gallery_, to emit messages on
      screen


> Note: for `galleryUtils.py` to be fully functional, `gallery` UPS product must
> be set up. That is not part of the usual MRB set up.



`cppUtils.py`: helpers for source management
---------------------------------------------

The `cppUtils.py` module provides a few utilities to manage the loading of C++
software into Python via ROOT.
This is normally achieved by sending the ROOT interpreter some commands that it
turns to its C++ compiler (based on Clang). ROOT can then bind the compiled
objects, types, functions etc. to the `ROOT` Python namespace.
These utilities simplify the procedure relying on the assumption that the code
is following some recommended naming practices and that it is set up in a way
compatible with Fermilab UPS.

> Although the examples in this section use a LArSoft data product
> (`raw::RawDigit`), this machinery is usually _not necessary for data products_
> because ROOT typically knows how to find those. For example, all is needed in
> a properly configured environment is: `import ROOT; d = ROOT.raw.RawDigit` to
> get a `raw::RawDigit` object. These utilities are instead useful when looking
> for algorithms and service providers.


`SourceCode`

:   this object manages the loading of C++ code into ROOT. Only the functions
    expected to be commonly used are documented here.
    
    `addIncPath`, `addIncPaths`
    
    :   adds the specified path(s) to the list of directories ROOT looks for C++
        headers into; it understands environment variables (e.g.
        `${LARDATAOBJ_INC}/RawData`; see `os.path.expandvars`). It won't add a
        path if it's already present, unless `force`d to.
        

    `addIncPathEnv`, `addIncPathEnvs`
    
    :   adds the path in the specified environment variable(s) to the list of
        directories ROOT looks for C++ headers into; example:
        `addIncPathEnv('LARDATAOBJ_INC')`.

    `loadHeader`, `loadHeaderFromUPS()` 
    
    :   instructs ROOT to read and compile the specified C++ header.
        `loadHeader` looks for the header, which is specified as a relative
         path, in the paths registered with `addIncPath` family of functions.
        In addition, `loadHeaderFromUPS()` also tries to look in the UPS product
        the header belongs to, whose name is deduced from the first element of
        the header path.
        Both return the full path of the specified header.
        Example:
        `SourceCode.loadHeader('lardataobj/RawData/RawDigit.h')` or,
        equivalent,
        `SourceCode.loadHeaderFromUPS('lardataobj/RawData/RawDigit.h')`.
        For the latter, there is no need to explain `SourceCode` where to find
        the header, since it will automatically consider `$LARDATAOBJ_INC`.
    
    `loadLibrary`
    
    :   instructs ROOT to load the specified compiled library into memory.
        Example: `SourceCode.loadLibrary('lardataobj_RawData')` (library name
        is expanded adding `lib` prefix and a proper suffix).
    
`readHeader`

:   makes ROOT compile the specified header; this is a necessary step to be able
    to use C++ code from precompiled libraries. Most often, using this function
    directly is not necessary if the `SourceCode` object is used to manage the
    code.




`ROOTutils.py`: common ROOT operations
---------------------------------------------

The module `ROOTutils.py` includes helpers for some often-needed tasks related
to ROOT. It also is a "replacement" for `ROOT` module: ROOT can be imported by:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
from ROOTutils import ROOT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The only reason to do that instead of plain `import ROOT` is that `ROOTutils.py`
attempts to load ROOT module in a way that does not interfere with command line
parsing (the most common symptom is ROOT intercepting command line options like
`--help` or `-b`).

Finally, the module also attaches string conversion functions to some commonly
used ROOT objects, so that they can be printed more easily:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
v = ROOT.TVector3(1.5, 0.0, -2.0)
print v
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

will show `Name: TVector3 Title: A 3D physics vector` with plain `ROOT` module,
and `( 1.5, 0, -2 )` if `ROOTutils` is loaded.

Other provided utilities are:

`splitROOTpath`

:   returns the specified path split into file path and ROOT directory path

`createROOTpath`

:   creates a complete ROOT directory path, including a new file if needed,
    and all `TDirectoryFile` objects in the path.

`getROOTclass`

:   returns the Python class for the specified class: for example,
    `ROOTutils.getROOTclass('geo::Point_t')` returns `ROOT.geo.Point_t`.
    The advantage of this function is that it works even if intermediate objects
    haven't been loaded yet. ROOT binds names to Python in a lazy way, only when
    they are explicitly requested. For example, immediately after loading
    `GeometryCore.h`, namespace `geo` is not bound to `ROOT.geo`. Asking members
    of namespace `geo` ends up being trouble in this situation.

`activateDirectory`

:   switches to a ROOT directory, and back when done. This is used as a context
    manager:
    
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.python}
    dir = outputFile.GetDirectory("plots")
    with activateDirectory(dir):
      for plot in plots: item.Write()
    # with
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    will switch to `plots` ROOT directory, write all plots into it and then,
    at the end of the `with` block, switch back to whatever ROOT directory was
    active before `activateDirectory` call. If an exception is thrown in that
    block, the active directory is restored as well.
    Note that the object acting as context manager,
    `ROOTutils.DirectoryChanger`, is also directly available if needed.



Execution environment
======================

This facility includes interactions between many diverse components: _gallery_,
ROOT, LArSoft, SBND code, Python... While it has been proven that this _can_
work, the setup of a working environment might be non-trivial.
The following recipes has proven to work on SBND SLF6 GPVM servers:

* use a binary distribution, setting it up with `setup icaruscode ...`
* use a MRB area with all of the following:
    * `icaruscode` checked out in the source directory, or explicitly set up;
      e.g. `( cd "$SOURCE_DIR" && mrb gitCheckout icaruscode ; )`
    * all checked out code properly compiled _and installed_: `mrb install -j4`
    * all local "products" set up: `mrbslp`

Both setup procedures should set up all is needed, include important parts like:

* UPS environment variables, like `$LARCOREALG_INC`, used to locate header files
  and compiled libraries
* ROOT machinery (`ROOTSYS` and more)
* `PYTHONPATH` environment variable, used to locate the Python modules described
  in this document

In addition, to use _gallery_, it must be explicitly set up. The easiest way to
get there is to set up `larsoftobj`:

* check which `larsoftobj` version is matching the LArSoft version in MRB
  working area in the [LArSoft release web page][LArRel]
* after the setup procedure described above, set up `larsoftobj` as well


[Markdown]: https://daringfireball.net/projects/markdown
[pandoc]: https://pandoc.org/MANUAL.html
[LArRel]: https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/LArSoft_release_list
