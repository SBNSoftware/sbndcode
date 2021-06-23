Here is a description of the gdml folder. 
In case of emergency, break the glass: gustavo.valdiviesso@unifal-mg.edu.br 

There are two main versions:
  - sbnd_v01_0X_base.gdml (legacy version)
  - sbnd_v02_0X_base.gdml ( new LArG4 compatible version)
  
 If you HAVE to update the geometry, do it in a _base file. The final useful files are generated from this one using the preparser. Compile a local version of the preparser (located at ./preparser/preparserGDML.cpp) and use it like so:

$ preparseGDML sbnd_v02_00_base.gdml -w

This will generate the two files required by larsoft:
  - bnd_v02_00.gdml
  - bnd_v02_00_nowires.gdml

If you need to ENABLE the overburden (it is disabled by default) you need to generate the geometries with the right setup:

$ preparseGDML sbnd_v02_00_base.gdml -w --setup DefaultWithShielding

You can visualize the geometries with the geoVis_sbnd.C root macro:

root 'geoVis_sbnd.C("sbnd_v02_00.gdml")'

Pay attention to the fact that the _base files is not compatible with ROOT and cannot be visualided like this. It is, however, fully compatible with Geant4. If you venture to write a G4 visualizer, let me know.
