import fhicl

pset = fhicl.make_pset('prodoverlay_corsika_cosmics_proton_genie_rockbox_sbnd.fcl')

#
# Extract flux and GDML used to make the GeomScan file
#
pset_geom_scan = pset['physics']['producers']['generator']['GeomScan']

geom_scan_flux_config = pset_geom_scan.split('flux')[1][0]
geom_scan_gdml_version = pset_geom_scan.split('gdml')[1][:6]

print('Using GeomScan file:', pset_geom_scan)
print('\t -> using flux configuration:', geom_scan_flux_config)
print('\t -> using gdml version:', geom_scan_gdml_version)


#
# Extract flux and GDML set in services and GENIE
#
pset_gdml = pset['services']['Geometry']['GDML']
pset_flux = pset['physics']['producers']['generator']['FluxSearchPaths']

flux_config = pset_flux.split('config')[1][0]
gdml_version = pset_gdml.split('sbnd_')[1][:6]

print('Using:')
print('\t -> gdml version:', gdml_version)
print('\t -> flux configuration:', flux_config)

#
# Check that the two are the same, if not, need to regenerate the GeomScan file
#
message = 'Please regenerate the GeomScan file. Instructions: link.'

assert geom_scan_flux_config == flux_config, f"Have you updated/changed the flux files? {message}"

assert geom_scan_gdml_version == gdml_version, f"Have you updated/changed the GDML file? {message}"


