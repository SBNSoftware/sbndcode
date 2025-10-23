# Updating Wire Geometry Files

When a gdml update changes any wire-related geometry (such as absolute wire coordinates), the input Wire-Cell geometry file needs to be updated as well. The wire geometry file is stored in `sbnd_data/<version>/WireCell` with filename `sbnd-wires-geometry-<gdml_version>.json.bz2`. 

There are two steps to generate the new `json.bz2` file.

1. Running `CTreeGeometry` via a `lar` command to generate `ChannelWireGeometry.txt`. It'll obtain the gdml from the SBND geometry service automatically. 
```bash
lar -n 1 -c ctree_geometry.fcl
```
2. Running the `wire-cell-python` (repo [here](https://github.com/WireCell/wire-cell-python)) convert function on the aforementioned text file: 

```bash
wirecell-util convert-multitpc-wires --face-style sbnd ChannelWireGeometry.txt sbnd-wires-geometry-<gdml_version>.json.bz2
```
To run WireCell with this new .json.bz2 file, `sbndcode/WireCell/cfg/pgrapher/experiment/sbnd/params.jsonnet` needs the new file name under base->files->wires, and the file needs to be visible in `$WIRECELL_PATH`.

To test locally, the directory containing this file can be append to `$WIRECELL_PATH`.
For this to be picked up by an official release, this .json.bz2 file should be added to `sbnd_data`.
