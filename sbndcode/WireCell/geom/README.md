# Updating Wire Geometry Files

When a gdml update changes any wire-related geometry, the input Wire-Cell geoemtry file needs to be updated as well. The wire geometry file is stored in `sbnd_data/<version>/WireCell` with filename `sbnd-wires-geometry-<gdml_version>.json.bz2`. 

There are two steps to generate the new `json.bz2` file.

1. Running `CTreeGeometry` via a `lar` command to generate `ChannelWireGeometry.txt`. It'll obtain the gdml from the SBND geometry service automatically. 
```bash
lar -n 1 -c ctree_geometry.fcl
```
2. Running the `wire-cell-python` (repo [here](https://github.com/WireCell/wire-cell-python)) convert function on the aforementioned text file: 

```bash
wirecell-util convert-multitpc-wires --face-style sbnd ChannelWireGeometry.txt sbnd-wires-geometry-<gdml_version>.json.bz2
```
