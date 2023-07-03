// This file provides configuration pnodes for various WCT components
// that perform file I/O.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";

{
    // Numpy file format (.npz) saves Numpy arrays which are useful
    // for quick processing in Python.  NOTE: if the file already
    // exists new data will be APPENDED.  User shoud delete prior file
    // if that is not desired.
    numpy: {

        // Inject this between two ports passing a frame to save it in a file.
        frame :: function(fname, tags=[], name="") g.pnode({
            type: "NumpyFrameSaver",
            name: name,
            data: {
                filename: fname,
                frame_tags: tags,
            }
        }, nin=1, nout=1),
                                                 
        // Inject this between two ports depos.
        depos :: function(fname, name="") g.pnode({
            type: "NumpyDepoSaver",
            name: name,
            data: {
                filename: fname
            }
        }, nin=1, nout=1),
                                                 
}
