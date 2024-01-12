// This defines some depo source nodes.  They make certain limiting
// choices such as volume 

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local v = import "vector.jsonnet";
local ar39spectrum = import "ar39.jsonnet";

function(params, tools)
{
    // Return a pnode that makes ar39 thoughout the given bounding box.
    ar39 :: function(bb = params.det.bounds, name="ar39depos") g.pnode({
        local vol = v.volume(v.frompoint(bb.tail), v.frompoint(bb.head)),
        type: "BlipSource",
        name: name,
        data: {
            rng: wc.tn(tools.random),
            charge: ar39spectrum,
            time: {
                type: "decay",
                start: params.daq.start_time,
                stop: params.daq.stop_time,
                activity: params.lar.ar39activity * params.lar.density * vol,
            },
            position: {
                type: "box",
                extent: bb,
            },
        },
    }, nin=0, nout=1, uses=[tools.random]),

    // Return a node that produces ideal track depos given track segments.
    tracks :: function(tracklist, name="trackdepos", step=1.0*wc.mm) g.pnode({
        type: "TrackDepos",
        name: name,
        data: {
            step_size: step,
            tracks: tracklist,
        }
    }, nin=0, nout=1),

    jsondepos :: function(name="jsondepos", file="noinputfile") g.pnode({
        type: 'JsonDepoSource',
        name: name,
        data : {
            filename: file,
            model: "electrons",  // take "n" from depo as already in number of electrons
            scale: 1.0,           // multiply by "n", there is a hard-coded -1.0 in cxx file
        }
    }, nin=0, nout=1),
}
