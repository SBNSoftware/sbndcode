// This defines nodes for Ar39 depo source

local wc = import "wirecell.jsonnet";
local ar39spectrum = import "ar39.jsonnet"

function(params, tools)
{
    local bb = params.det.bounds,
    local vol = v.volume(v.frompoint(bb.tail), v.frompoint(bb.head)),

    ar39fulldet : g.pnode({
        type: "BlipSource",
        name: "ar39fulldet",
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

}
