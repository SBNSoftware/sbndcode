// A slowdown when number of depos gets high is found in 


local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";

local io = import "pgrapher/common/fileio.jsonnet";

local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
// customized override -- moved to simparams.jsonnet
/* local params = default_params{ */
/*     files: super.files{ */
/*         chresp: null, */
/*     } */
/* }; */

local tools_maker = import "pgrapher/common/tools.jsonnet";

local tools = tools_maker(params);

local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";


local stubby = {
    tail: wc.point(1000, -1000, 5000.0, wc.mm),
    head: wc.point(1500, -1000, 6000.0, wc.mm),
};

local tracklist = [
    {
        time: 0.0*wc.ms,
        charge: -5000,          // negative means per step
        ray: stubby,
        //ray: params.det.bounds,
    },
    // {
    //     time: 2.0*wc.ms,
    //     charge: -5000,          // negative means per step
    //     ray: stubby,
    // },
];

local anode = tools.anodes[0];

local sim = sim_maker(params, tools);

local depos = sim.tracks(tracklist, step=0.1*wc.mm);

local drifter = sim.drifter;
local sink = g.pnode({type:'DumpDepos'}, nin=1, nout=0);

local graph = g.pipeline([depos, drifter, sink]);

local app = {
    type: "Pgrapher",
    data: {
        edges: g.edges(graph),
    },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSigProc", "WireCellImg"],
        apps: ["Pgrapher"]
    },
};


[cmdline] + g.uses(graph) + [app]
