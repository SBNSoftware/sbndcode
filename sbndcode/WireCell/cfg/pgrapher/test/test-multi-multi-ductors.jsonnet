// Test using either monolithic MultiDuctor or graph-based equivalent.
//
// This requires the '-V' option to either wire-cell or jsonnet in
// order to set the "multiductor" variable to either "component" or
// "graph".
// 
// Output is to a .npz file of the same name as this file.  Plots can
// be made to do some basic checks with "wirecell-gen plot-sim".

local multiductor = std.extVar("multiductor");

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";

local io = import "pgrapher/common/fileio.jsonnet";
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";

local tools = tools_maker(params);

local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";

local nf_maker = import "pgrapher/experiment/uboone/nf.jsonnet";
local chndb_maker = import "pgrapher/experiment/uboone/chndb.jsonnet";

local sp_maker = import "pgrapher/experiment/uboone/sp.jsonnet";

local stubby = {
    tail: wc.point(1000.0, 0.0, 5000.0, wc.mm),
    head: wc.point(1100.0, 0.0, 5100.0, wc.mm),
};

local tracklist = [
    {
        time: 1*wc.ms,
        charge: -5000,          // negative means per step
        ray: stubby,
        //ray: params.det.bounds,
    },
];
local output = "wct-sim-ideal-sn-nf-sp.npz";
    
local anode = tools.anodes[0];

local sim = sim_maker(params, tools);

//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
local depos = sim.tracks(tracklist);


local deposio = io.numpy.depos(output);

local drifter = sim.drifter;

local ductors = sim.make_anode_ductors(anode);
local md_chain = sim.multi_ductor_chain(ductors);
local md_component = sim.multi_ductor(anode, ductors, [md_chain]);
local md_pipes = sim.multi_ductor_pipes(ductors);
local md_graph = sim.multi_ductor_graph(anode, md_pipes, "mdg");
local ductor = if multiductor == "component" then md_component else md_graph;
    

local noise_model = sim.make_noise_model(anode, sim.miscfg_csdb);
local noise = sim.noise(anode, noise_model).return;

local digitizer = sim.digitizer(anode, tag="orig");
local sim_frameio = io.numpy.frames(output, "simframeio", tags="orig");

local sink = sim.frame_sink;

local graph = g.pipeline([depos, deposio, drifter, ductor, noise, digitizer,
                          sim_frameio, sink]);

local app = {
    type: "Pgrapher",
    data: {
        edges: g.eduges(graph),
    },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellImg"],
        apps: ["Pgrapher"]
    },
};


// Finally, the configuration sequence which is emitted.
[cmdline] + g.uses(graph) + [app]
