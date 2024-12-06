local params = import "../experiment/uboone/simparams.jsonnet";
local tools_maker = import "../common/tools.jsonnet";
local tools = tools_maker(params);
local sim_maker = import "../experiment/uboone/sim.jsonnet";
local sim = sim_maker(params, tools);

local wires_modes = [
    { name:"shorteduv",
      wires: sim.shorted_channels.uv,
      mode: "accept", },
    { name:"shortedvy",
      wires: sim.shorted_channels.vy,
      mode: "accept", },
    { name:"nominal",
      wires: sim.shorted_channels.uv + sim.shorted_channels.vy,
      mode: "reject", } ];

local wbdepos = [sim.make_wbdepo(tools.anode, wm.wires, wm.mode, name=wm.name) for wm in wires_modes];
local ductors = sim.make_anode_ductors(tools.anode);
local nductors = std.length(ductors);
local sanity = std.assertEqual(nductors, 3);
local tags = ['frame'+std.toString(n) for n in std.range(1,nductors)];
local pipes = [{
    wbdepos: wbdepos[n],
    ductor: ductors[n],
    tag: tags[n],
} for n in std.range(0, nductors-1)];
//pipes
sim.multi_ductor_graph(tools.anode, pipes).graph

