local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";

local tools = tools_maker(params);
local sim = sim_maker(params, tools);

local anode = tools.anodes[0];
local ductors = sim.make_anode_ductors(anode);
local md_pipes = sim.multi_ductor_pipes(ductors);
local ductor = sim.multi_ductor_graph(anode, md_pipes, "mdg");

//tools.fields
//tools.pirs
//ductors[0]
ductor
