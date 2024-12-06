local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local tools = tools_maker(params);

local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";

local tools = tools_maker(params);

local sim = sim_maker(params, tools);

local anode = tools.anodes[0];
local ductors = sim.make_anode_ductors(anode);

local md_chain = sim.multi_ductor_chain(ductors);
local md = sim.multi_ductor(anode, ductors, md_chain);

{
    uses: md.uses,
    edges: md.edges
}
