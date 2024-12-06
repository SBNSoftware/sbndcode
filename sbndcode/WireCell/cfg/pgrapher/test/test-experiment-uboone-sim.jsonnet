local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local sim_maker = import "pgrapher/experiment/uboone/sim.jsonnet";

local tools = tools_maker(params);

local sim = sim_maker(params, tools);

sim
