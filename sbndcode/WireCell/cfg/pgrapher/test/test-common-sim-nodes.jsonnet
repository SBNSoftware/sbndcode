

local params = import "../common/params.jsonnet";
local tools_maker = import "../common/tools.jsonnet";
local tools = tools_maker(params);
local nodes_maker = import "../common/sim/nodes.jsonnet";
local nodes = nodes_maker(params, tools);
nodes
