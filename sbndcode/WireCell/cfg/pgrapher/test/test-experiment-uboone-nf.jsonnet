local params = import "pgrapher/experiment/uboone/params.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local nf_maker = import "pgrapher/experiment/uboone/nf.jsonnet";
local chndb_maker = import "pgrapher/experiment/uboone/chndb.jsonnet";

local tools = tools_maker(params);
local anode = tools.anodes[0];

local chndbs = chndb_maker(params, tools.anodes[0], tools.field);
local nf = nf_maker(params, anode, chndbs.wct("before"));

nf
