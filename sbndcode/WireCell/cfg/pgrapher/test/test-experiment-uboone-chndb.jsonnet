local g = import "pgraph.jsonnet";
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local tools = tools_maker(params);
local chndb_maker = import "pgrapher/experiment/uboone/chndb.jsonnet";

local chndbs = chndb_maker(params, tools);

local temp = {
    wct: {
        before: chndbs.wct("before"),
        after: chndbs.wct("after"),
    },
    wcls: {
        before: chndbs.wcls("before"),
        after: chndbs.wcls("after"),
    },
    multi: chndbs.wcls_multi("multichndb"),
};


local one = chndbs.wct("before");

//std.objectHas(one, "uses")
//std.foldl(g.popuses, one.uses, [])

// std.mapWithKey(function (k,v) if k == "uses" then null else v, one)
// util.prune(std.mapWithKey(function (k,v) if k == "uses" then null else v, one))
// std.prune({a:42, uses:null, b: [[1,2,3],[4,5,6]]})

//g.popuses([],one);
//g.resolve_uses([chndbs.wct("before")])

temp


