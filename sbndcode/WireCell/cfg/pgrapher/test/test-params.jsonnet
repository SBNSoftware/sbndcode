local params = import "../common/params.jsonnet";

local want = std.split("lar,det,daq,adc,elec,sim,nf,files",",");
local has = [std.objectHas(params, f) for f in want];

[
    std.assertEqual(std.count(has, false),0),
    std.assertEqual(std.length(std.objectFields(params)), std.length(want)),
]
