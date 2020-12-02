// Functions originally provided by this file are moved into
// pgraph.jsonnet under the "fan." object.  Please avoid importing
// this file for new any configuration.  It is kept only to preserve
// backward compatibility for older files.
local g = import "pgraph.jsonnet";
{
    fanpipe :: g.fan.pipe,
    fansink :: g.fan.sink,
}
