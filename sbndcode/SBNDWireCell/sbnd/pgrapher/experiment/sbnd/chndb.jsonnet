// Some helpers for making channel noise "db" objects

// local perfect = import 'chndb-perfect.jsonnet';
local base = import 'chndb-base.jsonnet';

function(params, tools) {

    perfect(anode) :: {
        type:'OmniChannelNoiseDB',
        name: 'ocndbperfect-' + anode.name,
        data: base(params, anode, tools.field, anode.data.ident),
        uses: [anode, tools.field],
    },

}
