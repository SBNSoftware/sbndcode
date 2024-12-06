// This WCT configuration file defines configuration nodes which may
// be used in a WC/LS job.  They configure components defined not in
// WCT proper but in the larwirecell package of larsoft.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";


function(params, tools)
{
    // converters from data which is input to WCT
    input : {
        // Note: scale is -1 to correct a sign error in the SimDepoSource converter.
        depos : function(name="deposource", model="", scale=-1.0, art_tag="plopper:bogus", assn_art_tag="") g.pnode({
            type: 'wclsSimDepoSource',
            name: name,
            data: {
                model: model,
                scale: scale,
                art_tag: art_tag, //name of upstream art producer of depos "label:instance:processName"
                assn_art_tag: assn_art_tag,
            },
        }, nin=0, nout=1),      // fixme: should add model to uses?

        digits : function(arttag, name="adcdigits", tags=["orig"], nticks=params.nf.nsamples) g.pnode({
            type: 'wclsRawFrameSource',
            name: name,
            data: {
                // The art::Event inputTag that locates the input raw::RawDigit collection.
                art_tag: arttag,

                // These tags will be placed on the resulting frame.
                frame_tags: tags,

                // This component can pad/trunc traces as it converts
                nticks: nticks,
            },
        }, nin=0, nout=1),

    },

    // converters for data which is output from WCT
    output : {
        // Save a frame to raw::RawDigits
        digits : function(name="digitsaver", tags=["wct"], cmm=[]) g.pnode({
            type: "wclsFrameSaver",
            name: name, 
            data: {
                anode: wc.tn(tools.anode),
                digitize: true,         // true means save as RawDigit, else recob::Wire
                frame_tags: tags,
                nticks: params.daq.nticks,
                chanmaskmaps: cmm,
            },
        }, nin=1, nout=1, uses=[tools.anode]),

        // Save a frame to recob::Wires
        signals : function(tags=["wct"], name="signalsaver", cmm=[]) g.pnode({
            type: "wclsFrameSaver",
            name: name, 
            data: {
                anode: wc.tn(tools.anode),
                digitize: false,         // true means save as RawDigit, else recob::Wire
                frame_tags: tags,
                nticks: params.daq.nticks,
                chanmaskmaps: cmm,
            },
        },nin=1, nout=1, uses=[tools.anode]),

        // Save wiener RMS.  this used to be saved in threshold trace tag summary.
        thresholds: function(tags=["wiener"], name="thsaver", cmm=[]) g.pnode({
            type: "wclsFrameSaver",
            name: name,
            data: {
                anode: wc.tn(tools.anode),
                digitize: false,
                sparse: true,
                frame_tags: [],
                frame_scale: null,      // this turns off saving frames
                nticks: params.daq.nticks,
                summary_tags: tags,
                summary_scale: 1.0,
            },
        }, nin=1, nout=1, uses=[tools.anode]),
        
        // fixme: https://github.com/WireCell/larwirecell/issues/3

    },

}

