// This file provides a function which takes a params object (see
// ../params/) and returns a data structure with a number of
// sub-objects that may configure various WCT "tool" type componets
// which are not INodes.

// Some attributes are merely default and you may wish to override
// them.  For example, the default IDFT FftwDFT and to instead ues
// TorchDFT you may do something like:
//
// local default_tools = tools_maker(params)
// local tools = std.mergePatch(default_tools,
//   {dft: {type: "TorchDFT", data: {device: "gpu"}}});
// 

local wc = import "wirecell.jsonnet";

function(params)
{
    // The IRandom pRNG
    random : {
        type: "Random",
        data: {
            generator: "default",
            seeds: [0,1,2,3,4],
        }
    },
    // The IDFT FFT implementation 
    dft : {
        type: "FftwDFT",
    },

    // One FR per field file.
    fields : std.mapWithIndex(function (n, fname) {
        type: "FieldResponse",
        name: "field%d"%n,
        data: { filename: fname }
    }, params.files.fields),

    field: $.fields[0],         // the nominal field

    // The number of ticks for response waveforms.  For simulation
    // this must be enlarged beyond what is wanted in the output
    // readout in order to accept catch early activity into the first
    // few readout ticks.  When responses are used to deconvolve then
    // this should best be shortened.  Note: sigproc at time of
    // writing did not use components for Elec or RC so are not
    // currently subject to this config.
    local sim_response_binning = {
        tick: params.daq.tick,
        nticks: params.sim.ductor.nticks, // MUST match ductor
    },

    perchanresp : {
        type: "ParamsPerChannelResponse",
        data: {
            filename: params.files.chresp,
        }
    },
    // It's a little awkward to NOT use PerChannelResponse because you
    // need to give an empty string to a component that wants to use
    // it an an empty list to the "uses".  Here we package that little
    // "if" branch.  It's a general pattern not specific to this
    // object.  Might want to make wc.tn() convert "null" into "" and
    // g.uses() convert [null] into [].
    perchanresp_nameuses : if std.type(params.files.chresp) == 'null'
    then {name:"", uses:[]}
    else {name:wc.tn($.perchanresp), uses:[$.perchanresp]},


    wires : {
        type: "WireSchemaFile",
        data: { filename: params.files.wires }
    },

    elec_resp : {
        type: if std.objectHas(params.elec, "type")
              then params.elec.type
              else "ColdElecResponse", // default
        data: sim_response_binning {
            shaping: params.elec.shaping,
            gain: params.elec.gain,
            postgain: params.elec.postgain,
            filename: if std.objectHas(params.elec, "filename")
                      then params.elec.filename
                      else ""
        },
    },

    rc_resp : {
        type: "RCResponse",
        data: sim_response_binning {
            // width: 1.0*wc.ms,
            width: if std.objectHas(params, 'rc_resp')
            then params.rc_resp.width else 1.0*wc.ms,
        }
    },


    sys_resp : {
        type: "ResponseSys",
        data: sim_response_binning {
            start: params.sys_resp.start,
            magnitude: params.sys_resp.magnitude,
            time_smear: params.sys_resp.time_smear,
        }
    },

    // there is one trio of PIRs (one per wire plane in a face) for
    // each field response.  WARNING/fixme: this sets the default DFT
    // with no way to override!  This config structure needs a redo!
    pirs : std.mapWithIndex(function (n, fr) [
        {
            type: "PlaneImpactResponse",
            name : "PIR%splane%d" % [fr.name, plane],
            data : sim_response_binning {
                plane: plane,
                dft: wc.tn($.dft),
                field_response: wc.tn(fr),
                // note twice we give rc so we have rc^2 in the final convolution
                short_responses: if params.sys_status == false
                                    then [wc.tn($.elec_resp)]
                                    else [wc.tn($.elec_resp), wc.tn($.sys_resp)],
		overall_short_padding: if std.objectHas(params, 'overall_short_padding')
                                    then params.overall_short_padding
                                    else if params.sys_status == false
                                    then 0.1*wc.ms
                                    // cover the full time range of the convolved short responses
                                    else 0.1*wc.ms - params.sys_resp.start,
		// long_responses: [wc.tn($.rc_resp), wc.tn($.rc_resp)],
        long_responses: if std.objectHas(params, 'rc_resp')
        then std.makeArray(params.rc_resp.rc_layers, function(x) wc.tn($.rc_resp))
        else [wc.tn($.rc_resp), wc.tn($.rc_resp)],
		long_padding: 1.5*wc.ms,
	    },
            uses: [$.dft, fr, $.elec_resp, $.rc_resp, $.sys_resp],
        } for plane in [0,1,2]], $.fields),

    // One anode per detector "volume"
    anodes : [{
        type : "AnodePlane",
        name : vol.name,
        data : {
            // This ID is used to pick out which wires in the
            // WireSchema belong to this anode plane.
            ident : vol.wires,
            nimpacts: params.sim.nimpacts,
            // The wire schema file
            wire_schema: wc.tn($.wires),

            faces : vol.faces,
        },
        uses: [$.wires],
    } for vol in params.det.volumes],

    // Arbitrarily call out the first anode to make single-anode
    // detector config slightly cleaner.
    anode: $.anodes[0],

}
