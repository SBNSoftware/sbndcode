// This produces a function to configure DNN-ROI for one APA given
// anode and torch service (ts) objects.
// 
// The prefix is prepended to all internal node names if uniqueness
// beyond anode ID is needed.  The output_scale allows for an ad-hoc
// scaling of dnnroi output.  The U and W planes will go through
// dnnroi while hte W plane will be shunted.  What comes out will be a
// unified frame with frame tag "dnnspN" where "N" is the anode ID.

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";


function (anode, ts_p0, ts_p1, prefix="dnnroi", output_scale=1.0, nticks=3428, tick_per_slice=4, nchunks=1)
    local apaid = anode.data.ident;
    local prename = prefix + std.toString(apaid);
    local intags = ['loose_lf%d'%apaid, 'mp2_roi%d'%apaid,
                     'mp3_roi%d'%apaid];

    local dnnroi_u = pg.pnode({
        type: "DNNROIFinding",
        name: prename+"u",
        data: {
            anode: wc.tn(anode),
            plane: 0,
            intags: intags,
            decon_charge_tag: "decon_charge%d" %apaid,
            outtag: "dnnsp%du"%apaid,
            input_scale: 0.00025, // 1/4000
            output_scale: output_scale,
            forward: wc.tn(ts_p0),
            tick_per_slice: tick_per_slice,
            nticks: nticks,
            nchunks: nchunks
        }
    }, nin=1, nout=1, uses=[ts_p0, anode]);
    local dnnroi_v = pg.pnode({
        type: "DNNROIFinding",
        name: prename+"v",
        data: {
            anode: wc.tn(anode),
            plane: 1,
            intags: intags,
            decon_charge_tag: "decon_charge%d" %apaid,
            outtag: "dnnsp%dv"%apaid,
            input_scale: 0.00025, // 1/4000
            output_scale: output_scale,
            forward: wc.tn(ts_p1),
            tick_per_slice: tick_per_slice,
            nticks: nticks,
            nchunks: nchunks
        }
    }, nin=1, nout=1, uses=[ts_p1, anode]);
    local dnnroi_w = pg.pnode({
        type: "PlaneSelector",
        name: prename+"w",
        data: {
            anode: wc.tn(anode),
            plane: 2,
            tags: ["gauss%d"%apaid],
            tag_rules: [{
                frame: {".*":"DNNROIFinding"},
                trace: {["gauss%d"%apaid]:"dnnsp%dw"%apaid},
            }],
        }
    }, nin=1, nout=1, uses=[anode]);

    local dnnpipes = [dnnroi_u, dnnroi_v, dnnroi_w];
    local dnnfanout = pg.pnode({
        type: "FrameFanout",
        name: prename,
        data: {
            multiplicity: 3
        }
    }, nin=1, nout=3);

    local dnnfanin = pg.pnode({
        type: "FrameFanin",
        name: prename,
        data: {
            multiplicity: 3,
            tag_rules: [{
                frame: {".*": "dnnsp%d%s" % [apaid,plane]},
                trace: {".*": "dnnsp%d%s" % [apaid,plane]},
            } for plane in ["u", "v", "w"]]
        },
    }, nin=3, nout=1);

    local retagger = pg.pnode({
      type: "Retagger",
      name: 'dnnroi%d' % apaid,
      data: {
        // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
        tag_rules: [{
          // Retagger also handles "frame" and "trace" like fanin/fanout
          // merge separately all traces like gaussN to gauss.
          frame: {
            ".*": "dnnsp%d" % apaid
          },
          merge: {
            ".*": "dnnsp%d" % apaid
          },
        }],
      },
    }, nin=1, nout=1);
    
    pg.intern(innodes=[dnnfanout],
              outnodes=[retagger],
              centernodes=dnnpipes+[dnnfanin],
              edges=[pg.edge(dnnfanout, dnnpipes[ind], ind, 0) for ind in [0,1,2]] +
              [pg.edge(dnnpipes[ind], dnnfanin, 0, ind) for ind in [0,1,2]] +
              [pg.edge(dnnfanin, retagger, 0, 0)])
