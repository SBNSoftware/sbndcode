// This is a main entry point for configuring a wire-cell CLI job to
// simulate SBND.  It is simplest signal-only simulation with
// one set of nominal field response function.  

local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
local savetid = std.extVar("save_track_id");
local roi = std.extVar('roi');
local nchunks = std.extVar('nchunks');
local tick_per_slice = std.extVar('tick_per_slice');
local dnnroi_model_p0 = std.extVar('dnnroi_model_p0');
local dnnroi_model_p1 = std.extVar('dnnroi_model_p1');

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
local wc_device = std.extVar('wc_device');
local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

local data_params = import 'params.jsonnet';
local simu_params = import 'simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;

local base = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
local params = base {
  lar: super.lar { // <- super.lar overrides default values
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.s,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.s,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.ms,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
  sim: super.sim {
    // front porch size [us]
    tick0_time: std.extVar('tick0_time') * wc.us,
  }
};

local default_tools = tools_maker(params);
local tools = if wc_device == 'gpu' then std.mergePatch(default_tools,
    {dft: {type: "TorchDFT", data: {device: wc_device}}}) else default_tools;
local sim_maker = import 'pgrapher/experiment/sbnd/sim.jsonnet';
local sim = sim_maker(params, tools);
local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);

local wcls_maker = import "pgrapher/ui/wcls/nodes.jsonnet";
local wcls = wcls_maker(params, tools);

// added Ewerton 2023-03-14
local wcls_input_sim = {
    depos: wcls.input.depos(name="", art_tag=std.extVar('inputTag')),
    deposet: g.pnode({
            type: 'wclsSimDepoSetSource',
            name: "",
            data: {
                model: "",
                scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
                art_tag: std.extVar('inputTag'), //name of upstream art producer of depos "label:instance:processName"
                id_is_track: if (savetid == 'true') then false else true,
                assn_art_tag: "",
            },
        }, nin=0, nout=1),
};
// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};
local wcls_output_sim = {
  // ADC output from simulation
  // sim_digits: wcls.output.digits(name="simdigits", tags=["orig"]),
  sim_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'simdigits',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
     frame_tags: ['daq'],
      nticks: params.daq.nticks,
      pedestal_mean: 'native',
    },
  }, nin=1, nout=1, uses=[mega_anode]),
 // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: wcls.output.digits(name="nfdigits", tags=["raw"]),
  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: wcls.output.signals(name="spsignals", tags=["gauss", "wiener"]),
  // save "threshold" from normal decon for each channel noise
  // used in imaging
  sp_thresholds: wcls.output.thresholds(name="spthresholds", tags=["wiener"]),
};


local drifter = sim.drifter;

local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);

// signal plus noise pipelines
local sn_pipes = sim.splusn_pipelines;

local rng = tools.random;
local wcls_depoflux_writer = g.pnode({
  type: 'wclsDepoFluxWriter',
  name: 'postdrift',
  data: {
    anodes: [wc.tn(anode) for anode in tools.anodes],
    field_response: wc.tn(tools.field),
    tick: 0.5 * wc.us,
    window_start: 0.0 * wc.ms,
    window_duration: self.tick * params.daq.nticks,
    nsigma: 3.0,

    reference_time: -1700 * wc.us,

    //energy: 1, # equivalent to use_energy = true
    simchan_label: 'simpleSC',
    sed_label: if (savetid == 'true') then 'ionandscint' else '',
    sparse: false,
  },
}, nin=1, nout=1, uses=tools.anodes + [tools.field]);

local sp_maker = import 'pgrapher/experiment/sbnd/sp.jsonnet';

local sp_override = 
if roi == "dnn" then {
    sparse: true,
    use_roi_debug_mode: true,
    save_negative_charge: false, // TODO: no negative charge in gauss, default is false
    use_multi_plane_protection: true,
    do_not_mp_protect_traditional: false, // TODO: do_not_mp_protect_traditional to make a clear ref, defualt is false 
    mp_tick_resolution:4,
    tight_lf_tag: "",
    cleanup_roi_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
    //decon_charge_tag: "",
    gauss_tag: "",
    wiener_tag: "",
} 
else if roi == "both" then {
    sparse: true,
    use_roi_debug_mode: true,
    save_negative_charge: false, // TODO: no negative charge in gauss, default is false
    use_multi_plane_protection: true,
    do_not_mp_protect_traditional: false, // TODO: do_not_mp_protect_traditional to make a clear ref, defualt is false 
    mp_tick_resolution:4,
    tight_lf_tag: "",
    cleanup_roi_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
} 
else if roi == "trad" then {
    sparse: true,
};

local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local magoutput = 'sbnd-data-check.root';
local magnify = import 'pgrapher/experiment/sbnd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local base = import 'pgrapher/experiment/sbnd/chndb-base.jsonnet';
//local perfect = import 'pgrapher/experiment/sbnd/chndb-perfect.jsonnet';

local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  //data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/sbnd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local multipass1 = [
  g.pipeline([
               sn_pipes[n],
             ],
             'multipass%d' % n)
  for n in anode_iota
];

local multipass2 = [
  g.pipeline([
               sn_pipes[n],
               //sinks.orig_pipe[n],
               
               nf_pipes[n], // NEED to include this pipe for channelmaskmaps    
               //sinks.raw_pipe[n], 
               
               sp_pipes[n], 

               //sinks.decon_pipe[n],
               //sinks.threshold_pipe[n],
             ],
             'multipass%d' % n)
  for n in anode_iota
];


local f_sp = import 'pgrapher/experiment/sbnd/funcs.jsonnet';

local outtags = ['orig%d' % n for n in anode_iota];
local bi_manifold1 = f.fanpipe('DepoSetFanout', multipass1, 'FrameFanin', 'sn_mag_nf', outtags);
local bi_manifold2 = f_sp.fanpipe('DepoSetFanout', multipass2, 'FrameFanin', 'sn_mag_nf_mod2', outtags, "true");

local retagger_sim = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'orig',
      },
      merge: {
        'orig\\d': 'daq',
      },
    }],
  },
}, nin=1, nout=1);

local sink_sim = sim.frame_sink;

//===============================NF+SP============================================

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.

local wcls_output_sp = {
  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),


  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener'],

      // this may be needed to convert the decon charge [units:e-] to be consistent with the LArSoft default
      // for SBND, this scale is about ~50. Frame scale needed when using LArSoft producers reading in recob::Wire.
      frame_scale: [0.02, 0.02],
      nticks: params.daq.nticks,

      // uncomment the below configs to save summaries and cmm
       summary_tags: ['wiener'],
       summary_operator: {wiener: 'set'},
       summary_scale: [0.02], # summary scale should be the same as frame_scale
       chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),
 dnnsp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'dnnsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['dnnsp'],

      // this may be needed to convert the decon charge [units:e-] to be consistent with the LArSoft default ?unit? e.g. decon charge * 0.005 --> "charge value" to GaussHitFinder
      frame_scale: [0.02],
      nticks: params.daq.nticks,
      chanmaskmaps: [],
    },
  }, nin=1, nout=1, uses=[mega_anode]),

};
local dnnroi = import 'dnnroi.jsonnet';

local ts_p0 = {
    type: "TorchService",
    name: "dnnroi_p0",
    tick_per_slice: tick_per_slice, 
    data: {
        model: dnnroi_model_p0,
        device: wc_device,
        concurrency: 1,
    },
};

local ts_p1 = {
    type: "TorchService",
    name: "dnnroi_p1",
    tick_per_slice: tick_per_slice, 
    data: {
        model: dnnroi_model_p1,
        device: wc_device,
        concurrency: 1,
    },
};

local dnnroi_pipes = [ dnnroi(tools.anodes[n], ts_p0, ts_p1, output_scale=1, nchunks=nchunks) for n in std.range(0, std.length(tools.anodes) - 1) ];

local fanout = function (name, multiplicity=2)
  g.pnode({
    type: 'FrameFanout',
    name: name,
    data: {
        multiplicity: multiplicity
    },
  }, nin=1, nout=multiplicity);

local sp_fans = [fanout("sp_fan_%d" % n) for n in std.range(0, std.length(tools.anodes) - 1)];



local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(5638 * n, 5638 * (n + 1) - 1),
    },
  }, nin=1, nout=1)
  for n in anode_iota
];


local nfsp_pipes = 
if roi == "dnn" then
// oports: 0: dnnroi, 1: traditional sp
[
  g.intern(
    innodes=[chsel_pipes[n]],
    outnodes=[dnnroi_pipes[n],sp_fans[n]],
    centernodes=[nf_pipes[n], sp_pipes[n], sp_fans[n]],
    edges=[
      g.edge(chsel_pipes[n], nf_pipes[n], 0, 0),
      g.edge(nf_pipes[n], sp_pipes[n], 0, 0),
      g.edge(sp_pipes[n], sp_fans[n], 0, 0),
      g.edge(sp_fans[n], dnnroi_pipes[n], 0, 0),
    ],
    iports=chsel_pipes[n].iports,
    oports=dnnroi_pipes[n].oports+[sp_fans[n].oports[1]],
    name='nfsp_pipe_%d' % n,
  )
  for n in std.range(0, std.length(tools.anodes) - 1)
]
else if roi == "both" then
// oports: 0: dnnroi, 1: traditional sp
[
  g.intern(
    innodes=[chsel_pipes[n]],
    outnodes=[dnnroi_pipes[n],sp_fans[n]],
    centernodes=[nf_pipes[n], sp_pipes[n], sp_fans[n]],
    edges=[
      g.edge(chsel_pipes[n], nf_pipes[n], 0, 0),
      g.edge(nf_pipes[n], sp_pipes[n], 0, 0),
      g.edge(sp_pipes[n], sp_fans[n], 0, 0),
      g.edge(sp_fans[n], dnnroi_pipes[n], 0, 0),
    ],
    iports=chsel_pipes[n].iports,
    oports=dnnroi_pipes[n].oports+[sp_fans[n].oports[1]],
    name='nfsp_pipe_%d' % n,
  )
  for n in std.range(0, std.length(tools.anodes) - 1)
]
else if roi == "trad" then
[
  g.pipeline([
               chsel_pipes[n],
               //sinks.orig_pipe[n],

               nf_pipes[n], // NEED to include this pipe for channelmaskmaps 
               //sinks.raw_pipe[n],

               sp_pipes[n],
               //sinks.decon_pipe[n],
               //sinks.threshold_pipe[n],
               // sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
             ],
             'nfsp_pipe_%d' % n)
  for n in anode_iota
];

//local fanpipe = f_sp.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf'); // commented Ewerton 2023-05-24
local fanpipe = f_sp.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf_mod');   //added Ewerton 2023-05-24

local retagger = function(name) g.pnode({
  type: 'Retagger',
  name: name,
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
        'dnnsp\\d': 'dnnsp',
      },
    }],
  },
}, nin=1, nout=1);

local retagger_dnnroi = retagger("retagger_dnnroi");
local retagger_sp = retagger("retagger_sp");

local sink = function(name) g.pnode({ type: 'DumpFrames', name: name }, nin=1, nout=0);
local sink_dnnroi = sink("sink_dnnroi");
local sink_sp = sink("sink_sp");

local fanout_apa = g.pnode({
    type: 'FrameFanout',
    name: 'fanout_apa',
    data: {
        multiplicity: std.length(tools.anodes),
        "tag_rules": [
            {
               "frame": {
                  ".*": "orig%d" % n
               },
               "trace": { }
            }
            for n in std.range(0, std.length(tools.anodes) - 1)
        ]
        }},
    nin=1, nout=std.length(tools.anodes));

local framefanin = function(name) g.pnode({
    type: 'FrameFanin', 
    name: name,
    data: {
        multiplicity: std.length(tools.anodes),

         "tag_rules": [
            {     
                "frame": {
                  ".*": "framefanin"
                },
                trace: {
                  ['dnnsp%d' % n]: ['dnnsp%d' % n],
                  ['gauss%d' % n]: ['gauss%d' % n],
                  ['wiener%d' % n]: ['wiener%d' % n],
                  ['threshold%d' % n]: ['threshold%d' % n],
                },
            }
            for n in std.range(0, std.length(tools.anodes) - 1)
         ],    
         "tags": [ ]
    },
}, nin=std.length(tools.anodes), nout=1);
local fanin_apa_dnnroi = framefanin('fanin_apa_dnnroi');
local fanin_apa_sp = framefanin('fanin_apa_sp');

local graph1_trad = g.pipeline([
wcls_input_sim.deposet,         //sim
setdrifter,                     //sim
wcls_depoflux_writer,           //sim
bi_manifold1,                   //sim
retagger_sim,                   //sim
wcls_output_sim.sim_digits,     //sim
fanpipe,                        //sp
retagger_sp,                    //sp
wcls_output_sp.sp_signals,      //sp
sink_sp                         //sp
]);


local graph2_trad = g.pipeline([
wcls_input_sim.deposet,         //sim
setdrifter,                     //sim
wcls_depoflux_writer,           //sim
bi_manifold2,                   //sim
retagger_sp,                    //sp
wcls_output_sp.sp_signals,      //sp
sink_sp                         //sp
]);

local graph_sim = g.pipeline([
wcls_input_sim.deposet,         //sim
setdrifter,                     //sim
wcls_depoflux_writer,           //sim
bi_manifold1,                   //sim
]);

local graph_sp_dnn = 
g.intern(
  innodes=[retagger_sim],
  outnodes=[],
  centernodes=nfsp_pipes+[wcls_output_sim.sim_digits, fanout_apa, retagger_dnnroi, retagger_sp, fanin_apa_dnnroi, fanin_apa_sp, wcls_output_sp.dnnsp_signals, sink_dnnroi, sink_sp],
  edges=[
    g.edge(retagger_sim, wcls_output_sim.sim_digits, 0, 0),
    g.edge(wcls_output_sim.sim_digits, fanout_apa, 0, 0),
    g.edge(fanout_apa, nfsp_pipes[0], 0, 0),
    g.edge(fanout_apa, nfsp_pipes[1], 1, 0),
    g.edge(nfsp_pipes[0], fanin_apa_dnnroi, 0, 0),
    g.edge(nfsp_pipes[1], fanin_apa_dnnroi, 0, 1),
    g.edge(fanin_apa_dnnroi, retagger_dnnroi, 0, 0),
    g.edge(retagger_dnnroi, wcls_output_sp.dnnsp_signals, 0, 0),
    g.edge(wcls_output_sp.dnnsp_signals, sink_dnnroi, 0, 0),
    g.edge(nfsp_pipes[0], fanin_apa_sp, 1, 0),
    g.edge(nfsp_pipes[1], fanin_apa_sp, 1, 1),
    g.edge(fanin_apa_sp, retagger_sp, 0, 0),
    g.edge(retagger_sp, sink_sp, 0, 0),
  ]
);

local graph_sp_both = 
g.intern(
  innodes=[retagger_sim],
  outnodes=[],
  centernodes=nfsp_pipes+[wcls_output_sim.sim_digits, fanout_apa, retagger_dnnroi, retagger_sp, fanin_apa_dnnroi, fanin_apa_sp, wcls_output_sp.sp_signals, wcls_output_sp.dnnsp_signals, sink_dnnroi, sink_sp],
  edges=[
    g.edge(retagger_sim, wcls_output_sim.sim_digits, 0, 0),
    g.edge(wcls_output_sim.sim_digits, fanout_apa, 0, 0),
    g.edge(fanout_apa, nfsp_pipes[0], 0, 0),
    g.edge(fanout_apa, nfsp_pipes[1], 1, 0),
    g.edge(nfsp_pipes[0], fanin_apa_dnnroi, 0, 0),
    g.edge(nfsp_pipes[1], fanin_apa_dnnroi, 0, 1),
    g.edge(fanin_apa_dnnroi, retagger_dnnroi, 0, 0),
    g.edge(retagger_dnnroi, wcls_output_sp.dnnsp_signals, 0, 0),
    g.edge(wcls_output_sp.dnnsp_signals, sink_dnnroi, 0, 0),
    g.edge(nfsp_pipes[0], fanin_apa_sp, 1, 0),
    g.edge(nfsp_pipes[1], fanin_apa_sp, 1, 1),
    g.edge(fanin_apa_sp, retagger_sp, 0, 0),
    g.edge(retagger_sp, wcls_output_sp.sp_signals, 0, 0),
    g.edge(wcls_output_sp.sp_signals, sink_sp, 0, 0),
  ]
);

local graph_dnn = g.pipeline([
graph_sim,
graph_sp_dnn,
]);

local graph_both = g.pipeline([
graph_sim,
graph_sp_both,
]);

local save_simdigits = std.extVar('save_simdigits');

local graph = 
if roi == "dnn" then graph_dnn 
else if roi == "both" then graph_both 
else if save_simdigits == "true" then graph1_trad 
else graph2_trad;

local app = {
  type: 'TbbFlow',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]

