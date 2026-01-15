// This is a main entry point for configuring a wire-cell CLI job to
// simulate SBND.  It is simplest signal-only simulation with
// one set of nominal field response function.  

local epoch = std.extVar('epoch');  // eg "dynamic", "after", "before", "perfect"
local reality = std.extVar('reality');
// local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
local sigoutform = "sparse";  // eg "sparse" or "dense"
local savetid = std.extVar("save_track_id");

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
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
};

local tools = tools_maker(params);
local sim_maker = import 'pgrapher/experiment/sbnd/sim.jsonnet';
local sim = sim_maker(params, tools);
local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);

local wcls_maker = import "pgrapher/ui/wcls/nodes.jsonnet";
local wcls = wcls_maker(params, tools);

// added Ewerton 2023-03-14
local wcls_input_sim = {
    depos: wcls.input.depos(name="", art_tag=std.extVar('inputTag')), //commented Ewerton 2023-03-15
    deposet: g.pnode({
            type: 'wclsSimDepoSetSource',
            name: "",
            data: {
                model: "",
                scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
                // art_tag: std.extVar('inputTag'), //name of upstream art producer of depos "label:instance:processName"
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
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
     frame_tags: ['daq'],
      // nticks: params.daq.nticks,
      // chanmaskmaps: ['bad'],
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
  sp_thresholds: wcls.output.thresholds(name="spthresholds", tags=["threshold"]),
};


local drifter = sim.drifter;

// added Ewerton 2023-03-14
local setdrifter = g.pnode({
            type: 'DepoDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);

// signal plus noise pipelines
// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local hio_sp = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_sp%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['loose_lf%d' % n
        , 'tight_lf%d' % n
        , 'cleanup_roi%d' % n
        , 'break_roi_1st%d' % n
        , 'break_roi_2nd%d' % n
        , 'shrink_roi%d' % n
        , 'extend_roi%d' % n
        , 'mp3_roi%d' % n
        , 'mp2_roi%d' % n
        , 'decon_charge%d' % n
        , 'gauss%d' % n],
        filename: "g4-rec-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        tick0: 0,
        nticks: 3427,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];


local rng = tools.random;
local wcls_deposetsimchannel_sink = g.pnode({
  type: 'wclsSimChannelSink',
  name: 'postdrift',
  data: {
    artlabel: 'simpleSC',  // where to save in art::Event
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
    rng: wc.tn(rng),
    tick: 0.5 * wc.us,
    start_time: -0.205 * wc.ms,
    readout_time: self.tick * 3427,
    nsigma: 3.0,
    drift_speed: params.lar.drift_speed,
    u_to_rp: 100 * wc.mm,  // time to collection plane
    v_to_rp: 100 * wc.mm,  // time to collection plane
    y_to_rp: 100 * wc.mm,
    u_time_offset: 0.0 * wc.us,
    v_time_offset: 0.0 * wc.us,
    y_time_offset: 0.0 * wc.us,
    g4_ref_time: -1700 * wc.us,
    use_energy: true,
  },
}, nin=1, nout=1, uses=tools.anodes);

local hio_orig = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_orig%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['orig%d'%n],
        filename: "g4-rec-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local sp_maker = import 'pgrapher/experiment/sbnd/sp.jsonnet';
local sp = sp_maker(params, tools, { sparse: sigoutform == 'sparse' });
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local magoutput = 'sbnd-data-check.root';
local magnify = import 'pgrapher/experiment/sbnd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local perfect = import 'pgrapher/experiment/sbnd/chndb-perfect.jsonnet';
//local base = import 'chndb-base_sbnd.jsonnet';

local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  // data: base(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/sbnd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local depo_fanout = [g.pnode({
    type:'DepoFanout',
    name:'depo_fanout-%d'%n,
    data:{
        multiplicity:2,
        tags: [],
    }}, nin=1, nout=2) for n in anode_iota];

// local deposplats = [sim.make_ductor('splat%d'%n, tools.anodes[n], tools.pirs[0], 'DepoSplat', 'deposplat%d'%n) for n in anode_iota] ;

// The approximated sim+sigproc
local splat = function(params, tools, anode, name=null) {
    local apaid = anode.data.ident,
    local sufix = if std.type(name) == "null" then apaid else name,
    local bg = g.pnode({
        type:'DepoBagger',
        name: sufix,
        data: {
            gate: [params.sim.ductor.start_time,
                   params.sim.ductor.start_time+params.sim.ductor.readout_time],
        },
    }, nin=1, nout=1),
    local sp = g.pnode({
        type: 'DepoFluxSplat',
        name: sufix,
        data: {
            anode: wc.tn(anode),
            field_response: wc.tn(tools.field), // for speed and origin
            sparse: true,
            tick: 0.5 * wc.us,
            window_start: -0.205 * wc.ms, // TODO: there seems to be an offset between sim and rec?
            window_duration: self.tick * params.daq.nticks,
            reference_time: 0.0,
            "smear_long": [
                3.5,
                3.5,
                3.5
              ],
            "smear_tran": [
                0.4, //can use 0.3 if too strong, otherwise can use 0.4
                0.4,
                0.14
            ]
        },
    }, nin=1, nout=1, uses=[anode, tools.field]),
    local rt = g.pnode({
        type: 'Retagger',
        name: sufix,
        data: {
            // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
            tag_rules: [{
                // Retagger also handles "frame" and "trace" like fanin/fanout
                // merge separately all traces like gaussN to gauss.
                frame: {
                ".*": "deposplat%d" % apaid
                },
                merge: {
                ".*": "deposplat%d" % apaid
                },
            }],
        },
    }, nin=1, nout=1),
    ret: g.pipeline([bg, sp, rt],"%s-%s" % [bg.name, sp.name]),
}.ret;
    
local deposplats = [splat(params, tools, tools.anodes[n]) for n in anode_iota] ;

local hio_truth = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_truth%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['deposplat%d'%n],
        filename: "g4-tru-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        tick0: 0,
        nticks: 3427,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];


local multipass2 = [
  g.pipeline([
               deposplats[n],
               hio_truth[n],
             ],
             'multipass%d' % n)
  for n in anode_iota
];


local f_sp = import 'pgrapher/experiment/sbnd/funcs.jsonnet';
local outtags = ['orig%d' % n for n in anode_iota];
local bi_manifold2 = f_sp.fanpipe('DepoFanout', multipass2, 'FrameFanin');

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
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      // nticks: params.daq.nticks,
      chanmaskmaps: ['bad'],
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
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener'],

      // this may be needed to convert the decon charge [units:e-] to be consistent with the LArSoft default ?unit? e.g. decon charge * 0.005 --> "charge value" to GaussHitFinder
      frame_scale: [1.0, 1.0],
       nticks: params.daq.nticks,
      chanmaskmaps: [],
      //nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};


local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(5632 * n, 5632 * (n + 1) - 1),
      //tags: ['orig%d' % n], // traces tag
    },
  }, nin=1, nout=1)
  for n in anode_iota
];


local nfsp_pipes = [
  g.pipeline([
               chsel_pipes[n],
               //sinks.orig_pipe[n],

               //nf_pipes[n],
               //sinks.raw_pipe[n],

               sp_pipes[n],
              //  hio_sp[n],
               //sinks.decon_pipe[n],
               //sinks.threshold_pipe[n],
               // sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
             ],
             'nfsp_pipe_%d' % n)
  for n in anode_iota
];

local f_sp = import 'pgrapher/experiment/sbnd/funcs.jsonnet';
//local fanpipe = f_sp.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf'); // commented Ewerton 2023-05-24
local fanpipe = f_sp.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf_mod');   //added Ewerton 2023-05-24

local retagger_sp = g.pnode({
  type: 'Retagger',
  name: 'sp',  //added Ewerton 2023-05-24
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
      },
    }],
  },
}, nin=1, nout=1);

local sink_sp = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);

local graph = g.pipeline([
wcls_input_sim.depos,         //sim
drifter,                     //sim
wcls_deposetsimchannel_sink,    //sim
bi_manifold2,                   //sim // changed Ewerton 2023-05-28
sink_sp                         //sp
]);

// local save_simdigits = std.extVar('save_simdigits');
local save_simdigits = "false";

local app = {
  type: 'TbbFlow',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
