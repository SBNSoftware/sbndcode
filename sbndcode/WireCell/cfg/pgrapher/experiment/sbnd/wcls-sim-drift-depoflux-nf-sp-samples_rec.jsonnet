// This is a main entry point for configuring a wire-cell CLI job to
// simulate SBND.  It is simplest signal-only simulation with
// one set of nominal field response function.  

local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
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

local hio_sp = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_sp%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['loose_lf%d' % n
        #, 'tight_lf%d' % n
        #, 'cleanup_roi%d' % n
        #, 'break_roi_1st%d' % n
        #, 'break_roi_2nd%d' % n
        #, 'shrink_roi%d' % n
        #, 'extend_roi%d' % n
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

local sp_override = {
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
    decon_charge_tag: "",
    gauss_tag: "",
    wiener_tag: "",
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

local multipass2 = [
  g.pipeline([
               sn_pipes[n],
               //sinks.orig_pipe[n],
               
               nf_pipes[n], // NEED to include this pipe for channelmaskmaps    
               //sinks.raw_pipe[n], 
               
               sp_pipes[n], 
               hio_sp[n],

               //sinks.decon_pipe[n],
               //sinks.threshold_pipe[n],
             ],
             'multipass%d' % n)
  for n in anode_iota
];


local f_sp = import 'pgrapher/experiment/sbnd/funcs.jsonnet';

local outtags = ['orig%d' % n for n in anode_iota];
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

};


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


local nfsp_pipes = [
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
wcls_input_sim.deposet,         //sim
setdrifter,                     //sim
wcls_depoflux_writer,           //sim
bi_manifold2,                   //sim
retagger_sp,                    //sp
wcls_output_sp.sp_signals,      //sp
sink_sp                         //sp
]);

local save_simdigits = std.extVar('save_simdigits');

local app = {
  type: 'TbbFlow',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
  
