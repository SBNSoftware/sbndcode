// This is a main entry point to configure a WC/LS job that applies
// noise filtering and signal processing to existing RawDigits.  The
// FHiCL is expected to provide the following parameters as attributes
// in the "params" structure.
//
// epoch: the hardware noise fix expoch: "before", "after", "dynamic" or "perfect"
// reality: whether we are running on "data" or "sim"ulation.
// raw_input_label: the art::Event inputTag for the input RawDigit
//
// see the .fcl of the same name for an example
//
// Manual testing, eg:
//
// jsonnet -V reality=data -V epoch=dynamic -V raw_input_label=daq \\
//         -V signal_output_form=sparse \\
//         -J cfg cfg/pgrapher/experiment/uboone/wcls-nf-sp.jsonnet
//
// jsonnet -V reality=sim -V epoch=perfect -V raw_input_label=daq \\
//         -V signal_output_form=sparse \\
//         -J cfg cfg/pgrapher/experiment/uboone/wcls-nf-sp.jsonnet


//local epoch = std.extVar('epoch');  // eg "dynamic", "after", "before", "perfect"
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
local raw_input_label = std.extVar('raw_input_label');  // eg "daq"
local use_paramresp = std.extVar('use_paramresp');  // eg "true" or "false"
local use_dnnroi = std.extVar('use_dnnroi');
local nchunks = std.extVar('nchunks');
local tick_per_slice = std.extVar('tick_per_slice');
local dnnroi_model_p0 = std.extVar('dnnroi_model_p0');
local dnnroi_model_p1 = std.extVar('dnnroi_model_p1');

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/experiment/sbnd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

local data_params = import 'params.jsonnet';

local params = data_params {
    daq: super.daq { // <- super.daq overrides default values
      // Set the waveform sample length, eg, 6000, 15000, "auto"
      nticks: std.extVar('nticks'),
    },
    # provides filename for parametrized per channel electronics response correction
    # default is: use_parmresp = false
    files: super.files {
      chresp: if (use_paramresp == 'true') then "sbnd-params-chresp-ideal-v1.json.bz2" else null,
    },
};


local tools = tools_maker(params);


local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

// Collect the WC/LS input converters for use below.  Make sure the
// "name" argument matches what is used in the FHiCL that loads this
// file.  In particular if there is no ":" in the inputer then name
// must be the emtpy string.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      // nticks: params.daq.nticks,
    },
  }, nin=0, nout=1),

};


// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.

local wcls_output = {
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
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener'],

      // this may be needed to convert the decon charge [units:e-] to be consistent with the LArSoft default ?unit? e.g. decon charge * 0.005 --> "charge value" to GaussHitFinder
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
      frame_scale: [0.02, 0.02, 0.02],
      nticks: params.daq.nticks,

      chanmaskmaps: ['dnnspbad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

local base = import 'pgrapher/experiment/sbnd/chndb-base.jsonnet';

local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(5638 * n, 5638 * (n + 1) - 1),
      //tags: ['orig%d' % n], // traces tag
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local nf_maker = import 'pgrapher/experiment/sbnd/nf-data.jsonnet'; //added Ewerton 2024-08
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_maker = import 'pgrapher/experiment/sbnd/sp-lowthreshold.jsonnet';
local sp_override = if use_dnnroi then {
    sparse: true,
    use_roi_debug_mode: true,
    save_negtive_charge: true, // TODO: no negative charge in gauss, default is false
    use_multi_plane_protection: true,
    mp_tick_resolution: 4,
    tight_lf_tag: "",
    // loose_lf_tag: "",
    cleanup_roi_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
    // m_decon_charge_tag: "",
} else {
    sparse: true,
};
//local sp = sp_maker(params, tools, { sparse: sigoutform == 'sparse' });
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];
local dnnroi = import 'dnnroi.jsonnet';
local ts_p0 = {
    type: "TorchService",
    name: "dnnroi_p0",
    tick_per_slice: tick_per_slice, 
    data: {
        model: dnnroi_model_p0,
        device: "cpu",
        concurrency: 1,
    },
};

local ts_p1 = {
    type: "TorchService",
    name: "dnnroi_p1",
    tick_per_slice: tick_per_slice, 
    data: {
        model: dnnroi_model_p1,
        device: "cpu",
        concurrency: 1,
    },
};

local magoutput = 'sbnd-data-check.root';
local magnify = import 'pgrapher/experiment/sbnd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local fanout = function (name, multiplicity=2)
  g.pnode({
    type: 'FrameFanout',
    name: name,
    data: {
        multiplicity: multiplicity
    },
  }, nin=1, nout=multiplicity);

local sp_fans = [fanout("sp_fan_%d" % n) for n in std.range(0, std.length(tools.anodes) - 1)];
local dnnroi_pipes = [ dnnroi(tools.anodes[n], ts_p0, ts_p1, output_scale=1, nchunks=nchunks) for n in std.range(0, std.length(tools.anodes) - 1) ];
local nfsp_pipes = if use_dnnroi then 
[
  g.intern(
    innodes=[chsel_pipes[n]],
    outnodes=[dnnroi_pipes[n],sinks.decon_pipe[n]],
    centernodes=[nf_pipes[n], sp_pipes[n], sp_fans[n], sinks.decon_pipe[n]],
    edges=[
      g.edge(chsel_pipes[n], nf_pipes[n], 0, 0),
      g.edge(nf_pipes[n], sp_pipes[n], 0, 0),
      g.edge(sp_pipes[n], sp_fans[n], 0, 0),
      g.edge(sp_fans[n], dnnroi_pipes[n], 0, 0),
      g.edge(sp_fans[n], sinks.decon_pipe[n], 1, 0),
    ],
    iports=chsel_pipes[n].iports,
    oports=dnnroi_pipes[n].oports+sinks.decon_pipe[n].oports,
    name='nfsp_pipe_%d' % n,
  )
  for n in std.range(0, std.length(tools.anodes) - 1)
]
else
[
  g.pipeline([
               chsel_pipes[n],
               //sinks.orig_pipe[n],

               nf_pipes[n],
               //sinks.raw_pipe[n],

               sp_pipes[n],
               //sinks.decon_pipe[n],
               //sinks.threshold_pipe[n],
               // sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
             ],
             'nfsp_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local fanpipe = f.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf');

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

local retag_dnnroi = retagger("retag_dnnroi");
local retag_sp = retagger("retag_sp");


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

local graph = if use_dnnroi then
 g.intern(
  innodes=[wcls_input.adc_digits],
  outnodes=[],
  centernodes=nfsp_pipes+[fanout_apa, retag_dnnroi, retag_sp, fanin_apa_dnnroi, fanin_apa_sp, wcls_output.sp_signals, wcls_output.dnnsp_signals, sink_dnnroi, sink_sp],
  edges=[
    g.edge(wcls_input.adc_digits, fanout_apa, 0, 0),
    g.edge(fanout_apa, nfsp_pipes[0], 0, 0),
    g.edge(fanout_apa, nfsp_pipes[1], 1, 0),
    g.edge(nfsp_pipes[0], fanin_apa_dnnroi, 0, 0),
    g.edge(nfsp_pipes[1], fanin_apa_dnnroi, 0, 1),
    g.edge(fanin_apa_dnnroi, retag_dnnroi, 0, 0),
    g.edge(retag_dnnroi, wcls_output.dnnsp_signals, 0, 0),
    g.edge(wcls_output.dnnsp_signals, sink_dnnroi, 0, 0),
    g.edge(nfsp_pipes[0], fanin_apa_sp, 1, 0),
    g.edge(nfsp_pipes[1], fanin_apa_sp, 1, 1),
    g.edge(fanin_apa_sp, retag_sp, 0, 0),
    g.edge(retag_sp, wcls_output.sp_signals, 0, 0),
    g.edge(wcls_output.sp_signals, sink_sp, 0, 0),
  ]
)
else
g.pipeline([wcls_input.adc_digits, 
fanpipe, 
retag_sp, 
wcls_output.sp_signals, 
sink_sp]);

local app = {
  type: 'TbbFlow',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]

