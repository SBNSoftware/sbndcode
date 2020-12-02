// This provides signal processing related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// BIG FAT FIXME: we are taking from uboone.  If SBND needs tuning do
// four things: 0) read this comment, 1) cp this file into sbnd/, 2)
// fix the import and 3) delete this comment.
local spfilt = import 'pgrapher/experiment/sbnd/sp-filters.jsonnet';

function(params, tools, override = {}) {

  local pc = tools.perchanresp_nameuses,

  local resolution = params.adc.resolution,
  local fullscale = params.adc.fullscale[1] - params.adc.fullscale[0],
  local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale,

  // SBND needs a per-anode sigproc
  make_sigproc(anode, name=null):: g.pnode({
    type: 'OmnibusSigProc',
    name:
      if std.type(name) == 'null'
      then anode.name + 'sigproc%d' % anode.data.ident
      else name,

    data: {
      /**  
       *  Default SP parameters (till May 2019)
       */
      // anode: wc.tn(anode),
      // field_response: wc.tn(tools.field),
      // per_chan_resp: pc.name,
      // fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
      // postgain: 1,  // default 1.2
      // ADC_mV: 4096 / (1400.0 * wc.mV),  // default 4096/2000
      // r_fake_signal_low_th: 400,  // default 500
      // r_fake_signal_high_th: 800,  // default 1000
      // r_fake_signal_low_th_ind_factor: 1.5,  // default 1
      // r_fake_signal_high_th_ind_factor: 1.5,  // default 1
      // troi_col_th_factor: 5.0,  // default 5
      // troi_ind_th_factor: 3.5,  // default 3
      // r_th_factor: 3.5,  // default 3

      /**  
       *  Optimized SP parameters (May 2019)
       *  Associated tuning in sp-filters.jsonnet
       */
      anode: wc.tn(anode),
      field_response: wc.tn(tools.field),
      elecresponse: wc.tn(tools.elec_resp),
      ftoffset: 0.0, // default 0.0
      ctoffset: 1.0*wc.microsecond, // default -8.0
      per_chan_resp: pc.name,
      fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
      postgain: 1.0,  // default 1.2
      ADC_mV: ADC_mV_ratio, // 4096 / (1400.0 * wc.mV), 
      troi_col_th_factor: 5.0,  // default 5
      troi_ind_th_factor: 3.0,  // default 3
      lroi_rebin: 6, // default 6
      lroi_th_factor: 3.5, // default 3.5
      lroi_th_factor1: 0.7, // default 0.7
      lroi_jump_one_bin: 1, // default 0

      r_th_factor: 3.0,  // default 3
      r_fake_signal_low_th: 375,  // default 500
      r_fake_signal_high_th: 750,  // default 1000
      r_fake_signal_low_th_ind_factor: 1.0,  // default 1
      r_fake_signal_high_th_ind_factor: 1.0,  // default 1      
      r_th_peak: 3.0, // default 3.0
      r_sep_peak: 6.0, // default 6.0
      r_low_peak_sep_threshold_pre: 1200, // default 1200


      // frame tags
      wiener_tag: 'wiener%d' % anode.data.ident,
      wiener_threshold_tag: 'threshold%d' % anode.data.ident,
      decon_charge_tag: 'decon_charge%d' % anode.data.ident,
      gauss_tag: 'gauss%d' % anode.data.ident,

      use_roi_debug_mode: false,
      tight_lf_tag: 'tight_lf%d' % anode.data.ident,
      loose_lf_tag: 'loose_lf%d' % anode.data.ident,
      cleanup_roi_tag: 'cleanup_roi%d' % anode.data.ident,
      break_roi_loop1_tag: 'break_roi_1st%d' % anode.data.ident,
      break_roi_loop2_tag: 'break_roi_2nd%d' % anode.data.ident,
      shrink_roi_tag: 'shrink_roi%d' % anode.data.ident,
      extend_roi_tag: 'extend_roi%d' % anode.data.ident,

      use_multi_plane_protection: false,
      mp3_roi_tag: 'mp3_roi%d' % anode.data.ident,
      mp2_roi_tag: 'mp2_roi%d' % anode.data.ident,
      
      isWrapped: false,
      // process_planes: [0, 2],

    } + override,
  }, nin=1, nout=1, uses=[anode, tools.field, tools.elec_resp] + pc.uses + spfilt),

}
