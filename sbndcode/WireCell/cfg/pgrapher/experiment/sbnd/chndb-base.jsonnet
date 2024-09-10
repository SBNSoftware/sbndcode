// Base channel noise DB object configuration for microboone
// This does not include any run dependent RMS cuts.
// See chndb.jsonnet

local handmade = import 'chndb-resp.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params, anode, field, n, rms_cuts=[])
  {
    anode: wc.tn(anode),
    field_response: wc.tn(field),

    tick: params.daq.tick,

    // This sets the number of frequency-domain bins used in the noise
    // filtering.  It is not necessarily true that the time-domain
    // waveforms have the same number of ticks.  This must be non-zero.
    nsamples: params.nf.nsamples,

    // For MicroBooNE, channel groups is a 2D list.  Each element is
    // one group of channels which should be considered together for
    // coherent noise filtering.
    groups: [std.range(   0 + n * 5638 + g*32,    0 + n * 5638 + (g+1)*32 - 1) for g in std.range(0,149)] +
            [std.range(4806 + n * 5638 + g*32, 4806 + n * 5638 + (g+1)*32 - 1) for g in std.range(0,25)] ,
    

    // Externally determined "bad" channels.
    //
    // Dead channels: 3232:3263 (inclusive) (East V).   4160:4191 (East Y)
    // Shorted channels:  7169 (West U), 8378 (West V).
    // There are four physically missing wires ( = bad channels) due to combs, in the center of each 1/2 APA.
    // They are 4374 and 5231 (East Y), 10012 and 10869 (West Y).
    // So in total, there are 76 bad channels.
    // 
    //bad: [],
    bad: [546, 607] + std.range(3232, 3263) + std.range(4160, 4191) + [4374, 4800, 4801, 4802, 4803, 4804, 4805, 5060, 5231, 5636, 5637, 7169, 8378, 8574, 10012, 10869, 10438, 10439, 10440, 10441, 10442, 10443],

    // Overide defaults for specific channels.  If an info is
    // mentioned for a particular channel in multiple objects in this
    // list then last mention wins.
    channel_info: [

      // First entry provides default channel info across ALL
      // channels.  Subsequent entries override a subset of channels
      // with a subset of these entries.  There's no reason to
      // repeat values found here in subsequent entries unless you
      // wish to change them.
      {
        //channels: std.range(n * 2560, (n + 1) * 2560 - 1),
        channels: std.range(n * 5638, n * 5638 + 5637),
        nominal_baseline: 2001.0,  // adc count [879.5 mV]
        gain_correction: 1.0,  // unitless
        response_offset: 0.0,  // ticks?
        pad_window_front: 10,  // ticks?
        pad_window_back: 10,  // ticks?
        decon_limit: 0.02,
        decon_limit1: 0.09,
        adc_limit: 15,
        roi_min_max_ratio: 0.8, // default 0.8
        min_rms_cut: 1.0,  // units???
        max_rms_cut: 30.0,  // units???

        // parameter used to make "rcrc" spectrum
        rcrc: 1.1 * wc.millisecond, // 1.1 for collection, 3.3 for induction
        rc_layers: 1, // default 2

        // parameters used to make "config" spectrum
        reconfig: {},

        // list to make "noise" spectrum mask
        freqmasks: [],

        // field response waveform to make "response" spectrum.
        response: {},

      },


      {
        //channels: { wpid: wc.WirePlaneId(wc.Ulayer) },
//        channels: std.range(n * 2560, n * 2560 + 800- 1),
        channels: std.range(n * 5638, n * 5638 + 1984-1),
        freqmasks: [
//          { value: 1.0, lobin: 0, hibin: $.nsamples - 1 },
//          { value: 0.0, lobin: 169, hibin: 173 },
//          { value: 0.0, lobin: 513, hibin: 516 },
        ],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Ulayer) },
        /// this uses hard-coded waveform.
        response: { waveform: handmade.u_resp, waveformid: wc.Ulayer },
        response_offset: 125.6, // offset of the negative peak
        pad_window_front: 20,
        decon_limit: 0.02,
        decon_limit1: 0.07,
        roi_min_max_ratio: 3.0,
      },

      {
        //channels: { wpid: wc.WirePlaneId(wc.Vlayer) },
//        channels: std.range(n * 2560 + 800, n * 2560 + 1600- 1),
        channels: std.range(n * 5638 + 1984, n * 5638 + 3968-1),
        freqmasks: [
 //         { value: 1.0, lobin: 0, hibin: $.nsamples - 1 },
 //         { value: 0.0, lobin: 169, hibin: 173 },
 //         { value: 0.0, lobin: 513, hibin: 516 },
        ],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Vlayer) },
        /// this uses hard-coded waveform.
        response: { waveform: handmade.v_resp, waveformid: wc.Vlayer },
        response_offset: 129.5,
        decon_limit: 0.01,
        decon_limit1: 0.08,
        roi_min_max_ratio: 1.5,
      },

      local freqbinner = wc.freqbinner(params.daq.tick, params.nf.nsamples);
      local harmonic_freqs = [
        //f*wc.kilohertz for f in
        // [51.5, 102.8, 154.2, 205.5, 256.8, 308.2, 359.2, 410.5, 461.8, 513.2, 564.5, 615.8]
        //[51.5, 77.2, 102.8, 128.5, 154.2, 180.0, 205.5, 231.5, 256.8, 282.8, 308.2, 334.0, 359.2, 385.5, 410.5, 461.8, 513.2, 564.5, 615.8, 625.0]
      ];

      {
        //channels: { wpid: wc.WirePlaneId(wc.Wlayer) },
//        channels: std.range(n * 2560 + 1600, n * 2560 + 2560- 1),
        channels: std.range(n * 5638 + 3968, n * 5638 + 5638-1),
        nominal_baseline: 650,
        decon_limit: 0.05,
        decon_limit1: 0.08,
        freqmasks: freqbinner.freqmasks(harmonic_freqs, 5.0*wc.kilohertz),
      },
      
    ] + rms_cuts,
  }
