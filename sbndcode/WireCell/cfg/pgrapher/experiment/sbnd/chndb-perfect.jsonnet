// Perfect channel noise DB object configuration for microboone.


local wc = import 'wirecell.jsonnet';

function(params, anode, field, n)
  {
    anode: wc.tn(anode),
    field_response: wc.tn(field),

    tick: params.daq.tick,

    // This sets the number of frequency-domain bins used in the noise
    // filtering.  It is expected that time-domain waveforms have the
    // same number of samples.
    nsamples: params.nf.nsamples,

    // For MicroBooNE, channel groups is a 2D list.  Each element is
    // one group of channels which should be considered together for
    // coherent noise filtering.
    //groups: [std.range(g*48, (g+1)*48-1) for g in std.range(0,171)],
    groups: [std.range(n * 2560 + u * 40, n * 2560 + (u + 1) * 40 - 1) for u in std.range(0, 19)]
            + [std.range(n * 2560 + 800 + v * 40, n * 2560 + 800 + (v + 1) * 40 - 1) for v in std.range(0, 19)]
            + [std.range(n * 2560 + 1600 + w * 48, n * 2560 + 1600 + (w + 1) * 48 - 1) for w in std.range(0, 19)],

    // Externally determined "bad" channels.
    //bad: [],
    //bad: // shorted-U
    //     [296] + std.range(298, 315) + [317] + std.range(319,327) + std.range(336, 337)
    //     + std.range(343, 345) + std.range(348, 351) + std.range(376, 400) + std.range(410, 445)
    //     + std.range(447, 484) + std.range(501, 503) + std.range(505, 520) + std.range(522, 524)
    //     + std.range(536, 559) + std.range(561, 592) + std.range(595, 598) + std.range(600, 632)
    //     + std.range(634, 652) + [654] + std.range(656,671)
    //     // inverse "V" due to disconnected MB
    //     + std.range(864, 911)
    //     + std.range(3936,3983)
    //     // shorted-Y
    //     + std.range(7136, 7199) + std.range(7201, 7214) + std.range(7216, 7263),

    // Overide defaults for specific channels.  If an info is
    // mentioned for a particular channel in multiple objects in this
    // list then last mention wins.
    /*channel_info: [

        // First entry provides default channel info across ALL
        // channels.  Subsequent entries override a subset of channels
        // with a subset of these entries.  There's no reason to
        // repeat values found here in subsequent entries unless you
        // wish to change them.
        {
            channels: std.range(0, 2400 + 2400 + 3456 - 1),
            nominal_baseline: 2048.0,  // adc count
            gain_correction: 1.0,     // unitless
            response_offset: 0.0,      // ticks?
            pad_window_front: 10,     // ticks?
            pad_window_back: 10,      // ticks?
     decon_limit: 0.02,
     decon_limit1: 0.09,
     adc_limit: 15,
            min_rms_cut: 1.0,
            max_rms_cut: 5.0,

            // parameter used to make "rcrc" spectrum
            rcrc: 1.0*wc.millisecond,

            // parameters used to make "config" spectrum
            reconfig : {},

            // list to make "noise" spectrum mask
            freqmasks: [],

            // field response waveform to make "response" spectrum.
            response: {},

        },

        {
            channels: {wpid: wc.WirePlaneId(wc.Ulayer)},
            pad_window_front: 20,
     decon_limit: 0.02,
     decon_limit1: 0.09,
        },

        {
            channels: {wpid: wc.WirePlaneId(wc.Vlayer)},
     decon_limit: 0.01,
     decon_limit1: 0.08,
     },

        {
            channels: {wpid: wc.WirePlaneId(wc.Wlayer)},
            nominal_baseline: 400.0,
     decon_limit: 0.05,
     decon_limit1: 0.08,
        },

        {                       // these are before hardware fix
            channels: params.nf.misconfigured.channels,
            reconfig: {
                from: {gain:  params.nf.misconfigured.gain,
                       shaping: params.nf.misconfigured.shaping},
                to:   {gain: params.elec.gain,
                       shaping: params.elec.shaping},
            }
        },
    ],*/
  }
