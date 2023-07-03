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
    // groups: [std.range(g*48, (g+1)*48-1) for g in std.range(0,171)],
    groups: [std.range(n * 2560 + u * 40, n * 2560 + (u + 1) * 40 - 1) for u in std.range(0, 19)]
            + [std.range(n * 2560 + 800 + v * 40, n * 2560 + 800 + (v + 1) * 40 - 1) for v in std.range(0, 19)]
            + [std.range(n * 2560 + 1600 + w * 48, n * 2560 + 1600 + (w + 1) * 48 - 1) for w in std.range(0, 19)],


    // Externally determined "bad" channels.
    bad: [
    # CE group: Inactive FE
     4411,  # femb515x12
     4412,  # femb515x13
     9990,  # femb605x10
    11842,  # femb120x03
    # CE group: Broken connection
        1,  # femb311u39
      400,  # femb301u40
      401,  # femb301u39
      800,  # femb320v01
      801,  # femb320v02
      876,  # femb319v37
     1200,  # femb310v01
     2961,  # femb501u39
     5321,  # femb216u39
     5363,  # femb217u37
     6132,  # femb215v13
     7058,  # femb213x03
     7295,  # femb202x01
     7681,  # femb611u39
     8080,  # femb601u40
     8328,  # femb607u32
     8480,  # femb620v01
     9282,  # femb620x03
     9283,  # femb620x04
     9736,  # femb611x25
     9854,  # femb602x02
    10800,  # femb105u40
    11024,  # femb110u16
    11457,  # femb110v18
    11459,  # femb110v20
    11463,  # femb110v24
    11469,  # femb110v30: bad in 4875-185-1500, ok in 5803-76-3200
    11517,  # femb109v38
    11669,  # femb105v30
    11679,  # femb105v40
    12756,  # femb110x44
    12801,  # femb411u39
    13001,  # femb416u39
    13081,  # femb418u39
    # CE group: ENC > 2000e
     4410,  # femb515x11: High noise, no signal 5008-76
    #-----
    # CE group excessive sticky
    #femb318x
     1719,  # femb318x24
     5125,  # femb211u35
     7551,  # femb208x33
     7190,  # femb211x39
     7194,  # femb211x43
     7918,  # femb616u02, sticky pedestal (three peaks)
    #-----
    # CE group: good.
    # femb311
        2,  # femb311u38, no signal
        4,  # femb311u36, very sticky pedestal 5308-76
     1632,  # femb320x33, very sticky pedestal 5308-76
     2169,  # femb302x07, Mostly stuck on one bad code, 5308-76
     2450,  # femb308x14, Very noisy (1000 ADC) in run 5759 (20nov2019)
     3541,  # femb516v22, very sticky--signal near zero half the time (5308-81)
     3543,  # femb516v24, very sticky--signal near zero half the time (5308-81)
     3661,  # femb513v22, most signal near zero (5308-81)
     3663,  # femb513v24, most signal near zero (5308-81)
     4061,  # femb503v22, most signal near zero (5308-81)
     4063,  # femb503v24, most signal near zero (5308-81)
     4141,  # femb501v22, signal near zero half the time (5308-81)
     4143,  # femb501v24, signal sometimes near zero (5308-81)
     4377,  # femb516x26, very sticky pedestal
     4379,  # femb516x28, very sticky pedestal
     4381,  # femb516x30, very sticky pedestal
     4383,  # femb516x32, very sticky pedestal
     4385,  # femb516x34, very sticky pedestal
     4387,  # femb516x36, very sticky pedestal
     4521,  # femb513x26, very sticky pedestal
     4523,  # femb513x28, very sticky pedestal
     4525,  # femb513x30, very sticky pedestal
     4527,  # femb513x32, very sticky pedestal
     4529,  # femb513x34, very sticky pedestal
     4531,  # femb513x36, very sticky pedestal
     4652,  # femb501x36, very sticky pedestal
     4654,  # femb501x34, very sticky pedestal
     4656,  # femb501x32, very sticky pedestal
     4658,  # femb501x30, very sticky pedestal
     4660,  # femb501x28, very sticky pedestal
     4658,  # femb501x26, very sticky pedestal
     4748,  # femb503x36, very sticky pedestal
     4750,  # femb503x34, very sticky pedestal
     4752,  # femb503x32, very sticky pedestal
     4754,  # femb503x30, very sticky pedestal
     4756,  # femb503x28, very sticky pedestal
     4758,  # femb503x26, very sticky pedestal
     5361,  # femb217u39, no signal
     7680,  # femb611u40: No signal in 5308-76, end wire
     8501,  # femb620v22, very sticky pedestal
     8503,  # femb620v24, very sticky pedestal
     8821,  # femb612v22, very sticky pedestal
     8823,  # femb612v24, very sticky pedestal
     9261,  # femb601v22, very sticky pedestal
     9263,  # femb601v24, very sticky pedestal
     9305,  # femb620x26, very sticky pedestal
     9307,  # femb620x28, very sticky pedestal
     9309,  # femb620x30, very sticky pedestal
     9311,  # femb620x32, very sticky pedestal
     9313,  # femb620x34, very sticky pedestal
     9315,  # femb620x36, very sticky pedestal
     9689,  # femb612x26, very sticky pedestal
     9691,  # femb612x28, very sticky pedestal
     9693,  # femb612x30, very sticky pedestal
     9695,  # femb612x32, very sticky pedestal
     9697,  # femb612x34, very sticky pedestal
     9699,  # femb612x36, very sticky pedestal
     9772,  # femb601x26, very sticky pedestal
     9774,  # femb601x28, very sticky pedestal
     9776,  # femb601x30, very sticky pedestal
     9778,  # femb601x32, very sticky pedestal
     9780,  # femb601x34, very sticky pedestal
     9782,  # femb601x36, very sticky pedestal
    10102,  # femb608x42, mostly stuck on one code
    10189,  # femb609x03, mostly stuck on one code
    10697,  # femb102u23, mostly stuck on a few classic codes
    10907,  # femb107u13, mostly stuck on one code
    11203,  # femb116v04, stuck on many classic codes
    11270,  # femb115v31, stuck on many classic codes
    11902,  # femb119x15, stuck on two classic codes
    12324,  # femb101x44, stuck on many classic codes
    12333,  # femb101x35, stuck on many classic codes
    12744,  # femb109x08, stuck on many classic codes
    13363,  # femb405u37, very noisy, nosignal 5308-76-4800

    #-----
    # These 16 channels are an intermitently bad ASIC.
    # Matt W. 19oct2018.
    # femb316u
      200,   # femb316u40
      202,   # femb316u38
      204,   # femb316u36
      206,   # femb316u34
      208,   # femb316u32
    # femb316v
      991,   # femb316v32
      993,   # femb316v34
      995,   # femb316v36
      997,   # femb316v38
      999,   # femb316v40
    # femb316x
     1829,   # femb316x38
     1831,   # femb316x40
     1833,   # femb316x42
     1835,   # femb316x44
     1837,   # femb316x46
     1839    # femb316x48
    #-----
],

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
        channels: std.range(n * 2560, (n + 1) * 2560 - 1),
        nominal_baseline: 2048.0,  // adc count
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
	channels: std.range(n * 2560, n * 2560 + 800- 1),
	freqmasks: [
          { value: 1.0, lobin: 0, hibin: $.nsamples - 1 },
          { value: 0.0, lobin: 169, hibin: 173 },
          { value: 0.0, lobin: 513, hibin: 516 },
        ],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Ulayer) },
        /// this uses hard-coded waveform.
        response: { waveform: handmade.u_resp, waveformid: wc.Ulayer },
        response_offset: 120, // offset of the negative peak
        pad_window_front: 20,
        decon_limit: 0.02,
        decon_limit1: 0.07,
        roi_min_max_ratio: 3.0,
      },

      {
        //channels: { wpid: wc.WirePlaneId(wc.Vlayer) },
	channels: std.range(n * 2560 + 800, n * 2560 + 1600- 1),
        freqmasks: [
          { value: 1.0, lobin: 0, hibin: $.nsamples - 1 },
          { value: 0.0, lobin: 169, hibin: 173 },
          { value: 0.0, lobin: 513, hibin: 516 },
        ],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Vlayer) },
        /// this uses hard-coded waveform.
        response: { waveform: handmade.v_resp, waveformid: wc.Vlayer },
        response_offset: 124,
        decon_limit: 0.01,
        decon_limit1: 0.08,
        roi_min_max_ratio: 1.5,
      },

      local freqbinner = wc.freqbinner(params.daq.tick, params.nf.nsamples);
      local harmonic_freqs = [f*wc.kilohertz for f in
        // [51.5, 102.8, 154.2, 205.5, 256.8, 308.2, 359.2, 410.5, 461.8, 513.2, 564.5, 615.8]
        [51.5, 77.2, 102.8, 128.5, 154.2, 180.0, 205.5, 231.5, 256.8, 282.8, 308.2, 334.0, 359.2, 385.5, 410.5, 461.8, 513.2, 564.5, 615.8, 625.0]
      ];
      
      {
        //channels: { wpid: wc.WirePlaneId(wc.Wlayer) },
	channels: std.range(n * 2560 + 1600, n * 2560 + 2560- 1),
        nominal_baseline: 400.0,
        decon_limit: 0.05,
        decon_limit1: 0.08,
        freqmasks: freqbinner.freqmasks(harmonic_freqs, 5.0*wc.kilohertz),
      },

    ] + rms_cuts,
  }
