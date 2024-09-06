// This provides multiple MagnifySink for e.g. protoDUNE

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// multiple MagnifySink
// tagn (n = 0, 1, ... 5) for anode[n]
// FrameFanin tags configured in sim.jsonnet
function(tools, outputfile) {

  local nanodes = std.length(tools.anodes),

  local magorig = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magorig%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['orig%d' % n],
        trace_has_tag: false,   // traces from source have NO tag
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magraw = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magraw%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['raw%d' % n],
        trace_has_tag: true,
        cmmtree: [["noisy", "T_noisy%d"%n],
                  ["sticky", "T_stky%d"%n],
                  ["ledge", "T_ldg%d"%n],
                  ["harmonic", "T_hm%d"%n] ], // maskmap in nf.jsonnet 
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magdecon = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdecon%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['gauss%d' % n, 'wiener%d' % n],
        trace_has_tag: true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magdebug = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdebug%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['tight_lf%d' %n, 'loose_lf%d' %n, 'cleanup_roi%d' %n,
                 'break_roi_1st%d' %n, 'break_roi_2nd%d' %n,
                 'shrink_roi%d' %n, 'extend_roi%d' %n],
        trace_has_tag: true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magthr = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magthr%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        summaries: ['wiener%d' % n],  // note that if tag set, each apa should have a tag set for FrameFanin
        summary_operator: { ['wiener%d' % n]: 'set' },  // []: obj comprehension
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],


  return: {
    orig_pipe: [g.pipeline([magorig[n]], name='magorigpipe%d' % n) for n in std.range(0, nanodes - 1)],
    raw_pipe: [g.pipeline([magraw[n]], name='magrawpipe%d' % n) for n in std.range(0, nanodes - 1)],
    decon_pipe: [g.pipeline([magdecon[n]], name='magdeconpipe%d' % n) for n in std.range(0, nanodes - 1)],
    debug_pipe: [g.pipeline([magdebug[n]], name='magdebugpipe%d' % n) for n in std.range(0, nanodes - 1)],
    threshold_pipe: [g.pipeline([magthr[n]], name='magthrpipe%d' % n) for n in std.range(0, nanodes - 1)],
  },


}.return
