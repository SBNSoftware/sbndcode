// This provides some util functions.

local g = import 'pgraph.jsonnet';

{
  // Build a fanout-[pipelines]-fanin graph.  pipelines is a list of
  // pnode objects, one for each spine of the fan.
  fanpipe:: function(fout, pipelines, fin, name='fanpipe', outtags=[]) {

    local fanmult = std.length(pipelines),
    local fannums = std.range(0, fanmult - 1),

    local fanout = g.pnode({
      type: fout,
      name: name,
      data: {
        multiplicity: fanmult,
        tag_rules: [  // example in gen/test/test_fans.jsonnet
          {
            frame: {
              //'.*': 'number%d' % n,
              //'.*': 'gauss%d' % n,
              //'.*': 'framefanout%d ' % n,
              '.*': 'orig%d' % n,
            },
            trace: {
              // fake doing Nmult SP pipelines
              //orig: ['wiener', 'gauss'],
              //'.*': 'orig',
            },
          }
          for n in fannums
        ],
      },
    }, nin=1, nout=fanmult),


    local fanin = g.pnode({
      type: fin,
      name: name,
      data: {
        multiplicity: fanmult,
        tag_rules: [
          {
            frame: {
              //['number%d' % n]: ['output%d' % n, 'output'],
              '.*': 'framefanin',
            },
            trace: {
              //gauss: 'gauss%d' % n,
              //wiener: 'wiener%d' % n,
              ['gauss%d' % n]: ['gauss%d' % n],
              ['wiener%d' % n]: ['wiener%d' % n],
              ['threshold%d' % n]: ['threshold%d' % n],
            },

          }
          for n in fannums
        ],

        //tags: if outtags == [] then ['from-pipeline-%d' % n for n in fannums] else outtags,
        //tags: ['from-pipeline-%d' % n for n in fannums],
        tags: outtags,
      },
    }, nin=fanmult, nout=1),

    ret: g.intern(innodes=[fanout],
                  outnodes=[fanin],
                  centernodes=pipelines,
                  edges=
                  [g.edge(fanout, pipelines[n], n, 0) for n in std.range(0, fanmult - 1)] +
                  [g.edge(pipelines[n], fanin, 0, n) for n in std.range(0, fanmult - 1)],
                  name=name),
  }.ret,


}
