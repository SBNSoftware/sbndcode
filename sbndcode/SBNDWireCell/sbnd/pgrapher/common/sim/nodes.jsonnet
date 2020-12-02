// This file provides variety of simulation related Pnodes
// parameterized on tools and params.


local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local depos = import "pgrapher/common/sim/depos.jsonnet";

function(params, tools)
{
    // Create a drifter Pnode.
    drifter: g.pnode({
        local xregions = wc.unique_list(std.flattenArrays([v.faces for v in params.det.volumes])),

        type: "Drifter",
        data: params.lar {
            rng: wc.tn(tools.random),
            xregions: xregions,
            time_offset: params.sim.depo_toffset,

            drift_speed: params.lar.drift_speed,
            fluctuate: params.sim.fluctuate,

            DL: params.lar.DL,
            DT: params.lar.DT,
            lifetime: params.lar.lifetime,

        },
    }, nin=1, nout=1, uses=[tools.random]),

    // Implement "fixed" depo mode like LArG4 uses
    make_bagger :: function(name="bagger") g.pnode({
        type:'DepoBagger',
        name:name,
        data: {
            gate: [params.sim.ductor.start_time,
                   params.sim.ductor.start_time+params.sim.ductor.readout_time],
        },
    }, nin=1, nout=1),


    // The set of all ductors are formed as the "cross product" of all
    // anodes and all PIR trios.  For one element of that product this
    // function is called.  The name should be unique across all
    // anodes X PIR trios.
    make_ductor:: function(name, anode, pir_trio, type='Ductor') g.pnode({
        type: type,
        name: name,
        data: {
            rng: wc.tn(tools.random),
            anode: wc.tn(anode),
            pirs: std.map(function(pir) wc.tn(pir), pir_trio),

            continuous: params.sim.continuous,
            fixed: params.sim.fixed,

            fluctuate: params.sim.fluctuate,
            drift_speed: params.lar.drift_speed,
            first_frame_number: params.daq.first_frame_number,
            readout_time: params.sim.ductor.readout_time,
            start_time: params.sim.ductor.start_time,
            tick: params.daq.tick,
            nsigma: 3,
        },
    }, nin=1,nout=1,uses=[tools.random, anode] + pir_trio),
    
    make_depotransform :: function(name, anode, pirs) g.pnode({
        type:'DepoTransform',
        name:name,
        data: {
            rng: wc.tn(tools.random),
            anode: wc.tn(anode),
            pirs: std.map(function(pir) wc.tn(pir), pirs),
            fluctuate: params.sim.fluctuate,
            drift_speed: params.lar.drift_speed,
            first_frame_number: params.daq.first_frame_number,
            readout_time: params.sim.ductor.readout_time,
            start_time: params.sim.ductor.start_time,
            tick: params.daq.tick,
            nsigma: 3,
        },
    }, nin=1, nout=1, uses=[anode, tools.random] + pirs),

    // This may look similar to above but above is expected to diverge
    make_depozipper :: function(name, anode, pirs) g.pnode({
        type:'DepoZipper',
        name:name,
        data: {
            rng: wc.tn(tools.random),
            anode: wc.tn(anode),
            pirs: std.map(function(pir) wc.tn(pir), pirs),
            fluctuate: params.sim.fluctuate,
            drift_speed: params.lar.drift_speed,
            first_frame_number: params.daq.first_frame_number,
            readout_time: params.sim.ductor.readout_time,
            start_time: params.sim.ductor.start_time,
            tick: params.daq.tick,
            nsigma: 3,
        },
    }, nin=1, nout=1, uses=[anode, tools.random] + pirs),

    make_reframer :: function(name, anode, tags=[]) g.pnode({
        type: 'Reframer',
        name: name,
        data: {
            anode: wc.tn(anode),
            tags: tags,
            fill: 0.0,
            tbin: params.sim.reframer.tbin,
            toffset: params.sim.reframer.toffset,
            nticks: params.sim.reframer.nticks,
        },
    }, nin=1, nout=1, uses=[anode]),


    // make all ductors for given anode and for all PIR trios.
    make_anode_ductors:: function(anode)
    std.mapWithIndex(function(n, pir_trio)
                     $.make_ductor('ductor%d%s'%[n, anode.name], anode, pir_trio), tools.pirs),


    // Multi APA's are harder
    
    // make all ductors for a given PIR trio.  Basename should include
    // an identifier unique to the PIR trio.
    make_detector_ductors:: function(pirname, anodes, pir_trio)
    std.mapWithIndex(function (n, anode)
                     $.make_ductor(pirname + anode.name,
                                   anode, pir_trio), tools.anodes),
        

    // Map above function across all trio of PIRs.  Result is a 2D
    // array of ductor Pnodes indexed like: [ipir][ianode]
    ductors: std.mapWithIndex(function (n, pir_trio)
                              $.make_detector_ductors("ductor%d"%[n], tools.anodes, pir_trio), tools.pirs),


    // Make a WireBoundedDepos
    make_wbdepo:: function(anode, regions, mode="accept", name="")
    g.pnode({
        type: 'WireBoundedDepos',
        name: name,
        data: {
            anode: wc.tn(anode),
            regions: regions,
            mode: mode,
        },
    }, nin=1, nout=1, uses = [anode]),

    // Make a multiductor as a subgraph.
    // pipes = [{wbdepos:..., ductor:..., tag:...}, ...]
    multi_ductor_graph:: function(anode, pipes, name="") {
        size: std.length(pipes),
        fanout: g.pnode({
            type: 'DepoFanout',
            name: name + 'fanout',
            data: {
                multiplicity: std.length(pipes),
            },
        }, nin=1, nout=self.size),
        pipelines: [g.pipeline([pipes[n].wbdepos, pipes[n].ductor], name=name+'pipe'+std.toString(n))
                    for n in std.range(0, self.size-1)],
        fanin: g.pnode({
            type: 'FrameFanin',
            name: name+'fanin',
            data: {
                multiplicity: std.length(pipes),
                tags: [p.tag for p in pipes],
            },
        }, nin=self.size, nout=1),
        reframer: g.pnode({
            type: 'Reframer',
            name: name+'reframer',
            data: {
                anode: wc.tn(anode),
                tags: [p.tag for p in pipes],
                fill: 0.0,
                tbin: params.sim.reframer.tbin,
                toffset: params.sim.reframer.toffset,
                nticks: params.sim.reframer.nticks,
            },
        }, nin=1, nout=1),

        outpipe: g.pipeline([self.fanin, self.reframer], name=name+'outpipe'),

        // the thing you want
        //graph: g.pipeline([self.fanout, self.fanin, self.reframer], name+'graph'),
        graph: g.intern(innodes=[self.fanout],
                        outnodes=[self.outpipe],
                        centernodes=self.pipelines,
                        edges=
                        [g.edge(self.fanout, self.pipelines[n], n, 0)
                         for n in std.range(0,self.size-1)] +
                        [g.edge(self.pipelines[n], self.fanin, 0, n)
                         for n in std.range(0,self.size-1)],
                        name=name+'graph'),
    }.graph,
           
        

    // Make aone multiductor for a single anode from the primitive
    // ductors which are also featured in the given chain.  The chain
    // is left as an exercise to the caller.
    multi_ductor:: function(anode, ductors, chains, name="") g.pnode({
        type: "MultiDuctor",
        data : {
            anode: wc.tn(anode),
            continuous: params.sim.continuous,
            chains : chains,
            tick: params.daq.tick,
            start_time : params.daq.start_time,
            readout_time: params.daq.readout_time, 
            first_frame_number: params.daq.first_frame_number,
        }
    }, nin=1, nout=1, uses = [anode] + ductors),


    // This operates on all channels so needs a channel selector bypass.
    // Maybe fixme: there are tag-aware nodes inside.
    misconfigure:: function(params, chndbobj=null) {

        local split = g.pnode({
            type: "FrameSplitter",
            name: "misconsplit"
        }, nin=1, nout=2),

        local chsel_static = g.pnode({
            type: "ChannelSelector",
            name: "misconsel_static",
            data: {
                channels: params.nf.misconfigured.channels,
            }
        }, nin=1, nout=1),

        local chsel_dynamic = g.pnode({
            type: "DBChannelSelector",
            name: "misconsel_dynamic",
            data: {
                channelDB: wc.tn(chndbobj), 
            }
        }, nin=1, nout=1, uses = [chndbobj]),

        local chsel = if std.type(chndbobj)=='null'
                        then chsel_static
                        else chsel_dynamic,

        local miscon = g.pnode({
            type: "Misconfigure",
            name: "sigmisconfig",
            data: {
                // Must match what was actually used to start with
                from: {
                    gain: params.elec.gain,
                    shaping: params.elec.shaping,
                },
                to: {
                    gain: params.nf.misconfigured.gain,
                    shaping: params.nf.misconfigured.shaping,
                },
                tick: params.daq.tick, // sample period of the response

                // fixme: these should probably be set from params.
                nsamples: 50,   // number of samples of the response
                truncate:true // result is extended by nsamples, tuncate clips that off
            }
        }, nin=1, nout=1),

        local merge = g.pnode({
            type: "FrameMerger",
            name: "misconmerge",
            data: {
                rule: "replace",
                // note: the first two need to match the order of what data is
                // fed to ports 0 and 1 of this component in the pgraph below!
                mergemap: [
                ],
            }
        }, nin=2, nout=1),

        return : g.intern([split], [merge], [chsel, miscon],
                          edges=[
                              g.edge(split, chsel),
                              g.edge(chsel, miscon),
                              g.edge(miscon, merge, 0, 0),
                              g.edge(split, merge, 1, 1),
                          ],
                          name="misconfigure"),
    }.return,

    // Make a digitizer bound to an anode.
    digitizer:: function(anode, name="", tag="") g.pnode({
        type: "Digitizer",
        name: name,
        data : params.adc {
            anode: wc.tn(anode),
            frame_tag: tag
        }
    }, nin=1, nout=1, uses=[anode]),

    // Sinks cap off the end of the graph.
    frame_sink: g.pnode({ type: "DumpFrames" }, nin=1, nout=0),

    // Note, this is "noisy" and only useful for debugging/developing.
    depo_sink: g.pnode({ type: "DumpDepos" }, nin=1, nout=0),


    // drifter + ductor + digitizer = signal
    // signal: g.intern([self.drifter],[self.multi_ductor],
    //                  edges=[
    //                      g.edge(self.drifter, self.multi_ductor),
    //                  ]),

    
} + depos(params,tools)
    
