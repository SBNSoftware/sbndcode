local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local sim_maker = import "pgrapher/common/sim/nodes.jsonnet";


// return some nodes, includes base sim nodes.
function(params, tools) {
    local sim = sim_maker(params, tools),

    local nanodes = std.length(tools.anodes),

    // I rue the day that we must have an (anode) X (field) cross product!
    local ductors = sim.make_detector_ductors("nominal", tools.anodes, tools.pirs[0]),


    local zippers = [sim.make_depozipper("depozipper-"+tools.anodes[n].name, tools.anodes[n], tools.pirs[0])
                     for n in std.range(0, nanodes-1)],
    // One pRNG per anode.  The anode pipelines run concurrently under TbbFlow, so a pRNG shared between them
    // hands out its sequence in an order that depends on thread timing and the job is not reproducible.
    local rngs = [{
        type: "Random",
        name: "rng-" + tools.anodes[n].name,
        data: {
            generator: "default",
            seeds: [0, 1, 2, 3, 4, n + 1],
        },
    } for n in std.range(0, nanodes-1)],

    // As sim.make_depotransform, with the anode's own pRNG.
    local transforms = [g.pnode({
        type: 'DepoTransform',
        name: "depotransform-"+tools.anodes[n].name,
        data: {
            rng: wc.tn(rngs[n]),
            dft: wc.tn(tools.dft),
            anode: wc.tn(tools.anodes[n]),
            pirs: std.map(function(pir) wc.tn(pir), tools.pirs[0]),
            fluctuate: params.sim.fluctuate,
            drift_speed: params.lar.drift_speed,
            first_frame_number: params.daq.first_frame_number,
            readout_time: params.sim.ductor.readout_time,
            start_time: params.sim.ductor.start_time,
            tick: params.daq.tick,
            nsigma: 3,
        },
    }, nin=1, nout=1, uses=[tools.anodes[n], rngs[n], tools.dft] + tools.pirs[0])
                        for n in std.range(0, nanodes-1)],
    local depos2traces = transforms,
    //local depos2traces = zippers,

    local digitizers = [
        sim.digitizer(tools.anodes[n], name="digitizer-" + tools.anodes[n].name, tag="orig%d"%n)
        for n in std.range(0,nanodes-1)],

    local reframers = [
        g.pnode({
            type: 'Reframer',
            name: 'reframer-'+tools.anodes[n].name,
            data: {
                anode: wc.tn(tools.anodes[n]),
                tags: [],           // ?? what do?
                fill: 0.0,
                tbin: params.sim.reframer.tbin,
                toffset: 0,
                nticks: params.sim.reframer.nticks,
            },
        }, nin=1, nout=1) for n in std.range(0, nanodes-1)],


    // fixme: see https://github.com/WireCell/wire-cell-gen/issues/29
    local make_noise_model = function(anode, csdb=null) {
        type: "EmpiricalNoiseModel",
        name: "empericalnoise-" + anode.name,
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(tools.dft),
            chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
            spectra_file: params.files.noise,
            nsamples: params.daq.nticks,
            period: params.daq.tick,
            wire_length_scale: 1.0*wc.cm, // optimization binning
        },
        uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
    },
    local noise_models = [make_noise_model(anode) for anode in tools.anodes],


    local add_noise = function(model, rng) g.pnode({
        type: "AddNoise",
        name: "addnoise-" + model.name,
        data: {
            rng: wc.tn(rng),
            dft: wc.tn(tools.dft),
            model: wc.tn(model),
	    nsamples: params.daq.nticks,
            replacement_percentage: 0.02, // random optimization
        }}, nin=1, nout=1, uses=[rng, tools.dft, model]),

    local noises = [add_noise(noise_models[n], rngs[n]) for n in std.range(0, nanodes-1)],
    
    local outtags = ["orig%d"%n for n in std.range(0, nanodes-1)],

    ret : {

        signal_pipelines: [g.pipeline([depos2traces[n], reframers[n],  digitizers[n]],
                                      name="simsigpipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],

        splusn_pipelines:  [g.pipeline([depos2traces[n], reframers[n], noises[n], digitizers[n]],
                                       name="simsignoipipe-" + tools.anodes[n].name) for n in std.range(0, nanodes-1)],
    
        signal: f.fanpipe('DepoSetFanout', self.signal_pipelines, 'FrameFanin', "simsignalgraph", outtags),
        splusn: f.fanpipe('DepoSetFanout', self.splusn_pipelines, 'FrameFanin', "simsplusngraph", outtags),

    } + sim,                    // tack on base for user sugar.
}.ret
