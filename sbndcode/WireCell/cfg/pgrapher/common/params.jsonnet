// This file is part of wire-cell-toolkit/cfg/.
//
// This file provides a base data structure to define parameters that
// span all currently supported WCT functionality.  Not every
// parameter will be used and not every value here may be valid for
// your use and should be overridden.  The parameters are named and
// factored into sub-objects in order to be sympathetic to how the C++
// components are structured and name their configuration paramters.
// As such it's often possible to build a component configuration
// object by inheriting from one or more sub-objects in the parameter
// structure.  For most jobs, this structure should be derived and
// overriden before being passed to functions that produce other
// configuration structures.
//

local wc = import "wirecell.jsonnet";

{

    // Parameters relevant to the bulk liquid argon volume.
    lar : {
        // Longitudinal diffusion constant
        DL :  7.2 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 12.0 * wc.cm2/wc.s,
        // Electron lifetime
        lifetime : 8*wc.ms,
        // Electron drift speed, assumes a certain applied E-field
        drift_speed : 1.6*wc.mm/wc.us, // at 500 V/cm
        // LAr density
        density: 1.389*wc.g/wc.centimeter3,
        // Decay rate per mass for natural Ar39.
        ar39activity: 1*wc.Bq/wc.kg,
    },

    det: {

        // The detector volumes are defined as a set of planes
        // organized into a "front" and "back" face of an AnodePlane
        // and also used by the Drifter.  Each volume defined will map
        // to an AnodePlane which is given the same name.  The volumes
        // given here are several unrelated volumes and given just as
        // examples to show the structure.
        volumes : [
            {
                wires: 0,
                name: "onesided",
                faces: [ {anode: 0, response: 10*wc.cm, cathode: 2*wc.m},
                         null ],
            },
            {
                wires: 0,
                name: "twosided",
                faces: [ {anode: 0, response: +10*wc.cm, cathode: 2*wc.m},
                         {anode: 0, response: -10*wc.cm, cathode: -2*wc.m} ],
            },
            {
                wires: 0,
                name: "anothertwosided",
                faces: [ {anode: 4*wc.m, response: +10*wc.cm + 4*wc.m, cathode:  2*wc.m + 4*wc.m },
                         {anode: 4*wc.m, response: -10*wc.cm + 4*wc.m, cathode: -2*wc.m + 4*wc.m} ],
            },
        ],
    },
    
    // Parameters related to the DAQ
    daq : {
        // One digitization sampling period
        tick: 0.5*wc.us,

        // Number of ticks in one DAQ readout.  Note, some components
        // take an "nsamples".  This can be but need not be the same
        // as "nticks".  For example, NF will typicall differ.  Also,
        // in general this is not the number used for the Ductor for
        // simulation.
        nticks: 3415,

        // Readout period in units of time
        readout_time: self.tick*self.nticks,

        // In cases where a node limits its running based on number of
        // readouts (aka frames), this is how many are to be
        // processed.
        nreadouts: 1,

        // Where a node sets a readout (frame) time, this is the
        // starting time.
        start_time: 0.0*wc.s,

        // In case where a node limits running based on time, this is
        // when it will stop.
        stop_time: self.start_time + self.nreadouts*self.readout_time,

        // In a case where a node counts readouts (frames), this is
        // the first number.
        first_frame_number: 100,
    },

    // Parameters having to do with digitization.
    adc : {
        // A relative gain applied just prior to digitization.  This
        // is not FE gain, see elec for that.
        gain: 1.0,

        // Voltage baselines added to any input voltage signal listed
        // in a per plan (U,V,W) array.
        baselines: [900*wc.millivolt,900*wc.millivolt,200*wc.millivolt],

        // The resolution (bits) of the ADC
        resolution: 12,

        // The voltage range as [min,max] of the ADC, eg min voltage
        // counts 0 ADC, max counts 2^resolution-1.
        fullscale: [0*wc.volt, 2.0*wc.volt],
    },

    // Parameters having to do with the front end electronics
    elec : {
        // The FE amplifier gain in units of Voltage/Charge.
        gain : 14.0*wc.mV/wc.fC,

        // The shaping (aka peaking) time of the amplifier shaper.
        shaping : 2.0*wc.us,

        // An realtive gain applied after shaping and before ADC.
        // Don't use this to fix wrong sign depos.  If you neeed to
        // fix sign of larsoft depos, use:
        // wcls.input.depos([...], scale=-1)
        // or better, fix the C++ that sets the wrong sign to begin with.
        postgain: 1.0,

        fields : {
            
            // The start of the field response paths in X, measured from
            // the collection wires.  This is a positive, relative
            // distance.  It is set by the configuration that went into
            // the field response calcualtion (usually Garfield).  
            start_dx: 10*wc.cm,

            // The time it takes for "a" point deposition to drift
            // from teh reponse plane to a collection wire.  This
            // actually can vary over ~5us as a function of transverse
            // starting point (impact position) so it's only
            // approximated assuming nominal drift speed.  The nominal
            // 10cm and 1.1 mm/us gives 90.9us.  1.6mm/us is 62.5us.
            drift_dt: self.start_dx / $.lar.drift_speed,

            // The number of ticks that go by during the drift time
            // from the response plane to the collection plane.  Eg,
            // 182 and 125 for the example speeds above.
            nticks: wc.roundToInt(self.drift_dt / $.daq.tick),
        },

    },

    // Parameters related to simulation, not given elsewhere.
    sim : {

        // The number of impact bins per wire region gives the
        // granularity of the simulation convolution in the transverse
        // dimension.  Typically should match what the granularity at
        // which the field response functions are defined.
        nimpacts: 10,

        // if statistical fluctations should be applied
        fluctuate: true,
        // if continuous or discontinuous mode is is used.  See, eg
        // https://wirecell.github.io/news/posts/simulation-updates/
        continuous: true,

        // Fixed overrides continuous and simply makes frames at a fixed time.
        fixed: false,

        // A fixed time offset added to all drifted depo times which
        // can be useful if the origin depo source fails to provide
        // correct times.  
        depo_toffset: 0.0,


        // Default ductor parameters.  If your detector (eg MB) has
        // reason for a readout to start earlier, better override
        // this.
        ductor : {
            nticks: $.daq.nticks + $.elec.fields.nticks,
            readout_time: self.nticks * $.daq.tick,
            start_time: $.daq.start_time - $.elec.fields.drift_dt,
        },

        // If a ductor's time acceptance is increased then a Reframer
        // can be used to chop off the early excess to meet readout
        // assumptions.  Depending on the form of the ductor, the
        // reframer will likely need it's "tags" configured.
        reframer: {
            tbin: $.elec.fields.nticks,
            nticks: $.daq.nticks,
        }
    },

    // Parameters related to noise filtering.  
    nf : {                    

        // Number of frequency bins over which NF filters are applied.
        // Note, this likely differs from nticks so should be
        // overriden in the experiment specific parameters.  Where
        // they differ then truncation or extension of waveforms may
        // occur.  Although we set it to the default expected
        // daq.nticks, it may differ from number of ticks in the input
        // waveforms.
        nsamples: $.daq.nticks, 
        
    },    

    // Some configuration is too bulky to include directly and so are
    // expelicitly loaded by WCT from other files (typically as
    // compressed JSON).  The user is free to provide this files
    // themselves and there are wirecell-* Python programs that can
    // help generate them or convert from foreign formats.  However,
    // most common cases are provided for by files in the
    // wire-cell-data package.  Here, the attributes are given and a
    // more specfic parameter file should override and provide their
    // names.
    files : {
        
        // The "wire schema" file giving wire locations.
        // wirecell-util has generation and conersion commands to make
        // these.
        wires: null,

        // An array of field response files.  This array may span
        // alternative "universes" (such as "shorted wire regsion")
        // such as implemented in the ~MultiDuctor~ simulation
        // component.  The first field file is considered "nominal".
        // The info in these files are typically produced by Garfield
        // initially and then converted to WCT format using the
        // "wirecell-sigproc convert-garfield" command.  NOTE: you
        // must assure that elec.response.plane_dx is consistent with
        // what was used to generate the field files.
        fields: [],

        // A noise file provides a spectral lookup table.  This info
        // is usually provided by some analysis and converted.  One
        // such converter is "wirecell-sigproc convert-noise-spectra".
        noise: null,

        // This file gives per-channel calibrated responses.  See
        // "wirecell-sigproc channel-response" for one converter.
        chresp: null, 

        // The sigproc DNNROI model file. 
        dnnroi: "unet-l23-cosmic500-e50.ts",
    },
}


