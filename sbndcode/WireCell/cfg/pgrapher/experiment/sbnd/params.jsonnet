// ProtoDUNE-SP specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to SBND.

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

base {
    // This section will be overwritten in simparams.jsonnet
    det : {

        // The "faces" is consumed by, at least, the Drifter and
        // AnodePlane.  The "wires" number is used to set
        // AnodePlane.ident used to lookup the anode in WireSchema.
        // It corresponds to the anode number.
	// more positive X side is face 0 [this means the face.ident in geometry json.bz2 file]
	// in the below configuration of faces, if only one set of planes, 
	// just use the first element to configure the drift volume
        //
	// check wire geometry file to define these parameters
	local uplane_left = 201.45*wc.cm,  
	local wplane_left = 202.05*wc.cm,
	local uplane_right = -201.45*wc.cm, 
	local wplane_right = -202.05*wc.cm,
	local cpa_left = 0.45*wc.cm, // active volume min from GDML; x-positions of CPA surfaces 
	local cpa_right = -0.45*wc.cm, // active volume min from GDML; x-positions of CPA surfaces

        // The "response" plane is where the field response functions
        // start.  Garfield calcualtions start somewhere relative to
        // something, here's where that is made concrete.  This MUST
        // match what field response functions also used.
        response_plane: 10*wc.cm, // relative to collection wires
 	local res_plane = self.response_plane, // relative to the collection plane

        // The cathode plane is like the anode cut off plane.  Any
        // depo not between the two is dropped prior to drifting.
        //local cpa_plane = apa_cpa - 0.5*cpa_thick,


        // The volumes are then defined in terms of these above
        // numbers.  You can use "wirecell-util wires-info" or
        // "wirecell-util wires-volumes" or others to understand the
        // mapping of anode number to the 6 locations in X and Z.  For
        // Larsoft wires the numbering is column major starting at
        // small X and Z so the centerline is -/+/-/+/-/+.  Also
        // important is that the faces are listed "front" first.
        // Front is the one with the more positive X coordinates and
        // if we want to ignore a face it is made null.
        volumes: [
            {
                local sign = 2*(n%2)-1,
                //local centerline = sign*apa_cpa,
                wires: n,       // anode number
                name: "apa%d"%n,
                faces:
                // top, front face is against cryo wall
                if sign > 0
                then [ // tpc 1 face 1; left, west
                    {
                        anode: uplane_left, // x-position of induction plane
                        response: wplane_left - res_plane, // x-position of response plane
                        cathode: cpa_left,
                    },
                    null
                ]
                // bottom, back face is against cryo wall
                else [ // tpc 0 face 0; right, east
                    {
                        anode: uplane_right,
                        response: wplane_right + res_plane,
                        cathode: cpa_right,
                    },
                    null
                ],
            } for n in std.range(0,1)],

        // This describes some rough, overall bounding box.  It's not
        // directly needed but can be useful on the Jsonnet side, for
        // example when defining some simple kinematics.  It is
        // represented by a ray going from extreme corners of a
        // rectangular solid.  Again "wirecell-util wires-info" helps
        // to choose something.
        bounds : {
            tail: wc.point(-2.0, -2.0, 0.0, wc.m),
            head: wc.point(+2.0,  2.0, 5.0, wc.m),
        }
    },

    daq: super.daq {
        nticks: 3415,
    },

    adc: super.adc {
	// AC coupling in front of ADC to solve FE baseline distortion
	// induction plane: 6k Ohm vs 6k Ohm; collection: 6k Ohm vs 1.2k Ohm
	//Vref = 1.8 V
	// induce a RC filter  with tau = 1ms
        baselines: [879.5*wc.millivolt, 879.5*wc.millivolt, 286.0*wc.millivolt],

        // check this
        fullscale: [0*wc.volt, 1.8*wc.volt],
    },

    // This sets a relative gain at the input to the ADC.  Note, if
    // you are looking to fix SimDepoSource, you are in the wrong
    // place.  See the "scale" parameter of wcls.input.depos() defined
    // in pgrapher/common/ui/wcls/nodes.jsonnet.
    // also, see later overwriting in simparams.jsonnet
    elec: super.elec {
      postgain: 1.0, // pulser calibration: 41.649 ADC*tick/1ke
                       // theoretical elec resp (14mV/fC): 36.6475 ADC*tick/1ke
      shaping: 2.2 * wc.us,
    },

    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

        // The "absolute" time (ie, relative to trigger time?) that the lower edge
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = -200*wc.us,

        // Open the ductor's gate a bit early.
        local response_time_offset = $.det.response_plane / $.lar.drift_speed,
        local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),

        ductor : {
            nticks: $.daq.nticks + response_nticks,
            readout_time: self.nticks * $.daq.tick,
            start_time: tick0_time - response_time_offset,
        },

        // To counter the enlarged duration of the ductor, a Reframer
        // chops off the little early, extra time.  Note, tags depend on how 
        reframer: {
            tbin: response_nticks,
            nticks: $.daq.nticks,
        }
        
    },

    files: {
        wires: "sbnd-wires-geometry-v0202.json.bz2", // new SBND geometry

        fields: [ "garfield-sbnd-v1.json.bz2" ],

        // noise: "sbn_fd_incoherent_noise.json.bz2",
        noise: "sbnd-noise-spectra-v1.json.bz2", // Scaled from ProtoDUNE I measurement

        // coherent_noise: "sbn_fd_coherent_noise.json.bz2",

        chresp: null,
    },

}

