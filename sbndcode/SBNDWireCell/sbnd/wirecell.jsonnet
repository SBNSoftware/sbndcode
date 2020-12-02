/**
   Jsonnet version of Wire Cell System of units.  This is a one-time
   manual/Emacs conversion from the WireCellUtil/Units.h file.

   Below the units is some helper code.
*/

{
    // constants
    pi: 2*std.acos(0),
    twopi  : 2*self.pi,
    halfpi  : self.pi/2,
    pi2 : self.pi*self.pi,
    sqrtpi: std.sqrt(self.pi),

    //
    // Length [L]
    //
    millimeter  : 1.0,
    millimeter2 : self.millimeter*self.millimeter,
    millimeter3 : self.millimeter*self.millimeter*self.millimeter,

    centimeter  : 10.0*self.millimeter,
    centimeter2 : self.centimeter*self.centimeter,
    centimeter3 : self.centimeter*self.centimeter*self.centimeter,

    meter  : 1000.0*self.millimeter,
    meter2 : self.meter*self.meter,
    meter3 : self.meter*self.meter*self.meter,

    kilometer : 1000.0*self.meter,
    kilometer2 : self.kilometer*self.kilometer,
    kilometer3 : self.kilometer*self.kilometer*self.kilometer,

    parsec : 3.0856775807e+16*self.meter,

    micrometer : 1.0e-6 *self.meter,
    nanometer : 1.0e-9 *self.meter,
    angstrom  : 1.0e-10*self.meter,
    fermi     : 1.0e-15*self.meter,

    barn : 1.0e-28*self.meter2,
    millibarn : 1.0e-3 *self.barn,
    microbarn : 1.0e-6 *self.barn,
    nanobarn : 1.0e-9 *self.barn,
    picobarn : 1.0e-12*self.barn,

    // symbols
    nm  : self.nanometer,
    um  : self.micrometer,

    mm  : self.millimeter,
    mm2 : self.millimeter2,
    mm3 : self.millimeter3,

    cm  : self.centimeter,
    cm2 : self.centimeter2,
    cm3 : self.centimeter3,

    m  : self.meter,
    m2 : self.meter2,
    m3 : self.meter3,

    km  : self.kilometer,
    km2 : self.kilometer2,
    km3 : self.kilometer3,

    pc : self.parsec,

    //
    // Angle
    //
    radian      : 1.0,
    milliradian : 1.0e-3*self.radian,
    degree : (self.pi/180.0)*self.radian,

    steradian : 1.0,

    // symbols
    rad  : self.radian,
    mrad : self.milliradian,
    sr   : self.steradian,
    deg  : self.degree,

    //
    // Time [T]
    //
    nanosecond  : 1.0,
    second      : 1.0e+9 *self.nanosecond,
    millisecond : 1.0e-3 *self.second,
    microsecond : 1.0e-6 *self.second,
    picosecond : 1.0e-12*self.second,

    hertz : 1.0/self.second,
    kilohertz : 1.0e+3*self.hertz,
    megahertz : 1.0e+6*self.hertz,

    // symbols
    ns : self.nanosecond,
    s : self.second,
    ms : self.millisecond,
    us : self.microsecond,

    //
    // Electric charge [Q]
    //
    eplus : 1.0,// positron charge
    e_SI  : 1.602176487e-19,// positron charge in coulomb
    coulomb : self.eplus/self.e_SI,// coulomb : 6.24150 e+18 * eplus
    fC : 1.0e-15*self.coulomb,     // femtocoulomb

    //
    // Energy [E]
    //
    megaelectronvolt : 1.0,
    electronvolt : 1.0e-6*self.megaelectronvolt,
    kiloelectronvolt : 1.0e-3*self.megaelectronvolt,
    gigaelectronvolt : 1.0e+3*self.megaelectronvolt,
    teraelectronvolt : 1.0e+6*self.megaelectronvolt,
    petaelectronvolt : 1.0e+9*self.megaelectronvolt,

    joule : self.electronvolt/self.e_SI,// joule : 6.24150 e+12 * MeV

    // symbols
    MeV : self.megaelectronvolt,
    eV : self.electronvolt,
    keV : self.kiloelectronvolt,
    GeV : self.gigaelectronvolt,
    TeV : self.teraelectronvolt,
    PeV : self.petaelectronvolt,

    //
    // Mass [E][T^2][L^-2]
    //
    kilogram : self.joule*self.second*self.second/(self.meter*self.meter),
    gram : 1.0e-3*self.kilogram,
    milligram : 1.0e-3*self.gram,

    // symbols
    kg : self.kilogram,
    g : self.gram,
    mg : self.milligram,

    //
    // Power [E][T^-1]
    //
    watt : self.joule/self.second,// watt : 6.24150 e+3 * MeV/ns

    //
    // Force [E][L^-1]
    //
    newton : self.joule/self.meter,// newton : 6.24150 e+9 * MeV/mm

    //
    // Pressure [E][L^-3]
    //
    pascal : self.newton/self.m2,   // pascal : 6.24150 e+3 * MeV/mm3
    bar        : 100000*self.pascal, // bar    : 6.24150 e+8 * MeV/mm3
    atmosphere : 101325*self.pascal, // atm    : 6.32420 e+8 * MeV/mm3

    //
    // Electric current [Q][T^-1]
    //
    ampere : self.coulomb/self.second, // ampere : 6.24150 e+9 * eplus/ns
    milliampere : 1.0e-3*self.ampere,
    microampere : 1.0e-6*self.ampere,
    nanoampere : 1.0e-9*self.ampere,

    //
    // Electric potential [E][Q^-1]
    //
    megavolt : self.megaelectronvolt/self.eplus,
    kilovolt : 1.0e-3*self.megavolt,
    volt : 1.0e-6*self.megavolt,
    millivolt : 1.0e-3*self.volt,
    microvolt : 1.0e-6*self.volt,
    mV : self.millivolt,
    uV : self.microvolt,

    //
    // Electric resistance [E][T][Q^-2]
    //
    ohm : self.volt/self.ampere,// ohm : 1.60217e-16*(MeV/eplus)/(eplus/ns)

    //
    // Electric capacitance [Q^2][E^-1]
    //
    farad : self.coulomb/self.volt,// farad : 6.24150e+24 * eplus/Megavolt
    millifarad : 1.0e-3*self.farad,
    microfarad : 1.0e-6*self.farad,
    nanofarad : 1.0e-9*self.farad,
    picofarad : 1.0e-12*self.farad,

    //
    // Magnetic Flux [T][E][Q^-1]
    //
    weber : self.volt*self.second,// weber : 1000*megavolt*ns

    //
    // Magnetic Field [T][E][Q^-1][L^-2]
    //
    tesla     : self.volt*self.second/self.meter2,// tesla :0.001*megavolt*ns/mm2

    gauss     : 1.0e-4*self.tesla,
    kilogauss : 1.0e-1*self.tesla,

    //
    // Inductance [T^2][E][Q^-2]
    //
    henry : self.weber/self.ampere,// henry : 1.60217e-7*MeV*(ns/eplus)**2

    //
    // Temperature
    //
    kelvin : 1.0,

    //
    // Amount of substance
    //
    mole : 1.0,

    //
    // Activity [T^-1]
    //
    becquerel : 1.0/self.second ,
    curie : 3.7e+10 * self.becquerel,
    kilobecquerel : 1.0e+3*self.becquerel,
    megabecquerel : 1.0e+6*self.becquerel,
    gigabecquerel : 1.0e+9*self.becquerel,
    millicurie : 1.0e-3*self.curie,
    microcurie : 1.0e-6*self.curie,
    Bq : self.becquerel,
    kBq : self.kilobecquerel,
    MBq : self.megabecquerel,
    GBq : self.gigabecquerel,
    Ci : self.curie,
    mCi : self.millicurie,
    uCi : self.microcurie,

    //
    // Absorbed dose [L^2][T^-2]
    //
    gray : self.joule/self.kilogram ,
    kilogray : 1.0e+3*self.gray,
    milligray : 1.0e-3*self.gray,
    microgray : 1.0e-6*self.gray,

    //
    // Luminous intensity [I]
    //
    candela : 1.0,

    //
    // Luminous flux [I]
    //
    lumen : self.candela*self.steradian,

    //
    // Illuminance [I][L^-2]
    //
    lux : self.lumen/self.meter2,

    //
    // Miscellaneous
    //
    perCent     : 0.01 ,
    perThousand : 0.001,
    perMillion  : 0.000001,


    // some constants of nature.
    clight : 2.99792458e8*self.meter/self.second,

    //// Above are Wire Cell system of units.



    //// Below are some Jsonnet helpers


    // at 500 volts
    nominal_drift_velocity: 1.6*self.mm/self.us,

    // vectors
    point(x,y,z,u) :: {x:x*u, y:y*u, z:z*u},
    ray(p1,p2) :: {tail:p1, head:p2},

    Point :: {x:0,y:0,z:0},
    Ray :: {tail:self.Point,head:self.Point},
    Track :: { time:0.0, charge:-1, ray:self.Ray },

    // WirePlaneID is a packed integer.  WARNING, layer is NOT what
    // most people call "plane number".  It is a bit field.  For
    // 3-plane detectors the outer most wire plane layer is 1, then 2
    // and collection is 4 (not 3).  layer=0 is undefined.
    Ulayer:1<<0,
    Vlayer:1<<1,
    Wlayer:1<<2,
    WirePlaneId(layer, face=0, apa=0) :: (layer&7) | (face << 3) | (apa << 4),

    // Base class for a configurable.
    Component :: {
	type:"",
	name:"",
	data:{}
    },
    /// example usage: 
    TrackDepos :: self.Component + { type: "TrackDepos" },

    /// Return canonical "type:name" or just "type" if no name from a
    /// configuration object.  Use this instead of bare names to
    /// better guard against typos and changes in dependent
    /// configuration.  So instead of using:
    ///
    ///    anode: "Anode:myanode",
    ///
    /// use the more robust:
    ///
    ///    anode: wc.tn(myanode),
    ///
    tn(obj) :: if std.objectHas(obj, "name") && obj.name != "" then obj.type + ":" + obj.name else obj.type,


    // Return a new list where only the first occurrence of any object is kept.
    unique_helper(l, x):: if std.count(l,x) == 0 then l + [x] else l,
    unique_list(l):: std.foldl($.unique_helper, l, []),



    // Round a floating point to nearest integer.  It's a bit weird to
    // go through a format/parse.  Maybe there's a better way?
    roundToInt(x):: std.parseInt("%d" % (x+0.5)),

    freqbinner :: function(tick, nsamples) {
        nyquist : 0.5 / tick,
        hz_perbin : 1.0/(tick/$.second * nsamples),

        // A function to return the frequency bin holding the given
        // frequency.  The frequency is specified in the system of
        // units.
        bin :: function(frequency) std.floor((frequency/$.hertz) / self.hz_perbin),

        // Return a frequency bin range holding meanfreq +/- deltafreq,
        // both freqencies in system of units.
        bin_range :: function(meanfreq, deltafreq)  [
            self.bin(std.max(0, meanfreq-deltafreq)),
            self.bin(std.min(self.nyquist, meanfreq+deltafreq))
        ],
    
        
        // Return something suitable to set to chndb's
        // channel_info[].freqmasks.  The "meanfreqs" should be a list
        // of frequencies in the sytem of units and delta is a common
        // detla.  See bin_range.
        local _br = self.bin_range,
        freqmasks :: function(meanfreqs, delta) [
            { value: 1.0, lobin: 0, hibin: nsamples-1 }
        ] + [ {
            local br = _br(mf, delta),
            value: 0.0, lobin: br[0], hibin: br[1],
        } for mf in meanfreqs],

        testfreqs :: [f*$.kilohertz for f in [51.5, 102.8, 154.2, 205.5, 256.8, 308.2, 359.2, 410.5, 461.8, 513.2, 564.5, 615.8]],
        testmasks : self.freqmasks(self.testfreqs, 2*$.kilohertz),
        
    },

}


