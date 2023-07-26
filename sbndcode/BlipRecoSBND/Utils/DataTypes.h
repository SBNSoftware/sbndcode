
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <vector>

typedef std::vector<int>        vint_t;
typedef std::vector<bool>       vbool_t;
typedef std::vector<float>      vfloat_t;
typedef std::set<int>           si_t;
typedef std::map<int,float>     mif_t;

const int kNplanes  = 3;  

namespace blip {
  
  //###################################################
  //  Data structures
  //###################################################

  struct ParticleInfo {
    simb::MCParticle particle;
    int   trackId           = -9;
    int   index             = -9;
    int   isPrimary         = -9;
    int   numTrajPts          = -9;
    double depEnergy         = -9;
    int   depElectrons        = -9;
    double numElectrons      = -9;
    double mass              = -9;
    double E                 = -9;
    double endE              = -9;
    double KE                = -9;
    double endKE             = -9;
    double P                 = -9; 
    double Px                = -9; 
    double Py                = -9; 
    double Pz                = -9; 
    double pathLength        = -9;
    double time              = -9;
    double endtime           = -9;
    TVector3 startPoint;      
    TVector3 endPoint;
    TVector3 position;
  };
  
  // True energy depositions
  struct TrueBlip {
    int       ID            = -9;     // unique blip ID
    int       TPC           = -9;     // TPC ID
    float     Time          = -999e9; // time [us]
    float     Energy        = 0;      // energy dep [MeV]
    int       DepElectrons  = 0;      // deposited electrons
    int       NumElectrons  = 0;      // electrons reaching wires
    float     DriftTime     = -9;     // drift time [us]
    int       LeadG4ID      = -9;     // lead G4 track ID
    int       LeadG4Index   = -9;     // lead G4 track index
    int       LeadG4PDG     = -9;     // lead G4 PDG
    float     LeadCharge    = -9;     // lead G4 charge dep
    TVector3  Position;               // XYZ position
    mif_t     G4ChargeMap;          
    mif_t     G4PDGMap;
  };

  struct HitInfo {
    int   hitid         = -9;
    int   plane         = -9;
    int   tpc           = -9;
    int   wire          = -9;
    int   chan          = -9;
    float amp           = -9;
    float rms           = -9;
    int   trkid         = -9;
    int   shwrid        = -9;
    int   clustid       = -9;
    int   blipid        = -9;
    bool  ismatch       = false;
    float integralADC    = -999;     // [ADCs] from integral
    float sigmaintegral = -999;
    float sumADC        = -999;     // [ADCs] from sum 
    float charge        = -999;     // [e-]
    float peakTime      = -999999;
    float driftTime     = -999999;  // [tick]
    float gof           = -9;
    int   g4trkid       = -9;
    int   g4pdg         = -999;
    int   g4charge      = -999;     // [e-]
    float g4frac        = -99;      
    float g4energy      = -999;     // [MeV]
  };
  
  struct HitClust {
    int     ID              = -9;
    bool    isValid         = false;
    int     CenterChan      = -999;
    int     CenterWire      = -999;
    bool    isMerged        = false;
    bool    isMatched       = false;
    int     DeadWireSep     = 99;
    int     TPC             = -9;
    int     Plane           = -9;
    int     NHits           = -9;
    int     NWires          = -9;
    float   ADCs            = -999;
    float   Amplitude       = -999;
    float   Charge          = -999;
    float   SigmaCharge     = -999;
    float   Time            = -999;
    float   RMS             = -999;
    float   StartHitTime    = -999;
    float   EndHitTime      = -999;
    float   StartTime       = -999;
    float   EndTime         = -999;
    float   Timespan        = -999;
    int     StartWire       = -999;
    int     EndWire         = -999;
    int     NPulseTrainHits = -9;
    float   GoodnessOfFit   = -999;
    int     BlipID          = -9;
    int     EdepID          = -9;
    si_t    HitIDs;
    si_t    Wires;
    si_t    Chans;
    si_t    G4IDs;
    
    std::map<int,TVector3> IntersectLocations;
  };

  struct Blip {
    
    int       ID              = -9;         // Blip ID / index
    bool      isValid         = false;      // Blip passes basic checks
    int       TPC             = -9;         // TPC
    int       NPlanes         = -9;         // Num. matched planes
    int       MaxWireSpan     = -9;         // Maximum span of wires on any plane cluster
    float     Charge          = -9;         // Charge on calorimetry plane
    float     Energy          = -999;       // Energy (const dE/dx, fcl-configurable)
    float     EnergyESTAR     = -999;       // Energy (ESTAR method from ArgoNeuT)
    float     Time            = -999;       // Drift time [ticks]
    float     ProxTrkDist     = -9;         // Distance to cloest track
    int       ProxTrkID       = -9;         // ID of closest track
    bool      inCylinder      = false;      // Is it in a cone/cylinder region? 
    
    TVector3  Position;                     // 3D position TVector3
    float     SigmaYZ         = -9.;        // Uncertainty in YZ intersect [cm]
    float     dX              = -9;         // Equivalent length along drift direction [cm] 
    float     dYZ             = -9;         // Approximate length scale in YZ space [cm]

    // Plane/cluster-specific information
    blip::HitClust clusters[kNplanes];
    
    // Truth-matched energy deposition
    blip::TrueBlip truth;
    
    // Prototype getter functions
    double X() { return Position.X(); }
    double Y() { return Position.Y(); }
    double Z() { return Position.Z(); }
    
  };
  
}


