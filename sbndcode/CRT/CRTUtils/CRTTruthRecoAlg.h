#ifndef CRTTRUTHRECOALG_H_SEEN
#define CRTTRUTHRECOALG_H_SEEN


///////////////////////////////////////////////
// CRTTruthRecoAlg.h
//
// Functions for truth matching in CRT reconstruction
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include "TVector3.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TPolyLine3D.h"


namespace sbnd{

  class CRTTruthRecoAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
    };

    CRTTruthRecoAlg(const Config& config);

    CRTTruthRecoAlg();

    ~CRTTruthRecoAlg();

    void reconfigure(const Config& config);

    bool CrossesTagger(const simb::MCParticle& particle, int tag_i);

    bool CrossesStrip(const simb::MCParticle& particle, int tag_i);

    TVector3 TaggerCrossPoint(const simb::MCParticle& particle, int tag_i);

    bool IsThroughGoing(const simb::MCParticle& particle);

    void DrawCube(TCanvas *c1, double *rmin, double *rmax, int col);
    
    // Function to calculate the CRT crossing points of a true particle
    std::pair<TVector3, TVector3> TpcCrossPoints(simb::MCParticle const& particle);

    double TpcLength(simb::MCParticle const& particle);

    // Convert start time to CRT crossing point
    TVector3 T0ToXYZPosition(TVector3 position, TVector3 direction, std::string tagger, int tpc, double t0);

  private:

    detinfo::DetectorProperties const* fDetectorProperties;
    art::ServiceHandle<geo::AuxDetGeometry> fAuxDetGeoService;
    const geo::AuxDetGeometry* fAuxDetGeo;
    const geo::AuxDetGeometryCore* fAuxDetGeoCore;

    TPCGeoAlg fTpcGeo;

    // Positions of the CRT planes
    std::vector<double> crtPlanes = {-359.1, -357.3, 357.3, 359.1, -358.9, -357.1, 661.52, 663.32, 865.52, 867.32, -240.65, -238.85, 655.35, 657.15};
    std::vector<int> fixCoord   = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2}; // Fixed coordinate for each plane
    std::vector<int> widthCoord = {2, 1, 2, 1, 0, 2, 2, 0, 2, 0, 1, 0, 1, 0}; // Width direction for each plane
    std::vector<int> lenCoord   = {1, 2, 1, 2, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1}; // Length direction for each plane

    std::map<std::string, int> nameToInd;
    

  };

}

#endif
