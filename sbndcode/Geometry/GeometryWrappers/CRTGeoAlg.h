#ifndef CRTGEOALG_H_SEEN
#define CRTGEOALG_H_SEEN


///////////////////////////////////////////////
// CRTGeoAlg.h
//
// Wrapper for some awkward CRT geometry things
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
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// c++
#include <vector>

// ROOT
#include "TVector3.h"
#include "TGeoManager.h"

namespace sbnd{

  struct CRTSiPMGeo{
    CRTSiPMGeo(const std::string &_stripName, const uint32_t _channel, const double location[3])
    {
      stripName = _stripName;
      channel   = _channel;
      x         = location[0];
      y         = location[1];
      z         = location[2];

      null = false;
    }
    std::string stripName;
    uint16_t    channel;
    double      x;
    double      y;
    double      z;
    bool        null;
  };

  // CRT strip geometry struct contains dimensions and mother module
  struct CRTStripGeo{
    CRTStripGeo(const TGeoNode *stripNode, const geo::AuxDetSensitiveGeo &auxDetSensitive, 
		const uint16_t _adsID, const std::string &_moduleName,
		const uint16_t _channel0, const uint16_t _channel1)
    {
      name       = stripNode->GetName();
      moduleName = _moduleName;
      channel0   = _channel0;
      channel1   = _channel1;

      // Strip Dimensions
      double halfWidth  = auxDetSensitive.HalfWidth1();
      double halfHeight = auxDetSensitive.HalfHeight();
      double halfLength = auxDetSensitive.Length()/2.;

      // Find world coordinates for edges
      double limits[3] = {halfWidth, halfHeight, halfLength};
      double limitsWorld[3];
      auxDetSensitive.LocalToWorld(limits, limitsWorld);

      double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
      double limitsWorld2[3];
      auxDetSensitive.LocalToWorld(limits2, limitsWorld2);

      // Fill edges & width
      minX  = std::min(limitsWorld[0], limitsWorld2[0]);
      maxX  = std::max(limitsWorld[0], limitsWorld2[0]);
      minY  = std::min(limitsWorld[1], limitsWorld2[1]);
      maxY  = std::max(limitsWorld[1], limitsWorld2[1]);
      minZ  = std::min(limitsWorld[2], limitsWorld2[2]);
      maxZ  = std::max(limitsWorld[2], limitsWorld2[2]);
      width = halfHeight * 2.;
      
      adsID = _adsID;
      null  = false;
    }
    std::string name;
    std::string moduleName;
    uint16_t    channel0;
    uint16_t    channel1;
    double      minX;
    double      maxX;
    double      minY;
    double      maxY;
    double      minZ;
    double      maxZ;
    double      width;
    uint16_t    adsID;
    bool        null;
  };

  // CRT module geometry struct contains dimensions, daughter strips and mother tagger
  struct CRTModuleGeo{
    CRTModuleGeo(const TGeoNode *moduleNode, const geo::AuxDetGeo &auxDet, 
		 const uint16_t _adID, const std::string &_taggerName)
    {
      name       = moduleNode->GetName();
      taggerName = _taggerName;

      // Module Dimensions
      double halfWidth  = auxDet.HalfWidth1();
      double halfHeight = auxDet.HalfHeight();
      double halfLength = auxDet.Length()/2;

      // Find world coordinates for edges
      double limits[3] = {halfWidth, halfHeight, halfLength};
      double limitsWorld[3];
      auxDet.LocalToWorld(limits, limitsWorld);

      double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
      double limitsWorld2[3];
      auxDet.LocalToWorld(limits2, limitsWorld2);

      // Which plane within the tagger
      double origin[3] = {0, 0, 0};
      double modulePosMother[3];
      moduleNode->LocalToMaster(origin, modulePosMother);
      planeID = (modulePosMother[2] > 0);

      // Location of SiPMs - CRT BeamTelescope specific
      top = false;

      // Fill edges
      minX = std::min(limitsWorld[0], limitsWorld2[0]);
      maxX = std::max(limitsWorld[0], limitsWorld2[0]);
      minY = std::min(limitsWorld[1], limitsWorld2[1]);
      maxY = std::max(limitsWorld[1], limitsWorld2[1]);
      minZ = std::min(limitsWorld[2], limitsWorld2[2]);
      maxZ = std::max(limitsWorld[2], limitsWorld2[2]);

      adID = _adID;
      null = false;
    }
    std::string   name;
    std::string   taggerName;
    double        minX;
    double        maxX;
    double        minY;
    double        maxY;
    double        minZ;
    double        maxZ;
    size_t        planeID;
    bool          top;
    uint16_t      adID;
    bool          null;
  };

  // CRT tagger geometry struct contains dimensions and daughter modules
  struct CRTTaggerGeo{
    CRTTaggerGeo(const TGeoNode *taggerNode, const TGeoNode *detNode)
    {
      // Fill name
      name = taggerNode->GetName();

      // Tagger Dimensions
      double halfWidth  = ((TGeoBBox*)taggerNode->GetVolume()->GetShape())->GetDX();
      double halfHeight = ((TGeoBBox*)taggerNode->GetVolume()->GetShape())->GetDY();
      double halfLength = ((TGeoBBox*)taggerNode->GetVolume()->GetShape())->GetDZ()/2;

      // Find world coordinates for edges
      double limits[3] = {halfWidth, halfHeight, halfLength};
      double limitsDet[3];
      taggerNode->LocalToMaster(limits, limitsDet);
      double limitsWorld[3];
      detNode->LocalToMaster(limitsDet, limitsWorld);

      double limits2[3] = {-halfWidth, -halfHeight, -halfLength};
      double limitsDet2[3];
      taggerNode->LocalToMaster(limits2, limitsDet2);
      double limitsWorld2[3];
      detNode->LocalToMaster(limitsDet2, limitsWorld2);

      // Fill edges
      minX = std::min(limitsWorld[0], limitsWorld2[0]);
      maxX = std::max(limitsWorld[0], limitsWorld2[0]);
      minY = std::min(limitsWorld[1], limitsWorld2[1]);
      maxY = std::max(limitsWorld[1], limitsWorld2[1]);
      minZ = std::min(limitsWorld[2], limitsWorld2[2]);
      maxZ = std::max(limitsWorld[2], limitsWorld2[2]);

      null = false;
    }
    std::string name;
    double      minX;
    double      maxX;
    double      minY;
    double      maxY;
    double      minZ;
    double      maxZ;
    bool        null;
  };


  class CRTGeoAlg {
  public:

    CRTGeoAlg(geo::GeometryCore const *geometry, geo::AuxDetGeometryCore const *auxdet_geometry);
    CRTGeoAlg();

    ~CRTGeoAlg();

    std::vector<double> CRTLimits() const;

    size_t NumTaggers() const;

    size_t NumModules() const;

    size_t NumStrips() const;

    CRTTaggerGeo GetTagger(const std::string taggerName) const;

    CRTModuleGeo GetModule(const std::string moduleName) const;

    CRTStripGeo GetStrip(const std::string stripName) const;

    std::string GetTaggerName(const std::string name) const;

    std::string ChannelToStripName(const uint16_t channel) const;

    std::string ChannelToTaggerName(const uint16_t channel) const;

    size_t ChannelToPlaneID(const uint16_t channel) const;

    std::vector<double> StripLimitsWithChargeSharing(const std::string stripName, const double x, 
						     const double ex);

    geo::Point_t ChannelToSipmPosition(const uint16_t channel) const;

    std::pair<int, int> GetStripSipmChannels(const std::string stripName) const;

    double DistanceDownStrip(const geo::Point_t position, const std::string stripName) const;

    bool CheckOverlap(const CRTModuleGeo &module1, const CRTModuleGeo &module2);

    bool HasOverlap(const CRTModuleGeo &module);

    bool StripHasOverlap(const std::string stripName);

  private:

    std::map<std::string, CRTTaggerGeo> fTaggers;
    std::map<std::string, CRTModuleGeo> fModules;
    std::map<std::string, CRTStripGeo>  fStrips;
    std::map<uint16_t, CRTSiPMGeo>      fSiPMs;

    geo::GeometryCore const       *fGeometryService;
    const geo::AuxDetGeometryCore *fAuxDetGeoCore;
  };

}

#endif
