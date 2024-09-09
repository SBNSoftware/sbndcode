#ifndef CRTGEOALG_H_SEEN
#define CRTGEOALG_H_SEEN

///////////////////////////////////////////////
// CRTGeoAlg.h
//
// Wrapper for some awkward CRT geometry things
// T Brooks (tbrooks@fnal.gov), November 2018
// Edited heavily - Henry Lay, November 2022
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
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "TGeoManager.h"

// sbndcode
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

namespace sbnd::crt {

  struct CRTSiPMGeo{
    CRTSiPMGeo(const std::string &_stripName, const uint32_t _channel, const geo::Point_t location,
               const uint32_t _pedestal, const double _gain)
    {
      stripName = _stripName;
      channel   = _channel;
      x         = location.X();
      y         = location.Y();
      z         = location.Z();
      pedestal  = _pedestal;
      gain      = _gain;
      null      = false;
    }
    std::string stripName;
    uint16_t    channel;
    double      x;
    double      y;
    double      z;
    bool        null;
    uint32_t    pedestal;
    double      gain;
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
      geo::AuxDetSensitiveGeo::LocalPoint_t const limits{halfWidth, halfHeight, halfLength};
      geo::AuxDetSensitiveGeo::LocalPoint_t const limits2{-halfWidth, -halfHeight, -halfLength};

      auto const limitsWorld  = auxDetSensitive.toWorldCoords(limits);
      auto const limitsWorld2 = auxDetSensitive.toWorldCoords(limits2);

      // Fill edges & width
      minX  = std::min(limitsWorld.X(), limitsWorld2.X());
      maxX  = std::max(limitsWorld.X(), limitsWorld2.X());
      minY  = std::min(limitsWorld.Y(), limitsWorld2.Y());
      maxY  = std::max(limitsWorld.Y(), limitsWorld2.Y());
      minZ  = std::min(limitsWorld.Z(), limitsWorld2.Z());
      maxZ  = std::max(limitsWorld.Z(), limitsWorld2.Z());
      width = halfHeight * 2.;

      adsID = _adsID;
      null  = false;

      std::string volumeName = stripNode->GetVolume()->GetName();
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
    CRTModuleGeo()
    : name("")
    , taggerName("")
    , minX(-std::numeric_limits<double>::max())
    , maxX(std::numeric_limits<double>::max())
    , minY(-std::numeric_limits<double>::max())
    , maxY(std::numeric_limits<double>::max())
    , minZ(-std::numeric_limits<double>::max())
    , maxZ(std::numeric_limits<double>::max())
    , orientation(0)
    , top(false)
    , adID(std::numeric_limits<uint16_t>::max())
    , t0CableDelayCorrection(0)
    , t1CableDelayCorrection(0)
    , invertedOrdering(false)
    , minos(false)
    , null(false)
    {}

    CRTModuleGeo(const TGeoNode *moduleNode, const geo::AuxDetGeo &auxDet,
                 const uint16_t _adID, const std::string &_taggerName,
                 const int32_t _t0CableDelayCorrection,
                 const int32_t _t1CableDelayCorrection,
                 const bool _invertedOrdering,
                 const bool _minos)
    {
      name       = moduleNode->GetName();
      taggerName = _taggerName;

      // Module Dimensions
      double halfWidth  = auxDet.HalfWidth1();
      double halfHeight = auxDet.HalfHeight();
      double halfLength = auxDet.Length()/2;

      // Find world coordinates for edges
      geo::AuxDetGeo::LocalPoint_t const limits{halfWidth, halfHeight, halfLength};
      geo::AuxDetGeo::LocalPoint_t const limits2{-halfWidth, -halfHeight, -halfLength};

      auto const limitsWorld  = auxDet.toWorldCoords(limits);
      auto const limitsWorld2 = auxDet.toWorldCoords(limits2);

      // Which plane within the tagger
      double origin[3] = {0, 0, 0};
      double modulePosMother[3];
      moduleNode->LocalToMaster(origin, modulePosMother);

      if(_minos)
        orientation = (modulePosMother[2] < 0);
      else
        orientation = (modulePosMother[2] > 0);

      // Location of SiPMs
      if(CRTCommonUtils::GetTaggerEnum(taggerName) == kBottomTagger || CRTCommonUtils::GetTaggerEnum(taggerName) == kNorthTagger
         || CRTCommonUtils::GetTaggerEnum(taggerName) == kWestTagger || CRTCommonUtils::GetTaggerEnum(taggerName) == kEastTagger)
        top = (orientation == 1) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);
      else
        top = (orientation == 0) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);

      // Fill edges
      minX = std::min(limitsWorld.X(), limitsWorld2.X());
      maxX = std::max(limitsWorld.X(), limitsWorld2.X());
      minY = std::min(limitsWorld.Y(), limitsWorld2.Y());
      maxY = std::max(limitsWorld.Y(), limitsWorld2.Y());
      minZ = std::min(limitsWorld.Z(), limitsWorld2.Z());
      maxZ = std::max(limitsWorld.Z(), limitsWorld2.Z());

      t0CableDelayCorrection = _t0CableDelayCorrection;
      t1CableDelayCorrection = _t1CableDelayCorrection;

      invertedOrdering = _invertedOrdering;
      adID = _adID;
      minos = _minos;

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
    uint16_t      orientation;
    bool          top;
    uint16_t      adID;
    int32_t       t0CableDelayCorrection;
    int32_t       t1CableDelayCorrection;
    bool          invertedOrdering;
    bool          minos;
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

    CRTGeoAlg(fhicl::ParameterSet const &p, geo::GeometryCore const *geometry,
              geo::AuxDetGeometryCore const *auxdet_geometry);

    CRTGeoAlg(fhicl::ParameterSet const &p = fhicl::ParameterSet());

    ~CRTGeoAlg();

    std::vector<double> CRTLimits() const;

    size_t NumTaggers() const;

    size_t NumModules() const;

    size_t NumStrips() const;

    size_t NumSiPMs() const;

    std::map<std::string, CRTTaggerGeo> GetTaggers() const;

    std::map<std::string, CRTModuleGeo> GetModules() const;

    std::map<std::string, CRTStripGeo> GetStrips() const;

    std::map<uint16_t, CRTSiPMGeo> GetSiPMs() const;

    CRTTaggerGeo GetTagger(const std::string taggerName) const;

    CRTModuleGeo GetModule(const std::string moduleName) const;

    CRTModuleGeo GetModule(const uint16_t channel) const;

    CRTModuleGeo GetModuleByAuxDetIndex(const unsigned ad_i) const;

    CRTStripGeo GetStrip(const std::string stripName) const;

    CRTStripGeo GetStrip(const uint16_t channel) const;

    CRTStripGeo GetStripByAuxDetIndices(const unsigned ad_i, const unsigned ads_i) const;

    CRTSiPMGeo GetSiPM(const uint16_t channel) const;

    std::string GetTaggerName(const std::string name) const;

    std::string ChannelToStripName(const uint16_t channel) const;

    std::string ChannelToTaggerName(const uint16_t channel) const;

    enum CRTTagger ChannelToTaggerEnum(const uint16_t channel) const;

    size_t ChannelToOrientation(const uint16_t channel) const;

    std::array<double, 6> StripHit3DPos(const uint16_t channel, const double x, const double ex);

    std::vector<double> StripWorldToLocalPos(const CRTStripGeo &strip, const double x,
                                             const double y, const double z);

    std::vector<double> StripWorldToLocalPos(const uint16_t channel, const double x,
                                             const double y, const double z);

    std::array<double, 6> FEBWorldPos(const CRTModuleGeo &module);

    std::array<double, 6> FEBChannel0WorldPos(const CRTModuleGeo &module);

    geo::Point_t ChannelToSipmPosition(const uint16_t channel) const;

    std::pair<int, int> GetStripSipmChannels(const std::string stripName) const;

    double DistanceDownStrip(const geo::Point_t position, const std::string stripName) const;

    double DistanceDownStrip(const geo::Point_t position, const uint16_t channel) const;

    bool CheckOverlap(const CRTStripGeo &strip1, const CRTStripGeo &strip2, const double overlap_buffer = 0.);

    bool CheckOverlap(const uint16_t channel1, const uint16_t channel2, const double overlap_buffer = 0.);

    bool AdjacentStrips(const CRTStripGeo &strip1, const CRTStripGeo &strip2, const double overlap_buffer = 0.1);

    bool AdjacentStrips(const uint16_t channel1, const uint16_t channel2, const double overlap_buffer = 0.1);

    bool DifferentOrientations(const CRTStripGeo &strip1, const CRTStripGeo &strip2);

    enum CRTTagger WhichTagger(const double &x, const double &y, const double &z, const double &buffer = 1);

    enum CoordSet GlobalConstrainedCoordinates(const uint16_t channel);

    bool IsPointInsideCRTLimits(const geo::Point_t &point);

  private:

    std::map<std::string, CRTTaggerGeo> fTaggers;
    std::map<std::string, CRTModuleGeo> fModules;
    std::map<std::string, CRTStripGeo>  fStrips;
    std::map<uint16_t, CRTSiPMGeo>      fSiPMs;

    geo::GeometryCore const       *fGeometryService;
    const geo::AuxDetGeometryCore *fAuxDetGeoCore;

    std::vector<std::pair<unsigned, double>> fT0CableLengthCorrectionsVector;
    std::map<unsigned, double>               fT0CableLengthCorrections;
    std::vector<std::pair<unsigned, double>> fT1CableLengthCorrectionsVector;
    std::map<unsigned, double>               fT1CableLengthCorrections;
    double                                   fDefaultPedestal;
    std::vector<std::pair<unsigned, double>> fSiPMPedestalsVector;
    std::map<unsigned, double>               fSiPMPedestals;
    double                                   fDefaultGain;
    std::vector<std::pair<unsigned, double>> fSiPMGainsVector;
    std::map<unsigned, double>               fSiPMGains;
    std::vector<std::pair<unsigned, bool>>   fChannelInversionVector;
    std::map<unsigned, bool>                 fChannelInversion;
  };
}

#endif
