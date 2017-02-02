////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapsbndAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com // adapted for SBND by ariana.hackenburg@yale.edu
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELMAPSBNDALG_H
#define GEO_CHANNELMAPSBNDALG_H

#include <vector>
#include <set>
#include <stdint.h>

#include "larcore/Geometry/ChannelMapAlg.h"
#include "sbndcode/Geo/GeoObjectSortersbnd.h"
#include "fhiclcpp/ParameterSet.h"


namespace geo {

class ChannelMapsbndAlg : public ChannelMapAlg {

public:

  ChannelMapsbndAlg(fhicl::ParameterSet const& p);
  ~ChannelMapsbndAlg();

  void                     Initialize( GeometryData_t& geodata ) override;
  void                     Uninitialize();
  std::vector<WireID>      ChannelToWire(uint32_t channel)        const;
  uint32_t                 Nchannels()                            const;
  WireID                   NearestWireID(const TVector3& worldPos,
                                         unsigned int    PlaneNo,
                                         unsigned int    TPCNo,
                                         unsigned int    cstat)   const;
  double                   WireCoordinate(double Ypos, double Zpos,
                                          unsigned int PlaneNo,
                                          unsigned int TPCNo,
                                          unsigned int cstat) const;


  uint32_t                 PlaneWireToChannel(unsigned int plane, 
                                              unsigned int wire,
                                              unsigned int tpc, 
                                              unsigned int cstat) const;

  View_t                   View( uint32_t const channel )         const;
  SigType_t                SignalType( uint32_t const channel)    const;
  std::set<View_t>  const& Views()                                const;
  std::set<PlaneID> const& PlaneIDs()                             const;

private:

  unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
  uint32_t                                             fNchannels;      ///< number of channels in the detector
  uint32_t                                             fTopChannel;     ///< book keeping highest channel #
  std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat
  std::set<View_t>                                     fViews;          ///< vector of the views present in the detector
  std::set<PlaneID>                                    fPlaneIDs;       ///< vector of the PlaneIDs present in the detector

  std::vector< unsigned int >                          fWiresInPlane;
  unsigned int                                         fPlanesPerAPA;
  uint32_t                                             fChannelsPerAPA;
  std::vector<std::vector<std::vector<unsigned int>>>  nAnchoredWires;

  std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane
  ///< in the heirachy
  geo::GeoObjectSortersbnd                           fSorter;         ///< sorts geo::XXXGeo objects

  /// all data we need for each APA
  typedef struct {
    double fFirstWireCenterY;
    double fFirstWireCenterZ;
    /// +1 if the wire ID order follow z (larger z, or smaller intercept => larger wire ID); -1 otherwise
    float fWireSortingInZ;
  } PlaneData_t;
  ///< collects all data we need for each plane (indices: c t p)
  std::vector<std::vector<std::vector<PlaneData_t>>> fPlaneData;

  std::vector< double > fWirePitch;
  std::vector< double > fOrientation;
  std::vector< double > fSinOrientation; // to explore improving speed
  std::vector< double > fCosOrientation; // to explore improving speed



};

}
#endif // GEO_CHANNELMAPSBNDALG_H

