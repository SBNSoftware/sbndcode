////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAPAAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELlar1ndMAPALG_H
#define GEO_CHANNELlar1ndMAPALG_H

#include <vector>
#include <set>
#include <stdint.h>

#include "Geometry/ChannelMapAlg.h"
#include "lar1ndcode/Geo/GeoObjectSorterlar1nd.h"
#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class ChannelMaplar1ndAlg : public ChannelMapAlg{

  public:

    ChannelMaplar1ndAlg(fhicl::ParameterSet const& p);
    ~ChannelMaplar1ndAlg();
    
    void                     Initialize(std::vector<geo::CryostatGeo*> & cgeo);
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(uint32_t channel)        const;
    uint32_t                 Nchannels()                            const;
    WireID                   NearestWireID(const TVector3& worldPos,
					   unsigned int    PlaneNo,
					   unsigned int    TPCNo,
					   unsigned int    cstat)   const;
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

    std::vector< unsigned int >				 fWiresInPlane;
    unsigned int					 fPlanesPerAPA;   
    uint32_t					         fChannelsPerAPA;
    std::vector<std::vector<std::vector<unsigned int>>>	 nAnchoredWires;

    std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy
    geo::GeoObjectSorterlar1nd                               fSorter;         ///< sorts geo::XXXGeo objects
    
    std::vector<std::vector<std::vector<double>>> fFirstWireCenterY;
    std::vector<std::vector<std::vector<double>>> fFirstWireCenterZ;
    std::vector< double > fWirePitch;
    std::vector< double > fOrientation;
    std::vector< double > fTanOrientation; // to explore improving speed
    std::vector< double > fCosOrientation; // to explore improving speed



  };

}
#endif // GEO_CHANNELMAPlar1ndALG_H

