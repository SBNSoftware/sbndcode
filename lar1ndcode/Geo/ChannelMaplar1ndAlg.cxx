////////////////////////////////////////////////////////////////////////
/// \file  ChannelMaplar1ndAlg.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  tylerdalion@gmail.com
////////////////////////////////////////////////////////////////////////

#include "lar1ndcode/Geo/ChannelMaplar1ndAlg.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo{

  //----------------------------------------------------------------------------
  ChannelMaplar1ndAlg::ChannelMaplar1ndAlg(fhicl::ParameterSet const& p)
    : fSorter(geo::GeoObjectSorterlar1nd(p))
  {
  }

  //----------------------------------------------------------------------------
  ChannelMaplar1ndAlg::~ChannelMaplar1ndAlg()
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMaplar1ndAlg::Initialize(std::vector<geo::CryostatGeo*> & cgeo)
  {

    if(!fFirstChannelInThisPlane.empty() || !fFirstChannelInNextPlane.empty())
      {
	this->Uninitialize();
      }


    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMaplar1ndAlg") << "Sorting...";

    fSorter.SortCryostats(cgeo);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(fSorter);


    mf::LogInfo("ChannelMaplar1ndAlg") << "Initializing...";
      
    fNTPC.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    nAnchoredWires.resize(fNcryostat);
    fViews.clear();
    fPlaneIDs.clear();
    fPlanesPerAPA = cgeo[0]->TPC(0).Nplanes();

    fTopChannel = 0;

    // Size some vectors and initialize the FirstChannel vectors.
    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();

      nAnchoredWires[cs].resize(fNTPC[cs]);
      fWiresPerPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]/2);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]/2);

      for(unsigned int apa = 0; apa != fNTPC[cs]/2; ++apa){
	
        nAnchoredWires[cs][apa].resize(fPlanesPerAPA);
	fWiresPerPlane[cs][apa].resize(fPlanesPerAPA);
        fFirstChannelInThisPlane[cs][apa].resize(fPlanesPerAPA);
        fFirstChannelInNextPlane[cs][apa].resize(fPlanesPerAPA);

      }// end loop over apas
    }// end cryostats

    // Find the number of wires anchored to the frame
    for(unsigned int c = 0; c != fNcryostat; ++c){
      for(unsigned int a = 0; a != fNTPC[c]/2; ++a){
        for(unsigned int p = 0; p != fPlanesPerAPA; ++p){

          unsigned int t = 2*a;
          fWiresPerPlane[c][a][p] = cgeo[c]->TPC(t).Plane(p).Nwires();
          double xyz[3] = {0.};
          double xyz_next[3] = {0.};

	  fViews.emplace(cgeo[c]->TPC(t).Plane(p).View());

          for(unsigned int w = 0; w != fWiresPerPlane[c][a][p]; ++w){

	    // for vertical planes
	    if(cgeo[c]->TPC(t).Plane(p).View() == geo::kZ)   { 
	      nAnchoredWires[c][a][p] = fWiresPerPlane[c][a][p];      
	      break;
	    }

	    cgeo[c]->TPC(t).Plane(p).Wire(w).GetCenter(xyz);
	    cgeo[c]->TPC(t).Plane(p).Wire(w+1).GetCenter(xyz_next);

    	    if(xyz[2]==xyz_next[2]){
	      nAnchoredWires[c][a][p] = w-1;      
	      break;
	    }

          }
        }
      }
    }

    static uint32_t CurrentChannel = 0;
 
    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      for(unsigned int apa = 0; apa != fNTPC[cs]/2; ++apa){  
        for(unsigned int p = 0; p != fPlanesPerAPA; ++p){

          fFirstChannelInThisPlane[cs][apa][p] = CurrentChannel;
          CurrentChannel = CurrentChannel + 2*nAnchoredWires[cs][apa][p];
          fFirstChannelInNextPlane[cs][apa][p] = CurrentChannel;

        }// end plane loop
      }// end apa loop
    }// end cs


    // Save the number of channels
    fNchannels = CurrentChannel;

    // Save the number of channels
    fChannelsPerAPA = fFirstChannelInNextPlane[0][0][fPlanesPerAPA-1];

    //resize vectors
    fFirstWireCenterY.resize(fNcryostat);
    fFirstWireCenterZ.resize(fNcryostat);
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      fFirstWireCenterY[cs].resize(fNTPC[cs]);
      fFirstWireCenterZ[cs].resize(fNTPC[cs]);
      for (unsigned int tpc=0; tpc<fNTPC[cs]; tpc++){
        fFirstWireCenterY[cs][tpc].resize(fPlanesPerAPA);
        fFirstWireCenterZ[cs][tpc].resize(fPlanesPerAPA);
      }                                                                   
    }

    fWirePitch.resize(fPlanesPerAPA);
    fOrientation.resize(fPlanesPerAPA);
    fTanOrientation.resize(fPlanesPerAPA);
    fCosOrientation.resize(fPlanesPerAPA);


    //save data into fFirstWireCenterY and fFirstWireCenterZ
    for (unsigned int cs=0; cs<fNcryostat; cs++){
      for (unsigned int tpc=0; tpc<fNTPC[cs]; tpc++){
        for (unsigned int plane=0; plane<fPlanesPerAPA; plane++){
	  fPlaneIDs.emplace(PlaneID(cs, tpc, plane));
          double xyz[3]={0.0, 0.0, 0.0};
          cgeo[cs]->TPC(tpc).Plane(plane).Wire(0).GetCenter(xyz);
          fFirstWireCenterY[cs][tpc][plane]=xyz[1];
          fFirstWireCenterZ[cs][tpc][plane]=xyz[2];
        }
      }
    }

    //initialize fWirePitch and fOrientation
    for (unsigned int plane=0; plane<fPlanesPerAPA; plane++){
      fWirePitch[plane]=cgeo[0]->TPC(0).WirePitch(0,1,plane);
      fOrientation[plane]=cgeo[0]->TPC(0).Plane(plane).Wire(0).ThetaZ();
      fTanOrientation[plane] = tan(fOrientation[plane]);
      fCosOrientation[plane] = cos(fOrientation[plane]);
    }


    mf::LogVerbatim("GeometryTest") << "fNchannels = " << fNchannels ; 
    mf::LogVerbatim("GeometryTest") << "U channels per APA = " << 2*nAnchoredWires[0][0][0] ;
    mf::LogVerbatim("GeometryTest") << "V channels per APA = " << 2*nAnchoredWires[0][0][1] ;
    mf::LogVerbatim("GeometryTest") << "Z channels per APA side = " << nAnchoredWires[0][0][2] ;

    return;

  }
   
  //----------------------------------------------------------------------------
  void ChannelMaplar1ndAlg::Uninitialize()
  {

    std::vector< std::vector<std::vector<uint32_t> > >().swap(fFirstChannelInThisPlane);
    std::vector< std::vector<std::vector<uint32_t> > >().swap(fFirstChannelInNextPlane);

  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMaplar1ndAlg::ChannelToWire(uint32_t channel)  const
  {

    // first check if this channel ID is legal
    if(channel >= fNchannels )
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";

    std::vector< WireID > AllSegments;
    
    static unsigned int cstat;
    static unsigned int tpc;
    static unsigned int plane;
    static unsigned int wireThisPlane;
    static unsigned int NextPlane;
    static unsigned int ThisPlane;
    
    for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
      
      bool breakVariable = false;
      
      for(unsigned int apaloop = 0; apaloop != fNTPC[csloop]/2; ++apaloop){
	for(unsigned int planeloop = 0; planeloop != fPlanesPerAPA; ++planeloop){
	  
	  NextPlane = fFirstChannelInNextPlane[csloop][apaloop][planeloop];
       	  ThisPlane = fFirstChannelInThisPlane[csloop][apaloop][planeloop];
	  
	  if(channel < NextPlane){
	    
	    cstat = csloop;
	    tpc   = 2*apaloop;
	    plane = planeloop;
	    wireThisPlane  = channel - ThisPlane;
	    
	    breakVariable = true;
	    break;
	  }// end if break	  
	  if(breakVariable) break;
	  
	}// end plane loop	
	if(breakVariable) break;
	
      }// end apa loop      
      if(breakVariable) break;
      
    }// end cryostat loop
    

    int WrapDirection = 1; // go from tpc to (tpc+1) or tpc to (tpc-1)

    // find the lowest wire
    uint32_t ChannelGroup = std::floor( wireThisPlane/nAnchoredWires[cstat][tpc/2][plane] );
    unsigned int bottomwire = wireThisPlane-ChannelGroup*nAnchoredWires[cstat][tpc/2][plane];
    
    if(ChannelGroup%2==1){
      // start in the other TPC
      tpc += 1;
      WrapDirection  = -1;	 
    }
    
    for(unsigned int WireSegmentCount = 0; WireSegmentCount != 50; ++WireSegmentCount){
      
      tpc += WrapDirection*(WireSegmentCount%2);
      
      geo::WireID CodeWire(cstat, tpc, plane, bottomwire + WireSegmentCount*nAnchoredWires[cstat][std::floor(tpc/2)][plane]);
      
      AllSegments.push_back(CodeWire);
      
      // reset the tcp variable so it doesnt "accumulate value"
      tpc -= WrapDirection*(WireSegmentCount%2);
      
      if( bottomwire + (WireSegmentCount+1)*nAnchoredWires[cstat][std::floor(tpc/2)][plane] > 
	  fWiresPerPlane[cstat][std::floor(tpc/2)][plane]-1) break;
      
    } //end WireSegmentCount loop
    
    
    return AllSegments;
  }


  //----------------------------------------------------------------------------
  uint32_t ChannelMaplar1ndAlg::Nchannels() const
  {
    return fNchannels;
  }
  

  //----------------------------------------------------------------------------
  WireID  ChannelMaplar1ndAlg::NearestWireID(const TVector3& xyz,
					 unsigned int    plane,
					 unsigned int    tpc,
					 unsigned int    cryostat)     const
  {

    //get the position of first wire in a given cryostat, tpc and plane
    double firstxyz[3]={0.0, 0.0, 0.0};
    firstxyz[1]=fFirstWireCenterY[cryostat][tpc][plane];
    firstxyz[2]=fFirstWireCenterZ[cryostat][tpc][plane];

    double distance = 0.;

    //get the orientation angle of a given plane and calculate the distance between first wire
    //and a point projected in the plane
    int rotate = 1;
    if (tpc%2 == 1) rotate = -1;
 
    distance = std::abs( (xyz[1]-firstxyz[1] -rotate*fTanOrientation[plane]*(xyz[2]-firstxyz[2]))
                                * fCosOrientation[plane]);

    //if the distance between the wire and a given point is greater than the half of wirepitch,
    //then the point is closer to a i+1 wire thus add one
    //double res = distance/fWirePitch[plane] - int( distance/fWirePitch[plane] );
    //if (res > fWirePitch[plane]/2)	iwire+=1;

    // do it, but also check to see if we are on the edge

    double dwire=distance/fWirePitch[plane];
    uint32_t iwire=int(dwire);
    if (dwire-iwire>fWirePitch[plane]*0.5) ++iwire;
    uint32_t maxwireminus1=fWiresPerPlane[0][tpc/2][plane]-1;
    if(iwire>maxwireminus1) iwire=maxwireminus1;

    WireID wid(cryostat, tpc, plane, iwire);
    return wid;

  }
  
  //----------------------------------------------------------------------------
  uint32_t ChannelMaplar1ndAlg::PlaneWireToChannel(unsigned int plane,
					       unsigned int wire,
					       unsigned int tpc,
					       unsigned int cstat) const
  {

    unsigned int OtherSideWires = 0;

    uint32_t Channel = fFirstChannelInThisPlane[cstat][std::floor(tpc/2)][plane];

    // get number of wires starting on the first side of the APA if starting
    // on the other side TPC.
    OtherSideWires += (tpc%2)*nAnchoredWires[cstat][std::floor(tpc/2)][plane];
    
    // Lastly, account for the fact that channel number while moving up wire number in one
    // plane resets after 2 times the number of wires anchored -- one for each APA side.
    // At the same time, OtherSideWires accounts for the fact that if a channel starts on
    // the other side, it is offset by the number of wires on the first side.
    Channel += (OtherSideWires + wire)%(2*nAnchoredWires[cstat][std::floor(tpc/2)][plane]);

    return Channel;

  }


  //----------------------------------------------------------------------------
  SigType_t ChannelMaplar1ndAlg::SignalType( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    SigType_t sigt = kInduction;

    if(       chan <  fFirstChannelInThisPlane[0][0][2]     ){ sigt = kInduction;  }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ sigt = kCollection; }
    else{    mf::LogWarning("BadChannelSignalType") << "Channel " << channel 
						    << " (" << chan << ") not given signal type." << std::endl;         }
  
    return sigt;
  }

  //----------------------------------------------------------------------------
  View_t ChannelMaplar1ndAlg::View( uint32_t const channel )  const
  {
    uint32_t chan = channel % fChannelsPerAPA;
    View_t view = geo::kU;

    if(       chan <  fFirstChannelInNextPlane[0][0][0]     ){ view = geo::kU; }
    else if( (chan >= fFirstChannelInThisPlane[0][0][1]) &&
             (chan <  fFirstChannelInNextPlane[0][0][1])    ){ view = geo::kV; }
    else if( (chan >= fFirstChannelInThisPlane[0][0][2]) &&
             (chan <  fFirstChannelInNextPlane[0][0][2])    ){ view = geo::kZ; }
    else{    mf::LogWarning("BadChannelViewType") << "Channel " << channel 
						  << " (" << chan << ") not given view type.";}
    
    return view;
  }  
 
  //----------------------------------------------------------------------------
  std::set<View_t> const& ChannelMaplar1ndAlg::Views() const
  {
    return fViews;
  }

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& ChannelMaplar1ndAlg::PlaneIDs() const
  {
    return fPlaneIDs;
  }

} // namespace
