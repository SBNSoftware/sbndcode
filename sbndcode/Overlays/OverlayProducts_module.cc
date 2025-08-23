////////////////////////////////////////////////////////////////////////
// Class:       OverlayProducts
// Plugin Type: producer (Unknown Unknown)
// File:        OverlayProducts_module.cc
//
// Generated at Fri Sep 23 17:04:48 2022 by Bruce Howard using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

// NB: See https://github.com/SBNSoftware/icaruscode/blob/4ea3861350628402b7edccfe46efdd403e4c0432/icaruscode/TPC/SignalProcessing/RecoWire/ROIConverter_module.cc

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "cetlib_except/exception.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
//#include "icaruscode/IcarusObj/ChannelROI.h"
#include "sbnobj/ICARUS/TPC/ChannelROI.h"
// #include "icaruscode/IcarusObj/classes.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// PROBABLY TEMPORARY
//#include "art/Framework/Services/Registry/ServiceHandle.h" 
//#include "larcore/CoreUtils/ServiceUtil.h"
//#include "larcore/Geometry/Geometry.h"
/////////////////////

#include "sbnobj/Common/CRT/CRTHit.hh"

#include "art/Persistency/Common/PtrMaker.h"

#include <iostream>
#include <memory>
#include <string>

class OverlayProducts;


class OverlayProducts : public art::EDProducer {
public:
  explicit OverlayProducts(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OverlayProducts(OverlayProducts const&) = delete;
  OverlayProducts(OverlayProducts&&) = delete;
  OverlayProducts& operator=(OverlayProducts const&) = delete;
  OverlayProducts& operator=(OverlayProducts&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  struct opWaveformOverlay {
    raw::TimeStamp_t timestamp;
    std::vector< raw::ADC_Count_t > wvfm;
    float baseline; // TODO: does SBND have/plan to have such a product?
    double deltaT = 0.002; // TODO: automate? SBND?

    raw::TimeStamp_t endTime() const { return timestamp+(deltaT*wvfm.size()); }
    std::pair< raw::TimeStamp_t, raw::TimeStamp_t > timeBin( unsigned int binNumber ) { return std::make_pair(timestamp+(deltaT*binNumber),
													      timestamp+(deltaT*(binNumber+1))); } 
  };

private:

  // Declare member data here.
  // TPC
  bool fTPCOverlayRaw;   //< Set true if you should overlay at the RawDigit stage
  bool fTPCOverlayROI;   //< Set true if you should overlay at the recob::ChannelROI stage
  bool fTPCOverlayHits;  //< Set true if you should overlay at the recob::Hit stage
  bool fTPCHitsWireAssn; //< Set true if you want to make the hit<->wire association in HitCollectionCreator
  std::vector< art::InputTag > fTPCRawInputLabels;
  art::InputTag                fTPCRawOutputLabel;
  std::vector< art::InputTag > fTPCROIInputLabels;
  std::vector< art::InputTag > fTPCHitInputLabels;
  std::string fTPCHitCreatorInstanceName;

  // PMT
  bool fPMTOverlayRaw;
  bool fPMTOverlayHits;
  art::InputTag fPMTWaveDataLabel;
  art::InputTag fPMTWaveBaseLabel;
  art::InputTag fPMTWaveSimLabel;
  int fPMTWaveTestCh;
  std::vector< art::InputTag > fPMTHitInputLabels;

  // CRT
  bool fCRTOverlayHits;
  std::vector< art::InputTag > fCRTHitInputLabels;

  // PROBABLY TEMPORARY
  //geo::GeometryCore const* fGeom;
};


OverlayProducts::OverlayProducts(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTPCOverlayRaw      ( p.get< bool >("TPCOverlayRaw", false) ),
  fTPCOverlayROI      ( p.get< bool >("TPCOverlayROI", false) ),
  fTPCOverlayHits     ( p.get< bool >("TPCOverlayHits", false) ),
  fTPCHitsWireAssn    ( p.get< bool >("TPCHitsWireAssn", true) ),
  fTPCRawInputLabels  ( p.get< std::vector<art::InputTag> >("TPCRawInputLabels", {}) ),
  fTPCRawOutputLabel  ( p.get< art::InputTag >("TPCRawOutputLabel","") ),
  fTPCROIInputLabels  ( p.get< std::vector<art::InputTag> >("TPCROIInputLabels", {}) ),
  fTPCHitInputLabels  ( p.get< std::vector<art::InputTag> >("TPCHitInputLabels", {}) ),
  fTPCHitCreatorInstanceName ( p.get<std::string>("TPCHitCreatorInstanaceName","") ),
  fPMTOverlayRaw      ( p.get< bool >("PMTOverlayRaw", false) ),
  fPMTOverlayHits     ( p.get< bool >("PMTOverlayHits", false) ),
  fPMTWaveDataLabel   ( p.get< art::InputTag >("PMTWaveDataLabel", "") ),
  fPMTWaveBaseLabel   ( p.get< art::InputTag >("PMTWaveBaseLabel", "") ),
  fPMTWaveSimLabel    ( p.get< art::InputTag >("PMTWaveSimLabel", "") ),
  fPMTWaveTestCh      ( p.get< int >("PMTWaveTestCh", -1) ),
  fPMTHitInputLabels  ( p.get< std::vector<art::InputTag> >("PMTHitInputLabels", {}) ),
  fCRTOverlayHits     ( p.get< bool >("CRTOverlayHits", false) ),
  fCRTHitInputLabels  ( p.get< std::vector<art::InputTag> >("CRTHitInputLabels", {}) )
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  if ( !fTPCOverlayRaw && !fTPCOverlayROI && !fTPCOverlayHits && !fPMTOverlayRaw && !fPMTOverlayHits && !fCRTOverlayHits ) {
    throw cet::exception("OverlayProducts") << "Error... need to be doing SOME overlay." << std::endl;
  }

  if ( (fTPCOverlayRaw && fTPCRawInputLabels.size()<2) ||
       (fTPCOverlayROI && fTPCROIInputLabels.size()<2) ||
       (fTPCOverlayHits && fTPCHitInputLabels.size()<2) ||
       (fPMTOverlayRaw && fPMTWaveDataLabel==fPMTWaveSimLabel) ||
       (fPMTOverlayHits && fPMTHitInputLabels.size()<2) ||
       (fCRTOverlayHits && fCRTHitInputLabels.size()<2) )
  {
    throw cet::exception("OverlayProducts") << "Error... Need at least 2 input products to merge." << std::endl;
  }

  if ( fTPCOverlayRaw ) produces< std::vector<raw::RawDigit> >(fTPCRawOutputLabel.instance());
  if ( fTPCOverlayROI ) produces< std::vector<recob::ChannelROI> >();
  // basically from the GausHitFinder
  if ( fTPCOverlayHits ) recob::HitCollectionCreator::declare_products(producesCollector(), fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);

  if ( fPMTOverlayRaw ) {
    produces< std::vector<raw::OpDetWaveform> >();
    produces< std::vector<icarus::WaveformBaseline> >();
    produces< art::Assns<icarus::WaveformBaseline,raw::OpDetWaveform> >();
  }
  if ( fPMTOverlayHits ) produces< std::vector<recob::OpHit> >();

  if ( fCRTOverlayHits ) produces< std::vector<sbn::crt::CRTHit> >();
}

void OverlayProducts::produce(art::Event& e)
{
  // Implementation of required member function here.

  // PROBABLY TEMPORARY
  //fGeom = lar::providerFrom<geo::Geometry>();

  // TODO: Add CRT systems...

  // NOTE! IF OVERLAYING RAW DIGITS/WAVEFORMS, DO --DATA-- (OR JUST WHATEVER SAMPLE HAS THE NOISE AND THE BASE TIME) FIRST IN THE VECTOR OF LABELS...

  // TODO: recob::Hit mixing will just join the two products together for now. Probably need better version in future...
  //     NOTE: see also HitMerger_module.cc

  if ( fTPCOverlayRaw )
  {
    // Need to see what Assns exist for the raw::RawDigits but this is a start...
    std::map< raw::ChannelID_t, raw::RawDigit::ADCvector_t > rawdigitMap;
    std::map< raw::ChannelID_t, float > pedestalMap;
    std::map< raw::ChannelID_t, float > sigmaMap;
    std::map< raw::ChannelID_t, raw::Compress_t > compressMap;

    for ( auto const& iLabel : fTPCRawInputLabels ) {
      art::Handle< std::vector<raw::RawDigit> > digitsHandle;
      std::vector< art::Ptr<raw::RawDigit> > digits;
      if ( e.getByLabel(iLabel,digitsHandle) ) {
	art::fill_ptr_vector(digits,digitsHandle);
      }
      else{
	mf::LogWarning("OverlayProducts") << "Event failed to find raw::RawDigit with label " << iLabel << ".";
        return;
      }

      for ( auto const& iDigit : digits ) {
	// Since the overlaid product is pure signal (no noise) we will keep the same pedestal and sigma, so only write to the map if no previous obect...
	auto chID = iDigit->Channel();
	//std::vector<geo::WireID> wire = fGeom->ChannelToWire( chID );
	//std::cout << chID << " " << wire[0].Cryostat << " " << wire[0].TPC << " " << wire[0].Plane << " " << wire[0].Wire << std::endl;
	if ( pedestalMap.find( chID ) == pedestalMap.end() ) {
	  pedestalMap[ chID ] = iDigit->GetPedestal();
	  sigmaMap[ chID ] = iDigit->GetSigma();
	  compressMap[ chID ] = iDigit->Compression();
	  rawdigitMap[ chID ] = iDigit->ADCs();
	}
	else {
	  auto const& thisADC = iDigit->ADCs();
	  for ( unsigned int idx = 0; idx < rawdigitMap[ chID ].size(); ++idx ) {
	    rawdigitMap[ chID ][ idx ] += thisADC[ idx ];
	  }
	}
      } // loop RawDigits
    } // loop labels

    // Save the overlaid product
    std::unique_ptr< std::vector< raw::RawDigit > > rawDigitVec = std::make_unique< std::vector< raw::RawDigit > >();

    //unsigned int intCheck=0;
    for ( auto const &[chKey, adcVals] : rawdigitMap ) {
      const raw::ChannelID_t channel = chKey;
      const raw::RawDigit::ADCvector_t ADCs = adcVals;

      const float pedestal = pedestalMap[ chKey ];
      const float sigma = sigmaMap[ chKey ];
      const raw::Compress_t compression = compressMap[ chKey ];

      //std::cout << intCheck << " " << channel << " " << ADCs.size() << " " << compression << std::endl;

      raw::RawDigit rawDigit(channel, ADCs.size(), ADCs, compression);
      rawDigit.SetPedestal( pedestal, sigma );
      rawDigitVec->push_back(rawDigit);

      //std::cout << "    " << rawDigit.GetSigma() << std::endl;
      //intCheck+=1;
    }
    e.put( std::move(rawDigitVec) );
  }
  /*
  if ( fTPCOverlayROI ){
    std::map< raw::ChannelID_t, std::vector< std::vector<short int>::size_type > > chOffsetsMap;
    std::map< raw::ChannelID_t, std::vector< std::vector<short int> > > chDataMap;

    for ( auto const& iLabel : fTPCROIInputLabels ) {
      art::Handle< std::vector<recob::ChannelROI> > chROIsHandle;
      std::vector< art::Ptr<recob::ChannelROI> > chROIs;
      if ( e.getByLabel(iLabel,chROIsHandle) ) {
	art::fill_ptr_vector(chROIs,chROIsHandle);
      }
      else{
	mf::LogWarning("OverlayProducts") << "Event failed to find recob::ChannelROIs with label " << iLabel << ".";
        return;
      }

      for ( auto const& ichROI : chROIs ) {
	//std::cout << "New ChannelROI" << std::endl;
	const recob::ChannelROI::RegionsOfInterest_t& channelROIs = ichROI->SignalROI();
	const raw::ChannelID_t chid = ichROI->Channel();

	for ( auto const& range : channelROIs.get_ranges() ) {
	  //size_t startTick = range.begin_index();
	  //std::vector<float> dataVec(range.data().size());
	  //for(size_t binIdx = 0; binIdx < range.data().size(); binIdx++) dataVec[binIdx] = range.data()[binIdx];
	  //std::cout << "    range, start: " << startTick << ", " << dataVec.size() << " size -- [0] = " << dataVec[0] << std::endl;

	  std::vector<short int> dataVec(range.data().size());
	  for(size_t binIdx=0; binIdx < range.data().size(); ++binIdx) dataVec[binIdx] = range.data()[binIdx];

	  // If channel not in the map already then add it to the map. Otherwise, add the ROI to the map, checking if it overlaps with an already present one...
	  if ( chOffsetsMap.find( chid ) == chOffsetsMap.end() ) {
	    chOffsetsMap[ chid ] = std::vector< std::vector<short int>::size_type >{ range.begin_index() };
	    chDataMap[ chid ] = std::vector< std::vector<short int> >{ dataVec };
	  }
	  else {
	    // Check if this ROI overlaps with any other region for this channel
	    std::vector<short int>::size_type loRange = range.begin_index();
	    std::vector<short int>::size_type hiRange = range.begin_index() + dataVec.size();

	    std::vector< unsigned int > overlaps;

	    for( unsigned int iOffset=0; iOffset < chOffsetsMap[chid].size(); ++iOffset ) {	      
	      if ( ( loRange > chOffsetsMap[chid].at(iOffset) && loRange < chOffsetsMap[chid].at(iOffset)+chDataMap[chid].at(iOffset).size() ) ||
		   ( hiRange > chOffsetsMap[chid].at(iOffset) && hiRange < chOffsetsMap[chid].at(iOffset)+chDataMap[chid].at(iOffset).size() ) ) {
		overlaps.push_back( iOffset );
	      }
	    }

	    // If no overlaps, put the ROI into the map and move on
	    if ( overlaps.size() == 0 ) {
	      chOffsetsMap[ chid ].push_back( loRange );
	      chDataMap[ chid ].push_back( dataVec );
	      continue;
	    }

	    // If overlaps, find the min offset and max offset+size to appropriately pad out the vectors
	    unsigned int idxMinOffset = 0;
	    std::vector<short int>::size_type minOffset = loRange;
	    std::vector<short int>::size_type maxEnd = hiRange;
	    for ( auto const& overlap : overlaps ) {
	      if ( chOffsetsMap[chid].at(overlap) < minOffset ){
		minOffset = chOffsetsMap[chid].at(overlap);
		idxMinOffset = overlap;
	      }
	      if ( maxEnd < chOffsetsMap[chid].at(overlap)+chDataMap[chid].at(overlap).size() ) maxEnd = chOffsetsMap[chid].at(overlap)+chDataMap[chid].at(overlap).size();
	    }

	    // Sum up the ROIs
	    // 1. pad this data vector
	    for ( unsigned int iBegin=0; iBegin < loRange-minOffset; ++iBegin) dataVec.insert( dataVec.begin(), 0 );
	    for ( unsigned int iEnd=0; iEnd < maxEnd-hiRange; ++iEnd) dataVec.push_back( 0 );
	    // 2. now add in the other vectors
	    for ( auto const& overlap : overlaps ) {
	      std::vector<short int>::size_type delta = chOffsetsMap[chid].at(overlap) - minOffset;
	      for ( unsigned int idxData=0; idxData<chDataMap[chid].at(overlap).size(); ++idxData ) {
		dataVec[ delta+idxData ]+=chDataMap[chid].at(overlap).at(idxData);
	      }
	    }
	    // 3. Now delete the current entries in the maps and replace them with this new entry to the map...
	    for ( unsigned int iOverlap=0; iOverlap < overlaps.size(); ++iOverlap ) {
	      auto const& overlap = overlaps[ overlaps.size() - iOverlap];
	      chOffsetsMap[ chid ].erase( chOffsetsMap[ chid ].begin()+overlap );
	      chDataMap[ chid ].erase( chDataMap[ chid ].begin()+overlap );
	    }
	    chOffsetsMap[ chid ].insert( chOffsetsMap[ chid ].begin()+idxMinOffset, minOffset );
	    chDataMap[ chid ].insert( chDataMap[ chid ].begin()+idxMinOffset, dataVec );
	  }
	} // ROI ranges
      } // chROIs
    } // labels

    // Now put the products into the event
    std::unique_ptr< std::vector< recob::ChannelROI > > channelROIVec = std::make_unique< std::vector< recob::ChannelROI > >();

    for ( auto const &[chKey, offsetVec] : chOffsetsMap ) {
      recob::ChannelROI::RegionsOfInterest_t ROIVec;
      for ( unsigned int idxROI = 0; idxROI < offsetVec.size(); ++idxROI ) {
	ROIVec.add_range(offsetVec[idxROI], chDataMap[chKey].at(idxROI) );
      }

      recob::ChannelROI thisChROI(ROIVec,chKey);
      channelROIVec->push_back(thisChROI);
    }

    e.put( std::move(channelROIVec) );
  }*/
  
  if ( fTPCOverlayHits )
  {
    // As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc)
    recob::HitCollectionCreator hitCol(e, fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);

    // Loop through hit labels and put together the hit set
    for ( auto const& iLabel : fTPCHitInputLabels ) {
      art::Handle< std::vector<recob::Hit> > hitsHandle;
      std::vector< art::Ptr<recob::Hit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
	art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
	mf::LogWarning("OverlayProducts") << "Event failed to find recob::Hit with label " << iLabel << ".";
	return;
      }

      art::FindManyP<recob::Wire> fmwire(hitsHandle, e, iLabel);
      if( !fmwire.isValid() && fTPCHitsWireAssn ){
	mf::LogError("OverlayProducts") << "Error in validity of fmwire. Returning.";
	return;
      }

      for ( auto const& iHitPtr : hits ) {
	recob::Hit theHit = *iHitPtr;
	// And get the associated wire -- if we should
	if ( fTPCHitsWireAssn ) {
	  std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at(iHitPtr.key());
	  if ( hitWires.size() == 0 ) throw cet::exception("OverlayProducts") << "Hit found with no associated wires...\n";
	  else if ( hitWires.size() > 1 ) mf::LogWarning("OverlayProducts") << "Hit with >1 recob::Wire associated...";
	  hitCol.emplace_back(theHit,hitWires[0]);
	}
	else hitCol.emplace_back(theHit);
      }
    }

    // Put the hit into the event
    hitCol.put_into(e);
  }

  // PMT Overlays
  // Get the data waveforms and populate the map. Then get the sim waveforms and put them into the right waveforms.
  if ( fPMTOverlayRaw ) {
    // Produced vectors
    std::unique_ptr< std::vector< raw::OpDetWaveform > > opWaveformVec = std::make_unique< std::vector< raw::OpDetWaveform > >();
    std::unique_ptr< std::vector< icarus::WaveformBaseline > > baselineVec = std::make_unique< std::vector< icarus::WaveformBaseline > >();
    art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform> baselineToWaveforms;

    std::map< raw::Channel_t, std::vector< opWaveformOverlay > > opWaveformsMap;
    // Load in data waveforms
    art::Handle< std::vector<raw::OpDetWaveform> > dataWavesHandle;
    std::vector< art::Ptr<raw::OpDetWaveform> > dataWaves;
    if ( e.getByLabel(fPMTWaveDataLabel,dataWavesHandle) ) {
      art::fill_ptr_vector(dataWaves,dataWavesHandle);
    }
    else{
      mf::LogWarning("OverlayProducts") << "Event failed to find raw::OpDetWaveform with label " << fPMTWaveDataLabel << ".";
      return;
    }
    // ... and the baselines
    art::Handle< std::vector<icarus::WaveformBaseline> > baselinesHandle;
    std::vector< art::Ptr<icarus::WaveformBaseline> > baselines;
    if ( e.getByLabel(fPMTWaveBaseLabel,baselinesHandle) ) {
      art::fill_ptr_vector(baselines,baselinesHandle);
    }
    else{
      mf::LogWarning("OverlayProducts") << "Event failed to find raw::OpDetWaveform with label " << fPMTWaveBaseLabel << ".";
      return;
    }
    art::FindManyP<icarus::WaveformBaseline> fmbase( dataWavesHandle, e, fPMTWaveBaseLabel );
    if( !fmbase.isValid() ) {
      mf::LogError("OverlayProducts") << "Error in validity of fmbase. Returning.";
      return;
    }
    // TEST!!!! PRINT BASELINES ////////////////
    std::cout << "// ---- BASELINES ---- " << std::endl;
    for ( unsigned int idxBase=0; idxBase<baselines.size(); ++idxBase ) {
      std::cout << baselines[idxBase]->baseline() << " ";
    }
    std::cout << "---- BASELINES ---- //" << std::endl;
    std::cout << std::endl;
    ////////////////////////////////////////////
    // Now fill the map
    for ( auto const& dataWave : dataWaves ) {
      auto chID = dataWave->ChannelNumber();
      auto time = dataWave->TimeStamp();
      auto wvfm = dataWave->Waveform();

      if ( wvfm.size() == 0 ) continue;

      if ( opWaveformsMap.find( chID ) == opWaveformsMap.end() )
	opWaveformsMap[ chID ] = std::vector< opWaveformOverlay >();

      auto baselines = fmbase.at( dataWave.key() );
      if ( baselines.size() > 1 ) std::cout << "SHOULD NOT BE MORE THAN ONE BASELINE..." << std::endl;

      opWaveformOverlay dataOpWave;
      dataOpWave.timestamp = time;
      dataOpWave.wvfm = wvfm;
      dataOpWave.baseline = ( baselines.size()>0 ? baselines[0]->baseline() : wvfm[0] ); // TODO: better way of guarding?

      // Print
      if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	std::cout << "DATA WAVEFORM AT TIME " << std::setprecision(8) << time << ": [ ";
	for ( unsigned int i=0; i<wvfm.size(); ++i ) {
	  std::cout << wvfm[i];
	  if ( i==wvfm.size()-1 ) std::cout<< " ]" << std::endl;
	  else std::cout << ", ";
	}
      }

      opWaveformsMap[chID].push_back( dataOpWave );
    }

    // Load in the simulation waveforms and do overlay
    art::Handle< std::vector<raw::OpDetWaveform> > simWavesHandle;
    std::vector< art::Ptr<raw::OpDetWaveform> > simWaves;
    if ( e.getByLabel(fPMTWaveSimLabel,simWavesHandle) ) {
      art::fill_ptr_vector(simWaves,simWavesHandle);
    }
    else{
      mf::LogWarning("OverlayProducts") << "Event failed to find raw::OpDetWaveform with label " << fPMTWaveSimLabel << ".";
      return;
    }
    for ( auto const& simWave : simWaves ) {
      auto chID = simWave->ChannelNumber();
      auto time = simWave->TimeStamp();
      auto wvfm = simWave->Waveform();
      raw::ADC_Count_t simBaseline = 14999; // TODO: automate, parameter, not hard-code?

      if ( wvfm.size() == 0 ) continue;

      // PRINT
      if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	std::cout << "SIM WAVEFORM AT TIME " << std::setprecision(8) << time << ": [ ";
	for ( unsigned int i=0;i<wvfm.size(); ++i ) {
	  std::cout << wvfm[i];
          if ( i==wvfm.size()-1 ) std::cout<< " ]" << std::endl;
          else std::cout << ", ";
        }
      }

      // if this channel isn't in the map, then let's make a new struct to start the map
      if ( opWaveformsMap.find( chID ) == opWaveformsMap.end() ) {
	opWaveformOverlay simOpWave;
	simOpWave.timestamp = time;
	simOpWave.wvfm = wvfm;
	simOpWave.baseline = simBaseline;
	opWaveformsMap[chID] = {simOpWave};
      }
      // ... else let's look for overlaps
      else {
	std::vector< unsigned int > overlaps;
	if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	  std::cout << "Overlap idx vals: [ ";
	}
	for ( unsigned int idxWave=0; idxWave<opWaveformsMap[ chID ].size(); ++idxWave ) {
	  const opWaveformOverlay thisWaveformOverlay = opWaveformsMap[ chID ].at(idxWave);
	  double simEnd = time + wvfm.size()*thisWaveformOverlay.deltaT;

	  if ( ( time < thisWaveformOverlay.timestamp && simEnd > thisWaveformOverlay.timestamp ) ||
	       ( time < thisWaveformOverlay.endTime() && simEnd > thisWaveformOverlay.endTime() ) ) {
	    overlaps.push_back( idxWave );
	  }

	  if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	    std::cout << thisWaveformOverlay.timestamp << "-to-" << thisWaveformOverlay.endTime() << " ";
	  }
	} // iterate saved waveforms looking for overlap
	if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	  std::cout << " ]" << std::endl;
	}

	if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	  std::cout << std::endl;
	  std::cout << "OVERLAPS = " << overlaps.size() << std::endl;
	  std::cout << "  [ ";
	  for ( auto const& i : overlaps ) std::cout << i << " ";
	  std::cout << "]" << std::endl;
	  std::cout << std::endl;
	}

	// Case 1: No overlaps. Make a new entry in the map as above...
	if ( overlaps.size() == 0 ) {
	  opWaveformOverlay simOpWave;
	  simOpWave.timestamp = time;
	  simOpWave.wvfm = wvfm;
	  simOpWave.baseline = simBaseline;
	  opWaveformsMap[chID].push_back(simOpWave);
	}
	// Case 2: 1+ overlaps. Make a new struct and erase the overlaps
	if ( overlaps.size() >= 1 ) {
	  // Make a map from timestamp to idx within overlap to find the minimum time and then break these out into vectors
	  // TODO: better way?
	  std::map< raw::TimeStamp_t, unsigned int > timestampIdxMap;
	  for ( unsigned int idxOverlap=0; idxOverlap<overlaps.size(); ++idxOverlap ) {
	    timestampIdxMap[ opWaveformsMap[ chID ].at( overlaps[idxOverlap] ).timestamp ] = overlaps[idxOverlap];
	  }
	  std::vector< raw::TimeStamp_t > orderedTimes;
	  std::vector< unsigned int > orderedOverlaps;
	  for ( auto const &[time, idx] : timestampIdxMap ) {
	    orderedTimes.push_back(time);
	    orderedOverlaps.push_back(idx);
	  }

	  if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	    std::cout << "Sim Time: " << time << std::endl;
	    std::cout << "Times: [ ";
	    for ( unsigned int i=0; i<orderedTimes.size(); ++i ) {
	      std::cout << orderedTimes[i];
	      if ( i == orderedTimes.size()-1 ) std::cout << " ]" << std::endl;
	      else std::cout << ", ";
	    }
	  }

	  // New struct
	  opWaveformOverlay overlayOpWave;
	  overlayOpWave.timestamp = std::min( time, orderedTimes[0] );
	  overlayOpWave.baseline = opWaveformsMap[ chID ].at( orderedOverlaps[0] ).baseline;
	  std::vector< raw::ADC_Count_t > adcVec;
	  double deltaT = opWaveformsMap[ chID ].at( orderedOverlaps[0] ).deltaT;

	  unsigned int idxWvfmEntry=0;
	  double loTime=std::numeric_limits<double>::lowest();
	  for ( unsigned int iOverlap=0; iOverlap<orderedOverlaps.size(); ++iOverlap ) {
	    // If time < first overlay candidate, start adding before
	    while ( idxWvfmEntry<wvfm.size() && time+((double(idxWvfmEntry)+0.5)*deltaT) > loTime && time+((double(idxWvfmEntry)+0.5)*deltaT) < orderedTimes[iOverlap] ) {
	      adcVec.push_back( overlayOpWave.baseline + (wvfm[idxWvfmEntry] - simBaseline) );
	      idxWvfmEntry+=1;
	    }

	    opWaveformOverlay thisOverlap = opWaveformsMap[ chID ].at( orderedOverlaps[iOverlap] );

	    // If after, then start putting in entries before the overlapping time
	    unsigned int alreadyInOverlap = 0;
	    while ( idxWvfmEntry<wvfm.size() && orderedTimes[iOverlap]+((double(adcVec.size())+0.5)*deltaT) < time ) {
	      adcVec.push_back( thisOverlap.wvfm[ adcVec.size() ] );
	      alreadyInOverlap+=1;
	    }

	    for ( unsigned int idxInOverlap=alreadyInOverlap; idxInOverlap<thisOverlap.wvfm.size(); ++idxInOverlap ) {
	      raw::ADC_Count_t addTo = ( idxWvfmEntry < wvfm.size() ? (wvfm[idxWvfmEntry] - simBaseline) : 0 );
	      adcVec.push_back( thisOverlap.wvfm[idxInOverlap] + addTo );
	      idxWvfmEntry+=1;
	    }

	    loTime = thisOverlap.endTime();
	  } // loop overlaps

	  // and fill in any remaining bins left
	  while ( idxWvfmEntry < wvfm.size() ) {
	    adcVec.push_back( overlayOpWave.baseline + (wvfm[idxWvfmEntry] - simBaseline) );
	    idxWvfmEntry+=1;
	  }

	  // Now delete the overlapping waveforms
	  for ( unsigned int iOverlap=0; iOverlap<overlaps.size(); ++iOverlap ) {
	    unsigned int positionToErase = overlaps[ overlaps.size()-1-iOverlap ];

	    if ( int(chID)==fPMTWaveTestCh && wvfm.size()>10500 ) {
	      std::cout << "... We should now erase entry " << positionToErase << " which is STRUCT WITH:" << std::endl;
	      std::cout << "    time:     " << opWaveformsMap[ chID ].at( positionToErase+1 ).timestamp << std::endl;
	      std::cout << "    baseline: " << opWaveformsMap[ chID ].at( positionToErase+1 ).baseline << std::endl;
	      std::cout << "    wv size:  " << opWaveformsMap[ chID ].at( positionToErase+1 ).wvfm.size() << std::endl;
	      std::cout << "    end time: " << opWaveformsMap[ chID ].at( positionToErase+1 ).endTime() << std::endl;
	    }

	    opWaveformsMap[ chID ].erase( opWaveformsMap[ chID ].begin()+positionToErase );
	  }

	  // And save this new one
	  overlayOpWave.wvfm = adcVec;
	  opWaveformsMap[ chID ].push_back( overlayOpWave );
	} // deal with overlaps
      } // look for overlaps
    } // loop simulated waveforms

    // Now put all the waveforms into the collections
    art::PtrMaker<raw::OpDetWaveform> opdetwvfmPtrMaker(e);
    art::PtrMaker<icarus::WaveformBaseline> baselinePtrMaker(e);
    for ( auto const &[chID, wvfmStructVec] : opWaveformsMap ) {
      for ( auto const& wvfmStruct : wvfmStructVec ) {
	const raw::TimeStamp_t time = wvfmStruct.timestamp;
	const raw::Channel_t chan = chID;
	std::vector< raw::ADC_Count_t > vals = wvfmStruct.wvfm;
	const icarus::WaveformBaseline::Baseline_t base =  wvfmStruct.baseline;

	// Print
	if ( int(chID)==fPMTWaveTestCh ) {
	  std::cout << "Saving STRUCT WITH: " << std::endl;
	  std::cout << "    time:     " << wvfmStruct.timestamp << std::endl;
	  std::cout << "    baseline: " << wvfmStruct.baseline << std::endl;
	  std::cout << "    wv size:  " << wvfmStruct.wvfm.size() << std::endl;
	  std::cout << "    end time: " << wvfmStruct.endTime() << std::endl;
	}

	if ( int(chID)==fPMTWaveTestCh && vals.size()>10500 ) {
	  std::cout << "OVERLAID WAVEFORM AT TIME " << std::setprecision(8) << time << ": [ ";
	  for ( unsigned int i=0;i<vals.size(); ++i ) {
	    std::cout << vals[i];
	    if ( i==vals.size()-1 ) std::cout<< " ]" << std::endl;
	    else std::cout << ", ";
	  }
	}

	// TODO: why does constructor require unsigned short but then use type short???
	std::vector<uint16_t> valsUnsigned;
	for ( auto const& val : vals ) valsUnsigned.push_back( val );

	raw::OpDetWaveform        odw( time, chan, valsUnsigned );
	icarus::WaveformBaseline  bLine( base );

	opWaveformVec->push_back( odw );
	baselineVec->push_back( bLine );
	baselineToWaveforms.addSingle( baselinePtrMaker(baselineVec->size()-1), opdetwvfmPtrMaker(opWaveformVec->size()-1) );
      }
    }

    // ... and put collections into the event
    e.put( std::move( opWaveformVec ) );
    e.put( std::move( baselineVec ) );
    e.put( std::move( std::make_unique< art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform> >(std::move(baselineToWaveforms)) ) );
  }

  if ( fPMTOverlayHits ) {
    std::unique_ptr< std::vector< recob::OpHit > > opHitVec = std::make_unique< std::vector< recob::OpHit > >();

    for ( auto const& iLabel : fPMTHitInputLabels ) {
      art::Handle< std::vector<recob::OpHit> > hitsHandle;
      std::vector< art::Ptr<recob::OpHit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
	art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
	mf::LogWarning("OverlayProducts") << "Event failed to find recob::OpHit with label " << iLabel << ".";
        return;
      }

      for ( auto const& iHit : hits ) {
	recob::OpHit newOpHit( iHit->OpChannel(), iHit->PeakTime(), iHit->PeakTimeAbs(), iHit->Frame(), iHit->Width(), iHit->Area(), iHit->Amplitude(), iHit->PE(), iHit->FastToTotal() );
	opHitVec->push_back( newOpHit );
      } // hit
    } // label

    e.put( std::move(opHitVec) );
  }

  // CRT Overlays
  if ( fCRTOverlayHits ) {
    std::unique_ptr< std::vector< sbn::crt::CRTHit > > crtHitVec = std::make_unique< std::vector< sbn::crt::CRTHit > >();

    for ( auto const& iLabel : fCRTHitInputLabels ) {
      art::Handle< std::vector<sbn::crt::CRTHit> > hitsHandle;
      std::vector< art::Ptr<sbn::crt::CRTHit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
	art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
	mf::LogWarning("OverlayProducts") << "Event failed to find sbn::crt::CRTHit with label " << iLabel << ".";
        return;
      }

      for ( auto const& iHit : hits ) {
	sbn::crt::CRTHit newCRTHit;
	newCRTHit.feb_id      = iHit->feb_id;
	newCRTHit.pesmap      = iHit->pesmap;
	newCRTHit.peshit      = iHit->peshit;
	newCRTHit.ts0_s       = iHit->ts0_s;
	newCRTHit.ts0_s_corr  = iHit->ts0_s_corr;
	newCRTHit.ts0_ns      = iHit->ts0_ns;
	newCRTHit.ts0_ns_corr = iHit->ts0_ns_corr;
	newCRTHit.ts1_ns      = iHit->ts1_ns;
	newCRTHit.plane       = iHit->plane;
	newCRTHit.x_pos       = iHit->x_pos;
	newCRTHit.x_err       = iHit->x_err;
	newCRTHit.y_pos       = iHit->y_pos;
	newCRTHit.y_err       = iHit->y_err;
	newCRTHit.z_pos       = iHit->z_pos;
	newCRTHit.z_err       = iHit->z_err;
	newCRTHit.tagger      = iHit->tagger;

        crtHitVec->push_back( newCRTHit );
      } // hit
    } // label

    e.put( std::move(crtHitVec) );
  }
}

DEFINE_ART_MODULE(OverlayProducts)