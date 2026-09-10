/**
 * @file   sbndcode/TPC/NuGraph/SBNDHDF5Maker_module.cc
 * @brief  Implementation of `SBNDHDF5Maker` _art_ module.
 * @author Alex Antonakis (aantonakis@ucsb.edu) and Giuseppe Cerati (cerati@fnal.gov), V Hewes
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardataobj/RecoBase/Slice.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "hep_hpc/hdf5/make_ntuple.hpp"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include <set>
#include <optional>

#include "StitchingUtils.h"

class SBNDHDF5Maker : public art::EDAnalyzer {
public:
  explicit SBNDHDF5Maker(fhicl::ParameterSet const& p);

  SBNDHDF5Maker(SBNDHDF5Maker const&) = delete;
  SBNDHDF5Maker(SBNDHDF5Maker&&) = delete;
  SBNDHDF5Maker& operator=(SBNDHDF5Maker const&) = delete;
  SBNDHDF5Maker& operator=(SBNDHDF5Maker&&) = delete;

  void beginSubRun(art::SubRun const& sr) override;
  void endSubRun(art::SubRun const& sr) override;
  void analyze(art::Event const& e) override;

private:

  art::InputTag fTruthLabel;
  art::InputTag fHitLabel;
  art::InputTag fHitTruthLabel;
  art::InputTag fSPLabel;
  art::InputTag fOpHitLabel;
  art::InputTag fOpFlashLabel;

  bool fUseMap;
  std::string fOutputName;

  struct HDFDataFile {
    hep_hpc::hdf5::File file;  ///< output HDF5 file

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
                          hep_hpc::hdf5::Column<int, 1>,    // is cc
			  hep_hpc::hdf5::Column<int, 1>, // nu pdg
			  hep_hpc::hdf5::Column<float, 1>,  // nu energy
			  hep_hpc::hdf5::Column<float, 1>,  // lep energy
			  hep_hpc::hdf5::Column<float, 1>   // nu dir (x, y, z)
    > eventNtupleNu; ///< event ntuple with neutrino information

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
			  hep_hpc::hdf5::Column<int, 1>,    // spacepoint id
			  hep_hpc::hdf5::Column<float, 1>,  // 3d position (x, y, z)
			  hep_hpc::hdf5::Column<int, 1>,     // 2d hit (u, v, y)
			  hep_hpc::hdf5::Column<float, 1>     //ChiSquared of the hit 
    > spacePointNtuple; ///< spacepoint ntuple

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    			  hep_hpc::hdf5::Column<int, 1>,    // hit id
    			  hep_hpc::hdf5::Column<float, 1>,  // hit integral
    			  hep_hpc::hdf5::Column<float, 1>,  // hit rms
    			  hep_hpc::hdf5::Column<int, 1>,    // tpc id
    			  hep_hpc::hdf5::Column<int, 1>,    // global plane
    			  hep_hpc::hdf5::Column<float, 1>,  // global wire
    			  hep_hpc::hdf5::Column<float, 1>,  // global time
    			  hep_hpc::hdf5::Column<int, 1>,    // raw plane
    			  hep_hpc::hdf5::Column<float, 1>,  // raw wire
    			  hep_hpc::hdf5::Column<float, 1>,  // raw time
    			  hep_hpc::hdf5::Column<int, 1>    //cryostat
    > hitNtuple; ///< hit ntuple

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    			  hep_hpc::hdf5::Column<int, 1>,    // g4 id
    			  hep_hpc::hdf5::Column<int, 1>,    // particle type
    			  hep_hpc::hdf5::Column<int, 1>,    // parent g4 id
    			  hep_hpc::hdf5::Column<int, 1>,    // is from nu
    			  hep_hpc::hdf5::Column<float, 1>,  // momentum
    			  hep_hpc::hdf5::Column<float, 1>,  // start position (x, y, z)
    			  hep_hpc::hdf5::Column<float, 1>,  // end position (x, y, z)
    			  hep_hpc::hdf5::Column<std::string, 1>, // start process
    			  hep_hpc::hdf5::Column<std::string, 1>  // end process
    > particleNtuple; ///< particle ntuple

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    			  hep_hpc::hdf5::Column<int, 1>,    // hit id
    			  hep_hpc::hdf5::Column<int, 1>,    // g4 id
    			  hep_hpc::hdf5::Column<float, 1>,  // deposited energy [ MeV ]
    			  hep_hpc::hdf5::Column<float, 1>,  // x position
    			  hep_hpc::hdf5::Column<float, 1>,  // y position
    			  hep_hpc::hdf5::Column<float, 1>   // z position
    > energyDepNtuple; ///< energy deposition ntuple

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    			  hep_hpc::hdf5::Column<int, 1>,    // hit id
    			  hep_hpc::hdf5::Column<int, 1>,    // hit_channel
    			  hep_hpc::hdf5::Column<int, 1>,    // wire pos
    			  hep_hpc::hdf5::Column<float, 1>,  // peaktime
    			  hep_hpc::hdf5::Column<float, 1>,  // width
    			  hep_hpc::hdf5::Column<float, 1>,  // area
    			  hep_hpc::hdf5::Column<float, 1>,  // amplitude
    			  hep_hpc::hdf5::Column<float, 1>,  // pe
    			  hep_hpc::hdf5::Column<int, 1>     // usedInFlash
    > opHitNtuple; ///< PMT hit ntuple

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    			  hep_hpc::hdf5::Column<int, 1>,    // flash id
                          hep_hpc::hdf5::Column<int, 1>,    // wire pos
    			  hep_hpc::hdf5::Column<float, 1>,  // time
    			  hep_hpc::hdf5::Column<float, 1>,  // time width
    			  hep_hpc::hdf5::Column<float, 1>,  // Y center
    			  hep_hpc::hdf5::Column<float, 1>,  // Y width
    			  hep_hpc::hdf5::Column<float, 1>,  // Z center
    			  hep_hpc::hdf5::Column<float, 1>,  // Z width
    			  hep_hpc::hdf5::Column<float, 1>   // totalpe
    > opFlashNtuple; ///< Flash ntuple

    hep_hpc::hdf5::Ntuple<hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    			  hep_hpc::hdf5::Column<int, 1>,    // sumpe id
    			  hep_hpc::hdf5::Column<int, 1>,    // flash id
    			  hep_hpc::hdf5::Column<int, 1>,    // PMT channel
                          hep_hpc::hdf5::Column<float, 1>,  // YZ pos
    			  hep_hpc::hdf5::Column<float, 1>   // pe
    > opFlashSumPENtuple; ///< Flash SumPE ntuple

    static std::string makeOutputFileName(std::string const& outputName, art::SubRunID const& sr)
    {
      struct timeval now;
      gettimeofday(&now, NULL);
      std::ostringstream fileName;
      fileName << outputName
	       << "_r" << std::setfill('0') << std::setw(5) << sr.run()
	       << "_s" << std::setfill('0') << std::setw(5) << sr.subRun()
	       << "_ts" << std::setw(6) << now.tv_usec << ".h5";
      std::cout << fileName.str() << std::endl;
      return fileName.str();
    }
    
    HDFDataFile(std::string const& outputName, art::SubRunID const& sr)
      : file{ makeOutputFileName(outputName, sr), H5F_ACC_TRUNC }
      , eventNtupleNu{
      	hep_hpc::hdf5::make_ntuple({file, "event_table", 1000},
      	hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("is_cc"),
        hep_hpc::hdf5::make_scalar_column<int>("nu_pdg"),    
        hep_hpc::hdf5::make_scalar_column<float>("nu_energy"),
        hep_hpc::hdf5::make_scalar_column<float>("lep_energy"),
        hep_hpc::hdf5::make_column<float>("nu_dir", 3))
        }
      , spacePointNtuple{
        hep_hpc::hdf5::make_ntuple({file, "spacepoint_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("spacepoint_id"),
        hep_hpc::hdf5::make_column<float>("position", 3),
        hep_hpc::hdf5::make_column<int>("hit_id", 3),
        hep_hpc::hdf5::make_scalar_column<float>("Chi_squared"))
        }
      , hitNtuple{
        hep_hpc::hdf5::make_ntuple({file, "hit_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
        hep_hpc::hdf5::make_scalar_column<float>("integral"),
        hep_hpc::hdf5::make_scalar_column<float>("rms"),
        hep_hpc::hdf5::make_scalar_column<int>("tpc"),
        hep_hpc::hdf5::make_scalar_column<int>("global_plane"),
        hep_hpc::hdf5::make_scalar_column<float>("global_wire"),
        hep_hpc::hdf5::make_scalar_column<float>("global_time"),
        hep_hpc::hdf5::make_scalar_column<int>("local_plane"),
        hep_hpc::hdf5::make_scalar_column<float>("local_wire"),
        hep_hpc::hdf5::make_scalar_column<float>("local_time"),
      	hep_hpc::hdf5::make_scalar_column<int>("Cryostat"))
        }
      , particleNtuple{
        hep_hpc::hdf5::make_ntuple({file, "particle_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("g4_id"),
        hep_hpc::hdf5::make_scalar_column<int>("type"),
        hep_hpc::hdf5::make_scalar_column<int>("parent_id"),
        hep_hpc::hdf5::make_scalar_column<int>("from_nu"),
        hep_hpc::hdf5::make_scalar_column<float>("momentum"),
        hep_hpc::hdf5::make_column<float>("start_position", 3),
        hep_hpc::hdf5::make_column<float>("end_position", 3),
        hep_hpc::hdf5::make_scalar_column<std::string>("start_process"),
      	hep_hpc::hdf5::make_scalar_column<std::string>("end_process"))
        }
      , energyDepNtuple{
        hep_hpc::hdf5::make_ntuple({file, "edep_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
        hep_hpc::hdf5::make_scalar_column<int>("g4_id"),
        hep_hpc::hdf5::make_scalar_column<float>("energy"),
        hep_hpc::hdf5::make_scalar_column<float>("x_position"),
        hep_hpc::hdf5::make_scalar_column<float>("y_position"),
      	hep_hpc::hdf5::make_scalar_column<float>("z_position"))
        }
      , opHitNtuple{
        hep_hpc::hdf5::make_ntuple({file, "ophit_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("hit_id"),
        hep_hpc::hdf5::make_scalar_column<int>("hit_channel"),
        hep_hpc::hdf5::make_column<int>("wire_pos", 3),//3 views
        hep_hpc::hdf5::make_scalar_column<float>("peaktime"),
        hep_hpc::hdf5::make_scalar_column<float>("width"),
        hep_hpc::hdf5::make_scalar_column<float>("area"),
        hep_hpc::hdf5:: make_scalar_column<float>("amplitude"),
        hep_hpc::hdf5::make_scalar_column<float>("pe"),
      	hep_hpc::hdf5::make_scalar_column<int>("sumpe_id"))
        }
      , opFlashNtuple{
        hep_hpc::hdf5::make_ntuple({file, "opflash_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("flash_id"),
        hep_hpc::hdf5::make_column<int>("wire_pos", 3),//3 views
        hep_hpc::hdf5::make_scalar_column<float>("time"),
        hep_hpc::hdf5::make_scalar_column<float>("time_width"),
        hep_hpc::hdf5::make_scalar_column<float>("y_center"),
        hep_hpc::hdf5::make_scalar_column<float>("y_width"),
        hep_hpc::hdf5::make_scalar_column<float>("z_center"),
        hep_hpc::hdf5::make_scalar_column<float>("z_width"),
      	hep_hpc::hdf5::make_scalar_column<float>("totalpe"))
        }
      , opFlashSumPENtuple{
        hep_hpc::hdf5::make_ntuple({file, "opflashsumpe_table", 1000},
        hep_hpc::hdf5::make_column<int>("event_id", 3),
        hep_hpc::hdf5::make_scalar_column<int>("sumpe_id"),
        hep_hpc::hdf5::make_scalar_column<int>("flash_id"),
        hep_hpc::hdf5::make_scalar_column<int>("pmt_channel"),
        hep_hpc::hdf5::make_column<float>("yz_pos", 2),
      	hep_hpc::hdf5::make_scalar_column<float>("sumpe"))
        }
      { }
  };
  std::unique_ptr<HDFDataFile> fHDFData;

  int NearWire(const geo::PlaneGeo &plane, float x, float y, float z) const;
};

SBNDHDF5Maker::SBNDHDF5Maker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTruthLabel(p.get<art::InputTag>("TruthLabel")),
    fHitLabel(  p.get<art::InputTag>("HitLabel")),
    fHitTruthLabel(  p.get<art::InputTag>("HitTruthLabel","")),
    fSPLabel(   p.get<art::InputTag>("SPLabel")),
    fOpHitLabel(  p.get<art::InputTag>("OpHitLabel")),
    fOpFlashLabel(  p.get<art::InputTag>("OpFlashLabel")),
    fUseMap(    p.get<bool>("UseMap", false)),
    fOutputName(p.get<std::string>("OutputName"))
{ }

void SBNDHDF5Maker::analyze(art::Event const& e) {

  const cheat::BackTrackerService* bt = fUseMap? nullptr: art::ServiceHandle<cheat::BackTrackerService>().get();
  // TODO --> Not needed for opflash?
  //geo::WireReadoutGeom const& geom = art::ServiceHandle<geo::WireReadout const>()->Get();

  std::set<int> g4id;
  auto const* pi = lar::providerFrom<cheat::ParticleInventoryService>();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);
  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  // auto tpcgeo = art::ServiceHandle<geo::Geometry const>()->TPC();
  // std::cout << "tpc drift distance=" << tpcgeo.DriftDistance() << std::endl;
  // std::cout << "drift velocity=" << detProp.DriftVelocity() << std::endl;
  // std::cout << "drift time=" << tpcgeo.DriftDistance()/detProp.DriftVelocity() << " 2*dt/0.4=" << 2*tpcgeo.DriftDistance()/detProp.DriftVelocity()/clockData.TPCClock().TickPeriod() << std::endl;
  // std::cout << "tick period=" << clockData.TPCClock().TickPeriod() << std::endl;
  // std::cout << "tpcTime=" << clockData.TPCTime() << std::endl;
  // std::cout << "triggerOffsetTPC=" << clockData.TriggerOffsetTPC() << std::endl;
  // std::cout << "triggerTime=" << clockData.TriggerTime() << std::endl;
  // std::cout << "correction=" << 2*(tpcgeo.DriftDistance()/detProp.DriftVelocity()-clockData.TriggerOffsetTPC())/clockData.TPCClock().TickPeriod() << std::endl;

  // std::cout << "NEW EVENT: " << run <<"  "<< subrun << " "<< event<<std::endl;

  std::array<int, 3> evtID { run, subrun, event };
  // Fill event table
  // Get MC truth
  auto truthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fTruthLabel);
  std::cout << "truthHandle->size()=" << truthHandle->size() << std::endl;
  if (truthHandle->size() != 1) {
    std::cout << "truthHandle->size()=" << truthHandle->size() << "  skipping event" << std::endl;
    //avoid pile-up, which is not handled downstream
    return;
  }
  simb::MCNeutrino const& nutruth = truthHandle->at(0).GetNeutrino();

  auto up = nutruth.Nu().Momentum().Vect().Unit();
  std::array<float, 3> nuMomentum {(float)up.X(),(float)up.Y(),(float)up.Z()};
  std::cout << "Filling Neutrino info: "
      << "is cc? " << (nutruth.CCNC() == simb::kCC)
      << ", nu energy " << nutruth.Nu().E()
      << ", lepton energy " << nutruth.Lepton().E()
      << "\nnu momentum x " << nuMomentum[0] << ", y "
      << nuMomentum[1] << ", z " << nuMomentum[2]
      << std::endl;
  fHDFData->eventNtupleNu.insert( evtID.data(),
				  nutruth.CCNC() == simb::kCC,
				  nutruth.Nu().PdgCode(),
				  nutruth.Nu().E(),
				  nutruth.Lepton().E(),
				  nuMomentum.data()
				  );

  // for (int ip=0;ip<truthHandle->at(0).NParticles();ip++) {
  //   std::cout << "mcp tkid=" << truthHandle->at(0).GetParticle(ip).TrackId() << " pdg=" << truthHandle->at(0).GetParticle(ip).PdgCode() 
  // 		<< " mother=" << truthHandle->at(0).GetParticle(ip).Mother()
  // 		<< " vtx=" << truthHandle->at(0).GetParticle(ip).Vx() << " " << truthHandle->at(0).GetParticle(ip).Vy() << " " << truthHandle->at(0).GetParticle(ip).Vz()
  // 		<< std::endl;
  // }

  mf::LogDebug("SBNDHDF5Maker") << "Filling event table"
				  << "\nrun " << evtID[0] << ", subrun " << evtID[1]
				  << ", event " << evtID[2]
				  << "\nis cc? " << (nutruth.CCNC() == simb::kCC)
				  << ", nu energy " << nutruth.Nu().E()
				  << ", lepton energy " << nutruth.Lepton().E()
				  << "\nnu momentum x " << nuMomentum[0] << ", y "
				  << nuMomentum[1] << ", z " << nuMomentum[2];

  std::cout << "Grabbing SpacePoints" << std::endl;
  // Get spacepoints from the event record
  auto spListHandle = e.getValidHandle<std::vector<recob::SpacePoint>>(fSPLabel);
  std::vector<recob::SpacePoint> const& splist = *spListHandle;

  std::cout << "Grabbing Hits" << std::endl;
  // Get hits from the event record
  auto hitListHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
  std::vector<art::Ptr<recob::Hit>> hitlist;
  art::fill_ptr_vector(hitlist, hitListHandle);

  std::cout << "Grabbing Associations from SpacePoints to Hits" << std::endl;
  // Get assocations from spacepoints to hits
  art::FindManyP<recob::Hit> fmp(spListHandle, e, fSPLabel);

  std::cout << "Filling the SpacePoint table" << std::endl;
  // Fill spacepoint table
  for (size_t i = 0; i < spListHandle->size(); ++i) {

    std::vector<art::Ptr<recob::Hit>> const& associatedHits = fmp.at(i);
    //if (associatedHits.size() < 3) {
    //  throw art::Exception(art::errors::LogicError) << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits";
    //}

    std::array<float, 3> pos {(float)splist[i].XYZ()[0],(float)splist[i].XYZ()[1],(float)splist[i].XYZ()[2]};

    std::array<int, 3> hitID { -1, -1, -1 };
    for (art::Ptr<recob::Hit> const& hit: associatedHits) {
      int plane = util::stitchedPlane(hit->WireID());
      //std::cout << " tpc=" << hit->WireID().TPC << " plane=" << hit->WireID().Plane
      //          << " plane2=" << plane << " key=" <<  hit.key() << std::endl;
      hitID[plane] = hit.key();
    }
    fHDFData->spacePointNtuple.insert(evtID.data(),splist[i].ID(), pos.data(), hitID.data(),splist[i].Chisq() );

    mf::LogDebug("SBNDHDF5Maker") << "Filling spacepoint table"
				    << "\nrun " << evtID[0] << ", subrun " << evtID[1]
				    << ", event " << evtID[2]
				    << "\nspacepoint id " << splist[i].ID()
				    << "\nposition x " << pos[0] << ", y " << pos[1]
				    << ", z " << pos[2]
				    << "\nhit ids " << hitID[0] << ", " << hitID[1]
				    << ", " << hitID[2]
				    << "Chi_squared "  << splist[i].Chisq();
  } // for spacepoint

  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (fUseMap) {
    hittruth = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> >(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(hitListHandle, e, fHitTruthLabel));
  }

  std::cout << "Looping over hits" << std::endl;
  // Loop over hits
  for (art::Ptr<recob::Hit> hit : hitlist) {

    // Fill hit table
    geo::WireID wireid = hit->WireID();
    size_t plane = util::stitchedPlane(wireid);
    double time = util::stitchedTime(wireid,hit->PeakTime());
    size_t wire = util::stitchedWire(wireid);
    fHDFData->hitNtuple.insert(evtID.data(),
			       hit.key(), hit->Integral(), hit->RMS(), wireid.TPC,
			       plane, wire, time,
			       wireid.Plane, wireid.Wire, hit->PeakTime(),wireid.Cryostat
			       );

    mf::LogInfo("SBNDHDF5Maker") << "Filling hit table"
                             << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                             << ", event " << evtID[2]
                             << "\nhit id " << hit.key() << ", integral "
                             << hit->Integral() << ", RMS " << hit->RMS()
                             << ", TPC " << wireid.TPC
                             << "\nglobal plane " << plane << ", global wire "
                             << wire << ", global time " << time
                             << "\nlocal plane " << wireid.Plane
                             << ", local wire " << wireid.Wire
                             << ", local time " << hit->PeakTime()
                             << ", Cryostat " << wireid.Cryostat;

    // Fill energy deposit table
    if (fUseMap) {

    // Use the BackTracker/Hit truth map to store truth info per hit (microboone-style).
    if (!hittruth) {
      mf::LogWarning("SBNDHDF5Maker") << "Hit truth map requested but not available for label " << fHitTruthLabel;
    } else {
      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec = hittruth->data(hit.key());
      // loop over matched particles and record energy-fraction info
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
        int tid = particle_vec[i_p]->TrackId();
        g4id.insert(tid);
        float ideFrac = (match_vec[i_p] ? static_cast<float>(match_vec[i_p]->ideFraction) : 0.f);
        fHDFData->energyDepNtuple.insert(evtID.data(), hit.key(), tid, ideFrac, -1., -1., -1.);
        mf::LogInfo("SBNDHDF5Maker") << "Filling energy deposit table"
                 << "\nrun " << evtID[0] << ", subrun " << evtID[1]
                 << ", event " << evtID[2]
                 << "\nhit id " << hit.key() << ", g4 id "
                 << tid << ", energy fraction " << ideFrac << ", position (nan, nan, nan)";
      } // for particle
    }
    }  else {

      std::vector<sim::TrackIDE> const& h_ides = bt->ChannelToTrackIDEs(clockData, hit->Channel(), hit->StartTick(), hit->EndTick());
      for (auto const& tide : h_ides) {
      	int tid = tide.trackID;
	// if negative id, make sure there is a corresponding particle to look for before taking the abs.
	// This way negative means no associated particle (a convention that can be used in the code that processes the ntuple).
      	if (pi->TrackIdToParticle_P(abs(tid))) tid = abs(tid);
        g4id.insert(tid);
        fHDFData->energyDepNtuple.insert(evtID.data(),hit.key(), tid, tide.energyFrac, -1., -1., -1.);
      }

    } // if using microboone map method or not
  } // for hit

  std::unordered_map<int, simb::MCParticle const*> allIDs;

  // Add invisible particles to hierarchy
  for (int id : g4id) {
    const simb::MCParticle* p = pi->TrackIdToParticle_P(abs(id));
    allIDs.emplace(abs(id), p);
    while (p && p->Mother() != 0 ) {
      auto mid = abs(p->Mother());
      p = pi->TrackIdToParticle_P(mid);
      allIDs.emplace(mid, p);
    }
  }
  // Loop over true particles and fill table
  for (auto [ id, p ] : allIDs) {
    if (!p) {
      mf::LogError("SBNDHDF5Maker") << "Track ID=" << id << " not tracked back to any particle.";
      continue;
    }
    auto mct = pi->TrackIdToMCTruth_P(abs(id));
    int fromNu = (mct.isAvailable() ? mct->NeutrinoSet() : 0);
    // std::cout << "all mcp tkid=" << p->TrackId() << " pdg=" << p->PdgCode()
    // 	      << " mother=" << p->Mother()
    // 	      << " vtx=" << p->Vx() << " " << p->Vy() << " " << p->Vz()
    // 	      << " mctruth=" << mct->NeutrinoSet()
    // 	      << std::endl;
    std::array<float, 3> particleStart { (float)p->Vx(), (float)p->Vy(), (float)p->Vz() };
    std::array<float, 3> particleEnd { (float)p->EndX(), (float)p->EndY(), (float)p->EndZ() };
    fHDFData->particleNtuple.insert(evtID.data(),
				    abs(id), p->PdgCode(), p->Mother(), fromNu, (float)p->P(),
				    particleStart.data(), particleEnd.data(),
				    p->Process(), p->EndProcess()
				    );
    mf::LogDebug("SBNDHDF5Maker") << "Filling particle table"
				    << "\nrun " << evtID[0] << ", subrun " << evtID[1]
				    << ", event " << evtID[2]
				    << "\ng4 id " << abs(id) << ", pdg code "
				    << p->PdgCode() << ", parent " << p->Mother()
				    << ", momentum " << p->P()
				    << "\nparticle start x " << particleStart[0]
				    << ", y " << particleStart[1]
				    << ", z " << particleStart[2]
				    << "\nparticle end x " << particleEnd[0] << ", y "
				    << particleEnd[1] << ", z " << particleEnd[2]
				    << "\nstart process " << p->Process()
				    << ", end process " << p->EndProcess();
  }

  /*
  std::cout << "Grabbing OpFlashs" << std::endl;
  // Get OpFlashs from the event record
  art::Handle< std::vector< recob::OpFlash > > opFlashListHandle;
  std::optional<art::FindManyP<recob::OpHit> > assocFlashHit;
  std::vector<recob::OpFlash> const* opflashlist = nullptr;
  if (!fOpFlashLabel.empty()) {
    e.getByLabel(fOpFlashLabel, opFlashListHandle);
    opflashlist = opFlashListHandle.product();
    assocFlashHit.emplace(opFlashListHandle, e, fOpFlashLabel);
  }

  std::cout << "Grabbing Slice to OpFlash associations" << std::endl;
  // get the flash matching the slice we selected
  auto const& sliceOpFlashAssns = e.getProduct<art::Assns<recob::Slice, recob::OpFlash>>(fHitLabel);
  int flkey = sliceOpFlashAssns.size()==0 ? -1: sliceOpFlashAssns.at(0).second.key();

  std::vector<std::vector<int> > flashsumpepmtmap;//this has at most one element, due to the check on flkey

  // Pick only the flash we care about
  if (flkey>=0) {
    recob::OpFlash const& opflash = opflashlist->at(flkey);
    std::vector<float> pes;
    for (auto pe : opflash.PEs()) pes.push_back(pe);
    std::vector<int> nearwires;
    double xyz[3] = {0.,opflash.YCenter(),opflash.ZCenter()};
    for (auto const& p : geom.Iterate<geo::PlaneGeo>()) nearwires.push_back(NearWire(p,xyz[0],xyz[1],xyz[2]));

    // Fill opflash table
    fHDFData->opFlashNtuple.insert(evtID.data(),
				   flkey,
				   nearwires.data(),
				   opflash.Time(),opflash.TimeWidth(),
				   opflash.YCenter(),opflash.YWidth(),opflash.ZCenter(),opflash.ZWidth(),
				   opflash.TotalPE()
				   );
    size_t count = 0;
    std::vector<int> sumpepmtmap(art::ServiceHandle<geo::Geometry>()->NOpDets(),0);
    for (size_t ipmt=0;ipmt<pes.size();ipmt++) {
      if (pes[ipmt]<=0.) continue;
      auto xyz = geom.OpDetGeoFromOpChannel(ipmt).GetCenter();
      std::vector<float> yzpos = {float(xyz.Y()),float(xyz.Z())};
      fHDFData->opFlashSumPENtuple.insert(evtID.data(),count,flkey,ipmt,yzpos.data(),pes[ipmt]);
      sumpepmtmap[ipmt] = count;
      count++;
    }
    flashsumpepmtmap.push_back(sumpepmtmap);
    mf::LogDebug("SBNDHDF5Maker") << "Filling opflash table"
				    << "\nrun " << evtID[0] << ", subrun " << evtID[1]
				    << ", event " << evtID[2]
				    << "\nflash id " << flkey << ", Time " << opflash.Time()
				    << ", TotalPE " << opflash.TotalPE()//;
				    << "\nWireCenters size0 " << opflash.WireCenters().size()//;
				    << "\nYCenter " << opflash.YCenter()<< " ZCenter " << opflash.ZCenter()
				    << "\nYWidth " << opflash.YWidth()<< " ZWidth " << opflash.ZWidth()
				    << "\nInBeamFrame " << opflash.InBeamFrame()<< " OnBeamTime " << opflash.OnBeamTime()
				    << "\nAbsTime " << opflash.AbsTime() << " TimeWidth " << opflash.TimeWidth() << " FastToTotal " << opflash.FastToTotal();
  }

  std::cout << "Grabbing OpHits" << std::endl;
  // Get OpHits from the event record
  std::vector< art::Ptr< recob::OpHit > > ophitlist;
  if (!fOpHitLabel.empty()) {
    auto opHitListHandle = e.getValidHandle< std::vector< recob::OpHit > >(fOpHitLabel);
    art::fill_ptr_vector(ophitlist, opHitListHandle);
  }

  std::set<art::Ptr<recob::OpHit>> ophv;
  if (flkey >= 0) {
    auto const& flashHits = assocFlashHit->at(flkey);
    ophv.insert(begin(flashHits), end(flashHits));
  }
  // Loop over ophits
  for (art::Ptr< recob::OpHit > ophit : ophitlist) {
    std::vector<int> nearwires;
    auto xyz = geom.OpDetGeoFromOpChannel(ophit->OpChannel()).GetCenter();
    for (auto const& p : geom.Iterate<geo::PlaneGeo>()) nearwires.push_back(NearWire(p,xyz.X(),xyz.Y(),xyz.Z()));

    bool isInFlash = (ophv.count(ophit) > 0);

    // Fill ophit table
    fHDFData->opHitNtuple.insert(evtID.data(),ophit.key(),
				 ophit->OpChannel(),nearwires.data(),
				 ophit->PeakTime(),ophit->Width(),
				 ophit->Area(),ophit->Amplitude(),ophit->PE(),
				 (isInFlash ? flashsumpepmtmap[0][ophit->OpChannel()] : -1)
				 );
    mf::LogDebug("SBNDHDF5Maker") << "\nFilling ophit table"
				    << "\nrun " << evtID[0] << ", subrun " << evtID[1]
				    << ", event " << evtID[2]
				    << "\nhit id " << ophit.key() << ", channel "
				    << ophit->OpChannel() << ", PeakTime " << ophit->PeakTime()
				    << ", Width " << ophit->Width()
				    << "\narea " << ophit->Area() << ", amplitude "
				    << ophit->Amplitude() << ", PE " << ophit->PE();
  }
  */
  //End optical analyzer

} // function SBNDHDF5Maker::analyze

void SBNDHDF5Maker::beginSubRun(art::SubRun const& sr) {
  fHDFData = std::make_unique<HDFDataFile>(fOutputName, sr.id());
}

void SBNDHDF5Maker::endSubRun(art::SubRun const& sr) {
  fHDFData.reset();
}

int SBNDHDF5Maker::NearWire(const geo::PlaneGeo &plane, float x, float y, float z) const
{
  geo::WireID wireID;
  try {
    wireID = plane.NearestWireID(geo::Point_t(x,y,z));
  }
  catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane.ClosestWireID(e.suggestedWireID());
  }
  return wireID.Wire;
}

DEFINE_ART_MODULE(SBNDHDF5Maker)
