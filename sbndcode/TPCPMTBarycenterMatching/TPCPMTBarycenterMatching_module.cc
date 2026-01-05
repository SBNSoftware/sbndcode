////////////////////////////////////////////////////////////////////////
// Class:       TPCPMTBarycenterMatchProducer
// Plugin Type: producer (Unknown Unknown)
// File:        TPCPMTBarycenterMatchProducer_module.cc
//
// Generated at Sun Oct 22 14:43:16 2023 by John Smedley using cetskelgen
// from  version .
//
//  @file   icaruscode/PMT/OpReco/TPCPMTBarycenterMatchProducer_module.cc
//  @brief  Producer to match Pandora slices to their best match OpFlash by minimizing barycenter distance, as well as compare slices to the triggering OpFlash
//  @author Jack Smedley ( jsmedley@fnal.gov )
//
////////////////////////////////////////////////////////////////////////

// Ported to sbndcode by Alejandro Sánchez Castillo (asanchezcastillo@ugr.es) on March 2025. 

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"


//Data product includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"

//ROOT includes
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <cmath>
#include "TTree.h"


#include <cmath> // std::hypot(), std::abs(), std::sqrt()
#include <iostream>
#include <memory>
#include <string>
#include <utility> // std::move()
#include <vector>

/**
 * @brief Matches optical flashes and charge slices based on their location.
 * 
 * 
 * This algorithm associates slices of charge from the TPC (`recob::Slice`) to
 * reconstructed scintillation light flashes (`recob::OpFlash`).
 *
 * 
 * Each slice is associated independently to the closest among the flashes.
 * Distance is computed between the charge centroid of the slice and the one 
 * of the flash (either in 2D or 3D).
 *  * 
 * Matching of a slice immediately fails if no charge is associated to it. 
 * Whenever the candidate pool contains at least a flash, the matching is 
 * "successful", meaning that it yields a result.
 * It is the privilege of the analyzer to decide whether that
 * is a good match or not, based on the information in the matching data
 * product associated to the slice. If the matching did not succeed, matching
 * information will still be associated with the slice, but most or all of its
 * information will be placeholders.
 * 
 * It is possible that a flash ends being associated with multiple slices.
 * 
 * In addition, if a trigger time is available, the algorithm will attempt to
 * identify the triggering flash and, if found, will provide the distance of
 * each slice from that flash.
 * 
 * Content of the matching information:
 * 
 * * `ChargeT0` (microseconds): time associated to the slice via particle flow
 *     object; it is left in the original time reference. If no `anab::T0`
 *     object is associated, this variable is assigned the magic value `-9999`.
 * 
 * ### Definition of the centroids
 * 
 * The centroid of the slice is defined as the position of all the valid space
 * points associated with the slice, weighted by their charge.
 * Space points (`recob::SpacePoint`) are associated each to a TPC hit
 * (`recob::Hit`); the charge of the hit is used as a weight, with no further
 * calibration applied.
 * If `CollectionOnly` is specified, only space points associated to hits on
 * a collection plane are included.
 * 
 * The centroid of a flash is defined only in its projection on the PMT plane.
 * Its definition is delegated to the flash reconstruction algorithm
 * (`recob::OpFlash::XCenter()`, `recob::OpFlash::YCenter()` and `recob::OpFlash::ZCenter()`).
 * In particular note that there is no attempt to enforce a location based on
 * the earliest light. In fact, the standard flash reconstruction algorithm
 * builds the center of the flash from _all_ the associated hits no matter
 * their time.
 * 
 * ### Determination of the trigger flash
 * 
 * If a trigger time is found from the `TriggerLabel` data product, the module
 * will attempt to identify the triggering flash.
 * The identification is based solely on the time of the trigger and of the
 * flash.
 * In ICARUS data, the time reference is the trigger time, therefore the time
 * of the global trigger is a fixed value by definition. Nevertheless, to
 * accommodate relative delays between the trigger and the light systems, an
 * expected time difference between the reconstrcuted flash (`TriggerDelay`)
 * is accounted for with a cushion ('TriggerTolerance'). These parameters are
 * measured in microseconds and controlled directly by the module caller.
 * In ICARUS simulation the mechanism is the same; the times are stored with
 * respect to the simulation time, but the module will take care of
 * the appropriate conversion, so the configuration parameter holds the same
 * meaning and scale as in data.
 * 
 * 
 * Input
 * ------
 * 
 * * `std::vector<recob::OpFlash>` (based on `OpFlashLabel`): collection of
 *     reconstructed flashes to be associated to the charge; also their
 *     associations to `recob::OpHit` (and those hits themselves).
 * * `std::vector<recob::Slice>` (based on `PandoraLabel`): collection of
 *     reconstructed TPC charge slices to be associated to the light; also their
 *     associations to `recob::Hit` and `recob::PFParticle`, and
 *     `recob::SpacePoint` objects associated to those hits
 *     (and also the associated objects as well).
 * 
 * 
 * Output
 * -------
 * 
 * A single collection, merging all the input, is produced for each of the
 * following data products in the _art_/ROOT output:
 * 
 * * `std::vector<sbn::TPCPMTBarycenterMatch>`: collection of matching information;
 *     one matching information object is present for each slice in the input
 *     collections, in the same order as the input.
 * * `art::Assns<sbn::TPCPMTBarycenterMatch, recob::Slice>`: association of each
 *     matched slice with its matching information.
 * * `art::Assns<sbn::TPCPMTBarycenterMatch, recob::OpFlash>`: association of each
 *     matched light flash with its matching information.
 * 
 * Note that while there is currently no direct association between the slice
 * and the flash, the information contained in `sbn::TPCPMTBarycenterMatch` is very
 * detailed.
 * 
 * In addition, if `FillMatchTree` is set, a ROOT tree called `matchTree` will
 * be written via `TFileService`. The tree contains one entry per slice,
 * whether matched or not,
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * * `OpFlashLabel` (string, mandatory): base of the tag of the input flashes.
 * * `PandoraLabel` (string, mandatory): base of the tag of input slices, and
 *     their associations to hits, space points, particle flow objects etc.
 * * `CollectionOnly` (flag, default: `true`): if set, only hits from
 *     collection planes will contribute to the centroid of the TPC slice.
 * * `Verbose` (flag, default: `false`): enables verbose output directly to
 *     console standard output.
 * * `FillMatchTree` (flag, default: `false`): if set to `true`, a ROOT tree
 *     with detailed matching information called `"matchTree"` will be written
 *     via `TFileService`.
 * * `Do3DMatching` (bool, default: `true`): if set to true the matching is peformed
 *     using the 3 dimensions (XYZ). If false, the matching is performed using only
 *     the reconstruction on the detection plane (YZ).
 */
class TPCPMTBarycenterMatchProducer : public art::EDProducer {
public:
  explicit TPCPMTBarycenterMatchProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCPMTBarycenterMatchProducer(TPCPMTBarycenterMatchProducer const&) = delete;
  TPCPMTBarycenterMatchProducer(TPCPMTBarycenterMatchProducer&&) = delete;
  TPCPMTBarycenterMatchProducer& operator=(TPCPMTBarycenterMatchProducer const&) = delete;
  TPCPMTBarycenterMatchProducer& operator=(TPCPMTBarycenterMatchProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  void InitializeSlice();                                                                     ///< Re-initialize all slice-level data members
  void updateChargeVars(double sumCharge, TVector3 const& sumPos, TVector3 const& sumPosSqr, std::array<double, 3> const& triggerFlashCenter); ///< Update slice-level data members with charge and trigger match info
  void updateFlashVars(art::Ptr<recob::OpFlash> flash);                      ///< Update slice-level data members with best match info
  void updateMatchInfo(sbn::TPCPMTBarycenterMatch& matchInfo);                                      ///< Update match product with slice-level data members
  void GetPCA(std::vector<double> const& x, std::vector<double> const& y, std::vector<double> const& weight, std::vector<double>& PCA ); ///< Get the PCA
  double GetSliceCharge(const std::vector<art::Ptr<recob::Hit>> &tpcHitsVec, const detinfo::DetectorPropertiesData det_prop, int tpc);
  double GetFlashLight(double flash_pe, std::vector<double>& total_dir_visibility, std::vector<double>& total_ref_visibility);
  void CreateOpHitList( std::vector<art::Ptr<recob::OpHit>> ophitlist, std::vector<double>& ophit_z, std::vector<double>& ophit_y, std::vector<double>& ophit_weight);

  // Input parameters
  std::vector<std::string>  fInputTags;            ///< Suffix added onto fOpFlashLabel and fPandoraLabel, used by ICARUS for separate cryostat labels but could be empty
  std::vector<std::string>  fOpFlashesModuleLabel;         ///< Label for PMT reconstruction products
  std::string               fPandoraLabel;         ///< Label for Pandora output products
  std::string               fOpT0Label;
  bool                      fCollectionOnly;       ///< Only use TPC spacepoints from the collection plane
  double                    fDistanceCandidateFlashes; ///< Maximum distance between candidate flashes to be considered for matching (cm)
  std::vector<double>       fCalAreaConst;         /// Calibration area constants for wire plane
  std::vector<int>          fSkipChannelList;
  double                    fOpDetVUVEff;           // Efficiencies for PMT detection
  double                    fOpDetVISEff;           // Efficiencies for PMT detection
  bool                      fVerbose;              ///< Print extra info
  bool                      fFillMatchTree;        ///< Fill an output TTree in the supplemental file
  bool                      fDo3DMatching;         ///< Wether to perform the matching in 3D or 2D
  std::vector<double>       fLightChargeRatioBounds; ///< Vector to store the distance between the barycenter of the charge and the barycenter of the light for each slice
  double                    fXError;
  double                    fYError;
  double                    fZError;
  double                    fAngleError;
  std::vector<double>       fFlashVetoWindow;
  // Event-level data members
  int                       fRun;                  ///< Number of the run being processed
  int                       fEvent;                ///< Number of the event being processed
  int                       fTPC;                 ///< Cryostat this event occured in
  int                       fSliceNum;             ///< Number of slice in the event
  // Slice-level data members 
  double                    fChargeT0;             ///< Start time for cathode-crossing PFPs, not always available (us)
  double                    fChargeTotal;          ///< Total charge in slice (integrated ADC counts)
  double                    fChargeCenterX;  ///< Weighted mean X position of spacepoints (cm)
  double                    fChargeCenterY;        ///< Weighted mean Y position of spacepoints (cm)
  double                    fChargeCenterZ;        ///< Weighted mean Z position of spacepoints (cm)
  double                    fChargeWidthX;         ///< Weighted standard deviation of X position of spacepoints (cm)
  double                    fChargeWidthY;         ///< Weighted standard deviation of Y position of spacepoints (cm)
  double                    fChargeWidthZ;         ///< Weighted standard deviation of Z position of spacepoints (cm)
  double                    fFlashTime;            ///< Matched OpFlash time (us)
  double                    fFlashTimeDebug;         ///< Matched OpFlash time error (us)
  double                    fFlashPEs;             ///< Brightness of matched flash (photoelectrons)
  double                    fFlashCenterX;         ///< Flash position X obtained through eta_pmt curve (cm)
  double                    fFlashCenterY;         ///< Weighted mean Y postion of hit PMTs (cm)
  double                    fFlashCenterZ;         ///< Weighted mean Z postion of hit PMTs (cm)
  double                    fFlashWidthY;          ///< Weighted standard deviation of Y postion of hit PMTs (cm)
  double                    fFlashWidthZ;          ///< Weighted standard deviation of Z postion of hit PMTs (cm)
  double                    fDeltaT;               ///< | Matched flash time - charge T0 | when available (us)
  double                    fDeltaY;               ///< | Matched flash Y center - charge Y center | (cm)
  double                    fDeltaZ;               ///< | Matched flash Z center - charge Z center | (cm)
  double                    fRadius;               ///< Hypotenuse of DeltaY and DeltaZ *parameter minimized by matching* (cm)
  double                    fChi2;                 ///< Chi2 to be minimized when matching flash to slice (dimensionless)
  double                    fScore;                ///< Score to be maximized when matching flash to slice (dimensionless)
  double                    fAngle;                ///< Angle between charge PCA and light PCA (degrees)
  double                    fDeltaY_Trigger;       ///< | Triggering flash Y center - charge Y center | (cm)
  double                    fDeltaZ_Trigger;       ///< | Triggering flash Z center - charge Z center | (cm)
  double                    fRadius_Trigger;       ///< Hypotenuse of DeltaY_Trigger and DeltaZ_Trigger (cm)
  TTree*                    fMatchTree;            ///< Tree to store all match information
  // Geometry service
  geo::WireReadoutGeom const& fWireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  opdet::sbndPDMapAlg fPDSMap;

  //Vector for PMT position
  std::vector<double> fOpDetID;
  std::vector<int> fOpDetType;
  std::vector<double> fOpDetX;
  std::vector<double> fOpDetY;
  std::vector<double> fOpDetZ;

  // Semi-analytical model for VUV and VIS light
  std::unique_ptr<phot::SemiAnalyticalModel> _semi_model;
  fhicl::ParameterSet _vuv_params;
  fhicl::ParameterSet _vis_params;
  std::shared_ptr<phot::OpticalPath> _optical_path_tool;

};


TPCPMTBarycenterMatchProducer::TPCPMTBarycenterMatchProducer(fhicl::ParameterSet const& p)
  : EDProducer{p},
  // More initializers here.
  fOpFlashesModuleLabel(p.get<std::vector<std::string>>("OpFlashesModuleLabel")),
  fPandoraLabel(p.get<std::string>("PandoraLabel")),
  fCollectionOnly(p.get<bool>("CollectionOnly", true)),
  fDistanceCandidateFlashes(p.get<double>("DistanceCandidateFlashes")), // cm
  fCalAreaConst(p.get<std::vector<double>>("CalAreaConst")),
  fOpDetVUVEff (p.get<double>("OpDetVUVEff")),
  fOpDetVISEff (p.get<double>("OpDetVISEff")),
  fVerbose(p.get<bool>("Verbose")),
  fFillMatchTree(p.get<bool>("FillMatchTree")),
  fDo3DMatching(p.get<bool>("Do3DMatching")),
  fLightChargeRatioBounds(p.get<std::vector<double>>("LightChargeRatioBounds")),
  fXError(p.get<double>("XError")), // cm
  fYError(p.get<double>("YError")), // cm
  fZError(p.get<double>("ZError")), // cm
  fAngleError(p.get<double>("AngleError")), // deg
  fFlashVetoWindow(p.get<std::vector<double>>("FlashVetoWindow")) // us
  {
  // Call appropriate produces<>() functions here.

  produces< std::vector<sbn::TPCPMTBarycenterMatch> >();
  produces< art::Assns<sbn::TPCPMTBarycenterMatch, recob::Slice> >();
  produces< art::Assns<sbn::TPCPMTBarycenterMatch, recob::OpFlash> >();

  if(fOpFlashesModuleLabel.size()!=2){
    throw cet::exception("TPCPMTBarycenterMatching")
    << "Module has been missconfigured. Number of OpFlashesModuleLabel must be 2!";
  }

  art::ServiceHandle<art::TFileService> tfs;
  if ( fFillMatchTree ) {
      art::ServiceHandle<art::TFileService> tfs;
      fMatchTree = tfs->make<TTree>("matchTree","TPC Slice - OpFlash Matching Analysis");

    //Event Info
    fMatchTree->Branch("run",                 &fRun,                 "run/I"                );
    fMatchTree->Branch("event",               &fEvent,               "event/I"              );
    fMatchTree->Branch("cryo",                &fTPC,                "cryo/I"                );
    fMatchTree->Branch("sliceNum",            &fSliceNum,            "sliceNum/I"           );

    //Charge Info
    fMatchTree->Branch("chargeT0",            &fChargeT0,            "chargeT0/d"           );
    fMatchTree->Branch("chargeTotal",         &fChargeTotal,         "chargeTotal/d"        );
    fMatchTree->Branch("chargeCenterXGlobal", &fChargeCenterX,       "chargeCenterX/d");
    fMatchTree->Branch("chargeCenterY",       &fChargeCenterY,       "chargeCenterY/d"      );
    fMatchTree->Branch("chargeCenterZ",       &fChargeCenterZ,       "chargeCenterZ/d"      );
    fMatchTree->Branch("chargeWidthX",        &fChargeWidthX,        "chargeWidthX/d"       );
    fMatchTree->Branch("chargeWidthY",        &fChargeWidthY,        "chargeWidthY/d"       );
    fMatchTree->Branch("chargeWidthZ",        &fChargeWidthZ,        "chargeWidthZ/d"       );

    //Matched Flash Info
    fMatchTree->Branch("flashTime",           &fFlashTime,           "flashTime/d"          );
    fMatchTree->Branch("flashPEs",            &fFlashPEs,            "flashPEs/d"           );
    fMatchTree->Branch("flashCenterX",        &fFlashCenterX,        "flashCenterX/d"       );
    fMatchTree->Branch("flashCenterY",        &fFlashCenterY,        "flashCenterY/d"       );
    fMatchTree->Branch("flashCenterZ",        &fFlashCenterZ,        "flashCenterZ/d"       );
    fMatchTree->Branch("flashWidthY",         &fFlashWidthY,         "flashWidthY/d"        );
    fMatchTree->Branch("flashWidthZ",         &fFlashWidthZ,         "flashWidthZ/d"        );

    //Match Quality Info
    fMatchTree->Branch("deltaT",              &fDeltaT,              "deltaT/d"             );
    fMatchTree->Branch("deltaY",              &fDeltaY,              "deltaY/d"             );
    fMatchTree->Branch("deltaZ",              &fDeltaZ,              "deltaZ/d"             );
    fMatchTree->Branch("radius",              &fRadius,              "radius/d"             );
    fMatchTree->Branch("angle",               &fAngle,               "angle/d"              );
    fMatchTree->Branch("chi2",                &fChi2,                "chi2/d"               );
    fMatchTree->Branch("score",               &fScore,               "score/d"              );
    fMatchTree->Branch("deltaZ_Trigger",      &fDeltaZ_Trigger,      "deltaZ_Trigger/d"     );
    fMatchTree->Branch("deltaY_Trigger",      &fDeltaY_Trigger,      "deltaY_Trigger/d"     );
    fMatchTree->Branch("radius_Trigger",      &fRadius_Trigger,      "radius_Trigger/d"     );

  } //End MatchTree

  //Fill the OpDet positions
  for(unsigned int opch=0; opch<fWireReadout.NOpChannels(); opch++){
    auto pdCenter = fWireReadout.OpDetGeoFromOpChannel(opch).GetCenter();
    fOpDetID.push_back(opch);
    fOpDetX.push_back(pdCenter.X());
    fOpDetY.push_back(pdCenter.Y());
    fOpDetZ.push_back(pdCenter.Z());
    if(fPDSMap.pdType(opch)=="pmt_coated") fOpDetType.push_back(0);
    else if(fPDSMap.pdType(opch)=="pmt_uncoated") fOpDetType.push_back(1);
    else if(fPDSMap.pdType(opch)=="xarapuca_vuv") fOpDetType.push_back(2);
    else if(fPDSMap.pdType(opch)=="xarapuca_vis") fOpDetType.push_back(3);
    else fOpDetType.push_back(-1);
  }

  _vuv_params = p.get<fhicl::ParameterSet>("VUVHits");
  _vis_params = p.get<fhicl::ParameterSet>("VIVHits");
  _optical_path_tool = std::shared_ptr<phot::OpticalPath>(art::make_tool<phot::OpticalPath>(p.get<fhicl::ParameterSet>("OpticalPathTool")));
  _semi_model = std::make_unique<phot::SemiAnalyticalModel>(_vuv_params, _vis_params, _optical_path_tool, true, false);
}

void TPCPMTBarycenterMatchProducer::produce(art::Event& e)
{
  // Implementation of required member function here.
  fEvent  = e.id().event();
  fRun    = e.run();

  // Detector properties and clocks
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);


  //Initialize new data products
  auto matchInfoVector = std::make_unique< std::vector<sbn::TPCPMTBarycenterMatch> >();
  art::PtrMaker< sbn::TPCPMTBarycenterMatch > const makeInfoPtr(e); 
  auto sliceAssns = std::make_unique< art::Assns<sbn::TPCPMTBarycenterMatch, recob::Slice> >();
  auto flashAssns = std::make_unique< art::Assns<sbn::TPCPMTBarycenterMatch, recob::OpFlash> >();

  /* ~~~~~~~~~~~~~~~~~~~~ Flash Section
  *
  * Here we gather the OpFlashes found in this cryostat and their OpHits
  * We iterate through the flashes to identify a triggering flash
  */
    std::array<double, 3> triggerFlashCenter;
    for ( size_t tpc=0; tpc<2 ; tpc ++) {
      //Fetch the flashes and their associated hits, pointer vector needed for generating associations
      art::Handle flashHandle = e.getHandle<std::vector<recob::OpFlash>>(fOpFlashesModuleLabel[tpc]);
      int nFlashes = (*flashHandle).size();
      triggerFlashCenter = {-9999., -9999., -9999.};
      double minTime = 99999.;
      for (int i = 0; i < nFlashes; ++i) {
          const recob::OpFlash &flash = (*flashHandle)[i];
          double _flashtime = flash.AbsTime();
          if ( std::abs(_flashtime)< minTime) {
              minTime = std::abs(_flashtime);
              triggerFlashCenter = {flash.XCenter(), flash.YCenter(), flash.ZCenter()};
          }
      }
    }


  // TODO (acastill): evaluate if there are two flashes in time coincidence and treat them as one.

  /* ~~~~~~~~~~~~~~~~~~~~ TPC Section
  * Here we start by gathering the Slices in the event
  * For each slice, the charge centroid is first calculated
  * Then we iterate through flashes to identify the best match flash
  * If a triggering flash was found in this cyrostat, the barycenter distance to the triggering flash is also stored
  */

  //Fetch slices, TPC hits, and PFPs; pointer vector needed for generating associations
  ::art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fPandoraLabel, sliceHandle);

  // Slice to PFP assns
  unsigned nSlices = (*sliceHandle).size();
  ::art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPandoraLabel, pfpHandle);

  // PFP to metadata assns
  art::FindManyP<recob::PFParticle> slice_pfp_assns(sliceHandle, e, fPandoraLabel);

  // Hit to PFP assns
  art::FindManyP<recob::Cluster> pfp_cluster_assns(pfpHandle, e, fPandoraLabel);
  ::art::Handle<std::vector<recob::Cluster>> clusterHandle;
  e.getByLabel(fPandoraLabel, clusterHandle);

  // Cluster to Hit assns
  art::FindManyP<recob::Hit> cluster_hit_assns (clusterHandle, e, fPandoraLabel);

    

  //For slice...
  for ( unsigned j = 0; j < nSlices; j++ ) {
    //For each TPC
    std::vector<sbn::TPCPMTBarycenterMatch> sliceMatchInfoVector;
    std::vector<art::Ptr<sbn::TPCPMTBarycenterMatch>>  infoPtrVector;
    std::vector<art::Ptr<recob::OpFlash>> flashPtrVector;
    fSliceNum = j;
    const art::Ptr<recob::Slice> slicePtr { sliceHandle, j };
    //Vector for recob PFParticles
    std::vector<art::Ptr<recob::Hit>> tpcHitsVec;
    std::vector<art::Ptr<recob::PFParticle>> pfpVect = slice_pfp_assns.at(j);

    // Get the hits associated to the PFParticles in this slice
    for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
      std::vector<art::Ptr<recob::Cluster>> cluster_v = pfp_cluster_assns.at(pfp.key());
      for(size_t i=0; i<cluster_v.size(); i++){
        std::vector<art::Ptr<recob::Hit>> hitVect = cluster_hit_assns.at(cluster_v[i].key());
        tpcHitsVec.insert(tpcHitsVec.end(), hitVect.begin(), hitVect.end());
      }
    }

    int nPFPs = pfpVect.size();
    //Retrieve Pandora's T0 for this slice if available, same for every PFP in slice so we only need one
    if ( nPFPs != 0 ) {
      art::FindOne<anab::T0> f1T0( {pfpVect.at(0)}, e, fPandoraLabel);
      if ( f1T0.at(0).isValid() ) {
        fChargeT0 = f1T0.at(0).ref().Time() / 1e3;
      }
    }

    for ( size_t tpc=0; tpc<2 ; tpc ++) {
      fTPC = tpc;
      InitializeSlice();
      sbn::TPCPMTBarycenterMatch sliceMatchInfo;
      updateMatchInfo(sliceMatchInfo);
      art::FindOne<recob::SpacePoint> f1SpacePoint(tpcHitsVec, e, fPandoraLabel);

      std::vector<double> hit_z;
      std::vector<double> hit_y;
      std::vector<double> hit_weight;

      int nHits = tpcHitsVec.size();

      double thisCharge;
      double sumCharge = 0.;
      TVector3 sumPos {0.,0.,0.};
      TVector3 sumPosSqr {0.,0.,0.};

      size_t maxChargePlaneIdx=99999;

      // If we are not using collection plane only, we need to find the plane with the most charge
      if(!fCollectionOnly){
        std::vector<double> PlaneCharge(3, 0.); // Vector to store the charge for each plane
        //Loop to get the charge of each plane
        for ( int k = 0; k < nHits; k++ ) {
          if ( !f1SpacePoint.at(k).isValid() ) continue;
          // If hit does not belong to correct tpc then skip
          const art::Ptr<recob::Hit> &tpcHit = tpcHitsVec.at(k);
          const recob::SpacePoint point = f1SpacePoint.at(k).ref();
          TVector3 const thisPoint = point.XYZ();
          if ((tpc == 0) == (thisPoint.X() > 0)) continue; // Skip if the point is not in the TPC we are considering
          int plane = tpcHit->WireID().Plane;
          PlaneCharge[plane] += tpcHit->Integral();
        }
        // Idx to select the plane with the most charge
        maxChargePlaneIdx = std::distance(PlaneCharge.begin(), std::max_element(PlaneCharge.begin(), PlaneCharge.end()));
      }

      //For hit...
      for ( int k = 0; k < nHits; k++ ) {
        // If hit does not belong to correct tpc then skip 
        const art::Ptr<recob::Hit> &tpcHit = tpcHitsVec.at(k);

        //Only use hits with associated SpacePoints, and optionally only collection plane hits
        if ( !f1SpacePoint.at(k).isValid() ) continue;
        if ( fCollectionOnly && tpcHit->WireID().Plane != 2 ) continue;
        else if ( !fCollectionOnly && tpcHit->WireID().Plane != maxChargePlaneIdx) continue; // If not using collection plane only, skip hits not in the plane with the most charge
        const recob::SpacePoint point = f1SpacePoint.at(k).ref();
        TVector3 const thisPoint = point.XYZ();
        if ((tpc == 0) == (thisPoint.X() > 0)) continue; // Skip if the point is not in the TPC we are considering
        thisCharge = tpcHit->Integral();
        TVector3 const thisPointSqr {thisPoint.X()*thisPoint.X(), thisPoint.Y()*thisPoint.Y(), thisPoint.Z()*thisPoint.Z()};
        sumCharge += thisCharge;
        sumPos += thisPoint * thisCharge;
        sumPosSqr += thisPointSqr * thisCharge;
        // Store hit information for PCA
        hit_y.push_back(thisPoint.Y());
        hit_z.push_back(thisPoint.Z());
        hit_weight.push_back(thisCharge);
      } //End for hit

      //No charge found in slice...
      if ( sumCharge == 0. ) {
        if ( fFillMatchTree ) fMatchTree->Fill();
        if ( fVerbose ) std::cout << "No charge found in Event: " << fEvent << " Slice: " << j << "! Continuing..."  << std::endl;
        continue;
      }

      //Update charge variables
      updateChargeVars(sumCharge, sumPos, sumPosSqr, triggerFlashCenter);
      updateMatchInfo(sliceMatchInfo);
      
      std::vector<double> ChargePCA = {0.,0.};
      GetPCA(hit_z, hit_y, hit_weight, ChargePCA);

      //Get the charge of the slice 
      double sliceCharge = this->GetSliceCharge(tpcHitsVec, det_prop, tpc);
      //Get the visibility map of the charge barycenter
      geo::Point_t SliceXYZ = {fChargeCenterX, fChargeCenterY, fChargeCenterZ};
      std::vector<double> direct_visibility;
      std::vector<double> reflect_visibility;
      _semi_model->detectedDirectVisibilities(direct_visibility, SliceXYZ);
      _semi_model->detectedReflectedVisibilities(reflect_visibility, SliceXYZ);

      int matchIndex = -1;
      double minDistance = 1e6;
      double thisFlashCenterX, thisFlashCenterY, thisFlashCenterZ, thisDistance;

      // --- Read Recob OpFlash
      ::art::Handle<std::vector<recob::OpFlash>> flashHandle;
      e.getByLabel(fOpFlashesModuleLabel[tpc], flashHandle);
      std::vector< art::Ptr<recob::OpFlash> > flashVect;
      art::fill_ptr_vector(flashVect, flashHandle);

      //OpHit OpFlash assns
      art::FindManyP<recob::OpHit> flash_ophit_assns(flashHandle, e, fOpFlashesModuleLabel[tpc]);

      //Vector to store the idxs of the candidate flashes
      std::vector<int> candidateFlashIdxs;
      int nFlashes = flashVect.size();

      //For flash...
      for ( int m = 0; m < nFlashes; m++ ) {
        auto & flash = flashVect[m];

        double flashTime = flash->Time();

        if(flashTime<fFlashVetoWindow[0] || flashTime>fFlashVetoWindow[1]) continue;

        //Find index of flash that minimizes barycenter distance in XYZ place
        thisFlashCenterX = flash->XCenter();
        thisFlashCenterY = flash->YCenter();
        thisFlashCenterZ = flash->ZCenter();

        if(fDo3DMatching)
          thisDistance = std::hypot( (thisFlashCenterX - fChargeCenterX), (thisFlashCenterY - fChargeCenterY), (thisFlashCenterZ - fChargeCenterZ) );
        else thisDistance = std::hypot( (thisFlashCenterY - fChargeCenterY), (thisFlashCenterZ - fChargeCenterZ) );

        if ( thisDistance < minDistance ) {
          minDistance = thisDistance;
          matchIndex = m;
        }
        // Fill the vector with the idxs of the flashes
        if(thisDistance<fDistanceCandidateFlashes){
          candidateFlashIdxs.push_back(m);
        }
      } //End for flash

      //No valid match found...
      if ( matchIndex == -1 ) {
        if ( fFillMatchTree ) fMatchTree->Fill();
        if ( fVerbose ) std::cout << "No matching flash found for Event: " << fEvent << " Slice: " << j << "! Continuing..."  << std::endl;
        continue;
      }

      // -----------------> Here we would end with the first pool of candidate flashes <-----------------

      double minChi2 = 1e6;
      for(size_t i=0; i<candidateFlashIdxs.size(); i++)
      {
        fFlashTimeDebug=-99999;
        int idx = candidateFlashIdxs[i];
        auto & flash = flashVect[idx];
        double _currentFlashX = flash->XCenter();
        double _currentFlashY = flash->YCenter();
        double _currentFlashZ = flash->ZCenter();
        //Get the ophits associated to this opflash
        std::vector<art::Ptr<recob::OpHit>> ophitlist = flash_ophit_assns.at(flash.key());
        std::vector<double> ophit_z;
        std::vector<double> ophit_y;
        std::vector<double> ophit_weight;
        CreateOpHitList(ophitlist, ophit_z, ophit_y, ophit_weight);

        double FlashLight = this->GetFlashLight(flash->TotalPE(), direct_visibility, reflect_visibility);
        fFlashTimeDebug = flash->Time();
        double lightChargeRatio = (FlashLight/sliceCharge);
        std::vector<double> LightPCA = {0.,0.};
        GetPCA(ophit_z, ophit_y, ophit_weight, LightPCA);
        double LightChargeCosine = abs( ChargePCA[0]*LightPCA[0] + ChargePCA[1]*LightPCA[1] ) / ( std::hypot(ChargePCA[0], ChargePCA[1]) * std::hypot(LightPCA[0], LightPCA[1]) );
        double angle = std::acos(LightChargeCosine)*(180/M_PI);
        double chi2;
        if(fDo3DMatching) chi2 = std::pow(fChargeCenterX-_currentFlashX, 2)/std::pow(fXError, 2) + std::pow(fChargeCenterY-_currentFlashY, 2)/std::pow(fYError, 2) + std::pow(fChargeCenterZ-_currentFlashZ, 2)/std::pow(fZError, 2) + std::pow(angle, 2)/std::pow(fAngleError, 2);
        else
          chi2 = std::pow(fChargeCenterY-_currentFlashY, 2)/std::pow(fYError, 2) + std::pow(fChargeCenterZ-_currentFlashZ, 2)/std::pow(fZError, 2) + std::pow(angle, 2)/std::pow(fAngleError, 2);

        if(lightChargeRatio < fLightChargeRatioBounds[0] || lightChargeRatio > fLightChargeRatioBounds[1]) continue;
        if(chi2 < minChi2)
        {
          minChi2 = chi2;
          fAngle = angle;
          matchIndex = idx;
        }
      }

      //Best match flash pointer
      unsigned unsignedMatchIndex = matchIndex;
      const art::Ptr<recob::OpFlash> flashPtr { flashHandle, unsignedMatchIndex };
      //Update match info
      updateFlashVars(flashPtr);
      fChi2 = minChi2;
      fScore = 1./fChi2;
      updateMatchInfo(sliceMatchInfo);
      sliceMatchInfoVector.push_back(sliceMatchInfo);
      art::Ptr<sbn::TPCPMTBarycenterMatch> const infoPtr = makeInfoPtr(matchInfoVector->size());
      infoPtrVector.push_back(infoPtr);
      flashPtrVector.push_back(flashPtr);
      if ( fFillMatchTree ) fMatchTree->Fill();
    } //End for tpc

    int maxFlashIdx=-1;
    int minNPes = -1000000000;

    for(int i=0; i<static_cast<int>(flashPtrVector.size()); i++)
    {
      int nPEs = flashPtrVector[i]->TotalPE();
      if(nPEs>minNPes)
      {
        maxFlashIdx=i;
        minNPes = nPEs;
      }
    }

    if(maxFlashIdx!=-1) 
    {
      sliceAssns->addSingle(infoPtrVector[maxFlashIdx], slicePtr);
      flashAssns->addSingle(infoPtrVector[maxFlashIdx], flashPtrVector[maxFlashIdx]);
      matchInfoVector->push_back(std::move(sliceMatchInfoVector[maxFlashIdx]));
    }
  } //End for slice

  //Store new products at the end of the event
  e.put(std::move(matchInfoVector));
  e.put(std::move(sliceAssns));
  e.put(std::move(flashAssns));

} //End produce()

void TPCPMTBarycenterMatchProducer::InitializeSlice() {
  fChargeT0 = -9999.;
  fChargeTotal = -9999.;
  fChargeCenterX = -9999.;
  fChargeCenterY = -9999.;
  fChargeCenterZ = -9999.;
  fChargeWidthX = -9999.;
  fChargeWidthY = -9999.;
  fChargeWidthZ = -9999.;
  fFlashTime = -9999.;
  fFlashPEs = -9999.;
  fFlashCenterX = -9999.;
  fFlashCenterY = -9999.;
  fFlashCenterZ = -9999.;
  fFlashWidthY = -9999.;
  fFlashWidthZ = -9999.;
  fDeltaT = -9999.;
  fDeltaY = -9999.;
  fDeltaZ = -9999.;
  fRadius = -9999.;
  fChi2 = -9999.;
  fScore = -99999.;
  fAngle = -9999.;
  fDeltaZ_Trigger = -9999.;
  fDeltaY_Trigger = -9999.;
  fRadius_Trigger = -9999.;
} //End InitializeSlice()


void TPCPMTBarycenterMatchProducer::updateChargeVars(double sumCharge, TVector3 const& sumPos, TVector3 const& sumPosSqr, std::array<double, 3> const& triggerFlashCenter) {
  fChargeCenterX = sumPos[0] / sumCharge;
  fChargeCenterY = sumPos[1] / sumCharge;
  fChargeCenterZ = sumPos[2] / sumCharge;
  fChargeWidthX = std::sqrt( sumPosSqr[0]/sumCharge - (sumPos[0]/sumCharge)*(sumPos[0]/sumCharge) );
  fChargeWidthY = std::sqrt( sumPosSqr[1]/sumCharge - (sumPos[1]/sumCharge)*(sumPos[1]/sumCharge) );
  fChargeWidthZ = std::sqrt( sumPosSqr[2]/sumCharge - (sumPos[2]/sumCharge)*(sumPos[2]/sumCharge) );
  fChargeTotal = sumCharge;
  if ( triggerFlashCenter[1] != -9999. ) {
    fDeltaY_Trigger = abs(triggerFlashCenter[0] - fChargeCenterY);
    fDeltaZ_Trigger = abs(triggerFlashCenter[1] - fChargeCenterZ);
    fRadius_Trigger = std::hypot(fDeltaY_Trigger, fDeltaZ_Trigger);
  }
} //End updateChargeVars()


void TPCPMTBarycenterMatchProducer::updateFlashVars(art::Ptr<recob::OpFlash> flash) {
  double matchedTime = flash->Time();
  double matchedXCenter = flash->XCenter();
  double matchedYCenter = flash->YCenter();
  double matchedZCenter = flash->ZCenter();
  double matchedYWidth = flash->YWidth();
  double matchedZWidth = flash->ZWidth();

  fFlashTime = matchedTime;
  fFlashPEs =  flash->TotalPE();
  fFlashCenterX = matchedXCenter;
  fFlashCenterY = matchedYCenter;
  fFlashCenterZ = matchedZCenter;
  fFlashWidthY = matchedYWidth;
  fFlashWidthZ = matchedZWidth;
  if ( fChargeT0 != -9999 ) fDeltaT = abs(matchedTime - fChargeT0);
  fDeltaY = abs(matchedYCenter - fChargeCenterY);
  fDeltaZ = abs(matchedZCenter - fChargeCenterZ);
  fRadius = std::hypot(fDeltaY, fDeltaZ);
} //End updateFlashVars()


void TPCPMTBarycenterMatchProducer::updateMatchInfo(sbn::TPCPMTBarycenterMatch& matchInfo) {
  matchInfo.chargeTotal = fChargeTotal;
  matchInfo.chargeCenter = {fChargeCenterX, fChargeCenterY, fChargeCenterZ};
  matchInfo.chargeWidth = {fChargeWidthX, fChargeWidthY, fChargeWidthZ};
  matchInfo.flashTime = fFlashTime;
  matchInfo.flashPEs = fFlashPEs;
  matchInfo.flashCenter = {fFlashCenterX, fFlashCenterY, fFlashCenterZ};
  matchInfo.flashWidth = {-9999., fFlashWidthY, fFlashWidthZ};
  matchInfo.deltaT = fDeltaT;
  matchInfo.deltaY = fDeltaY;
  matchInfo.deltaZ = fDeltaZ;
  matchInfo.radius = fRadius;
  matchInfo.radius_Trigger = fRadius_Trigger;
  matchInfo.chi2 = fChi2;
  matchInfo.score = fScore;
  matchInfo.angle = fAngle;
} //End updateMatchInfo()

void TPCPMTBarycenterMatchProducer::GetPCA(std::vector<double> const&x, std::vector<double> const&y, std::vector<double> const&weight, std::vector<double> & PCA ) {

  double weight_sum = std::accumulate(weight.begin(), weight.end(), 0.0);
  //Create the matrix of points
  size_t n = x.size();
  Eigen::MatrixXd points(n, 2);
  for (size_t i = 0; i < n; ++i) {
      points(i, 0) = x[i];
      points(i, 1) = y[i];
  }

  //Get the weighted centroid
  Eigen::Vector2d weighted_sum(0.0, 0.0);
  for (size_t i = 0; i < n; ++i) {
      weighted_sum(0) += weight[i] * x[i];
      weighted_sum(1) += weight[i] * y[i];
  }

  //Center the points wrt centroid
  Eigen::Vector2d centroid = weighted_sum / weight_sum;
  Eigen::MatrixXd centered_points = points.rowwise() - centroid.transpose();

  //Get the weighted covariance matrix
  Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();
  for (size_t i = 0; i < n; ++i) {
      Eigen::Vector2d pt = centered_points.row(i);
      cov += std::abs(weight[i]) * (pt * pt.transpose());
  }
  cov /= weight_sum;

  // Get eigenvalues / eigenvectors 
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(cov);
  Eigen::VectorXd eigvals = solver.eigenvalues();
  Eigen::Matrix2d eigvecs = solver.eigenvectors();

  // Get Principal component (largest eigenvalue)
  int maxIndex;
  eigvals.maxCoeff(&maxIndex);
  Eigen::Vector2d principal_direction = eigvecs.col(maxIndex);

  PCA[0] = principal_direction(0);
  PCA[1] = principal_direction(1);

}



double TPCPMTBarycenterMatchProducer::GetSliceCharge(const std::vector<art::Ptr<recob::Hit>> &tpcHitsVec,  const detinfo::DetectorPropertiesData det_prop, int tpc) {

    std::vector<double> plane_charge{0.,0.,0.};
    std::vector<int>    plane_hits{0,0,0};
    for (size_t i=0; i < tpcHitsVec.size(); i++){
      auto hit = tpcHitsVec[i];
      if( hit->WireID().TPC != static_cast<unsigned>(tpc) ) continue; // Skip hits not in the TPC we are considering
      auto drift_time = (hit->PeakTime() - 500)*0.5; // assuming TPC beam readout starts at 500 ticks, conversion = 0.5 us/tick
      double atten_correction = std::exp(drift_time/det_prop.ElectronLifetime()); // exp(us/us)
      auto hit_plane = hit->View();
      plane_charge.at(hit_plane) += hit->Integral()*atten_correction*(1/fCalAreaConst.at(hit_plane));
      plane_hits.at(hit_plane)++; 
    }

    uint bestHits =  std::max_element(plane_hits.begin(), plane_hits.end()) - plane_hits.begin();
    double _comp_charge  = plane_charge.at(bestHits);

    return _comp_charge;
}

void TPCPMTBarycenterMatchProducer::CreateOpHitList( std::vector<art::Ptr<recob::OpHit>> ophitlist, std::vector<double>& ophit_z, std::vector<double>& ophit_y, std::vector<double>& ophit_weight) {

  std::map<std::pair<int, int>, double> footPrintMap;
  for (size_t i=0; i<ophitlist.size(); i++)
  {
    auto &ophit = ophitlist[i];
    int channel = ophit->OpChannel();
    int pos_z = static_cast<int>(fOpDetZ[channel]);
    int pos_y = static_cast<int>(fOpDetY[channel]);
    double weight = ophit->PE();
    footPrintMap[{pos_z, pos_y}] += weight;
  }

  for (const auto& opdet : footPrintMap) {
    int z = opdet.first.first;
    int y = opdet.first.second;
    int weight = opdet.second;
    ophit_z.push_back(z);
    ophit_y.push_back(y);
    ophit_weight.push_back(weight*weight);
  }

}

double TPCPMTBarycenterMatchProducer::GetFlashLight(double flash_pe, std::vector<double>& dir_visibility, std::vector<double>& ref_visibility) {

  double tot_visibility=0;

  for(size_t ch=0; ch<dir_visibility.size(); ch++){
    if (std::find(fSkipChannelList.begin(), fSkipChannelList.end(), ch) != fSkipChannelList.end()) continue;
    if(fOpDetType[ch]==0) tot_visibility += fOpDetVUVEff*dir_visibility[ch] + fOpDetVISEff*ref_visibility[ch];
    else if(fOpDetType[ch]==1) tot_visibility += fOpDetVISEff*ref_visibility[ch];
    else continue; // skip other types
  }
  //std::cout << " The number of PEs in the flash is " << flash_pe << " the total direct visibility is " << total_dir_visibility << " the total reflected visibility is " << total_ref_visibility << " with a VUV QE " << fOpDetVUVEff << " and a vis QE " << fOpDetVISEff << " so the total visibility is " << tot_visibility << std::endl;
  if((flash_pe == 0) || std::isinf(1/tot_visibility))
    return 0.0;
  // deposited light is inverse of visibility * PE count
  double total_gamma = (1/tot_visibility)*flash_pe;
  return total_gamma;
}


DEFINE_ART_MODULE(TPCPMTBarycenterMatchProducer)