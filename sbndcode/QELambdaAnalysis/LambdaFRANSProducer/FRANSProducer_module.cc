////////////////////////////////////////////////////////////////////////
// Class:       FRANSProducer
// Plugin Type: Producer
// File:        FRANSProducer_module.cc
////////////////////////////////////////////////////////////////////////

#include "FRANSProducer_module.hh"


namespace {
  template <typename T>
  bool fetchHandle(art::Event const& e, std::string const& label, art::Handle<T>& handle)
  {
    handle = e.getHandle<T>(label);
    if (!handle.isValid()) {
      mf::LogWarning("FRANSProducer")
        << "Product not found (label '" << label << "'), event "
        << e.id().run() << "_" << e.id().subRun() << "_" << e.id().event();
      return false;
    }
    return true;
  }
}


// -------- Constructor
qelambda::FRANSProducer::FRANSProducer(fhicl::ParameterSet const& p) 
  : EDProducer{p},
  fMCTruthLabel( p.get< std::string >("MCTruthLabel") ),
  fMCParticleLabel( p.get< std::string >("MCLabel") ),
  fHitLabel( p.get< std::string >("HitLabel") ),
  fVertexLabel( p.get< std::string >("VertexLabel") ),
  fRecoLabel( p.get< std::string >("RecoLabel") ),
  fTPCPMTBarycenterFMModuleLabel( p.get< std::string >("TPCPMTBarycenterFMModuleLabel") ),

  fTrackLabel(p.get<std::string>("TrackLabel")) ,
  fShowerLabel(p.get<std::string>("ShowerLabel")) ,
  fPFParticleMetadataLabel(p.get<std::string>("PFParticleMetadataLabel")) ,

  fView( p.get<int>("View") ),
  fBcFMScoreCut( p.get< float >("BcFMScoreCut") ),
  fTrkLenCut( p.get< float >("TrkLenCut") ),
  fShwrECut( p.get< float >("ShwrECut") ),
  fVtxDistCut( p.get< float >("VtxDistCut") ),
  fShwrNHits( p.get< float >("ShwrNHits") ),
  fnPrimPfpsCut1( p.get< float >("nPrimPfpsCut1") ),
  fnPrimPfpsCut2( p.get< float >("nPrimPfpsCut2") ),
  fnPrimTrksCut1( p.get< float >("nPrimTrksCut1") ),
  fnPrimTrksCut2( p.get< float >("nPrimTrksCut2") ),
  fnPrimShwrsCut1( p.get< float >("nPrimShwrsCut1") ),
  fnPrimShwrsCut2( p.get< float >("nPrimShwrsCut2") ),
  fnPrimTrksExitCut( p.get< float >("nPrimTrksExitCut") ),
  fminGapCut( p.get< float >("minGapCut") ),

  fFhiclPsetFRANS(  p.get<fhicl::ParameterSet>("ChargeDensityAlg") ),
  fFhiclPsetTrackFinder(  p.get<fhicl::ParameterSet>("TrackFinderAlg") ),
  fFhiclPsetVertexFinder(  p.get<fhicl::ParameterSet>("VertexFinderAlg") ),
  fFhiclPsetHough(  p.get<fhicl::ParameterSet>("HoughAlg") ),
  fFhiclPsetTPCLines(  p.get<fhicl::ParameterSet>("TPCLinesAlg") ),
  fFhiclPsetCalo(  p.get<fhicl::ParameterSet>("CaloAlg") ),
  fMichelPset(p.get<fhicl::ParameterSet>("MichelTagging"))
{

  fhicl::ParameterSet pt = fFhiclPsetFRANS;
  fPsetFRANS.ApplyRawSmoothing = pt.get<bool>("ApplyRawSmoothing");
  fPsetFRANS.ApplySmoothing = pt.get<bool>("ApplySmoothing");
  fPsetFRANS.ApplyCumulativeSmoothing = pt.get<bool>("ApplyCumulativeSmoothing");
  fPsetFRANS.NDriftPack = pt.get<unsigned int>("NDriftPack");
  fPsetFRANS.NWirePack = pt.get<unsigned int>("NWirePack");
  fPsetFRANS.ExpoAvSmoothPar = pt.get<double>("ExpoAvSmoothPar");
  fPsetFRANS.UnAvNeighbours = pt.get<int>("UnAvNeighbours");
  fPsetFRANS.CumulativeCut = pt.get<double>("CumulativeCut");
  fPsetFRANS.SlidingWindowN = pt.get<int>("SlidingWindowN");
  fPsetFRANS.NSamplesBeginSlope = pt.get<int>("NSamplesBeginSlope");
  fPsetFRANS.MaxRadius = pt.get<int>("MaxRadius");
  fPsetFRANS.UseHitWidth = pt.get<bool>("UseHitWidth");
  fPsetFRANS.Verbose = pt.get<int>("Verbose");
  fPsetFRANS.CalculateScore = pt.get<bool>("CalculateScore");
  fPsetFRANS.UseAlpha = pt.get<bool>("UseAlpha");
  fPsetFRANS.UseIota = pt.get<bool>("UseIota");
  fPsetFRANS.TMVAFilename = pt.get<std::string>("TMVAFilename");
  fPsetFRANS.OutputPath = pt.get<std::string>("OutputPath");

  // Find BDT weights files
  std::string file_name;
  cet::search_path sp("FW_SEARCH_PATH");
  if ( !sp.find_file(fPsetFRANS.TMVAFilename, file_name) )
    throw cet::exception("FRANSProducer") << "BDT file " <<
        fPsetFRANS.TMVAFilename << " not found in FW_SEARCH_PATH\n";
  // std::cout<<" BDT weights file name: "<<file_name<<std::endl;
  fPsetFRANS.TMVAFilename = file_name;

  pt = fFhiclPsetTrackFinder;
  fPsetTrackFinder.MaxDTube = pt.get<int>("MaxDTube");
  fPsetTrackFinder.MaxDCluster = pt.get<double>("MaxDCluster");
  fPsetTrackFinder.SingleWireMode = pt.get<bool>("SingleWireMode");
  fPsetTrackFinder.MinClusterHits = pt.get<int>("MinClusterHits");
  fPsetTrackFinder.DCleaning = pt.get<double>("DCleaning");
  fPsetTrackFinder.ClusterCompletenessCut = pt.get<double>("ClusterCompletenessCut");
  fPsetTrackFinder.ClusterAngleCut = pt.get<double>("ClusterAngleCut");
  fPsetTrackFinder.CaptureMissingHits = pt.get<bool>("CaptureMissingHits");
  fPsetTrackFinder.MinTrackHits = pt.get<int>("MinTrackHits");
  fPsetTrackFinder.HitDensityThreshold = pt.get<float>("HitDensityThreshold");
  fPsetTrackFinder.UseCompactness = pt.get<bool>("UseCompactness");
  fPsetTrackFinder.ConnectednessTol = pt.get<double>("ConnectednessTol");
  fPsetTrackFinder.ConnectednessWidthTol = pt.get<double>("ConnectednessWidthTol");
  fPsetTrackFinder.CompactnessTol = pt.get<double>("CompactnessTol");
  fPsetTrackFinder.Verbose = pt.get<int>("Verbose");

  pt = fFhiclPsetVertexFinder;
  fPsetVertexFinder.MaxDistToEdge = pt.get<double>("MaxDistToEdge");
  fPsetVertexFinder.RefineVertexIntersection = pt.get<bool>("RefineVertexIntersection");
  fPsetVertexFinder.SlidingWindowN = pt.get<int>("SlidingWindowN");
  fPsetVertexFinder.UseEdgesDiscard = pt.get<bool>("UseEdgesDiscard");
  fPsetVertexFinder.MaxTrackFractionInMain = pt.get<float>("MaxTrackFractionInMain");
  fPsetVertexFinder.DecideMainTrack = pt.get<bool>("DecideMainTrack");
  fPsetVertexFinder.AddCollinearLines = pt.get<bool>("AddCollinearLines");
  fPsetVertexFinder.VertexDistanceROI = pt.get<float>("VertexDistanceROI");
  fPsetVertexFinder.MinWires = pt.get<int>("MinWires");
  fPsetVertexFinder.AngleTolerance = pt.get<float>("AngleTolerance");
  fPsetVertexFinder.TriangleInequalityTol = pt.get<float>("TriangleInequalityTol");
  fPsetVertexFinder.VertexHitsTol = pt.get<float>("VertexHitsTol");
  fPsetVertexFinder.VertexHitsMinHits = pt.get<int>("VertexHitsMinHits");
  fPsetVertexFinder.VertexCompactnessTol = pt.get<float>("VertexCompactnessTol");
  fPsetVertexFinder.MinTrackOccupancy = pt.get<float>("MinTrackOccupancy");
  fPsetVertexFinder.MinTrackGoodness = pt.get<float>("MinTrackGoodness");
  fPsetVertexFinder.MakeCalorimetry = pt.get<bool>("MakeCalorimetry");
  fPsetVertexFinder.Verbose = pt.get<int>("Verbose");
  
  pt = fFhiclPsetHough;
  fPsetHough.MaxRadiusLineHypothesis = pt.get<double>("MaxRadiusLineHypothesis");
  fPsetHough.ThetaRes = pt.get<double>("ThetaRes");
  fPsetHough.MaxDistanceTube = pt.get<double>("MaxDistanceTube");
  fPsetHough.MinHoughHits = pt.get<int>("MinHoughHits");
  fPsetHough.Verbose = pt.get<int>("Verbose");

  pt = fFhiclPsetTPCLines;
  fPsetTPCLines.MaxRadius = pt.get<double>("MaxRadius");
  fPsetTPCLines.DriftConversion = pt.get<double>("DriftConversion");
  fPsetTPCLines.MaxMeanOccupation = pt.get<float>("MaxMeanOccupation");
  fPsetTPCLines.MaxHoughTracks = pt.get<int>("MaxHoughTracks");
  fPsetTPCLines.MinTrackHits = pt.get<int>("MinTrackHits");
  fPsetTPCLines.RemoveIsolatedHits = pt.get<bool>("RemoveIsolatedHits");
  fPsetTPCLines.MaxNeighbourDistance = pt.get<double>("MaxNeighbourDistance");
  fPsetTPCLines.MinNeighboursHits = pt.get<int>("MinNeighboursHits");
  fPsetTPCLines.MinTrackGoodness = pt.get<float>("MinTrackGoodness");
  fPsetTPCLines.CustomKinkPoint = pt.get<bool>("CustomKinkPoint");
  fPsetTPCLines.ParallelTracksMaxD = pt.get<int>("ParallelTracksMaxD");
  fPsetTPCLines.VertexAlgorithm = pt.get<int>("VertexAlgorithm");
  fPsetTPCLines.View = pt.get<int>("View");
  fPsetTPCLines.OutputPath = pt.get<std::string>("OutputPath");
  fPsetTPCLines.Verbose = pt.get<int>("Verbose");
  fPsetTPCLines.DebugMode = pt.get<int>("DebugMode");

  pt = fFhiclPsetCalo;
  fPsetCalo.MaxTrackHits = pt.get<int>("MaxTrackHits");
  fPsetCalo.MakeFit = pt.get<bool>("MakeFit");
  fPsetCalo.UseHitRMSInFit = pt.get<bool>("UseHitRMSInFit");
  fPsetCalo.HitWidthTolInFit = pt.get<double>("HitWidthTolInFit");
  fPsetCalo.Verbose = pt.get<int>("Verbose");

  fPsetTPCLines.HoughAlgorithmPset = fPsetHough;
  fPsetTPCLines.TrackFinderAlgorithmPset = fPsetTrackFinder;
  fPsetTPCLines.VertexFinderAlgorithmPset = fPsetVertexFinder;
  fPsetTPCLines.CaloAlgorithmPset = fPsetCalo;

  // Initialize the TPC lines algorithm
  fTPCLinesAlgo = new TPCLinesAlgo(fPsetTPCLines, fPsetFRANS);

  // Set the TPC constants
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

  fTriggerOffsetTPC = clockData.TriggerOffsetTPC(); //in us
  fTickPeriodTPC = clockData.TPCClock().TickPeriod(); //in us
  fReadoutWindow = detProp.ReadOutWindowSize();
  fDriftVelocity = detProp.DriftVelocity(); //in cm/us
  constexpr geo::TPCID tpcid{0, 0};
  fWirePlanePosition = std::abs( m_readout.Plane(geo::PlaneID{tpcid, 1}).GetCenter().X() );

  // Output data products
  produces<std::vector<recob::Hit>>();
  // Association with wires
  produces<art::Assns<recob::Wire, recob::Hit>>();
  // MCParticle for lambda
  produces<std::vector<simb::MCParticle>>();
  produces<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>();
  // Lambda products
  produces<qelambda::QELambdaReco>();


  // =========================
  //       MICHEL BDT
  // =========================

  // Find BDT weights files
  std::string fmichelModelConfig = fMichelPset.get<std::string> ("MichelBDTConfig");
  std::string michelModelPath;
  if ( !sp.find_file(fmichelModelConfig, michelModelPath) )
    throw cet::exception("FRANSProducer") << "Could not locate michel BDT configuration file: " <<
        fmichelModelConfig << " not found in FW_SEARCH_PATH\n";
  mf::LogInfo("FRANSProducer") << " Michel BDT weights file name: " << michelModelPath;

  fmichelBDT.Load(michelModelPath);
  
}


// --------- Function to get wire position
bool qelambda::FRANSProducer::GetWirePosition(int ch, double xyz_start[3], double xyz_end[3]){
  int maxCh = 11264;
  if (ch < 0 || ch >= maxCh) {   
    mf::LogError("FRANSProducer") << "Channel id, "<<  ch << " outside valid range [0, " << maxCh << ")";
    return false;
  }

  std::vector<geo::WireID> wireV = m_readout.ChannelToWire(ch);
  if (wireV.empty()) {           
    mf::LogWarning("FRANSProducer") << "No wire found for channel " << ch;
    return false;
  }
  for(size_t w=0; w<wireV.size(); w++){
    m_readout.WireEndPoints(wireV[w], xyz_start, xyz_end);
  }
  return true;
}


// -------- Function to get the recob hits in the ROI ball
std::vector<SHit> qelambda::FRANSProducer::GetFRANSHitsView(int _view, int _tpc, std::vector<art::Ptr<recob::Hit>> hitV, std::vector<int> hitClusterId, std::map<int, art::Ptr<recob::SpacePoint>> fhitToSpacePointMap)
{   
    // reset variables
    std::vector<SHit> hitList;
    
    // Get hits in the views
    for(size_t k=0; k<hitV.size(); k++){
      art::Ptr<recob::Hit> hit = hitV[k];
      if( (int)hit->WireID().TPC==_tpc && hit->View()==_view){

        float x = -999;
        float y = -999;
        float z = -999;

        if(fhitToSpacePointMap.find(hit.key())!=fhitToSpacePointMap.end()){
          x = fhitToSpacePointMap[hit.key()]->XYZ()[0];
          y = fhitToSpacePointMap[hit.key()]->XYZ()[1];
          z = fhitToSpacePointMap[hit.key()]->XYZ()[2];
        }
        
        hitList.push_back( SHit(-1,
                            (int)hit->Channel(),
                            hit->PeakTime(),
                            hit->RMS(),
                            hit->Integral(),
                            hit->StartTick(),
                            hit->EndTick(),
                            hit->GoodnessOfFit(),
                            hitClusterId[k],
                            x,
                            y,
                            z)
                          );
      }
    }

    return hitList;

}

bool qelambda::FRANSProducer::CheckClearCosmic(std::vector<art::Ptr<recob::PFParticle>> pfpVec, art::FindManyP<larpandoraobj::PFParticleMetadata> pfp2pfpmetaAssn) {
  bool isCosmic = false;
  
  for(const art::Ptr <recob::PFParticle> &pfp: pfpVec) {
    bool pfpIsCosmic = false;
    const std::vector<art::Ptr < larpandoraobj::PFParticleMetadata>> pfpMetaVec = pfp2pfpmetaAssn.at(pfp.key());
    if (pfpMetaVec.empty()) continue; 
    art::Ptr <larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();
    std::map<std::string, float> pfpPropMap = pfpMeta->GetPropertiesMap();

    if(pfpPropMap.find("IsClearCosmic") != pfpPropMap.end()) {
      pfpIsCosmic = pfpPropMap["IsClearCosmic"];
    }

    if (pfpIsCosmic) isCosmic = true;
  }

  return isCosmic;
}


// -------- Function to get the energy reconstructed as showers
bool qelambda::FRANSProducer::FindNuPfp(art::Ptr<recob::PFParticle> &nupfp,
  std::vector<art::Ptr<recob::PFParticle>> pfpVect){

  //for(size_t ix = 0; ix<pfpVect.size(); ix++ ){
  for(auto& pfp : pfpVect){
    if(pfp->IsPrimary() && ( std::abs(pfp->PdgCode())==12 || std::abs(pfp->PdgCode())==14 )  ){
      nupfp = pfp;
      break;
    }
  }

  return !nupfp.isNull();
}


// -------- Function to  convert vertex X in cm to vertex X in ticks
int qelambda::FRANSProducer::VertexToDriftTick(double vt, double vx){
  return int( ( vt/1000 + ( fWirePlanePosition-std::abs(vx) )/fDriftVelocity - fTriggerOffsetTPC)/fTickPeriodTPC );
}


// -------- Function to  convert drift tick to vertex X in cm
double qelambda::FRANSProducer::DriftTickToVertex(double tt, int tpc){
  double driftLength = (fTriggerOffsetTPC + tt * fTickPeriodTPC ) * fDriftVelocity ;
  if(tpc==0) return driftLength - fWirePlanePosition;
  else return fWirePlanePosition-driftLength;
}

VertexView qelambda::FRANSProducer::getPointInWires(const double& x, const double& y, const double& z, bool reco) {
  VertexView verView;

  //Get Vertex in wires
  int tpc = 0;
  if(x > 0)  tpc = 1;

  geo::PlaneID UPlane_ID(0,tpc,0);
  geo::PlaneID VPlane_ID(0,tpc,1);
  geo::PlaneID CPlane_ID(0,tpc,2);


  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  geo::PlaneGeo const& plane =  wireReadout.Plane(UPlane_ID);
  geo::WireID wireID;
  try {
    wireID = plane.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane.ClosestWireID(e.suggestedWireID());
  }
  //int recoVU = wireID.deepestIndex();
  int channelVU = wireReadout.PlaneWireToChannel(wireID);

  geo::PlaneGeo const& plane1 =  wireReadout.Plane(VPlane_ID);
  try {
    wireID = plane1.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane1.ClosestWireID(e.suggestedWireID());
  }
  //int recoVV = wireID.deepestIndex();
  int channelVV = wireReadout.PlaneWireToChannel(wireID);

  geo::PlaneGeo const& plane2 =  wireReadout.Plane(CPlane_ID);
  try {
    wireID = plane2.NearestWireID(geo::Point_t{x,y,z});
  } catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    wireID = plane2.ClosestWireID(e.suggestedWireID());
  }
  int channelVC = wireReadout.PlaneWireToChannel(wireID);

  int recoTimeTick = VertexToDriftTick(0, x);

  verView.x = x;
  verView.y = y;
  verView.z = z;
  verView.t = 0;
  verView.reco = reco;
  verView.U = channelVU;
  verView.V = channelVV;
  verView.C = channelVC;
  verView.TT = recoTimeTick;

  return verView;
}


// -------- Function to check if a point is in the fiducial volume
bool qelambda::FRANSProducer::PointInFV(double x, double y, double z){
  return ( std::abs(x)>fXFidCut1 && std::abs(x)<fXFidCut2 && std::abs(y)<fYFidCut && z>fZFidCut1 && z<fZFidCut2 );
}

bool qelambda::FRANSProducer::PointInFVNoCathode(double x, double y, double z){
  return (std::abs(x)<fXFidCut2 && std::abs(y)<fYFidCut && z>fZFidCut1 && z<fZFidCut2 );
}

// -------- Function to check if a point is in the active volume
bool qelambda::FRANSProducer::PointInAV(double x, double y, double z){
  return ( std::abs(x)<fXAVCut && std::abs(y)<fYAVCut && z>fZAVCut1 && z<fZAVCut2 );
}


// -------- Function to reset the variables
void qelambda::FRANSProducer::ResetVars()
{ 

  // Lambda preselection variables
  fnPrimPfps = 0;
  fnPrimTrks = 0;
  fnPrimShwrs = 0;
  fnTrksExits = 0;
  fmaxGap = 0;
  fFransReco = false;

  fHyperonVertexX = -999;
  fHyperonVertexY = -999;
  fHyperonVertexZ = -999;
  fPrimaryLeptonVertexX = -999;
  fPrimaryLeptonVertexY = -999;
  fPrimaryLeptonVertexZ = -999;
  
  fmuDecay = false;
  fmichelFeatures.clear();

}

std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<int>> qelambda::FRANSProducer::GetAssociatedHits(
  const art::Event &e,
  const std::vector<art::Ptr<recob::PFParticle>>& pfpVect,
  const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,      
  const art::Handle<std::vector<recob::Cluster>>& clusterHandle)      
{
  std::vector<art::Ptr<recob::Hit>> sliceRecobHits;
  std::vector<int> sliceRecobHitsClusterId;

  if (!pfpHandle.isValid() || !clusterHandle.isValid()) {
    mf::LogWarning("FRANSProducer") << "GetAssociatedHits: invalid PFParticle/Cluster handle";
    return {sliceRecobHits, sliceRecobHitsClusterId};
  }

  art::FindManyP<recob::Cluster> pfp_cluster_assns (pfpHandle, e, fRecoLabel);
  art::FindManyP<recob::Hit> cluster_hit_assns (clusterHandle, e, fRecoLabel);
  art::FindManyP<recob::Vertex> pfp_vertex_assns (pfpHandle, e, fRecoLabel);

  if (!pfp_cluster_assns.isValid() || !cluster_hit_assns.isValid()) {
    mf::LogWarning("FRANSProducer") << "GetAssociatedHits: PFParticle->Cluster or Cluster->Hit assns not available";
    return {sliceRecobHits, sliceRecobHitsClusterId};
  }
  const bool vertexAssnsOk = pfp_vertex_assns.isValid();   // only used for printing

  for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){
      //  Fill pfp id - pfp map
      fPfpIDToPfpMap[pfp->Self()] = pfp;
  
      // --- Read cluster and store hits
      std::vector<art::Ptr<recob::Cluster>> cluster_v = pfp_cluster_assns.at(pfp.key());
      for(size_t i=0; i<cluster_v.size(); i++){  
          std::vector<art::Ptr<recob::Hit>> hitVect = cluster_hit_assns.at(cluster_v[i].key());
          if((hitVect.size()>0)){
              sliceRecobHits.insert(sliceRecobHits.end(), hitVect.begin(), hitVect.end());
              sliceRecobHitsClusterId.insert(sliceRecobHitsClusterId.end(), hitVect.size(), pfp->Self());
          }
      }

      // --- Get the PFP vertex
      if (vertexAssnsOk) {
        std::vector<art::Ptr<recob::Vertex>> vertex_v = pfp_vertex_assns.at(pfp.key());
        if(vertex_v.size()>0){
            double pfp_vertex[3];
            vertex_v[0]->XYZ(pfp_vertex);
            // geo::Point_t clusterP={pfp_vertex[0], pfp_vertex[1], pfp_vertex[2]};
        }
      }

  }

  return {sliceRecobHits, sliceRecobHitsClusterId};
}



bool qelambda::FRANSProducer::ApplyLambdaPreselection(
  const art::Event &e,
  const art::Ptr<recob::Slice>& slice,
  const std::vector<art::Ptr<recob::PFParticle>>& slcPfpVec,
  const art::Handle<std::vector<recob::Slice>> &sliceHandle,
  const art::Handle<std::vector<recob::PFParticle>> &pfpHandle,
  const art::Handle<std::vector<recob::Shower>> &showerHandle,
  art::Ptr<recob::PFParticle>& nupfp)
{
  if (!sliceHandle.isValid() || !pfpHandle.isValid() || !showerHandle.isValid()) {
    mf::LogWarning("FRANSProducer") << "ApplyLambdaPreselection: invalid Slice/PFParticle/Shower handle";
    return false;
  }

  // Get Assns
  art::FindOneP<recob::Track> pfp2trackAssn(pfpHandle, e, fTrackLabel);
  art::FindOneP<recob::Shower> pfp2showerAssn(pfpHandle, e, fShowerLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns(pfpHandle, e, fRecoLabel);

  art::FindManyP<recob::Vertex> pfp_vertex_assns(pfpHandle, e, fRecoLabel);
  art::FindManyP<recob::Hit> shower2hitAssn(showerHandle, e, fShowerLabel);

  art::FindOneP<sbn::TPCPMTBarycenterMatch> slice2bcfmAssns(sliceHandle, e, fTPCPMTBarycenterFMModuleLabel);

  if (!pfp2trackAssn.isValid() || !pfp2showerAssn.isValid() || !pfp_pfpmeta_assns.isValid() ||
      !pfp_vertex_assns.isValid() || !shower2hitAssn.isValid()) {
    mf::LogWarning("FRANSProducer") << "ApplyLambdaPreselection: some PFParticle/Shower associations are not available";
    return false;
  }

  // Check clear cosmic
  fAnaTreeHandle.fIsClearCosmic = CheckClearCosmic(slcPfpVec, pfp_pfpmeta_assns);
  if (fAnaTreeHandle.fIsClearCosmic) return false;

  // // Get nu score
  if (!FindNuPfp(nupfp, slcPfpVec)) return false;

  // Get barycenter flash-matching (optional association: if missing, score = -1)
  float slcBcFMScore = -1;
  if (slice2bcfmAssns.isValid()) {
    const art::Ptr<sbn::TPCPMTBarycenterMatch> slcBcFM = slice2bcfmAssns.at(slice.key());
    if(slcBcFM) slcBcFMScore = slcBcFM->score;
  } else {
    mf::LogWarning("FRANSProducer") << "Barycenter flash-matching association not found (label '"
                                    << fTPCPMTBarycenterFMModuleLabel << "'), BcFMScore = -1";
  }
  fAnaTreeHandle.fBcFMScore = slcBcFMScore;
  if (slcBcFMScore < fBcFMScoreCut) return false;

  // ---- Slice in FV
  // --- Get the reco interaction vertex
  auto vertexVect = pfp_vertex_assns.at(nupfp.key());
  if (vertexVect.empty()) {    
    return false;
  }
  double slc_vertex[3];
  for(const art::Ptr<recob::Vertex> &ver : vertexVect){
      ver->XYZ(slc_vertex);
  }
  VertexView recoVertex = getPointInWires(slc_vertex[0], slc_vertex[1], slc_vertex[2], false);

  // --- Set the vertex and the TPC
  fRecoIsFiducial = PointInFVNoCathode(recoVertex.x, recoVertex.y, recoVertex.z);
  fRecoVertex = recoVertex;
  fTPC = (fRecoVertex.x<0) ? 0 : 1;
  
  // --- Reco vertex
  fAnaTreeHandle.fRecnuvX = fRecoVertex.x;
  fAnaTreeHandle.fRecnuvY = fRecoVertex.y;
  fAnaTreeHandle.fRecnuvZ = fRecoVertex.z;
  fAnaTreeHandle.fRecoIsFiducial = fRecoIsFiducial;

  if (!fRecoIsFiducial) return false;


  // Preselection cuts:
  // number of primary pfps 
  // number of primary tracks (> 0.5 track score) w lenght > trkLenCut
  // number of primary showers (< 0.5 track score) w n hits > ShwrNHitsCut
  // number of primary tracks exiting the FV
  // max gap between pairs of primary particles > fminGapCut

  std::vector<std::array<double,3>> tracksStartPos;
  for (const art::Ptr<recob::PFParticle> &pfp : slcPfpVec)
  {
    if (nupfp->Self() != pfp->Parent()) continue;

    // ---- Get track score
    double trkScore = -1;
    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfp_pfpmeta_assns.at(pfp.key());

    //If pfp_meta_vec is empty track score is -1
    if (!pfpMetaVec.empty()){
      std::map<std::string, float> pfpPropMap = pfpMetaVec.at(0)->GetPropertiesMap();
      if(!((pfpPropMap.find("TrackScore") == pfpPropMap.end()) || (pfpPropMap["TrackScore"] < 0))){
        trkScore = pfpPropMap["TrackScore"];
      }
    }

    if (trkScore < 0) continue;

    fnPrimPfps++;    

    const art::Ptr<recob::Track> track = pfp2trackAssn.at(pfp.key());
    const art::Ptr<recob::Shower> shower = pfp2showerAssn.at(pfp.key());

    if (trkScore > 0.5){ // try just > 0

      if (track.isNull()) {
        mf::LogWarning("FRANSProducer") << "\t LAMBDA PRESELECTION: Track NULL";
        continue;
      }

      double trkLen = track->Length();
      if (trkLen < fTrkLenCut) continue;
      
      fnPrimTrks++;

      // Checks if ends inside FV
      double endX = track->End().X();
      double endY = track->End().Y();
      double endZ = track->End().Z();
      if (!PointInFVNoCathode(endX, endY, endZ)) fnTrksExits++;


      // Add start position 
      double startX = track->Start().X();
      double startY = track->Start().Y();
      double startZ = track->Start().Z();

      if (startX < -200 || startX > 200 ||
          startY < -200 || startY > 200 ||
          startZ < 0 || startZ > 500) continue;

      tracksStartPos.push_back({startX, startY, startZ});


    } else {
      
      if (shower.isNull()) {
        mf::LogWarning("FRANSProducer") << "\t LAMBDA PRESELECTION: Shower NULL";
        continue;
      }

      std::map<int, int> shwrNHitMap;

      std::vector<art::Ptr<recob::Hit>> showerHitVec = shower2hitAssn.at(shower.key());

      for (const auto& hit : showerHitVec){
        shwrNHitMap[hit->WireID().Plane]++;
      }
      
      int maxHits = 0;

      if (!shwrNHitMap.empty()) {
          maxHits = std::max_element(
              shwrNHitMap.begin(),
              shwrNHitMap.end(),
              [](const auto& a, const auto& b) {
                  return a.second < b.second;
              })->second;
      }

      if (maxHits < fShwrNHits) continue;
      const int bp = shower->best_plane();
      if (bp < 0 || bp >= static_cast<int>(shower->Energy().size())) continue;
      if (shower->Energy()[bp] < fShwrECut) continue;

      fnPrimShwrs++;

      // Add starting position
      double startX = shower->ShowerStart().X();
      double startY = shower->ShowerStart().Y();
      double startZ = shower->ShowerStart().Z();

      if (startX < -200 || startX > 200 ||
          startY < -200 || startY > 200 ||
          startZ < 0 || startZ > 500) continue;

      tracksStartPos.push_back({startX, startY, startZ});
      
    }
  }

  if (fnPrimPfps < fnPrimPfpsCut1 || fnPrimPfps >= fnPrimPfpsCut2) return false;
  if (fnPrimTrks < fnPrimTrksCut1 || fnPrimTrks >= fnPrimTrksCut2) return false;
  if (fnPrimShwrs < fnPrimShwrsCut1 || fnPrimShwrs >= fnPrimShwrsCut2) return false;
  if (fnTrksExits >= fnPrimTrksExitCut) return false;


  fmaxGap = 0;

  if (tracksStartPos.size() >= 2)
  {
    for (size_t i = 0; i < tracksStartPos.size(); i++)
    {
      for (size_t j = i+1; j < tracksStartPos.size(); j++){
        double gap = std::hypot(
            tracksStartPos[i][0] - tracksStartPos[j][0],
            tracksStartPos[i][1] - tracksStartPos[j][1],
            tracksStartPos[i][2] - tracksStartPos[j][2]
        );

        if (gap > fmaxGap){
          fmaxGap = gap;
        }
      }
    }
  }

  if (fmaxGap < fminGapCut) return false;
  
  return true;
}


bool qelambda::FRANSProducer::getFRANSVertices(){
  // --- Vertex to be fed to PANDORA reprocessing 

  fFRANSslc.Triangle.GetVertexXYZ(fHyperonVertexX, fHyperonVertexY, fHyperonVertexZ);
  // Overwrite the X and Z depending on the view for the hyperon
  double xyzWire_start[3], xyzWire_end[3];
  if (!GetWirePosition(fFRANSslc.Triangle.GetMainVertex().X() + fTPCLinesAlgo->ShiftX(), xyzWire_start, xyzWire_end)) {
    return false;
  }
  double fHyperonVertexZCorrected = fHyperonVertexZ;
  if(fView==2)
    fHyperonVertexZCorrected = (xyzWire_start[2]+xyzWire_end[2])/2;
  double fHyperonVertexXCorrected = DriftTickToVertex(fFRANSslc.Triangle.GetMainVertex().Y() + fTPCLinesAlgo->ShiftY(), fTPC); 
  fHyperonVertexX = fHyperonVertexXCorrected;
  fHyperonVertexZ = fHyperonVertexZCorrected;
          
  // Primary lepton vertex
  fFRANSslc.Triangle.GetPrimaryLeptonVertexXYZ(fPrimaryLeptonVertexX, fPrimaryLeptonVertexY, fPrimaryLeptonVertexZ);

  return true;

}

art::Ptr<recob::PFParticle> qelambda::FRANSProducer::getMuonPfp(
  const art::Event &e,
  const art::Handle<std::vector<recob::PFParticle>> &pfpHandle)
{
  art::Ptr<recob::PFParticle> muon;

  if (!pfpHandle.isValid()) {
    mf::LogWarning("FRANSProducer") << "getMuonPfp: invalid PFParticle handle";
    return muon;
  }

  art::FindOneP<recob::Track> pfp_track_assns (pfpHandle, e, fTrackLabel);
  if (!pfp_track_assns.isValid()) {
    mf::LogWarning("FRANSProducer") << "getMuonPfp: PFParticle->Track assns not available";
    return muon;
  }

  std::vector<int> muPfpIDs = fFRANSslc.Triangle.GetMainTrackClusterIDs();

  // std::cout << "Main track custer id" << std::endl;
    // for (int id : muPfpIDs){
    //   // std::cout << " - " << id;
    // }
  // std::cout << std::endl;

  if (muPfpIDs.size() == 1) {
    // std::cout << "Muon found (only one pfp)" << std::endl;
    auto it = fPfpIDToPfpMap.find(muPfpIDs[0]);
    // std::cout << "Muon pfp found at fPfpIDToPfpMap: " << (it != fPfpIDToPfpMap.end()) << std::endl;

    if (it != fPfpIDToPfpMap.end()){
      fFRANSslc.MainTrackPfpID = muPfpIDs[0];
      muon = it->second;
    } 
  } else {

    double minDist = 999;
    int bestId = -1;

    for (const auto& muID : muPfpIDs){
      auto it = fPfpIDToPfpMap.find(muID);
      // std::cout << "Muon pfp found at fPfpIDToPfpMap: " << (it != fPfpIDToPfpMap.end()) << std::endl;
      if (it == fPfpIDToPfpMap.end()) continue;

      auto muPfp = it->second;

      art::Ptr<recob::Track> muTrack = pfp_track_assns.at(muPfp.key());

      if (!muTrack || muTrack.isNull()) continue;

      double dist = std::hypot(
        muTrack->Start().X() - fPrimaryLeptonVertexX,
        muTrack->Start().Y() - fPrimaryLeptonVertexY,
        muTrack->Start().Z() - fPrimaryLeptonVertexZ
      );

      if (dist < minDist){
        // std::cout << "Muon found" << std::endl;
        muon = muPfp;
        bestId = muID;
        minDist = dist;
        // std::cout << "Min muon dist: " << minDist << std::endl;
      }
    }
    fFRANSslc.MainTrackPfpID = bestId;
  }

  return muon;
}


void qelambda::FRANSProducer::getVPfps(
  const art::Event &e,
  const art::Handle<std::vector<recob::PFParticle>> &pfpHandle) 
{
  if (!pfpHandle.isValid()) return;

  art::FindOneP<recob::Track> pfp_track_assns(pfpHandle, e, fTrackLabel);
  if (!pfp_track_assns.isValid()) return;

  // Returns the PFParticle ID whose track starts closest to the primary lepton vertex
  auto selectClosestToHypVertex = [&](const std::vector<int>& pfpIDs) -> int {

    // Single candidate: nothing to select
    if (pfpIDs.size() == 1) return pfpIDs.front();

    double minDist2 = std::numeric_limits<double>::max();
    int bestId = -1;

    for (const int id : pfpIDs) {
      auto it = fPfpIDToPfpMap.find(id);
      if (it == fPfpIDToPfpMap.end()) continue;

      const art::Ptr<recob::Track> track = pfp_track_assns.at(it->second.key());
      if (track.isNull()) continue;

      // Squared distance: same ordering, no sqrt
      const auto& start = track->Start();
      const double dx = start.X() - fHyperonVertexX;
      const double dy = start.Y() - fHyperonVertexY;
      const double dz = start.Z() - fHyperonVertexZ;
      const double dist2 = dx*dx + dy*dy + dz*dz;

      if (dist2 < minDist2) {
        minDist2 = dist2;
        bestId = id;
      }
    }

    return bestId;
  };

  const std::vector<int>& track1PfpIDs = fFRANSslc.Triangle.GetTrack1ClusterIDs();
  const std::vector<int>& track2PfpIDs = fFRANSslc.Triangle.GetTrack2ClusterIDs();

  mf::LogDebug("FRANSProducer")
    << "Track1 cluster ids: " << track1PfpIDs.size()
    << ", Track2 cluster ids: " << track2PfpIDs.size();

  fFRANSslc.Track1PfpID = selectClosestToHypVertex(track1PfpIDs);
  fFRANSslc.Track2PfpID = selectClosestToHypVertex(track2PfpIDs);
}


void qelambda::FRANSProducer::updateLambdaInfo()
{

  fLambdaInfo.SliceId = fFRANSslc.SliceID;

  fLambdaInfo.FransNUnassocOrigins = fFRANSslc.LambdaAnaInfo.fNUnOrigins;
  fLambdaInfo.FransNMaxDirtUnassocHits = fFRANSslc.LambdaAnaInfo.fNMaxDirtUnassociatedHits;
  fLambdaInfo.FransNAngles = fFRANSslc.LambdaAnaInfo.fNAngles;

  fLambdaInfo.FransScore = fFRANSslc.LambdaAnaInfo.fAngleFRANSScore;
  fLambdaInfo.FransEta = fFRANSslc.LambdaAnaInfo.fAngleFRANSEta;
  fLambdaInfo.FransDelta = fFRANSslc.LambdaAnaInfo.fAngleFRANSDelta;
  fLambdaInfo.FransIota = fFRANSslc.LambdaAnaInfo.fAngleFRANSIota;
  fLambdaInfo.FransFitScore = fFRANSslc.LambdaAnaInfo.fAngleFRANSFitScore;

  fLambdaInfo.FransContainedPhi = fFRANSslc.LambdaAnaInfo.fAngleDecayContainedDiff;
  fLambdaInfo.FransPzsign = fFRANSslc.LambdaAnaInfo.fAnglePzSign;

  fLambdaInfo.FransTrack1ID = fFRANSslc.Track1PfpID;
  fLambdaInfo.FransTrack2ID = fFRANSslc.Track2PfpID;
  fLambdaInfo.FransMainTrackID = fFRANSslc.MainTrackPfpID;

  fLambdaInfo.HyperonVtxX = fbestHyperonVertexX;
  fLambdaInfo.HyperonVtxY = fbestHyperonVertexY;
  fLambdaInfo.HyperonVtxZ = fbestHyperonVertexZ;
  fLambdaInfo.PrimLeptonVtxX = fbestPrimaryLeptonVertexX;
  fLambdaInfo.PrimLeptonVtxY = fbestPrimaryLeptonVertexY;
  fLambdaInfo.PrimLeptonVtxZ = fbestPrimaryLeptonVertexZ;

  fLambdaInfo.MichelNMuonTracks = fmichelFeatures.nMuonTracks;
  fLambdaInfo.MichelRssRatio = fmichelFeatures.lowerChi2Ratio;
  fLambdaInfo.MichelAngleQ80 = fmichelFeatures.angleQ80;
  fLambdaInfo.MichelDropRadius = fmichelFeatures.dropRadius;
  fLambdaInfo.MichelNBlips = fmichelFeatures.nBlips;
  fLambdaInfo.MichelScore = fmichelFeatures.score;

}


void qelambda::FRANSProducer::setupHitMaps(
  const art::Event &e, 
  const art::Handle<std::vector<recob::SpacePoint>>& spHandle)
{
  fhitToSpacePointMap.clear();

  if (!spHandle.isValid()) {
    mf::LogWarning("FRANSProducer") << "setupHitMaps: invalid SpacePoint handle";
    return;
  }

  // Hit to sp map
  art::FindManyP<recob::Hit> SPToHitAssoc (spHandle, e, fRecoLabel);

  // Get all the sp from the event
  std::vector<art::Ptr<recob::SpacePoint>> eventSPVec;
  art::fill_ptr_vector(eventSPVec, spHandle);

  for (const art::Ptr<recob::SpacePoint> &SP: eventSPVec){
    std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());
    fhitToSpacePointMap[SPHit.at(0).key()] = SP;
  }
}

// -------- Main function
void qelambda::FRANSProducer::produce(art::Event &e)
{

  // --- Event General Info
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fEventID = e.id().event();
  std::string evLabel = std::to_string(fRunID)+"_"+std::to_string(fSubRunID)+"_"+std::to_string(fEventID);
  mf::LogInfo("FRANSProducer") << "Running FRANSProducer---"<<evLabel<<" View: "<<fView;

  // Clear FRANS slice
  fFRANSslc.clear();
  fLambdaHits.clear();
  fLambdaInfo = qelambda::QELambdaReco{};
  fLambdaInfo.SetDefault();

  fmichelFeatures.clear();

  fbestHyperonVertexX = -999;
  fbestHyperonVertexY = -999;
  fbestHyperonVertexZ = -999;
  fbestPrimaryLeptonVertexX = -999;
  fbestPrimaryLeptonVertexY = -999;
  fbestPrimaryLeptonVertexZ = -999;

  auto putEmptyProducts = [&e]() {
    qelambda::QELambdaReco emptyInfo;
    emptyInfo.SetDefault();
    e.put(std::make_unique<std::vector<recob::Hit>>());
    e.put(std::make_unique<art::Assns<recob::Wire, recob::Hit>>());
    e.put(std::make_unique<std::vector<simb::MCParticle>>());
    e.put(std::make_unique<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>());
    e.put(std::make_unique<qelambda::QELambdaReco>(emptyInfo));
  };


  // =========================
  //   Get handles
  // =========================
  art::Handle<std::vector<simb::MCTruth>>     mctruthHandle;
  art::Handle<std::vector<simb::MCParticle>>  mcparticleHandle;
  art::Handle<std::vector<recob::Slice>>      sliceHandle;
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  art::Handle<std::vector<recob::Shower>>     showerHandle;
  art::Handle<std::vector<recob::Hit>>        hitsHandle;
  art::Handle<std::vector<recob::Cluster>>    clusterHandle;
  art::Handle<std::vector<recob::SpacePoint>> spHandle;


  mctruthHandle = e.getHandle<std::vector<simb::MCTruth>>(fMCTruthLabel);
  mcparticleHandle = e.getHandle<std::vector<simb::MCParticle>>(fMCParticleLabel);
  sliceHandle = e.getHandle<std::vector<recob::Slice>>(fRecoLabel);
  pfpHandle = e.getHandle<std::vector<recob::PFParticle>>(fRecoLabel);
  showerHandle = e.getHandle<std::vector<recob::Shower>>(fShowerLabel);
  hitsHandle = e.getHandle<std::vector<recob::Hit>>(fHitLabel);
  clusterHandle = e.getHandle<std::vector<recob::Cluster>>(fRecoLabel);
  spHandle = e.getHandle<std::vector<recob::SpacePoint>>(fRecoLabel);

  if (!sliceHandle.isValid() || !pfpHandle.isValid() || !showerHandle.isValid() || !hitsHandle.isValid() || !clusterHandle.isValid() || !spHandle.isValid()) {
    mf::LogWarning("FRANSProducer") << "One or more required handles are invalid in event " << evLabel
                                    << ", writing empty products";
    putEmptyProducts();
    return;
  }

  // Associations used directly in produce()
  //Slice association for Hits
  art::FindManyP<recob::Hit> slice_hit_assns (sliceHandle, e, fRecoLabel);
  //Slice association for PFParticles
  art::FindManyP<recob::PFParticle> slice_pfp_assns (sliceHandle, e, fRecoLabel);
  //PF to track
  art::FindManyP<recob::Track> pfp_track_assns (pfpHandle, e, fTrackLabel);
  // Hit to Wires
  art::FindManyP<recob::Wire> hits_wires_assns(hitsHandle, e, fHitLabel);

  if (!slice_hit_assns.isValid() || !slice_pfp_assns.isValid() || !pfp_track_assns.isValid()) {
    mf::LogWarning("FRANSProducer") << "Slice->Hit, Slice->PFParticle or PFParticle->Track assns missing in event "
                                    << evLabel << ", writing empty products";
    putEmptyProducts();
    return;
  }
  const bool hitWireAssnsOk = hits_wires_assns.isValid();
  if (!hitWireAssnsOk) {
    mf::LogWarning("FRANSProducer") << "Hit->Wire assns missing in event " << evLabel
                                    << ", lambda hits will not be written";
  }

  std::vector<art::Ptr<recob::PFParticle>> pfpVect;


  // Setup hits map
  fPfpIDToPfpMap.clear();
  setupHitMaps(e, spHandle);


  // --- Get slice vector
  std::vector< art::Ptr<recob::Slice> > sliceVect;
  art::fill_ptr_vector(sliceVect, sliceHandle);

  fbestFRANSScore = -900;
  fbestSlice = -1;

  // --- Loop over slices
  for (auto& slice : sliceVect){
    fSliceID = slice->ID();

    // --- Reset the analysis tree
    fAnaTreeHandle.ResetVars();
    ResetVars();

    // Get slice hits
    std::vector<art::Ptr<recob::Hit>> slcHitVec = slice_hit_assns.at(slice.key());

    // Get slice pfps
    pfpVect = slice_pfp_assns.at(slice.key());
    art::Ptr<recob::PFParticle> nupfp;

    // Apply Lambda preselection to the slice to mimic the cuts done previously on Cafpyana
    if (!ApplyLambdaPreselection(e, slice, pfpVect, sliceHandle, pfpHandle, showerHandle, nupfp)) continue;

    // --- Get the associated hits to the clusters in the slice
    auto clusterHits = GetAssociatedHits(e, pfpVect, pfpHandle, clusterHandle);
    std::vector<art::Ptr<recob::Hit>> sliceRecobHits = clusterHits.first;
    std::vector<int> sliceRecobHitsClusterId = clusterHits.second; 

    // --- Get the SimpleHits for the given view and the SVertex
    std::vector<SHit> simpleHitsVect = GetFRANSHitsView(fView, fTPC, sliceRecobHits, sliceRecobHitsClusterId, fhitToSpacePointMap);
    SVertex simpleVertex( SPoint(fRecoVertex.C, fRecoVertex.TT), std::to_string( fView) );
    
    // --- Set the true and reco vertices
    std::vector<int> RecoVertexUVYT = {fRecoVertex.U, fRecoVertex.V, fRecoVertex.C, fRecoVertex.TT};
    std::vector<int> TrueVertexUVYT = RecoVertexUVYT; // REMOVE VERTEX TRUE {fTrueVertex.U, fTrueVertex.V, fTrueVertex.C, fTrueVertex.TT};

    // --- Fill with the selected hits and analyze
    bool filled = fTPCLinesAlgo->SetHitList(fView, RecoVertexUVYT, TrueVertexUVYT, simpleHitsVect);

    // Run TPCLines algorithm
    SEventId ev(fRunID, fSubRunID, fEventID);
    if(filled){
      fTPCLinesAlgo->AnaView(ev.Label());
    }


    // --- Fill FRANS variables 
    STriangle lambdaTriangleCandidate = fTPCLinesAlgo->FillLambdaAnaTree(fAnaTreeHandle);

    // --- Vertex to be fed to PANDORA reprocessing 
    // Hyperon vertex
    // std::cout << "  --- Slice FRANS score: " << fAnaTreeHandle.fAngleFRANSScore << std::endl;
    // std::cout << "  --- Previous best FRANS score: " << fbestFRANSScore << std::endl;

    fFransReco = (fAnaTreeHandle.fAngleFRANSScore > -2);

    // std::cout << "Saving FRANS slice info, sliceID: " << fSliceID << std::endl;

    // Update fFRANSslc
    fFRANSslc.SliceID = fSliceID;
    fFRANSslc.SlcHits = std::move(slcHitVec);
    fFRANSslc.SlcRecobHits = std::move(sliceRecobHits);
    fFRANSslc.SlcRecobHitsClusterId = std::move(sliceRecobHitsClusterId);
    fFRANSslc.Triangle = lambdaTriangleCandidate;
    fFRANSslc.LambdaAnaInfo = fAnaTreeHandle;

    // std::cout << "IS FRANS EXECUTED: " << fFransReco << std::endl;
    // std::cout << "GETTING SLICE TRUE INFO, FRANS VERTICES AND RUNNING MICHEL TAGGING" << std::endl;


    if (fFransReco){ 

      fvalidFransVertex = getFRANSVertices();

      if (fvalidFransVertex) {

        geo::Point_t fransVertex(fPrimaryLeptonVertexX, fPrimaryLeptonVertexY, fPrimaryLeptonVertexZ);


        // Get hits for second Pandora of best frans slice
        if (fFRANSslc.LambdaAnaInfo.fAngleFRANSScore > fbestFRANSScore){

          fbestSlice = fFRANSslc.SliceID;

          // --- Store Lambda decay pfps candidates hits for PANDORA
          std::vector<int> pfpIDs = fFRANSslc.Triangle.GetClusterIDs();
          fLambdaHits.clear();
          for (size_t kx = 0; kx < fFRANSslc.SlcRecobHits.size(); kx++){
            if (std::find(pfpIDs.begin(), pfpIDs.end(), fFRANSslc.SlcRecobHitsClusterId[kx]) != pfpIDs.end()){
              fLambdaHits.push_back(fFRANSslc.SlcRecobHits[kx]);
            }
          }

          fbestHyperonVertexX = fHyperonVertexX;
          fbestHyperonVertexY = fHyperonVertexY;
          fbestHyperonVertexZ = fHyperonVertexZ;
          fbestPrimaryLeptonVertexX = fPrimaryLeptonVertexX;
          fbestPrimaryLeptonVertexY = fPrimaryLeptonVertexY;
          fbestPrimaryLeptonVertexZ = fPrimaryLeptonVertexZ;
          // fNLambdaHits = fLambdaHits.size();
          
          // std::cout<<" NUMBERSTORELAMBDAHITS: "<<fLambdaHits.size()<<std::endl;
          fbestFRANSScore = fFRANSslc.LambdaAnaInfo.fAngleFRANSScore;

          // Get track 1 and 2 pfp ids
          getVPfps(e, pfpHandle);
          
          // --- Look for muon pfp - closest one to fPrimaryLeptonVertex
          art::Ptr<recob::PFParticle> muonPfp = getMuonPfp(e, pfpHandle);
          std::vector<art::Ptr<recob::Track>> muonTrackVec;
          if (muonPfp.isNonnull()) muonTrackVec = pfp_track_assns.at(muonPfp.key());

          if (muonTrackVec.size() == 0) {
            // std::cout << "* Muon track is empty - skipping Michel Tagging" << std::endl;
          } else{
            // std::cout << "* Muon track has lenght " << muonTrackVec.size() << ", taking first muon track" << std::endl;

            // auto muonTrack = muonTrackVec[0];

            // std::cout << "Muon PFP starting point: " << muonTrack->Start().X() << ", "<<  muonTrack->Start().Y() << ", "<<  muonTrack->Start().Z() << std::endl;
            // std::cout << "Muon PFP ending point: " << muonTrack->End().X() << ", "<<  muonTrack->End().Y() << ", "<<  muonTrack->End().Z() << std::endl;
            // std::cout << "FRANS muon starting point: " << fPrimaryLeptonVertexX << ", "<<  fPrimaryLeptonVertexY << ", "<<  fPrimaryLeptonVertexZ << std::endl;


            // ======== MICHEL BDT

            // std::cout << "*** Muon pfp input: " << muonPfp.key() << std::endl;

            auto michelExtractor = sbnd::MichelExtractor(fMichelPset, fmichelBDT);
            michelExtractor.Extract(e, fFRANSslc.SliceID, muonPfp, fransVertex);

            fmichelFeatures = michelExtractor.getMichelFeatures();
            
              // std::cout << "Michel features computed: " << fmichelFeatures.featuresComputed << std::endl;
              // std::cout << "Michel score: " << fmichelFeatures.score << std::endl;
              // std::cout << "True muon decay: " << fmuDecay << std::endl; 
          }

          updateLambdaInfo();

        }
      }
    }

  }


  // ===========================================================
  //              FILL ART EVENT
  // ===========================================================

  // --- Add hyperon hits and association between recob::Hit and recob::Wire for second stage reconstruction
  // std::cout << "  --- Adding hyperon hits and associations" << std::endl;

  std::unique_ptr<std::vector<recob::Hit>> lambdaHitsVec =
      std::make_unique<std::vector<recob::Hit>>();
  std::unique_ptr<art::Assns<recob::Wire, recob::Hit>> lambdaHitsWireAssns =
      std::make_unique<art::Assns<recob::Wire, recob::Hit>>();

  std::unique_ptr<qelambda::QELambdaReco> lambdaProducts = 
      std::make_unique<qelambda::QELambdaReco>();
  *lambdaProducts = fLambdaInfo;

  if (hitWireAssnsOk && fLambdaHits.size() > 0) {
    for (auto const& hit : fLambdaHits) {
      if (hit.key() >= hits_wires_assns.size()) {
        mf::LogWarning("FRANSProducer") << "Invalid hit key: " << hit.key();
        continue;
      }

      auto const& wire = hits_wires_assns.at(hit.key());

      if (wire.empty()) {
         mf::LogWarning("FRANSProducer") << "No wire associated to hit key: " << hit.key();
        continue;
      }

      lambdaHitsVec->push_back(*hit);
      util::CreateAssn(e, hit, wire[0], *lambdaHitsWireAssns);
    }
  }

  // --- Add reconstructed hyperon vertex for second stage reconstruction
  auto pPartColFirst = std::make_unique<std::vector<simb::MCParticle>>();
  auto tpassn = std::make_unique<art::Assns<simb::MCTruth, simb::MCParticle, sim::GeneratedParticleInfo>>();

  // Awful messy trick for now
  // std::cout << "  --- Second round of MCParticles" << std::endl;

  if (!mcparticleHandle.isValid()) {
    mf::LogWarning("FRANSProducer") << "Invalid MCParticle collection, writing empty MCParticle products";
  } else {
    auto partCol = std::make_unique<std::vector<simb::MCParticle>>();

    sim::MergeSimSourcesUtility MergeUtility{{0}};
    art::PtrMaker<simb::MCParticle> const makePartPtr{e};

    MergeUtility.MergeMCParticles(*partCol, *mcparticleHandle, 0);

    const auto& listMap = MergeUtility.GetMCParticleListMap();

    art::FindOneP<simb::MCTruth, sim::GeneratedParticleInfo> mctAssn(mcparticleHandle, e, fMCParticleLabel);

    if (partCol->empty() || listMap.empty() || listMap.at(0).empty() ||
        !mctAssn.isValid() || mctAssn.size() == 0) {
      mf::LogWarning("FRANSProducer") << "Could not build MCParticle/MCTruth products, writing empty ones";
    } else {
      const auto& assocVectorPrimitive = listMap.at(0);

      // Add just the first entry
      tpassn->addSingle(mctAssn.at(0), makePartPtr(assocVectorPrimitive[0]), mctAssn.data(0).ref());

      pPartColFirst->push_back(partCol->at(0));

      // Add the hyperon vertex candidate
      pPartColFirst->front().AddTrajectoryPoint(
          TLorentzVector(fbestHyperonVertexX, fbestHyperonVertexY, fbestHyperonVertexZ, 0),
          TLorentzVector(0, 0, 0, 0));

      // Print front particle start/end points
      // std::cout << "  --- Front Particle Start Point: "
      //           << pPartColFirst->front().Vx() << " "
      //           << pPartColFirst->front().Vy() << " "
      //           << pPartColFirst->front().Vz() << std::endl;
      // std::cout << "  --- Front Particle End Point: "
      //           << pPartColFirst->front().EndX() << " "
      //           << pPartColFirst->front().EndY() << " "
      //           << pPartColFirst->front().EndZ() << std::endl;
    }
  }

  // std::cout << "  *** lambdaProducts: " << std::endl;
  // std::cout << "  ---   SliceID: " << lambdaProducts->SliceId << std::endl;
  // std::cout << "  ---   FRANS Score: " << lambdaProducts->FransScore << std::endl;
  // std::cout << "  --- Michel Score: " << lambdaProducts->MichelScore << std::endl;
  

  // --- Fill the art event
  e.put(std::move(lambdaHitsVec));
  e.put(std::move(lambdaHitsWireAssns));
  e.put(std::move(pPartColFirst));
  e.put(std::move(tpassn));
  e.put(std::move(lambdaProducts));
}