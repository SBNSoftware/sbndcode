///////////////////////////////////////////////////////////////////////
// Class:       FlashPredict
// Plugin Type: producer (art v3_01_02)
// File:        FlashPredict_module.cc
//
// Generated at Mon May 20 07:11:36 2019 by David Caratelli using cetskelgen
// from cetlib version v3_05_01.
//
// Ported from MicroBooNE to SBND on 2019 Sept 19
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
//#include "larpandora/LArPandoraEventBuilding/LArPandoraSliceIdHelper.h"
//#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
//#include "larpandora/LArPandoraEventBuilding/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "OpT0FinderTypes.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"

#include <memory>
#include <string>

class FlashPredict;
class FlashPredict : public art::EDProducer {

public:
  explicit FlashPredict(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  // Plugins should not be copied or assigned.
  FlashPredict(FlashPredict const&) = delete;
  FlashPredict(FlashPredict&&) = delete;
  FlashPredict& operator=(FlashPredict const&) = delete;
  FlashPredict& operator=(FlashPredict&&) = delete;
  // Required functions.
  void produce(art::Event& e) override;
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.
  //  ::flashana::FlashMatchManager m_flashMatchManager; ///< The flash match manager
  // art::InputTag fFlashProducer;
  // art::InputTag fT0Producer; // producer for ACPT in-time anab::T0 <-> recob::Track assocaition
  std::string fPandoraProducer, fSpacePointProducer, fOpHitProducer, fInputFilename, fCaloProducer, fTrackProducer;
  float fBeamWindowEnd, fBeamWindowStart;
  float fLightWindowEnd, fLightWindowStart;
  float fMinFlashPE;
  float fPEscale;
  float fChargeToNPhotonsShower, fChargeToNPhotonsTrack;
  std::string fDetector; // SBND or ICARUS
  int fCryostat;  // =0 or =1 to match ICARUS reco chain selection
  bool fMakeTree,fSelectNeutrino, fUseCalo;
  std::vector<float> fPMTChannelCorrection;
  // geometry service
  const uint nMaxTPCs = 4;
  flashana::QCluster_t lightCluster[4];


  ::flashana::Flash_t GetFlashPESpectrum(const recob::OpFlash& opflash);
  void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                    const art::Ptr<recob::PFParticle> &particle,
                                    lar_pandora::PFParticleVector &downstreamPFParticles) const;
  void CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                    const lar_pandora::PFParticleVector &parentPFParticles,
                                    lar_pandora::PFParticleVector &downstreamPFParticles) const;
  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
                    std::vector<art::Ptr<recob::PFParticle> > &pfp_v);
  bool isPDInCryoTPC(float pd_x, int icryo,int itpc, std::string detector);
  bool isChargeInCryoTPC(float qp_x, int icryo,int itpc, std::string detector);

  // root stuff
  TTree* _flashmatch_acpt_tree;
  TTree* _flashmatch_nuslice_tree;
  TH1F *ophittime;
						
  // Tree variables
  std::vector<float> _pe_reco_v, _pe_hypo_v;
  float _trk_vtx_x, _trk_vtx_y, _trk_vtx_z, _trk_end_x, _trk_end_y, _trk_end_z;
  float _nuvtx_x, _nuvtx_y, _nuvtx_z, _nuvtx_q;
  float _flash_x,_flash_y,_flash_z,_flash_pe,_flash_unpe;
  float _flash_r, _score;
  int _evt, _run, _sub;
  float _flashtime;
  float _flashpe;
  // PFP map
  std::map<unsigned int, unsigned int> _pfpmap;

  std::vector<float> dysp,dzsp,rrsp,pesp,dymean,dzmean,rrmean,pemean;
  int rr_nbins,dy_nbins,dz_nbins,pe_nbins;

};

FlashPredict::FlashPredict(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
    // More initializers here.
{
  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::PFParticle, anab::T0> >();
  // fFlashProducer      = p.get<art::InputTag>("FlashProducer" );
  fOpHitProducer      = p.get<std::string>("OpHitProducer","ophit");
  fPandoraProducer    = p.get<std::string>("PandoraProducer" ,"pandora"  );
  fTrackProducer         = p.get<std::string>("TrackProducer","pandoraTrack"      );
  fCaloProducer         = p.get<std::string>("CaloProducer","pandoraCalo"      );
  fSpacePointProducer = p.get<std::string>("SpacePointProducer", "pandora" );
  fBeamWindowStart = p.get<float>("BeamWindowStart", 0.0);
  fBeamWindowEnd   = p.get<float>("BeamWindowEnd", 4000.0);  // in ns
  fMinFlashPE      = p.get<float>("MinFlashPE", 0.0);
  fChargeToNPhotonsShower   = p.get<float>("ChargeToNPhotonsShower", 1.0);  // ~40000/1600
  fChargeToNPhotonsTrack    = p.get<float>("ChargeToNPhotonsTrack", 1.0);   // ~40000/1600
  fInputFilename = p.get<std::string>("InputFileName","fmplots.root");  // root file with histograms for match score calc
  fMakeTree = p.get<bool>("MakeTree",false);  
  fUseCalo = p.get<bool>("UseCalo",false);  
  fSelectNeutrino = p.get<bool>("SelectNeutrino",true);  
  fLightWindowStart = p.get<float>("LightWindowStart", -0.010);  // in us w.r.t. flash time
  fLightWindowEnd   = p.get<float>("LightWindowEnd", 0.090);  // in us w.r.t flash time
  fDetector = p.get<std::string>("Detector", "SBND");
  fCryostat = p.get<int>("Cryostat",0); //set =0 ot =1 for ICARUS to match reco chain selection
  fPEscale = p.get<float>("PEscale",1.0);

  if (fDetector == "SBND" && fCryostat==1){
    throw cet::exception("FlashPredictSBND") << "SBND has only one cryostat. \n"
                                             << "Check Detector and Cryostat parameter." << std::endl;
  }

  art::ServiceHandle<art::TFileService> tfs;

  int nbins = 500*(fBeamWindowEnd-fBeamWindowStart);
  ophittime = tfs->make<TH1F>("ophittime","ophittime",nbins,fBeamWindowStart,fBeamWindowEnd); // in us
  ophittime->SetOption("HIST");

  if (fMakeTree) {
  _flashmatch_nuslice_tree = tfs->make<TTree>("nuslicetree","nu FlashPredict tree");
  _flashmatch_nuslice_tree->Branch("evt",&_evt,"evt/I");
  _flashmatch_nuslice_tree->Branch("run",&_run,"run/I");
  _flashmatch_nuslice_tree->Branch("sub",&_sub,"sub/I");
  _flashmatch_nuslice_tree->Branch("flashtime",&_flashtime,"flashtime/F");
  _flashmatch_nuslice_tree->Branch("flashpe"  ,&_flash_pe  ,"flashpe/F");
  _flashmatch_nuslice_tree->Branch("flash_x"  ,&_flash_x  ,"flash_x/F");
  _flashmatch_nuslice_tree->Branch("flash_y"  ,&_flash_y  ,"flash_y/F");
  _flashmatch_nuslice_tree->Branch("flash_z"  ,&_flash_z  ,"flash_z/F");
  _flashmatch_nuslice_tree->Branch("flash_r"  ,&_flash_r  ,"flash_r/F");
  _flashmatch_nuslice_tree->Branch("flash_unpe"  ,&_flash_unpe  ,"flash_unpe/F");
  _flashmatch_nuslice_tree->Branch("nuvtx_q",&_nuvtx_q,"nuvtx_q/F");
  _flashmatch_nuslice_tree->Branch("nuvtx_x",&_nuvtx_x,"nuvtx_x/F");
  _flashmatch_nuslice_tree->Branch("nuvtx_y",&_nuvtx_y,"nuvtx_y/F");
  _flashmatch_nuslice_tree->Branch("nuvtx_z",&_nuvtx_z,"nuvtx_z/F");
  _flashmatch_nuslice_tree->Branch("score",&_score,"score/F");
  }

  //read histograms and fill vectors for match score calculation
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fInputFilename, fname);
  //  std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
  TFile *infile = new TFile(fname.c_str(), "READ");
  if(!infile->IsOpen())
    {
      throw cet::exception("FlashPredictSBND") << "Could not find the light-charge match root file '" << fname << "'!\n";
    }
  //
  TH1 *temphisto = (TH1*)infile->Get("rrp1");
  rr_nbins = temphisto->GetNbinsX();
  if (rr_nbins<=0) {
    std::cout << " problem with input histos for rr " << rr_nbins << " bins " << std::endl;
    rr_nbins=1;
    rrmean.push_back(0);
    rrsp.push_back(0.001);
  }
  else {
  for (int ib=1;ib<=rr_nbins;++ib) {
    rrmean.push_back(temphisto->GetBinContent(ib));
    float tt = temphisto->GetBinError(ib);
    if (tt<=0) { tt=100.; std::cout << "zero value for bin spread in rr" << std::endl;}
    rrsp.push_back(tt);
  }
  }
  //
  temphisto = (TH1*)infile->Get("dyp1");
  dy_nbins = temphisto->GetNbinsX();
  if (dy_nbins<=0) {
    std::cout << " problem with input histos for dy " << dy_nbins << " bins " << std::endl;
    dy_nbins=1;
    dymean.push_back(0);
    dysp.push_back(0.001);
  }
  else {
  for (int ib=1;ib<=dy_nbins;++ib) {
    dymean.push_back(temphisto->GetBinContent(ib));
    float tt = temphisto->GetBinError(ib);
    if (tt<=0) { tt=100.; std::cout << "zero value for bin spread in dy" << std::endl;}
    dysp.push_back(tt);
  }
  }
  //
  temphisto = (TH1*)infile->Get("dzp1");
  dz_nbins = temphisto->GetNbinsX();
  if (dz_nbins<=0) {
    std::cout << " problem with input histos for dz " << dz_nbins << " bins " << std::endl;
    dz_nbins=1;
    dzmean.push_back(0);
    dzsp.push_back(0.001);
  }
  else {
  for (int ib=1;ib<=dz_nbins;++ib) {
    dzmean.push_back(temphisto->GetBinContent(ib));
    float tt = temphisto->GetBinError(ib);
    if (tt<=0) { tt=100.; std::cout << "zero value for bin spread in dz" << std::endl;}
    dzsp.push_back(tt);
  }
  }
  //
  if (fDetector == "SBND" ) {
    temphisto = (TH1*)infile->Get("pep1");
    pe_nbins = temphisto->GetNbinsX();
    if (pe_nbins<=0) {
      std::cout << " problem with input histos for pe " << pe_nbins << " bins " << std::endl;
      pe_nbins=1;
      pemean.push_back(0);
      pesp.push_back(0.001);
    }
    else {
      for (int ib=1;ib<=pe_nbins;++ib) {
	pemean.push_back(temphisto->GetBinContent(ib));
	float tt = temphisto->GetBinError(ib);
	if (tt<=0) { tt=100.; std::cout << "zero value for bin spread in pe" << std::endl;}
	pesp.push_back(tt);
      }
    }
  }
  else {
    pe_nbins=1;
    pemean.push_back(0);
    pesp.push_back(0.001);
  }
  //
  infile->Close();

  // Call appropriate produces<>() functions here.

  // Call appropriate consumes<>() for any products to be retrieved by this module.
} // FlashPredict::FlashPredict(fhicl::ParameterSet const& p)

void FlashPredict::produce(art::Event & e)
{

  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::PFParticle, anab::T0> > pfp_t0_assn_v( new art::Assns<recob::PFParticle, anab::T0>  );

  // SBND map for light detector type labels (pmt, uncoatedpmt, arapuca)
  opdet::sbndPDMapAlg map;

  // reset TTree variables
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  _flashtime = -9999.;
  _flashpe   = -9999.;
  _flash_unpe   = -9999.;
  _flash_r   = -9999.;
  _score = -9999.;

  // std::cout << "event no " << _evt << std::endl;

  const art::ServiceHandle<geo::Geometry> geometry;
  uint nTPCs(geometry->NTPC());
  if (nTPCs>2) {
    std::cout << "nTPC can't be larger than 2, resizing." << std::endl;
    nTPCs = 2;
  }
  geo::CryostatGeo geo_cryo = geometry->Cryostat(0);
  if (fDetector == "ICARUS" && fCryostat==1) geo_cryo = geometry->Cryostat(1);


  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraProducer);

  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPandoraProducer);  

  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointProducer);
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointProducer);

  // grab tracks associated with PFParticles
  // auto const& track_h = e.getValidHandle<std::vector<recob::Track> >(fTrackProducer);
  // art::FindManyP<recob::Track> pfp_track_assn_v(track_h, e, fTrackProducer);

  // grab calorimetry info for tracks
  // auto const& calo_h = e.getValidHandle<std::vector<anab::Calorimetry> >(fCaloProducer);
  // art::FindManyP<anab::Calorimetry>  track_calo_assn_v(calo_h, e, fCaloProducer);


  // load OpHits previously created
  art::Handle<std::vector<recob::OpHit> > ophit_h;
  e.getByLabel(fOpHitProducer,ophit_h);
  if(!ophit_h.isValid()){
    mf::LogWarning("FlashPredict")
      <<"No optical hits from producer module "<<fOpHitProducer;
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }
  std::vector<recob::OpHit> const& OpHitCollection(*ophit_h);

  _pfpmap.clear();
  for (unsigned int p=0; p < pfp_h->size(); p++)    _pfpmap[pfp_h->at(p).Self()] = p;

  // get flash time
  ophittime->Reset();
  for(size_t j = 0; j < OpHitCollection.size(); j++){
    recob::OpHit oph = OpHitCollection[j];
    double PMTxyz[3];
    geometry->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
    if ( fDetector == "SBND" && !map.pdType(oph.OpChannel(),"pmt")) continue; // use only uncoated PMTs for SBND for flashtime
    if (!geo_cryo.ContainsPosition(PMTxyz)) continue;   // use only PMTs in the specified cryostat for ICARUS
    //    std::cout << "op hit " << j << " channel " << oph.OpChannel() << " time " << oph.PeakTime() << " pe " << fPEscale*oph.PE() << std::endl;
    if ( (oph.PeakTime()<fBeamWindowStart) || (oph.PeakTime()> fBeamWindowEnd) ) continue;

       ophittime->Fill(oph.PeakTime(),fPEscale*oph.PE());
    // float thisPE = fPEscale*oph.PE();
    // if (thisPE>1) ophittime->Fill(oph.PeakTime(),thisPE);
  }

  if (ophittime->GetEntries()<=0 || ophittime->Integral() < fMinFlashPE) {
    e.put(std::move(T0_v));
    e.put(std::move(pfp_t0_assn_v));
    return;
  }    

  auto ibin =  ophittime->GetMaximumBin();
  float flashtime = (ibin*0.002)+fBeamWindowStart;  // in us
  float lowedge = flashtime+fLightWindowStart;
  float highedge = flashtime+fLightWindowEnd;
    std::cout << "light window " << lowedge << " " << highedge << std::endl;

  // Loop over pandora pfp particles
  for (unsigned int p=0; p < pfp_h->size(); p++){
    auto const& pfp = pfp_h->at(p);
    if (pfp.IsPrimary() == false) continue;
    if ( fSelectNeutrino && (abs(pfp.PdgCode())!=12) && (abs(pfp.PdgCode())!=14) )  continue; 

    lightCluster[0].clear();
    lightCluster[1].clear();

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;    
    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);

    //  loop over mothers and daughters, fill lightcluster for each tpc    
    for (size_t i=0; i < pfp_ptr_v.size(); i++) {
      const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
      auto key = pfp_ptr_v.at(i).key();
      recob::PFParticle pfp = *pfp_ptr_v.at(i);
      pfp_v.push_back(pfp);
      /*      
      if ( fUseCalo && lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr)) {
	// grab tracks associated with pfp particle
	auto const& track_ptr_v = pfp_track_assn_v.at(key);
	for (size_t tr=0; tr < track_ptr_v.size(); tr++) {
	  auto mytrack = track_ptr_v[tr];
	  auto const& trackkey = mytrack.key();
	  // grab calo objects associated with tracks
	  const std::vector< art::Ptr<anab::Calorimetry> > calo_ptr_v = track_calo_assn_v.at( trackkey );
	  for (size_t ca=0;  ca <  calo_ptr_v.size(); ca++) {
	    auto mycalo = calo_ptr_v.at( ca );
	    int npts = mycalo->dEdx().size();
	    for (int ip=0;ip<npts;++ip) {
	      Point_t pxyz=mycalo->fXYZ[ip];
	      float ds = mycalo->fTrkPitch[ip];
	      float dQdx = mycalo->fdQdx[ip];
	      float dEdx = mycalo->fdEdx[ip];
	      float alpha = dQdx/dEdx;
	      float charge = (1-alpha)*dEdx*ds;
	      // hardcode for now for SBND	      
	      float xpos = 0.0;
	      if (pxyz[0]<0) xpos = fabs(pxyz[0]+200.0);
	      else xpos = fabs(200.0-pxyz[0]);
	      lightCluster[tpcindex].emplace_back(xpos, position[1], position[2], charge);
	    }
	  }
	}
      }
      else { // this is a shower
      */
	auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
	std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
        // TODO: refactor this loop over spacepoints so that it's not
        // necessary to query the wire position every time.
        // There's just two different X wire positions on any given
        // cryostat for the collection wires.
	for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
	  auto SP = spacepoint_ptr_v[sp];
	  auto const& spkey = SP.key();
	  const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
	  for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
	    auto hit = this_hit_ptr_v.at( h );
	    // Only use hits from the collection plane
	    geo::WireID wid = hit->WireID();
	    if (geometry->SignalType(wid) != geo::kCollection) continue;
	    // Add the charged point to the vector
	    const auto &position(SP->XYZ());
            const auto tpcindex = wid.TPC;
            // throw the charge coming from another TPC
            if (!isChargeInCryoTPC(position[0], fCryostat, tpcindex, fDetector)) continue;
	    const auto charge(hit->Integral());
	    double Wxyz[3];
	    geometry->WireIDToWireGeo(wid).GetCenter(Wxyz);
            // xpos is the distance from the wire planes.
	    float xpos = fabs(position[0]-Wxyz[0]);
	    lightCluster[tpcindex].emplace_back(xpos, position[1], position[2], charge * (lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr) ? fChargeToNPhotonsTrack : fChargeToNPhotonsShower));
	  } // for all hits associated to this spacepoint
	} // for all spacepoints
	//      }  // if track or shower
    } // for all pfp pointers
    
    int icountPE=0;
    float mscore[2] ={0}; 
    float charge[2]={0};
    for (size_t itpc=0; itpc<nTPCs; ++itpc) {
      double xave=0.0; double yave=0.0; double zave=0.0; double norm=0.0;
      _nuvtx_q = 0;
      for (size_t i=0; i<lightCluster[itpc].size(); ++i) {
        flashana::QCluster_t this_cl = lightCluster[itpc];
        flashana::QPoint_t qp = this_cl[i];
        xave += 0.001*qp.q*qp.x;
        yave += 0.001*qp.q*qp.y;
        zave += 0.001*qp.q*qp.z;
        norm += 0.001*qp.q;
        _nuvtx_q += qp.q;
      }
      if (norm>0) {
	_nuvtx_x=xave/norm;
	_nuvtx_y=yave/norm;
	_nuvtx_z=zave/norm;
        charge[itpc] = _nuvtx_q;
	// store PMT photon counts in the tree as well
	double PMTxyz[3];
	double unpe_tot=0;
	double sumy=0; double sumz=0; double pnorm=0;
	double sum_Ay=0; double sum_Az=0;
	double sum_By=0; double sum_Bz=0;
	double sum_D=0;  double sum=0;
	double sum_Cy=0;double sum_Cz=0;
        // TODO: change this next loop, such that it only loops
        // through channels in the current fCryostat
	for(size_t j = 0; j < OpHitCollection.size(); j++){
	  recob::OpHit oph = OpHitCollection[j];
          std::string op_type = "pmt";
          if (fDetector == "SBND") op_type = map.pdName(oph.OpChannel());
	  geometry->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
          // check cryostat and tpc
          if (!isPDInCryoTPC(PMTxyz[0], fCryostat, itpc, fDetector)) continue;
	  // only use optical hits around the flash time
	  if ( (oph.PeakTime()<lowedge) || (oph.PeakTime()>highedge) || oph.PE()<=0 ) continue;
	  // only use PMTs for SBND
          if (op_type == "pmt") {
	    // Add up the position, weighting with PEs
	    _flash_x=PMTxyz[0];
	    pnorm+=oph.PE();
	    sum+=1.0;
	    sumy    += oph.PE()*PMTxyz[1];
	    sumz    += oph.PE()*PMTxyz[2];
	    sum_By  += PMTxyz[1];
	    sum_Bz  += PMTxyz[2];
	    sum_Ay  += oph.PE()*PMTxyz[1]*oph.PE()*PMTxyz[1];
	    sum_Az  += oph.PE()*PMTxyz[2]*oph.PE()*PMTxyz[2];
	    sum_D   +=oph.PE()*oph.PE();
	    sum_Cy  +=oph.PE()*oph.PE()*PMTxyz[1];
	    sum_Cz  +=oph.PE()*oph.PE()*PMTxyz[2];
	  }
          else if ( op_type == "barepmt") {
            unpe_tot+=oph.PE();
          }
          else if ( (op_type == "arapucaT1" || op_type == "arapucaT2") ) {
            //TODO: Use ARAPUCA
            // arape_tot+=oph.PE();
            continue;
          }
          else if ( op_type == "xarapucaT1" || op_type == "xarapucaT2")  {
            //TODO: Use XARAPUCA
            // xarape_tot+=oph.PE();
            continue;
          }
	}
	
	if (pnorm>0) {
	  _flashtime=flashtime;
	  _flash_pe=pnorm*fPEscale;
	  _flash_y=sum_Cy/sum_D;
	  _flash_z=sum_Cz/sum_D;
	  sum_By=_flash_y;        sum_Bz=_flash_z;
	  _flash_r=sqrt((sum_Ay-2.0*sum_By*sum_Cy+sum_By*sum_By*sum_D+sum_Az-2.0*sum_Bz*sum_Cz+sum_Bz*sum_Bz*sum_D)/sum_D);
	  _flash_unpe=unpe_tot*fPEscale;
	  icountPE  +=(int)(_flash_pe);
	}
	else { _flash_pe=0; _flash_y=0; _flash_z=0; _flash_unpe=0; _flash_r=0;}
	
	//      calculate match score here, put association on the event
	float slice=_nuvtx_x;
        float drift_distance = 200.0; // TODO: no hardcoded values
        if (fDetector == "ICARUS") {drift_distance=150.0;} // TODO: no hardcoded values
	_score = 0; int icount =0;
	int isl = int(dy_nbins*(slice/drift_distance));
        if (dysp[isl]>0) {_score+=abs(_flash_y-_nuvtx_y-dymean[isl])/dysp[isl];}
	icount++;
	isl = int(dz_nbins*(slice/drift_distance));
        if (dzsp[isl]>0) {_score+=abs(_flash_z-_nuvtx_z-dzmean[isl])/dzsp[isl];}
	icount++;
	isl = int(rr_nbins*(slice/drift_distance));
        if (rrsp[isl]>0) {_score+=abs(_flash_r-rrmean[isl])/rrsp[isl];}
	icount++;
	if (0) {  // this piece is broken?
	  //if (fDetector == "SBND") { // pe metric for sbnd only
	  isl = int(pe_nbins*(slice/drift_distance));	  
	  float myratio = 100.0*_flash_unpe;
          if (_flash_pe>1 && pesp[isl]>0) {
	    myratio/=_flash_pe;
	    _score+=abs(myratio-pemean[isl])/pesp[isl];
	    icount++;
	  }
	}
	//	_score/=icount;
	if (_flash_pe>0 ) {
	  mscore[itpc] = _score;
	  // fill tree
	  if (fMakeTree) _flashmatch_nuslice_tree->Fill();
	}
      } // if tpc charge>0
    }  // end loop over TPCs
    
    double this_score=0.0; int icount=0; double totc=0;
    for (size_t itpc=0; itpc<nTPCs; ++itpc) {
      this_score += mscore[itpc];
      totc += charge[itpc];
      if (mscore[itpc] > 0) icount++;
    }
    if (icount>0) { 
      this_score/=(icount*1.0);
      
      // create t0 and pfp-t0 association here
      T0_v->push_back(anab::T0( flashtime, icountPE, p, 0, this_score));
      //    util::CreateAssn(*this, e, *T0_v, pfp_h[p], *pfp_t0_assn_v);
      //    util::CreateAssn(*this, e, *T0_v, pfp, *pfp_t0_assn_v);
      util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);
    }
  } // over all PFparticles


  e.put(std::move(T0_v));
  e.put(std::move(pfp_t0_assn_v));



}// end of producer module

::flashana::Flash_t FlashPredict::GetFlashPESpectrum(const recob::OpFlash& opflash)
{
  // prepare container to store flash
  ::flashana::Flash_t flash;
  flash.time = opflash.Time();
  // geometry service
  const art::ServiceHandle<geo::Geometry> geometry;
  uint nOpDets(geometry->NOpDets());
  std::vector<float> PEspectrum;
  PEspectrum.resize(nOpDets);
  // apply gain to OpDets
  for (uint OpChannel = 0; OpChannel < nOpDets; ++OpChannel)
  {
    uint opdet = geometry->OpDetFromOpChannel(OpChannel);
    PEspectrum[opdet] = opflash.PEs().at(OpChannel);
  }
  _pe_reco_v = PEspectrum;

  // Reset variables
  flash.x = flash.y = flash.z = 0;
  flash.x_err = flash.y_err = flash.z_err = 0;
  float totalPE = 0.;
  float sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
  for (unsigned int opdet = 0; opdet < PEspectrum.size(); opdet++)
  {
    double PMTxyz[3];
    geometry->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
    // Add up the position, weighting with PEs
    sumy += PEspectrum[opdet] * PMTxyz[1];
    sumy2 += PEspectrum[opdet] * PMTxyz[1] * PMTxyz[1];
    sumz += PEspectrum[opdet] * PMTxyz[2];
    sumz2 += PEspectrum[opdet] * PMTxyz[2] * PMTxyz[2];
    totalPE += PEspectrum[opdet];
  }

  flash.y = sumy / totalPE;
  flash.z = sumz / totalPE;
  // This is just sqrt(<x^2> - <x>^2)
  if ((sumy2 * totalPE - sumy * sumy) > 0.)
    flash.y_err = std::sqrt(sumy2 * totalPE - sumy * sumy) / totalPE;
  if ((sumz2 * totalPE - sumz * sumz) > 0.)
    flash.z_err = std::sqrt(sumz2 * totalPE - sumz * sumz) / totalPE;
  // Set the flash properties
  flash.pe_v.resize(nOpDets);
  flash.pe_err_v.resize(nOpDets);
  // Fill the flash with the PE spectrum
  for (unsigned int i = 0; i < nOpDets; ++i)
  {
    const auto PE(PEspectrum.at(i));
    flash.pe_v.at(i) = PE;
    flash.pe_err_v.at(i) = std::sqrt(PE);
  }
  if (flash.pe_v.size() != nOpDets)
    throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;
  return flash;

}// ::flashana::Flash_t FlashPredict::GetFlashPESpectrum


void FlashPredict::CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                                const art::Ptr<recob::PFParticle> &particle,
                                                lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
    downstreamPFParticles.push_back(particle);
  for (const auto &daughterId : particle->Daughters())
  {
    const auto iter(pfParticleMap.find(daughterId));
    if (iter == pfParticleMap.end())
      throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;
    this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
  }
} // void FlashPredict::CollectDownstreamPFParticles

void FlashPredict::CollectDownstreamPFParticles(const lar_pandora::PFParticleMap &pfParticleMap,
                                                const lar_pandora::PFParticleVector &parentPFParticles,
                                                lar_pandora::PFParticleVector &downstreamPFParticles) const
{
  for (const auto &particle : parentPFParticles)
    this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
} // void FlashPredict::CollectDownstreamPFParticles

void FlashPredict::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v) {
  
  auto daughters = pfp_ptr->Daughters();
  pfp_v.push_back(pfp_ptr);
  for(auto const& daughterid : daughters) {
    if (_pfpmap.find(daughterid) == _pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }
    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, _pfpmap.at(daughterid) );
    AddDaughters(pfp_ptr, pfp_h, pfp_v);
  } // for all daughters
  return;
} // void FlashPredict::AddDaughters

// TODO: no hardcoding
// TODO: collapse with the next
bool FlashPredict::isPDInCryoTPC(float pd_x, int icryo, int itpc, std::string detector)
{
  // check whether this optical detector views the light inside this tpc.
  std::ostringstream lostPDMessage;
  lostPDMessage << "\nThere's an " << detector << "photo detector that belongs nowhere. \n"
                 << "icryo: " << icryo << "\n"
                 << "itpc: "  << itpc << "\n"
                 << "pd_x: " << pd_x << std::endl;

  if (detector == "ICARUS") {
    if (icryo == 0) {
      if (itpc == 0 && -400 < pd_x && pd_x < -300 ) return true;
      else if (itpc == 1 && -100 < pd_x && pd_x < 0) return true;
      // else {std::cout << lostPDMessage.str(); return false;}
    }
    else if (icryo == 1) {
      if (itpc == 0 && 0 < pd_x && pd_x < 100) return true;
      else if (itpc == 1 && 300 < pd_x && pd_x < 400) return true;
      // else {std::cout << lostPDMessage.str(); return false;}
    }
  }
  else if (detector == "SBND") {
    if ((itpc == 0 && -213. < pd_x && pd_x < 0) || (itpc == 1 && 0 < pd_x && pd_x < 213) ) return true;
    //    else {std::cout << lostPDMessage.str(); return false;}
    else { return false;}
  }
  return false;
}

// TODO: no hardcoding
// TODO: collapse with the previous
// TODO: figure out what to do with the charge that falls into the crevices 
bool FlashPredict::isChargeInCryoTPC(float qp_x, int icryo, int itpc, std::string detector)
{
  std::ostringstream lostChargeMessage;
  lostChargeMessage << "\nThere's " << detector << " charge that belongs nowhere. \n"
                 << "icryo: " << icryo << "\n"
                 << "itpc: "  << itpc << "\n"
                 << "qp_x: " << qp_x << std::endl;

  if (detector == "ICARUS") {
    if (icryo == 0) {
      if (itpc == 0 && -368.49 <= qp_x && qp_x <= -220.29 ) return true;
      else if (itpc == 1 && -220.14 <= qp_x && qp_x <= -71.94) return true;
      // else {std::cout << lostChargeMessage.str(); return false;}
    }
    else if (icryo == 1) {
      if (itpc == 0 && 71.94 <= qp_x && qp_x <= 220.14) return true;
      else if (itpc == 1 && 220.29 <= qp_x && qp_x <= 368.49) return true;
      // else {std::cout << lostChargeMessage.str(); return false;}
    }
  }
  else if (detector == "SBND") {
    if ((itpc == 0 && qp_x < 0) || (itpc == 1 && qp_x > 0) ) return true;
    else {return false;}
    //    else {std::cout << lostChargeMessage.str(); return false;}
  }
  return false;
}

void FlashPredict::beginJob()
{
  // Implementation of optional member function here.




}

void FlashPredict::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FlashPredict)
