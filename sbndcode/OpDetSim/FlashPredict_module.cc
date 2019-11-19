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

#include <iostream>
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
  std::string fPandoraProducer, fSpacePointProducer, fOpHitProducer, fInputFilename;
  float fBeamWindowEnd, fBeamWindowStart;
  float fLightWindowEnd, fLightWindowStart;
  float fMaxTotalPE;
  float fChargeToNPhotonsShower, fChargeToNPhotonsTrack;
  int fDetector; // ==0 for SBND, =1 ICARUS
  bool fMakeTree,fSelectNeutrino;
  std::vector<float> fPMTChannelCorrection;
  // geometry service
  //const uint nMaxTPCs = 4;

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

  // root stuff
  //TTree* _flashmatch_acpt_tree;
  TTree* _flashmatch_nuslice_tree;
  TH1F *ophittime;
						
  // Tree variables
  std::vector<float> _pe_reco_v, _pe_hypo_v;
  //float _trk_vtx_x, _trk_vtx_y, _trk_vtx_z, _trk_end_x, _trk_end_y, _trk_end_z;
  float _nuvtx_x, _nuvtx_y, _nuvtx_z, _nuvtx_q;
  float _flash_x,_flash_y,_flash_z,_flash_pe;
  float _flash_r, _score;
  int _evt, _run, _sub;
  float _flashtime;
  float _flashpe;
  // PFP map
  std::map<unsigned int, unsigned int> _pfpmap;

  std::vector<float> dysp,dzsp,rrsp,dymean,dzmean,rrmean;
  int rr_nbins,dy_nbins,dz_nbins;
  
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
  //   fT0Producer         = p.get<std::string>("T0Producer"      );
  fSpacePointProducer = p.get<std::string>("SpacePointProducer", "pandora" );
  fBeamWindowStart = p.get<float>("BeamWindowStart", 0.0);
  fBeamWindowEnd   = p.get<float>("BeamWindowEnd", 4000.0);  // in ns
  fMaxTotalPE      = p.get<float>("MaxTotalPE", 500000000.0);
  fChargeToNPhotonsShower   = p.get<float>("ChargeToNPhotonsShower", 1.0);  // ~40000/1600
  fChargeToNPhotonsTrack    = p.get<float>("ChargeToNPhotonsTrack", 1.0);   // ~40000/1600
  fInputFilename = p.get<std::string>("InputFileName","fmplots.root");  // root file with histograms for match score calc
  fMakeTree = p.get<bool>("MakeTree",false);  
  fSelectNeutrino = p.get<bool>("SelectNeutrino",true);  
  fLightWindowStart = p.get<float>("LightWindowStart", -0.010);  // in us w.r.t. flash time
  fLightWindowEnd   = p.get<float>("LightWindowEnd", 0.090);  // in us w.r.t flash time
  fDetector = p.get<int>("Detector",0); // ==0 for SBND, ==1 for ICARUS
  art::ServiceHandle<art::TFileService> tfs;

  ophittime = tfs->make<TH1F>("ophittime","ophittime",1100,-0.2,2.0); // in us

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
  TH1 *temphisto = (TH1*)infile->Get("rrp");
  rr_nbins = temphisto->GetNbinsX();
  if (rr_nbins<=0) {
    std::cout << " problem with input histos for rr " << rr_nbins << " bins " << std::endl;
    rr_nbins=1;
    rrmean.push_back(0);
    rrsp.push_back(0.001);
  }
  else {
  for (int ib=0;ib<rr_nbins;++ib) {
    rrmean.push_back(temphisto->GetBinContent(ib));
    rrsp.push_back(temphisto->GetBinError(ib));
  }
  }
  //
  temphisto = (TH1*)infile->Get("dyp");
  dy_nbins = temphisto->GetNbinsX();
  if (dy_nbins<=0) {
    std::cout << " problem with input histos for dy " << dy_nbins << " bins " << std::endl;
    dy_nbins=1;
    dymean.push_back(0);
    dysp.push_back(0.001);
  }
  else {
  for (int ib=0;ib<dy_nbins;++ib) {
    dymean.push_back(0.0);
    // not enough stats to believe mean values different than 0
    //    dymean.push_back(temphisto->GetBinContent(ib));
    dysp.push_back(temphisto->GetBinError(ib));
  }
  }
  //
  temphisto = (TH1*)infile->Get("dzp");
  dz_nbins = temphisto->GetNbinsX();
  if (dz_nbins<=0) {
    std::cout << " problem with input histos for dz " << dz_nbins << " bins " << std::endl;
    dz_nbins=1;
    dzmean.push_back(0);
    dzsp.push_back(0.001);
  }
  else {
  for (int ib=0;ib<dz_nbins;++ib) {
    dzmean.push_back(0.0);
    // not enough stats to believe mean values different than 0
    //    dzmean.push_back(temphisto->GetBinContent(ib));
    dzsp.push_back(temphisto->GetBinError(ib));
  }
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


  opdet::sbndPDMapAlg map;

  // reset TTree variables
  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  _flashtime = -9999.;
  _flashpe   = -9999.;

  // std::cout << "event no " << _evt << std::endl;

  const art::ServiceHandle<geo::Geometry> geometry;
  uint nTPCs(geometry->NTPC());
  //  int nOpDets(geometry->NOpDets());

  // std::cout << "nopdets " << nOpDets << std::endl;

  // _flashtime = beamflash.time;
  // _flashpe   = std::accumulate( beamflash.pe_v.begin(), beamflash.pe_v.end(), 0);


  // grab PFParticles in event
  auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPandoraProducer);

  // grab spacepoints associated with PFParticles
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfp_h, e, fPandoraProducer);

  // grab tracks associated with PFParticles
  art::FindManyP<recob::Track> pfp_track_assn_v(pfp_h, e, fPandoraProducer);

  // grab vertices associated with PFParticles
  // art::FindManyP<recob::Vertex> pfp_vertex_assn_v(pfp_h, e, fPandoraProducer);

  // grab associated metadata
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h, e, fPandoraProducer);
  auto const& spacepoint_h = e.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointProducer);
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, e, fSpacePointProducer);

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


  //get better access to the data
  std::vector<recob::OpHit> const& OpHitCollection(*ophit_h);
  // std::vector<recob::OpHit> FlashHits;


  _pfpmap.clear();
  for (unsigned int p=0; p < pfp_h->size(); p++)
    _pfpmap[pfp_h->at(p).Self()] = p;

  flashana::QCluster_t lightCluster[4];

  // get flash time
  ophittime->Reset();
  for(size_t j = 0; j < OpHitCollection.size(); j++){
    recob::OpHit oph = OpHitCollection[j];
    if ( (fDetector==0) && map.pdType(oph.OpChannel(),"pmt")) {
      if ( (oph.PeakTime()>fBeamWindowStart) && (oph.PeakTime()< fBeamWindowEnd) ) {
	ophittime->Fill(oph.PeakTime(),100*oph.PE());
      }
    }
    else if (fDetector == 1) {
      if ( (oph.PeakTime()>fBeamWindowStart) && (oph.PeakTime()< fBeamWindowEnd) ) {
        ophittime->Fill(oph.PeakTime(),100*oph.PE());
      }
    }
  }
  auto ibin =  ophittime->GetMaximumBin();
  float flashtime = (ibin*0.002)-0.2;  // in us
  float lowedge = flashtime-fLightWindowStart;
  float highedge = flashtime+fLightWindowEnd;
  

  // Loop over pandora pfp particles, select primary particles identified as neutrinos
  for (unsigned int p=0; p < pfp_h->size(); p++){
    auto const& pfp = pfp_h->at(p);
    if (pfp.IsPrimary() == false) continue;
    if (fSelectNeutrino && (pfp.PdgCode()!=12) && (pfp.PdgCode()!=14) ) continue; 

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, p);
    // now build vectors of PFParticles, space-points, and hits for this slice
    std::vector<recob::PFParticle> pfp_v;
    std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
    std::vector<std::vector<art::Ptr<recob::SpacePoint> > > spacepoint_v_v;
    std::vector<std::vector<art::Ptr<recob::Hit> > > hit_v_v;

    AddDaughters(pfp_ptr, pfp_h, pfp_ptr_v);

    // go through these pfparticles and fill info needed for matching
    for (size_t i=0; i < pfp_ptr_v.size(); i++) {
      auto key = pfp_ptr_v.at(i).key();
      recob::PFParticle pfp = *pfp_ptr_v.at(i);
      pfp_v.push_back(pfp);

      //auto const& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      const std::vector< art::Ptr<recob::SpacePoint> >& spacepoint_ptr_v = pfp_spacepoint_assn_v.at(key);
      std::vector< art::Ptr<recob::Hit> > hit_ptr_v;
      for (size_t sp=0; sp < spacepoint_ptr_v.size(); sp++) {
        auto SP = spacepoint_ptr_v[sp];
        auto const& spkey = SP.key();
        const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = spacepoint_hit_assn_v.at( spkey );
        for (size_t h=0; h < this_hit_ptr_v.size(); h++) {
          auto hit = this_hit_ptr_v.at( h );

          // Only use hits from the collection plane
          if (hit->View() != geo::kZ)
            continue;
          // Add the charged point to the vector
          const auto &position(SP->XYZ());
          const auto charge(hit->Integral());
          const auto tpcindex = (hit->WireID()).TPC;
          // const auto cryoindex = (hit->WireID()).Cryostat;
          // std::cout << " tpc index" << tpcindex << " x pos of sp " << position[0] << std::endl;
          lightCluster[tpcindex].emplace_back(position[0], position[1], position[2], charge * (lar_pandora::LArPandoraHelper::IsTrack(pfp_ptr) ? fChargeToNPhotonsTrack : fChargeToNPhotonsShower));
        } // for all hits associated to this spacepoint
      } // for all spacepoints
    } // for all pfp pointers

  // std::cout << "check " << lightCluster.size() << std::endl;
    float mscore[2] ={0}; float charge[2]={0};
    for (size_t it=0;it<nTPCs;++it) {
      double xave=0.0; double yave=0.0; double zave=0.0; double norm=0.0;
      _nuvtx_q=0;
      for (size_t i=0;i<lightCluster[it].size();++i) {
	flashana::QCluster_t this_cl = lightCluster[it];
	flashana::QPoint_t qp = this_cl[i];
	xave+=0.001*qp.q*qp.x;
	yave+=0.001*qp.q*qp.y;
	zave+=0.001*qp.q*qp.z;
	norm+=0.001*qp.q;
	_nuvtx_q+=qp.q;
      }
      if (norm>0) {
	_nuvtx_x=xave/norm;
	_nuvtx_y=yave/norm;
	_nuvtx_z=zave/norm;
	charge[it]=_nuvtx_q;
	// store PMT photon counts in the tree as well
	double PMTxyz[3];
	double sumy=0; double sumz=0; double pnorm=0;
	double sum_Ay=0; double sum_Az=0;
	double sum_By=0; double sum_Bz=0;
	double sum_D=0;  double sum=0;
	double sum_Cy=0;double sum_Cz=0;
	for(size_t j = 0; j < OpHitCollection.size(); j++){
	  recob::OpHit oph = OpHitCollection[j];
          if ( (fDetector==0) && map.pdType(oph.OpChannel(),"pmt")) {
	    if ( (oph.PeakTime()>lowedge) && (oph.PeakTime()< highedge) ) {
	      geometry->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
	      // std::cout << oph.OpChannel() << "   " << oph.PeakTime() << "   " << oph.PE()
	      //  <<  " PMT pos:   " << PMTxyz[0] <<"    " << PMTxyz[1] <<"    " << PMTxyz[2] << std::endl;
	      if ((it==0 && PMTxyz[0]<0) || (it==1 && PMTxyz[0]>0) ){
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
	    }
	  }
          else if (fDetector == 1) {
          if ( (oph.PeakTime()>lowedge) && (oph.PeakTime()< highedge) ) {
            geometry->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
            if ((it==0 && PMTxyz[0]<0) || (it==1 && PMTxyz[0]>0) ){ 
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
          }
        }
      }
	
	if (pnorm>0) {
	  _flashtime=flashtime;
	  _flash_pe=pnorm;
	  _flash_y=sum_Cy/sum_D;
	  _flash_z=sum_Cz/sum_D;
	  sum_By=_flash_y;        sum_Bz=_flash_z;
	  _flash_r=sqrt((sum_Ay-2.0*sum_By*sum_Cy+sum_By*sum_By*sum_D+sum_Az-2.0*sum_Bz*sum_Cz+sum_Bz*sum_Bz*sum_D)/sum_D);
	}
	else { _flash_pe=0; _flash_y=0; _flash_z=0;}
	
	//      calculate match score here, put association on the event
	float slice = (200.0-abs(_nuvtx_x));
	_score = 0;
	//      std::cout << "slice " << slice <<  " x pos " << _nuvtx_x << std::endl;
	int isl = int(dy_nbins*(slice/200.0));
	_score+=abs(_flash_y-_nuvtx_y-dymean[isl])/dysp[isl];
	//      std::cout << "bin " << isl <<  " " << dymean[isl] << " " << dysp[isl] << " score " << _score << std::endl;
	isl = int(dz_nbins*(slice/200.0));
	_score+=abs(_flash_z-_nuvtx_z-dzmean[isl])/dzsp[isl];
	//      std::cout << "bin " << isl <<  " " << dzmean[isl] << " " << dzsp[isl] << " score " << _score << std::endl;
	isl = int(rr_nbins*(slice/200.0));
	_score+=abs(_flash_r-rrmean[isl])/rrsp[isl];
	//      std::cout << "bin " << isl <<  " " << rrmean[isl] << " " << rrsp[isl] << " score " << _score << std::endl;
	mscore[it]=_score;
	if (fMakeTree) _flashmatch_nuslice_tree->Fill();
	
      } // if tpc charge>0
    }  // end loop over TPCs
    // fill tree
    
    // top line just to get code to compile
    double this_score=charge[0]*mscore[0]+charge[1]*mscore[1];
    this_score=mscore[0]+mscore[1];
    if (mscore[0]>0 && mscore[1]>0) this_score*=0.5;
    // create t0 and pfp-t0 association here

    T0_v->push_back(anab::T0( flashtime, 0, p, 0, this_score));
    //    util::CreateAssn(*this, e, *T0_v, pfp_h[p], *pfp_t0_assn_v);
    //    util::CreateAssn(*this, e, *T0_v, pfp, *pfp_t0_assn_v);
    util::CreateAssn(*this, e, *T0_v, pfp_ptr, *pfp_t0_assn_v);

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

void FlashPredict::beginJob()
{
  // Implementation of optional member function here.




}

void FlashPredict::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FlashPredict)
