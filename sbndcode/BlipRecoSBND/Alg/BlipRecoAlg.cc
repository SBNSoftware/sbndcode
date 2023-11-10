#include "sbndcode/BlipRecoSBND/Alg/BlipRecoAlg.h"

namespace blip {

  //###########################################################
  // Constructor
  //###########################################################
  BlipRecoAlg::BlipRecoAlg( fhicl::ParameterSet const& pset )
  : fGeom { *lar::providerFrom<geo::Geometry>() }
  {
    this->reconfigure(pset);
    
    auto const& detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    kLArDensity           = detProp.Density();
    kNominalEfield        = detProp.Efield();
    kDriftVelocity        = detProp.DriftVelocity(detProp.Efield(0),detProp.Temperature()); 
    kTickPeriod           = clockData.TPCClock().TickPeriod();
    kNominalRecombFactor  = ModBoxRecomb(fCalodEdx,kNominalEfield);
    kWion                 = 1000./util::kGeVToElectrons;
   
    // -------------------------------------------------------------------
    // Determine number cryostats, TPC, planes, wires.
    //
    // Also cache all the X tick offsets so we don't have to keep re-calculating them 
    // for every single hit. Note that 'detProp.GetXTicksOffset()' does not make intuitive 
    // sense for cases with wireplanes aren't at X~0 (i.e., SBND).  We will account 
    // for that here so that our calculated drift times for hits makes sense.
      

    // Loop over cryostats
    std::cout<<"NCryostats: "<<fGeom.Ncryostats()<<"\n";
    for(size_t cstat=0; cstat<fGeom.Ncryostats(); cstat++){
      auto const& cryoid = geo::CryostatID(cstat);

      // Loop TPCs in cryostat 'cstat'
      for(size_t tpc=0; tpc<fGeom.NTPC(cryoid); tpc++){
        auto const& tpcid = geo::TPCID(cryoid,tpc);

        // Loop planes in TPC 'tpc'
        for(size_t pl=0; pl<fGeom.Nplanes(tpcid); pl++){
          auto const& planeid = geo::PlaneID(cstat,tpc,pl);
            
          float offset = detProp.GetXTicksOffset(pl,tpc,cstat);
          std::cout<<"CRYOSTAT "<<cstat<<" / TPC "<<tpc<<" / PLANE "<<pl<<":  "<<fGeom.Nwires(planeid)<<" wires\n";
          std::cout<<"  XTicksOffset (from detProp): "<<offset<<"\n";
         
          kXTicksOffsets[cstat][tpc][pl] = fTimeOffset[pl];
          if( fApplyXTicksOffset ) {
            
            //kXTicksOffsets[cstat][tpc][pl] = offset;

            // subtract out the geometric time offset added to account for the 
            // distance between plane0 and X=0. This is based on code in
            // lardataalg/DetectorInfo/DetectorPropertiesStandard.cxx 
            // (as of lardataalg v9_15_01)
            auto const& cryostat  = fGeom.Cryostat(geo::CryostatID(cstat));
            auto const& tpcgeom   = cryostat.TPC(tpc);
            auto const xyz        = tpcgeom.Plane(0).GetCenter();
            const double dir((tpcgeom.DriftDirection() == geo::kNegX) ? +1.0 : -1.0);
            float x_ticks_coefficient = kDriftVelocity*kTickPeriod;
            
            float goofy_offset = -xyz.X() / (dir * x_ticks_coefficient);
            std::cout<<"  After geometric correction: "<<offset - goofy_offset<<"\n";

            kXTicksOffsets[cstat][tpc][pl] += offset - goofy_offset;
          }


        }
      }
    }

    /*
    std::cout<<"XticksOffset, Plane 0: "<<detProp.GetXTicksOffset(0,1,0)<<"\n";
    std::cout<<"XticksOffset, Plane 1: "<<detProp.GetXTicksOffset(1,1,0)<<"\n";
    std::cout<<"XticksOffset, Plane 2: "<<detProp.GetXTicksOffset(2,1,0)<<"\n";
    std::cout<<"XticksOffset, Plane 0: "<<detProp.GetXTicksOffset(0,0,0)<<"\n";
    std::cout<<"XticksOffset, Plane 1: "<<detProp.GetXTicksOffset(1,0,0)<<"\n";
    std::cout<<"XticksOffset, Plane 2: "<<detProp.GetXTicksOffset(2,0,0)<<"\n";
    */
    
    /*
    // initialize channel list
    fBadChanMask       .resize(8256,false);
    fBadChanMaskPerEvt = fBadChanMask;
    if( fBadChanFile != "" ) {
      cet::search_path sp("FW_SEARCH_PATH");
      std::string fullname;
      sp.find_file(fBadChanFile,fullname);
      if (fullname.empty()) {
        throw cet::exception("Bad channel list not found");
      } else {
        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        while (std::getline(inFile,line)) {
          if( line.find("#") != std::string::npos ) continue;
          std::istringstream ss(line);
          int ch1, ch2;
          ss >> ch1;
          if( !(ss >> ch2) ) ch2 = ch1;
          for(int i=ch1; i<=ch2; i++) fBadChanMask[i] = true;
        }
      }
    }
    int NBadChansFromFile     = std::count(fBadChanMask.begin(),fBadChanMask.end(),true);
    */

    EvtBadChanCount = 0;

    printf("******************************************\n");
    printf("Initializing BlipRecoAlg...\n");
    printf("  - Efield: %.4f kV/cm\n",kNominalEfield);
    printf("  - Drift velocity: %.4f cm/us\n",kDriftVelocity);
    printf("  - using dE/dx: %.2f MeV/cm\n",fCalodEdx);
    printf("  - equiv. recomb: %.4f\n",kNominalRecombFactor);
    //printf("  - custom bad chans: %i\n",NBadChansFromFile);
    printf("*******************************************\n");

    // create diagnostic histograms
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory hdir = tfs->mkdir("BlipRecoAlg");
   
    /*
    h_chanstatus     = hdir.make<TH1D>("chanstatus","Channel status for 'channels' list",5,0,5);
    h_chanstatus     ->GetXaxis()->SetBinLabel(1, "disconnected");
    h_chanstatus     ->GetXaxis()->SetBinLabel(2, "dead");
    h_chanstatus     ->GetXaxis()->SetBinLabel(3, "lownoise");
    h_chanstatus     ->GetXaxis()->SetBinLabel(4, "noisy");
    h_chanstatus     ->GetXaxis()->SetBinLabel(5, "good");
    
    h_hit_chanstatus     = hdir.make<TH1D>("hit_chanstatus","Channel status of hits",5,0,5);
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(1, "disconnected");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(2, "dead");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(3, "lownoise");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(4, "noisy");
    h_hit_chanstatus     ->GetXaxis()->SetBinLabel(5, "good");
    */

    //h_hit_times       = hdir.make<TH1D>("hit_peaktime","Hit peaktimes",500,-5000,5000);
    h_chan_nhits      = hdir.make<TH1D>("chan_nhits","Untracked hits;TPC readout channel;Total hits",8256,0,8256);
    h_chan_nclusts    = hdir.make<TH1D>("chan_nclusts","Untracked isolated hits;TPC readout channel;Total clusts",8256,0,8256);
    h_chan_bad        = hdir.make<TH1D>("chan_bad","Channels marked as bad;TPC readout channel",8256,0,8256);
    h_recomb          = hdir.make<TH1D>("recomb","Applied recombination factor",150,0.40,0.70);
    h_clust_nwires    = hdir.make<TH1D>("clust_nwires","Clusters (pre-cut);Wires in cluster",100,0,100);
    h_clust_timespan  = hdir.make<TH1D>("clust_timespan","Clusters (pre-cut);Time span [ticks]",300,0,300);

    int qbins = 200;
    float qmax = 100;
    //int wiresPerPlane[3]={2400,2400,3456};
    for(int i=0; i<kNplanes; i++) {
      //h_hit_maskfrac[i]       = dir_diag.make<TH1D>(Form("pl%i_hit_maskfrac",i),"",100,0,1.);
      //h_hit_maskfrac_true[i]  = dir_diag.make<TH1D>(Form("pl%i_hit_maskfrac_true",i),"",100,0,1.);
      //h_hit_mult[i]         = hdir.make<TH1D>(Form("pl%i_hit_mult",i),      Form("Plane %i;Num same-wire hits within +/- 50 ticks",i),20,0,20);
      if( i == fCaloPlane ) continue;
      //h_wire_nhits[i]       = hdir.make<TH1D>(Form("pl%i_wire_nhits",i),      Form("Plane %i untracked hits not plane-matched;TPC readout channel;Total hits",i),wiresPerPlane[i],0,wiresPerPlane[i]);
      h_clust_overlap[i]    = hdir.make<TH1D>(Form("pl%i_clust_overlap",i),   Form("Plane %i clusters;Overlap fraction",i),101,0,1.01);
      h_clust_dt[i]         = hdir.make<TH1D>(Form("pl%i_clust_dt",i),        Form("Plane %i clusters;dT [ticks]",i),200,-10,10);
      h_clust_dtfrac[i]     = hdir.make<TH1D>(Form("pl%i_clust_dtfrac",i),    Form("Plane %i clusters;Charge-weighted mean dT/RMS",i),120,-3,3);
      h_clust_q[i]     = hdir.make<TH2D>(Form("pl%i_clust_charge",i),  
        Form("Pre-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3} e-]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_q[i]    ->SetOption("colz");
      h_clust_q_cut[i]     = hdir.make<TH2D>(Form("pl%i_clust_charge_cut",i),  
        Form("Post-cut;Plane %i cluster charge [#times10^{3} e-];Plane %i cluster charge [#times10^{3}]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_q_cut[i]    ->SetOption("colz");
      
      h_clust_score[i]    = hdir.make<TH1D>(Form("pl%i_clust_matchscore",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);

      h_clust_picky_overlap[i]   = hdir.make<TH1D>(Form("pl%i_clust_picky_overlap",i),  Form("Plane %i clusters (3 planes, intersect #Delta cut);Overlap fraction",i),101,0,1.01);
      h_clust_picky_dt[i]        = hdir.make<TH1D>(Form("pl%i_clust_picky_dt",i),       Form("Plane %i clusters (3 planes, intersect #Delta cut);dT [ticks]",i),200,-10,10);
      h_clust_picky_dtfrac[i]      = hdir.make<TH1D>(Form("pl%i_clust_picky_dtfrac",i),Form("Plane %i clusters (3 planes, intersect #Delta cut);Charge-weighted mean dT/RMS",i),120,-3,3);
      h_clust_picky_q[i]  = hdir.make<TH2D>(Form("pl%i_clust_picky_charge",i),  
        Form("3 planes, intersect #Delta cut;Plane %i cluster charge [#times 10^{3} e-];Plane %i cluster charge [#times 10^{3} e-]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_picky_q[i]     ->SetOption("colz");
      
      h_clust_truematch_overlap[i]   = hdir.make<TH1D>(Form("pl%i_clust_truematch_overlap",i),  Form("Plane %i clusters (MC match);Overlap fraction",i),101,0,1.01);
      h_clust_truematch_dt[i]        = hdir.make<TH1D>(Form("pl%i_clust_truematch_dt",i),       Form("Plane %i clusters (MC match);dT [ticks]",i),200,-10,10);
      h_clust_truematch_dtfrac[i]      = hdir.make<TH1D>(Form("pl%i_clust_truematch_dtfrac",i),Form("Plane %i clusters (MC match);Charge-weighted mean dT/RMS",i),120,-3,3);
      h_clust_truematch_q[i]  = hdir.make<TH2D>(Form("pl%i_clust_truematch_charge",i),  
        Form("3 planes, intersect #Delta cut ;Plane %i cluster charge [#times 10^{3} e-];Plane %i cluster charge [#times 10^{3} e-]",fCaloPlane,i),
        qbins,0,qmax,qbins,0,qmax);
      h_clust_truematch_q[i]     ->SetOption("colz");
      h_clust_truematch_score[i]    = hdir.make<TH1D>(Form("pl%i_clust_truematch_matchscore",i),   Form("Plane %i clusters;Match score",i),101,0,1.01);
      
      
      h_nmatches[i]         = hdir.make<TH1D>(Form("pl%i_nmatches",i),Form("number of plane%i matches to single collection cluster",i),20,0,20);
    }
  
    // Efficiency as a function of energy deposited on a wire
    h_recoWireEff_denom = hdir.make<TH1D>("recoWireEff_trueCount","Collection plane;Electron energy deposited on wire [MeV];Count",150,0,1.5);
    h_recoWireEff_num   = hdir.make<TH1D>("recoWireEff","Collection plane;Electron energy deposited on wire [MeV];Hit reco efficiency",150,0,1.5);
    
    h_recoWireEffQ_denom = hdir.make<TH1D>("recoWireEffQ_trueCount","Collection plane;Charge deposited on wire [e-];Count",80,0,20000);
    h_recoWireEffQ_num   = hdir.make<TH1D>("recoWireEffQ","Collection plane;Charge deposited on wire [e-];Hit reco efficiency",80,0,20000);

  }
  
  //--------------------------------------------------------------
  //BlipRecoAlg::BlipRecoAlg( )
  //{
  //}
  
  //--------------------------------------------------------------  
  //Destructor
  BlipRecoAlg::~BlipRecoAlg()
  {
  }
  
  
  //###########################################################
  // Reconfigure fcl parameters
  //###########################################################
  void BlipRecoAlg::reconfigure( fhicl::ParameterSet const& pset ){
    
    fHitProducer        = pset.get<std::string>   ("HitProducer",       "gaushit");
    fTrkProducer        = pset.get<std::string>   ("TrkProducer",       "pandora");
    fGeantProducer      = pset.get<std::string>   ("GeantProducer",     "largeant");
    fSimDepProducer     = pset.get<std::string>   ("SimEDepProducer",   "ionization");
    fSimChanProducer    = pset.get<std::string>   ("SimChanProducer",   "driftWC:simpleSC");
    fSimGainFactor      = pset.get<float>         ("SimGainFactor",     -9);
    fTrueBlipMergeDist  = pset.get<float>         ("TrueBlipMergeDist", 0.3);
    fMaxHitTrkLength    = pset.get<float>               ("MaxHitTrkLength", 5);
    fDoHitFiltering     = pset.get<bool>                ("DoHitFiltering",  false);
    fMaxHitMult         = pset.get<int>                 ("MaxHitMult",      10);
    fMaxHitAmp          = pset.get<float>               ("MaxHitAmp",       200);  
    fMinHitAmp          = pset.get<std::vector<float>>  ("MinHitAmp",       {-99e9,-99e9,-99e9});
    fMaxHitRMS          = pset.get<std::vector<float>>  ("MaxHitRMS",       { 99e9, 99e9, 99e9});
    fMinHitRMS          = pset.get<std::vector<float>>  ("MinHitRMS",       {-99e9,-99e9,-99e9});
    fMaxHitRatio        = pset.get<std::vector<float>>  ("MaxHitRatio",     { 99e9, 99e9, 99e9});
    fMinHitRatio        = pset.get<std::vector<float>>  ("MinHitRatio",     {-99e9,-99e9,-99e9});
    fMaxHitGOF          = pset.get<std::vector<float>>  ("MaxHitGOF",       { 99e9, 99e9, 99e9});
    fMinHitGOF          = pset.get<std::vector<float>>  ("MinHitGOF",       {-99e9,-99e9,-99e9});
    
    fHitClustWidthFact  = pset.get<float>         ("HitClustWidthFact", 5.0);
    fHitClustWireRange  = pset.get<int>           ("HitClustWireRange", 1);
    fMaxWiresInCluster  = pset.get<int>           ("MaxWiresInCluster", 10);
    fMaxClusterSpan     = pset.get<float>         ("MaxClusterSpan",    30);
    fMinClusterCharge   = pset.get<float>         ("MinClusterCharge",  300);
    fMaxClusterCharge   = pset.get<float>         ("MaxClusterCharge",  12e6);
    
    fApplyXTicksOffset  = pset.get<bool>          ("ApplyXTicksOffset", true);
    fTimeOffset         = pset.get<std::vector<float>>("TimeOffset", {0.,0.,0.});
    fMatchMinOverlap    = pset.get<float>         ("ClustMatchMinOverlap",  0.5 );
    fMatchSigmaFact     = pset.get<float>         ("ClustMatchSigmaFact",   1.0);
    fMatchMaxTicks      = pset.get<float>         ("ClustMatchMaxTicks",    5.0 );
    fMatchQDiffLimit    = pset.get<float>         ("ClustMatchQDiffLimit",  15e3);
    fMatchMaxQRatio     = pset.get<float>         ("ClustMatchMaxQRatio",   4);
    
    fMinMatchedPlanes   = pset.get<int>           ("MinMatchedPlanes",    2);
    fPickyBlips         = pset.get<bool>          ("PickyBlips",          false);
    fApplyTrkCylinderCut= pset.get<bool>          ("ApplyTrkCylinderCut", false);
    fCylinderRadius     = pset.get<float>         ("CylinderRadius",      15);
    
    fCaloAlg            = new calo::CalorimetryAlg( pset.get<fhicl::ParameterSet>("CaloAlg") );
    fCaloPlane          = pset.get<int>           ("CaloPlane",           2);
    fCalodEdx           = pset.get<float>         ("CalodEdx",            2.8);
    fLifetimeCorr       = pset.get<bool>          ("LifetimeCorrection",  false);
    fSCECorr            = pset.get<bool>          ("SCECorrection",       false);
    fYZUniformityCorr   = pset.get<bool>          ("YZUniformityCorrection",true);
    fModBoxA            = pset.get<float>         ("ModBoxA",             0.93);
    fModBoxB            = pset.get<float>         ("ModBoxB",             0.212);
    
    fVetoBadChannels    = pset.get<bool>          ("VetoBadChannels",     true);
    fBadChanProducer    = pset.get<std::string>   ("BadChanProducer",     "nfspl1:badchannels");
    fBadChanFile        = pset.get<std::string>   ("BadChanFile",         "");
    fMinDeadWireGap     = pset.get<int>           ("MinDeadWireGap",      1);
    
    fKeepAllClusts[0] = pset.get<bool>          ("KeepAllClustersInd", false);
    fKeepAllClusts[1] = pset.get<bool>          ("KeepAllClustersInd", false);
    fKeepAllClusts[2] = pset.get<bool>          ("KeepAllClustersCol", true);
    
    keepAllClusts = true;
    for(auto& config : fKeepAllClusts ) {
      if( config == false ) {
        keepAllClusts = false; 
        break;
      }
    }
  }



  //###########################################################
  // Main reconstruction procedure.
  //
  // This function does EVERYTHING. The resulting collections of 
  // blip::HitClusts and blip::Blips can then be retrieved after
  // this function is run.
  //###########################################################
  void BlipRecoAlg::RunBlipReco( const art::Event& evt ) {
  
    //std::cout<<"\n"
    //<<"=========== BlipRecoAlg =========================\n"
    //<<"Event "<<evt.id().event()<<" / run "<<evt.id().run()<<"\n";
  
    //=======================================
    // Reset things
    //=======================================
    blips.clear();
    hitclust.clear();
    hitinfo.clear();
    pinfo.clear();
    trueblips.clear();
    EvtBadChanCount = 0;
    //map_plane_hitg4ids.clear();
   

  
    //=======================================
    // Get data products for this event
    //========================================
    
    // --- detector properties
    auto const& SCE_provider        = lar::providerFrom<spacecharge::SpaceChargeService>();
    auto const& chanFilt            = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    auto const& detProp             = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    auto const& clockData           = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    //auto const& detProp              = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    //auto const& lifetime_provider   = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
    //auto const& tpcCalib_provider   = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
    
  
    // -- geometry
    art::ServiceHandle<geo::Geometry> geom;

    // -- G4 particles
    art::Handle< std::vector<simb::MCParticle> > pHandle;
    std::vector<art::Ptr<simb::MCParticle> > plist;
    if (evt.getByLabel(fGeantProducer,pHandle))
      art::fill_ptr_vector(plist, pHandle);
 
    // -- SimEnergyDeposits
    art::Handle<std::vector<sim::SimEnergyDeposit> > sedHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit> > sedlist;
    if (evt.getByLabel(fSimDepProducer,sedHandle)) 
      art::fill_ptr_vector(sedlist, sedHandle);
    
    // -- SimChannels (usually dropped in reco)
    art::Handle<std::vector<sim::SimChannel> > simchanHandle;
    std::vector<art::Ptr<sim::SimChannel> > simchanlist;
    if (evt.getByLabel(fSimChanProducer,simchanHandle)) 
      art::fill_ptr_vector(simchanlist, simchanHandle);

    // -- hits (from input module, usually track-masked subset of gaushit)
    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitProducer,hitHandle))
      art::fill_ptr_vector(hitlist, hitHandle);

    // -- hits (from gaushit), these are used in truth-matching of hits
    art::Handle< std::vector<recob::Hit> > hitHandleGH;
    std::vector<art::Ptr<recob::Hit> > hitlistGH;
    if (evt.getByLabel("gaushit",hitHandleGH))
      art::fill_ptr_vector(hitlistGH, hitHandleGH);

    // -- tracks
    art::Handle< std::vector<recob::Track> > tracklistHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrkProducer,tracklistHandle))
      art::fill_ptr_vector(tracklist, tracklistHandle);
  
    // -- associations
    art::FindManyP<recob::Track> fmtrk(hitHandle,evt,fTrkProducer);
    art::FindManyP<recob::Track> fmtrkGH(hitHandleGH,evt,fTrkProducer);
    art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> fmhh(hitHandleGH,evt,"gaushitTruthMatch");
  
    /*
    //====================================================
    // Update map of bad channels for this event
    //====================================================
    if( fVetoBadChannels ) {
      fBadChanMaskPerEvt = fBadChanMask;
      if( fBadChanProducer != "" ) { 
        std::vector<int> badChans;
        art::Handle< std::vector<int>> badChanHandle;
        if( evt.getByLabel(fBadChanProducer, badChanHandle))
          badChans = *(badChanHandle);
        for(auto& ch : badChans ) {
          EvtBadChanCount++;
          fBadChanMaskPerEvt[ch] = true;
          h_chan_bad->Fill(ch);
        }
      }
    }
    */
   
    //====================================================
    // Prep the particle inventory service for MC+overlay
    //====================================================
    if( evt.isRealData() && plist.size() ) {
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      pi_serv->Rebuild(evt);
      pi_serv->provider()->PrepParticleList(evt);
    }
   
    //===============================================================
    // Map of each hit to its gaushit index (needed if the provided
    // hit collection is some filtered subset of gaushit, in order to
    // use gaushitTruthMatch later on)
    //===============================================================
    std::map< int, int > map_gh;
    // if input collection is already gaushit, this is trivial
    if( fHitProducer == "gaushit" ) {
      for(auto& h : hitlist ) map_gh[h.key()] = h.key(); 
    // ... but if not, find the matching gaushit. There's no convenient
    // hit ID, so we must loop through and compare channel/time (ugh)
    } else {
      std::map<int,std::vector<int>> map_chan_ghid;
      for(auto& gh : hitlistGH ) map_chan_ghid[gh->Channel()].push_back(gh.key());
      for(auto& h : hitlist ) {
        for(auto& igh : map_chan_ghid[h->Channel()]){
          if( hitlistGH[igh]->PeakTime() != h->PeakTime() ) continue;
          map_gh[h.key()] = igh;
          break;
        }
      }
    }
   
    //=====================================================
    // Record PDG for every G4 Track ID
    //=====================================================
    std::map<int,int> map_g4trkid_pdg;
    for(size_t i = 0; i<plist.size(); i++) map_g4trkid_pdg[plist[i]->TrackId()] = plist[i]->PdgCode();
  
    std::map<int, std::map<int,double> > map_g4trkid_chan_energy;
    std::map<int, std::map<int,double> > map_g4trkid_chan_charge;

    //======================================================
    // Use SimChannels to make a map of the collected charge
    // for every G4 particle, instead of relying on the TDC-tick
    // matching that's done by BackTracker's other functions
    //======================================================
    std::map<int,double> map_g4trkid_charge;
    for(auto const &chan : simchanlist ) {
      if( fGeom.View(chan->Channel()) != geo::kW ) continue;
      //std::map<int,double> map_g4trkid_perWireEnergyDep;
      for(auto const& tdcide : chan->TDCIDEMap() ) {
        for(auto const& ide : tdcide.second) {
          if( ide.trackID < 0 ) continue;
          double ne = ide.numElectrons;
          
          // ####################################################
          // # behavior in MicroBooNE as of Nov 2022            #
          // ####################################################
          // WireCell's detsim implements its gain "fudge factor" 
          // by scaling the SimChannel electrons (DocDB 31089)
          // instead of the electronics gain. So we need to correct 
          // for this effect to get accurate count of 'true' 
          // electrons collected on this channel.
          // ####################################################
          if( fSimGainFactor > 0 ) ne /= fSimGainFactor;
          map_g4trkid_charge[ide.trackID] += ne;
         
          // keep track of charge deposited per wire for efficiency plots
          // (coll plane only)
          if( chan->Channel() > 4800 ) {
            map_g4trkid_chan_charge[ide.trackID][chan->Channel()] += ne; 
            if( abs(map_g4trkid_pdg[ide.trackID]) == 11 ) 
              map_g4trkid_chan_energy[ide.trackID][chan->Channel()] += ide.energy;
          }
        
        }
      }
    
    }

    for(auto& m : map_g4trkid_chan_energy ) {
      for(auto& mm : m.second ) {
        if( mm.second > 0 ) h_recoWireEff_denom->Fill(mm.second);
      }
    }
    
    for(auto& m : map_g4trkid_chan_charge ) {
      for(auto& mm : m.second ) {
        if( mm.second > 0 ) h_recoWireEffQ_denom->Fill(mm.second);
      }
    }
   

    //==================================================
    // Use G4 information to determine the "true" blips in this event.
    //==================================================
    if( plist.size() ) {
      pinfo.resize(plist.size());
      for(size_t i = 0; i<plist.size(); i++){
        BlipUtils::FillParticleInfo( *plist[i], pinfo[i], sedlist, fCaloPlane);
        if( map_g4trkid_charge[pinfo[i].trackId] ) pinfo[i].numElectrons = (int)map_g4trkid_charge[pinfo[i].trackId];
        pinfo[i].index = i;
      }
      BlipUtils::MakeTrueBlips(pinfo, trueblips);
      BlipUtils::MergeTrueBlips(trueblips, fTrueBlipMergeDist);
    }


    //=======================================
    // Map track IDs to the index in the vector
    //=======================================
    //std::cout<<"Looping over tracks...\n";
    std::map<size_t,size_t> map_trkid_index;
    for(size_t i=0; i<tracklist.size(); i++) 
      map_trkid_index[tracklist.at(i)->ID()] = i;

    //=======================================
    // Fill vector of hit info
    //========================================
    hitinfo.resize(hitlist.size());
    std::map<int, std::map<int, std::map<int,std::vector<int> >>> cryo_tpc_plane_hitsMap;
    int nhits_untracked = 0;

    for(size_t i=0; i<hitlist.size(); i++){
      auto const& thisHit = hitlist[i];
      auto const& wireid  = thisHit->WireID();
      int   chan    = thisHit->Channel();
      int   cstat   = wireid.Cryostat;
      int   tpc     = wireid.TPC;
      int   plane   = wireid.Plane;
      int   wire    = wireid.Wire;
      
      /*
      const geo::TPCGeo& tpcgeom = fGeom.Cryostat(geo::CryostatID(cstat)).TPC(tpc);
      std::cout<<"Hit in cryo/TPC "<<cstat<<"/"<<tpc<<", drift direction "<<tpcgeom.DriftDirection()<<"\n";

      auto center = tpcgeom.GetCenter();
      auto planecenter = tpcgeom.Plane(0).GetBoxCenter();
      std::cout<<"Center of TPC: "<<center.X()<<","<<center.Y()<<","<<center.Z()<<"\n";
      std::cout<<"Center of planes: "<<planecenter.X()<<","<<planecenter.Y()<<","<<planecenter.Z()<<"\n";
      
      auto const driftVec = planecenter - center;
      std::cout<<"DriftVec: "<<driftVec.X()<<","<<driftVec.Y()<<","<<driftVec.Z()<<"\n";

      short int drift   = tpcgeom.DriftDirection();
      std::cout<<"Drift direction: "<<drift<<"\n";
      */
        
      hitinfo[i].hitid        = i;
      hitinfo[i].cryo         = cstat;
      hitinfo[i].tpc          = tpc;
      hitinfo[i].plane        = plane;
      hitinfo[i].chan         = chan;
      hitinfo[i].wire         = wire;
      hitinfo[i].amp          = thisHit->PeakAmplitude();
      hitinfo[i].rms          = thisHit->RMS();
      hitinfo[i].integralADC  = thisHit->Integral();
      hitinfo[i].sigmaintegral = thisHit->SigmaIntegral();
      hitinfo[i].sumADC       = thisHit->SummedADC();
      hitinfo[i].charge       = fCaloAlg->ElectronsFromADCArea(thisHit->Integral(),plane);
      hitinfo[i].gof          = thisHit->GoodnessOfFit() / thisHit->DegreesOfFreedom();
      hitinfo[i].peakTime     = thisHit->PeakTime();
      hitinfo[i].driftTime    = thisHit->PeakTime()-kXTicksOffsets[cstat][tpc][plane]; //detProp.GetXTicksOffset(wireid);

      //h_hit_times->Fill(thisHit->PeakTime());
      //h_hit_chanstatus->Fill( chanFilt.Status(chan) );

      if( plist.size() ) {
        
        //int truthid;
        //float truthidfrac, numElectrons, energy;
        //BlipUtils::HitTruth( thisHit, truthid, truthidfrac, energy, numElectrons);

        //--------------------------------------------------
        // since SimChannels aren't saved by default, the normal 
        // backtracker won't work, so instead the truth-matching metadata 
        // is stored in the event in the form of the "gaushitTruthMatch"
        // data association.
        //--------------------------------------------------
        int igh = map_gh[i];
        if( fmhh.at(igh).size() ) {
          std::vector<simb::MCParticle const*> pvec;
          std::vector<anab::BackTrackerHitMatchingData const*> btvec;
          fmhh.get(igh,pvec,btvec);
          hitinfo[i].g4energy = 0;
          hitinfo[i].g4charge = 0;
          float maxQ = -9;
          for(size_t j=0; j<pvec.size(); j++){
            hitinfo[i].g4energy += btvec.at(j)->energy;
            hitinfo[i].g4charge += btvec.at(j)->numElectrons;
            if( btvec.at(j)->numElectrons <= maxQ ) continue;
            maxQ = btvec.at(j)->numElectrons;
            hitinfo[i].g4trkid  = pvec.at(j)->TrackId();
            hitinfo[i].g4pdg    = pvec.at(j)->PdgCode();
            hitinfo[i].g4frac   = btvec.at(j)->ideNFraction;
          }
          
          // ###      uB behavior as of Nov 2022              ###
          // WireCell's detsim implements its gain "fudge factor" 
          // by scaling the SimChannel electrons for some reason.
          // So we need to correct for this effect to get accurate
          // count of 'true' electrons collected on channel.
          if( fSimGainFactor > 0 ) hitinfo[i].g4charge /= fSimGainFactor;
         
          if( map_g4trkid_chan_energy[hitinfo[i].g4trkid][chan] > 0 ) {
            double trueEnergyDep = map_g4trkid_chan_energy[hitinfo[i].g4trkid][chan];
            h_recoWireEff_num->Fill(trueEnergyDep);
          }
          if( map_g4trkid_chan_charge[hitinfo[i].g4trkid][chan] > 0 ) {
            double trueChargeDep = map_g4trkid_chan_charge[hitinfo[i].g4trkid][chan];
            h_recoWireEffQ_num->Fill(trueChargeDep);
          }

        }


      }//endif MC
      

      // find associated track
      if( fHitProducer == "gaushit" && fmtrk.isValid() ) {
        if(fmtrk.at(i).size()) hitinfo[i].trkid = fmtrk.at(i)[0]->ID();
      
      // if the hit collection didn't have associations made
      // to the tracks, try gaushit instead
      } else if ( fmtrkGH.isValid() && map_gh.size() ) {
        int gi = map_gh[i];
        if (fmtrkGH.at(gi).size()) hitinfo[i].trkid= fmtrkGH.at(gi)[0]->ID(); 
      }

      // add to the map
      //planehitsMap[plane].push_back(i);
      cryo_tpc_plane_hitsMap[cstat][tpc][plane].push_back(i);
      //tpc_plane_hitsMap[tpc][plane].push_back(i);
      if( hitinfo[i].trkid < 0 ) nhits_untracked++;
      //printf("  %lu   plane: %i,  wire: %i, time: %i\n",i,hitinfo[i].plane,hitinfo[i].wire,int(hitinfo[i].driftTime));

    }//endloop over hits
    
    //for(auto& a : tpc_plane_hitsMap ) {
      //for(auto& b : a.second ) 
        //std::cout<<"TPC "<<a.first<<", plane "<<b.first<<": "<<b.second.size()<<" hits\n";
    //}


    //=================================================================
    // Blip Reconstruction
    //================================================================
    //  
    //  Procedure
    //  [x] Look for hits that were not included in a track 
    //  [x] Filter hits based on hit width, etc
    //  [x] Merge together closely-spaced hits on same wires and adjacent wires
    //  [x] Plane-to-plane time matching
    //  [x] Wire intersection check to get XYZ
    //  [x] Create "blip" object

    // Create a series of masks that we'll update as we go along
    std::vector<bool> hitIsTracked(hitlist.size(),  false);
    std::vector<bool> hitIsGood(hitlist.size(),     true);
    std::vector<bool> hitIsClustered(hitlist.size(),false);
    
    
    // Basic track inclusion cut: exclude hits that were tracked
    for(size_t i=0; i<hitlist.size(); i++){
      if( hitinfo[i].trkid < 0 ) continue;
      auto it = map_trkid_index.find(hitinfo[i].trkid);
      if( it == map_trkid_index.end() ) continue;
      int trkindex = it->second;
      if( tracklist[trkindex]->Length() > fMaxHitTrkLength ) {
        hitIsTracked[i] = true;
        hitIsGood[i] = false;
      }
    }
        

    // Filter based on hit properties. For hits that are a part of
    // multi-gaussian fits (multiplicity > 1), need to re-think this.
    if( fDoHitFiltering ) {
      for(size_t i=0; i<hitlist.size(); i++){
        if( !hitIsGood[i] ) continue;
        hitIsGood[i] = false;
        auto& hit = hitlist[i];
        int plane = hit->WireID().Plane;
        if( hitinfo[i].gof        <= fMinHitGOF[plane] ) continue;
        if( hitinfo[i].gof        >= fMaxHitGOF[plane] ) continue;
        if( hit->RMS()            <= fMinHitRMS[plane] ) continue;
        if( hit->RMS()            >= fMaxHitRMS[plane] ) continue;
        if( hit->PeakAmplitude()  <= fMinHitAmp[plane] ) continue;
        if( hit->PeakAmplitude()  >= fMaxHitAmp )        continue;
        if( hit->Multiplicity()   >= fMaxHitMult )       continue;
        //float hit_ratio = hit->RMS() / hit->PeakAmplitude();
        //if( hit_ratio             < fMinHitRatio[plane] ) continue;
        //if( hit_ratio             > fMaxHitRatio[plane] ) continue;
        
        // we survived the gauntlet of cuts -- hit is good!
        hitIsGood[i] = true;
      }
    }

    //std::cout<<"Hit clustering\n";
    // ---------------------------------------------------
    // Hit clustering
    // ---------------------------------------------------
    std::map<int,std::map<int,std::vector<int>>> tpc_planeclustsMap;
   
    for(auto const& tpc_plane_hitsMap : cryo_tpc_plane_hitsMap ) {
    
    for(auto const& plane_hitsMap : tpc_plane_hitsMap.second ) {
      //std::cout<<"Looking at TPC "<<plane_hitsMap.first<<", which has hits appearing in "<<plane_hitsMap.second.size()<<" planes\n";

      for(auto const& planehits : plane_hitsMap.second){
        //std::cout<<"Looking at TPC "<<plane_hitsMap.first<<", plane "<<planehits.first<<", which has "<<planehits.second.size()<<" hits\n"; 

        for(auto const& hi : planehits.second ){
          
          //std::cout<<"hit "<<hi<<": good "<<hitIsGood[hi]<<", clustered "<<hitIsClustered[hi]<<"\n";

          // skip hits flagged as bad, or already clustered
          if( !hitIsGood[hi] || hitIsClustered[hi] ) continue; 
          
          // initialize a new cluster with this hit as seed
          std::vector<blip::HitInfo> hitinfoVec;
          std::set<int> hitIDs;
          
          hitinfoVec    .push_back(hitinfo[hi]);
          hitIDs        .insert(hi);
          int startWire = hitinfo[hi].wire;
          int endWire   = hitinfo[hi].wire;
          hitIsClustered[hi] = true;

          // see if we can add other hits to it; continue until 
          // no new hits can be lumped in with this clust
          int hitsAdded;
          do{
            hitsAdded = 0;  
            for(auto const& hj : planehits.second ) {
              
              if( !hitIsGood[hj] || hitIsClustered[hj] ) continue; 

              // skip hits outside overall cluster wire range
              int w1 = hitinfo[hj].wire - fHitClustWireRange;
              int w2 = hitinfo[hj].wire + fHitClustWireRange;
              if( w2 < startWire    || w1 > endWire ) continue;
              
              // check for proximity with every other hit added
              // to this cluster so far
              for(auto const& hii : hitIDs ) {

                if( hitinfo[hii].wire > w2 ) continue;
                if( hitinfo[hii].wire < w1 ) continue;
                
                float t1 = hitinfo[hj].peakTime;
                float t2 = hitinfo[hii].peakTime;
                float rms_sum = (hitinfo[hii].rms + hitinfo[hj].rms);
                if( fabs(t1-t2) > fHitClustWidthFact * rms_sum ) continue;

                hitinfoVec.push_back(hitinfo[hj]);
                startWire = std::min( hitinfo[hj].wire, startWire );
                endWire   = std::max( hitinfo[hj].wire, endWire );
                hitIDs.insert(hj);
                hitIsClustered[hj] = true;
                hitsAdded++;
                break;
              }
            

            }
          } while ( hitsAdded!=0 );
          
          blip::HitClust hc = BlipUtils::MakeHitClust(hitinfoVec);
          float span = hc.EndTime - hc.StartTime;
          h_clust_nwires->Fill(hc.NWires);
          h_clust_timespan->Fill(span);
          
          // basic cluster checks
          if( span      <= 0                )   continue;
          if( span      > fMaxClusterSpan   )   continue;
          if( hc.NWires > fMaxWiresInCluster )  continue;
          if( hc.Charge < fMinClusterCharge )   continue;
          if( hc.Charge > fMaxClusterCharge )   continue;
          
          //std::cout<<"Making a new cluster on plane "<<planehits.first<<"\n";
          //std::cout<<"span "<<span<<" ticks, "<<hc.NWires<<" wires, "<<hc.Charge<<" electrons\n";
         
          // Exclude cluster if it is *entirely* on bad channels
          if( fVetoBadChannels ) {
            int nbadchanhits = 0;
            for(auto const& hitID : hc.HitIDs ) {
              int chan = hitinfo[hitID].chan;
              if( chanFilt.Status(chan) < 4 ||
                fBadChanMaskPerEvt[chan] ) nbadchanhits++;
            }
            if( nbadchanhits == hc.NHits ) continue;
          }
          
          // measure wire separation to nearest dead region
          // (0 = directly adjacent)
          for(size_t dw=1; dw<=5; dw++){
            int  w1   = hc.StartWire-dw;
            int  w2   = hc.EndWire+dw;
            bool flag = false;
            // treat edges of wireplane as "dead"
            //if( w1 < 0 || w2 >= (int)fGeom.Nwires(hc.Plane) )
            if( w1 < 0 || w2 >= (int)fGeom.Nwires(geo::PlaneID(0,hc.TPC,hc.Plane)))
              flag=true;
            //otherwise, use channel filter service
            else {
              int ch1 = fGeom.PlaneWireToChannel(geo::WireID(0,hc.TPC,hc.Plane,w1));
              int ch2 = fGeom.PlaneWireToChannel(geo::WireID(0,hc.TPC,hc.Plane,w2));
              if( chanFilt.Status(ch1)<2 ) flag=true;
              if( chanFilt.Status(ch2)<2 ) flag=true;
            }
            if( flag ) { hc.DeadWireSep = dw-1; break; }
          }
          //std::cout<<"DeadWireSep "<<hc.DeadWireSep<<"\n";
         
          // veto this cluster if the gap between it and the
          // nearest dead wire (calculated above) isn't big enough
          if( fMinDeadWireGap > 0 && hc.DeadWireSep < fMinDeadWireGap ) continue;
          
          // **************************************
          // assign the ID, then go back and encode this 
          // cluster ID into the hit information
          // **************************************
          int idx = (int)hitclust.size();
          hc.ID = idx;
          tpc_planeclustsMap[hc.TPC][hc.Plane].push_back(idx);
          for(auto const& hitID : hc.HitIDs) hitinfo[hitID].clustid = hc.ID;
          // ... and find the associated truth-blip
          if( hc.G4IDs.size() ) {
            for(size_t j=0; j< trueblips.size(); j++){
              //std::cout<<"Looking at true blip "<<j<<" which has LeadG4ID of "<<trueblips[j].LeadG4ID<<"\n";
              if( hc.G4IDs.count(trueblips[j].LeadG4ID)) {
                hc.isTruthMatched = true;
                hc.EdepID         = trueblips[j].ID;
                break;
              }
            }
          }
          
          // finally, add the finished cluster to the stack
          hitclust.push_back(hc);

        }
      }//loop over planes
    }//loop over TPCs
    }//loop over cryostats
    //std::cout<<"All done with clustering\n";
    


    // =============================================================================
    // Plane matching and 3D blip formation
    // =============================================================================

    // --------------------------------------
    // Method 1A: Require match between calo plane ( typically collection) and
    //            1 or 2 induction planes. For every hitclust on the calo plane,
    //            do the following:
    //              1. Loop over hitclusts in one of the other planes (same TPC)
    //              3. Find closest-matched clust and add it to the histclust group
    //              4. Repeat for remaining plane(s)
    
    float _matchQDiffLimit= (fMatchQDiffLimit <= 0 ) ? std::numeric_limits<float>::max() : fMatchQDiffLimit;
    float _matchMaxQRatio = (fMatchMaxQRatio  <= 0 ) ? std::numeric_limits<float>::max() : fMatchMaxQRatio;
     
    for(auto& tpcMap : tpc_planeclustsMap ) { // loop on TPCs
      
      //std::cout
      //<<"Performing cluster matching in TPC "<<tpcMap.first<<", which has clusters in "<<tpcMap.second.size()<<" planes\n";
      auto& planeMap = tpcMap.second;
      if( planeMap.find(fCaloPlane) != planeMap.end() ){
        int   planeA              = fCaloPlane;
        auto&  hitclusts_planeA   = planeMap[planeA];
        //std::cout<<"using plane "<<fCaloPlane<<" as reference/calo plane ("<<planeMap[planeA].size()<<" clusts)\n";
        for(auto& i : hitclusts_planeA ) {
          auto& hcA = hitclust[i];
          
          // initiate hit-cluster group
          std::vector<blip::HitClust> hcGroup;
          hcGroup.push_back(hcA);

          // for each of the other planes, make a map of potential matches
          std::map<int, std::set<int>> cands;
         
          // map of cluster ID <--> match metrics
          std::map<int, float> map_clust_dtfrac;
          std::map<int, float> map_clust_dt;
          std::map<int, float> map_clust_overlap;
          std::map<int, float> map_clust_score;

          // ---------------------------------------------------
          // loop over other planes
          for(auto&  hitclusts_planeB : planeMap ) {
            int planeB = hitclusts_planeB.first;
            if( planeB == planeA ) continue;

            // Loop over all non-matched clusts on this plane
            for(auto const& j : hitclusts_planeB.second ) {
              auto& hcB = hitclust[j];
              if( hcB.isMatched ) continue;
              
              // *******************************************
              // Check that the two central wires intersect
              // *******************************************
              double y, z;
              int& chanA = hcA.CenterChan;
              int& chanB = hcB.CenterChan;
              if( !art::ServiceHandle<geo::Geometry>()
                ->ChannelsIntersect(chanA,chanB,y,z)) continue;
              // Save intersect location, so we don't have to
              // make another call to the Geometry service later
              TVector3 xloc(0,y,z);
              hcA.IntersectLocations[hcB.ID] = xloc;
              hcB.IntersectLocations[hcA.ID] = xloc;
              

              // ***********************************
              // Calculate the cluster overlap
              // ***********************************
              float overlapFrac = BlipUtils::CalcHitClustsOverlap(hcA,hcB);
              
              // *******************************************
              // Calculate time difference for start/end, and
              // check that Q-weighted means are comparable
              // *******************************************
              float dt_start  = (hcB.StartTime - hcA.StartTime);
              float dt_end    = (hcB.EndTime   - hcA.EndTime);
              float dt        = ( fabs(dt_start) < fabs(dt_end) ) ? dt_start : dt_end;
              float sigmaT    = std::sqrt(pow(hcA.RMS,2)+pow(hcB.RMS,2));
              float dtfrac    = (hcB.Time - hcA.Time) / sigmaT;
              
              // *******************************************
              // Check relative charge between clusters
              // *******************************************
              float qdiff     = fabs(hcB.Charge-hcA.Charge);
              float ratio     = std::max(hcA.Charge,hcB.Charge)/std::min(hcA.Charge,hcB.Charge);
              
              // *******************************************
              // Combine metrics into a consolidated score
              // *******************************************
              float score = overlapFrac * exp(-fabs(ratio-1.)) * exp(-fabs(dt)/float(fMatchMaxTicks));
              
              // If both clusters are matched to the same MC truth particle,
              // set flag to fill special diagnostic histograms...
              bool trueFlag = (hcA.isTruthMatched && (hcA.EdepID == hcB.EdepID)) ? true : false;
             
              // Diagnostic histograms
              h_clust_overlap[planeB] ->Fill(overlapFrac);
              h_clust_dt[planeB]      ->Fill(dt);
              h_clust_dtfrac[planeB]  ->Fill(dtfrac);
              h_clust_q[planeB]       ->Fill(0.001*hcA.Charge,0.001*hcB.Charge);
              if( score > 0 ) h_clust_score[planeB]   ->Fill(score);
              if( trueFlag ) {
                h_clust_truematch_overlap[planeB]->Fill(overlapFrac);
                h_clust_truematch_dt[planeB]     ->Fill(dt);
                h_clust_truematch_dtfrac[planeB] ->Fill(dtfrac);
                h_clust_truematch_q[planeB]      ->Fill(0.001*hcA.Charge,0.001*hcB.Charge);
                if( score > 0 ) h_clust_truematch_score[planeB]  ->Fill(score);
              }

              if( overlapFrac   < fMatchMinOverlap  ) continue;
              if( fabs(dt)      > fMatchMaxTicks    ) continue;
              if( fabs(dtfrac)  > fMatchSigmaFact   ) continue;
              if( qdiff         > _matchQDiffLimit 
               && ratio         > _matchMaxQRatio   ) continue;
              
              h_clust_q_cut[planeB]->Fill(0.001*hcA.Charge,0.001*hcB.Charge);
              
              // **************************************************
              // We made it through the cuts -- the match is good!
              // **************************************************
              map_clust_dt[j]       = dt;
              map_clust_dtfrac[j]   = dtfrac;
              map_clust_overlap[j]  = overlapFrac;
              map_clust_score[j]    = score;
              cands[planeB]         .insert(j);
            
            }
              
          }//endloop over other planes
          
          // ---------------------------------------------------
          // loop over the candidates found on each plane
          // and select the one with the largest score
          if( cands.size() ) {
            for(auto& c : cands ) {
              int plane = c.first;
              h_nmatches[plane]->Fill(c.second.size());
              float bestScore   = -9;
              int   bestID      = -9;
              for(auto cid : c.second) {
                if( map_clust_score[cid] > bestScore ) {
                  bestScore = map_clust_score[cid];
                  bestID = cid;
                }
              }
              if( bestID >= 0 ) hcGroup.push_back(hitclust[bestID]);
            }
            
            // ----------------------------------------
            // make our new blip, but if it isn't valid, forget it and move on
            blip::Blip newBlip = BlipUtils::MakeBlip(hcGroup,detProp,clockData);
            if( !newBlip.isValid ) continue;
            if( newBlip.NPlanes < fMinMatchedPlanes ) continue;
            
            // ---------------------------------------
            // does this qualify as a "picky" blip?
            bool picky = ( newBlip.NPlanes >= 3 && newBlip.SigmaYZ < 1. );

            // ----------------------------------------
            // save matching information
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].isMatched = true;
              for(auto hit : hitclust[hc.ID].HitIDs) hitinfo[hit].ismatch = true;
             
              // Diagnostic plots for successful 3-plane matches
              if( picky && hc.Plane != fCaloPlane ) {
                float q1 = (float)newBlip.clusters[fCaloPlane].Charge;
                float q2 = (float)newBlip.clusters[hc.Plane].Charge;
                h_clust_picky_overlap[hc.Plane]->Fill(map_clust_overlap[hc.ID]);
                h_clust_picky_dtfrac[hc.Plane] ->Fill(map_clust_dtfrac[hc.ID]);
                h_clust_picky_dt[hc.Plane]     ->Fill(map_clust_dt[hc.ID]);
                h_clust_picky_q[hc.Plane]  ->Fill(0.001*q1,0.001*q2);
              }//end diagnostic plots
              
            }//end loop over clusters

            if( fPickyBlips && !picky ) continue;

            
            // ----------------------------------------
            // apply cylinder cut 
            for(auto& trk : tracklist ){
              if( trk->Length() < fMaxHitTrkLength ) continue;
              auto& a = trk->Vertex();
              auto& b = trk->End();
              TVector3 p1(a.X(), a.Y(), a.Z() );
              TVector3 p2(b.X(), b.Y(), b.Z() );
              // TO-DO: if this track starts or ends at a TPC boundary, 
              // we should extend p1 or p2 to outside the AV to avoid blind spots
              TVector3 bp = newBlip.Position;
              float d = BlipUtils::DistToLine(p1,p2,bp);
              if( d > 0 ) {
                // update closest trkdist
                if( newBlip.ProxTrkDist < 0 || d < newBlip.ProxTrkDist ) {
                  newBlip.ProxTrkDist = d;
                  newBlip.ProxTrkID = trk->ID();
                }
                // need to do some math to figure out if this is in
                // the 45 degreee "cone" relative to the start/end 
                if( !newBlip.inCylinder && d < fCylinderRadius ) {
                  float angle1 = asin( d / (p1-bp).Mag() ) * 180./3.14159;
                  float angle2 = asin( d / (p2-bp).Mag() ) * 180./3.14159;
                  if( angle1 < 45. && angle2 < 45. ) newBlip.inCylinder = true;
                }
              }
            }//endloop over trks
           
            if( fApplyTrkCylinderCut && newBlip.inCylinder ) continue;
            
            // ----------------------------------------
            // if we made it this far, the blip is good!
            // associate this blip with the hits and clusters within it
            newBlip.ID = blips.size();
            blips.push_back(newBlip);
            for(auto& hc : hcGroup ) {
              hitclust[hc.ID].BlipID = newBlip.ID;
              for( auto& h : hc.HitIDs ) hitinfo[h].blipid = newBlip.ID;
            }

  
          }//endif ncands > 0
        }//endloop over caloplane ("Plane A") clusters
      }//endif calo plane has clusters
    }//endloop over TPCs

    // Re-index the clusters after removing unmatched
    if( !keepAllClusts ) {
      std::vector<blip::HitClust> hitclust_filt;
      for(size_t i=0; i<hitclust.size(); i++){
        auto& hc = hitclust[i];
        int blipID = hc.BlipID;
        if( fKeepAllClusts[hc.Plane] || blipID >= 0 ) {
          int idx = (int)hitclust_filt.size();
          hc.ID = idx;
          for( auto& h : hc.HitIDs ) hitinfo[h].clustid = hc.ID;
          if( blipID >= 0 ) blips[blipID].clusters[hc.Plane] = hc;
          hitclust_filt.push_back(hc);
        }
      }
      hitclust = hitclust_filt;
    }
    
    //std::cout<<"Found "<<hitclust.size()<<" clusters and "<<blips.size()<<" blips\n";
    

    for(size_t i=0; i<hitlist.size(); i++){
      if (hitinfo[i].trkid >= 0 ) continue;
      h_chan_nhits->Fill(fGeom.PlaneWireToChannel(geo::WireID(0,hitinfo[i].tpc,hitinfo[i].plane,hitinfo[i].wire)));
      
      int clustid = hitinfo[i].clustid;
      if( clustid >= 0 ) {
        if( hitclust[clustid].NWires > 1 ) continue;
        h_chan_nclusts->Fill(fGeom.PlaneWireToChannel(geo::WireID(0,hitinfo[i].tpc,hitinfo[i].plane,hitinfo[i].wire)));
      }
      //if( hitinfo[i].ismatch    ) continue;
      //if( hitclust[clustid].NWires > 1 ) continue;
      //h_chan_nclusts->Fill(fGeom.PlaneWireToChannel(hitinfo[i].plane,hitinfo[i].wire));
    }

    
    //*************************************************************************
    // Loop over the vector of blips and perform calorimetry calculations
    //*************************************************************************
    for(size_t i=0; i<blips.size(); i++){
      auto& blip = blips[i];
     
      blip.Charge = blip.clusters[fCaloPlane].Charge;
      //std::cout<<"blip "<<i<<": TPC "<<blip.TPC<<", XYZ "<<blip.X()<<","<<blip.Y()<<","<<blip.Z()<<", charge "<<blip.Charge<<"\n";
      
      // ***** MICROBOONE ************
      // --- YZ uniformity correction ---
      // Correct for charge-collection non-uniformity based on Y/Z position
      // (taken from CalibrationdEdx_module)
      //if( fYZUniformityCorr ) blip.Charge *= tpcCalib_provider.YZdqdxCorrection(fCaloPlane,blip.Position.Y(),blip.Position.Z());

      // ================================================================================
      // Calculate blip energy assuming T = T_beam (eventually can do more complex stuff
      // like associating blip with some nearby track/shower and using its tagged T0)
      //    Method 1: Assume a dE/dx for electrons, use that + local E-field to get recomb.
      //    Method 2: ESTAR lookup table method ala ArgoNeuT
      // ================================================================================
      float depEl   = std::max(0.0,(double)blip.Charge);
      float Efield  = kNominalEfield;

      // --- Lifetime correction ---
      // Ddisabled by default. Without knowing real T0 of a blip, attempting to 
      // apply this correction can do more harm than good! Note lifetime is in
      // units of 'ms', not microseconds, hence the 1E-3 conversion factor.
      if( fLifetimeCorr && blip.Time>0 ) depEl *= exp( 1e-3*blip.Time/detProp.ElectronLifetime());
      
      // --- SCE corrections ---
      geo::Point_t point( blip.Position.X(),blip.Position.Y(),blip.Position.Z() );
      if( fSCECorr ) {

        // 1) Spatial correction
        if( SCE_provider->EnableCalSpatialSCE() ) {
          // TODO: Deal with cases where X falls outside AV (diffuse out-of-time signal)
          //       For example, maybe re-assign to center of drift volume?
          geo::Vector_t loc_offset = SCE_provider->GetCalPosOffsets(point,blip.TPC);
          point.SetXYZ(point.X()-loc_offset.X(),point.Y()+loc_offset.Y(),point.Z()+loc_offset.Z());
        }
      
        // 2) E-field correction
        //
        // notes:
        //   - GetEfieldOffsets(xyz) and GetCalEfieldOffsets(xyz) return the exact
        //     same underlying E-field offset map; the only difference is the former
        //     is used in the simulation, and the latter in reconstruction (??).
        //   - The SpaceCharge service must have 'EnableCorSCE' and 'EnableCalEfieldSCE'
        //     enabled in order to use GetCalEfieldOffsets
        //   - Blips may appear to be outside the active volume if T0-corrections aren't 
        //     applied to the reconstructed 'X'. SCE map should return (0,0,0) in this case.
        if( SCE_provider->EnableCalEfieldSCE() ) {
          auto const field_offset = SCE_provider->GetCalEfieldOffsets(point,blip.TPC); 
          Efield = detProp.Efield()*std::hypot(1+field_offset.X(),field_offset.Y(),field_offset.Z());;
        }

      }
      
      // METHOD 1
      float recomb  = ModBoxRecomb(fCalodEdx,Efield);
      blip.Energy   = depEl * (1./recomb) * kWion;
      h_recomb      ->Fill(recomb);
      
      // METHOD 2 (TODO)
      //std::cout<<"Calculating ESTAR energy dep...  "<<depEl<<", "<<Efield<<"\n";
      //blips[i].EnergyESTAR = ESTAR->Interpolate(depEl, Efield); 
      
      // ================================================
      // Save the true blip into the object;
      // at least one cluster needs to match.
      // ================================================
      std::set<int> set_edepids;
      for(auto& hc : blip.clusters ) {
        if( !hc.isValid )   continue; 
        if( hc.EdepID < 0 ) continue;
        set_edepids.insert( hc.EdepID );
      }
      if( set_edepids.size() == 1 ) 
        blip.truth = trueblips[*set_edepids.begin()];
    
    }//endloop over blip vector

  }//End main blip reco function
 
  
  
  //###########################################################
  float BlipRecoAlg::ModBoxRecomb(float dEdx, float Efield) {
    float rho = kLArDensity;
    float Xi = fModBoxB * dEdx / ( Efield * rho );
    return log(fModBoxA+Xi)/Xi;
  }

  float BlipRecoAlg::dQdx_to_dEdx(float dQdx_e, float Efield){
    float rho = kLArDensity;
    float beta  = fModBoxB / (rho * Efield);
    float alpha = fModBoxA;
    return ( exp( beta * kWion * dQdx_e ) - alpha ) / beta;
  }
  
  float BlipRecoAlg::Q_to_E(float Q, float Efield){
    if( Efield != kNominalEfield )  return kWion * (Q/ModBoxRecomb(fCalodEdx,Efield));
    else                            return kWion * (Q/kNominalRecombFactor);
  }

  
  //###########################################################
  void BlipRecoAlg::PrintConfig() {
  
    printf("BlipRecoAlg Configurations\n\n");
    printf("  Input hit collection      : %s\n",          fHitProducer.c_str());
    printf("  Input trk collection      : %s\n",          fTrkProducer.c_str());
    printf("  Max wires per cluster     : %i\n",          fMaxWiresInCluster);
    printf("  Max cluster timespan      : %.1f ticks\n",    fMaxClusterSpan);
    printf("  Min cluster overlap       : %4.1f\n",       fMatchMinOverlap);
    printf("  Clust match sigma-factor  : %4.1f\n",       fMatchSigmaFact);
    printf("  Clust match max dT        : %4.1f ticks\n", fMatchMaxTicks);
    printf("  Charge diff limit         : %.1fe3\n",fMatchQDiffLimit/1000);
    printf("  Charge ratio maximum      : %.1f\n",fMatchMaxQRatio);      
    
    /*
    printf("  Min cluster overlap       : ");
    for(auto val : fMatchMinOverlap) { printf("%3.1f   ",val); } printf("\n");
    printf("  Clust match sigma factor  : ");
    for(auto val : fMatchSigmaFact) { printf("%3.1f   ",val); } printf("\n");
    printf("  Clust match max ticks     : ");
    for(auto val : fMatchMaxTicks) { printf("%3.1f   ",val); } printf("\n");
    */

    for(int i=0; i<kNplanes; i++){
    if( i == fCaloPlane ) continue;
    printf("  pl%i matches per cand      : %4.2f\n",       i,h_nmatches[i]->GetMean());
    }
    
    //printf("  Track-cylinder radius     : %.1f cm\n",       fCylinderRadius);
    //printf("  Applying cylinder cut?    : %i\n",          fApplyTrkCylinderCut);
    //printf("  Picky blip mode?          : %i\n",        fPickyBlips);
    printf("\n");
    
  }
  
}
