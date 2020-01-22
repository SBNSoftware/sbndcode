#include "FlashMatchAlg.h"

namespace sbnd{

FlashMatchAlg::FlashMatchAlg(const Config& config){

  this->reconfigure(config);

  fGeometryService = lar::providerFrom<geo::Geometry>();

}


FlashMatchAlg::FlashMatchAlg(){

  fGeometryService = lar::providerFrom<geo::Geometry>();

}


FlashMatchAlg::~FlashMatchAlg(){

}


void FlashMatchAlg::reconfigure(const Config& config){

  fSpacePointLabel = config.SpacePointLabel();
  fPandoraLabel = config.PandoraLabel();
  fInputFile = config.InputFile();
  fBeamFlashMin = config.BeamFlashMin();
  fBeamFlashMax = config.BeamFlashMax();
  fLightWindowMin = config.LightWindowMin();
  fLightWindowMax = config.LightWindowMax();
  fEventMin = config.EventMin();
  fEventMax = config.EventMax();
  fTimeResolution = config.TimeResolution();
  fFlashThreshold = config.FlashThreshold();
  fFlashTimeCorrection = config.FlashTimeCorrection();
  fFlashMatchCut = config.FlashMatchCut();

  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fInputFile, fname);
  TFile *infile = new TFile(fname.c_str(), "READ");

  TH1 *temphisto = (TH1*)infile->Get("rrp1");
  rr_nbins = temphisto->GetNbinsX();
  if (rr_nbins<=0) {
    rr_nbins=1;
    rrmean.push_back(0);
    rrsp.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= rr_nbins; ++ib) {
      rrmean.push_back(temphisto->GetBinContent(ib));
      float tt = temphisto->GetBinError(ib);
      if(tt<=0) tt = 100.;
      rrsp.push_back(tt);
    }
  }

  temphisto = (TH1*)infile->Get("dyp1");
  dy_nbins = temphisto->GetNbinsX();
  if (dy_nbins<=0) {
    dy_nbins=1;
    dymean.push_back(0);
    dysp.push_back(0.001);
  }
  else {
    for (int ib = 1;ib <= dy_nbins; ++ib) {
      dymean.push_back(temphisto->GetBinContent(ib));
      float tt = temphisto->GetBinError(ib);
      if(tt<=0) tt = 100.;
      dysp.push_back(tt);
    }
  }
  
  temphisto = (TH1*)infile->Get("dzp1");
  dz_nbins = temphisto->GetNbinsX();
  if (dz_nbins<=0) {
    dz_nbins=1;
    dzmean.push_back(0);
    dzsp.push_back(0.001);
  }
  else {
    for (int ib = 1; ib <= dz_nbins; ++ib) {
      dzmean.push_back(temphisto->GetBinContent(ib));
      float tt = temphisto->GetBinError(ib);
      if(tt<=0) tt = 100.;
      dzsp.push_back(tt);
    }
  }

  temphisto = (TH1*)infile->Get("pep1");
  pe_nbins = temphisto->GetNbinsX();
  if (pe_nbins<=0) {
    pe_nbins=1;
    pemean.push_back(0);
    pesp.push_back(0.001);
  }
  else {
    for (int ib=1; ib <= pe_nbins; ++ib) {
      pemean.push_back(temphisto->GetBinContent(ib));
      float tt = temphisto->GetBinError(ib);
      if(tt<=0) tt = 100.;
      pesp.push_back(tt);
    }
  }
  
  infile->Close();

  return;
}


// Find the biggest optical flash in the beam window
std::pair<double, double> FlashMatchAlg::BiggestBeamFlash(std::vector<recob::OpHit> ophits){

    int nbins = 1./fTimeResolution * (fBeamFlashMax - fBeamFlashMin);
    TH1F *ophittime = new TH1F("ophittime", "ophittime", nbins, fBeamFlashMin, fBeamFlashMax); // in us
    for(size_t i = 0; i < ophits.size(); i++){

      if ( !fChannelMap.pdType(ophits[i].OpChannel(),"pmt")) continue;

      if ( (ophits[i].PeakTime() > fBeamFlashMin) && (ophits[i].PeakTime() < fBeamFlashMax) ) {
	      ophittime->Fill(ophits[i].PeakTime(), ophits[i].PE());
      }

    }

    auto ibin =  ophittime->GetMaximumBin();
    float flashtime = (ibin * fTimeResolution) + fBeamFlashMin;  // in us
    float lowedge = flashtime + fLightWindowMin;
    float highedge = flashtime + fLightWindowMax;

    delete ophittime;
    return std::make_pair(lowedge, highedge);
}


// Find all optical flash times in the event
std::vector<double> FlashMatchAlg::OpFlashes(std::vector<recob::OpHit> ophits){

  std::vector<double> opflashes;

  // get flash time
  int nbins = 1./fTimeResolution * (fEventMax - fEventMin);
  TH1F *ophittimes = new TH1F("ophittime", "ophittime", nbins, fEventMin, fEventMax); // in us
  for(size_t i = 0; i < ophits.size(); i++){
    if ( !fChannelMap.pdType(ophits[i].OpChannel(), "pmt")) continue;
	  ophittimes->Fill(ophits[i].PeakTime(), ophits[i].PE());
  }

  // Find all flashes above threshold
  for(int i = 0; i < ophittimes->GetNbinsX(); i++){
    if(ophittimes->GetBinContent(i) > fFlashThreshold){
      double flash_time = ((double)i*fTimeResolution) + fEventMin + fFlashTimeCorrection;
      opflashes.push_back(flash_time);
      // Same flash if hit still above threshold or within 100 ns
      while(ophittimes->GetBinContent(i) > fFlashThreshold 
            || ((double)i*fTimeResolution)+fEventMin+fFlashTimeCorrection-flash_time < fLightWindowMax-fLightWindowMin) i++;
    }
  }

  delete ophittimes;
  return opflashes;
}


// Create real PDS optical flashes per TPC
std::pair<std::vector<double>, std::vector<double>> FlashMatchAlg::OpFlashes(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle){

  // Optical flash reconstruction for numuCC
  std::vector<recob::OpHit> ophits_tpc0;
  std::vector<recob::OpHit> ophits_tpc1;
  for(auto const& ophit : (*pdsHandle)){
    // Only look at PMTs
    if( fChannelMap.pdName(ophit.OpChannel()) != "pmt" ) continue;
    // Work out what TPC detector is in
    double PMTxyz[3];
	  fGeometryService->OpDetGeoFromOpChannel(ophit.OpChannel()).GetCenter(PMTxyz);
    if(PMTxyz[0] > 0) ophits_tpc1.push_back(ophit);
    else ophits_tpc0.push_back(ophit);
  }

  std::vector<double> opflashes_tpc0 = OpFlashes(ophits_tpc0);
  std::vector<double> opflashes_tpc1 = OpFlashes(ophits_tpc1);

  return std::make_pair(opflashes_tpc0, opflashes_tpc1);
}


// Determine if there is a PDS flash in time with the neutrino beam per TPC
std::pair<bool, bool> FlashMatchAlg::BeamFlash(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle){

  std::vector<recob::OpHit> ophits_tpc0;
  std::vector<recob::OpHit> ophits_tpc1;
  for(auto const& ophit : (*pdsHandle)){
    // Only look at PMTs
    if( fChannelMap.pdName(ophit.OpChannel()) != "pmt" ) continue;
    // Work out what TPC detector is in
    double PMTxyz[3];
	  fGeometryService->OpDetGeoFromOpChannel(ophit.OpChannel()).GetCenter(PMTxyz);
    if(ophit.PeakTime() < fBeamFlashMin || ophit.PeakTime() > fBeamFlashMax) continue;
    if(PMTxyz[0] > 0) ophits_tpc1.push_back(ophit);
    else ophits_tpc0.push_back(ophit);
  }

  std::vector<double> opflashes_tpc0 = OpFlashes(ophits_tpc0);
  std::vector<double> opflashes_tpc1 = OpFlashes(ophits_tpc1);

  bool tpc0Flash = opflashes_tpc0.size() > 0;
  bool tpc1Flash = opflashes_tpc1.size() > 0;

  return std::make_pair(tpc0Flash, tpc1Flash);
}

// Return the number of reconstructed PE in time with beam in TPC
std::pair<double, double> FlashMatchAlg::BeamPE(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle){

  double npe_tpc0 = 0;
  double npe_tpc1 = 0;
  for(auto const& ophit : (*pdsHandle)){
    if ( !fChannelMap.pdType(ophit.OpChannel(),"pmt")) continue;
    // Work out what TPC detector is in
    double PMTxyz[3];
	  fGeometryService->OpDetGeoFromOpChannel(ophit.OpChannel()).GetCenter(PMTxyz);
    if(ophit.PeakTime() < fBeamFlashMin || ophit.PeakTime() > fBeamFlashMax) continue;
    if(PMTxyz[0] > 0) npe_tpc1 += ophit.PE();
    else npe_tpc0 += ophit.PE();
  }

  return std::make_pair(npe_tpc0, npe_tpc1);

}

// Calculate the flash matching variables for optical hits in TPC and time window
std::vector<double> FlashMatchAlg::OpVariables(std::vector<recob::OpHit> ophits, int tpc, double start_t, double end_t){

  std::vector<double> variables {0, 0, -99999, -99999};

  double PMTxyz[3];
	double sum_Ay = 0; 
  double sum_Az = 0;
	double sum_Cy = 0; 
  double sum_Cz = 0;
	double sum_D = 0;
  double unpe_tot = 0;
  double pe_tot = 0;

	for(auto const& oph : ophits){

    // Time window of flash
	  if ( !((oph.PeakTime() > start_t) && (oph.PeakTime() < end_t)) ) continue;

    // Check if in correct TPC
	  fGeometryService->OpDetGeoFromOpChannel(oph.OpChannel()).GetCenter(PMTxyz);
	  if ((tpc==0 && PMTxyz[0]>0) || (tpc==1 && PMTxyz[0]<0) ) continue;

    // For TPB coated PMTs
	  if ( fChannelMap.pdType(oph.OpChannel(),"pmt")){
	    // Add up the position, weighting with PEs
      pe_tot  += oph.PE();
      variables[0] += oph.PE() * PMTxyz[1];
      variables[1] += oph.PE() * PMTxyz[2];
	    sum_Ay  += pow(oph.PE(),2.) * pow(PMTxyz[1],2.);
	    sum_Az  += pow(oph.PE(),2.) * pow(PMTxyz[2],2.);
	    sum_D   += pow(oph.PE(),2.);
	    sum_Cy  += pow(oph.PE(),2.) * PMTxyz[1];
	    sum_Cz  += pow(oph.PE(),2.) * PMTxyz[2];
    }
    // For uncoated PMTs
	  else if ( fChannelMap.pdType(oph.OpChannel(),"barepmt")){
      unpe_tot += oph.PE();
      variables[0] += oph.PE() * PMTxyz[1];
      variables[1] += oph.PE() * PMTxyz[2];
    }
	}
	
  // If no hits return null vector
  if(pe_tot <= 0) return variables;

  // Calculate charge weighted center
	variables[0] /= (pe_tot + unpe_tot);
	variables[1] /= (pe_tot + unpe_tot);
  // Calculate PE spread
	variables[2] = sqrt((sum_Ay - pow(sum_Cy,2.)/sum_D + sum_Az - pow(sum_Cz,2.)/sum_D)/sum_D);
  // Calculate PE ratio
  variables[3] = 100. * unpe_tot / pe_tot;

  return variables;

}


// Calculate the flash matching score from charge weighted mean and optical hit variables
double FlashMatchAlg::DyScore(double x, double y, std::vector<double> variables){
    
  //calculate match score here, put association on the event
	float slice = fTpcGeo.MaxX() - abs(x);
  float drift_distance = fTpcGeo.MaxX();

	int y_bin = int(dy_nbins * (slice / drift_distance));
  
  return abs(abs(variables[0] - y) - dymean[y_bin])/dysp[y_bin];
}


// Calculate the flash matching score from charge weighted mean and optical hit variables
double FlashMatchAlg::DzScore(double x, double z, std::vector<double> variables){
    
  //calculate match score here, put association on the event
	float slice = fTpcGeo.MaxX() - abs(x);
  float drift_distance = fTpcGeo.MaxX();

	int z_bin = int(dz_nbins * (slice / drift_distance));
  
  return abs(abs(variables[1] - z) - dzmean[z_bin])/dzsp[z_bin];
}


// Calculate the flash matching score from charge weighted mean and optical hit variables
double FlashMatchAlg::RrScore(double x, std::vector<double> variables){
    
  //calculate match score here, put association on the event
	float slice = fTpcGeo.MaxX() - abs(x);
  float drift_distance = fTpcGeo.MaxX();

	int r_bin = int(rr_nbins * (slice / drift_distance));
  
  return abs(variables[2] - rrmean[r_bin])/rrsp[r_bin];
}


// Calculate the flash matching score from charge weighted mean and optical hit variables
double FlashMatchAlg::PeScore(double x, std::vector<double> variables){
    
  //calculate match score here, put association on the event
	float slice = fTpcGeo.MaxX() - abs(x);
  float drift_distance = fTpcGeo.MaxX();

	int pe_bin = int(pe_nbins * (slice / drift_distance));
  
  return abs(variables[3] - pemean[pe_bin])/pesp[pe_bin];
}


// Calculate the flash matching score from charge weighted mean and optical hit variables
double FlashMatchAlg::FlashScore(double x, double y, double z, std::vector<double> variables, double w1, double w2, double w3, double w4){
    
  //calculate match score here, put association on the event
	float slice = fTpcGeo.MaxX() - abs(x);
  float drift_distance = fTpcGeo.MaxX();

	int y_bin = int(dy_nbins * (slice / drift_distance));
	int z_bin = int(dz_nbins * (slice / drift_distance));
	int r_bin = int(rr_nbins * (slice / drift_distance));
	int pe_bin = int(pe_nbins * (slice / drift_distance));
  
  double score = 0;
  if(dysp[y_bin] > 0) score += w1 * abs(abs(variables[0] - y) - dymean[y_bin])/dysp[y_bin];
  if(dzsp[z_bin] > 0) score += w2 * abs(abs(variables[1] - z) - dzmean[z_bin])/dzsp[z_bin];
  if(rrsp[r_bin] > 0) score += w3 * abs(variables[2] - rrmean[r_bin])/rrsp[r_bin];
  if(pesp[pe_bin] > 0) score += w4 * abs(variables[3] - pemean[pe_bin])/pesp[pe_bin];

	return score;
}



// Return flash score for a PFParticle
double FlashMatchAlg::FlashScore(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event, art::ValidHandle<std::vector<recob::OpHit>> pdsHandle){

  // Associations
  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  event.getByLabel(fPandoraLabel, pfParticleHandle);
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn_v(pfParticleHandle, event, fPandoraLabel);
  auto const& spacepoint_h = event.getValidHandle<std::vector<recob::SpacePoint> >(fSpacePointLabel);
  art::FindManyP<recob::Hit> spacepoint_hit_assn_v(spacepoint_h, event, fSpacePointLabel);

  // Get interesting ophits
  std::vector<recob::OpHit> ophs;
  for(auto const& oph : (*pdsHandle)){
    if ( fChannelMap.pdType(oph.OpChannel(),"pmt") || fChannelMap.pdType(oph.OpChannel(),"barepmt") ) {
      ophs.push_back(oph);
    }
  }

  // Reconstruct optical variables
  std::pair<double, double> beam_flash = BiggestBeamFlash(ophs);
  std::vector<std::vector<double>> opvars;
  for (size_t it=0; it<fGeometryService->NTPC(); ++it) {
    opvars.push_back(OpVariables(ophs, it, beam_flash.first, beam_flash.second));
  }

  // go through these pfparticles and fill info needed for matching
  int pfp_tpc = -1;
  double nuvtx_x = 0.;
  double nuvtx_y = 0.;
  double nuvtx_z = 0.;
  double norm = 0.;

  for (const size_t daughterId : pfparticle.Daughters()) {
    art::Ptr<recob::PFParticle> pDaughter = pfParticleMap.at(daughterId);

    const std::vector< art::Ptr<recob::SpacePoint> >& SPs = pfp_spacepoint_assn_v.at(pDaughter.key());
    for (size_t sp=0; sp < SPs.size(); sp++) {
      const std::vector< art::Ptr<recob::Hit> > hits = spacepoint_hit_assn_v.at( SPs[sp].key() );

      for (size_t h=0; h < hits.size(); h++) {

        // Only use hits from the collection plane
        if (fGeometryService->SignalType(hits[h]->WireID()) != geo::kCollection) continue;

        // Add the charged point to the vector
        const auto &position(SPs[sp]->XYZ());
        const auto charge(hits[h]->Integral());

        int tpcindex = (hits[h]->WireID()).TPC;
        if(pfp_tpc == -1) pfp_tpc = tpcindex;
        if(tpcindex != pfp_tpc) pfp_tpc = -2;

        nuvtx_x += charge * position[0];
        nuvtx_y += charge * position[1];
        nuvtx_z += charge * position[2];
        norm += charge;

      } // for all hits associated to this spacepoint
    } // for all spacepoints
  } // for all pfp pointers

  // No charge deposition in PFP, very strange, remove
  if(norm <= 0) return 99999;
  // PFP in two TPCs, should have t0 tag already, keep
  if(pfp_tpc <= -1) return 0;
  // No flash in TPC, remove
  if(opvars[pfp_tpc][2] == -99999) return 99999;

  nuvtx_x /= norm;
  nuvtx_y /= norm;
  nuvtx_z /= norm;

  return FlashScore(nuvtx_x, nuvtx_y, nuvtx_z, opvars[pfp_tpc], 2., 2., 4., 0.);

  if(FlashScore(nuvtx_x, nuvtx_y, nuvtx_z, opvars[pfp_tpc], 2., 2., 4., 0.) > fFlashMatchCut) return false;
  return true;
  
}


// Determine if PFParticle matches a beam flash
bool FlashMatchAlg::FlashMatch(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event, art::ValidHandle<std::vector<recob::OpHit>> pdsHandle){

  if(FlashScore(pfparticle, pfParticleMap, event, pdsHandle) > fFlashMatchCut) return false;
  return true;

}

}
