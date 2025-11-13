#include "sbndcode/OpDetSim/DigiPMTSBNDAlg.hh"

//------------------------------------------------------------------------------
//--- opdet::simpmtsbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet {


  DigiPMTSBNDAlg::DigiPMTSBNDAlg(ConfigurationParameters_t const& config)
    : fParams(config)
    , fSampling(fParams.frequency)
    , fPMTCoatedVUVEff_tpc0(fParams.PMTCoatedVUVEff_tpc0 / fParams.larProp->ScintPreScale())
    , fPMTCoatedVISEff_tpc0(fParams.PMTCoatedVISEff_tpc0 / fParams.larProp->ScintPreScale())
    , fPMTUncoatedEff_tpc0(fParams.PMTUncoatedEff_tpc0/ fParams.larProp->ScintPreScale())
    , fPMTCoatedVUVEff_tpc1(fParams.PMTCoatedVUVEff_tpc1 / fParams.larProp->ScintPreScale())
    , fPMTCoatedVISEff_tpc1(fParams.PMTCoatedVISEff_tpc1 / fParams.larProp->ScintPreScale())
    , fPMTUncoatedEff_tpc1(fParams.PMTUncoatedEff_tpc1/ fParams.larProp->ScintPreScale())
      //  , fSinglePEmodel(fParams.SinglePEmodel)
    , fUseDataNoise(fParams.UseDataNoise)
    , fOpDetNoiseFile(fParams.OpDetNoiseFile)
    , fEngine(fParams.engine)
    , fFlatGen(*fEngine)
    , fPoissonQGen(*fEngine)
    , fGaussQGen(*fEngine)
    , fExponentialGen(*fEngine)
  {

    mf::LogInfo("DigiPMTSBNDAlg") << "PMT corrected efficiencies TPC0 = "
                                  << fPMTCoatedVUVEff_tpc0 << " " << fPMTCoatedVISEff_tpc0 << " " << fPMTUncoatedEff_tpc0
                                  << "PMT corrected efficiencies TPC1 = "
                                  << fPMTCoatedVUVEff_tpc1 << " " << fPMTCoatedVISEff_tpc1 << " " << fPMTUncoatedEff_tpc1 <<"\n";

    if(fPMTCoatedVUVEff_tpc0 > 1.0001 || fPMTCoatedVISEff_tpc0 > 1.0001 || fPMTUncoatedEff_tpc0 > 1.0001 || fPMTCoatedVUVEff_tpc1 > 1.0001 || fPMTCoatedVISEff_tpc1 > 1.0001 || fPMTUncoatedEff_tpc1 > 1.0001)
      mf::LogWarning("DigiPMTSBNDAlg")
        << "Detection efficiencies set in fhicl file seem to be too large!\n"
        << "PMTCoatedVUVEff TPC0: " << fParams.PMTCoatedVUVEff_tpc0 << "\n"
        << "PMTCoatedVISEff TPC0: " << fParams.PMTCoatedVISEff_tpc0 << "\n"
        << "PMTUncoatedEff TPC0: " << fParams.PMTUncoatedEff_tpc0 << "\n"
        << "PMTCoatedVUVEff TPC1: " << fParams.PMTCoatedVUVEff_tpc1 << "\n"
        << "PMTCoatedVISEff TPC1: " << fParams.PMTCoatedVISEff_tpc1 << "\n"
        << "PMTUncoatedEff TPC1: " << fParams.PMTUncoatedEff_tpc1 << "\n"
        << "Final efficiency must be equal or smaller than the scintillation "
        << "pre scale applied at simulation time.\n"
        << "Please check this number (ScintPreScale): "
        << fParams.larProp->ScintPreScale();

    fSampling = fSampling / 1000.0; //in GHz, to cancel with ns
    fSamplingPeriod = 1./fSampling;

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fParams.PMTDataFile, fname);
    TFile* file = TFile::Open(fname.c_str(), "READ");

    std::vector<std::vector<double>>* SinglePEVec_p;
    std::vector<int>* fSinglePEChannels_p;
    std::vector<double>* fPeakAmplitude_p;

    file->GetObject("SERChannels", fSinglePEChannels_p);
    file->GetObject("SinglePEVec", SinglePEVec_p);
    file->GetObject("PeakAmplitude",  fPeakAmplitude_p);

    fSinglePEWaveVector = *SinglePEVec_p;
    fSinglePEChannels = *fSinglePEChannels_p;
    fPeakAmplitude = *fPeakAmplitude_p;

    // TPB emission time histogram for pmt_coated histogram
    std::vector<double>* timeTPB_p;
    file->GetObject("timeTPB", timeTPB_p);
    fTimeTPB = std::make_unique<CLHEP::RandGeneral>
      (*fEngine, timeTPB_p->data(), timeTPB_p->size());

    // PMT calibration database service
    fPMTCalibrationDatabaseService = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();
    fPMTHDOpticalWaveformsPtr = art::make_tool<opdet::HDOpticalWaveform>(fParams.HDOpticalWaveformParams);

    //Resize the SER vector to the number of channels
    fSinglePEWave_HD.resize(320);
    //shape of single pulse
    if (fParams.PMTSinglePEmodel=="testbench") {
      mf::LogDebug("DigiPMTSBNDAlg") << " using testbench pe response";
      std::vector<double>* SinglePEVec_p;
      file->GetObject("SinglePEVec_HD", SinglePEVec_p);
      fSinglePEWave = *SinglePEVec_p;
      // Prepare HD waveforms
      // Create the same SER for all channels
      for(size_t i=0; i<fSinglePEWave_HD.size(); i++){
        fPMTHDOpticalWaveformsPtr->produceSER_HD(fSinglePEWave_HD[i], fSinglePEWave);
      }
      pulsesize = fSinglePEWave_HD[0][0].size();
      mf::LogDebug("DigiPMTSBNDAlg")<<"HD wvfs size: "<<pulsesize;
    }
    else if (fParams.PMTSinglePEmodel=="ideal"){
      mf::LogDebug("DigiPMTSBNDAlg") << " using ideal pe response";
      //shape of single pulse
      sigma1 = fParams.PMTRiseTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      sigma2 = fParams.PMTFallTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      pulsesize = (int)((6 * sigma2 + fParams.TransitTime) * fSampling);
      fSinglePEWave.resize(pulsesize);
      Pulse1PE(fSinglePEWave);
    }
    else if (fParams.PMTSinglePEmodel=="data"){
      // Read the SER from the calibration database (channel independent)
      for(size_t i=0; i<fSinglePEWave_HD.size(); i++){
        if(!fPMTCalibrationDatabaseService->getReconstructChannel(i)) continue;
        fSinglePEWave = fPMTCalibrationDatabaseService->getSER(i);
        double SPEAmplitude =  fPMTCalibrationDatabaseService->getSPEAmplitude(i);
        double SPEPeakValue = *std::max_element(fSinglePEWave.begin(), fSinglePEWave.end(), [](double a, double b) {return std::abs(a) < std::abs(b);});
        double SinglePENormalization = std::abs(SPEAmplitude/SPEPeakValue);
        std::transform(fSinglePEWave.begin(), fSinglePEWave.end(), fSinglePEWave.begin(), [SinglePENormalization](double val) {return val * SinglePENormalization;});
        if(fSinglePEWave.size()==0) continue;
        if(fSinglePEWave.size()>0) pulsesize = fSinglePEWave.size();
        fPMTHDOpticalWaveformsPtr->produceSER_HD(fSinglePEWave_HD[i], fSinglePEWave);
      }
    }
    else{
      throw cet::exception("DigiPMTSBNDAlg") << "Wrong PMTSinglePEmodel configuration: " << fParams.PMTSinglePEmodel << std::endl;
    }

    if(fParams.MakeGainFluctuations){
      fPMTGainFluctuationsPtr = art::make_tool<opdet::PMTGainFluctuations>(fParams.GainFluctuationsParams);
    }
    if(fParams.SimulateNonLinearity){
      fPMTNonLinearityPtr = art::make_tool<opdet::PMTNonLinearity>(fParams.NonLinearityParams); 
      fPMTNonLinearityPtr->ConfigureNonLinearity();
    }

    // infer pulse polarity from SER peak sign
    double minADC_SinglePE = *min_element(fSinglePEWave.begin(), fSinglePEWave.end());
    double maxADC_SinglePE = *max_element(fSinglePEWave.begin(), fSinglePEWave.end());
    fPositivePolarity = std::abs(maxADC_SinglePE) > std::abs(minADC_SinglePE);
  
 
    // get ADC saturation value
    // currently assumes all dynamic range for PE (no overshoot)
    fADCSaturation = (fPositivePolarity ? fParams.PMTBaseline + fParams.PMTADCDynamicRange : fParams.PMTBaseline - fParams.PMTADCDynamicRange);
    file->Close();

    // Initialize noise file 
    std::string fname_noise;
    cet::search_path sp_noise("FW_SEARCH_PATH");
    if(fParams.UseDataNoise){
      std::cout << " Trying to open file " << fParams.OpDetNoiseFile << std::endl;
      sp_noise.find_file(fParams.OpDetNoiseFile, fname_noise);
      noise_file = TFile::Open(fname_noise.c_str(), "READ");
    }

  } // end constructor

  DigiPMTSBNDAlg::~DigiPMTSBNDAlg(){
    if(fParams.UseDataNoise) noise_file->Close();
  }


  void DigiPMTSBNDAlg::ConstructWaveformUncoatedPMT(
    int ch,
    sim::SimPhotons const& simphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformUncoatedPMT(simphotons, start_time, waves, ch, pdtype);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiPMTSBNDAlg::ConstructWaveformCoatedPMT(
    int ch,
    std::vector<short unsigned int>& waveform,
    std::unordered_map<int, sim::SimPhotons>& DirectPhotonsMap,
    std::unordered_map<int, sim::SimPhotons>& ReflectedPhotonsMap,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformCoatedPMT(ch, start_time, waves, DirectPhotonsMap, ReflectedPhotonsMap);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiPMTSBNDAlg::ConstructWaveformLiteUncoatedPMT(
    int ch,
    sim::SimPhotonsLite const& litesimphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformLiteUncoatedPMT(litesimphotons, start_time, waves, ch, pdtype);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiPMTSBNDAlg::ConstructWaveformLiteCoatedPMT(
    int ch,
    std::vector<short unsigned int>& waveform,
    std::unordered_map<int, sim::SimPhotonsLite>& DirectPhotonsMap,
    std::unordered_map<int, sim::SimPhotonsLite>& ReflectedPhotonsMap,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformLiteCoatedPMT(ch, start_time, waves, DirectPhotonsMap, ReflectedPhotonsMap);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiPMTSBNDAlg::CreatePDWaveformUncoatedPMT(
    sim::SimPhotons const& simphotons,
    double t_min,
    std::vector<double>& wave,
    int ch,
    std::string pdtype)
  {
    double ttsTime = 0;
    double tphoton;
    double ttpb=0;
    double _PMTUncoatedEff;
    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);
    _PMTUncoatedEff = (ch % 2 == 0) ? fPMTUncoatedEff_tpc0 : fPMTUncoatedEff_tpc1;
    for(size_t i = 0; i < simphotons.size(); i++) { //simphotons is here reflected light. To be added for all PMTs
      if(fFlatGen.fire(1.0) < _PMTUncoatedEff) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS);
        ttpb = fTimeTPB->fire(); //for including TPB emission time

        //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
        tphoton = ttsTime + simphotons[i].Time - t_min + ttpb + fParams.CableTime;

        // store the photon time if it's within the readout window
        if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++;
      }
    }

    for(size_t t=0; t<nPE_v.size(); t++){
      if(nPE_v[t] > 0) {
        if(fParams.SimulateNonLinearity){
          AddSPE(t, wave, ch, fPMTNonLinearityPtr->NObservedPE(ch, t, nPE_v) );
        }
        else{
          AddSPE(t, wave, ch, nPE_v[t]);
        }
      }
    }

    if(fParams.UseDataNoise) AddDataNoise(wave, ch);
    else
    {
      if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave, ch);
    }  

    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave, ch);

    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::CreatePDWaveformCoatedPMT(
    int ch,
    double t_min,
    std::vector<double>& wave,
    std::unordered_map<int, sim::SimPhotons>& DirectPhotonsMap,
    std::unordered_map<int, sim::SimPhotons>& ReflectedPhotonsMap)
  {

    double ttsTime = 0;
    double tphoton;
    double ttpb=0;
    double _PMTCoatedVUVEff;
    double _PMTCoatedVISEff;
    sim::SimPhotons auxphotons;

    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);

    //direct light
    _PMTCoatedVUVEff = (ch % 2 == 0) ? fPMTCoatedVUVEff_tpc0 : fPMTCoatedVUVEff_tpc1;
    if(auto it{ DirectPhotonsMap.find(ch) }; it != std::end(DirectPhotonsMap) )
    {auxphotons = it->second;}
    for(size_t j = 0; j < auxphotons.size(); j++) { //auxphotons is direct light
      if(fFlatGen.fire(1.0) < _PMTCoatedVUVEff) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
        ttpb = fTimeTPB->fire(); //for including TPB emission time

        //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
        tphoton = ttsTime + auxphotons[j].Time - t_min + ttpb + fParams.CableTime;

        // store the pgoton time if it's within the readout window
        if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++; 
      }
    }

    // reflected light
    _PMTCoatedVISEff = (ch % 2 == 0) ? fPMTCoatedVISEff_tpc0 : fPMTCoatedVISEff_tpc1;
    if(auto it{ ReflectedPhotonsMap.find(ch) }; it != std::end(ReflectedPhotonsMap) )
    {auxphotons = it->second;}
    for(size_t j = 0; j < auxphotons.size(); j++) { //auxphotons is now reflected light
      if(fFlatGen.fire(1.0) < _PMTCoatedVISEff) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
        ttpb = fTimeTPB->fire(); //for including TPB emission time

        //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
        tphoton = ttsTime + auxphotons[j].Time - t_min + ttpb + fParams.CableTime;
        
        // store the pgoton time if it's within the readout window
        if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++;
      }
    }

    for(size_t t=0; t<nPE_v.size(); t++){
      if(nPE_v[t] > 0) {
        if(fParams.SimulateNonLinearity){
          AddSPE(t, wave, ch, fPMTNonLinearityPtr->NObservedPE(ch, t, nPE_v) );
        }
        else{
          AddSPE(t, wave, ch, nPE_v[t]);
        }
      }
    }

    //Adding noise and saturation
    if(fParams.UseDataNoise) AddDataNoise(wave, ch);
    else
    {
      if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave, ch);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave, ch);

    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::CreatePDWaveformLiteUncoatedPMT(
    sim::SimPhotonsLite const& litesimphotons,
    double t_min,
    std::vector<double>& wave,
    int ch,
    std::string pdtype)
  {
    double mean_photons;
    size_t accepted_photons;
    double ttsTime = 0;
    double tphoton;
    double ttpb=0;
    double _PMTUncoatedEff;
    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);
    // here litesimphotons corresponds only to reflected light
    _PMTUncoatedEff = (ch % 2 == 0) ? fPMTUncoatedEff_tpc0 : fPMTUncoatedEff_tpc1;
    std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;
    for (auto const& reflectedPhotons : photonMap) {
      // TODO: check that this new approach of not using the last
      // (1-accepted_photons) doesn't introduce some bias. ~icaza
      mean_photons = reflectedPhotons.second*_PMTUncoatedEff;
      accepted_photons = fPoissonQGen.fire(mean_photons);
      for(size_t i = 0; i < accepted_photons; i++) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
        ttpb = fTimeTPB->fire(); // TPB emission time (cathode foils)

        //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
        tphoton = ttsTime + reflectedPhotons.first - t_min + ttpb + fParams.CableTime;
        
        // store the pgoton time if it's within the readout window
        if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++;
      }
    }

    for(size_t t=0; t<nPE_v.size(); t++){
      if(nPE_v[t] > 0) {
        if(fParams.SimulateNonLinearity){
          AddSPE(t, wave, ch, fPMTNonLinearityPtr->NObservedPE(ch, t, nPE_v) );
        }
        else{
          AddSPE(t, wave, ch, nPE_v[t]);
        }
      }
    }
    
    if(fParams.UseDataNoise) AddDataNoise(wave, ch);
    else
    {
      if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave, ch);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave, ch);

    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::CreatePDWaveformLiteCoatedPMT(
    int ch,
    double t_min,
    std::vector<double>& wave,
    std::unordered_map<int, sim::SimPhotonsLite>& DirectPhotonsMap,
    std::unordered_map<int, sim::SimPhotonsLite>& ReflectedPhotonsMap)
  {

    double mean_photons;
    size_t accepted_photons;
    double ttsTime = 0;
    double tphoton;
    double ttpb;
    double _PMTCoatedVUVEff;
    double _PMTCoatedVISEff;
    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);
    // direct light
    _PMTCoatedVUVEff = (ch % 2 == 0) ? fPMTCoatedVUVEff_tpc0 : fPMTCoatedVUVEff_tpc1;
    if ( auto it{ DirectPhotonsMap.find(ch) }; it != std::end(DirectPhotonsMap) ){
      for (auto& directPhotons : (it->second).DetectedPhotons) {
        // TODO: check that this new approach of not using the last
        // (1-accepted_photons) doesn't introduce some bias. ~icaza
        mean_photons = directPhotons.second*_PMTCoatedVUVEff;
        accepted_photons = fPoissonQGen.fire(mean_photons);
        for(size_t i = 0; i < accepted_photons; i++) {
          if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
          ttpb = fTimeTPB->fire(); // TPB emission time (PMT coating)

          //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
          tphoton = ttsTime + directPhotons.first - t_min + ttpb + fParams.CableTime;

          // store the pgoton time if it's within the readout window
          if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++;
        }
      }
    }

    // reflected light
    _PMTCoatedVISEff = (ch % 2 == 0) ? fPMTCoatedVISEff_tpc0 : fPMTCoatedVISEff_tpc1;
    if ( auto it{ ReflectedPhotonsMap.find(ch) }; it != std::end(ReflectedPhotonsMap) ){
      for (auto& reflectedPhotons : (it->second).DetectedPhotons) {
        // TODO: check that this new approach of not using the last
        // (1-accepted_photons) doesn't introduce some bias. ~icaza
        mean_photons = reflectedPhotons.second*_PMTCoatedVISEff;
        accepted_photons = fPoissonQGen.fire(mean_photons);
        for(size_t i = 0; i < accepted_photons; i++) {
          if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
          ttpb = fTimeTPB->fire(); // TPB emission time (in the cathode foils)
          
          //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
          tphoton = ttsTime + reflectedPhotons.first - t_min + ttpb + fParams.CableTime;

          // store the pgoton time if it's within the readout window
          if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++;
        }
      }
    }
    
    for(size_t t=0; t<nPE_v.size(); t++){
      if(nPE_v[t] > 0) {
        if(fParams.SimulateNonLinearity){
          AddSPE(t, wave, ch, fPMTNonLinearityPtr->NObservedPE(ch, t, nPE_v) );
        }
        else{
          AddSPE(t, wave, ch, nPE_v[t]);
        }
      }
    }

    //Adding noise and saturation
    if(fParams.UseDataNoise) AddDataNoise(wave, ch);
    else
    {
      if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave, ch);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave, ch);
    
    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::Pulse1PE(std::vector<double>& fSinglePEWave)//single pulse waveform
  {
    double time;
    double constT1 = fParams.PMTChargeToADC * fParams.PMTMeanAmplitude;
    double constT21 = 2.0 * sigma1 * sigma1;
    double constT22 = 2.0 * sigma2 * sigma2;
    for(size_t i = 0; i<fSinglePEWave.size(); i++) {
      time = static_cast<double>(i) / fSampling;
      if (time < fParams.TransitTime)
        fSinglePEWave[i] = constT1 * std::exp(-1.0 * std::pow(time - fParams.TransitTime, 2) / constT21);
      else
        fSinglePEWave[i] = constT1 * std::exp(-1.0 * std::pow(time - fParams.TransitTime, 2) / constT22);
    }
  }


  double DigiPMTSBNDAlg::Transittimespread(double fwhm)
  {
    double tts, sigma;
    sigma = fwhm / transitTimeSpread_frac;
    tts = fGaussQGen.fire(0., sigma);
    return tts;
  }


  void DigiPMTSBNDAlg::AddSPE(size_t time, std::vector<double>& wave, int ch, double npe)
  {
    if(!fPMTCalibrationDatabaseService->getReconstructChannel(ch)) return;
    // time bin HD (double precision)
    // used to gert the time-shifted SER
    double time_bin_hd = fSampling*time;
    size_t wvf_shift  = fPMTHDOpticalWaveformsPtr->TimeBinShift(time_bin_hd);
    
    // get actual time bin and waveform min/max iterators
    size_t time_bin=std::floor(time_bin_hd);
    size_t max = time_bin + pulsesize < wave.size() ? time_bin + pulsesize : wave.size();
    auto min_it = std::next(wave.begin(), time_bin);
    auto max_it = std::next(wave.begin(), max);
    
    // simulate gain fluctuations
    double npe_anode = npe;
    if(fParams.MakeGainFluctuations)
      //npe_anode=fPMTGainFluctuationsPtr->GainFluctuation(npe, fEngine);
      npe_anode=fPMTGainFluctuationsPtr->GainFluctuation(ch, npe, fEngine);
    // add SER to the waveform
    std::transform(min_it, max_it,
                     fSinglePEWave_HD[ch][wvf_shift].begin(), min_it,
                     [npe_anode](auto a, auto b) { return a+npe_anode*b; });
  }


  void DigiPMTSBNDAlg::CreateSaturation(std::vector<double>& wave)
  {
    if(fPositivePolarity)
      std::replace_if(wave.begin(), wave.end(),
                      [&](auto w){return w > fADCSaturation;}, fADCSaturation);
    else
      std::replace_if(wave.begin(), wave.end(),
                        [&](auto w){return w < fADCSaturation;}, fADCSaturation);
  }


  void DigiPMTSBNDAlg::AddLineNoise(std::vector<double>& wave, int ch)
  {
    // TODO: after running the profiler I can see that this is where
    // most cycles are being used.  Potentially some improvement could
    // be achieved with the below code but the array dynamic allocation
    // is not helping. The function would need to have it passed.
    // ~icaza
    //
    // double *array = new double[wave.size()]();
    // CLHEP::RandGaussQ::shootArray(fEngine, wave.size(), array, 0, fParams.PMTBaselineRMS);
    // for(size_t i = 0; i<wave.size(); i++) {
    //   wave[i] += array[i];
    // }
    // delete array;
    //
    std::transform(wave.begin(), wave.end(), wave.begin(),
                   [this](double w) -> double {
                     return w + fGaussQGen.fire(0., fParams.PMTBaselineRMS) ; });
  }


  void DigiPMTSBNDAlg::AddDataNoise(std::vector<double>& wave, int ch)
  {
    // Get the noise waveform from data
    // Get all the noise waveforms that correspond to the channel we are using
    // Select a random noise waveform
    // Compare the length of the noise and the simulation waveform:
    // if the simulation length is larger then choose another random noise waveform to fill up the missing items

    if (!noise_file || noise_file->IsZombie()) {
        throw cet::exception("DigiPMTSBNDAlg") << " PMT Noise file could not be opened " << std::endl;
        return;
    }

    TList *keys = noise_file->GetListOfKeys();
    if (!keys || keys->IsEmpty()) {
        throw cet::exception("DigiPMTSBNDAlg") << " PMT Noise file is empty " << std::endl;
        return;
    }

    std::string opChannelName = "opchannel_" + std::to_string(ch) + "_pmt";
    std::vector<TKey*> keylist;
    TIter nextkey(keys);
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        std::string channel_name = key->GetName();
        if (channel_name.find(opChannelName) != std::string::npos) {
            keylist.push_back(key);
        }
    }

    if (keylist.empty()) {
        throw cet::exception("DigiPMTSBNDAlg") << " PMT Noise file has no data for channel " << ch << std::endl;
        return;
    }

    std::vector<double> noise_wform;
    int waveBins = wave.size();
    int currentSize = 0;

    // Llenar el vector de ruido hasta igualar la longitud del waveform
    while (currentSize < waveBins) {
        int noiseWformIdx = static_cast<int>(fEngine->flat() * keylist.size());
        TH1 *noiseHist = (TH1*)keylist[noiseWformIdx]->ReadObj();

        for (int i = 1; i <= noiseHist->GetNbinsX(); i++) {
            noise_wform.push_back(noiseHist->GetBinContent(i));
            currentSize += 1;
            if (currentSize >= waveBins) break;
        }
        delete noiseHist;
    }

    // Agregar ruido al waveform
    for (size_t i = 0; i < wave.size(); i++) {
        wave[i] += noise_wform[i] + fParams.PMTBaseline;
    }
  }


  void DigiPMTSBNDAlg::AddDarkNoise(std::vector<double>& wave, int ch)
  {
    double timeBin;
    // Multiply by 10^9 since fParams.DarkNoiseRate is in Hz (conversion from s to ns)
    double mean =  1000000000.0 / fParams.PMTDarkNoiseRate;
    double darkNoiseTime = fExponentialGen.fire(mean);
    while(darkNoiseTime < wave.size()) {
      timeBin = std::round(darkNoiseTime);
      if(timeBin < wave.size()) {AddSPE(fSamplingPeriod*timeBin, wave, ch);}
      // Find next time to add dark noise
      darkNoiseTime += fExponentialGen.fire(mean);
    }
  }


  // TODO: this function is not being used anywhere! ~icaza
  double DigiPMTSBNDAlg::FindMinimumTime(
    sim::SimPhotons const& simphotons,
    int ch,
    std::string pdtype,
    std::unordered_map<int, sim::SimPhotons>& directPhotonsOnPMTS)
  {
    double t_min = 1e15;

    if(pdtype == "pmt_uncoated") {
      // TODO: use std::algorithm.  ~icaza
      for(size_t i = 0; i < simphotons.size(); i++) {
        if(simphotons[i].Time < t_min) t_min = simphotons[i].Time;
      }
    }
    else if(pdtype == "pmt_coated") {
      sim::SimPhotons auxphotons;
      if ( auto it{ directPhotonsOnPMTS.find(ch) }; it != std::end(directPhotonsOnPMTS) )
      {auxphotons = it->second;}
      auxphotons += (simphotons);
      // TODO: use std::algorithm.  ~icaza
      for(size_t i = 0; i < auxphotons.size(); i++) {
        if(auxphotons[i].Time < t_min) t_min = auxphotons[i].Time;
      }
    }
    else {
      throw cet::exception("DigiPMTSBNDAlg") << "Wrong pdtype: " << pdtype << std::endl;
    }
    return t_min;
  }

  // TODO: this function is not being used anywhere! ~icaza
  double DigiPMTSBNDAlg::FindMinimumTimeLite(
    sim::SimPhotonsLite const& litesimphotons,
    int ch,
    std::string pdtype,
    std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS)
  {

    if(pdtype == "pmt_uncoated") {
      std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;
      auto min = std::find_if(
        photonMap.begin(),
        photonMap.end(),
        [](const auto& pm) {return pm.second != 0; });
      if(min != photonMap.end()) return double(min->first);
    }
    else if(pdtype == "pmt_coated") {
      sim::SimPhotonsLite auxphotons;
      if ( auto it{ directPhotonsOnPMTS.find(ch) }; it != std::end(directPhotonsOnPMTS) )
      {auxphotons = it->second;}
      // TODO: this might be buggy:
      // potentially adding to a uninitialized object.  ~icaza
      auxphotons += (litesimphotons);
      std::map<int, int> const& auxphotonMap = auxphotons.DetectedPhotons;
      auto min = std::find_if(
        auxphotonMap.begin(),
        auxphotonMap.end(),
        [](const auto& pm) {return pm.second != 0; });
      if(min != auxphotonMap.end()) return double(min->first);
    }
    else {
      throw cet::exception("DigiPMTSBNDAlg") << "Wrong pdtype: " << pdtype << std::endl;
    }
    return 1e5;
  }

  // -----------------------------------------------------------------------------
  // ---  opdet::DigiPMTSBNDAlgMaker
  // -----------------------------------------------------------------------------

  DigiPMTSBNDAlgMaker::DigiPMTSBNDAlgMaker
  (Config const& config)
  {
    // settings
    fBaseConfig.PMTChargeToADC           = config.pmtchargeToADC();
    fBaseConfig.PMTBaseline              = config.pmtbaseline();
    fBaseConfig.PMTADCDynamicRange       = config.pmtADCDynamicRange();
    fBaseConfig.PMTCoatedVUVEff_tpc0     = config.pmtcoatedVUVEff_tpc0();
    fBaseConfig.PMTCoatedVISEff_tpc0     = config.pmtcoatedVISEff_tpc0();
    fBaseConfig.PMTUncoatedEff_tpc0      = config.pmtuncoatedEff_tpc0();
    fBaseConfig.PMTCoatedVUVEff_tpc1     = config.pmtcoatedVUVEff_tpc1();
    fBaseConfig.PMTCoatedVISEff_tpc1     = config.pmtcoatedVISEff_tpc1();
    fBaseConfig.PMTUncoatedEff_tpc1      = config.pmtuncoatedEff_tpc1();
    fBaseConfig.PMTSinglePEmodel         = config.PMTsinglePEmodel();
    fBaseConfig.PositivePolarity         = config.PositivePolarity();
    fBaseConfig.PMTRiseTime              = config.pmtriseTime();
    fBaseConfig.PMTFallTime              = config.pmtfallTime();
    fBaseConfig.PMTMeanAmplitude         = config.pmtmeanAmplitude();
    fBaseConfig.PMTDarkNoiseRate         = config.pmtdarkNoiseRate();
    fBaseConfig.PMTBaselineRMS           = config.pmtbaselineRMS();
    fBaseConfig.TransitTime              = config.transitTime();
    fBaseConfig.TTS                      = config.tts();
    fBaseConfig.CableTime                = config.cableTime();
    fBaseConfig.PMTDataFile              = config.pmtDataFile();
    fBaseConfig.UseDataNoise             = config.UseDataNoise();
    fBaseConfig.OpDetNoiseFile           = config.OpDetNoiseFile();
    fBaseConfig.MakeGainFluctuations = config.gainFluctuationsParams.get_if_present(fBaseConfig.GainFluctuationsParams);
    fBaseConfig.SimulateNonLinearity = config.nonLinearityParams.get_if_present(fBaseConfig.NonLinearityParams);
    config.hdOpticalWaveformParams.get_if_present(fBaseConfig.HDOpticalWaveformParams);
  }

  std::unique_ptr<DigiPMTSBNDAlg>
  DigiPMTSBNDAlgMaker::operator()(
    detinfo::LArProperties const& larProp,
    detinfo::DetectorClocksData const& clockData,
    CLHEP::HepRandomEngine* engine
    ) const
  {
    // set the configuration
    auto params = fBaseConfig;

    // set up parameters
    params.larProp = &larProp;
    params.frequency = clockData.OpticalClock().Frequency();
    params.engine = engine;

    return std::make_unique<DigiPMTSBNDAlg>(params);
  } // DigiPMTSBNDAlgMaker::create()

}
