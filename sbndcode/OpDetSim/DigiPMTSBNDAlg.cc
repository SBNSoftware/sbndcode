#include "sbndcode/OpDetSim/DigiPMTSBNDAlg.hh"

//------------------------------------------------------------------------------
//--- opdet::simpmtsbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet {


  DigiPMTSBNDAlg::DigiPMTSBNDAlg(ConfigurationParameters_t const& config)
    : fParams(config)
    , fSampling(fParams.frequency)
    , fPMTCoatedVUVEff(fParams.PMTCoatedVUVEff / fParams.larProp->ScintPreScale())
    , fPMTCoatedVISEff(fParams.PMTCoatedVISEff / fParams.larProp->ScintPreScale())
    , fPMTUncoatedEff(fParams.PMTUncoatedEff/ fParams.larProp->ScintPreScale())
    , fUseDataNoise(fParams.UseDataNoise)
    , fUseCustomSPEHeight(fParams.UseCustomSPEHeight)
    , fOpDetNoiseFile(fParams.OpDetNoiseFile)
    , fSPEPeakHeight(fParams.SPEPeakHeight)
      //  , fSinglePEmodel(fParams.SinglePEmodel)
    , fEngine(fParams.engine)
    , fFlatGen(*fEngine)
    , fPoissonQGen(*fEngine)
    , fGaussQGen(*fEngine)
    , fExponentialGen(*fEngine)
  {

    mf::LogInfo("DigiPMTSBNDAlg") << "PMT corrected efficiencies = "
                                  << fPMTCoatedVUVEff << " " << fPMTCoatedVISEff << " " << fPMTUncoatedEff <<"\n";

    if(fPMTCoatedVUVEff > 1.0001 || fPMTCoatedVISEff > 1.0001 || fPMTUncoatedEff > 1.0001)
      mf::LogWarning("DigiPMTSBNDAlg")
        << "Detection efficiencies set in fhicl file seem to be too large!\n"
        << "PMTCoatedVUVEff: " << fParams.PMTCoatedVUVEff << "\n"
        << "PMTCoatedVISEff: " << fParams.PMTCoatedVISEff << "\n"
        << "PMTUncoatedEff: " << fParams.PMTUncoatedEff << "\n"
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

    // TPB emission time histogram for pmt_coated histogram
    std::vector<double>* timeTPB_p;
    file->GetObject("timeTPB", timeTPB_p);
    fTimeTPB = std::make_unique<CLHEP::RandGeneral>
      (*fEngine, timeTPB_p->data(), timeTPB_p->size());



    //shape of single pulse
    if (fParams.PMTSinglePEmodel) {
      mf::LogDebug("DigiPMTSBNDAlg") << " using testbench pe response";
      std::vector<double>* SinglePEVec_p;
      file->GetObject("SinglePEVec_HD", SinglePEVec_p);
    fSinglePEWave = *SinglePEVec_p;
      if(fUseCustomSPEHeight)
      {
        double min_value = *std::min_element(fSinglePEWave.begin(), fSinglePEWave.end());
        double SPENormalization = std::abs(fSPEPeakHeight/min_value);
        // Transform the SPE response to have a given peak height
        std::transform(fSinglePEWave.begin(), fSinglePEWave.end(), fSinglePEWave.begin(),
                    [SPENormalization](double x) { return x * SPENormalization; });
      }

      // Prepare HD waveforms
      fPMTHDOpticalWaveformsPtr = art::make_tool<opdet::HDOpticalWaveform>(fParams.HDOpticalWaveformParams);
      fPMTHDOpticalWaveformsPtr->produceSER_HD(fSinglePEWave_HD,fSinglePEWave);

      pulsesize = fSinglePEWave_HD[0].size();
      mf::LogDebug("DigiPMTSBNDAlg")<<"HD wvfs size: "<<pulsesize;
    }
    else {
      mf::LogDebug("DigiPMTSBNDAlg") << " using ideal pe response";
      //shape of single pulse
      sigma1 = fParams.PMTRiseTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      sigma2 = fParams.PMTFallTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));

      pulsesize = (int)((6 * sigma2 + fParams.TransitTime) * fSampling);
      fSinglePEWave.resize(pulsesize);
      Pulse1PE(fSinglePEWave);
    }

    if(fParams.MakeGainFluctuations){
      fPMTGainFluctuationsPtr = art::make_tool<opdet::PMTGainFluctuations>(fParams.GainFluctuationsParams);
    }
    if(fParams.SimulateNonLinearity){
      fPMTNonLinearityPtr = art::make_tool<opdet::PMTNonLinearity>(fParams.NonLinearityParams);
    }  

    // infer pulse polarity from SER peak sign
    double minADC_SinglePE = *min_element(fSinglePEWave.begin(), fSinglePEWave.end());
    double maxADC_SinglePE = *max_element(fSinglePEWave.begin(), fSinglePEWave.end());
    fPositivePolarity = std::abs(maxADC_SinglePE) > std::abs(minADC_SinglePE); 
    
    // get ADC saturation value
    // currently assumes all dynamic range for PE (no overshoot)
    fADCSaturation = (fPositivePolarity ? fParams.PMTBaseline + fParams.PMTADCDynamicRange : fParams.PMTBaseline - fParams.PMTADCDynamicRange);


    file->Close();
  } // end constructor


  DigiPMTSBNDAlg::~DigiPMTSBNDAlg(){}


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

    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);

    for(size_t i = 0; i < simphotons.size(); i++) { //simphotons is here reflected light. To be added for all PMTs
      if(fFlatGen.fire(1.0) < fPMTUncoatedEff) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS);
        ttpb = fTimeTPB->fire(); //for including TPB emission time

        //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
        tphoton = ttsTime + simphotons[i].Time - t_min + ttpb + fParams.CableTime;

        // store the pgoton time if it's within the readout window
        if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++;
      }
    }

    for(size_t t=0; t<nPE_v.size(); t++){
      if(nPE_v[t] > 0) {
        if(fParams.SimulateNonLinearity){
          AddSPE(t, wave, fPMTNonLinearityPtr->NObservedPE(t, nPE_v) );
        }
        else{
          AddSPE(t, wave, nPE_v[t]);
        }
      }
    }

    if(fParams.PMTBaselineRMS > 0.0) 
    {
      if(fUseDataNoise)
      {
        std::cout << " Using data noise for channel " << ch << std::endl;
        AddDataNoise(wave, ch);
      } 
      else AddLineNoise(wave);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
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
    sim::SimPhotons auxphotons;

    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);

    //direct light
    if(auto it{ DirectPhotonsMap.find(ch) }; it != std::end(DirectPhotonsMap) )
    {auxphotons = it->second;}
    for(size_t j = 0; j < auxphotons.size(); j++) { //auxphotons is direct light
      if(fFlatGen.fire(1.0) < fPMTCoatedVUVEff) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
        ttpb = fTimeTPB->fire(); //for including TPB emission time

        //photon time in ns (w.r.t. the waveform start time a.k.a t_min)
        tphoton = ttsTime + auxphotons[j].Time - t_min + ttpb + fParams.CableTime;

        // store the pgoton time if it's within the readout window
        if(tphoton > 0 && tphoton < nPE_v.size()) nPE_v[(size_t)tphoton]++; 
      }
    }

    // reflected light
    if(auto it{ ReflectedPhotonsMap.find(ch) }; it != std::end(ReflectedPhotonsMap) )
    {auxphotons = it->second;}
    for(size_t j = 0; j < auxphotons.size(); j++) { //auxphotons is now reflected light
      if(fFlatGen.fire(1.0) < fPMTCoatedVISEff) {
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
          AddSPE(t, wave, fPMTNonLinearityPtr->NObservedPE(t, nPE_v) );
        }
        else{
          AddSPE(t, wave, nPE_v[t]);
        }
      }
    }

    //Adding noise and saturation
    if(fParams.PMTBaselineRMS > 0.0) 
    {
      if(fUseDataNoise) AddDataNoise(wave, ch);
      else AddLineNoise(wave);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
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

    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);

    // here litesimphotons corresponds only to reflected light
    std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;
    for (auto const& reflectedPhotons : photonMap) {
      // TODO: check that this new approach of not using the last
      // (1-accepted_photons) doesn't introduce some bias. ~icaza
      mean_photons = reflectedPhotons.second*fPMTUncoatedEff;
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
          AddSPE(t, wave, fPMTNonLinearityPtr->NObservedPE(t, nPE_v) );
        }
        else{
          AddSPE(t, wave, nPE_v[t]);
        }
      }
    }
    
    if(fParams.PMTBaselineRMS > 0.0) 
    {
      if(fUseDataNoise) AddDataNoise(wave, ch);
      else AddLineNoise(wave);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
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

    // we want to keep the 1 ns SimPhotonLite resolution
    // digitizer sampling period is 2 ns
    // create a PE accumulator vector with size x2 the waveform size
    std::vector<unsigned int> nPE_v( (size_t) fSamplingPeriod*wave.size(), 0);

    // direct light
    if ( auto it{ DirectPhotonsMap.find(ch) }; it != std::end(DirectPhotonsMap) ){
      for (auto& directPhotons : (it->second).DetectedPhotons) {
        // TODO: check that this new approach of not using the last
        // (1-accepted_photons) doesn't introduce some bias. ~icaza
        mean_photons = directPhotons.second*fPMTCoatedVUVEff;
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
    if ( auto it{ ReflectedPhotonsMap.find(ch) }; it != std::end(ReflectedPhotonsMap) ){
      for (auto& reflectedPhotons : (it->second).DetectedPhotons) {
        // TODO: check that this new approach of not using the last
        // (1-accepted_photons) doesn't introduce some bias. ~icaza
        mean_photons = reflectedPhotons.second*fPMTCoatedVISEff;
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
          AddSPE(t, wave, fPMTNonLinearityPtr->NObservedPE(t, nPE_v) );
        }
        else{
          AddSPE(t, wave, nPE_v[t]);
        }
      }
    }

    //Adding noise and saturation
    if(fParams.PMTBaselineRMS > 0.0) 
    {
      if(fUseDataNoise) AddDataNoise(wave, ch);
      else AddLineNoise(wave);
    }
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
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


  void DigiPMTSBNDAlg::AddSPE(size_t time, std::vector<double>& wave, double npe)
  {
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
      npe_anode=fPMTGainFluctuationsPtr->GainFluctuation(npe, fEngine);

    // add SER to the waveform
    std::transform(min_it, max_it,
                     fSinglePEWave_HD[wvf_shift].begin(), min_it,
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


  void DigiPMTSBNDAlg::AddLineNoise(std::vector<double>& wave)
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
    // Get all the noise waveforms that do correspond to the channel we are using 
    // Select a random noise waveform
    // Compare the length of the noise and the simulation waveform:
    // if the simulation length is larger then choose another random noise waveform to fill up the missing items

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fOpDetNoiseFile, fname);
    TFile* file = TFile::Open(fname.c_str(), "READ");

    TDirectory *dir = (TDirectory*)file->Get("flashana");
    if (!dir) {
        std::cerr << " File could not be opened" << std::endl;
        file->Close();
        return;
    }

    std::string opChannelName = "opchannel_" + std::to_string(ch)+"_pmt";
    std::vector<TKey*> keylist;
    for (auto it = dir->GetListOfKeys()->begin(); it != dir->GetListOfKeys()->end(); ++it) {
        TKey *key = (TKey*)(*it);
        std::string channel_name = key->GetName();
        if (channel_name.find(opChannelName.c_str()) != std::string::npos) {
          //Creating a list of keys that contains the channel number
            keylist.push_back(key);
        }
    }

    std::vector<double> noise_wform; 
    int waveBins = wave.size();

    int currentSize=0;
    while(currentSize<waveBins)
    {
      int noiseWformIdx = static_cast<int>(fEngine->flat() * keylist.size());
      TH1 *noiseHist = (TH1*)keylist[noiseWformIdx]->ReadObj();
      for(int i=1; i<=noiseHist->GetNbinsX(); i++)
      {
        noise_wform.push_back(noiseHist->GetBinContent(i));
        currentSize+=1;
        if(currentSize>=waveBins) break;
      }
      delete noiseHist;
    }

    for(size_t i=0; i<wave.size();i++)
    {
      wave[i]+=noise_wform[i];
    }

    //std::transform(wave.begin(), wave.end(), noise_wform.begin(), noise_wform.begin(), std::plus<double>());

    file->Close();
  }


  void DigiPMTSBNDAlg::AddDarkNoise(std::vector<double>& wave)
  {
    double timeBin;
    // Multiply by 10^9 since fParams.DarkNoiseRate is in Hz (conversion from s to ns)
    double mean =  1000000000.0 / fParams.PMTDarkNoiseRate;
    double darkNoiseTime = fExponentialGen.fire(mean);
    while(darkNoiseTime < wave.size()) {
      timeBin = std::round(darkNoiseTime);
      if(timeBin < wave.size()) {AddSPE(fSamplingPeriod*timeBin, wave);}
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
    fBaseConfig.PMTCoatedVUVEff          = config.pmtcoatedVUVEff();
    fBaseConfig.PMTCoatedVISEff          = config.pmtcoatedVISEff();
    fBaseConfig.PMTUncoatedEff           = config.pmtuncoatedEff();
    fBaseConfig.PMTSinglePEmodel         = config.PMTsinglePEmodel();
    fBaseConfig.PMTRiseTime              = config.pmtriseTime();
    fBaseConfig.PMTFallTime              = config.pmtfallTime();
    fBaseConfig.PMTMeanAmplitude         = config.pmtmeanAmplitude();
    fBaseConfig.PMTDarkNoiseRate         = config.pmtdarkNoiseRate();
    fBaseConfig.UseDataNoise             = config.UseDataNoise();
    fBaseConfig.UseCustomSPEHeight       = config.UseCustomSPEHeight();
    fBaseConfig.SPEPeakHeight            = config.SPEPeakHeight();
    fBaseConfig.OpDetNoiseFile           = config.OpDetNoiseFile();
    fBaseConfig.PMTBaselineRMS           = config.pmtbaselineRMS();
    fBaseConfig.TransitTime              = config.transitTime();
    fBaseConfig.TTS                      = config.tts();
    fBaseConfig.CableTime                = config.cableTime();
    fBaseConfig.PMTDataFile              = config.pmtDataFile();
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
