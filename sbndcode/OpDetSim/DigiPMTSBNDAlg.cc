#include "sbndcode/OpDetSim/DigiPMTSBNDAlg.hh"

//------------------------------------------------------------------------------
//--- opdet::simpmtsbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet {


  DigiPMTSBNDAlg::DigiPMTSBNDAlg(ConfigurationParameters_t const& config)
    : fParams(config)
    , fSampling(fParams.timeService->OpticalClock().Frequency())
    , fQEDirect(fParams.QEDirect / fParams.larProp->ScintPreScale())
    , fQERefl(fParams.QERefl / fParams.larProp->ScintPreScale())
      //  , fSinglePEmodel(fParams.SinglePEmodel)
    , fEngine(fParams.engine)
  {

    //art::ServiceHandle<rndm::NuRandomService> seedSvc;
    //fEngine = new CLHEP::HepJamesRandom;
    //seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), "DigiPMTSBNDAlg");

    mf::LogInfo("DigiPMTSBNDAlg") << "PMT corrected efficiencies = "
                                  << fQEDirect << " " << fQERefl;

    if(fQERefl > 1.0001 || fQEDirect > 1.0001)
      mf::LogWarning("DigiPMTSBNDAlg")
        << "Quantum efficiency set in fhicl file " << fParams.QERefl
        << " or " << fParams.QEDirect << " seems to be too large!\n"
        << "Final QE must be equal or smaller than the scintillation "
        << "pre scale applied at simulation time.\n"
        << "Please check this number (ScintPreScale): "
        << fParams.larProp->ScintPreScale();

    fSampling = fSampling / 1000.0; //in GHz, to cancel with ns

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fParams.PMTDataFile, fname);
    TFile* file = TFile::Open(fname.c_str());
    file->GetObject("timeTPB", timeTPB);

    //shape of single pulse
    if (fParams.SinglePEmodel) {
      mf::LogDebug("DigiPMTSBNDAlg") << " using testbench pe response";
      std::vector<double> *wsp_pointer;
      file->GetObject("wsp", wsp_pointer);
      wsp = *wsp_pointer;
      pulsesize = wsp.size();
    }
    else {
      mf::LogDebug("DigiPMTSBNDAlg") << " using ideal pe response";
      //shape of single pulse
      sigma1 = fParams.PMTRiseTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      sigma2 = fParams.PMTFallTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));

      pulsesize = (int)((6 * sigma2 + fParams.TransitTime) * fSampling);
      wsp.resize(pulsesize);
      Pulse1PE(wsp);
    }

    saturation = fParams.PMTBaseline + fParams.PMTSaturation * fParams.PMTChargeToADC * fParams.PMTMeanAmplitude;
  } // end constructor


  DigiPMTSBNDAlg::~DigiPMTSBNDAlg(){ }


  void DigiPMTSBNDAlg::ConstructWaveform(
    int ch,
    sim::SimPhotons const& simphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    std::unordered_map<int, sim::SimPhotons>& directPhotonsOnPMTS,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveform(simphotons, start_time, waves, ch, pdtype, directPhotonsOnPMTS);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiPMTSBNDAlg::ConstructWaveformLite(
    int ch,
    sim::SimPhotonsLite const& litesimphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformLite(litesimphotons, start_time, waves, ch, pdtype, directPhotonsOnPMTS);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiPMTSBNDAlg::CreatePDWaveform(
    sim::SimPhotons const& simphotons,
    double t_min,
    std::vector<double>& wave,
    int ch,
    std::string pdtype,
    std::unordered_map<int, sim::SimPhotons>& directPhotonsOnPMTS)
  {
    double ttsTime = 0;
    double tphoton;
    size_t timeBin;
    for(size_t i = 0; i < simphotons.size(); i++) { //simphotons is here reflected light. To be added for all PMTs
      if(CLHEP::RandFlat::shoot(fEngine, 1.0) < fQERefl) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS);
        tphoton = ttsTime + simphotons[i].Time - t_min + fParams.CableTime;
        if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
        timeBin = std::floor(tphoton*fSampling);
        if(timeBin < wave.size()) {AddSPE(timeBin, wave);}
      }
    }
    if(pdtype == "pmt_coated") { //To add direct light for TPB coated PMTs
      sim::SimPhotons auxphotons;
      double ttpb = 0;
      if(auto it{ directPhotonsOnPMTS.find(ch) }; it != std::end(directPhotonsOnPMTS) )
      {auxphotons = it->second;}
      for(size_t j = 0; j < auxphotons.size(); j++) { //auxphotons is direct light
        if(CLHEP::RandFlat::shoot(fEngine, 1.0) < fQEDirect) {
          if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
          // TODO: this uses root random machine!
          // use RandGeneral instead. ~icaza
          ttpb = timeTPB->GetRandom(); //for including TPB emission time
          tphoton = ttsTime + auxphotons[j].Time - t_min + ttpb + fParams.CableTime;
          if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
          timeBin = std::floor(tphoton*fSampling);
          if(timeBin < wave.size()) {AddSPE(timeBin, wave);}
        }
      }
    }
    if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::CreatePDWaveformLite(
    sim::SimPhotonsLite const& litesimphotons,
    double t_min,
    std::vector<double>& wave,
    int ch,
    std::string pdtype,
    std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS)
  {
    double mean_photons;
    size_t accepted_photons;
    double ttsTime = 0;
    double tphoton;
    size_t timeBin;
    // reflected light to be added to all PMTs
    std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;
    for (auto const& reflectedPhotons : photonMap) {
      // TODO: check that this new approach of not using the last
      // (1-accepted_photons) doesn't introduce some bias. ~icaza
      mean_photons = reflectedPhotons.second*fQEDirect;
      accepted_photons = CLHEP::RandPoissonQ::shoot(fEngine, mean_photons);
      for(size_t i = 0; i < accepted_photons; i++) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS);
        tphoton = ttsTime + reflectedPhotons.first - t_min + fParams.CableTime;
        if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
        timeBin = std::floor(tphoton*fSampling);
        if(timeBin < wave.size()) {AddSPE(timeBin, wave);}
      }
    }

    // direct light for TPB coated PMTs
    if(pdtype == "pmt_coated") {
      if ( auto it{ directPhotonsOnPMTS.find(ch) }; it != std::end(directPhotonsOnPMTS) ){
        double ttpb;
        for (auto& directPhotons : (it->second).DetectedPhotons) {
          // TODO: check that this new approach of not using the last
          // (1-accepted_photons) doesn't introduce some bias. ~icaza
          mean_photons = directPhotons.second*fQEDirect;
          accepted_photons = CLHEP::RandPoissonQ::shoot(fEngine, mean_photons);
          for(size_t i = 0; i < accepted_photons; i++) {
            if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
            // TODO: this uses root random machine!
            // use RandGeneral. ~icaza
            ttpb = timeTPB->GetRandom(); //for including TPB emission time
            tphoton = ttsTime + directPhotons.first - t_min + ttpb + fParams.CableTime;
            if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
            timeBin = std::floor(tphoton*fSampling);
            if(timeBin < wave.size()) {AddSPE(timeBin, wave);}
          }
        }
      }
    }

    if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::Pulse1PE(std::vector<double>& wsp)//single pulse waveform
  {
    double time;
    double constT1 = fParams.PMTChargeToADC * fParams.PMTMeanAmplitude;
    double constT21 = 2.0 * sigma1 * sigma1;
    double constT22 = 2.0 * sigma2 * sigma2;
    for(size_t i = 0; i<wsp.size(); i++) {
      time = static_cast<double>(i) / fSampling;
      if (time < fParams.TransitTime)
        wsp[i] = constT1 * std::exp(-1.0 * std::pow(time - fParams.TransitTime, 2) / constT21);
      else
        wsp[i] = constT1 * std::exp(-1.0 * std::pow(time - fParams.TransitTime, 2) / constT22);
    }
  }


  double DigiPMTSBNDAlg::Transittimespread(double fwhm)
  {
    double tts, sigma;
    sigma = fwhm / transitTimeSpread_frac;
    tts = CLHEP::RandGaussQ::shoot(fEngine, 0, sigma);
    return tts;
  }


  void DigiPMTSBNDAlg::AddSPE(size_t time_bin, std::vector<double>& wave)
  {
    size_t max = time_bin + pulsesize < wave.size() ? time_bin + pulsesize : wave.size();
    auto min_it = std::next(wave.begin(), time_bin);
    auto max_it = std::next(wave.begin(), max);
    std::transform(min_it, max_it,
                   wsp.begin(), min_it,
                   std::plus<double>( ));
  }


  void DigiPMTSBNDAlg::CreateSaturation(std::vector<double>& wave)
  {
    std::replace_if(wave.begin(), wave.end(),
                    [&](auto w){return w < saturation;}, saturation);
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
                     return w + CLHEP::RandGaussQ::shoot(fEngine, 0, fParams.PMTBaselineRMS) ; });
  }


  void DigiPMTSBNDAlg::AddDarkNoise(std::vector<double>& wave)
  {
    size_t timeBin;
    // Multiply by 10^9 since fParams.DarkNoiseRate is in Hz (conversion from s to ns)
    double mean =  1000000000.0 / fParams.PMTDarkNoiseRate;
    double darkNoiseTime = CLHEP::RandExponential::shoot(fEngine, mean);
    while(darkNoiseTime < wave.size()) {
      timeBin = std::round(darkNoiseTime);
      if(timeBin < wave.size()) {AddSPE(timeBin, wave);}
      // Find next time to add dark noise
      darkNoiseTime += CLHEP::RandExponential::shoot(fEngine, mean);
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
    fBaseConfig.PMTSaturation            = config.pmtsaturation();
    fBaseConfig.QEDirect                 = config.qEDirect();
    fBaseConfig.QERefl                   = config.qERefl();
    fBaseConfig.SinglePEmodel            = config.singlePEmodel();
    fBaseConfig.PMTRiseTime              = config.pmtriseTime();
    fBaseConfig.PMTFallTime              = config.pmtfallTime();
    fBaseConfig.PMTMeanAmplitude         = config.pmtmeanAmplitude();
    fBaseConfig.PMTDarkNoiseRate         = config.pmtdarkNoiseRate();
    fBaseConfig.PMTBaselineRMS           = config.pmtbaselineRMS();
    fBaseConfig.TransitTime              = config.transitTime();
    fBaseConfig.TTS                      = config.tts();
    fBaseConfig.CableTime                = config.cableTime();
    fBaseConfig.PMTDataFile              = config.pmtDataFile();
  }

  std::unique_ptr<DigiPMTSBNDAlg>
  DigiPMTSBNDAlgMaker::operator()(
    detinfo::LArProperties const& larProp,
    detinfo::DetectorClocks const& detClocks,
    CLHEP::HepRandomEngine* engine
  ) const
  {
    // set the configuration
    auto params = fBaseConfig;

    // set up parameters
    params.larProp = &larProp;
    params.timeService = &detClocks;
    params.engine = engine;

    return std::make_unique<DigiPMTSBNDAlg>(params);
  } // DigiPMTSBNDAlgMaker::create()

}
