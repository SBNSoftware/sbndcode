#include "sbndcode/OpDetSim/DigiArapucaSBNDAlg.hh"

//------------------------------------------------------------------------------
//--- opdet::simarapucasbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet {

  DigiArapucaSBNDAlg::DigiArapucaSBNDAlg(ConfigurationParameters_t const& config)
    : fParams(config)
    , fSampling(fParams.frequency)
    , fArapucaVUVEff(fParams.ArapucaVUVEff / fParams.larProp->ScintPreScale())
    , fArapucaVISEff(fParams.ArapucaVISEff / fParams.larProp->ScintPreScale())
    , fXArapucaVUVEff(fParams.XArapucaVUVEff / fParams.larProp->ScintPreScale())
    , fXArapucaVISEff(fParams.XArapucaVISEff / fParams.larProp->ScintPreScale())
    , fEngine(fParams.engine)
  {

    if(fArapucaVUVEff > 1.0001 || fArapucaVISEff > 1.0001 ||
       fXArapucaVUVEff > 1.0001 || fXArapucaVISEff > 1.0001)
      mf::LogWarning("DigiArapucaSBNDAlg")
        << "Quantum efficiency set in fhicl file "
        << fParams.ArapucaVUVEff << " or " << fParams.ArapucaVISEff << " or "
        << fParams.XArapucaVUVEff << " or " << fParams.XArapucaVISEff
        << " seems to be too large!\n"
        << "Final QE must be equal to or smaller than the scintillation "
        << "pre scale applied at simulation time.\n"
        << "Please check this number (ScintPreScale): "
        << fParams.larProp->ScintPreScale();

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fParams.ArapucaDataFile, fname);
    TFile* file = TFile::Open(fname.c_str());
    //Note: TPB time now implemented at digitization module for both coated pmts and (x)arapucas
    //OpDetSim/digi_arapuca_sbnd.root updated in sbnd_data (now including the TPB times vector)

    // TPB emission time histogram for visible (x)arapucas
    std::vector<double>* timeTPB_p;
    file->GetObject("timeTPB", timeTPB_p);
    fTimeTPB = std::make_unique<CLHEP::RandGeneral>
      (*fEngine, timeTPB_p->data(), timeTPB_p->size());

    std::vector<double>* TimeArapucaVUV_p;
    file->GetObject("TimeArapucaVUV", TimeArapucaVUV_p);
    fTimeArapucaVUV = std::make_unique<CLHEP::RandGeneral>
      (*fEngine, TimeArapucaVUV_p->data(), TimeArapucaVUV_p->size());
    std::vector<double>* TimeArapucaVIS_p;
    file->GetObject("TimeArapucaVIS", TimeArapucaVIS_p);
    fTimeArapucaVIS = std::make_unique<CLHEP::RandGeneral>
      (*fEngine, TimeArapucaVIS_p->data(), TimeArapucaVIS_p->size());
    std::vector<double>* TimeXArapucaVUV_p;
    file->GetObject("TimeXArapucaVUV", TimeXArapucaVUV_p);
    fTimeXArapucaVUV = std::make_unique<CLHEP::RandGeneral>
      (*fEngine, TimeXArapucaVUV_p->data(), TimeXArapucaVUV_p->size());

    fSampling = fSampling / 1000; //in GHz to cancel with ns
    pulsesize = fParams.PulseLength * fSampling;
    wsp.resize(pulsesize);

    Pulse1PE(wsp);
    saturation = fParams.Baseline + fParams.Saturation * fParams.ADC * fParams.MeanAmplitude;
  } // end constructor

  DigiArapucaSBNDAlg::~DigiArapucaSBNDAlg() {}


  void DigiArapucaSBNDAlg::ConstructWaveform(
    int ch,
    sim::SimPhotons const& simphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    double start_time,
    unsigned n_samples)
  {
    std::vector<double> waves(n_samples, fParams.Baseline);
    CreatePDWaveform(simphotons, start_time, waves, pdtype);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiArapucaSBNDAlg::ConstructWaveformLite(
    int ch,
    sim::SimPhotonsLite const& litesimphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    double start_time,
    unsigned n_samples)
  {
    std::vector<double> waves(n_samples, fParams.Baseline);
    std::map<int, int> const& photonMap = litesimphotons.DetectedPhotons;
    CreatePDWaveformLite(photonMap, start_time, waves, pdtype);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiArapucaSBNDAlg::CreatePDWaveform(
    sim::SimPhotons const& simphotons,
    double t_min,
    std::vector<double>& wave,
    std::string pdtype)
  {
    int nCT = 1;
    double tphoton;
    size_t timeBin;
    if(pdtype == "arapuca_vuv") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaVUVEff) { //Sample a random subset according to Arapuca's efficiency
          tphoton = (fTimeArapucaVUV->fire());
          tphoton += simphotons[i].Time - t_min;
          if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          timeBin = std::floor(tphoton * fSampling);
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else if(pdtype == "arapuca_vis") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaVISEff) { //Sample a random subset according to Arapuca's efficiency.
          tphoton = (fTimeArapucaVIS->fire());
          tphoton += simphotons[i].Time - t_min;
          tphoton +=fTimeTPB->fire();
          if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          timeBin = std::floor(tphoton * fSampling);
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else if(pdtype == "xarapuca_vuv") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fXArapucaVUVEff) {
          tphoton = (fTimeXArapucaVUV->fire());
          tphoton += simphotons[i].Time - t_min;
          if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          timeBin = std::floor(tphoton * fSampling);
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else if(pdtype == "xarapuca_vis") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fXArapucaVISEff) {
          tphoton = (CLHEP::RandExponential::shoot(fEngine, 8.5)); //decay time of EJ280 in ns
          tphoton += simphotons[i].Time - t_min;
          tphoton +=fTimeTPB->fire();
          if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          timeBin = std::floor(tphoton * fSampling);
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else{
      throw cet::exception("DigiARAPUCASBNDAlg") << "Wrong pdtype: " << pdtype << std::endl;
    }
    if(fParams.BaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  void DigiArapucaSBNDAlg::CreatePDWaveformLite(
    std::map<int, int> const& photonMap,
    double t_min,
    std::vector<double>& wave,
    std::string pdtype)
  {
    if(pdtype == "xarapuca_vuv"){
      SinglePDWaveformCreatorLite(fXArapucaVUVEff, fTimeXArapucaVUV, wave, photonMap, t_min);
    }
    else if(pdtype == "xarapuca_vis"){
      // creating the waveforms for xarapuca_vis is different than the rest
      // so there's an overload for that which lacks the timeHisto
      SinglePDWaveformCreatorLite(fXArapucaVISEff, wave, photonMap, t_min);
    }
    else if(pdtype == "arapuca_vuv"){
      SinglePDWaveformCreatorLite(fArapucaVUVEff, fTimeArapucaVUV, wave, photonMap, t_min);
    }
    else if(pdtype == "arapuca_vis"){
      SinglePDWaveformCreatorLite(fArapucaVISEff, fTimeArapucaVIS, wave, photonMap, t_min);
    }
    else{
      throw cet::exception("DigiARAPUCASBNDAlg") << "Wrong pdtype: " << pdtype << std::endl;
    }
    if(fParams.BaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  void DigiArapucaSBNDAlg::SinglePDWaveformCreatorLite(
    double effT,
    std::unique_ptr<CLHEP::RandGeneral>& timeHisto,
    std::vector<double>& wave,
    std::map<int, int> const& photonMap,
    double const& t_min
    )
  {
    // TODO: check that this new approach of not using the last
    // (1-accepted_photons) doesn't introduce some bias
    double meanPhotons;
    size_t acceptedPhotons;
    double tphoton;
    int nCT;
    size_t timeBin;
    for (auto const& photonMember : photonMap) {
      // TODO: check that this new approach of not using the last
      // (1-accepted_photons) doesn't introduce some bias
      meanPhotons = photonMember.second*effT;
      acceptedPhotons = CLHEP::RandPoissonQ::shoot(fEngine, meanPhotons);
      for(size_t i = 0; i < acceptedPhotons; i++) {
        tphoton = timeHisto->fire();
        tphoton += photonMember.first - t_min;
        if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
        if(fParams.CrossTalk > 0.0 &&
           (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
        else nCT = 1;
        timeBin = std::floor(tphoton * fSampling);
        if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
      }
    }
  }


  void DigiArapucaSBNDAlg::SinglePDWaveformCreatorLite(
    double effT,
    std::vector<double>& wave,
    std::map<int, int> const& photonMap,
    double const& t_min
    )
  {
    double meanPhotons;
    size_t acceptedPhotons;
    double tphoton;
    int nCT;
    size_t timeBin;
    for (auto const& photonMember : photonMap) {
      // TODO: check that this new approach of not using the last
      // (1-accepted_photons) doesn't introduce some bias
      meanPhotons = photonMember.second*effT;
      acceptedPhotons = CLHEP::RandPoissonQ::shoot(fEngine, meanPhotons);
      for(size_t i = 0; i < acceptedPhotons; i++) {
        tphoton = (CLHEP::RandExponential::shoot(fEngine, fParams.DecayTXArapucaVIS));
        tphoton += photonMember.first - t_min;
        if(tphoton < 0.) continue; // discard if it didn't made it to the acquisition
        if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
        else nCT = 1;
        timeBin = std::floor(tphoton * fSampling);
        if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
      }
    }
  }

  void DigiArapucaSBNDAlg::Pulse1PE(std::vector<double>& wsp)//single pulse waveform
  {
    double time;
    double constT1 = fParams.ADC * fParams.MeanAmplitude;
    for(size_t i = 0; i<wsp.size(); i++) {
      time = static_cast<double>(i) / fSampling;
      if (time < fParams.PeakTime)
        wsp[i] = constT1 * std::exp((time - fParams.PeakTime) / fParams.RiseTime);
      else
        wsp[i] = constT1 * std::exp(-(time - fParams.PeakTime) / fParams.FallTime);
    }
  }


  void DigiArapucaSBNDAlg::AddSPE(
    size_t time_bin,
    std::vector<double>& wave,
    int nphotons) //adding single pulse
  {
    // if(time_bin > wave.size()) return;
    size_t max = time_bin + pulsesize < wave.size() ? time_bin + pulsesize : wave.size();
    auto min_it = std::next(wave.begin(), time_bin);
    auto max_it = std::next(wave.begin(), max);
    std::transform(min_it, max_it,
                   wsp.begin(), min_it,
                   [nphotons](double w, double ws) -> double {
                     return w + ws*nphotons  ; });
  }


  void DigiArapucaSBNDAlg::CreateSaturation(std::vector<double>& wave)
  {
    std::replace_if(wave.begin(), wave.end(),
                    [&](auto w){return w > saturation;}, saturation);
  }


  void DigiArapucaSBNDAlg::AddLineNoise(std::vector< double >& wave)
  {
    // TODO: after running the profiler I can see that this is where
    // most cycles are being used.  Potentially some improvement could
    // be achieved with the below code but the array dynamic allocation
    // is not helping. The function would need to have it passed.
    // ~icaza
    //
    // double *array = new double[wave.size()]();
    // CLHEP::RandGaussQ::shootArray(fEngine, wave.size(), array, 0, fParams.BaselineRMS);
    // for(size_t i = 0; i<wave.size(); i++) {
    //   wave[i] += array[i];
    // }
    // delete array;
    //
    std::transform(wave.begin(), wave.end(), wave.begin(),
                   [this](double w) -> double {
                     return w + CLHEP::RandGaussQ::shoot(fEngine, 0, fParams.BaselineRMS) ; });
  }


  void DigiArapucaSBNDAlg::AddDarkNoise(std::vector<double>& wave)
  {
    int nCT;
    size_t timeBin;
    // Multiply by 10^9 since fDarkNoiseRate is in Hz (conversion from s to ns)
    double mean = 1000000000.0 / fParams.DarkNoiseRate;
    double darkNoiseTime = CLHEP::RandExponential::shoot(fEngine, mean);
    while(darkNoiseTime < wave.size()) {
      timeBin = std::round(darkNoiseTime);
      if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
      else nCT = 1;
      if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
      // Find next time to add dark noise
      darkNoiseTime += CLHEP::RandExponential::shoot(fEngine, mean);
    }
  }


  double DigiArapucaSBNDAlg::FindMinimumTime(sim::SimPhotons const& simphotons)
  {
    double t_min = 1e15;
    for(size_t i = 0; i < simphotons.size(); i++) {
      if(simphotons[i].Time < t_min) t_min = simphotons[i].Time;
    }
    return t_min;
  }


  double DigiArapucaSBNDAlg::FindMinimumTimeLite(std::map<int, int> const& photonMap)
  {
    for (auto const& mapMember : photonMap) {
      if(mapMember.second != 0) return (double)mapMember.first;
    }
    return 1e5;
  }


  // -----------------------------------------------------------------------------
  // // ---  opdet::DigiArapucaSBNDAlgMaker
  // // -----------------------------------------------------------------------------
  DigiArapucaSBNDAlgMaker::DigiArapucaSBNDAlgMaker
  (Config const& config)
  {
    // settings
    fBaseConfig.ADC               = config.voltageToADC();
    fBaseConfig.Baseline          = config.baseline();
    fBaseConfig.Saturation        = config.saturation();
    fBaseConfig.ArapucaVUVEff     = config.arapucaVUVEff();
    fBaseConfig.ArapucaVISEff     = config.arapucaVISEff();
    fBaseConfig.XArapucaVUVEff    = config.xArapucaVUVEff();
    fBaseConfig.XArapucaVISEff    = config.xArapucaVISEff();
    fBaseConfig.RiseTime          = config.riseTime();
    fBaseConfig.FallTime          = config.fallTime();
    fBaseConfig.MeanAmplitude     = config.meanAmplitude();
    fBaseConfig.DarkNoiseRate     = config.darkNoiseRate();
    fBaseConfig.BaselineRMS       = config.baselineRMS();
    fBaseConfig.CrossTalk         = config.crossTalk();
    fBaseConfig.PulseLength       = config.pulseLength();
    fBaseConfig.PeakTime          = config.peakTime();
    fBaseConfig.DecayTXArapucaVIS = config.decayTXArapucaVIS();
    fBaseConfig.ArapucaDataFile   = config.arapucaDataFile();

  }

  std::unique_ptr<DigiArapucaSBNDAlg> DigiArapucaSBNDAlgMaker::operator()(
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

    return std::make_unique<DigiArapucaSBNDAlg>(params);
  } // DigiArapucaSBNDAlgMaker::create()

}
