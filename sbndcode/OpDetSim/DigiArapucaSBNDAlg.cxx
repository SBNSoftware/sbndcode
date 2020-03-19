#include "sbndcode/OpDetSim/DigiArapucaSBNDAlg.h"

#ifndef SBND_OPDETSIM_DIGIARAPUCASBNDALG_CXX
#define SBND_OPDETSIM_DIGIARAPUCASBNDALG_CXX

//------------------------------------------------------------------------------
//--- opdet::simarapucasbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet {

  //DigiArapucaSBNDAlg::DigiArapucaSBNDAlg(fhicl::ParameterSet const& p)
  DigiArapucaSBNDAlg::DigiArapucaSBNDAlg(ConfigurationParameters_t const& config)
    : fParams(config)
    , fSampling(fParams.timeService->OpticalClock().Frequency())
    , fArapucaEffT1(fParams.ArapucaEffT1 / fParams.larProp->ScintPreScale())
    , fArapucaEffT2(fParams.ArapucaEffT2 / fParams.larProp->ScintPreScale())
    , fArapucaEffxT1(fParams.ArapucaEffxT1 / fParams.larProp->ScintPreScale())
    , fArapucaEffxT2(fParams.ArapucaEffxT2 / fParams.larProp->ScintPreScale())
    , fEngine(fParams.engine)
  {

    //art::ServiceHandle<rndm::NuRandomService> seedSvc;
    //fEngine = new CLHEP::HepJamesRandom;
    //seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(fEngine), "DigiArapucaSBNDAlg");

    //    std::cout << "arapucas corrected efficiencies = " << fArapucaEffT1 << ", " << fArapucaEffT2 << " and " << fArapucaEffx << std::endl;
    // std::cout << "optical clock = " << fSampling << std::endl;

    if(fArapucaEffT1 > 1.0001 || fArapucaEffT2 > 1.0001 || fArapucaEffxT1 > 1.0001 || fArapucaEffxT2 > 1.0001)
      std::cout << "WARNING: Quantum efficiency set in fhicl file " << fParams.ArapucaEffT1 << " or " << fParams.ArapucaEffT2 << " or " << fParams.ArapucaEffxT1 << " or " << fParams.ArapucaEffxT2 << " seems to be too large! Final QE must be equal to or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << fParams.larProp->ScintPreScale() << std::endl;

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fParams.ArapucaDataFile, fname);
    TFile* file = TFile::Open(fname.c_str());
    file->GetObject("TimeArapucaT1", TimeArapucaT1);
    file->GetObject("TimeArapucaT2", TimeArapucaT2);
    file->GetObject("TimeArapucaX", TimeArapucaX);

    fSampling = fSampling / 1000; //in GHz to cancel with ns
    pulsesize = fParams.PulseLength * fSampling;
    wsp.resize(pulsesize);

    for(int i = 0; i < pulsesize; i++)
      wsp[i] = (Pulse1PE(static_cast< double >(i) / fSampling));

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
    std::vector<double> waves(std::vector<double>(n_samples, fParams.Baseline));
    CreatePDWaveform(simphotons, start_time, waves, pdtype);
    waveform.resize(n_samples);
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
    std::vector<double> waves(std::vector<double>(n_samples, fParams.Baseline));
    std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
    CreatePDWaveformLite(photonMap, start_time, waves, pdtype);
    waveform.resize(n_samples);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  void DigiArapucaSBNDAlg::CreatePDWaveform(
    sim::SimPhotons const& simphotons,
    double t_min,
    std::vector<double>& wave,
    std::string pdtype)
  {
    int nCT = 1;
    double tphoton = 0;
    if(pdtype == "arapuca_vuv") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffT1) { //Sample a random subset according to Arapuca's efficiency
          tphoton = simphotons[i].Time;
          tphoton += (TimeArapucaT1->GetRandom());
          tphoton -= t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          double timeBin = tphoton * fSampling;
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else if(pdtype == "arapuca_vis") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffT2) { //Sample a random subset according to Arapuca's efficiency.
          tphoton = simphotons[i].Time;
          tphoton += (TimeArapucaT2->GetRandom());
          tphoton -= t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          double timeBin = tphoton * fSampling;
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else if(pdtype == "xarapuca_vuv") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffxT1) {
          tphoton = simphotons[i].Time;
          tphoton += (TimeArapucaX->GetRandom());
          tphoton -= t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          double timeBin = tphoton * fSampling;
          if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
        }
      }
    }
    else if(pdtype == "xarapuca_vis") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffxT2) {
          tphoton = simphotons[i].Time;
          tphoton += (CLHEP::RandExponential::shoot(fEngine, 8.5)); //decay time of EJ280 in ns
          tphoton -= t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          double timeBin = tphoton * fSampling;
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
    double tphoton = 0;
    int nCT = 1;
    for (auto const& mapMember : photonMap) {
      for(int i = 0; i < mapMember.second; i++) {
        double randFlat = CLHEP::RandFlat::shoot(fEngine, 1.0);
        if(pdtype == "arapuca_vuv"){
          if(randFlat < fArapucaEffT1) {
            tphoton = (TimeArapucaT1->GetRandom());
            tphoton += mapMember.first - t_min;
            if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
            else nCT = 1;
            double timeBin = tphoton * fSampling;
            if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
          }
        }
        else if(pdtype == "arapuca_vis"){
          if(randFlat < fArapucaEffT2) {
            tphoton = (TimeArapucaT2->GetRandom());
            tphoton += mapMember.first - t_min;
            if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
            else nCT = 1;
            double timeBin = tphoton * fSampling;
            if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
          }
        }
        else if(pdtype == "xarapuca_vuv"){
          if(randFlat < fArapucaEffxT1) {
            tphoton = (TimeArapucaX->GetRandom());
            tphoton += mapMember.first - t_min;
            if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
            else nCT = 1;
            double timeBin = tphoton * fSampling;
            if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
          }
        }
        else if(pdtype == "xarapuca_vis"){
          if( randFlat < fArapucaEffxT2) {
            tphoton = (CLHEP::RandExponential::shoot(fEngine, 8.5)); //decay time of EJ280 in ns
            tphoton += mapMember.first - t_min;
            if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
            else nCT = 1;
            double timeBin = tphoton * fSampling;
            if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
          }
        }
        else{
          throw cet::exception("DigiARAPUCASBNDAlg") << "Wrong pdtype: " << pdtype << std::endl;
        }
      }
    }
    if(fParams.BaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  double DigiArapucaSBNDAlg::Pulse1PE(double time) const//single pulse waveform
  {
    if (time < fParams.PeakTime) return (fParams.ADC * fParams.MeanAmplitude * std::exp((time - fParams.PeakTime) / fParams.RiseTime));
    else return (fParams.ADC * fParams.MeanAmplitude * std::exp(-(time - fParams.PeakTime) / fParams.FallTime));
  }


  void DigiArapucaSBNDAlg::AddSPE(
    size_t time_bin,
    std::vector<double>& wave,
    int nphotons) //adding single pulse
  {
    size_t max = time_bin + pulsesize < wave.size() ? time_bin + pulsesize : wave.size();
    auto min_it = std::next(wave.begin(), time_bin);
    auto max_it = std::next(wave.begin(), max);
    std::transform(min_it, max_it,
                   wsp.begin(), min_it,
                   [nphotons](double w, double ws) -> double {
                     return w + ws*nphotons  ; });
  }


  void DigiArapucaSBNDAlg::CreateSaturation(std::vector<double>& wave)  //Implementing saturation effects
  {
    for(size_t k = 0; k < wave.size(); k++) {
      if(wave[k] > saturation) wave[k] = saturation;
    }
  }


  void DigiArapucaSBNDAlg::AddLineNoise(std::vector< double >& wave)
  {
    double noise;
    for(size_t i = 0; i < wave.size(); i++) {
      //noise= gRandom->Gaus(0, fParams.BaselineRMS); //gaussian baseline noise
      noise = CLHEP::RandGauss::shoot(fEngine, 0, fParams.BaselineRMS); //gaussian baseline noise
      wave[i] += noise;
    }
  }

  void DigiArapucaSBNDAlg::AddDarkNoise(std::vector< double >& wave)
  {
    int nCT;
    // Multiply by 10^9 since fDarkNoiseRate is in Hz (conversion from s to ns)
    //double darkNoiseTime = static_cast< double >(gRandom->Exp((1.0/fParams.DarkNoiseRate)*1000000000.0));
    double darkNoiseTime = CLHEP::RandExponential::shoot(fEngine, (1.0 / fParams.DarkNoiseRate) * 1000000000.0);
    while (darkNoiseTime < wave.size()) {
      size_t timeBin = (darkNoiseTime);
      //if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
      if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
      else nCT = 1;
      if(timeBin < wave.size()) AddSPE(timeBin, wave, nCT);
      // Find next time to add dark noise
      //darkNoiseTime += static_cast< double >(gRandom->Exp((1.0/fParams.DarkNoiseRate)*1000000000.0));
      darkNoiseTime += CLHEP::RandExponential::shoot(fEngine, (1.0 / fParams.DarkNoiseRate) * 1000000000.0);
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
    fBaseConfig.ArapucaEffT1      = config.arapucaEffT1();
    fBaseConfig.ArapucaEffT2      = config.arapucaEffT2();
    fBaseConfig.ArapucaEffxT1     = config.arapucaEffxT1();
    fBaseConfig.ArapucaEffxT2     = config.arapucaEffxT2();
    fBaseConfig.RiseTime          = config.riseTime();
    fBaseConfig.FallTime          = config.fallTime();
    fBaseConfig.MeanAmplitude     = config.meanAmplitude();
    fBaseConfig.DarkNoiseRate     = config.darkNoiseRate();
    fBaseConfig.BaselineRMS       = config.baselineRMS();
    fBaseConfig.CrossTalk         = config.crossTalk();
    fBaseConfig.PulseLength       = config.pulseLength();
    fBaseConfig.PeakTime          = config.peakTime();
    fBaseConfig.ArapucaDataFile   = config.arapucaDataFile();

  }

  std::unique_ptr<DigiArapucaSBNDAlg> DigiArapucaSBNDAlgMaker::operator()(
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

    return std::make_unique<DigiArapucaSBNDAlg>(params);
  } // DigiArapucaSBNDAlgMaker::create()

}

#endif //SBND_OPDETSIM_DIGIARAPUCASBNDALG_CXX
