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

    if(fArapucaEffT1 > 1.0001 || fArapucaEffT2 > 1.0001 || fArapucaEffxT1 > 1.0001 || fArapucaEffxT2 > 1.0001)
      std::cout << "WARNING: Quantum efficiency set in fhicl file " << fParams.ArapucaEffT1 << " or " << fParams.ArapucaEffT2 << " or " << fParams.ArapucaEffxT1 << " or " << fParams.ArapucaEffxT2 << " seems to be too large! Final QE must be equal to or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << fParams.larProp->ScintPreScale() << std::endl;

    xvec.reserve(2000); //TODO: no hardcoded values

    xvec=AssignVector(fParams.ArapucaT1File);
    TimeArapucaT1 = new TH1D("Time Profile T1", "", xvec.size(), 0.0, 1.0*xvec.size());//histogram that stores the arrival time of photons at SiPM (t=0 is the time is reaches the outside of the optical window) for arapuca T1
    for(size_t i=1; i<=xvec.size(); i++)TimeArapucaT1->SetBinContent(i,xvec[i-1]);

    xvec=AssignVector(fParams.ArapucaT2File);
    TimeArapucaT2 = new TH1D("Time Profile T2", "", xvec.size(), 0.0, 1.0*xvec.size());//histogram that stores the arrival time of photons at SiPM (t=0 is the time is reaches the outside of the optical window) for arapuca T2
    for(size_t i=1; i<=xvec.size(); i++)TimeArapucaT2->SetBinContent(i,xvec[i-1]);

    xvec=AssignVector(fParams.ArapucaXT1File);
    TimeArapucaX = new TH1D("Time Profile X-Arapuca T1", "", xvec.size(), 0.0, 1.0*xvec.size());//histogram that stores the arrival time of photons at SiPM (t=0 is the time is reaches the outside of the optical window) for X-arapuca
    for(size_t i=1; i<=xvec.size(); i++)TimeArapucaX->SetBinContent(i,xvec[i-1]);

    fSampling=fSampling/1000; //in GHz to cancel with ns
    pulsesize=fParams.PulseLength*fSampling;
    wsp.resize(pulsesize);

    for(int i = 0; i < pulsesize; i++)
      wsp[i] = (Pulse1PE(static_cast< double >(i) / fSampling));
    //Random number engine initialization
    //int seed = time(NULL);
    //gRandom = new TRandom3(seed);
  }

  DigiArapucaSBNDAlg::~DigiArapucaSBNDAlg()
  {
  }

  void DigiArapucaSBNDAlg::ConstructWaveform(int ch, sim::SimPhotons const& simphotons, std::vector<short unsigned int>& waveform, std::string pdName, double start_time, unsigned n_samples)
  {

    std::vector<double> waves(std::vector<double>(n_samples, fParams.Baseline));
    CreatePDWaveform(simphotons, start_time, waves, pdName);
    waveform.resize(n_samples);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }

  void DigiArapucaSBNDAlg::ConstructWaveformLite(int ch, sim::SimPhotonsLite const& litesimphotons, std::vector<short unsigned int>& waveform, std::string pdName, double start_time, unsigned n_samples)
  {

    std::vector<double> waves(std::vector<double>(n_samples, fParams.Baseline));
    std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
    CreatePDWaveformLite(photonMap, start_time, waves, pdName);
    waveform.resize(n_samples);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }

  void DigiArapucaSBNDAlg::AddSPE(size_t time_bin, std::vector<double>& wave, int nphotons) //adding single pulse
  {
    size_t min = 0, max = 0;

    if(time_bin < wave.size()) {
      min = time_bin;
      max = time_bin + pulsesize < wave.size() ? time_bin + pulsesize : wave.size();
      for(size_t i = min; i < max; i++) {
        wave[i] += (wsp[i - min]) * (double)nphotons;
      }
    }
  }

  double DigiArapucaSBNDAlg::Pulse1PE(double time) const//single pulse waveform
  {
    if (time < fParams.PeakTime) return (fParams.ADC * fParams.MeanAmplitude * std::exp((time - fParams.PeakTime) / fParams.RiseTime));
    else return (fParams.ADC * fParams.MeanAmplitude * std::exp(-(time - fParams.PeakTime) / fParams.FallTime));
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
      AddSPE(timeBin, wave, nCT);
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

  void DigiArapucaSBNDAlg::CreatePDWaveform(sim::SimPhotons const& simphotons, double t_min, std::vector<double>& wave, std::string pdtype)
  {
    int nCT = 1;
    double tphoton = 0;
    if(pdtype == "arapucaT1") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        //if((gRandom->Uniform(1.0))<fArapucaEffT1){ //Sample a random subset according to Arapuca's efficiency
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffT1) { //Sample a random subset according to Arapuca's efficiency
          tphoton = simphotons[i].Time;
          tphoton += (TimeArapucaT1->GetRandom());
          tphoton -= t_min;
          //if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
      }
    }
    if(pdtype == "arapucaT2") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        //if((gRandom->Uniform(1.0))<fArapucaEffT2){ //Sample a random subset according to Arapuca's efficiency.
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffT2) { //Sample a random subset according to Arapuca's efficiency.
          tphoton = simphotons[i].Time;
          tphoton += (TimeArapucaT2->GetRandom());
          tphoton -= t_min;
          //if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
      }
    }
    if(pdtype == "xarapucaT1") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        //if((gRandom->Uniform(1.0))<fArapucaEffx){
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffxT1) {
          tphoton = simphotons[i].Time;
          tphoton += (TimeArapucaX->GetRandom());
          tphoton -= t_min;
          //if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
      }
    }
    if(pdtype == "xarapucaT2") {
      for(size_t i = 0; i < simphotons.size(); i++) {
        //if((gRandom->Uniform(1.0))<fArapucaEffx){
        if((CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffxT2) {
          tphoton = simphotons[i].Time;
          tphoton += (CLHEP::RandExponential::shoot(fEngine, 8.5)); //decay time of EJ280 in ns
          tphoton -= t_min;
          //if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
      }
    }

    if(fParams.BaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }

  void DigiArapucaSBNDAlg::CreateSaturation(std::vector<double>& wave)  //Implementing saturation effects
  {
    for(size_t k = 0; k < wave.size(); k++) {
      if(wave[k] > (fParams.Baseline + fParams.Saturation * fParams.ADC * fParams.MeanAmplitude))
        wave[k] = fParams.Baseline + fParams.Saturation * fParams.ADC * fParams.MeanAmplitude;
    }
  }

  void DigiArapucaSBNDAlg::CreatePDWaveformLite(std::map< int, int > const& photonMap, double t_min, std::vector<double>& wave, std::string pdtype)
  {
    double tphoton = 0;
    int nCT = 1;
    for (auto const& mapMember : photonMap) {
      for(int i = 0; i < mapMember.second; i++) {
        if(pdtype == "arapucaT1" && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffT1) {
          tphoton = (TimeArapucaT1->GetRandom());
          tphoton += mapMember.first - t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
        if(pdtype == "arapucaT2" && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffT2) {
          tphoton = (TimeArapucaT2->GetRandom());
          tphoton += mapMember.first - t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
        if(pdtype == "xarapucaT1" && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffxT1) {
          tphoton = (TimeArapucaX->GetRandom());
          tphoton += mapMember.first - t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
        if(pdtype == "xarapucaT2" && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fArapucaEffxT2) {
          tphoton = (CLHEP::RandExponential::shoot(fEngine, 8.5)); //decay time of EJ280 in ns
          tphoton += mapMember.first - t_min;
          if(fParams.CrossTalk > 0.0 && (CLHEP::RandFlat::shoot(fEngine, 1.0)) < fParams.CrossTalk) nCT = 2;
          else nCT = 1;
          AddSPE(tphoton * fSampling, wave, nCT);
        }
      }
    }
    if(fParams.BaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }

  double DigiArapucaSBNDAlg::FindMinimumTimeLite(std::map< int, int > const& photonMap)
  {
    for (auto const& mapMember : photonMap) {
      if(mapMember.second != 0) return (double)mapMember.first;
    }
    return 1e5;
  }

  std::vector<double> DigiArapucaSBNDAlg::AssignVector(std::string filename){

    cet::search_path sp("FW_SEARCH_PATH");
    std::string const& fname = sp.find_file(filename);

     // Read in vector from text file
    std::ifstream file (fname);
    std::string line;
 
    mf::LogInfo("DigiArapucaSBNDAlg")<<"DigiArapuca opening file " << filename.c_str();

    // Read in each line and place into vector
    std::vector<double> vec;
    if (file.is_open()){
      while ( file.good() ){
        getline (file, line);
        vec.push_back( strtod( line.c_str(), NULL ) );
      }
    }
   else throw cet::exception("DigiArapucaSBNDAlg") << "No File: Unable to open file " << fname << "\n";

   file.close();
   return(vec);
  
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
    fBaseConfig.ArapucaT1File     = config.arapucaT1File();
    fBaseConfig.ArapucaT2File     = config.arapucaT2File();
    fBaseConfig.ArapucaXT1File    = config.arapucaXT1File();
    fBaseConfig.RiseTime          = config.riseTime();
    fBaseConfig.FallTime          = config.fallTime();
    fBaseConfig.MeanAmplitude     = config.meanAmplitude();
    fBaseConfig.DarkNoiseRate     = config.darkNoiseRate();
    fBaseConfig.BaselineRMS       = config.baselineRMS();
    fBaseConfig.CrossTalk         = config.crossTalk();
    fBaseConfig.PulseLength       = config.pulseLength();
    fBaseConfig.PeakTime          = config.peakTime();

  }

  std::unique_ptr<DigiArapucaSBNDAlg>
  DigiArapucaSBNDAlgMaker::operator()(
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
