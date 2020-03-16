#include "sbndcode/OpDetSim/DigiPMTSBNDAlg.h"

#ifndef DIGIPMTSBNDALG_CXX
#define DIGIPMTSBNDALG_CXX

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

    std::cout << "PMT corrected efficiencies = " << fQEDirect << " " << fQERefl << std::endl;

    if(fQERefl > 1.0001 || fQEDirect > 1.0001)
      std::cout << "WARNING: Quantum efficiency set in fhicl file " << fParams.QERefl << " or " << fParams.QEDirect << " seems to be too large! Final QE must be equal or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << fParams.larProp->ScintPreScale() << std::endl;

    fSampling = fSampling / 1000.0; //in GHz, to cancel with ns

    //Random number engine initialization
    //int seed = time(NULL);
    //gRandom = new TRandom3(seed);

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fParams.PMTDataFile, fname);
    TFile* file = TFile::Open(fname.c_str());
    file->GetObject("timeTPB", timeTPB);

    //shape of single pulse
    if (fParams.SinglePEmodel) {
      std::cout << " using testbench pe response " << std::endl;
      std::vector<double> *wsp_pointer;
      file->GetObject("wsp", wsp_pointer);
      wsp = *wsp_pointer;
      pulsesize = wsp.size();
    }
    else {
      std::cout << " using ideal pe response " << std::endl;
      //shape of single pulse
      sigma1 = fParams.PMTRiseTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      sigma2 = fParams.PMTFallTime / (std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));

      pulsesize = (int)((6 * sigma2 + fParams.TransitTime) * fSampling);
      wsp.resize(pulsesize);
      for(int i = 0; i < pulsesize; i++) {
        wsp[i] = (Pulse1PE(static_cast< double >(i) / fSampling));
      }
    }
  } // end constructor

  DigiPMTSBNDAlg::~DigiPMTSBNDAlg()
  { }

  void DigiPMTSBNDAlg::ConstructWaveform(
    int ch,
    sim::SimPhotons const& simphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    std::map<int, sim::SimPhotons>& auxmap,
    double start_time,
    unsigned n_sample)
  {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveform(simphotons, start_time, waves, ch, pdtype, auxmap);
    waveform.resize(n_sample);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }

  void DigiPMTSBNDAlg::ConstructWaveformLite(
    int ch,
    sim::SimPhotonsLite const& litesimphotons,
    std::vector<short unsigned int>& waveform,
    std::string pdtype,
    std::map<int, sim::SimPhotonsLite>& auxmap,
    double start_time,
    unsigned n_sample)
  {

    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformLite(litesimphotons, start_time, waves, ch, pdtype, auxmap);
    waveform.resize(n_sample);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  double DigiPMTSBNDAlg::Pulse1PE(double time)//single pulse waveform
  {
    if (time < fParams.TransitTime) return (fParams.PMTChargeToADC * fParams.PMTMeanAmplitude * std::exp(-1.0 * pow(time - fParams.TransitTime, 2.0) / (2.0 * pow(sigma1, 2.0))));
    else return (fParams.PMTChargeToADC * fParams.PMTMeanAmplitude * std::exp(-1.0 * pow(time - fParams.TransitTime, 2.0) / (2.0 * pow(sigma2, 2.0))));
  }


  double DigiPMTSBNDAlg::Transittimespread(double fwhm)
  {

    double tts, sigma;

    sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
    //tts = gRandom->Gaus(0,sigma);
    tts = CLHEP::RandGauss::shoot(fEngine, 0, sigma);

    return tts;
  }


  void DigiPMTSBNDAlg::AddSPE(size_t time_bin, std::vector<double>& wave)
  {

    size_t min = 0;
    size_t max = 0;

    if(time_bin < wave.size()) {
      min = time_bin;
      max = time_bin + pulsesize < wave.size() ? time_bin + pulsesize : wave.size();
      for(size_t i = min; i < max; i++) {
        wave[i] += wsp[i - min];
      }
    }
  }


  void DigiPMTSBNDAlg::CreatePDWaveform(
    sim::SimPhotons const& simphotons,
    double t_min,
    std::vector<double>& wave,
    int ch,
    std::string pdtype,
    std::map<int, sim::SimPhotons>& auxmap)
  {

    double ttsTime = 0;
    for(size_t i = 0; i < simphotons.size(); i++) { //simphotons is here reflected light. To be added for all PMTs
      //if((gRandom->Uniform(1.0))<fQERefl){
      if(CLHEP::RandFlat::shoot(fEngine, 1.0) < fQERefl) {
        if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
        AddSPE((fParams.TransitTime + ttsTime + simphotons[i].Time - t_min)*fSampling, wave);
      }
    }
    if(pdtype == "pmt_coated") { //To add direct light for TPB coated PMTs
      sim::SimPhotons auxphotons;
      double ttpb = 0;
      // TODO: no need to loop only to overwrite member. ~icaza
      for (auto& mapMember : auxmap) {
        if(mapMember.first == ch) auxphotons = mapMember.second;
      }
      for(size_t j = 0; j < auxphotons.size(); j++) { //auxphotons is direct light
        //if((gRandom->Uniform(1.0))<fQEDirect){
        if(CLHEP::RandFlat::shoot(fEngine, 1.0) < fQEDirect) {
          if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
          ttpb = timeTPB->GetRandom(); //for including TPB emission time
          AddSPE((fParams.TransitTime + ttsTime + auxphotons[j].Time + ttpb - t_min)*fSampling, wave);
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
    std::map<int, sim::SimPhotonsLite>& auxmap)
  {

    double ttsTime = 0;
    std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
    for (auto const& mapMember : photonMap) { //including reflected light for all PMT channels
      for(int i = 0; i < mapMember.second; i++) {
        //if((gRandom->Uniform(1.0))<(fQERefl)){
        if(CLHEP::RandFlat::shoot(fEngine, 1.0) < (fQERefl)) {
          if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
          AddSPE((fParams.TransitTime + ttsTime + mapMember.first - t_min)*fSampling, wave);
        }
      }
    }
    if(pdtype == "pmt_coated") { //To add direct light for TPB coated PMTs
      double ttpb;
      sim::SimPhotonsLite auxphotons;
      // TODO: no need to loop only to overwrite member. ~icaza
      for (auto& mapMember : auxmap) //auxphotons is direct light
        if(mapMember.first == ch) auxphotons = mapMember.second;
      std::map< int, int > const& auxphotonMap = auxphotons.DetectedPhotons;
      for (auto& mapMember2 : auxphotonMap) {
        for(int i = 0; i < mapMember2.second; i++) {
          //if((gRandom->Uniform(1.0))<(fQEDirect)){
          if(CLHEP::RandFlat::shoot(fEngine, 1.0) < (fQEDirect)) {
            if(fParams.TTS > 0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
            ttpb = timeTPB->GetRandom(); //for including TPB emission time
            AddSPE((fParams.TransitTime + ttsTime + mapMember2.first + ttpb - t_min)*fSampling, wave);
          }
        }
      }
    }

    if(fParams.PMTBaselineRMS > 0.0) AddLineNoise(wave);
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::CreateSaturation(std::vector<double>& wave)  //Implementing saturation effects
  {

    for(size_t k = 0; k < wave.size(); k++) {
      if(wave[k] < (fParams.PMTBaseline + fParams.PMTSaturation * fParams.PMTChargeToADC * fParams.PMTMeanAmplitude))
        wave[k] = fParams.PMTBaseline + fParams.PMTSaturation * fParams.PMTChargeToADC * fParams.PMTMeanAmplitude;
    }
  }

  void DigiPMTSBNDAlg::AddLineNoise(std::vector<double>& wave)
  {

    double noise = 0.0;

    for(size_t i = 0; i < wave.size(); i++) {
      //noise = gRandom->Gaus(0,fParams.PMTBaselineRMS); //gaussian noise
      noise = CLHEP::RandGauss::shoot(fEngine, 0, fParams.PMTBaselineRMS); //gaussian noise
      wave[i] += noise;
    }
  }

  void DigiPMTSBNDAlg::AddDarkNoise(std::vector< double >& wave)
  {
    size_t timeBin = 0;

    // Multiply by 10^9 since fParams.DarkNoiseRate is in Hz (conversion from s to ns)
    //double darkNoiseTime = static_cast< double >(gRandom->Exp((1.0/fParams.PMTDarkNoiseRate)*1000000000.0));
    double darkNoiseTime = CLHEP::RandExponential::shoot(fEngine, (1.0 / fParams.PMTDarkNoiseRate) * 1000000000.0);
    while (darkNoiseTime < wave.size()) {
      timeBin = (darkNoiseTime);
      AddSPE(timeBin, wave);
      // Find next time to add dark noise
      //darkNoiseTime += static_cast< double >(gRandom->Exp((1.0/fParams.PMTDarkNoiseRate)*1000000000.0));
      darkNoiseTime += CLHEP::RandExponential::shoot(fEngine, (1.0 / fParams.PMTDarkNoiseRate) * 1000000000.0);
    }
  }

  double DigiPMTSBNDAlg::FindMinimumTime(
    sim::SimPhotons const& simphotons,
    int ch,
    std::string pdtype,
    std::map<int, sim::SimPhotons>& auxmap)
  {
    double t_min = 1e15;

    if(pdtype == "pmt_uncoated") { //TPB non-coated PMTs
      for(size_t i = 0; i < simphotons.size(); i++) {
        if(simphotons[i].Time < t_min) t_min = simphotons[i].Time;
      }
    }
    else {//for coated PMTs
      sim::SimPhotons auxphotons;
      // TODO: no need to loop only to overwrite member. ~icaza
      for (auto& mapMember : auxmap) {
        if(mapMember.first == ch) auxphotons = mapMember.second;
      }
      auxphotons += (simphotons);
      for(size_t i = 0; i < auxphotons.size(); i++) {
        if(auxphotons[i].Time < t_min) t_min = auxphotons[i].Time;
      }
    }
    return t_min;
  }

  double DigiPMTSBNDAlg::FindMinimumTimeLite(
    sim::SimPhotonsLite const& litesimphotons,
    int ch,
    std::string pdtype,
    std::map<int, sim::SimPhotonsLite>& auxmap)
  {

    if(pdtype == "pmt_uncoated") { //TPB non-coated PMTs
      std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
      for (auto const& mapMember : photonMap) {
        if(mapMember.second != 0) return (double)mapMember.first;
      }
    }
    else {
      sim::SimPhotonsLite auxphotons;
      // TODO: no need to loop only to overwrite member. ~icaza
      for (auto& mapMember : auxmap) {
        if(mapMember.first == ch) auxphotons = mapMember.second;
      }
      auxphotons += (litesimphotons);
      std::map< int, int > const& auxphotonMap = auxphotons.DetectedPhotons;
      for (auto & mapMember : auxphotonMap) {
        if(mapMember.second != 0) return (double)mapMember.first;
      }
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

#endif
