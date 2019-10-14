#include "DigiPMTSBNDAlg.h"

#ifndef DIGIPMTSBNDALG_CXX
#define DIGIPMTSBNDALG_CXX
 
//------------------------------------------------------------------------------
//--- opdet::simpmtsbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet{


  DigiPMTSBNDAlg::DigiPMTSBNDAlg(ConfigurationParameters_t const& config) 
  : fParams(config)
  , fSampling(fParams.timeService->OpticalClock().Frequency()) //This is number is wrong!!! Therefore, the hard coded value!!!
  , fQEDirect(fParams.QEDirect / fParams.larProp->ScintPreScale())
  , fQERefl(fParams.QERefl / fParams.larProp->ScintPreScale())
  {

    std::cout << "PMT corrected efficiencies = " << fQEDirect << " " << fQERefl << std::endl;

    if(fQERefl>1.0001 || fQEDirect>1.0001)
	std::cout << "WARNING: Quantum efficiency set in fhicl file " << fParams.QERefl << " or " << fParams.QEDirect << " seems to be too large! Final QE must be equal or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << fParams.larProp->ScintPreScale() << std::endl;

    fSampling=fSampling/1000.0; //in GHz, to cancel with ns
  
//Random number engine initialization
    int seed = time(NULL);
    gRandom = new TRandom3(seed);

    timeTPB = new TH1D("Time TPB", "", 1000, 0.0, 1000.0);//histogram that stores the emission time of photons converted by TPB

    double x[1000]={12321, 10239, 8303, 6975, 5684, 4667, 4031, 3446, 2791, 2485, 2062, 1724, 1419, 1367, 1111, 982, 974, 822, 732, 653, 665, 511, 500, 452, 411, 439, 409, 357, 342, 357, 302, 296, 316, 271, 286, 265, 260, 288, 279, 238, 214, 242, 232, 238, 251, 239, 200, 225, 182, 190, 206, 194, 188, 227, 210, 198, 170, 184, 158, 160, 170, 183, 168, 143, 158, 140, 167, 145, 154, 162, 155, 115, 143, 148, 124, 126, 133, 122, 91, 130, 90, 124, 135, 112, 94, 81, 107, 99, 109, 78, 83, 75, 68, 97, 69, 74, 91, 84, 84, 74, 68, 73, 71, 55, 68, 40, 55, 63, 71, 62, 63, 60, 71, 55, 62, 53, 54, 58, 63, 39, 42, 56, 44, 33, 36, 43, 60, 49, 50, 51, 52, 49, 47, 57, 39, 45, 41, 23, 41, 26, 29, 51, 23, 45, 26, 50, 39, 20, 44, 27, 14, 17, 13, 35, 20, 25, 26, 26, 29, 31, 20, 17, 28, 24, 28, 34, 22, 16, 17, 21, 23, 33, 15, 30, 8, 20, 15, 20, 14, 17, 18, 21, 16, 20, 22, 24, 14, 18, 25, 13, 10, 13, 11, 18, 9, 4, 13, 23, 10, 13, 15, 26, 21, 18, 15, 17, 6, 15, 9, 13, 14, 6, 13, 9, 9, 6, 8, 7, 13, 13, 11, 13, 8, 5, 8, 13, 7, 9, 6, 14, 11, 11, 9, 10, 13, 9, 4, 3, 17, 3, 5, 1, 5, 5, 5, 15, 4, 6, 3, 11, 3, 10, 8, 8, 7, 5, 8, 7, 13, 7, 7, 15, 5, 6, 9, 9, 7, 4, 9, 7, 5, 7, 5, 5, 6, 3, 8, 6, 4, 12, 7, 4, 4, 6, 7, 9, 3, 2, 3, 4, 4, 1, 9, 9, 2, 2, 2, 4, 3, 3, 1, 5, 1, 7, 4, 6, 4, 6, 7, 4, 4, 5, 2, 3, 2, 8, 4, 9, 4, 4, 8, 2, 2, 2, 0, 2, 14, 4, 3, 2, 3, 4, 5, 3, 7, 1, 4, 1, 1, 8, 3, 5, 2, 1, 7, 4, 5, 0, 5, 6, 4, 2, 6, 1, 4, 5, 0, 0, 4, 1, 4, 6, 2, 0, 4, 3, 4, 3, 3, 8, 4, 1, 2, 3, 2, 6, 7, 4, 2, 5, 6, 3, 2, 6, 5, 3, 1, 4, 6, 3, 0, 2, 2, 1, 0, 0, 5, 4, 3, 3, 3, 9, 0, 4, 2, 6, 0, 2, 6, 4, 6, 1, 0, 5, 3, 1, 1, 4, 0, 1, 1, 2, 2, 4, 5, 7, 5, 3, 7, 6, 3, 2, 1, 3, 0, 4, 4, 1, 2, 4, 6, 11, 7, 5, 5, 5, 4, 2, 5, 2, 2, 3, 0, 6, 3, 2, 3, 3, 8, 0, 0, 1, 2, 1, 0, 3, 6, 1, 6, 1, 4, 5, 0, 2, 6, 0, 3, 7, 0, 2, 5, 2, 6, 3, 5, 2, 2, 1, 5, 5, 0, 3, 3, 2, 3, 6, 0, 0, 1, 0, 1, 4, 4, 2, 3, 4, 3, 7, 1, 1, 3, 2, 2, 2, 2, 4, 9, 4, 8, 2, 2, 6, 5, 2, 6, 1, 2, 6, 7, 0, 5, 0, 4, 0, 1, 4, 1, 2, 2, 1, 0, 2, 4, 1, 0, 3, 3, 0, 6, 2, 0, 0, 3, 0, 2, 2, 3, 3, 2, 0, 2, 1, 3, 1, 2, 1, 1, 2, 4, 3, 0, 2, 4, 2, 3, 3, 3, 5, 5, 2, 2, 1, 0, 2, 0, 1, 0, 0, 5, 1, 3, 8, 4, 3, 2, 6, 4, 1, 3, 2, 0, 9, 4, 2, 7, 2, 0, 0, 2, 1, 3, 4, 2, 3, 3, 3, 2, 8, 6, 3, 1, 3, 3, 0, 0, 3, 0, 6, 1, 0, 2, 0, 0, 1, 2, 7, 2, 1, 0, 1, 6, 3, 2, 0, 1, 0, 2, 3, 5, 3, 6, 4, 1, 1, 0, 0, 7, 1, 1, 1, 1, 8, 1, 1, 0, 0, 0, 1, 2, 2, 7, 1, 1, 0, 1, 1, 3, 1, 1, 1, 3, 5, 2, 1, 0, 3, 2, 5, 0, 5, 4, 2, 5, 3, 3, 0, 0, 5, 0, 5, 1, 4, 0, 1, 6, 1, 6, 1, 2, 1, 2, 4, 0, 8, 3, 1, 7, 1, 2, 4, 4, 2, 3, 5, 0, 4, 5, 2, 1, 1, 5, 2, 0, 4, 2, 0, 2, 4, 4, 4, 4, 5, 0, 3, 0, 2, 3, 3, 0, 0, 6, 1, 1, 6, 10, 0, 2, 0, 1, 4, 1, 0, 1, 3, 2, 0, 1, 1, 0, 0, 8, 4, 0, 0, 3, 0, 0, 4, 0, 1, 4, 0, 1, 0, 2, 1, 5, 2, 2, 0, 0, 0, 2, 1, 0, 0, 4, 0, 3, 1, 1, 2, 1, 1, 2, 0, 1, 3, 3, 0, 2, 2, 3, 2, 0, 1, 2, 0, 0, 1, 0, 4, 2, 3, 0, 0, 1, 2, 2, 6, 2, 3, 4, 2, 3, 1, 2, 7, 3, 2, 3, 3, 1, 0, 0, 4, 2, 9, 0, 3, 2, 2, 0, 1, 1, 1, 8, 1, 4, 2, 0, 4, 0, 2, 0, 0, 0, 3, 0, 0, 4, 0, 2, 2, 5, 1, 3, 6, 0, 0, 0, 0, 0, 3, 2, 1, 0, 1, 1, 1, 1, 2, 0, 1, 4, 4, 1, 1, 0, 0, 9, 2, 1, 1, 2, 0, 2, 2, 2, 2, 0, 1, 7, 0, 7, 0, 5, 0, 5, 1, 1, 0, 1, 1, 1, 0, 1, 0, 2, 3, 2, 1, 1, 0, 0, 2, 4, 0, 0, 2, 4, 2, 5, 2, 2, 1, 1, 4, 0, 1, 2, 1, 3, 3, 0, 2, 1, 3, 0, 2, 2, 1, 0, 4, 5, 3, 0, 0, 2, 1, 3, 0, 2, 0, 3, 2, 2, 0, 3, 1, 0, 5, 2, 2, 4, 6, 3, 2, 2, 2, 1, 4, 6, 1, 2, 1, 2, 6, 1, 2};

  for(int i=1; i<=1000; i++)timeTPB->SetBinContent(i,x[i-1]);

  //shape of single pulse
    sigma1 = fParams.PMTRiseTime/(std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
    sigma2 = fParams.PMTFallTime/(std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));

    pulsesize=(int)((6*sigma2+fParams.TransitTime)*fSampling);
    wsp.resize(pulsesize);

    for(int i=0; i<pulsesize; i++){
	wsp[i]=(Pulse1PE(static_cast< double >(i)/fSampling));
    }
  }

  DigiPMTSBNDAlg::~DigiPMTSBNDAlg()
  { }

  void DigiPMTSBNDAlg::ConstructWaveform(int ch, sim::SimPhotons const& simphotons, std::vector<short unsigned int>& waveform, std::string pdtype, std::map<int,sim::SimPhotons> auxmap, double start_time, unsigned n_sample) {
    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveform(simphotons, start_time, waves, ch, pdtype, auxmap);
    waveform.resize(n_sample);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }

  void DigiPMTSBNDAlg::ConstructWaveformLite(int ch, sim::SimPhotonsLite const& litesimphotons, std::vector<short unsigned int>& waveform, std::string pdtype, std::map<int, sim::SimPhotonsLite> auxmap, double start_time, unsigned n_sample){	

    std::vector<double> waves(n_sample, fParams.PMTBaseline);
    CreatePDWaveformLite(litesimphotons, start_time, waves, ch, pdtype, auxmap);
    waveform.resize(n_sample);
    waveform = std::vector<short unsigned int> (waves.begin(), waves.end());
  }


  double DigiPMTSBNDAlg::Pulse1PE(double time)//single pulse waveform
  { 
    if (time < fParams.TransitTime) return (fParams.PMTChargeToADC*fParams.PMTMeanAmplitude*std::exp(-1.0*pow(time - fParams.TransitTime,2.0)/(2.0*pow(sigma1,2.0))));
    else return (fParams.PMTChargeToADC*fParams.PMTMeanAmplitude*std::exp(-1.0*pow(time - fParams.TransitTime,2.0)/(2.0*pow(sigma2,2.0))));
  }


  double DigiPMTSBNDAlg::Transittimespread(double fwhm){

    double tts,sigma;

    sigma = fwhm/(2.0*sqrt(2.0*log(2.0)));
    tts = gRandom->Gaus(0,sigma);	

    return tts;
  } 


  void DigiPMTSBNDAlg::AddSPE(size_t time_bin, std::vector<double>& wave){

    size_t min=0;
    size_t max=0;

    if(time_bin<wave.size()){
	min=time_bin;
	max=time_bin+pulsesize < wave.size() ? time_bin+pulsesize : wave.size();
	for(size_t i = min; i<= max; i++){
	  wave[i]+= wsp[i-min];	
	}	
    }
  }


  void DigiPMTSBNDAlg::CreatePDWaveform(sim::SimPhotons const& simphotons, double t_min, std::vector<double>& wave, int ch, std::string pdtype, std::map<int,sim::SimPhotons> auxmap){

    double ttsTime=0;
    for(size_t i=0; i<simphotons.size(); i++){//simphotons is here reflected light. To be added for all PMTs
      if((gRandom->Uniform(1.0))<fQERefl){ 
        if(fParams.TTS>0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
        AddSPE((fParams.TransitTime+ttsTime+simphotons[i].Time-t_min)*fSampling,wave);
      }
    }
    if(pdtype=="pmt"){ //To add direct light for TPB coated PMTs
      sim::SimPhotons auxphotons;
      double ttpb = 0; 
      for (auto& mapMember: auxmap){
 	 if(mapMember.first==ch) auxphotons=mapMember.second;}
      for(size_t j=0; j<auxphotons.size(); j++){//auxphotons is direct light
        if((gRandom->Uniform(1.0))<fQEDirect){ 
          if(fParams.TTS>0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
          ttpb = timeTPB->GetRandom(); //for including TPB emission time
          AddSPE((fParams.TransitTime+ttsTime+auxphotons[j].Time+ttpb-t_min)*fSampling,wave);
        }
      }
    }
    if(fParams.PMTBaselineRMS>0.0) AddLineNoise(wave);
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }

  void DigiPMTSBNDAlg::CreatePDWaveformLite(sim::SimPhotonsLite const& litesimphotons, double t_min, std::vector<double>& wave, int ch, std::string pdtype, std::map<int, sim::SimPhotonsLite> auxmap){
    
    double ttsTime=0;
    std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
    for (auto const& mapMember: photonMap){ //including reflected light for all PMT channels
      for(int i=0; i<mapMember.second; i++){
   	 if((gRandom->Uniform(1.0))<(fQERefl)){
           if(fParams.TTS>0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
 	   AddSPE((fParams.TransitTime+ttsTime+mapMember.first-t_min)*fSampling,wave); }
      }
    }
    if(pdtype=="pmt"){ //To add direct light for TPB coated PMTs
      double ttpb;
      sim::SimPhotonsLite auxphotons;
      for (auto& mapMember: auxmap) //auxphotons is direct light
 	 if(mapMember.first==ch) auxphotons=mapMember.second;
      std::map< int, int > const& auxphotonMap = auxphotons.DetectedPhotons;
      for (auto& mapMember2: auxphotonMap){ 
        for(int i=0; i<mapMember2.second; i++){
   	  if((gRandom->Uniform(1.0))<(fQEDirect)){
           if(fParams.TTS>0.0) ttsTime = Transittimespread(fParams.TTS); //implementing transit time spread
           ttpb = timeTPB->GetRandom(); //for including TPB emission time
           AddSPE((fParams.TransitTime+ttsTime+mapMember2.first+ttpb-t_min)*fSampling,wave);
          }
        }
      }
    }

    if(fParams.PMTBaselineRMS>0.0) AddLineNoise(wave);
    if(fParams.PMTDarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }


  void DigiPMTSBNDAlg::CreateSaturation(std::vector<double>& wave){ //Implementing saturation effects

    for(size_t k=0; k<wave.size(); k++){ 
	if(wave[k]<(fParams.PMTBaseline+fParams.PMTSaturation*fParams.PMTChargeToADC*fParams.PMTMeanAmplitude))
	  wave[k]=fParams.PMTBaseline+fParams.PMTSaturation*fParams.PMTChargeToADC*fParams.PMTMeanAmplitude;	  
    }
  }

  void DigiPMTSBNDAlg::AddLineNoise(std::vector<double>& wave){

    double noise = 0.0;

    for(size_t i = 0; i<wave.size(); i++){
	noise = gRandom->Gaus(0,fParams.PMTBaselineRMS); //gaussian noise
	wave[i] += noise;
    }
  }

  void DigiPMTSBNDAlg::AddDarkNoise(std::vector< double >& wave)
  {
    size_t timeBin=0;

    // Multiply by 10^9 since fParams.DarkNoiseRate is in Hz (conversion from s to ns)
    double darkNoiseTime = static_cast< double >(gRandom->Exp((1.0/fParams.PMTDarkNoiseRate)*1000000000.0));
    while (darkNoiseTime < wave.size()){
      timeBin = (darkNoiseTime);
      AddSPE(timeBin,wave);
      // Find next time to add dark noise
      darkNoiseTime += static_cast< double >(gRandom->Exp((1.0/fParams.PMTDarkNoiseRate)*1000000000.0));
    }
  }
  
  double DigiPMTSBNDAlg::FindMinimumTime(sim::SimPhotons const& simphotons, int ch, std::string pdtype, std::map<int,sim::SimPhotons> auxmap){
    double t_min=1e15;
 
    if(pdtype=="barepmt"){ //TPB non-coated PMTs
      for(size_t i=0; i<simphotons.size(); i++){	 	 
        if(simphotons[i].Time<t_min) t_min = simphotons[i].Time;
      }
    }else{//for coated PMTs
      sim::SimPhotons auxphotons;
      for (auto& mapMember: auxmap){
        if(mapMember.first==ch) auxphotons=mapMember.second;
      }
      auxphotons+=(simphotons); 
      for(size_t i=0; i<auxphotons.size(); i++){	 	 
        if(auxphotons[i].Time<t_min) t_min = auxphotons[i].Time;
      }
    }
    return t_min;
  }
  
  double DigiPMTSBNDAlg::FindMinimumTimeLite(sim::SimPhotonsLite const& litesimphotons, int ch, std::string pdtype, std::map<int,sim::SimPhotonsLite>auxmap){
  
    if(pdtype=="barepmt"){ //TPB non-coated PMTs
      std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
      for (auto const& mapMember: photonMap){
        if(mapMember.second!=0) return (double)mapMember.first;
      }
    }else{
      sim::SimPhotonsLite auxphotons;
      for (auto& mapMember: auxmap){
        if(mapMember.first==ch) auxphotons=mapMember.second;
      }
      auxphotons+=(litesimphotons);
      std::map< int, int > const& auxphotonMap = auxphotons.DetectedPhotons;
      for (auto & mapMember: auxphotonMap){
        if(mapMember.second!=0) return (double)mapMember.first;
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
    fBaseConfig.PMTRiseTime              = config.pmtriseTime();
    fBaseConfig.PMTFallTime              = config.pmtfallTime();
    fBaseConfig.PMTMeanAmplitude         = config.pmtmeanAmplitude();
    fBaseConfig.PMTDarkNoiseRate         = config.pmtdarkNoiseRate();
    fBaseConfig.PMTBaselineRMS           = config.pmtbaselineRMS();
    fBaseConfig.TransitTime              = config.transitTime();
    fBaseConfig.TTS                      = config.tts();
  } 

  std::unique_ptr<DigiPMTSBNDAlg>
  DigiPMTSBNDAlgMaker::operator()(
    detinfo::LArProperties const& larProp,
    detinfo::DetectorClocks const& detClocks
    ) const
  {
    // set the configuration 
    auto params = fBaseConfig;

    // set up parameters
    params.larProp = &larProp;
    params.timeService = &detClocks;
                  
    return std::make_unique<DigiPMTSBNDAlg>(params);
  } // DigiPMTSBNDAlgMaker::create()
    
}

#endif
