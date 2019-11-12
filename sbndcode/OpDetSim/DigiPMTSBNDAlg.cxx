#include "DigiPMTSBNDAlg.h"

#ifndef DIGIPMTSBNDALG_CXX
#define DIGIPMTSBNDALG_CXX
 
//------------------------------------------------------------------------------
//--- opdet::simpmtsbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet{


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

    if(fQERefl>1.0001 || fQEDirect>1.0001)
	std::cout << "WARNING: Quantum efficiency set in fhicl file " << fParams.QERefl << " or " << fParams.QEDirect << " seems to be too large! Final QE must be equal or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << fParams.larProp->ScintPreScale() << std::endl;

    fSampling=fSampling/1000.0; //in GHz, to cancel with ns
  
//Random number engine initialization
    //int seed = time(NULL);
    //gRandom = new TRandom3(seed);

    timeTPB = new TH1D("Time TPB", "", 1000, 0.0, 1000.0);//histogram that stores the emission time of photons converted by TPB

    double x[1000]={12321, 10239, 8303, 6975, 5684, 4667, 4031, 3446, 2791, 2485, 2062, 1724, 1419, 1367, 1111, 982, 974, 822, 732, 653, 665, 511, 500, 452, 411, 439, 409, 357, 342, 357, 302, 296, 316, 271, 286, 265, 260, 288, 279, 238, 214, 242, 232, 238, 251, 239, 200, 225, 182, 190, 206, 194, 188, 227, 210, 198, 170, 184, 158, 160, 170, 183, 168, 143, 158, 140, 167, 145, 154, 162, 155, 115, 143, 148, 124, 126, 133, 122, 91, 130, 90, 124, 135, 112, 94, 81, 107, 99, 109, 78, 83, 75, 68, 97, 69, 74, 91, 84, 84, 74, 68, 73, 71, 55, 68, 40, 55, 63, 71, 62, 63, 60, 71, 55, 62, 53, 54, 58, 63, 39, 42, 56, 44, 33, 36, 43, 60, 49, 50, 51, 52, 49, 47, 57, 39, 45, 41, 23, 41, 26, 29, 51, 23, 45, 26, 50, 39, 20, 44, 27, 14, 17, 13, 35, 20, 25, 26, 26, 29, 31, 20, 17, 28, 24, 28, 34, 22, 16, 17, 21, 23, 33, 15, 30, 8, 20, 15, 20, 14, 17, 18, 21, 16, 20, 22, 24, 14, 18, 25, 13, 10, 13, 11, 18, 9, 4, 13, 23, 10, 13, 15, 26, 21, 18, 15, 17, 6, 15, 9, 13, 14, 6, 13, 9, 9, 6, 8, 7, 13, 13, 11, 13, 8, 5, 8, 13, 7, 9, 6, 14, 11, 11, 9, 10, 13, 9, 4, 3, 17, 3, 5, 1, 5, 5, 5, 15, 4, 6, 3, 11, 3, 10, 8, 8, 7, 5, 8, 7, 13, 7, 7, 15, 5, 6, 9, 9, 7, 4, 9, 7, 5, 7, 5, 5, 6, 3, 8, 6, 4, 12, 7, 4, 4, 6, 7, 9, 3, 2, 3, 4, 4, 1, 9, 9, 2, 2, 2, 4, 3, 3, 1, 5, 1, 7, 4, 6, 4, 6, 7, 4, 4, 5, 2, 3, 2, 8, 4, 9, 4, 4, 8, 2, 2, 2, 0, 2, 14, 4, 3, 2, 3, 4, 5, 3, 7, 1, 4, 1, 1, 8, 3, 5, 2, 1, 7, 4, 5, 0, 5, 6, 4, 2, 6, 1, 4, 5, 0, 0, 4, 1, 4, 6, 2, 0, 4, 3, 4, 3, 3, 8, 4, 1, 2, 3, 2, 6, 7, 4, 2, 5, 6, 3, 2, 6, 5, 3, 1, 4, 6, 3, 0, 2, 2, 1, 0, 0, 5, 4, 3, 3, 3, 9, 0, 4, 2, 6, 0, 2, 6, 4, 6, 1, 0, 5, 3, 1, 1, 4, 0, 1, 1, 2, 2, 4, 5, 7, 5, 3, 7, 6, 3, 2, 1, 3, 0, 4, 4, 1, 2, 4, 6, 11, 7, 5, 5, 5, 4, 2, 5, 2, 2, 3, 0, 6, 3, 2, 3, 3, 8, 0, 0, 1, 2, 1, 0, 3, 6, 1, 6, 1, 4, 5, 0, 2, 6, 0, 3, 7, 0, 2, 5, 2, 6, 3, 5, 2, 2, 1, 5, 5, 0, 3, 3, 2, 3, 6, 0, 0, 1, 0, 1, 4, 4, 2, 3, 4, 3, 7, 1, 1, 3, 2, 2, 2, 2, 4, 9, 4, 8, 2, 2, 6, 5, 2, 6, 1, 2, 6, 7, 0, 5, 0, 4, 0, 1, 4, 1, 2, 2, 1, 0, 2, 4, 1, 0, 3, 3, 0, 6, 2, 0, 0, 3, 0, 2, 2, 3, 3, 2, 0, 2, 1, 3, 1, 2, 1, 1, 2, 4, 3, 0, 2, 4, 2, 3, 3, 3, 5, 5, 2, 2, 1, 0, 2, 0, 1, 0, 0, 5, 1, 3, 8, 4, 3, 2, 6, 4, 1, 3, 2, 0, 9, 4, 2, 7, 2, 0, 0, 2, 1, 3, 4, 2, 3, 3, 3, 2, 8, 6, 3, 1, 3, 3, 0, 0, 3, 0, 6, 1, 0, 2, 0, 0, 1, 2, 7, 2, 1, 0, 1, 6, 3, 2, 0, 1, 0, 2, 3, 5, 3, 6, 4, 1, 1, 0, 0, 7, 1, 1, 1, 1, 8, 1, 1, 0, 0, 0, 1, 2, 2, 7, 1, 1, 0, 1, 1, 3, 1, 1, 1, 3, 5, 2, 1, 0, 3, 2, 5, 0, 5, 4, 2, 5, 3, 3, 0, 0, 5, 0, 5, 1, 4, 0, 1, 6, 1, 6, 1, 2, 1, 2, 4, 0, 8, 3, 1, 7, 1, 2, 4, 4, 2, 3, 5, 0, 4, 5, 2, 1, 1, 5, 2, 0, 4, 2, 0, 2, 4, 4, 4, 4, 5, 0, 3, 0, 2, 3, 3, 0, 0, 6, 1, 1, 6, 10, 0, 2, 0, 1, 4, 1, 0, 1, 3, 2, 0, 1, 1, 0, 0, 8, 4, 0, 0, 3, 0, 0, 4, 0, 1, 4, 0, 1, 0, 2, 1, 5, 2, 2, 0, 0, 0, 2, 1, 0, 0, 4, 0, 3, 1, 1, 2, 1, 1, 2, 0, 1, 3, 3, 0, 2, 2, 3, 2, 0, 1, 2, 0, 0, 1, 0, 4, 2, 3, 0, 0, 1, 2, 2, 6, 2, 3, 4, 2, 3, 1, 2, 7, 3, 2, 3, 3, 1, 0, 0, 4, 2, 9, 0, 3, 2, 2, 0, 1, 1, 1, 8, 1, 4, 2, 0, 4, 0, 2, 0, 0, 0, 3, 0, 0, 4, 0, 2, 2, 5, 1, 3, 6, 0, 0, 0, 0, 0, 3, 2, 1, 0, 1, 1, 1, 1, 2, 0, 1, 4, 4, 1, 1, 0, 0, 9, 2, 1, 1, 2, 0, 2, 2, 2, 2, 0, 1, 7, 0, 7, 0, 5, 0, 5, 1, 1, 0, 1, 1, 1, 0, 1, 0, 2, 3, 2, 1, 1, 0, 0, 2, 4, 0, 0, 2, 4, 2, 5, 2, 2, 1, 1, 4, 0, 1, 2, 1, 3, 3, 0, 2, 1, 3, 0, 2, 2, 1, 0, 4, 5, 3, 0, 0, 2, 1, 3, 0, 2, 0, 3, 2, 2, 0, 3, 1, 0, 5, 2, 2, 4, 6, 3, 2, 2, 2, 1, 4, 6, 1, 2, 1, 2, 6, 1, 2};


    for(int i=1; i<=1000; i++)timeTPB->SetBinContent(i,x[i-1]);

    //shape of single pulse
    if (fParams.SinglePEmodel == 1  ) {
      
           std::cout << " using testbench pe response " << std::endl;
      pulsesize=450;
      wsp={0, 0, 0, 0, 0, 0, -0.0165309, -0.01653, -0.01653, -0.01653, -0.01653, -0.03398, -0.03398, -0.03398, -0.03398, -0.03398, -0.03398, -0.123, -0.123, -0.2137, -0.2556, -0.5209, -1.51, -7.384, -20.064, -35.448, -47.957, -51.942, -45.182, -33.576, -22.9686, -16.41, -11.1456, -7.4829, -5.568, -3.217, -1.889, -1.889, -0.9656, -0.451, -0.325, -0.252, 0.0357, 0.09189, 0.12, 0.1738, 0.29, 0.318, 0.487, 0.5, 0.52, 0.5329, 0.56, 0.598, 0.61, 0.65, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.781, 0.787, 0.793, 0.799, 0.805, 0.811, 0.817, 0.823, 0.829, 0.83, 0.841, 0.847, 0.853, 0.859, 0.865, 0.87, 0.875, 0.88, 0.885, 0.89, 0.895, 0.9, 0.912, 0.917, 0.92, 0.923, 0.924, 0.925, 0.926, 0.927, 0.928, 0.929, 0.93, 0.931, 0.932, 0.931, 0.93, 0.9285, 0.9272, 0.9242, 0.9212, 0.9184, 0.9156, 0.9128, 0.91, 0.9072, 0.9044, 0.9016, 0.8988, 0.896, 0.8932, 0.8904, 0.8876, 0.8848, 0.882, 0.8792, 0.8764, 0.8736, 0.8708, 0.868, 0.8652, 0.8624, 0.8596, 0.8568, 0.854, 0.8512, 0.8484, 0.8456, 0.8428, 0.84, 0.8372, 0.8344, 0.8316, 0.8288, 0.826, 0.8232, 0.8204, 0.8176, 0.8148, 0.812, 0.8092, 0.8064, 0.8036, 0.8008, 0.798, 0.7952, 0.7924, 0.7896, 0.7868, 0.784, 0.7812, 0.7784, 0.7756, 0.7728, 0.77, 0.7672, 0.7644, 0.7616, 0.7588, 0.756, 0.7532, 0.7504, 0.7476, 0.7448, 0.742, 0.7392, 0.7364, 0.7336, 0.7308, 0.728, 0.7252, 0.7224, 0.7196, 0.7168, 0.714, 0.7112, 0.7084, 0.7056, 0.7028, 0.7, 0.6972, 0.6944, 0.6916, 0.6888, 0.686, 0.6832, 0.6804, 0.6776, 0.6748, 0.672, 0.6692, 0.6664, 0.6636, 0.6608, 0.658, 0.6552, 0.6524, 0.6496, 0.6468, 0.644, 0.6412, 0.6384, 0.6356, 0.6328, 0.63, 0.6272, 0.6244, 0.6216, 0.6188, 0.616, 0.6132, 0.6104, 0.6076, 0.6048, 0.602, 0.5992, 0.5964, 0.5936, 0.5908, 0.588, 0.5852, 0.5824, 0.5796, 0.5768, 0.5740, 0.5712, 0.5684, 0.5656, 0.5628, 0.56, 0.5572, 0.5544, 0.5516, 0.5488, 0.546, 0.5432, 0.5404, 0.5376, 0.5348, 0.5320, 0.5292, 0.5264, 0.5236, 0.5208, 0.5180, 0.5152, 0.5124, 0.5096, 0.5068, 0.504, 0.5012, 0.4984, 0.4956, 0.4928, 0.49, 0.4872, 0.4844, 0.4816, 0.4788, 0.476, 0.4732, 0.4704, 0.4676, 0.4648, 0.462, 0.4592, 0.4564, 0.4536, 0.4508, 0.448, 0.4452, 0.4424, 0.4396, 0.4368, 0.434, 0.4312, 0.4284, 0.4256, 0.4228, 0.42, 0.4172, 0.4144, 0.4116, 0.4088, 0.406, 0.4032, 0.4004, 0.3976, 0.3948, 0.392, 0.3892, 0.3864, 0.3836, 0.3808, 0.378, 0.3752, 0.3724, 0.3696, 0.3668, 0.364, 0.3612, 0.3584, 0.3556, 0.3528, 0.35, 0.3472, 0.3444, 0.3416, 0.3388, 0.336, 0.3332, 0.3304, 0.3276, 0.3248, 0.322, 0.3192, 0.3164, 0.3136, 0.3108, 0.308, 0.3052, 0.3024, 0.2996, 0.2968, 0.2940, 0.2912, 0.2884, 0.2856, 0.2828, 0.28, 0.2772, 0.2744, 0.2716, 0.2688, 0.266, 0.2632, 0.2604, 0.2576, 0.2548, 0.252, 0.2492, 0.2464, 0.2436, 0.2408, 0.238, 0.2352, 0.2324, 0.2296, 0.2276, 0.2256, 0.22446, 0.22246, 0.22046, 0.21846, 0.21646, 0.21447, 0.21247, 0.21047, 0.20847, 0.20647, 0.20447, 0.20247, 0.20047, 0.19847, 0.19647, 0.19447, 0.19247, 0.19047, 0.18847, 0.18647, 0.18447, 0.18247, 0.18047, 0.17847, 0.17647, 0.17447, 0.17247, 0.17047, 0.16847, 0.16647, 0.16447, 0.16247, 0.16047, 0.15847, 0.15647, 0.15447, 0.15247, 0.15047, 0.14847, 0.14647, 0.14447, 0.14247, 0.14047, 0.13847, 0.13647, 0.13447, 0.13247, 0.13047, 0.12847, 0.12647, 0.12447, 0.12247, 0.12047, 0.11847, 0.11647, 0.11447, 0.11247, 0.11047, 0.10847, 0.10647, 0.10447, 0.10247, 0.10047, 0.09847, 0.09647, 0.09447, 0.09247, 0.09047, 0.08847, 0.08647, 0.08447, 0.08, 0.075, 0.07, 0.065, 0.06, 0.054, 0.05148, 0.04948, 0.04748, 0.04548, 0.04348, 0.04148, 0.03948, 0.03748, 0.03548, 0.03348, 0.03148, 0.02948, 0.02748, 0.02548, 0.02348, 0.02148, 0.01948, 0.01748, 0.01548, 0.01348, 0.01078};
    }
    else {
      std::cout << " using ideal pe response " << std::endl;
      //shape of single pulse
      sigma1 = fParams.PMTRiseTime/(std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      sigma2 = fParams.PMTFallTime/(std::sqrt(2.0) * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9))));
      
      pulsesize=(int)((6*sigma2+fParams.TransitTime)*fSampling);
      wsp.resize(pulsesize);
      for(int i=0; i<pulsesize; i++){
	wsp[i]=(Pulse1PE(static_cast< double >(i)/fSampling));
      }
    }
  } // end constructor

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
    //tts = gRandom->Gaus(0,sigma);	
    tts = CLHEP::RandGauss::shoot(fEngine, 0, sigma);

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
      //if((gRandom->Uniform(1.0))<fQERefl){ 
      if(CLHEP::RandFlat::shoot(fEngine, 1.0)<fQERefl){ 
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
        //if((gRandom->Uniform(1.0))<fQEDirect){ 
        if(CLHEP::RandFlat::shoot(fEngine, 1.0)<fQEDirect){ 
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
   	 //if((gRandom->Uniform(1.0))<(fQERefl)){
   	 if(CLHEP::RandFlat::shoot(fEngine, 1.0)<(fQERefl)){
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
   	  //if((gRandom->Uniform(1.0))<(fQEDirect)){
   	  if(CLHEP::RandFlat::shoot(fEngine, 1.0)<(fQEDirect)){
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
	//noise = gRandom->Gaus(0,fParams.PMTBaselineRMS); //gaussian noise
	noise = CLHEP::RandGauss::shoot(fEngine, 0, fParams.PMTBaselineRMS); //gaussian noise
	wave[i] += noise;
    }
  }

  void DigiPMTSBNDAlg::AddDarkNoise(std::vector< double >& wave)
  {
    size_t timeBin=0;

    // Multiply by 10^9 since fParams.DarkNoiseRate is in Hz (conversion from s to ns)
    //double darkNoiseTime = static_cast< double >(gRandom->Exp((1.0/fParams.PMTDarkNoiseRate)*1000000000.0));
    double darkNoiseTime = CLHEP::RandExponential::shoot(fEngine, (1.0/fParams.PMTDarkNoiseRate)*1000000000.0);
    while (darkNoiseTime < wave.size()){
      timeBin = (darkNoiseTime);
      AddSPE(timeBin,wave);
      // Find next time to add dark noise
      //darkNoiseTime += static_cast< double >(gRandom->Exp((1.0/fParams.PMTDarkNoiseRate)*1000000000.0));
      darkNoiseTime += CLHEP::RandExponential::shoot(fEngine, (1.0/fParams.PMTDarkNoiseRate)*1000000000.0);
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
    fBaseConfig.SinglePEmodel            = config.singlePEmodel();
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
