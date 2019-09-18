#include "DigiArapucaSBNDAlg.h" 

#ifndef DIGIARAPUCASBNDALG_CXX
#define DIGIARAPUCASBNDALG_CXX
 
//------------------------------------------------------------------------------
//--- opdet::simarapucasbndAlg implementation
//------------------------------------------------------------------------------

namespace opdet{

//DigiArapucaSBNDAlg::DigiArapucaSBNDAlg(fhicl::ParameterSet const& p) 
DigiArapucaSBNDAlg::DigiArapucaSBNDAlg(ConfigurationParameters_t const& config) 
  : fParams(config)
  , fSampling(fParams.timeService->OpticalClock().Frequency()) 
  , fArapucaEffT1(fParams.ArapucaEffT1 / fParams.larProp->ScintPreScale())
  , fArapucaEffT2(fParams.ArapucaEffT2 / fParams.larProp->ScintPreScale())
  , fArapucaEffx(fParams.ArapucaEffx / fParams.larProp->ScintPreScale())
  {

    std::cout << "arapucas corrected efficiencies = " << fArapucaEffT1 << ", " << fArapucaEffT2 << " and " << fArapucaEffx << std::endl;
   // std::cout << "optical clock = " << fSampling << std::endl;

    if(fArapucaEffT1>1.0001 || fArapucaEffT2>1.0001 || fArapucaEffx>1.0001)
	std::cout << "WARNING: Quantum efficiency set in fhicl file " << fParams.ArapucaEffT1 << " or " << fParams.ArapucaEffT2 << " or " << fParams.ArapucaEffx << " seems to be too large! Final QE must be equal to or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << fParams.larProp->ScintPreScale() << std::endl;

    TimeArapucaT1 = new TH1D("Time Profile T1", "", 150, 0.0, 150.0);//histogram that stores the arrival time of photons at SiPM (t=0 is the time is reaches the outside of the optical window) for arapuca T1

    double x[150] = {374, 1455, 2002, 2230, 2336, 2296, 2087, 2006, 1831, 1716, 1623, 1553, 1437, 1327, 1348, 1260, 1237, 1234, 1164, 1122, 1070, 1092, 993, 1002, 892, 969, 951, 907, 909, 961, 863, 857, 915, 900, 853, 842, 790, 780, 779, 808, 781, 749, 723, 737, 723, 755, 732, 680, 724, 631, 656, 693, 669, 632, 636, 643, 632, 640, 608, 615, 597, 633, 602, 545, 591, 595, 551, 574, 567, 507, 545, 535, 552, 519, 537, 563, 523, 461, 550, 510, 514, 469, 517, 493, 466, 460, 488, 446, 474, 516, 451, 451, 457, 465, 450, 456, 493, 441, 441, 475, 433, 419, 435, 405, 392, 410, 430, 404, 392, 407, 435, 411, 383, 422, 394, 397, 413, 366, 389, 376, 366, 372, 375, 345, 370, 368, 370, 390, 351, 382, 373, 380, 377, 339, 372, 371, 351, 360, 338, 365, 309, 187, 95, 41, 19, 7, 2, 3, 0, 0};
 
    for(size_t i=1; i<=150; i++)TimeArapucaT1->SetBinContent(i,x[i-1]);

    TimeArapucaT2 = new TH1D("Time Profile T1", "", 90, 0.0, 90.0);//histogram that stores the arrival time of photons at SiPM (t=0 is the time is reaches the outside of the optical window) for arapuca T2
 
    double x2[90]={5051, 8791, 9054, 8777, 8045, 7009, 6304, 5637, 4828, 4320, 3821, 3333, 2968, 2629, 2364, 2060, 1786, 1624, 1368, 1174, 1115, 935, 820, 725, 627, 535, 500, 453, 386, 355, 315, 279, 221, 196, 198, 181, 120, 115, 128, 109, 79, 67, 77, 59, 48, 39, 40, 37, 32, 39, 22, 25, 18, 20, 14, 19, 8, 9, 6, 11, 13, 4, 11, 4, 4, 5, 3, 2, 4, 4, 5, 2, 6, 0, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 1, 0, 1, 0, 0};

    for(size_t i=1; i<=90; i++)TimeArapucaT2->SetBinContent(i,x2[i-1]);

    TimeArapucaX = new TH1D("Time Profile X-Arapuca", "", 20, 0.0, 20.0);//histogram that stores the arrival time of photons at SiPM (t=0 is the time is reaches the outside of the optical window) for X-arapuca
 
    double x3[20]= {0.0639064, 0.205121, 0.225523, 0.184618, 0.121112, 0.0859086, 0.049705, 0.0268027, 0.0171017, 0.010401, 0.00410041, 0.00290029, 0.00190019, 0.00050005, 0.00010001, 0, 0.00020002, 0.00010001, 0, 0};

    for(size_t i=1; i<=20; i++)TimeArapucaX->SetBinContent(i,x3[i-1]);

    fSampling=fSampling/1000; //in GHz to cancel with ns
    pulsesize=fParams.PulseLength*fSampling;
    wsp.resize(pulsesize);

    for(int i=0; i<pulsesize; i++)
	wsp[i]=(Pulse1PE(static_cast< double >(i)/fSampling));
//Random number engine initialization
    int seed = time(NULL);
    gRandom = new TRandom3(seed);
  }

  DigiArapucaSBNDAlg::~DigiArapucaSBNDAlg()
  { }

  void DigiArapucaSBNDAlg::ConstructWaveform(int ch, sim::SimPhotons const& simphotons, std::vector<std::vector<short unsigned int>>& waveforms, std::string pdName, double start_time, unsigned n_samples){	

    std::vector<double> waves(std::vector<double>(n_samples,fParams.Baseline));
    CreatePDWaveform(simphotons, start_time, waves, pdName);
    waveforms[ch].resize(n_samples);
    waveforms[ch] = std::vector<short unsigned int> (waves.begin(), waves.end());
  }

  void DigiArapucaSBNDAlg::ConstructWaveformLite(int ch, sim::SimPhotonsLite const& litesimphotons, std::vector<std::vector<short unsigned int>>& waveforms, std::string pdName, double start_time, unsigned n_samples){	

    std::vector<double> waves(std::vector<double>(n_samples,fParams.Baseline));
    std::map< int, int > const& photonMap = litesimphotons.DetectedPhotons;
    CreatePDWaveformLite(photonMap, start_time, waves, pdName);
    waveforms[ch].resize(n_samples);
    waveforms[ch] = std::vector<short unsigned int> (waves.begin(), waves.end());
  }

  void DigiArapucaSBNDAlg::AddSPE(size_t time_bin, std::vector<double>& wave, int nphotons){//adding single pulse
    size_t min=0, max=0;

    if(time_bin<fNsamples){
	min=time_bin;
	max=time_bin+pulsesize < fNsamples ? time_bin+pulsesize : fNsamples;
	for(size_t i = min; i<= max; i++){
		wave[i]+= (wsp[i-min])*(double)nphotons;	
	}		
    }
  }

  double DigiArapucaSBNDAlg::Pulse1PE(double time) const//single pulse waveform
  {
    if (time < fParams.PeakTime) return (fParams.ADC*fParams.MeanAmplitude*std::exp((time - fParams.PeakTime)/fParams.RiseTime));
    else return (fParams.ADC*fParams.MeanAmplitude*std::exp(-(time - fParams.PeakTime)/fParams.FallTime));
  }

  void DigiArapucaSBNDAlg::AddLineNoise(std::vector< double >& wave)
  {
    double noise;
    for(size_t i = 0; i < wave.size(); i++){
        noise= gRandom->Gaus(0, fParams.BaselineRMS); //gaussian baseline noise  
        wave[i] += noise; 
    }
  }

  void DigiArapucaSBNDAlg::AddDarkNoise(std::vector< double >& wave)
  {
    int nCT;
    // Multiply by 10^9 since fDarkNoiseRate is in Hz (conversion from s to ns)
    double darkNoiseTime = static_cast< double >(gRandom->Exp((1.0/fParams.DarkNoiseRate)*1000000000.0));
    while (darkNoiseTime < wave.size()){
	size_t timeBin = (darkNoiseTime);
	  if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	  else nCT=1;
	  AddSPE(timeBin,wave,nCT);
        // Find next time to add dark noise
        darkNoiseTime += static_cast< double >(gRandom->Exp((1.0/fParams.DarkNoiseRate)*1000000000.0));
    }
  }

  double DigiArapucaSBNDAlg::FindMinimumTime(sim::SimPhotons const& simphotons){
    double t_min=1e15;
      for(size_t i=0; i<simphotons.size(); i++){	 	 
      	if(simphotons[i].Time<t_min) t_min = simphotons[i].Time;
      }
    return t_min;
  }

  void DigiArapucaSBNDAlg::CreatePDWaveform(sim::SimPhotons const& simphotons, double t_min, std::vector<double>& wave, std::string pdtype){
    int nCT=1;
    double tphoton=0;
    if(pdtype=="arapucaT1"){
	for(size_t i=0; i<simphotons.size(); i++){
	  if((gRandom->Uniform(1.0))<fArapucaEffT1){ //Sample a random subset according to Arapuca's efficiency
	    tphoton=simphotons[i].Time;
	    tphoton+=(TimeArapucaT1->GetRandom());
            tphoton-=t_min;
 	    if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	    else nCT=1;
	    AddSPE(tphoton*fSampling,wave,nCT);
	  }
	}
    }
    if(pdtype=="arapucaT2"){   
	for(size_t i=0; i<simphotons.size(); i++){
 	  if((gRandom->Uniform(1.0))<fArapucaEffT2){ //Sample a random subset according to Arapuca's efficiency.
	    tphoton=simphotons[i].Time;
  	    tphoton+=(TimeArapucaT2->GetRandom());
            tphoton-=t_min;
 	    if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	    else nCT=1;
	    AddSPE(tphoton*fSampling,wave,nCT);
	  }
	}
    }
    if(pdtype=="xarapucaprime"){   
	for(size_t i=0; i<simphotons.size(); i++){
 	  if((gRandom->Uniform(1.0))<fArapucaEffx){ 
	    tphoton=simphotons[i].Time;
 	    tphoton+=(TimeArapucaX->GetRandom());
            tphoton-=t_min;
 	    if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	    else nCT=1;
	    AddSPE(tphoton*fSampling,wave,nCT);
	  }
	}
    }

    if(fParams.BaselineRMS>0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }

  void DigiArapucaSBNDAlg::CreateSaturation(std::vector<double>& wave){ //Implementing saturation effects
    for(size_t k=0; k<fNsamples; k++){
	if(wave[k]>(fParams.Baseline+fParams.Saturation*fParams.ADC*fParams.MeanAmplitude))
	  wave[k]=fParams.Baseline+fParams.Saturation*fParams.ADC*fParams.MeanAmplitude;	  
    }
  }

  void DigiArapucaSBNDAlg::CreatePDWaveformLite(std::map< int, int > const& photonMap, double t_min, std::vector<double>& wave, std::string pdtype){
    double tphoton=0;
    int nCT=1;
    for (auto const& mapMember: photonMap){
      for(int i=0; i<mapMember.second; i++){
         if(pdtype=="arapucaT1" && (gRandom->Uniform(1.0))<fArapucaEffT1){
  	    tphoton=(TimeArapucaT1->GetRandom());
 	    tphoton+=mapMember.first-t_min;
 	    if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	    else nCT=1;
	    AddSPE(tphoton*fSampling,wave,nCT);
         }
         if(pdtype=="arapucaT2" && (gRandom->Uniform(1.0))<fArapucaEffT2){
  	    tphoton=(TimeArapucaT2->GetRandom());
 	    tphoton+=mapMember.first-t_min;
 	    if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	    else nCT=1;
	    AddSPE(tphoton*fSampling,wave,nCT);
         }
	if(pdtype=="xarapucaprime" && (gRandom->Uniform(1.0))<fArapucaEffx){
	   tphoton=(TimeArapucaX->GetRandom());
	   tphoton+=mapMember.first-t_min;
 	   if(fParams.CrossTalk>0.0 && (gRandom->Uniform(1.0))<fParams.CrossTalk) nCT=2;
	   else nCT=1;
	   AddSPE(tphoton*fSampling,wave,nCT);
	}
      }
    }
    if(fParams.BaselineRMS>0.0) AddLineNoise(wave);
    if(fParams.DarkNoiseRate > 0.0) AddDarkNoise(wave);
    CreateSaturation(wave);
  }

  double DigiArapucaSBNDAlg::FindMinimumTimeLite(std::map< int, int > const& photonMap){
    for (auto const& mapMember: photonMap){
 	 if(mapMember.second!=0) return (double)mapMember.first;
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
    fBaseConfig.ArapucaEffx       = config.arapucaEffx();
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
    detinfo::DetectorClocks const& detClocks
    ) const
  {
    // set the configuration 
    auto params = fBaseConfig;

    // set up parameters
    params.larProp = &larProp;
    params.timeService = &detClocks;
 
    return std::make_unique<DigiArapucaSBNDAlg>(params);
  } // DigiArapucaSBNDAlgMaker::create()

}

#endif
