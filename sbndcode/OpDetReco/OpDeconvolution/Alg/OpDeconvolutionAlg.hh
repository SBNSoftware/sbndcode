////////////////////////////////////////////////////////////////////////
// Class:       OpDeconvolutionAlg
// File:        OpDeconvolutionAlg.hh
//
// Generated at Fri Jun 11 05:27:39 2021 by Francisco Nicolas-Arnaldos
////////////////////////////////////////////////////////////////////////

#ifndef SBNDOPDECONVOLUTION_H
#define SBNDOPDECONVOLUTION_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/Simulation/sim.h"

#include <memory>

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/Utilities/LArFFT.h"
#include "TFile.h"

#include <cmath>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TComplex.h"

#include "nurandom/RandomUtils/NuRandomService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


namespace opdet {
  class OpDeconvolutionAlg;
}


class opdet::OpDeconvolutionAlg {
public:
  explicit OpDeconvolutionAlg(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpDeconvolutionAlg(OpDeconvolutionAlg const&) = delete;
  OpDeconvolutionAlg(OpDeconvolutionAlg&&) = delete;
  OpDeconvolutionAlg& operator=(OpDeconvolutionAlg const&) = delete;
  OpDeconvolutionAlg& operator=(OpDeconvolutionAlg&&) = delete;

  // Required functions.
  std::vector<raw::OpDetWaveform> RunDeconvolution(std::vector<raw::OpDetWaveform> const& wfHandle);


private:
  bool fDebug;
  int fMaxFFTSizePow;
  std::vector<double> fSinglePEWave;
  bool fApplyExpoAvSmooth;
  bool fApplyUnAvSmooth;
  float fExpoAvSmoothPar;
  short unsigned int fUnAvNeighbours;
  double fHypoSignalTimeWindow;
  bool fHypoSignalCustom;
  std::vector<double> fHypoSignalTimeConsts;
  std::vector<double> fHypoSignalWeights;
  double fHypoSignalScale;
  double fPMTChargeToADC;
  double fDecoWaveformPrecision;
  short unsigned int fBaselineSample;
  std::string fOpDetDataFile;
  std::string fFilter;
  bool fScaleHypoSignal;
  bool fUseParamFilter;
  std::vector<double> fFilterParams;

  double fNormUnAvSmooth;
  double fSamplingFreq;
  size_t MaxBinsFFT;

  unsigned int NDecoWf;

  //Random Numbers
  std::unique_ptr<CLHEP::RandGauss>  fCLHEPEGauss;

  TF1 *fFilterTF1;
  std::vector<double> fSignalHypothesis;
  std::vector<double> fNoiseHypothesis;

  // Declare member data here.

  // Declare member functions
  void ApplyExpoAvSmoothing(std::vector<double>& wf);
  void ApplyUnAvSmoothing(std::vector<double>& wf);
  size_t WfSizeFFT(size_t n);
  std::vector<double> ScintArrivalTimesShape(size_t n, detinfo::LArProperties const& lar_prop);
  void SubtractBaseline(std::vector<double> &wf, double baseline);
  void EstimateBaselineStdDev(std::vector<double> &wf, double &_mean, double &_stddev);
  std::vector<TComplex> DeconvolutionKernel(size_t size, double baseline_stddev, double snr_scaling);

  //Load TFileService serrvice
  art::ServiceHandle<art::TFileService> tfs;
  //Load FFT serrvice
  art::ServiceHandle<util::LArFFT> fft_service;
};


opdet::OpDeconvolutionAlg::OpDeconvolutionAlg(fhicl::ParameterSet const& p)
{
  //read fhicl paramters
  fDebug = p.get< bool >("Debug");
  fMaxFFTSizePow = p.get< int >("MaxFFTSizePow", 15);
  fApplyExpoAvSmooth   = p.get< bool >("ApplyExpoAvSmooth");
  fApplyUnAvSmooth   = p.get< bool >("ApplyUnAvSmooth");
  fExpoAvSmoothPar = p.get< float >("ExpoAvSmoothPar");
  fUnAvNeighbours = p.get< short unsigned int >("UnAvNeighbours");
  fHypoSignalTimeWindow = p.get< double >("HypoSignalTimeWindow");
  fHypoSignalCustom = p.get< bool >("HypoSignalCustom", false);
  if(fHypoSignalCustom){
    fHypoSignalTimeConsts = p.get< std::vector<double> >("HypoSignalTimeConsts");
    fHypoSignalWeights = p.get< std::vector<double> >("HypoSignalWeights");
  }
  fHypoSignalScale = p.get< double >("HypoSignalScale");
  fPMTChargeToADC = p.get< double >("PMTChargeToADC");
  fDecoWaveformPrecision = p.get< double >("DecoWaveformPrecision");
  fBaselineSample = p.get< short unsigned int >("BaselineSample");
  fOpDetDataFile = p.get< std::string >("OpDetDataFile");
  fFilter = p.get< std::string >("Filter");
  fScaleHypoSignal = p.get< bool >("ScaleHypoSignal");
  fUseParamFilter = p.get< bool >("UseParamFilter");
  fFilterParams = p.get< std::vector<double> >("FilterParams");


  fNormUnAvSmooth=1./(2*fUnAvNeighbours+1);
  NDecoWf=0;
  MaxBinsFFT=std::pow(2, fMaxFFTSizePow);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSamplingFreq=clockData.OpticalClock().Frequency()/1000.;//in GHz
  auto const* lar_prop = lar::providerFrom<detinfo::LArPropertiesService>();

  //Initizalize random numbers
  art::ServiceHandle<rndm::NuRandomService> seedSvc;
  CLHEP::HepJamesRandom *engine = new CLHEP::HepJamesRandom;
  seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(engine), "opdeconvolutionSBND");
  CLHEP::RandGauss GaussGen(engine);
  fCLHEPEGauss = std::make_unique<CLHEP::RandGauss>(engine);

  //Load SER
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fOpDetDataFile, fname);
  TFile* file = TFile::Open(fname.c_str(), "READ");
  std::vector<double>* SinglePEVec_p;
  file->GetObject("SinglePEVec", SinglePEVec_p);
  fSinglePEWave = *SinglePEVec_p;
  //for(size_t k=0; k<fSinglePEWave.size(); k++) std::cout<<"  "<<k<<":"<<fSinglePEWave[k];std::cout<<std::endl;
  mf::LogInfo("OpDeconvolutionAlg")<<"Loaded SER from "<<fOpDetDataFile<<"... size="<<fSinglePEWave.size()<<std::endl;
  fSinglePEWave.resize(MaxBinsFFT, 0);


  if(fUseParamFilter){
    fFilterTF1 = new TF1("FilterTemplate", fFilter.c_str());
    for(size_t k=0; k<fFilterParams.size(); k++)
      fFilterTF1->SetParameter(k, fFilterParams[k]);
    mf::LogInfo("OpDeconvolutionAlg")<<"Creating parametrized filter... TF1:"<<fFilter<<std::endl;
  }
  else{
    //Create light signal hypothesis for "on-fly" Wiener filter
    fSignalHypothesis.resize(MaxBinsFFT, 0);
    if(fFilter=="Wiener")
      fSignalHypothesis = ScintArrivalTimesShape(MaxBinsFFT, *lar_prop);
    else if(fFilter=="Wiener1PE")
      fSignalHypothesis[0]=1;
    mf::LogInfo("OpDeconvolutionAlg")<<"Built light signal hypothesis... L="<<fFilter<<" size"<<fSignalHypothesis.size()<<std::endl;
  }
}


std::vector<raw::OpDetWaveform> opdet::OpDeconvolutionAlg::RunDeconvolution(std::vector<raw::OpDetWaveform> const& wfVector)
{
  std::vector<raw::OpDetWaveform> wfDeco;
  wfDeco.reserve(wfVector.size());

  for(auto const& wf : wfVector)
  {
    //Read waveform
    size_t wfsize=wf.Waveform().size();
    if(wfsize>MaxBinsFFT){
      mf::LogWarning("OpDeconvolutionAlg")<<"Skipping waveform...waveform size is"<<wfsize<<"...maximum allowed FFT size is="<<MaxBinsFFT<<std::endl;
      continue;
    }
    mf::LogInfo("OpDeconvolutionAlg")<<"Deconvolving waveform:"<<NDecoWf<<"...size="<<wfsize<<std::endl;
    size_t wfsizefft=WfSizeFFT(wfsize);

    std::vector<double> wave;
    wave.reserve(wfsizefft);
    wave.assign(wf.Waveform().begin(), wf.Waveform().end());

    //Reserve minimum ADCC value for Wiener filter
    double minADC=*min_element(wave.begin(), wave.end());

    //Apply waveform smoothing
    if(fApplyExpoAvSmooth)
      ApplyExpoAvSmoothing(wave);
    if(fApplyUnAvSmooth)
      ApplyUnAvSmoothing(wave);

    //Estimate baseline standrd deviation
    double baseline_mean=0., baseline_stddev=1.;
    EstimateBaselineStdDev(wave, baseline_mean, baseline_stddev);
    double wfPeakPE=fHypoSignalScale*(baseline_mean-minADC)/fPMTChargeToADC;
    SubtractBaseline(wave, baseline_mean);

    //Create deconvolution kernel
    wave.resize(wfsizefft, 0);
    std::vector<TComplex> fDeconvolutionKernel=DeconvolutionKernel(wfsize, baseline_stddev, wfPeakPE);

    //Deconvolve raw signal (covolve with kernel)
    fft_service->ReinitializeFFT(wfsizefft, "", 20);
    fft_service->Convolute(wave, fDeconvolutionKernel);
    wave.resize(wfsize);

    //Set deconvlved waveform precision and restore baseline before saving
    EstimateBaselineStdDev(wave, baseline_mean, baseline_stddev);
    SubtractBaseline(wave, baseline_mean);
    double fDecoWfScaleFactor=1./fDecoWaveformPrecision;
    std::transform(wave.begin(), wave.end(), wave.begin(), [fDecoWfScaleFactor](double &dec){ return fDecoWfScaleFactor*dec; } );

    //Debbuging and save wf in hist file
    if(fDebug){
      std::cout<<".....Debbuging.....\n";
      auto minADC_ix=min_element(wave.begin(), wave.end());
      std::cout<<"Stamp="<<wf.TimeStamp()<<" OpCh"<<wf.ChannelNumber()<<" MinADC="<<minADC<<" (";
      std::cout<<minADC_ix-wave.begin()<<") Size="<<wf.Waveform().size()<<" ScFactor="<<wfPeakPE<<"\n\n";

      std::string name="h_deco"+std::to_string(NDecoWf)+"_"+std::to_string(wf.ChannelNumber())+"_"+std::to_string(wf.TimeStamp());
      TH1F * h_deco = tfs->make< TH1F >(name.c_str(),";Bin;#PE", MaxBinsFFT, 0, MaxBinsFFT);
      for(size_t k=0; k<wave.size(); k++){
        //if(fDebug) std::cout<<k<<":"<<wave[k]<<":"<<rawsignal[k]<<"  ";
        h_deco->Fill(k, wave[k]);
      }

      name="h_raw"+std::to_string(NDecoWf)+"_"+std::to_string(wf.ChannelNumber())+"_"+std::to_string(wf.TimeStamp());
      TH1F * h_raw = tfs->make< TH1F >(name.c_str(),";Bin;ADC", MaxBinsFFT, 0, MaxBinsFFT);
      for(size_t k=0; k<wf.Waveform().size(); k++){
        //if(fDebug) std::cout<<k<<":"<<wave[k]<<":"<<rawsignal[k]<<"  ";
        h_raw->Fill(k, wf.Waveform()[k]);
      }
    }

    //raw::OpDetWaveform decowf(wf.TimeStamp(), wf.ChannelNumber(), std::vector<short unsigned int> (wave.begin(),  std::next(wave.begin(), wf.Waveform().size()) ) );
    raw::OpDetWaveform decowf( wf.TimeStamp(), wf.ChannelNumber(), std::vector<short unsigned int> (wave.begin(),  wave.end()) );
    wfDeco.push_back(decowf);
    NDecoWf++;
  }

  return wfDeco;
}


void opdet::OpDeconvolutionAlg::ApplyExpoAvSmoothing(std::vector<double>& wf){
  std::transform (std::next(wf.begin(), 1), wf.end(), wf.begin(), std::next(wf.begin(), 1),
    [&](double _x, double _y) { return  fExpoAvSmoothPar*_x+ (1. - fExpoAvSmoothPar)*_y; }  );
}


void opdet::OpDeconvolutionAlg::ApplyUnAvSmoothing(std::vector<double>& wf){
  std::vector<double> wf_aux(wf.begin(), wf.end());
  for(size_t bin=fUnAvNeighbours; bin<wf.size()-fUnAvNeighbours; bin++){
    double sum=0.;
    for(size_t nbin=bin-fUnAvNeighbours; nbin<=bin+fUnAvNeighbours; nbin++)
      sum+=wf_aux[nbin];
    //std::cout<<bin<<" "<<sum<<" "<<sum*fNormUnAvSmooth<<std::endl;
    wf[bin]=sum*fNormUnAvSmooth;
  }
}


size_t opdet::OpDeconvolutionAlg::WfSizeFFT(size_t n){
  if (n && !(n & (n - 1)))
       return n;
  size_t cont=0;
  while(n>0){
    cont++;
    n=(n>>1);
  }
  return 1<<cont;
}


std::vector<double> opdet::OpDeconvolutionAlg::ScintArrivalTimesShape(size_t n, detinfo::LArProperties const& lar_prop){
  if(fHypoSignalCustom && fHypoSignalWeights.size()!=fHypoSignalTimeConsts.size()){
    mf::LogWarning("OpDeconvolutionAlg")<<"HypoSignalWeights and HypoSignalTimeConsts must have the same size"<<
      "...Switching to 2-exponential (default values in LArProperties)"<<std::endl;
  }
  if(!fHypoSignalCustom || fHypoSignalWeights.size()!=fHypoSignalTimeConsts.size()){
    fHypoSignalTimeConsts={lar_prop.ScintFastTimeConst(), lar_prop.ScintSlowTimeConst()};
    fHypoSignalWeights={lar_prop.ScintYieldRatio(), 1.-lar_prop.ScintYieldRatio()};
  }

  mf::LogInfo("OpDeconvolutionAlg")<<"Creating hypothesis signal for Wiener filter..."<<std::endl;
  mf::LogInfo("OpDeconvolutionAlg")<<"Using a "<<fHypoSignalTimeConsts.size()<<"-exponential"<<std::endl;
  mf::LogInfo("OpDeconvolutionAlg")<<"Time constants [ns]:";
  for(auto &tau:fHypoSignalTimeConsts) mf::LogInfo("OpDeconvolutionAlg")<<tau<<"  ";
  mf::LogInfo("OpDeconvolutionAlg")<<"\nWeights:";
  for(auto &w:fHypoSignalWeights) mf::LogInfo("OpDeconvolutionAlg")<<w<<"  ";
  mf::LogInfo("OpDeconvolutionAlg")<<"\n";

  std::vector<double> v(n, 0);
  size_t maxbin=fSamplingFreq*fHypoSignalTimeWindow;
  if(n<maxbin)
    maxbin=n;

  for(size_t k=0; k<maxbin; k++){
    double t = (double)(k) / fSamplingFreq; //in ns
    double value=0;
    for(size_t ix=0; ix<fHypoSignalTimeConsts.size(); ix++){
      value+=fHypoSignalWeights[ix]*std::exp(-1.*t/fHypoSignalTimeConsts[ix]);
    }
    v[k]=value;
  }
  return v;
}


void opdet::OpDeconvolutionAlg::SubtractBaseline(std::vector<double> &wf, double baseline){
  std::transform (wf.begin(), wf.end(), wf.begin(), [&](double _x) { return _x-baseline; }  );
  return;
}


void opdet::OpDeconvolutionAlg::EstimateBaselineStdDev(std::vector<double> &wf, double &_mean, double &_stddev){
  double minADC=*min_element(wf.begin(), wf.end());
  double maxADC=*max_element(wf.begin(), wf.end());
  unsigned nbins=25*ceil(maxADC-minADC);
  TH1F h_std = TH1F("",";;", nbins, 0, (maxADC-minADC)/2);
  TH1F h_mean = TH1F("",";;", nbins, minADC, maxADC);

  for(size_t ix=0; ix<wf.size()-fBaselineSample; ix++){
    double sum2=0, sum=0;
    for(size_t jx=ix; jx<ix+fBaselineSample; jx++){
      sum = sum + wf.at(jx);
    }
    sum/=fBaselineSample;
    for(size_t jx=ix; jx<ix+fBaselineSample; jx++){
      sum2 = sum2 + pow( wf.at(jx)-sum, 2 );
    }
    h_std.Fill( std::sqrt(sum2/fBaselineSample) );
    h_mean.Fill( sum );
    //std::cout<<ix<<":"<<wf[ix]<<":"<<sum<<" ";
  }

  _stddev=h_std.GetXaxis()->GetBinCenter(h_std.GetMaximumBin());
  _mean=h_mean.GetXaxis()->GetBinCenter(h_mean.GetMaximumBin());

  if(fDebug){
    std::cout<<"   -- Estimating baseline...StdDev: "<<_stddev<<" Bias="<<_mean<<std::endl;
    std::cout<<"      .. "<<minADC<<" "<<maxADC<<" "<<maxADC-minADC<<" "<<ceil(log10(maxADC-minADC))<<" "<<nbins<<"\n";

    std::string name="h_baselinestddev_"+std::to_string(NDecoWf)+std::to_string(_mean);
    TH1F * hs_std = tfs->make< TH1F > (name.c_str(),"Baseline StdDev;ADC;# entries",
      h_std.GetNbinsX(), h_std.GetXaxis()->GetXmin(), h_std.GetXaxis()->GetXmax());
    for(int k=1; k<=h_std.GetNbinsX(); k++)
      hs_std->SetBinContent(k, h_std.GetBinContent(k));

    name="h_baselinemean_"+std::to_string(NDecoWf)+std::to_string(_mean);
    TH1F * hs_mean = tfs->make< TH1F >(name.c_str(),"Baseline Mean;ADC;# entries",
      h_mean.GetNbinsX(), h_mean.GetXaxis()->GetXmin(), h_mean.GetXaxis()->GetXmax());
    for(int k=1; k<=h_mean.GetNbinsX(); k++)
      hs_mean->SetBinContent(k, h_mean.GetBinContent(k));
  }

  return;
}


std::vector<TComplex> opdet::OpDeconvolutionAlg::DeconvolutionKernel(size_t wfsize, double baseline_stddev, double snr_scaling){
  //Initizalize kernel
  size_t size=WfSizeFFT(wfsize);
  TComplex kerinit(0,0,false);
  std::vector<TComplex> kernel(size, kerinit);
  fft_service->ReinitializeFFT (size, "", 20);

  //Prepare detector response FFT
  std::vector<double> ser( fSinglePEWave.begin(), std::next(fSinglePEWave.begin(), size) );
  std::vector<TComplex> serfft;
  serfft.resize(size);
  fft_service->DoFFT(ser, serfft);

  if(fUseParamFilter){
    double freq_step=fSamplingFreq/size;
    for(size_t k=0; k<size/2; k++){
      kernel[k]= fFilterTF1->Eval(k*freq_step) / serfft[k] ;
      //std::cout<<k<<" "<<freq<<" "<<filval<<" "<<fSamplingFreq<<"\n";
    }
  }
  else{
    //Build Wiener filter kernel: G = Conj(R) / ( |R|^2 + |N|^2/|L|^2)
    //R=Detector resopnse FFT
    //N=Noise mean spectral power
    //L=True signal mean spectral power

    //Prepare L
    std::vector<double> hypo( fSignalHypothesis.begin(), std::next(fSignalHypothesis.begin(), size) );
    std::vector<TComplex> hypofft;
    hypofft.resize(size);
    fft_service->DoFFT(hypo, hypofft);

    //Prepare Noise Spectral Power
    double noise_power=wfsize*baseline_stddev*baseline_stddev;
    if(fScaleHypoSignal){
      noise_power/=pow(snr_scaling, 2);
    }

    for(size_t k=0; k<size/2; k++){
      double den = pow(TComplex::Abs(serfft[k]), 2) + noise_power / pow(TComplex::Abs(hypofft[k]), 2) ;
      kernel[k]= TComplex::Conjugate( serfft[k] ) / den;
      //std::cout<<k<<":"<<serfft[k]<<":"<<hypofft[k]<<":"<<den << " ";
    }
  }


  if(fDebug){
    std::string name="h_wienerfilter_"+std::to_string(NDecoWf);
    TH1F * hs_wiener = tfs->make< TH1F >
      (name.c_str(),"Wiener Filter;Frequency Bin;Magnitude",size/2, 0, size/2);
    for(size_t k=0; k<size/2; k++)
      hs_wiener->SetBinContent(k, TComplex::Abs( kernel[k]*serfft[k] ) );
  }

  return kernel;
}


#endif
