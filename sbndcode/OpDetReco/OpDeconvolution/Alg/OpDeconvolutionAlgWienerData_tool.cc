////////////////////////////////////////////////////////////////////////
// Specific class tool for OpDeconvolutionAlg
// File: OpDeconvolutionAlgWienerData_tool.hh
// Base class:        OpDeconvolutionAlg.hh
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include <memory>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/Utilities/LArFFT.h"
#include "TFile.h"

#include <cmath>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TComplex.h"

#include "sbndcode/OpDetReco/OpDeconvolution/Alg/OpDeconvolutionAlg.hh"

#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"

namespace opdet {
  class OpDeconvolutionAlgWiener;
}


class opdet::OpDeconvolutionAlgWiener : opdet::OpDeconvolutionAlg {
public:
  explicit OpDeconvolutionAlgWiener(fhicl::ParameterSet const& p);

  ~OpDeconvolutionAlgWiener() {delete fFilterTF1;}

  // Required functions.
  std::vector<raw::OpDetWaveform> RunDeconvolution(std::vector<raw::OpDetWaveform> const& wfHandle) override;

private:
  bool fDebug;
  int fMaxFFTSizePow;
  std::vector<std::vector<double>> fSinglePEWaveVector;
  std::vector<double> fSinglePEWave;
  std::vector<int> fSinglePEChannels;
  std::vector<double> fPeakAmplitude;
  std::vector<int> fSkipChannelList;
  bool fPositivePolarity;
  bool fUseSaturated;
  int fADCSaturationValue;
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
  std::string fFilter;
  std::string fElectronics;
  bool fScaleHypoSignal;
  bool fUseParamFilter;
  std::vector<double> fFilterParams;
  std::vector<std::vector<double>> fFilterParamVector;
  std::vector<double> fFilterParamChannel;
  double fNormUnAvSmooth;
  double fSamplingFreq;
  double fDaphne_Freq;
  size_t MaxBinsFFT;
  unsigned int NDecoWf;

  TF1 *fFilterTF1;

  std::vector<double> fSignalHypothesis;
  std::vector<double> fNoiseHypothesis;

  bool fUseParamFilterInidividualChannel;
  bool fCorrectBaselineOscillations;
  short unsigned int fBaseSampleBins;
  double fBaseVarCut;

  sbndDB::PMTCalibrationDatabase const* fPMTCalibrationDatabaseService;

  // Declare member functions
  void ApplyExpoAvSmoothing(std::vector<double>& wf);
  void ApplyUnAvSmoothing(std::vector<double>& wf);
  size_t WfSizeFFT(size_t n);
  std::vector<double> ScintArrivalTimesShape(size_t n, detinfo::LArProperties const& lar_prop);
  void SubtractBaseline(std::vector<double> &wf, double baseline);
  void EstimateBaselineStdDev(std::vector<double> &wf, double &_mean, double &_stddev);
  void CorrectBaselineOscillations(std::vector<double> &wf);
  std::vector<TComplex> DeconvolutionKernel(size_t size, double baseline_stddev, double snr_scaling);

  //Load TFileService serrvice
  art::ServiceHandle<art::TFileService> tfs;
  //Load FFT serrvice
  art::ServiceHandle<util::LArFFT> fft_service;
};


opdet::OpDeconvolutionAlgWiener::OpDeconvolutionAlgWiener(fhicl::ParameterSet const& p)
{
  //read fhicl paramters
  fDebug = p.get< bool >("Debug");
  fMaxFFTSizePow = p.get< int >("MaxFFTSizePow");
  fPositivePolarity = p.get< bool >("PositivePolarity");
  fUseSaturated = p.get< bool >("UseSaturated");
  fADCSaturationValue = p.get< int >("ADCSaturationValue");
  fApplyExpoAvSmooth   = p.get< bool >("ApplyExpoAvSmooth");
  fApplyUnAvSmooth   = p.get< bool >("ApplyUnAvSmooth");
  fExpoAvSmoothPar = p.get< float >("ExpoAvSmoothPar");
  fUnAvNeighbours = p.get< short unsigned int >("UnAvNeighbours");
  fHypoSignalTimeWindow = p.get< double >("HypoSignalTimeWindow");
  fHypoSignalCustom = p.get< bool >("HypoSignalCustom");
  fHypoSignalTimeConsts = p.get< std::vector<double> >("HypoSignalTimeConsts");
  fHypoSignalWeights = p.get< std::vector<double> >("HypoSignalWeights");
  fHypoSignalScale = p.get< double >("HypoSignalScale");
  fPMTChargeToADC = p.get< double >("PMTChargeToADC");
  fDecoWaveformPrecision = p.get< double >("DecoWaveformPrecision");
  fBaselineSample = p.get< short unsigned int >("BaselineSample");
  fSkipChannelList = p.get< std::vector<int>>("SkipChannelList");
  fFilter = p.get< std::string >("Filter");
  fElectronics = p.get< std::string >("Electronics");
  fDaphne_Freq  = p.get< double >("DaphneFreq");
  fScaleHypoSignal = p.get< bool >("ScaleHypoSignal");
  fUseParamFilter = p.get< bool >("UseParamFilter");
  fUseParamFilterInidividualChannel = p.get< bool >("UseParamFilterInidividualChannel");
  fCorrectBaselineOscillations = p.get< bool >("CorrectBaselineOscillations");
  fBaseSampleBins = p.get< short unsigned int >("BaseSampleBins");
  fBaseVarCut = p.get< double >("BaseVarCut");
  fFilterParams = p.get< std::vector<double> >("FilterParams");

  fFilterTF1 = new TF1("FilterTemplate", fFilter.c_str());

  fNormUnAvSmooth=1./(2*fUnAvNeighbours+1);
  NDecoWf=0;
  MaxBinsFFT=std::pow(2, fMaxFFTSizePow);

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
  fSamplingFreq=clockData.OpticalClock().Frequency()/1000.;//in GHz
  if (fElectronics=="Daphne") fSamplingFreq=fDaphne_Freq/1000.;//in GHz
  auto const* lar_prop = lar::providerFrom<detinfo::LArPropertiesService>();

  //Load PMT Calibration Database
  fPMTCalibrationDatabaseService = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();

  if (!fUseParamFilter){
    //Create light signal hypothesis for "on-fly" Wiener filter
    fSignalHypothesis.resize(MaxBinsFFT, 0);
    if(fFilter=="Wiener")
      fSignalHypothesis = ScintArrivalTimesShape(MaxBinsFFT, *lar_prop);
    else if(fFilter=="Wiener1PE")
      fSignalHypothesis[0]=1;
    else
      throw std::runtime_error("Deconvolution not chosen correctly ");
    mf::LogInfo("OpDeconvolutionAlg")<<"Built light signal hypothesis... L="<<fFilter<<" size"<<fSignalHypothesis.size()<<std::endl;
  }
}


std::vector<raw::OpDetWaveform> opdet::OpDeconvolutionAlgWiener::RunDeconvolution(std::vector<raw::OpDetWaveform> const& wfVector)
{
  std::vector<raw::OpDetWaveform> wfDeco;
  wfDeco.reserve(wfVector.size());
  for(auto const& wf : wfVector)
  {
    int channelNumber = wf.ChannelNumber();
    auto it = std::find(fSkipChannelList.begin(), fSkipChannelList.end(), channelNumber);
    if (it == fSkipChannelList.end()) {
      fSinglePEWave = fPMTCalibrationDatabaseService->getSER(channelNumber);
      double SPEAmplitude =  fPMTCalibrationDatabaseService->getSPEAmplitude(channelNumber);
      double SPEPeakValue = *std::max_element(fSinglePEWave.begin(), fSinglePEWave.end(), [](double a, double b) {return std::abs(a) < std::abs(b);});
      double SinglePENormalization = std::abs(SPEAmplitude/SPEPeakValue);
      std::transform(fSinglePEWave.begin(), fSinglePEWave.end(), fSinglePEWave.begin(), [SinglePENormalization](double val) {return val * SinglePENormalization;});
      fSinglePEWave.resize(MaxBinsFFT, 0);
      // If use channel dependent param filter
      if(fUseParamFilter)
      {
        double GaussFilterPower = fPMTCalibrationDatabaseService->getGaussFilterPower(channelNumber);
        double GaussFilterWC = fPMTCalibrationDatabaseService->getGaussFilterWC(channelNumber);
        fFilterParamChannel = {GaussFilterWC, GaussFilterPower};
        for(size_t k=0; k<fFilterParamChannel.size(); k++)
        {
          fFilterTF1->SetParameter(k, fFilterParamChannel[k]);
        }
          mf::LogInfo("OpDeconvolutionAlg")<<"Creating parametrized filter... TF1:"<<fFilter << " for channel " << channelNumber <<std::endl;
      }
    }
    //Read waveform
    size_t wfsize=wf.Waveform().size();
    if(wfsize>MaxBinsFFT){
      mf::LogWarning("OpDeconvolutionAlg")<<"Skipping waveform...waveform size is"<<wfsize<<"...maximum allowed FFT size is="<<MaxBinsFFT<<std::endl;
      continue;
    }
    size_t wfsizefft=WfSizeFFT(wfsize);

    std::vector<double> wave;
    wave.reserve(wfsizefft);
    wave.assign(wf.Waveform().begin(), wf.Waveform().end());

    //Get peak ADC value
    double wfPeakADC;
    bool saturated=false;
    if(fPositivePolarity) {
      wfPeakADC = *max_element(wave.begin(), wave.end());
      saturated = wfPeakADC>=fADCSaturationValue;
    }
    else{
      wfPeakADC = *min_element(wave.begin(), wave.end());
      saturated = wfPeakADC<=fADCSaturationValue;
    }

    if(!fUseSaturated && saturated ){
      mf::LogWarning("OpDeconvolutionAlg")<<"Skip saturated waveform @ OpCh "<< wf.ChannelNumber()<<" with time stamp "<<wf.TimeStamp()<<"\n";
      continue;
    }

    //Apply waveform smoothing
    if(fApplyExpoAvSmooth)
      ApplyExpoAvSmoothing(wave);
    if(fApplyUnAvSmooth)
      ApplyUnAvSmoothing(wave);

    //Estimate baseline standard deviation
    double baseline_mean=0., baseline_stddev=1.;
    EstimateBaselineStdDev(wave, baseline_mean, baseline_stddev);
    double wfPeakPE;
    if(fPositivePolarity) wfPeakPE = fHypoSignalScale*(wfPeakADC-baseline_mean)/fPMTChargeToADC;
    else wfPeakPE = fHypoSignalScale*(baseline_mean-wfPeakADC)/fPMTChargeToADC;
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
    if(fCorrectBaselineOscillations) CorrectBaselineOscillations(wave);

    //Debbuging and save wf in hist file
    if(fDebug){
      std::string name="h_deco"+std::to_string(NDecoWf)+"_"+std::to_string(wf.ChannelNumber())+"_"+std::to_string(wf.TimeStamp());
      TH1F * h_deco = tfs->make< TH1F >(name.c_str(),";Bin;#PE", MaxBinsFFT, 0, MaxBinsFFT);
      for(size_t k=0; k<wave.size(); k++){
        h_deco->Fill(k, wave[k]);
      }

      name="h_raw"+std::to_string(NDecoWf)+"_"+std::to_string(wf.ChannelNumber())+"_"+std::to_string(wf.TimeStamp());
      TH1F * h_raw = tfs->make< TH1F >(name.c_str(),";Bin;ADC", MaxBinsFFT, 0, MaxBinsFFT);
      for(size_t k=0; k<wf.Waveform().size(); k++){
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


void opdet::OpDeconvolutionAlgWiener::ApplyExpoAvSmoothing(std::vector<double>& wf){
  std::transform (std::next(wf.begin(), 1), wf.end(), wf.begin(), std::next(wf.begin(), 1),
    [&](double _x, double _y) { return  fExpoAvSmoothPar*_x+ (1. - fExpoAvSmoothPar)*_y; }  );
}


void opdet::OpDeconvolutionAlgWiener::ApplyUnAvSmoothing(std::vector<double>& wf){
  std::vector<double> wf_aux(wf.begin(), wf.end());
  for(size_t bin=fUnAvNeighbours; bin<wf.size()-fUnAvNeighbours; bin++){
    double sum=0.;
    for(size_t nbin=bin-fUnAvNeighbours; nbin<=bin+fUnAvNeighbours; nbin++)
      sum+=wf_aux[nbin];
    wf[bin]=sum*fNormUnAvSmooth;
  }
}


size_t opdet::OpDeconvolutionAlgWiener::WfSizeFFT(size_t n){
  if (n && !(n & (n - 1)))
       return n;
  size_t cont=0;
  while(n>0){
    cont++;
    n=(n>>1);
  }
  return 1<<cont;
}


std::vector<double> opdet::OpDeconvolutionAlgWiener::ScintArrivalTimesShape(size_t n, detinfo::LArProperties const& lar_prop){
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


void opdet::OpDeconvolutionAlgWiener::SubtractBaseline(std::vector<double> &wf, double baseline){
  std::transform (wf.begin(), wf.end(), wf.begin(), [&](double _x) { return _x-baseline; }  );
  return;
}


void opdet::OpDeconvolutionAlgWiener::EstimateBaselineStdDev(std::vector<double> &wf, double &_mean, double &_stddev){
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
  }

  _stddev=h_std.GetXaxis()->GetBinCenter(h_std.GetMaximumBin());
  _mean=h_mean.GetXaxis()->GetBinCenter(h_mean.GetMaximumBin());

  if(fDebug){

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


void opdet::OpDeconvolutionAlgWiener::CorrectBaselineOscillations(std::vector<double>& wave) 
{
  // subtract baseline using linear interpolation between regions defined
  // by the datasize and baseSampleBins
  if(fBaseSampleBins <= 0) return;
  // number of points to characterize the baseline
  unsigned short nBasePts = wave.size() / fBaseSampleBins;
  // the baseline offset vector
  std::vector<double> base(nBasePts, 0.);
  // find the average value in each region, using values that are
  // similar
  double fbins = fBaseSampleBins;
  unsigned short nfilld = 0;

  for(unsigned short ii = 0; ii < nBasePts; ++ii) {
    unsigned short loBin = ii * fBaseSampleBins;
    unsigned short hiBin = loBin + fBaseSampleBins;
    double ave = 0.;
    double sum = 0.;
    for(unsigned short bin = loBin; bin < hiBin; ++bin) {
      ave += wave[bin];
      sum += wave[bin] * wave[bin];
    } // bin
    ave = ave / fbins;
    double var = (sum - fbins * ave * ave) / (fbins - 1.);
    // Set the baseline for this region if the variance is small
    if(var < fBaseVarCut) {
      base[ii] = ave;
      ++nfilld;
    }
  } // ii
    // interpolate and subtract
    // find the first non-zero sample
  unsigned short thisRegion = 0;
  unsigned short nextRegion = 1;
  for(nextRegion = 1; nextRegion < base.size(); ++nextRegion) if(base[nextRegion] != 0) break;
  double nBins = (nextRegion - thisRegion) * fBaseSampleBins;
  double slp = (base[nextRegion] - base[thisRegion]) / nBins;
  for(unsigned short bin = 0; bin < wave.size(); ++bin) {
    // determine which region we are in
    short region = bin / fBaseSampleBins;
    if(region == nextRegion) {
      // update the slope for the next region
      unsigned short nnr = 0;
      for(nnr = nextRegion + 1; nnr < base.size(); ++nnr) if(base[nnr] != 0) break;
      if(nnr == base.size()) nnr = base.size() - 1;
      unsigned short previousRegion = thisRegion;
      thisRegion = nextRegion;
      nextRegion = nnr;
      if(thisRegion == nextRegion) {
        for (unsigned short kk = thisRegion * fBaseSampleBins; kk < thisRegion * fBaseSampleBins + fBaseSampleBins; ++ kk) {
          double bs = base[previousRegion] + (double)(kk - previousRegion * fBaseSampleBins) * slp;
            wave[kk] -= bs;
        }
        return;
        }
      nBins = (nextRegion - thisRegion) * fBaseSampleBins;
        slp = (base[nextRegion] - base[thisRegion]) / nBins;
    }
    double bs = base[thisRegion] + (double)(bin - thisRegion * fBaseSampleBins) * slp;
    wave[bin] -= bs;
  } // bin
  return;
}


std::vector<TComplex> opdet::OpDeconvolutionAlgWiener::DeconvolutionKernel(size_t wfsize, double baseline_stddev, double snr_scaling){
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

DEFINE_ART_CLASS_TOOL(opdet::OpDeconvolutionAlgWiener)
