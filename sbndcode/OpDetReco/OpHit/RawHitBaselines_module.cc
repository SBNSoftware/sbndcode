////////////////////////////////////////////////////////////////////////
// Class:       RawHitBaselines
// Plugin Type: producer (Unknown Unknown)
// File:        RawHitBaselines_module.cc
//
// Generated at Tue Oct 22 09:27:52 2024 by Jacob McLaughlin using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////
//Must be run after summed waveform production
//Takes in OpDetWaveform from PMTs and performs a running average and guassian smoothing
//Subtracts baseline from OpDetWaveform and saves the baseline-subtracted copy for 
//raw hit finding modules

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"

#include <memory>
#include <numeric>
#include <chrono>


namespace opdet {
  class RawHitBaselines;
}


class opdet::RawHitBaselines : public art::EDProducer {
public:
  explicit RawHitBaselines(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RawHitBaselines(RawHitBaselines const&) = delete;
  RawHitBaselines(RawHitBaselines&&) = delete;
  RawHitBaselines& operator=(RawHitBaselines const&) = delete;
  RawHitBaselines& operator=(RawHitBaselines&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void CalcRunningAvg(std::vector<int> &Baseline, raw::OpDetWaveform wvf);
  void GaussianSmoothing(std::vector<int> &Baseline);
  double GetWaveformMedian(raw::OpDetWaveform wvf);
private:
  int fRunningAvgSampleWidth;
  int fGuassianConvlSize;
  double fGaussianConvlWidth;
  std::string fWaveformName;
  std::string fch_instance_name;
  // Declare member data here.

};


opdet::RawHitBaselines::RawHitBaselines(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fRunningAvgSampleWidth = p.get< int >("RunningAvgWidth", 50 );
  fGuassianConvlSize = p.get< int >("GaussianConvlSize", 500 );
  fGaussianConvlWidth = p.get< double >("GaussianConvlWidth", 100 );
  fWaveformName = p.get< std::string >("ChannelWaveformName" );
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fch_instance_name = p.get<std::string>("pmtBackSubChannel","pmtBackSubChannel");
  produces< std::vector< raw::OpDetWaveform > >(fch_instance_name);
}

void opdet::RawHitBaselines::produce(art::Event& e)
{
  //Load in summed waveform and channel waveform vectors
  std::cout << "My module on event #" << e.id().event() << std::endl;
  art::Handle< std::vector< raw::OpDetWaveform > > chWaveHandle; //User handle for vector of OpDetWaveforms
  e.getByLabel(fWaveformName, chWaveHandle);
  //Make holders for the baseline subtracted waveform
  std::unique_ptr< std::vector< raw::OpDetWaveform > > channelWvfmVec(std::make_unique< std::vector< raw::OpDetWaveform > > ());
  //Now we can generate a baseline for each waveform
  int Index=0;
  for(auto const& wvf : (*chWaveHandle)) 
  {
    Index=Index+1;
    //auto start = std::chrono::high_resolution_clock::now();
    std::vector<uint16_t> chWvfms( wvf.size() );
    std::vector<int> Baseline( wvf.size() );
    CalcRunningAvg( Baseline,  wvf);
    //auto TimeRunningAvg = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
    //Do smoothing
    //GaussianSmoothing(Baseline);
    //auto TimeGauss =  std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() - TimeRunningAvg;
    //Record baseline subtracted waveform
    for(int i=0; i<int(wvf.size()); i++) chWvfms[i] = wvf[i] - Baseline[i];
    //Convert chWvfms object into a raw::opdetwaveform
    raw::TimeStamp_t tempTimeStamp = wvf.TimeStamp();
    unsigned int tempChannel = wvf.ChannelNumber();
    channelWvfmVec->push_back(raw::OpDetWaveform(tempTimeStamp, tempChannel, chWvfms  ) ); 
    //auto TimeCopy =  std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() - TimeGauss;
    //std::cout << "Time for running avg " <<  TimeRunningAvg << " Time for gauss " <<  TimeGauss << " Time for copy " << TimeCopy <<std::endl;
  }
  e.put(std::move(channelWvfmVec),fch_instance_name); 
  
}

void opdet::RawHitBaselines::CalcRunningAvg(std::vector<int> &Baseline, raw::OpDetWaveform wvf)
{
  //int index =0;
  double Median = GetWaveformMedian(wvf);
  //Bury a nested loop in a thing Im not sure is any faster
  for(int i =0 ; i<int(wvf.size())-fRunningAvgSampleWidth; i++)
  {
    int EndIndex = i+fRunningAvgSampleWidth;
    double sum=0;
    for(int j=i; j<EndIndex; j++)
    {
      if(wvf[j] > Median-40) sum = sum + wvf[j]; //mask out pulses
      else sum = sum + Median;
    }
    Baseline[i] = sum/(EndIndex - i);
  }
  //std::for_each(Baseline.begin(), Baseline.begin()+(wvf.size()-fRunningAvgSampleWidth), [wvf, &index, this](int& Val) 
  //          { 
  //          Val = std::accumulate(wvf.begin()+index, wvf.begin()+index+fRunningAvgSampleWidth, 0)/fRunningAvgSampleWidth;
  //          index = index+1;
  //          } );
  for(int i =int(wvf.size())-fRunningAvgSampleWidth ; i<int(wvf.size()); i++)
  {
    int EndIndex = i+fRunningAvgSampleWidth;
    if(EndIndex>int(wvf.size())) EndIndex = int(wvf.size());
    double sum=0;
    for(int j=i; j<EndIndex; j++)
    {
      if(wvf[j] > Median-40) sum = sum + wvf[j]; //mask out pulses
      else sum = sum + Median;
    }
    Baseline[i] = sum/(EndIndex - i);
  }
  //Baseline is all updated an can return
}

double opdet::RawHitBaselines::GetWaveformMedian(raw::OpDetWaveform wvf)
  {
    std::sort(wvf.begin(), wvf.end());
    int MedianIndex = int(wvf.size()/2);
    return wvf[MedianIndex];
  }

void opdet::RawHitBaselines::GaussianSmoothing(std::vector<int> &Baseline)
{
  std::vector<int> Out(Baseline.size());
  std::vector<int> X(fGuassianConvlSize*2+1);
  std::iota(X.begin(), X.end(), -fGuassianConvlSize);
  std::vector<double> Weights(fGuassianConvlSize*2+1);
  double sum=0;
  for(int i=0; i<int(Weights.size()); i++)
  {
    Weights[i] = TMath::Exp( - TMath::Power(double(X[i]), 2.0) / (2*TMath::Power(double(fGaussianConvlWidth), 2.0)) );
    sum += Weights[i];
  }
  //Now do the convolution
  for(int i=fGuassianConvlSize+1; i<int(Baseline.size()) - (fGuassianConvlSize+1); i++)
  {
    double PointSum=0;
    std::for_each(X.begin(), X.end(), [&PointSum, i, Weights, Baseline, this](int Index){ PointSum += Baseline[i+Index]*Weights[fGuassianConvlSize + Index]; });
    Out[i] = PointSum;
  }
  //Handle edges properly
  for(int i=0; i<fGuassianConvlSize+1; i++)
  {
    double PointSum=0;
    std::for_each(X.begin()+fGuassianConvlSize-i, X.end(), [&PointSum, i, Weights, Baseline, this](int Index){ PointSum += Baseline[i+Index]*Weights[fGuassianConvlSize + Index]; });
    Out[i] = PointSum;
  }
  for(int i=int(Baseline.size()) - (fGuassianConvlSize+1); i<int(Baseline.size()); i++)
  {
    double PointSum=0;
    std::for_each(X.begin(), X.begin()+fGuassianConvlSize-(i-int(Baseline.size())), [&PointSum, i, Weights, Baseline, this](int Index){ PointSum += Baseline[i+Index]*Weights[fGuassianConvlSize + Index]; });
    Out[i] = PointSum;
  }
  //Finally copy out to baseline
  for(int i=0; i<int(Baseline.size()); i++) Baseline[i] = Out[i]/sum;
}
DEFINE_ART_MODULE(opdet::RawHitBaselines)

