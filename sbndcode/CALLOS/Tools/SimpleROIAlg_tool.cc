////////////////////////////////////////////////////////////////////////
// Specific class tool for ROIFindetAlg tool
// File: SimpleROI_tool.cc
// Base class:        ROIFINDERALG.hh
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "sbndcode/CALLOS/Tools/ROIFinderAlg.hh"

namespace callos {
  class SimpleROIAlg;
}


class callos::SimpleROIAlg : callos::ROIFINDERALG {
public:
  explicit SimpleROIAlg(fhicl::ParameterSet const& p);

  ~SimpleROIAlg() {}

  // Required functions.
  bool ProcessWaveform(std::vector<float> const& wvf ,std::vector<SimpleROI> & ROI) ;

private:

  // Declare member data here.
  bool fDebug;
  float fConstantBaseline;
  float fThresholdLow;
  float fThresholdHigh;
  float fThresholdUndershoot;
  int fBinsBefore;

  // Declare member functions
  void SubtractBaseline(std::vector<float> &wf, float baseline);
  void SubtractBaseline(std::vector<float> &wf, std::vector<float> &baseline);  //overload for non-constant baselines
  float EstimateBaseline(std::vector<double> &wf);

  // //Load TFileService serrvice
  // art::ServiceHandle<art::TFileService> tfs;
  // //Load FFT serrvice
  // art::ServiceHandle<util::LArFFT> fft_service;
};

callos::SimpleROIAlg::SimpleROIAlg(fhicl::ParameterSet const& p)
{
  fDebug = p.get< bool >("Debug",true);
  fConstantBaseline = p.get< float >( "ConstantBaseline", 0.0);
  fThresholdLow = p.get< float >( "ThresholdLow", 0.0);
  fThresholdHigh = p.get< float >( "ThresholdHigh", 0.0);
  fBinsBefore = p.get< int >( "BinsBefore", 0);
  fThresholdUndershoot = p.get< float >( "ThresholdUndershoot", -4);

  if (fDebug) 
  {
    std::cout << "SimpleROIAlg fhicl congif::"<<std::endl;
    std::cout << "ConstantBaseline: "<<fConstantBaseline<<std::endl;
    std::cout << "ThresholdLow: "<<fThresholdLow<<std::endl;
    std::cout << "ThresholdHigh: "<<fThresholdHigh<<std::endl;
  }

}

// constant baseline substractor
void callos::SimpleROIAlg::SubtractBaseline(std::vector<float> &wf, float baseline){
  std::transform (wf.begin(), wf.end(), wf.begin(), [&](float _x) { return _x-baseline; }  );
  return;
}

// // variable baseline subtractor 
// void callos::SimpleROIAlg::SubtractBaseline(std::vector<float> &wf, std::vector<float> &baseline){
//   std::transform (wf.begin(), wf.end(), baseline.begin(), wf.begin(), std::minus<float>());
//   return;
// }

bool callos::SimpleROIAlg::ProcessWaveform(std::vector<float> const& wvf ,std::vector<SimpleROI> & ROI)
{
  // constant for now, will update 
  std::vector<float> wave;
  // look at the second half of the waveform (scintillation tail, after the prompt light is gone)
  int half=wvf.size()/2;
  wave.reserve(half);
  wave.assign(wvf.begin()+half, wvf.end());

  SubtractBaseline(wave,fConstantBaseline);
  auto wvf_max = *std::max_element(wave.begin(), wave.end());
  auto wvf_min = *std::min_element(wave.begin(), wave.end());

  // std::cout << "SimpleROIAlg::ProcessWaveform: max of the waveform is "<<wvf_max<<std::endl;
  //check max of the waveform is below ThresholdHigh
  if (wvf_max > fThresholdHigh ) return false;
  if (wvf_min < fThresholdUndershoot ) return false;

  // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: passed max threshold requisite"<<std::endl;
  // Signal could be empty or contain single PEs. Read the wvf, look for  peaks above a threshold,
  // compute the charge untill the waveform is below fConstantBaseline again to left and right
  // store the ROI in the ROI vector

  for (size_t i = 0; i < wave.size(); i++) {
    

    if (wave[i] > fThresholdLow) { //found a peak
      SimpleROI roi(10, 1);//this should be updated with something more meaningful later on
      roi.SetCharge(0);
      float aux_charge = 0;
      int j = i;
      // Add the bins after the threshold until the waveform crosses the baseline
      while ( (wave[j] > 0) && (j< (int) wave.size()) ) {
        // roi[j-i] = wave[j];
        aux_charge += wave[j];
        j++;
      }
      // Add the bins before the threshold
      for (int k = 0; k < fBinsBefore; k++) {
        if ((i-k) < 0) break;

        if (wave[i-k] > fConstantBaseline) {
          // roi[j-i+k] = wave[i-k];
          aux_charge += wave[i-k];
        }
      }
      if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: found a peak at "<<i<<" with charge "<<aux_charge<<std::endl;
      roi.SetCharge(aux_charge);
      ROI.push_back(roi);
      i = j;//skip the waveform, find the next peak
    }
  }
  
  //  if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform"<<std::endl;
  return true;

}

DEFINE_ART_CLASS_TOOL(callos::SimpleROIAlg)

