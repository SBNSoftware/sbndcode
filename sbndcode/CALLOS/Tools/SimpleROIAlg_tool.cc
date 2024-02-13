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
  bool ProcessWaveform(std::vector<float> const& wvf ,std::vector<float> & ROI, std::vector<float> & charges) ;

private:

  // Declare member data here.
  bool fDebug;

  // Declare member functions
  void SubtractBaseline(std::vector<double> &wf, double baseline);
  float EstimateBaseline(std::vector<double> &wf);

  // //Load TFileService serrvice
  // art::ServiceHandle<art::TFileService> tfs;
  // //Load FFT serrvice
  // art::ServiceHandle<util::LArFFT> fft_service;
};

callos::SimpleROIAlg::SimpleROIAlg(fhicl::ParameterSet const& p)
{
  fDebug = p.get< bool >("Debug",true);
}

void callos::SimpleROIAlg::SubtractBaseline(std::vector<double> &wf, double baseline){
  std::transform (wf.begin(), wf.end(), wf.begin(), [&](double _x) { return _x-baseline; }  );
  return;
}

bool callos::SimpleROIAlg::ProcessWaveform(std::vector<float> const& wvf ,std::vector<float> & ROI, std::vector<float> & charges)
{
 if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform"<<std::endl;
  return true;

}

DEFINE_ART_CLASS_TOOL(callos::SimpleROIAlg)

