#ifndef SERPULSEFINDER_H
#define SERPULSEFINDER_H

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

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>
#include "TVector2.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"


namespace opdet {
  class SERPulseFinder;
}

class opdet::SERPulseFinder{
public:

    SERPulseFinder(fhicl::ParameterSet const& p);
    ~SERPulseFinder();

    void RunSERCalibration(std::vector<raw::OpDetWaveform> const&  , std::vector<TH1D> * );
    std::vector<int>* find_peaks(std::vector<double> * , int , double );
    double calculateProminence(std::vector<double> * data, int peakIndex);
    
    void setProminence(double Prominence){ fProminence = Prominence; return;}
    void setWindowLength(int WindowLenght){ fWindowLenght = WindowLenght; return;}
    void setWindowInitialTime(int WindowInitialTime){ fWindowInitialTime = WindowInitialTime; return;}
    void setPeakWidhtL(int PeakWidthL){ fPeakWidthL = PeakWidthL; return;}
    void setPeakWidhtR(int PeakWidthR){ fPeakWidthR = PeakWidthR; return;}
    void setInputWaveform(TH1D* InputWaveform){ fInputWaveform = InputWaveform; return;}
    void SubstractBaseline();
    void CalculateFFT();
    void AddWaveforms(TH1D& , int);
    bool IsIsolatedPeak(int );
    int GetPeakValue(int );

    std::vector<int>* GetTriggerIdx();
    std::vector<int>* GetTriggerIdxFFT();

private:

    int fWindowLenght;
    int fSmoothingWidth;
    int fPeakWidthL;
    int fPeakWidthR;
    int fProminence;
    int fWindowInitialTime;
    size_t fDecreasingInterval;
    int fWaveformMax;
    int fProminenceLeftLimit;
    int fProminenceRightLimit;
    int fSummedWaveformLength;
    int fProminenceLowerBound;
    int fBaselineSample;
    
    double fTriggerTimeDifference;

    TH1D * TimeSignal;
    TH1D * CumulativeTimeSignal;

    TH1D* fInputWaveform;
    TH1D* fInputWaveformFFT;
    std::vector<int> * TriggerIdx;

    int fSERStart;
    int fSEREnd;

    std::vector<int> DetectedPeaksPerChannel;
};

#endif