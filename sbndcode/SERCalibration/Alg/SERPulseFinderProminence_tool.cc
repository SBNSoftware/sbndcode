
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

#include "SERPulseFinderBase.hh"


namespace opdet {
  class SERPulseFinderProminence;
}

class opdet::SERPulseFinderProminence : opdet::SERPulseFinderBase{

public:

    explicit  SERPulseFinderProminence(fhicl::ParameterSet const& p);
    
    virtual void RunSERCalibration(std::vector<raw::OpDetWaveform> const&  , std::vector<TH1D>& ) override;

    void find_peaks(std::vector<double> * , int , double );
    double calculateProminence(std::vector<double> * data, int peakIndex);
    
    void setProminence(double Prominence){ fProminence = Prominence; return;}
    void setWindowLength(int WindowLenght){ fWindowLenght = WindowLenght; return;}
    void setPeakWidhtL(int PeakWidthL){ fPeakWidthL = PeakWidthL; return;}
    void setPeakWidhtR(int PeakWidthR){ fPeakWidthR = PeakWidthR; return;}
    void SubstractBaseline();
    void CalculateFFT();
    void AddWaveforms(TH1D& , int);
    bool IsIsolatedPeak(int );
    int GetPeakValue(int );
    std::vector<std::vector<int>> GetPeakAmplitudeVector() override {return fPeakAmplitudeVector;}
    std::vector<int> GetNumberOfPeaksVector() override {return fNumberOfPeaksVector;}

    void GetTriggerIdx();
    void GetTriggerIdxFFT();

private:

    int fWindowLenght;
    int fSmoothingWidth;
    int fPeakWidthL;
    int fPeakWidthR;
    int fProminence;
    size_t fDecreasingInterval;
    int fWaveformMax;
    int fProminenceLeftLimit;
    int fProminenceRightLimit;
    int fProminenceLowerBound;
    int fBaselineSample;
    int fLowerSERAmplitude;
    int fHigherSERAmplitude;
    bool fUseFrequencyFilter;
    size_t fMaxNumPeaks;
    TH1D * TimeSignal;
    TH1D * CumulativeTimeSignal;

    TH1D* fInputWaveform= nullptr;
    TH1D* fInputWaveformFFT = nullptr;
    std::vector<int> TriggerIdx;
    std::vector<std::vector<int>> fPeakAmplitudeVector;

    int fSERStart;
    int fSEREnd;

    std::vector<int> fNumberOfPeaksVector;

    //Load FFT serrvice
    art::ServiceHandle<util::LArFFT> fft_service;


};

opdet::SERPulseFinderProminence::SERPulseFinderProminence(fhicl::ParameterSet const& p) 
{
    fSmoothingWidth = p.get< int >("SmoothingWidth");
    fPeakWidthL = p.get< int >("PeakWidthL");
    fPeakWidthR = p.get< int >("PeakWidthR");
    fProminence = p.get< int >("Prominence");
    fDecreasingInterval = p.get< size_t >("DecreasingInterval", 0);
    fProminenceLeftLimit = p.get< int >("ProminenceLeftLimit");
    fProminenceRightLimit = p.get< int >("ProminenceRightLimit");
    fProminenceLowerBound = p.get< int >("ProminenceLowerBound", -9999);
    fBaselineSample = p.get< int >("BaselineSample");
    fSERStart = p.get< int >("SERStart");
    fSEREnd = p.get< int >("SEREnd");
    fLowerSERAmplitude = p.get< int >("LowerSERAmplitude");
    fHigherSERAmplitude = p.get< int >("HigherSERAmplitude");
    fMaxNumPeaks = p.get< size_t >("MaxNumPeaks");
    fUseFrequencyFilter = p.get< bool >("UseFrequencyFilter");
    fNumberOfPeaksVector.resize(320);
    fPeakAmplitudeVector.resize(320);
}

void opdet::SERPulseFinderProminence::GetTriggerIdx()
{
    Int_t b_max = fInputWaveform->GetMaximumBin();
    fWaveformMax = fInputWaveform->GetBinContent(b_max);

    //Create the auxiliary vector to pass to find_peaks
    std::vector<double> * y = new std::vector<double>();
    for(int i=0; i<fInputWaveform->GetNbinsX(); i++)
    {
        y->push_back(fInputWaveform->GetBinContent(i));
    }

    find_peaks(y, fSmoothingWidth, fProminence);

    delete y;
}

void opdet::SERPulseFinderProminence::RunSERCalibration(std::vector<raw::OpDetWaveform> const& wfVector , std::vector<TH1D>& calibratedSER_v)
{
    std::cout << " Entering the run SER calibration with a wf vector of size " << wfVector.size() << std::endl;
    //Loop over all waveforms and find the peaks 
    for(auto const& wf : wfVector)
    {
        int ChNumber = wf.ChannelNumber();
        fInputWaveform = new TH1D("Input WForm", "" , wf.size(), 0, double(wf.size()));

        for(size_t i = 0; i < wf.size(); i++) {
            fInputWaveform->SetBinContent(i + 1, (double)wf[i]);
        }
        // Substract waveform baseline
        SubstractBaseline();
        // Calculate waveform FFT
        if(fUseFrequencyFilter)
        {
            CalculateFFT();
            GetTriggerIdxFFT();
        }
        else GetTriggerIdx();

        // Loop over the found peaks to see if we add it to the waveform
        for(size_t k=0; k<TriggerIdx.size(); k++)
        {
            int PeakValue = this->GetPeakValue(TriggerIdx[k]);
            bool isIsolated = this->IsIsolatedPeak(TriggerIdx[k]);
            fPeakAmplitudeVector.at(ChNumber).push_back(PeakValue);
            if(PeakValue>fHigherSERAmplitude ||PeakValue<fLowerSERAmplitude) continue;
            if(TriggerIdx[k]>4900 || TriggerIdx[k]<100) continue; // Substitue by function that tells wether it is within waveform limits 
            if(TriggerIdx.size()>fMaxNumPeaks || !isIsolated) continue;
            AddWaveforms(calibratedSER_v.at(ChNumber), TriggerIdx[k]);
            fNumberOfPeaksVector.at(ChNumber)+=1;
            //myWaveformPeakFinder.PlotPeak(PeakIdx->at(k));
        }
        delete fInputWaveform;
        if (fInputWaveformFFT != nullptr) {
            delete fInputWaveformFFT;
            fInputWaveformFFT = nullptr;  // No es estrictamente necesario, pero es buena práctica
        }
        TriggerIdx.clear();
    }
}


void opdet::SERPulseFinderProminence::GetTriggerIdxFFT()
{
    Int_t b_max = fInputWaveformFFT->GetMaximumBin();
    fWaveformMax = fInputWaveformFFT->GetBinContent(b_max);

    //Create the auxiliary vector to pass to find_peaks
    std::vector<double> * y = new std::vector<double>();
    for(int i=0; i<fInputWaveformFFT->GetNbinsX(); i++)
    {
        y->push_back(fInputWaveformFFT->GetBinContent(i));
    }

    find_peaks(y, fSmoothingWidth, fProminence);

    delete y;
}


void opdet::SERPulseFinderProminence::find_peaks(std::vector<double> * y, int width, double prominence) {
        
    
    int MaxWidth = std::max(fPeakWidthL, fPeakWidthR);

    for (size_t i = MaxWidth; i < y->size() - MaxWidth; ++i) {
        bool is_peak = true;
        // Check that i is a peak
        for (size_t j = i - fPeakWidthL; j <= i + fPeakWidthR; j++) {
            if(j != i && y->at(i) < y->at(j)) {
                is_peak = false;
                break;
            }
        }

        //Check also that the bins surrounding the peak are monotonically increasing/descreasing
        if (is_peak) {
            for (size_t j = i - fDecreasingInterval; j <= i + fDecreasingInterval; j++) {
                if ( j<i && y->at(j) > y->at(j + 1)) { // If the bin is to the left it should increase
                    is_peak = false;
                    break;
                }
                else if( j>i && y->at(j) < y->at(j + 1)) { // If the bin is to the right it should decrease
                    is_peak = false;
                    break;
                }
            }
        }

        // Check that the prominene is high enough
        if (is_peak) {
            double current_prominence = calculateProminence(y,  i);
            if(current_prominence>=prominence) 
            {
                TriggerIdx.push_back(i);
            }
        }
    }
}

double opdet::SERPulseFinderProminence::calculateProminence(std::vector<double> * data, int peakIndex) {
   
    int dataSize = data->size();
    // Verificar si el índice del pico está dentro de los límites del vector
    if (peakIndex < 0 || peakIndex >= dataSize) {
        std::cerr << "Índice de pico fuera de los límites." << std::endl;
        return 0.0;
    }

    int leftIndex = peakIndex - 1;
    int rightIndex = peakIndex + 1;
    while (leftIndex > 0 && data->at(leftIndex) < data->at(peakIndex) && (leftIndex>=peakIndex-fProminenceLeftLimit)) {
        leftIndex--;
    }

    while (rightIndex < dataSize && data->at(rightIndex) < data->at(peakIndex) && (rightIndex>= peakIndex+fProminenceRightLimit)) {
        rightIndex++;
    }


    //Find the min of the histogram from the peakIndex to the left
    double peakHeight = data->at(peakIndex);
    double leftMin = peakHeight;
    for(int i=leftIndex; i<peakIndex; i++)
    {
        if(data->at(i)<leftMin && data->at(i)!=0) 
        {
            leftMin=data->at(i);
            leftIndex=i;
        }
    }

    //Find the min of the histogram from the peakIndex to the right
    double rightMin = peakHeight;
    for(int i=peakIndex; i<rightIndex; i++)
    {
        if(data->at(i)<rightMin && data->at(i)!=0)
        {
            rightMin=data->at(i);
            rightIndex=i;
        }
    }

    double prominence = peakHeight - std::min(leftMin, rightMin);

    if(std::min(leftMin, rightMin)<fProminenceLowerBound) return 0;
    else return prominence;
}

/*
void opdet::SERPulseFinderProminence::CalculateFFT()
{
    int N = fInputWaveform->GetNbinsX();

    TComplex wforminit(0,0,false);
    std::vector<TComplex> wformfft(N, wforminit);
    fft_service->ReinitializeFFT (N, "", 20);

    std::vector<double> wform;
    for(int i=0; i<N; i++)
    {
        wform.push_back(fInputWaveform->GetBinContent(i));
    }
   
    wformfft.resize(N);
    fft_service->DoFFT(wform, wformfft);

    int lowerFreqCut=1500;

    for(int i=lowerFreqCut; i<N; i++)
    {
        wformfft[i]=0;
    }

    for(int i=N/2; i< N-lowerFreqCut; i++)
    {
        wformfft[i]=0;
    }

    std::vector<double> filtered_wform;
    filtered_wform.resize(N);
    fft_service->DoInvFFT(wformfft, filtered_wform);

    fInputWaveformFFT = new TH1D("", "", N, 0, N);
    for(int i=0; i<N; i++)
    {
        fInputWaveformFFT->SetBinContent(i, filtered_wform[i]);
    }

}
*/


void opdet::SERPulseFinderProminence::CalculateFFT()
{
    int N = fInputWaveform->GetNbinsX();

    double* signal_array = new double[N];
    for (int i = 0; i < N; i++) {
        signal_array[i] = fInputWaveform->GetBinContent(i + 1);
    }

    TH1 * h_fft = 0;
    TVirtualFFT::SetTransform(0);
    h_fft = fInputWaveform->FFT(h_fft, "MAG");
    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    Double_t *re_full = new Double_t[N];
    Double_t *im_full = new Double_t[N];

    fft->GetPointsComplex(re_full,im_full);

    int lowerFreqCut=1500;
    for(int i=lowerFreqCut; i<N; i++)
    {
        re_full[i] = 0;
        im_full[i] = 0;
    }

    for(int i=N/2; i< N-lowerFreqCut; i++)
    {
        re_full[i] = 0;
        im_full[i] = 0;
    }

    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R M K");
    fft_back->SetPointsComplex(re_full, im_full);
    fft_back->Transform();

    fInputWaveformFFT= new TH1D("", "", N , 0 , N);
    fInputWaveformFFT = dynamic_cast<TH1D*>(TH1D::TransformHisto(fft_back, fInputWaveformFFT, "Re"));



    for(int i=0; i<N; i++)
    {
        fInputWaveformFFT->SetBinContent(i+1, fInputWaveformFFT->GetBinContent(i+1)/N);
    }

    delete fft_back;
    delete fft;
    delete h_fft;
    delete [] re_full;
    delete [] im_full;
    delete [] signal_array;
}


void opdet::SERPulseFinderProminence::AddWaveforms(TH1D& SERWaveform, int PeakTime)
{
    int SERLength = fSEREnd-fSERStart;
    for(int j=0; j<SERLength; j++)
    {
        if(PeakTime-int(SERLength/2)<0) break;
        double updatedValue = SERWaveform.GetBinContent(j)+fInputWaveform->GetBinContent(PeakTime-int(SERLength/2)+j);
        //std::cout << " Setting content of the calibrated ser in bin " << j << " to " << updatedValue << std::endl;
        SERWaveform.SetBinContent(j, updatedValue);
    }
}

int opdet::SERPulseFinderProminence::GetPeakValue(int PeakIdx)
{
    return fInputWaveform->GetBinContent(PeakIdx);
}

bool opdet::SERPulseFinderProminence::IsIsolatedPeak(int PeakIdx)
{
    int TimeDifferenceNext;
    int TimeDifferencePrevious;
    int NextTriggerIdx;
    int PreviousTriggerIdx;
    unsigned int PeakPositionOnVector = std::distance(TriggerIdx.begin(), std::find(TriggerIdx.begin(), TriggerIdx.end(), PeakIdx));
    
    if(PeakPositionOnVector<TriggerIdx.size()-1)
    {
        NextTriggerIdx = TriggerIdx.at(PeakPositionOnVector+1);
    }
    if(PeakPositionOnVector>0)
    {
        PreviousTriggerIdx = TriggerIdx.at(PeakPositionOnVector-1);
    }
    
    TimeDifferenceNext = abs(NextTriggerIdx - PeakIdx);
    TimeDifferencePrevious = abs(PeakIdx - PreviousTriggerIdx);

    if(std::min(TimeDifferenceNext, TimeDifferencePrevious)>550) return true;
    else return false;

    return false;
}

void opdet::SERPulseFinderProminence::SubstractBaseline(){

    double minADC= fInputWaveform->GetMinimum();
    double maxADC=fInputWaveform->GetMaximum();
    unsigned nbins=25*ceil(maxADC-minADC);

    if(nbins==0) nbins=100;

    TH1F h_mean = TH1F("",";;", nbins, minADC, maxADC);

    for(int ix=0; ix<fInputWaveform->GetNbinsX()-fBaselineSample; ix++){
        double sum=0;
        for(int jx=ix; jx<ix+fBaselineSample; jx++){
        sum = sum + fInputWaveform->GetBinContent(jx);
        }
        sum/=fBaselineSample;
        h_mean.Fill( sum );
    }

    double baseline_mean=h_mean.GetXaxis()->GetBinCenter(h_mean.GetMaximumBin());

    for(int i=1; i<=fInputWaveform->GetNbinsX(); i++)
    {
        fInputWaveform->SetBinContent(i, -(fInputWaveform->GetBinContent(i)-baseline_mean));
    }

    return;
}

DEFINE_ART_CLASS_TOOL(opdet::SERPulseFinderProminence)