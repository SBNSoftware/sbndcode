#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>

#include "sbndcode/SERCalibration/Alg/SERPulseFinder_tool.h"


opdet::SERPulseFinder::SERPulseFinder(fhicl::ParameterSet const& p) 
{
    fSmoothingWidth = 20;
    fPeakWidthL=10;
    fPeakWidthR=10;
    fProminence = 21;
    //fProminence = 18;
    fWindowInitialTime=0;
    fTriggerTimeDifference=0.;
    TriggerIdx = nullptr;
    fInputWaveformFFT=0;
    fDecreasingInterval=0;
    fProminenceLeftLimit=200; // Prev. 50
    fProminenceRightLimit=50; // Prev. 50
    fSummedWaveformLength=1000;
    fProminenceLowerBound=-99999;
    fBaselineSample=15;
    fSERStart = -200;
    fSEREnd = 300;
    DetectedPeaksPerChannel.resize(320);
}

opdet::SERPulseFinder::~SERPulseFinder()
{}

std::vector<int>* opdet::SERPulseFinder::GetTriggerIdx()
{

    Int_t b_max = fInputWaveform->GetMaximumBin();
    fWaveformMax = fInputWaveform->GetBinContent(b_max);

    //Create the auxiliary vector to pass to find_peaks
    std::vector<double> * y = new std::vector<double>();
    for(int i=0; i<fInputWaveform->GetNbinsX(); i++)
    {
        y->push_back(fInputWaveform->GetBinContent(i));
    }

    TriggerIdx = find_peaks(y, fSmoothingWidth, fProminence);

    delete y;
    return TriggerIdx;
}


void opdet::SERPulseFinder::RunSERCalibration(std::vector<raw::OpDetWaveform> const& wfVector , std::vector<TH1D> * calibratedSER_v)
{
    //Loop over all waveforms and find the peaks 
    for(auto const& wf : wfVector)
    {
        int ChNumber = wf.ChannelNumber();

        TH1D* currentWform = new TH1D("Current WForm", "" , wf.size(), 0, double(wf.size()));

        for(size_t i = 0; i < wf.size(); i++) {
            currentWform->SetBinContent(i + 1, (double)wf[i]);
        }

        //Set the inpuf waveform for the pulse finder
        setInputWaveform(currentWform);
        // Substract waveform baseline
        SubstractBaseline();
        // Calculate waveform FFT
        CalculateFFT();
        std::vector<int> *PeakIdx = GetTriggerIdxFFT();
        // Loop over the found peaks to see if we add it to the waveform
        for(size_t k=0; k<PeakIdx->size(); k++)
        {
            int PeakValue = this->GetPeakValue(PeakIdx->at(k));
            bool isIsolated = this->IsIsolatedPeak(PeakIdx->at(k));
            if(PeakValue>40 ||PeakValue<7) continue;
            if(PeakIdx->at(k)>4900 || PeakIdx->at(k)<100) continue;
            if(PeakIdx->size()>40 || !isIsolated) continue;
            AddWaveforms(calibratedSER_v->at(ChNumber), PeakIdx->at(k));
            DetectedPeaksPerChannel.at(ChNumber)+=1;
            //myWaveformPeakFinder.PlotPeak(PeakIdx->at(k));
        }
    }
}


std::vector<int>* opdet::SERPulseFinder::GetTriggerIdxFFT()
{
    Int_t b_max = fInputWaveformFFT->GetMaximumBin();
    fWaveformMax = fInputWaveformFFT->GetBinContent(b_max);

    //Create the auxiliary vector to pass to find_peaks
    std::vector<double> * y = new std::vector<double>();
    for(int i=0; i<fInputWaveformFFT->GetNbinsX(); i++)
    {
        y->push_back(fInputWaveformFFT->GetBinContent(i));
    }

    TriggerIdx = find_peaks(y, fSmoothingWidth, fProminence);

    delete y;
    return TriggerIdx;
}


std::vector<int>* opdet::SERPulseFinder::find_peaks(std::vector<double> * y, int width, double prominence) {
        
    std::vector<int> * peaks = new std::vector<int>();
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
                peaks->push_back(i);
            }
        }
    }

    return peaks;
}

double opdet::SERPulseFinder::calculateProminence(std::vector<double> * data, int peakIndex) {
   
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

    //return prominence;

}

void opdet::SERPulseFinder::CalculateFFT()
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

    int lowerFreqCut=2000;
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
    fInputWaveformFFT = dynamic_cast<TH1D*>(TH1D::TransformHisto(fft_back, fInputWaveformFFT, "Re"));
    for(int i=0; i<N; i++)
    {
        fInputWaveformFFT->SetBinContent(i+1, fInputWaveformFFT->GetBinContent(i+1)/N);
    }

    fft_back=0;
    delete fft_back;
    delete [] re_full;
    delete [] im_full;
    delete [] signal_array;
    delete h_fft;
}


void opdet::SERPulseFinder::AddWaveforms(TH1D& SERWaveform, int PeakTime)
{
    int SERLength = fSEREnd-fSERStart;
    for(int j=0; j<SERLength; j++)
    {
        if(PeakTime-int(SERLength/2)<0) break;
        double updatedValue = SERWaveform.GetBinContent(j)+fInputWaveform->GetBinContent(PeakTime-int(SERLength/2)+j);
        SERWaveform.SetBinContent(j, updatedValue);
    }
}

int opdet::SERPulseFinder::GetPeakValue(int PeakIdx)
{
    return fInputWaveform->GetBinContent(PeakIdx);
}

bool opdet::SERPulseFinder::IsIsolatedPeak(int PeakIdx)
{
    int TimeDifferenceNext;
    int TimeDifferencePrevious;
    int NextTriggerIdx;
    int PreviousTriggerIdx;
    unsigned int PeakPositionOnVector = std::distance(TriggerIdx->begin(), std::find(TriggerIdx->begin(), TriggerIdx->end(), PeakIdx));
    
    if(PeakPositionOnVector<TriggerIdx->size()-1)
    {
        NextTriggerIdx = TriggerIdx->at(PeakPositionOnVector+1);
    }
    if(PeakPositionOnVector>0)
    {
        PreviousTriggerIdx = TriggerIdx->at(PeakPositionOnVector-1);
    }
    
    TimeDifferenceNext = abs(NextTriggerIdx - PeakIdx);
    TimeDifferencePrevious = abs(PeakIdx - PreviousTriggerIdx);

    if(std::min(TimeDifferenceNext, TimeDifferencePrevious)>550) return true;
    else return false;

    return false;
}

void opdet::SERPulseFinder::SubstractBaseline(){

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

DEFINE_ART_CLASS_TOOL(opdet::SERPulseFinder)