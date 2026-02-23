// Simple class used to average the waveforms of the same channel.
// Also stores the number of waveforms added to the average.


#ifndef CALLOS_AV_WVF_H
#define CALLOS_AV_WVF_H

#include <vector>
#include "lardataobj/RawData/OpDetWaveform.h"


class AverageWaveform {
  public:
      AverageWaveform(int size) : size_(size), wvf_count_(1) 
      {
        array_.resize(size, 0);
      }

      void addToAverage(const std::vector<float>& inputArray) {
        for (int i = 0; i < size_; i++) {
          // divide by wvf_count_ to avoid overflow, maybe for large number of waveforms is better to stack a bunch and divide after
          array_[i] += inputArray[i] / wvf_count_;
        }
        wvf_count_++;
      }

      int getWvfCount() const { return wvf_count_; }

  private:
    std::vector<float> array_;
    int size_;
    int wvf_count_;
};

//ROI finder class placeholder, we need floats as we'll be subtracting baselines
class SimpleROI :  public std::vector<float> 
{
  public:
    SimpleROI(int size, int channel) : fROISamples(size), fChannel(channel) 
    {
      this->resize(size, 0);
    }
    void SetCharge(double charge) {fCharge = charge;}
    double Charge() {return fCharge;}
    void SetPedestalCharge(double pedestalCharge) {fPedestalCharge = pedestalCharge;}
    double PedestalCharge() {return fPedestalCharge;}
    void SetStartTick(unsigned int startTick) {fStartTick = startTick;}
    int StartTick() {return fStartTick;}
    void SetBaseline(double baseline) {fBaseline = baseline;}
    double Baseline() {return fBaseline;}
    void SetBaselineSTD(double baselineSTD) {fBaselineSTD = baselineSTD;}
    double BaselineSTD() {return fBaselineSTD;}
    void SetEndTick(unsigned int endTick) {fEndTick = endTick;}
    int EndTick() {return fEndTick;}
    // void SetWaveform(const std::vector<float>& inputArray) {this->assign(inputArray.begin(), inputArray.end());}
    // std::vector<float>& Waveform() { return *this; }
    void SetPeak(double peak) {fPeak = peak;}
    double Peak() {return fPeak;}
    void SetPeakTick(unsigned int peakTick) {fPeakTick = peakTick;}
    int PeakTick() {return fPeakTick;}
    // void SetFFT_ROI(const std::vector<float>& inputArray) {this->assign(inputArray.begin(), inputArray.end());}
    // std::vector<float>& FFT_ROI() { return *this; }

    void SetWaveform(const std::vector<float>& waveform) {
        this->waveform = waveform;
    }

    void SetFFT_ROI(const std::vector<float>& fft_roi) {
        this->fft_roi = fft_roi;
    }

    const std::vector<float>& Waveform() const {
        return waveform;
    }

    const std::vector<float>& FFT_ROI() const {
        return fft_roi;
    }
    // void SetDerivateROI(const std::vector<float>& inputArray) {this->assign(inputArray.begin(), inputArray.end());}
    // std::vector<float>& DerivateROI() { return *this; }
    // std::vector<float>& AverageWaveform_2() { return *this; }


    // // Functions included for backwards compatability with previous data types
    // std::vector<float>& Waveform() { return *this; }

    // Functions included for backwards compatability with previous data types
    // std::vector<float> const& Waveform() const { return *this; }

    // Destructor
    ~SimpleROI() {}

private:
  int fROISamples;
  unsigned int fChannel;
  double fCharge;
  unsigned int fStartTick;
  double fBaseline;
  double fBaselineSTD;
  unsigned int fEndTick;
  double fPeak;
  unsigned int fPeakTick;
  std::vector<float> waveform;
  std::vector<float> fft_roi;
  double fPedestalCharge;



};

#endif