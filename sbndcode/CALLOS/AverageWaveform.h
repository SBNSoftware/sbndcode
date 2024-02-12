// Simple class used to average the waveforms of the same channel.
// Also stores the number of waveforms added to the average.


#include <vector>

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
