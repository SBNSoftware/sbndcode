// Simple class used to average the waveforms of the same channel.
// Also stores the number of waveforms added to the average.

class AverageWaveform {
public:
  AverageWaveform(int size) : size_(size){
    array_ = new float[size_];
    wvf_count_ = 1;
  }

  ~AverageWaveform() {
    delete[] array_;
  }

  void addToAverage(float* inputArray) {
    for (int i = 0; i < size_; i++) {
      // divide by wvf_count_ to avoid overflow, maybe for large number of waveforms is better to stack a bunch and divide after
      array_[i] += inputArray[i]/wvf_count_;
    }
    wvf_count_++;
  }

private:
  float* array_;
  int size_;
  int wvf_count_;
};
