///////////////////////////////////////////////////////////////////////
///
/// Interface class for PMTNonLinearity tool
///
////////////////////////////////////////////////////////////////////////

#ifndef SBND_PMTNonLinearity_H
#define SBND_PMTNonLinearity_H



namespace opdet {
  class PMTNonLinearity;
}

//Base class
class opdet::PMTNonLinearity {
public:
  //Constructor
  virtual ~PMTNonLinearity() noexcept = default;
  
  //Configure non-linearity tool
  virtual void ConfigureNonLinearity() = 0;

  //Returns rescaled number of PE
  virtual double NObservedPE(size_t bin, std::vector<unsigned int> & pe_vector){
    return pe_vector[bin];
  }
  virtual double NObservedPE(int opch, size_t bin, std::vector<unsigned int> & pe_vector){
    // Default implementation calls the non-channel-dependent version 
    return NObservedPE(bin, pe_vector);
  }

};

#endif
