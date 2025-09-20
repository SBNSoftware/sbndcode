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

  //Returns rescaled number of PE
  virtual double NObservedPE(size_t bin, std::vector<unsigned int> & pe_vector) = 0;
};

#endif
