#ifndef NCPIZEROSTRUCTS_H_SEEN
#define NCPIZEROSTRUCTS_H_SEEN

constexpr double kPiZeroMass = 134.9769;

std::vector<double> kKFCovMatrix = {11964.3, 0, 0, 0, 3179.01, 0, 0, 0, 0.315488};

constexpr int n_fluxweight_univs   = 1000;
constexpr int n_genieweight_univs  = 500;
constexpr int n_geant4weight_univs = 1000;

const std::vector<std::string> flux_weight_names = { "expskin_Flux",
                                                     "horncurrent_Flux",
                                                     "kminus_Flux",
                                                     "kplus_Flux",
                                                     "kzero_Flux",
                                                     "nucleoninexsec_Flux",
                                                     "nucleonqexsec_Flux",
                                                     "nucleontotxsec_Flux",
                                                     "piminus_Flux",
                                                     "pioninexsec_Flux",
                                                     "pionqexsec_Flux",
                                                     "piontotxsec_Flux",
                                                     "piplus_Flux"
};

const std::vector<std::string> genie_weight_names = { "GENIEReWeight_SBND_v5_multisigma_CoulombCCQE",
                                                      "GENIEReWeight_SBND_v5_multisigma_DecayAngMEC",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarnCC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarnCC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarnNC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarnNC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarpCC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarpCC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarpNC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvbarpNC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvnCC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvnCC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvnNC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvnNC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvpCC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvpCC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvpNC1pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NonRESBGvpNC2pi",
                                                      "GENIEReWeight_SBND_v5_multisigma_NormCCMEC",
                                                      "GENIEReWeight_SBND_v5_multisigma_NormNCMEC",
                                                      "GENIEReWeight_SBND_v5_multisigma_NormCCCOH",
                                                      "GENIEReWeight_SBND_v5_multisigma_NormNCCOH",
                                                      "GENIEReWeight_SBND_v5_multisigma_RDecBR1eta",
                                                      "GENIEReWeight_SBND_v5_multisigma_RDecBR1gamma",
                                                      "GENIEReWeight_SBND_v5_multisigma_RPA_CCQE",
                                                      "GENIEReWeight_SBND_v5_multisigma_ThetaDelta2NRad",
                                                      "GENIEReWeight_SBND_v5_multisigma_Theta_Delta2Npi",
                                                      "GENIEReWeight_SBND_v5_multisigma_VecFFCCQEshape",
                                                      "GENIEReWeight_SBND_v5_multisim_CCRESVariationResponse",
                                                      "GENIEReWeight_SBND_v5_multisim_DISBYVariationResponse",
                                                      "GENIEReWeight_SBND_v5_multisim_FSI_N_VariationResponse",
                                                      "GENIEReWeight_SBND_v5_multisim_FSI_pi_VariationResponse",
                                                      "GENIEReWeight_SBND_v5_multisim_NCELVariationResponse",
                                                      "GENIEReWeight_SBND_v5_multisim_NCRESVariationResponse",
                                                      "GENIEReWeight_SBND_v5_multisim_ZExpAVariationResponse"
};

const std::vector<std::string> geant4_weight_names = { "reinteractions_Geant4",
};

namespace NC {
  enum EventType
  {
    kSignalNCPiZero,
    kOtherNCPiZero,
    kOtherNC,
    kCCNuMu,
    kCCNuE,
    kDirt,
    kNonFV,
    kCosmic,
    kFailedTruthMatch,
    kUnknownEv = -1
  };
}

namespace CC {
  enum EventType
  {
    kSignalCCPiZero,
    kOtherCCPiZero,
    kNC,
    kOtherCCNuMu,
    kCCNuE,
    kDirt,
    kNonFV,
    kCosmic,
    kFailedTruthMatch,
    kUnknownEv = -1
  };
}

enum VarType
  {
    kBool,
    kInt,
    kUInt,
    kFloat,
    kDouble,
    kUnknownVar = -1
  };

enum VecType
  {
    kOneD,
    kTwoD,
    kUnknownVec = -1
  };
    
class VecVar
{
  std::string name;

  public:
  
  VecVar(std::string n = "")
    : name(n)
  {}

  std::string Name() const
  {
    return name;
  }

  virtual void Clear() {}

  virtual VarType IdentifyVar() const
  {
    return kUnknownVar;
  }

  virtual VecType IdentifyVec() const
  {
    return kUnknownVec;
  }

  virtual void Resize(const int size) = 0;

  template<typename T>
  void Assign(const int size, const T val) {}

  template<typename T>
  void SetVal(const int pos, const T val) {}

  template<typename T>
  void SetVal(const int pos, const std::vector<T> &val) {}

  virtual void Resize(const int pos, const int size) = 0;

  template<typename T>
  void Assign(const int pos, const int size, const T val) {}

  template<typename T>
  void SetVal(const int posA, const int posB, const T val) {}
};

template<typename T>
class InhVecVar;

template<typename T>
class InhVecVar : public VecVar
{
  std::vector<T> var;

  public:

  InhVecVar(std::string n, std::vector<T> v)
    : VecVar(n)
    , var(v)
    {}

  InhVecVar(std::string n = "")
    : VecVar(n)
    , var(std::vector<T>())
    {}

  std::vector<T>& Var()
  {
    return var;
  }

  void Clear()
  {
    var.clear();
  }

  virtual VarType IdentifyVar() const;

  VecType IdentifyVec() const
  {
    return kOneD;
  }

  void Resize(const int size)
  {
    var.resize(size);
  }

  void Assign(const int size, const T value)
  {
    var.assign(size, value);
  }

  void SetVal(const int pos, const T value)
  {
    var[pos] = value;
  }

  T GetVal(const int pos)
  {
    return var[pos];
  }

  void Resize(const int pos, const int size) {};

  void Assign(const int pos, const int size, const T val) {}

  void SetVal(const int pos, const std::vector<T> value) {}

  void SetVal(const int posA, const int posB, const T val) {}
};

template<>
VarType InhVecVar<bool>::IdentifyVar() const { return kBool; }

template<>
VarType InhVecVar<int>::IdentifyVar() const { return kInt; }

template<>
VarType InhVecVar<size_t>::IdentifyVar() const { return kUInt; }

template<>
VarType InhVecVar<float>::IdentifyVar() const { return kFloat; }

template<>
VarType InhVecVar<double>::IdentifyVar() const { return kDouble; }

template<typename T>
class InhVecVecVar;

template<typename T>
class InhVecVecVar : public VecVar
{
  std::vector<std::vector<T>> var;

  public:

 InhVecVecVar(std::string n, std::vector<std::vector<T>> v)
    : VecVar(n)
    , var(v)
    {}

  InhVecVecVar(std::string n = "")
    : VecVar(n)
    , var(std::vector<std::vector<T>>())
    {}

  std::vector<std::vector<T>>& Var()
  {
    return var;
  }

  void Clear()
  {
    var.clear();
  }

  virtual VarType IdentifyVar() const;

  VecType IdentifyVec() const
  {
    return kTwoD;
  }

  void Resize(const int size)
  {
    var.resize(size);
  }

  void Resize(const int pos, const int size)
  {
    var[pos].resize(size);
  }

  void Assign(const int pos, const int size, const T value)
  {
    var[pos].assign(size, value);
  }

  void SetVal(const int pos, const std::vector<T> value)
  {
    var[pos] = value;
  }

  void SetVal(const int posA, const int posB, const T value)
  {
    var[posA][posB] = value;
  }

  std::vector<T> GetVecVal(const int pos)
  {
    return var[pos];
  }

  T GetVal(const int posA, const int posB)
  {
    return var[posA][posB];
  }
};

template<>
VarType InhVecVecVar<bool>::IdentifyVar() const { return kBool; }

template<>
VarType InhVecVecVar<int>::IdentifyVar() const { return kInt; }

template<>
VarType InhVecVecVar<size_t>::IdentifyVar() const { return kUInt; }

template<>
VarType InhVecVecVar<float>::IdentifyVar() const { return kFloat; }

template<>
VarType InhVecVecVar<double>::IdentifyVar() const { return kDouble; }

typedef std::map<std::string, VecVar*> VecVarMap;

template <typename T,
  typename TIter = decltype(std::begin(std::declval<T>())),
  typename = decltype(std::end(std::declval<T>()))>
  constexpr auto enumerate(T && iterable)
{
  struct iterator
  {
    size_t i;
    TIter iter;
    bool operator != (const iterator & other) const { return iter != other.iter; }
    void operator ++ () { ++i; ++iter; }
    auto operator * () const { return std::tie(i, *iter); }
  };
  struct iterable_wrapper
  {
    T iterable;
    auto begin() { return iterator{ 0, std::begin(iterable) }; }
    auto end() { return iterator{ 0, std::end(iterable) }; }
  };
  return iterable_wrapper{ std::forward<T>(iterable) };
}

#endif
