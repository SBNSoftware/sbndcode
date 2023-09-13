enum EventType
  {
    kNCPiZero,
    kOtherNC,
    kCCNuMu,
    kCCNuE,
    kDirt,
    kNonFV,
    kCosmic,
    kBadRecoSignal,
    kUnknownEv = -1
  };

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

  void SetVal(const int posA, const int posB, const T value)
  {
    var[posA][posB] = value;
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
