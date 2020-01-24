//Holder class for the pointer to the branch. 
#ifndef BranchType_HH
#define BranchType_HH

namespace optimiser{
  class BranchTypeBase;
  template <class T> class BranchType;
}

class optimiser::BranchTypeBase{
public: 
  virtual ~BranchTypeBase() = default;
};

template <class T>
class optimiser::BranchType : public optimiser::BranchTypeBase{
  
public:
  
  BranchType(){
    Branch = 0;
  }
  
  T*& GetPointer(){return Branch;}
  
private:
  
  T* Branch;
  
};
  
#endif
