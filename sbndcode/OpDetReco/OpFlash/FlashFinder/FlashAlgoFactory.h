/**
 * \file FlashAlgoFactory.h
 * 
 * \ingroup FlashFinder
 *  
 * \brief Class def header for a class FlashAlgoFactory
 * 
 * @author kazuhiro
 */

/** \addtogroup FlashAlgoFinder

    @{*/
#ifndef FLASHALGOFACTORY_H
#define FLASHALGOFACTORY_H

#include <iostream>
#include <map>
#include "FlashAlgoBase.h"

namespace lightana {
  /**
      \class FlashFactoryBase
      \brief Abstract base class for factory (to be implemented per flash)
  */

  class FlashAlgoFactoryBase {

  public:
    /// Default ctor
    FlashAlgoFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~FlashAlgoFactoryBase(){}
    /// Abstract constructor method
    virtual FlashAlgoBase* create(const std::string instance_name) = 0;
  };

    /**
      \class FlashAlgoFactory
      \brief Factory class for instantiating flash algorithm instance
    */
    class FlashAlgoFactory {
    private:
    /// Default ctor, shouldn't be used
      FlashAlgoFactory() {}
    public:
    /// Default dtor
      ~FlashAlgoFactory() {_factory_map.clear();}
    /// Static sharable instance getter
      static FlashAlgoFactory& get()
      { if(!_me) _me = new FlashAlgoFactory; return *_me; }
      /// Factory registration method (should be called by global factory instance in algorithm header)
      void add_factory(const std::string name, lightana::FlashAlgoFactoryBase* factory)
      { _factory_map[name] = factory; }
      /// Factory creation method (should be called by clients, possibly you!)
      FlashAlgoBase* create(const std::string name, const std::string instance_name) {
        auto iter = _factory_map.find(name);
        if(iter == _factory_map.end() || !((*iter).second)) {
          std::cerr << "Found no registered class " << name << std::endl;
          return nullptr;
        }
        auto ptr = (*iter).second->create(instance_name);
        return ptr;
      }
    
    private:
    /// Static factory container
    std::map<std::string,lightana::FlashAlgoFactoryBase*> _factory_map;
    /// Static self
    static FlashAlgoFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group
