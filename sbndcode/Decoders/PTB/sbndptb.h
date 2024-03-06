////////////////////////////////////////////////////////////////////////
//
// An offline-friendly packaging of the Penn Trigger Board (PTB)'s data for SBND
// Tom Junk, November 10, 2023
// modeled after pdspctb.h from ProtoDUNE-SP
//
////////////////////////////////////////////////////////////////////////

#ifndef  sbndptb_H
#define  sbndptb_H

#include "RtypesCore.h"
#include <stdint.h>
#include <vector>

namespace raw {

  namespace ptb {

    struct Trigger {
      ULong64_t timestamp;
      ULong64_t trigger_word;
      uint32_t word_type;
    };

    struct ChStatus {
      ULong64_t timestamp;
      uint32_t beam;
      uint32_t crt;
      uint32_t pds;
      uint32_t mtca;
      uint32_t nim;
      uint32_t auxpds;
      uint32_t word_type;
    };

    struct Feedback {
      ULong64_t timestamp;
      uint32_t code;
      uint32_t source;
      uint64_t payload;
      uint32_t word_type;
    };

    struct Misc {
      ULong64_t timestamp;
      ULong64_t payload;
      uint32_t word_type;
    };

    struct WordIndex {
      uint32_t word_type;
      uint32_t index;
    };

    class sbndptb
    {

    public:

      sbndptb() {}; // Constructor of an emtpy data product

      // constructor with all of the vectors set

      sbndptb(std::vector<raw::ptb::Trigger> &HLtrigs,
	    std::vector<raw::ptb::Trigger> &LLtrigs,
	    std::vector<raw::ptb::ChStatus> &chstats,
	    std::vector<raw::ptb::Feedback> &fbs,
	    std::vector<raw::ptb::Misc> &m,
	    std::vector<raw::ptb::WordIndex> &wordindexes) : 
      fHLTriggers(HLtrigs),
      fLLTriggers(LLtrigs),
      fChStatuses(chstats),
      fFeedbacks(fbs),
      fMiscs(m),
      fIndexes(wordindexes) {};

      const std::vector<raw::ptb::Trigger>&     GetHLTriggers() const;   
      const std::vector<raw::ptb::Trigger>&     GetLLTriggers() const;   
      const std::vector<raw::ptb::ChStatus>&    GetChStatuses() const; 
      const std::vector<raw::ptb::Feedback>&    GetFeedbacks() const;  
      const std::vector<raw::ptb::Misc>&        GetMiscs() const;
      const std::vector<raw::ptb::WordIndex>&   GetIndexes() const;

      size_t  GetNTriggers() const;   
      size_t  GetNHLTriggers() const;   
      size_t  GetNLLTriggers() const;   
      size_t  GetNChStatuses() const; 
      size_t  GetNFeedbacks() const;  
      size_t  GetNMiscs() const;      
      size_t  GetNIndexes() const;      

      const raw::ptb::Trigger&    GetHLTrigger(size_t i) const;   
      const raw::ptb::Trigger&    GetLLTrigger(size_t i) const;   
      const raw::ptb::ChStatus&   GetChStatuse(size_t i) const; 
      const raw::ptb::Feedback&   GetFeedback(size_t i) const;  
      const raw::ptb::Misc&       GetMisc(size_t i) const;      
      const raw::ptb::WordIndex&  GetIndex(size_t i) const;      

    private:

      std::vector<raw::ptb::Trigger> fHLTriggers;
      std::vector<raw::ptb::Trigger> fLLTriggers;
      std::vector<raw::ptb::ChStatus> fChStatuses;
      std::vector<raw::ptb::Feedback> fFeedbacks;
      std::vector<raw::ptb::Misc> fMiscs;
      std::vector<raw::ptb::WordIndex> fIndexes;
    };
  } // namespace ptb
} // namespace raw

// accessors

const std::vector<raw::ptb::Trigger>&       raw::ptb::sbndptb::GetHLTriggers()   const { return fHLTriggers; }
const std::vector<raw::ptb::Trigger>&       raw::ptb::sbndptb::GetLLTriggers()   const { return fLLTriggers; }
const std::vector<raw::ptb::ChStatus>&      raw::ptb::sbndptb::GetChStatuses() const { return fChStatuses; }
const std::vector<raw::ptb::Feedback>&      raw::ptb::sbndptb::GetFeedbacks()  const { return fFeedbacks; }
const std::vector<raw::ptb::Misc>&          raw::ptb::sbndptb::GetMiscs()      const { return fMiscs; }
const std::vector<raw::ptb::WordIndex>&     raw::ptb::sbndptb::GetIndexes()    const { return fIndexes; }

size_t  raw::ptb::sbndptb::GetNTriggers()   const { return fHLTriggers.size() + fLLTriggers.size(); }
size_t  raw::ptb::sbndptb::GetNHLTriggers()   const { return fHLTriggers.size(); }
size_t  raw::ptb::sbndptb::GetNLLTriggers()   const { return fLLTriggers.size(); }
size_t  raw::ptb::sbndptb::GetNChStatuses() const { return fChStatuses.size(); }
size_t  raw::ptb::sbndptb::GetNFeedbacks()  const { return fFeedbacks.size(); }
size_t  raw::ptb::sbndptb::GetNMiscs()      const { return fMiscs.size(); }
size_t  raw::ptb::sbndptb::GetNIndexes()    const { return fIndexes.size(); }

const raw::ptb::Trigger&     raw::ptb::sbndptb::GetHLTrigger(size_t i)   const { return fHLTriggers.at(i); }
const raw::ptb::Trigger&     raw::ptb::sbndptb::GetLLTrigger(size_t i)   const { return fLLTriggers.at(i); }
const raw::ptb::ChStatus&    raw::ptb::sbndptb::GetChStatuse(size_t i) const { return fChStatuses.at(i); }
const raw::ptb::Feedback&    raw::ptb::sbndptb::GetFeedback(size_t i)  const { return fFeedbacks.at(i); }
const raw::ptb::Misc&        raw::ptb::sbndptb::GetMisc(size_t i)      const { return fMiscs.at(i); }
const raw::ptb::WordIndex&   raw::ptb::sbndptb::GetIndex(size_t i)     const { return fIndexes.at(i); }

#endif // sbndptb_H
