/// \file    SBNDPTBRawUtils.cxx
/// \brief   SBND PTB raw data utilities
/// \author  trj@fnal.gov

#include "SBNDPTBRawUtils.h"
#include "sbndaq-artdaq-core/Overlays/SBND/PTB_content.h"

namespace raw {
  namespace ptb {
    const std::vector<raw::ptb::ChStatus>  GetChStatusBeforeHLTs(const raw::ptb::sbndptb &pdata)
    {
      std::vector<raw::ptb::ChStatus> chs;
      raw::ptb::ChStatus emptychstat;
      emptychstat.timestamp = 0;
      emptychstat.beam = 0;
      emptychstat.crt = 0;
      emptychstat.pds = 0;
      emptychstat.mtca = 0;
      emptychstat.nim = 0;
      emptychstat.auxpds = 0;
      emptychstat.word_type = 0;

      // Find the last CHStatus before each HLT

      const auto &hlts = pdata.GetHLTriggers();
      const auto &idxs = pdata.GetIndexes();
      const auto &chst = pdata.GetChStatuses();
      
      for (size_t i=0; i<hlts.size(); ++i)
	{
	  for (size_t j=0; j<idxs.size(); ++j)
	    {
	      if (idxs.at(j).word_type == (uint32_t) ::ptb::content::word::t_gt && idxs.at(j).index == i)
		{
		  size_t kstatindex = j;
		  if (kstatindex > 0)
		    {
		      kstatindex --;  // it's the word before the HLT that has the chstat
		      if (idxs.at(kstatindex).word_type == (uint32_t) ::ptb::content::word::t_ch)
			{
			  size_t kstat = idxs.at(kstatindex).index;
		      
			  if (kstat < chst.size())
			    {
			      chs.push_back(chst.at(kstat));
			    }
			  else
			    {
			      chs.push_back(emptychstat);
			    }
			}
		      else
			{
			  chs.push_back(emptychstat);
			}
		    }
		  else
		    {
		      chs.push_back(emptychstat);
		    }
		}
	    }
	}  
      return chs;
    }
  }
}
