///////////////////////////////////////////////////////////////////////
///
/// \file   IROIFinder.h
///
/// \brief  This provides an interface for tools which are tasked with
///         finding ROI's in input waveforms. This allows different
///         approaches to be tried interchangeably
///
/// \author T. Usher
///
/// Ported from ICARUS to SBND by A. Scarff
////////////////////////////////////////////////////////////////////////
#ifndef IROIFinder_H
#define IROIFinder_H
#include "fhiclcpp/ParameterSet.h"
namespace art
{
  class TFileDirectory;
}
namespace sbnd_tool
{
    class IROIFinder
    {
    public:
      virtual ~IROIFinder() noexcept = default;
        
      virtual void   configure(const fhicl::ParameterSet& pset)                = 0;
      virtual void   initializeHistograms(art::TFileDirectory&)          const = 0;
      virtual size_t plane()                                             const = 0;
        
      // Define the waveform container
      using Waveform        = std::vector<float>;
        
      // Define the ROI and its container
      using CandidateROI    = std::pair<size_t, size_t>;
      using CandidateROIVec = std::vector<CandidateROI>;
        
      // Find the ROI's
      virtual void FindROIs(const Waveform&, size_t, CandidateROIVec&) const = 0;
    };
}
#endif
