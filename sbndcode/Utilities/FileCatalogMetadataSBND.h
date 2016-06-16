////////////////////////////////////////////////////////////////////////
//
// Name:  FileCatalogMetadataSBND.h.  
//
// Purpose:  Art service adds microboone-specific per-job sam metadata.
//
//           FCL parameters:
//
//           FCLName        - FCL file name.
//           FCLVersion     - FCL file version.
//           ProjectName    - Project name.
//           ProjectStage   - Project stage.
//           ProjectVersion - Project version.
//
//           Above values will be added in internal metadata of artroot
//           output files whenever this service is included in job
//           configuration (service does not need to be called).  The
//           public interface of this service consists of accessors that
//           allow other code to discover above metadata parameters.
//
// Created:  28-Oct-2014,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#ifndef FILECATALOGMETADATASBND_H
#define FILECATALOGMETADATASBND_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace util {

  class FileCatalogMetadataSBND
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataSBND(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataSBND() = default;

    // Accessors.

    const std::string& FCLName() const {return fFCLName;}
    const std::string& FCLVersion() const {return fFCLVersion;}
    const std::string& ProjectName() const {return fProjectName;}
    const std::string& ProjectStage() const {return fProjectStage;}
    const std::string& ProjectVersion() const {return fProjectVersion;}

  private:

    // Callbacks.

    void postBeginJob();

    // Data members.

    std::string fFCLName;
    std::string fFCLVersion;
    std::string fProjectName;
    std::string fProjectStage;
    std::string fProjectVersion;
  };

} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataSBND, LEGACY)

#endif
