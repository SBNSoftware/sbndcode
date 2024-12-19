#include "sbndcode/SBNDCVN/module_helpers/SBNDICVNZlibMaker.h"
#include "lardataobj/RecoBase/Slice.h"

#include "canvas/Persistency/Common/FindManyP.h"

namespace fs = boost::filesystem;

namespace lcvn {

  SBNDICVNZlibMaker::SBNDICVNZlibMaker(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  SBNDICVNZlibMaker::~SBNDICVNZlibMaker() {}

  //......................................................................
  void SBNDICVNZlibMaker::reconfigure(const fhicl::ParameterSet& pset)
  {
    fOutputDir = pset.get<std::string>("OutputDir", "");
    fPixelMapInput = pset.get<std::string>("PixelMapInput");
    fSetLog = pset.get<bool>("SetLog");
    fReverseViews = pset.get<std::vector<bool>>("ReverseViews");

    fPlaneLimit = pset.get<unsigned int>("PlaneLimit");
    fTDCLimit = pset.get<unsigned int>("TDCLimit");
    fverbose = pset.get<bool>("verbose");
    fUseSlice = pset.get<bool>("UseSlice");
    fSliceLabel = pset.get<std::string>("SliceLabel");
  }

  void SBNDICVNZlibMaker::beginJob()
  {
    fImage = CVNImageUtils(fPlaneLimit, fTDCLimit, 3);
    fImage.SetLogScale(fSetLog);
    fImage.SetViewReversal(fReverseViews);

    if (fOutputDir != "")
      out_dir = fOutputDir;

    else
      out_dir = ".";

    // Throw an error if the specified output directory doesn't exist
    if (!fs::exists(out_dir))
      throw art::Exception(art::errors::FileOpenError)
        << "Output directory " << out_dir << " does not exist!" << std::endl;
  }

}
