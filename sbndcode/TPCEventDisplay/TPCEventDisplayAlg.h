#ifndef TPCEVENTDISPLAYALG_H_SEEN
#define TPCEVENTDISPLAYALG_H_SEEN

///////////////////////////////////////////////
// TPCEventDisplayAlg.h
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

//larsoft
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// ROOT
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TView3D.h"
#include "TAxis3D.h"

namespace detinfo { class DetectorClocksData; }

namespace sbnd::crt {
  
  class TPCEventDisplayAlg {
  public:
    
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> SliceLabel {
        Name("SliceLabel")
          };

      fhicl::Atom<bool> DrawTPC {
        Name("DrawTPC")
          };
      fhicl::Atom<bool> DrawSlices {
        Name("DrawSlices")
          };

      fhicl::Atom<int> TPCColour {
        Name("TPCColour")
          };
      fhicl::Atom<int> SliceStartingColour {
        Name("SliceStartingColour")
          };
      fhicl::Atom<int> SliceColourInterval {
        Name("SliceColourInterval")
          };

      fhicl::Atom<double> LineWidth {
        Name("LineWidth")
          };
    };
    
    TPCEventDisplayAlg(const Config& config);
    
    TPCEventDisplayAlg(const fhicl::ParameterSet& pset) :
      TPCEventDisplayAlg(fhicl::Table<Config>(pset, {})()) {}
    
    TPCEventDisplayAlg();

    ~TPCEventDisplayAlg();

    void reconfigure(const Config& config);

    void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour, int lineWidth = -1);

    void Draw(detinfo::DetectorClocksData const& clockData, const art::Event& event,
              const TString& saveName);

  private:

    TPCGeoAlg fTPCGeoAlg;

    art::InputTag fSliceLabel;

    bool fDrawTPC;
    bool fDrawSlices;

    int fTPCColour;
    int fSliceStartingColour;
    int fSliceColourInterval;

    double fLineWidth;
  };
}

#endif
