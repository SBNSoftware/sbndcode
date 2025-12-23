
// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/Ptr.h"

// Tool include
#include "larreco/Calorimetry/INormalizeCharge.h"

// Services
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// ROOT includes
#include "TFile.h"
#include "TH2F.h"

// C++ includes
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>



namespace sbnd {
  namespace calo {
    
    class NormalizeYZ : public INormalizeCharge{
    public:
      explicit NormalizeYZ(fhicl::ParameterSet const& pset);
      
      void configure(const fhicl::ParameterSet& pset) override;
      
      double Normalize(double dQdx,
                       const art::Event& evt,
                       const recob::Hit& hit,
                       const geo::Point_t& location,
                       const geo::Vector_t& direction,
                       double t0) override;
      
    private:
      void reconfigure(const fhicl::ParameterSet& pset);
      
      std::map<std::string, TH2F*> fCorrHists;
      std::string fFileName;
      bool fVerbose;
    };
    
    NormalizeYZ::NormalizeYZ(fhicl::ParameterSet const& pset){
      reconfigure(pset);  // delegate config to reusable function
    }
    
    void NormalizeYZ::configure(const fhicl::ParameterSet& pset){
      reconfigure(pset); 
    }
    
    void NormalizeYZ::reconfigure(const fhicl::ParameterSet& pset){
      fFileName = pset.get<std::string>("FileName");
      fVerbose = pset.get<bool>("Verbose", false);

      std::string fname;
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(fFileName, fname);  
    
      TFile* f = TFile::Open(fname.c_str(), "READ");
      if (!f || f->IsZombie()) {
        throw cet::exception("NormalizeYZ") << "Failed to open correction map file: " << fFileName;
      }
      
      for (int plane = 0; plane < 3; plane++){   // planes : 2 inductions (idx : 0, 1) and 1 collection (idx : 2)
        for (int tpc = 0; tpc < 2; tpc++){     // tpc : east (idx : 0) and west (idx : 1) TPCs
          std::string histname = Form("CzyHist_%d_%d", plane, tpc);
          TH2F* h = (TH2F*)f->Get(histname.c_str());
          if (!h) {
            throw cet::exception("NormalizeYZ") << "Missing histogram in file: " << histname;
          }
	  
          fCorrHists[Form("plane%d_%d", plane, tpc)] = (TH2F*)h->Clone();
          fCorrHists[Form("plane%d_%d", plane, tpc)]->SetDirectory(nullptr);
        }
      }
      f->Close();
    }

    double NormalizeYZ::Normalize(double dQdx,
				  const art::Event& evt,
				  const recob::Hit& hit,
				  const geo::Point_t& location,
				  const geo::Vector_t&,
				  double){
      
      int plane = hit.WireID().Plane;
      //int tpc = hit.WireID().TPC;
   
      // just to be sure and consistent with the logic of input YZ map
      // seems like current setup of directly calling hit.WireID().TPC has a tolerance of 30cm
      // meaning - an approx region of -30<x<30 lies in both the TPCs    
      int tpc;
      if(location.X()<0){
	tpc = 0;
      } else tpc = 1;
     

      std::string key = Form("plane%d_%d", plane, tpc);
      
      auto it = fCorrHists.find(key);
      if (it == fCorrHists.end()) {
        mf::LogWarning("NormalizeYZ") << "No correction histogram for " << key << ". Returning uncorrected dQdx";
        return dQdx;
      }
      
      TH2F* hCorr = it->second;
      int binX = hCorr->GetXaxis()->FindBin(location.Z());
      int binY = hCorr->GetYaxis()->FindBin(location.Y());
      double scale = hCorr->GetBinContent(binX, binY);

      if (fVerbose){
        std::cout << "[NormalizeYZ] Plane: " << plane << ", TPC: " << tpc
                  << ", Y: " << location.Y() << ", Z: " << location.Z()
                  << ", Scale: " << scale << ", dQdx (raw): " << dQdx
                  << ", dQdx (corrected): " << dQdx / scale << std::endl;
      }
      
      if (scale < 1e-3){
        mf::LogWarning("NormalizeYZ") << "Invalid scale at (Y,Z)=(" << location.Y() << "," << location.Z() << "). Returning uncorrected dQdx";
        return dQdx;
      }
      
      return dQdx*scale;
    }
    
    DEFINE_ART_CLASS_TOOL(NormalizeYZ)
    
  } // namespace calo
} // namespace sbnd

