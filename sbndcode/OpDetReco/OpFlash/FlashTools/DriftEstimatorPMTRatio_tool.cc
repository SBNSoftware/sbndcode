///////////////////////////////////////////////////////////////////////
/// File: DriftEstimatorPMTRatio_tool.cc
///
/// Base class: DriftEstimatorBase.hh
///
/// Tool description: this tool estimates the drift coordinate
/// from the ratio between the #PE reconstructed for the
/// uncoated/coated PMTs. It requires a calibration curve
/// (speficied in the CalibrationFile fhicl parameter).
/// Once the drift has been estimated, the photon propagation
/// time is calculated using the VUV and VIS light group velocities
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/types/Atom.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

#include <vector>
#include <string>
#include <map>

#include "DriftEstimatorBase.hh"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

//ROOT includes
#include "TProfile.h"
#include "TFile.h"

namespace lightana
{
  class DriftEstimatorPMTRatio : DriftEstimatorBase{

  public:

    //Configuration parameters
    struct Config {

      fhicl::Atom<std::string> CalibrationFile {
        fhicl::Name("CalibrationFile"),
        fhicl::Comment("Filepath to the calibration ROOT file")
      };

      fhicl::Atom<double> VGroupVUV {
        fhicl::Name("VGroupVUV"),
        fhicl::Comment("Group velocity for VUV photons")
      };

      fhicl::Atom<double> VGroupVIS {
        fhicl::Name("VGroupVIS"),
        fhicl::Comment("Group velocity for VIS photons")
      };

    };

    // Default constructor
    explicit DriftEstimatorPMTRatio(art::ToolConfigTable<Config> const& config);

    // Method giving the estimated drift coordinate
    double GetDriftPosition(std::vector<double> PE_v) override;

    // Method giving the photon propagation
    double GetPropagationTime(double drift) override;

    // Method giving the photon propagation from PE vector
    double PEToPropagationTime(std::vector<double> PE_v);

  private:
    double Interpolate(double val);

    // Input filepah with calibration curve
    std::string fCalibrationFile;

    // Scintillation light group velocities
    double fVGroupVUV;
    double fVGroupVIS;
    double fVGroupVUV_I;

    //Geo properties
    double fDriftDistance;
    double fVISLightPropTime;
    double fKinkDistance;

    // Vectors with calibration variables
    std::vector<double> fPMTRatioCal;
    std::vector<double> fDriftCal;
    int fNCalBins;
    double fPMTRatio_MinVal;
    double fPMTRatio_MaxVal;

    // PDS mapping
    opdet::sbndPDMapAlg fPDSMap;

    std::set<int> fPDSBoxIDs;

  };

  DriftEstimatorPMTRatio::DriftEstimatorPMTRatio(art::ToolConfigTable<Config> const& config)
    : fCalibrationFile { config().CalibrationFile() },
    fVGroupVUV { config().VGroupVUV() },
    fVGroupVIS { config().VGroupVIS() }
  {

    // Open input file
    std::string file_name;
    cet::search_path sp("FW_SEARCH_PATH");
    if ( !sp.find_file(fCalibrationFile, file_name) )
      throw cet::exception("DriftEstimatorPMTRatio") << "Calibration file " <<
          fCalibrationFile << " not found in FW_SEARCH_PATH\n";

    TFile* input_file = TFile::Open(file_name.c_str(), "READ");
    TProfile * hProf_Calibration = (TProfile*)input_file->Get("PMTRatioCalibrationProfile");

    //Fill calibration variables
    fNCalBins = hProf_Calibration->GetNbinsX();
    for (int ix=1; ix<=fNCalBins; ix++){
      fPMTRatioCal.push_back( hProf_Calibration->GetBinCenter(ix) );
      fDriftCal.push_back( hProf_Calibration->GetBinContent(ix) );
    }
    fPMTRatio_MinVal = hProf_Calibration->GetBinCenter(1);
    fPMTRatio_MaxVal = hProf_Calibration->GetBinCenter(fNCalBins);

    input_file->Close();

    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

    fDriftDistance = geom.TPC().DriftDistance();
    fVISLightPropTime = fDriftDistance/fVGroupVIS;
    fKinkDistance = 0.5*fDriftDistance*(1-fVGroupVUV/fVGroupVIS);

    fVGroupVUV_I = 1./fVGroupVUV;

    for(size_t oc=0; oc<fPDSMap.size(); oc++){
      fPDSBoxIDs.insert( fPDSMap.pdBox(oc) );
    }
  }

  double DriftEstimatorPMTRatio::GetDriftPosition(std::vector<double> PE_v){

    std::map<int, double> BoxMap_PECoated;
    std::map<int, double> BoxMap_PEUncoated;
    std::map<int, int> BoxMap_NCoatedCh;
    std::map<int, int> BoxMap_NUncoatedCh;

    for(size_t oc=0; oc<PE_v.size(); oc++){
      // skip 0 PE channels
      if(PE_v[oc]==0) continue;

      std::string pd_type=fPDSMap.pdType(oc);
      // exclude xarapucas by now
      if(pd_type=="xarapuca_vuv" || pd_type=="xarapuca_vis") continue;

      // get PDS box
      int box_id = fPDSMap.pdBox(oc);

      // we store the pe in each box per PMT flavour
      // and the number of "triggered" PMTs
      if(pd_type=="pmt_coated") {
        BoxMap_PECoated[box_id]+=PE_v[oc];
        BoxMap_NCoatedCh[box_id]+=1;
      }
      else if(pd_type=="pmt_uncoated") {
        BoxMap_PEUncoated[box_id]+=PE_v[oc];
        BoxMap_NUncoatedCh[box_id]+=1;
      }
    }

    // compute PMTRatio metric
    double PECoated=0, PEUncoated=0;
    for(size_t boxID=0; boxID<fPDSBoxIDs.size(); boxID++){
      //we need the uncoated PMT in each window and at least one coated
      if( BoxMap_NUncoatedCh[boxID]==1 && BoxMap_NCoatedCh[boxID]>=1){
        double CoWeight = 1./BoxMap_NCoatedCh[boxID];
        PECoated+=CoWeight * BoxMap_PECoated[boxID];
        PEUncoated+=BoxMap_PEUncoated[boxID];
      }
    }

    if(PECoated!=0){
      double pmtratio = PEUncoated/PECoated;

      double drift_distance;
      if(pmtratio<=fPMTRatioCal[0])
        drift_distance=fDriftCal[0];
      else if(pmtratio>=fPMTRatioCal[fNCalBins-1])
        drift_distance=fDriftCal[fNCalBins-1];
      else
        drift_distance=Interpolate(pmtratio);

      return drift_distance;
    }
    else return fDriftCal[fNCalBins-1]; 
  }

  double DriftEstimatorPMTRatio::GetPropagationTime(double drift){

    // drift is here the X coordinate
    // cathode: x=0 cm, PDS: x=200 cm
    if(std::abs(drift) > fKinkDistance)
      return (fDriftDistance-std::abs(drift)) * fVGroupVUV_I ;
    else
      return std::abs(drift) * fVGroupVUV_I + fVISLightPropTime;
  }

  double DriftEstimatorPMTRatio::PEToPropagationTime(std::vector<double> PE_v){

    double _drift = GetDriftPosition(PE_v);

    return GetPropagationTime(_drift);
  }

  double DriftEstimatorPMTRatio::Interpolate(double val){

    size_t upix = std::upper_bound(fPMTRatioCal.begin(), fPMTRatioCal.end(), val)-fPMTRatioCal.begin();

    double slope = ( fDriftCal[upix]-fDriftCal[upix] ) / ( fPMTRatioCal[upix]-fPMTRatioCal[upix-1] );

    return fDriftCal[upix-1] + slope * ( val - fPMTRatioCal[upix-1] );
  }

}

DEFINE_ART_CLASS_TOOL(lightana::DriftEstimatorPMTRatio)
