///////////////////////////////////////////////////////////////////////
/// File: FlashGeoThreshold_tool.cc
///
/// Base class: FlashGeoBase.hh
///
/// It computes the PMTs Threshold
/// weighted by the reconstructed number of PE
///
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Utilities/ToolConfigTable.h"

#include "FlashGeoBase.hh"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbndcode/Calibration/PDSDatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/Calibration/PDSDatabaseInterface/IPMTCalibrationDatabaseService.h"

#include <string>
#include <map>

namespace lightana{

  class FlashGeoThreshold : FlashGeoBase
  {

  public:

    //Configuration parameters
     struct Config {

       fhicl::Atom<float> ThresholdY {
         fhicl::Name("ThresholdY"),
         fhicl::Comment("Use channels above this threshold (Y coord)")
       };


       fhicl::Atom<float> ThresholdZ {
         fhicl::Name("ThresholdZ"),
         fhicl::Comment("Use channels above this threshold (Z coord)")
       };

       fhicl::Sequence<std::string> PDTypes {
         fhicl::Name("PDTypes"),
         fhicl::Comment("Only use PD types specified in this list")
       };

       fhicl::Atom<bool> NormalizeByPDType {
         fhicl::Name("NormalizeByPDType"),
         fhicl::Comment("Normalize number of PE in each ZY position by PD type")
       };

       fhicl::Atom<unsigned int> WeightExp {
         fhicl::Name("WeightExp"),
         fhicl::Comment("Weight exponent for YZ barycenters")
       };  

     };


    // Default constructor
    explicit FlashGeoThreshold(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash geometry
    void GetFlashLocation(std::vector<double> pePerOpChannel,
                            double& Ycenter, double& Zcenter,
                            double& Ywidth, double& Zwidth) override;
    void InitializeFlashGeoAlgo() override;


  private:

    void ResetVars();
    void GetCenter(std::map<int, double> PEAcc, double& center, double& width, double& threshold);
    // Fhicl configuration parameters
    std::vector<std::string> fPDTypes;
    double fThresholdY;
    double fThresholdZ;

    bool fNormalizeByPDType;
    unsigned int fWeightExp;

    // Store #PE at each YZ position
    std::map<int, double> fYPEAccumulator;
    std::map<int, double> fZPEAccumulator;

    // Store PD type at each YZ position
    std::map<int, std::string> fPDTypeByY;
    std::map<int, std::string> fPDTypeByZ;

    // Store normalization factor for each PD type
    std::map<std::string, double> fNormFactorsY;
    std::map<std::string, double> fNormFactorsZ;

    // Store the number of ON PMTs at each YZ position
    std::map<int, unsigned int> fNumONPMTsY_tpc0;
    std::map<int, unsigned int> fNumONPMTsZ_tpc0;
    std::map<int, unsigned int> fNumONPMTsY_tpc1;
    std::map<int, unsigned int> fNumONPMTsZ_tpc1;

    // PDS mapping
    opdet::sbndPDMapAlg fPDSMap;
    
    // Array to store PD xyz position
    double fPDxyz[3];

    static constexpr double fDefCenterValue = -999.;

    // Calibration database service
    sbndDB::PMTCalibrationDatabase const* fPMTCalibrationDatabaseService;

  };


  FlashGeoThreshold::FlashGeoThreshold(art::ToolConfigTable<Config> const& config)
    : fPDTypes{ config().PDTypes() }
    , fThresholdY{ config().ThresholdY() }
    , fThresholdZ{ config().ThresholdZ() }
    , fNormalizeByPDType{ config().NormalizeByPDType() }
    , fWeightExp{ config().WeightExp() }
  {
    //Load PMT Calibration Database
    fPMTCalibrationDatabaseService = lar::providerFrom<sbndDB::IPMTCalibrationDatabaseService const>();
  }

  void FlashGeoThreshold::InitializeFlashGeoAlgo()
  {
    // Initialize YZ map
    for(size_t opch=0; opch<::lightana::NOpDets(); opch++){
      ::lightana::OpDetCenterFromOpChannel(opch, fPDxyz);
      if( std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(opch)) == fPDTypes.end() ) continue;
      fYPEAccumulator[ (int) fPDxyz[1] ] = 0;
      fZPEAccumulator[ (int) fPDxyz[2] ] = 0;

      fNumONPMTsY_tpc0[(int) fPDxyz[1]] = 0;
      fNumONPMTsY_tpc1[(int) fPDxyz[1]] = 0;
      fNumONPMTsZ_tpc0[(int) fPDxyz[2]] = 0;
      fNumONPMTsZ_tpc1[(int) fPDxyz[2]] = 0;

      if(fNormalizeByPDType){
        if( std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(opch)) == fPDTypes.end() ) continue;
        fPDTypeByY[ (int) fPDxyz[1] ] = fPDSMap.pdType(opch);
        fPDTypeByZ[ (int) fPDxyz[2] ] = fPDSMap.pdType(opch);
      }
    }

    // Initialize accumulators for ON PMTs
    for(size_t opch=0; opch<::lightana::NOpDets(); opch++){
      ::lightana::OpDetCenterFromOpChannel(opch, fPDxyz);
      if(!fPMTCalibrationDatabaseService->getReconstructChannel(opch)) continue;
      if(fPDxyz[0]<0 && (fPDSMap.pdType(opch)=="pmt_coated" || fPDSMap.pdType(opch)=="pmt_uncoated")) {
        fNumONPMTsY_tpc0[(int) fPDxyz[1]]++;
        fNumONPMTsZ_tpc0[(int) fPDxyz[2]]++;
      }
      else if(fPDxyz[0]>0 && (fPDSMap.pdType(opch)=="pmt_coated" || fPDSMap.pdType(opch)=="pmt_uncoated")) {
        fNumONPMTsY_tpc1[(int) fPDxyz[1]]++;
        fNumONPMTsZ_tpc1[(int) fPDxyz[2]]++;
      }
    }

    // Initialize normalization factors by PD type
    if(fNormalizeByPDType){
      for(auto & pd:fPDTypes) {
        fNormFactorsY[pd]=0.;
        fNormFactorsZ[pd]=0.;
      }
    }
  
  }

  void FlashGeoThreshold::GetFlashLocation(std::vector<double> pePerOpChannel,
                                            double& Ycenter, double& Zcenter,
                                            double& Ywidth, double& Zwidth)
  {
    // Reset variables
    ResetVars();
    Ycenter = Zcenter = fDefCenterValue;
    Ywidth  = Zwidth  = fDefCenterValue;

    // Fill PE accumulators
    for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {
      if( std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(opch)) == fPDTypes.end() ) continue;
      // Get physical detector location for this opChannel
      ::lightana::OpDetCenterFromOpChannel(opch, fPDxyz);
      if(opch%2==0){
        fYPEAccumulator[(int) fPDxyz[1] ]+=pePerOpChannel[opch]/fNumONPMTsY_tpc0[(int) fPDxyz[1]];
        fZPEAccumulator[(int) fPDxyz[2] ]+=pePerOpChannel[opch]/fNumONPMTsZ_tpc0[(int) fPDxyz[2]];
      }
      else{
        fYPEAccumulator[(int) fPDxyz[1] ]+=pePerOpChannel[opch]/fNumONPMTsY_tpc1[(int) fPDxyz[1]];
        fZPEAccumulator[(int) fPDxyz[2] ]+=pePerOpChannel[opch]/fNumONPMTsZ_tpc1[(int) fPDxyz[2]];
      }
    }

    // Normalize PE accumulators
    if(fNormalizeByPDType){
      
      // Get normalization values for each PD type (Z)
      for(auto & z: fZPEAccumulator){
        if(z.second>fNormFactorsZ[ fPDTypeByZ[z.first] ]) 
          fNormFactorsZ[ fPDTypeByZ[z.first] ] = z.second;
      }

      // Get normalization values for each PD type (Y)
      for(auto & y: fYPEAccumulator){
        if(y.second>fNormFactorsY[ fPDTypeByY[y.first] ]) {
          fNormFactorsY[ fPDTypeByY[y.first] ] = y.second;
        }
      }

      for(auto & y: fYPEAccumulator)
        fYPEAccumulator[y.first] = y.second/ fNormFactorsY[ fPDTypeByY[y.first] ];

      for(auto & z: fZPEAccumulator)
        fZPEAccumulator[z.first] = z.second/ fNormFactorsZ[ fPDTypeByZ[z.first] ];

    }

    else{

      double normY = std::max_element(fYPEAccumulator.begin(), fYPEAccumulator.end(), 
                      [](const std::pair<int, double> &pe1, const std::pair<int, double> &pe2) 
                      {return pe1.second < pe2.second;})->second;
      
      double normZ =  std::max_element(fZPEAccumulator.begin(), fZPEAccumulator.end(), 
                      [](const std::pair<int, double> &pe1, const std::pair<int, double> &pe2) 
                      {return pe1.second < pe2.second;})->second;

      for(auto & y: fYPEAccumulator) fYPEAccumulator[y.first] = y.second / normY; 
      for(auto & z: fZPEAccumulator) fZPEAccumulator[z.first] = z.second / normZ;
 
    }

    double maxVal = 0.0;
    for (const auto& pair : fYPEAccumulator) {
        if (pair.second > maxVal)
            maxVal = pair.second;
    }
    if (maxVal > 0.0) {
      for (auto& pair : fYPEAccumulator) {
            pair.second /= maxVal;
      }
    }

    maxVal = 0.0;
    for (const auto& pair : fZPEAccumulator) {
        if (pair.second > maxVal)
            maxVal = pair.second;
    }
    if (maxVal > 0.0) {
      for (auto& pair : fZPEAccumulator) {
            pair.second /= maxVal;
      }
    }


    // Get YZ position of selected channels (above threshold)
    GetCenter(fYPEAccumulator, Ycenter, Ywidth, fThresholdY);
    GetCenter(fZPEAccumulator, Zcenter, Zwidth, fThresholdZ);

  }

  void FlashGeoThreshold::ResetVars(){
    // Reset PE accumulators
    for(auto & y: fYPEAccumulator) fYPEAccumulator[y.first]=0;
    for(auto & z: fZPEAccumulator) fZPEAccumulator[z.first]=0;
    // Reset normalization
    for(auto & pd:fPDTypes) {
      fNormFactorsY[pd]=0.;
      fNormFactorsZ[pd]=0.;
    }
  }

  void FlashGeoThreshold::GetCenter(std::map<int, double> PEAcc, double& center, double& width, double& threshold){

    // set variables
    center = 0.;
    double weight;
    double weightNormCenter = 0.;
    double weightNormWidth = 0.;
    double sum=0;
    double sum2=0;

    for(auto & pe: PEAcc){

      // Get weight for this channel
      if(fWeightExp==1) weight = pe.second;
      else if(fWeightExp==2) weight = pe.second*pe.second;
      else weight = std::pow(pe.second, fWeightExp);
      
      // For OpFlash width consider all channels
      weightNormWidth+=weight;
      sum+=weight*pe.first;
      sum2+=weight*pe.first*pe.first;

      // For OpFlash center consider only channels above ceertain threshold
      if(pe.second>threshold){  
        center+=weight*pe.first;
        weightNormCenter+=weight;
      }

    }

    center = center/weightNormCenter;
    width = std::sqrt(sum2*weightNormWidth - sum*sum)/weightNormWidth;

  }

}

DEFINE_ART_CLASS_TOOL(lightana::FlashGeoThreshold)
