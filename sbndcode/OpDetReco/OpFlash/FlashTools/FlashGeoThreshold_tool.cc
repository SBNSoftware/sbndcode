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

#include <string>
#include <map>

namespace lightana{

  class FlashGeoThreshold : FlashGeoBase
  {

  public:

    //Configuration parameters
     struct Config {

       fhicl::Atom<float> Threshold {
         fhicl::Name("Threshold"),
         fhicl::Comment("Use channels abnove this thershold")
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

    // Method to calculate the OpFlash t0
    void GetFlashLocation(std::vector<double> pePerOpChannel,
                            double& Ycenter, double& Zcenter,
                            double& Ywidth, double& Zwidth) override;

  private:

    void ResetVars();
    void GetCenter(std::map<int, double> PEAcc, double& center, double& width);

    // Fhicl configuration parameters
    std::vector<std::string> fPDTypes;
    float fThreshold;
    bool fNormalizeByPDType;
    unsigned int fWeightExp;

    // Store #PE at each YZ position
    std::map<int, double> fYPEAccumulator;
    std::map<int, double> fZPEAccumulator;

    // Store PD type ar each YZ position
    std::map<int, std::string> fPDTypeByY;
    std::map<int, std::string> fPDTypeByZ;

    // Store normalization factor for each PD type
    std::map<std::string, double> fNormFactorsY;
    std::map<std::string, double> fNormFactorsZ;

    // PDS mapping
    opdet::sbndPDMapAlg fPDSMap;
    
    // Array to store PD xyz position
    double fPDxyz[3];

    static constexpr double fDefCenterValue = -999.;

  };


  FlashGeoThreshold::FlashGeoThreshold(art::ToolConfigTable<Config> const& config)
    : fPDTypes{ config().PDTypes() }
    , fThreshold{ config().Threshold() }
    , fNormalizeByPDType{ config().NormalizeByPDType() }
    , fWeightExp{ config().WeightExp() }
  {

    // Initialize YZ map
    for(size_t opch=0; opch<::lightana::NOpDets(); opch++){
      
      if( std::find(fPDTypes.begin(), fPDTypes.end(), fPDSMap.pdType(opch)) == fPDTypes.end() ) continue;
      
      ::lightana::OpDetCenterFromOpChannel(opch, fPDxyz);
      
      fYPEAccumulator[ (int) fPDxyz[1] ] = 0;
      fZPEAccumulator[ (int) fPDxyz[2] ] = 0;
      
      if(fNormalizeByPDType){
        fPDTypeByY[ (int) fPDxyz[1] ] = fPDSMap.pdType(opch);
        fPDTypeByZ[ (int) fPDxyz[2] ] = fPDSMap.pdType(opch);
      }

    }

    // Initialize normalixation factors by PD type
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
      fYPEAccumulator[(int) fPDxyz[1] ]+=pePerOpChannel[opch];
      fZPEAccumulator[(int) fPDxyz[2] ]+=pePerOpChannel[opch];
    }
    
    std::cout<<" YYY \n";
    for(auto & y: fYPEAccumulator){
      std::cout<<y.first<<":"<<y.second<<" ";
    }
    std::cout<<std::endl;

    std::cout<<" ZZZ \n";
    for(auto & z: fZPEAccumulator){
      std::cout<<z.first<<":"<<z.second<<" ";
    }
    std::cout<<std::endl;

    // Normalize PE accumulators
    if(fNormalizeByPDType){
      std::cout<<"Normalizing by PDType\n";

      // Get normalization values for each PD type
      for(auto & z: fZPEAccumulator){
        if(z.second>fNormFactorsZ[ fPDTypeByZ[z.first] ]) 
          fNormFactorsZ[ fPDTypeByZ[z.first] ] = z.second;
      }

      for(auto & y: fYPEAccumulator){
        if(y.second>fNormFactorsY[ fPDTypeByY[y.first] ]) 
          fNormFactorsY[ fPDTypeByY[y.first] ] = y.second;
      }

      std::cout<<" Norms are: \n";
      for(auto & pd:fPDTypes){
        std::cout<<"PDTypes: "<<pd<<" NormY: "<<fNormFactorsY[pd]<<" NormZ: "<<fNormFactorsZ[pd]<<std::endl;
      }

      for(auto & y: fYPEAccumulator)
        fYPEAccumulator[y.first] = y.second/ fNormFactorsY[ fPDTypeByY[y.first] ]; 
      
      for(auto & z: fZPEAccumulator)
        fZPEAccumulator[z.first] = z.second/ fNormFactorsZ[ fPDTypeByZ[z.first] ]; 

    }
    else{
      std::cout<<"Not normalizing by PDType\n";

      double normY = std::max_element(fYPEAccumulator.begin(), fYPEAccumulator.end(), 
                      [](const std::pair<int, double> &pe1, const std::pair<int, double> &pe2) 
                      {return pe1.second < pe2.second;})->second;
      
      double normZ =  std::max_element(fZPEAccumulator.begin(), fZPEAccumulator.end(), 
                      [](const std::pair<int, double> &pe1, const std::pair<int, double> &pe2) 
                      {return pe1.second < pe2.second;})->second;
      
        
      std::cout<<" Norms are: Y=" <<normY<<" Z="<<normZ<<"\n";

      for(auto & y: fYPEAccumulator) fYPEAccumulator[y.first] = y.second / normY; 
      for(auto & z: fZPEAccumulator) fZPEAccumulator[z.first] = z.second / normZ;
 
    }

    std::cout<<" After norm...\n YYY \n";
    for(auto & y: fYPEAccumulator){
      std::cout<<y.first<<":"<<y.second<<" ";
    }
    std::cout<<std::endl;
    std::cout<<" ZZZ \n";
    for(auto & z: fZPEAccumulator){
      std::cout<<z.first<<":"<<z.second<<" ";
    }
    std::cout<<std::endl;

    // Get YZ position of selected channels (above threshold)
    GetCenter(fYPEAccumulator, Ycenter, Ywidth);
    GetCenter(fZPEAccumulator, Zcenter, Zwidth);

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

  void FlashGeoThreshold::GetCenter(std::map<int, double> PEAcc, double& center, double& width){

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
      if(pe.second>fThreshold){  
        center+=weight*pe.first;
        weightNormCenter+=weight;
      }

    }

    center = center/weightNormCenter;
    width = std::sqrt(sum2*weightNormWidth - sum*sum)/weightNormWidth;

  }

}

DEFINE_ART_CLASS_TOOL(lightana::FlashGeoThreshold)
