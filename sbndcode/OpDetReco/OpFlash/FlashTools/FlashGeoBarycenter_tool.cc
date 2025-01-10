///////////////////////////////////////////////////////////////////////
/// File: FlashGeoBarycenter_tool.cc
///
/// Base class: FlashGeoBase.hh
///
/// It computes the PMTs barycenter
/// weighted by the reconstructed number of PE
///
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Utilities/ToolConfigTable.h"

#include "FlashGeoBase.hh"

namespace lightana{

  class FlashGeoBarycenter : FlashGeoBase
  {

  public:

    //Configuration parameters
     struct Config {

       fhicl::Atom<unsigned int> WeightExp {
         fhicl::Name("WeightExp"),
         fhicl::Comment("Weight exponent for YZ barycenters")
       };

     };


    // Default constructor
    explicit FlashGeoBarycenter(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    void GetFlashLocation(std::vector<double> pePerOpChannel,
                            double& Ycenter, double& Zcenter,
                            double& Ywidth, double& Zwidth) override;

  private:

    unsigned int fWeightExp;

    static constexpr double fDefCenterValue = -999.;

  };


  FlashGeoBarycenter::FlashGeoBarycenter(art::ToolConfigTable<Config> const& config)
    : fWeightExp{ config().WeightExp() }
  {
  }


  void FlashGeoBarycenter::GetFlashLocation(std::vector<double> pePerOpChannel,
                                            double& Ycenter, double& Zcenter,
                                            double& Ywidth, double& Zwidth)
  {

    // Reset variables
    Ycenter = Zcenter = fDefCenterValue;
    Ywidth  = Zwidth  = fDefCenterValue;
    double totalPE = 0.;
    double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;
    double weight =1.;
    for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {
      // Get physical detector location for this opChannel
      double PMTxyz[3];
      ::lightana::OpDetCenterFromOpChannel(opch, PMTxyz);

      // Get weight for this channel
      if(fWeightExp==1) weight = pePerOpChannel[opch];
      else if(fWeightExp==2) weight = pePerOpChannel[opch]*pePerOpChannel[opch];
      else weight = std::pow(pePerOpChannel[opch], fWeightExp);

      // Add up the position, weighting with PEs
      sumy    += weight*PMTxyz[1];
      sumy2   += weight*PMTxyz[1]*PMTxyz[1];
      sumz    += weight*PMTxyz[2];
      sumz2   += weight*PMTxyz[2]*PMTxyz[2];
      totalPE += weight;
    }

    Ycenter = sumy/totalPE;
    Zcenter = sumz/totalPE;

    // This is just sqrt(<x^2> - <x>^2)
    if ( (sumy2*totalPE - sumy*sumy) > 0. )
      Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

    if ( (sumz2*totalPE - sumz*sumz) > 0. )
      Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
  }

}

DEFINE_ART_CLASS_TOOL(lightana::FlashGeoBarycenter)
