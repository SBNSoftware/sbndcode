////////////////////////////////////////////////////////////////////////
//// File:        DigiPMTSBNDAlg.hh
////
//// This algorithm is used for the electronic response of PMTs
//// Created by L. Paulucci, F. Marinho, and I.L. de Icaza
//// Based on OpDetDigitizerDUNE_module.cc
//////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_DIGIPMTSBNDALG_HH
#define SBND_OPDETSIM_DIGIPMTSBNDALG_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandExponential.h"
#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "sbndcode/OpDetSim/PMTAlg/PMTGainFluctuations.hh"
#include "sbndcode/OpDetSim/PMTAlg/PMTNonLinearity.hh"
#include "sbndcode/OpDetSim/HDWvf/HDOpticalWaveforms.hh"

#include "TFile.h"

namespace opdet {

  class DigiPMTSBNDAlg {

  public:

    struct ConfigurationParameters_t {
      double TransitTime; //ns
      double TTS; //ns
      double CableTime; //ns
      double PMTChargeToADC; //voltage to ADC conversion scale
      double PMTBaseline; //waveform baseline in ADC
      double PMTFallTime; //pulse decaying time constant (exponential) in ns
      double PMTRiseTime; //pulse rising time constant (exponential) in ns
      double PMTMeanAmplitude; //mean amplitude for single pe in pC
      double PMTBaselineRMS; //Pedestal RMS in ADC counts
      double PMTDarkNoiseRate; //in Hz
      double PMTADCDynamicRange; //ADC dynbamic range
      double PMTCoatedVUVEff; //PMT (coated) efficiency for direct (VUV) light
      double PMTCoatedVISEff; //PMT (coated) efficiency for reflected (VIS) light
      double PMTUncoatedEff; //PMT (uncoated) efficiency
      std::string PMTDataFile; //File containing timing emission structure for TPB, and single PE profile from data
      bool PMTSinglePEmodel; //Model for single pe response, false for ideal, true for test bench meas
      bool MakeGainFluctuations; //Fluctuate PMT gain
      fhicl::ParameterSet GainFluctuationsParams;
      bool SimulateNonLinearity; //Fluctuate PMT gain
      fhicl::ParameterSet NonLinearityParams;
      
      fhicl::ParameterSet HDOpticalWaveformParams;

      detinfo::LArProperties const* larProp = nullptr; //< LarProperties service provider.
      double frequency;       //wave sampling frequency (GHz)
      CLHEP::HepRandomEngine* engine = nullptr;
    };// ConfigurationParameters_t

    //Default constructor
    DigiPMTSBNDAlg(ConfigurationParameters_t const& config);
    //Default destructor
    ~DigiPMTSBNDAlg();

    void ConstructWaveformUncoatedPMT(
      int ch,
      sim::SimPhotons const& simphotons,
      std::vector<short unsigned int>& waveform,
      std::string pdtype,
      double start_time,
      unsigned n_sample);

    void ConstructWaveformCoatedPMT(
      int ch,
      std::vector<short unsigned int>& waveform,
      std::unordered_map<int, sim::SimPhotons>& DirectPhotonsMap,
      std::unordered_map<int, sim::SimPhotons>& ReflectedPhotonsMap,
      double start_time,
      unsigned n_sample);

    void ConstructWaveformLiteUncoatedPMT(
      int ch,
      sim::SimPhotonsLite const& litesimphotons,
      std::vector<short unsigned int>& waveform,
      std::string pdtype,
      double start_time,
      unsigned n_sample);

    void ConstructWaveformLiteCoatedPMT(
      int ch,
      std::vector<short unsigned int>& waveform,
      std::unordered_map<int, sim::SimPhotonsLite>& DirectPhotonsMap,
      std::unordered_map<int, sim::SimPhotonsLite>& ReflectedPhotonsMap,
      double start_time,
      unsigned n_sample);

    double Baseline()
      {
        return fParams.PMTBaseline;
      }

  private:

    ConfigurationParameters_t fParams;
    // Declare member data here.
    double fSampling;       //wave sampling frequency (GHz)
    double fSamplingPeriod; //wave sampling period (ns)
    double fPMTCoatedVUVEff;
    double fPMTCoatedVISEff;
    double fPMTUncoatedEff;
    bool fPositivePolarity;
    int fADCSaturation;

    double sigma1;
    double sigma2;

    const double transitTimeSpread_frac = 2.0 * std::sqrt(2.0 * std::log(2.0));

    CLHEP::HepRandomEngine* fEngine; //!< Reference to art-managed random-number engine
    CLHEP::RandFlat fFlatGen;
    CLHEP::RandPoissonQ fPoissonQGen;
    CLHEP::RandGaussQ fGaussQGen;
    CLHEP::RandExponential fExponentialGen;
    std::unique_ptr<CLHEP::RandGeneral> fTimeTPB; // histogram for getting the TPB emission time for coated PMTs

    //PMTFluctuationsAlg
    std::unique_ptr<opdet::PMTGainFluctuations> fPMTGainFluctuationsPtr;
    //HDWaveforms
    std::unique_ptr<opdet::HDOpticalWaveform> fPMTHDOpticalWaveformsPtr;

    //PMTNonLinearity
    std::unique_ptr<opdet::PMTNonLinearity> fPMTNonLinearityPtr;

    void AddSPE(size_t time, std::vector<double>& wave, double npe = 1); // add single pulse to auxiliary waveform
    void Pulse1PE(std::vector<double>& wave);
    double Transittimespread(double fwhm);

    std::vector<double> fSinglePEWave; // single photon pulse vector
    std::vector<std::vector<double>> fSinglePEWave_HD; // single photon pulse vector
    int pulsesize; //size of 1PE waveform
    std::unordered_map< raw::Channel_t, std::vector<double> > fFullWaveforms;

    void CreatePDWaveformUncoatedPMT(
      sim::SimPhotons const& SimPhotons,
      double t_min,
      std::vector<double>& wave,
      int ch,
      std::string pdtype);
    void CreatePDWaveformCoatedPMT(
      int ch,
      double t_min,
      std::vector<double>& wave,
      std::unordered_map<int, sim::SimPhotons>& DirectPhotonsMap,
      std::unordered_map<int, sim::SimPhotons>& ReflectedPhotonsMap);
    void CreatePDWaveformLiteUncoatedPMT(
      sim::SimPhotonsLite const& litesimphotons,
      double t_min,
      std::vector<double>& wave,
      int ch,
      std::string pdtype);
    void CreatePDWaveformLiteCoatedPMT(
      int ch,
      double t_min,
      std::vector<double>& wave,
      std::unordered_map<int, sim::SimPhotonsLite>& DirectPhotonsMap,
      std::unordered_map<int, sim::SimPhotonsLite>& ReflectedPhotonsMap);
    void CreateSaturation(std::vector<double>& wave);//Including saturation effects (dynamic range)
    void AddLineNoise(std::vector<double>& wave); //add noise to baseline
    void AddDarkNoise(std::vector<double>& wave); //add dark noise
    double FindMinimumTime(
      sim::SimPhotons const&,
      int ch,
      std::string pdtype,
      std::unordered_map<int, sim::SimPhotons>& directPhotonsOnPMTS);
    double FindMinimumTimeLite(
      sim::SimPhotonsLite const& litesimphotons,
      int ch,
      std::string pdtype,
      std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS);
  };//class DigiPMTSBNDAlg

  class DigiPMTSBNDAlgMaker {

  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> transitTime {
        Name("TransitTime"),
        Comment("Single pe: Time of maximum amplitude in the PMT pulse in ns")
      };

      fhicl::Atom<double> tts {
        Name("TTS"),
        Comment("Single pe: Transit time spread in ns")
      };

      fhicl::Atom<double> cableTime {
        Name("CableTime"),
        Comment("Time delay of the 30 m long readout cable in ns")
      };

      fhicl::Atom<double> pmtmeanAmplitude {
        Name("PMTMeanAmplitude"),
        Comment("Single pe: mean amplitude of PMT pulse in pC")
      };

      fhicl::Atom<double> pmtriseTime {
        Name("PMTRiseTime"),
        Comment("Single pe: Pulse rise time constant (exponential), from 0.1 to 0.9 of maximum amplitude")
      };

      fhicl::Atom<double> pmtfallTime {
        Name("PMTFallTime"),
        Comment("Single pe: Pulse decay time constant (exponential), from 0.1 to 0.9 of maximum amplitude")
      };

      fhicl::Atom<double> pmtchargeToADC {
        Name("PMTChargeToADC"),
        Comment("Charge to ADC convertion factor")
      };

      fhicl::Atom<double> pmtbaseline {
        Name("PMTBaseline"),
        Comment("Waveform baseline in ADC")
      };

      fhicl::Atom<double> pmtbaselineRMS {
        Name("PMTBaselineRMS"),
        Comment("RMS of the electronics noise fluctuations in ADC counts")
      };

      fhicl::Atom<double> pmtdarkNoiseRate {
        Name("PMTDarkNoiseRate"),
        Comment("Dark noise rate in Hz")
      };

      fhicl::Atom<double> pmtADCDynamicRange {
        Name("PMTADCDynamicRange"),
        Comment("Saturation in number of ADCs")
      };

      fhicl::Atom<double> pmtcoatedVUVEff {
        Name("PMTCoatedVUVEff"),
        Comment("PMT (coated) detection efficiency for direct (VUV) light")
      };

      fhicl::Atom<double> pmtcoatedVISEff {
        Name("PMTCoatedVISEff"),
        Comment("PMT (coated) detection efficiency for reflected (VIS) light")
      };

      fhicl::Atom<double> pmtuncoatedEff {
        Name("PMTUncoatedEff"),
        Comment("PMT (uncoated) detection efficiency")
      };

      fhicl::Atom<bool> PMTsinglePEmodel {
        Name("PMTSinglePEmodel"),
        Comment("Model used for single PE response of PMT. =0 is ideal, =1 is testbench")
      };

      fhicl::Atom<std::string> pmtDataFile {
        Name("PMTDataFile"),
        Comment("File containing timing emission distribution for TPB and single pe pulse from data")
      };

      fhicl::OptionalDelegatedParameter gainFluctuationsParams {
        Name("GainFluctuationsParams"),
        Comment("Parameters used for SinglePE response fluctuations")
      };

      fhicl::OptionalDelegatedParameter nonLinearityParams {
        Name("NonLinearityParams"),
        Comment("Parameters used for simulating PMT non linear effects")
      };
      
      fhicl::OptionalDelegatedParameter hdOpticalWaveformParams {
        Name("HDOpticalWaveformParamsPMT"),
        Comment("Parameters used for high definition waveform")
      };

    };    //struct Config

    DigiPMTSBNDAlgMaker(Config const& config); //Constructor

    std::unique_ptr<DigiPMTSBNDAlg> operator()(
      detinfo::LArProperties const& larProp,
      detinfo::DetectorClocksData const& clockData,
      CLHEP::HepRandomEngine* engine
      ) const;

  private:
    // Part of the configuration learned from configuration files.
    DigiPMTSBNDAlg::ConfigurationParameters_t fBaseConfig;
  }; //class DigiPMTSBNDAlgMaker

} // namespace opdet

#endif //SBND_OPDETSIM_DIGIPMTSBNDALG_HH
