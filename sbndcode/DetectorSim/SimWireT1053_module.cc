////////////////////////////////////////////////////////////////////////
// $Id: SimWireT1053.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireT1053 class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
// - ported from the MicroBooNE class by A.Szlec
////////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "larsim/RandomUtils/LArSeedService.h"


#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/TriggerData.h"
// #include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "sbndcode/Utilities/SignalShapingServiceT1053.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

///Detector simulation of raw signals on wires
namespace detsim {

// Base class for creation of raw signals on wires.
class SimWireT1053 : public art::EDProducer {

public:

  explicit SimWireT1053(fhicl::ParameterSet const& pset);
  virtual ~SimWireT1053();

  // read/write access to event
  void produce (art::Event& evt);
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const& p);

private:

  void GenNoiseInTime(std::vector<float> &noise);
  void GenNoiseInFreq(std::vector<float> &noise);

  std::string            fDriftEModuleLabel;///< module making the ionization electrons
  raw::Compress_t        fCompression;      ///< compression type to use

  double                 fNoiseFact;        ///< noise scale factor (ADC RMS for gaussian noise)
  double                 fNoiseWidth;       ///< exponential noise width (kHz)
  double                 fNoiseRand;        ///< fraction of random "wiggle" in noise in freq. spectrum
  double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
  size_t                 fNTicks;           ///< number of ticks of the clock
  double                 fSampleRate;       ///< sampling rate in ns
  unsigned int           fNTimeSamples;     ///< number of ADC readout samples in all readout frames (per event)
  float                  fCollectionPed;    ///< ADC value of baseline for collection plane
  float                  fInductionPed;     ///< ADC value of baseline for induction plane
  float                  fBaselineRMS;      ///< ADC value of baseline RMS within each channel
  TH1D*                  fNoiseDist;        ///< distribution of noise counts
  bool fGetNoiseFromHisto;                  ///< if True -> Noise from Histogram of Freq. spectrum
  bool fGenNoiseInTime;                     ///< if True -> Noise with Gaussian dsitribution in Time-domain
  bool fGenNoise;                           ///< if True -> Gen Noise. if False -> Skip noise generation entierly
  std::string fNoiseFileFname;
  std::string fNoiseHistoName;
  TH1D*                  fNoiseHist = nullptr; ///< distribution of noise counts

  std::string fTrigModName;                 ///< Trigger data product producer name
  //define max ADC value - if one wishes this can
  //be made a fcl parameter but not likely to ever change
  const float adcsaturation = 4095;

  ::detinfo::ElecClock fClock; ///< TPC electronics clock

}; // class SimWireT1053

DEFINE_ART_MODULE(SimWireT1053)

//-------------------------------------------------
SimWireT1053::SimWireT1053(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);

  produces< std::vector<raw::RawDigit>   >();

  fCompression = raw::kNone;
  TString compression(pset.get< std::string >("CompressionType"));
  if (compression.Contains("Huffman", TString::kIgnoreCase)) fCompression = raw::kHuffman;

  // create a default random engine; obtain the random seed from LArSeedService,
  // unless overridden in configuration with key "Seed" and "SeedPedestal"
  art::ServiceHandle<sim::LArSeedService> Seeds;
  Seeds->createEngine(*this, "HepJamesRandom", "noise", pset, "Seed");
  Seeds->createEngine(*this, "HepJamesRandom", "pedestal", pset, "SeedPedestal");


}

//-------------------------------------------------
SimWireT1053::~SimWireT1053()
{
  delete fNoiseHist;
}

//-------------------------------------------------
void SimWireT1053::reconfigure(fhicl::ParameterSet const& p)
{
  fDriftEModuleLabel = p.get< std::string         >("DriftEModuleLabel");
  fNoiseFact        = p.get< double              >("NoiseFact");
  fNoiseWidth       = p.get< double              >("NoiseWidth");
  fNoiseRand        = p.get< double              >("NoiseRand");
  fLowCutoff        = p.get< double              >("LowCutoff");
  fGetNoiseFromHisto = p.get< bool                >("GetNoiseFromHisto");
  fGenNoiseInTime   = p.get< bool                >("GenNoiseInTime");
  fGenNoise         = p.get< bool                >("GenNoise");
  fCollectionPed    = p.get< float               >("CollectionPed");
  fInductionPed     = p.get< float               >("InductionPed");
  fBaselineRMS      = p.get< float               >("BaselineRMS");

  fTrigModName      = p.get< std::string         >("TrigModName");

  if (fGetNoiseFromHisto)
  {
    fNoiseHistoName = p.get< std::string         >("NoiseHistoName");

    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(p.get<std::string>("NoiseFileFname"), fNoiseFileFname);

    TFile in(fNoiseFileFname.c_str(), "READ");
    if (!in.IsOpen()) {
      throw art::Exception(art::errors::FileOpenError)
        << "Could not find Root file '" << fNoiseHistoName
        << "' for noise histogram\n";
    }
    fNoiseHist = (TH1D*) in.Get(fNoiseHistoName.c_str());

    if (!fNoiseHist) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find noise histogram '" << fNoiseHistoName
        << "' in Root file '" << in.GetPath() << "'\n";
    }
    // release the histogram from its file; now we own it
    fNoiseHist->SetDirectory(nullptr);
    in.Close();
  }
  
  //detector properties information
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fNTimeSamples  = detprop->NumberTimeSamples();    



  return;
}

//-------------------------------------------------
void SimWireT1053::beginJob()
{

  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fNoiseDist  = tfs->make<TH1D>("Noise", ";Noise  (ADC);", 1000,   -10., 10.);

  art::ServiceHandle<util::LArFFT> fFFT;
  fNTicks = fFFT->FFTSize();

  if ( fNTicks % 2 != 0 )
    LOG_DEBUG("SimWireT1053") << "Warning: FFTSize not a power of 2. "
                              << "May cause issues in (de)convolution.\n";

  if ( fNTimeSamples > fNTicks )
    mf::LogError("SimWireT1053") << "Cannot have number of readout samples "
                                 << "greater than FFTSize!";

  return;

}

//-------------------------------------------------
void SimWireT1053::endJob()
{}

void SimWireT1053::produce(art::Event& evt)
{

  // auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
  art::ServiceHandle<detinfo::DetectorClocksServiceStandard> tss;
  // In case trigger simulation is run in the same job...
  tss->preProcessEvent(evt);
  auto const* ts = tss->provider();

  // get the geometry to be able to figure out signal types and chan -> plane mappings
  art::ServiceHandle<geo::Geometry> geo;
  //unsigned int signalSize = fNTicks;

  std::vector<const sim::SimChannel*> chanHandle;
  evt.getView(fDriftEModuleLabel, chanHandle);

  //Get fIndShape and fColShape from SignalShapingService, on the fly
  art::ServiceHandle<util::SignalShapingServiceT1053> sss;

  // make a vector of const sim::SimChannel* that has same number
  // of entries as the number of channels in the detector
  // and set the entries for the channels that have signal on them
  // using the chanHandle
  std::vector<const sim::SimChannel*> channels(geo->Nchannels(), nullptr);
  for (size_t c = 0; c < chanHandle.size(); ++c) {
    channels.at(chanHandle.at(c)->Channel()) = chanHandle.at(c);
  }

  const auto NChannels = geo->Nchannels();

  // vectors for working
  std::vector<short>    adcvec(fNTimeSamples, 0);
  std::vector<double>   chargeWork(fNTicks, 0.);



  // make a unique_ptr of sim::SimDigits that allows ownership of the produced
  // digits to be transferred to the art::Event after the put statement below
  std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
  digcol->reserve(NChannels);

  unsigned int chan = 0;
  art::ServiceHandle<util::LArFFT> fFFT;

  //LOOP OVER ALL CHANNELS
  std::map<int, double>::iterator mapIter;
  for (chan = 0; chan < geo->Nchannels(); chan++) {

    // get the sim::SimChannel for this channel
    const sim::SimChannel* sc = channels.at(chan);
    std::fill(chargeWork.begin(), chargeWork.end(), 0.);
    if ( sc ) {

      // loop over the tdcs and grab the number of electrons for each
      for (int t = 0; t < (int)(chargeWork.size()); ++t) {

        int tdc = ts->TPCTick2TDC(t);

        // continue if tdc < 0
        if ( tdc < 0 ) continue;

        chargeWork.at(t) = sc->Charge(tdc);

      }

      // Convolve charge with appropriate response function
      sss->Convolute(chan, chargeWork);

    }

    //Generate Noise:
    std::vector<float> noisetmp(fNTicks, 0.);
    if (fGenNoise) {
      if (fGenNoiseInTime)
        GenNoiseInTime(noisetmp);
      else
        GenNoiseInFreq(noisetmp);
    }

    //Pedestal determination
    float ped_mean = fCollectionPed;
    geo::SigType_t sigtype = geo->SignalType(chan);
    if (sigtype == geo::kInduction)
      ped_mean = fInductionPed;
    else if (sigtype == geo::kCollection)
      ped_mean = fCollectionPed;
    //slight variation on ped on order of RMS of baselien variation
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("pedestal");
    CLHEP::RandGaussQ rGaussPed(engine, 0.0, fBaselineRMS);
    ped_mean += rGaussPed.fire();

    for (unsigned int i = 0; i < fNTimeSamples; ++i) {
      float adcval = noisetmp.at(i) + chargeWork.at(i) + ped_mean;

      //Add Noise to NoiseDist Histogram
      if (i % 100 == 0)
        fNoiseDist->Fill(noisetmp.at(i));

      //allow for ADC saturation
      if ( adcval > adcsaturation )
        adcval = adcsaturation;
      //don't allow for "negative" saturation
      if ( adcval < 0 )
        adcval = 0;

      adcvec.at(i) = (unsigned short)(adcval);

    }// end loop over signal size

    // resize the adcvec to be the correct number of time samples,
    // just drop the extra samples
    //adcvec.resize(fNTimeSamples);

    // compress the adc vector using the desired compression scheme,
    // if raw::kNone is selected nothing happens to adcvec
    // This shrinks adcvec, if fCompression is not kNone.
    raw::Compress(adcvec, fCompression);

    // add this digit to the collection
    raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
    rd.SetPedestal(ped_mean);
    digcol->push_back(rd);

  }// end loop over channels

  evt.put(std::move(digcol));
  return;
}

//-------------------------------------------------
void SimWireT1053::GenNoiseInTime(std::vector<float> &noise)
{
  //ART random number service
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
  CLHEP::RandGaussQ rGauss(engine, 0.0, fNoiseFact);

  //In this case fNoiseFact is a value in ADC counts
  //It is going to be the Noise RMS
  //loop over all bins in "noise" vector
  //and insert random noise value
  for (unsigned int i = 0; i < noise.size(); i++)
    noise.at(i) = rGauss.fire();
}


//-------------------------------------------------
void SimWireT1053::GenNoiseInFreq(std::vector<float> &noise)
{
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
  CLHEP::RandFlat flat(engine, -1, 1);

  if (noise.size() != fNTicks)
    throw cet::exception("SimWireT1053")
        << "\033[93m"
        << "Frequency noise vector length must match fNTicks (FFT size)"
        << " ... " << noise.size() << " != " << fNTicks
        << "\033[00m"
        << std::endl;

  // noise in frequency space
  std::vector<TComplex> noiseFrequency(fNTicks / 2 + 1, 0.);

  double pval = 0.;
  double lofilter = 0.;
  double phase = 0.;
  double rnd[2] = {0.};

  // width of frequencyBin in kHz
  double binWidth = 1.0 / (fNTicks * fSampleRate * 1.0e-6);

  for (size_t i = 0; i < fNTicks / 2 + 1; ++i) {
    // exponential noise spectrum
    flat.fireArray(2, rnd, 0, 1);
    //if not from histo or in time --> then hardcoded freq. spectrum
    if ( !fGetNoiseFromHisto )
    {
      pval = fNoiseFact * exp(-(double)i * binWidth / fNoiseWidth);
      // low frequency cutoff
      lofilter = 1.0 / (1.0 + exp(-(i - fLowCutoff / binWidth) / 0.5));
      // randomize 10%

      pval *= lofilter * ((1 - fNoiseRand) + 2 * fNoiseRand * rnd[0]);
    }


    else
    {

      pval = fNoiseHist->GetBinContent(i) * ((1 - fNoiseRand) + 2 * fNoiseRand * rnd[0]) * fNoiseFact;
      //mf::LogInfo("SimWireT1053")  << " pval: " << pval;
    }

    phase = rnd[1] * 2.*TMath::Pi();
    TComplex tc(pval * cos(phase), pval * sin(phase));
    noiseFrequency.at(i) += tc;
  }


  // mf::LogInfo("SimWireT1053") << "filled noise freq";

  // inverse FFT MCSignal
  art::ServiceHandle<util::LArFFT> fFFT;
  fFFT->DoInvFFT(noiseFrequency, noise);

  noiseFrequency.clear();

  // multiply each noise value by fNTicks as the InvFFT
  // divides each bin by fNTicks assuming that a forward FFT
  // has already been done.
  for (unsigned int i = 0; i < noise.size(); ++i) noise.at(i) *= 1.*fNTicks;

}


}

