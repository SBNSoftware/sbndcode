// SBNDThermalNoiseServiceInFreq.cxx
// Andrew Scarff
// July 2019
// Based upon SPhaseChannelNoiseService.cxx developed by Jingbo Wang for ProtoDUNE.

#include "sbndcode/DetectorSim/Services/SBNDThermalNoiseServiceInFreq.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

SBNDThermalNoiseServiceInFreq::
SBNDThermalNoiseServiceInFreq(fhicl::ParameterSet const& pset)
  : fRandomSeed(0), fLogLevel(1),  m_pran(nullptr), fNoiseEngine(nullptr)
{
  const string myname = "SBNDThermalNoiseServiceInFreq::ctor: ";
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints");
  bool haveSeed      = pset.get_if_present<int>("RandomSeed", fRandomSeed);


  fShapingTimeOrder = { {0.5, 0 }, {1.0, 1}, {2.0, 2}, {3.0, 3} };
  fNoiseWidth        = pset.get< double              >("NoiseWidth");
  fNoiseRand         = pset.get< double              >("NoiseRand");
  fLowCutoff         = pset.get< double              >("LowCutoff");


  if ( fRandomSeed == 0 ) haveSeed = false;
  pset.get_if_present<int>("LogLevel", fLogLevel);
  int seed = fRandomSeed;
  
  string rname = "SBNDThermalNoiseServiceInFreq";
  if ( haveSeed ) {
    if ( fLogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(seed);
  } else {
    if ( fLogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( fLogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( fLogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;
  generateNoise();
  if ( fLogLevel > 1 ) print() << endl;
}

//**********************************************************************

SBNDThermalNoiseServiceInFreq::
SBNDThermalNoiseServiceInFreq(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: SBNDThermalNoiseServiceInFreq(pset) { }

//**********************************************************************

SBNDThermalNoiseServiceInFreq::~SBNDThermalNoiseServiceInFreq() {
  const string myname = "SBNDThermalNoiseServiceInFreq::dtor: ";
  if ( fLogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int SBNDThermalNoiseServiceInFreq::addNoise(Channel chan, AdcSignalVector& sigs) const {

  //Get services.
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::SignalShapingServiceSBND> sss;
  art::ServiceHandle<util::LArFFT> fFFT;
  
  //Generate Noise:
  size_t view = (size_t)geo->View(chan);
  
  double noise_factor;
  auto tempNoiseVec = sss->GetNoiseFactVec();
  double shapingTime = 2.0; //sss->GetShapingTime(chan);
  double asicGain = sss->GetASICGain(chan);

  size_t fNTicks = fFFT->FFTSize();

  
  if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() ) {
    noise_factor = tempNoiseVec[view].at( fShapingTimeOrder.find( shapingTime )->second );
    noise_factor *= asicGain/4.7;
  }
  else {
    throw cet::exception("SBNDThermalNoiseServiceInFreq_service.cc")
      << "\033[93m"
      << "Shaping Time recieved from signalshapingservices_sbnd.fcl is not one of the allowed values"
      << std::endl
      << "Allowed values: 0.5, 1.0, 2.0, 3.0 us"
      << "\033[00m"
      << std::endl;
  }

  
  CLHEP::RandFlat flat(*fNoiseEngine, -1, 1);


  if (sigs.size() != fNTicks)
    throw cet::exception("SBNDThermalNoiseServiceInFreq_service.cc")
        << "\033[93m"
        << "Frequency noise vector length must match fNTicks (FFT size)"
        << " ... " << sigs.size() << " != " << fNTicks
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
    std::cout << "rnd[0]: " << rnd[0] << std::endl;
    std::cout << "rnd[1]: " << rnd[1] << std::endl;

    pval = noise_factor * exp(-(double)i * binWidth / fNoiseWidth);
    // low frequency cutoff
    lofilter = 1.0 / (1.0 + exp(-(i - fLowCutoff / binWidth) / 0.5));
    // randomize 10%

    pval *= lofilter * ((1 - fNoiseRand) + 2 * fNoiseRand * rnd[0]);

    /*
    //if not from histo or in time --> then hardcoded freq. spectrum
    if ( !fGetNoiseFromHisto )
    {
      pval = noise_factor * exp(-(double)i * binWidth / fNoiseWidth);
      // low frequency cutoff
      lofilter = 1.0 / (1.0 + exp(-(i - fLowCutoff / binWidth) / 0.5));
      // randomize 10%

      pval *= lofilter * ((1 - fNoiseRand) + 2 * fNoiseRand * rnd[0]);
    }


    else
    {

      pval = fNoiseHist->GetBinContent(i) * ((1 - fNoiseRand) + 2 * fNoiseRand * rnd[0]) * noise_factor;
    }
    */
    phase = rnd[1] * 2.*TMath::Pi();
    TComplex tc(pval * cos(phase), pval * sin(phase));
    noiseFrequency.at(i) += tc;
  }


  // inverse FFT MCSignal
  fFFT->DoInvFFT(noiseFrequency, sigs);

  noiseFrequency.clear();

  // multiply each noise value by fNTicks as the InvFFT
   // divides each bin by fNTicks assuming that a forward FFT
  // has already been done.
  for (unsigned int i = 0; i < sigs.size(); ++i) {
    sigs.at(i) *= 1.*fNTicks;
  }
  
  return 0;
}


//**********************************************************************

ostream& SBNDThermalNoiseServiceInFreq::print(ostream& out, string prefix) const {
  out << prefix << "SBNDThermalNoiseServiceInFreq: " << endl;
  
  out << prefix << "          LogLevel: " <<  fLogLevel << endl;
  out << prefix << "        RandomSeed: " <<  fRandomSeed << endl;
  out << prefix << "  NoiseArrayPoints: " << fNoiseArrayPoints << endl;
  
  return out;
}



//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(SBNDThermalNoiseServiceInFreq, ChannelNoiseService)

//**********************************************************************
