
////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceSBND_service.cc
/// \author H. Greenlee   (adapted to LArIAT by A. Szelc)
////////////////////////////////////////////////////////////////////////

#include "sbndcode/Utilities/SignalShapingServiceSBND.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/LArFFT.h"
#include "TFile.h"

//----------------------------------------------------------------------
// Constructor.
util::SignalShapingServiceSBND::SignalShapingServiceSBND(const fhicl::ParameterSet& pset,
								    art::ActivityRegistry& /* reg */) 
  : fInit(false)
{
  reconfigure(pset);
}


//----------------------------------------------------------------------
// Destructor.
util::SignalShapingServiceSBND::~SignalShapingServiceSBND()
{}


//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShapingServiceSBND::reconfigure(const fhicl::ParameterSet& pset)
{
  // Reset initialization flag.

  fInit = false;

  // Reset kernels.

  fColSignalShaping.Reset();
  fIndUSignalShaping.Reset();
  fIndVSignalShaping.Reset();

  // Fetch fcl parameters.

  fDeconNorm = pset.get<double>("DeconNorm");
  fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
  fASICGainInMVPerFC = pset.get<std::vector<double> >("ASICGainInMVPerFC");

  fNFieldBins = pset.get<int>("FieldBins");
  fCol3DCorrection = pset.get<double>("Col3DCorrection");
  fInd3DCorrection = pset.get<double>("Ind3DCorrection");
  fColFieldRespAmp = pset.get<double>("ColFieldRespAmp");
  fIndUFieldRespAmp = pset.get<double>("IndUFieldRespAmp");
  fIndVFieldRespAmp = pset.get<double>("IndVFieldRespAmp");
  fShapeTimeConst = pset.get<std::vector<double> >("ShapeTimeConst");
  fNoiseFactVec = pset.get<std::vector<DoubleVec> >("NoiseFactVec");

  fInputFieldRespSamplingPeriod = pset.get<double>("InputFieldRespSamplingPeriod"); 
  fFieldResponseTOffset = pset.get<std::vector<double> >("FieldResponseTOffset");  

  fUseFunctionFieldShape= pset.get<bool>("UseFunctionFieldShape");
  fUseSimpleFieldShape = pset.get<bool>("UseSimpleFieldShape");
  fUseHistogramFieldShape = pset.get<bool>("UseHistogramFieldShape");

  if(fUseSimpleFieldShape) {
    fNFieldBins = 300;
  }
  fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");
  
  // Construct parameterized collection filter function.
  if(!fGetFilterFromHisto) {

    mf::LogInfo("SignalShapingServiceSBND") << "Getting Filter from .fcl file" ;
    std::string colFilt = pset.get<std::string>("ColFilter");
    std::vector<double> colFiltParams = pset.get<std::vector<double> >("ColFilterParams");
    fColFilterFunc = new TF1("colFilter", colFilt.c_str());
    for(unsigned int i=0; i<colFiltParams.size(); ++i)
      fColFilterFunc->SetParameter(i, colFiltParams[i]);
    
    // Construct parameterized induction filter function.

    std::string indUFilt = pset.get<std::string>("IndUFilter");
    std::vector<double> indUFiltParams = pset.get<std::vector<double> >("IndUFilterParams");
    fIndUFilterFunc = new TF1("indUFilter", indUFilt.c_str());
    for(unsigned int i=0; i<indUFiltParams.size(); ++i)
      fIndUFilterFunc->SetParameter(i, indUFiltParams[i]);

    std::string indVFilt = pset.get<std::string>("IndVFilter");
    std::vector<double> indVFiltParams = pset.get<std::vector<double> >("IndVFilterParams");
    fIndVFilterFunc = new TF1("indVFilter", indVFilt.c_str());
    for(unsigned int i=0; i<indVFiltParams.size(); ++i)
      fIndVFilterFunc->SetParameter(i, indVFiltParams[i]);
  } else {
    constexpr unsigned int NPlanes = 3;
    
    std::string histoname = pset.get<std::string>("FilterHistoName");
    mf::LogInfo("SignalShapingServiceSBND") << " using filter from .root file " ;
   
    // constructor decides if initialized value is a path or an environment variable
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(pset.get<std::string>("FilterFunctionFname"), fname);
    
    TFile in(fname.c_str(), "READ");
    if (!in.IsOpen()) {
      throw cet::exception("SignalShapingServiceSBND")
        << "Can't open filter function file '" << fname << "'!\n";
    }
    mf::LogInfo("SignalShapingServiceSBND")
      << "Reading filter histograms from '" << fname << "'";
    for(unsigned int i = 0; i < NPlanes; ++i) {
      auto pHist = dynamic_cast<TH1*>(in.Get(Form(histoname.c_str(),i)));
      if (!pHist) {
        // this happens also if there is an object but it's not a TH1
        throw cet::exception("SignalShapingServiceSBND")
          << "Can't find filter histogram '" << histoname << "' for plane #" << i
          << " in '" << fname << "'!\n";
      }
      pHist->SetDirectory(nullptr); // detach the histogram from its source file
      fFilterHist[i] = pHist;
    }
    in.Close();
  }
 
  /////////////////////////////////////
  if(fUseFunctionFieldShape) {

    std::string colField = pset.get<std::string>("ColFieldShape");
    std::vector<double> colFieldParams = pset.get<std::vector<double> >("ColFieldParams");
    fColFieldFunc = new TF1("colField", colField.c_str());
    for(unsigned int i=0; i<colFieldParams.size(); ++i)
      fColFieldFunc->SetParameter(i, colFieldParams[i]);

    // Construct parameterized induction filter function.
    
    std::string indUField = pset.get<std::string>("IndUFieldShape");
    std::vector<double> indUFieldParams = pset.get<std::vector<double> >("IndUFieldParams");
    fIndUFieldFunc = new TF1("indUField", indUField.c_str());
    for(unsigned int i=0; i<indUFieldParams.size(); ++i)
      fIndUFieldFunc->SetParameter(i, indUFieldParams[i]);
    // Warning, last parameter needs to be multiplied by the FFTSize, in current version of the code,
    
    std::string indVField = pset.get<std::string>("IndVFieldShape");
    std::vector<double> indVFieldParams = pset.get<std::vector<double> >("IndVFieldParams");
    fIndVFieldFunc = new TF1("indVField", indVField.c_str());
    for(unsigned int i=0; i<indVFieldParams.size(); ++i)
      fIndVFieldFunc->SetParameter(i, indVFieldParams[i]);
    // Warning, last parameter needs to be multiplied by the FFTSize, in current version of the code,

  } else if (fUseHistogramFieldShape){
    constexpr unsigned int NPlanes = 3;

    //constructor decides if initialized value is a path or an environment variable
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file( pset.get<std::string>("FieldResponseFname"), fname);
    std::string histoname = pset.get<std::string>("FieldResponseHistoName");
    
    mf::LogInfo("SignalShapingServiceSBND")
      << "Using the field response provided from '" << fname
      << "' (histograms '" << histoname << "_*')";

    TFile fin(fname.c_str(), "READ");
    if ( !fin.IsOpen() ) {
      throw cet::exception("SignalShapingServiceSBND")
        << "Could not find the field response file '" << fname << "'!\n";
    }

    const std::string iPlane[3] = {"U", "V", "Y"};

    for(unsigned int i = 0; i < NPlanes; ++i) {
      std::string PlaneHistoName = histoname + "_" + iPlane[i];
      LOG_DEBUG("SignalShapingServiceSBND")
        << "Field Response " << i << ": " << PlaneHistoName;
      
      auto pHist = dynamic_cast<TH1*>(fin.Get(PlaneHistoName.c_str()));
      if (!pHist) {
        throw cet::exception("SignalShapingServiceSBND")
          << "Could not find the field response histogram '" << PlaneHistoName
          << "' in file '" << fname << "'\n";
      } 
      if (pHist->GetNbinsX() > fNFieldBins) {
        throw art::Exception( art::errors::Configuration ) << "FieldBins (" << fNFieldBins
          << ") should always be larger than or equal to the number of the bins in the input histogram ("
          << pHist->GetNbinsX() << " in '" << PlaneHistoName << "')!\n";
      }
      pHist->SetDirectory(nullptr); // detach the histogram from his source file
      
      fFieldResponseHist[i] = pHist;
      LOG_DEBUG("SignalShapingServiceSBND")
        << "RESPONSE HISTOGRAM " << iPlane[i] << ": " << pHist->GetEntries() << " entries in "
        << pHist->GetNbinsX() << " bins (" << pHist->GetBinLowEdge(1)
        << " to " << pHist->GetBinLowEdge(pHist->GetNbinsX() + 1);
    }

    fin.Close();
  }

}


//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceSBND::SignalShaping(unsigned int channel) const
{
  if(!fInit)
    init();

  // Figure out plane type.

  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  // we need to distiguish the U and V planes
  geo::View_t view = geom->View(channel);

  // Return appropriate shaper.
  //geo::SigType_t sigtype = geom->SignalType(channel);

  if (view == geo::kU)
    return fIndUSignalShaping;
  else if (view == geo::kV)
    return fIndVSignalShaping;
  else if (view == geo::kZ)
    return fColSignalShaping;
  else
    throw cet::exception("SignalShapingServiceSBND")<< "1 can't determine"
                                                          << " SignalType\n";  
  

  return fColSignalShaping;
}

//---Give Gain Settings to SimWire ---//
double util::SignalShapingServiceSBND::GetASICGain(unsigned int const channel) const
{
  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  // we need to distiguish the U and V planes
  geo::View_t view = geom->View(channel);

  double gain = 0.0;
  if(view == geo::kU)
    gain = fASICGainInMVPerFC.at(0);
  else if(view == geo::kV)
     gain = fASICGainInMVPerFC.at(1);
  else if(view == geo::kZ)
    gain = fASICGainInMVPerFC.at(2);
  else
    throw cet::exception("SignalShapingServiceSBND")<< "2 can't determine"
                                                    << " SignalType\n";
  return gain;
} 

// //---Give Shaping time Settings to SimWire ---//
// double util::SignalShapingServiceSBND::GetShapingTime(unsigned int const channel) const
// {
//   art::ServiceHandle<geo::Geometry> geom;
//   //geo::SigType_t sigtype = geom->SignalType(channel);

//   // we need to distiguish the U and V planes
//   geo::View_t view = geom->View(channel);

//   double shaping_time = 0;
//   if(view == geo::kU)
//     shaping_time = fShapeTimeConst.at(0);
//   if(view == geo::kV)
//     shaping_time = fShapeTimeConst.at(1);
//   else if(view == geo::kZ)
//     shaping_time = fShapeTimeConst.at(2);
//   else
//     throw cet::exception("SignalShapingServiceSBND")<< "3 can't determine"
//                                                     << " SignalType\n";
//   return shaping_time;
// } 

double util::SignalShapingServiceSBND::GetRawNoise(unsigned int const channel) const
{
  unsigned int plane;
  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  // we need to distiguish the U and V planes
  geo::View_t view = geom->View(channel);

  if(view == geo::kU)
    plane = 0;
  else if(view == geo::kV)
    plane = 1;
  else if(view == geo::kZ)
    plane = 2;
  else
    throw cet::exception("SignalShapingServiceSBND")<< "4 can't determine"
                                                    << " SignalType\n";

  double shapingtime = fShapeTimeConst.at(plane);
  double gain = fASICGainInMVPerFC.at(plane);
  int temp;
  if (shapingtime == 0.5){
    temp = 0;
  }else if (shapingtime == 1.0){
    temp = 1;
  }else if (shapingtime == 2.0){
    temp = 2;
  }else{
    temp = 3;
  }
  double rawNoise;
  
  auto tempNoise = fNoiseFactVec.at(plane);
  rawNoise = tempNoise.at(temp);

  rawNoise *= gain/4.7;
  return rawNoise;
}

double util::SignalShapingServiceSBND::GetDeconNoise(unsigned int const channel) const
{
  unsigned int plane;
  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  // we need to distiguish the U and V planes
  geo::View_t view = geom->View(channel);

  if(view == geo::kU)
    plane = 0;
  else if(view == geo::kV)
    plane = 1;
  else if(view == geo::kZ)
    plane = 2;
  else
    throw cet::exception("SignalShapingServiceSBND")<< "5 can't determine"
                                                    << " SignalType\n";

  double shapingtime = fShapeTimeConst.at(plane);
  int temp;
  if (shapingtime == 0.5){
    temp = 0;
  }else if (shapingtime == 1.0){
    temp = 1;
  }else if (shapingtime == 2.0){
    temp = 2;
  }else{
    temp = 3;
  }
  double deconNoise;
  
  auto tempNoise = fNoiseFactVec.at(plane);
  deconNoise = tempNoise.at(temp);

  // replaced 2000 with fADCPerPCAtLowestASICGain/4.7 because 2000 V/ADC is specific to MicroBooNE
  deconNoise = deconNoise /4096.*(fADCPerPCAtLowestASICGain/4.7/4.7) *6.241*1000/fDeconNorm;
  return deconNoise;
}



//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceSBND::init()
{
  if(!fInit) {
    fInit = true;

    // Do microboone-specific configuration of SignalShaping by providing
    // microboone response and filter functions.

    // Calculate field and electronics response functions.

    SetFieldResponse();
    SetElectResponse(fShapeTimeConst.at(2),fASICGainInMVPerFC.at(2));

    // Configure convolution kernels.

    fColSignalShaping.AddResponseFunction(fColFieldResponse);
    fColSignalShaping.AddResponseFunction(fElectResponse);
    fColSignalShaping.save_response();
    fColSignalShaping.set_normflag(false);
    //fColSignalShaping.SetPeakResponseTime(0.);

    SetElectResponse(fShapeTimeConst.at(0),fASICGainInMVPerFC.at(0));

    fIndUSignalShaping.AddResponseFunction(fIndUFieldResponse);
    fIndUSignalShaping.AddResponseFunction(fElectResponse);
    fIndUSignalShaping.save_response();
    fIndUSignalShaping.set_normflag(false);
    //fIndUSignalShaping.SetPeakResponseTime(0.);

    SetElectResponse(fShapeTimeConst.at(1),fASICGainInMVPerFC.at(1));

    fIndVSignalShaping.AddResponseFunction(fIndVFieldResponse);
    fIndVSignalShaping.AddResponseFunction(fElectResponse);
    fIndVSignalShaping.save_response();
    fIndVSignalShaping.set_normflag(false);
    //fIndVSignalShaping.SetPeakResponseTime(0.);

    SetResponseSampling();

    // Calculate filter functions.

    SetFilters();

    // Configure deconvolution kernels.

    fColSignalShaping.AddFilterFunction(fColFilter);
    fColSignalShaping.CalculateDeconvKernel();

    fIndUSignalShaping.AddFilterFunction(fIndUFilter);
    fIndUSignalShaping.CalculateDeconvKernel();

    fIndVSignalShaping.AddFilterFunction(fIndVFilter);
    fIndVSignalShaping.CalculateDeconvKernel();
  }
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceSBND::SetFieldResponse()
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  auto const* detprop = lar::providerFrom<::detinfo::DetectorPropertiesService>();

  // Get plane pitch.
 
  double xyz1[3] = {0.};
  double xyz2[3] = {0.};
  double xyzl[3] = {0.};
  // should always have at least 2 planes
  geo->Plane(0).LocalToWorld(xyzl, xyz1);
  geo->Plane(1).LocalToWorld(xyzl, xyz2);

  // this assumes all planes are equidistant from each other,
  // probably not a bad assumption
  double pitch = xyz2[0] - xyz1[0]; ///in cm

  fColFieldResponse.resize(fNFieldBins, 0.);
  fIndUFieldResponse.resize(fNFieldBins, 0.);
  fIndVFieldResponse.resize(fNFieldBins, 0.);

  // set the response for the collection plane first
  // the first entry is 0

  double driftvelocity=detprop->DriftVelocity()/1000.;  
  double integral = 0.;
  ////////////////////////////////////////////////////
  if(fUseFunctionFieldShape) {

    art::ServiceHandle<util::LArFFT> fft;
    int signalSize = fft->FFTSize();
    std::vector<double> ramp(signalSize);
    // TComplex kernBin;
    // int size = signalSize/2;
    // int bin=0;
    //std::vector<TComplex> freqSig(size+1);
    std::vector<double> bipolar(signalSize);    
    
    fColFieldResponse.resize(signalSize, 0.);
    fIndUFieldResponse.resize(signalSize, 0.);
    fIndVFieldResponse.resize(signalSize, 0.);
   
    // Hardcoding. Bad. Temporary hopefully.
    fIndUFieldFunc->SetParameter(4,fIndUFieldFunc->GetParameter(4)*signalSize);
    fIndVFieldFunc->SetParameter(4,fIndVFieldFunc->GetParameter(4)*signalSize);

    for(int i = 0; i < signalSize; i++) {
      ramp[i]=fColFieldFunc->Eval(i);
      fColFieldResponse[i]=ramp[i];
      integral += fColFieldResponse[i];
      // rampc->Fill(i,ramp[i]);
      bipolar[i]=fIndUFieldFunc->Eval(i);
      fIndUFieldResponse[i]=bipolar[i];
      bipolar[i]=fIndVFieldFunc->Eval(i);
      fIndVFieldResponse[i]=bipolar[i];
      // bipol->Fill(i,bipolar[i]);
    }
     
    for(int i = 0; i < signalSize; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
    }
      
    //this might be not necessary if the function definition is not defined in the middle of the signal range  
    fft->ShiftData(fIndUFieldResponse,signalSize/2.0);
    fft->ShiftData(fIndVFieldResponse,signalSize/2.0);

  }else if (fUseHistogramFieldShape){ 
    //Ticks in nanoseconds
    //Calculate the normalization of the collection plane

    for(int ibin=1; ibin<=fFieldResponseHist[2]->GetNbinsX(); ibin++)
      integral += fFieldResponseHist[2]->GetBinContent(ibin);

    //Induction U Plane
    for(int ibin=1; ibin<=fFieldResponseHist[0]->GetNbinsX(); ibin++)
      fIndUFieldResponse[ibin-1] = fIndUFieldRespAmp*fFieldResponseHist[0]->GetBinContent(ibin)/integral;
    
    //Induction V Plane
    for(int ibin=1; ibin<=fFieldResponseHist[1]->GetNbinsX(); ibin++)
      fIndVFieldResponse[ibin-1] = fIndVFieldRespAmp*fFieldResponseHist[1]->GetBinContent(ibin)/integral;

    //Collection Plane
    for(int ibin=1; ibin<=fFieldResponseHist[2]->GetNbinsX(); ibin++)
      fColFieldResponse[ibin-1] = fColFieldRespAmp*fFieldResponseHist[2]->GetBinContent(ibin)/integral;

  } else if (fUseSimpleFieldShape) {
   
    mf::LogInfo("SignalShapingServiceSBND") << " using try-2 hard-coded field shapes " ;

    const int nbincPlane = 16;
    double cPlaneResponse[nbincPlane] = {
      0,               0,               0,   0.02620087336,   0.02620087336, 
      0.04366812227,    0.1310043668,    0.1659388646,    0.1397379913,    0.3711790393, 
      0.06550218341,    0.0480349345,  -0.01310043668, -0.004366812227,               0, 
      0
    };

    for(int i = 1; i < nbincPlane; ++i){
      fColFieldResponse[i] = cPlaneResponse[i];
      integral += fColFieldResponse[i];
    }

    for(int i = 0; i < nbincPlane; ++i){
      //fColFieldResponse[i] *= fColFieldRespAmp/integral;
      fColFieldResponse[i] /= integral;
    }

    //const int nbiniOld = 6;
    const int nbinuPlane = 228;
    // now induction plane 0 ("U")
    // this response function has a very long (first) positive lobe, ~ 100 usec
    // So for starters, we us the single-lobe filter

    double uPlaneResponse[nbinuPlane] = {
      0, 0.0001881008778, 0.0003762017556, 0.0005643026334, 0.0007524035112, 
      0.000940504389,  0.001128605267,  0.001316706145,  0.001504807022,    0.0016929079, 
      0.001881008778,  0.002069109656,  0.002257210534,  0.002445311411,  0.002633412289, 
      0.002821513167,  0.003009614045,  0.003197714923,    0.0033858158,  0.003573916678, 
      0.003762017556,  0.003950118434,  0.004138219312,  0.004326320189,  0.004514421067, 
      0.004702521945,  0.004890622823,  0.005078723701,  0.005266824579,  0.005454925456, 
      0.005643026334,  0.005831127212,   0.00601922809,  0.006207328968,  0.006395429845, 
      0.006583530723,  0.006771631601,  0.006959732479,  0.007147833357,  0.007335934234, 
      0.007524035112,   0.00771213599,  0.007900236868,  0.008088337746,  0.008276438623, 
      0.008464539501,  0.008652640379,  0.008840741257,  0.009028842135,  0.009216943012, 
      0.00940504389,  0.009593144768,  0.009781245646,  0.009969346524,    0.0101574474, 
      0.01053554828,   0.01053364916,   0.01072175003,   0.01090985091,   0.01109795179, 
      0.01128605267,   0.01147415355,   0.01166225442,    0.0118503553,   0.01203845618, 
      0.01222655706,   0.01241465794,   0.01260275881,   0.01279085969,   0.01297896057, 
      0.01316706145,   0.01335516232,    0.0135432632,   0.01373136408,   0.01391946496, 
      0.01410756584,   0.01429566671,   0.01448376759,   0.01467186847,   0.01485996935, 
      0.01504807022,    0.0152361711,   0.01542427198,   0.01561237286,   0.01580047374, 
      0.01598857461,   0.01617667549,   0.01636477637,   0.01655287725,   0.01674097812, 
      0.016929079,   0.01711717988,   0.01730528076,   0.01749338164,   0.01768148251, 
      0.01786958339,   0.01805768427,   0.01824578515,   0.01843388602,    0.0186219869, 
      0.01881008778,   0.01899818866,   0.01918628954,   0.01937439041,   0.01956249129, 
      0.01975059217,   0.01993869305,   0.02012679393,    0.0203148948,   0.02050299568, 
      0.02069109656,   0.02087919744,   0.02106729831,   0.02125539919,   0.02144350007, 
      0.02163160095,   0.02181970183,    0.0220078027,   0.02219590358,   0.02238400446, 
      0.02257210534,   0.02302354744,   0.02347498955,   0.02392643166,   0.02437787376, 
      0.02482931587,   0.02528075798,   0.02573220008,   0.02618364219,    0.0266350843, 
      0.0270865264,   0.02753796851,   0.02798941062,   0.02844085272,   0.02889229483, 
      0.02934373694,   0.02979517904,   0.03024662115,   0.03069806326,   0.03114950536, 
      0.03160094747,    0.0321652501,   0.03272955274,   0.03329385537,     0.033858158, 
      0.03442246064,   0.03498676327,    0.0355510659,   0.03611536854,   0.03667967117, 
      0.03724397381,   0.03780827644,   0.03837257907,   0.03893688171,   0.03950118434, 
      0.04006548697,   0.04062978961,   0.04119409224,   0.04175839487,   0.04232269751, 
      0.04288700014,    0.0435641633,   0.04424132646,   0.04491848962,   0.04559565278, 
      0.04627281594,    0.0469499791,   0.04762714226,   0.04830430542,   0.04898146858, 
      0.04965863174,    0.0503357949,   0.05101295806,   0.05169012122,   0.05236728438, 
      0.05304444754,    0.0537216107,   0.05439877386,   0.05507593702,   0.05575310018, 
      0.05643026334,    0.0572579072,   0.05808555107,   0.05891319493,   0.05974083879, 
      0.06056848265,   0.06139612652,   0.06222377038,   0.06305141424,    0.0638790581, 
      0.06470670196,   0.06553434583,   0.06636198969,   0.06718963355,   0.06801727741, 
      0.06884492128,   0.06967256514,     0.070500209,   0.07132785286,   0.07215549673, 
      0.07298314059,   0.07305838094,   0.07313362129,   0.07320886164,   0.07328410199, 
      0.07335934234,   0.07524035112,   0.07524035112,   0.07524035112,   0.07524035112, 
      0.07524035112,   0.07524035112,   0.07524035112,   0.07565835307,   0.06688031211, 
      -1.508982036,    -1.401197605,    -1.293413174,   -0.5748502994,   -0.3233532934, 
      -0.2694610778,   -0.2694610778,   -0.1796407186,   -0.1437125749,  -0.03592814371, 
      0,               0,               0
    };

    for(int i = 0; i < nbinuPlane; ++i){
      //fIndUFieldResponse[i] = fIndUFieldRespAmp*uPlaneResponse[i]/(nbiniOld);
      fIndUFieldResponse[i] = uPlaneResponse[i]/integral;
    }
   
    const int nbinvPlane = 20;
    double vPlaneResponse[nbinvPlane] = {
      0,               0,   0.01090909091,   0.01090909091,   0.01090909091, 
      0.02181818182,   0.03272727273,    0.7636363636,     2.018181818,            2.04, 
      1.090909091,     -1.03861518,    -1.757656458,    -1.757656458,    -1.118508655, 
      -0.2396804261,  -0.07989347537, -0.007989347537,               0,               0
    };

    for (int i = 0; i < nbinvPlane; ++i) {
      //fIndVFieldResponse[i] = vPlaneResponse[i]*fIndVFieldRespAmp/(nbiniOld);
      fIndVFieldResponse[i] = vPlaneResponse[i]/integral;
    }

   } else {

    //////////////////////////////////////////////////
    mf::LogInfo("SignalShapingServiceSBND") << " using the old field shape " ;
    int nbinc = TMath::Nint(fCol3DCorrection*(std::abs(pitch))/(driftvelocity*detprop->SamplingRate())); ///number of bins //KP
    
    double integral = 0.;
    for(int i = 1; i < nbinc; ++i){
      fColFieldResponse[i] = fColFieldResponse[i-1] + 1.0;
      integral += fColFieldResponse[i];
    }

    for(int i = 0; i < nbinc; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
    }

    // now the induction plane
    
    int nbini = TMath::Nint(fInd3DCorrection*(std::abs(pitch))/(driftvelocity*detprop->SamplingRate()));//KP
    for(int i = 0; i < nbini; ++i){
      fIndUFieldResponse[i] = fIndUFieldRespAmp/(1.*nbini);
      fIndUFieldResponse[nbini+i] = -fIndUFieldRespAmp/(1.*nbini);
    }

    for(int i = 0; i < nbini; ++i){
      fIndVFieldResponse[i] = fIndVFieldRespAmp/(1.*nbini);
      fIndVFieldResponse[nbini+i] = -fIndVFieldRespAmp/(1.*nbini);
    }

  }
  
  return;
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceSBND::SetElectResponse(double shapingtime, double gain)
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  //auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingSBND") << "Setting SBND electronics response function...";

  int nticks = fft->FFTSize();
  fElectResponse.resize(nticks, 0.);
  std::vector<double> time(nticks,0.);

  //Gain and shaping time variables from fcl file:    
  double Ao = 1.0;  //gain
  double To = shapingtime;  //peaking time
    
  // this is actually sampling time, in ns
  // mf::LogInfo("SignalShapingSBND") << "Check sampling intervals: " 
  //                                  << fSampleRate << " ns" 
  //                                  << "Check number of samples: " << fNTicks;

  // The following sets the microboone electronics response function in 
  // time-space. Function comes from BNL SPICE simulation of SBND 
  // electronics. SPICE gives the electronics transfer function in 
  // frequency-space. The inverse laplace transform of that function 
  // (in time-space) was calculated in Mathematica and is what is being 
  // used below. Parameters Ao and To are cumulative gain/timing parameters 
  // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain. 
  // They have been adjusted to make the SPICE simulation to match the 
  // actual electronics response. Default params are Ao=1.4, To=0.5us. 
  double max=0.;
  
  for(size_t i = 0; i < fElectResponse.size(); ++i){

    //convert time to microseconds, to match fElectResponse[i] definition
    time[i] = (1.*i)*fInputFieldRespSamplingPeriod*1e-3; 
    fElectResponse[i] = 
      4.31054*exp(-2.94809*time[i]/To)*Ao - 2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*Ao
      -2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*cos(2.38722*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*cos(5.18561*time[i]/To)*Ao
      +0.762456*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*Ao
      -0.762456*exp(-2.82833*time[i]/To)*cos(2.38722*time[i]/To)*sin(1.19361*time[i]/To)*Ao
      +0.762456*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao
      -2.6202*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao 
      -0.327684*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*Ao + 
      +0.327684*exp(-2.40318*time[i]/To)*cos(5.18561*time[i]/To)*sin(2.5928*time[i]/To)*Ao
      -0.327684*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao;

      if(fElectResponse[i] > max) max = fElectResponse[i];
  }// end loop over time buckets
    
  LOG_DEBUG("SignalShapingSBND") << " Done.";

  // normalize fElectResponse[i], before the convolution
  // Put in overall normalization in a pedantic way:
  // first put in the pulse area per eleectron at the lowest gain setting,
  // then normalize by the actual ASIC gain setting used.
  // This code is executed only during initialization of service,
  // so don't worry about code inefficiencies here.
  for(auto& element : fElectResponse){
     element /= max;
     element *= fADCPerPCAtLowestASICGain*1.60217657e-7;
     element *= gain/4.7;
   }
   
  return;

}


//----------------------------------------------------------------------
// Calculate microboone filter functions.
void util::SignalShapingServiceSBND::SetFilters()
{
  // Get services.

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<util::LArFFT> fft;

  double ts = detprop->SamplingRate();
  int n = fft->FFTSize() / 2;

  // Calculate collection filter.

  fColFilter.resize(n+1);
  fIndUFilter.resize(n+1);
  fIndVFilter.resize(n+1);
  
  if(!fGetFilterFromHisto)
  {
  fColFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fColFilterFunc->Eval(freq);
    fColFilter[i] = TComplex(f, 0.);
  }
  

  // Calculate induction filters.
 
  fIndUFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fIndUFilterFunc->Eval(freq);
    fIndUFilter[i] = TComplex(f, 0.);
    }

   fIndVFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fIndVFilterFunc->Eval(freq);
    fIndVFilter[i] = TComplex(f, 0.);
    }
  
  }
  else
  {
    
    for(int i=0; i<=n; ++i) {
      double f = fFilterHist[2]->GetBinContent(i);  // hardcoded plane numbers. Bad. To change later.
      fColFilter[i] = TComplex(f, 0.);
      double g = fFilterHist[1]->GetBinContent(i);
      fIndVFilter[i] = TComplex(g, 0.);
      double h = fFilterHist[0]->GetBinContent(i);
      fIndUFilter[i] = TComplex(h, 0.);
      
    }
  }
  
  fIndUSignalShaping.AddFilterFunction(fIndUFilter);
  fIndVSignalShaping.AddFilterFunction(fIndVFilter);
  fColSignalShaping.AddFilterFunction(fColFilter);
  
}

//---------------------------------------------------------------
// Sample SBND response (the convoluted field and electronics
// response)
void util::SignalShapingServiceSBND::SetResponseSampling()
{
  //Get services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::LArFFT> fft;
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  //Operations permitted only if output of rebinning has a larger bin size
  if( fInputFieldRespSamplingPeriod > detprop->SamplingRate() )
    throw cet::exception("SignalShapingServiceSBND") << "\033[93m"
                                     << "Invalid operation: cannot rebin to a more finely binned vector!"    
                                     << "\033[00m" << std::endl;
  
  int nticks = fft->FFTSize();
  std::vector<double> SamplingTime(nticks,0.);
  for(int itime = 0; itime < nticks; itime++) {
    SamplingTime[itime] = (1.*itime) * detprop->SamplingRate();
  }

  //Sampling
  for(int iplane = 0; iplane <= 2; iplane++) {
    const std::vector<double>* pResp;
    switch( iplane ) {
    case 0: pResp = &(fIndUSignalShaping.Response_save()); break;
    case 1: pResp = &(fIndVSignalShaping.Response_save()); break;
    default: pResp = &(fColSignalShaping.Response_save()); break;          
    }
  
  std::vector<double> SamplingResp(nticks, 0.);
  
  int nticks_input = pResp->size();
  std::vector<double> InputTime(nticks_input, 0.);
  for(int itime = 0; itime < nticks_input; itime++) {
    InputTime[itime] = (1.*itime) * fInputFieldRespSamplingPeriod;
  }

  /*
    Much more sophisticated approach using a linear (trapezoidal) interpolation
    current deafult!
  */
  int SamplingCount = 0;
  for(int itime = 0; itime < nticks; itime++) {
    int low = -1, up = -1;
    for(int jtime = 0; jtime < nticks; jtime++) {
      if(InputTime[jtime] == SamplingTime[itime]) {
	SamplingResp[itime] = (*pResp)[jtime];
	SamplingCount++;
	break;
      } else if(InputTime[jtime] > SamplingTime[itime]) {
	low = jtime - 1;
	up = jtime;
	SamplingResp[itime] = (*pResp)[low] + (SamplingTime[itime] - InputTime[low]) * ( (*pResp)[up] - (*pResp)[low]) / (InputTime[up] - InputTime[low] );
	SamplingCount++;
	break;
      } else {
	SamplingResp[itime] = 0;
      }
    }// for(int jtime = 0; jtime < nticks; jtime++)
    
  }// for(int itime = 0; itime < nticks; itime++)

  SamplingResp.resize(SamplingCount, 0.);

  switch( iplane ) {
  case 0: fIndUSignalShaping.AddResponseFunction(SamplingResp, true); break;
  case 1: fIndVSignalShaping.AddResponseFunction(SamplingResp, true);; break;
  default: fColSignalShaping.AddResponseFunction(SamplingResp, true);; break;
  }

  }// for(int iplane = 0; iplane < 2; iplane++)

  return;
}

int util::SignalShapingServiceSBND::FieldResponseTOffset(unsigned int const channel) const
{
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  double time_offset = 0;
  if(view == geo::kU)
    time_offset = fFieldResponseTOffset.at(0);
  else if(view == geo::kV)
    time_offset = fFieldResponseTOffset.at(1);
  else if(view == geo::kZ)
    time_offset = fFieldResponseTOffset.at(2);
  else
    throw cet::exception("SignalShapingServiceSBND")<< "6 can't determine"
                                                    << " SignalType\n";
// std::cout << "TIME OFFSET" << 	time_offset << " " << view << std::endl;

  //auto tpc_clock = art::ServiceHandle<util::TimeService>()->TPCCLOCK();
  auto tpc_clock = lar::providerFrom<detinfo::DetectorClocksService>()->TPCClock();
  return tpc_clock.Ticks(time_offset/1.e3);
}


namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceSBND)

}
