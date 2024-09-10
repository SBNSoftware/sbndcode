#include "sbndcode/TPC1DSignalProcessing/tools/IHitEfficiencyHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// Eigen
#include <Eigen/Dense>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

#include <cmath>
#include <algorithm>

namespace TrackHitEfficiencyAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       TrackHitEfficiencyAnalysis
    // Module Type: producer
    // File:        TrackHitEfficiencyAnalysis.cc
    //
    //              The intent of this module is to provide methods for
    //              "analyzing" hits on waveforms
    //
    // Configuration parameters:
    //
    // TruncMeanFraction     - the fraction of waveform bins to discard when
    //
    // Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
    //
    ////////////////////////////////////////////////////////////////////////
    
// The following typedefs will, obviously, be useful
using HitPtrVec       = std::vector<art::Ptr<recob::Hit>>;
using ViewHitMap      = std::map<size_t,HitPtrVec>;
using TrackViewHitMap = std::map<int,ViewHitMap>;

class TrackHitEfficiencyAnalysis : virtual public IHitEfficiencyHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TrackHitEfficiencyAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~TrackHitEfficiencyAnalysis();
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&) override;

    /**
     *  @brief Interface for initializing the tuple variables
     *
     *  @param TTree          pointer to a TTree object to which to add variables
     */
    void initializeTuple(TTree*) override;

    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const art::Event&)  const override;
    
private:
    
    // Clear mutable variables
    void clear() const;
    
    // Fcl parameters.
    std::vector<art::InputTag>  fRawDigitProducerLabelVec;
    std::vector<art::InputTag>  fWireProducerLabelVec;
    std::vector<art::InputTag>  fHitProducerLabelVec;
    art::InputTag               fMCParticleProducerLabel;
    art::InputTag               fSimChannelProducerLabel;
    art::InputTag               fBadChannelProducerLabel;
    bool                        fUseBadChannelDB;
    std::string                 fLocalDirName;           ///< Fraction for truncated mean
    std::vector<int>            fOffsetVec;              ///< Allow offsets for each plane
    std::vector<float>          fSigmaVec;               ///< Window size for matching to SimChannels
    int                         fMinAllowedChanStatus;   ///< Don't consider channels with lower status

    // Pointers to the histograms we'll create.
    std::vector<TH1F*>          fTotalElectronsHistVec;
    std::vector<TH1F*>          fMaxElectronsHistVec;
    std::vector<TH1F*>          fHitElectronsVec;
    std::vector<TH1F*>          fHitSumADCVec;
    std::vector<TH1F*>          fHitIntegralHistVec;
    std::vector<TH1F*>          fHitPulseHeightVec;
    std::vector<TH1F*>          fHitPulseWidthVec;
    std::vector<TH1F*>          fSimNumTDCVec;
    std::vector<TH1F*>          fHitNumTDCVec;
    std::vector<TH1F*>          fSnippetLenVec;
    std::vector<TH1F*>          fNMatchedHitVec;
    std::vector<TH1F*>          fDeltaMidTDCVec;
    std::vector<TProfile*>      fWireEfficVec;
    std::vector<TProfile*>      fWireEfficPHVec;
    std::vector<TProfile*>      fHitEfficVec;
    std::vector<TProfile*>      fHitEfficPHVec;
    std::vector<TProfile*>      fHitEfficXZVec;
    std::vector<TProfile*>      fHitEfficRMSVec;
    std::vector<TProfile*>      fCosXZvRMSVec;
    std::vector<TProfile2D*>    fHitENEvXZVec;
    std::vector<TH1F*>          fSimDivHitChgVec;
    std::vector<TH1F*>          fSimDivHitChg1Vec;
    std::vector<TH2F*>          fHitVsSimChgVec;
    std::vector<TH2F*>          fHitVsSimIntVec;
    std::vector<TH2F*>          fToteVHitEIntVec;

    std::vector<TH1F*>          fNSimChannelHitsVec;
    std::vector<TH1F*>          fNRecobHitVec;
 //std::vector<TH1F*>          fNRejectedHitVec;
    std::vector<TH1F*>          fHitEfficiencyVec;
    
 std::vector<TH1F*>          fNFakeHitVec;

    // TTree variables
    mutable TTree*             fTree;
    
    mutable std::vector<int>   fTPCVec;
    mutable std::vector<int>   fCryoVec;
    mutable std::vector<int>   fPlaneVec;
    mutable std::vector<int>   fWireVec;
    
    mutable std::vector<float> fTotalElectronsVec;
    mutable std::vector<float> fMaxElectronsVec;
    mutable std::vector<int>   fStartTickVec;
    mutable std::vector<int>   fStopTickVec;
    mutable std::vector<int>   fMaxETickVec;
    mutable std::vector<float> fPartDirX;
    mutable std::vector<float> fPartDirY;
    mutable std::vector<float> fPartDirZ;

    mutable std::vector<int>   fNMatchedWires;
    mutable std::vector<int>   fNMatchedHits;
    
    mutable std::vector<float> fHitPeakTimeVec;
    mutable std::vector<float> fHitPeakAmpVec;
    mutable std::vector<float> fHitPeakRMSVec;
    mutable std::vector<float> fHitBaselinevec;
    mutable std::vector<float> fHitSummedADCVec;
    mutable std::vector<float> fHitIntegralVec;
    mutable std::vector<int>   fHitStartTickVec;
    mutable std::vector<int>   fHitStopTickVec;
    mutable std::vector<int>   fHitMultiplicityVec;
    mutable std::vector<int>   fHitLocalIndexVec;
    mutable std::vector<float> fHitGoodnessVec;
    mutable std::vector<int>   fNumDegreesVec;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
TrackHitEfficiencyAnalysis::TrackHitEfficiencyAnalysis(fhicl::ParameterSet const & pset) : fTree(nullptr)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("TrackHitEfficiencyAnalysis") << "TrackHitEfficiencyAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
TrackHitEfficiencyAnalysis::~TrackHitEfficiencyAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TrackHitEfficiencyAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fRawDigitProducerLabelVec = pset.get< std::vector<art::InputTag>>("RawDigitLabelVec",   std::vector<art::InputTag>() = {"rawdigitfilter"});
    fWireProducerLabelVec     = pset.get< std::vector<art::InputTag>>("WireModuleLabelVec", std::vector<art::InputTag>() = {"decon1droi"});
    fHitProducerLabelVec      = pset.get< std::vector<art::InputTag>>("HitModuleLabelVec",  std::vector<art::InputTag>() = {"gauss"});
    fMCParticleProducerLabel  = pset.get< art::InputTag             >("MCParticleLabel",    "largeant");
    fSimChannelProducerLabel  = pset.get< art::InputTag             >("SimChannelLabel",    "largeant");
    fBadChannelProducerLabel  = pset.get< art::InputTag             >("BadChannelLabel",    "simnfspl1:badchannels");
    fUseBadChannelDB          = pset.get< bool                      >("UseBadChannelDB",    true);
    fLocalDirName             = pset.get<std::string                >("LocalDirName",       std::string("wow"));
    fOffsetVec                = pset.get<std::vector<int>           >("OffsetVec",          std::vector<int>()={0,0,0});
    fSigmaVec                 = pset.get<std::vector<float>         >("SigmaVec",           std::vector<float>()={1.,1.,1.});
    fMinAllowedChanStatus     = pset.get< int                       >("MinAllowedChannelStatus");
}

//----------------------------------------------------------------------------
/// Begin job method.
void TrackHitEfficiencyAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    fTotalElectronsHistVec.resize(fGeometry->Nplanes());
    fMaxElectronsHistVec.resize(fGeometry->Nplanes());
    fHitElectronsVec.resize(fGeometry->Nplanes());
    fHitSumADCVec.resize(fGeometry->Nplanes());
    fHitIntegralHistVec.resize(fGeometry->Nplanes());
    fHitPulseHeightVec.resize(fGeometry->Nplanes());
    fHitPulseWidthVec.resize(fGeometry->Nplanes());
    fSimNumTDCVec.resize(fGeometry->Nplanes());
    fHitNumTDCVec.resize(fGeometry->Nplanes());
    fSnippetLenVec.resize(fGeometry->Nplanes());
    fNMatchedHitVec.resize(fGeometry->Nplanes());
    fDeltaMidTDCVec.resize(fGeometry->Nplanes());
    fHitVsSimChgVec.resize(fGeometry->Nplanes());
    fHitVsSimIntVec.resize(fGeometry->Nplanes());
    fCosXZvRMSVec.resize(fGeometry->Nplanes());
    fToteVHitEIntVec.resize(fGeometry->Nplanes());
    fHitENEvXZVec.resize(fGeometry->Nplanes());
    fNSimChannelHitsVec.resize(fGeometry->Nplanes());
    fNRecobHitVec.resize(fGeometry->Nplanes());
    fNFakeHitVec.resize(fGeometry->Nplanes());
    fHitEfficiencyVec.resize(fGeometry->Nplanes());
    fSimDivHitChgVec.resize(fGeometry->Nplanes());
    fSimDivHitChg1Vec.resize(fGeometry->Nplanes());

    fWireEfficVec.resize(fGeometry->Nplanes());
    fWireEfficPHVec.resize(fGeometry->Nplanes());

    fHitEfficVec.resize(fGeometry->Nplanes());
    fHitEfficPHVec.resize(fGeometry->Nplanes());
    fHitEfficXZVec.resize(fGeometry->Nplanes());
    fHitEfficRMSVec.resize(fGeometry->Nplanes());

    for(size_t plane = 0; plane < fGeometry->Nplanes(); plane++)
    {
        fTotalElectronsHistVec.at(plane) = dir.make<TH1F>(("TotalElecs"  + std::to_string(plane)).c_str(), ";Total # electrons",     250,    0.,  100000.);
        fMaxElectronsHistVec.at(plane)   = dir.make<TH1F>(("MaxElecs"    + std::to_string(plane)).c_str(), ";Max # electrons",       250,    0.,  20000.);
        fHitElectronsVec.at(plane)       = dir.make<TH1F>(("HitElecs"    + std::to_string(plane)).c_str(), ";# e- in Hit Range",     250,    0.,  100000.);
        fHitSumADCVec.at(plane)          = dir.make<TH1F>(("SumADC"      + std::to_string(plane)).c_str(), "Hit Sum ADC",            200,    0.,  1000.);
        fHitIntegralHistVec.at(plane)    = dir.make<TH1F>(("Integral"    + std::to_string(plane)).c_str(), "Hit Integral",           200,    0.,  1000.);
        fHitPulseHeightVec.at(plane)     = dir.make<TH1F>(("PulseHeight" + std::to_string(plane)).c_str(), "Hit PH (ADC)",           200,    0.,  150.);
        fHitPulseWidthVec.at(plane)      = dir.make<TH1F>(("PulseWidth"  + std::to_string(plane)).c_str(), "Hit PulseWidth;Hit RMS", 100,    0.,  20.);
        fSimNumTDCVec.at(plane)          = dir.make<TH1F>(("SimNumTDC"   + std::to_string(plane)).c_str(), ";TDC ticks",             100,    0.,  100.);
        fHitNumTDCVec.at(plane)          = dir.make<TH1F>(("HitNumTDC"   + std::to_string(plane)).c_str(), ";ticks",                 100,    0.,  100.);
        fSnippetLenVec.at(plane)         = dir.make<TH1F>(("SnippetLen"  + std::to_string(plane)).c_str(), ";ticks",                 100,    0.,  100.);
        fNMatchedHitVec.at(plane)        = dir.make<TH1F>(("NMatched"    + std::to_string(plane)).c_str(), ";# hits",                 20,    0.,  20.);
        fDeltaMidTDCVec.at(plane)        = dir.make<TH1F>(("DeltaMid"    + std::to_string(plane)).c_str(), ";# hits",                 50,  -25.,  25.);
        fNSimChannelHitsVec.at(plane)    = dir.make<TH1F>(("NSimChan"    + std::to_string(plane)).c_str(), ";# hits",                100,    0.,  1200.);
        fNRecobHitVec.at(plane)          = dir.make<TH1F>(("NRecobHit"   + std::to_string(plane)).c_str(), ";# hits",                100,    0.,  1200.);
        fNFakeHitVec.at(plane)           = dir.make<TH1F>(("NFakeHit"    + std::to_string(plane)).c_str(), ";# hits",                100,    0.,  50.);
        fHitEfficiencyVec.at(plane)      = dir.make<TH1F>(("PlnEffic"    + std::to_string(plane)).c_str(), ";# hits",                101,    0.,  1.01);
        fSimDivHitChgVec.at(plane)       = dir.make<TH1F>(("SimDivHit"   + std::to_string(plane)).c_str(), ";# e / SummedADC",       200,    0.,  200.);
        fSimDivHitChg1Vec.at(plane)      = dir.make<TH1F>(("SimDivHit1"  + std::to_string(plane)).c_str(), ";# e / Integral",        200,    0.,  200.);

        fHitVsSimChgVec.at(plane)        = dir.make<TH2F>(("HitVSimQ"  + std::to_string(plane)).c_str(), "# e vs Hit SumADC;SumADC;# e",     50, 0.,   1000., 50, 0., 100000.);
        fHitVsSimIntVec.at(plane)        = dir.make<TH2F>(("HitVSimI"  + std::to_string(plane)).c_str(), "# e vs Hit Integral;Integral;# e", 50, 0.,   1000., 50, 0., 100000.);
        fToteVHitEIntVec.at(plane)       = dir.make<TH2F>(("ToteVHite" + std::to_string(plane)).c_str(), "Tot e vs Hit e;Total #e;Hit #e",   50, 0., 100000., 50, 0., 100000.);
        
        fWireEfficVec.at(plane)          = dir.make<TProfile>  (("WireEffic"   + std::to_string(plane)).c_str(), "Wire Efficiency;# electrons",      50, 0., 100000.,  0., 1.);
        fWireEfficPHVec.at(plane)        = dir.make<TProfile>  (("WireEfficPH" + std::to_string(plane)).c_str(), "Wire Efficiency;# electrons",      50, 0.,  20000.,  0., 1.);

        fHitEfficVec.at(plane)           = dir.make<TProfile>  (("HitEffic"    + std::to_string(plane)).c_str(), "Hit Efficiency;Total # e-",        50,  0., 100000., 0., 1.);
        fHitEfficPHVec.at(plane)         = dir.make<TProfile>  (("HitEfficPH"  + std::to_string(plane)).c_str(), "Hit Efficiency;Max # e-",          50,  0.,  20000., 0., 1.);
        fHitEfficXZVec.at(plane)         = dir.make<TProfile>  (("HitEfficXZ"  + std::to_string(plane)).c_str(), "Hit Efficiency;cos(thetaXZ)",      50, -1.,      1., 0., 1.);
        fHitEfficRMSVec.at(plane)        = dir.make<TProfile>  (("HitEfficRMS" + std::to_string(plane)).c_str(), "Hit Efficiency;MC Length (ticks)", 50,  0.,    100., 0., 1.);
        fCosXZvRMSVec.at(plane)          = dir.make<TProfile>  (("CosXZVRMS"   + std::to_string(plane)).c_str(), "CosXZ v. RMS;Cos(XZ);Ticks",       50, -1.,      1., 0., 100.);
        fHitENEvXZVec.at(plane)          = dir.make<TProfile2D>(("HitENEvXZ"   + std::to_string(plane)).c_str(), "XZ v # e;cos(thetaXZ);#electrons", 50, -1.,      1.,
                                                                50, 0., 100000., 0., 1.);
    }
    
    return;
}
    
void TrackHitEfficiencyAnalysis::initializeTuple(TTree* tree)
{
    fTree = tree;
    
    fTree->Branch("CryostataVec",      "std::vector<int>",   &fCryoVec);
    fTree->Branch("TPCVec",            "std::vector<int>",   &fTPCVec);
    fTree->Branch("PlaneVec",          "std::vector<int>",   &fPlaneVec);
    fTree->Branch("WireVec",           "std::vector<int>",   &fWireVec);

    fTree->Branch("TotalElectronsVec", "std::vector<float>", &fTotalElectronsVec);
    fTree->Branch("MaxElectronsVec",   "std::vector<float>", &fMaxElectronsVec);
    fTree->Branch("StartTick",         "std::vector<int>",   &fStartTickVec);
    fTree->Branch("StopTick",          "std::vector<int>",   &fStopTickVec);
    fTree->Branch("MaxETick",          "std::vector<int>",   &fMaxETickVec);
    fTree->Branch("PartDirX",          "std::vector<float>", &fPartDirX);
    fTree->Branch("PartDirY",          "std::vector<float>", &fPartDirY);
    fTree->Branch("PartDirZ",          "std::vector<float>", &fPartDirZ);

    fTree->Branch("NMatchedWires",     "std::vector<int>",   &fNMatchedWires);
    fTree->Branch("NMatchedHits",      "std::vector<int>",   &fNMatchedHits);
    
    fTree->Branch("HitPeakTimeVec",    "std::vector<float>", &fHitPeakTimeVec);
    fTree->Branch("HitPeakAmpVec",     "std::vector<float>", &fHitPeakAmpVec);
    fTree->Branch("HitPeakRMSVec",     "std::vector<float>", &fHitPeakRMSVec);
    fTree->Branch("HitBaselineVec",    "std::vector<float>", &fHitBaselinevec);
    fTree->Branch("HitSummedADCVec",   "std::vector<float>", &fHitSummedADCVec);
    fTree->Branch("HitIntegralVec",    "std::vector<float>", &fHitIntegralVec);
    fTree->Branch("HitStartTickVec",   "std::vector<int>",   &fHitStartTickVec);
    fTree->Branch("HitStopTickVec",    "std::vector<int>",   &fHitStopTickVec);
    fTree->Branch("HitMultiplicity",   "std::vector<int>",   &fHitMultiplicityVec);
    fTree->Branch("HitLocalIndex",     "std::vector<int>",   &fHitLocalIndexVec);
    fTree->Branch("HitGoodness",       "std::vector<float>", &fHitGoodnessVec);
    fTree->Branch("HitNumDegrees",     "std::vector<int>",   &fNumDegreesVec);

    clear();

    return;
}
    
void TrackHitEfficiencyAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();
    fWireVec.clear();

    fTotalElectronsVec.clear();
    fMaxElectronsVec.clear();
    fStartTickVec.clear();
    fStopTickVec.clear();
    fMaxETickVec.clear();
    fPartDirX.clear();
    fPartDirY.clear();
    fPartDirZ.clear();
    
    fNMatchedWires.clear();
    fNMatchedHits.clear();

    fHitPeakTimeVec.clear();
    fHitPeakAmpVec.clear();
    fHitPeakRMSVec.clear();
    fHitBaselinevec.clear();
    fHitSummedADCVec.clear();
    fHitIntegralVec.clear();
    fHitStartTickVec.clear();
    fHitStopTickVec.clear();
    fHitMultiplicityVec.clear();
    fHitLocalIndexVec.clear();
    fHitGoodnessVec.clear();
    fNumDegreesVec.clear();
    
    return;
}

void TrackHitEfficiencyAnalysis::fillHistograms(const art::Event& event) const
{
    std::cout << " filling histos " << std::endl;
    // Basic assumption is that the producer label vecs for RawDigits and Wire data are
    // all the same length and in the same order. Here we just check for length
    if (fRawDigitProducerLabelVec.size() != fWireProducerLabelVec.size()) return;
    std::cout<<fRawDigitProducerLabelVec.size()<<std::endl;
    // Always clear the tuple
    clear();
    
    art::Handle< std::vector<sim::SimChannel>> simChannelHandle;
    event.getByLabel(fSimChannelProducerLabel, simChannelHandle);
    
    art::Handle< std::vector<simb::MCParticle>> mcParticleHandle;
    event.getByLabel(fMCParticleProducerLabel, mcParticleHandle);

    // If there is no sim channel informaton then exit
    if (!simChannelHandle.isValid() || simChannelHandle->empty() || !mcParticleHandle.isValid()) 
      {
	std::cout<<"No sim channel information"<<std::endl;
	return;
      }
    std::cout<<"There is sim channel information"<<std::endl;
    
    
    // There are several things going on here... for each channel we have particles (track id's) depositing energy in a range to ticks
    // So... for each channel we want to build a structure that relates particles to tdc ranges and deposited energy (or electrons)
    // Here is a complicated structure:
    using TDCToIDEMap             = std::map<unsigned short, sim::IDE>; // We need this one in order
    using ChanToTDCToIDEMap       = std::map<raw::ChannelID_t, TDCToIDEMap>;
    using PartToChanToTDCToIDEMap = std::unordered_map<int, ChanToTDCToIDEMap>;
    
    PartToChanToTDCToIDEMap partToChanToTDCToIDEMap;
    
    // Build out the above data structure
    for(const auto& simChannel : *simChannelHandle)
    {
        for(const auto& tdcide : simChannel.TDCIDEMap())
        {
            for(const auto& ide : tdcide.second) partToChanToTDCToIDEMap[ide.trackID][simChannel.Channel()][tdcide.first] = ide;
        }
    }
    
    // what needs to be done?
    // First we define a straightforward channel to Wire map so we can look up a given
    // channel's Wire data as we loop over SimChannels.
    using ChanToWireMap = std::unordered_map<raw::ChannelID_t,const recob::Wire*>;
    
    ChanToWireMap channelToWireMap;
    
    // We will use the presence of a RawDigit as an indicator of a good channel... So
    // we want a mapping between channel and RawDigit
    using ChanToRawDigitMap = std::unordered_map<raw::ChannelID_t,const raw::RawDigit*>;
    
    ChanToRawDigitMap chanToRawDigitMap;
    
    // Look up the list of bad channels
    art::Handle< std::vector<int>> badChannelHandle;
    event.getByLabel(fBadChannelProducerLabel, badChannelHandle);

    //    std::cout<<"Pre TPC loop (L455ish)"<<std::endl;

    // Now start a loop over the individual TPCs to build out the structures for RawDigits and Wires
    for(size_t tpcID = 0; tpcID < fRawDigitProducerLabelVec.size(); tpcID++)
    {
      //std::cout<<"starting TPC loop"<<std::endl;
        art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
        event.getByLabel(fRawDigitProducerLabelVec[tpcID], rawDigitHandle);

	//std::cout << "rawDigitHandle: " << rawDigitHandle << std::endl;
	art::Handle< std::vector<recob::Wire> > wireHandle;;
        event.getByLabel(fWireProducerLabelVec[tpcID], wireHandle);

	//std::cout << "wireHandle" << wireHandle << std::endl;

        if (!rawDigitHandle.isValid() || !wireHandle.isValid()) 
	{
	  //	  std::cout<<"Inside TPC loop if statement (L470ish)"<<std::endl;
	  return;
        }
        for(const auto& wire : *wireHandle) channelToWireMap[wire.Channel()] = &wire;
        
        for(const auto& rawDigit : *rawDigitHandle) chanToRawDigitMap[rawDigit.Channel()] = &rawDigit;
    }
    //    std::cout<<"Post TPC loop (L470ish)"<<std::endl;

    // Now we create a data structure to relate hits to their channel ID
    using ChanToHitVecMap = std::unordered_map<raw::ChannelID_t,std::vector<const recob::Hit*>>;
    
    ChanToHitVecMap channelToHitVec;
    
    // And now fill it
    for(const auto& hitLabel : fHitProducerLabelVec)
    {
        art::Handle< std::vector<recob::Hit> > hitHandle;
        event.getByLabel(hitLabel, hitHandle);

        for(const auto& hit : *hitHandle) channelToHitVec[hit.Channel()].push_back(&hit);
    }
    
    // It is useful to create a mapping between trackID and MCParticle
    using TrackIDToMCParticleMap = std::unordered_map<int, const simb::MCParticle*>;
    
    TrackIDToMCParticleMap trackIDToMCParticleMap;
    
    for(const auto& mcParticle : *mcParticleHandle)
      trackIDToMCParticleMap[mcParticle.TrackId()] = &mcParticle;
    
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    
    std::vector<int> nSimChannelHitVec  = {0,0,0};
    std::vector<int> nRecobHitVec       = {0,0,0};
    std::vector<int> nFakeHitVec        = {0,0,0};
    std::vector<int> nSimulatedWiresVec = {0,0,0};
    
    unsigned int lastwire=-1;
    //std::cout<<"Just before eternal loop of loops (L502)"<<std::endl;
    
    for(const auto& partToChanInfo : partToChanToTDCToIDEMap)
    {
      //std::cout<<"Starting eternal loop of loops"<<std::endl;
      TrackIDToMCParticleMap::const_iterator trackIDToMCPartItr = trackIDToMCParticleMap.find(partToChanInfo.first);
        
        if (trackIDToMCPartItr == trackIDToMCParticleMap.end()) continue;
        
        int         trackPDGCode = trackIDToMCPartItr->second->PdgCode();
        std::string processName  = trackIDToMCPartItr->second->Process();

        // Looking for primary muons (e.g. CR Tracks)
        if (fabs(trackPDGCode) != 13 || processName != "primary") continue;

        // Recover particle position and angle information
        Eigen::Vector3f partStartPos(trackIDToMCPartItr->second->Vx(),trackIDToMCPartItr->second->Vy(),trackIDToMCPartItr->second->Vz());
        Eigen::Vector3f partStartDir(trackIDToMCPartItr->second->Px(),trackIDToMCPartItr->second->Py(),trackIDToMCPartItr->second->Pz());
        
        partStartDir.normalize();
        
        Eigen::Vector2f partStartDirVecXZ(partStartDir[0],partStartDir[2]);
        
        partStartDirVecXZ.normalize();
        
        // Assuming the SimChannels contain position information (currently not true for WC produced SimChannels)
        // then we want to keep a running position
        std::vector<Eigen::Vector3f> lastPositionVec = {partStartPos,partStartPos,partStartPos};

        for(const auto& chanToTDCToIDEMap : partToChanInfo.second)
        {
            // skip bad channels
            if (fUseBadChannelDB)
            {
                // This is the "correct" way to check and remove bad channels...
                if( chanFilt.Status(chanToTDCToIDEMap.first) < fMinAllowedChanStatus)
                {
                std::vector<geo::WireID> wids = fGeometry->ChannelToWire(chanToTDCToIDEMap.first);
                std::cout << "*** skipping bad channel with status: " << chanFilt.Status(chanToTDCToIDEMap.first) << " for channel: " << chanToTDCToIDEMap.first << ", plane: " << wids[0].Plane << ", wire: " << wids[0].Wire    << std::endl;
                    continue;
                }
            }
        
            // Was a list made available?
            // If so then we try that
            if (badChannelHandle.isValid())
            {
                // Here we query the input list from the wirecell processing
                std::vector<int>::const_iterator badItr = std::find(badChannelHandle->begin(),badChannelHandle->end(),chanToTDCToIDEMap.first);
    
                if (badItr != badChannelHandle->end()) continue;
                //            {
                //                ChanToRawDigitMap::const_iterator rawDigitItr = chanToRawDigitMap.find(chanToTDCToIDEMap.first);
                //
                //                if (rawDigitItr != chanToRawDigitMap.end())
                //                {
                //                    float nSig(3.);
                //                    float mean(0.);
                //                    float rmsFull(0.);
                //                    float rmsTrunc(0.);
                //                    int   nTrunc(0);
                //
                //                    getTruncatedMeanRMS(rawDigitItr->second->ADCs(), nSig, mean, rmsFull, rmsTrunc, nTrunc);
                //
                //                    std::cout << "--> Rejecting channel: " << chanToTDCToIDEMap.first << " from bad channel list, rms: " << rmsFull << std::endl;
                //                }
                //
                //                continue;
                //            }
            }
        
            TDCToIDEMap    tdcToIDEMap = chanToTDCToIDEMap.second;
            float          totalElectrons(0.);
            float          maxElectrons(0.);
            unsigned short maxElectronsTDC(0);
            int            nMatchedWires(0);
            int            nMatchedHits(0);
        
            // The below try-catch block may no longer be necessary
            // Decode the channel and make sure we have a valid one
            std::vector<geo::WireID> wids = fGeometry->ChannelToWire(chanToTDCToIDEMap.first);
        
            // Recover plane and wire in the plane
            unsigned int plane = wids[0].Plane;
            unsigned int wire  = wids[0].Wire;
            
            Eigen::Vector3f avePosition(0.,0.,0.);

            if(wire!=lastwire) nSimulatedWiresVec[plane]++;
            lastwire=wire;
    
            for(const auto& ideVal : tdcToIDEMap)
            {
                totalElectrons += ideVal.second.numElectrons;
        
                if (maxElectrons < ideVal.second.numElectrons)
                {
                    maxElectrons    = ideVal.second.numElectrons;
                    maxElectronsTDC = ideVal.first;
                }
        
                avePosition += Eigen::Vector3f(ideVal.second.x,ideVal.second.y,ideVal.second.z);
            }
        
            // Get local track direction by using the average position of deposited charge as the current position
            // and then subtracting the last position
            avePosition /= float(tdcToIDEMap.size());
        
            Eigen::Vector3f partDirVec = avePosition - lastPositionVec[plane];
        
            partDirVec.normalize();
        
            lastPositionVec[plane] = avePosition;
        
            // Now what we want is the projected angle in the XZ plane
            Eigen::Vector2f projPairDirVec(partDirVec[0],partDirVec[2]);
        
            projPairDirVec.normalize();
        
            float cosThetaXZ = projPairDirVec[0];
    
            totalElectrons = std::min(totalElectrons, float(99900.));
        
            fTotalElectronsHistVec[plane]->Fill(totalElectrons, 1.);
            fMaxElectronsHistVec[plane]->Fill(maxElectrons, 1.);
        
            nSimChannelHitVec[plane]++;
    
            unsigned short startTDC = tdcToIDEMap.begin()->first;
            unsigned short stopTDC  = tdcToIDEMap.rbegin()->first;

            // Convert to ticks to get in same units as hits
            unsigned short startTick = clockData.TPCTDC2Tick(startTDC)        + fOffsetVec[plane];
            unsigned short stopTick  = clockData.TPCTDC2Tick(stopTDC)         + fOffsetVec[plane];
            unsigned short maxETick  = clockData.TPCTDC2Tick(maxElectronsTDC) + fOffsetVec[plane];
    
            fSimNumTDCVec[plane]->Fill(stopTick - startTick, 1.);
        
            // Set up to extract the "best" parameters in the event of more than one hit for this pulse train
            float          nElectronsTotalBest(0.);
            float          hitSummedADCBest(0.);
            float          hitIntegralBest(0.);
            float          hitPeakTimeBest(0.);
            float          hitPeakAmpBest(-100.);
            float          hitRMSBest(0.);
            int            hitMultiplicityBest(0);
            int            hitLocalIndexBest(0);
            float          hitGoodnessBest(0.);
            int            hitNumDegreesBest(0);
            float          hitBaselineBest(0.);
            float          hitSnippetLenBest(0.);
            unsigned short hitStopTickBest(0);
            unsigned short hitStartTickBest(0);
        
            // Start by recovering the Wire associated to this channel
            ChanToWireMap::const_iterator wireItr = channelToWireMap.find(chanToTDCToIDEMap.first);
            
            if (wireItr != channelToWireMap.end())
            {
                const recob::Wire::RegionsOfInterest_t&       signalROI = wireItr->second->SignalROI();
                const lar::sparse_vector<float>::datarange_t* wireRangePtr(NULL);
                
                // Here we need to match the range of the ROI's on the given Wire with the tick range from the SimChannel
                for(const auto& range : signalROI.get_ranges())
                {
                    // #################################################
                    // ### Getting a vector of signals for this wire ###
                    // #################################################
                    //std::vector<float> signal(wire->Signal());
                    
                    const std::vector<float>& signal          = range.data();
                    raw::TDCtick_t            roiFirstBinTick = range.begin_index();
                    raw::TDCtick_t            roiLastBinTick  = roiFirstBinTick + signal.size();
                    
                    // If no overlap then go to next
                    if (roiFirstBinTick > stopTick || roiLastBinTick < startTick) continue;
                    
                    wireRangePtr = &range;
                    break;
                }
                
                // Check that we have found the wire range
                // Note that if we have not matched an ROI then we can't have a hit either so skip search for that...
                if (wireRangePtr)
                {
                    const recob::Hit* rejectedHit = 0;
                    const recob::Hit* bestHit     = 0;
                    
                    nMatchedWires++;
                    
                    // The next mission is to recover the hits associated to this Wire
                    // The easiest way to do this is to simply look up all the hits on this channel and then match
                    ChanToHitVecMap::iterator hitIter = channelToHitVec.find(chanToTDCToIDEMap.first);
                    
                    if (hitIter != channelToHitVec.end())
                    {
                        // Loop through the hits for this channel and look for matches
                        // In the event of more than one hit associated to the sim channel range, keep only
                        // the best match (assuming the nearby hits are "extra")
                        // Note that assumption breaks down for long pulse trains but worry about that later
                        for(const auto& hit : hitIter->second)
                        {
                            unsigned short hitStartTick = hit->PeakTime() - fSigmaVec[plane] * hit->RMS();
                            unsigned short hitStopTick  = hit->PeakTime() + fSigmaVec[plane] * hit->RMS();
                    
                            // If hit is out of range then skip, it is not related to this particle
                            if (hitStartTick > stopTick || hitStopTick < startTick)
                            {
                                 nFakeHitVec[plane]++;
                                 rejectedHit = hit;
                                continue;
                            }
                    
                            float hitHeight = hit->PeakAmplitude();
                    
                            // Use the hit with the largest pulse height as the "best"
                            if (hitHeight < hitPeakAmpBest) continue;
                    
                            hitPeakAmpBest   = hitHeight;
                            bestHit          = hit;
                            hitStartTickBest = hitStartTick;
                            hitStopTickBest  = hitStopTick;
                        }
                    
                        // Find a match?
                        if (bestHit)
                        {
			  //			    std::cout << "found a bestHit (L720)" << std::endl;
			    nElectronsTotalBest = 0.;
                            hitPeakTimeBest     = bestHit->PeakTime();
                            hitIntegralBest     = bestHit->Integral();
                            hitSummedADCBest    = bestHit->SummedADC();
                            hitRMSBest          = bestHit->RMS();
                            hitMultiplicityBest = bestHit->Multiplicity();
                            hitLocalIndexBest   = bestHit->LocalIndex();
                            hitGoodnessBest     = bestHit->GoodnessOfFit();
                            hitNumDegreesBest   = bestHit->DegreesOfFreedom();
                            hitSnippetLenBest   = bestHit->EndTick() - bestHit->StartTick();
                            hitBaselineBest     = 0.;  // To do...
                    
                            nMatchedHits++;
                    
                            // Get the number of electrons
                            for(unsigned short tick = hitStartTickBest; tick <= hitStopTickBest; tick++)
                            {
                                unsigned short hitTDC = clockData.TPCTick2TDC(tick - fOffsetVec[plane]);
                    
                                TDCToIDEMap::iterator ideIterator = tdcToIDEMap.find(hitTDC);
                    
                                if (ideIterator != tdcToIDEMap.end()) nElectronsTotalBest += ideIterator->second.numElectrons;
                            }
                        }
        
                        if (nMatchedHits > 0)
                        {
                            float chgRatioADC = hitSummedADCBest > 0. ? totalElectrons / hitSummedADCBest : 0.;
                            float chgRatioInt = hitIntegralBest  > 0. ? totalElectrons / hitIntegralBest  : 0.;
        
                            fHitSumADCVec[plane]->Fill(hitSummedADCBest, 1.);
                            fHitIntegralHistVec[plane]->Fill(hitIntegralBest, 1.);
                            fSimDivHitChgVec[plane]->Fill(chgRatioADC, 1.);
                            fSimDivHitChg1Vec[plane]->Fill(chgRatioInt, 1.);
                            fHitVsSimChgVec[plane]->Fill(std::min(hitSummedADCBest,float(999.)), totalElectrons, 1.);
                            fHitVsSimIntVec[plane]->Fill(std::min(hitIntegralBest,float(999.)), totalElectrons, 1.);
                            fToteVHitEIntVec[plane]->Fill(std::min(totalElectrons,float(99999.)),std::min(nElectronsTotalBest,float(99999.)),1.);
                            fHitPulseHeightVec[plane]->Fill(std::min(hitPeakAmpBest,float(149.5)), 1.);
                            fHitPulseWidthVec[plane]->Fill(std::min(hitRMSBest,float(19.8)), 1.);
                            fHitElectronsVec[plane]->Fill(nElectronsTotalBest, 1.);
                            fHitNumTDCVec[plane]->Fill(std::min(float(hitStopTickBest - hitStartTickBest),float(99.5)), 1.);
                            fSnippetLenVec[plane]->Fill(std::min(hitSnippetLenBest, float(99.5)), 1.);
                            fDeltaMidTDCVec[plane]->Fill(hitPeakTimeBest - maxETick, 1.);
                    
                            nRecobHitVec[plane]++;

                        }
                        else if (rejectedHit)
                        { 
                            unsigned short hitStartTick = rejectedHit->PeakTime() - fSigmaVec[plane] * rejectedHit->RMS();
                            unsigned short hitStopTick  = rejectedHit->PeakTime() + fSigmaVec[plane] * rejectedHit->RMS();
        
                            mf::LogDebug("TrackHitEfficiencyAnalysis") << "**> TPC: " << rejectedHit->WireID().TPC << ", Plane " << rejectedHit->WireID().Plane << ", wire: " << rejectedHit->WireID().Wire << ", hit startstop            tick: " << hitStartTick << "/" << hitStopTick << ", start/stop ticks: " << startTick << "/" << stopTick << std::endl;
                            mf::LogDebug("TrackHitEfficiencyAnalysis") << "    TPC/Plane/Wire: " << wids[0].TPC << "/" << plane << "/" << wids[0].Wire << ", Track # hits: " << partToChanInfo.second.size() << ", # hits: "<<         hitIter->second.size() << ", # electrons: " << totalElectrons << ", pulse Height: " << rejectedHit->PeakAmplitude() << ", charge: " << rejectedHit->Integral()      << ", " <<rejectedHit->SummedADC() << std::endl;
                        }
                        else
                        {
                            mf::LogDebug("TrackHitEfficiencyAnalysis") << "==> No match, TPC/Plane/Wire: " << "/" << wids[0].TPC << "/" << wids[0].Plane << "/" << wids[0].Wire << ", # electrons: " << totalElectrons << ",           startTick: " << startTick << ", stopTick: " << stopTick << std::endl;
                        }
                    }
                }
            }

            fWireEfficVec.at(plane)->Fill(totalElectrons, std::min(nMatchedWires,1), 1.);
            fWireEfficPHVec.at(plane)->Fill(maxElectrons, std::min(nMatchedWires,1), 1.);
            
            float matchHit   = std::min(nMatchedHits,1);
            float snippetLen = std::min(float(stopTick - startTick),float(99.5));
            
            fNMatchedHitVec[plane]->Fill(nMatchedHits, 1.);
            fHitEfficVec[plane]->Fill(totalElectrons,   matchHit,   1.);
            fHitEfficPHVec[plane]->Fill(maxElectrons,   matchHit,   1.);
            fHitEfficXZVec[plane]->Fill(cosThetaXZ,     matchHit,   1.);
            fCosXZvRMSVec[plane]->Fill(cosThetaXZ,      snippetLen, 1.);
            fHitEfficRMSVec[plane]->Fill(snippetLen,    matchHit,   1.);
            
            fHitENEvXZVec[plane]->Fill(cosThetaXZ, totalElectrons, matchHit);
            
            // Store tuple variables
            fTPCVec.push_back(wids[0].TPC);
            fCryoVec.push_back(wids[0].Cryostat);
            fPlaneVec.push_back(wids[0].Plane);
            fWireVec.push_back(wids[0].Wire);
            
            fTotalElectronsVec.push_back(totalElectrons);
            fMaxElectronsVec.push_back(maxElectrons);
            fStartTickVec.push_back(startTick);
            fStopTickVec.push_back(stopTick);
            fMaxETickVec.push_back(maxETick);
            fPartDirX.push_back(partDirVec[0]);
            fPartDirY.push_back(partDirVec[1]);
            fPartDirZ.push_back(partDirVec[2]);
            
            fNMatchedWires.push_back(nMatchedWires);
            fNMatchedHits.push_back(nMatchedHits);
            
            fHitPeakTimeVec.push_back(hitPeakTimeBest);
            fHitPeakAmpVec.push_back(hitPeakAmpBest);
            fHitPeakRMSVec.push_back(hitRMSBest);
            fHitBaselinevec.push_back(hitBaselineBest);
            fHitSummedADCVec.push_back(hitSummedADCBest);
            fHitIntegralVec.push_back(hitIntegralBest);
            fHitStartTickVec.push_back(hitStartTickBest);
            fHitStopTickVec.push_back(hitStopTickBest);
            fHitMultiplicityVec.push_back(hitMultiplicityBest);
            fHitLocalIndexVec.push_back(hitLocalIndexBest);
            fHitGoodnessVec.push_back(hitGoodnessBest);
            fNumDegreesVec.push_back(hitNumDegreesBest);
        }
    }

    for(size_t idx = 0; idx < fGeometry->Nplanes();idx++)
    {
        if (nSimChannelHitVec[idx] > 10)
        {
            float hitEfficiency = float(nRecobHitVec[idx]) / float(nSimChannelHitVec[idx]);
        
            fNSimChannelHitsVec[idx]->Fill(std::min(nSimChannelHitVec[idx],1999),1.);
            fNRecobHitVec[idx]->Fill(std::min(nRecobHitVec[idx],1999), 1.);
            fNFakeHitVec[idx]->Fill(nFakeHitVec[idx]/(float)nSimulatedWiresVec[idx],1.);
            fHitEfficiencyVec[idx]->Fill(hitEfficiency, 1.);
        }
    }

    return;
}
    
// Useful for normalizing histograms
void TrackHitEfficiencyAnalysis::endJob(int numEvents)
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(TrackHitEfficiencyAnalysis)
}
