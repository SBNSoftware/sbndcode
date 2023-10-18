#include <functional>

//ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"

//LArSoft includes
// #include "larcore/Geometry/Geometry.h"
// #include "larcorealg/Geometry/LocalTransformationGeo.h"
// #include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SupernovaTruth.h"
// #include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// #include "lardata/DetectorInfoServices/DetectorClocksService.h"
// #include "larsim/MCCheater/BackTrackerService.h"
// #include "larsim/MCCheater/ParticleInventoryService.h"
// #include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//ART includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

enum PType {kUnknown=0, kGen};
std::map<PType, std::string> PTypeString {{kUnknown, "Unknown"}, {kGen    , "Gen"    }};

class FullGeoAna : public art::EDAnalyzer {

public:
  explicit FullGeoAna(fhicl::ParameterSet const & p);

  FullGeoAna(FullGeoAna const &) = delete;
  FullGeoAna(FullGeoAna      &&) = delete;
  FullGeoAna & operator = (FullGeoAna const &) = delete;
  FullGeoAna & operator = (FullGeoAna      &&) = delete;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob();

private:

  // --- Private member variables ---
  void  ResetVariables();

  // --- LArSoft FHiCL inputs ---
  std::string fname              ;
  std::string fRawDigitLabel     ;
  std::string fHitLabel          ;
  std::string fCalDataModuleLabel;
  std::string fOpHitModuleLabel  ;

  std::string fGEANTLabel;
  std::string fGenLabel  ; std::map<int, simb::MCParticle> GenParts;

  std::map<int, const simb::MCParticle*> truthmap;
  std::map<int, PType> trkIDToPType;

  int firstCatch ;
  int secondCatch;
  int thirdCatch ;

  bool fSaveIDEs ;
  bool fSaveTruth;
  bool fSaveTPC  ;
  bool fSavePDS  ;

  TTree* fFullGeoAnaTree;

  int Run   ;
  int SubRun;
  int Event ;

  int NTotHit ;
  int NColHit ;
  int NIndHit ;
  int NHitNoBT;

  // --- Truth information ---
  std::vector<int>         True_Bck_PDG			;
  std::vector<int>         True_Bck_Mother		;
  std::vector<int>         True_Bck_EndProcess		;
  std::vector<int>         True_Bck_ID			;
  std::vector<double>      True_Bck_VertX		;
  std::vector<double>      True_Bck_VertY		;
  std::vector<double>      True_Bck_VertZ		;
  std::vector<double>      True_Bck_Time		;
  std::vector<double>      True_Bck_Energy		;
  std::vector<double>      True_Bck_EndE		;
  std::vector<double>      True_Bck_EndX		;
  std::vector<double>      True_Bck_EndY		;
  std::vector<double>      True_Bck_EndZ		;
  std::vector<double>      True_Bck_EndT		;
  std::vector<std::string> True_Bck_StartMaterial	;
  std::vector<std::string> True_Bck_EndMaterial		;

  // --- PDS information ---
  std::vector<int>                  PDS_OpHit_OpChannel      ;
  std::vector<double>               PDS_OpHit_X              ;
  std::vector<double>               PDS_OpHit_Y              ;
  std::vector<double>               PDS_OpHit_Z              ;
  std::vector<double>               PDS_OpHit_PeakTimeAbs    ;
  std::vector<double>               PDS_OpHit_PeakTime       ;
  std::vector<unsigned short>       PDS_OpHit_Frame          ;
  std::vector<double>               PDS_OpHit_Width          ;
  std::vector<double>               PDS_OpHit_Area           ;
  std::vector<double>               PDS_OpHit_Amplitude      ;
  std::vector<double>               PDS_OpHit_PE             ;
  std::vector<double>               PDS_OpHit_FastToTotal    ;
  std::vector<int>                  PDS_OpHit_True_GenType   ;
  std::vector<int>                  PDS_OpHit_True_Index     ;
  std::vector<double>               PDS_OpHit_True_Energy    ;
  std::vector<int>                  PDS_OpHit_True_TrackID   ;
  std::vector<int>                  PDS_OpHit_True_GenTypeAll;
  std::vector<double>               PDS_OpHit_True_EnergyAll ;
  std::vector<int>                  PDS_OpHit_True_TrackIDAll;
  std::vector<int>                  PDS_OpHit_True_IndexAll  ;

  // --- TPC information ---
  std::vector<int>                  Hit_View                 ;
  std::vector<int>                  Hit_Size                 ;
  std::vector<int>                  Hit_TPC                  ;
  std::vector<int>                  Hit_Chan                 ;
  std::vector<double>               Hit_X_start              ;
  std::vector<double>               Hit_Y_start              ;
  std::vector<double>               Hit_Z_start              ;
  std::vector<double>               Hit_X_end                ;
  std::vector<double>               Hit_Y_end                ;
  std::vector<double>               Hit_Z_end                ;
  std::vector<float>                Hit_Time                 ;
  std::vector<float>                Hit_RMS                  ;
  std::vector<float>                Hit_SADC                 ;
  std::vector<float>                Hit_Int                  ;
  std::vector<float>                Hit_Peak                 ;
  std::vector<int>                  Hit_True_GenType         ;
  std::vector<int>                  Hit_True_MainTrID        ;
  std::vector<int>                  Hit_True_TrackID         ;
  std::vector<float>                Hit_True_EvEnergy        ;
  std::vector<int>                  Hit_True_MarleyIndex     ;
  std::vector<float>                Hit_True_X               ;
  std::vector<float>                Hit_True_Y               ;
  std::vector<float>                Hit_True_Z               ;
  std::vector<float>                Hit_True_Energy          ;
  std::vector<float>                Hit_True_nElec           ;
  std::vector<int>                  Hit_True_nIDEs           ;

  // --- IDE information ---
  int   NTotIDEs;
  std::vector<int>                  IDEChannel               ;
  std::vector<int>                  IDEStartTime             ;
  std::vector<int>                  IDEEndTime               ;
  std::vector<float>                IDEEnergy                ;
  std::vector<float>                IDEElectrons             ;
  std::vector<int>                  IDEParticle              ;

  int TotGen_Gen;

  // art::ServiceHandle<geo::Geometry> geo;
  // art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  // art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

};

FullGeoAna::FullGeoAna(fhicl::ParameterSet const & p) 
  : EDAnalyzer(p), fname("FullGeoAna_module") {
  this->reconfigure(p);
}

void FullGeoAna::reconfigure(fhicl::ParameterSet const & p) {

  fRawDigitLabel      = p.get<std::string>("RawDigitLabel"     );
  fHitLabel           = p.get<std::string>("HitLabel"          );
  fCalDataModuleLabel = p.get<std::string>("CalDataModuleLabel");
  fOpHitModuleLabel   = p.get<std::string>("OpHitModuleLabel"  );

  fGEANTLabel = p.get<std::string>("GEANT4Label"   );
  fGenLabel   = p.get<std::string>("GeneratorLabel");

  fSaveIDEs  = p.get<bool>("SaveIDEs" , 0);
  fSaveTruth = p.get<bool>("SaveTruth", 0);
  fSaveTPC   = p.get<bool>("SaveTPC"  , 1);
  fSavePDS   = p.get<bool>("SavePDS"  , 0);

  mf::LogInfo(fname) << "Reconfigured " << this->processName()  << " with:\n"
		     << "SaveIDEs:  "   << std::boolalpha       << fSaveIDEs  << "\n"
    		     << "SaveTruth: "   << std::boolalpha       << fSaveTruth << "\n"
		     << "SaveTPC:   "   << std::boolalpha       << fSaveTPC   << "\n"
		     << "SavePDS:   "   << std::boolalpha       << fSavePDS   << std::endl;
}

void FullGeoAna::ResetVariables() {

  trkIDToPType.clear();
  GenParts    .clear();

  Run = SubRun = Event = -1;
  TotGen_Gen           =  0;

  NTotHit    = 0;
  NColHit    = 0;
  NIndHit    = 0;
  NHitNoBT   = 0;

  // --- Truth information ---
  True_Bck_PDG			.clear();
  True_Bck_Mother		.clear();
  True_Bck_EndProcess		.clear();
  True_Bck_ID			.clear();
  True_Bck_VertX		.clear();
  True_Bck_VertY		.clear();
  True_Bck_VertZ		.clear();
  True_Bck_Time			.clear();
  True_Bck_Energy		.clear();
  True_Bck_EndE			.clear();
  True_Bck_EndX			.clear();
  True_Bck_EndY			.clear();
  True_Bck_EndZ			.clear();
  True_Bck_EndT			.clear();
  True_Bck_StartMaterial	.clear();
  True_Bck_EndMaterial		.clear();

  // --- TPC information ---
  Hit_View			.clear();
  Hit_Size			.clear();
  Hit_TPC			.clear();
  Hit_Chan			.clear();
  Hit_X_start			.clear();
  Hit_Y_start			.clear();
  Hit_Z_start			.clear();
  Hit_X_end			.clear();
  Hit_Y_end			.clear();
  Hit_Z_end			.clear();
  Hit_Time			.clear();
  Hit_RMS			.clear();
  Hit_SADC			.clear();
  Hit_Int			.clear();
  Hit_Peak			.clear();
  Hit_True_GenType		.clear();
  Hit_True_MainTrID		.clear();
  Hit_True_TrackID		.clear();
  Hit_True_EvEnergy		.clear();
  Hit_True_MarleyIndex		.clear();
  Hit_True_X			.clear();
  Hit_True_Y			.clear();
  Hit_True_Z			.clear();
  Hit_True_Energy		.clear();
  Hit_True_nElec		.clear();
  Hit_True_nIDEs		.clear();

  // --- PDS information ---
  PDS_OpHit_OpChannel		.clear();
  PDS_OpHit_X			.clear();
  PDS_OpHit_Y			.clear();
  PDS_OpHit_Z			.clear();
  PDS_OpHit_PeakTimeAbs		.clear();
  PDS_OpHit_PeakTime		.clear();
  PDS_OpHit_Frame		.clear();
  PDS_OpHit_Width		.clear();
  PDS_OpHit_Area		.clear();
  PDS_OpHit_Amplitude		.clear();
  PDS_OpHit_PE			.clear();
  PDS_OpHit_FastToTotal		.clear();
  PDS_OpHit_True_GenType	.clear();
  PDS_OpHit_True_Index		.clear();
  PDS_OpHit_True_Energy		.clear();
  PDS_OpHit_True_TrackID	.clear();
  PDS_OpHit_True_GenTypeAll	.clear();
  PDS_OpHit_True_EnergyAll	.clear();
  PDS_OpHit_True_TrackIDAll	.clear();
  PDS_OpHit_True_IndexAll	.clear();

  // --- IDE information ---
  NTotIDEs=0;
  IDEChannel			.clear();
  IDEStartTime			.clear();
  IDEEndTime			.clear();
  IDEEnergy			.clear();
  IDEElectrons			.clear();
  IDEParticle			.clear();

}

void FullGeoAna::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  fFullGeoAnaTree = tfs->make<TTree>("FullGeoAna", "Analysis TTree for Full Geometry simulations");

  fFullGeoAnaTree->Branch("Run"     , &Run     , "Run/I"     );
  fFullGeoAnaTree->Branch("SubRun"  , &SubRun  , "SubRun/I"  );
  fFullGeoAnaTree->Branch("Event"   , &Event   , "Event/I"   );
  fFullGeoAnaTree->Branch("NTotHit" , &NTotHit , "NTotHits/I");
  fFullGeoAnaTree->Branch("NColHit" , &NColHit , "NColHits/I");
  fFullGeoAnaTree->Branch("NIndHit" , &NIndHit , "NIndHits/I");
  fFullGeoAnaTree->Branch("NHitNoBT", &NHitNoBT, "NHitNoBT/I");

  if (fSaveTruth) {
    fFullGeoAnaTree->Branch("True_Bck_PDG"          , &True_Bck_PDG          );
    fFullGeoAnaTree->Branch("True_Bck_Mother"       , &True_Bck_Mother       );
    fFullGeoAnaTree->Branch("True_Bck_EndProcess"   , &True_Bck_EndProcess   );
    fFullGeoAnaTree->Branch("True_Bck_ID"	    , &True_Bck_ID           );
    fFullGeoAnaTree->Branch("True_Bck_VertX"        , &True_Bck_VertX        );
    fFullGeoAnaTree->Branch("True_Bck_VertY"        , &True_Bck_VertY        );
    fFullGeoAnaTree->Branch("True_Bck_VertZ"        , &True_Bck_VertZ        );
    fFullGeoAnaTree->Branch("True_Bck_Time"         , &True_Bck_Time         );
    fFullGeoAnaTree->Branch("True_Bck_Energy"       , &True_Bck_Energy       );
    fFullGeoAnaTree->Branch("True_Bck_EndE"         , &True_Bck_EndE         );
    fFullGeoAnaTree->Branch("True_Bck_EndX"         , &True_Bck_EndX         );
    fFullGeoAnaTree->Branch("True_Bck_EndY"         , &True_Bck_EndY         );
    fFullGeoAnaTree->Branch("True_Bck_EndZ"         , &True_Bck_EndZ         );
    fFullGeoAnaTree->Branch("True_Bck_EndT"         , &True_Bck_EndT         );
    fFullGeoAnaTree->Branch("True_Bck_StartMaterial", &True_Bck_StartMaterial);
    fFullGeoAnaTree->Branch("True_Bck_EndMaterial"  , &True_Bck_EndMaterial  );
  }
 
  if (fSaveTPC) {
    fFullGeoAnaTree->Branch("Hit_View"                 , &Hit_View                 );
    fFullGeoAnaTree->Branch("Hit_Size"                 , &Hit_Size                 );
    fFullGeoAnaTree->Branch("Hit_TPC"                  , &Hit_TPC                  );
    fFullGeoAnaTree->Branch("Hit_Chan"                 , &Hit_Chan                 );
    fFullGeoAnaTree->Branch("Hit_X_start"              , &Hit_X_start              );
    fFullGeoAnaTree->Branch("Hit_Y_start"              , &Hit_Y_start              );
    fFullGeoAnaTree->Branch("Hit_Z_start"              , &Hit_Z_start              );
    fFullGeoAnaTree->Branch("Hit_X_end"                , &Hit_X_end                );
    fFullGeoAnaTree->Branch("Hit_Y_end"                , &Hit_Y_end                );
    fFullGeoAnaTree->Branch("Hit_Z_end"                , &Hit_Z_end                );
    fFullGeoAnaTree->Branch("Hit_Time"                 , &Hit_Time                 );
    fFullGeoAnaTree->Branch("Hit_RMS"                  , &Hit_RMS                  );
    fFullGeoAnaTree->Branch("Hit_SADC"                 , &Hit_SADC                 );
    fFullGeoAnaTree->Branch("Hit_Int"                  , &Hit_Int                  );
    fFullGeoAnaTree->Branch("Hit_Peak"                 , &Hit_Peak                 );
    fFullGeoAnaTree->Branch("Hit_True_GenType"         , &Hit_True_GenType         );
    fFullGeoAnaTree->Branch("Hit_True_MainTrID"        , &Hit_True_MainTrID        );
    fFullGeoAnaTree->Branch("Hit_True_TrackID"         , &Hit_True_TrackID         );
    fFullGeoAnaTree->Branch("Hit_True_EvEnergy"        , &Hit_True_EvEnergy        );
    fFullGeoAnaTree->Branch("Hit_True_MarleyIndex"     , &Hit_True_MarleyIndex     );
    fFullGeoAnaTree->Branch("Hit_True_X"               , &Hit_True_X               );
    fFullGeoAnaTree->Branch("Hit_True_Y"               , &Hit_True_Y               );
    fFullGeoAnaTree->Branch("Hit_True_Z"               , &Hit_True_Z               );
    fFullGeoAnaTree->Branch("Hit_True_Energy"          , &Hit_True_Energy          );
    fFullGeoAnaTree->Branch("Hit_True_nElec"           , &Hit_True_nElec           );
    fFullGeoAnaTree->Branch("Hit_True_nIDEs"           , &Hit_True_nIDEs           );
  }

  if (fSavePDS) {
    fFullGeoAnaTree->Branch("PDS_OpHit_OpChannel"      , &PDS_OpHit_OpChannel      );
    fFullGeoAnaTree->Branch("PDS_OpHit_X"              , &PDS_OpHit_X              );
    fFullGeoAnaTree->Branch("PDS_OpHit_Y"              , &PDS_OpHit_Y              );
    fFullGeoAnaTree->Branch("PDS_OpHit_Z"              , &PDS_OpHit_Z              );
    fFullGeoAnaTree->Branch("PDS_OpHit_PeakTimeAbs"    , &PDS_OpHit_PeakTimeAbs    );
    fFullGeoAnaTree->Branch("PDS_OpHit_PeakTime"       , &PDS_OpHit_PeakTime       );
    fFullGeoAnaTree->Branch("PDS_OpHit_Frame"          , &PDS_OpHit_Frame          );
    fFullGeoAnaTree->Branch("PDS_OpHit_Width"          , &PDS_OpHit_Width          );
    fFullGeoAnaTree->Branch("PDS_OpHit_Area"           , &PDS_OpHit_Area           );
    fFullGeoAnaTree->Branch("PDS_OpHit_Amplitude"      , &PDS_OpHit_Amplitude      );
    fFullGeoAnaTree->Branch("PDS_OpHit_PE"             , &PDS_OpHit_PE             );
    fFullGeoAnaTree->Branch("PDS_OpHit_FastToTotal"    , &PDS_OpHit_FastToTotal    );
    fFullGeoAnaTree->Branch("PDS_OpHit_True_GenType"   , &PDS_OpHit_True_GenType   );
    fFullGeoAnaTree->Branch("PDS_OpHit_True_Energy"    , &PDS_OpHit_True_Energy    );
    fFullGeoAnaTree->Branch("PDS_OpHit_True_TrackID"   , &PDS_OpHit_True_TrackID   );
    fFullGeoAnaTree->Branch("PDS_OpHit_True_GenTypeAll", &PDS_OpHit_True_GenTypeAll);
    fFullGeoAnaTree->Branch("PDS_OpHit_True_EnergyAll" , &PDS_OpHit_True_EnergyAll );
    fFullGeoAnaTree->Branch("PDS_OpHit_True_TrackIDAll", &PDS_OpHit_True_TrackIDAll);
    fFullGeoAnaTree->Branch("PDS_OpHit_True_IndexAll"  , &PDS_OpHit_True_IndexAll  );
  }

  if (fSaveIDEs) {
    fFullGeoAnaTree->Branch("NTotIDEs"                 , &NTotIDEs  , "NTotIDEs/I" );
    fFullGeoAnaTree->Branch("IDEChannel"               , &IDEChannel               );
    fFullGeoAnaTree->Branch("IDEStartTime"             , &IDEStartTime             );
    fFullGeoAnaTree->Branch("IDEEndTime"               , &IDEEndTime               );
    fFullGeoAnaTree->Branch("IDEEnergy"                , &IDEEnergy                );
    fFullGeoAnaTree->Branch("IDEElectrons"             , &IDEElectrons             );
    fFullGeoAnaTree->Branch("IDEParticle"              , &IDEParticle              );
  }

  fFullGeoAnaTree->Branch("TotGen_Gen", &TotGen_Gen, "TotGen_Gen/I");

}

void FullGeoAna::analyze(art::Event const & evt) {

  ResetVariables();

  Run    = evt.run   ();
  SubRun = evt.subRun();
  Event  = evt.event ();

  /*
    |    _____                   _______         _   _      
    |   / ____|                 |__   __|       | | | |     
    |  | (___   __ ___   _____     | |_ __ _   _| |_| |__   
    |   \___ \ / _` \ \ / / _ \    | | '__| | | | __| '_ \  
    |   ____) | (_| |\ V /  __/    | | |  | |_| | |_| | | | 
    |  |_____/ \__,_| \_/ \___|    |_|_|   \__,_|\__|_| |_| 
    |                                                   
  */
  art::Handle< std::vector<simb::MCParticle> > particleHandle;

  if (!evt.getByLabel(fGEANTLabel, particleHandle)) {
    throw cet::exception("AnalysisExample") 
      << " No simb::MCParticle objects in this event - " 
      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  for (size_t it=0; it<particleHandle->size(); it++) {
    const simb::MCParticle& part = (*particleHandle)[it];

    True_Bck_PDG       .push_back(part.PdgCode());
    True_Bck_Mother    .push_back(part.Mother ());
    True_Bck_ID        .push_back(part.TrackId());
    True_Bck_VertX     .push_back(part.Vx     ());
    True_Bck_VertY     .push_back(part.Vy     ());
    True_Bck_VertZ     .push_back(part.Vz     ());
    True_Bck_Time      .push_back(part.T      ());
    True_Bck_Energy    .push_back(part.E      ());
    True_Bck_EndE      .push_back(part.EndE   ());
    True_Bck_EndX      .push_back(part.EndX   ());
    True_Bck_EndY      .push_back(part.EndY   ());
    True_Bck_EndZ      .push_back(part.EndZ   ());
    True_Bck_EndT      .push_back(part.EndT   ());


    if(part.EndProcess() == "hadElastic"){
      True_Bck_EndProcess.push_back(1);}
    else if(part.EndProcess() == "protonInelastic"){
      True_Bck_EndProcess.push_back(2);}
    else if(part.EndProcess() == "neutronInelastic"){
      True_Bck_EndProcess.push_back(3);}
    else if(part.EndProcess() == "dInelastic"){
      True_Bck_EndProcess.push_back(4);}
    else if(part.EndProcess() == "tInelastic"){
      True_Bck_EndProcess.push_back(5);}
    else if(part.EndProcess() == "hBertiniCaptureAtRest"){
      True_Bck_EndProcess.push_back(6);}
    else if(part.EndProcess() == "alphaInelastic"){
      True_Bck_EndProcess.push_back(7);}
    else if(part.EndProcess() == "FastScintillation"){
      True_Bck_EndProcess.push_back(8);}
    else if(part.EndProcess() == "phot"){
      True_Bck_EndProcess.push_back(9);}
    else if(part.EndProcess() == "nCapture"){
      True_Bck_EndProcess.push_back(10);}
    else{
      True_Bck_EndProcess.push_back(0);
    }
  }

  // -------------------------------------
  
  fFullGeoAnaTree->Fill();
}

void FullGeoAna::endJob() {}

DEFINE_ART_MODULE(FullGeoAna)
