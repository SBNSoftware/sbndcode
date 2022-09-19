// See header file for documentation

//LArSoft
#include "larcoreobj/SummaryData/POTSummary.h"

//ART, ...
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCFlux.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/NuChoice.h"

#include "sbnobj/SBND/CRT/FEBData.hh"

#include "lardata/Utilities/AssociationUtil.h"

#include "ImportCRTData.h"

namespace crt {

  ImportCRTData::ImportCRTData(fhicl::ParameterSet const & ps,
                               art::ProductRegistryHelper &helper,
                               art::SourceHelper const &pm)
    :
    fSourceHelper(pm),
    fSubRunID(),
    fEventCounter(0),
    fInputType(ps.get<std::string>("inputType")),
    fSelfIncrementRuns(ps.get<bool>("SelfIncrementRun")),
    fIncrement(1),
    fTLfebdata{helper.reconstitutes<std::vector<sbnd::crt::FEBData>, art::InEvent>("data")}
    // fTLmctruth{helper.reconstitutes<std::vector<simb::MCTruth>, art::InEvent>("flux")},
    // fTLmcflux{helper.reconstitutes<std::vector<simb::MCFlux>, art::InEvent>("flux")},
    // fTLdk2nu{fInputType=="dk2nu" ? helper.reconstitutes<std::vector<bsim::Dk2Nu>, art::InEvent>("flux"):fTLmctruth},
    // fTLnuchoice{fInputType=="dk2nu" ? helper.reconstitutes<std::vector<bsim::NuChoice>, art::InEvent>("flux"):fTLmctruth}
  {

    fPOT = fCurrentPOT = 0;

    // helper.reconstitutes<sumdata::POTSummary, art::InSubRun >("flux");

    // if (fInputType=="dk2nu") {
    //   helper.reconstitutes<art::Assns<simb::MCTruth, bsim::Dk2Nu>, art::InEvent>("flux");
    //   helper.reconstitutes<art::Assns<simb::MCTruth, bsim::NuChoice>, art::InEvent>("flux");
    //   helper.reconstitutes<art::Assns<simb::MCTruth, simb::MCFlux>, art::InEvent>("flux");

    //   fConfigPS=ps.get<fhicl::ParameterSet>(ps.get<std::string>("dk2nuConfig"));
    // }

  }

  void ImportCRTData::closeCurrentFile()
  {
    //mf::LogInfo(__FUNCTION__)<<"File boundary (processed "<<fEventCounter<<" events)"<<std::endl;
    fSubRunID.flushSubRun();
    fEventCounter=0;
    fCRTInputFile->Close();
    delete fCRTInputFile;
    fSkipEvents=0; //if crossing file boundary don't skip events in next file
    fEntry=0;
  }

  void ImportCRTData::readFile(std::string const &name,
                               art::FileBlock* &fb)
  {
    std::cout << "readFile " << name << std::endl;
    // Fill and return a new Fileblock.
    fb = new art::FileBlock(art::FileFormatVersion(1, "ImportCRTData"), name);

    fCRTInputFile=new TFile(name.c_str());
    if (fCRTInputFile->IsZombie()) {
      //throw cet::exception(__PRETTY_FUNCTION__) << "Failed to open "<<fCRTInputFile<<std::endl;
    }

    fTree = (TTree*)(fCRTInputFile->Get("t"));
    fTree->SetBranchAddress("hit1_feb", &fHit1Feb);
    fTree->SetBranchAddress("hit1_t0", &fHit1T0);
    fTree->SetBranchAddress("hit1_t1", &fHit1T1);
    fTree->SetBranchAddress("hit1_adc", &fHit1Adc);

    fTree->SetBranchAddress("hit2_feb", &fHit2Feb);
    fTree->SetBranchAddress("hit2_t0", &fHit2T0);
    fTree->SetBranchAddress("hit2_t1", &fHit2T1);
    fTree->SetBranchAddress("hit2_adc", &fHit2Adc);



    TTree * aux = (TTree*)(fCRTInputFile->Get("aux"));

    double pot = 0.;
    aux->SetBranchAddress("pot", &pot);
    aux->GetEntry(0);

    // fTree = (TTree*)(fCRTInputFile->Get("tree"));

    // fTree->SetBranchAddress("mac5", &fMac5);
    // fTree->SetBranchAddress("flags", &flags);
    // fTree->SetBranchAddress("lostcpu", &lostcpu);
    // fTree->SetBranchAddress("lostfpga", &lostfpga);
    // fTree->SetBranchAddress("ts0", &ts0);
    // fTree->SetBranchAddress("ts1", &ts1);
    // fTree->SetBranchAddress("adc", &adc);
    // fTree->SetBranchAddress("s", &s);
    // fTree->SetBranchAddress("ms", &ms);
    // fTree->SetBranchAddress("cable", &cable);
    // fTree->GetEntry(0);
    // std::cout << "--------> mac5 " << fMac5 << std::endl;


    // if (fInputType=="gsimple") {
      // fFluxDriver=new GSimpleInterface();
      // ((GSimpleInterface*)fFluxDriver)->SetRootFile(fCRTInputFile);
    // } else if (fInputType=="dk2nu") {
    //   fFluxDriver=new DK2NuInterface();
    //   ((DK2NuInterface*)fFluxDriver)->SetRootFile(fCRTInputFile);
    //   ((DK2NuInterface*)fFluxDriver)->Init(fConfigPS);
    // } else if (fInputType=="boone") {
    //   fFluxDriver=new BooNEInterface();
    //   ((BooNEInterface*)fFluxDriver)->SetRootFile(fCRTInputFile);
    // } else {
    //   throw cet::exception(__PRETTY_FUNCTION__) << "Ntuple format "<<fInputType<<" not supported"<<std::endl;
    // }
    std::cout << "Reading in file " << name << std::endl;
    std::cout << "POT = " << pot << std::endl;
    fPOT+=pot;
    fCurrentPOT = pot;

    // if (fSelfIncrementRuns) {
    //   fIncrement ++;
    //   fFluxDriver->SetRun(fFluxDriver->GetRun() + fIncrement);
    // }
  }


  bool ImportCRTData::readNext(art::RunPrincipal* const &/*inR*/,
                               art::SubRunPrincipal* const &/*inSR*/,
                               art::RunPrincipal* &outR,
                               art::SubRunPrincipal* &outSR,
                               art::EventPrincipal* &outE)
  {
    if (fMaxEvents > 0 && fEventCounter == unsigned(fMaxEvents))
      return false;

    // Create empty result, then fill it from current file:
    std::unique_ptr< std::vector<sbnd::crt::FEBData>  > febdata_v(new std::vector<sbnd::crt::FEBData>);


    fTree->GetEntry(fEventCounter);

    std::array<uint16_t, 32> adc;

    for (size_t j = 0; j < fHit1Feb.size(); j++)
    {
      for (int s = 0; s < 32; s++) { adc[s] = fHit1Adc[j][s]; }
      sbnd::crt::FEBData feb_data_1(fHit1Feb[j],
                                    fHit1T0[j],
                                    fHit1T1[j],
                                    adc,
                                    0);

      febdata_v->push_back(feb_data_1);

      for (int s = 0; s < 32; s++) { adc[s] = fHit2Adc[j][s]; }
      sbnd::crt::FEBData feb_data_2(fHit2Feb[j],
                                    fHit2T0[j],
                                    fHit2T1[j],
                                    adc,
                                    0);

      febdata_v->push_back(feb_data_2);
    }

    fEventCounter++;
    fEntry++;

    art::RunNumber_t rn = 1000;
    if (rn==0) rn=999999;
    art::Timestamp tstamp(time(0));

    art::SubRunID newID(rn, 0); //subrun not used in flux files, so set to 0
    if (fSubRunID.runID() != newID.runID()) { // New Run
      std::cout << "Making new run" << std::endl;
      outR = fSourceHelper.makeRunPrincipal(rn, tstamp);
    }
    if (fSubRunID != newID) { // New SubRun
      std::cout << "Making new subrun" << std::endl;
      outSR = fSourceHelper.makeSubRunPrincipal(rn,0,tstamp);
      std::unique_ptr<sumdata::POTSummary> pot(new sumdata::POTSummary);
      pot->totpot = fPOT;
      pot->totgoodpot = fPOT;

      std::cout << "run " << outSR->run() << ", subrun " << outSR->subRun() << std::endl;

      art::put_product_in_principal(std::move(pot),
                                    *outSR,
                                    "crtdata");

      fSubRunID = newID;
    }


    outE = fSourceHelper.makeEventPrincipal(fSubRunID.run(),
                                            fSubRunID.subRun(),
                                            fEventCounter,
                                            tstamp);

    // Put products in the event.
    art::put_product_in_principal(std::move(febdata_v),
                                  *outE,
                                  "crtdata"); // Module label


    return true;
  }

}
