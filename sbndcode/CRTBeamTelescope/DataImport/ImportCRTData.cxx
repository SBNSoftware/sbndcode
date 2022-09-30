// See header file for documentation

// ART
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardata/Utilities/AssociationUtil.h"

// SBN
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "ImportCRTData.h"


namespace crt {

  ImportCRTData::ImportCRTData(fhicl::ParameterSet const & ps,
                               art::ProductRegistryHelper &helper,
                               art::SourceHelper const &pm)
    :
    fSourceHelper(pm),
    fSubRunID(),
    fTLfebdata{helper.reconstitutes<std::vector<sbnd::crt::FEBData>, art::InEvent>("crtdata")}
  {

    fPOT = fCurrentPOT = 0;
    helper.reconstitutes<sumdata::POTSummary, art::InSubRun >("crtdata");

    fEventCounter = 0;
    fMaxEvents = ps.get<int>("MaxEvents", -1);
    fVerbose   = ps.get<bool>("Verbose", false);

    fCableLengthCorrectionsVector = ps.get<std::vector<std::pair<unsigned, double>>>("CableLengthCorrectionsMap");
    fCableLengthCorrections       = std::map<unsigned, double>(fCableLengthCorrectionsVector.begin(), fCableLengthCorrectionsVector.end());

    fT1Offset = ps.get<unsigned>("T1Offset",100000);
  }

  void ImportCRTData::closeCurrentFile()
  {
    //mf::LogInfo(__FUNCTION__)<<"File boundary (processed "<<fEventCounter<<" events)"<<std::endl;
    fSubRunID.flushSubRun();
    fEventCounter=0;
    fCRTInputFile->Close();
    delete fCRTInputFile;
  }

  void ImportCRTData::readFile(std::string const &name,
                               art::FileBlock* &fb)
  {
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

    std::cout << "Reading in file " << name << std::endl;
    std::cout << "POT = " << pot << std::endl;
    fPOT+=pot;
    fCurrentPOT = pot;

    fTotalTreeEvents = fTree->GetEntries();
  }


  bool ImportCRTData::readNext(art::RunPrincipal* const &/*inR*/,
                               art::SubRunPrincipal* const &/*inSR*/,
                               art::RunPrincipal* &outR,
                               art::SubRunPrincipal* &outSR,
                               art::EventPrincipal* &outE)
  {
    if ((fMaxEvents > 0 && fEventCounter == unsigned(fMaxEvents)) || fEventCounter == fTotalTreeEvents) {
      return false;
    }

    // Create empty result, then fill it from current file:
    std::unique_ptr< std::vector<sbnd::crt::FEBData>  > febdata_v(new std::vector<sbnd::crt::FEBData>);

    if(fVerbose)
      std::cout << "==================================\n"
		<< "Event: " << fEventCounter << " / " << fTotalTreeEvents << '\n'
		<< "==================================" << std::endl;
    
    fTree->GetEntry(fEventCounter);

    std::array<uint16_t, 32> adc;

    for (size_t j = 0; j < fHit1Feb->size(); j++)
    {
      for (int s = 0; s < 32; s++) { adc[s] = fHit1Adc->at(j)[s]; }
      sbnd::crt::FEBData feb_data_1(fHit1Feb->at(j),
                                    fHit1T0->at(j) + fCableLengthCorrections[fHit1Feb->at(j)],
                                    fHit1T1->at(j) + fCableLengthCorrections[fHit1Feb->at(j)] + fT1Offset,
                                    adc,
                                    0);

      febdata_v->push_back(feb_data_1);

      for (int s = 0; s < 32; s++) { adc[s] = fHit2Adc->at(j)[s]; }
      sbnd::crt::FEBData feb_data_2(fHit2Feb->at(j),
                                    fHit2T0->at(j) + fCableLengthCorrections[fHit2Feb->at(j)],
                                    fHit2T1->at(j) + fCableLengthCorrections[fHit2Feb->at(j)] + fT1Offset,
                                    adc,
                                    0);

      febdata_v->push_back(feb_data_2);

      if(fVerbose)
	{
	  std::cout << "---------- FEB 1 ---------\n"
		    << "Mac5: " << feb_data_1.Mac5() << '\n'
		    << "T0:   " << feb_data_1.Ts0() << '\n'
		    << "T1:   " << feb_data_1.Ts1() << '\n'
		    << "ADC:  [";
	  for(auto const &adc : feb_data_1.ADC()) 
	    { std::cout << adc << ", ";}
	  std::cout << "]\n" << std::endl;

	  std::cout << "---------- FEB 2 ---------\n"
		    << "Mac5: " << feb_data_2.Mac5() << '\n'
		    << "T0:   " << feb_data_2.Ts0() << '\n'
		    << "T1:   " << feb_data_2.Ts1() << '\n'
		    << "ADC:  [";
	  for(auto const &adc : feb_data_2.ADC()) 
	    { std::cout << adc << ", ";}
	  std::cout << "]\n" << std::endl;
	}
    }
    fEventCounter++;

    art::RunNumber_t rn = 1000;
    if (rn==0) rn=999999;
    art::Timestamp tstamp(time(0));

    art::SubRunID newID(rn, 0);
    if (fSubRunID.runID() != newID.runID()) { // New Run
      outR = fSourceHelper.makeRunPrincipal(rn, tstamp);
      outSR = fSourceHelper.makeSubRunPrincipal(rn,0,tstamp);
      std::cout << "Made new run " << outSR->run() << ", subrun " << outSR->subRun() << std::endl;

      // Save POTs
      std::unique_ptr<sumdata::POTSummary> pot(new sumdata::POTSummary);
      pot->totpot = fPOT;
      pot->totgoodpot = fPOT;

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
