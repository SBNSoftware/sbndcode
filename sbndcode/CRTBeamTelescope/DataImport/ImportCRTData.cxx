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
    fSpills = fCurrentSpills = 0;
    helper.reconstitutes<sumdata::POTSummary, art::InSubRun >("crtdata");

    fEventCounter = 0;
    fMaxEvents = ps.get<int>("MaxEvents", -1);
    fVerbose   = ps.get<bool>("Verbose", false);

    fMac5ToGeoIDVec = ps.get<std::vector<std::pair<unsigned, unsigned>>>("FEBMac5ToGeometryIDMap");
    fMac5ToGeoID    = std::map<unsigned, unsigned>(fMac5ToGeoIDVec.begin(), fMac5ToGeoIDVec.end());
    
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

    //fCRTInputFile=new TFile(name.c_str());
    fCRTInputFile = TFile::Open(name.c_str());
    fCRTInputFileName = name;

    // the name will be 
    if (fCRTInputFile->IsZombie()) {
      throw cet::exception(__PRETTY_FUNCTION__) << "Failed to open "<<name<<std::endl;
    }else{
      fTree = (TTree*)(fCRTInputFile->Get("t"));

      fTree->SetBranchAddress("s", &fUnixS);

      fTree->SetBranchAddress("hit1_feb", &fHit1Feb);
      fTree->SetBranchAddress("hit1_flags", &fHit1Flags);
      fTree->SetBranchAddress("hit1_t0", &fHit1T0);
      fTree->SetBranchAddress("hit1_t1", &fHit1T1);
      fTree->SetBranchAddress("hit1_adc", &fHit1Adc);

      fTree->SetBranchAddress("hit2_feb", &fHit2Feb);
      fTree->SetBranchAddress("hit2_flags", &fHit2Flags);
      fTree->SetBranchAddress("hit2_t0", &fHit2T0);
      fTree->SetBranchAddress("hit2_t1", &fHit2T1);
      fTree->SetBranchAddress("hit2_adc", &fHit2Adc);


      TTree * aux = (TTree*)(fCRTInputFile->Get("aux"));
      double pot = 0., spills = 0.;
      aux->SetBranchAddress("pot", &pot);
      aux->SetBranchAddress("spills", &spills);
      aux->GetEntry(0);

      std::cout << "Reading in file " << name << std::endl;
      std::cout << "POT = " << pot << std::endl;
      std::cout << "Spills = " << spills << std::endl;
      fPOT           += pot;
      fCurrentPOT    =  pot;
      fSpills        += spills;
      fCurrentSpills =  spills;

      fTotalTreeEvents = fTree->GetEntries();
    }
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
      int max_s = 999999; uint16_t max_adc = 0;
      for (int s = 0; s < 32; s++) 
      { 
        adc[s] = fHit1Adc->at(j)[s];
        if(adc[s] > max_adc) { max_adc = adc[s]; max_s = s;}
      }

      sbnd::crt::FEBData feb_data_1(fMac5ToGeoID[fHit1Feb->at(j)],
                                    fHit1Flags->at(j),
                                    fHit1T0->at(j),
                                    fHit1T1->at(j) + fT1Offset,
                                    fUnixS->at(j),
                                    adc,
                                    max_s);

      febdata_v->push_back(feb_data_1);

      max_s = 999999; max_adc = 0;

      for (int s = 0; s < 32; s++) 
      { 
        adc[s] = fHit2Adc->at(j)[s];
        if(adc[s] > max_adc) { max_adc = adc[s]; max_s = s;}
      }
      sbnd::crt::FEBData feb_data_2(fMac5ToGeoID[fHit2Feb->at(j)],
                                    fHit2Flags->at(j),
                                    fHit2T0->at(j),
                                    fHit2T1->at(j) + fT1Offset,
                                    fUnixS->at(j),
                                    adc,
                                    max_s);

      febdata_v->push_back(feb_data_2);

      if(fVerbose)
      {
        std::cout << "---------- FEB 1 ---------\n"
                  << "Mac5:  " << feb_data_1.Mac5() << '\n'
                  << "Flags: " << feb_data_1.Flags() << '\n'
                  << "T0:    " << feb_data_1.Ts0() << '\n'
                  << "T1:    " << feb_data_1.Ts1() << '\n'
                  << "UnixS: " << feb_data_1.UnixS() << '\n'
                  << "ADC:   [";
        for(auto const &adc : feb_data_1.ADC()) 
          { std::cout << adc << ", ";}
        std::cout << "]\n" << std::endl;

        std::cout << "---------- FEB 2 ---------\n"
                  << "Mac5:  " << feb_data_2.Mac5() << '\n'
                  << "Flags: " << feb_data_2.Flags() << '\n'
                  << "T0:    " << feb_data_2.Ts0() << '\n'
                  << "T1:    " << feb_data_2.Ts1() << '\n'
                  << "UnixS: " << feb_data_2.UnixS() << '\n'
                  << "ADC:   [";
        for(auto const &adc : feb_data_2.ADC()) 
          { std::cout << adc << ", ";}
        std::cout << "]\n" << std::endl;
      }
    }
    fEventCounter++;

    /*art::RunNumber_t rn = 1000;
    if (rn==0) rn=999999;*/

    art::RunNumber_t runNumber; int subrunNumber; std::string dateTimeStr;
    convertFileNameToDateTimeRunSubRun(fCRTInputFileName, dateTimeStr, runNumber, subrunNumber);

    uint64_t timeStamp_thisFile = dateToUnixEpochNs(dateTimeStr);

    //art::Timestamp tstamp(time(0));
    art::Timestamp tstamp(timeStamp_thisFile);
    art::SubRunID newID(runNumber, subrunNumber);

    if (fVerbose){
      std::cout<<"Date: "<<dateTimeStr<<std::endl;
      std::cout<<runNumber<<" "<<tstamp.timeHigh()<<" "<<tstamp.timeLow()<<" "<<newID<<std::endl;
      std::cout<<"timestamp: "<<timeStamp_thisFile<<", in human-readable format: "<<unixEpochNanosecondsToDateTime(timeStamp_thisFile)<<std::endl;
    }
    
    if (fSubRunID.runID() != newID.runID()) { // New Run
      outR = fSourceHelper.makeRunPrincipal(runNumber, tstamp);
    }

    if (fSubRunID != newID) { // New SubRun
      outSR = fSourceHelper.makeSubRunPrincipal(runNumber, subrunNumber, tstamp);
      std::cout << "Made new run " << outSR->run() << ", subrun " << outSR->subRun() << std::endl;

      // Save POTs
      std::unique_ptr<sumdata::POTSummary> pot(new sumdata::POTSummary);
      pot->totpot     = fPOT;
      pot->totgoodpot = fPOT;
      pot->totspills  = fSpills;
      pot->goodspills = fSpills;

      fSubRunID = newID;
      // this means products of sumdata::POTSummary is being written. 

      art::put_product_in_principal(std::move(pot), //
                                    *outSR,
                                    "crtdata");
      //fSubRunID = newID;

    }
    
    std::sort(febdata_v->begin(), febdata_v->end(),
              [](const sbnd::crt::FEBData &a, const sbnd::crt::FEBData &b) -> bool {
                return a.Mac5() < b.Mac5() || (a.Mac5() == b.Mac5() && a.Ts1() < b.Ts1());
              });

    febdata_v->erase(std::unique(febdata_v->begin(), febdata_v->end(),
				 [](const sbnd::crt::FEBData &a, const sbnd::crt::FEBData &b) -> bool {
				   return a.Mac5() == b.Mac5() && a.Ts1() == b.Ts1();
				 }), febdata_v->end());

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


  uint64_t ImportCRTData::dateToUnixEpochNs(const std::string &dateStr) {
    // Parse the date string into a tm struct
    std::tm timeinfo = {};
    std::istringstream ss(dateStr);
    ss >> std::get_time(&timeinfo, "%Y-%m-%d %H:%M:%S");

    // Convert the tm struct into a time_t value using timegm instead of mktime
    time_t unixTime = timegm(&timeinfo);

    // Convert the time_t value into nanoseconds
    return unixTime * 1000000000;
  }

  void ImportCRTData::convertFileNameToDateTimeRunSubRun(const std::string& fileName, std::string& dateTimeStr, art::RunNumber_t & runNumber, int& subrunNumber) {
      std::string dateStr = fileName.substr(7, 8); // Extract the date part from the file name
      // Extract the year, month, and day from the date string
      std::string year  = dateStr.substr(0, 4);
      std::string month = dateStr.substr(4, 2);
      std::string day   = dateStr.substr(6, 2);
      std::string formattedDate = year + "-" + month + "-" + day;

      // Extract the hour, minute, and second from the time string
      std::string timeStr = fileName.substr(16, 6); // Extract the time part from the file name
      std::string hour    = timeStr.substr(0, 2);
      std::string minute  = timeStr.substr(2, 2);
      std::string second  = timeStr.substr(4, 2);

      // Convert the year, month, day, hour, minute, and second to integers
      int yearInt   = std::stoi(year);
      int monthInt  = std::stoi(month);
      int dayInt    = std::stoi(day);
      int hourInt   = std::stoi(hour);
      int minuteInt = std::stoi(minute);
      int secondInt = std::stoi(second);

      // form dataTimeStr
      dateTimeStr = formattedDate + " " + hour + ":" + minute + ":" + second;

      // Calculate the run number based on the day
      runNumber = (yearInt-2017+1) * 1000 + monthInt * 100 + dayInt;

      // Calculate the subrun number based on the time
      subrunNumber = hourInt * 10 + minuteInt + secondInt;
  }


  std::string ImportCRTData::unixEpochNanosecondsToDateTime(uint64_t unixEpochNanoseconds) {
      time_t unixTime = unixEpochNanoseconds / 1000000000; // Convert nanoseconds to seconds
      struct tm timeinfo;
      gmtime_r(&unixTime, &timeinfo);

      char dateTimeStr[20];
      strftime(dateTimeStr, sizeof(dateTimeStr), "%Y-%m-%d %H:%M:%S", &timeinfo);

      return dateTimeStr;
  }
}