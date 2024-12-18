////////////////////////////////////////////////////////////////////////
// Class:       PTBAnalysis
// Plugin Type: analyzer
// File:        PTBAnalysis_module.cc
// Author:      Max Dubnowski (madubnowski@gmail.com)
// Purpose:     Analyze the outputs from the decoder
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"
#include <bitset>

#include "sbndcode/Decoders/PTB/sbndptb.h"

namespace sbnd::ptb {
  class PTBAnalysis;
}

class sbnd::ptb::PTBAnalysis : public art::EDAnalyzer {
public:
  explicit PTBAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PTBAnalysis(PTBAnalysis const&) = delete;
  PTBAnalysis(PTBAnalysis&&) = delete;
  PTBAnalysis& operator=(PTBAnalysis const&) = delete;
  PTBAnalysis& operator=(PTBAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec, int eventNum);


private:

  //CRTGeoAlg fCRTGeoAlg;
  //TPCGeoAlg fTPCGeoAlg;
  //CRTBackTrackerAlg fCRTBackTrackerAlg;

  std::string fMCParticleModuleLabel, fSimDepositModuleLabel, fFEBDataModuleLabel, fCRTStripHitModuleLabel,
    fCRTClusterModuleLabel, fCRTSpacePointModuleLabel, fCRTTrackModuleLabel, fTPCTrackModuleLabel,
    fCRTSpacePointMatchingModuleLabel, fCRTTrackMatchingModuleLabel, fPFPModuleLabel, fPTBModuleLabel,
    fTDCModuleLabel;
  bool fDebug, fDataMode, fNoTPC, fHasPTB, fHasTDC;

  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;


  std::vector<uint64_t> _ptb_hlt_trigger;
  std::vector<uint64_t> _ptb_hlt_timestamp;
  std::vector<int> _ptb_hlt_trunmask;
  std::vector<uint64_t> _ptb_llt_trigger;
  std::vector<uint64_t> _ptb_llt_timestamp;
  std::vector<int> _ptb_llt_trunmask;

  std::vector<uint64_t> _ptb_chStatus_timestamp;
  std::vector<uint64_t> _ptb_chStatus_beam;
  std::vector<uint64_t> _ptb_chStatus_crt;
  std::vector<uint64_t> _ptb_chStatus_pds;
  std::vector<uint64_t> _ptb_chStatus_mtca;
  std::vector<uint64_t> _ptb_chStatus_nim;
  std::vector<uint64_t> _ptb_chStatus_auxpds;

  std::vector<int> _ptb_crt_unmask;
  std::vector<int> _ptb_beam_unmask;
  std::vector<int> _ptb_pds_unmask;
  std::vector<int> _ptb_mtca_unmask;
  std::vector<int> _ptb_nim_unmask;
  std::vector<int> _ptb_auxpds_unmask;

  int _ptb_hlt_perevent;
  int _ptb_llt_perevent;
  
  int _ptb_crt_perevent;
  int _ptb_beam_perevent;
  int _ptb_pds_perevent;
  int _ptb_mtca_perevent;
  int _ptb_nim_perevent;
  int _ptb_auxpds_perevent;

};

sbnd::ptb::PTBAnalysis::PTBAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  
  fPTBModuleLabel                   = p.get<std::string>("PTBModuleLabel", "ptbdecoder");
  
  fDebug                            = p.get<bool>("Debug", false);
  fDataMode                         = p.get<bool>("DataMode", false);
  fNoTPC                            = p.get<bool>("NoTPC", false);
  fHasPTB                           = p.get<bool>("HasPTB", false);
  fHasTDC                           = p.get<bool>("HasTDC", false);
  
  art::ServiceHandle<art::TFileService> fs;

  fTree = fs->make<TTree>("tree","");
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);
  fTree->Branch("event", &_event);
  
  if(fHasPTB)
    {
      fTree->Branch("ptb_hlt_trigger", "std::vector<uint64_t>", &_ptb_hlt_trigger);
      fTree->Branch("ptb_hlt_trunmask", "std::vector<int>", &_ptb_hlt_trunmask);
      fTree->Branch("ptb_hlt_timestamp", "std::vector<uint64_t>", &_ptb_hlt_timestamp);
      fTree->Branch("ptb_llt_trigger", "std::vector<uint64_t>", &_ptb_llt_trigger);
      fTree->Branch("ptb_chStatus_timestamp", "std::vector<uint64_t>", &_ptb_chStatus_timestamp);
      fTree->Branch("ptb_chStatus_beam", "std::vector<uint64_t>", &_ptb_chStatus_beam);
      fTree->Branch("ptb_chStatus_crt", "std::vector<uint64_t>", &_ptb_chStatus_crt);
      fTree->Branch("ptb_chStatus_pds", "std::vector<uint64_t>", &_ptb_chStatus_pds);
      fTree->Branch("ptb_chStatus_mtca", "std::vector<uint64_t>", &_ptb_chStatus_mtca);
      fTree->Branch("ptb_chStatus_nim", "std::vector<uint64_t>", &_ptb_chStatus_nim);
      fTree->Branch("ptb_chStatus_auxpds", "std::vector<uint64_t>", &_ptb_chStatus_auxpds);
      fTree->Branch("ptb_llt_trunmask", "std::vector<int>", &_ptb_llt_trunmask);
      fTree->Branch("ptb_llt_perevent", &_ptb_llt_perevent);
      fTree->Branch("ptb_hlt_perevent", &_ptb_hlt_perevent);
      fTree->Branch("ptb_llt_timestamp", "std::vector<uint64_t>", &_ptb_llt_timestamp);
      fTree->Branch("ptb_crt_unmask", "std::vector<int>", &_ptb_crt_unmask);
      fTree->Branch("ptb_beam_unmask", "std::vector<int>", &_ptb_beam_unmask);
      fTree->Branch("ptb_pds_unmask", "std::vector<int>", &_ptb_pds_unmask);
      fTree->Branch("ptb_nim_unmask", "std::vector<int>", &_ptb_nim_unmask);
      fTree->Branch("ptb_mtca_unmask", "std::vector<int>", &_ptb_mtca_unmask);
      fTree->Branch("ptb_auxpds_unmask", "std::vector<int>", &_ptb_auxpds_unmask);
      fTree->Branch("ptb_crt_perevent", &_ptb_crt_perevent);
      fTree->Branch("ptb_beam_perevent", &_ptb_beam_perevent);
      fTree->Branch("ptb_pds_perevent", &_ptb_pds_perevent);
      fTree->Branch("ptb_mtca_perevent", &_ptb_mtca_perevent);
      fTree->Branch("ptb_nim_perevent", &_ptb_nim_perevent);
      fTree->Branch("ptb_auxpds_perevent", &_ptb_auxpds_perevent);
    
    }
}


void sbnd::ptb::PTBAnalysis::analyze(art::Event const& e)
{
 
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

 

  if(fHasPTB)
    {
      // Get PTBs
      art::Handle<std::vector<raw::ptb::sbndptb>> PTBHandle;
      e.getByLabel(fPTBModuleLabel, PTBHandle);
      if(!PTBHandle.isValid()){
        std::cout << "PTB product " << fPTBModuleLabel << " not found..." << std::endl;
        throw std::exception();
      }
      std::vector<art::Ptr<raw::ptb::sbndptb>> PTBVec;
      art::fill_ptr_vector(PTBVec, PTBHandle);

      // Fill PTB variables
      std::cout << "Run: " << _run  << "  SubRun: " << _subrun  << "  Event: " << _event  << std::endl; 
      AnalysePTBs(PTBVec, _event);
    }
 
  // Fill the Tree
  fTree->Fill();
}


void sbnd::ptb::PTBAnalysis::AnalysePTBs(std::vector<art::Ptr<raw::ptb::sbndptb>> &PTBVec, int eventNum)
{
  /*

  //######################################################

  //Channel Status words

   unsigned nChWs = 0;
   for (auto const& ptb : PTBVec){
     nChWs += ptb->GetNChStatuses();
   }

   _ptb_chStatus_timestamp.resize(nChWs);


   //Running through the Channel word sizes to resize properly and fill the bitmasked channel statuses
    
   int chw_i = 0; int c_i =0;  int b_i=0;  int p_i =0; int m_i =0; int n_i=0; int a_i=0; //These are the number of channel words that are given for each type of word
   for (auto const& ptb : PTBVec)
     {
       for (unsigned i =0; i < ptb->GetNChStatuses(); i++)
	 {
	  
   	  _ptb_chStatus_timestamp[chw_i] = ptb->GetChStatuse(i).timestamp;  //typo in the header file
  
	  if (ptb->GetChStatuse(i).crt != 0){
	    _ptb_chStatus_crt.resize(c_i+1);
	    _ptb_chStatus_crt[c_i] = ptb->GetChStatuse(i).crt;
	    c_i++;
	  }

	  
	  if (ptb->GetChStatuse(i).beam != 0){
	    _ptb_chStatus_beam.resize(b_i+1);
	    _ptb_chStatus_beam[b_i] = ptb->GetChStatuse(i).beam;
	    b_i++;
	  }


	  if (ptb->GetChStatuse(i).pds != 0){
	    _ptb_chStatus_pds.resize(p_i+1);
	    _ptb_chStatus_pds[p_i] = ptb->GetChStatuse(i).pds;
	    p_i++;
	  }


	  if (ptb->GetChStatuse(i).mtca != 0){
	    _ptb_chStatus_mtca.resize(m_i+1);
	    _ptb_chStatus_mtca[m_i] = ptb->GetChStatuse(i).mtca;
	    m_i++;
	  }

	  if (ptb->GetChStatuse(i).nim != 0){
	    _ptb_chStatus_nim.resize(n_i+1);
	    _ptb_chStatus_nim[n_i] = ptb->GetChStatuse(i).nim;
	    n_i++;
	  }

	  if (ptb->GetChStatuse(i).auxpds != 0){
	    _ptb_chStatus_auxpds.resize(a_i+1);
	    _ptb_chStatus_auxpds[a_i] = ptb->GetChStatuse(i).auxpds;
	    a_i++;
	  }
	  
	  chw_i++;
	  }//End of Channel statuses loop
     } //End of PTBVec loop
   
   //Resize with Proper channel size
   _ptb_crt_unmask.resize(c_i);
   _ptb_beam_unmask.resize(b_i);
   _ptb_pds_unmask.resize(p_i);
   _ptb_mtca_unmask.resize(m_i);
   _ptb_nim_unmask.resize(n_i);
   _ptb_auxpds_unmask.resize(a_i);
   

   //Looking at Unmasking CRT Channel Statuses
   
   unsigned crt_i =0; 
   for(auto const& ptb : PTBVec)
     {
      for(unsigned i = 0; i < ptb->GetNChStatuses(); ++i)
        {

	  int val = ptb->GetChStatuse(i).crt;
	  int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	  int numOfTrig =0;
	  for(int i=0; i<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = i;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	    
	  if (numOfTrig ==1){
	    _ptb_crt_unmask[crt_i] = upBit[0];
	    crt_i++;
	  }
	  
	  else if (numOfTrig > 1){
	    c_i += (numOfTrig -1);
	    _ptb_crt_unmask.resize(c_i);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_crt_unmask[crt_i] = upBit[mult];
	      crt_i++;
	    } //End of loop over numOfTrig 
	  } //End of else statement for multiple trigs
	} //End of loop over all channel statuses
    } //End of loop over elements in container PTBVec




   // Looking at Beam Channel Statuses

   unsigned beam_i =0; 
   for(auto const& ptb : PTBVec)
     {
      for(unsigned i = 0; i < ptb->GetNChStatuses(); ++i)
        {

	  int val = ptb->GetChStatuse(i).beam;
	  int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	  int numOfTrig =0;
	  for(int i=0; i<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = i;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	    
	  if (numOfTrig ==1){
	    _ptb_beam_unmask[beam_i] = upBit[0];
	    beam_i++;
	  }
	  
	  else if (numOfTrig > 1){
	    b_i += (numOfTrig -1);
	    _ptb_beam_unmask.resize(b_i);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_beam_unmask[beam_i] = upBit[mult];
	      beam_i++;
	    } //End of loop over numOfTrig 
	  } //End of else statement for multiple trigs
	} //End of loop over all channel statuses
    } //End of loop over elements in container PTBVec



   //Looking at Unmasking PDS Channel Statuses
   
   unsigned pds_i =0; 
   for(auto const& ptb : PTBVec)
     {
      for(unsigned i = 0; i < ptb->GetNChStatuses(); ++i)
        {

	  int val = ptb->GetChStatuse(i).pds;
	  int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	  int numOfTrig =0;
	  for(int i=0; i<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = i;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	    
	  if (numOfTrig ==1){
	    _ptb_pds_unmask[pds_i] = upBit[0];
	    pds_i++;
	  }
	  
	  else if (numOfTrig > 1){
	    p_i += (numOfTrig -1);
	    _ptb_pds_unmask.resize(p_i);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_pds_unmask[pds_i] = upBit[mult];
	      pds_i++;
	    } //End of loop over numOfTrig 
	  } //End of else statement for multiple trigs
	} //End of loop over all channel statuses
    } //End of loop over elements in container PTBVec



   //Looking at Unmasking MTCA Channel Statuses
   
   unsigned mtca_i =0; 
   for(auto const& ptb : PTBVec)
     {
      for(unsigned i = 0; i < ptb->GetNChStatuses(); ++i)
        {

	  int val = ptb->GetChStatuse(i).mtca;
	  int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	  int numOfTrig =0;
	  for(int i=0; i<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = i;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	    
	  if (numOfTrig ==1){
	    _ptb_mtca_unmask[mtca_i] = upBit[0];
	    mtca_i++;
	  }
	  
	  else if (numOfTrig > 1){
	    m_i += (numOfTrig -1);
	    _ptb_mtca_unmask.resize(m_i);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_mtca_unmask[mtca_i] = upBit[mult];
	      mtca_i++;
	    } //End of loop over numOfTrig 
	  } //End of else statement for multiple trigs
	} //End of loop over all channel statuses
    } //End of loop over elements in container PTBVec




   //Looking at Unmasking NIM Channel Statuses
   
   unsigned nim_i =0; 
   for(auto const& ptb : PTBVec)
     {
      for(unsigned i = 0; i < ptb->GetNChStatuses(); ++i)
        {

	  int val = ptb->GetChStatuse(i).nim;
	  int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	  int numOfTrig =0;
	  for(int i=0; i<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = i;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	    
	  if (numOfTrig ==1){
	    _ptb_nim_unmask[nim_i] = upBit[0];
	    nim_i++;
	  }
	  
	  else if (numOfTrig > 1){
	    n_i += (numOfTrig -1);
	    _ptb_nim_unmask.resize(n_i);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_nim_unmask[nim_i] = upBit[mult];
	      nim_i++;
	    } //End of loop over numOfTrig 
	  } //End of else statement for multiple trigs
	} //End of loop over all channel statuses
    } //End of loop over elements in container PTBVec




   //Looking at Unmasking AuxPDS Channel Statuses
   
   unsigned auxpds_i =0; 
   for(auto const& ptb : PTBVec)
     {
      for(unsigned i = 0; i < ptb->GetNChStatuses(); ++i)
        {

	  int val = ptb->GetChStatuse(i).auxpds;
	  int upBit[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	  int numOfTrig =0;
	  for(int i=0; i<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = i;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	    
	  if (numOfTrig ==1){
	    _ptb_auxpds_unmask[auxpds_i] = upBit[0];
	    auxpds_i++;
	  }
	  
	  else if (numOfTrig > 1){
	    a_i += (numOfTrig -1);
	    _ptb_auxpds_unmask.resize(a_i);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_auxpds_unmask[a_i] = upBit[mult];
	      auxpds_i++;
	    } //End of loop over numOfTrig 
	  } //End of else statement for multiple trigs
	} //End of loop over all channel statuses
    } //End of loop over elements in container PTBVec

   _ptb_crt_perevent = c_i; _ptb_beam_perevent = b_i; _ptb_pds_perevent = p_i; 
   _ptb_mtca_perevent = m_i; _ptb_nim_perevent = n_i; _ptb_auxpds_perevent = a_i;



  */

   // #################################################################

  
  // HLT Words
  unsigned nHLTs = 0;

  for(auto const& ptb : PTBVec)
    nHLTs += ptb->GetNHLTriggers();
  
  _ptb_hlt_trigger.resize(nHLTs);
  _ptb_hlt_timestamp.resize(nHLTs);
  _ptb_hlt_trunmask.resize(nHLTs);
  
  
  unsigned hlt_i = 0; //For multiple upbits in trigger words for unmasking
  unsigned h_i = 0; //For trigger with bitmask
  for(auto const& ptb : PTBVec)
    {
      for(unsigned i = 0; i < ptb->GetNHLTriggers(); ++i)
        {
	  _ptb_hlt_trigger[h_i] = ptb->GetHLTrigger(i).trigger_word;
	  _ptb_hlt_timestamp[h_i] = ptb->GetHLTrigger(i).timestamp; //Units can be found in the Decoder Module 
	  h_i++;

	  int val = ptb->GetHLTrigger(i).trigger_word;
	  int upBit[32];

	  //Set all default values of upBit to -1, in case there are multiple HLTs in the same PTB word
	  for (int u=0; u< 32; u++){
	    upBit[u] = -1;
	  }
	  
	  int numOfTrig =0;
	  for(int b=0; b<32;i++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = b;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }
	  
	  if (numOfTrig ==1){
	    _ptb_hlt_trunmask[hlt_i] = upBit[0];
	    hlt_i++;
	  }//End of if statement for single upbit
	  
	  else if (numOfTrig > 1){
	    nHLTs += (numOfTrig -1);
	    _ptb_hlt_timestamp.resize(nHLTs);
	    _ptb_hlt_trunmask.resize(nHLTs);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_hlt_trunmask[hlt_i] = upBit[mult];
	      hlt_i++;
	    } //End of loop over multiple upbits
	  } //End of else statement for multiple triggers
        } //End of loop over nHLTriggers
    } //End of loop over ptb in PTBVec
  
  _ptb_hlt_perevent = hlt_i;






  // ######################################################

  //LLT Words
  
  unsigned nLLTs = 0;

  for(auto const& ptb : PTBVec)
    nLLTs += ptb->GetNLLTriggers();

  _ptb_llt_trigger.resize(nLLTs);
  _ptb_llt_timestamp.resize(nLLTs);
  _ptb_llt_trunmask.resize(nLLTs);
  
  unsigned llt_i = 0; //For ptb_llt_trunmask and perevent, includes the multiples with multiple upbits per word
  unsigned l_i=0; //For ptb_llt_trigger, records only unique words
  for(auto const& ptb : PTBVec)
    {
      for(unsigned i = 0; i < ptb->GetNLLTriggers(); ++i)
        {
          _ptb_llt_trigger[l_i]   = ptb->GetLLTrigger(i).trigger_word;
	  _ptb_llt_timestamp[l_i] = ptb->GetLLTrigger(i).timestamp; //Units can be found in the Decoder Module    

	  l_i++;
	 
	  int val = ptb->GetLLTrigger(i).trigger_word;
	  
	  int upBit[32];
	  
	  for (int u =0; u< 32; u++){
	    upBit[u] = -1;
	  }

	  int numOfTrig =0;
	  for(int b=0; b<32;b++){
	    if ((val & 0x01) ==1){
	      upBit[numOfTrig] = b;
	      numOfTrig++;
	    }
	    val = val >> 1;
	  }//End of loop over bits
	 

	  if (numOfTrig ==1){
	    _ptb_llt_trunmask[llt_i] = upBit[0];
	    llt_i++;
	  } //End of if statement for single triggers
	  
	  else if (numOfTrig > 1){
	    nLLTs += (numOfTrig -1);
	    _ptb_llt_timestamp.resize(nLLTs);
	    _ptb_llt_trunmask.resize(nLLTs);

	    for (int mult =0; mult < numOfTrig; mult++){ 
	      _ptb_llt_trunmask[llt_i] = upBit[mult];
	      llt_i++;
	    } //End of loop over multiple upBits
	  } //End of else statement for multiple triggers
	} //End of loop over nLLTriggers
    } //End of Loop over ptb in PTBVec

  _ptb_llt_perevent = llt_i;
 
  std::cout << std::endl;

}


DEFINE_ART_MODULE(sbnd::ptb::PTBAnalysis)
