////////////////////////////////////////////////////////////////////////
// \file    CVNZlibMaker_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
//           - wrote the zlib code used in this module
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <cstdlib>

#include "boost/filesystem.hpp"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Data products
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "sbndcode/RecoUtils/RecoUtils.h"
// #include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/LArTrainingData.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/PixelMap.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDICVNZlibMaker.h"

// Compression
#include "zlib.h"
#include "math.h"

#include "TH1.h"
#include "TTree.h"
#include "larcoreobj/SummaryData/POTSummary.h"


namespace lcvn {

  class SBNDCVNZlibMaker : public SBNDICVNZlibMaker {
  public:

    explicit SBNDCVNZlibMaker(fhicl::ParameterSet const& pset);
    ~SBNDCVNZlibMaker();

    void beginJob();
    void endSubRun(art::SubRun const &sr);
    void analyze(const art::Event& evt);
    void reconfigure(const fhicl::ParameterSet& pset);

  private:

    unsigned int fTopologyHitsCut;
    std::string fGenieGenModuleLabel;
    bool fVerbose;
    bool fUseBackTrackInfo;
    bool fUseNuContainment;
    double fCosEfrac;
    double fNuEfrac;
    double fVolCut;
    bool fSaveTree;
    std::string fPFParticleModuleLabel;
    std::string fHitModuleLabel;
    std::string fT0Label;
    
    // ParticleInventoryService
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    // BackTrackerService
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    void write_files(LArTrainingNuData td, std::string evtid);
    //bool Is_in_TPC(double stX, double stY, double stZ);
    bool Is_in_Fiducial_Vol(double stX, double stY, double stZ);
    void HitTruth(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& truthid);
    bool HitTruthId(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& mcid);
    bool TrackIdToMCTruth(Int_t const trkID, art::Ptr<simb::MCTruth>& mctruth);
    double HitEfrmTrkID(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t truthid);
    double HitTotE(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit);
    bool Is_Slice_Nu(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce);
    art::Ptr<recob::PFParticle> Get_Nu_like_PFP(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce);
    double Get_Slice_Score(art::FindManyP<larpandoraobj::PFParticleMetadata> const& assn, art::Ptr<recob::PFParticle> const& pfp);
    std::vector<int> Get_Best_Slice_ID_PFP_pdg(art::FindManyP<recob::PFParticle> const& pfp_assn, art::FindManyP<larpandoraobj::PFParticleMetadata> const& metadata_assn, std::vector<art::Ptr<recob::Slice>> const& slice_vec);
    int Get_Slice_PFP_ID(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce);
    int Get_Hit_Count(unsigned int tpc_no, unsigned int pln_no, std::vector<art::Ptr<recob::Hit>> const& hit_vec);
    double Get_Tot_nu_E(detinfo::DetectorClocksData const& clockData,art::Ptr<simb::MCTruth> const& mcneutrino, std::vector<art::Ptr<recob::Hit>> const& hit_vec);
    
    double HitTotE_frm_given_orgin(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, std::string origin_name);
    double HitTotE_frm_mctruth(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, art::Ptr<simb::MCTruth> const& my_mctruth);
    double TotE_from_mctruth(detinfo::DetectorClocksData const& clockData, std::vector<art::Ptr<recob::Hit>> const& hit_vec, art::Ptr<simb::MCTruth> const& my_mctruth);

    void Clear();

    TH1D* hPOT;
    double fPOT;
    int fRun;
    int fSubRun; 
    
    TTree *fEventTree;
    int frun;
    int fsubrun;
    int fevent;
    int fsliceID;
    double ftotsliceE;
    double ftotsliceNuE;
    double ftotsliceCosE;
    double ftotsliceOthE;
    double fsliceNuEfrac;
    double fsliceCosEfrac;
    double fsliceOthEfrac;	    
    bool fIsSliceNu;
    double fSliceScore;	
    bool fIsbestSlice;
    int fbestpfppdg;
    int fbestsliceID;
    int fpfppdg;
    bool fNuIsContained;
    int fNuID;
    int fNhits_tpc_0_pl_0;
    int fNhits_tpc_0_pl_1;
    int fNhits_tpc_0_pl_2;
    int fNhits_tpc_1_pl_0;
    int fNhits_tpc_1_pl_1;
    int fNhits_tpc_1_pl_2;
    int fNhits_total;
    double fNuSlice_purity;
    double fNuSlice_completeness;
    double fT0;
  };

  //......................................................................
  SBNDCVNZlibMaker::SBNDCVNZlibMaker(fhicl::ParameterSet const& pset)
    : SBNDICVNZlibMaker(pset)
  {
    reconfigure(pset);
  }

  //......................................................................
  SBNDCVNZlibMaker::~SBNDCVNZlibMaker()
  {  }

  //......................................................................
  void SBNDCVNZlibMaker::reconfigure(const fhicl::ParameterSet& pset)
  {
    SBNDICVNZlibMaker::reconfigure(pset);
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");
    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
    fVerbose = pset.get<bool>("Verbose");
    fUseBackTrackInfo = pset.get<bool>("UseBackTrackInfo");
    fUseNuContainment = pset.get<bool>("UseNuContainment");
    fCosEfrac = pset.get<double>("CosEfrac");
    fNuEfrac = pset.get<double>("NuEfrac");
    fVolCut = pset.get<double>("VolCut");
    fSaveTree = pset.get<bool>("SaveTree");
    fPFParticleModuleLabel = pset.get<std::string>("PFParticleModuleLabel");
    fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
    fT0Label = pset.get<std::string>("T0Label");
  }

  //......................................................................
  void SBNDCVNZlibMaker::endSubRun(const art::SubRun & sr){

    std::string fPOTModuleLabel = "generator";
    fRun = sr.run();
    fSubRun = sr.subRun();

    art::Handle< sumdata::POTSummary > potListHandle;
    if(sr.getByLabel(fPOTModuleLabel,potListHandle))
      fPOT = potListHandle->totpot;
    else
      fPOT = 0.;
    if(hPOT) hPOT->Fill(0.5, fPOT);
  }

  //......................................................................
  
  void SBNDCVNZlibMaker::beginJob()
  {
    SBNDICVNZlibMaker::beginJob();

    art::ServiceHandle<art::TFileService> tfs;
    hPOT = tfs->make<TH1D>("TotalPOT", "Total POT;; POT", 1, 0, 1);
  }

  //......................................................................
  
  void SBNDCVNZlibMaker::analyze(const art::Event& evt)
  {	  
    std::vector<art::Ptr<lcvn::PixelMap>> pixelmaps;
    art::InputTag itag1(fPixelMapInput, "cvnmap");
    auto h_pixelmaps = evt.getHandle<std::vector<lcvn::PixelMap>>(itag1);
    if (h_pixelmaps)
      art::fill_ptr_vector(pixelmaps, h_pixelmaps); 	  
	  
    if (pixelmaps.size() == 0) return;

    // Get associated slice for each pixel map
    art::FindOneP<recob::Slice> findOneSlice(h_pixelmaps, evt, itag1);
           	  
    AssignLabels labels;
    TDNuInfo info;
    double event_weight =1;
    InteractionType interaction = kOther;
    bool Nu_is_contained = false;
       
    std::vector<art::Ptr<simb::MCTruth>> mctruth_list;
    auto h_mctruth = evt.getHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
    if (h_mctruth)
      art::fill_ptr_vector(mctruth_list, h_mctruth);
       
    for(auto const& mctruth : mctruth_list){
      if(mctruth->Origin() == simb::kBeamNeutrino){
        simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();
	      
        if(fUseNuContainment){
          if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())){ 
            interaction = labels.GetInteractionType(true_neutrino);
            labels.GetTopology(mctruth, fTopologyHitsCut);
            float nu_energy = true_neutrino.Nu().E();
            float lep_energy = true_neutrino.Lepton().E();
            fNuID = true_neutrino.Nu().TrackId();
            info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
            info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), labels.GetNPizeros(), labels.GetNNeutrons(), labels.GetTopologyType(), labels.GetTopologyTypeAlt());
            Nu_is_contained = true;
            break;
          }
        }
	      
        else{
          interaction = labels.GetInteractionType(true_neutrino);
          labels.GetTopology(mctruth, fTopologyHitsCut);
          float nu_energy = true_neutrino.Nu().E();
          float lep_energy = true_neutrino.Lepton().E();
          fNuID = true_neutrino.Nu().TrackId();
          info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
          info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), labels.GetNPizeros(), labels.GetNNeutrons(), labels.GetTopologyType(), labels.GetTopologyTypeAlt());
          if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())) Nu_is_contained = true;
          break;
        }
	      
      } // select neutrino interactions
		   
		   
      else if(mctruth->Origin() == simb::kCosmicRay){
        interaction = kCosmic;
        break;
      } // select cosmic interactions
    }
       
       
    /////////////////////// use slice section ////////////////////////////////////////
       
    if(fUseSlice){
      std::cout << "***************** Used slice method ****************\n";
	  
      art::Handle< std::vector<recob::Slice> > SliceListHandle;
      std::vector< art::Ptr<recob::Slice> > SliceList;
      if( evt.getByLabel(fSliceLabel,SliceListHandle))
        art::fill_ptr_vector(SliceList,SliceListHandle);
	  
      art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
      std::vector< art::Ptr<recob::PFParticle> > PFPList;
      if( evt.getByLabel(fPFParticleModuleLabel,PFPListHandle))
        art::fill_ptr_vector(PFPList,PFPListHandle);
	  
      art::Handle< std::vector<recob::Hit> > HitListHandle;
      std::vector< art::Ptr<recob::Hit> > HitList;
      if( evt.getByLabel(fHitModuleLabel,HitListHandle))
        art::fill_ptr_vector(HitList,HitListHandle);
	  
	  
      art::FindManyP<recob::Hit> findManyHits(SliceListHandle, evt, fSliceLabel);
      art::FindManyP<recob::PFParticle> findManyPFPs(SliceListHandle, evt, fPFParticleModuleLabel);
      art::FindManyP<larpandoraobj::PFParticleMetadata> fm_pfpmd(PFPListHandle, evt, fPFParticleModuleLabel);
      art::FindManyP<anab::T0> findManyT0s(PFPListHandle, evt, fT0Label);
	  
      if(fUseBackTrackInfo){
        std::cout << "***************** Used backtracker information to get neutrino information ****************\n";
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
	     
        for(unsigned int i=0; i<pixelmaps.size(); i++){
          Clear();
          frun = evt.run();
          fsubrun = evt.subRun();
          fevent = evt.id().event();

          if (findOneSlice.isValid()){
            auto const & slice = findOneSlice.at(i);
		     
            std::cout << evt.run() << "  " << evt.subRun() << "  " << evt.event() << "  " << slice->ID() << "\n";
            AssignLabels labels;
            TDNuInfo info;
            double event_weight =1;
            InteractionType interaction = kOther;

            fsliceID = slice->ID();
            if(findManyHits.isValid()){
              std::vector<art::Ptr<recob::Hit>> slice_hits = findManyHits.at(slice.key());
              double tot_slice_eng = 0;
              double tot_slice_nu_eng = 0;
              double tot_slice_cos_eng = 0;
			   
              fNhits_tpc_0_pl_0 = Get_Hit_Count(0,0,slice_hits);
              fNhits_tpc_0_pl_1 = Get_Hit_Count(0,1,slice_hits);
              fNhits_tpc_0_pl_2 = Get_Hit_Count(0,2,slice_hits);
              fNhits_tpc_1_pl_0 = Get_Hit_Count(1,0,slice_hits);
              fNhits_tpc_1_pl_1 = Get_Hit_Count(1,1,slice_hits);
              fNhits_tpc_1_pl_2 = Get_Hit_Count(1,2,slice_hits);
              fNhits_total = slice_hits.size();
			   
              std::vector<double> mc_truth_eng(mctruth_list.size(),0.);
              for(auto const hit : slice_hits){
                //Int_t trkId; // newly deleted
                tot_slice_eng += HitTotE(clockData,hit);
                tot_slice_nu_eng += HitTotE_frm_given_orgin(clockData,hit,"nu"); // newly added
                tot_slice_cos_eng += HitTotE_frm_given_orgin(clockData,hit,"cos"); // newly added
			       
                // newly added section
			       
                if(mctruth_list.size()){
                  for(auto const truth : mctruth_list){
                    if(truth->Origin() == simb::kBeamNeutrino){
                      mc_truth_eng[truth.key()] += HitTotE_frm_mctruth(clockData,hit,truth);
                    }
                  }
                }
              } // loop over hits in the selected slice
			    
              ftotsliceE = tot_slice_eng;
    	
              if(tot_slice_eng > 0){
                std::cout << "Total energy : " << tot_slice_eng << "\n";
                std::cout << "Total cosmic energy : " << tot_slice_cos_eng << "\n";
                std::cout << "Total nu energy : " << tot_slice_nu_eng << "\n";
                std::cout << "Cosmic fraction : " << double(tot_slice_cos_eng)/tot_slice_eng << "\n";
                std::cout << "Neutrino fraction : " << double(tot_slice_nu_eng)/tot_slice_eng << "\n";
			      
                ftotsliceNuE = tot_slice_nu_eng;
                ftotsliceCosE = tot_slice_cos_eng;
                ftotsliceOthE = ftotsliceE - ftotsliceNuE - ftotsliceCosE;
                fsliceNuEfrac = double(ftotsliceNuE)/ftotsliceE;
                fsliceCosEfrac = double(ftotsliceCosE)/ftotsliceE;
                fsliceOthEfrac = double(ftotsliceOthE)/ftotsliceE;
			      
                ///////////////////////////////////////////////////// Check whether this has pandora T0 ///////////////////////////////
			      
                std::vector<float> pfp_T0_vec;
                if(findManyPFPs.isValid()){
                  std::vector<art::Ptr<recob::PFParticle>> slicePFPs = findManyPFPs.at(slice.key());
                  if(slicePFPs.size()){
                    for(auto const &pfp : slicePFPs){
                      if(findManyT0s.isValid()){
                        std::vector<art::Ptr<anab::T0>> T0_vec = findManyT0s.at(pfp.key());
                        if(T0_vec.size()){
                          for(auto const& T0 : T0_vec){
                            pfp_T0_vec.push_back(T0->Time());
                          }
                        }
                      }
                    }
                  }
                }
			      
                if(pfp_T0_vec.size()){
                  fT0 = *min_element(pfp_T0_vec.begin(), pfp_T0_vec.end());
                } 
			      
                ////////////////////////////////////////////////////// ////////////////////////////////////////////////////////////////
			      
                if((double(tot_slice_cos_eng)/tot_slice_eng) >= fCosEfrac){
                  interaction = kCosmic;
                } // select cosmic
                else if((double(tot_slice_nu_eng)/tot_slice_eng) >= fNuEfrac){
                  int index = -1;
                  double min_E = 0.; 
                  std::cout << "Size of the truth energy list : " << mc_truth_eng.size() << "\n";
                  for(unsigned int k=0; k<mc_truth_eng.size(); k++){
                    if(mc_truth_eng[k] > min_E){
                      //std::cout << "MC truth energy : " << mc_truth_eng[k] << "\n";	    
                      art::Ptr<simb::MCTruth> mctruth = mctruth_list[k];
                      simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();
                      //std::cout << "Neutrino vertext : " << true_neutrino.Nu().Vx() << "  " << true_neutrino.Nu().Vy() << "  " << true_neutrino.Nu().Vz() << "\n";
				       
                      if(fUseNuContainment){
                        if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())){ 
                          //std::cout << "================ Checking containment cut =================\n";
                          index = k;
                          min_E = mc_truth_eng[k];
                          fNuIsContained = true;
                        }
                      }
				       
                      else{
                        //std::cout << "================ Not Checking containment cut =================\n";
                        index = k;
                        min_E = mc_truth_eng[k];
                        if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())) fNuIsContained = true;
                        else fNuIsContained = false;
                      }
				       
                      /*if(fUseNuFullContainment){
                        if(Is_in_Fiducial_Vol(true_neutrino.Nu().Vx(),true_neutrino.Nu().Vy(),true_neutrino.Nu().Vz())){
                        index = k;
                        min_E = mc_truth_eng[k];
                        }
                        }*/
				       
				       
                    }
                  }
			      
                  if(index != -1){
                    simb::MCNeutrino true_neutrino = mctruth_list[index]->GetNeutrino();
                    interaction = labels.GetInteractionType(true_neutrino);
                    labels.GetTopology(mctruth_list[index], fTopologyHitsCut);
                    float nu_energy = true_neutrino.Nu().E();
                    float lep_energy = true_neutrino.Lepton().E();
                    fNuID = true_neutrino.Nu().TrackId();
                    info.SetTruthInfo(nu_energy, lep_energy, 0., event_weight);
                    info.SetTopologyInformation(labels.GetPDG(), labels.GetNProtons(), labels.GetNPions(), labels.GetNPizeros(), labels.GetNNeutrons(), labels.GetTopologyType(), labels.GetTopologyTypeAlt());
                    //fNuSlice_purity = double(Get_Tot_nu_E(clockData,mctruth_list[index],slice_hits))/tot_slice_eng; // newly deleted
                    fNuSlice_purity = double(mc_truth_eng[index])/tot_slice_eng; // newly added
                    //if(Get_Tot_nu_E(clockData,mctruth_list[index],HitList)>0)fNuSlice_completeness = double(Get_Tot_nu_E(clockData,mctruth_list[index],slice_hits))/Get_Tot_nu_E(clockData,mctruth_list[index],HitList); // newly deleted
                    if(TotE_from_mctruth(clockData, HitList, mctruth_list[index])) fNuSlice_completeness = double(mc_truth_eng[index])/TotE_from_mctruth(clockData, HitList, mctruth_list[index]); // newly added
                    std::cout << "Nu purity : " << fNuSlice_purity << " Nu completeness : " << fNuSlice_completeness << "\n";
                  }
                } // select nutrino
              } // total slice energy is greater that 0
            } // valid hit association
		     
            if(findManyPFPs.isValid()){
              fIsSliceNu=Is_Slice_Nu(findManyPFPs,slice);
              if(fIsSliceNu){
                if(fm_pfpmd.isValid()){
                  fSliceScore = Get_Slice_Score(fm_pfpmd,Get_Nu_like_PFP(findManyPFPs,slice));
                  if(slice->ID() == Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,SliceList)[0]) fIsbestSlice = true;
                  fbestpfppdg = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,SliceList)[1];
                  fbestsliceID = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,SliceList)[0];
                  fpfppdg = Get_Slice_PFP_ID(findManyPFPs, slice);
                }
              }
            }
		     
            LArTrainingNuData train(interaction, *pixelmaps[i], info);
            std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_sl"+std::to_string(slice->ID())+"_h"+std::to_string(time(0));
            write_files(train, evtid);
            break;
            //} // found the matching slice
          } // find associated slice
        } // loop over pixel maps
      } // use backtrack inoformation
	  
	  
      else{
        std::cout << "***************** Used truth level information to get neutrino information ****************\n";
        for(unsigned int i=0; i<pixelmaps.size(); i++){
          Clear();
          frun = evt.run();
          fsubrun = evt.subRun();
          fevent = evt.id().event();
          if(Nu_is_contained) fNuIsContained = true;
		      
          LArTrainingNuData train(interaction, *pixelmaps[i], info);
          if (findOneSlice.isValid()){
            auto const & slice = findOneSlice.at(i);
            
            if(findManyHits.isValid()){
              std::vector<art::Ptr<recob::Hit>> slice_hits = findManyHits.at(slice.key());
              fNhits_tpc_0_pl_0 = Get_Hit_Count(0,0,slice_hits);
              fNhits_tpc_0_pl_1 = Get_Hit_Count(0,1,slice_hits);
              fNhits_tpc_0_pl_2 = Get_Hit_Count(0,2,slice_hits);
              fNhits_tpc_1_pl_0 = Get_Hit_Count(1,0,slice_hits);
              fNhits_tpc_1_pl_1 = Get_Hit_Count(1,1,slice_hits);
              fNhits_tpc_1_pl_2 = Get_Hit_Count(1,2,slice_hits);
              fNhits_total = slice_hits.size();
            }      
			      
            if(findManyPFPs.isValid()){
              fIsSliceNu=Is_Slice_Nu(findManyPFPs,slice);
			    
              ////////////////////////////////////////// pandora T0 //////////////////////////////////////////////////////
			    
              std::vector<float> pfp_T0_vec;
              if(findManyPFPs.isValid()){
                std::vector<art::Ptr<recob::PFParticle>> slicePFPs = findManyPFPs.at(slice.key());
                if(slicePFPs.size()){
                  for(auto const &pfp : slicePFPs){
                    if(findManyT0s.isValid()){
                      std::vector<art::Ptr<anab::T0>> T0_vec = findManyT0s.at(pfp.key());
                      if(T0_vec.size()){
                        for(auto const& T0 : T0_vec){
                          pfp_T0_vec.push_back(T0->Time());
                        }
                      }
                    }
                  }
                }
              }
			      
              if(pfp_T0_vec.size()){
                fT0 = *min_element(pfp_T0_vec.begin(), pfp_T0_vec.end());
              } 
			    
              ///////////////////////////////////////////////////////////////////////////////////////////////////////////
			    
              if(fIsSliceNu){			   
                if(fm_pfpmd.isValid()){
                  fSliceScore = Get_Slice_Score(fm_pfpmd,Get_Nu_like_PFP(findManyPFPs,slice));
                  if(slice->ID() == Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,SliceList)[0]) fIsbestSlice = true;
                  fbestpfppdg = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,SliceList)[1];
                  fbestsliceID = Get_Best_Slice_ID_PFP_pdg(findManyPFPs,fm_pfpmd,SliceList)[0];
                  fpfppdg = Get_Slice_PFP_ID(findManyPFPs, slice);
                }
              }
            }
            break;
            //}
            std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_sl"+std::to_string(slice->ID())+"_h"+std::to_string(time(0));
            write_files(train, evtid);
          }
		  
        }
      } // don't use slice information to extrac truth info
    } // use slcie to make images
	
    else{
      std::cout << "***************** Used entire event ****************\n";
      LArTrainingNuData train(interaction, *pixelmaps[0], info);
      std::string evtid = "r"+std::to_string(evt.run())+"_s"+std::to_string(evt.subRun())+"_e"+std::to_string(evt.event())+"_h"+std::to_string(time(0));
      write_files(train, evtid);
    } // use event to make images
  }

  //......................................................................
  
  void SBNDCVNZlibMaker::write_files(LArTrainingNuData td, std::string evtid)
  {
    // cropped from 2880 x 500 to 500 x 500 here 
    std::vector<unsigned char> pixel_array(3 * fPlaneLimit * fTDCLimit);
    //fImage.DisableRegionSelection();
    fImage.SetPixelMapSize(td.fPMap.NWire(), td.fPMap.NTdc());
    fImage.ConvertPixelMapToPixelArray(td.fPMap, pixel_array);

    ulong src_len = 3 * fPlaneLimit * fTDCLimit; // pixelArray length
    ulong dest_len = compressBound(src_len);     // calculate size of the compressed data               
    char* ostream = (char *) malloc(dest_len);  // allocate memory for the compressed data

    int res = compress((Bytef *) ostream, &dest_len, (Bytef *) &pixel_array[0], src_len);
    
    // Buffer error

    if (res == Z_BUF_ERROR)
      std::cout << "Buffer too small!" << std::endl;

    // Memory error
    else if (res ==  Z_MEM_ERROR)
      std::cout << "Not enough memory for compression!" << std::endl;

    // Compression ok 
    else {
      // Create output files 
      std::string image_file_name = out_dir + "/event_" + evtid + ".gz";
      std::string info_file_name = out_dir + "/event_" +  evtid + ".info";

      std::ofstream image_file (image_file_name, std::ofstream::binary);
      std::ofstream info_file  (info_file_name);

      if(image_file.is_open() && info_file.is_open()) {

        // Write compressed data to file

        image_file.write(ostream, dest_len);

        image_file.close(); // close file

        // Write records to file

        // Category
	
	info_file << td.fInt << std::endl;
        //info_file << td.fInfo << std::endl;
	info_file << td.fInfo;
        info_file << td.fPMap.GetTotHits() << "," << fNhits_tpc_0_pl_0 << "," << fNhits_tpc_0_pl_1 << "," << fNhits_tpc_0_pl_2 << "," << fNhits_tpc_1_pl_0 << "," << fNhits_tpc_1_pl_1 << "," << fNhits_tpc_1_pl_2 << "," << fNhits_total <<  std::endl;   
	info_file << frun << "," << fsubrun << "," << fevent << "," << fsliceID << std::endl;   
	info_file << ftotsliceE << "," << ftotsliceNuE << "," << ftotsliceCosE << "," << ftotsliceOthE << "," << fsliceNuEfrac << "," << fsliceCosEfrac << "," << fsliceOthEfrac << "," << fNuSlice_purity << "," << fNuSlice_completeness <<  std::endl;
	info_file << fIsSliceNu << "," << fSliceScore << "," << fIsbestSlice << "," << fbestpfppdg << "," << fbestsliceID << "," << fpfppdg << std::endl;
	info_file << fNuIsContained << "," << fNuID << "  " << fT0 << std::endl;
	info_file.close(); // close file
      }
      else {

        if (image_file.is_open())
          image_file.close();
        else 
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << image_file_name << "!" << std::endl;

        if (info_file.is_open())
          info_file.close();
        else
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << info_file_name << "!" << std::endl;
      }
    }
    
    free(ostream);  // free allocated memory

  } // lcvn::CVNZlibMaker::write_files
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool SBNDCVNZlibMaker::Is_in_Fiducial_Vol(double stX, double stY, double stZ)
  {
    double X_bound = 196.5 - fVolCut;
    double Y_bound = 200. - fVolCut;
    double Z_up_bound = 500. - fVolCut;
    double Z_low_bound = fVolCut;
	
    if(TMath::Abs(stX) <= X_bound && TMath::Abs(stY) <= Y_bound && (stZ >=Z_low_bound && stZ <= Z_up_bound)) return true;
    else return false;
  }  

  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  void  SBNDCVNZlibMaker::HitTruth( detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& truthid){
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    if( !trackIDEs.size() ) return;
    Float_t maxe = 0;
    Int_t bestid = 0;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      if( trackIDEs[i].energy > maxe ) {
        maxe = trackIDEs[i].energy;
        bestid = trackIDEs[i].trackID;
      }
    }
    truthid = bestid;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool  SBNDCVNZlibMaker::HitTruthId(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t& mcid){
    mcid = std::numeric_limits<int>::lowest();
    HitTruth(clockData,hit,mcid);
    if( mcid > std::numeric_limits<int>::lowest() ) return true;
    else return false;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool  SBNDCVNZlibMaker::TrackIdToMCTruth(Int_t const trkID, art::Ptr<simb::MCTruth>& mctruth){
    bool matchFound = false;
    try {
      mctruth = pi_serv->TrackIdToMCTruth_P(trkID);
      matchFound = true;
    } catch(...) {
      std::cout<<"Exception thrown matching TrackID "<<trkID<<" to MCTruth\n";
      matchFound = false;
    }
    return matchFound;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::HitEfrmTrkID(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, Int_t truthid){
    double hit_E = 0.;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    if( !trackIDEs.size() ) return hit_E;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      if(trackIDEs[i].trackID == truthid) hit_E += trackIDEs[i].energy;
    }
    return hit_E;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::HitTotE(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit){
    double Tot_E = 0.;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    if( !trackIDEs.size() ) return Tot_E;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      Tot_E += trackIDEs[i].energy;
    }
    return Tot_E;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool SBNDCVNZlibMaker::Is_Slice_Nu(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce){
    if(assn.isValid()){
      std::vector<art::Ptr<recob::PFParticle>> slicePFPs = assn.at(slce.key());
      for(auto const &pfp : slicePFPs){
        if((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)) return true;
      }
    }
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::Get_Slice_Score(art::FindManyP<larpandoraobj::PFParticleMetadata> const& assn, art::Ptr<recob::PFParticle> const& pfp){
    double best_score = -9999;
    if(assn.isValid() && pfp.isNonnull()){
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = assn.at(pfp->Self());
      if(pfpMetaVec.size()){
        for (auto const pfpMeta : pfpMetaVec){
          larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
          double score = propertiesMap.at("NuScore");
          if (score > best_score) best_score = score;
        }
      }
    }
    return best_score;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  art::Ptr<recob::PFParticle> SBNDCVNZlibMaker::Get_Nu_like_PFP(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce){
    art::Ptr<recob::PFParticle> result;
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs = assn.at(slce.key());
    for(auto const &pfp : slicePFPs){
      if((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){ 
        result = pfp;
        break;
      }
    }
    return result;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<int> SBNDCVNZlibMaker::Get_Best_Slice_ID_PFP_pdg(art::FindManyP<recob::PFParticle> const& pfp_assn, art::FindManyP<larpandoraobj::PFParticleMetadata> const& metadata_assn,
                                                               std::vector<art::Ptr<recob::Slice>> const& slice_vec){
    std::vector<int> Sl_ID_pfp_pdf = {-999,-999};
    std::vector<art::Ptr<recob::PFParticle>> nuPfpVec;
    std::vector<art::Ptr<recob::Slice>> nuSliceVec;
    for(auto const &slice : slice_vec){
      std::vector<art::Ptr<recob::PFParticle>> slicePFPs = pfp_assn.at(slice.key());
      for (auto const &pfp : slicePFPs){
        if ((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){
          nuPfpVec.push_back(pfp);
          nuSliceVec.push_back(slice);
          break;
        }
      }
    }
    
    if(nuPfpVec.size()){
      float bestScore = -999.;
      art::Ptr<recob::Slice> bestNuSlice;
      art::Ptr<recob::PFParticle> bestNuPfp;
      for(size_t i=0; i<nuPfpVec.size(); i++){
        float score = -999.;
        const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMetaVec = metadata_assn.at(nuPfpVec.at(i)->Self());
        for (auto const pfpMeta : pfpMetaVec) {
          larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
          score = propertiesMap.at("NuScore");
          if( score > bestScore ) {
            bestScore = score;
            bestNuPfp = nuPfpVec.at(i);
            bestNuSlice = nuSliceVec.at(i);
          }
        }
      }
       
      Sl_ID_pfp_pdf[0] = bestNuSlice->ID();
      Sl_ID_pfp_pdf[1] = bestNuPfp->PdgCode();
    }
    return Sl_ID_pfp_pdf;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int SBNDCVNZlibMaker::Get_Slice_PFP_ID(art::FindManyP<recob::PFParticle> const& assn, art::Ptr<recob::Slice> const& slce){
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs = assn.at(slce.key());
    int result = -9999;
    for(auto const &pfp : slicePFPs){
      if((pfp->PdgCode() == 12 || pfp->PdgCode() == 14)){ 
        result = pfp->PdgCode();
        break;
      }
    }
    return result;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int SBNDCVNZlibMaker::Get_Hit_Count(unsigned int tpc_no, unsigned int pln_no, std::vector<art::Ptr<recob::Hit>> const& hit_vec){
    int n_hits = 0;
    if(hit_vec.size()){
      for(auto const hit : hit_vec){
        if(hit->WireID().TPC == tpc_no){
          if(hit->WireID().Plane == pln_no){
            n_hits++;
          }
        }
      }
    }
    
    return n_hits;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::Get_Tot_nu_E(detinfo::DetectorClocksData const& clockData,art::Ptr<simb::MCTruth> const& mcneutrino, std::vector<art::Ptr<recob::Hit>> const& hit_vec){
    double tot_nu_E = 0;
    for(auto const hit : hit_vec){
      Int_t trkId;
      if(HitTruthId(clockData,hit,trkId)){
        art::Ptr<simb::MCTruth> mctruth;
        if(TrackIdToMCTruth(trkId,mctruth)){
          if(mctruth->Origin() == simb::kBeamNeutrino){
            if(mctruth == mcneutrino){
              tot_nu_E += HitEfrmTrkID(clockData,hit,trkId);
            }
          }
        }
      }
    }
    return tot_nu_E;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::HitTotE_frm_given_orgin(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, std::string origin_name){
    double Tot_E = 0.;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    if( !trackIDEs.size() ) return Tot_E;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      art::Ptr<simb::MCTruth> mctruth;
      if(TrackIdToMCTruth(trackIDEs[i].trackID,mctruth)){
        if(origin_name == "nu"){
          if(mctruth->Origin() == simb::kBeamNeutrino) Tot_E += trackIDEs[i].energy;
        }
        else if(origin_name == "cos"){
          if(mctruth->Origin() == simb::kCosmicRay) Tot_E += trackIDEs[i].energy;
        }
      }
           
    }
    return Tot_E;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::HitTotE_frm_mctruth(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit> const& hit, art::Ptr<simb::MCTruth> const& my_mctruth){
    double Tot_E = 0.;
    std::vector<sim::TrackIDE> trackIDEs = bt_serv->HitToTrackIDEs(clockData, hit);
    if( !trackIDEs.size() ) return Tot_E;
    for(size_t i = 0; i < trackIDEs.size(); ++i){
      art::Ptr<simb::MCTruth> mctruth;
      if(TrackIdToMCTruth(trackIDEs[i].trackID,mctruth)){
        if(my_mctruth.isNonnull()){
          if(mctruth == my_mctruth){
            Tot_E += trackIDEs[i].energy;
          }
        }
      }
    }
    return Tot_E;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double SBNDCVNZlibMaker::TotE_from_mctruth(detinfo::DetectorClocksData const& clockData, std::vector<art::Ptr<recob::Hit>> const& hit_vec, art::Ptr<simb::MCTruth> const& my_mctruth){
    double tot_E = 0;
    if(hit_vec.size() && my_mctruth.isNonnull()){
      for(auto const hit : hit_vec){
        tot_E += HitTotE_frm_mctruth(clockData,hit,my_mctruth);
      }
    }
    return tot_E;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void SBNDCVNZlibMaker::Clear(){
    frun = -9999;
    fsubrun = -9999;
    fevent = -9999;
    fsliceID = -9999;
    ftotsliceE = -9999;
    ftotsliceNuE = -9999;
    ftotsliceCosE = -9999;
    ftotsliceOthE = -9999;
    fsliceNuEfrac = -9999;
    fsliceCosEfrac = -9999;
    fsliceOthEfrac = -9999;	    
    fIsSliceNu = false;
    fSliceScore = -9999;
    fIsbestSlice = false;
    fbestpfppdg = -9999;
    fbestsliceID = -9999;
    fpfppdg = -9999;
    fNuIsContained = false;
    fNuID = -9999;
    fNhits_tpc_0_pl_0 = -9999;
    fNhits_tpc_0_pl_1 = -9999;
    fNhits_tpc_0_pl_2 = -9999;
    fNhits_tpc_1_pl_0 = -9999;
    fNhits_tpc_1_pl_1 = -9999;
    fNhits_tpc_1_pl_2 = -9999;
    fNhits_total = -9999;
    fNuSlice_purity = -9999.;
    fNuSlice_completeness = -9999.;
    fT0 = 0.;
  }

  DEFINE_ART_MODULE(SBNDCVNZlibMaker)
} // namespace cvn
