//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  7 10:22:34 2024 by ROOT version 6.26/07
// from TTree tree/
// found on file: /exp/sbnd/data/users/hlay/crt/clustering/merge_checks_jan2024/production/ana/crtana_sbnd.root
//////////////////////////////////////////////////////////

#ifndef tree_h
#define tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   vector<short>   *mc_trackid;
   vector<short>   *mc_pdg;
   vector<short>   *mc_status;
   vector<unsigned short> *mc_ndaughters;
   vector<double>  *mc_vx;
   vector<double>  *mc_vy;
   vector<double>  *mc_vz;
   vector<double>  *mc_vt;
   vector<double>  *mc_endx;
   vector<double>  *mc_endy;
   vector<double>  *mc_endz;
   vector<double>  *mc_endt;
   vector<double>  *mc_vpx;
   vector<double>  *mc_vpy;
   vector<double>  *mc_vpz;
   vector<double>  *mc_ve;
   vector<double>  *mc_endpx;
   vector<double>  *mc_endpy;
   vector<double>  *mc_endpz;
   vector<double>  *mc_ende;
   vector<short>   *ide_trackid;
   vector<float>   *ide_e;
   vector<float>   *ide_entryx;
   vector<float>   *ide_entryy;
   vector<float>   *ide_entryz;
   vector<float>   *ide_entryt;
   vector<float>   *ide_exitx;
   vector<float>   *ide_exity;
   vector<float>   *ide_exitz;
   vector<float>   *ide_exitt;
   vector<unsigned short> *feb_mac5;
   vector<unsigned short> *feb_flags;
   vector<unsigned int> *feb_ts0;
   vector<unsigned int> *feb_ts1;
   vector<unsigned int> *feb_unixs;
   vector<vector<unsigned short> > *feb_adc;
   vector<unsigned int> *feb_coinc;
   vector<unsigned int> *sh_channel;
   vector<unsigned int> *sh_ts0;
   vector<unsigned int> *sh_ts1;
   vector<unsigned int> *sh_unixs;
   vector<double>  *sh_pos;
   vector<double>  *sh_err;
   vector<unsigned short> *sh_adc1;
   vector<unsigned short> *sh_adc2;
   vector<bool>    *sh_saturated1;
   vector<bool>    *sh_saturated2;
   vector<int>     *sh_truth_trackid;
   vector<double>  *sh_truth_completeness;
   vector<double>  *sh_truth_purity;
   vector<double>  *sh_truth_pos;
   vector<double>  *sh_truth_energy;
   vector<double>  *sh_truth_time;
   vector<unsigned int> *cl_ts0;
   vector<unsigned int> *cl_ts1;
   vector<unsigned int> *cl_unixs;
   vector<unsigned short> *cl_nhits;
   vector<short>   *cl_tagger;
   vector<unsigned char> *cl_composition;
   vector<int>     *cl_truth_trackid;
   vector<double>  *cl_truth_completeness;
   vector<double>  *cl_truth_purity;
   vector<double>  *cl_truth_hit_completeness;
   vector<double>  *cl_truth_hit_purity;
   vector<int>     *cl_truth_pdg;
   vector<double>  *cl_truth_energy;
   vector<double>  *cl_truth_time;
   vector<double>  *cl_truth_x;
   vector<double>  *cl_truth_y;
   vector<double>  *cl_truth_z;
   vector<double>  *cl_truth_core_energy;
   vector<double>  *cl_truth_core_time;
   vector<double>  *cl_truth_core_x;
   vector<double>  *cl_truth_core_y;
   vector<double>  *cl_truth_core_z;
   vector<bool>    *cl_has_sp;
   vector<double>  *cl_sp_x;
   vector<double>  *cl_sp_ex;
   vector<double>  *cl_sp_y;
   vector<double>  *cl_sp_ey;
   vector<double>  *cl_sp_z;
   vector<double>  *cl_sp_ez;
   vector<double>  *cl_sp_pe;
   vector<double>  *cl_sp_time;
   vector<double>  *cl_sp_etime;
   vector<bool>    *cl_sp_complete;
   vector<int>     *td_tag_trackid;
   vector<int>     *td_tag_pdg;
   vector<short>   *td_tag_tagger;
   vector<double>  *td_tag_energy;
   vector<double>  *td_tag_time;
   vector<double>  *td_tag_x;
   vector<double>  *td_tag_y;
   vector<double>  *td_tag_z;
   vector<bool>    *td_tag_reco_status;
   vector<int>     *td_trackid;
   vector<int>     *td_pdg;
   vector<double>  *td_energy;
   vector<double>  *td_time;
   vector<bool>    *td_reconstructable;
   vector<bool>    *td_reco_status;
   vector<bool>    *td_reco_triple;
   vector<double>  *tr_start_x;
   vector<double>  *tr_start_y;
   vector<double>  *tr_start_z;
   vector<double>  *tr_end_x;
   vector<double>  *tr_end_y;
   vector<double>  *tr_end_z;
   vector<double>  *tr_dir_x;
   vector<double>  *tr_dir_y;
   vector<double>  *tr_dir_z;
   vector<double>  *tr_time;
   vector<double>  *tr_etime;
   vector<double>  *tr_pe;
   vector<double>  *tr_length;
   vector<double>  *tr_tof;
   vector<double>  *tr_theta;
   vector<double>  *tr_phi;
   vector<bool>    *tr_triple;
   vector<int>     *tr_truth_trackid;
   vector<double>  *tr_truth_completeness;
   vector<double>  *tr_truth_purity;
   vector<int>     *tr_truth_pdg;
   vector<double>  *tr_truth_energy;
   vector<double>  *tr_truth_time;
   vector<double>  *tr_truth_start_x;
   vector<double>  *tr_truth_start_y;
   vector<double>  *tr_truth_start_z;
   vector<double>  *tr_truth_end_x;
   vector<double>  *tr_truth_end_y;
   vector<double>  *tr_truth_end_z;
   vector<double>  *tr_truth_dir_x;
   vector<double>  *tr_truth_dir_y;
   vector<double>  *tr_truth_dir_z;
   vector<double>  *tr_truth_particle_energy;
   vector<double>  *tr_truth_length;
   vector<double>  *tr_truth_tof;
   vector<double>  *tr_truth_theta;
   vector<double>  *tr_truth_phi;
   vector<double>  *tpc_start_x;
   vector<double>  *tpc_start_y;
   vector<double>  *tpc_start_z;
   vector<double>  *tpc_end_x;
   vector<double>  *tpc_end_y;
   vector<double>  *tpc_end_z;
   vector<double>  *tpc_dir_x;
   vector<double>  *tpc_dir_y;
   vector<double>  *tpc_dir_z;
   vector<double>  *tpc_length;
   vector<double>  *tpc_track_score;
   vector<int>     *tpc_truth_trackid;
   vector<int>     *tpc_truth_pdg;
   vector<double>  *tpc_truth_energy;
   vector<double>  *tpc_truth_time;
   vector<bool>    *tpc_sp_matchable;
   vector<bool>    *tpc_sp_matched;
   vector<bool>    *tpc_sp_good_match;
   vector<double>  *tpc_sp_time;
   vector<double>  *tpc_sp_score;
   vector<bool>    *tpc_tr_matchable;
   vector<bool>    *tpc_tr_matched;
   vector<bool>    *tpc_tr_good_match;
   vector<double>  *tpc_tr_time;
   vector<double>  *tpc_tr_score;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_mc_trackid;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_status;   //!
   TBranch        *b_mc_ndaughters;   //!
   TBranch        *b_mc_vx;   //!
   TBranch        *b_mc_vy;   //!
   TBranch        *b_mc_vz;   //!
   TBranch        *b_mc_vt;   //!
   TBranch        *b_mc_endx;   //!
   TBranch        *b_mc_endy;   //!
   TBranch        *b_mc_endz;   //!
   TBranch        *b_mc_endt;   //!
   TBranch        *b_mc_vpx;   //!
   TBranch        *b_mc_vpy;   //!
   TBranch        *b_mc_vpz;   //!
   TBranch        *b_mc_ve;   //!
   TBranch        *b_mc_endpx;   //!
   TBranch        *b_mc_endpy;   //!
   TBranch        *b_mc_endpz;   //!
   TBranch        *b_mc_ende;   //!
   TBranch        *b_ide_trackid;   //!
   TBranch        *b_ide_e;   //!
   TBranch        *b_ide_entryx;   //!
   TBranch        *b_ide_entryy;   //!
   TBranch        *b_ide_entryz;   //!
   TBranch        *b_ide_entryt;   //!
   TBranch        *b_ide_exitx;   //!
   TBranch        *b_ide_exity;   //!
   TBranch        *b_ide_exitz;   //!
   TBranch        *b_ide_exitt;   //!
   TBranch        *b_feb_mac5;   //!
   TBranch        *b_feb_flags;   //!
   TBranch        *b_feb_ts0;   //!
   TBranch        *b_feb_ts1;   //!
   TBranch        *b_feb_unixs;   //!
   TBranch        *b_feb_adc;   //!
   TBranch        *b_feb_coinc;   //!
   TBranch        *b_sh_channel;   //!
   TBranch        *b_sh_ts0;   //!
   TBranch        *b_sh_ts1;   //!
   TBranch        *b_sh_unixs;   //!
   TBranch        *b_sh_pos;   //!
   TBranch        *b_sh_err;   //!
   TBranch        *b_sh_adc1;   //!
   TBranch        *b_sh_adc2;   //!
   TBranch        *b_sh_saturated1;   //!
   TBranch        *b_sh_saturated2;   //!
   TBranch        *b_sh_truth_trackid;   //!
   TBranch        *b_sh_truth_completeness;   //!
   TBranch        *b_sh_truth_purity;   //!
   TBranch        *b_sh_truth_pos;   //!
   TBranch        *b_sh_truth_energy;   //!
   TBranch        *b_sh_truth_time;   //!
   TBranch        *b_cl_ts0;   //!
   TBranch        *b_cl_ts1;   //!
   TBranch        *b_cl_unixs;   //!
   TBranch        *b_cl_nhits;   //!
   TBranch        *b_cl_tagger;   //!
   TBranch        *b_cl_composition;   //!
   TBranch        *b_cl_truth_trackid;   //!
   TBranch        *b_cl_truth_completeness;   //!
   TBranch        *b_cl_truth_purity;   //!
   TBranch        *b_cl_truth_hit_completeness;   //!
   TBranch        *b_cl_truth_hit_purity;   //!
   TBranch        *b_cl_truth_pdg;   //!
   TBranch        *b_cl_truth_energy;   //!
   TBranch        *b_cl_truth_time;   //!
   TBranch        *b_cl_truth_x;   //!
   TBranch        *b_cl_truth_y;   //!
   TBranch        *b_cl_truth_z;   //!
   TBranch        *b_cl_truth_core_energy;   //!
   TBranch        *b_cl_truth_core_time;   //!
   TBranch        *b_cl_truth_core_x;   //!
   TBranch        *b_cl_truth_core_y;   //!
   TBranch        *b_cl_truth_core_z;   //!
   TBranch        *b_cl_has_sp;   //!
   TBranch        *b_cl_sp_x;   //!
   TBranch        *b_cl_sp_ex;   //!
   TBranch        *b_cl_sp_y;   //!
   TBranch        *b_cl_sp_ey;   //!
   TBranch        *b_cl_sp_z;   //!
   TBranch        *b_cl_sp_ez;   //!
   TBranch        *b_cl_sp_pe;   //!
   TBranch        *b_cl_sp_time;   //!
   TBranch        *b_cl_sp_etime;   //!
   TBranch        *b_cl_sp_complete;   //!
   TBranch        *b_td_tag_trackid;   //!
   TBranch        *b_td_tag_pdg;   //!
   TBranch        *b_td_tag_tagger;   //!
   TBranch        *b_td_tag_energy;   //!
   TBranch        *b_td_tag_time;   //!
   TBranch        *b_td_tag_x;   //!
   TBranch        *b_td_tag_y;   //!
   TBranch        *b_td_tag_z;   //!
   TBranch        *b_td_tag_reco_status;   //!
   TBranch        *b_td_trackid;   //!
   TBranch        *b_td_pdg;   //!
   TBranch        *b_td_energy;   //!
   TBranch        *b_td_time;   //!
   TBranch        *b_td_reconstructable;   //!
   TBranch        *b_td_reco_status;   //!
   TBranch        *b_td_reco_triple;   //!
   TBranch        *b_tr_start_x;   //!
   TBranch        *b_tr_start_y;   //!
   TBranch        *b_tr_start_z;   //!
   TBranch        *b_tr_end_x;   //!
   TBranch        *b_tr_end_y;   //!
   TBranch        *b_tr_end_z;   //!
   TBranch        *b_tr_dir_x;   //!
   TBranch        *b_tr_dir_y;   //!
   TBranch        *b_tr_dir_z;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_tr_etime;   //!
   TBranch        *b_tr_pe;   //!
   TBranch        *b_tr_length;   //!
   TBranch        *b_tr_tof;   //!
   TBranch        *b_tr_theta;   //!
   TBranch        *b_tr_phi;   //!
   TBranch        *b_tr_triple;   //!
   TBranch        *b_tr_truth_trackid;   //!
   TBranch        *b_tr_truth_completeness;   //!
   TBranch        *b_tr_truth_purity;   //!
   TBranch        *b_tr_truth_pdg;   //!
   TBranch        *b_tr_truth_energy;   //!
   TBranch        *b_tr_truth_time;   //!
   TBranch        *b_tr_truth_start_x;   //!
   TBranch        *b_tr_truth_start_y;   //!
   TBranch        *b_tr_truth_start_z;   //!
   TBranch        *b_tr_truth_end_x;   //!
   TBranch        *b_tr_truth_end_y;   //!
   TBranch        *b_tr_truth_end_z;   //!
   TBranch        *b_tr_truth_dir_x;   //!
   TBranch        *b_tr_truth_dir_y;   //!
   TBranch        *b_tr_truth_dir_z;   //!
   TBranch        *b_tr_truth_particle_energy;   //!
   TBranch        *b_tr_truth_length;   //!
   TBranch        *b_tr_truth_tof;   //!
   TBranch        *b_tr_truth_theta;   //!
   TBranch        *b_tr_truth_phi;   //!
   TBranch        *b_tpc_start_x;   //!
   TBranch        *b_tpc_start_y;   //!
   TBranch        *b_tpc_start_z;   //!
   TBranch        *b_tpc_end_x;   //!
   TBranch        *b_tpc_end_y;   //!
   TBranch        *b_tpc_end_z;   //!
   TBranch        *b_tpc_dir_x;   //!
   TBranch        *b_tpc_dir_y;   //!
   TBranch        *b_tpc_dir_z;   //!
   TBranch        *b_tpc_length;   //!
   TBranch        *b_tpc_track_score;   //!
   TBranch        *b_tpc_truth_trackid;   //!
   TBranch        *b_tpc_truth_pdg;   //!
   TBranch        *b_tpc_truth_energy;   //!
   TBranch        *b_tpc_truth_time;   //!
   TBranch        *b_tpc_sp_matchable;   //!
   TBranch        *b_tpc_sp_matched;   //!
   TBranch        *b_tpc_sp_good_match;   //!
   TBranch        *b_tpc_sp_time;   //!
   TBranch        *b_tpc_sp_score;   //!
   TBranch        *b_tpc_tr_matchable;   //!
   TBranch        *b_tpc_tr_matched;   //!
   TBranch        *b_tpc_tr_good_match;   //!
   TBranch        *b_tpc_tr_time;   //!
   TBranch        *b_tpc_tr_score;   //!

   tree(TTree *tree=0);
   virtual ~tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tree_cxx
tree::tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/exp/sbnd/data/users/hlay/crt/clustering/merge_checks_jan2024/production/ana/crtana_sbnd.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/exp/sbnd/data/users/hlay/crt/clustering/merge_checks_jan2024/production/ana/crtana_sbnd.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/exp/sbnd/data/users/hlay/crt/clustering/merge_checks_jan2024/production/ana/crtana_sbnd.root:/crtana");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

tree::~tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mc_trackid = 0;
   mc_pdg = 0;
   mc_status = 0;
   mc_ndaughters = 0;
   mc_vx = 0;
   mc_vy = 0;
   mc_vz = 0;
   mc_vt = 0;
   mc_endx = 0;
   mc_endy = 0;
   mc_endz = 0;
   mc_endt = 0;
   mc_vpx = 0;
   mc_vpy = 0;
   mc_vpz = 0;
   mc_ve = 0;
   mc_endpx = 0;
   mc_endpy = 0;
   mc_endpz = 0;
   mc_ende = 0;
   ide_trackid = 0;
   ide_e = 0;
   ide_entryx = 0;
   ide_entryy = 0;
   ide_entryz = 0;
   ide_entryt = 0;
   ide_exitx = 0;
   ide_exity = 0;
   ide_exitz = 0;
   ide_exitt = 0;
   feb_mac5 = 0;
   feb_flags = 0;
   feb_ts0 = 0;
   feb_ts1 = 0;
   feb_unixs = 0;
   feb_adc = 0;
   feb_coinc = 0;
   sh_channel = 0;
   sh_ts0 = 0;
   sh_ts1 = 0;
   sh_unixs = 0;
   sh_pos = 0;
   sh_err = 0;
   sh_adc1 = 0;
   sh_adc2 = 0;
   sh_saturated1 = 0;
   sh_saturated2 = 0;
   sh_truth_trackid = 0;
   sh_truth_completeness = 0;
   sh_truth_purity = 0;
   sh_truth_pos = 0;
   sh_truth_energy = 0;
   sh_truth_time = 0;
   cl_ts0 = 0;
   cl_ts1 = 0;
   cl_unixs = 0;
   cl_nhits = 0;
   cl_tagger = 0;
   cl_composition = 0;
   cl_truth_trackid = 0;
   cl_truth_completeness = 0;
   cl_truth_purity = 0;
   cl_truth_hit_completeness = 0;
   cl_truth_hit_purity = 0;
   cl_truth_pdg = 0;
   cl_truth_energy = 0;
   cl_truth_time = 0;
   cl_truth_x = 0;
   cl_truth_y = 0;
   cl_truth_z = 0;
   cl_truth_core_energy = 0;
   cl_truth_core_time = 0;
   cl_truth_core_x = 0;
   cl_truth_core_y = 0;
   cl_truth_core_z = 0;
   cl_has_sp = 0;
   cl_sp_x = 0;
   cl_sp_ex = 0;
   cl_sp_y = 0;
   cl_sp_ey = 0;
   cl_sp_z = 0;
   cl_sp_ez = 0;
   cl_sp_pe = 0;
   cl_sp_time = 0;
   cl_sp_etime = 0;
   cl_sp_complete = 0;
   td_tag_trackid = 0;
   td_tag_pdg = 0;
   td_tag_tagger = 0;
   td_tag_energy = 0;
   td_tag_time = 0;
   td_tag_x = 0;
   td_tag_y = 0;
   td_tag_z = 0;
   td_tag_reco_status = 0;
   td_trackid = 0;
   td_pdg = 0;
   td_energy = 0;
   td_time = 0;
   td_reconstructable = 0;
   td_reco_status = 0;
   td_reco_triple = 0;
   tr_start_x = 0;
   tr_start_y = 0;
   tr_start_z = 0;
   tr_end_x = 0;
   tr_end_y = 0;
   tr_end_z = 0;
   tr_dir_x = 0;
   tr_dir_y = 0;
   tr_dir_z = 0;
   tr_time = 0;
   tr_etime = 0;
   tr_pe = 0;
   tr_length = 0;
   tr_tof = 0;
   tr_theta = 0;
   tr_phi = 0;
   tr_triple = 0;
   tr_truth_trackid = 0;
   tr_truth_completeness = 0;
   tr_truth_purity = 0;
   tr_truth_pdg = 0;
   tr_truth_energy = 0;
   tr_truth_time = 0;
   tr_truth_start_x = 0;
   tr_truth_start_y = 0;
   tr_truth_start_z = 0;
   tr_truth_end_x = 0;
   tr_truth_end_y = 0;
   tr_truth_end_z = 0;
   tr_truth_dir_x = 0;
   tr_truth_dir_y = 0;
   tr_truth_dir_z = 0;
   tr_truth_particle_energy = 0;
   tr_truth_length = 0;
   tr_truth_tof = 0;
   tr_truth_theta = 0;
   tr_truth_phi = 0;
   tpc_start_x = 0;
   tpc_start_y = 0;
   tpc_start_z = 0;
   tpc_end_x = 0;
   tpc_end_y = 0;
   tpc_end_z = 0;
   tpc_dir_x = 0;
   tpc_dir_y = 0;
   tpc_dir_z = 0;
   tpc_length = 0;
   tpc_track_score = 0;
   tpc_truth_trackid = 0;
   tpc_truth_pdg = 0;
   tpc_truth_energy = 0;
   tpc_truth_time = 0;
   tpc_sp_matchable = 0;
   tpc_sp_matched = 0;
   tpc_sp_good_match = 0;
   tpc_sp_time = 0;
   tpc_sp_score = 0;
   tpc_tr_matchable = 0;
   tpc_tr_matched = 0;
   tpc_tr_good_match = 0;
   tpc_tr_time = 0;
   tpc_tr_score = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("mc_trackid", &mc_trackid, &b_mc_trackid);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
   fChain->SetBranchAddress("mc_ndaughters", &mc_ndaughters, &b_mc_ndaughters);
   fChain->SetBranchAddress("mc_vx", &mc_vx, &b_mc_vx);
   fChain->SetBranchAddress("mc_vy", &mc_vy, &b_mc_vy);
   fChain->SetBranchAddress("mc_vz", &mc_vz, &b_mc_vz);
   fChain->SetBranchAddress("mc_vt", &mc_vt, &b_mc_vt);
   fChain->SetBranchAddress("mc_endx", &mc_endx, &b_mc_endx);
   fChain->SetBranchAddress("mc_endy", &mc_endy, &b_mc_endy);
   fChain->SetBranchAddress("mc_endz", &mc_endz, &b_mc_endz);
   fChain->SetBranchAddress("mc_endt", &mc_endt, &b_mc_endt);
   fChain->SetBranchAddress("mc_vpx", &mc_vpx, &b_mc_vpx);
   fChain->SetBranchAddress("mc_vpy", &mc_vpy, &b_mc_vpy);
   fChain->SetBranchAddress("mc_vpz", &mc_vpz, &b_mc_vpz);
   fChain->SetBranchAddress("mc_ve", &mc_ve, &b_mc_ve);
   fChain->SetBranchAddress("mc_endpx", &mc_endpx, &b_mc_endpx);
   fChain->SetBranchAddress("mc_endpy", &mc_endpy, &b_mc_endpy);
   fChain->SetBranchAddress("mc_endpz", &mc_endpz, &b_mc_endpz);
   fChain->SetBranchAddress("mc_ende", &mc_ende, &b_mc_ende);
   fChain->SetBranchAddress("ide_trackid", &ide_trackid, &b_ide_trackid);
   fChain->SetBranchAddress("ide_e", &ide_e, &b_ide_e);
   fChain->SetBranchAddress("ide_entryx", &ide_entryx, &b_ide_entryx);
   fChain->SetBranchAddress("ide_entryy", &ide_entryy, &b_ide_entryy);
   fChain->SetBranchAddress("ide_entryz", &ide_entryz, &b_ide_entryz);
   fChain->SetBranchAddress("ide_entryt", &ide_entryt, &b_ide_entryt);
   fChain->SetBranchAddress("ide_exitx", &ide_exitx, &b_ide_exitx);
   fChain->SetBranchAddress("ide_exity", &ide_exity, &b_ide_exity);
   fChain->SetBranchAddress("ide_exitz", &ide_exitz, &b_ide_exitz);
   fChain->SetBranchAddress("ide_exitt", &ide_exitt, &b_ide_exitt);
   fChain->SetBranchAddress("feb_mac5", &feb_mac5, &b_feb_mac5);
   fChain->SetBranchAddress("feb_flags", &feb_flags, &b_feb_flags);
   fChain->SetBranchAddress("feb_ts0", &feb_ts0, &b_feb_ts0);
   fChain->SetBranchAddress("feb_ts1", &feb_ts1, &b_feb_ts1);
   fChain->SetBranchAddress("feb_unixs", &feb_unixs, &b_feb_unixs);
   fChain->SetBranchAddress("feb_adc", &feb_adc, &b_feb_adc);
   fChain->SetBranchAddress("feb_coinc", &feb_coinc, &b_feb_coinc);
   fChain->SetBranchAddress("sh_channel", &sh_channel, &b_sh_channel);
   fChain->SetBranchAddress("sh_ts0", &sh_ts0, &b_sh_ts0);
   fChain->SetBranchAddress("sh_ts1", &sh_ts1, &b_sh_ts1);
   fChain->SetBranchAddress("sh_unixs", &sh_unixs, &b_sh_unixs);
   fChain->SetBranchAddress("sh_pos", &sh_pos, &b_sh_pos);
   fChain->SetBranchAddress("sh_err", &sh_err, &b_sh_err);
   fChain->SetBranchAddress("sh_adc1", &sh_adc1, &b_sh_adc1);
   fChain->SetBranchAddress("sh_adc2", &sh_adc2, &b_sh_adc2);
   fChain->SetBranchAddress("sh_saturated1", &sh_saturated1, &b_sh_saturated1);
   fChain->SetBranchAddress("sh_saturated2", &sh_saturated2, &b_sh_saturated2);
   fChain->SetBranchAddress("sh_truth_trackid", &sh_truth_trackid, &b_sh_truth_trackid);
   fChain->SetBranchAddress("sh_truth_completeness", &sh_truth_completeness, &b_sh_truth_completeness);
   fChain->SetBranchAddress("sh_truth_purity", &sh_truth_purity, &b_sh_truth_purity);
   fChain->SetBranchAddress("sh_truth_pos", &sh_truth_pos, &b_sh_truth_pos);
   fChain->SetBranchAddress("sh_truth_energy", &sh_truth_energy, &b_sh_truth_energy);
   fChain->SetBranchAddress("sh_truth_time", &sh_truth_time, &b_sh_truth_time);
   fChain->SetBranchAddress("cl_ts0", &cl_ts0, &b_cl_ts0);
   fChain->SetBranchAddress("cl_ts1", &cl_ts1, &b_cl_ts1);
   fChain->SetBranchAddress("cl_unixs", &cl_unixs, &b_cl_unixs);
   fChain->SetBranchAddress("cl_nhits", &cl_nhits, &b_cl_nhits);
   fChain->SetBranchAddress("cl_tagger", &cl_tagger, &b_cl_tagger);
   fChain->SetBranchAddress("cl_composition", &cl_composition, &b_cl_composition);
   fChain->SetBranchAddress("cl_truth_trackid", &cl_truth_trackid, &b_cl_truth_trackid);
   fChain->SetBranchAddress("cl_truth_completeness", &cl_truth_completeness, &b_cl_truth_completeness);
   fChain->SetBranchAddress("cl_truth_purity", &cl_truth_purity, &b_cl_truth_purity);
   fChain->SetBranchAddress("cl_truth_hit_completeness", &cl_truth_hit_completeness, &b_cl_truth_hit_completeness);
   fChain->SetBranchAddress("cl_truth_hit_purity", &cl_truth_hit_purity, &b_cl_truth_hit_purity);
   fChain->SetBranchAddress("cl_truth_pdg", &cl_truth_pdg, &b_cl_truth_pdg);
   fChain->SetBranchAddress("cl_truth_energy", &cl_truth_energy, &b_cl_truth_energy);
   fChain->SetBranchAddress("cl_truth_time", &cl_truth_time, &b_cl_truth_time);
   fChain->SetBranchAddress("cl_truth_x", &cl_truth_x, &b_cl_truth_x);
   fChain->SetBranchAddress("cl_truth_y", &cl_truth_y, &b_cl_truth_y);
   fChain->SetBranchAddress("cl_truth_z", &cl_truth_z, &b_cl_truth_z);
   fChain->SetBranchAddress("cl_truth_core_energy", &cl_truth_core_energy, &b_cl_truth_core_energy);
   fChain->SetBranchAddress("cl_truth_core_time", &cl_truth_core_time, &b_cl_truth_core_time);
   fChain->SetBranchAddress("cl_truth_core_x", &cl_truth_core_x, &b_cl_truth_core_x);
   fChain->SetBranchAddress("cl_truth_core_y", &cl_truth_core_y, &b_cl_truth_core_y);
   fChain->SetBranchAddress("cl_truth_core_z", &cl_truth_core_z, &b_cl_truth_core_z);
   fChain->SetBranchAddress("cl_has_sp", &cl_has_sp, &b_cl_has_sp);
   fChain->SetBranchAddress("cl_sp_x", &cl_sp_x, &b_cl_sp_x);
   fChain->SetBranchAddress("cl_sp_ex", &cl_sp_ex, &b_cl_sp_ex);
   fChain->SetBranchAddress("cl_sp_y", &cl_sp_y, &b_cl_sp_y);
   fChain->SetBranchAddress("cl_sp_ey", &cl_sp_ey, &b_cl_sp_ey);
   fChain->SetBranchAddress("cl_sp_z", &cl_sp_z, &b_cl_sp_z);
   fChain->SetBranchAddress("cl_sp_ez", &cl_sp_ez, &b_cl_sp_ez);
   fChain->SetBranchAddress("cl_sp_pe", &cl_sp_pe, &b_cl_sp_pe);
   fChain->SetBranchAddress("cl_sp_time", &cl_sp_time, &b_cl_sp_time);
   fChain->SetBranchAddress("cl_sp_etime", &cl_sp_etime, &b_cl_sp_etime);
   fChain->SetBranchAddress("cl_sp_complete", &cl_sp_complete, &b_cl_sp_complete);
   fChain->SetBranchAddress("td_tag_trackid", &td_tag_trackid, &b_td_tag_trackid);
   fChain->SetBranchAddress("td_tag_pdg", &td_tag_pdg, &b_td_tag_pdg);
   fChain->SetBranchAddress("td_tag_tagger", &td_tag_tagger, &b_td_tag_tagger);
   fChain->SetBranchAddress("td_tag_energy", &td_tag_energy, &b_td_tag_energy);
   fChain->SetBranchAddress("td_tag_time", &td_tag_time, &b_td_tag_time);
   fChain->SetBranchAddress("td_tag_x", &td_tag_x, &b_td_tag_x);
   fChain->SetBranchAddress("td_tag_y", &td_tag_y, &b_td_tag_y);
   fChain->SetBranchAddress("td_tag_z", &td_tag_z, &b_td_tag_z);
   fChain->SetBranchAddress("td_tag_reco_status", &td_tag_reco_status, &b_td_tag_reco_status);
   fChain->SetBranchAddress("td_trackid", &td_trackid, &b_td_trackid);
   fChain->SetBranchAddress("td_pdg", &td_pdg, &b_td_pdg);
   fChain->SetBranchAddress("td_energy", &td_energy, &b_td_energy);
   fChain->SetBranchAddress("td_time", &td_time, &b_td_time);
   fChain->SetBranchAddress("td_reconstructable", &td_reconstructable, &b_td_reconstructable);
   fChain->SetBranchAddress("td_reco_status", &td_reco_status, &b_td_reco_status);
   fChain->SetBranchAddress("td_reco_triple", &td_reco_triple, &b_td_reco_triple);
   fChain->SetBranchAddress("tr_start_x", &tr_start_x, &b_tr_start_x);
   fChain->SetBranchAddress("tr_start_y", &tr_start_y, &b_tr_start_y);
   fChain->SetBranchAddress("tr_start_z", &tr_start_z, &b_tr_start_z);
   fChain->SetBranchAddress("tr_end_x", &tr_end_x, &b_tr_end_x);
   fChain->SetBranchAddress("tr_end_y", &tr_end_y, &b_tr_end_y);
   fChain->SetBranchAddress("tr_end_z", &tr_end_z, &b_tr_end_z);
   fChain->SetBranchAddress("tr_dir_x", &tr_dir_x, &b_tr_dir_x);
   fChain->SetBranchAddress("tr_dir_y", &tr_dir_y, &b_tr_dir_y);
   fChain->SetBranchAddress("tr_dir_z", &tr_dir_z, &b_tr_dir_z);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("tr_etime", &tr_etime, &b_tr_etime);
   fChain->SetBranchAddress("tr_pe", &tr_pe, &b_tr_pe);
   fChain->SetBranchAddress("tr_length", &tr_length, &b_tr_length);
   fChain->SetBranchAddress("tr_tof", &tr_tof, &b_tr_tof);
   fChain->SetBranchAddress("tr_theta", &tr_theta, &b_tr_theta);
   fChain->SetBranchAddress("tr_phi", &tr_phi, &b_tr_phi);
   fChain->SetBranchAddress("tr_triple", &tr_triple, &b_tr_triple);
   fChain->SetBranchAddress("tr_truth_trackid", &tr_truth_trackid, &b_tr_truth_trackid);
   fChain->SetBranchAddress("tr_truth_completeness", &tr_truth_completeness, &b_tr_truth_completeness);
   fChain->SetBranchAddress("tr_truth_purity", &tr_truth_purity, &b_tr_truth_purity);
   fChain->SetBranchAddress("tr_truth_pdg", &tr_truth_pdg, &b_tr_truth_pdg);
   fChain->SetBranchAddress("tr_truth_energy", &tr_truth_energy, &b_tr_truth_energy);
   fChain->SetBranchAddress("tr_truth_time", &tr_truth_time, &b_tr_truth_time);
   fChain->SetBranchAddress("tr_truth_start_x", &tr_truth_start_x, &b_tr_truth_start_x);
   fChain->SetBranchAddress("tr_truth_start_y", &tr_truth_start_y, &b_tr_truth_start_y);
   fChain->SetBranchAddress("tr_truth_start_z", &tr_truth_start_z, &b_tr_truth_start_z);
   fChain->SetBranchAddress("tr_truth_end_x", &tr_truth_end_x, &b_tr_truth_end_x);
   fChain->SetBranchAddress("tr_truth_end_y", &tr_truth_end_y, &b_tr_truth_end_y);
   fChain->SetBranchAddress("tr_truth_end_z", &tr_truth_end_z, &b_tr_truth_end_z);
   fChain->SetBranchAddress("tr_truth_dir_x", &tr_truth_dir_x, &b_tr_truth_dir_x);
   fChain->SetBranchAddress("tr_truth_dir_y", &tr_truth_dir_y, &b_tr_truth_dir_y);
   fChain->SetBranchAddress("tr_truth_dir_z", &tr_truth_dir_z, &b_tr_truth_dir_z);
   fChain->SetBranchAddress("tr_truth_particle_energy", &tr_truth_particle_energy, &b_tr_truth_particle_energy);
   fChain->SetBranchAddress("tr_truth_length", &tr_truth_length, &b_tr_truth_length);
   fChain->SetBranchAddress("tr_truth_tof", &tr_truth_tof, &b_tr_truth_tof);
   fChain->SetBranchAddress("tr_truth_theta", &tr_truth_theta, &b_tr_truth_theta);
   fChain->SetBranchAddress("tr_truth_phi", &tr_truth_phi, &b_tr_truth_phi);
   fChain->SetBranchAddress("tpc_start_x", &tpc_start_x, &b_tpc_start_x);
   fChain->SetBranchAddress("tpc_start_y", &tpc_start_y, &b_tpc_start_y);
   fChain->SetBranchAddress("tpc_start_z", &tpc_start_z, &b_tpc_start_z);
   fChain->SetBranchAddress("tpc_end_x", &tpc_end_x, &b_tpc_end_x);
   fChain->SetBranchAddress("tpc_end_y", &tpc_end_y, &b_tpc_end_y);
   fChain->SetBranchAddress("tpc_end_z", &tpc_end_z, &b_tpc_end_z);
   fChain->SetBranchAddress("tpc_dir_x", &tpc_dir_x, &b_tpc_dir_x);
   fChain->SetBranchAddress("tpc_dir_y", &tpc_dir_y, &b_tpc_dir_y);
   fChain->SetBranchAddress("tpc_dir_z", &tpc_dir_z, &b_tpc_dir_z);
   fChain->SetBranchAddress("tpc_length", &tpc_length, &b_tpc_length);
   fChain->SetBranchAddress("tpc_track_score", &tpc_track_score, &b_tpc_track_score);
   fChain->SetBranchAddress("tpc_truth_trackid", &tpc_truth_trackid, &b_tpc_truth_trackid);
   fChain->SetBranchAddress("tpc_truth_pdg", &tpc_truth_pdg, &b_tpc_truth_pdg);
   fChain->SetBranchAddress("tpc_truth_energy", &tpc_truth_energy, &b_tpc_truth_energy);
   fChain->SetBranchAddress("tpc_truth_time", &tpc_truth_time, &b_tpc_truth_time);
   fChain->SetBranchAddress("tpc_sp_matchable", &tpc_sp_matchable, &b_tpc_sp_matchable);
   fChain->SetBranchAddress("tpc_sp_matched", &tpc_sp_matched, &b_tpc_sp_matched);
   fChain->SetBranchAddress("tpc_sp_good_match", &tpc_sp_good_match, &b_tpc_sp_good_match);
   fChain->SetBranchAddress("tpc_sp_time", &tpc_sp_time, &b_tpc_sp_time);
   fChain->SetBranchAddress("tpc_sp_score", &tpc_sp_score, &b_tpc_sp_score);
   fChain->SetBranchAddress("tpc_tr_matchable", &tpc_tr_matchable, &b_tpc_tr_matchable);
   fChain->SetBranchAddress("tpc_tr_matched", &tpc_tr_matched, &b_tpc_tr_matched);
   fChain->SetBranchAddress("tpc_tr_good_match", &tpc_tr_good_match, &b_tpc_tr_good_match);
   fChain->SetBranchAddress("tpc_tr_time", &tpc_tr_time, &b_tpc_tr_time);
   fChain->SetBranchAddress("tpc_tr_score", &tpc_tr_score, &b_tpc_tr_score);
   Notify();
}

Bool_t tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_cxx
