//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 13 13:59:29 2021 by ROOT version 6.22/08
// from TTree ncrecoTree/Tree with NC Reco information
// found on file: ncrecoTree.root
//////////////////////////////////////////////////////////

#ifndef ncreco_h
#define ncreco_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class ncreco {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           SubRun;
   Int_t           Event;
   vector<double>  *E;
   vector<double>  *P;
   vector<double>  *Px;
   vector<double>  *Py;
   vector<double>  *Pz;
   vector<int>     *Status;
   Int_t           Mode;
   Int_t           CCNC;
   double	   NuE;
   Int_t           NuPdg;
   vector<int>     *Pdg;
   Int_t           NParticles;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_E;   //!
   TBranch        *b_P;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_Status;   //!
   TBranch        *b_Mode;   //!
   TBranch        *b_NuPdg;   //!
   TBranch        *b_CCNC;   //!
   TBranch        *b_NuE;   //!
   TBranch        *b_Pdg;   //!
   TBranch        *b_NParticles;   //!

   ncreco(TTree *tree=0);
   virtual ~ncreco();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ncreco_cxx
ncreco::ncreco(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ncrecoTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ncrecoTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ncrecoTree.root:/ncreco");
      dir->GetObject("ncrecoTree",tree);

   }
   Init(tree);
}

ncreco::~ncreco()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ncreco::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ncreco::LoadTree(Long64_t entry)
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

void ncreco::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   E = 0;
   P = 0;
   Px = 0;
   Py = 0;
   Pz = 0;
   Pdg = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("P", &P, &b_P);
   fChain->SetBranchAddress("Px", &Px, &b_Px);
   fChain->SetBranchAddress("Py", &Py, &b_Py);
   fChain->SetBranchAddress("Pz", &Pz, &b_Pz);
   fChain->SetBranchAddress("Status", &Status, &b_Status);
   fChain->SetBranchAddress("Mode", &Mode, &b_Mode);
   fChain->SetBranchAddress("NuE", &NuE, &b_NuE);
   fChain->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
   fChain->SetBranchAddress("NuPdg", &NuPdg, &b_NuPdg);
   fChain->SetBranchAddress("Pdg", &Pdg, &b_Pdg);
   fChain->SetBranchAddress("NParticles", &NParticles, &b_NParticles);
   Notify();
}

Bool_t ncreco::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ncreco::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ncreco::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ncreco_cxx
