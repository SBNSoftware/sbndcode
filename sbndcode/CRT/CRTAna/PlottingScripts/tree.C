#define tree_cxx
#include "tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L tree.C
//      root> tree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("td_tag_trackid",1);
   fChain->SetBranchStatus("td_tag_pdg",1);
   fChain->SetBranchStatus("td_tag_tagger",1);
   fChain->SetBranchStatus("td_tag_time",1);
   fChain->SetBranchStatus("td_tag_energy",1);

   int n = 0, ntot = 0, nev = 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      bool gotOneEv = false;

      for(int td_i = 0; td_i < td_tag_trackid->size(); ++td_i)
	{
	  if(td_tag_energy->at(td_i) < 0.003)
	    continue;

	  bool gotOne = false;

	  for(int td_j = 0; td_j < td_tag_trackid->size(); ++td_j)
	    {
	      if(td_i == td_j)
		continue;

	      if(td_tag_energy->at(td_j) < 0.003)
		continue;

	      if(td_tag_tagger->at(td_i) == td_tag_tagger->at(td_j)
		 && td_tag_trackid->at(td_i) != td_tag_trackid->at(td_j)
		 && abs(td_tag_time->at(td_i) - td_tag_time->at(td_j)) < 50.
		 && abs(td_tag_pdg->at(td_i)) == 13 && abs(td_tag_pdg->at(td_j)) == 13)
		{
		  gotOne = true;
		  gotOneEv = true;
		  std::cout << td_tag_pdg->at(td_i) <<  " & " << td_tag_pdg->at(td_j) << std::endl;
		}
	    }
	  if(gotOne)
	    ++n;
	  ++ntot;
	}
      if(gotOneEv)
	++nev;
   }

   std::cout << "\n\nn: " << n << " / " << ntot << " (" << n * 100. / ntot << ")" << std::endl;
   std::cout << "\n\nnev: " << nev << " / " << nentries << " (" << nev * 100. / nentries << ")" << std::endl;
}
