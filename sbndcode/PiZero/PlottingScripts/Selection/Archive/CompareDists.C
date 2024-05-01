#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"

void SimpleCompare(TCanvas *c, const TString &fileNameA, const TString &fileNameB,
		   const TString &treeName, const Plot &plot, const TString &saveDir);

void CompareDists()
{
  const std::vector<Plot> plots = { { "ncpizero_incl_pizero_mom", "pizero_mom", ";p_{#pi^{0}} (MeV);", 10, 0, 2000, kBlack, false, "event_type_incl==0" },
				    { "ncpizero_incl_cos_theta_pizero", "cos_theta_pizero", ";cos(#theta_{#pi^{0}});", 10, -1, 1, kBlack, false, "event_type_incl==0" } };

  for(auto const &plot : plots)
    {
      TCanvas *c = new TCanvas("c" + plot.name, "c" + plot.name);
      SimpleCompare(c, "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroBv3/NCPiZeroBv3_rockbox.root",
		    "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroBv3/NCPiZeroBv3_ncpizero.root",
		    "ncpizeroxsectrees/neutrinos", plot, "");
    }
}

void SimpleCompare(TCanvas *c, const TString &fileNameA, const TString &fileNameB,
		   const TString &treeName, const Plot &plot, const TString &saveDir)
{
  c->cd();

  TChain *treeA = new TChain(treeName);
  TChain *treeB = new TChain(treeName);

  treeA->Add(fileNameA);
  treeB->Add(fileNameB);

  TH1D* histA = new TH1D(plot.name + "_A", "", plot.nbins, plot.xlow, plot.xhigh);
  TH1D* histB = new TH1D(plot.name + "_B", "", plot.nbins, plot.xlow, plot.xhigh);

  treeA->Draw(plot.var + ">>" + plot.name + "_A", plot.req);
  treeB->Draw(plot.var + ">>" + plot.name + "_B", plot.req);

  histA->DrawNormalized("histe");
  histB->SetLineColor(kRed+2);
  histB->DrawNormalized("histesame");
}
