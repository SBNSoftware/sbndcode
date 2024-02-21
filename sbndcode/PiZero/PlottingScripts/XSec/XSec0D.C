#include "XSecCommon.C"
#include "Common.C"

void XSec0D(const TString productionVersion, const TString saveDirExt)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_zero_d/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  TChain *rockboxNus = new TChain("ncpizeroxsectrees/neutrinos");
  rockboxNus->Add(rockboxFile);
  TChain *rockboxSlices = new TChain("ncpizeroxsectrees/slices");
  rockboxSlices->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroxsectrees/subruns");
  rockboxSubruns->Add(rockboxFile);

  TChain *intimeNus = new TChain("ncpizeroxsectrees/neutrinos");
  intimeNus->Add(intimeFile);
  TChain *intimeSlices = new TChain("ncpizeroxsectrees/slices");
  intimeSlices->Add(intimeFile);
  TChain *intimeSubruns = new TChain("ncpizeroxsectrees/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  XSecSamples samples = { { rockboxNus, rockboxSlices, rockboxScaling },
                          { intimeNus, intimeSlices, intimeScaling }
  };

  XSecPlot *xsec_incl  = new XSecPlot("xsec_incl", "#sigma (cm^{2}/nucleon)", 1);
  XSecPlot *xsec_0p0pi = new XSecPlot("xsec_0p0pi", "#sigma (cm^{2}/nucleon)", 1);
  XSecPlot *xsec_Np0pi = new XSecPlot("xsec_Np0pi", "#sigma (cm^{2}/nucleon)", 1);

  Selection ncpizero_incl  = { "ncpizero_incl", "sel_incl", "event_type_incl==0", xsec_incl };
  Selection ncpizero_0p0pi = { "ncpizero_0p0pi", "sel_0p0pi", "event_type_0p0pi==0", xsec_0p0pi };
  Selection ncpizero_Np0pi = { "ncpizero_Np0pi", "sel_Np0pi", "event_type_Np0pi==0", xsec_Np0pi };

  Selections selections = { ncpizero_incl, ncpizero_0p0pi, ncpizero_Np0pi };

  for(auto const& selection : selections)
    {
      CalculateNominal(samples, selection);

      Bin *bin = selection.plot->GetBin(0);
      bin->CalculateXSecPurity(nTargets, intFlux);
      bin->Print();
    }
}
