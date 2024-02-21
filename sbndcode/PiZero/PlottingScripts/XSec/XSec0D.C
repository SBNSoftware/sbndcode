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

  Selection incl = { "ncpizero_incl", "sel_incl", "event_type_incl==0" };

  XSecPlot *inclPlot = new XSecPlot("incl0D", "#sigma (cm^{2}/nucleon)", 1);

  CalculateNominalBin(samples, inclPlot, 0, incl);

  Bin& bin = inclPlot->GetBin(0);
  bin.CalculateXSecPurity(nTargets, intFlux);
  bin.Print();
}
