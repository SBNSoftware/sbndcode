#include "XSecPlot.h"
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

  XSecPlot *inclPlot = new XSecPlot("incl0D", "#sigma (cm^{2}/nucleon)", 1);

  int rockboxTrueSig      = rockboxNus->GetEntries("event_type_incl==0");
  int rockboxSelSlices    = rockboxSlices->GetEntries("sel_incl");
  int intimeSelSlices     = intimeSlices->GetEntries("sel_incl");
  int rockboxSelSigSlices = rockboxSlices->GetEntries("sel_incl && event_type_incl==0 && comp>.5");

  double rawCount    = rockboxSelSlices + intimeSelSlices * (intimeScaling / rockboxScaling);
  double scaledCount = rawCount * rockboxScaling;

  double bkgdCount       = (rockboxSelSlices - rockboxSelSigSlices)
    + intimeSelSlices * (intimeScaling / rockboxScaling);
  double scaledBkgdCount = bkgdCount * rockboxScaling;

  double purity     = (rawCount - bkgdCount) / rawCount;
  double efficiency = rockboxSelSigSlices / rockboxTrueSig;

  Bin& bin = inclPlot->GetBin(0);

  bin.SetRawCount(rawCount);
  bin.SetScaledCount(scaledCount);
  bin.SetBackgroundCount(bkgdCount);
  bin.SetScaledBackgroundCount(scaledBkgdCount);
  bin.SetPurity(purity);
  bin.SetEfficiency(efficiency);
  bin.CalculateXSecPurity(nTargets, intFlux);

  bin.Print();
}
