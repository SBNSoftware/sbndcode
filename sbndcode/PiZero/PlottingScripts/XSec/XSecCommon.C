#include "XSecPlot.h"
#include "XSecSample.h"
#include "Selection.h"

void CalculateNominalBin(XSecSamples samples, Bin *bin, const Selection selection)
{
  double trueSig = 0., rawCount = 0., sigCount = 0.;

  const double normalisation = samples[0].scaling;

  for(auto const& sample : samples)
    {
      trueSig += sample.nutree->GetEntries(selection.signal) * (sample.scaling / normalisation);
      rawCount += sample.slicetree->GetEntries(selection.cut) * (sample.scaling / normalisation);
      sigCount += sample.slicetree->GetEntries(selection.cut + " && " + selection.signal + " && comp>.5")
        * (sample.scaling / normalisation);
    }

  double bkgdCount = rawCount - sigCount;

  double scaledCount     = rawCount * normalisation;
  double scaledBkgdCount = bkgdCount * normalisation;

  double purity     = (rawCount - bkgdCount) / rawCount;
  double efficiency = sigCount / trueSig;

  double fracStatErr = std::sqrt(rawCount) / rawCount;

  bin->SetRawCount(rawCount);
  bin->SetScaledCount(scaledCount);
  bin->SetBackgroundCount(bkgdCount);
  bin->SetScaledBackgroundCount(scaledBkgdCount);
  bin->SetPurity(purity);
  bin->SetEfficiency(efficiency);
  bin->SetFracStatErr(fracStatErr);
}

void CalculateNominal(XSecSamples samples, const Selection selection)
{
  XSecPlot *plot = selection.plot;

  for(Bin* bin : plot->GetBins())
    {
      CalculateNominalBin(samples, bin, selection);
    }
}

void CalculateNominalXSecPurity(const Selection selection, const double nTargets, const double intFlux)
{
  XSecPlot *plot = selection.plot;

  for(Bin* bin : plot->GetBins())
    {
      bin->CalculateXSecPurity(nTargets, intFlux);
    }
}
