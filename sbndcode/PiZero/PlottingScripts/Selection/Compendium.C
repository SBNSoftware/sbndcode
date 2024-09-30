#include "Observables2DPlots.C"
#include "Observables.C"
#include "ObservablesResolution.C"
#include "ObservablesResolution2.C"
#include "ObservablesResolution3.C"
#include "PhotonRecoEff.C"
#include "RecoEff.C"
#include "Selection.C"
#include "SelectionEff.C"
#include "TrueEventMode2DPlots.C"
#include "TrueEventModePlots.C"

void Compendium(const TString &productionVersion)
{
  for(auto const& selection : selections)
    {
      Selection(productionVersion, selection, no_plots);
      SelectionEff(productionVersion, selection, true_slc_observables);
      Observables2DPlots(productionVersion, observables_twod_sets, selection);
      Observables(productionVersion, selection, observables);
      Observables(productionVersion, selection, observables_corr);
    }

  ObservablesResolution(productionVersion);
  ObservablesResolution2(productionVersion);
  ObservablesResolution3(productionVersion);
  PhotonRecoEff(productionVersion);
  RecoEff(productionVersion);
  TrueEventMode2DPlots(productionVersion, true_observables_twod_sets, true_categories);
  TrueEventModePlots(productionVersion, true_observables, true_categories);
}
