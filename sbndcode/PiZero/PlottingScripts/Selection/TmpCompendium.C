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

void TmpCompendium(const TString &productionVersion)
{
  for(auto const& selection : {ncpizero_incl})
    {
      Selection(productionVersion, selection, selection_plots);
      Observables(productionVersion, selection, observables);
      Observables(productionVersion, selection, observables_corr);
    }
}
