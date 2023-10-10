#include "Enumerate.h"
#include "Categories.h"

#include <numeric>

void MutuallyExclusiveAndComplete(std::vector<Cut> categories)
{
  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv3/NCPiZeroAv3_rockbox.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv3/NCPiZeroAv3_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add(rockboxFile);
  events->Add(intimeFile);

  const int N = events->Draw("slc_true_e", "");

  std::vector<int> n(categories.size(), 0);

  std::cout << "\nTotal slices: " << N << '\n' << std::endl;
  TCut totalcuts = "";

  for(auto&& [i, category] : enumerate(categories))
    {
      totalcuts = totalcuts || category.cut;

      n[i] = events->Draw("slc_true_e", category.cut);
      std::cout << "\tCategory: " << setw(15) << category.name 
                << right << setw(10) << n[i] << std::endl;

      for(auto&& [j, categorytwo] : enumerate(categories))
        {
          if(i==j)
            continue;

          const int overlap = events->Draw("slc_true_e", category.cut && categorytwo.cut);
          if(overlap != 0)
            std::cout << "\n\t\tUh Oh! Overlap of " << overlap << " events between categories "
                      << category.name << " and " << categorytwo.name << std::endl;
        }
    }

  int total = std::reduce(n.begin(), n.end());

  std::cout << "\nTotal: " << total << " vs. True: " << N << " (" << N - total << ")\n" << std::endl;

  if((N - total) > 0)
    {
      std::cout << "Some slices not accounted for..." << std::endl;
      std::cout << events->Draw("slc_true_e", !totalcuts) << std::endl;
      events->Scan("slc_true_event_type:slc_n_pfps:slc_is_clear_cosmic:slc_comp:slc_pur:slc_true_mctruth_id", !totalcuts);
    }
}
