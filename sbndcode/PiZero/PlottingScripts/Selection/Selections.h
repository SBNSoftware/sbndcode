#include "Cuts.h"
#include "Categories.h"

SelectionParams ncpizero_incl = { "ncpizero_incl",
                                  ncpizero_incl_cuts,
                                  ncpizero_incl_broad_cuts,
                                  ncpizero_incl_categories,
                                  true_ncpizero_incl_cut
};

SelectionParams ncpizero_0p0pi = { "ncpizero_0p0pi",
                                   ncpizero_0p0pi_cuts,
                                   ncpizero_0p0pi_broad_cuts,
                                   ncpizero_0p0pi_categories,
                                   true_ncpizero_0p0pi_cut
};

SelectionParams ncpizero_Np0pi = { "ncpizero_Np0pi",
                                   ncpizero_Np0pi_cuts,
                                   ncpizero_Np0pi_broad_cuts,
                                   ncpizero_Np0pi_categories,
                                   true_ncpizero_Np0pi_cut
};
