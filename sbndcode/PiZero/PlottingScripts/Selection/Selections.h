#pragma once

#include "Cuts.h"
#include "Categories.h"

SelectionParams ncpizero_incl = { "ncpizero_incl",
                                  ncpizero_incl_cuts,
                                  ncpizero_incl_broad_cuts,
                                  ncpizero_incl_categories,
                                  true_ncpizero_incl_cut,
                                  { 0, 8 },
                                  { 1, 2, 3, 4, 6, 7 }
};

SelectionParams ncpizero_0p0pi = { "ncpizero_0p0pi",
                                   ncpizero_0p0pi_cuts,
                                   ncpizero_0p0pi_broad_cuts,
                                   ncpizero_0p0pi_categories,
                                   true_ncpizero_0p0pi_cut,
                                   { 0, 1, 8 },
                                   { 2, 3, 4, 5, 6, 7 }
};

SelectionParams ncpizero_Np0pi = { "ncpizero_Np0pi",
                                   ncpizero_Np0pi_cuts,
                                   ncpizero_Np0pi_broad_cuts,
                                   ncpizero_Np0pi_categories,
                                   true_ncpizero_Np0pi_cut,
                                   { 0, 1, 8 },
                                   { 2, 3, 4, 5, 6, 7 }
};

SelectionParams ccpizero = { "ccpizero",
                             ccpizero_cuts,
                             ccpizero_broad_cuts,
                             ccpizero_categories,
                             true_ccpizero_cut
};

std::vector<SelectionParams> selections = { ncpizero_incl,
                                            ncpizero_0p0pi,
                                            ncpizero_Np0pi
};
