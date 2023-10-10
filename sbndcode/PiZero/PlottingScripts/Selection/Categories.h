#include "/sbnd/app/users/hlay/plotting_utils/Structs.h"

const std::vector<Cut> ncpizero_categories = {
  { "Signal", "slc_true_event_type==0 && slc_comp>.5", "Signal (NC #pi^{0})", kMagenta+2 },
  { "NC", "slc_true_event_type==1", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type==2", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type==3", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type==4", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type==5", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type==6", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_event_type==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};

const std::vector<Cut> selection_categories = ncpizero_categories;

const std::vector<Cut> ncpizero_0p0pi_categories = {
  { "Signal", "slc_true_ccnc && slc_true_n_neutral_pions==1 && slc_true_n_protons==0 && slc_true_n_charged_pions==0 && slc_comp>.5", "Signal (NC 1#pi^{0}0p0#pi^{1pm})", kMagenta+2 },
  { "NC", "slc_true_ccnc && ((slc_true_n_neutral_pions==1 && (slc_true_n_protons!=0 || slc_true_n_charged_pions!=0)) || slc_true_n_neutral_pions>1)", "Other NC#pi^{0}", kBlue+1 },
  { "NC", "slc_true_ccnc && slc_true_n_neutral_pions==0", "Other NC", kOrange+2 },
  { "CCNuMu", "slc_true_event_type==2", "CC #nu_{#mu}", kGreen+2 },
  { "CCNuE", "slc_true_event_type==3", "CC #nu_{e}", kCyan+2 },
  { "Dirt", "slc_true_event_type==4", "Dirt", kOrange+3 },
  { "NonFVNu", "slc_true_event_type==5", "Non-FV #nu", kGray+2 },
  { "Cosmic", "slc_true_event_type==6", "Cosmic", kRed+1 },
  { "BadRecoSignal", "slc_true_ccnc && slc_true_n_neutral_pions==1 && slc_true_n_protons==0 && slc_true_n_charged_pions==0 && slc_comp<=.5", "Bad Reco Signal", kBlack }
};
