struct Selection {
  TString   name;
  TString   cut;
  TString   signal;
  XSecPlot *plot;
};

typedef std::vector<Selection> Selections;

Selection ncpizero_incl  = { "ncpizero_incl", "sel_incl", "event_type_incl" };
Selection ncpizero_0p0pi = { "ncpizero_0p0pi", "sel_0p0pi", "event_type_0p0pi" };
Selection ncpizero_Np0pi = { "ncpizero_Np0pi", "sel_Np0pi", "event_type_Np0pi" };

Selections selections = { ncpizero_incl, ncpizero_0p0pi, ncpizero_Np0pi };
