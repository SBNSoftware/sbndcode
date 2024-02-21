struct Selection {
  TString   name;
  TString   cut;
  TString   signal;
  XSecPlot *plot;
};

typedef std::vector<Selection> Selections;
