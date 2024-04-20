struct XSecSample {
  TString       name;
  TChain       *nutree;
  TChain       *slicetree;
  double        scaling;
  std::set<int> mask;
};

typedef std::vector<XSecSample> XSecSamples;
