struct XSecSample {
  TString name;
  TChain *nutree;
  TChain *slicetree;
  double  scaling;
};

typedef std::vector<XSecSample> XSecSamples;
