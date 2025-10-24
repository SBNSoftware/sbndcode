const SpillCut kSingleNeutrinoCandidateSlice([](const caf::SRSpillProxy *sp) {
    return kNNeutrinoCandidateSlices(sp) == 1;
  });
