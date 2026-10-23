const SpillVar kTotalSlices([](const caf::SRSpillProxy* sp) -> double {
    return sp->slc.size();
  });

const SpillVar kNuSlices([](const caf::SRSpillProxy* sp) -> double {
    int n = 0;
    for(auto const& slc : sp->slc) {
      if(!slc.is_clear_cosmic) n++;
    }
    return n;
  });

const Var kCompleteness([](const caf::SRSliceProxy *slc) -> double {
    return slc->tmatch.eff;
  });
