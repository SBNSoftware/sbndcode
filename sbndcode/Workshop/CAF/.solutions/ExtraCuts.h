const SpillCut kExactlyOneNuSlice([](const caf::SRSpillProxy* sp) {
    int nu_slices = 0;
    for(auto const& slc : sp->slc) {
      if(!slc.is_clear_cosmic) nu_slices++;
    }
    return nu_slices == 1;
  });
