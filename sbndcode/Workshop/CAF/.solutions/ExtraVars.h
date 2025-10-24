const SpillMultiVar kNeutrinoEnergy([](const caf::SRSpillProxy *sp) -> std::vector<double> {
    std::vector<double> energies;

    for(auto const& nu : sp->mc.nu)
      energies.push_back(nu.E);

    return energies;
  });

const SpillVar kNSlices([](const caf::SRSpillProxy *sp) -> int {
    return sp->nslc;
  });

const SpillVar kNNeutrinoCandidateSlices([](const caf::SRSpillProxy *sp) -> int {
    int n = 0;

    for(auto const& slc : sp->slc)
      {
        if(!slc.is_clear_cosmic)
          ++n;
      }

    return n;
  });

const Var kSliceCompleteness([](const caf::SRSliceProxy *slc) -> double {
    if(slc->tmatch.index == -999)
      return 0.;

    return slc->tmatch.eff;
  });
