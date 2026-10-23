const Cut kIsNeutrinoCandidateSlice([](const caf::SRSliceProxy *slc) {
    return !slc->is_clear_cosmic;
  });
