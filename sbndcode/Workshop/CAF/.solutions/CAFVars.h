const MultiVar kChildTrackLengths([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> trackLengths;
    for(auto const& pfp : slc->reco.pfp) {
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue; // 2026 Bug Fix: remove showers
      trackLengths.push_back(pfp.trk.len);
    }
    return trackLengths;
  });

const MultiVar kChildTrackdEdx([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> dEdx;
    for(auto const &pfp : slc->reco.pfp) {
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue;
      for(auto const& point : pfp.trk.calo[2].points) {
	dEdx.push_back(point.dedx);
      }
    }
    return dEdx;
  });

const MultiVar kChildTrackResRange([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> rr;
    for(auto const &pfp : slc->reco.pfp) {
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue;
      for(auto const& point : pfp.trk.calo[2].points) {
	rr.push_back(point.rr);
      }
    }
    return rr;
  });

const Var kLongestTrack([](const caf::SRSliceProxy *slc) -> double {
    double maxLength = -1.0;
    int maxLengthID  = -1;
    int i = -1;
    for(auto const& pfp : slc->reco.pfp) {
      ++i;
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue;
      if(pfp.trk.len > maxLength) {
	maxLength   = pfp.trk.len;
	maxLengthID = i;
      }
    }
    return maxLengthID; // Returns as double for CAFAna Var
  });

const Var kChildTrackLengthLongestTrack([](const caf::SRSliceProxy *slc) -> double {
    int maxLengthID = (int)kLongestTrack(slc);
    if(maxLengthID < 0) return -999.0;
    return slc->reco.pfp[maxLengthID].trk.len;
  });

const MultiVar kChildTrackLengthOtherTracks([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> trackLengths;
    int maxLengthID = (int)kLongestTrack(slc);
    int i = -1;
    for(auto const& pfp : slc->reco.pfp) {
      ++i;
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue;
      if(maxLengthID == i) continue;
      trackLengths.push_back(pfp.trk.len);
    }
    return trackLengths;
  });

const MultiVar kChildTrackdEdxLongestTrack([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> dEdx;
    int maxLengthID = (int)kLongestTrack(slc);
    if(maxLengthID < 0) return dEdx;
    for(auto const& point : slc->reco.pfp[maxLengthID].trk.calo[2].points) {
      dEdx.push_back(point.dedx);
    }
    return dEdx;
  });

const MultiVar kChildTrackdEdxOtherTracks([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> dEdx;
    int maxLengthID = (int)kLongestTrack(slc);
    int i = -1;
    for(auto const& pfp : slc->reco.pfp) {
      ++i;
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue;
      if(maxLengthID == i) continue;
      for(auto const& point : pfp.trk.calo[2].points) {
	dEdx.push_back(point.dedx);
      }
    }
    return dEdx;
  });

const MultiVar kChildTrackResRangeLongestTrack([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> rr;
    int maxLengthID = (int)kLongestTrack(slc);
    if(maxLengthID < 0) return rr;
    for(auto const& point : slc->reco.pfp[maxLengthID].trk.calo[2].points) {
      rr.push_back(point.rr);
    }
    return rr;
  });

const MultiVar kChildTrackResRangeOtherTracks([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> rr;
    int maxLengthID = (int)kLongestTrack(slc);
    int i = -1;
    for(auto const& pfp : slc->reco.pfp) {
      ++i;
      if(!pfp.parent_is_primary || pfp.parent == -1) continue;
      if(pfp.trackScore < 0.5) continue;
      if(maxLengthID == i) continue;
      for(auto const& point : pfp.trk.calo[2].points) {
	rr.push_back(point.rr);
      }
    }
    return rr;
  });

const Var kOpT0Time([](const caf::SRSliceProxy *slc) -> double {
    return slc->opt0.time;
  });
