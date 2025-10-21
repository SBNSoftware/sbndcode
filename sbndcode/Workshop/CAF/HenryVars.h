const MultiVar kChildTrackLengths([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> trackLengths;

    for(auto const& pfp : slc->reco.pfp)
      {
	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	trackLengths.push_back(pfp.trk.len);
      }

    return trackLengths;
  });

const MultiVar kChildTrackdEdx([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> dEdx;

    for(auto const &pfp : slc->reco.pfp)
      {
	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	// Only interested in the collection plane (2)
	for(auto const& point : pfp.trk.calo[2].points)
	  dEdx.push_back(point.dedx);
      }

    return dEdx;
  });

const MultiVar kChildTrackResRange([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> rr;

    for(auto const &pfp : slc->reco.pfp)
      {
	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	// Only interested in the collection plane (2)
	for(auto const& point : pfp.trk.calo[2].points)
	  rr.push_back(point.rr);
      }

    return rr;
  });

const Var kLongestTrack([](const caf::SRSliceProxy *slc) -> int {
    double maxLength = std::numeric_limits<double>::lowest();
    int maxLengthID  = std::numeric_limits<int>::signaling_NaN();

    int i = -1;

    for(auto const& pfp : slc->reco.pfp)
      {
	++i;

	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	if(pfp.trk.len > maxLength)
	  {
	    maxLength   = pfp.trk.len;
	    maxLengthID = i;
	  }
      }

    return maxLengthID;
  });

const Var kChildTrackLengthLongestTrack([](const caf::SRSliceProxy *slc) -> double {
    int maxLengthID = kLongestTrack(slc);

    if(maxLengthID == std::numeric_limits<int>::signaling_NaN())
      return std::numeric_limits<double>::signaling_NaN();

    return slc->reco.pfp[maxLengthID].trk.len;
  });

const MultiVar kChildTrackLengthOtherTracks([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> trackLengths;

    int maxLengthID = kLongestTrack(slc);

    int i = -1;

    for(auto const& pfp : slc->reco.pfp)
      {
	++i;

	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	if(maxLengthID == i)
	  continue;

	trackLengths.push_back(pfp.trk.len);
      }

    return trackLengths;
  });

const MultiVar kChildTrackdEdxLongestTrack([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> dEdx;

    int maxLengthID = kLongestTrack(slc);

    if(maxLengthID == std::numeric_limits<int>::signaling_NaN())
      return dEdx;

    // Only interested in the collection plane (2)
    for(auto const& point : slc->reco.pfp[maxLengthID].trk.calo[2].points)
      dEdx.push_back(point.dedx);

    return dEdx;
  });

const MultiVar kChildTrackdEdxOtherTracks([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> dEdx;

    int maxLengthID = kLongestTrack(slc);

    int i = -1;

    for(auto const& pfp : slc->reco.pfp)
      {
	++i;

	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	if(maxLengthID == i)
	  continue;

	for(auto const& point : pfp.trk.calo[2].points)
	  dEdx.push_back(point.dedx);
      }

    return dEdx;
  });

const MultiVar kChildTrackResRangeLongestTrack([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> rr;

    int maxLengthID = kLongestTrack(slc);

    if(maxLengthID == std::numeric_limits<int>::signaling_NaN())
      return rr;

    // Only interested in the collection plane (2)
    for(auto const& point : slc->reco.pfp[maxLengthID].trk.calo[2].points)
      rr.push_back(point.rr);

    return rr;
  });

const MultiVar kChildTrackResRangeOtherTracks([](const caf::SRSliceProxy *slc) -> std::vector<double> {
    std::vector<double> rr;

    int maxLengthID = kLongestTrack(slc);

    int i = -1;

    for(auto const& pfp : slc->reco.pfp)
      {
	++i;

	if(!pfp.parent_is_primary || pfp.parent == -1)
	  continue;

	if(maxLengthID == i)
	  continue;

	for(auto const& point : pfp.trk.calo[2].points)
	  rr.push_back(point.rr);
      }

    return rr;
  });

const Var kOpT0Time([](const caf::SRSliceProxy *slc) -> double {
    return slc->opt0.time;
  });
