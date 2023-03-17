void MatchingBasicPerformance()
{
  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt/clustering/crtana_v09_67_00_no_duplicates_no_showers.root");

  const TCut sp_cut = "tpc_length > 20";
  const TCut sp_cut_reco = sp_cut + "tpc_sp_score < 30";
  const TCut tr_cut = "tpc_length > 10";
  const TCut tr_cut_reco = tr_cut + "tpc_tr_score < 100";

  const int n_tracks = tree->Draw("tpc_truth_energy", "");
  const int n_tracks_sp = tree->Draw("tpc_truth_energy", sp_cut);
  const int n_tracks_sp_matchable = tree->Draw("tpc_truth_energy", "tpc_sp_matchable" + sp_cut);
  const int n_tracks_sp_matched = tree->Draw("tpc_truth_energy", "tpc_sp_matched" + sp_cut_reco);
  const int n_tracks_sp_matchable_matched = tree->Draw("tpc_truth_energy", "tpc_sp_matchable && tpc_sp_matched" + sp_cut_reco);
  const int n_tracks_sp_matched_good = tree->Draw("tpc_truth_energy", "tpc_sp_matched && tpc_sp_good_match" + sp_cut_reco);
  const int n_tracks_tr = tree->Draw("tpc_truth_energy", tr_cut);
  const int n_tracks_tr_matchable = tree->Draw("tpc_truth_energy", "tpc_tr_matchable" + tr_cut);
  const int n_tracks_tr_matched = tree->Draw("tpc_truth_energy", "tpc_tr_matched" + tr_cut_reco);
  const int n_tracks_tr_matchable_matched = tree->Draw("tpc_truth_energy", "tpc_tr_matchable && tpc_tr_matched" + tr_cut_reco);
  const int n_tracks_tr_matched_good = tree->Draw("tpc_truth_energy", "tpc_tr_matched && tpc_tr_good_match" + tr_cut_reco);

  std::cout << "\n----------------------------------------\n"
	    << "N Tracks: " << n_tracks << '\n'
	    << "N Tracks SP: " << n_tracks_sp << '\n'
	    << "N Tracks TR: " << n_tracks_tr << '\n'
	    << "----------------------------------------\n"
	    << "SP Matchable:           " << n_tracks_sp_matchable << Form(" (%.3f%%)", (100. * n_tracks_sp_matchable) / n_tracks_sp) << '\n'
	    << "SP Matched:             " << n_tracks_sp_matched << Form(" (%.3f%%)", (100. * n_tracks_sp_matched) / n_tracks_sp) << '\n'
	    << "SP Matchable & Matched: " << n_tracks_sp_matchable_matched << Form(" (%.3f%%)", (100. * n_tracks_sp_matchable_matched) / n_tracks_sp_matchable) << '\n'
	    << "SP Matched & Good:      " << n_tracks_sp_matched_good << Form(" (%.3f%%)", (100. * n_tracks_sp_matched_good) / n_tracks_sp_matched) << '\n'
	    << "----------------------------------------\n"
	    << "TR Matchable:           " << n_tracks_tr_matchable << Form(" (%.3f%%)", (100. * n_tracks_tr_matchable) / n_tracks_tr) << '\n'
	    << "TR Matched:             " << n_tracks_tr_matched << Form(" (%.3f%%)", (100. * n_tracks_tr_matched) / n_tracks_tr) << '\n'
	    << "TR Matchable & Matched: " << n_tracks_tr_matchable_matched << Form(" (%.3f%%)", (100. * n_tracks_tr_matchable_matched) / n_tracks_tr_matchable) << '\n'
	    << "TR Matched & Good:      " << n_tracks_tr_matched_good << Form(" (%.3f%%)", (100. * n_tracks_tr_matched_good) / n_tracks_tr_matched) << '\n'
	    << "----------------------------------------\n"
	    << std::endl;
}
