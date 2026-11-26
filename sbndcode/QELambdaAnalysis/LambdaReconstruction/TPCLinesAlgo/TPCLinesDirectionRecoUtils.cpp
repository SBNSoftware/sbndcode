////////////////////////////////////////////////////////////////////////////
//
// \file DirectionRecoUtils.h
//
// \brief Definition of direction reconstruction
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesDirectionRecoUtils.h"

namespace TPCLinesDirectionUtils{

    std::vector<SLinearCluster> SlopeTrackMerger(std::vector<SLinearCluster> trackList, double distTh, double slopeTh, int verbose){
        if(verbose>=1) std::cout << "\n\n+-+-+-+-+-+-+- Slope track merger +-+-+-+-+-+-+-\n";
        if(verbose>=1) std::cout << "NTracks: " << trackList.size() << '\n';

        std::vector<SLinearCluster> mergedTracks;
        bool foundNewMerge = false;

        for (size_t trkix = 0; trkix < trackList.size(); ++trkix) {
            bool merge_flag = false;
            SLinearCluster trk = trackList[trkix];
            if(verbose>=1) std::cout << trkix << " " << trk.GetMinX() << " " << trk.GetMaxX() << '\n';

            for (size_t mix = 0; mix < mergedTracks.size(); ++mix) {
                SLinearCluster mergedTrk = mergedTracks[mix];
                double minXDist = std::min(std::abs(mergedTrk.GetMinX() - trk.GetMaxX()), std::abs(mergedTrk.GetMaxX() - trk.GetMinX()));
                if(verbose>=1) std::cout << "  " << mergedTrk.GetMinX() << " " << mergedTrk.GetMaxX() << '\n';

                if (minXDist <= distTh) {
                    double dC1C2 = TPCLinesDistanceUtils::GetClusterConnectedness(trk.GetHitCluster(), mergedTrk.GetHitCluster());
                    if(verbose>=1) std::cout << "*** Merging  " << mergedTrk.GetMinX() << " " << mergedTrk.GetMaxX() << "  " << trk.GetMinX() << " " << trk.GetMaxX() << "   d=" << dC1C2 << "  MergedTrackComp " << mergedTrk.GetConnectedness() << '\n';

                    if (dC1C2 < distTh * mergedTrk.GetConnectedness()) {
                        double m1 = 0.0;
                        double m2 = 100.0;

                        if (trk.GetMaxX() <= mergedTrk.GetMinX()) {
                            m1 = trk.GetTrackEquationEnd().Slope();
                            m2 = mergedTrk.GetTrackEquationStart().Slope();
                        } else if (trk.GetMinX() >= mergedTrk.GetMaxX()) {
                            m1 = trk.GetTrackEquationStart().Slope();
                            m2 = mergedTrk.GetTrackEquationEnd().Slope();
                        }

                        double angle1 = std::atan(m1) * 180.0 / M_PI;
                        double angle2 = std::atan(m2) * 180.0 / M_PI;
                        double angle_between = std::abs(angle1 - angle2);
                        if(verbose>=1) std::cout << "            a1=" << angle1 << " a2=" << angle2 << " angle_diff=" << angle_between << '\n';

                        if (angle_between < slopeTh) {
                            if(verbose>=1) std::cout << "                  Merged!\n";

                            std::vector<SHit> newTrackHits = trk.GetHits();
                            std::vector<SHit> hitsToMerge = mergedTrk.GetHits();
                            newTrackHits.insert(newTrackHits.end(), hitsToMerge.begin(), hitsToMerge.end());

                            SLinearCluster newTrack(newTrackHits);
                            mergedTracks[mix] = newTrack;

                            merge_flag = true;
                            foundNewMerge = true;
                            break;
                        }
                    }
                }
            }

            if (!merge_flag) {
                mergedTracks.push_back(trk);
            }
        }

        if(verbose>=1) std::cout << "End loop continue " << foundNewMerge << '\n';

        if (foundNewMerge) {
            return SlopeTrackMerger(mergedTracks, distTh, slopeTh);
        }

        return mergedTracks;
    }


    bool GetLineHypoDistance(SLinearCluster mergeTrack, SLinearCluster track, int tol, int verbose) {
        
        LineEquation trackEq;
        SHit trackEdgeHit;
        LineEquation mergeTrackEq;
        SHit mergeTrackEdgeHit;

        if (mergeTrack.GetMeanX() > track.GetMeanX()) {
            trackEq = track.GetTrackEquationEnd();
            trackEdgeHit = track.GetHits().back();
            mergeTrackEq = mergeTrack.GetTrackEquationStart();
            mergeTrackEdgeHit = mergeTrack.GetHits().front();
        }
        else{
            trackEq = track.GetTrackEquationStart();
            trackEdgeHit = track.GetHits().front();
            mergeTrackEq = mergeTrack.GetTrackEquationEnd();
            mergeTrackEdgeHit = mergeTrack.GetHits().back();
        }

        double yHypoEdgeTrack = mergeTrackEq.Slope() * trackEdgeHit.X() + mergeTrackEq.Intercept();
        double yHypoEdgeMergeTrack = trackEq.Slope() * mergeTrackEdgeHit.X() + trackEq.Intercept();

        if(verbose>=1){
            std::cout << "     FreeSeg hit " << trackEdgeHit;
            std::cout << "     OverEdge hit " << mergeTrackEdgeHit;
            std::cout << "      YEdge=" << trackEdgeHit.Y() << " YWidth=" << trackEdgeHit.Width() << " HypoY=" << yHypoEdgeTrack << std::endl;
            std::cout << "      YEdgeMergeTrack=" << mergeTrackEdgeHit.Y() << " YWidth=" << mergeTrackEdgeHit.Width() << " HypoY=" << yHypoEdgeMergeTrack << std::endl;
        }
        
        return (std::abs(trackEdgeHit.Y() - yHypoEdgeTrack) < tol * trackEdgeHit.Width()) || (std::abs(mergeTrackEdgeHit.Y() - yHypoEdgeMergeTrack) < tol * mergeTrackEdgeHit.Width());
    }


    bool FullTrackContained(SLinearCluster mergeTrack, SLinearCluster track, double tol) {

        LineEquation mergeTrackEq = (mergeTrack.GetMeanX() > track.GetMeanX()) ?
                                            mergeTrack.GetTrackEquationStart() :
                                            mergeTrack.GetTrackEquationEnd();

        int nhits = 0;
        for (SHit hit : track.GetHits()) {
            double yHypo = mergeTrackEq.Slope() * hit.X() + mergeTrackEq.Intercept();
            if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
                nhits++;
            }
        }

        return nhits == track.NHits();
    }


    std::pair<bool, int> GetNHitsInHypo(SLinearCluster track1, SLinearCluster track2, int fNHits, double fWidthTol) {
        int nhits1 = std::min(fNHits, track1.NHits());
        std::vector<SHit> trk1Hits = std::vector<SHit>(track1.GetHits().end() - nhits1, track1.GetHits().end());
        SLinearCluster trk1(trk1Hits);
        //trk1.FitTrack();
        
        int nhits2 = std::min(fNHits, track2.NHits());
        std::vector<SHit> trk2Hits = std::vector<SHit>(track2.GetHits().begin(), track2.GetHits().begin() + nhits2);
        SLinearCluster trk2(trk2Hits);
        //trk2.FitTrack();
        
        std::cout << "Edge tracks " << trk1.GetYatMinX() << " " << trk1.GetMaxX() << " " << trk2.GetMinX() << " " << trk2.GetMaxX() << std::endl;
        
        int nhits1In = 0;
        for (SHit &hit : track1.GetHits()) {
            double yHypo = trk2.GetTrackEquation().Slope() * hit.X() + trk2.GetTrackEquation().Intercept();
            if (std::abs(yHypo - hit.Y()) < fWidthTol * hit.Width()) {
                nhits1In += 1;
            }
        }
        
        int nhits2In = 0;
        for (SHit &hit : track2.GetHits()) {
            double yHypo = trk1.GetTrackEquation().Slope() * hit.X() + trk1.GetTrackEquation( ).Intercept();
            if (std::abs(yHypo - hit.Y()) < fWidthTol * hit.Width()) {
                nhits2In += 1;
            }
        }
        
        int maxIn = std::max(nhits1In, nhits2In);
        return std::make_pair(maxIn > fNHits, maxIn);
    }


    bool check_overlap_and_find_region(int a1, int a2, int b1, int b2, int & overlapStart, int & overlapEnd) {
        overlapEnd=0;
        overlapStart=0;
        //if ((a1 <= b1 && b1 <= a2) || (a1 <= b2 && b2 <= a2)) {
        if( (a1 <= b1 && b1 <= a2) || (b1 <= a1 && a1 <= b2)) {
            int overlap_start = std::max(a1, b1);
            int overlap_end = std::min(a2, b2);
            
            if (overlap_start < overlap_end) {
                overlapStart = overlap_start;
                overlapEnd = overlap_end;
                return true;
            }
        }
        return false;
    }


    struct CompareSLinearClustersByMinX {
        bool operator()(SLinearCluster trk1, SLinearCluster trk2) const {
            return trk1.GetMinX() < trk2.GetMinX();
        }
    };


    std::vector<SLinearCluster> GetVertexTracks(
        std::vector<SLinearCluster> trackList,
        std::map<int, std::vector<int>>& shortToLongMap,
        std::map<int, int> & connectionsMap,
        int maxHitsShortTrack,
        float dTol,
        float connTolEps,
        int verbose
    ) {
        if(verbose>=1) std::cout << "\n+-+-+-+-+-+-+- Analyzing short tracks +-+-+-+-+-+-+-" << std::endl;
        
        // Get short tracks
        std::vector<SLinearCluster> shortTracks;
        for (SLinearCluster& track : trackList) {
            if (track.NHits() <= maxHitsShortTrack) {
                shortTracks.push_back(track);
            }
        }

        std::vector<SLinearCluster> shortTracks2;
        std::map<int, std::vector<int>> shortToLongDict;

        // Get short tracks with posterior long tracks
        for (SLinearCluster & sTrack : shortTracks) {
            std::vector<SLinearCluster> posteriorTracks;
            float sTrackMaxX = sTrack.GetMaxX();
            
            if(verbose>=1) std::cout << " Short Track Candidate ID: "<<sTrack.GetId() << std::endl;

            for (SLinearCluster lTrack : trackList) {
                if (lTrack.GetId() == sTrack.GetId()) {
                    continue;
                }

                if (lTrack.GetMinX() - sTrackMaxX <= dTol && lTrack.GetMinX() - sTrackMaxX > 0) {
                    float connTol = connTolEps * (sTrack.GetConnectedness() + lTrack.GetConnectedness()) / 2;
                    float conn = TPCLinesDistanceUtils::GetClusterConnectedness(sTrack.GetHitCluster(), lTrack.GetHitCluster());

                    if(verbose>=1) std::cout << "    cand long " << lTrack.GetId() << " " << lTrack.GetMinX() - sTrackMaxX << " Conn " << conn << " Tol " << connTol << std::endl;

                    if (conn < connTol) {
                        posteriorTracks.push_back(lTrack);
                    }
                }
            }

            if (posteriorTracks.size() == 2) {
                shortTracks2.push_back(sTrack);
                shortToLongDict[sTrack.GetId()] = { posteriorTracks[0].GetId(), posteriorTracks[1].GetId() };
            }
        }

        
        std::vector<SLinearCluster> finalShortTracks;
        shortToLongMap.clear();
        connectionsMap.clear();

        // Check slope and y position
        for (SLinearCluster& sTrack : shortTracks2) {

            SLinearCluster lTrack1 = trackList[shortToLongDict[sTrack.GetId()][0]];
            SLinearCluster lTrack2 = trackList[shortToLongDict[sTrack.GetId()][1]];
            
            float minSlope = std::min(lTrack1.GetTrackEquation().Slope(), lTrack2.GetTrackEquation().Slope());
            float maxSlope = std::max(lTrack1.GetTrackEquation().Slope(), lTrack2.GetTrackEquation().Slope());

            float minY = std::min(lTrack1.GetMeanY(), lTrack2.GetMeanY());
            float maxY = std::max(lTrack1.GetMeanY(), lTrack2.GetMeanY());
            float sTrackSlope = sTrack.GetTrackEquation().Slope();
            
            if(verbose>=1) std::cout << " Short Track Candidate ID: "<<sTrack.GetId() <<" , slope analysis: "<< std::endl;
            if(verbose>=1) std::cout << "minSlope " << minSlope << " maxSlope " << maxSlope << " shortSlope " << sTrackSlope << " Min/Max Y:"<< minY<<" "<<maxY << std::endl;
            //if(verbose>=1) std::cout << "minY " << minY << " maxY " << maxY << " shortSlope " << sTrack.GetMeanY() << std::endl;
            
            if (minSlope < sTrackSlope && sTrackSlope < maxSlope && minY){// < sTrack.GetMeanY() && sTrack.GetMeanY() < maxY) {
                finalShortTracks.push_back(sTrack);
                shortToLongMap[sTrack.GetId()] = { lTrack1.GetId(), lTrack2.GetId() };
                connectionsMap[lTrack1.GetId()] = lTrack2.GetId();
                connectionsMap[lTrack2.GetId()] = lTrack1.GetId();
            }
        }

        return finalShortTracks;
    }



    std::vector<std::vector<SLinearCluster>> GetParallelTracks(
        std::vector<SLinearCluster>& trackList, double fDist1DTh, double fAngleTh, double fMaxDWires, int verbose) {
        
        if(verbose>=1) std::cout << "\n+-+-+-+-+-+-+- Parallel track finder +-+-+-+-+-+-+-\n";

        if(trackList.size()<=1){
            return {trackList};
        }
        // Step 1: sort tracks by start point
        std::vector<SLinearCluster> sortedTracks = trackList;
        std::sort(sortedTracks.begin(), sortedTracks.end(), CompareSLinearClustersByMinX());

        // first give a socre for each pair of tracks
        // if one of the track is fully contained: score is 1/the angle between tracks
        // if not fully contained, score is the -angle between the tracks
        // the highest the socre, the better match
        std::map<double, std::vector<int>> scorePairs;
        std::map<int, int> trackIdToIndexDict;
        std::map<int, std::vector<int>> candidatesPerTrackId;
        //initialize
        for (size_t trkix = 0; trkix < sortedTracks.size(); ++trkix) {
            trackIdToIndexDict[sortedTracks[trkix].GetId()] = trkix;
            candidatesPerTrackId[sortedTracks[trkix].GetId()] = {};
        }

        // loop
        for (size_t trkix = 0; trkix < sortedTracks.size(); ++trkix) {
            
            SLinearCluster trk = sortedTracks[trkix];
            if(verbose>=1) std::cout << "\n *** Parallel tracks function... Track: " << trk.GetId()
                    << " min/max X = " << trk.GetMinX() << "/" << trk.GetMaxX() << "\n";
                    
            for (size_t trkjx = trkix + 1; trkjx < sortedTracks.size(); ++trkjx) {
                SLinearCluster trk2 = sortedTracks[trkjx];
                if(verbose>=1) std::cout << "    -- candidate track:" << trk2.GetId() 
                        << " min/max X:" << trk2.GetMinX() << " " << trk2.GetMaxX() 
                        << " Xdiff:" << trk2.GetMinX() - trk.GetMaxX() << "\n";

                // Consider the track if it starts after the track1
                if (trk2.GetMinX() - trk.GetMaxX() > fDist1DTh) {
                    
                    // Close enough in X
                    if (std::abs(trk2.GetMinX() - trk.GetMaxX()) > fMaxDWires) {
                        continue;
                    }
     
                    // Check the tracks match
                    bool hypoConnected = GetLineHypoDistance(trk, trk2);
                    if (!hypoConnected) {
                        continue;
                    }

                    // Check if track is fully contained
                    bool fullContained = false;
                    if (std::abs(trk2.GetMinX() - trk.GetMaxX()) <= 10) {
                        bool fullContained1 = FullTrackContained(trk, trk2);
                        bool fullContained2 = FullTrackContained(trk2, trk);
                        if(verbose>=1) std::cout << "       Full contained " << fullContained1 << " " << fullContained2 << "\n";
                        fullContained = fullContained1 || fullContained2;
                    }

                    // Check the angle compatibility
                    LineEquation trackEq1, trackEq2;
                    
                    // Get closest edges
                    if (trk2.GetMeanX() > trk.GetMeanX()) {
                        trackEq2 = trk2.GetTrackEquationStart();
                        trackEq1 = trk.GetTrackEquationEnd();
                    } else {
                        trackEq2 = trk2.GetTrackEquationEnd();
                        trackEq1 = trk.GetTrackEquationStart();
                    }

                    double angle1 = std::atan(trackEq1.Slope()) * 180.0 / M_PI;
                    double angle2 = std::atan(trackEq2.Slope()) * 180.0 / M_PI;
                    double slpCoM = (trk2.GetCoMPoint().Y() - trk.GetCoMPoint().Y()) / (trk2.GetCoMPoint().X() - trk.GetCoMPoint().X());
                    double angleCoM = std::atan(slpCoM) * 180.0 / M_PI;
                    double angle_between = std::abs(angle1 - angle2);
                    bool slopesCompatible = angle_between < fAngleTh && std::abs(angleCoM - angle1) < fAngleTh && std::abs(angleCoM - angle2) < fAngleTh;
                    
                    if(verbose>=1){
                        std::cout << "        angle1: " << angle1 << " angle2: " << angle2 
                            << " angleCOM: " << angleCoM << " angle_diff= " << angle_between << "\n";
                        std::cout << "        Slopes Compatible " << slopesCompatible << "\n";
                    }

                    if (fullContained) {
                        scorePairs[1.0 / angle_between] = {trk.GetId(), trk2.GetId()};
                        candidatesPerTrackId[trk.GetId()].push_back(trk2.GetId());
                        candidatesPerTrackId[trk2.GetId()].push_back(trk.GetId());
                    } else if (slopesCompatible) {
                        scorePairs[-angle_between] = {trk.GetId(), trk2.GetId()};
                        candidatesPerTrackId[trk.GetId()].push_back(trk2.GetId());
                        candidatesPerTrackId[trk2.GetId()].push_back(trk.GetId());
                    }
                }
            }
        }

        // Sort by max score
        std::vector<std::pair<double, std::vector<int>>> sortedScorePairs(scorePairs.begin(), scorePairs.end());
        std::sort(sortedScorePairs.begin(), sortedScorePairs.end(),
                [](const std::pair<double, std::vector<int>>& lhs,
                    const std::pair<double, std::vector<int>>& rhs) {
                    return lhs.first > rhs.first;
                });
        if(verbose>=1){
            std::cout << "Candidare parallel merge SCORES\n";
            for (const auto& pair : sortedScorePairs)
                std::cout << "Score " << pair.first << " Pair [" << pair.second[0] << ", " << pair.second[1] << "]\n";
        }

        // Perform the merging
        std::vector<std::vector<int>> finalTrackCluster;
        for (const auto& pair : sortedScorePairs){
            double score = pair.first;
            const std::vector<int>& ids = pair.second;
            if(verbose>=1) std::cout << "\n Score " << score << " Pair [" << ids[0] << ", " << ids[1] << "]\n";
            int id1 = ids[0];
            int id2 = ids[1];

            int id1_Ix = -1, id2_Ix = -1;
            for (size_t cix = 0; cix < finalTrackCluster.size(); ++cix) {
                const std::vector<int>& cluster = finalTrackCluster[cix];
                if (std::find(cluster.begin(), cluster.end(), id1) != cluster.end()) {
                    id1_Ix = static_cast<int>(cix);
                }
                if (std::find(cluster.begin(), cluster.end(), id2) != cluster.end()) {
                    id2_Ix = static_cast<int>(cix);
                }
            }

            std::vector<int> mergingTracksIx;
            int mergingClusterIx = -1;

            if (id1_Ix != -1 && id2_Ix != -1) {
                continue;
            } else if (id1_Ix != -1 && id2_Ix == -1) {
                mergingTracksIx.push_back(id2);
                mergingClusterIx = id1_Ix;
            } else if (id1_Ix == -1 && id2_Ix != -1) {
                mergingTracksIx.push_back(id1);
                mergingClusterIx = id2_Ix;
            } else {
                for (size_t cix = 0; cix < finalTrackCluster.size(); ++cix) {
                    const std::vector<int>& cluster = finalTrackCluster[cix];
                    for (int trkIx : cluster) {
                        if (std::find(candidatesPerTrackId[id1].begin(), candidatesPerTrackId[id1].end(), trkIx) != candidatesPerTrackId[id1].end() ||
                            std::find(candidatesPerTrackId[id2].begin(), candidatesPerTrackId[id2].end(), trkIx) != candidatesPerTrackId[id2].end()) {
                            mergingTracksIx.push_back(id1);
                            mergingTracksIx.push_back(id2);
                            mergingClusterIx = static_cast<int>(cix);
                        }
                    }
                }
            }

            // If no merging cluster found, create a new one
            if (mergingClusterIx == -1) {
                if(verbose>=1) std::cout << "  No previous cluster found, creating new cluster with IDs=" << id1 << " " << id2 << "\n";
                finalTrackCluster.push_back({id1, id2});
            }
            // If merging cluster found, check if the tracks overlap with the tracks in the cluster
            else {
                if(verbose>=1) std::cout << "  Previous cluster found, checking overlap region\n";
                std::vector<SLinearCluster> mergingTracks;
                for (int trkIx : mergingTracksIx) {
                    mergingTracks.push_back(sortedTracks[trackIdToIndexDict[trkIx]]);
                }

                bool overlaps = false;
                const std::vector<int>& cluster = finalTrackCluster[mergingClusterIx];
                for (int trackId : cluster) {
                    SLinearCluster trackInCluster = sortedTracks[trackIdToIndexDict[trackId]];
                    for (SLinearCluster mTrack : mergingTracks) {
                        int overlap_start, overlap_end;
                        bool overlap_exists = check_overlap_and_find_region(trackInCluster.GetMinX(), trackInCluster.GetMaxX(), mTrack.GetMinX(), mTrack.GetMaxX(), overlap_start, overlap_end);

                        if(verbose>=1) std::cout << "    Overlap with cluster track " << trackInCluster.GetId() << " for merging track " << mTrack.GetId() << " start/end " << overlap_start << " " << overlap_end << "\n";
                        if (overlap_exists) {
                            overlaps = true;
                            break;
                        }

                    }
                    if (overlaps) {
                        break;
                    }
                }

                if (!overlaps) {
                    
                    // if the score is negative, we have to check if the tracks are slopes compatible 
                    bool hasToMerge = score > 0;
                    if(score<0){
                        // cluster hits
                        std::vector<SHit> hitsInCluster;
                        for (int trackId : cluster) {
                            SLinearCluster trackInCluster = sortedTracks[trackIdToIndexDict[trackId]];
                            std::vector<SHit> auxHits = trackInCluster.GetHits();
                            hitsInCluster.insert(hitsInCluster.end(), auxHits.begin(), auxHits.end());
                        }
                        //potential merging hits
                        std::vector<SHit> mergingHits;
                        for (SLinearCluster mTrack : mergingTracks) {
                            std::vector<SHit> auxHits = mTrack.GetHits();
                            mergingHits.insert(mergingHits.end(), auxHits.begin(), auxHits.end());
                        }

                        SLinearCluster mergingTrack(mergingHits);
                        SLinearCluster mergingCluster(hitsInCluster);       
                        // Check the angle compatibility
                        LineEquation mergingTrackEq, clusterTrackEq;
                        
                        // Get closest edges
                        if (mergingTrack.GetMeanX() > mergingCluster.GetMeanX()) {
                            mergingTrackEq = mergingTrack.GetTrackEquationStart();
                            clusterTrackEq = mergingCluster.GetTrackEquationEnd();
                        } else {
                            mergingTrackEq = mergingTrack.GetTrackEquationEnd();
                            clusterTrackEq = mergingCluster.GetTrackEquationStart();
                        }

                        double angle1 = std::atan(mergingTrackEq.Slope()) * 180.0 / M_PI;
                        double angle2 = std::atan(clusterTrackEq.Slope()) * 180.0 / M_PI;
                        double slpCoM = (mergingCluster.GetCoMPoint().Y() - mergingTrack.GetCoMPoint().Y()) / (mergingCluster.GetCoMPoint().X() - mergingTrack.GetCoMPoint().X());
                        double angleCoM = std::atan(slpCoM) * 180.0 / M_PI;
                        double angle_between = std::abs(angle1 - angle2);
                        bool slopesCompatible = angle_between < fAngleTh && std::abs(angleCoM - angle1) < fAngleTh && std::abs(angleCoM - angle2) < fAngleTh;

                        hasToMerge = slopesCompatible;
                    }

                    if(hasToMerge){
                        if(verbose>=1){
                            std::cout << "         Merging tracks ";
                            for (int trkIx : mergingTracksIx) std::cout << trkIx << " ";
                            std::cout << "in cluster ";
                            for (int trackId : cluster) std::cout << trackId << " ";
                            std::cout << "\n";
                        }

                        finalTrackCluster[mergingClusterIx].insert(finalTrackCluster[mergingClusterIx].end(), mergingTracksIx.begin(), mergingTracksIx.end());
                    }

                } else {
                    finalTrackCluster.push_back(mergingTracksIx);
                }
            }
        }

        // Fill the vectors to return
        // matched parallel tracks
        std::vector<std::vector<SLinearCluster>> finalTrackClusterList;
        std::vector<std::vector<int>> finalTrackClusterIndexes = finalTrackCluster;
        for (const std::vector<int>& cluster : finalTrackCluster) {
            std::vector<SLinearCluster> trackList;
            for (int trkIx : cluster) {
                trackList.push_back(sortedTracks[trackIdToIndexDict[trkIx]]);
            }
            finalTrackClusterList.push_back(trackList);
        }
        // unmatched tracks
        std::vector<int> allParallelTracks;
        for (const std::vector<int>& sublist : finalTrackCluster) {
            allParallelTracks.insert(allParallelTracks.end(), sublist.begin(), sublist.end());
        }
        for (SLinearCluster track : sortedTracks) {
            if (std::find(allParallelTracks.begin(), allParallelTracks.end(), track.GetId()) == allParallelTracks.end()) {
                finalTrackClusterList.push_back({sortedTracks[trackIdToIndexDict[track.GetId()]]});
            }
        }

        return finalTrackClusterList;

    }



    SLinearCluster GetMainDirectionMaxHits(std::vector<std::vector<SLinearCluster>> trackClusterList, std::vector<SLinearCluster> & selectedTracksList, std::vector<SLinearCluster> & freeTracksList ) {
        size_t maxIndex = 0;
        int maxNHits = 0;
        for (size_t k = 0; k < trackClusterList.size(); ++k) {
            int nhits = 0;
            for (SLinearCluster track : trackClusterList[k]) {
                nhits += track.NHits();
            }
            if (nhits > maxNHits) {
                maxNHits = nhits;
                maxIndex = k;
            }
        }

        selectedTracksList.clear();
        freeTracksList.clear();

        std::vector<SHit> mainDirectionHits;
        for (SLinearCluster & track : trackClusterList[maxIndex]) {
            std::vector<SHit> hitsToAdd = track.GetHits();
            if(!mainDirectionHits.empty())
                mainDirectionHits.insert(mainDirectionHits.end(), hitsToAdd.begin(), hitsToAdd.end());
            else
                mainDirectionHits = std::vector<SHit>(hitsToAdd);

            selectedTracksList.push_back(track);
        }

        std::vector<SLinearCluster> freeTrackList;
        for (size_t k = 0; k < trackClusterList.size(); ++k) {
            if (k == maxIndex) continue;
            for (SLinearCluster& track : trackClusterList[k]) {
                freeTracksList.push_back(track);
            }
        }

        SLinearCluster mainDirection(mainDirectionHits);

        return mainDirection;
    }


    SLinearCluster GetMainDirectionLongest(std::vector<std::vector<SLinearCluster>> trackClusterList, std::vector<SLinearCluster> & selectedTracksList, std::vector<SLinearCluster> & freeTracksList, int minHits) {
        
        size_t maxIndex = 0;
        double maxTrackLength = 0.0;

        for (size_t k = 0; k < trackClusterList.size(); ++k) {
            int nhits = 0;
            for(SLinearCluster track : trackClusterList[k]) {
                nhits += track.NHits();
            }

            if (nhits > minHits) {
                double xStart = 1e6;
                double xEnd = -1e6;
                double yStart = 1e6;
                double yEnd = -1e6;

                for (SLinearCluster track : trackClusterList[k]) {
                    if (track.GetMinX() < xStart) {
                        xStart = track.GetMinX();
                    }
                    if (track.GetMaxX() > xEnd) {
                        xEnd = track.GetMaxX();
                    }
                    if (track.GetMinY() < yStart) {
                        yStart = track.GetMinY();
                    }
                    if (track.GetMaxY() > yEnd) {
                        yEnd = track.GetMaxY();
                    }
                }

                double trackLength = std::sqrt(0.3 * 0.3 * std::pow(xStart - xEnd, 2) + 0.08 * 0.08 * std::pow(yStart - yEnd, 2));
                
                if (trackLength > maxTrackLength) {
                    maxTrackLength = trackLength;
                    maxIndex = k;
                }
            }
        }

        selectedTracksList.clear();
        freeTracksList.clear();

        std::vector<SHit> mainDirectionHits;
        for (SLinearCluster & track : trackClusterList[maxIndex]) {
            std::vector<SHit> hitsToAdd = track.GetHits();
            if(!mainDirectionHits.empty())
                mainDirectionHits.insert(mainDirectionHits.end(), hitsToAdd.begin(), hitsToAdd.end());
            else
                mainDirectionHits = std::vector<SHit>(hitsToAdd);

            selectedTracksList.push_back(track);
        }

        std::vector<SLinearCluster> freeTrackList;
        for (size_t k = 0; k < trackClusterList.size(); ++k) {
            if (k == maxIndex) continue;
            for (SLinearCluster& track : trackClusterList[k]) {
                freeTracksList.push_back(track);
            }
        }

        SLinearCluster mainDirection(mainDirectionHits);

        return mainDirection;
    }


    SLinearCluster GetMainDirectionDownstream(std::vector<std::vector<SLinearCluster>> trackClusterList, std::vector<SLinearCluster> & selectedTracksList, std::vector<SLinearCluster> & freeTracksList, int minHits) {
        
        size_t maxIndex = 0;
        double minXStart = 1e6;

        for (size_t k = 0; k < trackClusterList.size(); ++k) {
            int nhits = 0;
            for (SLinearCluster track : trackClusterList[k]) {
                nhits += track.NHits();
            }

            if (nhits > minHits) {
                double xStart = 1e6;
                for (SLinearCluster track : trackClusterList[k]) {
                    if (track.GetMinX() < xStart) {
                        xStart = track.GetMinX();
                    }
                }
                
                if (xStart < minXStart) {
                    minXStart = xStart;
                    maxIndex = k;
                }
            }
        }

        selectedTracksList.clear();
        freeTracksList.clear();

        std::vector<SHit> mainDirectionHits;
        for (SLinearCluster & track : trackClusterList[maxIndex]) {
            std::vector<SHit> hitsToAdd = track.GetHits();
            if(!mainDirectionHits.empty())
                mainDirectionHits.insert(mainDirectionHits.end(), hitsToAdd.begin(), hitsToAdd.end());
            else
                mainDirectionHits = std::vector<SHit>(hitsToAdd);

            selectedTracksList.push_back(track);
        }

        std::vector<SLinearCluster> freeTrackList;
        for (size_t k = 0; k < trackClusterList.size(); ++k) {
            if (k == maxIndex) continue;
            for (SLinearCluster& track : trackClusterList[k]) {
                freeTracksList.push_back(track);
            }
        }

        SLinearCluster mainDirection(mainDirectionHits);

        return mainDirection;
    }

}

