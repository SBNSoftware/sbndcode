#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void AnaDump()
{
  TChain *tree = new TChain("crtana/tree");
  tree->Add("/sbnd/data/users/hlay/crt_top_high/production/crttophighana.root");

  TCut baseCut = "abs(pdg)==13";

  std::cout << "======================================" << std::endl;
  std::cout << "ANADUMP FOR CRT TOP HIGH INVESTIGATION" << std::endl;
  std::cout << "======================================" << std::endl << std::endl;

  std::cout << "Base Cut: " << baseCut << std::endl << std::endl;

  int totalTrue = tree->GetEntries(baseCut);
  int totalTrueEnterExitTPC = tree->GetEntries(baseCut + "entersOrExits");
  int totalTrueEnterExitTPCCRT = tree->GetEntries(baseCut + "entersOrExits && entersOrExitsCRT");

  TCut passesThrough = baseCut + "entersOrExits && entersOrExitsCRT";

  int totalTrueEnterExitTPCCRTHaveTrack = tree->GetEntries(passesThrough + "hasTPCTrack");
  int totalTrueEnterExitTPCCRTHaveLongTrack = tree->GetEntries(passesThrough + "hasTPCLongTrack");

  int intersectsAny = tree->GetEntries(passesThrough + "intersectsAny");
  int intersectsAnyStrips = tree->GetEntries(passesThrough + "intersectsAnyStrips");
  int intersectsAnyNoBottom = tree->GetEntries(passesThrough + "intersectsAnyNoBottom");
  int intersectsAnyStripsNoBottom = tree->GetEntries(passesThrough + "intersectsAnyStripsNoBottom");

  TCut otherCRTStrips = "intersectsBottomStrips || intersectsEastStrips || intersectsWestStrips || intersectsNorthStrips || intersectsSouthStrips";
  
  int intersectsOtherCRTStrips = tree->GetEntries(passesThrough + otherCRTStrips);
  int intersectsOtherPlusNormal = tree->GetEntries(passesThrough + (otherCRTStrips || "intersectsTopLowStrips"));
  int intersectsOtherPlusExtended = tree->GetEntries(passesThrough + (otherCRTStrips || "intersectsExtendedTopLow"));
  int intersectsOtherPlusRaised = tree->GetEntries(passesThrough + (otherCRTStrips || "intersectsRaisedTopLow"));
  int intersectsOtherPlusExtendedRaised = tree->GetEntries(passesThrough + (otherCRTStrips || "intersectsExtendedRaisedTopLow"));
  int intersectsOtherPlusExtendedRaisedFullNorth = tree->GetEntries(passesThrough + (otherCRTStrips || "intersectsExtendedRaisedTopLowFullNorth"));
  int intersectsOtherPlusExtendedRaisedFullSouth = tree->GetEntries(passesThrough + (otherCRTStrips || "intersectsExtendedRaisedTopLowFullSouth"));

  TCut createsTrack = passesThrough + "hasTPCLongTrack";
  TCut hitable = createsTrack + "intersectsBottomStrips || intersectsEastStrips || intersectsWestStrips || intersectsNorthStrips || intersectsSouthStrips || intersectsTopLowStrips || intersectsTopHighStrips";
  int isHitable = tree->GetEntries(hitable);
  int hasHitMatch = tree->GetEntries(hitable + "hasHitMatch");
  int hasGoodHitMatch = tree->GetEntries(hitable + "hasGoodHitMatch");
  int hitMatchTopHigh = tree->GetEntries(hitable + "hasGoodHitMatch && hitMatchTagger==6");

  TCut trackable = createsTrack + "(intersectsBottomStrips + intersectsEastStrips + intersectsWestStrips + intersectsNorthStrips + intersectsSouthStrips + intersectsTopLowStrips + intersectsTopHighStrips) > 1";
  int isTrackable = tree->GetEntries(trackable);
  int hasTrackMatch = tree->GetEntries(trackable + "hasTrackMatch");
  int hasGoodTrackMatch = tree->GetEntries(trackable + "hasGoodTrackMatch");
  int trackInAbsenceOfHit = tree->GetEntries(trackable + "hasTrackMatch && !hasHitMatch");
  int goodTrackInAbsenceOfGoodHit = tree->GetEntries(trackable + "hasGoodTrackMatch && !hasGoodHitMatch");

  TCut telescope = trackable + "hasGoodTrackMatch && (trackMatchTagger1 == 5 || trackMatchTagger2 == 5 || trackMatchTagger3 == 5) && (trackMatchTagger1 == 6 || trackMatchTagger2 == 6 || trackMatchTagger3 == 6)";
  int telescoped = tree->GetEntries(telescope);
  TCut telescopeOnly = telescope + "(trackMatchTagger1 == -1|| trackMatchTagger2 == -1 || trackMatchTagger3 == -1)";
  int telescopedOnly = tree->GetEntries(telescopeOnly);
  
  std::cout << "Total number of true particles considered:      " << totalTrue << '\n'
            << "of which enter or exit TPC:                     " << totalTrueEnterExitTPC << '\n'
            << "of which enter or exit CRT:                     " << totalTrueEnterExitTPCCRT << '\n'
            << "of which have a truth matched TPC track:        " << totalTrueEnterExitTPCCRTHaveTrack << '\n'
            << "of which are > 5cm:                             " << totalTrueEnterExitTPCCRTHaveLongTrack << '\n'
            << '\n'
            << "of which intersects any CRT simple rectangle:   " << intersectsAny << " (" <<  intersectsAny * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << "       - no bottom:                             " << intersectsAnyNoBottom << " (" <<  intersectsAnyNoBottom * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << "of which intersects any CRT strip areas:        " << intersectsAnyStrips << " (" <<  intersectsAnyStrips * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << "       - no bottom:                             " << intersectsAnyStripsNoBottom << " (" <<  intersectsAnyStripsNoBottom * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << '\n'
            << "of which intersects 'other' CRT strips (BNESW): " << intersectsOtherCRTStrips << " (" <<  intersectsOtherCRTStrips * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << " + normal top low:                              " << intersectsOtherPlusNormal << " (" <<  intersectsOtherPlusNormal * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << " + extended:                                    " << intersectsOtherPlusExtended << " (" <<  intersectsOtherPlusExtended * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << " + raised:                                      " << intersectsOtherPlusRaised << " (" <<  intersectsOtherPlusRaised * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << " + extended + raised:                           " << intersectsOtherPlusExtendedRaised << " (" <<  intersectsOtherPlusExtendedRaised * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << " + extended + raised (full north):              " << intersectsOtherPlusExtendedRaisedFullNorth << " (" <<  intersectsOtherPlusExtendedRaisedFullNorth * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << " + extended + raised (full south):              " << intersectsOtherPlusExtendedRaisedFullSouth << " (" <<  intersectsOtherPlusExtendedRaisedFullSouth * 100.f / totalTrueEnterExitTPCCRT << "%)" << '\n'
            << '\n'
            << std::endl;

  std::cout << "Total number of true particles considered:     " << totalTrue << '\n'
            << "of which enter or exit TPC:                    " << totalTrueEnterExitTPC << '\n'
            << "of which enter or exit CRT:                    " << totalTrueEnterExitTPCCRT << '\n'
            << "of which have a truth matched TPC track:       " << totalTrueEnterExitTPCCRTHaveTrack << '\n'
            << "of which are > 5cm:                            " << totalTrueEnterExitTPCCRTHaveLongTrack << '\n'
            << "of which are hitable:                          " << isHitable << '\n'
            << "of which have a hit match:                     " << hasHitMatch << " (" << hasHitMatch * 100.f / isHitable << "%)" << '\n'
            << " - good:                                       " << hasGoodHitMatch << " (" << hasGoodHitMatch * 100.f / isHitable << "%)" << '\n'
            << " - good using top high:                        " << hitMatchTopHigh << " (" << hitMatchTopHigh * 100.f / isHitable << "%)" << '\n'
            << "of which are trackable:                        " << isTrackable << '\n'
            << "of which have a track match:                   " << hasTrackMatch << " (" << hasTrackMatch * 100.f / isTrackable << "%)" << '\n'
            << " - good:                                       " << hasGoodTrackMatch << " (" << hasGoodTrackMatch * 100.f / isTrackable << "%)" << '\n'
            << " - good using telescope:                       " << telescoped << " (" << telescoped * 100.f / isTrackable << "%)" << '\n'
            << " - good using telescope ONLY:                  " << telescopedOnly << " (" << telescopedOnly * 100.f / isTrackable << "%)" << '\n'
            << "track in absence of hit:                       " << trackInAbsenceOfHit << " (" << trackInAbsenceOfHit * 100.f / isTrackable << "%)" << '\n'
            << "good track in absence of good hit:             " << goodTrackInAbsenceOfGoodHit << " (" << goodTrackInAbsenceOfGoodHit * 100.f / isTrackable << "%)" << '\n'
            << std::endl;
}
