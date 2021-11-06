/**
 * @file SBND_QGSP_BERT_NNC
 *
 * @brief A Geant4 physics list. Same as QGSP_BERT,
 * but without neutron cut.
 *
 * @author Marco Del Tutto (mdeltutt@fnal.gov)
 *
 * A Geant4 physics list based on:
 * geant4/source/physics_lists/lists/src/QGSP_BERT.*
 * and removing the G4NeutronTrackingCut.
 *
 */

#ifndef SBND_QGSP_BERT_NNC_h
#define SBND_QGSP_BERT_NNC_h 1

#include "Geant4/globals.hh"
#include "Geant4/G4VModularPhysicsList.hh"

class SBND_QGSP_BERT_NNC: public G4VModularPhysicsList
{
public:
  SBND_QGSP_BERT_NNC(G4int ver = 1);
  virtual ~SBND_QGSP_BERT_NNC()=default;

  SBND_QGSP_BERT_NNC(const SBND_QGSP_BERT_NNC &) = delete;
  SBND_QGSP_BERT_NNC & operator=(const SBND_QGSP_BERT_NNC &)=delete;

};

#endif
