//-----------------------------------------------------------------------------
//
// PenelopeEMProcess.hh - version 2021.05.07
//
// Geant4 physics constructor for defining Penelope physics
//
//-----------------------------------------------------------------------------

#ifndef PenelopeEMPhysics_h
#define PenelopeEMPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class PenelopeEMPhysicsMessenger;

class PenelopeEMPhysics : public G4VPhysicsConstructor
{
public:

  PenelopeEMPhysics(const G4String& name = "PenelopeEM");
  virtual ~PenelopeEMPhysics();

public:

  virtual void ConstructParticle();
  virtual void ConstructProcess();

public:  // Set methods
  void SetEthr(G4double value)     { fEthr = value; }
  void SetEpmax(G4double value)    { fEpmax = value; }

private:

   // hide assignment operator
  PenelopeEMPhysics & operator=(const PenelopeEMPhysics &right);
  PenelopeEMPhysics(const PenelopeEMPhysics&);

private: // Data members

  PenelopeEMPhysicsMessenger* fMessenger;

  G4double fEthr;
  G4double fEpmax;
};

#endif
