//-----------------------------------------------------------------------------
//
// PenElectron.hh
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#ifndef PenElectron_h
#define PenElectron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class PenElectron : public G4ParticleDefinition
{
 private:
   static PenElectron* theInstance;
   PenElectron(){}
   ~PenElectron(){}

 public:
   static PenElectron* Definition();
   static PenElectron* ElectronDefinition();
   static PenElectron* Electron();
};

#endif











