//-----------------------------------------------------------------------------
//
// PenPositron.hh
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#ifndef PenPositron_h
#define PenPositron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class PenPositron : public G4ParticleDefinition
{
 private:
   static PenPositron* theInstance;
   PenPositron(){}
   ~PenPositron(){}

 public:
   static PenPositron* Definition();
   static PenPositron* PositronDefinition();
   static PenPositron* Positron();
};


#endif



