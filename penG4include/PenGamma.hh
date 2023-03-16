//-----------------------------------------------------------------------------
//
// PenGamma.hh
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#ifndef PenGamma_h
#define PenGamma_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

class PenGamma : public G4ParticleDefinition
{
 private:
   static PenGamma* theInstance;

 private: // hide constructor as private
   PenGamma(){}

 public:
   ~PenGamma(){}

   static PenGamma* Definition();
   static PenGamma* GammaDefinition();
   static PenGamma* Gamma();
};


#endif


