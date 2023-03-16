
#ifndef PenelopeEMPhysicsMessenger_h
#define PenelopeEMPhysicsMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PenelopeEMPhysics;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

// ----------------------------------------------------------------------------

class PenelopeEMPhysicsMessenger: public G4UImessenger
{

public:

  PenelopeEMPhysicsMessenger( PenelopeEMPhysics* );
  ~PenelopeEMPhysicsMessenger();

  void SetNewValue( G4UIcommand*, G4String );

private:

  PenelopeEMPhysics*          fPenEMPhys;
  G4UIdirectory*              fPenEMPhysDir;

  G4UIcmdWithADoubleAndUnit*  fEthrCmd;
  G4UIcmdWithADoubleAndUnit*  fEpmaxCmd;
};

// ----------------------------------------------------------------------------

#endif
