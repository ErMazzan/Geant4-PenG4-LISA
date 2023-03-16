//MACG 2019-04-26: Messenger class for PenPrimaryGeneratorAction
//                 It is meant to follow PENELOPE beam parametrization.


#ifndef PenPrimaryGeneratorMessenger_h
#define PenPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PenPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4UIcommand;

class PenPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  PenPrimaryGeneratorMessenger(PenPrimaryGeneratorAction*);
  ~PenPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);

private:
  void IonCommand(G4String);

 private:
  PenPrimaryGeneratorAction* fPenPrimary;
  // Directory
  G4UIdirectory* fPenSourceDir;
  // General commands
  G4UIcmdWithAnInteger* fSourceCmd;
  G4UIcmdWithAString*   fParticleCmd;
  G4UIcommand*          fIonCmd;
  // Energy commands
  G4UIcmdWithADoubleAndUnit*  fE0Cmd;
  // Position commands
  G4UIcmdWith3VectorAndUnit*  fPositCmd;
  // Angular commands
  G4UIcmdWith3VectorAndUnit*  fConeCmd;

};

#endif
