
#ifndef PenOutputMessenger_h
#define PenOutputMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PenOutputPrinter;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

// ----------------------------------------------------------------------------

class PenOutputMessenger: public G4UImessenger
{

public:

  PenOutputMessenger( PenOutputPrinter* );
  ~PenOutputMessenger();

  void SetNewValue( G4UIcommand*, G4String );

private:

  PenOutputPrinter*          fOutput;
  G4UIdirectory*              fOutputDir;

  G4UIcmdWithABool*           fRZDoseMapCmd;
  G4UIcmdWithADoubleAndUnit*  fGridZminCmd;
  G4UIcmdWithADoubleAndUnit*  fGridZmaxCmd;
  G4UIcmdWithAnInteger*       fGridZnoBinsCmd;
  G4UIcmdWithADoubleAndUnit*  fGridRmaxCmd;
  G4UIcmdWithAnInteger*       fGridRnoBinsCmd;

  G4UIcmdWithABool*           fSpcEnddetCmd;
  G4UIcmdWithAString*         fSpcEnddetNameCmd;
  G4UIcommand*                fSpcEnddetEndetcCmd;

  G4UIcmdWithABool*           fEmergingCmd;
  G4UIcmdWithAnInteger*       fPolarAngleNbThCmd;
  G4UIcmdWithADoubleAndUnit*  fEmergElCmd;
  G4UIcmdWithADoubleAndUnit*  fEmergEuCmd;
  G4UIcmdWithAnInteger*       fEmergNbECmd;

  // G4UIcmdWithABool*           fTrkCmd;

};

// ----------------------------------------------------------------------------

#endif
