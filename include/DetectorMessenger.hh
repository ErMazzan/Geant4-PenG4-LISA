
#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

// ----------------------------------------------------------------------------

class DetectorMessenger: public G4UImessenger
{

public:

  DetectorMessenger( DetectorConstruction* );
  ~DetectorMessenger();

  void SetNewValue( G4UIcommand*, G4String );

private:

  DetectorConstruction*   fGeom;
  G4UIdirectory*             fGeomDir;
  G4UIdirectory*             fPenMatDir;

  //G4UIcmdWithADoubleAndUnit*  fCubeLCmd;
  //G4UIcmdWithADoubleAndUnit*  fDiscZCmd;
  //G4UIcommand*                fDiscDSMAXCmd;
  G4UIcommand*        fPenMatRegCmd;
  G4UIcommand*        fPenMatMSIMPACmd;
};

// ----------------------------------------------------------------------------

#endif
