
#include "globals.hh"

#include "PenelopeEMPhysicsMessenger.hh"
#include "PenelopeEMPhysics.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


//-------1---------2---------3---------4---------5---------6---------7---------8
// ----------------------------------------------------------------------------

PenelopeEMPhysicsMessenger::
PenelopeEMPhysicsMessenger( PenelopeEMPhysics* penEMPhys )
  : fPenEMPhys( penEMPhys )
{
  fPenEMPhysDir = new G4UIdirectory( "/penEMPhys/" );
  fPenEMPhysDir->SetGuidance("Configuration of PenelopeEMPhysics.");

  fEthrCmd = new G4UIcmdWithADoubleAndUnit("/penEMPhys/Ethr", this);
  fEthrCmd->SetGuidance("Set the energy threshold below which a Geant4");
  fEthrCmd->SetGuidance("particle is converted to a Penelope particle:");
  fEthrCmd->SetGuidance("   e-    --> pe-");
  fEthrCmd->SetGuidance("   gamma --> pgamma");
  fEthrCmd->SetGuidance("   e+    --> pe+");
  fEthrCmd->SetParameterName("Ethr", false);
  fEthrCmd->SetUnitCategory("Energy");
  fEthrCmd->SetRange("Ethr > 0.0");
  fEthrCmd->SetUnitCandidates("eV keV MeV GeV");
  fEthrCmd->AvailableForStates(G4State_PreInit);

  fEpmaxCmd = new G4UIcmdWithADoubleAndUnit("/penEMPhys/Epmax", this);
  fEpmaxCmd->SetGuidance("Set the value of maximum particle energy for ");
  fEpmaxCmd->SetGuidance("interpolation of material tables. For primary...");
  fEpmaxCmd->SetGuidance("   * e- or gamma, Epmax must be >= beam Emax.");
  fEpmaxCmd->SetGuidance("   * e, Epmax must be >= 1.21*(Emax+me*c^2), which");
  fEpmaxCmd->SetGuidance("     is the maximum E of annihilation photons.");
  fEpmaxCmd->SetParameterName("Epmax", false);
  fEpmaxCmd->SetUnitCategory("Energy");
  fEpmaxCmd->SetRange("Epmax > 0.0");
  fEpmaxCmd->SetUnitCandidates("eV keV MeV GeV");
  fEpmaxCmd->AvailableForStates(G4State_PreInit);

}

// ----------------------------------------------------------------------------

PenelopeEMPhysicsMessenger::~PenelopeEMPhysicsMessenger()
{
  delete fPenEMPhysDir;

  delete fEthrCmd;
  delete fEpmaxCmd;
}

// ----------------------------------------------------------------------------

void PenelopeEMPhysicsMessenger::SetNewValue(G4UIcommand* command,
					     G4String newValue)
{ 
  if (command == fEthrCmd) {
    fPenEMPhys->SetEthr(fEthrCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fEpmaxCmd) {
    fPenEMPhys->SetEpmax(fEpmaxCmd->GetNewDoubleValue(newValue));
  }
}
