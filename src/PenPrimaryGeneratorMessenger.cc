//MACG 2019-04-26: Messenger class for PenPrimaryGeneratorAction
//                 It is meant to follow PENELOPE beam parametrization.


#include "PenPrimaryGeneratorMessenger.hh"

#include "PenPrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"


PenPrimaryGeneratorMessenger::
PenPrimaryGeneratorMessenger(PenPrimaryGeneratorAction* penPrimary)
:fPenPrimary(penPrimary)
{ 
  //
  // Definition of the interactive commands to modify the parameters of the
  // generation of primary particles
  // 
  fPenSourceDir = new G4UIdirectory("/penSource/");
  fPenSourceDir->SetGuidance("Set properties of the beam with a logic");
  fPenSourceDir->SetGuidance("similar to that of PENMAIN.");

  // =================
  // GENERAL COMMANDS
  // =================

  fSourceCmd = new G4UIcmdWithAnInteger("/penSource/source", this);
  fSourceCmd->SetGuidance("Set the type of radiation source based on a ");
  fSourceCmd->SetGuidance("3-digit number:");
  fSourceCmd->SetGuidance(" - 1st digit (hundreds): Energy Distribution.");
  fSourceCmd->SetGuidance("    0 = mono or gaussian.");
  fSourceCmd->SetGuidance(" - 2nd digit (tens): Position Distribution.");
  fSourceCmd->SetGuidance("    0 = point or 2D-Gaussian distribution.");
  fSourceCmd->SetGuidance(" - 3rd digit (units): Angular Distribution.");
  fSourceCmd->SetGuidance("    0 = uniform conical distribution.");
  fSourceCmd->SetGuidance(" Examples:");
  fSourceCmd->SetGuidance(" - Default value = 000");
  fSourceCmd->SetParameterName("choice", false);
  fSourceCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


  fParticleCmd = new G4UIcmdWithAString("/penSource/particle", this);
  fParticleCmd -> SetGuidance("Type of primary particle");
  fParticleCmd -> SetParameterName("choice", false);
  fParticleCmd -> AvailableForStates(G4State_PreInit, G4State_Idle);

  fIonCmd = new G4UIcommand("/penSource/ion",this);
  fIonCmd->SetGuidance("Set an ion as beam particle.");
  fIonCmd->SetGuidance("[usage] /penSource/ion Z A E");
  fIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  fIonCmd->SetGuidance("        A:(int) AtomicMass");
  fIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");
  G4UIparameter* param;
  param = new G4UIparameter("Z",'i',false);
  param->SetDefaultValue("1");
  fIonCmd->SetParameter(param);
  param = new G4UIparameter("A",'i',false);
  param->SetDefaultValue("1");
  fIonCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  fIonCmd->SetParameter(param);
  fIonCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  //=================
  // ENERGY COMMANDS
  //=================

  fE0Cmd = new G4UIcmdWithADoubleAndUnit("/penSource/energ", this);
  fE0Cmd->SetGuidance("Set the source mean kinetic energy");
  fE0Cmd->SetUnitCandidates("MeV keV eV");
  fE0Cmd->SetParameterName("E0", false);
  fE0Cmd->SetRange("E0>0.0");
  fE0Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  
  //===================
  // POSITION COMMANDS
  //===================

  fPositCmd = new G4UIcmdWith3VectorAndUnit("/penSource/posit", this);
  fPositCmd->SetGuidance("Set the coordinates of the source");
  fPositCmd->SetUnitCandidates("km m cm mm um");
  fPositCmd->SetParameterName("X0","Y0","Z0",false);
  fPositCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


  //==================
  // ANGULAR COMMANDS
  //==================

  fConeCmd = new G4UIcmdWith3VectorAndUnit("/penSource/cone", this);
  fConeCmd->SetGuidance("Set uniform conical beam around (theta,phi)");
  fConeCmd->SetGuidance("direction with aperture alpha.");
  fConeCmd->SetUnitCandidates("deg rad");
  fConeCmd->SetParameterName("theta","phi","alpha",false);
  fConeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


}


PenPrimaryGeneratorMessenger::~PenPrimaryGeneratorMessenger()
{
  // Directories
  delete fPenSourceDir;
  // General
  delete fSourceCmd;
  delete fParticleCmd;
  delete fIonCmd;
  // Energy
  delete fE0Cmd;
  // Position
  delete fPositCmd;
  // Angular
  delete fConeCmd;
}  


void PenPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* cmd,
					       G4String newValue)
{
  //=================
  // GENERAL COMMANDS
  //=================

  if ( cmd == fSourceCmd )
    { fPenPrimary->SetSource(fSourceCmd->GetNewIntValue(newValue)); }

  else if ( cmd == fParticleCmd )
    { fPenPrimary->SetParticle(newValue); }

  else if ( cmd == fIonCmd )
    { IonCommand(newValue); }

  //=================
  // ENERGY COMMANDS
  //=================

  else if ( cmd == fE0Cmd )
    { fPenPrimary->SetMeanKinEne(fE0Cmd->GetNewDoubleValue(newValue)); }


  //===================
  // POSITION COMMANDS
  //===================

  else if ( cmd == fPositCmd )
    { fPenPrimary->SetPosition(fPositCmd->GetNew3VectorValue(newValue)); }


  //==================
  // ANGULAR COMMANDS
  //==================

  else if ( cmd == fConeCmd )
    { fPenPrimary->SetCone(fConeCmd->GetNew3VectorValue(newValue)); }

}



void PenPrimaryGeneratorMessenger::IonCommand(G4String newValues)
{
  // Read Command
  // -------------

  G4Tokenizer next( newValues );
  // check argument
  G4int Z = StoI(next());
  G4int A = StoI(next());
  G4double eE = 0.0;

  G4String sE = next();
  if (!(sE.isNull()))
    eE = StoD(sE) * keV;

  // Call the suitable method
  // -------------------------
  fPenPrimary->SetIon(Z, A, eE);
}
