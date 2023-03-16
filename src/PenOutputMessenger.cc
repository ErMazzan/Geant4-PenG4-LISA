
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "PenOutputMessenger.hh"
#include "PenOutputPrinter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"


//-------1---------2---------3---------4---------5---------6---------7---------8
// ----------------------------------------------------------------------------

PenOutputMessenger::PenOutputMessenger( PenOutputPrinter* myOutput )
  : fOutput( myOutput )
{ 
  fOutputDir = new G4UIdirectory( "/penOutput/" );
  fOutputDir->SetGuidance("Controls for output files.");

  // IMPORTANT NOTE: Commands must be available in Idle state in order
  // to be issued by worker threads!!
  // If they are not available in Idle state, command is read but not issued!
  
  fRZDoseMapCmd = new G4UIcmdWithABool("/penOutput/createRZDoseMapFile", this);
  fRZDoseMapCmd->SetGuidance("Flag to create \"depth-dose.dat\" file");
  fRZDoseMapCmd->SetParameterName("choice", false);
  fRZDoseMapCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGridZminCmd = new G4UIcmdWithADoubleAndUnit("/penOutput/gridZmin", this);
  fGridZminCmd->SetGuidance("Set the grid-Z lower limit.");
  fGridZminCmd->SetParameterName("gridZmin", false);
  fGridZminCmd->SetUnitCategory("Length");
  fGridZminCmd->SetUnitCandidates("um mm cm m");
  fGridZminCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fGridZmaxCmd = new G4UIcmdWithADoubleAndUnit("/penOutput/gridZmax", this);
  fGridZmaxCmd->SetGuidance("Set the grid-Z upper limit.");
  fGridZmaxCmd->SetParameterName("gridZmax", false);
  fGridZmaxCmd->SetUnitCategory("Length");
  fGridZmaxCmd->SetUnitCandidates("um mm cm m");
  fGridZmaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fGridZnoBinsCmd = new G4UIcmdWithAnInteger("/penOutput/gridZnoBins", this);
  fGridZnoBinsCmd->SetGuidance("Set the grid-Z number of bins.");
  fGridZnoBinsCmd->SetParameterName("gridZnoBins", false);
  fGridZnoBinsCmd->SetRange("gridZnoBins > 0");
  fGridZnoBinsCmd->AvailableForStates(G4State_PreInit, G4State_Idle); 

  fGridRmaxCmd = new G4UIcmdWithADoubleAndUnit("/penOutput/gridRmax", this);
  fGridRmaxCmd->SetGuidance("Set the grid-R upper limit.");
  fGridRmaxCmd->SetParameterName("gridRmax", false);
  fGridRmaxCmd->SetUnitCategory("Length");
  fGridRmaxCmd->SetUnitCandidates("um mm cm m");
  fGridRmaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fGridRnoBinsCmd = new G4UIcmdWithAnInteger("/penOutput/gridRnoBins", this);
  fGridRnoBinsCmd->SetGuidance("Set the grid-R number of bins.");
  fGridRnoBinsCmd->SetParameterName("gridRnoBins", false);
  fGridRnoBinsCmd->SetRange("gridRnoBins > 0");
  fGridRnoBinsCmd->AvailableForStates(G4State_PreInit, G4State_Idle); 


  fSpcEnddetCmd = new G4UIcmdWithABool("/penOutput/create-spc-enddet", this);
  fSpcEnddetCmd->SetGuidance("Flag to create \"spc-enddet-##.dat\" file");
  fSpcEnddetCmd->SetParameterName("choice", false);
  fSpcEnddetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSpcEnddetNameCmd = new G4UIcmdWithAString("/penOutput/spc-enddet-PV", this);
  fSpcEnddetNameCmd->SetGuidance("Name of the PV to score Edep spectrum");
  fSpcEnddetNameCmd->SetParameterName("volName", false);
  fSpcEnddetNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSpcEnddetEndetcCmd = new G4UIcommand("/penOutput/spc-enddet-ENDETC", this);
  fSpcEnddetEndetcCmd->SetGuidance("set EL[eV] EU[eV] and Nbins");
  G4UIparameter* el =
    new G4UIparameter("Histogram Elow [eV]", 'd', false);
  fSpcEnddetEndetcCmd->SetParameter(el);
  //
  G4UIparameter* eu =
    new G4UIparameter("Histogram Eup [eV]", 'd', false);
  fSpcEnddetEndetcCmd->SetParameter(eu);
  //
  G4UIparameter* nBins =
    new G4UIparameter("Histogram Nbins", 'i', false);
  fSpcEnddetEndetcCmd->SetParameter(nBins);
  //
  fSpcEnddetEndetcCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


  fEmergingCmd = new G4UIcmdWithABool("/penOutput/create-emerging", this);
  fEmergingCmd->SetGuidance("Flag to create emerging particle files:");
  fEmergingCmd->SetGuidance(" - \"polar-angle.dat\"");
  fEmergingCmd->SetGuidance(" - \"energy-up.dat\" & \"energy-down.dat\"");
  fEmergingCmd->SetParameterName("choice", false);
  fEmergingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fPolarAngleNbThCmd = new G4UIcmdWithAnInteger("/penOutput/polar-angle-NBTH",
						this);
  fPolarAngleNbThCmd->SetGuidance("Number of bins in theta angle.");
  fPolarAngleNbThCmd->SetParameterName("NBTH", false);
  fPolarAngleNbThCmd->SetRange("NBTH > 0");
  fPolarAngleNbThCmd->AvailableForStates(G4State_PreInit, G4State_Idle); 
  
  fEmergElCmd = new G4UIcmdWithADoubleAndUnit("/penOutput/emerg-EL", this);
  fEmergElCmd->SetGuidance("Set E_min in \"energy-down\" & \"energy-up\".");
  fEmergElCmd->SetParameterName("Emin", false);
  fEmergElCmd->SetUnitCategory("Energy");
  fEmergElCmd->SetUnitCandidates("eV keV MeV");
  fEmergElCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fEmergEuCmd = new G4UIcmdWithADoubleAndUnit("/penOutput/emerg-EU", this);
  fEmergEuCmd->SetGuidance("Set E_max in \"energy-down\" & \"energy-up\".");
  fEmergEuCmd->SetParameterName("Emax", false);
  fEmergEuCmd->SetUnitCategory("Energy");
  fEmergEuCmd->SetUnitCandidates("eV keV MeV");
  fEmergEuCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEmergNbECmd = new G4UIcmdWithAnInteger("/penOutput/emerg-NBE", this);
  fEmergNbECmd->SetGuidance("No.bins in \"energy-down\" & \"energy-up\".");
  fEmergNbECmd->SetParameterName("NBE", false);
  fEmergNbECmd->SetRange("NBE > 0");
  fEmergNbECmd->AvailableForStates(G4State_PreInit, G4State_Idle); 
  
  // fTrkCmd = new G4UIcmdWithABool("/penOutput/createTrkFiles", this);
  // fTrkCmd->SetGuidance("Flag to create \".trk\" files");
  // fTrkCmd->SetParameterName("choice", false);
  // fTrkCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

// ----------------------------------------------------------------------------

PenOutputMessenger::~PenOutputMessenger()
{
  delete fOutputDir;

  delete fRZDoseMapCmd;
  delete fGridZminCmd;
  delete fGridZmaxCmd;
  delete fGridZnoBinsCmd;
  delete fGridRmaxCmd;
  delete fGridRnoBinsCmd;

  delete fSpcEnddetCmd;
  delete fSpcEnddetNameCmd;
  delete fSpcEnddetEndetcCmd;

  delete fEmergingCmd;
  delete fPolarAngleNbThCmd;
  delete fEmergElCmd;
  delete fEmergEuCmd;
  delete fEmergNbECmd;

  // delete fTrkCmd;
}

// ----------------------------------------------------------------------------

void PenOutputMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fRZDoseMapCmd) {
    fOutput->SetRZDoseMapFlag(fRZDoseMapCmd->GetNewBoolValue(newValue));
  }
  else if (command == fGridZminCmd) {
    fOutput->SetGridZmin(fGridZminCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fGridZmaxCmd) {
    fOutput->SetGridZmax(fGridZmaxCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fGridZnoBinsCmd) {
    fOutput->SetGridZnoBins(fGridZnoBinsCmd->GetNewIntValue(newValue));
  }
  else if (command == fGridRmaxCmd) {
    fOutput->SetGridRmax(fGridRmaxCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fGridRnoBinsCmd) {
    fOutput->SetGridRnoBins(fGridRnoBinsCmd->GetNewIntValue(newValue));
  }

  else if (command == fSpcEnddetCmd) {
    fOutput->SetSpcEnddetFlag(fSpcEnddetCmd->GetNewBoolValue(newValue));
  }
  else if (command == fSpcEnddetNameCmd) {
    fOutput->AddSpcEnddetPVName(newValue);
  }
  else if (command == fSpcEnddetEndetcCmd) {
    const char* strvaluelist= newValue.c_str();
    std::istringstream is(strvaluelist);
    G4double el, eu;
    G4int nbins;
    is >> el >> eu >> nbins;
    fOutput->AddSpcEnddetMin(el*eV);
    fOutput->AddSpcEnddetMax(eu*eV);
    fOutput->AddSpcEnddetNoBins(nbins);    
  }

  else if (command == fEmergingCmd) {
    fOutput->SetEmergingFlag( fEmergingCmd->GetNewBoolValue(newValue));
  }
  else if (command == fPolarAngleNbThCmd) {
    fOutput->SetPolarAngleNoBins(fPolarAngleNbThCmd->GetNewIntValue(newValue));
  }
  else if (command == fEmergElCmd) {
    fOutput->SetEmergEneMin(fEmergElCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fEmergEuCmd) {
    fOutput->SetEmergEneMax(fEmergEuCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fEmergNbECmd) {
    fOutput->SetEmergEneNoBins(fEmergNbECmd->GetNewIntValue(newValue));
  }

  // if (command == fTrkCmd) {
  //   fOutput->SetTrkFlag(fTrkCmd->GetNewBoolValue(newValue));
  // }
}
