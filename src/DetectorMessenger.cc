
#include "globals.hh"

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


//-------1---------2---------3---------4---------5---------6---------7---------8
// ----------------------------------------------------------------------------

DetectorMessenger::DetectorMessenger( DetectorConstruction* myDet )
  : fGeom( myDet )
{
  fGeomDir = new G4UIdirectory( "/mygeom/" );
  fGeomDir->SetGuidance("Geometry parameters.");

  fPenMatDir = new G4UIdirectory( "/mygeom/peng4mat/" );
  fPenMatDir->SetGuidance("Setup for Penelope materials.");

  fPenMatRegCmd = new G4UIcommand("/mygeom/peng4mat/register", this);
  fPenMatRegCmd->SetGuidance("Register Penelope material for Disc volume.");
  fPenMatRegCmd->SetGuidance("[See PENELOPE manual for details]");
  //
  G4UIparameter* matName =
    new G4UIparameter("Material name (only NIST and user registered ones)",
		      's', false);
  fPenMatRegCmd->SetParameter(matName);
  //
  G4UIparameter* matPenID =
    new G4UIparameter("Material ID in Penelope database (0 for a G4 material)",
		      'i', false);
  fPenMatRegCmd->SetParameter(matPenID);
  //
  fPenMatRegCmd->AvailableForStates(G4State_PreInit);

  fPenMatMSIMPACmd = new G4UIcommand("/mygeom/peng4mat/MSIMPA", this);
  fPenMatMSIMPACmd->SetGuidance("Set MSIMPA parameters for a PenG4 mat. ");
  fPenMatMSIMPACmd->SetGuidance("Order: EABS[1-3] C1 C2 WCC WCR. ");  
  fPenMatMSIMPACmd->SetGuidance("This command should be issued as many times ");
  fPenMatMSIMPACmd->SetGuidance("as /mygeom/peng4mat/register. ");
  fPenMatMSIMPACmd->SetGuidance("Correspondence between PenG4 mat and MSIMPA ");
  fPenMatMSIMPACmd->SetGuidance("is done by order of command issuance.");
  //
  G4UIparameter* eabs1 =
    new G4UIparameter("EABS[1] value (Eabs for e- in eV)", 'd', false);
  fPenMatMSIMPACmd->SetParameter(eabs1);
  //
  G4UIparameter* eabs2 =
    new G4UIparameter("EABS[2] value (Eabs for gamma in eV)", 'd', false);
  fPenMatMSIMPACmd->SetParameter(eabs2);
  //
  G4UIparameter* eabs3 =
    new G4UIparameter("EABS[3] value (Eabs for e+ in eV)", 'd', false);
  fPenMatMSIMPACmd->SetParameter(eabs3);
  //
  G4UIparameter* c1 = new G4UIparameter("C1 value", 'd', false);
  fPenMatMSIMPACmd->SetParameter(c1);
  //
  G4UIparameter* c2 = new G4UIparameter("C2 value", 'd', false);
  fPenMatMSIMPACmd->SetParameter(c2);
  //
  G4UIparameter* wcc = new G4UIparameter("WCC value (in eV)", 'd', false);
  fPenMatMSIMPACmd->SetParameter(wcc);
  //
  G4UIparameter* wcr = new G4UIparameter("WCR value (in eV)", 'd', false);
  fPenMatMSIMPACmd->SetParameter(wcr);
  //
  fPenMatMSIMPACmd->AvailableForStates(G4State_PreInit);

}

// ----------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{
  delete fGeomDir;
  delete fPenMatDir;

  delete fPenMatRegCmd;
  delete fPenMatMSIMPACmd;
}

// ----------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 

  if (command == fPenMatRegCmd) {
    const char* strvaluelist= newValue.c_str();
    std::istringstream is(strvaluelist);
    G4String matName;
    G4int matID;
    is >> matName >> matID;
    fGeom->RegisterPenG4Mat(matName, matID);
  }
  else if (command == fPenMatMSIMPACmd) {
    const char* strvaluelist= newValue.c_str();
    std::istringstream is(strvaluelist);
    G4double eabs1, eabs2, eabs3, c1, c2, wcc, wcr;
    is >> eabs1 >> eabs2 >> eabs3 >> c1 >> c2 >> wcc >> wcr;
    fGeom->SetPenG4MatMSIMPA(eabs1, eabs2, eabs3, c1, c2, wcc, wcr);
  }
}
