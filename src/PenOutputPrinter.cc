 
#include "PenOutputPrinter.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4ExceptionHandler.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Navigator.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"

#include "DetectorConstruction.hh"
#include "PenOutputMessenger.hh"
#include "Run.hh"

#include "PenElectron.hh"
#include "PenGamma.hh"
#include "PenPositron.hh"

#include <cstring>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>


//-------1---------2---------3---------4---------5---------6---------7---------8
PenOutputPrinter::PenOutputPrinter()
{
  G4cout << "void PenOutputPrinter is not meant to be used!" << G4endl;
}


//-------1---------2---------3---------4---------5---------6---------7---------8
PenOutputPrinter::PenOutputPrinter(const DetectorConstruction *detConst)
{
  fDetConst = detConst;

  fTrkFlag = false;

  // 2d-dose-map.dat file parameters
  fRZDoseMapFlag = false;
  fGridZmin = 0.000*cm;
  fGridZmax = 0.005*cm;
  fGridZnoBins = 100;
  fGridRmax = 0.01*cm;
  fGridRnoBins = 50;
  fGridBinDZ = (fGridZmax-fGridZmin)/fGridZnoBins;
  fGridBinDR = (fGridRmax)/fGridRnoBins;
  fGridRZDoseEvt = 0;
  fDepthDoseEvt = 0;

  // spc-enddet-##.dat file parameters
  // null is default
  fSpcEnddetFlag = false;
  fNoSpcEnddet = 0;
  fSpcEnddetFileVec = new std::vector<std::ofstream*>;
  fSpcEnddetNameVec = new std::vector<G4String>;
  fSpcEnddetEvtVec = new std::vector<G4double>;
  fSpcEnddetMinVec = new std::vector<G4double>;
  fSpcEnddetMaxVec = new std::vector<G4double>;
  fSpcEnddetNoBinsVec = new std::vector<G4int>;

  // files for emerging particles
  // ----------------------------
  fEmergingFlag = true;
  // polar-angle.dat file parameters
  fPolarAngleNoBins = 90;
  fPolarAngleBinDT = (180.*deg)/fPolarAngleNoBins;
  fPolarAngleElectEvt = 0;
  fPolarAngleGammaEvt = 0;
  fPolarAnglePositEvt = 0;
  // energy-up.dat & energy-down.dat parameters
  fEmergEneMin =  0.0*eV;
  fEmergEneMax = 40.000001*keV;
  fEmergEneNoBins = 100;
  fEmergEneBinDE = (fEmergEneMax-fEmergEneMin)/fEmergEneNoBins;
  fEmergEneUpElectEvt = 0;
  fEmergEneUpGammaEvt = 0;
  fEmergEneUpPositEvt = 0;
  fEmergEneDownElectEvt = 0;
  fEmergEneDownGammaEvt = 0;
  fEmergEneDownPositEvt = 0;

  fMessenger = new PenOutputMessenger(this);
  G4cout << "PenOutputPrinter object created" << G4endl;
}


PenOutputPrinter::~PenOutputPrinter()
{
  delete fMessenger;

 if (fGridRZDoseEvt) {
    fGridRZDoseEvt->clear();
    fDepthDoseEvt->clear();
    delete fGridRZDoseEvt;
    delete fDepthDoseEvt;
  }

  if (fSpcEnddetFlag) {
    for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
      fSpcEnddetFileVec->pop_back();
    }
    fSpcEnddetFileVec->clear();
    fSpcEnddetNameVec->clear();
    fSpcEnddetMinVec->clear();
    fSpcEnddetMaxVec->clear();
    fSpcEnddetNoBinsVec->clear();
    fSpcEnddetEvtVec->clear();
    delete fSpcEnddetFileVec;
    delete fSpcEnddetNameVec;
    delete fSpcEnddetMinVec;
    delete fSpcEnddetMaxVec;
    delete fSpcEnddetNoBinsVec;
    delete fSpcEnddetEvtVec;
  }

  if (fEmergingFlag) {
    fPolarAngleElectEvt->clear();
    fPolarAngleGammaEvt->clear();
    fPolarAnglePositEvt->clear();
    delete fPolarAngleElectEvt;
    delete fPolarAngleGammaEvt;
    delete fPolarAnglePositEvt;
    fEmergEneUpElectEvt->clear();
    fEmergEneUpGammaEvt->clear();
    fEmergEneUpPositEvt->clear();
    fEmergEneDownElectEvt->clear();
    fEmergEneDownGammaEvt->clear();
    fEmergEneDownPositEvt->clear();
    delete fEmergEneUpElectEvt;
    delete fEmergEneUpGammaEvt;
    delete fEmergEneUpPositEvt;
    delete fEmergEneDownElectEvt;
    delete fEmergEneDownGammaEvt;
    delete fEmergEneDownPositEvt;
  }
  G4cout << "PenOutputPrinter object destroyed" << G4endl;
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CreateFiles()
{
  if (fRZDoseMapFlag) CreateRZDoseMapFile();
  if (fSpcEnddetFlag) CreateSpcEnddetFile();
  if (fEmergingFlag) {
    CreatePolarAngleFile();
    CreateEmergEneFiles();
  }
  if (fTrkFlag) CreateTrkFiles();
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CreateRZDoseMapFile()
{
  // Check fGridZmin & fGridZmax are correct
  if (fGridZmin > fGridZmax) {
    G4Exception("PenOutputPrinter::CreateRZDoseMapFile()",
		"OutputPrinter0001", FatalException,
		"fGridZmin > fGridZmax");
  }

  fRZDoseMapFile.open("2d-dose-map.dat");
  if (fRZDoseMapFile.is_open()) {
    G4cout << "Opening \"2d-dose-map.dat\" file..." << G4endl;
    fRZDoseMapFile << std::scientific;
  }
  fDoseZFile.open("z-dose.dat");
  if (fDoseZFile.is_open()) {
    G4cout << "Opening \"z-dose.dat\" file..." << G4endl;
    fDoseZFile << std::scientific;
  }
  fDepthDoseFile.open("depth-dose.dat");
  if (fDepthDoseFile.is_open()) {
    G4cout << "Opening \"depth-dose.dat\" file..." << G4endl;
    fDepthDoseFile << std::scientific;
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CreateSpcEnddetFile()
{

  for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
    // Check (*fSpcEnddetMinVec)[ii] & (*fSpcEnddetMaxVec)[ii] are correct
    if ((*fSpcEnddetMinVec)[ii] > (*fSpcEnddetMaxVec)[ii]) {
      G4cout << "fSpcEnddetMin > fSpcEnddetMax at vector element " << ii << "."
	     << G4endl;
      G4Exception("PenOutputPrinter::CreateSpcEnddetFile()",
		  "OutputPrinter0002", FatalException,
		  "fSpcEnddetMin > fSpcEnddetMax");
    }

    std::stringstream ss;
    G4String enddetOrder;
    ss << std::setfill('0') << std::setw(2);
    ss << (ii+1);
    ss >> enddetOrder;
    G4String spcEnddetFileName;
    spcEnddetFileName = "spc-enddet-" + enddetOrder + ".dat";

    std::ofstream* spcEnddetFile = new std::ofstream;
    spcEnddetFile->open(spcEnddetFileName.c_str());
    if (spcEnddetFile->is_open()) {
      G4cout << "Opening \"" << spcEnddetFileName <<"\" file..." << G4endl;
      (*spcEnddetFile) << std::scientific;
    }

    fSpcEnddetFileVec->push_back(spcEnddetFile);
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CreatePolarAngleFile()
{
  fPolarAngleFile.open("polar-angle.dat");
  if (fPolarAngleFile.is_open()) {
    G4cout << "Opening \"polar-angle.dat\" file..." << G4endl;
    fPolarAngleFile << std::scientific;
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CreateEmergEneFiles()
{
  fEmergEneUpFile.open("energy-up.dat");
  if (fEmergEneUpFile.is_open()) {
    G4cout << "Opening \"energy-up.dat\" file..." << G4endl;
    fEmergEneUpFile << std::scientific;
  }
  fEmergEneDownFile.open("energy-down.dat");
  if (fEmergEneDownFile.is_open()) {
    G4cout << "Opening \"energy-down.dat\" file..." << G4endl;
    fEmergEneDownFile << std::scientific;
  }
}




//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CreateTrkFiles()
{
  fElectronTrkFile.open("electron.trk");
  fElectronTrkFile << std::scientific;
  if (fElectronTrkFile.is_open())
    G4cout << "Opening \"electron.trk\" file..." << G4endl;

  fGammaTrkFile.open("gamma.trk");
  fGammaTrkFile << std::scientific;
  if (fGammaTrkFile.is_open())
    G4cout << "Opening \"gamma.trk\" file..." << G4endl;

  fPositronTrkFile.open("positron.trk");
  fPositronTrkFile << std::scientific;
  if (fPositronTrkFile.is_open())
    G4cout << "Opening \"positron.trk\" file..." << G4endl;
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::BookForMaster()
{
  G4cout << "PenOutputPrinter::BookForMaster() called" << G4endl;

  if (fRZDoseMapFlag) {
    // Update fGridBinDZ and fGridBinDR with values passed through UI command
    fGridBinDZ = (fGridZmax-fGridZmin)/fGridZnoBins; 
    fGridBinDR = (fGridRmax)/fGridRnoBins;

    fGridRZDoseEvt = new std::vector<G4double>;
    fGridRZDoseEvt->reserve( fGridRnoBins*fGridZnoBins );
    for (G4int ii = 0; ii < fGridRnoBins; ii++) {
      for (G4int jj = 0; jj < fGridZnoBins; jj++) {
	fGridRZDoseEvt->push_back(0.0);
      }
    }
    fDepthDoseEvt = new std::vector<G4double>;
    fDepthDoseEvt->reserve( fGridZnoBins );
    for (G4int ii = 0; ii < fGridZnoBins; ii++) {
      fDepthDoseEvt->push_back(0.0);
    }
  }

  if (fSpcEnddetFlag) {
    // Push back a zero for each spc-enddep detector
    for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
      fSpcEnddetEvtVec->push_back(0.0);
    }
  }

  if (fEmergingFlag) {
    // polar-angle
    fPolarAngleElectEvt = new std::vector<G4double>;
    fPolarAngleGammaEvt = new std::vector<G4double>;
    fPolarAnglePositEvt = new std::vector<G4double>;
    fPolarAngleElectEvt->reserve(fPolarAngleNoBins);
    fPolarAngleGammaEvt->reserve(fPolarAngleNoBins);
    fPolarAnglePositEvt->reserve(fPolarAngleNoBins);
    for (G4int ii = 0; ii < fPolarAngleNoBins; ii++) {
      fPolarAngleElectEvt->push_back(0.0);
      fPolarAngleGammaEvt->push_back(0.0);
      fPolarAnglePositEvt->push_back(0.0);
    }
    // Update fPolarAngleBinDT with values passed through UI command
    fPolarAngleBinDT = (180.*deg)/fPolarAngleNoBins;

    // energy-up & energy-down
    fEmergEneUpElectEvt = new std::vector<G4double>;
    fEmergEneUpGammaEvt = new std::vector<G4double>;
    fEmergEneUpPositEvt = new std::vector<G4double>;
    fEmergEneUpElectEvt->reserve(fEmergEneNoBins);
    fEmergEneUpGammaEvt->reserve(fEmergEneNoBins);
    fEmergEneUpPositEvt->reserve(fEmergEneNoBins);
    fEmergEneDownElectEvt = new std::vector<G4double>;
    fEmergEneDownGammaEvt = new std::vector<G4double>;
    fEmergEneDownPositEvt = new std::vector<G4double>;
    fEmergEneDownElectEvt->reserve(fEmergEneNoBins);
    fEmergEneDownGammaEvt->reserve(fEmergEneNoBins);
    fEmergEneDownPositEvt->reserve(fEmergEneNoBins);
    for (G4int ii = 0; ii < fEmergEneNoBins; ii++) {
      fEmergEneUpElectEvt->push_back(0.0);
      fEmergEneUpGammaEvt->push_back(0.0);
      fEmergEneUpPositEvt->push_back(0.0);
      fEmergEneDownElectEvt->push_back(0.0);
      fEmergEneDownGammaEvt->push_back(0.0);
      fEmergEneDownPositEvt->push_back(0.0);
    }
    // Update fEmergEneBinDE using values passed through UI command
    fEmergEneBinDE = (fEmergEneMax-fEmergEneMin)/fEmergEneNoBins;
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::BookForWorker()
{
  G4cout << "PenOutputPrinter::BookForWorker() called" << G4endl;

  if (fRZDoseMapFlag) {
    // Update fGridBinDZ and fGridBinDR with values passed through UI command
    fGridBinDZ = (fGridZmax-fGridZmin)/fGridZnoBins;
    fGridBinDR = (fGridRmax)/fGridRnoBins;

    fGridRZDoseEvt = new std::vector<G4double>;
    fGridRZDoseEvt->reserve( fGridRnoBins*fGridZnoBins );
    for (G4int ii = 0; ii < fGridRnoBins; ii++) {
      for (G4int jj = 0; jj < fGridZnoBins; jj++) {
	fGridRZDoseEvt->push_back(0.0);
      }
    }
    fDepthDoseEvt = new std::vector<G4double>;
    fDepthDoseEvt->reserve( fGridZnoBins );
    for (G4int ii = 0; ii < fGridZnoBins; ii++) {
      fDepthDoseEvt->push_back(0.0);
    }
  }

  if (fSpcEnddetFlag) {
    // Push back a zero for each spc-enddep detector
    for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
      fSpcEnddetEvtVec->push_back(0.0);
    }
  }

  if (fEmergingFlag) {
    // polar-angle
    fPolarAngleElectEvt = new std::vector<G4double>;
    fPolarAngleGammaEvt = new std::vector<G4double>;
    fPolarAnglePositEvt = new std::vector<G4double>;
    fPolarAngleElectEvt->reserve(fPolarAngleNoBins);
    fPolarAngleGammaEvt->reserve(fPolarAngleNoBins);
    fPolarAnglePositEvt->reserve(fPolarAngleNoBins);
    for (G4int ii = 0; ii < fPolarAngleNoBins; ii++) {
      fPolarAngleElectEvt->push_back(0.0);
      fPolarAngleGammaEvt->push_back(0.0);
      fPolarAnglePositEvt->push_back(0.0);
    }
    // Update fPolarAngleBinDT with values passed through UI command
    fPolarAngleBinDT = (180.*deg)/fPolarAngleNoBins;

    // energy-up & energy-down
    fEmergEneUpElectEvt = new std::vector<G4double>;
    fEmergEneUpGammaEvt = new std::vector<G4double>;
    fEmergEneUpPositEvt = new std::vector<G4double>;
    fEmergEneUpElectEvt->reserve(fEmergEneNoBins);
    fEmergEneUpGammaEvt->reserve(fEmergEneNoBins);
    fEmergEneUpPositEvt->reserve(fEmergEneNoBins);
    fEmergEneDownElectEvt = new std::vector<G4double>;
    fEmergEneDownGammaEvt = new std::vector<G4double>;
    fEmergEneDownPositEvt = new std::vector<G4double>;
    fEmergEneDownElectEvt->reserve(fEmergEneNoBins);
    fEmergEneDownGammaEvt->reserve(fEmergEneNoBins);
    fEmergEneDownPositEvt->reserve(fEmergEneNoBins);
    for (G4int ii = 0; ii < fEmergEneNoBins; ii++) {
      fEmergEneUpElectEvt->push_back(0.0);
      fEmergEneUpGammaEvt->push_back(0.0);
      fEmergEneUpPositEvt->push_back(0.0);
      fEmergEneDownElectEvt->push_back(0.0);
      fEmergEneDownGammaEvt->push_back(0.0);
      fEmergEneDownPositEvt->push_back(0.0);
    }
    // Update fEmergEneBinDE using values passed through UI command
    fEmergEneBinDE = (fEmergEneMax-fEmergEneMin)/fEmergEneNoBins;
  }
}




//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::PrintPreTrackInfo(const G4Track* aTrack)
{
  if (fTrkFlag) {

    G4ParticleDefinition* pDef = aTrack->GetDefinition();

    // Use initial track information to introduce
    // blank lines into the track output files of each PENELOPE particle
    // (Note: This information is ONLY accesible at PreUserTrackingAction
    // ------------------------------------------------------------------

    // pe-
    if (pDef == PenElectron::Definition()) {
      if (aTrack->GetParentID() == 0) {
	// single blank if a primary track starts (maybe double in future)
	fElectronTrkFile /*<< std::endl*/ << std::endl;
      } else {
	// single blank if a secondary track starts
	fElectronTrkFile << std::endl;
      }
      // G4ThreeVector pos = aTrack->GetPosition();
      // G4double ene = aTrack->GetKineticEnergy();
      // G4cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << ene
      // 	     << G4endl;
    }
    // pgamma
    else if (pDef == PenGamma::Definition()) {
      if (aTrack->GetParentID() == 0) {
	// single blank if a primary track starts (maybe double in future)
	fGammaTrkFile /*<< std::endl*/ << std::endl;
      } else {
	// single blank if a secondary track starts
	fGammaTrkFile << std::endl;
      }
      // G4ThreeVector pos = aTrack->GetPosition();
      // G4double ene = aTrack->GetKineticEnergy();
      // G4cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << ene
      // 	     << G4endl;
    }
    // pe+
    else if (pDef == PenPositron::Definition()) {
      if (aTrack->GetParentID() == 0) {
	// single blank if a primary track starts (maybe double in future)
	fPositronTrkFile /*<< std::endl*/ << std::endl;
      } else {
	// single blank if a secondary track starts
	fPositronTrkFile << std::endl;
      }
      // G4ThreeVector pos = aTrack->GetPosition();
      // G4double ene = aTrack->GetKineticEnergy();
      // G4cout << pos.x() << " " << pos.y() << " " << pos.z() << " " << ene
      // 	     << G4endl;
    }
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::PrintStepInfo(const G4Step* aStep)
{
  if (fTrkFlag) {

    // Dump the step position for each PENELOPE particle
    // -------------------------------------------------

    const G4Track* aTrack = aStep->GetTrack();
    const G4ParticleDefinition* pDef = aTrack->GetDefinition();

    // pe-
    if (pDef == PenElectron::Definition()) {
      G4ThreeVector pos = aTrack->GetPosition();
      G4double ene = aTrack->GetKineticEnergy();
      fElectronTrkFile << std::setw(13) << std::setprecision(5)
		       << pos.x()/cm << " "
		       << std::setw(13) << std::setprecision(5)
		       << pos.y()/cm << " "
		       << std::setw(13) << std::setprecision(5)
		       << pos.z()/cm << " "
		       << std::setw(13) << std::setprecision(5)
		       << ene/eV
		       << std::endl;
    }
    // pgamma
    if (pDef == PenGamma::Definition()) {
      G4ThreeVector pos = aTrack->GetPosition();
      G4double ene = aTrack->GetKineticEnergy();
      fGammaTrkFile << std::setw(13) << std::setprecision(5)
		    << pos.x()/cm << " "
		    << std::setw(13) << std::setprecision(5)
		    << pos.y()/cm << " "
		    << std::setw(13) << std::setprecision(5)
		    << pos.z()/cm << " "
		    << std::setw(13) << std::setprecision(5)
		    << ene/eV
		    << std::endl;
    }
    // pe-
    if (pDef == PenPositron::Definition()) {
      G4ThreeVector pos = aTrack->GetPosition();
      G4double ene = aTrack->GetKineticEnergy();
      fPositronTrkFile << std::setw(13) << std::setprecision(5)
		       << pos.x()/cm << " "
		       << std::setw(13) << std::setprecision(5)
		       << pos.y()/cm << " "
		       << std::setw(13) << std::setprecision(5)
		       << pos.z()/cm << " "
		       << std::setw(13) << std::setprecision(5)
		       << ene/eV
		       << std::endl;
    }
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::StoreStepInfo(const G4Step* aStep)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4StepPoint* preP = aStep->GetPreStepPoint();
  G4StepPoint* postP = aStep->GetPostStepPoint();

  // edep is NOT tallied if postP is at world volume
  // This is needed because, with very very tiny probability,
  // a particle can interact at world volume with.
  if ( edep > 0.0 &&
       postP->GetTouchableHandle()->GetHistoryDepth() != 0 ) {

    G4ThreeVector point = postP->GetPosition();
    // G4cout << "postPos = " << point << G4endl;

    //MACG Soft Edep is deposited at hinges (method NOT USED)
    // if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0) {
    //   G4ThreeVector prePos = aStep->GetPreStepPoint()->GetPosition();
    //   point = prePos + G4UniformRand()*(point-prePos);
    //   // G4cout << "prePos = " << prePos << G4endl;
    //   // G4cout << "point = " << point << G4endl;
    // }

    // depth-dose & 2d-dose-map files
    if (fRZDoseMapFlag) {
      G4double zPos = point.z();
      G4double rPos = point.rho();
      if (zPos < fGridZmax && zPos > fGridZmin) {
	// Edep is only deposited at Post-step point, as
	// all processes are ''discrete'' in current version.
	G4double dens = postP->GetMaterial()->GetDensity();
        // if (dens < 0.001*g/cm3)  {
	//   G4cout << "Density[g/cm3]= " << dens/(g/cm3)
	// 	 << "  EvtID= "
	// 	 << G4EventManager::GetEventManager()->GetConstCurrentEvent()
	//     ->GetEventID()
	// 	 << "  edep[eV]= " << edep/eV << G4endl;
	//   G4cout << "\tPOST: zPos[cm]= " << zPos/cm << "  rPos[cm]= "
	// 	 << rPos/cm << G4endl;
	//   G4cout << "\tPRE: zPos[cm]= " << (preP->GetPosition().z())/cm
	// 	 << "  rPos[cm]= " << (preP->GetPosition().rho())/cm
	// 	 << G4endl;
	// }
	G4int binZ = static_cast<G4int>( (zPos-fGridZmin)/fGridBinDZ);
	(*fDepthDoseEvt)[binZ] += (edep/dens);
	if (rPos < fGridRmax) {
	  G4int binR = static_cast<G4int>( rPos/fGridBinDR );
	  (*fGridRZDoseEvt)[binR*fGridZnoBins+binZ] += (edep/dens);
	  // DEBUG!!
	  // if (rPos < 2.*fGridRmax/fGridRnoBins) {
	  //   G4cout << "\trPos[cm]=" << rPos/cm;
	  //   G4cout << "\tzPos[cm]=" << zPos/cm;
	  //   G4cout << "\tedep[eV]= " << edep/eV;
	  //   G4cout << "\tat bin= " << binR*fGridZnoBins+binZ << G4endl;
	  // }
	}
      }
    }

    // spc-enddet-## files
    if (fSpcEnddetFlag) {
      G4String volName = aStep->GetTrack()->GetVolume()->GetName();
      // G4cout << "edep[eV] at " << volName << ": " << edep/eV << G4endl;
      for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
	if ( volName == (*fSpcEnddetNameVec)[ii] ) {
	  (*fSpcEnddetEvtVec)[ii] += edep;
	  // G4cout << "\tfSpcEnddetEvtVec[" << ii << "]= " << (*fSpcEnddetEvtVec)[ii]
	  // 	 << G4endl;
	}
      }
    }
  }

  // emerging particles
  // ------------------
  
  if (fEmergingFlag) {

    // Store the theta angle and energy of the particles emerging
    // ouf of the enclosure (a.k.a. world volume).
    // NOTE: To avoid the fact of not defining a vacuum material in Interface
    // particle's theta is recorded once it escapes from geometry to world

    G4int preDepth = preP->GetTouchableHandle()->GetHistoryDepth();
    G4int postDepth = postP->GetTouchableHandle()->GetHistoryDepth();
    // G4cout << "\tpreDepth = " << preDepth;
    // G4cout << "\tpostDepth = " << postDepth << G4endl;

    if (postDepth == 0 && preDepth > 0) {
      // Particle coming back to the enclosure.
      // Get type and momentum theta and kinetic energy.
      
      G4ParticleDefinition* pDef = aStep->GetTrack()->GetDefinition();
      // G4cout << "ParticleName = " << pDef->GetParticleName() << G4endl;

      G4double thetaP = postP->GetMomentumDirection().theta();
      // G4cout << "\tthetaP[deg] = " << thetaP/deg << G4endl;
 
      // polar-angle.dat
      // ---------------
      G4int binT = static_cast<G4int>((thetaP)/fPolarAngleBinDT);
      // G4cout << "\tbinT = " << binT << G4endl;

      if (pDef == PenElectron::Definition() ||
	  pDef == G4Electron::Definition() )
	(*fPolarAngleElectEvt)[binT] += 1.0;
      else if (pDef == PenGamma::Definition() ||
	       pDef == G4Gamma::Definition() )
	(*fPolarAngleGammaEvt)[binT] += 1.0;
      else if (pDef == PenPositron::Definition() ||
	       pDef == G4Positron::Definition() )
	(*fPolarAnglePositEvt)[binT] += 1.0;

      // energy-up.dat & energy-down.dat
      // -------------------------------
      G4double kinE = postP->GetKineticEnergy();
      // G4cout << "\tkinE[eV] = " << kinE/eV << G4endl;

      if (!(kinE < fEmergEneMin || kinE > fEmergEneMax)) {
	G4int binE = static_cast<G4int>((kinE-fEmergEneMin)/fEmergEneBinDE);
	// G4cout << "\tbinE = " << binE << G4endl;

	if (thetaP < M_PI/2.) { // energy-up.dat
	  if (pDef == PenElectron::Definition() ||
	      pDef == G4Electron::Definition() )
	    (*fEmergEneUpElectEvt)[binE] += 1.0;
	  else if (pDef == PenGamma::Definition() ||
		   pDef == G4Gamma::Definition() )
	    (*fEmergEneUpGammaEvt)[binE] += 1.0;
	  else if (pDef == PenPositron::Definition() ||
		   pDef == G4Positron::Definition() )
	    (*fEmergEneUpPositEvt)[binE] += 1.0;
	}
	else { // energy-down.dat
	  if (pDef == PenElectron::Definition() ||
	      pDef == G4Electron::Definition() )
	    (*fEmergEneDownElectEvt)[binE] += 1.0;
	  else if (pDef == PenGamma::Definition() ||
		   pDef == G4Gamma::Definition() )
	    (*fEmergEneDownGammaEvt)[binE] += 1.0;
	  else if (pDef == PenPositron::Definition() ||
		   pDef == G4Positron::Definition() )
	    (*fEmergEneDownPositEvt)[binE] += 1.0;
	}
      }
    }
  }

}




//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::ResetEventInfo()
{
  // This is called from BeginOfEventAction
  // ---------------------------------------
  // Reset variables persistent during one event
  if (fRZDoseMapFlag) {
    G4int rzBins = this->GetGridRZnoBins();
    for (G4int ii = 0; ii < rzBins; ii++) {
      (*fGridRZDoseEvt)[ii] = 0.0;
    }
    G4int zBins = this->GetGridZnoBins();
    for (G4int ii = 0; ii < zBins; ii++) {
      (*fDepthDoseEvt)[ii] = 0.0;
    }    
  }
  if (fSpcEnddetFlag) {
    for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
      (*fSpcEnddetEvtVec)[ii] = 0.0;
    }
  }
  if (fEmergingFlag) {
    for (G4int ii = 0; ii < fPolarAngleNoBins; ii++) {
      (*fPolarAngleElectEvt)[ii] = 0.;
      (*fPolarAngleGammaEvt)[ii] = 0.;
      (*fPolarAnglePositEvt)[ii] = 0.;
    }
    for (G4int ii = 0; ii < fEmergEneNoBins; ii++) {
      (*fEmergEneUpElectEvt)[ii] = 0.;
      (*fEmergEneUpGammaEvt)[ii] = 0.;
      (*fEmergEneUpPositEvt)[ii] = 0.;
      (*fEmergEneDownElectEvt)[ii] = 0.;
      (*fEmergEneDownGammaEvt)[ii] = 0.;
      (*fEmergEneDownPositEvt)[ii] = 0.;
    }
  }
}




//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::DumpAll(const Run* masterRun)
{
  G4cout << "PenOutputPrinter::DumpAll() called" << G4endl;
  // 2d-dose-map & z-dose file
  if (fRZDoseMapFlag) {
    DumpRZDoseMap(masterRun);
    DumpDepthDoseMap(masterRun);
  }
  // spc-enddet-## files
  if (fSpcEnddetFlag) DumpSpcEnddet(masterRun);
  // emerging particle files
  if (fEmergingFlag) {
    // polar-angle file
    DumpPolarAngle(masterRun);
    // energy-up & energy-down files
    DumpEmergEne(masterRun);
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::DumpRZDoseMap(const Run* masterRun)
{
  G4double nEvt = static_cast<G4double>
    (G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent());
  // G4cout << "nEvt = " << nEvt << G4endl;

  // // Creating a navigator to locate points in geometry and get densities
  // G4Navigator aNavigator;
  // aNavigator.SetWorldVolume(fDetConst->GetPhysWorld() );

  // =======================
  // 2d-dose-map.dat HEADER
  // =======================
  fRZDoseMapFile
    << " #  Results from PEN-GEANT4. Dose distribution.\n"
    << " #  Dose-map cylinder:         RU =  " << fGridRmax/cm << " cm\n"
    << " #     ZL =  " << fGridZmin/cm << " cm,  ZU =  " << fGridZmax/cm 
    << " cm\n"
    << " #  Numbers of bins:     NBR = " << fGridRnoBins << ", NBZ = " 
    << fGridZnoBins << std::endl;
  fRZDoseMapFile
    << " #\n"
    << " #  columns 1 and 2: coordinates R,Z of the bin centres\n"
    << " #  3rd column: dose (eV/g).\n"
    << " #  4th column: statistical uncertainty (3 sigma)." << std::endl;

  // ===================
  // z-dose.dat HEADER
  // ===================
  fDoseZFile
    << " #  Results from PEN-GEANT4.\n"
    << " #  Dose distribution along the central Z axis.\n"
    << " #  1st column: z (cm).\n"
    << " #  2nd column: dose (eV/g).\n"
    << " #  3rd column: statistical uncertainty (3 sigma)." << std::endl;

  // retrieve Dose RZ maps
  std::vector<G4double>* gridRZDoseRun = masterRun->GetGridRZDoseRun();
  std::vector<G4double>* gridRZDoseRun2 = masterRun->GetGridRZDoseRun2();

  for (G4int ii = 0; ii < fGridRnoBins; ii++) {
    G4double binR = (0.5+ii)*fGridBinDR;
    G4double binVol = M_PI*(2*ii+1)*fGridBinDR*fGridBinDR*fGridBinDZ;

    for (G4int jj = 0; jj < fGridZnoBins; jj++) {
      G4double binZ = (0.5+jj)*fGridBinDZ + fGridZmin;
      // G4ThreeVector point(binR, 0., binZ); // Using XZ plane as reference
      // G4double binDensity = aNavigator.LocateGlobalPointAndSetup(point)
      // 	->GetLogicalVolume()->GetMaterial()->GetDensity();
      // G4double binMass = binDensity*binVol;
      // // G4cout << "binR[cm]= " << binR/cm << "\tbinZ[cm]= " << binZ/cm
      // // 	     << "\tdens[g/cm3]= " << binDensity/(g/cm3)
      // // 	     << "\tvol[cm3]= " << binVol/cm3
      // // 	     << "\tmass[g]= " << binMass/g << G4endl;
      G4int bin = ii*fGridZnoBins+jj;
      G4double dose = ((*gridRZDoseRun)[bin]/binVol)/nEvt;
      G4double varDose = 0.0;
      if (dose > 0) {
	varDose = ((*gridRZDoseRun2)[bin]/binVol/binVol)/nEvt - dose*dose;
	if (varDose/dose/dose < 1.e-14) {
	  // This may happen due to rounding errors (A-B=epsilon)
	  varDose = 0.0;
	}
      }
      G4double doseErr = std::sqrt(varDose/nEvt); // approx (nEvt-1)
      doseErr *= 3.0;

      fRZDoseMapFile << std::setw(14) << std::setprecision(6)
		     << binR/cm
		     << std::setw(14) << std::setprecision(6)
		     << binZ/cm
		     << std::setw(14) << std::setprecision(6)
		     << dose/(eV/g)
		     << std::setw(10) << std::setprecision(2)
		     << doseErr/(eV/g)
		     << std::endl;
      if ( ii == 0 ) {  // Dump dose on axis to z-dose.dat
	fDoseZFile << std::setw(14) << std::setprecision(6)
		   << binZ/cm
		   << std::setw(14) << std::setprecision(6)
		   << dose/(eV/g)
		   << std::setw(10) << std::setprecision(2)
		   << doseErr/(eV/g)
		   << std::endl;
      }
    }
    fRZDoseMapFile << std::endl;  // A single-blank for GNUplot Grid format
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::DumpDepthDoseMap(const Run* masterRun)
{
  G4double nEvt = static_cast<G4double>
    (G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent());
  // G4cout << "nEvt = " << nEvt << G4endl;

  // // Creating a navigator to locate points in geometry and get densities
  // G4Navigator aNavigator;
  // aNavigator.SetWorldVolume(fDetConst->GetPhysWorld() );

  // ======================
  // depth-dose.dat HEADER
  // ======================
  fDepthDoseFile
    << " #  Results from PEN-GEANT4. Depth-dose distribution.\n"
    << " #  (integrated over X and Y within the volume of the material system).\n"
    << " #  1st column: z coordinate (cm).\n"
    << " #  2nd column: depth-dose (eV/(g/cm**2)).\n"
    << " #  3rd column: statistical uncertainty (3 sigma).\n"
    << " #\n"
    << " #  NOTE: The calculated dose distribution is correct only when the\n"
    << " #        Z bins have uniform mass density.\n"
    << " #"
    << std::endl;

  // retrieve depth-dose maps
  std::vector<G4double>* depthDoseRun = masterRun->GetDepthDoseRun();
  std::vector<G4double>* depthDoseRun2 = masterRun->GetDepthDoseRun2();

  for (G4int jj = 0; jj < fGridZnoBins; jj++) {
    G4double binZ = (0.5+jj)*fGridBinDZ + fGridZmin;
    // G4ThreeVector point(0., 0., binZ); // Using z-axis as reference
    // G4double binDensity = aNavigator.LocateGlobalPointAndSetup(point)
    //   ->GetLogicalVolume()->GetMaterial()->GetDensity();
    // // G4cout << "binZ[cm]= " << binZ/cm
    // // 	   << "\tdens[g/cm3]= " << binDensity/(g/cm3) << G4endl;
    G4double dose = ((*depthDoseRun)[jj]/fGridBinDZ)/nEvt;
    G4double varDose = 0.0;
    if (dose > 0) {
      varDose = ((*depthDoseRun2)[jj]/fGridBinDZ/fGridBinDZ)/nEvt - dose*dose;
      if (varDose/dose/dose < 1.e-14) {
	// This may happen due to rounding errors (A-B=epsilon)
	varDose = 0.0;
      }
    }
    G4double doseErr = std::sqrt(varDose/nEvt); // approx (nEvt-1)
    doseErr *= 3.0;
    
    fDepthDoseFile << std::setw(14) << std::setprecision(6)
		   << binZ/cm
		   << std::setw(14) << std::setprecision(6)
		   << dose/(eV/(g/cm2))
		   << std::setw(10) << std::setprecision(2)
		   << doseErr/(eV/(g/cm2))
		   << std::endl;
  }

}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::DumpSpcEnddet(const Run* masterRun)
{
  // Get number of events simulated  
  G4double nEvt = static_cast<G4double>
    (G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent());

  for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
    // Write file header
    (*(*fSpcEnddetFileVec)[ii])
      << " #  Results from PEN-GEANT4. Output from energy-deposition detector # "
      << ii+1 << "\n"
      << " #  Deposited energy spectrum.\n"
      << " #  WARNING: May be strongly biased if interaction forcing is used!\n"
      << " #  1st column: deposited energy (eV).\n"
      << " #  2nd column: probability density (1/(eV*particle)).\n"
      << " #  3rd column: statistical uncertainty (3 sigma)." << std::endl;

    G4double binDE = ((*fSpcEnddetMaxVec)[ii]-(*fSpcEnddetMinVec)[ii]);
    binDE /= (*fSpcEnddetNoBinsVec)[ii];

    // Retrieve results from Run and write to file
    std::vector<G4double>* spcEnddetRun = masterRun->GetSpcEnddetRun(ii);
    std::vector<G4double>* spcEnddetRun2 = masterRun->GetSpcEnddetRun2(ii);

    for (G4int jj = 0; jj < (*fSpcEnddetNoBinsVec)[ii]; jj++) {
      G4double bin = (0.5+jj)*binDE + (*fSpcEnddetMinVec)[ii];

      G4double w = ((*spcEnddetRun)[jj])/binDE/nEvt; // per unitE and nEvt
      G4double wErr = std::sqrt( (*spcEnddetRun2)[jj] )/binDE/nEvt;
      wErr *= 3.0;
    
      (*(*fSpcEnddetFileVec)[ii]) << std::setw(14) << std::setprecision(6)
				  << bin/eV
				  << std::setw(14) << std::setprecision(6)
				  << w/(1/eV)
				  << std::setw(10) << std::setprecision(2)
				  << wErr/(1/eV)
				  << std::endl;
    }
  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::DumpPolarAngle(const Run* masterRun)
{
  // Write file header
  fPolarAngleFile
    << " #  Results from PEN-GEANT4.\n"
    << " #  Angular distributions of emerging particles.\n"
    << " #  1st column: THETA (deg).\n"
    << " #  2nd and 3rd columns: PDF and STU for electrons.\n"
    << " #  4th and 5th columns: PDF and STU for photons.\n"
    << " #  6th and 7th columns: PDF and STU for positrons\n"
    << " #  PDFs and STUs in units of 1/(sr*primary_particle)." << std::endl;

  // Get number of events simulated  
  G4double nEvt = static_cast<G4double>
    (G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent());

  // Retrieve results from Run and write to file
  
  std::vector<G4double>* polarAngleE = masterRun->GetPolarAngleElectRun();
  std::vector<G4double>* polarAngleE2 = masterRun->GetPolarAngleElectRun2();
  std::vector<G4double>* polarAngleG = masterRun->GetPolarAngleGammaRun();
  std::vector<G4double>* polarAngleG2 = masterRun->GetPolarAngleGammaRun2();
  std::vector<G4double>* polarAngleP = masterRun->GetPolarAnglePositRun();
  std::vector<G4double>* polarAngleP2 = masterRun->GetPolarAnglePositRun2();

  for (G4int ii = 0; ii < fPolarAngleNoBins; ii++) {
    G4double theta = (0.5+ii)*fPolarAngleBinDT;
    G4double sAng = 2*M_PI*(std::cos(theta-0.5*fPolarAngleBinDT)
			    - std::cos(theta+0.5*fPolarAngleBinDT));
				     
    G4double wE = ((*polarAngleE)[ii])/sAng/ nEvt;
    G4double wE2 = std::sqrt( (*polarAngleE2)[ii] )/sAng/nEvt;
    wE2 *= 3.0;
    G4double wG = ((*polarAngleG)[ii])/sAng/ nEvt;
    G4double wG2 = std::sqrt( (*polarAngleG2)[ii] )/sAng/nEvt;
    wG2 *= 3.0;
    G4double wP = ((*polarAngleP)[ii])/sAng/ nEvt;
    G4double wP2 = std::sqrt( (*polarAngleP2)[ii] )/sAng/nEvt;
    wP2 *= 3.0;

    fPolarAngleFile << std::setw(14) << std::setprecision(6)
		    << theta/deg
		    << std::setw(14) << std::setprecision(6)
		    << wE/(1/sr)
		    << std::setw(10) << std::setprecision(2)
		    << wE2/(1/sr)
		    << std::setw(14) << std::setprecision(6)
		    << wG/(1/sr)
		    << std::setw(10) << std::setprecision(2)
		    << wG2/(1/sr)
		    << std::setw(14) << std::setprecision(6)
		    << wP/(1/sr)
		    << std::setw(10) << std::setprecision(2)
		    << wP2/(1/sr)
		    << std::endl;

  }
}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::DumpEmergEne(const Run* masterRun)
{
  // Write file headers
  fEmergEneUpFile
    << " #  Results from PEN-GEANT4.\n"
    << " #  Energy distributions of upbound particles.\n"
    << " #  1st column: E (eV).\n"
    << " #  2nd and 3rd columns: PDF and STU for electrons.\n"
    << " #  4th and 5th columns: PDF and STU for photons.\n"
    << " #  6th and 7th columns: PDF and STU for positrons.\n"
    << " #  PDFs and STUs in units of 1/(eV*primary_particle)." << std::endl;
  fEmergEneDownFile
    << " #  Results from PEN-GEANT4.\n"
    << " #  Energy distributions of downbound particles.\n"
    << " #  1st column: E (eV).\n"
    << " #  2nd and 3rd columns: PDF and STU for electrons.\n"
    << " #  4th and 5th columns: PDF and STU for photons.\n"
    << " #  6th and 7th columns: PDF and STU for positrons.\n"
    << " #  PDFs and STUs in units of 1/(eV*primary_particle)." << std::endl;

  // Get number of events simulated  
  G4double nEvt = static_cast<G4double>
    (G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent());

  // Retrieve energy-up results from Run and write to file
  // ---------------------------------------------------------

  std::vector<G4double>* emergEneUpE = masterRun->GetEmergEneUpElectRun();
  std::vector<G4double>* emergEneUpE2 = masterRun->GetEmergEneUpElectRun2();
  std::vector<G4double>* emergEneUpG = masterRun->GetEmergEneUpGammaRun();
  std::vector<G4double>* emergEneUpG2 = masterRun->GetEmergEneUpGammaRun2();
  std::vector<G4double>* emergEneUpP = masterRun->GetEmergEneUpPositRun();
  std::vector<G4double>* emergEneUpP2 = masterRun->GetEmergEneUpPositRun2();

  for (G4int ii = 0; ii < fEmergEneNoBins; ii++) {
    G4double ene = (0.5+ii)*fEmergEneBinDE + fEmergEneMin;

    G4double wE = ((*emergEneUpE)[ii])/fEmergEneBinDE/ nEvt;
    G4double wE2 = std::sqrt( (*emergEneUpE2)[ii] )/fEmergEneBinDE/nEvt;
    wE2 *= 3.0;
    G4double wG = ((*emergEneUpG)[ii])/fEmergEneBinDE/ nEvt;
    G4double wG2 = std::sqrt( (*emergEneUpG2)[ii] )/fEmergEneBinDE/nEvt;
    wG2 *= 3.0;
    G4double wP = ((*emergEneUpP)[ii])/fEmergEneBinDE/ nEvt;
    G4double wP2 = std::sqrt( (*emergEneUpP2)[ii] )/fEmergEneBinDE/nEvt;
    wP2 *= 3.0;

    fEmergEneUpFile << std::setw(14) << std::setprecision(6)
		    << ene/eV
		    << std::setw(14) << std::setprecision(6)
		    << wE/(1/eV)
		    << std::setw(10) << std::setprecision(2)
		    << wE2/(1/eV)
		    << std::setw(14) << std::setprecision(6)
		    << wG/(1/eV)
		    << std::setw(10) << std::setprecision(2)
		    << wG2/(1/eV)
		    << std::setw(14) << std::setprecision(6)
		    << wP/(1/eV)
		    << std::setw(10) << std::setprecision(2)
		    << wP2/(1/eV)
		    << std::endl;
  }

  // Retrieve energy-down results from Run and write to file
  // ---------------------------------------------------------

  std::vector<G4double>* emergEneDownE = masterRun->GetEmergEneDownElectRun();
  std::vector<G4double>* emergEneDownE2 = masterRun->GetEmergEneDownElectRun2();
  std::vector<G4double>* emergEneDownG = masterRun->GetEmergEneDownGammaRun();
  std::vector<G4double>* emergEneDownG2 = masterRun->GetEmergEneDownGammaRun2();
  std::vector<G4double>* emergEneDownP = masterRun->GetEmergEneDownPositRun();
  std::vector<G4double>* emergEneDownP2 = masterRun->GetEmergEneDownPositRun2();

  for (G4int ii = 0; ii < fEmergEneNoBins; ii++) {
    G4double ene = (0.5+ii)*fEmergEneBinDE + fEmergEneMin;

    G4double wE = ((*emergEneDownE)[ii])/fEmergEneBinDE/ nEvt;
    G4double wE2 = std::sqrt( (*emergEneDownE2)[ii] )/fEmergEneBinDE/nEvt;
    wE2 *= 3.0;
    G4double wG = ((*emergEneDownG)[ii])/fEmergEneBinDE/ nEvt;
    G4double wG2 = std::sqrt( (*emergEneDownG2)[ii] )/fEmergEneBinDE/nEvt;
    wG2 *= 3.0;
    G4double wP = ((*emergEneDownP)[ii])/fEmergEneBinDE/ nEvt;
    G4double wP2 = std::sqrt( (*emergEneDownP2)[ii] )/fEmergEneBinDE/nEvt;
    wP2 *= 3.0;

    fEmergEneDownFile << std::setw(14) << std::setprecision(6)
		      << ene/eV
		      << std::setw(14) << std::setprecision(6)
		      << wE/(1/eV)
		      << std::setw(10) << std::setprecision(2)
		      << wE2/(1/eV)
		      << std::setw(14) << std::setprecision(6)
		      << wG/(1/eV)
		      << std::setw(10) << std::setprecision(2)
		      << wG2/(1/eV)
		      << std::setw(14) << std::setprecision(6)
		      << wP/(1/eV)
		      << std::setw(10) << std::setprecision(2)
		      << wP2/(1/eV)
		      << std::endl;
  }

}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::CloseAllFiles()
{
  // 2d-dose-map, z-dose & depth-dose
  if (fRZDoseMapFlag) {
    fRZDoseMapFile.close();
    if ( !(fRZDoseMapFile.is_open()) )
      G4cout << "\"2d-dose-map.dat\" file closed successfully!" << G4endl;
    fDoseZFile.close();
    if ( !(fDoseZFile.is_open()) )
      G4cout << "\"z-dose.dat\" file closed successfully!" << G4endl;
    fDepthDoseFile.close();
    if ( !(fDepthDoseFile.is_open()) )
      G4cout << "\"depth-dose.dat\" file closed successfully!" << G4endl;
  }

  // spc-enddet-##
  if (fSpcEnddetFlag) {
    for (size_t ii = 0; ii < fNoSpcEnddet; ii++) {
      (*fSpcEnddetFileVec)[ii]->close();
      if ( !( (*fSpcEnddetFileVec)[ii]->is_open() ))
	G4cout << "\"spc-enddet-"
	       << std::setfill('0') << std::setw(2) << ii+1
	       << ".dat\" file closed successfully!" << G4endl;
    }
  }

  // emerging particles
  if (fEmergingFlag) {
    // polar-angle
    fPolarAngleFile.close();
    if ( !(fPolarAngleFile.is_open()) )
      G4cout << "\"polar-angle.dat\" file closed successfully!" << G4endl;
    // energy-up & energy-down
    fEmergEneUpFile.close();
    if ( !(fEmergEneUpFile.is_open()) )
      G4cout << "\"energy-up.dat\" file closed successfully!" << G4endl;
    fEmergEneDownFile.close();
    if ( !(fEmergEneDownFile.is_open()) )
      G4cout << "\"energy-down.dat\" file closed successfully!" << G4endl;
  }

  // tracking files
  if (fTrkFlag) {
    fElectronTrkFile.close();
    if ( !(fElectronTrkFile.is_open()) )
      G4cout << "\"electron.trk\" file closed successfully!" << G4endl;
    
    fGammaTrkFile.close();
    if ( !(fGammaTrkFile.is_open()) )
      G4cout << "\"gamma.trk\" file closed successfully!" << G4endl;

    fPositronTrkFile.close();
    if ( !(fPositronTrkFile.is_open()) )
      G4cout << "\"positron.trk\" file closed successfully!" << G4endl;
  }

}



//-------1---------2---------3---------4---------5---------6---------7---------8
void PenOutputPrinter::AddSpcEnddetPVName(G4String pvname)
{
  fSpcEnddetNameVec->push_back(pvname);

  // Set automatically fNoSpcEnddet and fSpcEnddetFlag if a PV is registered.
  fNoSpcEnddet = fSpcEnddetNameVec->size();
  if (!fSpcEnddetFlag) fSpcEnddetFlag = true;
}
