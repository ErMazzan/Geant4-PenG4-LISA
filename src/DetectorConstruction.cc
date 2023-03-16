
#include "DetectorConstruction.hh"

#include <array>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "DetectorMessenger.hh"

#include "PenInterface.hh"

DetectorConstruction::DetectorConstruction()
  :fConstructed(false)
{
  fPhysWorld = 0;
  fWorldMat = 0;
  fWorldSize=100.*m;
  
  fSCSize = 88.*mm;
  fSCThick = 4.*mm;
  
  fHouSize = 70.*mm;
  fHouThick = 8.*mm;

  fShellThick = 1 *um; //by now, both for Mo shielding and TM
  fTMSize = 44.*mm; //Total test mass area = fTMSize+2*fShellThick
  

  fMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if(!fConstructed)
  {
    fConstructed = true;
    DefineMaterials();
    SetupGeometry();
  }
  return fPhysWorld;
}

void DetectorConstruction::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
 
  // Define vacuum for world volume
  fWorldMat = nistManager->FindOrBuildMaterial("G4_Galactic");
  //
  // 1. Al fo the SpaceCraft; PEN:  13  Aluminium
  //G4Material* matPtr = nistManager->FindOrBuildMaterial(matName);
  fSCMat=nistManager->FindOrBuildMaterial("G4_Al");
  //
  // 2. Mo fo the Housing; PEN: 42  Molybdenum
  fHouMat=nistManager->FindOrBuildMaterial("G4_Mo");
  fHouShellMat=nistManager->FindOrBuildMaterial("G4_Mo");
  //
  // 3. Create Pt-Al mixture for the TestMass
  G4String name; //, symbol;             //a=mass of a mole;
  G4double density; //a, z, density;            //z=mean number of protons;
  G4int ncomponents; //, natoms;
  G4double fractionmass;

  density = 19.807*g/cm3; // from PEN sims
  fTMMat = new G4Material(name="Pt-Au"  , density, ncomponents=2);
  fTMMat->AddElement(nistManager->FindOrBuildElement("Au"), fractionmass=0.73);
  fTMMat->AddElement(nistManager->FindOrBuildElement("Pt"), fractionmass=0.27);
  
  // Create a replica of that material for describing a thin shell of the TM
  fTMshellMat = new G4Material(name="Pt-Au_shell"  , density, ncomponents=2);
  fTMshellMat->AddElement(nistManager->FindOrBuildElement("Au"), fractionmass=0.73);
  fTMshellMat->AddElement(nistManager->FindOrBuildElement("Pt"), fractionmass=0.27);


  //**********************************************
  // Registering G4Material to Penelope interface
  //**********************************************

    G4int mid = -1;  // Integer to get material index in PENELOPE execution.

    // Get PenInterface singleton
    PenInterface* penInterface = PenInterface::GetInstance();
    penInterface->SetVerbose(2);

    // Option A:
    // Use material properties defined in PENELOPE pre-defined material list.
    // See PENELOPE material list for further details.
    // It is important that G4Material* and PENELOPE material to be consistent,
    // so that G4 and PENELOPE particles can see the same material.
    //
    mid = penInterface->RegisterMaterial(fSCMat,13);
    penInterface->SetMSIMPA(mid,             // mat reg id
                5.e+4, 5.e+4, 5.e+4, //EABS(EGP)
                0.1, 0.1, // C1, C2
                5.e+4, 5.e+4); // WCC, WCR
    mid = penInterface->RegisterMaterial(fHouMat,42);
    penInterface->SetMSIMPA(mid,             // mat reg id
                5.e+4, 5.e+4, 5.e+4, //EABS(EGP)
                0.1, 0.1, // C1, C2
                5.e+4, 5.e+4); // WCC, WCR
    mid = penInterface->RegisterMaterial(fHouShellMat,42);
    penInterface->SetMSIMPA(mid,             // mat reg id
                5.e+1, 5.e+1, 5.e+1, //EABS(EGP)
                0.1, 0.1, // C1, C2
                5.e+1, 5.e+1); // WCC, WCR

    // Option B:
    // Build PENELOPE material using material properties defined in G4Material.
    // This is the only choice in case the material is not in PENELOPE database.
    //
    mid = penInterface->RegisterMaterial(fWorldMat);
    // Use minimal EABS for vacuum
    penInterface->SetMSIMPA(mid,             // mat reg id
                5.e+1, 5.e+1, 5.e+1, //EABS(EGP)
                0.05, 0.05, // C1, C2
                5.e+1, 5.e+1); // WCC, WCR
    mid = penInterface->RegisterMaterial(fTMMat); // from David's sims
    penInterface->SetMSIMPA(mid,             // mat reg id
                5.e+4, 5.e+4, 5.e+4, //EABS(EGP)
                0.1, 0.1, // C1, C2
                5.e+4, 5.e+4); // WCC, WCR
                
    mid = penInterface->RegisterMaterial(fTMshellMat); //
    penInterface->SetMSIMPA(mid,             // mat reg id
                5.e+1, 5.e+1, 5.e+1, //EABS(EGP)
                0.1, 0.1, // C1, C2
                5.e+1, 5.e+1); // WCC, WCR

  // In both cases, transport parameters MUST BE SET HERE, as done above,
  // as PenInterface::Initialize() has not been issued at this stage yet.
  //
  // NOTE: In case the same material needs to be used with *different* transport
  // parameters, then the G4Material* must point to a *different* objects,
  // so that they can be distinguished in PenInterface.
  // This means that G4NistMaterial::FindOrBuildMaterial() is NOT suitable
  // for such a purpose, as second issuance of this method passing same material
  // just retrieves the same G4Material object.
  // Thus, if such a case is needed (similar situation as if one would use
  // prod cuts per region), then it must be ensured that new G4Material* is
  // created with a different name. The simplest way is to use
  // the constructor that creates a material from base material.
  // Example assuming we need to copy 'mat1' properties:
  // G4Material* mat2 = new G4Material("newName", mat1->GetDensity(), mat1);

  // Reset verbosity to default
  penInterface->SetVerbose(1);
}

// Dani: changed the geometry to the standard Geant4 configuration (with daughter volumes)
void DetectorConstruction::SetupGeometry()
{
  // World
  G4VSolid* sWorld = new G4Box("World", fWorldSize,fWorldSize,fWorldSize);
  G4LogicalVolume* lWorld = new G4LogicalVolume(sWorld, fWorldMat, "World");
  fPhysWorld = new G4PVPlacement(0,G4ThreeVector(),lWorld,"World",0,false,0);

  //
  // Cubes created following Matrioshka logic and daughter volumes
  //
    
  // 1. Space craft
  G4VSolid* sSpaceCraft = new G4Box("sSpaceCraft",fSCSize/2.,fSCSize/2.,fSCSize/2.);
  G4LogicalVolume* lSpaceCraft = new G4LogicalVolume(sSpaceCraft, fSCMat,"SpaceCraft");
  G4ThreeVector positionSpaceCraft = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pSpaceCraft = new G4PVPlacement(0,positionSpaceCraft,"pSpaceCraft",lSpaceCraft,fPhysWorld,false,0);
  
  // 2. First Vacuum gap
  G4VSolid* sVacGap0 = new G4Box("sVacGap0",fSCSize/2.-fSCThick,fSCSize/2.-fSCThick,fSCSize/2.-fSCThick);
  G4LogicalVolume* lVacGap0 = new G4LogicalVolume(sVacGap0, fWorldMat,"VacGap0");
  G4ThreeVector positionVacGap0 = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pVacGap0 = new G4PVPlacement(0,positionVacGap0,"pVacGap0",lVacGap0,pSpaceCraft,false,0);
  
  // 3. Electrode housing
  G4VSolid* sHousing = new G4Box("sHousing",fHouSize/2.,fHouSize/2.,fHouSize/2.);
  G4LogicalVolume* lHousing = new G4LogicalVolume(sHousing, fHouMat,"Housing");
  G4ThreeVector positionHousing = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pHousing = new G4PVPlacement(0,positionHousing,"pHousing",lHousing,pVacGap0,false,0);

  // 4. Electrode housing inner shell
  G4VSolid* sHousingShell = new G4Box("sHousing",fHouSize/2.-fHouThick+fShellThick,fHouSize/2.-fHouThick+fShellThick,fHouSize/2.-fHouThick+fShellThick);
  G4LogicalVolume* lHousingShell = new G4LogicalVolume(sHousingShell, fHouShellMat,"HousingShell");
  G4ThreeVector positionHousingShell = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pHousingShell = new G4PVPlacement(0,positionHousingShell,"pHousingShell",lHousingShell,pHousing,false,0);

  // 5. Second Vacuum gap
  G4VSolid* sVacGap1 = new G4Box("sVacGap1",fHouSize/2.-fHouThick,fHouSize/2.-fHouThick,fHouSize/2.-fHouThick);
  G4LogicalVolume* lVacGap1 = new G4LogicalVolume(sVacGap1, fWorldMat,"VacGap1");
  G4ThreeVector positionVacGap1 = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pVacGap1 = new G4PVPlacement(0,positionVacGap1,"pVacGap1",lVacGap1,pHousingShell,false,0);
  
  // 6. Test Mass shell
  G4VSolid* sTestMassShell = new G4Box("sTestMassShell",fTMSize/2.,fTMSize/2.,fTMSize/2.);
  G4LogicalVolume* lTestMassShell = new G4LogicalVolume(sTestMassShell, fTMshellMat,"TestMassShell");
  G4ThreeVector positionTestMassShell = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pTestMassShell = new G4PVPlacement(0,positionTestMassShell,"pTestMassShell",lTestMassShell,pVacGap1,false,0);
  
  // 7. Test Mass
  G4VSolid* sTestMass = new G4Box("sTestMass",fTMSize/2.-fShellThick,fTMSize/2.-fShellThick,fTMSize/2.-fShellThick);
  G4LogicalVolume* lTestMass = new G4LogicalVolume(sTestMass, fTMMat,"TestMass");
  G4ThreeVector positionTestMass = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pTestMass = new G4PVPlacement(0,positionTestMass,"pTestMass",lTestMass,pTestMassShell,false,0);
  
  //
  // Spherical source using GPS
  //
    
  // David: No mother volume associated since its the source.
  //        World Material (G4_Galactic) associated to avoid possible interactions of the particles with the material.
  
  // G4Sphere* sHemisphere = new G4Sphere("sHemisphere", fSCSize+19.9*mm, fSCSize+20.*mm,
    //-90.*deg, 180.*deg, 0.*deg, 180.*deg);
    
  G4Sphere* sHemisphere = new G4Sphere("sHemisphere", 110.*mm, 110.1*mm,0.*deg, 360.*deg, 0.*deg, 180.*deg);
  G4LogicalVolume* lHemisphere = new G4LogicalVolume(sHemisphere, fWorldMat, "Hemisphere");
  G4ThreeVector positionSphereSource = G4ThreeVector(0.,0.,0.);
  G4VPhysicalVolume* pHemisphere = new G4PVPlacement(0, positionSphereSource, "pHemisphere", lHemisphere, fPhysWorld, false, 0);
  
  // Visualization
  G4VisAttributes* Att_Green = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  G4VisAttributes* Att_Gray= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;  // yellow(1.0,1.0,0.0));
  G4VisAttributes* Att_Blue= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* Att_White= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* Att_Red= new G4VisAttributes(G4Colour(1.0,0.,0.));
  G4VisAttributes* Att_Yellow= new G4VisAttributes(G4Colour(1.0,1.0,0.));
  G4VisAttributes* Att_Extra= new G4VisAttributes(G4Colour(0.,0.7,0.3));
  
  lWorld->SetVisAttributes(G4VisAttributes::Invisible);
  lVacGap0->SetVisAttributes(Att_Gray);
  lVacGap1->SetVisAttributes(Att_Gray);
  lSpaceCraft->SetVisAttributes(Att_Blue);
  lHousing->SetVisAttributes(Att_Red);
  lHousing->SetVisAttributes(Att_Red);
  lTestMassShell->SetVisAttributes(Att_Extra);
  lTestMass->SetVisAttributes(Att_Yellow);
  
  // Dani: some extra comments from the old code
  //
  // DSMAX (g4double) is the maximum step length for a given material.
  // It can be inserted as PenInterface::GetInstance()->RegisterDSMAX(logicVolumeName, dsmax);

}


//-------1---------2---------3---------4---------5---------6---------7---------8
//
void DetectorConstruction::RegisterPenG4Mat(G4String name, G4int penID)
{
  fPenMatVec.push_back( std::pair<G4String,G4int>(name, penID) );
}


//-------1---------2---------3---------4---------5---------6---------7---------8
//
void DetectorConstruction::
SetPenG4MatMSIMPA(G4double eabs1, G4double eabs2, G4double eabs3,
          G4double c1, G4double c2, G4double wcc, G4double wcr)
{
  // This method should be invoked in macro files of this application
  // right after RegisterPenG4Mat, to ensure that position of
  // registered PenG4Mat coincides with this MSIMPA parameters
  std::array<G4double,7> msimpa;
  msimpa[0] = eabs1;  // EABS[1]
  msimpa[1] = eabs2;  // EABS[2]
  msimpa[2] = eabs3;  // EABS[3]
  msimpa[3] = c1;  // C1
  msimpa[4] = c2;  // C2
  msimpa[5] = wcc;  // WCC
  msimpa[6] = wcr;  // WCR

  fPenMatMSIMPAVec.push_back(msimpa);
}
