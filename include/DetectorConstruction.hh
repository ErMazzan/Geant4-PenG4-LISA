
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <array>
#include <map>
#include <vector>
#include <utility> // std::pair

class G4String;
class G4Material;
class G4VPhysicalVolume;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();

// function inherited from PenG4 example used to register materials
  void RegisterPenG4Mat(G4String name, G4int =0);
  void SetPenG4MatMSIMPA(G4double, G4double, G4double, G4double,
			 G4double, G4double, G4double);

  G4VPhysicalVolume* GetPhysWorld() const  { return fPhysWorld; }

private:
  void DefineMaterials();
  void SetupGeometry();

private:
  DetectorMessenger* fMessenger;
  // Include each G4Material* registered for PenelopeEMPhysics
  G4Material* fWorldMat;
  G4VPhysicalVolume* fPhysWorld;

  std::vector<G4Material*> fMatVec;
  std::vector< std::pair<G4String,G4int> > fPenMatVec;
  std::vector< std::array<G4double, 7> > fPenMatMSIMPAVec;
  //EBAS[1-3],C1,C2,WCC,WCR

  G4bool fConstructed;
  
  G4double fWorldSize;
  
  G4double fSCSize;
  G4Material* fSCMat;
  G4double fSCThick;
  
  G4double fHouSize;
  G4Material* fHouMat;
  G4Material* fHouShellMat;
  G4double fHouThick;
  
  G4Material* fTMshellMat;
  G4double fShellThick;
  
  G4double fTMSize;
  G4Material* fTMMat;
  
};

#endif
