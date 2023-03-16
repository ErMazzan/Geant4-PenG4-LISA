//-----------------------------------------------------------------------------
//
// PenInterface.hh - version 2021.05.06
//
// This is the class that encapsulates Penelope material and initialization
// functions being invoked during the initialization phase.
//
//-----------------------------------------------------------------------------
//
// Methods to be used by the user
//
// static PenInterface* GetInstance();
//
//  Returns the pointer to the PenInterface object. PenInterface must be a
//  singleton and the object is instantiated at the first use of this method.
//
// void RegisterDSMAX(const G4LogicalVolume*, G4double);
//
//  Register the value of DSMAX for a given logical volume. This method is
//  meant to be used in Detector Construction class. If none is set,
//  then default is assumed (1.0e+20 cm).
//
// G4double GetDSMAX(const G4LogicalVolume* );
//
//  Returns the value of DSMAX (length) for a given logical volume. This is
//  used by PenPhys class whenever a particle does its 1st step in a volume.
// 
// G4int RegisterMaterial(const G4Material* g4mat,G4int penmat=0);
//
//  Register the Geant4 material to be used in Penelope. "penmat" is the
//  pre-defined Penelope material number (1-280). See page 226 of Penelope2008
//  manual for the list of pre-defined materials. If "penmat" is not given,
//  Penelope material is created with the properties in G4Materal. In this
//  case, G4Material is assumed to consist of mixture of G4Element.
//
//  The return value is the material index number used in Penelope.
//  This index number should be used to set the material specific parameters
//  of Penelope.
//
// G4int PenMat(const G4Material* g4mat, G4bool abortIfUnknownMaterial=true);
//
//  The return value is the pre-defined Penelope material number (1-280).
//  See page 226 of Penelope2008 manual for the list of pre-defined materials.
//  If the material is created with G4Material properties, 0 is returned.
//  If given Geant4 material is not registered, -1 is returned if
//  "abortIfUnknownMaterial" is set to false. Otherwise it causes G4Exception.
//
// G4int PenMatID(const G4Material* g4mat, G4bool abortIfUnknownMaterial=true);
//
//  The return value is the material index number used in Penelope.
//  This index number should be used to set the material specific parameters
//  of Penelope. (1<=index<=MAXMAT)
//  If given Geant4 material is not registered, -1 is returned if
//  "abortIfUnknownMaterial" is set to false. Otherwise it causes G4Exception.
//
// G4int PenMatID(G4String& g4matname, G4bool abortIfUnknownMaterial=true);
//
//  Same as previous method, except the argument is now the name of the
//  material defined for G4Material.
//  If given Geant4 material is not registered, -1 is returned if
//  "abortIfUnknownMaterial" is set to false. Otherwise it causes G4Exception.
//
// void SetVerbose(G4int val);
//
//  Set verbosity of this interface with regard to material registrations.
//  If val>1 information is more detailed.
//
// void SetEPMAX(G4double emx);
// //MACG void SetEMin(G4double emn);  //MACG This is not set from outside.
//
//  Set the maximum/minimum kinetic energy (in eV) to be used in Penelope.
//  By default they are set to 1 GeV and 50 eV.
//  This must be set prior to Initialize() method.
//
// void SetMSIMPA(G4int matID,G4double eabs1,G4double eabs2,G4double eabs3,
//         G4double c1,G4double c2,G4double wcc,G4double wcr);
//
//  Set the material specific parameters of Penelope.
//  See details on page 281 of the Penelope2014 manual.
//    matID;    material index number (1-MAXMAT)
//    eabs1;    cut-off energy of e- in eV
//    eabs2;    cut-off energy of gamma in eV
//    eabs3;    cut-off energy of e+ in eV
//    c1;       average angular deflection
//    c2;       maximum fraction of energy loss
//    wcc;      cutoff energy loss in eV for hard inelastic collision
//    wcr;      cutoff energy loss in eV for hard brems emission
//  This must be set prior to Initialize() method. Ideally, right after
//  material registration.
//
//-----------------------------------------------------------------------------
//

#ifndef PenInterface_h
#define PenInterface_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <stdio.h>
#include <math.h>
#include <map>
class G4Material;
class G4LogicalVolume;

#include "PenelopeDefines.hh"

class PenInterface
{
  public:
      static PenInterface* GetInstance();

  private:
      static PenInterface* pObject;

  private:
      PenInterface();

  public:
      ~PenInterface();

  public:
      void RegisterDSMAX(const G4LogicalVolume*, G4double);
      G4double GetDSMAX(const G4LogicalVolume* );
      G4int RegisterMaterial(const G4Material*,G4int=0);

  private:
      G4int PenMat(const G4Material*,G4bool=true);

  public:
      G4int PenMatID(const G4Material*,G4bool=true);
      G4int PenMatID(G4String&,G4bool=true);

  public:
      void Initialize();
      void Initialize(G4String*);

  public:
      inline void SetVerbose(G4int val) { verbose = val; }
      void SetMSIMPA(G4int,G4double,G4double,G4double,
                     G4double,G4double,G4double,G4double);
      G4double GetEABS(G4int,G4int);

  private:
      G4bool initialized;
      G4int verbose;
      std::map<const G4Material*,G4int> matMap;
      std::map<const G4LogicalVolume*, G4double> dsmaxMap;

  public: // Penelope parameters. They must be set before Initialize()
      //MACG  void SetEMin(G4double emn); // min kinetic energy in eV
      void SetEPMAX(G4double emx); // max kinetic energy in eV

  private: // Penelope parameters
      //MACG  G4double epmin;   // min kinetic energy in eV
      G4double epmax;   // max kinetic energy in eV
      G4int    maxnmat; // maximum number of materials

  private:
      CEGRID CEGRID_;

  public:
      const CEGRID& GetGlobalCEGRID();

  public: // Penelope functions
#include "share.h"

};

#endif



