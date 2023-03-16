//-------------------------------------------------------------------------------
//
// PenInterface.cpp for Penelope-2018
//
//  version 2018.09.30
//   -- Note: in this version, this source code file is included into PenPhys.cc
//      to avoid duplicated instantiation of global-scope namespace objects.
//  version 2021.05.06
//   -- Changes done to account properly EPMAX and DSMAX.
//
//-------------------------------------------------------------------------------

#include "PenInterface.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"

#include <vector>
#include <stdlib.h>
#include <unistd.h>

// Penelope methods

PenInterface* PenInterface::pObject = 0;

PenInterface* PenInterface::GetInstance()
{
  if(!pObject)
  { pObject = new PenInterface(); }
  return pObject;
}

PenInterface::PenInterface()
:initialized(false),verbose(1)
{
  //MACG epmin = 50.;  // default EMIN in eV - not needed at this scope
  epmax = 1.e9; // default EPMAX in eV
  maxnmat = MAXMAT; // maximum number of materials
  // Set default transport parameters. To modify them for a given  material
  // use SetMSIMPA() method when material is registered.
  for(int imat=0;imat<maxnmat;imat++)
  {
    PENELOPE_mod::EABS[0][imat] = 50.; // absorption energy for e- in eV
    PENELOPE_mod::EABS[1][imat] = 50.; // absorption energy for gamma in eV
    PENELOPE_mod::EABS[2][imat] = 50.; // absorption energy for e+ in eV
    PENELOPE_mod::C1[imat] = 0.01;     // average angular deflection
    PENELOPE_mod::C2[imat] = 0.01;     // maximum average fractional energy loss
    PENELOPE_mod::WCC[imat] = 50.;    // cutoff energy loss in eV for hard inelastic collision
    PENELOPE_mod::WCR[imat] = 50.;    // cutoff energy loss in eV for hard brems emission
    PENELOPE_mod::DEN[imat] = 1.0;     // material densities
    PENELOPE_mod::RDEN[imat] = 1.0;    // material density reciprocals
  }
}

PenInterface::~PenInterface()
{ ; }

void PenInterface::SetMSIMPA(G4int nmt,G4double e1,G4double e2,G4double e3,
              G4double cc1,G4double cc2,G4double wcc1,G4double wcr1)
{
  if(initialized)
  { G4Exception("PenInterface::SetMSIMPA","PenG4_0012",FatalException,
                "Penelope has already been initialized."); }
  if(nmt<1||nmt>maxnmat)
  { G4Exception("PenInterface::SetMSIMPA","PenG4_0002",FatalException,
                "material index must be [1,MAXMAT]."); }
  //MACG 2021.03.31 Bug fixed; these indexes are in [0,MAXMAT-1] (C-style)
  // and not in [1,MAXMAT].
  nmt--;
  //
  PENELOPE_mod::EABS[0][nmt] = e1; // absorption energy for e- in eV
  PENELOPE_mod::EABS[1][nmt] = e2; // absorption energy for gamma in eV
  PENELOPE_mod::EABS[2][nmt] = e3; // absorption energy for e+ in eV
  PENELOPE_mod::C1[nmt] = cc1;     // average angular deflection
  PENELOPE_mod::C2[nmt] = cc2;     // maximum average fractional energy loss
  PENELOPE_mod::WCC[nmt] = wcc1;   // cutoff energy loss in eV for hard inelastic collision
  PENELOPE_mod::WCR[nmt] = wcr1;   // cutoff energy loss in eV for hard brems emission
}

G4double PenInterface::GetEABS(G4int iParticle,G4int iMaterial)
{ return PENELOPE_mod::EABS[iParticle-1][iMaterial-1]; }


void PenInterface::RegisterDSMAX(const G4LogicalVolume* lv, G4double dsmax)
{
  dsmaxMap[lv] = dsmax;
  if (verbose > 0) {
    G4cout << "DSMAX[cm]= " << dsmax/cm << " registered for LV \""
	   << lv->GetName() << "\"." << G4endl;
  }
}


G4double PenInterface::GetDSMAX(const G4LogicalVolume* aLV)
{
  for (auto lvdsmax : dsmaxMap) {
    if ( lvdsmax.first == aLV )
      return lvdsmax.second;
  }
  G4double dsmax = 1.0e+20*cm;
  return dsmax;
}


G4int PenInterface::RegisterMaterial(const G4Material* mat,G4int matID)
{
  if((G4int)(matMap.size())==maxnmat)
  { G4Exception("PenInterface::RegisterMaterial","PenG4_0001",FatalException,
                "PenInterface::RegisterMaterial -- Too many materials!"); }
  matMap[mat] = matID;
  auto tmpid = PenMatID(mat); // returns the Penelope material index (1<=index<=maxnmap)
  if(verbose>1)
  { G4cout << "Material " << matMap.size() << " : " << mat->GetName()
           << " - Penelope material index = " << matID << " index " << tmpid << G4endl; }
  return tmpid; // returns the Penelope material index (1<=index<=maxnmap)
}

G4int PenInterface::PenMat(const G4Material* g4mat, G4bool abortIfUnknownMaterial)
{
  for(auto mat : matMap)
  { if(mat.first==g4mat) return mat.second; }
  if(abortIfUnknownMaterial)
  {
    G4ExceptionDescription msg;
    msg <<  "Material <" << g4mat->GetName() << "> is not registered.";
    G4Exception("PenInterface::PenMat","PenG4_0002",FatalException,msg);
  }
  return -1;
}

G4int PenInterface::PenMatID(const G4Material* g4mat, G4bool abortIfUnknownMaterial)
{
  G4int matID = 0;
  for(auto mat : matMap)
  {
    matID++; // Note: Penelope material ID starts from 1, not 0.
    if(mat.first==g4mat) return matID;
  }
  if(abortIfUnknownMaterial)
  {
    G4ExceptionDescription msg;
    msg <<  "Material <" << g4mat->GetName() << "> is not registered.";
    G4Exception("PenInterface::PenMatID","PenG4_0003",FatalException,msg);
  }
  return -1;
}

G4int PenInterface::PenMatID(G4String& g4matname, G4bool abortIfUnknownMaterial)
{
  G4int matID = 0;
  for(auto mat : matMap)
  {
    matID++; // Note: Penelope material ID starts with 1, not 0.
    if(mat.first->GetName()==g4matname) return matID;
  }
  if(abortIfUnknownMaterial)
  {
    G4ExceptionDescription msg;
    msg <<  "Material <" << g4matname << "> is not registered.";
    G4Exception("PenInterface::PenMatID","PenG4_0004",FatalException,msg);
  }
  return -1;
}

void PenInterface::SetEPMAX(G4double emx)
{
  if(initialized)
  {
    G4ExceptionDescription msg;
    msg <<  "Penelope has already been initialized with Emax = "<<epmax<<" (eV).";
    G4Exception("PenInterface::SetEPMax","PenG4_0014",FatalException,msg);
  }
  epmax = emx;
}

//MACG void PenInterface::SetEMin(G4double emn)
//MACG {
//MACG   if(initialized)
//MACG   {
//MACG     G4ExceptionDescription msg;
//MACG     msg <<  "Penelope has already been initialized with Emin = "<<epmin<<" (eV).";
//MACG     G4Exception("PenInterface::SetEMin","PenG4_0014",FatalException,msg);
//MACG   }
//MACG   epmin = emn;
//MACG }

void PenInterface::Initialize()
{
  if(initialized) return;

  G4String matFile = ".mat";
  G4String pwds = getenv("PWD");
  if(matMap.size()<1)
  { G4Exception("PenInterface::Initialize","PenG4_0005",FatalException,
                "No material is registered for Penelope."); }
  G4int matVec[MAXMAT];
  for(int ixix=0;ixix<MAXMAT;ixix++)
  { matVec[ixix] = -1; }
  G4int matVecIdx = 0;

  if(!getenv("PENDBASE_DIR"))
  { G4Exception("PenInterface::Initialize","PenG4_0006",FatalException,
    "Environment variable PENDBASE_DIR is not defined. Program abort..."); }
  G4String pendb = getenv("PENDBASE_DIR");
  G4cout << G4endl << "Creating Penelope material input data from " << pendb << G4endl;
  chdir(pendb.c_str());

  if(verbose>1) G4cout << "Number of materials : " << matMap.size() << G4endl;
  for(auto mat : matMap)
  {
    auto midx = mat.second;
    if(verbose>1) 
    { G4cout << "Material index : " << midx << " " << mat.first->GetName() << G4endl; }
    matVec[matVecIdx++] = midx;
    G4String pmdt = pwds;
    G4String pmfn = mat.first->GetName() + matFile;
    pmdt += "/";
    pmdt += pmfn;
    char fndt[80];
    G4cout << " Penelope Material File : " << pmdt << G4endl;
    std::strcpy(fndt,pmdt.c_str());
    FILE *IWR = fopen(fndt,"w");  //// MATERIAL DEFINITION FILE
    if(midx>0)
    { // Material pre-defined by PENELOPE
      PEMATZ(midx,IWR);
    }
    else if(midx==0)
    { // Material copied from G4
      auto g4mat = mat.first;
      G4String matnm = g4mat->GetName();
      char matnam[62] = "                                                             ";
      std::strcpy(matnam,matnm.c_str());
      int iatmn[30];
      double frmas[30];
      for(int iq=0;iq<30;iq++)
      {
        iatmn[iq]=0;
        frmas[iq]=0.0;
      }
      G4int nele = g4mat->GetNumberOfElements();
      if(nele>30)
      {
        G4cerr << "Definition of G4Material <" << matnm << "> has " << nele << " elements." << G4endl;
        G4Exception("PenInterface::Initialize","PenG4_0007",FatalException,
                    "Too many elements.");
      }
      auto frMassVec = g4mat->GetFractionVector();
      for(int iele=0;iele<nele;iele++)
      {
        auto elm = g4mat->GetElement(iele);
        iatmn[iele] = int(elm->GetZ());
        frmas[iele] = frMassVec[iele];
      }
      double dnst = (g4mat->GetDensity())/(g/cm3);
      PEMATC(matnam,nele,iatmn,frmas,dnst,IWR);
    }
    fclose(IWR); //// MATERIAL DEFINITION FILE
  }
  chdir(pwds.c_str());

  Initialize(&matFile);
}

void PenInterface::Initialize(G4String* matData)
{
  if(initialized) return;

  G4int nmat = matMap.size();        // number of materials
  G4int info = 3;        // PEINIT verbose level
  char pmfiles[MAXMAT][81];
  unsigned int imatidx = 0;
  for(auto mat : matMap)
  {
    G4String pmfn = mat.first->GetName();
    pmfn += *matData;
    std::strcpy(pmfiles[imatidx++],pmfn.c_str());
  }
  char fnwr[80] = "material.dat";
  FILE *IWR = fopen(fnwr,"w");  //// OUTPUT FILE
  PEINIT(epmax,nmat,IWR,info,pmfiles);

  G4cout << "PEINIT successfully called." << G4endl;
  initialized = true;
}

const CEGRID& PenInterface::GetGlobalCEGRID() 
{
  if(!initialized)
  { G4Exception("PenInterface::GetGlobalCEGRID()","PenInt0000",
    FatalException,"GetGlobalCEGRID() is called before initialized."); }
  return CEGRID_; 
}

/*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    PPPPP   EEEEEE  N    N  EEEEEE  L        OOOO   PPPPP   EEEEEE    C
C    P    P  E       NN   N  E       L       O    O  P    P  E         C
C    P    P  E       N N  N  E       L       O    O  P    P  E         C
C    PPPPP   EEEE    N  N N  EEEE    L       O    O  PPPPP   EEEE      C
C    P       E       N   NN  E       L       O    O  P       E         C
C    P       EEEEEE  N    N  EEEEEE  LLLLLL   OOOO   P       EEEEEE    C
C                                                                      C
C                                                   (version 2018).    C
C                                                                      C
C  Subroutine package for Monte Carlo simulation of coupled electron-  C
C  photon transport in homogeneous media.                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2018)                                     C
C  Copyright (c) 2001-2018                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided "as is" without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The convention used to name the interaction simulation subroutines is
C  the following:
C  - The first letter indicates the particle (E for electrons, P for
C    positrons, G for photons).
C  - The second and third letters denote the interaction mechanism
C    (EL for elastic, IN for inelastic, SI for inner-shell ionisation,
C    BR for bremsstrahlung, AN for annihilation, RA for Rayleigh, CO for
C    Compton, PH for photoelectric and PP for pair production).
C  - The fourth (lowercase) letter indicates the theoretical model
C    used to describe the interactions. This serves to distinguish the
C    default model (denoted by the letter 'a') from alternative models.
C  - The random sampling routines have four-letter names. Auxiliary
C    routines, which perform specific calculations, have longer names,
C    with the fifth and subsequent letters and/or numbers indicating
C    the kind of calculation (T for total x-section, D for differen-
C    tial x-section) or action (W for write data in a file, R for read
C    data from a file, I for initialisation of simulation algorithm).
C
C  The present subroutines may print warning and error messages in
C  the I/O unit 26. This is the default output unit in the example main
C  programs PENCYL and PENMAIN.
C
C  Subroutine PEMATW connects files to the I/O units 3 (input) and 7
C  (output). However, this does not conflict with the main program,
C  because PEMATW is not invoked during simulation. This subroutine is
C  used only by the program MATERIAL, to generate material data files.
C
C  Subroutine PEINIT connects material definition files to the input
C  unit 3; this unit is closed before returning to the main program.
C
*/

//  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//  *********************************************************************
//  *********************************************************************
//  **                      SHARED SUBROUTINES                         **
//  *********************************************************************
//  *********************************************************************

//  *********************************************************************
//                       SUBROUTINE PEINIT
//  *********************************************************************
#include "shareSubs/PEINIT.cpp"

//  *********************************************************************
//                       SUBROUTINE EGRID
//  *********************************************************************
#include "shareSubs/EGRID.cpp"

//  *********************************************************************
//                       SUBROUTINE PEMATR
//  *********************************************************************
#include "shareSubs/PEMATR.cpp"

//  *********************************************************************
//                       SUBROUTINE PEMATW
//  *********************************************************************
#include "shareSubs/PEMATW.cpp"
#include "shareSubs/pematc.cpp"
#include "shareSubs/pematz.cpp"

//  *********************************************************************
//                  FUNCTION PRANGE
//  *********************************************************************
#include "shareSubs/PRANGE.cpp"

//  *********************************************************************
//                  FUNCTION AVNCOL
//  *********************************************************************
#include "shareSubs/AVNCOL.cpp"

//  *********************************************************************
//                  FUNCTION PHMFP
//  *********************************************************************
#include "shareSubs/PHMFP.cpp"

//  *********************************************************************
//                  FUNCTION RYIELD
//  *********************************************************************
#include "shareSubs/RYIELD.cpp"

//  *********************************************************************
//                       SUBROUTINE EELaR
//  *********************************************************************
#include "shareSubs/EELaR.cpp"

//  *********************************************************************
//                       SUBROUTINE EELa0
//  *********************************************************************
#include "shareSubs/EELa0.cpp"

//  *********************************************************************
//                       SUBROUTINE EELaW
//  *********************************************************************
#include "shareSubs/EELaW.cpp"

//  *********************************************************************
//                       SUBROUTINE EINaT
//  *********************************************************************
#include "shareSubs/EINaT.cpp"

//  *********************************************************************
//                       SUBROUTINE EINaT1
//  *********************************************************************
#include "shareSubs/EINaT1.cpp"

//  *********************************************************************
//                       FUNCTION EINaDS
//  *********************************************************************
#include "shareSubs/EINaDS.cpp"

//  *********************************************************************
//                       SUBROUTINE PINaT
//  *********************************************************************
#include "shareSubs/PINaT.cpp"

//  *********************************************************************
//                       SUBROUTINE PINaT1
//  *********************************************************************
#include "shareSubs/PINaT1.cpp"

//  *********************************************************************
//                       FUNCTION PINaDS
//  *********************************************************************
#include "shareSubs/PINaDS.cpp"

//  *********************************************************************
//                       SUBROUTINE ESIaR
//  *********************************************************************
#include "shareSubs/ESIaR.cpp"

//  *********************************************************************
//                       SUBROUTINE ESIa0
//  *********************************************************************
#include "shareSubs/ESIa0.cpp"

//  *********************************************************************
//                       SUBROUTINE ESIaW
//  *********************************************************************
#include "shareSubs/ESIaW.cpp"

//  *********************************************************************
//                       SUBROUTINE PSIaR
//  *********************************************************************
#include "shareSubs/PSIaR.cpp"

//  *********************************************************************
//                       SUBROUTINE PSIa0
//  *********************************************************************
#include "shareSubs/PSIa0.cpp"

//  *********************************************************************
//                       SUBROUTINE PSIaW
//  *********************************************************************
#include "shareSubs/PSIaW.cpp"

//  *********************************************************************
//                       SUBROUTINE EBRaR
//  *********************************************************************
#include "shareSubs/EBRaR.cpp"

//  *********************************************************************
//                       SUBROUTINE EBRaW
//  *********************************************************************
#include "shareSubs/EBRaW.cpp"

//  *********************************************************************
//                       SUBROUTINE EBRaT
//  *********************************************************************
#include "shareSubs/EBRaT.cpp"

//  *********************************************************************
//                       SUBROUTINE PBRaT
//  *********************************************************************
#include "shareSubs/PBRaT.cpp"

//  *********************************************************************
//                       FUNCTION RLMOM
//  *********************************************************************
#include "shareSubs/RLMOM.cpp"

//  *********************************************************************
//                       SUBROUTINE RLPAC
//  *********************************************************************
#include "shareSubs/RLPAC.cpp"

//  *********************************************************************
//                       SUBROUTINE BRaAR
//  *********************************************************************
#include "shareSubs/BRaAR.cpp"

//  *********************************************************************
//                       SUBROUTINE BRaAW
//  *********************************************************************
#include "shareSubs/BRaAW.cpp"

//  *********************************************************************
//                       SUBROUTINE PANaT
//  *********************************************************************
#include "shareSubs/PANaT.cpp"

//  *********************************************************************
//                       SUBROUTINE GRAbT
//  *********************************************************************
#include "shareSubs/GRAbT.cpp"

//  *********************************************************************
//                       FUNCTION GRAaD
//  *********************************************************************
#include "shareSubs/GRAaD.cpp"

//  *********************************************************************
//                       FUNCTION GRAaF2
//  *********************************************************************
#include "shareSubs/GRAaF2.cpp"

//  *********************************************************************
//                       SUBROUTINE GRAaTI
//  *********************************************************************
#include "shareSubs/GRAaTI.cpp"

//  *********************************************************************
//                       SUBROUTINE GRAaR
//  *********************************************************************
#include "shareSubs/GRAaR.cpp"

//  *********************************************************************
//                       SUBROUTINE GRAaW
//  *********************************************************************
#include "shareSubs/GRAaW.cpp"

//  *********************************************************************
//                       SUBROUTINE GCOaT
//  *********************************************************************
#include "shareSubs/GCOaT.cpp"

//  *********************************************************************
//                        FUNCTION GCOaD
//  *********************************************************************
#include "shareSubs/GCOaD.cpp"

//  *********************************************************************
//                       SUBROUTINE GPHaT
//  *********************************************************************
#include "shareSubs/GPHaT.cpp"

//  *********************************************************************
//                       SUBROUTINE GPHaR
//  *********************************************************************
#include "shareSubs/GPHaR.cpp"

//  *********************************************************************
//                       SUBROUTINE GPHa0
//  *********************************************************************
#include "shareSubs/GPHa0.cpp"

//  *********************************************************************
//                       SUBROUTINE GPHaW
//  *********************************************************************
#include "shareSubs/GPHaW.cpp"

//  *********************************************************************
//                       SUBROUTINE GPPa0
//  *********************************************************************
#include "shareSubs/GPPa0.cpp"

//  *********************************************************************
//                       SUBROUTINE GPPaW
//  *********************************************************************
#include "shareSubs/GPPaW.cpp"

//  *********************************************************************
//                       SUBROUTINE PEILB4
//  *********************************************************************
#include "shareSubs/PEILB4.cpp"

//  *********************************************************************
//                       SUBROUTINE RELAXE
//  *********************************************************************
#include "shareSubs/RELAXE.cpp"

//  *********************************************************************
//                       SUBROUTINE RELAX0
//  *********************************************************************
#include "shareSubs/RELAX0.cpp"

//  *********************************************************************
//                       SUBROUTINE RELAXR
//  *********************************************************************
#include "shareSubs/RELAXR.cpp"

//  *********************************************************************
//                       SUBROUTINE RELAXW
//  *********************************************************************
#include "shareSubs/RELAXW.cpp"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//     *********************************
//     *  SUBROUTINE PACKAGE PENELAST  *
//     *********************************
//
//     The following subroutines perform class II simulation of elastic
//  scattering of electrons and positrons using the numerical cross
//  sections of the ELSEPA database, which covers the energy range from
//  50 eV to 100 MeV.
//
//   ****  DO NOT USE THESE SUBROUTINES FOR ENERGIES ABOVE 100 MeV  ****
//
//                                         Francesc Salvat. October 2004.
//
//  *********************************************************************
//                       SUBROUTINE EELdW
//  *********************************************************************
#include "shareSubs/EELdW.cpp"

//  *********************************************************************
//                       SUBROUTINE EELdR
//  *********************************************************************
#include "shareSubs/EELdR.cpp"

//  *********************************************************************
//                 SUBROUTINES ELINIT, DCSEL0 AND DCSEL
//  *********************************************************************
//
//     These subroutines read the ELSEPA database for elastic scattering
//  of electrons and positrons by neutral atoms, and generate the molecu-
//  lar DCS of a compound for arbitrary values of the energy (from 50 eV
//  to 100 MeV) and the angular deflection RMU=(1-COS(THETA))/2.
//
//  Other subprograms needed: subroutine SPLINE,
//                            function FINDI.
//
//  --> Subroutine ELINIT reads atomic elastic DCSs for electrons and
//  positrons from the database files and determines a table of the mole-
//  cular DCS, for the standard grid of energies and angular deflections,
//  as the incoherent sum of atomic DCSs. It is assumed that the database
//  files are in the same directory as the binary executable module. If
//  you wish to keep the database files on a separate directory, you can
//  edit the present source file and change the string FILE1, which con-
//  tains the names of the database files of the element, to include the
//  directory path.
//
//  --> Subroutine DCSEL0(E,IELEC) initializes the calculation of DCSs
//  for electrons (IELEC=-1) or positrons (IELEC=+1) with kinetic energy
//  E (eV). It builds a table of DCS values for the standard grid of
//  angular deflections RMU using log-log cubic spline interpolation in E
//  of the tables prepared by subroutine ELINIT.
//
//  --> Function DCSEL(RMU) computes the DCS in (cm**2/sr) at RMU by
//  linear log-log interpolation of the values tabulated by subroutine
//  DCSEL0. Notice that the delivered DCSEL value corresponds to the
//  particle and energy defined in the last call to subroutine DCSEL0.
//
//  EXAMPLE: To calculate cross sections for water, our main program must
//  contain the following definitions and calls:
//
//     CALL ELINIT(IZ,STF,NELEM) with NELEM=2
//                                    IZ(1)=1, STF(1)=2   <-- 2 H atoms
//                                    IZ(2)=8, STF(2)=1   <-- 1 O atom
//  (This sets the standard tabulation of cross sections for water.)
//
//  Now, to obtain the DCS for electrons with a kinetic energy of 10 keV
//  at RMU=0.0D0 (forward scattering) we simply insert the following two
//  statements in the main program,
//
//     CALL DCSEL0(1.0D4,-1)
//     DCS=DCSEL(0.0D0)

//  *********************************************************************
//                       SUBROUTINE ELINIT
//  *********************************************************************
#include "shareSubs/ELINIT.cpp"

//  *********************************************************************
//                       SUBROUTINE DCSEL0
//  *********************************************************************
#include "shareSubs/DCSEL0.cpp"

//  *********************************************************************
//                       FUNCTION DCSEL
//  *********************************************************************
#include "shareSubs/DCSEL.cpp"

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//                                                                      C
//  NUMERICAL TOOLS     (Francesc Salvat. Barcelona. February, 2012.)   C
//                                                                      C
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

//  *********************************************************************
//                       SUBROUTINE MERGE2
//  *********************************************************************
#include "shareSubs/MERGE2.cpp"

//  *********************************************************************
//                       SUBROUTINE SORT2
//  *********************************************************************
#include "shareSubs/SORT2.cpp"

//  *********************************************************************
//                       SUBROUTINE SPLINE
//  *********************************************************************
#include "shareSubs/SPLINE.cpp"

//  *********************************************************************
//                       SUBROUTINE FINDI
//  *********************************************************************
#include "shareSubs/FINDI.cpp"

//  *********************************************************************
//                       SUBROUTINE SINTEG
//  *********************************************************************
#include "shareSubs/SINTEG.cpp"

//  *********************************************************************
//                       SUBROUTINE SLAG6
//  *********************************************************************
#include "shareSubs/SLAG6.cpp"

//  *********************************************************************
//                       FUNCTION SUMGA
//  *********************************************************************
#include "shareSubs/SUMGA.cpp"

//  *********************************************************************
//                       FUNCTION RMOMX
//  *********************************************************************
#include "shareSubs/RMOMX.cpp"

//  *********************************************************************
//                       FUNCTION RNDG30
//  *********************************************************************
#include "shareSubs/RNDG30.cpp"

//  *********************************************************************
//                       FUNCTION RNDG3F
//  *********************************************************************
#include "shareSubs/RNDG3F.cpp"

//  *********************************************************************
//                       SUBROUTINE IRND0
//  *********************************************************************
#include "shareSubs/IRND0.cpp"

//  *********************************************************************
//                       SUBROUTINE RITA0
//  *********************************************************************
#include "shareSubs/RITA0.cpp"

//  *********************************************************************
//                       SUBROUTINE RITAI0
//  *********************************************************************
#include "shareSubs/RITAI0.cpp"

//  *********************************************************************
//                       SUBROUTINE RITAV
//  *********************************************************************
#include "shareSubs/RITAV.cpp"

//  *********************************************************************
//                       SUBROUTINE RITAM
//  *********************************************************************
#include "shareSubs/RITAM.cpp"

//  *********************************************************************
//                       SUBROUTINE ErrorFunction
//  *********************************************************************
#include "shareSubs/errors.cpp"



