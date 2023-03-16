//-------------------------------------------------------------------------------
//
// PenPhys.cc for Penelope-2018
//
//  version 2021.05.08
//   -- Note: in this version, the source code file of PenInterface.cpp is
//      included at the bottom of this file to avoid duplicated instantiation
//      of global-scope namespace objects.
//   -- Note: The value of dsmax is set from the logical volume and passed to
//      JUMP subroutine.
//   -- Note: The soft energy loss could be considered along step instead of at
//      hinges. This is key for very thin slabs. This code is commented out.
//
//-------------------------------------------------------------------------------

#include "globals.hh"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "PenPhys.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "PenInterface.hh"

#include "PenGamma.hh"
#include "PenElectron.hh"
#include "PenPositron.hh"

#include "G4Threading.hh"

//==================================================================
//
// PENELOPE namespaced global variables
//
#include "common-share.h"
//
//==================================================================

PenPhys::PenPhys()
  :verbose(1),sp(0.),dsmax(1.0e+20),carryOn(0.),prevMat(-1)
{
  pElectron = PenElectron::Electron();
  pGamma = PenGamma::Gamma();
  pPositron = PenPositron::Positron();
  pPenelope = PenInterface::GetInstance();

// Following line is needed if the original Penelope RAND() is used.
//  RAND0(G4Threading::G4GetThreadId()+1);

  IPOL = 0;
  LAGE = false;
  PAGE = 0.0;
}

PenPhys::~PenPhys()
{;}

void PenPhys::NewTrack(const G4Track*)
{
  prevMat = -1;
  carryOn = 0.;
  secVec.clear();
  CLEANS();
}

void PenPhys::NewVolume(const G4Track*)
{
  START();
}

G4double PenPhys::StepLength(const G4Track* aTrack, G4double /*previousStepSize*/)
{
  SetTrack(aTrack);

  G4double ds = 0.;

  // Check value of dsmax when entering into a new volume or doing 1st step
  if (aTrack->GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ||
      aTrack->GetCurrentStepNumber() == 1) {
    dsmax = (pPenelope->GetDSMAX( aTrack->GetVolume()->GetLogicalVolume() ))/cm;
    if (verbose > 1) {
      G4cout << "dsmax[cm] = " << dsmax << " : lv name= "
	     << aTrack->GetVolume()->GetLogicalVolume()->GetName()
	     << " : STEP No.: " << aTrack->GetCurrentStepNumber() << G4endl;
    }
  }
  
  JUMP(dsmax,ds);

  sp =  SSOFT*eV/cm;
  if(verbose>1)
  {
    G4cout << "++++ Step length proposed " << ds << " [cm], "
           << "effective stopping power " << sp/eV*cm << " [eV/cm]" << G4endl;
  }

  return ds*cm;
}

G4bool PenPhys::WorthToTrack(const G4Track* aTrack)
{
  G4int ipar = 0;
  auto pp = aTrack->GetDefinition();
  if(pp==pElectron)
  { ipar = 1; }
  else if(pp==pGamma)
  { ipar = 2; }
  else if(pp==pPositron)
  { ipar = 3; }
  else
  { G4Exception("PenPhys::SetTrack","PenG4_0008",FatalException,"Unknown particle type"); }

  auto imat = pPenelope->PenMatID(aTrack->GetMaterial());

  return (aTrack->GetKineticEnergy() > (pPenelope->GetEABS(ipar,imat))*eV);
}

void PenPhys::SetTrack(const G4Track* aTrack)
{
  E = (aTrack->GetKineticEnergy())/eV;
  X = (aTrack->GetPosition().x())/cm;
  Y = (aTrack->GetPosition().y())/cm;
  Z = (aTrack->GetPosition().z())/cm;
  U = aTrack->GetMomentumDirection().x();
  V = aTrack->GetMomentumDirection().y();
  W = aTrack->GetMomentumDirection().z();
  WGHT = aTrack->GetWeight();

  G4int ipar = 0;
  G4ParticleDefinition* pp = aTrack->GetDefinition();
  if(pp==pElectron)
  { ipar = 1; }
  else if(pp==pGamma)
  { ipar = 2; }
  else if(pp==pPositron)
  { ipar = 3; }
  else
  { G4Exception("PenPhys::SetTrack","PenG4_0008",FatalException,"Unknown particle type"); }

  KPAR = ipar;
  IBODY = 1;
  MAT = pPenelope->PenMatID(aTrack->GetMaterial());
  if(MAT!=prevMat)
  {
    NewVolume(aTrack);
    prevMat =  MAT;
  }
  ILB[0] = 1;
  ILB[1] = 0;
  ILB[2] = 0;
  ILB[3] = 0;
  ILB[4] = 0;

  if(verbose>1)
  {
    G4cout << "SetTrack > kpar = " <<  KPAR << " in material " <<  MAT
           << " -- ekin " <<  E << G4endl
           << "      at " << G4ThreeVector(X, Y, Z)
           << " dir " << G4ThreeVector(U, V, W) << G4endl;
  }
}


//MACG
// void PenPhys::AlongStepDoIt(const G4Track* aTrack, const G4Step* aStep, G4ParticleChange* aPChange)
// {
//   aPChange->Initialize(*aTrack);
  
//   // Account the soft energy loss along step instead of at the hinge
//   if (aTrack->GetDefinition() != PenGamma::Definition() ) {
//     G4double dE = (this->GetSSoft())*(aStep->GetStepLength());
//     if (dE > 0) {
//       aPChange->ProposeLocalEnergyDeposit(dE);
//       E -= dE/eV;    // To pass proper E value to KNOCK
//       aPChange->ProposeEnergy(E*eV);
//     }
//     if ( verbose > 1) {
//       G4cout << "+++ PenPhys::AlongStepDoIt: " << G4endl;
//       G4cout << "\tTrueStepLength[cm]= " << aStep->GetStepLength()/cm
// 	     << "\tSSoft[eV/cm]= " << this->GetSSoft()/eV*cm << G4endl;
//       G4cout << "\tdE[eV]= " << dE/eV << G4endl;
//     }
//   }
// }


void PenPhys::PostStepDoIt(const G4Track* aTrack, G4ParticleChange* aPChange)
{
  aPChange->Initialize(*aTrack);
  if (verbose > 1) {
    G4cout << "\tE priorPostDoIt[eV]= " << aPChange->GetEnergy()/eV << G4endl;
  }

  UpdateTrack(aPChange);
  //MACG UpdateTrack(aPChange, aTrack->GetDefinition());
  PushSecondaries(aTrack,aPChange);
  G4double de = aPChange->GetLocalEnergyDeposit() - carryOn;
  if(de>0.)
  {
     aPChange->ProposeLocalEnergyDeposit(de);
     carryOn = 0.;
  }
  else
  {
     aPChange->ProposeLocalEnergyDeposit(0.);
     carryOn = -de;
  }
  if (verbose > 1) {
    G4cout << "\tPostStep TrueEdep[eV]= "
	   << aPChange->GetLocalEnergyDeposit()/eV << G4endl
	   << "\tTrue StepLength[cm]= "
	   << aTrack->GetStep()->GetStepLength()/cm << G4endl;
  }
}

//MACG void PenPhys::UpdateTrack(G4ParticleChange* aPChange,
//MACG 			  const G4ParticleDefinition* pDef)
void PenPhys::UpdateTrack(G4ParticleChange* aPChange)
{
  using namespace PENELOPE_mod;

  G4double de = 0.;
  G4int icol = 0;
  KNOCK(de,icol);
  if(verbose>1) PrintProcessOccured(de,icol);

  // if (pDef != PenGamma::Definition() && icol == 1) {
  //   // hinge Edep is AlongStep for e-/e+
  //   E += de;
  //   de = 0;
  // }

  if( E > EABS[KPAR-1][MAT-1])
  { // track survives
    aPChange->ProposeEnergy(  E*eV );
    aPChange->ProposeLocalEnergyDeposit( de*eV );
    aPChange->ProposeMomentumDirection(  U,  V,  W );
    aPChange->ProposeWeight(  WGHT );
  }
  else
  { // track is dead
    aPChange->ProposeTrackStatus(fStopAndKill);
    aPChange->ProposeEnergy(0.);
    aPChange->ProposeLocalEnergyDeposit( (de + E)*eV );
  }

  if(verbose>1)
  {
    G4cout << "UpdateTrack > kpar = " <<  KPAR
           << " -- ekin " <<  E << G4endl
           << "      at " << G4ThreeVector( X, Y, Z)
           << " dir " << G4ThreeVector( U, V, W) << G4endl;
  }
}

void PenPhys::PushSecondaries(const G4Track* aTrack,G4ParticleChange* aPChange)
{
  using namespace PENELOPE_mod;

  G4int left = 0;
  SECPAR(left);
  G4double de = aPChange->GetLocalEnergyDeposit();

  if(verbose>1)
  { G4cout << " PenPhys::PushSecondaries -- " << left << " secondaries created." << G4endl; }
  while(left)
  {
    if( E > EABS[KPAR-1][MAT-1])
    {
      G4ParticleDefinition* pd = 0;
      switch( KPAR)
      {
       case 1:
        pd = pElectron; break;
       case 2:
        pd = pGamma; break;
       case 3:
        pd = pPositron; break;
      }
      G4DynamicParticle* dp
        = new G4DynamicParticle(pd,G4ThreeVector( U, V, W), E*eV);
      G4Track* sec
        = new G4Track(dp,aPChange->GetGlobalTime(),*(aPChange->GetPosition()));
      sec->SetTouchableHandle(aTrack->GetTouchableHandle());
      sec->SetWeight( WGHT);
      de -=  E*eV;
      secVec.push_back(sec);
      if(verbose>1)
      {
        G4cout << "Secondary > kpar = " <<  KPAR
               << " -- ekin " <<  E << G4endl
               << "      at " << G4ThreeVector( X, Y, Z)
               << " dir " << G4ThreeVector( U, V, Z) << G4endl;
      }
    }
    else
    {
      if(verbose>1)
      {
        G4cout << "Secondary NOT GENERATED > kpar = " <<  KPAR
               << " -- ekin " <<  E << G4endl
               << "      at " << G4ThreeVector( X, Y, Z)
               << " dir " << G4ThreeVector( U, V, Z) << G4endl;
      }
    }
    SECPAR(left);
  }

  aPChange->SetNumberOfSecondaries(secVec.size());
  for(size_t i=0;i<secVec.size();i++)
  { aPChange->AddSecondary(secVec[i]); }
  aPChange->ProposeLocalEnergyDeposit(de);
  secVec.clear();
  CLEANS();
}

void PenPhys::PrintProcessOccured(G4double de,G4int icol)
{
  static G4String eInt[8] = {
   "artificial soft event (random hinge)","hard elastic collision","hard inelastic collision",
   "hard bremsstrahlung emission","inner-shell impact ionization","n/a",
   "delta interaction","auxiliary interaction" };
  static G4String gInt[8] = {
   "coherent (Rayleigh) scattering","incoherent (Compton) scattering","photoelectric absorption",
   "electron-positron pair production","n/a","n/a",
   "delta interaction","auxiliary interaction" };
  static G4String pInt[8] = {
   "artificial soft event (random hinge)","hard elastic collision","hard inelastic collision",
   "hard bremsstrahlung emission","inner-shell impact ionization","annihilation",
   "delta interaction","auxiliary interaction" };

  G4cout << G4endl;
  G4cout << "PenPhys::PrintProcessOccured() - " << de << " [eV] : " << icol;
  switch ( KPAR)
  {
   case 1:
    G4cout << " : electron " << eInt[icol-1]; break;
   case 2:
    G4cout << " : gamma " << gInt[icol-1]; break;
   case 3:
    G4cout << " : positron " << pInt[icol-1]; break;
   default:
    G4cerr << " kpar = " <<  KPAR << G4endl;
    G4Exception("PenPhys::PenMat","PenG4_0002",FatalException,
                "PenPhys::PrintProcessOccured() --- Unknown particle type");
  }
  //MACG
  // if (icol == 1 && KPAR != 2)
  //   G4cout << " : Soft energy loss counted along step, not at hinge!";

  G4cout << G4endl << G4endl;
}

void PenPhys::InitializeCEGRID()
{
  const CEGRID CEGRID_global = pPenelope->GetGlobalCEGRID();
  CEGRID_.EMIN = CEGRID_global.EMIN;
  CEGRID_.EL = CEGRID_global.EL;
  CEGRID_.EU = CEGRID_global.EU;
  CEGRID_.DLEMP1 = CEGRID_global.DLEMP1;
  CEGRID_.DLFC = CEGRID_global.DLFC;
  for(int i=0; i<NEGP; i++)
  {
    CEGRID_.ET[i] = CEGRID_global.ET[i];
    CEGRID_.DLEMP[i] = CEGRID_global.DLEMP[i];
  }
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

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//                                                                      X
//    TRANSPORT ROUTINES (Francesc Salvat. Barcelona. May, 2018.)    X
//                                                                      X
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

//  *********************************************************************
//                       SUBROUTINE JUMP
//  *********************************************************************
#include "localSubs/JUMP.cpp"

//  *********************************************************************
//                       SUBROUTINE KNOCK
//  *********************************************************************
#include "localSubs/KNOCK.cpp"

//  *********************************************************************
//                       SUBROUTINE EELa
//  *********************************************************************
#include "localSubs/EELa.cpp"

//  *********************************************************************
//                       SUBROUTINE EINa
//  *********************************************************************
#include "localSubs/EINa.cpp"

//  *********************************************************************
//                       SUBROUTINE PINa
//  *********************************************************************
#include "localSubs/PINa.cpp"

//  *********************************************************************
//                       SUBROUTINE ESIa
//  *********************************************************************
#include "localSubs/ESIa.cpp"

//  *********************************************************************
//                       SUBROUTINE PSIa
//  *********************************************************************
#include "localSubs/PSIa.cpp"

//  *********************************************************************
//                       SUBROUTINE EBRa
//  *********************************************************************
#include "localSubs/EBRa.cpp"

//  *********************************************************************
//                       SUBROUTINE EBRaA
//  *********************************************************************
#include "localSubs/EBRaA.cpp"

//  *********************************************************************
//                       SUBROUTINE PANaR 
//  *********************************************************************
#include "localSubs/PANaR.cpp"

//  *********************************************************************
//                       SUBROUTINE PANa
//  *********************************************************************
#include "localSubs/PANa.cpp"

//  *********************************************************************
//                       SUBROUTINE GRAa
//  *********************************************************************
#include "localSubs/GRAa.cpp"

//  *********************************************************************
//                       SUBROUTINE GCOa
//  *********************************************************************
#include "localSubs/GCOa.cpp"

//  *********************************************************************
//                       SUBROUTINE GPHa
//  *********************************************************************
#include "localSubs/GPHa.cpp"

//  *********************************************************************
//                       SUBROUTINE SAUTER
//  *********************************************************************
#include "localSubs/SAUTER.cpp"

//  *********************************************************************
//                       SUBROUTINE GPPa
//  *********************************************************************
#include "localSubs/GPPa.cpp"

//  *********************************************************************
//                       SUBROUTINE SCHIFF
//  *********************************************************************
#include "localSubs/SCHIFF.cpp"

//  *********************************************************************
//                       SUBROUTINE RELAX
//  *********************************************************************
#include "localSubs/RELAX.cpp"

//  *********************************************************************
//                       SUBROUTINE PELd
//  *********************************************************************
#include "localSubs/PELd.cpp"

//  *********************************************************************
//                       FUNCTION RNDG3
//  *********************************************************************
#include "localSubs/RNDG3.cpp"

//  *********************************************************************
//                       SUBROUTINE CLEANS
//  *********************************************************************
#include "localSubs/CLEANS.cpp"

//  *********************************************************************
//                       SUBROUTINE START
//  *********************************************************************
#include "localSubs/START.cpp"

//  *********************************************************************
//                       SUBROUTINE DIRECT
//  *********************************************************************
#include "localSubs/DIRECT.cpp"

//  *********************************************************************
//                       SUBROUTINE DIRPOL
//  *********************************************************************
#include "localSubs/DIRPOL.cpp"

//  *********************************************************************
//                       SUBROUTINE STORES
//  *********************************************************************
#include "localSubs/STORES.cpp"

//  *********************************************************************
//                       SUBROUTINE SECPAR
//  *********************************************************************
#include "localSubs/SECPAR.cpp"

//  *********************************************************************
//                       SUBROUTINE EIMFP
//  *********************************************************************
#include "localSubs/EIMFP.cpp"

//  *********************************************************************
//                       SUBROUTINE PIMFP
//  *********************************************************************
#include "localSubs/PIMFP.cpp"

//  *********************************************************************
//                       SUBROUTINE GIMFP
//  *********************************************************************
#include "localSubs/GIMFP.cpp"

//  *********************************************************************
//                       SUBROUTINE EAUX
//  *********************************************************************
void PenPhys::EAUX()
{
  //  Auxiliary interaction mechanism for electrons, definable by the user.
  //  Usually it is not active.

///////  using namespace PENELOPE_mod;
///////  using namespace TRACK_mod;
  
  printf(" Warning: Subroutine EAUX has been entered.\n");
}

//  *********************************************************************
//                       SUBROUTINE PAUX
//  *********************************************************************
void PenPhys::PAUX()
{
  //  Auxiliary interaction mechanism for positrons, definable by the user.
  //  Usually it is not active.

///////  using namespace PENELOPE_mod;
///////  using namespace TRACK_mod;
  
  printf(" Warning: Subroutine PAUX has been entered.\n");
}

//  *********************************************************************
//                       SUBROUTINE GAUX
//  *********************************************************************
void PenPhys::GAUX()
{
  //  Auxiliary interaction mechanism for photons, definable by the user.
  //  Usually it is not active.

///////  using namespace PENELOPE_mod;
///////  using namespace TRACK_mod;

  printf(" Warning: Subroutine GAUX has been entered.\n");
}

//  *********************************************************************
//                       SUBROUTINE EELd
//  *********************************************************************
#include "localSubs/EELd.cpp"

//  *********************************************************************
//                         SUBROUTINE RAND0
//  *********************************************************************
#include "localSubs/RAND0.cpp"

//  *********************************************************************
//                         FUNCTION RAND (Random number generator)
//  *********************************************************************
#include "localSubs/RAND.cpp"

#include "PenInterface.cpp"

