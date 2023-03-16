//-----------------------------------------------------------------------------
//
// PenPhys.hh - version 2021.05.08
//
// This is the class that encapsulates Penelope physics
// functions being invoked during the tracking phase.
// All public methods of this class are supposed to be
// invoked from PenEMProcess class.
//
//-----------------------------------------------------------------------------

#ifndef PENEPHYS_H
#define PENEPHYS_H 1

class G4Step;
class G4Track;
class G4ParticleChange;
class G4ParticleDefinition;
class PenInterface;

#include "globals.hh"
#include "G4ios.hh"
#include <vector>     //MACG needed for gcc 7.5.0

#include "PenelopeDefines.hh"

class PenPhys
{
 public:
  PenPhys();
  ~PenPhys();

 public: 
  // methods invoked from PenEMProcess
      void NewTrack(const G4Track*);
      G4bool WorthToTrack(const G4Track*);
      G4double StepLength(const G4Track*,G4double);
     //MACG void AlongStepDoIt(const G4Track*,const G4Step*, G4ParticleChange*);
      void PostStepDoIt(const G4Track*, G4ParticleChange*);
      void InitializeCEGRID();
     //MACG G4double GetSSoft() const { return sp;}    //MACG

 private:
  // methods internally used to invoke Penelope methods
      void SetTrack(const G4Track*);
      void NewVolume(const G4Track*);
      void UpdateTrack(G4ParticleChange*);
     //MACG void UpdateTrack(G4ParticleChange*, const G4ParticleDefinition* );
      void PushSecondaries(const G4Track*,G4ParticleChange*);
      void PrintProcessOccured(G4double,G4int);

 private:
  // Penelope methods used within this class
#include "local.h"

 private:
      G4ParticleDefinition* pElectron;
      G4ParticleDefinition* pGamma;
      G4ParticleDefinition* pPositron;
      PenInterface* pPenelope;
      G4int verbose;
      G4double sp;  // effective stopping power = DESOFT/step_length
      G4double dsmax;  // current value of DSMAX in cm (depends on volume)
      G4double carryOn;
      G4int prevMat;
      std::vector<G4Track*> secVec;

 private: // variables used to be defined in namespace TRACK_mod

  //  ****  Particle TRACK variables (to be initialised before calling
  //        subroutine START).  
  //  ----  Energy, position, direction, and weight.
      double E, X, Y, Z, U, V, W, WGHT;
  //  ----  Particle type, current body, and material.
      int KPAR, IBODY, MAT;
  //  ----  Particle history flags.
      int ILB[5];

  //  ****  Photon polarisation.
  //  ----  Polarised photons if IPOL=1, otherwise unpolarised photons.
      int IPOL /*= 0*/;  // Unpolarized photons by default.
  //  ----  Stokes parameters.
      double SP1, SP2, SP3;

  //  ****  The particle age (time elapsed since the start of the shower)
  //        is recorded when LAGE=.TRUE.
      bool LAGE /*= false*/;
      double PAGE /*= 0.0*/;

  //  ****  Random-hinge slowing-down parameters (output from subroutines
  //        JUMP and KNOCK).
  //
  //            <----------------- step ----------------->
  //           hard event -------- hinge ------- hard event
  //            <----- segment -----><----- segment ----->
  //
  //        MHINGE ... labels the two segments (jumps) of a step between
  //                   two hard events;
  //                   = 0 (1) for the segment before (after) the hinge,
  //        E0STEP ... energy at the beginning of the segment,
  //        DESOFT ... energy loss along the step,
  //        SSOFT .... effective stopping power, = DESOFT/step_length.
      int MHINGE;
      double E0STEP, DESOFT, SSOFT;

  private: // structs converted from namespaces
      CEGRID CEGRID_; // struct CEGRID is defined in PenelopeDefines.hh.
      struct CHIST
      {
        int ILBA[5];  
      } CHIST_;
      struct SECST
      {
        double ES[NMS], XS[NMS], YS[NMS], ZS[NMS], US[NMS], VS[NMS], WS[NMS],
               WGHTS[NMS], SP1S[NMS], SP2S[NMS], SP3S[NMS], PAGES[NMS];
        int KPS[NMS], IBODYS[NMS], MS[NMS], ILBS[5][NMS], IPOLS[NMS], NSEC;
      } SECST_;
      struct CJUMP1
      {
        double ELAST1, ELAST2;
        int KSOFTE, KSOFTI, KDELTA;
      } CJUMP1_;
      struct CJUMP0
      {
        double P[8], ST, DST, DSR, W1, W2, T1, T2;  
      } CJUMP0_;
      struct RSEED
      {
        int ISEED1, ISEED2;
        int* seed1=&ISEED1; 
        int* seed2=&ISEED2; 
      } RSEED_;


};

#endif
