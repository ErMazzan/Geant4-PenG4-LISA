//-----------------------------------------------------------------------------
//
// PenEMProcess.hh - version 2021.05.08
//
// A wrapper class for Penelope physics.
// Code added to consider soft energy loss of charged particles along step.
// So far, this is commented out - Thus soft edep is deposited at hinges.
//
//-----------------------------------------------------------------------------

#ifndef PenEMProcess_h
#define PenEMProcess_h 1

class PenInterface;
class PenPhys;
class G4Track;
class G4Step;
class G4VParticleChange;
class G4ParticleChange;
class G4ParticleDefinition;
#include "globals.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

class PenEMProcess : public G4VProcess
{
  public:     
      PenEMProcess(const G4String& processName = "PenelopeEM",
                    G4ProcessType theType = fParameterisation);
      virtual ~PenEMProcess();

  public:
      virtual void BuildPhysicsTable(const G4ParticleDefinition&);
      virtual void StartTracking(G4Track*);
      virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track&, G4double, G4ForceCondition*);
      virtual G4VParticleChange* PostStepDoIt( const G4Track&, const G4Step&);
      virtual G4bool IsApplicable(const G4ParticleDefinition&);

  public: // Methods not used
      virtual G4double AlongStepGetPhysicalInteractionLength(
                         const G4Track&, G4double, G4double, G4double&, G4GPILSelection*);
      virtual G4double AtRestGetPhysicalInteractionLength(
                         const G4Track&, G4ForceCondition*);
      virtual G4VParticleChange* AlongStepDoIt(
                         const G4Track&, const G4Step&);
      virtual G4VParticleChange* AtRestDoIt(
                         const G4Track&, const G4Step&);

  private:
      G4bool pWorthToTrack;
      PenInterface* pPenelope;
      PenPhys* pPenPhys;
      G4ParticleChange* pPChange;
};

#endif


