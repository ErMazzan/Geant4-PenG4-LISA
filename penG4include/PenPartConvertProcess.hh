//-----------------------------------------------------------------------------
//
// PenPartConvertProcess.hh
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#ifndef PenPartConvertProcess_h
#define PenPartConvertProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"
class PenInterface;

class PenPartConvertProcess : public G4VDiscreteProcess 
{
  public:
      PenPartConvertProcess(G4double te,const G4String& processName ="PenPartSwitch");
      virtual ~PenPartConvertProcess();

  private:
      PenPartConvertProcess(const PenPartConvertProcess &);
      PenPartConvertProcess & operator=(const PenPartConvertProcess &);

  public:
     virtual G4VParticleChange* PostStepDoIt(const G4Track&,const G4Step&);
     virtual void BuildPhysicsTable(const G4ParticleDefinition&); 
     virtual G4bool IsApplicable(const G4ParticleDefinition&);
     virtual G4double PostStepGetPhysicalInteractionLength(
        const G4Track& track,G4double previousStepSize,G4ForceCondition* condition);
     virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

  public:
     void SetVerboseLevel(G4int value)
     { verboseLevel = value; }
     G4int GetVerboseLevel() const
     { return verboseLevel; }
     void SetThresholdEnergy(G4double te)
     { threshE = te; }
     G4double GetThresholdEnergy() const
     { return threshE; }

  private:
     G4int verboseLevel;
     G4double threshE;
     PenInterface* pPenInterface;

  private:
     G4ParticleChangeForDecay fParticleChangeForDecay;
};

#endif

