
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <fstream>

class RunAction;
class PenOutputPrinter;

class G4Event;

class EventAction : public G4UserEventAction
{
private:
  EventAction();

public:
  EventAction(RunAction* );
  ~EventAction();
    
public:
  //G4double fKinEnergy;
  //void SetKinEnergy(G4double ene)     { fKinEnergy = ene; }
  //G4double GetKinEnergy()         { return fKinEnergy; }
    
    
public:
  // virtual methods from G4UserEventAction

  void BeginOfEventAction(const G4Event* );
  void EndOfEventAction(const G4Event* );
    
private:

  // Data members
  PenOutputPrinter* fOutputPrinter;
  RunAction* fRunAction;
  G4int genParticles = 1;

};
#endif
