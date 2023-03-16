
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "EventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Step;
class G4Track;
class G4Event;
class RunAction;
class EventAction;
class PenOutputPrinter;


class SteppingAction : public G4UserSteppingAction
{
private:
  SteppingAction();
public:
  SteppingAction(RunAction* );
  ~SteppingAction();

  void UserSteppingAction(const G4Step*);
  void WritingFileCondition(const G4Event*);
    
public:
    
    G4int write = 0;
    G4int EnteringCharge;
    G4int ExitCharge;
    G4int SQCharge;
    G4int RootSQCharge;
private:
  // Method which classify tracks to fKill if they fulfill conditions set.
  void KillTrackCondition(G4Track*);

private: //DATA MEMBERS
    
  PenOutputPrinter* fOutputPrinter;
  RunAction* fRunAct;
  EventAction* fEventAct;
  G4double birthE;

};
#endif
