
#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"


ActionInitialization
::ActionInitialization(DetectorConstruction* aDet)
{
  fDetector = aDet;
}


ActionInitialization::~ActionInitialization()
{;}

void ActionInitialization::BuildForMaster() const
{
  RunAction* runAction = new RunAction(fDetector);
  SetUserAction(runAction);
}

#include "PenPrimaryGeneratorAction.hh"

void ActionInitialization::Build() const
{
  SetUserAction(new PenPrimaryGeneratorAction);

  RunAction* runAction = new RunAction(fDetector);
  EventAction* evtAction = new EventAction(runAction);
  SteppingAction* stpAction = new SteppingAction(runAction);
  SetUserAction(runAction);
  SetUserAction(evtAction);
  SetUserAction(stpAction);
}
