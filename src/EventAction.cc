
#include "EventAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "g4root.hh"

#include "RunAction.hh"
#include "PenOutputPrinter.hh"

/*
EventAction::EventAction()
{
  G4cout << "void EventAction constructor is not meant to be used!."
	 << G4endl;
}
*/

EventAction::EventAction(RunAction* runAction) : fRunAction(runAction)
{
  fOutputPrinter = runAction->GetOutputPrinter();
}


EventAction::~EventAction()
{
  // fOutputPrinter is deleted by RunAction!
}


//-----------------------------------------------------------------------------
void EventAction::BeginOfEventAction(const G4Event* evt)
{
  if (evt->GetEventID()%100 == 0) G4cout << "<<< Begin of Event  " << evt->GetEventID() << G4endl;

  fOutputPrinter->ResetEventInfo();
  fRunAction->ResetChargeCounter();
  fRunAction->ResetEnteringCharge();
  fRunAction->ResetExitingCharge();
  
  // fKinEnergy = 0.;
    
  /*
    G4int runID = fRunAction->GetRunID();
   
  // This allow us to reset the number of generated particles in each run

  if (runID == 0) {
     genParticles = 1;
     runID = 1;
  }
  
  G4cout << "Particles generated: " << genParticles << G4endl;
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(0,1);
  analysisManager->AddNtupleRow();    
  */
   //std::cout << "here"<<std::endl;
    
  
}


//-----------------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event* evt)
{
	
    genParticles += 1;
    
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    G4int ChargePerEvent = fRunAction->GetDepositedCharge();
    G4double KinEnergy = fRunAction->GetKinEnergy();
    
    G4double EnteringCharge = fRunAction->GetEnteringCharge();
    G4double ExitingCharge = fRunAction->GetExitingCharge();
    
    /*
    G4int SquareChargePerEvent = ChargePerEvent*ChargePerEvent;
    G4int SquareRootChargePerEvent = sqrt(SquareChargePerEvent);
    
    */
    
	analysisManager->FillH1(0,ChargePerEvent);
    
    G4cout << "Deposited Charge in TM per event: " << ChargePerEvent << G4endl;
    
    if (ChargePerEvent != 0) {
        analysisManager->FillH2(4, KinEnergy, EnteringCharge);
        analysisManager->FillH2(5, KinEnergy, ExitingCharge);
        analysisManager->FillH1(11, ChargePerEvent);
        analysisManager->FillH1(12, KinEnergy);
    }
    
    fRunAction->UpdateChargeCounterPerRun(ChargePerEvent);
    
    /*
    fRunAction->UpdateSquareCharge(SquareChargePerEvent);
    fRunAction->UpdateSquareRootCharge(SquareRootChargePerEvent);
    */
    
    // fRunAction->SetRunID(1);
}
