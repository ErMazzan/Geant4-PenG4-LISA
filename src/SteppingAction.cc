
#include "SteppingAction.hh"
#include "EventAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"

#include "G4RunManager.hh"
#include "g4root.hh"

#include "RunAction.hh"
#include "PenOutputPrinter.hh"

const std::string txtred("\033[0;31m");
const std::string txtgreen("\033[0;32m");
const std::string txtyellow("\033[0;33m");
const std::string txtblue("\033[0;34m");
const std::string txtreset("\033[0m");


SteppingAction::SteppingAction()
{
  G4cout << "void SteppingAction constructor is not meant to be used!."
	 << G4endl;
}


SteppingAction::SteppingAction(RunAction* runAct):  fRunAct(runAct)
{
  fOutputPrinter = runAct->GetOutputPrinter();
  //birthE = 0.;
}


SteppingAction::~SteppingAction()
{
  // fOutputPrinter is deleted within RunAction!
}

//-----------------------------------------------------------------------------
void SteppingAction::WritingFileCondition(const G4Event* evt)
{
    G4int EventID = evt->GetEventID();
    if (EventID%100 == 0) write = 1;
    else if (EventID%100 == 1) write = 0;
    G4cout << write << G4endl;
}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //  G4Track* aTrack = aStep->GetTrack();
  
  // Dump step info into output files and temporal variables
  // TO DECIDE MT-APPROACH
  // fOutputPrinter->PrintStepInfo(aStep);
    
  fOutputPrinter->StoreStepInfo(aStep);
 
  // Some definitions
  const G4String CriticalVolume = "pTestMassShell"; //Volume in which we will look for entering/exiting particles
  
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4String currentPhysicalName = prePoint->GetPhysicalVolume()->GetName();
	//const G4String ProcessName = endPoint->GetProcessDefinedStep()->GetProcessName(); // this may give problems if nex process is not defined
  
  G4Track* theTrack = aStep->GetTrack();
	G4int track_id = theTrack->GetTrackID();
    G4int parent_id = theTrack->GetParentID();
	G4double trackLength  = theTrack->GetTrackLength();
	G4double kinE  = theTrack->GetKineticEnergy();
	const G4ParticleDefinition* particle = theTrack->GetDefinition();
		const G4String particleName= particle->GetParticleName();
    
  
  //  auto posx = aStep->GetPreStepPoint()->GetPosition();
    
  //  G4cout << "X: " << posx.x() << G4endl;
    

  if (endPoint->GetPhysicalVolume()->GetName() == "pHemisphere") theTrack->SetTrackStatus(fStopAndKill); //artificially kills particles that go back to the world
  
  // You may uncomment next line for debugging
  // std::cout<<txtyellow<<"I am a "<<particleName<<" that now is in "<<currentPhysicalName<<" towards "<<endPoint->GetPhysicalVolume()->GetName()<<txtreset<<std::endl;
    
  const G4String nextPhysicalName = endPoint->GetPhysicalVolume()->GetName();
  // for debugging
  /*
  if (particleName=="pe-" || particleName=="e-" ){
	  std::cout<<txtyellow<<"Track "<<track_id<<"; Step: "<<theTrack->GetCurrentStepNumber()<<": I am a "<<particleName<<" that now is in "<<currentPhysicalName<<" towards "<<endPoint->GetPhysicalVolume()->GetName();
	  std::cout<<". Track length = "<<trackLength/um<<" um"<<"; kinE ="<<kinE/eV<<" eV."<<txtreset<<std::endl;
	  if (theTrack->GetTrackStatus()!= fAlive) std::cout<<txtyellow<<particleName<<": Track "<<track_id<<" died after running for "<<trackLength/um<<" um"<<txtreset<<std::endl;
  }
  */	
	
     if (currentPhysicalName == "World" && nextPhysicalName == "pSpaceCraft" && parent_id == 0) {
        fRunAct->SetKinEnergy(kinE);
    }
    
  // Look at new generated charged particles
  if (particle->GetPDGCharge()!=0 && theTrack->GetCurrentStepNumber() == 1 ){
	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  
	  if (particleName=="pe-" || particleName=="e-"){
		  
          G4int exitFlag = 0;
          
         // if (write == 1) {
              analysisManager->FillH1(1,log10f(kinE/eV));
		  
              // G4cout << "Me escribo :D" << G4endl;
              if (currentPhysicalName == "pHousing") analysisManager->FillH1(5,log10f(kinE/eV));
              else if (currentPhysicalName == "pHousingShell") analysisManager->FillH1(6,log10f(kinE/eV));
              else if (currentPhysicalName == "pTestMass") analysisManager->FillH1(7,log10f(kinE/eV));
              else if (currentPhysicalName == "pTestMassShell") analysisManager->FillH1(8,log10f(kinE/eV));
              
              if (kinE/keV >= 1) analysisManager->FillH2(1,theTrack->GetPosition().x()/mm,theTrack->GetPosition().z()/mm);
              else analysisManager->FillH2(0,theTrack->GetPosition().x()/mm,theTrack->GetPosition().z()/mm);
       //   }
          if (currentPhysicalName == "pTestMassShell" && nextPhysicalName == "pVacGap1") exitFlag = 1;
          if (currentPhysicalName == "pVacGap1" && nextPhysicalName == "pTesMassShell" && exitFlag == 1) {
              G4RunManager::GetRunManager()->rndmSaveThisEvent();
              exitFlag = 0;
          }
	  }
	  
	  if (particleName=="pe+" || particleName=="e+"){
		  		
        //  if (write == 1) {
              // G4cout << "Me escribo :D" << G4endl;
		      analysisManager->FillH1(2,log10f(kinE/eV));
        
              if (currentPhysicalName == "pHousing") analysisManager->FillH1(5,log10f(kinE/eV));
              else if (currentPhysicalName == "pHousingShell") analysisManager->FillH1(6,log10f(kinE/eV));
              else if (currentPhysicalName == "pTestMass") analysisManager->FillH1(7,log10f(kinE/eV));
              else if (currentPhysicalName == "pTestMassShell") analysisManager->FillH1(8,log10f(kinE/eV));
              
              if (kinE/keV >= 1) analysisManager->FillH2(1,theTrack->GetPosition().x()/mm,theTrack->GetPosition().z()/mm);
              else analysisManager->FillH2(0,theTrack->GetPosition().x()/mm,theTrack->GetPosition().z()/mm);
       //   }
	  }
	  
	  fRunAct->SetBirthE(kinE);
	  fRunAct->SetBirthDist(theTrack->GetPosition().mag());
	  fRunAct->SetBirthPos(theTrack->GetPosition());
  }
  
  if (endPoint->GetPhysicalVolume()->GetName() != "World"){
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      
	  // Look for particles entering the TM
	  if (currentPhysicalName == "pVacGap1" && endPoint->GetPhysicalVolume()->GetName() == CriticalVolume){
		  
		  if (particle->GetPDGCharge() != 0){
			  fRunAct->UpdateChargeCounter(particle->GetPDGCharge());
              
              EnteringCharge = particle->GetPDGCharge();
              SQCharge = EnteringCharge*EnteringCharge;
              RootSQCharge = sqrt(SQCharge);
              
              fRunAct->UpdateCharge(EnteringCharge);
              fRunAct->UpdateSQCharge(SQCharge);
              fRunAct->UpdateRootSQCharge(SQCharge);
              
              fRunAct->UpdateEnteringCharge(EnteringCharge);
              
              // analysisManager->FillH2(4,kinE,EnteringCharge);
			  
			  if (particleName=="pe-" || particleName=="e-" || particleName=="pe+" || particleName=="e+"){
				  analysisManager->FillH1(4,log10f(kinE/eV)); //Also includes positrons since 2022.08.24
				  
				  //analysisManager->FillH2(1,log10f(fRunAct->GetBirthE()/eV),log10f(trackLength/m)); //Removed since 2022.08.24
				  analysisManager->FillH2(3,fRunAct->GetBirthPos().x()/mm,fRunAct->GetBirthPos().z()/mm);
				  
			  }
		  }
		  
		  /*//for debugging only
		  if (particleName=="pe-" || particleName=="e-") {
			  std::cout<<txtblue<<"e- entering the TM\n"<<txtreset;
			  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
			  std::cout<<"Track id = "<<track_id<<std::endl;
			  std::cout<<"Track length = "<<trackLength/um<<" um"<<std::endl;
			  std::cout<<"Kinetic Energy = "<<kinE/eV<<" eV"<<std::endl;
		  }
		  if (particleName=="pe+" || particleName=="e+"){
			  std::cout<<txtyellow<<"e+ entering the TM\n"<<txtreset;
			  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
			  std::cout<<"Track id = "<<track_id<<std::endl;
			  std::cout<<"Track length = "<<trackLength/um<<" um"<<std::endl;
			  std::cout<<"Kinetic Energy = "<<kinE/eV<<" eV"<<std::endl;
		  }
		  if (particleName=="proton") {
			  std::cout<<txtred<<"proton entering the TM\n"<<txtreset;
			  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
			  std::cout<<"Track id = "<<track_id<<std::endl;
			  std::cout<<"Track length = "<<trackLength/cm<<" cm"<<std::endl;
		  }*/
	  }
	  
	  // particle leaving the TM
	  if (currentPhysicalName == CriticalVolume && endPoint->GetPhysicalVolume()->GetName() == "pVacGap1"){
          G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
          
		  if (particle->GetPDGCharge() != 0){
			  fRunAct->UpdateChargeCounter(-1*particle->GetPDGCharge());
              
              ExitCharge = -1*particle->GetPDGCharge();
              SQCharge = ExitCharge*ExitCharge;
              RootSQCharge = sqrt(SQCharge);
              
              fRunAct->UpdateCharge(ExitCharge);
              fRunAct->UpdateSQCharge(SQCharge);
              fRunAct->UpdateRootSQCharge(RootSQCharge);
              
              fRunAct->UpdateExitingCharge(ExitCharge);
              
              
              // analysisManager->FillH2(5,kinE,ExitCharge);
              
			  if (particleName=="pe-" || particleName=="e-" || particleName=="pe+" || particleName=="e+"){
                  
				  analysisManager->FillH1(3,log10f(kinE/eV)); //Also includes positrons since 2022.08.24
				  
				  //analysisManager->FillH2(0,log10f(fRunAct->GetBirthE()/eV),log10f(trackLength/m)); //Removed since 2022.08.24
				  analysisManager->FillH2(2,fRunAct->GetBirthPos().x()/mm,fRunAct->GetBirthPos().z()/mm);
				  
				  
			  }
		  }
		  /*
		  //for debugging
		  if (particle->GetPDGCharge() != 0 && track_id > 1){
			  if (particleName=="pe-" || particleName=="e-") {
				  std::cout<<txtblue<<"e- leaving the TM\n"<<txtreset;
				  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
				  std::cout<<"Track id = "<<track_id<<std::endl;
				  std::cout<<"Track length = "<<trackLength/um<<" um"<<std::endl;
				  std::cout<<"Kinetic Energy = "<<kinE/eV<<" eV"<<std::endl;
				  std::cout<<"Next volume: "<<endPoint->GetPhysicalVolume()->GetName()<<std::endl;
			  }
			  else if (particleName=="pe+" || particleName=="e+") {
				  std::cout<<txtyellow<<"e+ leaving the TM\n"<<txtreset;
				  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
				  std::cout<<"Track id = "<<track_id<<std::endl;
				  std::cout<<"Track length = "<<trackLength/um<<" um"<<std::endl;
				  std::cout<<"Kinetic Energy = "<<kinE/eV<<" eV"<<std::endl;
				  std::cout<<"Next volume: "<<endPoint->GetPhysicalVolume()->GetName()<<std::endl;
			  }
			  /*else if (particleName=="proton") {
				  std::cout<<txtred<<"proton leaving the TM\n"<<txtreset;
				  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
				  std::cout<<"Track id = "<<track_id<<std::endl;
				  std::cout<<"Track length = "<<trackLength/cm<<" cm"<<std::endl;
				  std::cout<<"Kinetic Energy = "<<kinE/eV<<" eV"<<std::endl;
				  std::cout<<"Next volume: "<<endPoint->GetPhysicalVolume()->GetName()<<std::endl;
			  }
			  else{
				  std::cout<<txtgreen<<particleName<<" leaving the TM\n"<<txtreset;
				  std::cout<<"Particle charge = "<<particle->GetPDGCharge()<<std::endl;
				  std::cout<<"Track id = "<<track_id<<std::endl;
				  std::cout<<"Track length = "<<trackLength/cm<<" cm"<<std::endl;
				  std::cout<<"Kinetic Energy = "<<trackLength/eV<<" eV"<<std::endl;
				  std::cout<<"Next volume: "<<endPoint->GetPhysicalVolume()->GetName()<<std::endl;
			  }
		  }*/
	  }
  }

  // Apply conditions to kill the track if needed
  // This must be applied at the end in order to score the last step
  //  KillTrackCondition(aTrack);
}


//-----------------------------------------------------------------------------
void SteppingAction::KillTrackCondition(G4Track* aTrack)
{
  
  // Kill tracks when they escape from the disc
  const G4Step* aStep = aTrack->GetStep();
  G4int preDepth = aStep->GetPreStepPoint()
    ->GetTouchableHandle()->GetHistoryDepth();
  G4int postDepth = aStep->GetPostStepPoint()
    ->GetTouchableHandle()->GetHistoryDepth();
  // G4cout << "\tpreDepth = " << preDepth;
  // G4cout << "\tpostDepth = " << postDepth << G4endl;

  //following lines should be commented?
    if (postDepth == 0 && preDepth > 0) {
    // Particle coming back to the enclosure, kill it
    aTrack->SetTrackStatus(fStopAndKill);
    }
}
