
#include "PenelopeEMPhysics.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "Hadr06PhysicsList.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"
// #include "G4UIterminal.hh"
// #include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4VModularPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"


int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Dani: setting the seeds here is essential for sending multiple jobs!
  G4long seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);
  G4Random::setTheSeed(seed);
	
 // Construct the run manager
 //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  // runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
  // runManager->SetNumberOfThreads(1);
#else
  G4RunManager* runManager = new G4RunManager;
#endif
  G4ScoringManager::GetScoringManager();

  // Set detector construction
  //
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // Set physics list
  //
  /*
  G4VModularPhysicsList* physics = new FTFP_BERT;
  physics->RegisterPhysics(new G4StepLimiterPhysics);
  //physics->RegisterPhysics(new PenelopeEMPhysics);
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(new Hadr06PhysicsList());
  */
  G4VModularPhysicsList* physics = new Hadr06PhysicsList();
  //physics->RegisterPhysics(new G4StepLimiterPhysics);
  runManager->SetUserInitialization(physics);

  

  // Set user action classes
  //
  G4VUserActionInitialization* user_action
    = new ActionInitialization(detector);
  runManager->SetUserInitialization(user_action);
  //

  // Initialize G4 kernel
  //
  // runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}
