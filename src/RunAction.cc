    
#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "g4root.hh"

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PenOutputPrinter.hh"


//-------1---------2---------3---------4---------5---------6---------7---------8
RunAction::RunAction(const DetectorConstruction* detConst)
{
  G4cout << "RunAction object created" << G4endl;
  // Create the printer
  
  ////////////////////////////////////////////
  // 1D histograms 
  ////////////////////////////////////////////
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // 0. Deposited charge in TM per event
  analysisManager->CreateH1("hQdep","Deposited Charge in TM per event", 1001, -500.5, 500.5); // -> Para comprobar el loop
  
  // 1. log(Energy) of the created electrons
  analysisManager->CreateH1("eleKinE","Kinetic energy of new secondary electrons", 40, 1, 9); //
  
  // 2. log(Energy) of the created positrons
  analysisManager->CreateH1("posKinE","Kinetic energy of new secondary positrons", 40, 1, 9); //
  
  // 3. log(Energy) of the escaping electrons/positrons
  analysisManager->CreateH1("EscKinE","Kinetic energy of the escaping electrons/positrons", 40, 1, 9);
  
  // 4. log(Energy) of the entering electrons/positrons
  analysisManager->CreateH1("EntKinE","Kinetic energy of the entering electrons/positrons", 40, 1, 9);
  
  // 5. log(Energy) of the created electrons/positrons in the electron housing
  analysisManager->CreateH1("HouAlleleKinE","Kinetic energy of new secondary electrons", 40, 1, 9); //
  
  // 6. log(Energy) of the created electrons/positrons in the electron housing inner shell
  analysisManager->CreateH1("HouShAlleleKinE","Kinetic energy of new secondary positrons", 40, 1, 9); //
  
  // 7. log(Energy) of the created electrons/positrons in the TM
  analysisManager->CreateH1("TMAlleleKinE","Kinetic energy of new secondary electrons", 40, 1, 9); //
  
  // 8. log(Energy) of the created electrons/positrons in the TM shell
  analysisManager->CreateH1("TMShAlleleKinE","Kinetic energy of new secondary positrons", 40, 1, 9); //
    
  // 9. Deposited charge in TM per Run
  analysisManager->CreateH1("hQdepPerRun", "Deposited Charge in TM per run", 1001, -500.5, 500.5);
      
  // 10. Square charge for each entering/exiting particle
  analysisManager->CreateH1("SquareChargeCounter", "Square charge for each entering/exiting particle", 1000, 0.0, 5000);

  // 11. Event rate for each net event charge
  analysisManager->CreateH1("EventRate", "Event rate for each net event charge", 80, -40, 40);
    
  // 12. Charging rate for energy value
  analysisManager->CreateH1("ChargingRate", "Charging rate for energy value", 50, 10, 100000, "MeV", "none", "log");
    
  ////////////////////////////////////////////
  // 2D histograms 
  ////////////////////////////////////////////
  
  // Removed: log(Initial Energy) vs range of the escaping electrons/positrons
  //analysisManager->CreateH2("rangeEsc","Range vs E for escaping e-/e+", 40, 1, 9, 18, -10, -1);
  
  // Removed: log(Initial Energy) vs range of the entering electrons/positrons
  //analysisManager->CreateH2("rangeEnt","Range vs E for entering e-/e+", 40, 1, 9, 18, -10, -1);
  
  // Removed: log (Initial Energy) vs distance to the center
  //analysisManager->CreateH2("eleKinEdist","Distance to the center of the generated e-", 40, 1, 9, 100, 0, 100);

  // 0. x,z coordinates of the generated electrons (e<1 keV)
  analysisManager->CreateH2("distrLE","x,y coordinates of the generated e-", 201, -45.5, 45.5, 201, -45.5, 45.5); //

  // 1. x,z coordinates of the generated electrons (e>1 keV)
  analysisManager->CreateH2("distrHE","x,y coordinates of the generated e-", 201, -45.5, 45.5, 201, -45.5, 45.5); //
  
  // 2. x,z coordinates  electrons/positrons escaping the TM
  analysisManager->CreateH2("EscDistr","x,z coordinates of the generated e-/e+", 201, -45.5, 45.5, 201, -45.5, 45.5);
  
  // 3. x,z coordinates for electrons/positrons entering the TM
  analysisManager->CreateH2("EntDistr","x,z coordinates of the entering e-/e+", 201, -45.5, 45.5, 201, -45.5, 45.5);
    
  // 4. Distribution of entering charges of the TM as a function of the initial energy
  // analysisManager->CreateH2("EntCharge","net charge for entering particles", 100, 10, 1000000, 85, -70, 15, "MeV", "none", "none", "none", "log", "linear");
  
  // 5. Distribution of exiting charges of the TM as a function of the initial energy
  // analysisManager->CreateH2("ExtCharge","net charge for exiting particles", 100, 10, 1000000, 85, -70, 15, "MeV", "none", "none", "none", "log", "linear");
    
  // 4. Distribution of entering charges of the TM as a function of the initial energy per event
  analysisManager->CreateH2("EntChargePerEvent","net charge for entering particles", 100, 10, 1000000, 81, -60.5, 20.5, "MeV", "none", "none", "none", "log", "linear");
    
  // 5. Distribution of exiting charges of the TM as a function of the initial energy per event
  analysisManager->CreateH2("ExtChargePerEvent","net charge for exiting particles", 100, 10, 1000000, 81, -60.5, 20.5, "MeV", "none", "none", "none", "log", "linear");

  ////////////////////////////////////////////
  // 1D NTuples
  ////////////////////////////////////////////
    
  // Particles generated per unit of time (second) for an integrated flux
    
  /*
  analysisManager->CreateNtuple("LISA", "Particles per second");
    
  // 0. Particles generated for each run
  analysisManager->CreateNtupleDColumn("Particles");
    
  // 1. Time step for each run
  analysisManager->CreateNtupleDColumn("Time");
  
  analysisManager->FinishNtuple();
  */
  fOutputPrinter = new PenOutputPrinter(detConst);
}


//-------1---------2---------3---------4---------5---------6---------7---------8
RunAction::RunAction()
{
  G4cout << "void RunAction is not meant to be used!" << G4endl;
}


RunAction::~RunAction()
{
  // delete the printer
  if (fOutputPrinter) delete fOutputPrinter;
  G4cout << "RunAction object destroyed" << G4endl;
}


//-----------------------------------------------------------------------------
G4Run* RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated to store run-persistent data.
  return new Run(fOutputPrinter);
}


//-----------------------------------------------------------------------------
void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4cout << "RunAction::BeginOfRunAction() " << G4endl;
  
  fDepositedCharge = 0;
  fDepositedChargePerRun = 0;
  fSQCharge = 0;
    
  /*
  fSquareDepositedCharge = 0;
  fSquareRootDepositedCharge = 0;
  */
  
  fBirthE = 0.;
  fBirthDist = 0.;
  fRunID = 0;
    
  //fBirthPos = G4ThreeVector(0.,0.,0.);

  
  if (IsMaster()) {
    fOutputPrinter->CreateFiles();
    fOutputPrinter->BookForMaster();
  }
  else {
    fOutputPrinter->BookForWorker();
  }
  
  // Open an output file
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4String fileName = "Output_LISA.root";
    
  if (global_count == 0) analysisManager->CloseFile(fileName);
  analysisManager->OpenFile(fileName);
  
   
  G4cout << "Run number: " << global_count+1 << G4endl;
  
}


//-----------------------------------------------------------------------------
void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "RunAction::EndOfRunAction() " << G4endl;
  
  // global_count allows us to store data between multiple runs
  global_count += 1;
    
  if (IsMaster()) {
    const Run* masterRun = static_cast<const Run*>(aRun);
    fOutputPrinter->DumpAll(masterRun);
    fOutputPrinter->CloseAllFiles();
  }
    
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  /*
  analysisManager->FillH1(9,fDepositedChargePerRun);
  analysisManager->FillH1(10,fSquareDepositedCharge);
  analysisManager->FillH1(11,fSquareRootDepositedCharge);
  */
    
  analysisManager->FillH1(9, fDepositedChargePerRun);
  analysisManager->FillH1(10, fSQCharge);
  // analysisManager->FillH1(11, fRootSQCharge);
    
  analysisManager->Write();
  /*
  G4cout << "Deposited Charge in TM per run: " << fDepositedChargePerRun << G4endl;
  G4cout << "Squared value of deposited Charge in TM per run: " << fSquareDepositedCharge << G4endl;
  G4cout << "Root value of squared deposited Charge in TM per run: " << fSquareRootDepositedCharge << G4endl;
  */
    
   static std::ofstream stuff("charge_data.csv");
   static bool first = true;
   if (first) {
     first = false;
     stuff << "Dep. Charge Run" << ";" << " Square Charge Run" << std::endl;
   }
   stuff << fDepositedChargePerRun << "; " << fSQCharge << "; " << std::endl;

    
  // Not closing the File allows us to overwrite data from multiple runs
  // It will provide just a warning message
    
  // analysisManager->CloseFile(fileName);
}
