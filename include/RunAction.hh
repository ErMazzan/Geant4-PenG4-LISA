
#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <fstream>

class DetectorConstruction;
class PenOutputPrinter;

class G4Run;


class RunAction : public G4UserRunAction
{
private:
  RunAction();

public:
  RunAction(const DetectorConstruction*);
  ~RunAction();

public:
  // virtual methods from G4UserRunAction

  G4Run* GenerateRun();
  // A derived G4Run is needed to score quantities during the run

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run* );

  // 'Get' methods
  PenOutputPrinter* GetOutputPrinter() const { return fOutputPrinter; }
  
  G4int fDepositedCharge;
  G4int GetDepositedCharge() const		{ return fDepositedCharge; }
  void UpdateChargeCounter(G4int count)		{ fDepositedCharge = fDepositedCharge+count; }
    
  G4int fCharge;
  G4int GetCharge() const       { return fCharge; }
  void UpdateCharge(G4int count)      { fCharge = fCharge+count; }
    
  G4int fSQCharge;
  void UpdateSQCharge(G4int count)        { fSQCharge = fSQCharge+count; }
    
  G4int fRootSQCharge;
  void UpdateRootSQCharge(G4int count)        { fRootSQCharge = fRootSQCharge+count; }
    
  G4int fDepositedChargePerRun;
  void UpdateChargeCounterPerRun(G4int count)       { fDepositedChargePerRun = fDepositedChargePerRun+count; }
     
  /*
  G4int fSquareDepositedCharge;
  void UpdateSquareCharge(G4int scount)       { fSquareDepositedCharge = fSquareDepositedCharge+scount; }
    
  G4int fSquareRootDepositedCharge;
  void UpdateSquareRootCharge(G4int srcount)       { fSquareRootDepositedCharge = fSquareRootDepositedCharge+srcount; }
  */
    
  G4double fKinEnergy;
  void SetKinEnergy(G4double Ene)        { fKinEnergy = Ene; }
  G4double GetKinEnergy()        { return fKinEnergy; }
    
  G4double fTotalEnteringCharge;
  void UpdateEnteringCharge(G4double charge)        { fTotalEnteringCharge += charge; }
  G4double GetEnteringCharge()        { return fTotalEnteringCharge; }
    
  G4double fTotalExitingCharge;
  void UpdateExitingCharge(G4double charge)        { fTotalExitingCharge += charge; }
  G4double GetExitingCharge()        { return fTotalExitingCharge; }
    
  G4double fBirthE;
  void SetBirthE(G4double Ene)		{ fBirthE = Ene; }
  G4double GetBirthE()		{ return fBirthE; }
    
  G4double fBirthDist;
  void SetBirthDist(G4double val)		{ fBirthDist = val; }
  G4double GetBirthDist()		{ return fBirthDist; }

  G4ThreeVector fBirthPos;
  void SetBirthPos(G4ThreeVector vec)		{ fBirthPos = vec; }
  G4ThreeVector GetBirthPos()		{ return fBirthPos; }
  
  void ResetChargeCounter() { fDepositedCharge = 0; }
    
  void ResetEnteringCharge() { fTotalEnteringCharge = 0; }
    
  void ResetExitingCharge() { fTotalExitingCharge = 0; }
    
  G4int fRunID;
  void SetRunID(G4int value) { fRunID = value; }
  G4int GetRunID() { return fRunID; }


private:

  // Data members
  PenOutputPrinter* fOutputPrinter;
    
  // File count
  G4int global_count = 0;

};

#endif
