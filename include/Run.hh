
#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"

#include "G4THitsMap.hh"
#include <vector>

class PenOutputPrinter;


class Run : public G4Run
{
public:
  // constructor and destructor.
  Run(PenOutputPrinter* );
  virtual ~Run();

private:
  Run();

public:
  // virtual method from G4Run. 
  // The method is overriden in this class for scoring.
  virtual void RecordEvent(const G4Event*);

  // virtual method form G4Run.
  // To merge local G4Run object into the global G4Run object.
  virtual void Merge(const G4Run*);

public:  // Get/Set methods
  inline std::vector<G4double>* GetGridRZDoseRun() const
  { return fGridRZDoseRun; }
  inline std::vector<G4double>* GetGridRZDoseRun2() const
  { return fGridRZDoseRun2; }
  inline std::vector<G4double>* GetDepthDoseRun() const
  { return fDepthDoseRun; }
  inline std::vector<G4double>* GetDepthDoseRun2() const
  { return fDepthDoseRun2; }
  
  inline std::vector<G4double>* GetSpcEnddetRun(size_t ii) const
  { return (*fSpcEnddetRun)[ii]; }
  inline std::vector<G4double>* GetSpcEnddetRun2(size_t ii) const
  { return (*fSpcEnddetRun2)[ii]; }
  inline std::vector<G4double>* GetPolarAngleElectRun() const
  { return fPolarAngleElectRun; }
  inline std::vector<G4double>* GetPolarAngleElectRun2() const
  { return fPolarAngleElectRun2; }
  inline std::vector<G4double>* GetPolarAngleGammaRun() const
  { return fPolarAngleGammaRun; }
  inline std::vector<G4double>* GetPolarAngleGammaRun2() const
  { return fPolarAngleGammaRun2; }
  inline std::vector<G4double>* GetPolarAnglePositRun() const
  { return fPolarAnglePositRun; }
  inline std::vector<G4double>* GetPolarAnglePositRun2() const
  { return fPolarAnglePositRun2; }
  inline std::vector<G4double>* GetEmergEneUpElectRun() const
  { return fEmergEneUpElectRun; }
  inline std::vector<G4double>* GetEmergEneUpElectRun2() const
  { return fEmergEneUpElectRun2; }
  inline std::vector<G4double>* GetEmergEneUpGammaRun() const
  { return fEmergEneUpGammaRun; }
  inline std::vector<G4double>* GetEmergEneUpGammaRun2() const
  { return fEmergEneUpGammaRun2; }
  inline std::vector<G4double>* GetEmergEneUpPositRun() const
  { return fEmergEneUpPositRun; }
  inline std::vector<G4double>* GetEmergEneUpPositRun2() const
  { return fEmergEneUpPositRun2; }
  inline std::vector<G4double>* GetEmergEneDownElectRun() const
  { return fEmergEneDownElectRun; }
  inline std::vector<G4double>* GetEmergEneDownElectRun2() const
  { return fEmergEneDownElectRun2; }
  inline std::vector<G4double>* GetEmergEneDownGammaRun() const
  { return fEmergEneDownGammaRun; }
  inline std::vector<G4double>* GetEmergEneDownGammaRun2() const
  { return fEmergEneDownGammaRun2; }
  inline std::vector<G4double>* GetEmergEneDownPositRun() const
  { return fEmergEneDownPositRun; }
  inline std::vector<G4double>* GetEmergEneDownPositRun2() const
  { return fEmergEneDownPositRun2; }


private:
  
  PenOutputPrinter* fOutputPrinter;
  
  std::vector<G4double>* fGridRZDoseRun;   // Edep/rho for 2d-dose-map
  std::vector<G4double>* fGridRZDoseRun2;  // Sum(Edep/rho)^2 for 2d-dose-map
  std::vector<G4double>* fDepthDoseRun;   // Edep/rho for depth-dose
  std::vector<G4double>* fDepthDoseRun2;  // Sum(Edep/rho)^2 for depth-dose

  // Edep spectrum in volume
  std::vector< std::vector<G4double>* >* fSpcEnddetRun;
  // Sum of w^2 Edep spectrum in volume
  std::vector< std::vector<G4double>* >* fSpcEnddetRun2;
  // BinDE values, to avoid its calculation at RecordEvent
  std::vector<G4double>* fSpcEnddetBinDEVec;

  // polar-angle.dat
  std::vector<G4double>* fPolarAngleElectRun;
  std::vector<G4double>* fPolarAngleElectRun2;
  std::vector<G4double>* fPolarAngleGammaRun;
  std::vector<G4double>* fPolarAngleGammaRun2;
  std::vector<G4double>* fPolarAnglePositRun;
  std::vector<G4double>* fPolarAnglePositRun2;
  // energy-up.dat
  std::vector<G4double>* fEmergEneUpElectRun;
  std::vector<G4double>* fEmergEneUpElectRun2;
  std::vector<G4double>* fEmergEneUpGammaRun;
  std::vector<G4double>* fEmergEneUpGammaRun2;
  std::vector<G4double>* fEmergEneUpPositRun;
  std::vector<G4double>* fEmergEneUpPositRun2;
  // energy-down.dat
  std::vector<G4double>* fEmergEneDownElectRun;
  std::vector<G4double>* fEmergEneDownElectRun2;
  std::vector<G4double>* fEmergEneDownGammaRun;
  std::vector<G4double>* fEmergEneDownGammaRun2;
  std::vector<G4double>* fEmergEneDownPositRun;
  std::vector<G4double>* fEmergEneDownPositRun2;

};

#endif
