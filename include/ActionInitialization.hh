
#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization(DetectorConstruction* );
  virtual ~ActionInitialization();

public:
  virtual void BuildForMaster() const;
  virtual void Build() const;


public:
  DetectorConstruction* fDetector;

};


#endif


