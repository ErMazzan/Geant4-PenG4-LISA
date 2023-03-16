//MACG 2019-04-25: Multi-purpose PenPrimaryGeneratorAction class.


#ifndef PenPrimaryGeneratorAction_h
#define PenPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Event;
// class G4ParticleGun;
class G4GeneralParticleSource;

class PenPrimaryGeneratorMessenger;


class PenPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PenPrimaryGeneratorAction();    
  ~PenPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  
private:
  void SetDefaultParameters();
  void SetDefaultPrimaryParticle();

public:

  // =================
  // MODIFIER METHODS
  // =================

  // Method to set the type of beam used
  void SetSource(G4int source);

  // Methods to change the parameters of the default radiation source
  // interactively
  void SetParticle(G4String);
  void SetIon(G4int, G4int, G4double );  // only nucleus properties

  inline void SetMeanKinEne(G4double val)  { fMeanEne = val; }
  inline void SetSigmaKinEne(G4double val) { fSigmaEne = val; }

  inline void SetPosition(G4ThreeVector pos) {
    fX0 = pos.x();
    fY0 = pos.y();
    fZ0 = pos.z();
  }
  inline void SetSigmaX(G4double val) { fSigmaX = val; }
  inline void SetSigmaY(G4double val) { fSigmaY = val; }
  inline void SetCone(G4ThreeVector coneDir) {
    fTheta = coneDir.x();
    fPhi = coneDir.y();
    fAlpha = coneDir.z();
    fUx = std::sin(fTheta)*std::cos(fPhi);
    fUy = std::sin(fTheta)*std::sin(fPhi);
    fUz = std::cos(fTheta);
  }
  

  // ===============
  // ACCESS METHODS
  // ===============

  inline G4String GetSourceID() const     {return fSourceID;}


private:

  // 3-digit integer to choose the radiation source. The values are:
  // - 1st digit (hundreds): Energy Distribution.
  //    0 = mono or gaussian.
  // - 2nd digit (tens): Position Distribution.
  //    0 = point or 2D Gaussian distribution.
  // - 3rd digit (units): Angular Distribution.
  //    0 = conical beam around (theta,phi) direction.
  // Examples:
  // - Default value = 000;


  G4int fSourceID;
  G4int fEneType;   // These are obtained from fSourceID
  G4int fPosType;   // automatically.
  G4int fAngType;   //

  // Particle Gun
  // G4ParticleGun* fParticleGun;
  G4GeneralParticleSource* fGPSGun;


  // Data members for each kind of distribution
  //
  // Energy
  G4double fMeanEne;
  G4double fSigmaEne;
  // Position
  G4double fX0;
  G4double fY0;
  G4double fZ0;
  G4double fSigmaX;
  G4double fSigmaY;
  // Direction
  G4double fTheta;
  G4double fPhi;
  G4double fAlpha;
  G4double fUx;
  G4double fUy;
  G4double fUz;

  // Messenger class
  PenPrimaryGeneratorMessenger* fMessenger;

};

#endif
