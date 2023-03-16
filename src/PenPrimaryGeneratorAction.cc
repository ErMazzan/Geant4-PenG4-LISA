//MACG 2019-04-25: Multi-purpose PenPrimaryGeneratorAction class.


#include "PenPrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Event.hh"
// #include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include "PenPrimaryGeneratorMessenger.hh"

#include <vector>


//-------1---------2---------3---------4---------5---------6---------7---------8
//-----------------------------------------------------------------------------
PenPrimaryGeneratorAction::PenPrimaryGeneratorAction()
{
  // Define the messenger
  fMessenger = new PenPrimaryGeneratorMessenger(this);

  // Define particle generator
  // G4int n_particle = 1;
  // fParticleGun  = new G4GeneralParticleSource();
  fGPSGun = new G4GeneralParticleSource();

  // Default beam parameters
  SetDefaultParameters();
  SetDefaultPrimaryParticle();
}  


//-----------------------------------------------------------------------------
PenPrimaryGeneratorAction::~PenPrimaryGeneratorAction()
{
  delete fGPSGun;
  // delete fMessenger;
}


//-----------------------------------------------------------------------------
void PenPrimaryGeneratorAction::SetDefaultParameters()
{

  // Default source ID
  fSourceID = 0;
  fEneType = 0;
  fPosType = 0;
  fAngType = 0;

  // Default energy distribution: 10 keV monoenergetic
  fMeanEne = 40.*keV;
  fSigmaEne = 0.*keV;

  // Default beam shape: Point source at origin of coordinates
  // or Gaussian over surface perpendicular to +Z axis.
  fX0 = 0.0;
  fY0 = 0.0;
  fZ0 = 0.0;
  fSigmaX = 0.0;
  fSigmaY = 0.0;

  // Default momentum direction distribution: only +Z axis.
  fTheta = 0.*deg;
  fPhi =   0.*deg;
  fAlpha = 5.*deg;  // aperture of uniform conical beam
  // Main component auto-calculation to avoid sin & cos calculations per event. 
  fUx = std::sin(fTheta)*std::cos(fPhi);
  fUy = std::sin(fTheta)*std::sin(fPhi);
  fUz = std::cos(fTheta);

}
  

//-----------------------------------------------------------------------------
void PenPrimaryGeneratorAction::SetDefaultPrimaryParticle()
{
  // Define primary particle: pe-
  //-----------------------------------
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("pe-");


  // General Particle Source default setup:
  //  Monoenergetic and plane wave (source ID 000)
  //-----------------------------------------------
    fGPSGun->SetParticleDefinition(particle);
}


//-----------------------------------------------------------------------------
void PenPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // ===================
  // TIME DISTRIBUTION
  // ===================
    
    // fGPSGun->GeneratePrimaryVertex(anEvent);
    G4double t0 = 0.0;
    G4double tf = 1.0;
    G4double time = (tf+t0)*G4UniformRand();
    fGPSGun->SetParticleTime(time*s);
    
  // ====================
  // ENERGY DISTRIBUTION
  // ====================

  if (fEneType == 0) {
    // 
    // Default energy distribution:
    // Set the beam energy and Gaussian energy spread
    // -----------------------------------------------

    G4double kineticEnergy = fMeanEne;

    if (fSigmaEne > 0.0)  {
      G4double randEne = 0.0;
      // This condition avoids too large deviations
      do {
	randEne = G4RandGauss::shoot();
      } while ( randEne < -5. || randEne > +5. );

      randEne *= fSigmaEne;
      kineticEnergy += randEne;
    }
      
    // fParticleGun->SetParticleEnergy( kineticEnergy );
  }

  // ======================
  // POSITION DISTRIBUTION
  // ======================

  if ( fPosType == 0 ) {
    //
    // Default position distribution:
    // Set the beam spot size (Gaussian shape)
    // ----------------------------------------
      
    // Set the position of the primary particles
    G4double x = fX0;
    G4double y = fY0;
    G4double z = fZ0;
      
    if ( fSigmaX > 0.0 ) {
      G4double randX = 0.0;
      // This condition avoids too large deviations
      do {
	randX = G4RandGauss::shoot();
      } while ( randX < -5. || randX > +5. );
	
      randX *= fSigmaX;
      x += randX;
    }

    if ( fSigmaY > 0.0 ) {
      G4double randY = 0.0;
      // This condition avoids too large deviations
      do {
	randY = G4RandGauss::shoot();
      } while ( randY < -5. || randY > +5. );

      randY *= fSigmaY;
      y += randY;
    }

      fGPSGun->SetParticlePosition(G4ThreeVector(x, y, z));
  }


  // =====================
  // ANGULAR DISTRIBUTION
  // =====================

  if (fAngType == 0) {
    //
    // Default angular distribution:
    // Conical beam
    // ------------------------------

    // Beam direction components
    G4ThreeVector dir(fUx, fUy, fUz);
    // G4cout << "dir = "<< dir << G4endl;
    
    if (fAlpha > 0.) {
      // Calculate random change of vector as if the original one were +Uz.
      G4double cosThetaMin = std::cos(fAlpha);
      G4double rndPhi = G4UniformRand()*2.*M_PI*rad;
      G4double rndUz = 1.-G4UniformRand()*(1.-cosThetaMin);  // cos(theta)
      G4ThreeVector rndDelta(std::sqrt(1.-rndUz*rndUz)*std::cos(rndPhi),
			     std::sqrt(1.-rndUz*rndUz)*std::sin(rndPhi),
			     rndUz-1.);
      // G4cout << "rndDelta = " << rndDelta << G4endl;

      // Rotate rndDelta according to tranformation +Uz -> beam direction
      rndDelta.rotateUz(dir);
      // G4cout << "rndDelta = " << rndDelta << G4endl;
      
      // Sum change to beam direction
      dir += rndDelta;
    }
    // G4cout << "dir = " << dir
    // 	   << "\tdir.mag2() = " << dir.mag2() << G4endl;

    // fParticleGun->SetParticleMomentumDirectionZ(dir);
  }

  // ===================
  // GENERATE PRIMARIES
  // ===================

    fGPSGun->GeneratePrimaryVertex( anEvent );
  
}



//-----------------------------------------------------------------------------
void PenPrimaryGeneratorAction::SetParticle(G4String name)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(name);
    fGPSGun->SetParticleDefinition(particle);
}


//-----------------------------------------------------------------------------
void PenPrimaryGeneratorAction::SetIon(G4int Z, G4int A, G4double E = 0.0)
{
  G4ParticleDefinition* ion = 0;
  G4IonTable* ionTable = G4IonTable::GetIonTable();
  ion = ionTable->GetIon(Z, A, E);

  if (ion == 0)  {
    G4cout << "Ion with Z = " << Z << ", A = " << "is not defined!"; 
    G4cout << G4endl;
    G4Exception("PenPrimaryGeneratorAction::SetIon()",
		"BadIonDef", FatalException,
		"Unknown definition of ion.");
  }
  else  {
      fGPSGun->SetParticleDefinition(ion);
  }
}



//-----------------------------------------------------------------------------
void PenPrimaryGeneratorAction::SetSource(G4int source)
{
  fSourceID = source;

  // Get the units, tens and hundreds of fSourceID
  // and store into the corresponding data members
  //-------------------------------------------------
  fAngType = (source % 10);
  fPosType = ((source/10) % 10);
  fEneType = source/100;          // 3-digit integer
}
