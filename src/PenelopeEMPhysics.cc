
#include "PenelopeEMPhysics.hh"

#include "PenelopeEMPhysicsMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(PenelopeEMPhysics);

PenelopeEMPhysics::PenelopeEMPhysics(const G4String& name)
: G4VPhysicsConstructor(name)
{
  // Probar con 500 keV y con 50 eV
  fEthr = 500.*MeV; //*keV;
  fEpmax = 1.0*GeV;    // This is maximum value in penelope2018

  fMessenger = new PenelopeEMPhysicsMessenger(this);
}

PenelopeEMPhysics::~PenelopeEMPhysics()
{
  delete fMessenger;
}

#include "PenGamma.hh"
#include "PenElectron.hh"
#include "PenPositron.hh"

void PenelopeEMPhysics::ConstructParticle()
{
  PenGamma::GammaDefinition();
  PenElectron::ElectronDefinition();
  PenPositron::PositronDefinition();
}

#include "PenInterface.hh"
#include "PenPartConvertProcess.hh"
#include "PenEMProcess.hh"

void PenelopeEMPhysics::ConstructProcess()
{
  static G4bool first = true;
  if(first) {
    // Maximum kinetic energy in eV to be used in Penelope
    PenInterface* penIF = PenInterface::GetInstance();
    penIF->SetEPMAX(fEpmax/eV);
    // Make sure this limitation is respected in G4
    if(fEthr > fEpmax) {
      G4Exception("PenelopeEMPhysics::ConstructProcess","PenEM0001",
                  JustWarning,"fEthr is pulled down to fEpmax.");
      fEthr = 0.999*fEpmax;
    }
    first = false;
    G4cout << "++++ PenelopeEMPhysics::ConstructProcess() ++++" << G4endl;
    G4cout << "  EPMAX[eV]= " << fEpmax/eV << G4endl;
    G4cout << "  Ethr[eV]= " << fEthr/eV << G4endl;
  }

  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma" ||
        particleName == "e-" ||
        particleName == "e+") {
      pmanager->AddDiscreteProcess(new PenPartConvertProcess(fEthr));
    }

    else if( particleName == "pgamma" ||
             particleName == "pe-" ||
             particleName == "pe+" ) {
      // gamma, electron, positron for Penelope
      pmanager->AddDiscreteProcess(new PenEMProcess);
      // pmanager->AddProcess(new PenEMProcess, -1, 1, 1);
    }
  }
}
