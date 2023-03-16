//-----------------------------------------------------------------------------
//
// PenEMProcess.cc
//
//  version 2021.05.08
//
//-----------------------------------------------------------------------------

#include "PenEMProcess.hh"
#include "PenInterface.hh"
#include "PenPhys.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "PenGamma.hh"
#include "PenElectron.hh"
#include "PenPositron.hh"

#include "G4SystemOfUnits.hh"

#include "G4AutoLock.hh"
namespace {
 G4Mutex penelopePhysInitializerMutex = G4MUTEX_INITIALIZER;
}

PenEMProcess::PenEMProcess(const G4String& processName, G4ProcessType theType)
:G4VProcess(processName,theType),
 pWorthToTrack(false)
{
  pPenelope = PenInterface::GetInstance();
  pPenPhys = new PenPhys();
  pParticleChange = pPChange = new G4ParticleChange();
}

PenEMProcess::~PenEMProcess()
{
  delete pPenPhys;
  delete pPChange;
}

G4bool PenEMProcess::IsApplicable(const G4ParticleDefinition& p)
{ return (&p==PenGamma::Gamma()||&p==PenElectron::Electron()||&p==PenPositron::Positron()); }

void PenEMProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
  G4AutoLock l(&penelopePhysInitializerMutex);
  pPenelope->Initialize();
  pPenPhys->InitializeCEGRID();
}

void PenEMProcess::StartTracking(G4Track* aTrack)
{
  pPenPhys->NewTrack(aTrack);
}

G4double PenEMProcess::AlongStepGetPhysicalInteractionLength(
  const G4Track&, G4double, G4double, G4double&, G4GPILSelection*)
{
  return 0.;
  // return DBL_MAX;
}

G4double PenEMProcess::PostStepGetPhysicalInteractionLength(
  const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition)
{
  G4double proposingStepLength = DBL_MIN;
  pWorthToTrack = pPenPhys->WorthToTrack(&aTrack);
  if(pWorthToTrack)
  { proposingStepLength = pPenPhys->StepLength(&aTrack,previousStepSize); }
  *condition = NotForced;
  // G4cout << "PenEMProcess::PostStepGPIL" << "\tProposedStepLength[cm]= "
  // 	 << proposingStepLength/cm << G4endl;
  return proposingStepLength;
}

G4VParticleChange* PenEMProcess::AlongStepDoIt(const G4Track& /*aTrack*/,
					       const G4Step& /*aStep*/)
{
  return 0;
  // if(pWorthToTrack)
  // { pPenPhys->AlongStepDoIt(&aTrack,&aStep,pPChange); }
  // // G4cout << "\tEnergy[eV]= " << pPChange->GetEnergy()/eV << G4endl;
  // // G4cout << "\tAlongStepEdep[eV]= "
  // // 	 << pPChange->GetLocalEnergyDeposit()/eV << G4endl;
  // return pPChange;
}

G4VParticleChange* PenEMProcess::PostStepDoIt(const G4Track& aTrack,
					      const G4Step& /*aStep*/)
{
  if(pWorthToTrack)
  { pPenPhys->PostStepDoIt(&aTrack,pPChange); }
  else
  {
    pPChange->Initialize(aTrack);
    pPChange->ProposeTrackStatus(fStopAndKill);
    pPChange->ProposeEnergy(0.);
    pPChange->ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy());
  }
  // G4cout << "\tTrueStepLength[cm]= " << aStep.GetStepLength()/cm << G4endl;
  // G4cout << "\tEnergy [eV] = " << pPChange->GetEnergy()/eV << G4endl;
  // G4cout << "\tPostStepEdep[eV] = " << pPChange->GetLocalEnergyDeposit()/eV
  //        << G4endl;
  return pPChange;
}

G4double PenEMProcess::AtRestGetPhysicalInteractionLength(
                             const G4Track&, G4ForceCondition*)
{ return 0.; }

G4VParticleChange* PenEMProcess::AtRestDoIt(const G4Track&, const G4Step&)
{ return 0; }

