//-----------------------------------------------------------------------------
//
// PenPartConvertProcess.cc
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#include "PenPartConvertProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "PenGamma.hh"
#include "PenElectron.hh"
#include "PenPositron.hh"
#include "PenInterface.hh"

PenPartConvertProcess::PenPartConvertProcess(G4double te,const G4String& processName)
:G4VDiscreteProcess(processName,fDecay),verboseLevel(0),threshE(te)
{
  pParticleChange = &fParticleChangeForDecay;
  pPenInterface = PenInterface::GetInstance();
}

PenPartConvertProcess::~PenPartConvertProcess()
{;}

PenPartConvertProcess::PenPartConvertProcess(const PenPartConvertProcess &right)
:G4VDiscreteProcess(right.GetProcessName(),fDecay),verboseLevel(1),threshE(right.threshE)
{;}

PenPartConvertProcess& PenPartConvertProcess::operator=(const PenPartConvertProcess & /*right*/)
{ return *this; }

G4bool PenPartConvertProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if(aParticleType.GetParticleName()=="gamma") return true;
  if(aParticleType.GetParticleName()=="e-") return true;
  if(aParticleType.GetParticleName()=="e+") return true;
  return false;
}

G4double PenPartConvertProcess::PostStepGetPhysicalInteractionLength(
    const G4Track& track,G4double /*previousStepSize*/,G4ForceCondition* condition)
{
  G4double remainder = DBL_MAX;
  if(track.GetKineticEnergy()<threshE &&
     pPenInterface->PenMatID(track.GetMaterial(),false)>0)
  { remainder = DBL_MIN; }
  *condition = NotForced;
  return remainder;
}

G4double PenPartConvertProcess::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*)
{
  return DBL_MAX;
}

void PenPartConvertProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
  return;
}

G4VParticleChange* PenPartConvertProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& )
{
  fParticleChangeForDecay.Initialize(aTrack);

  const G4ParticleDefinition* aType = aTrack.GetDefinition();
  G4ParticleDefinition* pd = 0;
  if(aType==G4Gamma::Gamma())
  { pd = PenGamma::Gamma(); }
  else if(aType==G4Electron::Electron())
  { pd = PenElectron::Electron(); }
  else if(aType==G4Positron::Positron())
  { pd = PenPositron::Positron(); }
  else
  { G4Exception("PenPartConvertProcess::PostStepDoIt","PenG4_0010",FatalException,
                "PenPartConvertProcess::PostStepDoIt - Unknown particle type."); }

  G4DynamicParticle* dp
    = new G4DynamicParticle(pd,aTrack.GetMomentumDirection(),aTrack.GetKineticEnergy());

  fParticleChangeForDecay.SetNumberOfSecondaries(1);
  G4Track* secondary = new G4Track(dp,aTrack.GetGlobalTime(),aTrack.GetPosition());
  secondary->SetGoodForTrackingFlag();
  secondary->SetTouchableHandle(aTrack.GetTouchableHandle());
  secondary->SetWeight(aTrack.GetWeight());
  fParticleChangeForDecay.AddSecondary(secondary);

  fParticleChangeForDecay.ProposeTrackStatus( fStopAndKill ) ;
  fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.);

  return &fParticleChangeForDecay ;
}

