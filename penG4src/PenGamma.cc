//-----------------------------------------------------------------------------
//
// PenGamma.cc
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#include "PenGamma.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

PenGamma* PenGamma::theInstance = 0;

PenGamma*  PenGamma::Definition() 
{
  if (theInstance !=0) return theInstance;

  const G4String name = "pgamma";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
  {
  // create particle
  //      
  //    Arguments for constructor are as follows 
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 
  //             shortlived      subType    anti_encoding
   anInstance = new G4ParticleDefinition(
	         name,         0.0*MeV,       0.0*MeV,         0.0, 
		    2,              -1,            -1,          
		    0,               0,             0,             
	      "gamma",               0,             0,          22,
	      true,                0.0,          NULL,
             false,           "photon",          22
	      );
  }
  theInstance = reinterpret_cast<PenGamma*>(anInstance);
  return theInstance;
}

PenGamma*  PenGamma::GammaDefinition() 
{
  return Definition();
}

PenGamma*  PenGamma::Gamma() 
{
  return Definition();
}
