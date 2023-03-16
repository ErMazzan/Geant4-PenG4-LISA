//-----------------------------------------------------------------------------
//
// PenPositron.cc
//
//  version 2018.06.05
//
//-----------------------------------------------------------------------------

#include "PenPositron.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
    
PenPositron* PenPositron::theInstance = 0;

PenPositron* PenPositron::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "pe+";
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
                 name,  0.51099906*MeV,       0.0*MeV,    +1.*eplus, 
		    1,               0,             0,          
		    0,               0,             0,             
	     "lepton",              -1,             0,          -11,
		 true,            -1.0,          NULL,
                false,             "e"
              );
 
    // Bohr Magnetron
   G4double muB =  0.5*eplus*CLHEP::hbar_Planck/(0.51099906*MeV/CLHEP::c_squared) ;
   
   anInstance->SetPDGMagneticMoment( muB * 1.0011596521859 );

  }
  theInstance = reinterpret_cast<PenPositron*>(anInstance);
  return theInstance;
}

PenPositron*  PenPositron::PositronDefinition()
{
  return Definition();
}

PenPositron*  PenPositron::Positron()
{
  return Definition();
}


