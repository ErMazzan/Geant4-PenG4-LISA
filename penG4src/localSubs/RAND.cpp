
//  *********************************************************************
//                         FUNCTION RAND (Random number generator)
//  *********************************************************************
#include "Randomize.hh"

double PenPhys::RAND(double /*DUMMY*/)
{
  //  This is an adapted version of subroutine RANECU written by F. James
  //  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
  //  give a single random number at each call.
  //
  //  The 'seeds' ISEED1 and ISEED2 must be initialised in the main program
  //  and transferred through the named common block /RSEED/.
  //
  //  Some compilers incorporate an intrinsic random number generator with
  //  the same name (but with different argument lists). To avoid conflict,
  //  it is advisable to declare RAND as an external function in all sub-
  //  programs that call it.

////////////////  using namespace RSEED;
/**********************************************************  
  const double USCALE = 1.0E0/2.147483563E9;
  
  int I1 = RSEED_.ISEED1/53668;
  RSEED_.ISEED1 = 40014*(RSEED_.ISEED1-I1*53668)-I1*12211;
  if(RSEED_.ISEED1 < 0){ RSEED_.ISEED1 = RSEED_.ISEED1+2147483563;}

  int I2 = RSEED_.ISEED2/52774;
  RSEED_.ISEED2 = 40692*(RSEED_.ISEED2-I2*52774)-I2*3791;
  if(RSEED_.ISEED2 < 0){ RSEED_.ISEED2 = RSEED_.ISEED2+2147483399;}

  int IZ = RSEED_.ISEED1-RSEED_.ISEED2;
  if(IZ < 1){ IZ = IZ+2147483562;}
  double RAND_RETURN = IZ*USCALE;

  return RAND_RETURN;
************************************************************/
  return G4UniformRand();
}
