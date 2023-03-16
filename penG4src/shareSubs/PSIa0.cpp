
//  *********************************************************************
//                       SUBROUTINE PSIa0
//  *********************************************************************
void PenInterface::PSIa0()
{
  //  This subroutine sets all variables in common /CPSI0/ to zero.
  //
  //  It has to be invoked before reading the first material definition
  //  file.

  using namespace CPSI0;

  for(int I = 0; I < 99; I++)
  {
    IPSIF[I] = 0;
    IPSIL[I] = 0;
    NSPSI[I] = 0;
  }

  for(int I = 0; I < NRP; I++)
  {
    for(int J = 0; J < 16; J++)
	  {
      XPSI[I][J] = -80.6;
    }
  }
  NCURP = 0;
}

