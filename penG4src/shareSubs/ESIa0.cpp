
//  *********************************************************************
//                       SUBROUTINE ESIa0
//  *********************************************************************
void PenInterface::ESIa0()
{
  //  This subroutine sets all variables in common /CESI0/ to zero.
  //
  //  It has to be invoked before reading the first material definition
  //  file.

  using namespace CESI0;

  for(int I = 0; I < 99; I++)
  {
    IESIF[I] = 0;
    IESIL[I] = 0;
    NSESI[I] = 0;
  }

  for(int I = 0; I < NRP; I++)
  {
    for(int J = 0; J < 16; J++)
    {
      XESI[I][J] = -80.6;
    }
  }
  NCURE = 0;
}

