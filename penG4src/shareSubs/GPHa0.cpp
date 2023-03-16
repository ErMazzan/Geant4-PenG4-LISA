
//  *********************************************************************
//                       SUBROUTINE GPHa0
//  *********************************************************************
void PenInterface::GPHa0()
{
  //  This subroutine sets all variables in common /CGPH00/ to zero.
  //  It has to be invoked before reading the first material definition
  //  file.

  using namespace CADATA;
  using namespace CGPH00;
  
  for(int I = 0; I < 99; I++)
  {
    NPHS[I] = 0;
    IPHF[I] = 0;
    IPHL[I] = 0;
  }

  for(int I = 0; I < NTP; I++)
  {
    EPH[I] = 0.0;
    for(int J = 0; J < 17; J++)
    {
      XPH[I][J] = 1.0E-35;
    }
  }
  NCUR = 0;
}

