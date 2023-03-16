
//  *********************************************************************
//                       SUBROUTINE RELAX0
//  *********************************************************************
void PenInterface::RELAX0()
{
  //  This subroutine sets all variables in common /CRELAX/ to zero.

  using namespace CADATA;
  using namespace CRELAX;
  
  for(int I = 0; I < 99; I++)
  {
    NSHT[I] = 0;
    for(int J = 0; J < 16; J++)
    {
      IFIRST[I][J] = 0;
      ILAST[I][J] = 0;
    }
    for(int J = 0; J < 30; J++)
    {
      IFI[I][J] = 0;
      EB[I][J] = 0.0;
      IKS[I][J] = 0;
    }
  }

  for(int I = 0; I < NRX; I++)
  {
    IS0[I] = 0;
    IS1[I] = 0;
    IS2[I] = 0;
    ET[I] = 0.0;
    P[I] = 0.0;
    F[I] = 0.0;
  }
  NCUR = 0;
}

