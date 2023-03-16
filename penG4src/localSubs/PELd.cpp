
//  *********************************************************************
//                       SUBROUTINE PELd
//  *********************************************************************
void PenPhys::PELd(double &RNDC, double &RMU)
{
  //     Simulation of positron hard elastic events. Cross sections from
  //  the ELSEPA numerical database.

  //  Argument value:
  //    RNDC ... cutoff value of the uniform random number
  //             (only hard events are simulated).
  //    RMU .... sampled angular deflection, =(1-CDT)/2.

  using namespace PENELOPE_mod;
//////////////  using namespace TRACK_mod;

//////////////  using namespace CEGRID;
  using namespace CPELDB;

  const unsigned int NP_P=128;
  const unsigned int NPM1=NP_P-1;
  //  ****  Energy grid point.
  double PK = (CEGRID_.XEL-CEGRID_.DLEMP[CEGRID_.KE-1])*CEGRID_.DLFC;
  int JE;
  if(RAND(1.0) < PK)
  {
    JE = CEGRID_.KE+1;
  }
  else
  {
    JE = CEGRID_.KE;
  }
  //  ****  Pointer.
  double RU = RNDC+RAND(2.0)*(1.0-RNDC);
  //  ****  Selection of the interval (binary search in a restricted
  //        interval).
  int ITN = (int)(RU*NPM1)+1;
  int I = ITLP[ITN-1][JE-1][MAT-1];
  int J = ITUP[ITN-1][JE-1][MAT-1];
  if(J-I < 2){}
  else
  {
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      int K=(I+J)/2;
      if(RU > PSP[K-1][JE-1][MAT-1])
      {
        I = K;
      }
      else
      {
        J = K;
      }
      if(J-I > 1){ brkIt = false; continue;}
    }
  }
  //  ****  Sampling from the rational inverse cumulative distribution.
  double PP = PSP[I-1][JE-1][MAT-1];
  double RR = RU-PP;
  if(RR > 1.0E-16)
  {
    double XX = XSP[I-1][JE-1][MAT-1];
    double AA = ASP[I-1][JE-1][MAT-1];
    double BB = BSP[I-1][JE-1][MAT-1];
    double D = PSP[I+1-1][JE-1][MAT-1]-PP;
    RMU = XX+((1.0+AA+BB)*D*RR/(D*D+(AA*D+BB*RR)*RR))*(XSP[I+1-1][JE-1][MAT-1]-XX);
  }
  else
  {
    RMU = XSP[I-1][JE-1][MAT-1];
  }
}

