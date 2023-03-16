
//  *********************************************************************
//                       SUBROUTINE GRAa
//  *********************************************************************
void PenPhys::GRAa(double &E_, double &CDT, int &IEFF, int &M)
{
  //  Random sampling of coherent (Rayleigh) scattering.

  using namespace PENELOPE_mod;

  using namespace COMPOS;
  using namespace CADATA;
////////////  using namespace CEGRID;
  using namespace CGIMFP;
  using namespace CGRA01;
  using namespace CGRA03;
  
  const unsigned int NPM1=NP-1;

  const double RREV = 1.0/REV;

  //  ****  Binary search.

  int II = IED[CEGRID_.KE-1];
  int IU = IEU[CEGRID_.KE-1];
  int IT;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    IT = (II+IU)/2;
    if(CEGRID_.XEL > ERA[IT-1])
    {
      II = IT;
    }
    else
    {
      IU = IT;
    }
    if(IU-II > 1){ brkIt = false; continue;}
  }
  double XSE = exp(XSRA[M-1][II-1]+(XSRA[M-1][II+1-1]-XSRA[M-1][II-1])*(CEGRID_.XEL-ERA[II-1])/(ERA[II+1-1]-ERA[II-1]));
  if(RAND(1.0)*SGRA[M-1][CEGRID_.KE-1] > XSE)
  {
    IEFF = 0;
    CDT = 1.0;
    return;
  }

  IEFF = 1;
  double QMAX = 2.0*E_*RREV;
  if(QMAX < 1.0E-10)
  {
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      CDT = 1.0-2.0*RAND(1.0);
      double G = 0.5*(1.0+CDT*CDT);
      if(RAND(2.0) > G){ brkIt = false; continue;}
    }
    return;
  }
  double Q2MAX = QMAX*QMAX;
  if(Q2MAX > QRA[NP-1][M-1]){ Q2MAX = QRA[NP-1][M-1];}

  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    double RU = RAND(3.0)*PMAX[CEGRID_.KE+1-1][M-1];
  
      //  ****  Selection of the interval
      //        (binary search within pre-calculated limits).

    int ITN = (int)(RU*NPM1+1);
    int I = ITLRA[ITN-1][M-1];
    int J = ITURA[ITN-1][M-1];
    if(J-I < 2){}
    else
    {
      bool brkIt2 = false;
      while(!brkIt2)
      {
        brkIt2 = true;
        int K = (I+J)/2;
        if(RU > PRA[K-1][M-1])
        {
          I = K;
        }
        else
        {
          J = K;
        }
        if(J-I > 1){ brkIt2 = false; continue;}
      }
    
      //  ****  Sampling from the rational inverse cumulative distribution.
    }

    double XX;
    double RR = RU-PRA[I-1][M-1];
    if(RR > 1.0E-16)
    {
      double D = DPRA[I-1][M-1];
      XX = QRA[I-1][M-1]+((1.0+ARA[I-1][M-1]+BRA[I-1][M-1])*D*RR/(D*D+(ARA[I-1][M-1]*D+BRA[I-1][M-1]*RR)*RR))*(QRA[I+1-1][M-1]-QRA[I-1][M-1]);
    }
    else
    {
      XX = QRA[I-1][M-1];
    }
    if(XX > Q2MAX){ brkIt = false; continue;}
    CDT = 1.0-2.0*XX/Q2MAX;
    //  ****  Rejection.
    double G = 0.5*(1.0+CDT*CDT);
    if(RAND(4.0) > G){ brkIt = false; continue;}
  }
}

