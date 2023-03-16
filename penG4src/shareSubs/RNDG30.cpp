
//  *********************************************************************
//                       FUNCTION RNDG30
//  *********************************************************************
void PenInterface::RNDG30()
{
  //  Initialisation of the RNDG3 sampling function.

  using namespace PENERROR_mod;
  using namespace CRITAA;
  using namespace CRNDG3;

  int N = NR;
  int NU = N/4;
  double XMIN = -3.0;
  double XMAX = +3.0;
  double ERRM;
  int IER;
  RITA0(RNDG3F,XMIN,XMAX,N,NU,ERRM,0,IER);
  if(N != NR || IER != 0){ ErrorFunction(1912); return;}
  if(ERRM > 1.0E-6){ ErrorFunction(1913); return;}

  for(int I = 0; I < NR; I++)
  {
    X_[I] = XA[I];
    A[I] = AA[I];
    B[I] = BA[I];
    F[I] = FA[I];
    KA[I] = IA[I];
  }
  NPM1 = NPM1A;
}

