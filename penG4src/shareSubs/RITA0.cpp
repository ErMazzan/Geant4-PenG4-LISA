
//  *********************************************************************
//                       SUBROUTINE RITA0
//  *********************************************************************
void PenInterface::RITA0(PDF_RITA PDF, double XLOW, double XHIGH, int &N, int &NU, double &ERRM, int fileidx, int &IER)
{
  //  Initialisation of the RITA algorithm for random sampling of a
  //  continuous random variable X from a probability distribution function
  //  PDF(X) defined in the interval (XLOW,XHIGH). The external function
  //  PDF(X) --not necessarily normalised-- must be provided by the user.
  //  N is the number of points in the sampling grid. These points are
  //  determined by means of an adaptive strategy that minimizes local
  //  interpolation errors. The first NU grid points are uniformly spaced
  //  in (XLOW,XHIGH); when NU is negative, the initial grid consists of
  //  -NU points logarithmically spaced (in this case, XLOW must be
  //  nonnegative).
  //
  //  ERRM is a measure of the interpolation error (the largest value of
  //  the absolute error of the rational interpolation integrated over each
  //  grid interval).
  //
  //  IER is an error code. An output value different from zero indicates
  //  that the initialization procedure has not succeeded; the reason is
  //  written in the standard output unit (6).
  //
  //  ****  Interpolation coefficients and PDF tables are printed on
  //        separate files (UNIT=IWR) if IWR is greater than zero.
  //
  //  Other subprograms needed: EXTERNAL function PDF,
  //                            subroutines RITAI0 and IRND0.
  //

  using namespace CRITA;
  using namespace CRITAA;
  using namespace CRITAN;

  //  ****  Initialisation of the RITA algorithm.
  IER=0;
  RITAI0(PDF,XLOW,XHIGH,N,NU,ERRM,fileidx,IER);
  if(IER != 0){ return;}
  //  ****  Walker's aliasing; cutoff and alias values.
  NPM1A = NPM1;
  for(int I = 0; I < NPM1; I++)
  {
    XA[I] = XT[I];
    AA[I] = A[I];
    BA[I] = B[I];
  }
  XA[NPM1+1-1] = XT[NPM1+1-1];
  double SAVE = CNORM;
  IRND0(DPAC,FA,IA,NPM1A);
  CNORM = SAVE;
  FA[NPM1+1-1] = 1.0;
  IA[NPM1+1-1] = NPM1A;
}

