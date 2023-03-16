
//  *********************************************************************
//                       SUBROUTINE GPHaT
//  *********************************************************************
void PenInterface::GPHaT(double &E, double &XS, int &M)
{
  //  Delivers the photoelectric cross section XS (in cm**2) for photons of
  //  energy E in material M.

  using namespace PENELOPE_mod;

  using namespace COMPOS;
  using namespace CGPH00;

  double XEL = log(E);
  XS = 0.0;
  double DEE, PCSL;
  int IZZ, I, IU, IT;
  for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
  {
    IZZ = IZ[M-1][IEL];
    I = IPHF[IZZ-1];
    IU = IPHL[IZZ-1];
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      IT = (I+IU)/2;
      if(XEL > EPH[IT-1])
      {
        I = IT;
      }
      else
      {
        IU = IT;
      }
      if(IU-I > 1){ brkIt = false; continue;}
    }
    DEE = EPH[I+1-1]-EPH[I-1];
    if(DEE > 1.0E-15)
    {
      PCSL = XPH[I-1][0]+(XPH[I+1-1][0]-XPH[I-1][0])*(XEL-EPH[I-1])/DEE;
    }
    else
    {
      PCSL = XPH[I-1][0];
    }
    XS = XS+STF[M-1][IEL]*exp(PCSL);
  }
}

