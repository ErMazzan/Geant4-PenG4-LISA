
//  *********************************************************************
//                  FUNCTION PRANGE
//  *********************************************************************
double PenInterface::PRANGE(double E, int KPAR, int M)
{
//  This function computes the range (in cm) of particles of type KPAR
//  and energy E in material M. For electrons and positrons, the output
//  is the CSDA range. For photons, the delivered value is the mean free
//  path (=inverse attenuation coefficient).

  using namespace PENELOPE_mod;

  using namespace COMPOS;
////////////  using namespace CEGRID;
  using namespace CRANGE;
  using namespace CGIMFP;
  using namespace CGRA01;
  using namespace CGPH00;
  
  double PRANGE_RETURN; 
  double EE;
  if(E < CEGRID_.EL)
  {
    EE = CEGRID_.EL;
  }
  else if(E > CEGRID_.EU)
  {
    EE = CEGRID_.EU;
  }
  else
  {
    EE = E;
  }
  CEGRID_.XEL = log(EE);
  CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
  CEGRID_.KE = (int)CEGRID_.XE;
  if(CEGRID_.KE < 1){ CEGRID_.KE = 1;}
  if(CEGRID_.KE >= NEGP){ CEGRID_.KE = NEGP-1;}
  CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;

  //  ************  Electrons and positrons. (electron KPAR = 0, fotons KPAR = 1, positrons KPAR = 2)

  if(KPAR != 1)
  {
    PRANGE_RETURN = exp(RANGEL[KPAR-1][M-1][CEGRID_.KE-1]+(RANGEL[KPAR-1][M-1][CEGRID_.KE+1-1]-RANGEL[KPAR-1][M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    return PRANGE_RETURN;
  }

  //  ************  Photons.

  //  ****  Rayleigh scattering.
  int II = IED[CEGRID_.KE];
  int IU = IEU[CEGRID_.KE];
  int IT;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    IT = (II+IU)/2;
    if(CEGRID_.XEL > ERA[IT])
    {
      II = IT;
    }
    else
    {
      IU = IT;
    }
    if(IU-II > 1){ brkIt = false; continue;}
  }
  double HMF1 = exp(XSRA[M-1][II-1]+(XSRA[M-1][II+1-1]-XSRA[M-1][II-1])*(CEGRID_.XEL-ERA[II-1])/(ERA[II+1-1]-ERA[II-1]));
  //  ****  Compton scattering.
  double HMF2 = exp(SGCO[M-1][CEGRID_.KE-1]+(SGCO[M-1][CEGRID_.KE+1-1]-SGCO[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
  //  ****  Photoelectric absorption.
  double XS;
  GPHaT(EE,XS,M);
  double HMF3 = XS*VMOL[M-1];
  double HMF4;
  //  ****  Pair production.
  if(E > 1.022E6)
  {
    HMF4 = exp(SGPP[M-1][CEGRID_.KE-1]+(SGPP[M-1][CEGRID_.KE+1-1]-SGPP[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  else
  {
    HMF4 = 0.0;
  }
  double HMF_TOTAL = HMF1+HMF2+HMF3+HMF4;
  if(HMF_TOTAL > 1.0E-35){ PRANGE_RETURN = 1.0/HMF_TOTAL;}
  else{ PRANGE_RETURN = 1.0/1.0E-35;}
      
  return PRANGE_RETURN;
}

