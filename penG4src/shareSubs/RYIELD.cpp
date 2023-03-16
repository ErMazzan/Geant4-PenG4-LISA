
//  *********************************************************************
//                  FUNCTION RYIELD
//  *********************************************************************
double PenInterface::RYIELD(double E, int KPAR, int M)
{
  //  This function calculates the radiative (bremsstrahlung) yields of
  //  electrons and positrons.

  using namespace PENELOPE_mod;

  using namespace COMPOS;
///////////  using namespace CEGRID;
  using namespace CBRYLD;
  
  double RYIELD_RETURN;
  double EE; 
  if(KPAR == 2) 
  {
    RYIELD_RETURN = 0.0;
    return RYIELD_RETURN;
  }

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

  if(KPAR == 1)
  {
    RYIELD_RETURN = exp(EBRY[M-1][CEGRID_.KE-1]+(EBRY[M-1][CEGRID_.KE+1-1]-EBRY[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  else
  {
    RYIELD_RETURN = exp(PBRY[M-1][CEGRID_.KE-1]+(PBRY[M-1][CEGRID_.KE+1-1]-PBRY[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  
  return RYIELD_RETURN;
}

