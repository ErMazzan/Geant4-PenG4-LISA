
//  *********************************************************************
//                       FUNCTION DCSEL
//  *********************************************************************
double DCSEL(double RMU)
{
  //  This function computes the DCS in (cm**2/sr) by linear log-log inter-
  //  polation in RMU=(1-cos(theta))/2. It is initialized by subroutine
  //  DCSEL0, which must be called before using DCSEL.

  using namespace CDCSEP;

  double RMUL;
  if(RMU < 1.0E-35){ RMUL = 1.0E-35;}
  else{ RMUL = RMU;}

  if(RMUL > 0.999999999999){ RMUL = 0.999999999999;}
  
  RMUL = log(RMUL);
  int I;
  FINDI(XMUL,RMUL,NA,I);
  double DCSEL_RETURN = exp(DCSIL[I-1]+(DCSIL[I+1-1]-DCSIL[I-1])*((RMUL-XMUL[I-1])/(XMUL[I+1-1]-XMUL[I-1]))); //Variable que torna la funcio
  //printf("I=%d DCSIL=%12.5E XMUL=%12.5E RMUL=%12.5E DCSEL=%12.5E\n",I,DCSIL[I-1],XMUL[I-1],RMUL,DCSEL_RETURN);
  return DCSEL_RETURN;
}

