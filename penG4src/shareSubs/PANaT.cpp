
//  *********************************************************************
//                       SUBROUTINE PANaT
//  *********************************************************************
void PenInterface::PANaT(double &E, double &TXS)
{
  //  Total cross section (per electron) for annihilation of positrons with
  //  kinetic energy E. Computed from Heitler's dcs formula for annihila-
  //  tion with free electrons at rest.

  //  Output argument:
  //    XST ... total annihilation cross section (cm**2).


  const double PIELR2 = PI*ELRAD*ELRAD;

  double GAM;
  if(E < 1.0){ GAM = 1.0+1.0/REV;}
  else{ GAM = 1.0+E/REV;}
  
  double GAM2 = GAM*GAM;
  double F2 = GAM2-1.0;
  double F1 = sqrt(F2);
  TXS = PIELR2*((GAM2+4.0*GAM+1.0)*log(GAM+F1)/F2-(GAM+3.0)/F1)/(GAM+1.0);
}

