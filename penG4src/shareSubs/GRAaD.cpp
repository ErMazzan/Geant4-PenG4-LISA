
//  *********************************************************************
//                       FUNCTION GRAaD
//  *********************************************************************
double GRAaD(double CDT)
{
  //  Differential x-section for Rayleigh scattering. Form-factor approx.

  using namespace PENELOPE_mod;

  using namespace CGRA00;
  using namespace CGRA02;

  double GRAaD_RETURN;
  double Q2 = FACTE*(1.0-CDT);
  if(Q2 > QQM)
  {
    GRAaD_RETURN = 0.0;
    return GRAaD_RETURN;
  }
  GRAaD_RETURN = (1.0+CDT*CDT)*GRAaF2(Q2);
  if(MOM > 0){ GRAaD_RETURN = GRAaD_RETURN*pow((1.0-CDT)*0.5,MOM);}
  return GRAaD_RETURN;
}

