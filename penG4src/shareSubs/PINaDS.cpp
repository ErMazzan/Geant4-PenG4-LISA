
//  *********************************************************************
//                       FUNCTION PINaDS
//  *********************************************************************
double PINaDS(double RMU)
{
  //  Angular differential cross section for soft close inelastic colli-
  //  sions of positrons.

  using namespace CPIN01;

  double AUX = 2.0*RMU*(1.0-RMU);
  double DENOM = EI*AUX+REV;
  double W = CPS*AUX/DENOM;
  double DWDMU = CPS*REV*(2.0-4.0*RMU)/(DENOM*DENOM);
  double WE = W/EI;
  double PINaDS_RETURN = (1.0-WE*(BHA1-WE*(BHA2-WE*(BHA3-WE*BHA4))))*DWDMU*pow(RMU,MOM)/(W*W);
  return PINaDS_RETURN;
}

