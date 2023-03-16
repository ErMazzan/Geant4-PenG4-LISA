
//  *********************************************************************
//                       FUNCTION EINaDS
//  *********************************************************************
double EINaDS(double RMU)
{
  //  Angular differential cross section for soft close inelastic colli-
  //  sions of electrons.
  //

  using namespace CEIN01;

  double EINaDS_RETURN;
  double AUX = 2.0*RMU*(1.0-RMU);
  double DENOM = EI*AUX+REV;
  double W = CPS*AUX/DENOM;
  double DWDMU = CPS*REV*(2.0-4.0*RMU)/(DENOM*DENOM);
  EINaDS_RETURN = (1.0+((W/(EE-W))*(W/(EE-W)))-(1.0-AMOL)*(W/(EE-W))+AMOL*((W/EE)*(W/EE)))*DWDMU*pow(RMU,MOM)/(W*W);
  return EINaDS_RETURN;
}

