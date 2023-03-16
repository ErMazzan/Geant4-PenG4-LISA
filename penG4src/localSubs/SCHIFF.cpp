
//  *********************************************************************
//                       SUBROUTINE SCHIFF
//  *********************************************************************
void PenPhys::SCHIFF(double &B, double &G1, double &G2)
{
  //  Screening functions F1(B) and F2(B) in the Bethe-Heitler differential
  //  cross section for pair production.

  const double TWOPI = PI+PI;
  double B2 = B*B;
  double F1 = 2.0-2.0*log(1.0+B2);
  double F2 = F1-6.666666666666666E-1;
  if(B < 1.0E-10)
  {
    F1 = F1-TWOPI*B;
  }
  else
  {
    double A0 = 4.0*B*atan2(1.0,B);
    F1 = F1-A0;
    F2 = F2+2.0*B2*(4.0-A0-3.0*log((1.0+B2)/B2));
  }
  G1 = 0.5*(3.0*F1-F2);
  G2 = 0.25*(3.0*F1+F2);
}

