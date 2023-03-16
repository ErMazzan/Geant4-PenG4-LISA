
//  *********************************************************************
//                       FUNCTION GRAaF2
//  *********************************************************************
double GRAaF2( double Q2)
{
  //  Squared molecular form factor, as a function of (Q*SL/REV)**2.

  using namespace PENELOPE_mod;

  using namespace CGRA00;
  using namespace CGRA02;

  double GRAaF2_RETURN; //Substitueix a la variable GRAaF2 en fortran. Es el valor que torna la funcio
  if(Q2 < 1.0E-9)
  {
    GRAaF2_RETURN = FF0[MM-1];
  }
  else if(Q2 > QQM)
  {
    GRAaF2_RETURN = 0.0;
  }
  else
  {
    double QL = log(Q2);
    int I;
    FINDI(QQ,QL,NQ,I);
    double F2 = AR[MM-1][I-1]+QL*(BR[MM-1][I-1]+QL*(CR[MM-1][I-1]+QL*DR[MM-1][I-1]));
    GRAaF2_RETURN = exp(F2);
  }
  return GRAaF2_RETURN;
}

