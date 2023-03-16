
//  *********************************************************************
//                       SUBROUTINE EBRaT
//  *********************************************************************
void PenInterface::EBRaT(double &E, double &WCRM, double &XH0, double &XH1, double &XH2,double &XS1, double &XS2, int M)
{
  //  Integrated cross sections for bremss emission by electrons of energy
  //  E in material M, restricted to energy losses larger than and less
  //  than the cutoff energy WCRM.

  //  Output arguments:
  //    XH0 ... total cross section for hard emission (cm**2).
  //    XH1 ... stopping cross section for hard emission (eV*cm**2).
  //    XH2 ... straggling cross section for hard emission (eV**2*cm**2).
  //    XS1 ... stopping cross section for soft emission (eV*cm**2).
  //    XS2 ... straggling cross section for soft emission (eV**2*cm**2).

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

/////////////////  using namespace CEGRID;
  using namespace CEBR;
  using namespace CEBR01;
  using namespace CEBR02;

  const double TREV = 2.0*REV;

  CEGRID_.XEL = log(E);
  if(CEGRID_.XEL < CEGRID_.DLEMP1){ CEGRID_.XEL = CEGRID_.DLEMP1;}
  
  CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
  CEGRID_.KE = (int)CEGRID_.XE;
  CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
  //  ****  Global x-section factor.
  double FACT = ZBR2[M-1]*( ((E+REV)*(E+REV))/(E*(E+TREV)) )*1.0E-27;

  //  ****  Moments of the scaled bremss x-section.

  double WCRE = WCRM/E;
  for(int IW = 0; IW < NBW; IW++)
  {
    X[IW] = WB[IW];
    Y[IW] = P0[M-1][CEGRID_.KE-1][IW];
  }
  double XH0A = RLMOM(X,Y,X[NBW-1],NBW,-1)-RLMOM(X,Y,WCRE,NBW,-1);
  if(IRETRN != 0){ return;}
  double XS1A = RLMOM(X,Y,WCRE,NBW,0);
  if(IRETRN != 0){ return;}
  double XS2A = RLMOM(X,Y,WCRE,NBW,1);
  if(IRETRN != 0){ return;}
  double XH1A = RLMOM(X,Y,X[NBW-1],NBW,0)-XS1A;
  if(IRETRN != 0){ return;}
  double XH2A = RLMOM(X,Y,X[NBW-1],NBW,1)-XS2A;
  if(IRETRN != 0){ return;}
  for(int IW = 0; IW < NBW; IW++)
  {
    if(CEGRID_.KE+1 < NEGP){ Y[IW] = P0[M-1][CEGRID_.KE+1-1][IW];}
    else{ Y[IW] = P0[M-1][NEGP-1][IW];}
  }
  double XH0B = RLMOM(X,Y,X[NBW-1],NBW,-1)-RLMOM(X,Y,WCRE,NBW,-1);
  if(IRETRN != 0){ return;}
  double XS1B = RLMOM(X,Y,WCRE,NBW,0);
  if(IRETRN != 0){ return;}
  double XS2B = RLMOM(X,Y,WCRE,NBW,1);
  if(IRETRN != 0){ return;}
  double XH1B = RLMOM(X,Y,X[NBW-1],NBW,0)-XS1B;
  if(IRETRN != 0){ return;}
  double XH2B = RLMOM(X,Y,X[NBW-1],NBW,1)-XS2B;
  if(IRETRN != 0){ return;}

  XH0 = ((1.0-CEGRID_.XEK)*XH0A+CEGRID_.XEK*XH0B)*FACT;
  XS1 = ((1.0-CEGRID_.XEK)*XS1A+CEGRID_.XEK*XS1B)*FACT*E;
  XH1 = ((1.0-CEGRID_.XEK)*XH1A+CEGRID_.XEK*XH1B)*FACT*E;
  XS2 = ((1.0-CEGRID_.XEK)*XS2A+CEGRID_.XEK*XS2B)*FACT*E*E;
  XH2 = ((1.0-CEGRID_.XEK)*XH2A+CEGRID_.XEK*XH2B)*FACT*E*E;
}

