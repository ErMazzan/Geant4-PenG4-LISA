
//  *********************************************************************
//                       SUBROUTINE EELa0
//  *********************************************************************
void PenInterface::EELa0(double &XS0, double &XS1, double &XS2, double &XS0H, double &A, double &B, double &RNDC, double &XS1S, double &XS2S)
{
  using namespace PENERROR_mod;

  //     This subroutine determines the parameters of the MW model for
  //  elastic scattering of electrons and positrons and initializes the
  //  mixed simulation algorithm (for particles with a given energy).

  //  Input arguments:
  //    XS0 .... total x-section (cm**2).
  //    XS1 .... 1st transport x-section (cm**2).
  //    XS2 .... 2nd transport x-section (cm**2).
  //    XS0H ... suggested x-section for hard events (cm**2).

  //  Output values:
  //    A, B ... angular distribution parameters.
  //    RNDC ... cutoff probability.
  //    XS0H ... adopted x-section for hard events (cm**2).
  //    XS1S ... 1st transport x-section for soft events (cm**2).
  //    XS2S ... 2nd transport x-section for soft events (cm**2). 


  double RMU1, RMU2;
  if(XS0 < 0.0){ ErrorFunction(1301); return;}
  RMU1 = XS1/(2.0*XS0);
  if(RMU1 > 0.48){ RMU1 = 0.48;} // Ensures numerical consistency.
  
  RMU2 = (3.0*XS1-XS2)/(6.0*XS0);
  if(RMU2 > 0.32){ RMU2 = 0.32;}
  
  if(RMU1 < 0.0 || RMU1 < RMU2)
  {
    printf("\n\n   *** The arguments in subroutine EELa0 are inconsistent.\n       XS0 = %14.7E, XS1 = %14.7E", XS0, XS1);
    ErrorFunction(1302); return;
  }

  //  ****  Wentzel screening parameter.

  double TST;
  double AU, AL;
  A = 1.0;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    A = A+A;
    TST = A*(A+1.0)*log((1.0+A)/A)-A-RMU1;
    if(TST < 0.0){ brkIt = false; continue;}
  }
  AU = A;
  AL = 0.0;
  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    A = 0.5*(AL+AU);
    TST = A*(A+1.0)*log((1.0+A)/A)-A-RMU1;
    if(TST > 0.0)
    {
      AU = A;
    }
    else
    {
      AL = A;
    }
    if(fabs(TST) > 1.0E-15){ brkIt = false; continue;}
  }
  //  ****  At high energies, when truncation errors in the input tables
  //  are significant, we use delta scattering.
  if(RMU2-RMU1*RMU1 < 1.0E-12 || A < 1.0E-9)
  {
    B = 1.0;
    RNDC = 1.0-XS0H/XS0;
    if(RNDC < 1.0E-14){ RNDC = 0.0;}
    XS1S = XS1*RNDC;
    XS2S = XS2*RNDC;
    return;
  }

  double RMU1W, RMU2W;
  RMU1W = A*(A+1.0)*log((1.0+A)/A)-A;
  RMU2W = A*(1.0-2.0*RMU1W);
  B = (RMU2W-RMU2)/(RMU2W-RMU1W*RMU1W);

  //  ****  Case I.
  double RMUC;
  if(B > 0.0)
  {
    RNDC = 1.0-XS0H/XS0;
    if(RNDC < 1.0E-6)
    {
      RNDC = 0.0;
      XS0H = XS0;
      XS1S = 0.0;
      XS2S = 0.0;
      return;
    }

    double A1 = A+1.0;
    double B1 = 1.0-B;
    double RND0 = B1*A1*RMU1/(A+RMU1);
    RNDC = 1.0-XS0H/XS0;
    if(RNDC < RND0)
    {
      RMUC = RNDC*A/(B1*A1-RNDC);
      XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
      XS2S = B1*(A*A1*(RMUC*RMUC)/(A+RMUC))-2.0*A*XS1S;
    }
    else if(RNDC > RND0+B)
    {
      double RNDMB = RNDC-B;
      RMUC = RNDMB*A/(B1*A1-RNDMB);
      XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
      XS2S = B1*(A*A1*(RMUC*RMUC)/(A+RMUC))-2.0*A*XS1S;
      XS1S = XS1S+B*RMU1;
      XS2S = XS2S+B*(RMU1*RMU1);
    }
    else
    {
      RMUC = RMU1;
      double WB = RNDC-RND0;
      XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
      XS2S = B1*(A*A1*(RMUC*RMUC)/(A+RMUC))-2.0*A*XS1S;
      XS1S = XS1S+WB*RMU1;
      XS2S = XS2S+WB*(RMU1*RMU1);
    }
    XS2S = 6.0*XS0*(XS1S-XS2S);
    XS1S = 2.0*XS0*XS1S;
    return;
  }
  if(B > -1.0E-12)
  {
    B = 0.0;
    RNDC = 1.0-XS0H/XS0;
    double A1 = A+1.0;
    RMUC = RNDC*A/(A1-RNDC);
    XS1S = A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
    XS2S = (A*A1*(RMUC*RMUC)/(A+RMUC))-2.0*A*XS1S;
    XS2S = 6.0*XS0*(XS1S-XS2S);
    XS1S = 2.0*XS0*XS1S;
    return;
  }

  //  ****  Case II. 

  double D1, D2;
  
  double C1_Aux = 8.333333333333333E-1;
  double C2_Aux = 7.083333333333333E-1;
  D1 = C2_Aux-RMU2;
  D2 = C1_Aux-RMU1;
  double D3 = C2_Aux*RMU1-C1_Aux*RMU2;
  AL = 1.0E-24;
  AU = A;
  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    A = 0.5*(AL+AU);
    RMU1W = A*(A+1.0)*log((1.0+A)/A)-A;
    RMU2W = A*(1.0-2.0*RMU1W);
    double F = D1*RMU1W-D2*RMU2W-D3;
    if(F < 0.0)
    {
      AL = A;
    }
    else
    {
      AU = A;
    }
    if(AU-AL > 1.0E-14*A){ brkIt = false; continue;}
  }
  B = (RMU1W-RMU1)/(C1_Aux-RMU1W);
  
  RNDC = 1.0-XS0H/XS0;
  if(RNDC < 1.0E-10)
  {
    RNDC = 0.0;
    XS0H = XS0;
    XS1S = 0.0;
    XS2S = 0.0;
    return;
  }
  double A1 = A+1.0;
  double B1 = 1.0+B;
  double RNDCM = B1*A1*0.5/(A+0.5);
  if(RNDC > RNDCM)
  {
    RNDC = RNDCM;
    XS0H = XS0*(1.0-RNDC);
  }
  RMUC = RNDC*A/(B1*A1-RNDC);
  XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
  XS2S = B1*(A*A1*(RMUC*RMUC)/(A+RMUC))-2.0*A*XS1S;
  XS2S = 6.0*XS0*(XS1S-XS2S);
  XS1S = 2.0*XS0*XS1S;
}

