
//  *********************************************************************
//                       SUBROUTINE GPPa
//  *********************************************************************
void PenPhys::GPPa(double &EE, double &CDTE, double &EP, double &CDTP, int &IZZ,int &ISH)
{
  //  Random sampling of electron-positron pair and triplet production by
  //  photons. Bethe-Heitler differential cross section.

  //  Output values:
  //    EE .....  kinetic energy of the electron.
  //    CDTE ...  polar direction cosine of the electron.
  //    EP .....  kinetic energy of the positron.
  //    CDTP ...  polar direction cosine of the positron.
  //    IZZ  ... atomic number of the atom where absorption has occurred.
  //    ISH .... atomic electron shell that has been ionized.

  using namespace PENELOPE_mod;
//////////////  using namespace TRACK_mod;

//////////////  using namespace CEGRID;
  using namespace CECUTR;
  using namespace CGCO;
  using namespace CGPP00;
  using namespace CGPP01;

  const double TREV = 2.0*REV;

  double EKI = REV/E;
  double EPS, ALZ, T, F00, G0, BMIN, G1, G2, G1MIN, G2MIN, A1, P1, XR, RU2M1;
  if(E < 1.1E6)
  {
    EPS = EKI+(1.0-2.0*EKI)*RAND(1.0);
  }
  else
  {
    //  ****  Low-energy and Coulomb corrections.
    ALZ = ZEQPP[MAT-1]/SL;
    T = sqrt(2.0*EKI);
    double TT= T*T;
    F00 = (-1.774-1.210E1*ALZ+1.118E1*ALZ*ALZ)*T+(8.523+7.326E1*ALZ-4.441E1*ALZ*ALZ)*(TT)-(1.352E1+1.211E2*ALZ-9.641E1*ALZ*ALZ)*(TT*T)+(8.946+6.205E1*ALZ-6.341E1*ALZ*ALZ)*(TT*TT);
    G0 = F0[MAT-1][1]+F00;
    BMIN = 4.0*EKI/BCB[MAT-1];
    SCHIFF(BMIN,G1,G2);
    G1MIN = G1+G0;
    G2MIN = G2+G0;
    XR = 0.5-EKI;
    A1 = 6.666666666666666E-1*G1MIN*(XR*XR);
    P1 = A1/(A1+G2MIN);
    //  ****  Random sampling of EPS.
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      if(RAND(2.0) > P1){}
      else
      {
        RU2M1 = 2.0*RAND(3.0)-1.0;
        if(RU2M1 < 0.0)
        {
          EPS = 0.5-XR*pow(fabs(RU2M1),3.333333333333333E-1);
        }
        else
        {
          EPS = 0.5+XR*pow(RU2M1,3.333333333333333E-1);
        }
        double B = EKI/(BCB[MAT-1]*EPS*(1.0-EPS));
        SCHIFF(B,G1,G2);
        G1 += G0;
        if(G1 < 0.0){ G1 = 0.0;}
        if(RAND(4.0)*G1MIN > G1){ brkIt = false; continue;}

        break;
      }
      EPS = EKI+2.0*XR*RAND(5.0);
      double B = EKI/(BCB[MAT-1]*EPS*(1.0-EPS));
      SCHIFF(B,G1,G2);
      G2 += G0;
      if(G2 < 0.0){ G2 = 0.0;}
      if(RAND(6.0)*G2MIN > G2){ brkIt = false; continue;}
    }
  }
  //  ****  Electron.
  EE = EPS*E-REV;
  CDTE = 2.0*RAND(7.0)-1.0;
  A1 = EE+REV;
  double A2 = sqrt(EE*(EE+TREV));
  CDTE = (CDTE*A1+A2)/(A1+CDTE*A2);
  //  ****  Positron.
  EP = (1.0-EPS)*E-REV;
  CDTP = 2.0*RAND(8.0)-1.0;
  A1 = EP+REV;
  A2 = sqrt(EP*(EP+TREV));
  CDTP = (CDTP*A1+A2)/(A1+CDTP*A2);

  //  ****  Triplet production.

  int ISHELL;
  double TRIPL = TRIP[MAT-1][CEGRID_.KE-1]+(TRIP[MAT-1][CEGRID_.KE+1-1]-TRIP[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
  IZZ = 0;
  ISH = 30;
  if(TRIPL < 1.0E-5){ return;}
  if(RAND(9.0) > TRIPL){ return;}
  double TST = RAND(10.0);
  //  ****  Binary search.
  if(TST < PTRSH[MAT-1][0])
  {
    ISHELL = 1;
  }
  else
  {
    ISHELL = 1;
    int JO = NOSCCO[MAT-1]+1;
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      int I = (ISHELL+JO)/2;
      if(TST > PTRSH[MAT-1][I-1])
      {
        ISHELL = I;
      }
      else
      {
        JO = I;
      }
      if(JO-ISHELL > 1){ brkIt = false; continue;}
    }
    ISHELL = ISHELL+1;
  }
  IZZ = KZCO[MAT-1][ISHELL-1];
  ISH = KSCO[MAT-1][ISHELL-1];
}

